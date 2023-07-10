import datetime
start_time = datetime.datetime.now()

L = 1
W = 1
mu = 1
rho = 1
delta = W/L
gamma = 0.4*delta**2
beta = 1.25
lambda_ = beta
g = gamma

import numpy as np
import ufl

from mpi4py import MPI
from petsc4py.PETSc import ScalarType
import dolfinx.io
from dolfinx import mesh, fem, plot, io

# mesh reading
from xdmf_utils import XDMFReader
geometry = XDMFReader("MeshDir/Matrix250")
domain, subdomains, facet_tags = geometry.getAll()

t_imap = domain.topology.index_map(domain.topology.dim)
num_cells = t_imap.size_local + t_imap.num_ghosts
print("Total Number of Cells: ", num_cells)

V = fem.VectorFunctionSpace(domain, ("CG", 1))

u_bottom = np.array([0,0,0], dtype=ScalarType)

def impose_bc(mesh, facet_tags, u,tag,V):
    fdim = domain.topology.dim - 1
    boundary_facets = np.array(facet_tags.indices[facet_tags.values == tag])
    bc = fem.dirichletbc(u, fem.locate_dofs_topological(V, fdim, boundary_facets), V)
    return bc

u_top = np.array([0,-25*0.001,0], dtype=ScalarType)

bc1 = impose_bc(mesh, facet_tags, u_top,1,V)
bc2 = impose_bc(mesh, facet_tags, u_bottom,2,V)

bcs = [bc1,bc2]

T = fem.Constant(domain, ScalarType((0, 0, 0)))

ds = ufl.Measure("ds", domain=domain)

def epsilon(u):
    return ufl.sym(ufl.grad(u)) # Equivalent to 0.5*(ufl.nabla_grad(u) + ufl.nabla_grad(u).T)
def sigma(u):
    return lambda_ * ufl.nabla_div(u) * ufl.Identity(u.geometric_dimension()) + 2*mu*epsilon(u)

u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)
f = fem.Constant(domain, ScalarType((0, -rho*g, 0)))
a = ufl.inner(sigma(u), epsilon(v)) * ufl.dx
L = ufl.dot(f, v) * ufl.dx + ufl.dot(T, v) * ds

print("Solver Started..")

problem = fem.petsc.LinearProblem(a, L, bcs=bcs, petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
uh = problem.solve()

with io.XDMFFile(domain.comm, "ResultsDir/deformation.xdmf", "w") as xdmf:
    xdmf.write_mesh(domain)
    uh.name = "Deformation"
    xdmf.write_function(uh)

## Stress computation

s = sigma(uh) -1./3*ufl.tr(sigma(uh))*ufl.Identity(uh.geometric_dimension())
von_Mises = ufl.sqrt(3./2*ufl.inner(s, s))

# %%
V_von_mises = fem.FunctionSpace(domain, ("DG", 0))
stress_expr = fem.Expression(von_Mises, V_von_mises.element.interpolation_points)
stresses = fem.Function(V_von_mises)
stresses.interpolate(stress_expr)

with io.XDMFFile(domain.comm, "ResultsDir/stresses.xdmf", "w") as xdmf:
    xdmf.write_mesh(domain)
    uh.name = "Stresses"
    xdmf.write_function(stresses)

print("Total Execution Time: ", datetime.datetime.now()-start_time)
