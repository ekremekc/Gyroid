import numpy as np
import ufl
import pyvista
import matplotlib.pyplot as plt

from petsc4py import PETSc
from mpi4py import MPI
from dolfinx import fem, mesh, plot
from dolfinx import *

#### Defining domain and function spaces
domain = mesh.create_rectangle(MPI.COMM_WORLD,[[0.0,0.0],[1,1]],[10,10],mesh.CellType.triangle)
V = fem.VectorFunctionSpace(domain,("CG",2))    # function space for the displacements and test functions
Vsig = fem.TensorFunctionSpace(domain,("DG",0)) # function space for the stress tensor

# gravity vector
b = -5
K = -0.1 # boundary height

B = fem.Constant(domain,PETSc.ScalarType((0,b)))

#### Getting facets of the bottom edge that will come in contact ####
def bottom(x):
    return np.isclose(x[1],0)

fdim = domain.topology.dim -1 
bottom_facets = mesh.locate_entities_boundary(domain,fdim,bottom)

marked_facets = np.hstack([bottom_facets])
marked_values = np.hstack([np.full(len(bottom_facets),3,dtype=np.int32)])
sorted_facets = np.argsort(marked_facets)
facet_tag = mesh.meshtags(domain,fdim,marked_facets[sorted_facets],marked_values[sorted_facets])


#### Defining elasticity constants ####
E = 100000.0
nu = 0.3
mu = fem.Constant(domain,E / (2.0*(1+nu)))
lmbda = fem.Constant(domain,E*nu / ((1.0+nu)*(1.0-2.0*nu)))

k_pen = 1000000

rho = 1.0
eta_m = 0.01
eta_k = 0.01
#### Defining parameters associated with Generalized alpha method and timestepping ####
am = 0.2
af = 0.4
g = 0.5 + af - am

alpha_m = fem.Constant(domain,am)
alpha_f = fem.Constant(domain,af)
gamma = fem.Constant(domain,g)
beta = fem.Constant(domain,(g+0.5)**2 / 4.0)

T = 3
Nsteps = 10000
dt = fem.Constant(domain,T/Nsteps)

#### Defining functions of interest ####
u_ = ufl.TestFunction(V)
u = fem.Function(V,name='Displacement')

u_old = fem.Function(V)
v_old = fem.Function(V)
a_old = fem.Function(V)


##########################################
###### Independent functions for the #####
###### values related to penetration #####
##########################################
penetrate = fem.Function(V) # value of the penetration at each element
abs_penetrate = fem.Function(V)
penalty_f = fem.Function(V)


X0 = fem.Function(V) # Reference Configuration of the domain
X0.interpolate(lambda x: [x[0],x[1]])
X = fem.Function(V) # Current Configuration of the domain
X.interpolate(X0)

NN = len(X.x.array)
# array indices of the x- and y-components of the vector-valued function arrays
x_inds = np.arange(0,NN,2)
y_inds = np.arange(1,NN+1,2)

#v_old.x.array[x_inds] = 5
#v_old.x.array[y_inds] = 5

#### Defining measures ####
metadata = {"quadrature_degree":4}
ds = ufl.Measure('ds',domain=domain,subdomain_data=facet_tag,metadata=metadata)
dx = ufl.Measure("dx",domain=domain,metadata=metadata)


#### defining elements of GA method ####
def sigma(r):
    return lmbda*ufl.nabla_div(r)*ufl.Identity(len(r)) + 2*mu*ufl.sym(ufl.grad(r))
    
def m(u,u_):
    return rho*ufl.inner(u,u_)*dx

def k(u,u_):
    return ufl.inner(sigma(u),ufl.grad(u_))*dx

def c(u,u_):
    return eta_m*m(u,u_) + eta_k*k(u,u_)


##########################################
##### Independent function to update #####
##### the penalty function, used in  #####
##### in the update_fields function. #####
##########################################

def update_penalty(u,X,X0,penetrate):
    X.x.array[:] = u.x.array + X0.x.array
    penetrate.x.array[y_inds] = K - X.x.array[y_inds] # Boundary at y=K 
    penetrate.x.array[x_inds] = 0
    return 0.5*k_pen*(penetrate.x.array[:] + np.abs(penetrate.x.array[:]))

##########################################
##### Calculates the work associated #####
##### with the penalty as a Form to  #####
##### incorporate with the residual. #####
##########################################

def penalty(penalty_f,u_):
    return ufl.inner(penalty_f,u_)*ds


def update_a(u, u_old, v_old, a_old, ufl_=True):
    if ufl_:
        dt_ = dt
        beta_ = beta
    else:
        dt_ = dt
        beta_ = beta.value
    return (u-u_old-dt_*v_old)/beta_/dt_**2 - (1-2*beta_)/2/beta_*a_old


def update_v(a,u_old,v_old,a_old,ufl_=True):
    if ufl_:
        dt_ = dt
        gamma_ = gamma
    else:
        dt_ = dt
        gamma_ = gamma.value
    return v_old + dt_*((1-gamma_)*a_old + gamma_*a)

def update_fields(u,u_old,v_old,a_old):
    u_vec,u0_vec = u.x.array[:],u_old.x.array[:]
    v0_vec,a0_vec = v_old.x.array[:],a_old.x.array[:]
    
    a_vec = update_a(u_vec,u0_vec,v0_vec,a0_vec,ufl_=False)
    v_vec = update_v(a_vec,u0_vec,v0_vec,a0_vec,ufl_=False)
    v_old.x.array[:] = v_vec
    a_old.x.array[:] = a_vec
    u_old.x.array[:] = u_vec
    ########################################
    ##### Update the penalty function. #####
    ########################################
    penalty_f.x.array[:] = update_penalty(u,X,X0,penetrate)
    
def avg(x_old,x_new,alpha):
    return alpha*x_old + (1-alpha)*x_new

a_new = update_a(u,u_old,v_old,a_old,ufl_=True)
v_new = update_v(a_new,u_old,v_old,a_old,ufl_=True)

##########################################
##### Initialize the penalty function ####
##########################################
penalty_f.x.array[:] = update_penalty(u,X,X0,penetrate)

#### Defining the problem and solver ####
res =  m(avg(a_old,a_new,alpha_m),u_) + c(avg(v_old,v_new,alpha_f),u_) + k(avg(u_old,u,alpha_f),u_) - rho*ufl.inner(B,u_)*dx - penalty(penalty_f,u_) #residual to solve
problem = fem.petsc.NonlinearProblem(res,u)#,bcs=bcs)

from dolfinx import nls
solver = nls.petsc.NewtonSolver(domain.comm, problem)

solver.atol = 1e-8
solver.rtol = 1e-8
solver.convergence_criterion = "incremental"

#### Quatities to track ####
time = np.linspace(0, T, Nsteps+1)

V0 = 1

# center of mass positions
Mx_ = np.zeros((Nsteps+1,))
My_ = np.zeros((Nsteps+1,))
Mx = fem.form((1/V0)*(u[0]+X0[0])*dx)
My = fem.form((1/V0)*(u[1]+X0[1])*dx)
Mx_[0] = fem.assemble_scalar(Mx)
My_[0] = fem.assemble_scalar(My)
max_pen = np.zeros((Nsteps+1,))

from dolfinx import log
#log.set_log_level(log.LogLevel.INFO)
#### Iterating through time ####
for (i, dt) in enumerate(np.diff(time)):
    t = time[i+1]   
    num_its,converged = solver.solve(u) # solve the current time step
    assert(converged)
    u.x.scatter_forward()
    
    # Update old fields with new quantities
    update_fields(u, u_old, v_old, a_old)
        
    Mx_[i+1] = fem.assemble_scalar(Mx)
    My_[i+1] = fem.assemble_scalar(My)
    max_pen[i+1] = max(penalty_f.x.array[y_inds])/k_pen ## Maximum penetration of the cube

# Plot center of mass evolution
plt.figure()
plt.plot(time,Mx_,time,My_)
plt.legend(("Mx","My"))
plt.xlabel("Time")
plt.ylabel("Center of Mass")
plt.show()   

# Plot of the maximum penetration
plt.figure()
plt.plot(time,max_pen)
plt.legend(("Max Penetration"))
plt.xlabel("Time")
plt.ylabel("Penetration")
plt.show()
