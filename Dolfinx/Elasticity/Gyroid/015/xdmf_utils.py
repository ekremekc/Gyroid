from mpi4py import MPI

rank = MPI.COMM_WORLD.rank 

import meshio
import dolfinx.io
from dolfinx.fem import Function, FunctionSpace, locate_dofs_topological
from numpy import array
def create_mesh(mesh, cell_type, prune_z):
    
    cells = mesh.get_cells_type(cell_type)
    cell_data = mesh.get_cell_data("gmsh:physical", cell_type)
    points = mesh.points[:,:2] if prune_z else mesh.points
    out_mesh = meshio.Mesh(points=points, cells={cell_type: cells}, cell_data={"name_to_read":[cell_data]})
    return out_mesh

def write_xdmf_mesh(name, dimension):

    if rank == 0:
        # Read in mesh
        msh_name = name + ".msh"
        msh = meshio.read(msh_name)

        if dimension == 2:
            prune_z = True
            volume_mesh = create_mesh(msh, "triangle",prune_z)
            tag_mesh = create_mesh(msh, "line",prune_z)

        elif dimension == 3:
            prune_z = False
            volume_mesh = create_mesh(msh, "tetra",prune_z)
            tag_mesh = create_mesh(msh, "triangle",prune_z)
            
        # Create and save one file for the mesh, and one file for the facets 
        
        xdmf_name = name + ".xdmf"
        xdmf_tags_name = name + "_tags.xdmf"
        meshio.write(xdmf_name, volume_mesh)
        meshio.write(xdmf_tags_name, tag_mesh)
    print(str(dimension)+"D XDMF mesh is generated.")

def load_xdmf_mesh(name):
    mesh_loader_name = name + ".xdmf"
    tag_loader_name = name + "_tags.xdmf"
    with dolfinx.io.XDMFFile(MPI.COMM_WORLD, mesh_loader_name, "r") as xdmf:
        mesh = xdmf.read_mesh(name="Grid")
        cell_tags = xdmf.read_meshtags(mesh, name="Grid")
    tdim = mesh.topology.dim
    mesh.topology.create_connectivity(tdim-1, tdim)
    with dolfinx.io.XDMFFile(MPI.COMM_WORLD, tag_loader_name, "r") as xdmf:
        facet_tags = xdmf.read_meshtags(mesh, name="Grid")
    if MPI.COMM_WORLD.rank == 0:
        print("XDMF Mesh is loaded.")
    return mesh, cell_tags, facet_tags

class XDMFReader:
    def __init__(self, name):
        self.name = name
        self._mesh = None
        self._cell_tags = None
        self._facet_tags = None
        self._gdim = None
        mesh_loader_name = name + ".xdmf"
        tag_loader_name = name + "_tags.xdmf"
        with dolfinx.io.XDMFFile(MPI.COMM_WORLD, mesh_loader_name, "r") as xdmf:
            self._mesh = xdmf.read_mesh(name="Grid")
            self._cell_tags = xdmf.read_meshtags(self.mesh, name="Grid")
        self.mesh.topology.create_connectivity(self.mesh.topology.dim, self.mesh.topology.dim-1)
        with dolfinx.io.XDMFFile(MPI.COMM_WORLD, tag_loader_name, "r") as xdmf:
            self._facet_tags = xdmf.read_meshtags(self.mesh, name="Grid")
        if MPI.COMM_WORLD.rank == 0:
            print("XDMF Mesh is loaded.")

    @property
    def mesh(self):
        return self._mesh
    @property
    def subdomains(self):
        return self._cell_tags
    @property
    def facet_tags(self):
        return self._facet_tags   
    @property
    def dimension(self):
        return self._mesh.topology.dim    

    def getAll(self):
        return self.mesh, self.subdomains, self.facet_tags
    
    def getNumberofCells(self):
        t_imap = self.mesh.topology.index_map(self.mesh.topology.dim)
        num_cells = t_imap.size_local + t_imap.num_ghosts
        return num_cells
        # print("Number of cells: ", num_cells)


