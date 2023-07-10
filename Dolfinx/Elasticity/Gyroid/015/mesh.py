import datetime
start_time = datetime.datetime.now()
import gmsh
import os
import sys
dir_path = os.path.dirname(os.path.realpath(__file__))
print(dir_path)
def fltk_options():

    # Type of entity label (0: description,
    #                       1: elementary entity tag,
    #                       2: physical group tag)
    gmsh.option.setNumber("Geometry.LabelType", 2)

    gmsh.option.setNumber("Geometry.PointNumbers", 0)
    gmsh.option.setNumber("Geometry.LineNumbers", 0)
    gmsh.option.setNumber("Geometry.SurfaceNumbers", 2)
    gmsh.option.setNumber("Geometry.VolumeNumbers", 2)

    # Mesh coloring(0: by element type, 1: by elementary entity,
    #                                   2: by physical group,
    #                                   3: by mesh partition)
    gmsh.option.setNumber("Mesh.ColorCarousel", 0)

    gmsh.option.setNumber("Mesh.Lines", 0)
    gmsh.option.setNumber("Mesh.SurfaceEdges", 0)
    gmsh.option.setNumber("Mesh.SurfaceFaces", 0) # CHANGE THIS FLAG TO 0 TO SEE LABELS

    gmsh.option.setNumber("Mesh.VolumeEdges", 2)
    gmsh.option.setNumber("Mesh.VolumeFaces", 2)

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)

gmsh.model.add("Gyroid")
gmsh.option.setString('Geometry.OCCTargetUnit', 'M')

path = os.path.dirname(os.path.abspath(__file__))


gmsh.model.occ.importShapes(os.path.join(path, 'MeshDir/250_lz1.step'))
gmsh.model.occ.synchronize()

import numpy as np

surfaces = gmsh.model.occ.getEntities(dim=2)
print("Total surface: ", len(surfaces))
top, top_mark       = [], 1
bottom, bottom_mark = [], 2

y_start = 2.5*0.001
y_end = 52.5*0.001

for surface in surfaces:
    com = gmsh.model.occ.getCenterOfMass(surface[0], surface[1])
    #print(com)
    if np.isclose(com[1], [y_end]): # TOP TAG1
        print("top worked")
        top.append(surface[1])
        
    elif np.isclose(com[1], [y_start]): #  BOTTOM TAG2
        print("bottom worked")
        bottom.append(surface[1])

gmsh.model.addPhysicalGroup(2, top, top_mark)
gmsh.model.addPhysicalGroup(2, bottom, bottom_mark)

vol_tags=gmsh.model.getEntities(dim=3)
print("Total Number of Volumes: ", len(vol_tags))
gmsh.model.addPhysicalGroup(3, [vol_tags[0][1]], tag=1)

#gmsh.option.setNumber("Mesh.MeshSizeMin", 0.010)
#gmsh.option.setNumber("Mesh.MeshSizeMax", 0.0006)

print("Meshing starting..")
gmsh.model.mesh.generate(3)

gmsh.write("{}.msh".format(os.path.join(path,"MeshDir/Matrix250")))

#if '-nopopup' not in sys.argv:
#    fltk_options()
#    #gmsh.fltk.run()

gmsh.finalize()

print("Total Meshing Time: ", datetime.datetime.now()-start_time)

from xdmf_utils import  write_xdmf_mesh

write_xdmf_mesh(os.path.join(path,"MeshDir/Matrix250"),dimension=3)
