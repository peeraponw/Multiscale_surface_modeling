from abaqus import *
from abaqusConstants import *
from caeModules import *
from regionToolset import *
import numpy as np

# ------------------------------------------------------------------------------------------------------------

# myModel = mdb.models['Model-1']
# myPart = myModel.parts['Part-1']
# myAssembly = myModel.rootAssembly
# cutPart = myAssembly.instances['Cut_Part-1']

meshSize = 10
localMeshSeed = 2
surfMeshSize = 10 

eps = 1e-5

# ------------------------------------------------------------------------------------------------------------

# # Meshing
myAssembly.deleteMesh(regions=cutPart.cells)
myAssembly.deleteSeeds(regions=cutPart.edges)

myAssembly.setElementType(elemTypes = (ElemType(elemCode=C3D8R, elemLibrary=EXPLICIT, secondOrderAccuracy=OFF, 
    hourglassControl=DEFAULT, distortionControl=DEFAULT),), regions = Region(cells = cutPart.cells))
myAssembly.setMeshControls(elemShape = HEX_DOMINATED, regions = cutPart.cells)

myAssembly.seedPartInstance(size = meshSize, regions = (cutPart, ))
myAssembly.seedEdgeByNumber(number = localMeshSeed, edges = cutPart.edges.getByBoundingBox(
                        xMin = -boxsize,    xMax = boxsize,
                        yMin = partitionH, yMax = boxsize,
                        zMin = -boxsize,    zMax = boxsize))
myAssembly.seedEdgeBySize(size = surfMeshSize, edges = cutPart.edges.getByBoundingBox(
                        xMin = -boxsize,    xMax = boxsize,
                        yMin = minH+eps, yMax = boxsize,
                        zMin = -boxsize,    zMax = boxsize),
						constraint=FINER)						
myAssembly.generateMesh(regions = (cutPart, ))

