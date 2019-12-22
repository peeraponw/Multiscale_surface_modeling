from abaqus import *
from abaqusConstants import *
from caeModules import *
from regionToolset import *
from mesh import *
from driverUtils import executeOnCaeStartup
import math
import os
import re
import random
import platform
import time
import shutil
import datetime
import part
import material
import section
import assembly
import step
import load
import mesh
import sketch
import interaction
import numpy as np

from odbAccess import *
from odbSection import *

# ------------------------------------------------------------------------------------------------------------

boxsize = 450
surfFile = 'recon2D_Drill2_1_450um_19pnts_trim350-800_waviness.csv'

boxH_Ratio = 'full'            # 'full' for cubic box or number for ratio box
partition_ratio = 0.95         # local mesh: ratio to (minH - boxY_min)
partition_ratio_trans = 0.80    # trans mesh: ration to (minH - boxY_min)

if boxH_Ratio == 'full':
    boxY_min = -boxsize/2
else:
    boxY_min = boxH_Ratio *(boxsize/2)

eps = 1e-5 

# ------------------------------------------------------------------------------------------------------------

path = 'X:/HET_submodel3D/HET_component_model/'

# # 150um component
# odbName = 'HET_2_upgraded'  # without .odb
# elemLocalName = [22738]
# nodesTopLabel = [8373, 8301, 8300, 8372] # order must be ccw around the surface's normal

# # 100um component
# odbName = 'HET_1_upgraded'  # without .odb
# elemLocalName = [80975]
# nodesTopLabel = [19276, 19167, 19166, 19275] # order must be ccw around the surface's normal

# # 50um component
# odbName = 'HET_33_upgraded'  # without .odb
# elemLocalName = [70602]
# nodesTopLabel = [14324, 14301, 14300, 14323] # order must be ccw around the surface's normal

# #30um component
# odbName = 'HET_49_re_upgraded'  # without .odb
# elemLocalName = [164034]
# nodesTopLabel = [28140, 28121, 28120, 28139] # order must be ccw around the surface's normal

# reference 450um submodel to 150um component model
odbName = 'HET_2_upgraded'  # without .odb
elemLocalName = [999]                    # dymmy, not relevant to submodel
nodesTopLabel = [8518, 8302, 8299, 8515] # order must be ccw around the surface's normal


instanceName = 'PLATE'
stepName = 'move'

extPlane = 1 # 0 or 1, define the orientation of the extruded roughness techniques, will no longer used soon

# ------------------------------------------------------------------------------------------------------------

# meshing parameters
surfMeshSize = 8
localMeshSeed = 4
transMeshSeed = 3
globalMeshSize = boxsize/10 

# ------------------------------------------------------------------------------------------------------------

# # #  auxiliary functions

def updateNodeLabelToIdx(label, values):
    '''
    This function receive node label and its field values variable.
    Return idx which refers to corresponding node.
    It also update global array to speed up idx searching process.
    '''
    if label in AllNodesLabel: 
        # if this label has been searched before, the label exists in AllNodesLabel
        # return idx accordingly
        return AllNodesIdx[AllNodesLabel.index(label)]
    if label == values[label-1].nodeLabel: 
        # if label has not yet been searched, assume firstly idx = label-1
        # update AllNodesLabel and AllNodesIdx and return idx = label-1
        AllNodesLabel[label-1] = label
        AllNodesIdx[label-1] = label-1
        return label-1
    for i in range(len(AllNodesIdx)):
        # if idx != label-1, search it from the beginning and find the correct one.
        # update AllNodesLabel and AllNodesIdx and return idx accordingly 
        if AllNodesIdx[i] != -1:
            continue # continue if this idx is searched already
        else:
            AllNodesLabel[i] = values[i].nodeLabel
            AllNodesIdx[i] = i
            if AllNodesLabel[i] == label:
                return i

def updateNodeCoordLabelToIdx(label, instance):
    '''
    This function receive node label and its instance.
    Return idx which refers to corresponding node.
    It also update global array to speed up idx searching process.
    '''
    if label in AllNodesCoordLabel: 
        # if this label has been searched before, the label exists in AllNodesLabel
        # return idx accordingly
        return AllNodesCoordIdx[AllNodesCoordLabel.index(label)]
    if label == instance.nodes[label-1].label: 
        # if label has not yet been searched, assume firstly idx = label-1
        # update AllNodesLabel and AllNodesIdx and return idx = label-1
        AllNodesCoordLabel[label-1] = label
        AllNodesCoordIdx[label-1] = label-1
        return label-1
    for i in range(len(AllNodesCoordIdx)):
        # if idx != label-1, search it from the beginning and find the correct one.
        # update AllNodesLabel and AllNodesIdx and return idx accordingly 
        if AllNodesCoordIdx[i] != -1:
            continue # continue if this idx is searched already
        else:
            AllNodesCoordLabel[i] = instance.nodes[i].label
            AllNodesCoordIdx[i] = i
            if AllNodesCoordLabel[i] == label:
                return i

def updateElemLabelToIdx(label, values):
    '''
    This function receive element label and its field values variable.
    Return idx which refers to corresponding element.
    It also update global array to speed up idx searching process.
    '''
    if label in AllElemsLabel:
        # if this label has been searched before, the label exists in AllNodesLabel
        # return idx accordingly
        return AllElemsIdx[AllElemsLabel.index(label)]
    if label == values[label-1].elementLabel:
        # if label has not yet been searched, assume firstly idx = label-1
        # update AllNodesLabel and AllNodesIdx and return idx = label-1
        AllElemsLabel[label-1] = label
        AllElemsIdx[label-1] = label-1
        return label-1
    for i in range(len(AllElemsIdx)):
        # if idx != label-1, search it from the beginning and find the correct one.
        # update AllNodesLabel and AllNodesIdx and return idx accordingly 
        if AllElemsIdx[i] != -1:
            continue # continue if this idx is searched already
        else:
            AllElemsLabel[i] = values[i].elementLabel
            AllElemsIdx[i] = i
            if AllElemsLabel[i] == label:
                return i

# ------------------------------------------------------------------------------------------------------------
# # # --------------
# # # load odb file
# # # --------------
myOdb = openOdb(path = path+odbName + '.odb')
odbAsm = myOdb.rootAssembly
odbInstance = odbAsm.instances[instanceName]
odbStep = myOdb.steps[stepName]
nElems = len(odbInstance.elements)
nNodes = len(odbInstance.nodes)

# # initialize nodeLabel and nodeIdx array as well as elemLabel and ElemIdx
# The concept is due to inconsistency in node/element label and index.
# Example:
# For model with small amount of elements, one can refer to node 1's reaction force (label=1) by:
## myFrame.fieldOutputs['RF'].values[1-1].data
# However, in case of huge model, the reaction force obtained by
## myFrame.fieldOutputs['RF'].values[10000-1].data
# may not be the reaction force of node 10000.
# The number 1 and 10000 in this example are called "label" which are shown as node number in ABAQUS/CAE,
# while 1-1 and 10000-1 in the brackets are called "index" or "idx". 
# Therefore, the functions "updateNodeLabelToIdx" and "updateElemLabelToIdx" are made 
# to convert label into idx.
global AllNodesLabel, AllNodesIdx, AllElemsLabel, AllElemsIdx, AllNodesCoordIdx
AllNodesLabel = [-1]*nNodes
AllNodesIdx = [-1]*nNodes
AllNodesCoordLabel = [-1]*nNodes
AllNodesCoordIdx = [-1]*nNodes
AllElemsLabel = [-1]*nElems
AllElemsIdx = [-1]*nElems

# # # -----------
# # # get local element
# # # -----------
if type(elemLocalName) is str: # elemLocalName is entered as string
    # obtain elements in the LOCAL set
    nLocalSetsElems = len(odbAsm.elementSets[elemLocalName].elements[0])
    # obtain elements label in the set
    localElemInSetsLabel = [odbAsm.elementSets[elemLocalName].elements[0][i].label for i in range(nLocalSetsElems)]
else: # elemLocaName is entered as list of integer
    # count number of elements
    nLocalSetsElems = len(elemLocalName)
    # element label can be assigned direcly to localElemInSetsLabel
    localElemInSetsLabel = elemLocalName
localElemInSetsIdx = [-1]*len(localElemInSetsLabel) # initilaize with negative values
# convert label to idx
for i, lbl in enumerate(localElemInSetsLabel):
    localElemInSetsIdx[i] = int(updateElemLabelToIdx(lbl, odbStep.frames[0].fieldOutputs['S'].values))

# # # -----------
# # # get local node
# # # -----------
def getNodeOriginCoord(nodesLabel):
    nodesCoord = []
    for i, lbl in enumerate(nodesLabel):
        nodeIdx = updateNodeLabelToIdx(lbl, odbStep.frames[0].fieldOutputs['U'].values)
        nodesCoord.append(odbInstance.nodes[nodeIdx].coordinates)
    return nodesCoord
    
nodesTopCoord = getNodeOriginCoord(nodesTopLabel)

# ------------------------------------------------------------------------------------------------------------

# # # -------------------
# # # create part
# # # -------------------

# Start modelling
executeOnCaeStartup()
myModel = mdb.models['Model-1']

# ------------------------------------------------------------------------------------------------------------

# Import data points
rough_data = np.genfromtxt(surfFile, delimiter=',', skip_header=1)
samplingPnt = int(rough_data.item((-1, 0)))
pntInterval = int(rough_data.item((1, 1))) - int(rough_data.item((0, 1)))
planeNo = int(boxsize/pntInterval)

if(samplingPnt != boxsize):
    print('wrong input file by boxsize')

## Seperate points from each plane to list
pt = []
for i in range(0, boxsize+pntInterval, pntInterval):
    pt.append([[x-boxsize/2, y+boxsize/2] for [z,x,y] in rough_data if z == i])

# find minimum height
h = [y+boxsize/2 for [z,x,y] in rough_data]
minH = np.amin(h)
partitionH = boxY_min + partition_ratio*(minH - boxY_min)
partitionH_trans = boxY_min + partition_ratio_trans*(minH - boxY_min)

# ------------------------------------------------------------------------------------------------------------

# Create box
## Draw first 3D wire as base on first plane
s = myModel.ConstrainedSketch(name = '__profile__', sheetSize = boxsize*2)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.sketchOptions.setValues(decimalPlaces=6)

topSurfPt = pt[0]
leftEnd = topSurfPt[0][1]
rightEnd = topSurfPt[-1][1]

s.Spline(points = topSurfPt)
s.Line(point1 = (-0.5*boxsize, boxY_min), point2 = (0.5*boxsize, boxY_min))
s.Line(point1 = (0.5*boxsize,  boxY_min), point2 = (0.5*boxsize, rightEnd))
s.Line(point1 = (-0.5*boxsize, boxY_min), point2 = (-0.5*boxsize, leftEnd))

myPart = myModel.Part(dimensionality=THREE_D, name='Part-1', type=DEFORMABLE_BODY)
myPart.BaseWire(sketch=s)

del myModel.sketches['__profile__']

# ------------------------------------------------------------------------------------------------------------

## Loop to create sketch on plane

for i in range(pntInterval, boxsize+pntInterval, pntInterval):

    ## Create datum plane
    myPart.DatumPlaneByPrincipalPlane(offset=i, principalPlane=XYPLANE)
    key_xyDatum = myPart.datums.keys()[-1]

    ## Draw sketch
    s = myModel.ConstrainedSketch(name='__profile__', sheetSize=boxsize*2, transform=
        myPart.MakeSketchTransform(sketchPlane=myPart.datums[key_xyDatum], origin=(0,0,i)))
    s.sketchOptions.setValues(decimalPlaces=6)

    topSurfPt = pt[i/pntInterval]

    a = topSurfPt[0]
    leftEnd = a[1]
    b = topSurfPt[-1]
    rightEnd = b[1]

    s.Spline(points = topSurfPt)
    s.Line(point1 = (-0.5*boxsize, boxY_min), point2 = (0.5*boxsize, boxY_min))
    s.Line(point1 = (0.5*boxsize,  boxY_min), point2 = (0.5*boxsize, rightEnd))
    s.Line(point1 = (-0.5*boxsize, boxY_min), point2 = (-0.5*boxsize, leftEnd))

    myPart.Wire(sketch=s, sketchPlane=myPart.datums[key_xyDatum])

    del myModel.sketches['__profile__']

# ------------------------------------------------------------------------------------------------------------

## Solid loft

loftPlane_list = []
for i in range(0, boxsize+pntInterval, pntInterval):
    edge_i = (myPart.edges.getByBoundingBox(xMin=-boxsize, xMax=boxsize, yMin=-boxsize, yMax=boxsize, zMin=-0.5+i, zMax=0.5+i))
    loftPlane_list.append(edge_i)

try:
    (mdb.models['Model-1'].parts['Part-1'].SolidLoft(endCondition=NONE, loftsections=loftPlane_list, startCondition=NONE))
except:
    print('Solid loft error')
    
# ------------------------------------------------------------------------------------------------------------

## Create YZ-datum plane for partition
key_yzDatum = myPart.datums.keys()[-1]
for i in range(pntInterval, boxsize, pntInterval):
    myPart.DatumPlaneByPrincipalPlane(offset=i-boxsize/2, principalPlane=YZPLANE)

## Create XZ-datum plane for local mesh partition
myPart.DatumPlaneByPrincipalPlane(offset=partitionH, principalPlane=XZPLANE)
key_xzDatum = myPart.datums.keys()[-1]

## Create XZ-datum plane for trans mesh partition
myPart.DatumPlaneByPrincipalPlane(offset=partitionH_trans, principalPlane=XZPLANE)
key_xzDatum_trans = myPart.datums.keys()[-1]

## count no. of datum plane
### xy-datum: 2,4,6,...,boxsize*2 = 20 planes (not included most left)
### yz-datum: boxsize*2+2+1, boxsize*2+2+2, ..., boxsize*2+2+(boxsize-1) = 19 planes (not inluded both end sides)
### xz-datum: boxsize*2+2+(boxsize-1)+1 = 1 plane

# key_yzDatum = boxsize*2+2
# key_xzDatum = boxsize*2+2+(boxsize-1)+1

# ------------------------------------------------------------------------------------------------------------

# Partitioning
## XZ-plane partition
myPart.PartitionCellByDatumPlane(cells=myPart.cells, datumPlane=myPart.datums[key_xzDatum])

## Create set for upperBox
upperBox = myPart.Set(name='upperBox', cells=myPart.cells.getByBoundingBox(
    xMin=-boxsize, xMax=boxsize, yMin=partitionH, yMax=boxsize, zMin=-boxsize, zMax=boxsize))

## XY-plane partition
for i in range(1, planeNo, 1):
    myPart.PartitionCellByDatumPlane(cells=upperBox.cells, datumPlane=myPart.datums[2*i])
    print('xy-partition success: ' + str(i) + '/' + str(planeNo-1))

## YZ-plane partition
for i in range(key_yzDatum+3, key_yzDatum+3+planeNo-1, 1):
    myPart.PartitionCellByDatumPlane(cells=upperBox.cells, datumPlane=myPart.datums[i])
    print('yz-partition success: ' + str(i-2-key_yzDatum) + '/' + str(planeNo-1))

## XZ-plane partition for trans mesh
myPart.PartitionCellByDatumPlane(cells=myPart.cells, datumPlane=myPart.datums[key_xzDatum_trans])

# ------------------------------------------------------------------------------------------------------------

# create reference to component model
## not use this variables since boxsize is manaully assigned, not referenced to macro medel
xLocal = np.linalg.norm(nodesTopCoord[extPlane+1]-nodesTopCoord[extPlane])   # mesh width at surface
yLocal = np.linalg.norm(nodesTopCoord[extPlane-1]-nodesTopCoord[extPlane])   # mesh height at surface

# ------------------------------------------------------------------------------------------------------------

# create surface sets
boxFaces_list = []

# back face
face_i = (myPart.faces.getByBoundingBox(
    xMin=-boxsize, xMax=boxsize, 
    yMin=-boxsize, yMax=boxsize,
    zMin=-boxsize, zMax=0))
boxFaces_list.append(face_i)
# right face
face_i = (myPart.faces.getByBoundingBox(
    xMin=boxsize/2., xMax=boxsize, 
    yMin=-boxsize, yMax=boxsize,
    zMin=-boxsize, zMax=boxsize))
boxFaces_list.append(face_i)
# left face
face_i = (myPart.faces.getByBoundingBox(
    xMin=-boxsize, xMax=-boxsize/2., 
    yMin=-boxsize, yMax=boxsize,
    zMin=-boxsize, zMax=boxsize))
boxFaces_list.append(face_i)
# bottom face
face_i = (myPart.faces.getByBoundingBox(
    xMin=-boxsize, xMax=boxsize, 
    yMin=-boxsize, yMax=-boxsize/2.,
    zMin=-boxsize, zMax=boxsize))
boxFaces_list.append(face_i)
# front face
face_i = (myPart.faces.getByBoundingBox(
    xMin=-boxsize, xMax=boxsize, 
    yMin=-boxsize, yMax=boxsize,
    zMin=boxsize, zMax=2*boxsize))
boxFaces_list.append(face_i)

myPart.Set(name='allFaces-pre', faces=boxFaces_list)

# ------------------------------------------------------------------------------------------------------------

# create local csys for transformation in assembly

myPart.DatumCsysByThreePoints(coordSysType=CARTESIAN, name='partCsys', 
    origin=(boxsize/2,boxsize/2,0),
    point1=(-boxsize/2,boxsize/2,0), 
    point2=(boxsize/2,boxsize/2,boxsize)
    )
partCsysIdx = myPart.datums.keys()[-1]
session.viewports['Viewport: 1'].setValues(displayedObject=myPart)

# ------------------------------------------------------------------------------------------------------------

# # Meshing
myPart.deleteMesh(regions=myPart.cells)
myPart.deleteSeeds(regions=myPart.edges)

myPart.setElementType(elemTypes = (ElemType(elemCode=C3D8R, elemLibrary=EXPLICIT, secondOrderAccuracy=OFF, 
    hourglassControl=DEFAULT, distortionControl=DEFAULT),), regions = Region(cells = myPart.cells))
myPart.setMeshControls(elemShape = HEX_DOMINATED, regions = myPart.cells)

myPart.seedPart(size = globalMeshSize)
myPart.seedEdgeByNumber(number = transMeshSeed, edges = myPart.edges.getByBoundingBox(
                        xMin = -boxsize,    xMax = boxsize,
                        yMin = partitionH_trans, yMax = partitionH,
                        zMin = -boxsize,    zMax = boxsize))
myPart.seedEdgeByNumber(number = localMeshSeed, edges = myPart.edges.getByBoundingBox(
                        xMin = -boxsize,    xMax = boxsize,
                        yMin = partitionH, yMax = boxsize,
                        zMin = -boxsize,    zMax = boxsize))
myPart.seedEdgeBySize(size = surfMeshSize, edges = myPart.edges.getByBoundingBox(
                        xMin = -boxsize,    xMax = boxsize,
                        yMin = minH+eps, yMax = boxsize,
                        zMin = -boxsize,    zMax = boxsize),
                        constraint=FINER)                       
myPart.generateMesh(regions = (myPart.cells, ))

# ------------------------------------------------------------------------------------------------------------

# # # -------------------
# # # assembly
# # # -------------------
myAsm = myModel.rootAssembly
myAsm.DatumCsysByDefault(CARTESIAN)
ad = myAsm.datums
myAsm.Instance(dependent=ON, name='subPart', part=myPart)
partAsm = myAsm.instances['subPart']
# myAsm.makeIndependent(instances=(partAsm, ))
pad = partAsm.datums
globalCsys = ad[1]

# # # ------------------
# # # create csys for transformation
# # # ------------------
myAsm.DatumPointByCoordinate(coords=nodesTopCoord[extPlane])
myAsm.DatumPointByCoordinate(coords=nodesTopCoord[extPlane+1])
myAsm.DatumPointByCoordinate(coords=nodesTopCoord[extPlane-1])

myAsm.DatumCsysByThreePoints(coordSysType=CARTESIAN, name='localCsys',
    origin=ad[ad.keys()[-3]], point1=ad[ad.keys()[-2]], point2=ad[ad.keys()[-1]])
localCsys = ad[ad.keys()[-1]]

# # # --------------------
# # # transform the coordinate
# # # --------------------
# minus (boxsize/2., boxsize/2., 0) is to correct origin coordinate of part which is different to P'Warm setup
partAsm.translate(vector=nodesTopCoord[extPlane] - (boxsize/2., boxsize/2., 0))
myAsm.ParallelCsys(fixedCsys=localCsys, movableCsys=partAsm.datums[partCsysIdx])

# # # --------------------
# # # set view point and close .odb
# # # --------------------
session.viewports['Viewport: 1'].setValues(displayedObject=myPart)
myOdb.close()

