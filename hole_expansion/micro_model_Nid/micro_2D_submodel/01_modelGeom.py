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
from odbAccess import *
from odbSection import *

# ------------------------------------------------------------------------------------------------------------

unitFactor = 1000. # 1 for mm-unit, 1000 for um-unit

boxsize = 30

dimension = '3D'

surfFile = 'recon_1D_flat_elemSize30.csv'

# ------------------------------------------------------------------------------------------------------------

path = 'X:/HET_submodel/HET_component_model/'

# data of the macro model
odbName = 'HET_49_re_upgraded'                  # without .odb
elemLocalName = [999]                           # dummy (can put any number)
nodesTopLabel = [28599, 28409, 28395, 28585]    # order must be ccw around the surface's normal

# ------------------------------------------------------------------------------------------------------------

instanceName = 'PLATE'
stepName = 'move'

extPlane = 1 # 0 or 1, define the orientation of the extruded roughness techniques, will no longer used soon

# ------------------------------------------------------------------------------------------------------------

# # # meshing parameters
meshSizeLocal = 1
meshSizeTrans = meshSizeLocal * 2
meshSizeGlobal = meshSizeLocal * 3

localPartitionH_ratio = 0.8
transPartitionH_ratio = 0.5

eps = 1e-5 

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
# # # get local element
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

# initialize data point
rough_data = np.genfromtxt(surfFile, delimiter=',', skip_header=1)
samplingPt = len(rough_data)

a = np.zeros((samplingPt, 2))
for i in range(0, samplingPt):
    a[(i,0)] = -boxsize/2
    a[(i,1)] = boxsize/2
topSurfPt = rough_data + a

leftEnd = topSurfPt[(0,1)]
rightEnd = topSurfPt[(-1,1)]

# ------------------------------------------------------------------------------------------------------------

# generate part
p = myModel.Part(dimensionality=THREE_D, name='subPart', type=DEFORMABLE_BODY)
f, e, d = p.faces, p.edges, p.datums

# Create roughness box
s = myModel.ConstrainedSketch(name = '__profile__', sheetSize = 200.0,)
s.sketchOptions.setValues(decimalPlaces=6)

# sketch part
s.Spline(points = topSurfPt)
s.Line(point1 = (-0.5*boxsize, -0.5*boxsize), point2 = (0.5*boxsize, -0.5*boxsize))
s.Line(point1 = (-0.5*boxsize, -0.5*boxsize), point2 = (-0.5*boxsize, leftEnd))
s.Line(point1 = (0.5*boxsize, -0.5*boxsize), point2 = (0.5*boxsize, rightEnd))

p.BaseSolidExtrude(sketch=s, depth=boxsize)

del myModel.sketches['__profile__']

# ------------------------------------------------------------------------------------------------------------

# create reference to component model
## not use this variables since boxsize is manaully assigned, not referenced to macro medel
xLocal = np.linalg.norm(nodesTopCoord[extPlane+1]-nodesTopCoord[extPlane])   # mesh width at surface
yLocal = np.linalg.norm(nodesTopCoord[extPlane-1]-nodesTopCoord[extPlane])   # mesh height at surface

# ------------------------------------------------------------------------------------------------------------

# # # -----------------------
# # # create surface sets
# # # -----------------------

boxFaces_list = []

# back face
face_i = (p.faces.getByBoundingBox(
    xMin=-boxsize, xMax=boxsize, 
    yMin=-boxsize, yMax=boxsize,
    zMin=-boxsize, zMax=0))
boxFaces_list.append(face_i)
# right face
face_i = (p.faces.getByBoundingBox(
    xMin=boxsize/2., xMax=boxsize, 
    yMin=-boxsize, yMax=boxsize,
    zMin=-boxsize, zMax=boxsize))
boxFaces_list.append(face_i)
# left face
face_i = (p.faces.getByBoundingBox(
    xMin=-boxsize, xMax=-boxsize/2., 
    yMin=-boxsize, yMax=boxsize,
    zMin=-boxsize, zMax=boxsize))
boxFaces_list.append(face_i)
# bottom face
face_i = (p.faces.getByBoundingBox(
    xMin=-boxsize, xMax=boxsize, 
    yMin=-boxsize, yMax=-boxsize/2.,
    zMin=-boxsize, zMax=boxsize))
boxFaces_list.append(face_i)
# front face
face_i = (p.faces.getByBoundingBox(
    xMin=-boxsize, xMax=boxsize, 
    yMin=-boxsize, yMax=boxsize,
    zMin=boxsize, zMax=2*boxsize))
boxFaces_list.append(face_i)

p.Set(name='allFaces-pre', faces=boxFaces_list)

# # # -------------------------------------------------
# # # create local csys for transformation in assembly
# # # -------------------------------------------------

p.DatumCsysByThreePoints(coordSysType=CARTESIAN, name='partCsys', 
    origin=(boxsize/2,boxsize/2,0),
    point1=(-boxsize/2,boxsize/2,0), 
    point2=(boxsize/2,boxsize/2,boxsize)
    )
partCsysIdx = d.keys()[-1]
session.viewports['Viewport: 1'].setValues(displayedObject=p)

# # # --------------------------
# # # partition part for meshing
# # # --------------------------

a = [y+boxsize/2. for [x,y] in rough_data]
minH = np.amin(a)

localPartitionH = localPartitionH_ratio * minH
transPartitionH = transPartitionH_ratio * minH

p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=localPartitionH)
p.PartitionCellByDatumPlane(cells=p.cells, datumPlane=d[d.keys()[-1]])
p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=transPartitionH)
p.PartitionCellByDatumPlane(cells=p.cells, datumPlane=d[d.keys()[-1]])

# # # -------------------
# # # assembly
# # # -------------------
myAsm = myModel.rootAssembly
myAsm.DatumCsysByDefault(CARTESIAN)
ad = myAsm.datums
myAsm.Instance(dependent=ON, name='subPart', part=p)
partAsm = myAsm.instances['subPart']
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

# # # ------------------
# # # Mesh
# # # ------------------
# global mesh
p.seedPart(size = meshSizeGlobal)
# transition zone
p.seedEdgeBySize(size=meshSizeTrans, edges=p.edges.getByBoundingBox(
    xMin=-boxsize, xMax=boxsize, 
    yMin=-transPartitionH, yMax=boxsize,
    zMin=-boxsize, zMax=boxsize),
constraint=FINER)
# local zone
p.seedEdgeBySize(size=meshSizeLocal, edges=p.edges.getByBoundingBox(
    xMin=-boxsize, xMax=boxsize, 
    yMin=-localPartitionH, yMax=boxsize,
    zMin=-boxsize, zMax=boxsize), 
constraint=FINER)
# generate mesh
p.generateMesh(regions=(p.cells,))
session.viewports['Viewport: 1'].setValues(displayedObject=p)
myOdb.close()
