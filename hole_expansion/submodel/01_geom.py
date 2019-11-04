# Import all auxilliary modules
from abaqus import *
from abaqusConstants import *
from part import *
from material import *
from section import *
from assembly import *
from step import *
from regionToolset import *
from interaction import *
from load import *
from mesh import *
import numpy as np

unitFactor = 1000. # 1 for mm-unit, 1000 for um-unit

# # # odb parameters
'''
odbName = 'het1tipH37'
elemLocalName = [10228]
instanceName = 'PLATE'
stepName = 'move'
nodesTopLabel = [4916, 4917, 4863, 4862] # order must be ccw around the surface's normal
'''
'''
odbName = 'het1tipH37_f'
elemLocalName = [81019]
instanceName = 'PLATE'
stepName = 'move'
nodesTopLabel = [19319, 19320, 19211, 19210] # order must be ccw around the surface's normal
'''
odbName = 'het1tipH37_f'
elemLocalName = [81019]
instanceName = 'PLATE'
stepName = 'move'
nodesTopLabel = [19208, 19206, 19424, 19426] # order must be ccw around the surface's normal
# # # sub-model parameters
modelHeightRatio = 0.5 # extrusion depth ratio, proportional to diagonal distance between surface nodes
modelWidthRatio = 1.0
extPlane = 1 # 0 or 1, define the orientation of the extruded roughness techniques, will no longer used soon

# # # geometrical parameters
roughnessMean = 0   # in micrometer
roughnessSD = 4.9438 # dimensionless # D=2.8644(rs50), M=2.2196(rs77), WJ=4.9438(rs11), WC=1.5469(rs24)
rdmseed = 11 # seed for random generator
randomPtNum = 49 # number of points in the middle to random depth
order = 2      # spline equation's order (cannot be bigger than randomPtNum+1)
if order > randomPtNum+1: order = randomPtNum+1

# # # meshing parameters
meshSizeLocal = 1.0
meshSizeTrans = meshSizeLocal * 2
meshSizeGlobal = meshSizeLocal * 3

# # #  auxiliary functions
eps = 1e-5 
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


# # # --------------
# # # load odb file
# # # --------------
myOdb = openOdb(path = odbName + '.odb')
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

# # # -------------------
# # # create part
# # # -------------------
myModel = mdb.models['Model-1']
myModel.Part(name='subPart', dimensionality=THREE_D, type=DEFORMABLE_BODY)
p = myModel.parts['subPart']
f, e, d = p.faces, p.edges, p.datums
p.DatumAxisByTwoPoint(point1 = (0,0,0), point2 = (1,1,1))
p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=0)
# for coord in nodesTopCoord:
    # p.DatumPointByCoordinate(coords = coord)
xLocal = np.linalg.norm(nodesTopCoord[extPlane+1]-nodesTopCoord[extPlane])
yLocal = np.linalg.norm(nodesTopCoord[extPlane-1]-nodesTopCoord[extPlane])
modelHeight = modelHeightRatio*np.linalg.norm(nodesTopCoord[extPlane]- nodesTopCoord[extPlane+2])

s = myModel.ConstrainedSketch(name='__profile__', sheetSize=2*modelHeight)
s.rectangle(point1=(0,0), point2=(xLocal, modelWidthRatio*yLocal))
p.BaseSolidExtrude(sketch=s, depth=modelHeight)
p.SolidExtrude(sketchPlane=d[2], sketchUpEdge=p.edges.findAt((0,modelWidthRatio*yLocal/2.,0)),
    sketchPlaneSide=SIDE1, sketch=s, sketchOrientation=LEFT,
    depth=modelHeight, flipExtrudeDirection=ON)

p.DatumPointByCoordinate((0,0,0))
p.DatumPointByCoordinate((xLocal,0,0))
p.DatumPointByCoordinate((0,yLocal,0))

# # # -----------------------
# # # create surface sets
# # # -----------------------
p.Set(name='allFaces-pre', faces=p.faces)

# # # -----------------------
# # # sketch roughness region
# # # -----------------------
t = p.MakeSketchTransform(sketchPlane=p.faces.findAt((xLocal/2.,0,0)),
    sketchUpEdge=p.edges.findAt((xLocal,0,0)),
    sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0,0,0))
boxerr = 1
while boxerr:
    
    s = myModel.ConstrainedSketch(name='__profile__', sheetSize=xLocal, transform=t)
    if roughnessMean == 0 and roughnessSD == 0:
        s.rectangle(point1=(0,0), point2=(xLocal, modelHeight*3))
    else:
        np.random.seed(rdmseed)
        # sampling point in local x-axis
        samplingPt = np.column_stack(np.arange(xLocal/(randomPtNum+1), xLocal, xLocal/(randomPtNum+1)))
        samplingPt = np.row_stack((0, samplingPt.T, xLocal))
        # height in local y-axis
        height = np.random.normal(roughnessMean, roughnessSD, randomPtNum)
        height = np.reshape(height, (randomPtNum, 1))
        height = np.row_stack((0, height, 0))
        # concatenate `samplingPt` and `height`
        topSurfPt = np.column_stack((samplingPt, height)).tolist()
        # sketch rough surface
        for i in range(0, len(topSurfPt)-order+1, order):
            s.Spline(points = topSurfPt[i:i+order+1])
        # sketch other lines
        s.Line(point1=(0,0), point2=(-2*xLocal,0))
        s.Line(point1=(-2*xLocal,0), point2=(-2*xLocal,2*xLocal))
        s.Line(point1=(-2*xLocal,2*xLocal), point2=(2*xLocal,2*xLocal))
        s.Line(point1=(2*xLocal,2*xLocal), point2=(2*xLocal,0))
        s.Line(point1=(2*xLocal,0), point2=(xLocal,0))
        
    try:
        p.CutExtrude(sketch=s, depth=2*xLocal, 
            sketchPlane=p.faces.findAt((xLocal/2.,0,0)), sketchPlaneSide=SIDE1,
            sketchUpEdge = p.edges.findAt((xLocal,0,0)), sketchOrientation=RIGHT,
            flipExtrudeDirection=OFF)
        boxerr=0
        print("Create box successfully")
    except:
        boxerr=1
        rdmseed += 1
        print("Resketch the box with random seed = ", rdmseed)
# # # -------------------------------------------------
# # # create local csys for transformation in assembly
# # # -------------------------------------------------
p.DatumCsysByThreePoints(coordSysType=CARTESIAN, name='partCsys', origin=(0,0,0),
    point1=p.vertices.findAt((xLocal,0,0)), point2=p.vertices.findAt((0,modelWidthRatio*yLocal,0)))
partCsysIdx = d.keys()[-1]
session.viewports['Viewport: 1'].setValues(displayedObject=p)
# # # --------------------------
# # # partition part for meshing
# # # --------------------------
internalDepth = modelHeight+np.min(height)
p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=np.min(height)-0.5*internalDepth)
p.PartitionCellByDatumPlane(cells=p.cells.findAt((0,0,0)), datumPlane=d[d.keys()[-1]])
p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=np.min(height)-0.1*internalDepth)
p.PartitionCellByDatumPlane(cells=p.cells.findAt((0,0,0)), datumPlane=d[d.keys()[-1]])

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
myAsm.ParallelCsys(fixedCsys=localCsys, movableCsys=partAsm.datums[partCsysIdx])
partAsm.translate(vector=nodesTopCoord[extPlane])

# # # ------
# # # Mesh
# # # ------
# global mesh
p.seedPart(size = meshSizeGlobal)
# local mesh
# transition zone
p.seedEdgeBySize(size=meshSizeTrans, edges=p.edges.getByBoundingBox(
    xMin=-10*modelHeight, xMax=10*modelHeight, 
    yMin=-10*modelHeight, yMax=10*modelHeight,
    zMin=np.min(height)-0.5*internalDepth, zMax=10*modelHeight), constraint=FINER)
# local zone
p.seedEdgeBySize(size=meshSizeLocal, edges=p.edges.getByBoundingBox(
    xMin=-10*modelHeight, xMax=10*modelHeight, 
    yMin=-10*modelHeight, yMax=10*modelHeight,
    zMin=np.min(height)-0.1*internalDepth, zMax=10*modelHeight), constraint=FINER)
# generate mesh
p.generateMesh(regions=(p.cells,))
session.viewports['Viewport: 1'].setValues(displayedObject=p)
myOdb.close()
