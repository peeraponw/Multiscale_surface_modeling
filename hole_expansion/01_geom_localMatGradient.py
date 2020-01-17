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

# # # geometrical parameters
unitFactor = 1000. # 1 for mm-unit, 1000 for um-unit
symmFac = 1 # 4 for full model, 2 for half model, 1 for quarter model, 0 to use revolve angle
revolveAngle = 0 # Revolution angle if symmFac == 0 (in deg)
# plate dimension
d = 10              *unitFactor # hole diameter
w = 110         /2.  *unitFactor # plate width
t = 1.5             *unitFactor # plate thickness
# holder dimension
holderD = 90        *unitFactor # fixing diameter
# punch dimension
punchD = 50         *unitFactor # punch diameter
punchH = 20         *unitFactor # punch length (does not matter)
tipD = 6            *unitFactor
tipH = 37           *unitFactor


# # # meshing parameters
meshSizeLocal = 0.1 *unitFactor
meshSizeTrans = meshSizeLocal * 6
meshSizeGlobal = meshSizeLocal * 6*3
fineMeshRegionFac = 0.05
transMeshRegionFac = 0.2

# # # geometry calculation # # # do not touch this part
tr = w - (2-np.sqrt(2))*w # width after edge trimming
eps = 1e-5 *unitFactor
localMeshD = fineMeshRegionFac*(holderD-d) + d
transMeshD = transMeshRegionFac*(holderD-localMeshD) + localMeshD
revolveAngle= np.deg2rad(revolveAngle)




# # # start modeling
myModel = mdb.models['Model-1']
# # sketch plate
myModel.ConstrainedSketch(name='__profile__', sheetSize=w*2)
mySketch = myModel.sketches['__profile__']
ptProfile = [   (0,  0),
                (0,  w),
                (tr, w),
                (w,  tr),
                (w,  0),
                (0,  0)]
for idx in range(1, len(ptProfile)):
    mySketch.Line(point1 = ptProfile[idx-1], point2 = ptProfile[idx])
myModel.Part(dimensionality=THREE_D, name='plate', type=DEFORMABLE_BODY)
platePart = myModel.parts['plate']
platePart.BaseSolidExtrude(depth=t, sketch=mySketch)
del myModel.sketches['__profile__']

# # mirror and partition part
ptRef = (0,0,0)
if symmFac >= 2:
    platePart.Mirror(keepOriginal=ON, mirrorPlane=platePart.faces.findAt((w/4.,0,t/2.)))
if symmFac >= 4:
    platePart.Mirror(keepOriginal=ON, mirrorPlane=platePart.faces.findAt((0,0,t/2.)))
# holder region
myModel.ConstrainedSketch(name='__profile__', sheetSize=w*2,
                            transform=platePart.MakeSketchTransform(sketchPlane=platePart.faces.findAt((0,0,t)),
                            sketchPlaneSide=SIDE1, sketchUpEdge=platePart.edges.findAt((w,0,0)),
                            sketchOrientation=RIGHT, origin=(0,0,t)))
mySketch = myModel.sketches['__profile__']
mySketch.CircleByCenterPerimeter(center=(0, 0), point1=(0, holderD/2.))
platePart.PartitionFaceBySketch(faces=platePart.faces.getByBoundingBox(
            xMin=-w, xMax=w, yMin=-w, yMax=w, zMin=t, zMax=2*t),
            sketch=mySketch, sketchUpEdge=platePart.edges.findAt((w,0,0)))
platePart.DatumAxisByPrincipalAxis(ZAXIS)
zaxisKey = platePart.datums.keys()[-1]
platePart.PartitionCellByExtrudeEdge(edges=
        platePart.edges.getByBoundingCylinder(center1=(0,0,t), center2=(0,0,t/2.), radius=holderD/2.)[0],
        cells=platePart.cells, line=platePart.datums[zaxisKey], sense=REVERSE)
del myModel.sketches['__profile__']
# auxBC region
myModel.ConstrainedSketch(name='__profile__', sheetSize=w*2,
                            transform=platePart.MakeSketchTransform(sketchPlane=platePart.faces.findAt((0,0,t)),
                            sketchPlaneSide=SIDE1, sketchUpEdge=platePart.edges.findAt((w,0,0)),
                            sketchOrientation=RIGHT, origin=(0,0,t)))
mySketch = myModel.sketches['__profile__']
mySketch.CircleByCenterPerimeter(center=(0, 0), point1=(0, holderD/2.- meshSizeGlobal))
platePart.PartitionFaceBySketch(faces=platePart.faces.getByBoundingBox(
            xMin=-w, xMax=w, yMin=-w, yMax=w, zMin=t, zMax=2*t),
            sketch=mySketch, sketchUpEdge=platePart.edges.findAt((w,0,0)))
platePart.DatumAxisByPrincipalAxis(ZAXIS)
zaxisKey = platePart.datums.keys()[-1]
platePart.PartitionCellByExtrudeEdge(edges=
        platePart.edges.getByBoundingCylinder(center1=(0,0,t), center2=(0,0,t/2.), radius=holderD/2. - 2*eps)[0],
        cells=platePart.cells, line=platePart.datums[zaxisKey], sense=REVERSE)
del myModel.sketches['__profile__']
# transition mesh region
myModel.ConstrainedSketch(name='__profile__', sheetSize=w*2,
                            transform=platePart.MakeSketchTransform(sketchPlane=platePart.faces.findAt((0,0,t)),
                            sketchPlaneSide=SIDE1, sketchUpEdge=platePart.edges.findAt((w,0,0)),
                            sketchOrientation=RIGHT, origin=(0,0,t)))
mySketch = myModel.sketches['__profile__']
mySketch.CircleByCenterPerimeter(center=(0, 0), point1=(0, transMeshD/2.))
platePart.PartitionFaceBySketch(faces=platePart.faces.getByBoundingBox(
            xMin=-w, xMax=w, yMin=-w, yMax=w, zMin=t, zMax=2*t),
            sketch=mySketch, sketchUpEdge=platePart.edges.findAt((w,0,0)))
platePart.PartitionCellByExtrudeEdge(edges=
        platePart.edges.getByBoundingCylinder(center1=(0,0,t), center2=(0,0,t/2.), radius=transMeshD/2.)[0],
        cells=platePart.cells,line=platePart.datums[zaxisKey], sense=REVERSE)
del myModel.sketches['__profile__']
# local mesh region
myModel.ConstrainedSketch(name='__profile__', sheetSize=w*2,
                            transform=platePart.MakeSketchTransform(sketchPlane=platePart.faces.findAt((0,0,t)),
                            sketchPlaneSide=SIDE1, sketchUpEdge=platePart.edges.findAt((w,0,0)),
                            sketchOrientation=RIGHT, origin=(0,0,t)))
mySketch = myModel.sketches['__profile__']
mySketch.CircleByCenterPerimeter(center=(0, 0), point1=(0, localMeshD/2.))
platePart.PartitionFaceBySketch(faces=platePart.faces.getByBoundingBox(
            xMin=-w, xMax=w, yMin=-w, yMax=w, zMin=t, zMax=2*t),
            sketch=mySketch, sketchUpEdge=platePart.edges.findAt((w,0,0)))
platePart.PartitionCellByExtrudeEdge(edges=
        platePart.edges.getByBoundingCylinder(center1=(0,0,t), center2=(0,0,t/2.), radius=localMeshD/2.)[0],
        cells=platePart.cells, line=platePart.datums[zaxisKey], sense=REVERSE)
del myModel.sketches['__profile__']
# local material region 3
myModel.ConstrainedSketch(name='__profile__', sheetSize=w*2,
                            transform=platePart.MakeSketchTransform(sketchPlane=platePart.faces.findAt((0,0,t)),
                            sketchPlaneSide=SIDE1, sketchUpEdge=platePart.edges.findAt((w,0,0)),
                            sketchOrientation=RIGHT, origin=(0,0,t)))
mySketch = myModel.sketches['__profile__']
mySketch.CircleByCenterPerimeter(center=(0, 0), point1=(0, d/2.+ 3*meshSizeLocal))
platePart.PartitionFaceBySketch(faces=platePart.faces.getByBoundingBox(
            xMin=-w, xMax=w, yMin=-w, yMax=w, zMin=t, zMax=2*t),
            sketch=mySketch, sketchUpEdge=platePart.edges.findAt((w,0,0)))
platePart.PartitionCellByExtrudeEdge(edges=
        platePart.edges.getByBoundingCylinder(center1=(0,0,t), center2=(0,0,t/2.), radius=d/2.+ 3*meshSizeLocal)[0],
        cells=platePart.cells, line=platePart.datums[zaxisKey], sense=REVERSE)
del myModel.sketches['__profile__']
# local material region 2
myModel.ConstrainedSketch(name='__profile__', sheetSize=w*2,
                            transform=platePart.MakeSketchTransform(sketchPlane=platePart.faces.findAt((0,0,t)),
                            sketchPlaneSide=SIDE1, sketchUpEdge=platePart.edges.findAt((w,0,0)),
                            sketchOrientation=RIGHT, origin=(0,0,t)))
mySketch = myModel.sketches['__profile__']
mySketch.CircleByCenterPerimeter(center=(0, 0), point1=(0, d/2.+ 2*meshSizeLocal))
platePart.PartitionFaceBySketch(faces=platePart.faces.getByBoundingBox(
            xMin=-w, xMax=w, yMin=-w, yMax=w, zMin=t, zMax=2*t),
            sketch=mySketch, sketchUpEdge=platePart.edges.findAt((w,0,0)))
platePart.PartitionCellByExtrudeEdge(edges=
        platePart.edges.getByBoundingCylinder(center1=(0,0,t), center2=(0,0,t/2.), radius=d/2.+ 2*meshSizeLocal)[0],
        cells=platePart.cells, line=platePart.datums[zaxisKey], sense=REVERSE)
del myModel.sketches['__profile__']
# local material region 1
myModel.ConstrainedSketch(name='__profile__', sheetSize=w*2,
                            transform=platePart.MakeSketchTransform(sketchPlane=platePart.faces.findAt((0,0,t)),
                            sketchPlaneSide=SIDE1, sketchUpEdge=platePart.edges.findAt((w,0,0)),
                            sketchOrientation=RIGHT, origin=(0,0,t)))
mySketch = myModel.sketches['__profile__']
mySketch.CircleByCenterPerimeter(center=(0, 0), point1=(0, d/2.+ meshSizeLocal))
platePart.PartitionFaceBySketch(faces=platePart.faces.getByBoundingBox(
            xMin=-w, xMax=w, yMin=-w, yMax=w, zMin=t, zMax=2*t),
            sketch=mySketch, sketchUpEdge=platePart.edges.findAt((w,0,0)))
platePart.PartitionCellByExtrudeEdge(edges=
        platePart.edges.getByBoundingCylinder(center1=(0,0,t), center2=(0,0,t/2.), radius=d/2.+ meshSizeLocal)[0],
        cells=platePart.cells, line=platePart.datums[zaxisKey], sense=REVERSE)
del myModel.sketches['__profile__']

## for radial setup
if symmFac == 0:
    # angular cut
    myModel.ConstrainedSketch(name='__profile__', sheetSize=w*2,
                            transform=platePart.MakeSketchTransform(sketchPlane=platePart.faces.findAt((0,0,t)),
                            sketchPlaneSide=SIDE1, sketchUpEdge=platePart.edges.findAt((w,t,0)),
                            sketchOrientation=RIGHT, origin=(0,0,t)))
    mySketch = myModel.sketches['__profile__']
    mySketch.Line(point1=(0,0), point2=(2*w, 2*w*np.tan(revolveAngle)))
    mySketch.Line(point1=(2*w,2*w*np.tan(revolveAngle)), point2=(0,2*w))
    mySketch.Line(point1=(0,2*w), point2=(0,0))    
    
    platePart.CutExtrude(sketch=mySketch, sketchOrientation=RIGHT, 
                    sketchPlane=platePart.faces.findAt((0,0,t)),
                    sketchPlaneSide=SIDE1, sketchUpEdge=platePart.edges.findAt((w,t,0)))
    del myModel.sketches['__profile__']
    # platePart.DatumCsysByThreePoints(coordSysType=
        # CARTESIAN, name='Aux-CSYS', 
        # origin=platePart.vertices.findAt((d*np.cos(revolveAngle), d*np.sin(revolveAngle), 0)),
        # point1=platePart.vertices.findAt((d*np.cos(revolveAngle), d*np.sin(revolveAngle), t)),
        # point2=platePart.vertices.findAt((holderD*np.cos(revolveAngle), holderD*np.sin(revolveAngle), 0)))
if symmFac >= 2:
    # horizontal cut
    platePart.DatumPlaneByPrincipalPlane(offset = 0, principalPlane = XZPLANE)
    platePart.PartitionCellByDatumPlane(cells = platePart.cells, datumPlane = platePart.datums[len(platePart.features)])
if symmFac >= 4:
    # vertical cut
    platePart.DatumPlaneByPrincipalPlane(offset = 0, principalPlane = YZPLANE)
    platePart.PartitionCellByDatumPlane(cells = platePart.cells, datumPlane = platePart.datums[len(platePart.features)])

# # make hole
myModel.ConstrainedSketch(name='__profile__', sheetSize=w*2,
                            transform=platePart.MakeSketchTransform(sketchPlane=platePart.faces.findAt((0,0,t)),
                            sketchPlaneSide=SIDE1, sketchUpEdge=platePart.edges.findAt((w,0,0)),
                            sketchOrientation=RIGHT, origin=(0,0,t)))
mySketch = myModel.sketches['__profile__']
mySketch.CircleByCenterPerimeter(center=(0, 0), point1=(0, d/2.))
platePart.CutExtrude(sketch=mySketch, sketchOrientation=RIGHT, 
                    sketchPlane=platePart.faces.findAt((0,0,t)),
                    sketchPlaneSide=SIDE1, sketchUpEdge=platePart.edges.findAt((w,0,0)))
del myModel.sketches['__profile__']

# # local coord for bcs if symmFac == 0
if symmFac == 0:
    platePart.DatumCsysByThreePoints(coordSysType=CARTESIAN, name='Aux-CSYS',         
        origin=platePart.vertices.findAt((d/2.*np.cos(revolveAngle), d/2.*np.sin(revolveAngle), 0)),
        point1=platePart.vertices.findAt((d/2.*np.cos(revolveAngle), d/2.*np.sin(revolveAngle), t)),
        point2=platePart.vertices.findAt((holderD/2.*np.cos(revolveAngle), holderD/2.*np.sin(revolveAngle), 0)))
    csysKey = platePart.datums.keys()[-1]
    

## create set for local material
# global mat
platePart.Set(name='wholePlate', cells=platePart.cells)
platePart.Set(name='localMat1-2-3', cells=platePart.cells.getByBoundingCylinder(
    center1=(0,0,0), center2=(0,0,t), radius=d/2.+3*meshSizeLocal))
platePart.SetByBoolean(name='globalMat', 
    sets=(platePart.sets['wholePlate'], platePart.sets['localMat1-2-3']),
    operation=DIFFERENCE)
# local mat 3
platePart.Set(name='localMat1-2', cells=platePart.cells.getByBoundingCylinder(
    center1=(0,0,0), center2=(0,0,t), radius=d/2.+2*meshSizeLocal))
platePart.SetByBoolean(name='localMat3', 
    sets=(platePart.sets['localMat1-2-3'], platePart.sets['localMat1-2']),
    operation=DIFFERENCE)
# local mat 1 and 2
platePart.Set(name='localMat1', cells=platePart.cells.getByBoundingCylinder(
    center1=(0,0,0), center2=(0,0,t), radius=d/2.+meshSizeLocal))
platePart.SetByBoolean(name='localMat2', 
    sets=(platePart.sets['localMat1-2'], platePart.sets['localMat1']),
    operation=DIFFERENCE)

    
# # sketch punch
myModel.Part(dimensionality=THREE_D, name='punch', type=ANALYTIC_RIGID_SURFACE)
punchPart = myModel.parts['punch']
myModel.ConstrainedSketch(name='__profile__', sheetSize=w*2)
mySketch = myModel.sketches['__profile__']
mySketch.ConstructionLine(point1=(0, 0), point2=(0,1))
ptProfile = [   (0      ,  0),
                (tipD/2.,  0),
                (punchD/2., -tipH) ,
                (punchD/2., -tipH-punchH)]
for idx in range(1, len(ptProfile)):
    mySketch.Line(point1 = ptProfile[idx-1], point2 = ptProfile[idx])
punchPart.AnalyticRigidSurfRevolve(sketch=mySketch)
del myModel.sketches['__profile__']
punchPart.ReferencePoint(point=(0,0,0))

# # # Assembly
myAsm = myModel.rootAssembly
myAsm.DatumCsysByDefault(CARTESIAN)
# import plate part
myAsm.Instance(dependent=OFF, name='plate', part=platePart)
plateAsm = myAsm.instances['plate']
myAsm.rotate(angle=90.0, axisDirection=(1.0, 0.0, 
    0.0), axisPoint=(0.0, 0.0, 0.0), instanceList=('plate', ))
myAsm.translate(instanceList=('plate', ), vector=(
    0.0, t, 0.0))
# # create set for fixed BC
# top
myAsm.Set(name='top-surf', 
    faces=plateAsm.faces.getByBoundingBox(
    xMin=-w, xMax=w, yMin=t, yMax=2*t, zMin=-w, zMax=w))
myAsm.Set(name='top-nonfix', 
    faces=plateAsm.faces.getByBoundingCylinder(
    center1=(0,t,0), center2=(0,t/2.,0), radius=holderD/2.))
myAsm.SetByBoolean(name='top-fix', sets=(myAsm.sets['top-surf'], myAsm.sets['top-nonfix']), 
        operation=DIFFERENCE)
# bot
myAsm.Set(name='bot-surf', 
    faces=plateAsm.faces.getByBoundingBox(
    xMin=-w, xMax=w, yMin=-t, yMax=0, zMin=-w, zMax=w))
myAsm.Set(name='bot-nonfix', 
    faces=plateAsm.faces.getByBoundingCylinder(
    center1=(0,-t,0), center2=(0,0,0), radius=holderD/2.))
myAsm.SetByBoolean(name='bot-fix', sets=(myAsm.sets['bot-surf'], myAsm.sets['bot-nonfix']), 
        operation=DIFFERENCE)
# inner hole
myAsm.Set(name='hole',
    faces=plateAsm.faces.getByBoundingCylinder(
    center1=(0,t,0), center2=(0,0,0), radius=d/2.))
myAsm.SetByBoolean(name='contact', sets=(myAsm.sets['bot-nonfix'], myAsm.sets['hole']),
        operation=UNION)
        
# # create set for symFac == 0 BC
if symmFac == 0:
    asmAngle = np.pi/2. - revolveAngle
    myAsm.Set('auxBC',
        faces=plateAsm.faces.getByBoundingCylinder(
        center1=(d/2.*np.sin(asmAngle),t/2.,d/2.*np.cos(asmAngle)), 
        center2=((holderD/2-meshSizeGlobal)*np.sin(asmAngle),t/2.,(holderD/2-meshSizeGlobal)*np.cos(asmAngle)),
        radius=t/2.))

# import punch part
myAsm.Instance(dependent=OFF, name='punch', part=punchPart)
punchAsm = myAsm.instances['punch']
punchOffset = tipH * (1- (punchD-d)/(punchD-tipD) )
myAsm.translate(instanceList=('punch', ), vector=(
    0.0, punchOffset, 0.0))
myAsm.Set(name='refPunch', referencePoints=(myAsm.instances['punch'].referencePoints[2], ))
# # # Mesh
# global mesh
myAsm.seedPartInstance(size = meshSizeGlobal, regions = (plateAsm, ))
# local mesh
# transition zone
myAsm.seedEdgeBySize(size=meshSizeTrans, edges=plateAsm.edges.getByBoundingCylinder(
    center1=(0,0,0), center2=(0,t+eps,0), radius=transMeshD/2.), constraint=FINER)
myAsm.seedEdgeBySize(size=meshSizeLocal, edges=plateAsm.edges.getByBoundingCylinder(
    center1=(0,0,0), center2=(0,t+eps,0), radius=localMeshD/2.), constraint=FINER)
# mesh control
myAsm.setMeshControls(allowMapped=True, regions=plateAsm.cells, 
    elemShape=HEX_DOMINATED, technique=SWEEP, algorithm=ADVANCING_FRONT)
myAsm.setMeshControls(allowMapped=True, regions=plateAsm.cells.getByBoundingCylinder(
            center1=(0,0,0), center2=(0,t,0), radius=holderD/2.
            ), elemShape=HEX, technique=STRUCTURED)
    
# generate mesh
myAsm.generateMesh(regions = (plateAsm, ))
print('Total nodes: ' + str(len(plateAsm.nodes)))
print('Total elements: ' + str(len(plateAsm.elements)))
