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
unitFactor = 1000 # 1 for mm-unit, 1000 for um-unit
# plate dimension
w = 20            /2. *unitFactor # width
l = 60                *unitFactor # length
t = 1.47               *unitFactor  # thickness
# support dimension
R = 8              *unitFactor     # support radius
supportDist = 24      *unitFactor  
# punch dimension
#rtratio = 7          # r/t ratio
r = 0.5             *unitFactor# punch radius
h = 20              *unitFactor# punch height

# # # meshing parameters
meshSizeLocal = 0.05 *unitFactor *2
meshSizeTrans = meshSizeLocal * 3
meshSizeGlobal = meshSizeLocal * 6

# # # # Start modeling
myModel = mdb.models['Model-1']
# # sketch plate
myModel.ConstrainedSketch(name = '__profile__', sheetSize = l*2)
mySketch = myModel.sketches['__profile__']
ptProfile = [ (-l/2., 0),
              (-l/2., t),
              ( l/2., t),
              ( l/2., 0),
              (-l/2., 0)] # draw a rectangle
for idx in range(1, len(ptProfile)):
    mySketch.Line(point1 = ptProfile[idx-1], point2 = ptProfile[idx])
myModel.Part(dimensionality=THREE_D, name='plate', type=DEFORMABLE_BODY)
platePart = myModel.parts['plate']
platePart.BaseSolidExtrude(depth=w, sketch=mySketch)
del myModel.sketches['__profile__']
# # partition part
ptRef = (0, 0, 0)
# vertical cut at left support
platePart.DatumPlaneByPrincipalPlane(offset = -(supportDist+R)/2., principalPlane = YZPLANE)
platePart.PartitionCellByDatumPlane(cells = platePart.cells.findAt(ptRef), datumPlane = platePart.datums[len(platePart.features)])
platePart.DatumPlaneByPrincipalPlane(offset = -(supportDist)/2., principalPlane = YZPLANE)
platePart.PartitionCellByDatumPlane(cells = platePart.cells.findAt(ptRef), datumPlane = platePart.datums[len(platePart.features)])
platePart.DatumPlaneByPrincipalPlane(offset = -(supportDist-R)/2., principalPlane = YZPLANE)
platePart.PartitionCellByDatumPlane(cells = platePart.cells.findAt(ptRef), datumPlane = platePart.datums[len(platePart.features)])
# vertical cut at right support
platePart.DatumPlaneByPrincipalPlane(offset = (supportDist+R)/2., principalPlane = YZPLANE)
platePart.PartitionCellByDatumPlane(cells = platePart.cells.findAt(ptRef), datumPlane = platePart.datums[len(platePart.features)])
platePart.DatumPlaneByPrincipalPlane(offset = (supportDist)/2., principalPlane = YZPLANE)
platePart.PartitionCellByDatumPlane(cells = platePart.cells.findAt(ptRef), datumPlane = platePart.datums[len(platePart.features)])
platePart.DatumPlaneByPrincipalPlane(offset = (supportDist-R)/2., principalPlane = YZPLANE)
platePart.PartitionCellByDatumPlane(cells = platePart.cells.findAt(ptRef), datumPlane = platePart.datums[len(platePart.features)])
# vertical cut at the middle
platePart.DatumPlaneByPrincipalPlane(offset = -2*r, principalPlane = YZPLANE)
platePart.PartitionCellByDatumPlane(cells = platePart.cells.findAt(ptRef), datumPlane = platePart.datums[len(platePart.features)])
platePart.DatumPlaneByPrincipalPlane(offset = 2*r, principalPlane = YZPLANE)
platePart.PartitionCellByDatumPlane(cells = platePart.cells.findAt(ptRef), datumPlane = platePart.datums[len(platePart.features)])
platePart.DatumPlaneByPrincipalPlane(offset = 0, principalPlane = YZPLANE)
platePart.PartitionCellByDatumPlane(cells = platePart.cells.findAt(ptRef), datumPlane = platePart.datums[len(platePart.features)])
# horizontal cut at left, right and middle
ptRefLeft = (-supportDist, 0, 0)
ptRefRight = (supportDist, 0, 0)
platePart.DatumPlaneByPrincipalPlane(offset = 0.2*t, principalPlane = XZPLANE)
horPlaneSup = platePart.datums[len(platePart.features)]
platePart.PartitionCellByDatumPlane(cells = platePart.cells.getByBoundingBox(
    xMin= -(supportDist+R)/2., xMax= -(supportDist-R)/2., yMin=0, yMax=t, zMin=0, zMax=w), datumPlane = horPlaneSup)
platePart.PartitionCellByDatumPlane(cells = platePart.cells.getByBoundingBox(
    xMin= (supportDist-R)/2., xMax= (supportDist+R)/2., yMin=0, yMax=t, zMin=0, zMax=w), datumPlane = horPlaneSup)
platePart.PartitionCellByDatumPlane(cells = platePart.cells.getByBoundingBox(
    xMin= -2*r, xMax=  2*r, yMin=0, yMax=t, zMin=0, zMax=w), datumPlane = horPlaneSup)
platePart.DatumPlaneByPrincipalPlane(offset = 0.8*t, principalPlane = XZPLANE)
horPlanePunch = platePart.datums[len(platePart.features)]
platePart.PartitionCellByDatumPlane(cells = platePart.cells.getByBoundingBox(
    xMin= -2*r, xMax= 2*r, yMin=0, yMax=t, zMin=0, zMax=w), datumPlane = horPlanePunch)

# # sketch support
myModel.ConstrainedSketch(name='__profile__', sheetSize = w*2)
mySketch = myModel.sketches['__profile__']
mySketch.ConstructionLine(point1 = (0,0), point2 = (0,1))
mySketch.Line(point1 = (R, 0), point2 = (R, w*2))
myModel.Part(dimensionality=THREE_D, name='support', type=ANALYTIC_RIGID_SURFACE)
supportPart = myModel.parts['support']
supportPart.AnalyticRigidSurfRevolve(sketch=mySketch)
del myModel.sketches['__profile__']
supportPart.ReferencePoint(point=(0,0,0))
# # sketch punch
myModel.ConstrainedSketch(name='__profile__', sheetSize = w*2)
mySketch = myModel.sketches['__profile__']
beta=np.deg2rad(2)
# arc at center
mySketch.ArcByCenterEnds(center=(0,0), 
                         point1=(-r*np.cos(beta), -r*np.sin(beta)),
                         point2=( r*np.cos(beta), -r*np.sin(beta)))
# lines on left side
mySketch.Line(  point1=(-r*np.cos(beta), -r*np.sin(beta)),
                point2=(-r*np.cos(beta)-0.5*h*np.tan(beta), -r*np.sin(beta)+0.5*h))
mySketch.Line(  point1=(-r*np.cos(beta)-0.5*h*np.tan(beta), -r*np.sin(beta)+0.5*h),
                point2=(-r*np.cos(beta)-0.5*h*np.tan(beta), -r*np.sin(beta)+1.0*h))
# lines on right side
mySketch.Line(  point1=( r*np.cos(beta), -r*np.sin(beta)),
                point2=( r*np.cos(beta)+0.5*h*np.tan(beta), -r*np.sin(beta)+0.5*h))
mySketch.Line(  point1=( r*np.cos(beta)+0.5*h*np.tan(beta), -r*np.sin(beta)+0.5*h),
                point2=( r*np.cos(beta)+0.5*h*np.tan(beta), -r*np.sin(beta)+1.0*h))
myModel.Part(dimensionality=THREE_D, name='punch', type=ANALYTIC_RIGID_SURFACE)
punchPart = myModel.parts['punch']
punchPart.AnalyticRigidSurfExtrude(depth=w*2, sketch=mySketch)
del myModel.sketches['__profile__']
punchPart.ReferencePoint(point=(0,0,0))
# # # Assembly
myAsm = myModel.rootAssembly
myAsm.DatumCsysByDefault(CARTESIAN)
# import plate part
myAsm.Instance(dependent=OFF, name='plate', part=platePart)
plateAsm = myAsm.instances['plate']
# import left support
myAsm.Instance(dependent=OFF, name='support-left', part=supportPart)
myAsm.rotate(angle=90, axisPoint=(0,0,0), axisDirection=(1,0,0), instanceList=('support-left', ))
myAsm.translate(vector=(-supportDist/2., -R, -w/2.), instanceList=('support-left', ))
myAsm.Set(name='refLeftSup', referencePoints=(myAsm.instances['support-left'].referencePoints[2], ))
# import right support
myAsm.Instance(dependent=OFF, name='support-right', part=supportPart)
myAsm.rotate(angle=90, axisPoint=(0,0,0), axisDirection=(1,0,0), instanceList=('support-right', ))
myAsm.translate(vector=(supportDist/2., -R, -w/2.), instanceList=('support-right', ))
myAsm.Set(name='refRightSup', referencePoints=(myAsm.instances['support-right'].referencePoints[2], ))
# import punch
myAsm.Instance(dependent=OFF, name='punch', part=punchPart)
myAsm.translate(vector=(0, r+t, w/2.), instanceList=('punch', ))
myAsm.Set(name='refPunch', referencePoints=(myAsm.instances['punch'].referencePoints[2], ))
# # # Mesh
# global mesh
myAsm.seedPartInstance(size = meshSizeGlobal, regions = (plateAsm, ))
# local mesh
# left support contact area
# z=0
myAsm.seedEdgeBySize(size=meshSizeTrans, edges=plateAsm.edges.getByBoundingBox(
    xMin= -(supportDist+R)/2., xMax= -(supportDist-R)/2., yMin=0.0*t, yMax=t, zMin=0, zMax=0), constraint=FINER) # transition area
# myAsm.seedEdgeBySize(size=meshSizeLocal, edges=plateAsm.edges.getByBoundingBox(
    # xMin= -(supportDist+R)/2., xMax= -(supportDist-R)/2., yMin=0, yMax=0.2*t, zMin=0, zMax=0), constraint=FINER) # contact area
# z=w
myAsm.seedEdgeBySize(size=meshSizeTrans, edges=plateAsm.edges.getByBoundingBox(
    xMin= -(supportDist+R)/2., xMax= -(supportDist-R)/2., yMin=0.0*t, yMax=t, zMin=w, zMax=w), constraint=FINER) # transition area
# myAsm.seedEdgeBySize(size=meshSizeLocal, edges=plateAsm.edges.getByBoundingBox(
    # xMin= -(supportDist+R)/2., xMax= -(supportDist-R)/2., yMin=0, yMax=0.2*t, zMin=w, zMax=w), constraint=FINER) # contact area
# right support contact area
# z=0
myAsm.seedEdgeBySize(size=meshSizeTrans, edges=plateAsm.edges.getByBoundingBox(
    xMin= (supportDist-R)/2., xMax= (supportDist+R)/2., yMin=0.0*t, yMax=t, zMin=0, zMax=0), constraint=FINER) # transition area
# myAsm.seedEdgeBySize(size=meshSizeLocal, edges=plateAsm.edges.getByBoundingBox(
    # xMin= (supportDist-R)/2., xMax= (supportDist+R)/2., yMin=0, yMax=0.2*t, zMin=0, zMax=0), constraint=FINER) # contact area
# z=w
myAsm.seedEdgeBySize(size=meshSizeTrans, edges=plateAsm.edges.getByBoundingBox(
    xMin= (supportDist-R)/2., xMax= (supportDist+R)/2., yMin=0.0*t, yMax=t, zMin=w, zMax=w), constraint=FINER) # transition area
# myAsm.seedEdgeBySize(size=meshSizeLocal, edges=plateAsm.edges.getByBoundingBox(
    # xMin= (supportDist-R)/2., xMax= (supportDist+R)/2., yMin=0, yMax=0.2*t, zMin=w, zMax=w), constraint=FINER) # contact area
# punch contact area
# z=0
myAsm.seedEdgeBySize(size=meshSizeTrans, edges=plateAsm.edges.getByBoundingBox(
    xMin= -2*r, xMax= 2*r, yMin=0, yMax=0.8*t, zMin=0, zMax=0), constraint=FINER) # transition area     
myAsm.seedEdgeBySize(size=meshSizeLocal, edges=plateAsm.edges.getByBoundingBox(
    xMin= -2*r, xMax= 2*r, yMin=0.8*t, yMax=t, zMin=0, zMax=0), constraint=FINER) # contact area
myAsm.seedEdgeBySize(size=meshSizeLocal, edges=plateAsm.edges.getByBoundingBox(
    xMin= -2*r, xMax= 2*r, yMin=0, yMax=0.2*t, zMin=0, zMax=0), constraint=FINER) # critical area
# z=w
myAsm.seedEdgeBySize(size=meshSizeTrans, edges=plateAsm.edges.getByBoundingBox(
    xMin= -2*r, xMax= 2*r, yMin=0, yMax=0.8*t, zMin=w, zMax=w), constraint=FINER) # transition area     
myAsm.seedEdgeBySize(size=meshSizeLocal, edges=plateAsm.edges.getByBoundingBox(
    xMin= -2*r, xMax= 2*r, yMin=0.8*t, yMax=t, zMin=w, zMax=w), constraint=FINER) # contact area
myAsm.seedEdgeBySize(size=meshSizeLocal, edges=plateAsm.edges.getByBoundingBox(
    xMin= -2*r, xMax= 2*r, yMin=0, yMax=0.2*t, zMin=w, zMax=w), constraint=FINER) # critical area
# mesh control
myAsm.setMeshControls(technique=SWEEP, algorithm=ADVANCING_FRONT, elemShape=HEX_DOMINATED,
    regions=plateAsm.cells, allowMapped=False)
# set mapped mesh
myAsm.setMeshControls(allowMapped=True, regions=plateAsm.cells.getByBoundingBox(
    xMin= -(supportDist+R)/2., xMax= -(supportDist-R)/2., yMin=0, yMax=0.2*t, zMin=0, zMax=w))
myAsm.setMeshControls(allowMapped=True, regions=plateAsm.cells.getByBoundingBox(
    xMin= (supportDist-R)/2., xMax= (supportDist+R)/2., yMin=0, yMax=0.2*t, zMin=0, zMax=w))
myAsm.setMeshControls(allowMapped=True, regions=plateAsm.cells.getByBoundingBox(
    xMin= -r, xMax= r, yMin=0, yMax=t, zMin=0, zMax=w))
# generate mesh
myAsm.generateMesh(regions = (plateAsm, ))
print('Total nodes: ' + str(len(plateAsm.nodes)))
print('Total elements: ' + str(len(plateAsm.elements)))