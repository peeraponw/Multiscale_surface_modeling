import numpy as np
from mesh import *
from abaqus import *
from abaqusConstants import *
from caeModules import *
from regionToolset import *

nIntervals = 80
myModel = mdb.models['Model-1']
platePart = myModel.parts['plate']
myAsm = myModel.rootAssembly
plateAsm = myAsm.instances['plate']

jobName = 'bend_sm_f'

materialName = 'DP1000'
materialFile = 'MaterialData_CB_um_mod.inp'
isMBW = 1

moveDistance = 16 *unitFactor
simTime = 0.001

supportFric = 0.20
punchFric = 0.10

varList = ('S', 'PEEQ', 'U', 'RF', 'SDV', 'STATUS') # if apply MBW, do not forget SDV

# # material params
# do not care this if use mbw input file
finalStrain = 5.0
resoStrain = 0.01
youngMod = 210000
poissonRatio = 0.3
density = 7.85e-9 / (unitFactor)**3

if materialName[-1] == '1': # if material == 149-1
    eps0 = 0.0001376
    A = 1586
    n = 0.08684
    k0 = 1100
    Q = 311.1
    beta = 266
    alpha = 0.3
elif materialName[-1] == '2': # if material == 149-2
    eps0 = 0.000007054
    A = 2116
    n = 0.1118
    k0 = 1354
    Q = 571.9
    beta = 277.5
    alpha = 0.45
    
# # # aux functions    
def SwiftFn(strain, eps0, A, n):
    '''
    Swift hardening law: power function
    '''
    stress = A * (strain + eps0)**n
    return stress
def VoceFn(strain, k0, Q, beta):
    '''
    Voce hardening law: saturation function
    '''
    stress = k0 + Q * (1-np.exp(-beta*strain))
    return stress
def AlphaFn(strain, eps0, A, n, k0, Q, beta, alpha):
    '''
    Alpha hardening law: combination of Swift and Voce
    '''
    swiftStrs = SwiftFn(strain, eps0, A, n)
    voceStrs = VoceFn(strain, k0, Q, beta)
    stress = alpha*swiftStrs + (1-alpha)*voceStrs
    return stress
def readFlowcurve(filename, mat):
    '''
    Read flow curve from txt file which composed of only strain-stress
    '''
    with open(filename, 'rU') as matfile:
        lines = matfile.read().rsplit('\n')
        lines.remove(lines[-1]) # remove the last blank line
        flow = []
        for i, line in enumerate(lines):
            line = line.rsplit(',')
            flow.append([float(line[1]), float(line[0])])
    return flow
# # # Create material
myModel.Material(name = materialName)
myMaterial = myModel.materials[materialName]
if materialFile == '':
    if materialName[-1] == '1' or materialName[-1] == '2':
        strain = np.arange(0.0, finalStrain+resoStrain, resoStrain)
        stress = AlphaFn(strain, eps0, A, n, k0, Q, beta, alpha)
        flow = np.c_[stress, strain]
    else:
        readFlowCurve(materialName + '.csv', myMaterial)
    myMaterial.Plastic(table=flow.tolist())
myMaterial.Elastic(table=((youngMod, poissonRatio), ))
myMaterial.Density(table=((density, ), ))
# # # Section Assignment
myModel.HomogeneousSolidSection(material = materialName, name = 'Section-1', thickness = None)
platePart.SectionAssignment(sectionName = 'Section-1', region = Region(cells = platePart.cells))
# # # Define BCs
# symm BC
myModel.ZsymmBC(createStepName = 'Initial', name = 'zsymm',
    region = Region(faces=plateAsm.faces.getByBoundingBox(
    xMin=-l, xMax=l, yMin=0, yMax=t, zMin=0, zMax=0)))
# fixed BC - left
myModel.EncastreBC(createStepName = 'Initial', name = 'supLeft', region = myAsm.sets['refLeftSup'])
# fixed BC - right
myModel.EncastreBC(createStepName = 'Initial', name = 'supRight', region = myAsm.sets['refRightSup'])
# fixed BC - mid
myModel.DisplacementBC(createStepName = 'Initial', name = 'punch', region = myAsm.sets['refPunch'],
    amplitude=UNSET, u1=SET, u2=UNSET, u3=SET, ur1=SET, ur2=SET, ur3=SET)
# move BC
moveAmp = myModel.TabularAmplitude(name='ramp', data=((0, 0), (simTime, 1),))
myModel.ExplicitDynamicsStep(name = 'move', previous = 'Initial', timePeriod = simTime)
myModel.DisplacementBC(createStepName = 'move', name = 'move', 
    region=myAsm.sets['refPunch'], u2=-moveDistance, amplitude='ramp')

# # # interaction proerties
# support 
myModel.ContactProperty('support')
if supportFric == 0:
    myModel.interactionProperties['support'].TangentialBehavior(
        formulation=FRICTIONLESS)
else:
    myModel.interactionProperties['support'].TangentialBehavior(
        formulation=PENALTY, fraction=0.005, maximumElasticSlip=FRACTION, 
        table=((supportFric, ), ))
# punch        
myModel.ContactProperty('punch')
if punchFric == 0:
    myModel.interactionProperties['punch'].TangentialBehavior(
        formulation=FRICTIONLESS)
else:
    myModel.interactionProperties['punch'].TangentialBehavior(
        formulation=PENALTY, fraction=0.005, maximumElasticSlip=FRACTION, 
        table=((punchFric, ), ))
# contact
# left
myAsm.Surface(name='supLeftSurf', side2Faces=
    myAsm.instances['support-left'].faces.getByBoundingCylinder(
    center1=(-supportDist/2., -R, -w/2.), center2=(-supportDist/2., -R, -w/2.+2*w), radius=R+1e-5*unitFactor))
myAsm.Surface(name='plateLeftSurf', side1Faces=
    myAsm.instances['plate'].faces.getByBoundingBox(
    xMin= -(supportDist+R)/2.   , xMax= -(supportDist-R)/2., 
    yMin= 0     , yMax= 0,
    zMin= 0                     , zMax= w
    ))
myModel.SurfaceToSurfaceContactExp(name='support-left', createStepName='Initial', 
    master=myAsm.surfaces['supLeftSurf'], slave=myAsm.surfaces['plateLeftSurf'],
    sliding=FINITE, interactionProperty='support')
# right
myAsm.Surface(name='supRightSurf', side2Faces=
    myAsm.instances['support-right'].faces.getByBoundingCylinder(
    center1=(supportDist/2., -R, -w/2.), center2=(supportDist/2., -R, -w/2.+2*w), radius=R+1e-5*unitFactor))
myAsm.Surface(name='plateRightSurf', side1Faces=
    myAsm.instances['plate'].faces.getByBoundingBox(
    xMin= (supportDist-R)/2.    , xMax= (supportDist+R)/2., 
    yMin= 0     , yMax= 0,
    zMin= 0                     , zMax= w
    ))
myModel.SurfaceToSurfaceContactExp(name='support-right', createStepName='Initial', 
    master=myAsm.surfaces['supRightSurf'], slave=myAsm.surfaces['plateRightSurf'],
    sliding=FINITE, interactionProperty='support')
# mid
myAsm.Surface(name='punchSurf', side2Faces=
    myAsm.instances['punch'].faces.getByBoundingBox(
    xMin=-2*r, xMax=2*r, yMin=t, yMax=h+r, zMin=-w/2., zMax=-w/2.+2*w))
myAsm.Surface(name='plateMidSurf', side1Faces=
    myAsm.instances['plate'].faces.getByBoundingBox(
    xMin= -2*r                  , xMax= 2*r, 
    yMin= t     , yMax= t,
    zMin= 0                     , zMax= w
    ))
myModel.SurfaceToSurfaceContactExp(name='punch', createStepName='Initial', 
    master=myAsm.surfaces['punchSurf'], slave=myAsm.surfaces['plateMidSurf'],
    sliding=FINITE, interactionProperty='punch')
# # # RigidBody Constraint
myModel.RigidBody(name='supLeftConstr', refPointRegion=myAsm.sets['refLeftSup'], 
    surfaceRegion=myAsm.surfaces['supLeftSurf'])
myModel.RigidBody(name='supRightConstr', refPointRegion=myAsm.sets['refRightSup'], 
    surfaceRegion=myAsm.surfaces['supRightSurf'])
myModel.RigidBody(name='punchConstr', refPointRegion=myAsm.sets['refPunch'], 
    surfaceRegion=myAsm.surfaces['punchSurf'])

# output request
myModel.FieldOutputRequest(name = 'F-Output-1', createStepName = 'move',
        region = MODEL, variables = varList, numIntervals = nIntervals)
myModel.FieldOutputRequest(name = 'ForceDispReq', createStepName = 'move', 
        region = myAsm.sets['refPunch'], variables = ('U', 'RF'), numIntervals = nIntervals)
# crete set for critical element
refNodes = plateAsm.nodes.getByBoundingSphere(center=(0, 0, 0), radius=1e-5*unitFactor)
localElems = [refNodes[0].getElements()[idx].label for idx in range(0, len(refNodes[0].getElements()))]
myAsm.Set(name='LOCAL', elements=plateAsm.elements.sequenceFromLabels(localElems))
# create set for angle measurement
leftEndCoord =  (-l/2., 0, -w/2.)
rightEndCoord = ( l/2., 0, -w/2.)
myAsm.Set(name='leftNodes', nodes=plateAsm.nodes.getByBoundingCylinder(
            center1=(-l/2., 0, w), center2=(-l/2.+1.5*meshSizeGlobal, 0, w), radius=1e-5*unitFactor))
myAsm.Set(name='rightNodes', nodes=plateAsm.nodes.getByBoundingCylinder(
            center1=( l/2., 0, w), center2=( l/2.-1.5*meshSizeGlobal, 0, w), radius=1e-5*unitFactor))

# # # create job
myAsm.regenerate()
job = mdb.Job(name=jobName, model=myModel)
job.writeInput()
if isMBW == 1:
    with open(jobName+'.inp', 'rU') as inpfile:
        lines = inpfile.read().rsplit('\n')
    with open(jobName+'.inp', 'w') as outfile:
        skiplines = 7 # the predefined material data generates 7 lines in input file
        counter = 0
        for line in lines:
            if line == '** MATERIALS':
                outfile.write('* INCLUDE, input='+materialFile+'\n')
                counter += 1
            if counter == 0 or counter > skiplines:
                outfile.write(line+'\n')
                continue
            counter += 1
            