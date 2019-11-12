import numpy as np
from mesh import *
from abaqus import *
from abaqusConstants import *
from caeModules import *
from regionToolset import *


boxsize = 50
nIntervals = 100
myModel = mdb.models['Model-1']
myAssembly = myModel.rootAssembly
cutPart = myAssembly.instances['Cut_Part-1']

jobName = 'micro_bending_0'

materialName = 'DP1000'
materialFile = 'MaterialData_CB_um_mod.inp'

moveDistance = 0.2*boxsize
simTime = 0.01
simScheme = 'EXPLICIT'
meshSize = 1
localMeshSize = 0.1

isMBW = 1

varList = ('S', 'PEEQ', 'U', 'RF', 'EVOL', 'SDV', 'STATUS') # if apply MBW, do not forget SDV 

# do not care this if use mbw input file
youngMod = 210000
poissonRatio = 0.3
density = 7.85e-9 / 1000**3
'''
def readMaterialFromFile(filename, mat):
    with open(filename, 'rU') as matfile:
        isElasPlas = 0
        isDIL = 0
        isFlow = 0
        isDensity = 0
        flow = []
        for idx, line in enumerate(matfile):
            if isElasPlas == 1:
                line = line.rsplit(',')
                youngMod = float(line[0])
                poissonRatio = float(line[1])
                isElasPlas = 0         
            if isDIL == 1:
                line = line.rsplit(',')
                D1 = float(line[0])
                D2 = float(line[1])
                D3 = float(line[2])
                D4 = float(line[3])
                D5 = float(line[4])
                D6 = float(line[5])
                Gf = float(line[6])
                Dcrit = float(line[7])
                isDIL = 0
            if isFlow == 1:
                if '*' not in line:
                    line = line.rsplit(',')
                    for num in range(0, len(line), 2):
                        flow.append([float(line[num]), float(line[num+1])])
                else:
                    isFlow = 0
            if isDensity == 1:
                isDensity = 0
                density = float(line.rsplit(',')[0])
            if '** E,NU,ceta,eta0,cthetas,cthetat,cthetac,m' in line:
                isElasPlas = 1
            if '** D1,D2,D3,D4,D5,D6,Gf,Dcrit' in line:
                isDIL = 1
            if '** FLOW CURVE' in line:
                isFlow = 1
            if '*Density' in line:
                isDensity = 1
    myMaterial.Elastic(table = ((youngMod, poissonRatio), ))
    myMaterial.Density(table = ((density, ), ))
    myMaterial.Plastic(table = flow)
'''    
# # Create material
myModel.Material(name = materialName)
myMaterial = myModel.materials[materialName]
myMaterial.Elastic(table=((youngMod, poissonRatio), ))
myMaterial.Density(table=((density, ), ))

if dimension == '2D':
    # # Section Assignment
    myModel.HomogeneousSolidSection(material = materialName, name = materialName, thickness = None)
    myModel.parts['Cut_Part'].SectionAssignment(sectionName = materialName, region = Region(faces = myModel.parts['Cut_Part'].faces))
    # # BCs
    myModel.XsymmBC(createStepName = 'Initial', name = 'xsymm',
                    region = Region(edges = cutPart.edges.getByBoundingBox(
                        xMin = -boxsize,    xMax = -0.5*boxsize,
                        yMin = -boxsize,    yMax = boxsize,
                        zMin = -boxsize,    zMax = boxsize )))
    myModel.YsymmBC(createStepName = 'Initial', name = 'ysymm',
                    region = Region(edges = cutPart.edges.getByBoundingBox(
                        xMin = -boxsize,    xMax = boxsize,
                        yMin = -boxsize,    yMax = -0.5*boxsize,
                        zMin = -boxsize,    zMax = boxsize )))
    moveAmplitude = myModel.TabularAmplitude(name = 'ramp', data = ((0, 0), (simTime, 1),))
    if simScheme == 'EXPLICIT':
        myModel.ExplicitDynamicsStep(name = 'move', previous = 'Initial', timePeriod=simTime)
        myModel.DisplacementBC(name = 'move', createStepName = 'move', distributionType = UNIFORM,
               region = Region(edges = cutPart.edges.getByBoundingBox(
                    xMin = 0.5*boxsize,     xMax = boxsize,
                    yMin = -boxsize,        yMax = boxsize,
                    zMin = -boxsize,        zMax = boxsize)),
               u1 = moveDistance, amplitude = 'ramp')
    elif simScheme == 'IMPLICIT':
        myModel.StaticStep(name = 'move', previous = 'Initial', nlgeom = ON, initialInc = 1e-6, minInc = 1e-15, maxNumInc=100000)
        myModel.DisplacementBC(name = 'move', createStepName = 'move', distributionType = UNIFORM,
                           region = Region(edges = cutPart.edges.getByBoundingBox(
                            xMin = 0.5*boxsize,     xMax = boxsize,
                            yMin = -boxsize,        yMax = boxsize,
                            zMin = -boxsize,        zMax = boxsize)),
                           u1 = moveDistance)
    # # # Meshing
    myAssembly.setElementType(elemTypes = (ElemType(elemCode=CPE4R, elemLibrary=STANDARD, secondOrderAccuracy=OFF, 
        hourglassControl=DEFAULT, distortionControl=DEFAULT), ElemType(elemCode=CPE3, elemLibrary = STANDARD)), regions = Region(faces = cutPart.faces))
    myAssembly.setMeshControls(elemShape = QUAD_DOMINATED, regions = cutPart.cells)
    myAssembly.seedPartInstance(size = meshSize, regions = (cutPart, ) )
    myAssembly.seedEdgeBySize(size = localMeshSize, edges = cutPart.edges.getByBoundingBox(
                            xMin = -boxsize,    xMax = boxsize,
                            yMin = 0.4*boxsize, yMax = boxsize,
                            zMin = -boxsize,    zMax = boxsize))
    myAssembly.generateMesh(regions = (cutPart, ))
else:
    # # Section Assignment
    myModel.HomogeneousSolidSection(material = materialName, name = materialName, thickness = None)
    myModel.parts['Cut_Part'].SectionAssignment(sectionName = materialName, region = Region(cells = myModel.parts['Cut_Part'].cells))
    # # BCs
    myModel.XsymmBC(createStepName = 'Initial', name = 'xsymm',
                    region = Region(faces = cutPart.faces.getByBoundingBox(
                        xMin = -boxsize,    xMax = -0.5*boxsize,
                        yMin = -boxsize,    yMax = boxsize,
                        zMin = -boxsize,    zMax = boxsize )))
    myModel.YsymmBC(createStepName = 'Initial', name = 'ysymm',
                    region = Region(faces = cutPart.faces.getByBoundingBox(
                        xMin = -boxsize,    xMax = boxsize,
                        yMin = -boxsize,    yMax = -0.5*boxsize,
                        zMin = -boxsize,    zMax = boxsize )))
    myModel.DisplacementBC(name = 'plane_strain', createStepName = 'Initial', distributionType = UNIFORM,
                    region = Region(cells = cutPart.cells.getByBoundingBox(
                        xMin = -boxsize,        xMax = boxsize,
                        yMin = -boxsize,        yMax = boxsize,
                        zMin = -boxsize,        zMax = boxsize)),
                    u3 = 0)
    moveAmplitude = myModel.TabularAmplitude(name = 'ramp', data = ((0, 0), (0.01, 1),))
    if simScheme == 'EXPLICIT':
        myModel.ExplicitDynamicsStep(name = 'move', previous = 'Initial', timePeriod=0.01)
        myModel.DisplacementBC(name = 'move', createStepName = 'move', distributionType = UNIFORM,
               region = Region(faces = cutPart.faces.getByBoundingBox(
                    xMin = 0.5*boxsize,     xMax = boxsize,
                    yMin = -boxsize,        yMax = boxsize,
                    zMin = -boxsize,        zMax = boxsize)),
               u1 = moveDistance, amplitude = 'ramp')
    elif simScheme == 'IMPLICIT':
        myModel.StaticStep(name = 'move', previous = 'Initial', nlgeom = ON, initialInc = 1e-6, minInc = 1e-15)
        myModel.DisplacementBC(name = 'move', createStepName = 'move', distributionType = UNIFORM,
                           region = Region(faces = cutPart.faces.getByBoundingBox(
                            xMin = 0.5*boxsize,     xMax = boxsize,
                            yMin = -boxsize,        yMax = boxsize,
                            zMin = -boxsize,        zMax = boxsize)),
                           u1 = moveDistance)
    # # # Meshing
    myAssembly.setElementType(elemTypes = (ElemType(elemCode=C3D8R, elemLibrary=EXPLICIT, secondOrderAccuracy=OFF, 
        hourglassControl=DEFAULT, distortionControl=DEFAULT),), regions = Region(cells = cutPart.cells))
    myAssembly.setMeshControls(elemShape = HEX_DOMINATED, regions = cutPart.cells)
    myAssembly.seedPartInstance(size = meshSize, regions = (cutPart, ) )
    myAssembly.seedEdgeBySize(size = localMeshSize, edges = cutPart.edges.getByBoundingBox(
                            xMin = -boxsize,    xMax = boxsize,
                            yMin = 0.80*0.5*boxsize, yMax = 0.5*boxsize+roughnessMean,
                            zMin = -boxsize,    zMax = 0))
    myAssembly.generateMesh(regions = (cutPart, ))

myModel.fieldOutputRequests['F-Output-1'].setValues(variables = varList, numIntervals = nIntervals)
# # # create job
myAssembly.regenerate()
job = mdb.Job(name=jobName, model=myModel)
job.writeInput()
if isMBW == 1:
    nMat = len(myModel.materials)
    with open(jobName+'.inp', 'rU') as inpfile:
        lines = inpfile.read().rsplit('\n')
    with open(jobName+'.inp', 'w') as outfile:
        skiplines = 2+5*nMat # the predefined material data generates 7 lines in input file
        counter = 0
        for line in lines:
            if line == '** MATERIALS':
                outfile.write('* INCLUDE, input='+materialFile+'\n')
                counter += 1
            if counter == 0 or counter > skiplines:
                outfile.write(line+'\n')
                continue
            counter += 1