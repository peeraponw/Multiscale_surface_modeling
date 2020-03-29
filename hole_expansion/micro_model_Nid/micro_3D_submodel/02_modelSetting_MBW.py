import numpy as np
from mesh import *
from abaqus import *
from abaqusConstants import *
from caeModules import *
from regionToolset import *


nIntervals = 80

jobName='HET_submodel3D_1_MBW'
materialName='DP1000'
materialFile='DP1000M.inp'

simTime = 1e-5

isMBW = 1
varList = ('S', 'PEEQ', 'U', 'RF', 'EVOL', 'SDV', 'STATUS') # if apply MBW, do not forget SDV 

# ------------------------------------------------------------------------------------------------------------

# do not care this if use mbw input file
youngMod = 210000
poissonRatio = 0.3
density = 7.85e-9 / 1000**3

def readMaterialFromFile(filename):
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
    myMaterial.Plastic(table = flow)

# ------------------------------------------------------------------------------------------------------------
   
# # Create material
myModel.Material(name = materialName)
myMaterial = myModel.materials[materialName]
myMaterial.Elastic(table=((youngMod, poissonRatio), ))
myMaterial.Density(table=((density, ), ))
if isMBW == 0:
    readMaterialFromFile(materialFile)


# # # Section Assignment
myModel.HomogeneousSolidSection(material = materialName, name = 'Section-1', thickness = None)
myPart.Set(name='subPart', cells=myPart.cells)
myPart.SectionAssignment(sectionName = 'Section-1', region = myPart.sets['subPart'])


# # # Create step
myModel.ExplicitDynamicsStep(name = 'move', previous = 'Initial', timePeriod = simTime)


# # # output request
myModel.FieldOutputRequest(name = 'F-Output-1', createStepName = 'move',
        region = MODEL, variables = varList, numIntervals = nIntervals)


# # # apply submodel BC
myModel.SubmodelBC(name='submodelBC', createStepName='move', dof=(1,2,3), 
        globalStep='1', timeScale=ON, region=partAsm.sets['allFaces-pre'])
myModel.setValues(globalJob=odbName)

# ------------------------------------------------------------------------------------------------------------

# # # create job
myAsm.regenerate()
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
