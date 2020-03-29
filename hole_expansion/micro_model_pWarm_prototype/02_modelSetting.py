import numpy as np
from mesh import *
from abaqus import *
from abaqusConstants import *
from caeModules import *
from regionToolset import *

nIntervals = 80

jobName='sub_dir0_2x'
materialName='DP1000'
materialFile='DP1000M.inp'

isMBW = 1

simTime = 0.01
varList = ('S', 'PEEQ', 'U', 'RF', 'SDV', 'STATUS') # if apply MBW, do not forget SDV

# # material params
# do not care this if use mbw input file
youngMod = 210000
poissonRatio = 0.3
density = 7.85e-9 / (unitFactor)**3
# # # Create material
myModel.Material(name = materialName)
myMaterial = myModel.materials[materialName]
myMaterial.Elastic(table=((youngMod, poissonRatio), ))
myMaterial.Density(table=((density, ), ))
# # # Section Assignment
myModel.HomogeneousSolidSection(material = materialName, name = 'Section-1', thickness = None)
p.Set(name='subPart', cells=p.cells)
p.SectionAssignment(sectionName = 'Section-1', region = p.sets['subPart'])
# # # Create step
myModel.ExplicitDynamicsStep(name = 'move', previous = 'Initial', timePeriod = simTime)
# # # output request
myModel.FieldOutputRequest(name = 'F-Output-1', createStepName = 'move',
        region = MODEL, variables = varList, numIntervals = nIntervals)

# # # apply submodel BC
myModel.SubmodelBC(name='submodelBC', createStepName='move', dof=(1,2,3), 
        globalStep='1', timeScale=ON, region=partAsm.sets['allFaces-pre'])
myModel.setValues(globalJob=odbName)
        
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
                if nMat > 1:
                    outfile.write('* INCLUDE, input='+materialFile_local+'\n')
                counter += 1
            if counter == 0 or counter > skiplines:
                outfile.write(line+'\n')
                continue
            counter += 1