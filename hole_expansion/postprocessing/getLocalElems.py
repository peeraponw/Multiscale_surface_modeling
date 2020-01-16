from abaqus import *
from abaqusConstants import *
from odbAccess import *
from odbAccess import *
from odbSection import *

path = ''
caeName = 'HET_submodel_1_re' # without .cae
if caeName[-4:] == '.cae': caeName = caeName[:-4]
# boxsize = 300
# localPoint:   149-1, rd96: (-39.140285, 45.293041)
#               149-1, rd48: ( 39.065548, 44.745098)
localPoint = (4997.73, 44, 180.929)
localbox = 3
localRatioY = 0.9
# # # Define name
instanceName = 'subPart'
stepName = 'move'
# # Start extraction
x, y, z = localPoint
mdb = openMdb(path + caeName + '.cae')
caeAsm = mdb.models['Model-1'].rootAssembly.instances['subPart']
elemLocal = caeAsm.elements.getByBoundingBox(xMin = x-0.5*localbox, xMax = x+0.5*localbox,
                                             yMin = y-localRatioY*localbox, yMax = y+(1-localRatioY)*localbox,
                                             zMin = z-0.5*localbox, zMax = z+0.5*localbox)
nLocalSetsElems = len(elemLocal)
localElemInSetsLabel = [elemLocal[i].label for i in range(0, nLocalSetsElems)]
mdb.close()
with open('localElems'+'.csv', 'w') as out:
    for e in localElemInSetsLabel:
        out.write(str(e)+'\n')
    