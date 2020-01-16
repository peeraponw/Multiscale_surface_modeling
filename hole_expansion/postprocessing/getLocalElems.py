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
localPoint = (4997.49,32,182.451)
localRadius = 5
shiftRatio = 0.1
# # # Define name
instanceName = 'subPart'
stepName = 'move'
# # Start extraction
x, y, z = localPoint
shiftVect = [x, z]/np.linalg.norm([x, z])
x_shift, z_shift = np.array([x, z]) + shiftVect*localRadius*shiftRatio
localPoint = (x_shift, y, z_shift)
mdb = openMdb(path + caeName + '.cae')
caeAsm = mdb.models['Model-1'].rootAssembly.instances['subPart']
elemLocal = caeAsm.elements.getByBoundingSphere(center=localPoint, radius=localRadius)
nLocalSetsElems = len(elemLocal)
localElemInSetsLabel = [elemLocal[i].label for i in range(0, nLocalSetsElems)]
mdb.close()
with open('localElems'+'.csv', 'w') as out:
    for e in localElemInSetsLabel:
        out.write(str(e)+'\n')
    