from abaqus import *
from abaqusConstants import *
from odbAccess import *
from odbAccess import *
from odbSection import *
import numpy as np

path = ''
caeName = 'HET_submodel3D_1' # without .cae
if caeName[-4:] == '.cae': caeName = caeName[:-4]

localPoint = (4997.49, 32, 182.451)
localRadius = 2.5
shiftRatio = 0.1    # 1 = whole sphere in submodel box
		    # -1 = whole sphere out of submodel box
		    # 0.8 = comparable to localRatioY = 0.9 of getByBoundingBox

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
with open('localElems_pnt2'+'.csv', 'w') as out:
    for e in localElemInSetsLabel:
        out.write(str(e)+'\n')
    