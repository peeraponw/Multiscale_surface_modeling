# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 11:35:23 2016

@author: peeraponw
"""

from abaqus import *
from abaqusConstants import *
from caeModules import *
from regionToolset import *
from mesh import *
from driverUtils import executeOnCaeStartup
import displayGroupMdbToolset
import numpy as np

name = 'rve03'

myModel = mdb.models['Model-1']
p = myModel.parts['Cut_Part']
myAssembly = myModel.rootAssembly
a = myAssembly.instances['Cut_Part-1']
# Point param
x = -2
y = 46.369054
box = 10

selectElems = a.elements.getByBoundingBox(xMin = x - 0.5*box, xMax = x + 0.5*box,
                                          yMin = y - 0.8*box, yMax = y + 0.2*box,
                                          zMin = 0, zMax = 0)
idx = [selectElems[i].label for i in range(0, len(selectElems))]

f = open(name+'_rect_label.csv', 'w')
for i in range(0, len(idx)):
    s = str(idx[i]) + '\n'
f.close()
