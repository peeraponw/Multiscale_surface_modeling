from abaqus import *
from abaqusConstants import *
from caeModules import *
from regionToolset import *
from mesh import *
from driverUtils import executeOnCaeStartup
import math
import os
import re
import random
import platform
import time
import shutil
import datetime
import part
import material
import section
import assembly
import step
import load
import mesh
import sketch
import numpy as np

# ------------------------------------------------------------------------------------------------------------

boxsize = 30
isTopHalf = 0    # 1 for top-half mesh
                 # 0 for full cubic mesh
surfFile = 'recon_2D_149-1A1_30pnts_domAmp50.csv'

# ------------------------------------------------------------------------------------------------------------

# Start modelling
executeOnCaeStartup()
myModel = mdb.models['Model-1']


# Create sketch
s1 = myModel.ConstrainedSketch(name='profile_1', sheetSize=0.05)
s1.sketchOptions.setValues(decimalPlaces=6)
s1.rectangle(point1=(-2.0*boxsize, -2.0*boxsize), point2=(2.0*boxsize, 2.0*boxsize)


# Create part
myPart = myModel.Part(name='elem', dimensionality=THREE_D, type=DEFORMABLE_BODY)


# Create datum plane
plane_1 = myPart.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=0)


# Create wire
myPart.Wire(sketchPlane=plane_1, sketchPlaneSide=SIDE1, sketchUpEdge=YAXIS, sketch=s1)
# profile_2 = myModel.myPart.Wire(sketchPlane=plane_2, sketchPlaneSide=SIDE1, sketchUpEdge=YAXIS, sketch=s2)


'''
Path

mdb.models[name].parts[name].Wire


Required arguments
sketchPlane
A Datum plane object or a planar Face object specifying the plane on which to sketch.
sketchPlaneSide
A SymbolicConstant specifying the direction of feature creation. Possible values are SIDE1 and SIDE2.
sketchUpEdge
An Edge object or a Datum axis object specifying the vertical (Y) direction of the sketch.
sketch
A ConstrainedSketch object specifying the planar sketch to be revolved.
'''

'''
# Solid loft


Path

mdb.models[name].parts[name].SolidLoft


Required arguments
loftsections
A sequence of sequences of edges specifying the cross-sections to be lofted. Each outer sequence specifies a section through which Abaqus will pass the loft feature. Each outer sequence can be defined as a sequence of edges or as an EdgeArray. The edges specifying a section must form a simple closed profile and must not contain multiple loops.


import part
mdb.models[name].parts[name].allInternalSets[name].edges[i]
mdb.models[name].parts[name].allInternalSurfaces[name].edges[i]
mdb.models[name].parts[name].allSets[name].edges[i]
mdb.models[name].parts[name].allSurfaces[name].edges[i]
mdb.models[name].parts[name].edges[i]
mdb.models[name].parts[name].sets[name].edges[i]
mdb.models[name].parts[name].surfaces[name].edges[i]


# rough_data = np.genfromtxt(surfFile, delimiter=',', skip_header=1)
# samplingPt = len(rough_data)

# topSurfPt = 

# s.Spline(points = topSurfPt)
# if isTopHalf == 0:
#   s.Line(point1 = (-0.5*boxsize, -0.5*boxsize), point2 = (0.5*boxsize, -0.5*boxsize))
#   s.Line(point1 = (0.5*boxsize,  -0.5*boxsize), point2 = (0.5*boxsize, topSurfPt[(-1,1)]))
#   s.Line(point1 = (-0.5*boxsize, -0.5*boxsize), point2 = (-0.5*boxsize, topSurfPt[(0,1)]))
# if isTopHalf == 1:
#   s.Line(point1 = (-0.5*boxsize, 0*boxsize), point2 = (0.5*boxsize, 0*boxsize))
#   s.Line(point1 = (0.5*boxsize,  0*boxsize), point2 = (0.5*boxsize, topSurfPt[(-1,1)]))
#   s.Line(point1 = (-0.5*boxsize, 0*boxsize), point2 = (-0.5*boxsize, topSurfPt[(0,1)]))


# ## Draw outer square
# s.rectangle(point1=(-2.0*boxsize, -2.0*boxsize), point2=(2.0*boxsize, 2.0*boxsize))


# ## Extrude part
# pbox = myModel.Part(name='box', dimensionality=THREE_D, type=DEFORMABLE_BODY)
# try:
#   pbox.BaseSolidExtrude(sketch=s, depth=boxsize)
#   boxerr = 0
#   print('Create box successfully')
# except:
#   boxerr = 1
#   print('Recreate box')
# del myModel.sketches['__profile__']
  
# ------------------------------------------------------------------------------------------------------------

# Create a plate to be cut
## Sketch part
r = myModel.ConstrainedSketch(name = '__profile__', sheetSize = boxsize*2)
r.rectangle(point1=(-0.7*boxsize, -0.7*boxsize), point2=(0.7*boxsize, 0.7*boxsize))

## Extrude part
p = myModel.Part(dimensionality=THREE_D, name='plate', type=DEFORMABLE_BODY)
p.BaseSolidExtrude(depth=boxsize, sketch=r)

del myModel.sketches['__profile__']

# ------------------------------------------------------------------------------------------------------------

# Assembly
## Declare assembly
a = mdb.models['Model-1'].rootAssembly

## Instance box part
a.Instance(name='box-0', part=pbox, dependent=ON)
## Instance plate part
a.Instance(name='Plate_beforeCut', part=p, dependent=ON)

## Cut part by box
a.InstanceFromBooleanCut(name='Cut_Part', 
    instanceToBeCut=mdb.models['Model-1'].rootAssembly.instances['Plate_beforeCut'], 
    cuttingInstances=(a.instances['box-0'], ), 
    originalInstances=DELETE)
a.makeIndependent(instances=(a.instances['Cut_Part-1'], ))

'''