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
import VORONOI
import numpy as np


# rough: SD=1.1504, rmdseed=25
# grind: SD=0.6293, rdmseed=32
# smooth: SD=0.0429, rdmseed=29
roughnessMean = 0 #3.8017 # in micrometer
roughnessSD = 1.1504 # dimensionless
rdmseed = 25 # seed for random generator
randomPtNum = 49 # number of points in the middle to random depth
order = 5      # spline equation's order (cannot be bigger than randomPtNum+1)
if order > randomPtNum+1: order = randomPtNum+1
roughnessPeriodic = 0   # 1 for 4-sided periodic roughness modeling
                        # 0 for top surface roughness modeling
dimension = '3D'
executeOnCaeStartup()
myModel = mdb.models['Model-1']

boxsize = 50
isTopHalf = 1   # 1 for only top-half
                # 0 for full-plate

# ------------------------------------------------------------------------------------------------------------

# Create a hollowed square part to cut out
boxerr = 1
while boxerr:
  ## Sketch part
  s = myModel.ConstrainedSketch(name = '__profile__', sheetSize = 0.05)
  g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
  s.sketchOptions.setValues(decimalPlaces=6)
  if roughnessMean == 0 and roughnessSD == 0:
    s.rectangle(point1=(-0.5*boxsize, -0.5*boxsize), point2=(0.5*boxsize, 0.5*boxsize)) # Create inner square
  else:
    np.random.seed(rdmseed)
    roughH = np.random.normal(roughnessMean, roughnessSD, randomPtNum)
    roughH = np.reshape(roughH, (randomPtNum, 1))
    roughH = np.row_stack((0, roughH, 0))
    
    roughV = np.random.normal(roughnessMean, roughnessSD, randomPtNum)
    roughV = np.reshape(roughV, (randomPtNum, 1))
    roughV = np.row_stack((0, roughV, 0))
    
    samplingPt = np.column_stack(np.arange(boxsize/(randomPtNum+1), boxsize, boxsize/(randomPtNum+1)))
    samplingPt = np.row_stack((0, samplingPt.T, boxsize)) - 0.5*boxsize
    
    botSurfPt = np.column_stack((samplingPt, roughH - 0.5*boxsize)).tolist()
    topSurfPt = np.column_stack((samplingPt, roughH + 0.5*boxsize)).tolist()
    
    rightSurfPt = np.column_stack((roughV + 0.5*boxsize, samplingPt)).tolist()
    leftSurfPt = np.column_stack((roughV - 0.5*boxsize, samplingPt)).tolist()
    
    for i in range(0, len(topSurfPt)-order+1, order):
      s.Spline(points = topSurfPt[i:i+order+1])
      if roughnessPeriodic == 1:
        s.Spline(points = botSurfPt[i:i+order+1])
        s.Spline(points = rightSurfPt[i:i+order+1])
        s.Spline(points = leftSurfPt[i:i+order+1])
    if roughnessPeriodic == 0:
      if isTopHalf == 0:
        s.Line(point1 = (-0.5*boxsize, -0.5*boxsize), point2 = (0.5*boxsize, -0.5*boxsize))
        s.Line(point1 = (0.5*boxsize,  -0.5*boxsize), point2 = (0.5*boxsize, 0.5*boxsize))
        s.Line(point1 = (-0.5*boxsize, -0.5*boxsize), point2 = (-0.5*boxsize,  0.5*boxsize))
      if isTopHalf == 1:
        s.Line(point1 = (-0.5*boxsize, 0*boxsize), point2 = (0.5*boxsize, 0*boxsize))
        s.Line(point1 = (0.5*boxsize,  0*boxsize), point2 = (0.5*boxsize, 0.5*boxsize))
        s.Line(point1 = (-0.5*boxsize, 0*boxsize), point2 = (-0.5*boxsize,  0.5*boxsize))
  
  ## Extrude part
  s.rectangle(point1=(-2.0*boxsize, -2.0*boxsize), point2=(2.0*boxsize, 2.0*boxsize)) # Create outer square
  if dimension == '2D':
    pbox = myModel.Part(name='box', dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)
    try:
      pbox.BaseShell(sketch=s)
      boxerr = 0
      print('Create box successfully')
    except:
      boxerr = 1
      print('Recreate box')
  else:
    pbox = myModel.Part(name='box', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    try:
      pbox.BaseSolidExtrude(sketch=s, depth=boxsize)
      boxerr = 0
      print('Create box successfully')
    except:
      boxerr = 1
      print('Recreate box')
  del myModel.sketches['__profile__']
  
# ------------------------------------------------------------------------------------------------------------

# Create a plate to be cut
## Sketch part
r = myModel.ConstrainedSketch(name = '__profile__', sheetSize = boxsize*2)
r.rectangle(point1=(-0.7*boxsize, -0.7*boxsize), point2=(0.7*boxsize, 0.7*boxsize))

## Extrude part
p = myModel.Part(dimensionality=THREE_D, name='plate', type=DEFORMABLE_BODY)
p.BaseSolidExtrude(depth=boxsize/100.0, sketch=r)
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

