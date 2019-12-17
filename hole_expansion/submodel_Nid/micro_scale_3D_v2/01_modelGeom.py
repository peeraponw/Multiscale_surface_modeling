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
import interaction
import load
import mesh
import optimization
import sketch
import job
import visualization
import connectorBehavior
import numpy as np
import sys

# ------------------------------------------------------------------------------------------------------------

boxsize = 300
surfFile = 'recon2D_Mill4_1_300um_21pnts_trim400-1000_waviness.csv'

boxH_Ratio = 1/3.
partition_ratio = 0.9

boxY_min = boxH_Ratio *(boxsize/2)

# ------------------------------------------------------------------------------------------------------------

# Import data points
rough_data = np.genfromtxt(surfFile, delimiter=',', skip_header=1)
samplingPnt = int(rough_data.item((-1, 0)))
pntInterval = int(rough_data.item((1, 1))) - int(rough_data.item((0, 1)))
planeNo = int(boxsize/pntInterval)

if(samplingPnt != boxsize):
    print('wrong input file by boxsize')

## Seperate points from each plane to list
pt = []
for i in range(0, boxsize+pntInterval, pntInterval):
    pt.append([[x-boxsize/2, y+boxsize/2] for [z,x,y] in rough_data if z == i])

# find minimum height
h = [y+boxsize/2 for [z,x,y] in rough_data]
minH = np.amin(h)
partitionH = boxY_min + partition_ratio*(minH - boxY_min)

# ------------------------------------------------------------------------------------------------------------

# Start modelling
executeOnCaeStartup()
myModel = mdb.models['Model-1']

# ------------------------------------------------------------------------------------------------------------

# Create surface box
## Draw first 3D wire as base on first plane
s = myModel.ConstrainedSketch(name = '__profile__', sheetSize = boxsize*2)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.sketchOptions.setValues(decimalPlaces=6)

topSurfPt = pt[0]
leftEnd = topSurfPt[0][1]
rightEnd = topSurfPt[-1][1]

s.Spline(points = topSurfPt)
s.Line(point1 = (-0.5*boxsize, boxY_min), point2 = (0.5*boxsize, boxY_min))
s.Line(point1 = (0.5*boxsize,  boxY_min), point2 = (0.5*boxsize, rightEnd))
s.Line(point1 = (-0.5*boxsize, boxY_min), point2 = (-0.5*boxsize, leftEnd))

myPart = myModel.Part(dimensionality=THREE_D, name='Part-1', type=DEFORMABLE_BODY)
myPart.BaseWire(sketch=s)

del myModel.sketches['__profile__']

# ------------------------------------------------------------------------------------------------------------

# Loop to create sketch on plane

for i in range(pntInterval, boxsize+pntInterval, pntInterval):

	## Create datum plane
	myPart.DatumPlaneByPrincipalPlane(offset=i, principalPlane=XYPLANE)
	key_xyDatum = myPart.datums.keys()[-1]

	## Draw sketch
	s = myModel.ConstrainedSketch(name='__profile__', sheetSize=boxsize*2, transform=
		myPart.MakeSketchTransform(sketchPlane=myPart.datums[key_xyDatum], origin=(0,0,i)))
	s.sketchOptions.setValues(decimalPlaces=6)

	topSurfPt = pt[i/pntInterval]

	a = topSurfPt[0]
	leftEnd = a[1]
	b = topSurfPt[-1]
	rightEnd = b[1]

	s.Spline(points = topSurfPt)
	s.Line(point1 = (-0.5*boxsize, boxY_min), point2 = (0.5*boxsize, boxY_min))
	s.Line(point1 = (0.5*boxsize,  boxY_min), point2 = (0.5*boxsize, rightEnd))
	s.Line(point1 = (-0.5*boxsize, boxY_min), point2 = (-0.5*boxsize, leftEnd))

	myPart.Wire(sketch=s, sketchPlane=myPart.datums[key_xyDatum])

	del myModel.sketches['__profile__']

# ------------------------------------------------------------------------------------------------------------

# Solid loft

loftPlane_list = []
for i in range(0, boxsize+pntInterval, pntInterval):
    edge_i = (myPart.edges.getByBoundingBox(xMin=-boxsize, xMax=boxsize, yMin=-boxsize, yMax=boxsize, zMin=-0.5+i, zMax=0.5+i))
    loftPlane_list.append(edge_i)

try:
    (mdb.models['Model-1'].parts['Part-1'].SolidLoft(endCondition=NONE, loftsections=loftPlane_list, startCondition=NONE))
except:
    print('Solid loft error')
    
# ------------------------------------------------------------------------------------------------------------

## Create YZ-datum plane for partition
key_yzDatum = myPart.datums.keys()[-1]
for i in range(pntInterval, boxsize, pntInterval):
    myPart.DatumPlaneByPrincipalPlane(offset=i-boxsize/2, principalPlane=YZPLANE)

## Create XZ-datum plane for partition
myPart.DatumPlaneByPrincipalPlane(offset=partitionH, principalPlane=XZPLANE)
key_xzDatum = myPart.datums.keys()[-1]

## count no. of datum plane
### xy-datum: 2,4,6,...,boxsize*2 = 20 planes (not included most left)
### yz-datum: boxsize*2+2+1, boxsize*2+2+2, ..., boxsize*2+2+(boxsize-1) = 19 planes (not inluded both end sides)
### xz-datum: boxsize*2+2+(boxsize-1)+1 = 1 plane

# key_yzDatum = boxsize*2+2
# key_xzDatum = boxsize*2+2+(boxsize-1)+1

# ------------------------------------------------------------------------------------------------------------

# Partitioning
## XZ-plane partition
myPart.PartitionCellByDatumPlane(cells=myPart.cells, datumPlane=myPart.datums[key_xzDatum])

## Create set for upperBox
upperBox = myPart.Set(name='upperBox', cells=myPart.cells.getByBoundingBox(
    xMin=-boxsize, xMax=boxsize, yMin=partitionH, yMax=boxsize, zMin=-boxsize, zMax=boxsize))

## XY-plane partition
for i in range(1, planeNo, 1):
    myPart.PartitionCellByDatumPlane(cells=upperBox.cells, datumPlane=myPart.datums[2*i])
    print('xy-partition success: ' + str(i) + '/' + str(planeNo-1))

## YZ-plane partition
for i in range(key_yzDatum+3, key_yzDatum+3+planeNo-1, 1):
    myPart.PartitionCellByDatumPlane(cells=upperBox.cells, datumPlane=myPart.datums[i])
    print('yz-partition success: ' + str(i-2-key_yzDatum) + '/' + str(planeNo-1))

# ------------------------------------------------------------------------------------------------------------

# Assembly
myAssembly = mdb.models['Model-1'].rootAssembly
myAssembly.Instance(name='Cut_Part-1', part=myPart, dependent=ON)
myAssembly.makeIndependent(instances=(myAssembly.instances['Cut_Part-1'], ))
cutPart = myAssembly.instances['Cut_Part-1']