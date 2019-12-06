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

boxsize = 30
surfFile = 'recon2D_Drill3_1_30pnts.csv'
partition_ratio = 0.8

# ------------------------------------------------------------------------------------------------------------

# Import data points
rough_data = np.genfromtxt(surfFile, delimiter=',', skip_header=1)
samplingPt = int(rough_data.item((-1, 0)))

if(samplingPt != boxsize):
    print('wrong input file by boxsize')

## Seperate points from each plane to list
pt = []
for i in range(boxsize + 1):
    pt.append([[x-boxsize/2, y+boxsize/2] for [z,x,y] in rough_data if z == i])

# find minimum height
h = [y+boxsize/2 for [z,x,y] in rough_data]
minH = np.amin(h)
partitionH = partition_ratio * minH

# ------------------------------------------------------------------------------------------------------------

# Start modelling
executeOnCaeStartup()
myModel = mdb.models['Model-1']

# ------------------------------------------------------------------------------------------------------------

# Create surface box
## Draw first 3D wire as base on first plane
s = myModel.ConstrainedSketch(name = '__profile__', sheetSize = 200.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.sketchOptions.setValues(decimalPlaces=6)

topSurfPt = pt[0]
leftEnd = topSurfPt[0][1]
rightEnd = topSurfPt[-1][1]

s.Spline(points = topSurfPt)
s.Line(point1 = (-0.5*boxsize, -0.5*boxsize), point2 = (0.5*boxsize, -0.5*boxsize))
s.Line(point1 = (0.5*boxsize,  -0.5*boxsize), point2 = (0.5*boxsize, rightEnd))
s.Line(point1 = (-0.5*boxsize, -0.5*boxsize), point2 = (-0.5*boxsize, leftEnd))

myPart = myModel.Part(dimensionality=THREE_D, name='Part-1', type=DEFORMABLE_BODY)
myPart.BaseWire(sketch=s)

del myModel.sketches['__profile__']

# ------------------------------------------------------------------------------------------------------------

# Loop to create sketch on plane

for i in range(1, boxsize+1):

	## Create datum plane
	myPart.DatumPlaneByPrincipalPlane(offset=i, principalPlane=XYPLANE)

	## Draw sketch
	s = myModel.ConstrainedSketch(name='__profile__', sheetSize=200.0, transform=
		myPart.MakeSketchTransform(sketchPlane=myPart.datums[2*i], origin=(0,0,i)))
	s.sketchOptions.setValues(decimalPlaces=6)

	topSurfPt = pt[i]

	a = topSurfPt[0]
	leftEnd = a[1]
	b = topSurfPt[-1]
	rightEnd = b[1]

	s.Spline(points = topSurfPt)
	s.Line(point1 = (-0.5*boxsize, -0.5*boxsize), point2 = (0.5*boxsize, -0.5*boxsize))
	s.Line(point1 = (0.5*boxsize,  -0.5*boxsize), point2 = (0.5*boxsize, rightEnd))
	s.Line(point1 = (-0.5*boxsize, -0.5*boxsize), point2 = (-0.5*boxsize, leftEnd))

	myPart.Wire(sketch=s, sketchPlane=myPart.datums[2*i])

	del myModel.sketches['__profile__']

# ------------------------------------------------------------------------------------------------------------

# Solid loft

loftPlane_list = []
for i in range(0, boxsize+1):
    edge_i = (myPart.edges.getByBoundingBox(xMin=-boxsize, xMax=boxsize, yMin=-boxsize, yMax=boxsize, zMin=-0.5+i, zMax=0.5+i))
    loftPlane_list.append(edge_i)

try:
    (mdb.models['Model-1'].parts['Part-1'].SolidLoft(endCondition=NONE, loftsections=loftPlane_list, startCondition=NONE))
except:
    print('Solid loft error')
    sys.exit()

# ------------------------------------------------------------------------------------------------------------

## Create YZ-datum plane for partition
for i in range(1, boxsize):
    myPart.DatumPlaneByPrincipalPlane(offset=i-boxsize/2, principalPlane=YZPLANE)

## Create XZ-datum plane for partition
myPart.DatumPlaneByPrincipalPlane(offset=partitionH, principalPlane=XZPLANE)

## count no. of datum plane
### xy-datum: 2,4,6,...,boxsize*2 = 20 planes (not included most left)
### yz-datum: boxsize*2+2+1, boxsize*2+2+2, ..., boxsize*2+2+(boxsize-1) = 19 planes (not inluded both end sides)
### xz-datum: boxsize*2+2+(boxsize-1)+1 = 1 plane

start_yzDatum = boxsize*2+2
start_xzDatum = boxsize*2+2+(boxsize-1)+1

# ------------------------------------------------------------------------------------------------------------

# Partitioning
## XZ-plane partition
myPart.PartitionCellByDatumPlane(cells=myPart.cells, datumPlane=myPart.datums[start_xzDatum])

## Create set for upperBox
upperBox = myPart.Set(name='upperBox', cells=myPart.cells.getByBoundingBox(
    xMin=-boxsize, xMax=boxsize, yMin=0.2*boxsize, yMax=boxsize, zMin=-boxsize, zMax=boxsize))

## XY-plane partition
for i in range(1, boxsize):
    myPart.PartitionCellByDatumPlane(cells=upperBox.cells, datumPlane=myPart.datums[2*i])

## YZ-plane partition
for i in range(1, boxsize):
    myPart.PartitionCellByDatumPlane(cells=upperBox.cells, datumPlane=myPart.datums[start_yzDatum+i])

# ------------------------------------------------------------------------------------------------------------

# Assembly
myAssembly = mdb.models['Model-1'].rootAssembly
myAssembly.Instance(name='Cut_Part-1', part=myPart, dependent=ON)
myAssembly.makeIndependent(instances=(myAssembly.instances['Cut_Part-1'], ))
