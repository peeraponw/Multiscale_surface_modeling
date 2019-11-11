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

# Import nGrains.txt
f = open('nGrains.txt','rU') # this file MUST be encoded in ANSI
nGrains = [int(i) for i in f.read().rsplit()]
f.close()
boxsize = float(nGrains[0])
nGrains.remove(boxsize)
# Import grainDia.txt
g = open('grainDia.txt','rU')
lines = g.read().rsplit()
g.close()
grainDia = []
for i in range(0,len(lines)):
  line = [int(word) if word.isdigit() else float(word) for word in lines[i].rsplit(',')]
  grainDia.append(line)
# Import arcData.txt
f = open('arcData.csv','rU')
lines = f.read().rsplit()
f.close()
# Each geometrical line contains 3 points to construct
pt1 = []
pt2 = []
pt3 = []
for i in range(0,len(lines)):
  line = [int(word) if word.isdigit() else float(word) for word in lines[i].rsplit(',')]
  if i%3 == 0:
    pt1.append(line)
  elif i%3 == 1:
    pt2.append(line)
  elif i%3 == 2:
    pt3.append(line)

nPart = pt3[-1][0] # number of parts in total

# Create a hollowed square part to cut out after periodicity
boxerr = 1
while boxerr:
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
      s.Line(point1 = (-0.5*boxsize, -0.5*boxsize), point2 = ( 0.5*boxsize, -0.5*boxsize))
      s.Line(point1 = (-0.5*boxsize,  0.5*boxsize), point2 = (-0.5*boxsize, -0.5*boxsize))
      s.Line(point1 = ( 0.5*boxsize, -0.5*boxsize), point2 = ( 0.5*boxsize,  0.5*boxsize))
  
  s.rectangle(point1=(-2.0*boxsize, -2.0*boxsize), point2=(2.0*boxsize, 2.0*boxsize)) # Create outer square
  if dimension == '2D':
    p = myModel.Part(name='box', dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)
    try:
      p.BaseShell(sketch=s)
      boxerr = 0
      print('Create box successfully')
    except:
      boxerr = 1
      print('Recreate box')
  else:
    p = myModel.Part(name='box', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    try:
      p.BaseSolidExtrude(sketch=s, depth=boxsize)
      boxerr = 0
      print('Create box successfully')
    except:
      boxerr = 1
      print('Recreate box')
  del myModel.sketches['__profile__']
  
session.viewports['Viewport: 1'].setValues(displayedObject=p)



# # Create parts by connecting points in pt1 to pt2 and pt2 to pt3
for i in range(0,int(nPart)):
  s = myModel.ConstrainedSketch(name = '__profile__', sheetSize = 1.2)
  g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
  s.sketchOptions.setValues(decimalPlaces=6)
  s.setPrimaryObject(option=STANDALONE)
  for k in range(0,len(pt1)):
    if pt1[k][0] == i+1 and pt2[k][0] == i+1 and pt3[k][0] == i+1 :
      sketchPt1 = (pt1[k][1], pt1[k][2])
      sketchPt2 = (pt2[k][1], pt2[k][2])
      sketchPt3 = (pt3[k][1], pt3[k][2])
      if sketchPt1 != sketchPt3:
        if sketchPt2 != sketchPt1:
          s.Line(point1 = (pt1[k][1], pt1[k][2]), point2 = (pt2[k][1], pt2[k][2]))
        if sketchPt2 != sketchPt3:
          s.Line(point1 = (pt2[k][1], pt2[k][2]), point2 = (pt3[k][1], pt3[k][2]))
        ori = int(pt1[k][3])
  
  
  if dimension == '2D':
    p = myModel.Part(name = 'Part-'+str(i+1)+'_'+str(ori), dimensionality = TWO_D_PLANAR, type = DEFORMABLE_BODY)
    p.BaseShell(sketch = s) 
  else:
    p = myModel.Part(name = 'Part-'+str(i+1)+'_'+str(ori), dimensionality = THREE_D, type = DEFORMABLE_BODY)
    p.BaseSolidExtrude(sketch=s, depth=boxsize/100.)
  s.unsetPrimaryObject()
  session.viewports['Viewport: 1'].setValues(displayedObject = p)
  del mdb.models['Model-1'].sketches['__profile__']
  
##n1 = nGrains[0]	# Number of grains 1
##n2 = nGrains[1]	# Number of grains 2

# Assembly
myAssembly = myModel.rootAssembly
myAssembly.DatumCsysByDefault(CARTESIAN)
##listVol1 = []
##listVol2 = []
##sVol1 = 0
##sVol2 = 0
sVol = [0] * len(nGrains) # Initiate n zero-volume elements
listVol = [[] for kk in range(0, len(nGrains))]
for i in range(0,int(nPart)):
  ori = VORONOI.getOriNum(i+1, pt1)
  try:
    p = myModel.parts['Part-'+str(i+1)+'_'+str(ori)]
    myAssembly.Instance(name='Cell-'+str(i+1)+'_'+str(ori), part=p, dependent=ON) # assemble parts
    # Check total volume
    if i < sum(nGrains): # Count only non-periodic parts
####        
        for kk in range(0, len(nGrains)):
            if i < sum(nGrains[0:kk+1]):
                tVol = p.getVolume()
                sVol[kk] += tVol
                listVol[kk].append(tVol/1) 
##      if i < n1: # Phase1
##        tVol = p.getVolume()
##        sVol1 += tVol
##        listVol1.append(tVol/1) # 1 is part extrusion thickness
##      else: # Phase2
##        tVol = p.getVolume()
##        sVol2 += tVol
##        listVol2.append(tVol/1) # 1 is part extrusion thickness
  except:
    pass
if 0:
  # Define material
  # matName = ['phase1', 'phase2', 'phase3', 'phase4']
  matName = ['phase'+str(x+1) for x in range(0, len(nGrains))]
  if len(matName) != len(nGrains):
      print 'Error: Number of material name is not equal to material type. Please respecify the material name.'
  
  for i in range(0,int(nPart)):
    ori = VORONOI.getOriNum(i+1, pt1)
    if i <= sum(nGrains):
      for kk in range(0, len(nGrains)-1):
          if ori <= sum(nGrains[0:kk+1]):
            myModel.Material(name = matName[kk]+'_'+str(i+1))
            myModel.materials[matName[kk]+'_'+str(i+1)].Depvar(n=176)
            myModel.materials[matName[kk]+'_'+str(i+1)].UserMaterial(mechanicalConstants=(i+1, 3.0))
            myModel.HomogeneousSolidSection(name='Section-'+str(i+1),material=matName[kk]+'_'+str(i+1), thickness=None)
          else:
            myModel.Material(name='Martensite')
            myModel.materials['Martensite'].Elastic(table=((0.21, 0.3), ))
            myModel.HomogeneousSolidSection(name='Section-'+str(i+1), material='Martensite', thickness=None)
    try:
      p = myModel.parts['Part-'+str(i+1)+'_'+str(ori)]
      f=p.faces
      for cc in range(0,len(f)):
        faces = f.findAt(f[cc].pointOn)
        region = regionToolset.Region(faces=faces)
        p.SectionAssignment(region=region, sectionName='Section-'+str(i+1), offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='')
    except:
      pass
  # Write graindata.inp
  f = open('graindata.inp', 'w')
  f.write('!MMM Crystal Plasticity Input File'+'\n'+'\n')
  for i in range(0,int(nPart)):
    ori = VORONOI.getOriNum(i+1, pt1)
    if i <= sum(nGrains) and ori <= sum(nGrains[0:-1]):
      f.write('Grain : '+str(i+1)+' : '+str(random.randint(0,360))+' : '+str(random.randint(0,180))+' : '+str(random.randint(0,360))+' : '+str(grainDia[i][1])+'\n')
  
  f.close()
# To be continues in MWMerge2D.py
#region = regionToolset.Region(faces=faces)
#p.Set(faces=faces, name='Set-'+str(i+1))