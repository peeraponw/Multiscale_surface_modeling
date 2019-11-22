from abaqus import *
from abaqusConstants import *
from caeModules import *
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
  
def getOriNum(curNum, ptList):
  for k in range(0, len(ptList)):
    if ptList[k][0] == curNum:
      break
  oriNum = int(ptList[k][3])
  return oriNum
def assignSection(c, ori):
  for cc in range(0,len(c)):
    cells = c.findAt(c[cc].pointOn)
    region = regionToolset.Region(cells=cells)
    p.SectionAssignment(region=region, sectionName='Section-'+str(ori), offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='')
  return
