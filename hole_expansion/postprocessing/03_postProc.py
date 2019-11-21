from abaqus import *
from abaqusConstants import *
from odbAccess import *
from odbAccess import *
from odbSection import *
import numpy as np
from time import time
import sys
import os

odbName = 'HET_1' # without .odb
odbDir = 'P:/abaqus_result/trial_hole_expansion/macro_scale/HET_1/'

fullpath = os.getcwd()

dirName = fullpath.split('\\')
if len(dirName) == 1: 
    # if it's Linux path
    dirName = fullpath.split('/')
odbName = dirName[-1]

elemLocalName = 'LOCAL'
instanceName = 'PLATE'
stepName = 'move'
isEcho = True


# symmFac = float(odbName[3])   # extract from file name
symmFac = 1                     # manual input
if symmFac == 0:
    revolveAngle = 3            # please fill in if symmFac==0


class elem:
    '''
    Element object which contains local variables.
    This object has nothing more than variable container.
    '''
    def __init__(self, idx, myStep):
        # Get number of frames
        nFrames = len(myStep.frames)-1
        # evaluate stress-related variables
        p = np.empty((nFrames))
        q = np.empty((nFrames))
        r = np.empty((nFrames))
        peeq = np.empty((nFrames))
        v = np.empty((nFrames))
        for j in range(0, nFrames):
            myFrame = myStep.frames[j+1]
            fieldSvalues = myFrame.fieldOutputs['S'].values[idx]
            try:
                fieldPEEQvalues = myFrame.fieldOutputs['PEEQ'].values[idx]
            except:
                fieldPEEQvalues = myFrame.fieldOutputs['SDV1'].values[idx]
            p[j] = fieldSvalues.press
            q[j] = fieldSvalues.mises
            r[j] = fieldSvalues.inv3
            peeq[j] = fieldPEEQvalues.data
            try:
                v[j] = myFrame.fieldOutputs['EVOL'].values[idx].data 
            except:
                v[j] = 1
        p[np.abs(p) < 1e-6] = 1e-6
        q[np.abs(q) < 1e-6] = 1e-6
        eta = -p/q
        J2 = 1/3. * q**2
        arg = (r/q)**3
        arg = arg.clip(-1, 1)
        theta = 1/3. * np.arccos(arg)
        lode = 1 - (6*theta/pi)
        # assign values
        self.press = p
        self.mises = q
        self.triax = eta
        self.lode = lode
        self.peeq = peeq
        self.label = myFrame.fieldOutputs['S'].values[idx].elementLabel
        self.volume = v
def updateNodeLabelToIdx(label, values):
    '''
    This function receive node label and its field values variable.
    Return idx which refers to corresponding node.
    It also update global array to speed up idx searching process.
    '''
    if label in AllNodesLabel: 
        # if this label has been searched before, the label exists in AllNodesLabel
        # return idx accordingly
        return AllNodesIdx[AllNodesLabel.index(label)]
    if label == values[label-1].nodeLabel: 
        # if label has not yet been searched, assume firstly idx = label-1
        # update AllNodesLabel and AllNodesIdx and return idx = label-1
        AllNodesLabel[label-1] = label
        AllNodesIdx[label-1] = label-1
        return label-1
    for i in range(len(AllNodesIdx)):
        # if idx != label-1, search it from the beginning and find the correct one.
        # update AllNodesLabel and AllNodesIdx and return idx accordingly 
        if AllNodesIdx[i] != -1:
            continue # continue if this idx is searched already
        else:
            AllNodesLabel[i] = values[i].nodeLabel
            AllNodesIdx[i] = i
            if AllNodesLabel[i] == label:
                return i
def updateNodeCoordLabelToIdx(label, instance):
    '''
    This function receive node label and its instance.
    Return idx which refers to corresponding node.
    It also update global array to speed up idx searching process.
    '''
    if label in AllNodesCoordLabel: 
        # if this label has been searched before, the label exists in AllNodesLabel
        # return idx accordingly
        return AllNodesCoordIdx[AllNodesCoordLabel.index(label)]
    if label == instance.nodes[label-1].label: 
        # if label has not yet been searched, assume firstly idx = label-1
        # update AllNodesLabel and AllNodesIdx and return idx = label-1
        AllNodesCoordLabel[label-1] = label
        AllNodesCoordIdx[label-1] = label-1
        return label-1
    for i in range(len(AllNodesCoordIdx)):
        # if idx != label-1, search it from the beginning and find the correct one.
        # update AllNodesLabel and AllNodesIdx and return idx accordingly 
        if AllNodesCoordIdx[i] != -1:
            continue # continue if this idx is searched already
        else:
            AllNodesCoordLabel[i] = instance.nodes[i].label
            AllNodesCoordIdx[i] = i
            if AllNodesCoordLabel[i] == label:
                return i
def updateElemLabelToIdx(label, values):
    '''
    This function receive element label and its field values variable.
    Return idx which refers to corresponding element.
    It also update global array to speed up idx searching process.
    '''
    if label in AllElemsLabel:
        # if this label has been searched before, the label exists in AllNodesLabel
        # return idx accordingly
        return AllElemsIdx[AllElemsLabel.index(label)]
    if label == values[label-1].elementLabel:
        # if label has not yet been searched, assume firstly idx = label-1
        # update AllNodesLabel and AllNodesIdx and return idx = label-1
        AllElemsLabel[label-1] = label
        AllElemsIdx[label-1] = label-1
        return label-1
    for i in range(len(AllElemsIdx)):
        # if idx != label-1, search it from the beginning and find the correct one.
        # update AllNodesLabel and AllNodesIdx and return idx accordingly 
        if AllElemsIdx[i] != -1:
            continue # continue if this idx is searched already
        else:
            AllElemsLabel[i] = values[i].elementLabel
            AllElemsIdx[i] = i
            if AllElemsLabel[i] == label:
                return i
def findNearest(array, value):
    '''
    Find a point in the "array" which is the nearest to "value".
    Return index "idx" and value "array[idx]"
    '''
    idx = (np.abs(array - value)).argmin()
    return (idx, array[idx])
def extractFD(filename):
    '''
    Read file "filename" and extract its "force" and "disp"
    Displacement and force must be at the first and second column correspondingly
    '''
    f = open(filename, 'rU')
    lines = f.read().rsplit('\n')
    f.close()
    lines = lines[1:-1]
    force = np.zeros((len(lines), 1))
    disp = np.zeros((len(lines), 1))
    for idx, line in enumerate(lines):
        words = line.rsplit(';')
        disp[idx] = float(words[0])
        force[idx] = float(words[1])
    return (disp, force)
def adjustFD(refD, refF, simD, simF):
    '''
    To compare data point from reference displacement "refD" and
    simulation displacement "simD". Due to the nature of "refD" which should be
    finer than "simD", "simD" and corresponding simulation force "simF" shall be
    interpolated to achieve equal number of data points for reference and 
    simulation array. This function process on data from displacement = 0 to 
    displacement = lastD.
    Return adjusted simulation disp "adjD" and adjusted simulation force "adjF"
    which have as many sampling points as refD and refF respectively.    
    '''
    if refD[-1] < simD[-1]: # if experiment breaks before, shorten the simulation result
        # The last simulated displacement is the first element of a set of 
        # simulated displacement which is bigger than refD
        lastSimD = simD[simD >= refD[-1]][0]
        # Shorten simD to the last element of refD
        simD = simD[simD < refD[-1]]
        # Append simD with another point to avoid error
        simD = np.r_[simD, lastSimD]
    else: # if simulation is shorter, trail zeros behind
        trailD = refD[refD > simD[-1]]
        if np.shape(trailD[0] > 1):
            trailD = np.reshape(trailD, (np.shape(trailD)[0], 1))
        trailF = np.zeros_like(trailD)
        simD = np.r_[simD, trailD]
        simF = np.r_[simF, trailF]
    
    adjD = np.zeros_like(refD)
    adjF = np.zeros_like(refD)
    for idx, disp in enumerate(refD):
        nearestIdx, nearestDisp = findNearest(simD, disp)
        if nearestDisp > disp:
            upperDisp = nearestDisp
            lowerDisp = simD[nearestIdx-1]
            upperForce = simF[nearestIdx]
            lowerForce = simF[nearestIdx-1]
        elif nearestDisp < disp:
            lowerDisp = nearestDisp
            upperDisp = simD[nearestIdx+1]
            lowerForce = simF[nearestIdx]
            upperForce = simF[nearestIdx+1]
        else:
            adjD[idx] = disp
            adjF[idx] = simF[nearestIdx]
            continue
        slope = (upperForce - lowerForce)/(upperDisp - lowerDisp)
        adjD[idx] = disp
        adjF[idx] = lowerForce + slope * (adjD[idx] - lowerDisp)
    return adjD, adjF
def shortenRegion(disp, force, lastD):
    d = disp[disp<=lastD]
    idx = len(d)
    f = force[:idx]
    return d, f
def getVolumeAvgVar(elemList):
    nElemsSet = len(elemList)
    nFrames = len(elemList[0].triax)-1
    
    volume = np.empty((nFrames, nElemsSet))
    triax = np.empty((nFrames, nElemsSet))
    lode = np.empty((nFrames, nElemsSet))
    peeq = np.empty((nFrames, nElemsSet))
    for idx, e in enumerate(elemList):
        volume[:, idx] = e.volume
        triax[:, idx] = e.triax
        lode[:, idx] = e.lode
        peeq[:, idx] = e.peeq
    sumvolume = np.sum(volume, axis=1)
    avgTriax = np.sum(triax*volume, axis=1)/sumvolume
    avgLode = np.sum(lode*volume, axis=1)/sumvolume
    avgPeeq = np.sum(peeq*volume, axis=1)/sumvolume
    return avgPeeq, avgTriax, avgLode
def echo(s):
    with open(odbName+'.log', 'a') as logfile: # append the log file after initiation previously.
        logfile.write(s)
    print(s)
def getElapseTime(sec):
    m, s = divmod(sec, 60)
    h, m = divmod(m, 60)
    return int(h), int(m), round(s, 3)
    
    
def getNodeCoord(nIdx, cIdx, myFrame):
    value = myFrame.fieldOutputs['U'].values[nIdx]
    disp = value.data
    coord0 = value.instance.nodes[cIdx].coordinates
    coord_new = coord0+disp
    return coord_new[:2]
    
def findIntersection(coordL0, coordL1, coordR0, coordR1):
    coordCen = (0,0)
    # find m (slope) and c (intercept) of left end
    slopeL = (coordL0[1] - coordL1[1])/(coordL0[0] - coordL1[0]) # m = (y0-y1)/(x0-x1)
    if abs(slopeL) <= 1e-6:
        interceptL = coordL0[1]
    else:
        interceptL = coordL0[1]-(slopeL * coordL0[0]) # c = y0-(m*x0)
    # find m (slope) and c (intercept) of left end
    slopeR = (coordR1[1] - coordR0[1])/(coordR1[0] - coordR0[0]) # m = (y0-y1)/(x0-x1)
    if abs(slopeR) <= 1e-6:
        interceptR = coordR0[1]
    else:
        interceptR = coordR0[1]-(slopeR * coordR0[0]) # c = y0-(m*x0)
    # # 
    if abs(slopeR - slopeL) <= 1e-6:
        x = 0
        y = 0
    else:
        x = (interceptL - interceptR)/(slopeR - slopeL)
        y = slopeL * x + interceptL
    coordCen = (x, y)
    return coordCen

    
def getBendingAngle(nodesL, nodesR, nodesCoordL, nodesCoordR, myFrame):    
    coordL0 = getNodeCoord(nodesL[0], nodesCoordL[0], myFrame)
    coordL1 = getNodeCoord(nodesL[1], nodesCoordL[1], myFrame)
    coordR0 = getNodeCoord(nodesR[0], nodesCoordR[0], myFrame)
    coordR1 = getNodeCoord(nodesR[1], nodesCoordR[1], myFrame)
    coordCen = findIntersection(coordL0, coordL1, coordR0, coordR1)
    leftDist = np.linalg.norm(coordCen-coordL0)
    rightDist = np.linalg.norm(coordCen-coordR0)
    lrDist = np.linalg.norm(coordL0 - coordR0)
    angle = np.arccos((lrDist**2 - leftDist**2 - rightDist**2)/(-2*leftDist*rightDist))
    return np.rad2deg(angle)
        
# # # # Open ODB and initialize variables
# myOdb = openOdb(path = odbName + '.odb')

odbPath = odbDir+odbName # concatenate odb file and dir
if os.path.isfile(odbPath+'_upgraded.odb'): 
    # if it was uphraded before, use the upgraded one
    odbPath = odbPath+'_upgraded' 
if odbAccess.isUpgradeRequiredForOdb(odbPath + '.odb'):
    # if odb needs to be upgraded, create a new one with '_upgraded' suffix
	odbAccess.upgradeOdb(existingOdbPath=odbPath + '.odb',
		upgradedOdbPath=odbPath + '_upgraded' + '.odb')
	myOdb = openOdb(path = odbPath + '_upgraded'+'.odb')
else:
	myOdb = openOdb(path = odbPath + '.odb')

myAsm = myOdb.rootAssembly
myInstance = myAsm.instances[instanceName]
myStep = myOdb.steps[stepName]
nFrames =  len(myStep.frames)-1
nElems = len(myInstance.elements)
nNodes = len(myInstance.nodes)
prevTime = time()

# # initialize nodeLabel and nodeIdx array as well as elemLabel and ElemIdx
# The concept is due to inconsistency in node/element label and index.
# Example:
# For model with small amount of elements, one can refer to node 1's reaction force (label=1) by:
## myFrame.fieldOutputs['RF'].values[1-1].data
# However, in case of huge model, the reaction force obtained by
## myFrame.fieldOutputs['RF'].values[10000-1].data
# may not be the reaction force of node 10000.
# The number 1 and 10000 in this example are called "label" which are shown as node number in ABAQUS/CAE,
# while 1-1 and 10000-1 in the brackets are called "index" or "idx". 
# Therefore, the functions "updateNodeLabelToIdx" and "updateElemLabelToIdx" are made 
# to convert label into idx.
global AllNodesLabel, AllNodesIdx, AllElemsLabel, AllElemsIdx, AllNodesCoordIdx
AllNodesLabel = [-1]*nNodes
AllNodesIdx = [-1]*nNodes
AllNodesCoordLabel = [-1]*nNodes
AllNodesCoordIdx = [-1]*nNodes
AllElemsLabel = [-1]*nElems
AllElemsIdx = [-1]*nElems

# # # determine force-disp
yForce = np.empty((nFrames))
yDisp = np.empty((nFrames))
for fr in range(0, nFrames):
    myFrame = myStep.frames[fr+1]
    yDisp[fr] = myFrame.fieldOutputs['U'].getSubset(region=myAsm.nodeSets['REFPUNCH']).values[0].data[1]/1000
    if symmFac != 0:
        yForce[fr] = 4./symmFac *myFrame.fieldOutputs['RF'].getSubset(region=myAsm.nodeSets['REFPUNCH']).values[0].data[1]/1e9
    else:
        yForce[fr] = 360./revolveAngle *myFrame.fieldOutputs['RF'].getSubset(region=myAsm.nodeSets['REFPUNCH']).values[0].data[1]/1e9
FD = np.c_[yDisp, yForce] # prepare to be written in the output file
session.XYData(data=FD.tolist(), name='force-disp')
FDTitle = 'Disp,Force' # header to be filled in the output file
# # # analyze element
if type(elemLocalName) is str: # elemLocalName is entered as string
    # obtain elements in the LOCAL set
    nLocalSetsElems = len(myAsm.elementSets[elemLocalName].elements[0])
    # obtain elements label in the set
    localElemInSetsLabel = [myAsm.elementSets[elemLocalName].elements[0][i].label for i in range(nLocalSetsElems)]
else: # elemLocaName is entered as list of integer
    # count number of elements
    nLocalSetsElems = len(elemLocalName)
    # element label can be assigned direcly to localElemInSetsLabel
    localElemInSetsLabel = elemLocalName
# initialize element idx array before conversion
localElemInSetsIdx = [-1]*len(localElemInSetsLabel) # initilaize with negative values
# convert label to idx
for i, lbl in enumerate(localElemInSetsLabel):
    localElemInSetsIdx[i] = int(updateElemLabelToIdx(lbl, myStep.frames[0].fieldOutputs['S'].values))
# calculate variables and assign to elem objects
localElems = [elem(idx, myStep) for idx in localElemInSetsIdx]
localVal = np.empty((nFrames, nLocalSetsElems*3)) # 3 columns contain peeq, eta, thetabar
localTitle = '' # header to be filled in the output file
for idx, el in enumerate(localElems):
    # assign variables to local value array
    localVal[:, 3*idx:3*idx+3] = np.c_[el.peeq, el.triax, el.lode]
    localTitle = localTitle + 'PEEQ-' + str(localElemInSetsLabel[idx]) + ','
    localTitle = localTitle + 'eta-' + str(localElemInSetsLabel[idx]) + ','
    localTitle = localTitle + 'thetabar-' + str(localElemInSetsLabel[idx]) + ','
    if isEcho == True:
        pltTriax = np.c_[el.triax, el.peeq].tolist()
        pltLode = np.c_[el.lode, el.peeq].tolist()
        session.XYData(data = pltTriax[1:], name = 'Local Triax vs PEEQ' + '@elem_' + str(localElemInSetsLabel[idx]))
        session.XYData(data = pltLode[1:], name = 'Local Lode vs PEEQ'  + '@elem_' + str(localElemInSetsLabel[idx]))
localTitle = localTitle[:-2] # remove the last two characters which are '; '
# # # Writing results to output file
# if localVal exists, merge localVal with FD
if 'localVal' in locals():
    listVal = np.c_[FD, localVal].tolist()
    title = FDTitle + ',' + localTitle
else:
    listVal = FD.tolist()
    title = FDTitle
with open(odbName + '.csv', 'w') as f:
    f.write( title + '\n')
    for val in listVal:
        s = str(val)[1:-1] # remove '[' and ']'
        f.write( s + '\n')

myOdb.close()