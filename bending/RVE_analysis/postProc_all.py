from abaqus import *
from abaqusConstants import *
from odbAccess import *
from odbAccess import *
from odbSection import *
import odbSection
import odbAccess
import multiprocessing
import numpy as np
from time import time
import sys

# # # Parameter section
path = 'postprocessing/'
odbName = 'grind-mbw06' # without .odb
features = ['peeq', 'mises', 'triax', 'lode', 'volume']
if odbName[-4:] == '.inp': odbName = odbName[:-4]
boxsize = 50 
# localPoint:   smooth (-20.8932, 24.8898)
#               grind (10.9539, 23.8227)
#               rough (-15.9562, 22.5945)
localPoint = (10.9539, 23.8227)
localbox = 10
localRatioY = 0.9
isEcho = True
# # # Define name
instanceName = 'CUT_PART-1'
stepName = 'move'

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
    nFrames = len(elemList[0].triax)
    
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
# # # Open ODB and initialize variables
myOdb = openOdb(path = path + odbName + '.odb')
myAsm = myOdb.rootAssembly
myInstance = myAsm.instances[instanceName]
myStep = myOdb.steps[stepName]
nFrames = len(myStep.frames)-1
nElems = len(myInstance.elements)
nNodes = len(myInstance.nodes)
prevTime = time()
initTime = time()
if isEcho == True:
    '''
    isEcho is a variable that indicates if the script will report its status during post processing or not.
    '''
    s = ''' 
        Opening file: {odbName}.odb
        The database contains {nElems} elements with {nFrames} frames.
        '''.format(odbName=odbName, nElems=nElems, nFrames=nFrames)
    with open(odbName+'.log', 'w') as logfile: # Initiate the log file by 'w'
        logfile.write(s)
    print(s)
# # # Initialize arrays for label-idx transformation
global AllNodesLabel, AllNodesIdx, AllElemsLabel, AllElemsIdx
AllNodesLabel = [-1]*nNodes
AllNodesIdx = [-1]*nNodes
AllElemsLabel = [-1]*nElems
AllElemsIdx = [-1]*nElems
# # # Evaluate all element for global variables
globalElems = []
for idx in range(0, nElems):
    globalElems.append(elem(idx, myStep))
    if isEcho == True:
        hr, mn, sc = getElapseTime(time()-initTime)
        s = '''Processing element: {elem}/{nElems}. Elapsed time: {hr} hr, {mn} min, {sc} sec.
            '''.format(elem=idx, nElems=nElems, hr=hr, mn=mn, sc=sc)
        echo(s)
        prevTime = time()

# # # rearrange the data
nFeatures = len(features)
title = ['el']*nElems
for feat in features:
    # initialize feature variables
    exec(feat + '= np.empty((nFrames, nElems))')
for e in globalElems:
    # put element's feature values to their corresponding columns
    lbl = e.label
    title[lbl-1] = title[lbl-1]+'-{:d}'.format(lbl)
    for feat in features:
        exec(feat + '[:, lbl-1] = e.'+feat)
if isEcho == True:
    hr, mn, sc = getElapseTime(time()-initTime)
    s = '''Writing data to file... Elapsed time: {hr} hr, {mn} min, {sc} sec.
        '''.format(hr=hr, mn=mn, sc=sc)
# # # write everything to csv file
import csv
for feat in features:
    with open(feat+'.csv', 'w') as f:
        t = ','.join(title)
        f.write(t+'\n')
        writer = csv.writer(f, delimiter=',', lineterminator='\n')
        exec('writer.writerows('+feat+'.tolist())')      
myOdb.close()