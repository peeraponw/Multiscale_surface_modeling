# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 16:00:44 2016

@author: peeraponw
"""

import numpy as np

filename = 'rve03'

def readFile(filename):
    f = open(filename, 'rU')
    lines = f.read().rsplit('\n')
    dat = []
    for line in lines[:-1]:
        line = line.rsplit(',')
        buf = [float(e) for e in line]    
        dat.append(buf)
    return np.array(dat)
def writeFile(filename, data):
    f = open(filename, 'w')
    for idx in range(0, len(data)):    
        s = str(data[idx, :].tolist())[1:-1] + '\n'
        f.write(s)
    f.close()


label = readFile(filename + '_label.csv')
triax = readFile(filename + '_triax.csv')
lode = readFile(filename + '_lode.csv')
peeq = readFile(filename + '_peeq.csv')
volume = readFile(filename + '_volume.csv')

idx = [int(lbl)-1 for lbl in label]
# check if there are identical elements in idx
if len(idx) != len(set(idx)): 
    print "Error there are identical Element IDs"
    exit()
# Init new array
triaxNew = np.zeros_like(triax)
lodeNew = np.zeros_like(lode)
peeqNew = np.zeros_like(peeq)
volumeNew = np.zeros_like(volume)
#
# Sort array
for i in range(0, len(idx)):
    triaxNew[idx[i], :] = triax[i, :]
    lodeNew[idx[i], :] = lode[i, :]
    peeqNew[idx[i], :] = peeq[i, :]
    volumeNew[idx[i], :] = volume[i, :]

writeFile(filename + '_triax.csv', triaxNew)
writeFile(filename + '_lode.csv', lodeNew)
writeFile(filename + '_peeq.csv', peeqNew)
writeFile(filename + '_volume.csv', volumeNew)