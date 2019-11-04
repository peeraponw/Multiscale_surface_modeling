import matplotlib.pyplot as plt
from matplotlib import rc
import os, os.path
import numpy as np
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'], 'size':'14'})
rc('text', usetex=True)

def readFile(filename):
    f = open(filename, 'rU')
    lines = f.read().rsplit('\n')
    dat = []
    for line in lines[:-1]:
        line = line.rsplit(',')
        buf = [float(e) for e in line]    
        dat.append(buf)
    return np.array(dat)
def volumeAvg(val, volume):
    wval = val*volume # weighted value
    out = np.sum(wval, axis=0)/np.sum(volume, axis=0)
    return out
    
class result:
    def __init__(self, name):        
        peeq = readFile(name + '_peeq_local4.csv')
        triax = readFile(name + '_triax_local4.csv')
        lode = readFile(name + '_lode_local4.csv')
        volume = readFile(name + '_volume_local4.csv')
        
        self.avgPeeq = volumeAvg(peeq, volume)
        self.avgTriax = volumeAvg(triax, volume)
        self.avgLode = volumeAvg(lode, volume)
        self.peeq = peeq
        self.triax = triax
        self.lode = lode
        self.volume = volume


mesh10 = result('rve10')
mesh06 = result('rve06')
mesh03 = result('rve03')
mesh01 = result('rve01')

elemIdx10 = 1 - 1# 6129 -1
elemIdx06 = 1 - 1# 7386 -1
elemIdx03 = 1 - 1# 12181 -1
elemIdx01 = 1 - 1# 55437 -1

# # Triaxiality Figures 
fig = plt.figure(figsize = (15,5))
axTriax1 = fig.add_subplot(121)
axTriax1.plot(mesh10.peeq[elemIdx10], mesh10.triax[elemIdx10], 'b', label = r'Mesh size = 1.0')
axTriax1.plot(mesh06.peeq[elemIdx06], mesh06.triax[elemIdx06], 'y', label = r'Mesh size = 0.6')
axTriax1.plot(mesh03.peeq[elemIdx03], mesh03.triax[elemIdx03], 'g', label = r'Mesh size = 0.3')
axTriax1.plot(mesh01.peeq[elemIdx01], mesh01.triax[elemIdx01], 'r', label = r'Mesh size = 0.1')
axTriax1.axis([0, 2.5, 0, 1.1])
axTriax1.set_title('Local')
axTriax1.set_xlabel(r'$\bar{\varepsilon}_p$', fontsize = 16)
axTriax1.set_ylabel(r'$\eta$', fontsize = 16)
axTriax1.legend(loc = 4, fontsize = 12)

axTriax2 = fig.add_subplot(122)
axTriax2.plot(mesh10.avgPeeq, mesh10.avgTriax, 'b', label = r'Mesh size = 1.0')
axTriax2.plot(mesh06.avgPeeq, mesh06.avgTriax, 'y', label = r'Mesh size = 0.6')
axTriax2.plot(mesh03.avgPeeq, mesh03.avgTriax, 'g', label = r'Mesh size = 0.3')
axTriax2.plot(mesh01.avgPeeq, mesh01.avgTriax, 'r', label = r'Mesh size = 0.1')
axTriax2.axis([0, 0.4, 0, 0.9])
axTriax2.set_title('Local')
axTriax2.set_xlabel(r'$\bar{\varepsilon}_p$', fontsize = 16)
axTriax2.set_ylabel(r'$\eta$', fontsize = 16)
axTriax2.legend(loc = 4, fontsize = 12)

# # Lode angle parameter Figures
fig2 = plt.figure(figsize = (15,5))

axLode1 = fig2.add_subplot(121)
axLode1.plot(mesh10.peeq[elemIdx10], mesh10.lode[elemIdx10], 'b', label = r'Mesh size = 1.0')
axLode1.plot(mesh06.peeq[elemIdx06], mesh06.lode[elemIdx06], 'y', label = r'Mesh size = 0.6')
axLode1.plot(mesh03.peeq[elemIdx03], mesh03.lode[elemIdx03], 'g', label = r'Mesh size = 0.3')
axLode1.plot(mesh01.peeq[elemIdx01], mesh01.lode[elemIdx01], 'r', label = r'Mesh size = 0.1')
axLode1.axis([0, 2.4, -0.1, 0.5])
axLode1.set_title('Local')
axLode1.set_xlabel(r'$\bar{\varepsilon}_p$', fontsize = 16)
axLode1.set_ylabel(r'$\bar{\theta}$', fontsize = 16)
axLode1.legend(fontsize = 12)

axLode2 = fig2.add_subplot(122)
axLode2.plot(mesh10.avgPeeq, mesh10.avgLode, 'b', label = r'Mesh size = 1.0')
axLode2.plot(mesh06.avgPeeq, mesh06.avgLode, 'y', label = r'Mesh size = 0.6')
axLode2.plot(mesh03.avgPeeq, mesh03.avgLode, 'g', label = r'Mesh size = 0.3')
axLode2.plot(mesh01.avgPeeq, mesh01.avgLode, 'r', label = r'Mesh size = 0.1')
axLode2.axis([0, 0.4, -0.1, 0.5])
axLode2.set_title('Local')
axLode2.set_xlabel(r'$\bar{\varepsilon}_p$', fontsize = 16)
axLode2.set_ylabel(r'$\bar{\theta}$', fontsize = 16)
axLode2.legend(fontsize = 12)

lodeavg01 = mesh01.avgLode
lodeavg03 = mesh03.avgLode
lodeavg06 = mesh06.avgLode
lodeavg10 = mesh10.avgLode 

peeqavg01 = mesh01.avgPeeq
peeqavg03 = mesh03.avgPeeq
peeqavg06 = mesh06.avgPeeq
peeqavg10 = mesh10.avgPeeq

etaavg01 = mesh01.avgTriax
etaavg03 = mesh03.avgTriax
etaavg06 = mesh06.avgTriax
etaavg10 = mesh10.avgTriax