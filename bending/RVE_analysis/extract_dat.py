filename = 'rve03'

# # # extract variables
from time import time
#allElems = [elem(idx, myStep) for idx in range(0, nElems)]
print 'Extracting variables from ' + str(nElems) + ' elements.'
t0 = time()
allElems = []
for idx in range(0, nElems):
    allE = elem(idx, myStep)
    allElems.append(allE)
    print 'Element ' + str(idx+1) + ' processed out of ' + str(nElems) + '. (' + str(round(float(idx+1)*100/nElems), 2) + '%)'
    print 'Elapsed time ' + str(round(time() - t0, 2))
ftriax = open(filename + '_triax.csv', 'w')
flode = open(filename + '_lode.csv', 'w')
fpeeq = open(filename + '_peeq.csv', 'w')
# fmises = open(filename + '_mises.csv', 'w')
flabel = open(filename + '_label.csv', 'w')
fvolume = open(filename + '_volume.csv', 'w')
for i in range(0, nElems):
    s = str(allElems[i].triax.tolist())[1:-1] + '\n'
    ftriax.write(s)
    s = str(allElems[i].lode.tolist())[1:-1] + '\n'
    flode.write(s)
    s = str(allElems[i].peeq.tolist())[1:-1] + '\n'
    fpeeq.write(s)
    # s = str(allElems[i].mises.tolist())[1:-1] + '\n'
    # fmises.write(s)
    s = str(allElems[i].label) + '\n'
    flabel.write(s)
    s = str(allElems[i].volume) + '\n'
    fvolume.write(s)
ftriax.close()
flode.close()
fpeeq.close()
# fmises.close()
flabel.close()
fvolume.close()