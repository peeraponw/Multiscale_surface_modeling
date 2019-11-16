a = mdb.models['Model-1'].rootAssembly
myAssembly = mdb.models['Model-1'].rootAssembly
a.InstanceFromBooleanMerge(name='Merged_Part-0', instances=myAssembly.instances.values(), 
    keepIntersections=ON, originalInstances=DELETE, domain=GEOMETRY)
p = myModel.parts['Merged_Part-0']
pbox = mdb.models['Model-1'].parts['box']
myAssembly.Instance(name='Merged_Part-TOP', part=p, dependent=ON) # assemble parts
myAssembly.Instance(name='Merged_Part-BOT', part=p, dependent=ON) # assemble parts
myAssembly.Instance(name='Merged_Part-RIGHT', part=p, dependent=ON) # assemble parts
myAssembly.Instance(name='Merged_Part-LEFT', part=p, dependent=ON) # assemble parts
myAssembly.Instance(name='Merged_Part-TOPLEFT', part=p, dependent=ON) # assemble parts
myAssembly.Instance(name='Merged_Part-TOPRIGHT', part=p, dependent=ON) # assemble parts
myAssembly.Instance(name='Merged_Part-BOTLEFT', part=p, dependent=ON) # assemble parts
myAssembly.Instance(name='Merged_Part-BOTRIGHT', part=p, dependent=ON) # assemble parts

myAssembly.translate(instanceList = ['Merged_Part-TOP'], vector = (0., boxsize, 0.))
a.Instance(name='box-1', part=pbox, dependent=ON)
a.InstanceFromBooleanCut(name = 'Merged_Part-TOP_Cut', instanceToBeCut = a.instances['Merged_Part-TOP'],cuttingInstances = (a.instances['box-1'], ),originalInstances = DELETE)

myAssembly.translate(instanceList = ['Merged_Part-BOT'], vector = (0., -0.5*boxsize, 0.))
a.Instance(name='box-2', part=pbox, dependent=ON)
a.InstanceFromBooleanCut(name = 'Merged_Part-BOT_Cut', instanceToBeCut = a.instances['Merged_Part-BOT'],cuttingInstances = (a.instances['box-2'], ),originalInstances = DELETE)

myAssembly.translate(instanceList = ['Merged_Part-RIGHT'], vector = (boxsize, 0., 0.))
a.Instance(name='box-3', part=pbox, dependent=ON)
a.InstanceFromBooleanCut(name = 'Merged_Part-RIGHT_Cut', instanceToBeCut = a.instances['Merged_Part-RIGHT'],cuttingInstances = (a.instances['box-3'], ),originalInstances = DELETE)

myAssembly.translate(instanceList = ['Merged_Part-LEFT'], vector = (-boxsize, 0., 0.))
a.Instance(name='box-4', part=pbox, dependent=ON)
a.InstanceFromBooleanCut(name = 'Merged_Part-LEFT_Cut', instanceToBeCut = a.instances['Merged_Part-LEFT'],cuttingInstances = (a.instances['box-4'], ),originalInstances = DELETE)

myAssembly.translate(instanceList = ['Merged_Part-TOPLEFT'], vector = (-boxsize, boxsize, 0.))
a.Instance(name='box-5', part=pbox, dependent=ON)
a.InstanceFromBooleanCut(name = 'Merged_Part-TOPLEFT_Cut', instanceToBeCut = a.instances['Merged_Part-TOPLEFT'],cuttingInstances = (a.instances['box-5'], ),originalInstances = DELETE)

myAssembly.translate(instanceList = ['Merged_Part-TOPRIGHT'], vector = (boxsize, boxsize, 0.))
a.Instance(name='box-6', part=pbox, dependent=ON)
a.InstanceFromBooleanCut(name = 'Merged_Part-TOPRIGHT_Cut', instanceToBeCut = a.instances['Merged_Part-TOPRIGHT'],cuttingInstances = (a.instances['box-6'], ),originalInstances = DELETE)

myAssembly.translate(instanceList = ['Merged_Part-BOTLEFT'], vector = (-boxsize, -0.5*boxsize, 0.))
a.Instance(name='box-7', part=pbox, dependent=ON)
a.InstanceFromBooleanCut(name = 'Merged_Part-BOTLEFT_Cut', instanceToBeCut = a.instances['Merged_Part-BOTLEFT'],cuttingInstances = (a.instances['box-7'], ),originalInstances = DELETE)

myAssembly.translate(instanceList = ['Merged_Part-BOTRIGHT'], vector = (boxsize, -0.5*boxsize, 0.))
a.Instance(name='box-8', part=pbox, dependent=ON)
a.InstanceFromBooleanCut(name = 'Merged_Part-BOTRIGHT_Cut', instanceToBeCut = a.instances['Merged_Part-BOTRIGHT'],cuttingInstances = (a.instances['box-8'], ),originalInstances = DELETE)

a.InstanceFromBooleanMerge(name='Merged_Part-ALL', instances=myAssembly.instances.values(), 
    keepIntersections=ON, originalInstances=DELETE, domain=GEOMETRY)
# Periodic cut
# Box_Import and cut
a.Instance(name='box-0', part=pbox, dependent=ON)
a = mdb.models['Model-1'].rootAssembly
a.InstanceFromBooleanCut(name='Cut_Part', 
    instanceToBeCut=mdb.models['Model-1'].rootAssembly.instances['Merged_Part-ALL-1'], 
    cuttingInstances=(a.instances['box-0'], ), 
    originalInstances=DELETE)
a.makeIndependent(instances=(a.instances['Cut_Part-1'], ))


# volumeFraction1 = sVol1/myAssembly.getVolume()
# volumeFraction2 = sVol2/myAssembly.getVolume()

# Write volume list to file
# f = open('volume1.csv', 'w')
# for eachVol1 in listVol1:
#   f.write(str(eachVol1)+'\n')
# f.close()

# f = open('volume2.csv', 'w')
# for eachVol2 in listVol2:
#   f.write(str(eachVol2)+'\n')
# f.close()

