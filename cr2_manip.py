
#from PIL import Image

import rawpy
import os, sys
import numpy as np
from matplotlib import pyplot as plt
from PIL import Image
import pickle

print np.__version__

rawdir = sys.argv[1]
print rawdir

filelist = [a for a in os.listdir(rawdir) if a[-4:]==".cr2"]
print "Importing "+str(len(filelist))+" files."

#rawlist = []
#vislist = []

##### need to redo this so that I'm not load all images into memory at once
## first pass: store sum and sumofsquares
## calculate avg and stddev from sum and sumofsquares for each pixel
## second pass, compute normalized residual for each pixel

firstvis = rawpy.imread(rawdir+filelist[0]).raw_image_visible

collector = np.zeros(firstvis.shape)
runSum = np.zeros(firstvis.shape)
runSum2 = np.zeros(firstvis.shape)
runNum = 0
assumedPixelAverage = np.mean(firstvis)		#guess average pixel value from first image
print "avg pixel value: ", assumedPixelAverage 

fileLim = None

for i,f in enumerate(filelist[:fileLim]):
	raw = rawpy.imread(rawdir+f)
	vis = raw.raw_image_visible
	adjVal = np.subtract(vis, assumedPixelAverage)

	runNum = i+1
	#runSum = np.add(runSum, np.subtract(vis, assumedPixelAverage))
	runSum = np.add(runSum, adjVal)
	runSum2 = np.add(runSum2, np.square(adjVal))

sigmaLim = np.sqrt(runNum) + 1.0
bbox = 10
sbox = range(-bbox,bbox+1)
marker = range(-bbox,bbox+1)[:5]+range(-bbox,bbox+1)[-5:]

mean = np.divide(runSum, runNum)
print "min mean = ", np.amin(mean)
print "max mean = ", np.amax(mean)
mean2 = np.divide(runSum2, runNum)
print "min mean2 = ", np.amin(mean2)
print "max mean2 = ", np.amax(mean2)
var = np.clip(np.subtract(mean2,np.square(mean)), 0.00001, None)
print "min var = ", np.amin(var)
print "max var = ", np.amax(var)
stddev = np.sqrt(var)
print "min std = ", np.amin(stddev)
print "max std = ", np.amax(stddev)
stdclip = np.clip(stddev, 0.00001, None)
print "min stdclip = ", np.amin(stdclip)
print "max stdclip = ", np.amax(stdclip)

with open("./5secV/mean.pickle","wb") as pickleFile:
	pickle.dump(mean, pickleFile, protocol=pickle.HIGHEST_PROTOCOL)
with open("./5secV/mean2.pickle","wb") as pickleFile:
	pickle.dump(mean2, pickleFile, protocol=pickle.HIGHEST_PROTOCOL)
# temporarily saved for later: how to read the pickeled data 
#	with open("mean.pickle","rb") as pickleFile:
#		mean = pickle.load(pickleFile)

for i,f in enumerate(filelist[:fileLim]):
	raw = rawpy.imread(rawdir+f)
	vis = np.subtract(raw.raw_image_visible, assumedPixelAverage)

	subAvgOnly = np.subtract(vis,mean)
	SAOmin = np.amin(subAvgOnly)
	SAOmax = np.amax(subAvgOnly)
	print "Image "+str(i)+" max pixel value: "+str(SAOmax)
	thresh = np.where(subAvgOnly>1020.0,1,0)
	for ii,a in enumerate(thresh):
		for jj,b in enumerate(a):
			if b==1 :
				print "Image "+str(i)+": super-bright pixel found! Indices: "+str(ii)+" "+str(jj)
				for iii in sbox:
					for jjj in sbox:
						if ii+iii>=0 and ii+iii<thresh.shape[0] and jj+jjj>=0 and jj+jjj<thresh.shape[1]:
							collector[ii+iii][jj+jjj] = subAvgOnly[ii+iii][jj+jjj]
							if (iii in marker) and (jjj in marker):
								subAvgOnly[ii+iii][jj+jjj] = 2000.0
								collector[ii+iii][jj+jjj] = 2000.0
	#SAOscaled = np.clip(np.multiply(np.add(SAOmin,subAvgOnly),255./(SAOmin+SAOmax)),0,255).astype('uint8')
	SAOscaled = np.clip(np.divide(subAvgOnly,4.0),0,255).astype('uint8')
	SAOim = Image.fromarray(SAOscaled)
	SAOim.save("./5secV/SAO"+str(i)+".png")
	continue

	darkSub = np.divide(np.subtract(vis,mean),stdclip)
	x = np.linspace(-sigmaLim,sigmaLim,100)

	test = darkSub.reshape(mean.shape[0]*mean.shape[1])
	print "max anomaly in image "+str(i)+": "+str(np.amax(test))

	#plt.hist(test,bins=x)
	hist, edges = np.histogram(test, x)
	centers = [(a+b)/2.0 for a,b in zip(edges[:-1], edges[1:])]
	plt.plot(centers, hist, 'ob')
	plt.semilogy()
	plt.savefig("./5secV/hist"+str(i)+".png")#plt.show()
	plt.gcf().clear()
	m = np.clip(np.multiply(np.add(darkSub,sigmaLim),255./(2.0*sigmaLim)),0,255).astype('uint8')
	im = Image.fromarray(np.clip(m,0,None))
	im.save("./5secV/test"+str(i)+".png")


collectorScaled = np.clip(np.divide(collector,4.0),0,255).astype('uint8')
collectorIM = Image.fromarray(collectorScaled)
collectorIM.save("./5secV/all_hits.png")

sys.exit(0)
##old


for f in filelist[:30]:
	rawlist.append(rawpy.imread(rawdir+f))
for f in rawlist:
	vislist.append(f.raw_image_visible)

#print vislist[0]

visStack = np.dstack(vislist)
#print visStack
#print visStack.shape
#print "dstack", visStack[0][0][0]

#median = np.median(visStack, axis=2).astype(int)
median = np.mean(visStack, axis=2)
stddev = np.std(visStack, axis=2)
stdclip = np.clip(stddev,1.0, None)
print np.amin(stdclip)
#print median.shape
#print median



#darkSubList = [np.subtract(a,median) for a in vislist]
darkSubList = [np.divide(np.subtract(a,median),stdclip) for a in vislist]
#print darkSubList[0]

#x = np.linspace(-400.,400.,100)
x = np.linspace(-10.,10.,100)
#print x


for i,l in enumerate(darkSubList):
	test = l.reshape(median.shape[0]*median.shape[1])
	print "max anomaly in image "+str(i)+": "+str(np.amax(l))

	#plt.hist(test,bins=x)
	hist, edges = np.histogram(test, x)
	centers = [(a+b)/2.0 for a,b in zip(edges[:-1], edges[1:])]
	plt.plot(centers, hist, 'ob')
	plt.semilogy()
	plt.savefig("./5secV/hist"+str(i)+".png")#plt.show()
	plt.gcf().clear()
	m = np.clip(np.multiply(np.add(l,7),255./14.),0,255).astype('uint8')
	im = Image.fromarray(np.clip(m,0,None))
	im.save("./5secV/test"+str(i)+".png")







