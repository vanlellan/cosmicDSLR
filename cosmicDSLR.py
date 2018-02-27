#
# ToDo
#DONE	-pickle sum, sum2, and num instead of mean and mean2 (and save as just one picklejar)
#DONE	-store everything in source directory by default (create new cosmic-search directory within source dir), with option to store in a different directory
#DONE	-add command-line options
#DONE	-make this code a function, add wrapper script for live-updating search during exposure session
#DONE	-save name of each added image to list, check against list before each add (save this list to the picklejar as well)
#	-add option to delete raw images after import (useful for a long-running live display, to save disk space)


import rawpy
import os, sys
import numpy as np
from matplotlib import pyplot as plt
from PIL import Image
import pickle
import subprocess
import time
import signal
from optparse import OptionParser

def process(rawdir, fileLim, batchMode):
	
	
	print np.__version__
	
	print "Importing from:", rawdir
	
	filelist = [a for a in os.listdir(rawdir) if a[-4:]==".cr2"]
	print "Trying "+str(len(filelist[:fileLim]))+" of "+str(len(filelist))+" found files."
	
	#check for existence of cosmic pickle and load or create data storage variables
	try:
		with open(rawdir+"cosmic.pickle","rb") as pickleFile:
			runNum, runSum, runSum2, assumedPixelAverage, doneFileList, collector = pickle.load(pickleFile)
		print "Loaded Data from "+str(len(doneFileList))+" Files Successfully!"
	except:
		print "Load Unsuccessful, initializing data storage"
		firstvis = rawpy.imread(rawdir+filelist[0]).raw_image_visible	#fails here if directory is empty
		collector = np.zeros(firstvis.shape)
		runSum = np.zeros(firstvis.shape)
		runSum2 = np.zeros(firstvis.shape)
		runNum = 0
		doneFileList = []
		assumedPixelAverage = np.mean(firstvis)		#guess average pixel value from first image
	
	newFileList = []
	for i,f in enumerate(filelist[:fileLim]):
		if f not in doneFileList:
			print "Adding file "+str(len(doneFileList))+" to picklejar." 
			raw = rawpy.imread(rawdir+f)
			vis = raw.raw_image_visible
			adjVal = np.subtract(vis, assumedPixelAverage)
			runNum = i+1
			runSum = np.add(runSum, adjVal)
			runSum2 = np.add(runSum2, np.square(adjVal))
			doneFileList.append(f)
			newFileList.append(f)
	
	sigmaLim = np.sqrt(runNum) + 1.0
	bbox = 10
	sbox = range(-bbox,bbox+1)
	marker = range(-bbox,bbox+1)[:5]+range(-bbox,bbox+1)[-5:]
	
	mean = np.divide(runSum, runNum)
	mean2 = np.divide(runSum2, runNum)
	var = np.clip(np.subtract(mean2,np.square(mean)), 0.00001, None)
	#stddev = np.sqrt(var)
	#stdclip = np.clip(stddev, 0.00001, None)
	
	displayList = []
	for i,f in enumerate(doneFileList):
		if batchMode or f in newFileList:
			cosmicFlag = False
			print "Analyzing file "+str(i)
			raw = rawpy.imread(rawdir+f)
			vis = np.subtract(raw.raw_image_visible, assumedPixelAverage)
	
			subAvgOnly = np.subtract(vis,mean)
			SAOmin = np.amin(subAvgOnly)
			SAOmax = np.amax(subAvgOnly)
			print "Image "+str(i)+" max pixel value: "+str(SAOmax)
			thresh = np.where(subAvgOnly>1020.0,1,0)	#change threshold to scale with average pixel brightness?
			for ii,a in enumerate(thresh):
				for jj,b in enumerate(a):
					if b==1 :
						print "Image "+str(i)+": cosmic candidate found! Indices: "+str(ii)+" "+str(jj)
						cosmicFlag = True
						for iii in sbox:
							for jjj in sbox:
								if ii+iii>=0 and ii+iii<thresh.shape[0] and jj+jjj>=0 and jj+jjj<thresh.shape[1]:
									collector[ii+iii][jj+jjj] = subAvgOnly[ii+iii][jj+jjj]
									if (iii in marker) and (jjj in marker):
										subAvgOnly[ii+iii][jj+jjj] = 2000.0
										collector[ii+iii][jj+jjj] = 2000.0
			SAOscaled = np.clip(np.divide(subAvgOnly,4.0),0,255).astype('uint8')	#change denominator to scale with average brightness
			SAOim = Image.fromarray(SAOscaled)
			saveFileName = rawdir+"SAO"+str(i)+".png"
			SAOim.save(saveFileName)
			if cosmicFlag and not batchMode:
				displayList.append(subprocess.Popen(["eog",saveFileName]))
	
	#	darkSub = np.divide(np.subtract(vis,mean),stdclip)
	#	x = np.linspace(-sigmaLim,sigmaLim,100)
	#
	#	test = darkSub.reshape(mean.shape[0]*mean.shape[1])
	#	print "max anomaly in image "+str(i)+": "+str(np.amax(test))
	#
	#	#plt.hist(test,bins=x)
	#	hist, edges = np.histogram(test, x)
	#	centers = [(a+b)/2.0 for a,b in zip(edges[:-1], edges[1:])]
	#	plt.plot(centers, hist, 'ob')
	#	plt.semilogy()
	#	plt.savefig(rawdir+"hist"+str(i)+".png")#plt.show()
	#	plt.gcf().clear()
	#	m = np.clip(np.multiply(np.add(darkSub,sigmaLim),255./(2.0*sigmaLim)),0,255).astype('uint8')
	#	im = Image.fromarray(np.clip(m,0,None))
	#	im.save(rawdir+"test"+str(i)+".png")
	
	with open(rawdir+"cosmic.pickle","wb") as pickleFile:
		pickle.dump([runNum, runSum, runSum2, assumedPixelAverage, doneFileList, collector], pickleFile, protocol=pickle.HIGHEST_PROTOCOL)
	
	collectorScaled = np.clip(np.divide(collector,4.0),0,255).astype('uint8')
	collectorIM = Image.fromarray(collectorScaled)
	collectorIM.save(rawdir+"all_hits.png")

if __name__ == "__main__":

        ## Define the command line input
        parser = OptionParser('Usage: %prog [options] raw-dir-path')
        parser.add_option('-l', '--file-lim', type='int', dest = 'fileLim', help = 'Limit on number of files to read in.', default = None, metavar="LIMIT")
        parser.add_option('-b', '--batch', action="store_true", dest = 'batch', help = 'Enable batch mode.', default = False)
        (options, args) = parser.parse_args()


	while True:
		rawdir = sys.args[0]	#possibly make this handle a list of input directories
		fileLim = options.fileLim
		batchMode = options.batch
		process(rawdir, fileLim, batchMode)
		if batchMode:
			print "Batch Analysis Completed. Exiting."
			sys.exit(0)
		try:
			print "Sleeping 10 seconds..."
			time.sleep(10)
		except(KeyboardInterrupt):
			print "Ctrl-C Caught... Aborting."
			sys.exit(0)
