#!/usr/bin/env python
# Samuel Shepard - doSVM - 9.2011
# Tests unknown sequences for the genotyping system.
# The load_matrix function is given below and taken
# from the SHOGUN toolkit help files.
from numpy import double, fromfile
from sys import argv, exit
from os.path import basename
from shogun.Features import RealFeatures, Labels
from shogun.Kernel import PolyKernel
from shogun.Classifier import GMNPSVM
from shogun.IO import SerializableAsciiFile

# Courtesy the Shogun Toolkit 1.0.0
# Soeren Sonnenburg et al. 2010
# http://www.shogun-toolbox.org/
# load.py in examples/undocumented/python_static/tools
# modified for arbritrary lines in the file
def load_numbers(filename):
	matrix=fromfile(filename, sep=' ')
	# whole matrix is 1-dim now, so reshape
	L=len(open(filename).readlines())
	matrix=matrix.reshape(L, len(matrix)/L)
	return matrix

def load_labels(filename):
	return fromfile(filename, dtype=double, sep=' ')

# ARG LIST:
# doSVM.py GROUP TEST RESULT TRAIN LABEL
# doSVM.py GROUP TEST RESULT TRAIN LABEL CLASSIFIER


# Instantiate variables.
size_cache=20
degree=20
inhomogene=True
use_normalization=True
alist=argv

# Load Data
testdata=load_numbers( alist[2] )
traindata=load_numbers( alist[4] )
labeldata=load_labels( alist[5] )

feats_test=RealFeatures(testdata)
feats_train=RealFeatures(traindata)
labels=Labels(labeldata)
kernel=PolyKernel(size_cache, degree, inhomogene)

# Load or train SVM
if ( len(alist) == 6 ):
	svm=GMNPSVM(1, kernel, labels)
	svm.train(feats_train)
	#######
	saveFile= "%s_classifier.dat" % (alist[1])
	saveFile=basename(saveFile)
	fstream = SerializableAsciiFile(saveFile, "w")
	status = svm.save_serializable(fstream)
	#######
else:
	svm=GMNPSVM(1, kernel, labels)
	fstream = SerializableAsciiFile( alist[6], "r")
	status = svm.load_serializable(fstream)
	svm.train()

# Classify
kernel.init(feats_train, feats_test)
R=svm.apply(feats_test).get_labels()

# Output results.
fp=open( alist[3], "a")
for i in range(len(R)):
	fp.write("%d\t%d\t%s\n" % ((i+1), R[i], alist[1]) )
fp.close()
