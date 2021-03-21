#!/usr/bin/env python
# Samuel Shepard - doSVM - 9.2011
# Tests unknown sequences for the genotyping system.
# The load_matrix function is given below and taken
# from the SHOGUN toolkit help files.
from numpy import double, fromfile
from sg import sg
from sys import argv, exit
from os.path import basename

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
	print matrix
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
sg('new_classifier', 'GMNPSVM')
sg('set_kernel', 'POLY', 'REAL', size_cache, degree, inhomogene, use_normalization)
traindata=load_numbers( alist[4] )
labels=load_labels( alist[5] )
sg('set_features', 'TRAIN', traindata)
sg('set_labels', 'TRAIN', labels)

# Load or train SVM
if ( len(alist) == 6 ):
	sg('train_classifier')
	#######
	saveFile= "%s_classifier.dat" % (alist[1])
	saveFile=basename(saveFile)
	sg('save_classifier',saveFile)
	#######
else:
	sg('load_classifier', alist[6], 'GMNPSVM' )


# Perform classification
testdata=load_numbers( alist[2] )
sg('set_features', 'TEST', testdata)
R=sg('classify')
#kernel_matrix = sg('get_kernel_matrix', 'TEST')

# Output results.
fp=open( alist[3], "a")
for i in range(len(R)):
	fp.write("%d\t%d\t%s\n" % ((i+1), R[i], alist[1]) )
fp.close()
