#!/bin/bash
# INSTALLATION TOOL for LABEL.
# Sam Shepard, 2012

bin=binaries_and_licenses
shogun=shogun1.1.0_cmdline_static
sam=sam3.5/bin
muscle=muscle3.8.31
fasttree=fasttree2.1.4

ln -s $bin/$fasttree/FastTreeMP_darwin64 FastTreeMP_Darwin
ln -s $bin/$muscle/muscle3.8.31_darwin64 muscle_Darwin
ln -s $bin/$sam/hmmscore_darwin64 hmmscore_Darwin
ln -s $bin/$sam/modelfromalign_darwin64 modelfromalign_Darwin
ln -s $bin/$shogun/shogun_darwin64 shogun_Darwin
ln -s $bin/$sam/align2model_darwin64 align2model_Darwin

ln -s $bin/$fasttree/FastTreeMP_linux64 FastTreeMP_Linux
ln -s $bin/$muscle/muscle3.8.31_linux64 muscle_Linux
ln -s $bin/$sam/hmmscore_linux64 hmmscore_Linux
ln -s $bin/$sam/modelfromalign_linux64 modelfromalign_Linux
ln -s $bin/$shogun/shogun_linux64 shogun_Linux
ln -s $bin/$sam/align2model_linux64 align2model_Linux

ln -s $bin/$fasttree/FastTreeMP_linux32 FastTreeMP_Linux32
ln -s $bin/$muscle/muscle3.8.31_linux32 muscle_Linux32
ln -s $bin/$sam/hmmscore_linux32 hmmscore_Linux32
ln -s $bin/$sam/modelfromalign_linux32 modelfromalign_Linux32
ln -s $bin/$shogun/shogun_linux32 shogun_Linux32
ln -s $bin/$sam/align2model_linux32 align2model_Linux32
