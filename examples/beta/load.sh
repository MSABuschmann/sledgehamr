#!/bin/bash

ml gpu
ml pytorch
ml cray-mpich
ml cudatoolkit/11.7
ml cray-fftw
ml cray-hdf5-parallel
export LIBRARY_PATH=$LD_LIBRARY_PATH
