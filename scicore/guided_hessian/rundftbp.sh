#!/bin/bash

export OMP_NUM_THREADS=1
## path to dftb+ executeble. The first path is on scicore.
/kernph/hubhan00/dftbplus-20.2.1/bin/dftb+ &> dftbp.out
#/home/hannes/miniconda3/bin/dftb+ &> dftbp.out
#/home/hannes/Downloads/dftbplus-20.1/bin/dftb+ &> dftbp.out
