#!/bin/bash
#   Post this job with `qsub -V $0'
#PBS -l nodes=2:ppn=3
#PBS -l walltime=6000
#PBS -N TEST
#PBS -j oe
#PBS -W umask=02
echo "mpiexec is at: $(which mpiexec)"
eups list -s
echo $PATH
mpiexec --verbose python $HSCPIPE_DIR/bin/test.py suprimecam price 126969
