#!/bin/bash
#   Post this job with `qsub -V $0'
#PBS -l nodes=1:ppn=1
#PBS -l walltime=100
#PBS -N TEST
#PBS -j oe
#PBS -W umask=02
echo "mpiexec is at: $(which mpiexec)"
eups list -s
echo $PATH
mpiexec --verbose python $HSCPIPE_DIR/bin/test.py foobar
