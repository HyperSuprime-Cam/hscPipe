#!/bin/bash

# Configuration parameters
HERE=$(pwd)  # Starting directory (where the package is checked out)
WORKDIR=/tigress/pprice/hscPipe8  # Working directory
DISTRIB_SRC=$WORKDIR/distrib/src  # Distribution directory for sources
DISTRIB_BIN=$WORKDIR/distrib/Linux64  # Distribution directory for binary tarballs
export SCONSFLAGS="-j 4"  # SCons build flags
export EUPSPKG_NJOBS=4  # Number of build jobs

# Need these on tiger to get the right environment
. /etc/profile  # Get "module"
module load rh/devtoolset/8  # Get modern compiler
module load git  # Get git-lfs

set -ev

# Set parameters from Jenkins envvars
GIT_TAG=$(git describe --tags --always)  # Tag to build
VERSION=$(echo "$GIT_TAG" | sed 's|[/ ]|_|g')  # Version to call it
env

# Build the stack
pushd $WORKDIR
. $WORKDIR/lsstsw/bin/setup.sh
rebuild -u -r $VERSION -r ${VERSION}-hsc -t current hscPipe

# Build the source distribution
eups distrib create --server-dir=$DISTRIB_SRC -S REPOSITORY_PATH='git://github.com/HyperSuprime-Cam/$PRODUCT.git' -f generic -d eupspkg hscPipe $VERSION
eups distrib declare --server-dir=$DISTRIB_SRC -t current

# Build the binary distribution
eups distrib create --server-dir=$DISTRIB_BIN -d tarball hscPipe $VERSION
eups distrib declare --server-dir=$DISTRIB_BIN -t current

# Clean up
popd
