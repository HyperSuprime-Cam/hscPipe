FROM centos:7

MAINTAINER paprice
LABEL description="Hyper Suprime-Cam pipeline" \
      name="paprice/hscPipe"

ARG DEVTOOLSET=7
ARG INSTALL_PFS_ARGS
ARG INSTALL_DIR=/opt/hsc
ARG HSCPIPE_VERSION=7.0

# Set to something unique, like $(date +%s), to force rebuild
# https://stackoverflow.com/a/47017887
ARG FORCE=XXX

ENV OMP_NUM_THREADS=1
ENV SCONSFLAGS="-j 2"
ENV MAKEFLAGS="-j 1"

SHELL ["/bin/bash", "-c"]

USER root

# Upgrade compiler using devtoolset
RUN yum -y update \
 && yum install -y centos-release-scl \
 && yum install -y \
        "devtoolset-${DEVTOOLSET}-gcc" \
        "devtoolset-${DEVTOOLSET}-gcc-c++" \
        "devtoolset-${DEVTOOLSET}-gcc-gfortran" \
 && echo ". /opt/rh/devtoolset-${DEVTOOLSET}/enable" > "/etc/profile.d/devtoolset-${DEVTOOLSET}.sh" \
 && yum clean all -y

# LSST stack prerequisites
RUN yum -y update \
 && yum install -y bison curl curl-devel blas bzip2-devel bzip2 flex fontconfig \
        freetype-devel git libuuid-devel \
        libXext libXrender libXt-devel make openssl-devel patch perl \
        readline-devel tar zlib-devel ncurses-devel cmake glib2-devel \
        java-1.8.0-openjdk gettext perl-ExtUtils-MakeMaker \
        which \
 && yum clean all -y

# Conveniences
RUN yum install -y screen tree \
 && yum clean all -y

RUN useradd -m hsc
RUN mkdir -p ${INSTALL_DIR} && chown hsc:hsc ${INSTALL_DIR}

USER hsc
WORKDIR ${INSTALL_DIR}

RUN . /opt/rh/devtoolset-${DEVTOOLSET}/enable \
 && git clone https://github.com/HyperSuprime-Cam/lsstsw.git \
 && lsstsw/bin/deploy

# Strips out anything excess after the install
RUN . /opt/rh/devtoolset-${DEVTOOLSET}/enable \
 && . lsstsw/bin/setup.sh \
 && rebuild -u hscPipe \
 && ( find ${INSTALL_DIR} -exec strip --strip-unneeded --preserve-dates {} + \
      > /dev/null 2>&1 || true ) \
 && ( find ${INSTALL_DIR} -maxdepth 5 -name tests -type d -exec rm -rf {} + \
      > /dev/null 2>&1 || true ) \
 && ( find ${INSTALL_DIR} -maxdepth 5 -path "*doc/html" -type d -exec rm -rf {} + \
      > /dev/null 2>&1 || true ) \
 && ( find ${INSTALL_DIR} -maxdepth 5 -name src -type d -exec rm -rf {} + \
      > /dev/null 2>&1 || true )

WORKDIR /home/hsc
RUN echo $'. '${INSTALL_DIR}'/loadLSST.bash\n\
setup hscPipe\n' >> /home/hsc/.bashrc

CMD /bin/bash
