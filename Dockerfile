FROM ubuntu:latest
MAINTAINER Caoxiang Zhu <czhu@pppl.gov> & STELLOPT developers

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get -q update && apt-get -y install \
    gfortran g++ libopenmpi-dev openmpi-bin \
    libnetcdf-dev libnetcdff-dev libhdf5-openmpi-dev hdf5-tools \
    libblas-dev liblapack-dev libscalapack-openmpi-dev \
    python3 python3-numpy python3-h5py make curl \
    git-all


# Build and install the following libraries with OpenMP support:
# Adapted from https://hub.docker.com/r/berkeleygw/bgw-docker
# 1. OpenBLAS
# 2. FFTW
# 3. LAPACK
# 4. ScaLAPCK

WORKDIR /home/STELLOPT

# ENV OB_VER=0.2.19
# ENV BLAS=OpenBLAS-${OB_VER}
# ENV OB_FLAGS="CC=gcc FC=gfortran NO_SHARED=1 NO_CBLAS=1 NO_LAPACK=1"
# RUN mkdir ${BLAS}
# RUN /bin/bash -l -c '\
#     curl -sSL "http://github.com/xianyi/OpenBLAS/archive/v${OB_VER}.tar.gz" | \
#     tar xz && cd ${BLAS} && \
#     ( make ${OB_FLAGS} USE_OPENMP=1 USE_THREADS=1 MAKE_NB_JOBS=8 || \
#       make ${OB_FLAGS} USE_OPENMP=0 USE_THREADS=0 MAKE_NB_JOBS=1 ) && \
#     make ${OB_FLAGS} PREFIX=/opt/${BLAS} install && \
#     cd ../ && rm -rf ${BLAS}'

# ENV FFTW=fftw-3.3.5
# RUN mkdir ${FFTW}
# RUN /bin/bash -l -c '\
#     curl -sSL "ftp://ftp.fftw.org/pub/fftw/${FFTW}.tar.gz" | \
#     tar xz && cd ${FFTW} && \
#     ./configure --enable-openmp --prefix=/opt/${FFTW} F77=gfortran F90=gfortran && \
#     make -j 8 && make install && \
#     cd ../ && rm -rf ${FFTW}'

# ENV LAPACK=lapack-3.6.1
# RUN mkdir ${LAPACK}
# COPY make.inc ./${LAPACK}/
# RUN /bin/bash -l -c '\
#     curl -sSL "http://www.netlib.org/lapack/${LAPACK}.tgz" | \
#     tar xz && cd ${LAPACK} && \
#     make lapacklib -j 8 && cp liblapack.a /lib64 && \
#     cd ../ && rm -rf ${LAPACK}'

# ENV SCALAPACK=scalapack-2.0.2
# RUN mkdir ${SCALAPACK}
# COPY SLmake.inc ./${SCALAPACK}/
# RUN /bin/bash -l -c '\
#     curl -sSL "http://www.netlib.org/scalapack/${SCALAPACK}.tgz" | \
#     tar xz && cd ${SCALAPACK} && \
#     make lib && cp libscalapack.a /lib64 && \
#     cd ../ && rm -rf ${SCALAPACK}'

# Get source code
#RUN ls -al
#RUN git clone https://github.com/PrincetonUniversity/STELLOPT.git
COPY . /home/STELLOPT
# Compile STELLOPT
ENV MACHINE="docker"
ENV STELLOPT_PATH=/home/STELLOPT
RUN echo $STELLOPT_PATH
RUN cd $STELLOPT_PATH  && ./build_all

# Add paths
ENV PATH="${STELLOPT_PATH}/bin:${PATH}"
CMD ["/bin/bash"]
