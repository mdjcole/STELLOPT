FROM ubuntu:latest
MAINTAINER Caoxiang Zhu <czhu@pppl.gov> & STELLOPT developers

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get -q update && apt-get -y install \
    gfortran g++ libopenmpi-dev openmpi-bin \
    libnetcdf-dev libnetcdff-dev libhdf5-openmpi-dev hdf5-tools \
    libblas-dev liblapack-dev libscalapack-openmpi-dev \
    python3 python3-numpy python3-h5py make curl \
    git-all

WORKDIR /home/STELLOPT

COPY . /home/STELLOPT
# Compile STELLOPT
ENV MACHINE="docker"
ENV STELLOPT_PATH=/home/STELLOPT
RUN echo $STELLOPT_PATH
RUN cd $STELLOPT_PATH  && ./build_all 2>&1 | tee log.build 
RUN chmod -R 755 ${STELLOPT_PATH}/BENCHMARKS
RUN cp -RP ${STELLOPT_PATH}/bin/* /usr/local/bin/

# add user
WORKDIR /home/visitor
RUN groupadd -r visitor -g 433 && \
useradd -u 431 -r -g visitor -d /home/visitor -s /sbin/nologin -c "Docker image user" visitor && \
chown -R visitor:visitor /home/visitor
USER visitor

# Set commands
CMD ["/bin/bash"]