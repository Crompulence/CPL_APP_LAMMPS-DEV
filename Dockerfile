# start from cpl library base
FROM cpllibrary/cpl-library
MAINTAINER Edward Smith <edward.smith05@imperial.ac.uk>

#Number of processes to use in build
ENV NPROCS=1

#Get LAMMPS
RUN git clone https://github.com/lammps/lammps.git /lammps &&  \
    git clone https://github.com/Crompulence/CPL_APP_LAMMPS-DEV.git /CPL_APP_LAMMPS-DEV

#Build LAMMPS with USER-CPL package from APP 
WORKDIR /CPL_APP_LAMMPS-DEV
RUN echo "/lammps" > /CPL_APP_LAMMPS-DEV/CODE_INST_DIR && \
    echo granular >> config/lammps_packages.in && \
    cd config && \
    sh ./enable-packages.sh make && \
    cd ../ && \
    make patch-lammps

RUN make -j $NPROCS

ENV PATH="/CPL_APP_LAMMPS-DEV/bin:${PATH}"
