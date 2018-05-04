LAMMPS_DIR=`cat CODE_INST_DIR`
LAMMPS_SRC_DIR=$(LAMMPS_DIR)/src
LAMMPSVERSION=`cat $(LAMMPS_SRC_DIR)/version.h`

.PHONY: all test clean clean-all yes-all patch-lammps

all:
	cp -R src/USER-CPL $(LAMMPS_SRC_DIR)
	cp ./config/Makefile.cpl $(LAMMPS_SRC_DIR)/MAKE
	cd $(LAMMPS_SRC_DIR) && $(MAKE) no-USER-CPL && $(MAKE) yes-USER-CPL
	#cd $(LAMMPS_SRC_DIR) && $(MAKE) yes-granular
	cd $(LAMMPS_SRC_DIR) && $(MAKE) cpl
	rm -rf bin > /dev/null
	mkdir bin
	cp -f $(LAMMPS_SRC_DIR)/lmp_cpl ./bin

patch-lammps:
	echo "Specify patch-lammps-Aug17 or patch-lammps-Oct17"
	echo $(LAMMPSVERSION)

patch-lammps-Aug17:
	cp ./config/mpmd_Aug2017.patch $(LAMMPS_DIR)
	cd $(LAMMPS_DIR) && patch -N -p1 < mpmd_Aug2017.patch

patch-lammps-Oct17:
	cp ./config/mpmd_Oct2017.patch $(LAMMPS_DIR)
	cd $(LAMMPS_DIR) && patch -N -p1 < mpmd_Oct2017.patch

unpatch-lammps-Aug17:
	cp ./config/mpmd_Oct2017.patch $(LAMMPS_DIR)
	cd $(LAMMPS_DIR) && patch -R -p1 < mpmd_Aug2017.patch

unpatch-lammps-Oct17:
	cp ./config/mpmd_Oct2017.patch $(LAMMPS_DIR)
	cd $(LAMMPS_DIR) && patch -R -p1 < mpmd_Oct2017.patch

yes-all:
	bash config/enable-packages.sh $(MAKE)

clean-tests:
	cd test/forceC-P/debug && ./clean.sh
	cd test/velocityP-C/debug && ./clean.sh

clean: clean-tests
	rm -rf bin

clean-all: clean
	cd $(LAMMPS_SRC_DIR) && $(MAKE) clean-all

test:
	py.test -v ./test

test-single:
	py.test -v -s ./test/constant_force

test-simwrap:
	cd ./test/Example_simwrap
	py.test -v -s ./
