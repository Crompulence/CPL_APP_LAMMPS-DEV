LAMMPS_DIR=`cat CODE_INST_DIR`
LAMMPS_SRC_DIR=$(LAMMPS_DIR)/src
LAMMPSVERSION=`cat $(LAMMPS_SRC_DIR)/version.h`

.PHONY: all test clean clean-all yes-all patch-lammps unpatch-lammps

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
	python ./config/patch_main.py $(LAMMPS_SRC_DIR)

unpatch-lammps:
	cp ./config/main_prepatch.cpp $(LAMMPS_SRC_DIR)/main.cpp

patch-lammps-old:
	python ./config/get_patch.py $(LAMMPS_SRC_DIR)
	cp ./config/mpmd.patch $(LAMMPS_DIR)
	cd $(LAMMPS_DIR) && patch -N -p1 < mpmd.patch

unpatch-lammps-old:
	python ./config/get_patch.py $(LAMMPS_SRC_DIR)
	cp ./config/mpmd.patch $(LAMMPS_DIR)
	cd $(LAMMPS_DIR) && patch -N -R -p1 < mpmd.patch

patch-lammps-Oct2017:
	cp ./config/mpmd_Oct2017.patch $(LAMMPS_DIR)
	cd $(LAMMPS_DIR) && patch -N -p1 < mpmd_Oct2017.patch

patch-lammps-Apr2018:
	cp ./config/mpmd_Apr2018.patch $(LAMMPS_DIR)
	cd $(LAMMPS_DIR) && patch -p2 < mpmd_Apr2018.patch

unpatch-lammps-Aug2017:
	cp ./config/mpmd_Oct2017.patch $(LAMMPS_DIR)
	cd $(LAMMPS_DIR) && patch -R -p1 < mpmd_Aug2017.patch

unpatch-lammps-Oct2017:
	cp ./config/mpmd_Oct2017.patch $(LAMMPS_DIR)
	cd $(LAMMPS_DIR) && patch -R -p1 < mpmd_Oct2017.patch

unpatch-lammps-Apr2018:
	cp ./config/mpmd_Apr2018.patch $(LAMMPS_DIR)
	cd $(LAMMPS_DIR) && patch -R -p1 < mpmd_Apr2018.patch

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
	pytest -v ./test

test-single:
	pytest -vs ./test/constant_force

test-simwrap:
	pytest -vs ./test/Example_simwrap/

test-couette:
	pytest -vs ./test/Couette_coupled/Partial_overlap/

test-granular:
	cd test/granular/single/terminal_velocity && pytest -v ./test_terminal.py
	cd test/granular/single/constant_velocity && pytest -v ./test_constant.py
	cd test/granular/single/resting_wall && pytest -v ./test_resting.py
	cd test/granular/column && pytest -v ./test_column.py
	cd test/granular/suzuki && pytest -v ./test_suzuki.py
