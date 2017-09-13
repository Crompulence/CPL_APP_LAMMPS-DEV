LAMMPS_DIR=`cat CODE_INST_DIR`
LAMMPS_SRC_DIR=$(LAMMPS_DIR)/src

.PHONY: all test clean clean-all yes-all patch-lammps

all:
	cp -R src/USER-CPL $(LAMMPS_SRC_DIR)
	cp ./config/Makefile.cpl $(LAMMPS_SRC_DIR)/MAKE
	cd $(LAMMPS_SRC_DIR) && $(MAKE) no-USER-CPL && $(MAKE) yes-USER-CPL
	cd $(LAMMPS_SRC_DIR) && $(MAKE) cpl
	rm -rf bin > /dev/null
	mkdir bin
	cp -f $(LAMMPS_SRC_DIR)/lmp_cpl ./bin

patch-lammps:
	cp ./config/lammps_cpl.patch $(LAMMPS_DIR)
	cd $(LAMMPS_DIR) && patch -p1 < lammps_cpl.patch

yes-all:
	bash config/enable-packages.sh $(MAKE)

clean: clean-tests
	rm -rf bin

clean-tests:
	cd test/forceC-P/debug && ./clean.sh
	cd test/velocityP-C/debug && ./clean.sh

clean-all: clean
	cd $(LAMMPS_SRC_DIR) && $(MAKE) clean-all

test:
	py.test -v ./test
