LAMMPS_SRC_DIR=`cat CODE_INST_DIR`/src

.PHONY: all test clean yes-all

all:
	cp -R src/USER-CPL $(LAMMPS_SRC_DIR)
	cp ./config/Makefile.cpl $(LAMMPS_SRC_DIR)/MAKE
	cd $(LAMMPS_SRC_DIR) && $(MAKE) yes-USER-CPL
	cd $(LAMMPS_SRC_DIR) && $(MAKE) cpl
	cp $(LAMMPS_SRC_DIR)/lmp_cpl ./bin

yes-all:
	bash config/enable-packages.sh $(MAKE)

clean:
	cd $(LAMMPS_SRC_DIR) && $(MAKE) clean-all

test:
	py.test -v ./test
