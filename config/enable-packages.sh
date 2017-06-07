MAKE=$1
APP_DIR=$(pwd)
LAMMPS_SRC_DIR=$(cat CODE_INST_DIR)/src
cd $LAMMPS_SRC_DIR
while read p; do
  $MAKE "no-$p"
  $MAKE "yes-$p"
done <$APP_DIR/config/lammps_packages.in
