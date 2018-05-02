MAKE=$1
APP_DIR=$(pwd)
LAMMPS_SRC_DIR=$(cat ../CODE_INST_DIR)/src
cd $LAMMPS_SRC_DIR
$MAKE "no-all"
$MAKE ps
while read p; do
  $MAKE "yes-$p"
done <$APP_DIR/lammps_packages.in
$MAKE ps
