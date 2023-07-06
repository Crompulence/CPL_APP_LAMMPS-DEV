CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
export PATH=$CWD/bin:$PATH

if [ -f CODE_INST_DIR ]; then
	export LAMMPS_PATH=$(cat $CWD/CODE_INST_DIR)
else
	printf "ERROR:\n\tFile CODE_INST_DIR not found. Please provide one as described in the README.md file.\n"
	return 1
fi

