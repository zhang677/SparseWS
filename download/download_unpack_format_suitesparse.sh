#!/bin/bash
# Command: export SUITESPARSE_PATH=/scratch/zgh23/sparse_mat i.e. <path_to_suitesparse>
# Command: ./scripts/formatting/download_unpack_format_suitesparse.sh <tensor_names.txt>

basedir=$(pwd)
download_script=download_suitesparse_partial.sh

[ -e $download_script ] && rm $download_script
echo "mkdir -p ${SUITESPARSE_PATH}" >> $download_script
echo "pushd ." >> $download_script
echo "cd ${SUITESPARSE_PATH}" >> $download_script
grep -F -f $1 download_suitesparse.sh >> $download_script 
echo "popd" >> $download_script

# Make it an executable
chmod ugo+x $download_script

# Call the download_script (created above)
./$download_script

# Unpack the downloaded suitesparse files since they come in .tar format
./unpack_suitesparse.sh $(realpath $1)