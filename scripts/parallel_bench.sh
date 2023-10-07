#!/bin/bash
# exec=/home/zgh23/code/SparseWS/bench-coord-cf
exec=$1
name=$2
matrix_root=/home/zgh23/code/SparseWS/data/origin/
trans_root=/home/zgh23/code/SparseWS/data/shift/
# get the binary name from exec
# e.g. bench-coord-cf
exec_name=$(echo $exec | cut -d '/' -f 6)
output=/home/zgh23/code/SparseWS/data/ablation/CSC_CSR_T/all/$exec_name.csv
# name is the name of the matrix, e.g. powerlaw_5000_10_0.1.mtx
# Store the "powerlaw" part in the variable "type"
type=$(echo $name | cut -d '_' -f 1)
$exec $matrix_root/$type/$name.mtx $trans_root/$type/$name-st.mtx 10 0 0 $output

# /home/zgh23/spack/opt/spack/linux-ubuntu20.04-broadwell/gcc-10.5.0/parallel-20220522-cycypvpmliyzsy7pmdlgu7hqwjdz5om3/bin