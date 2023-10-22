#!/bin/bash
input=$1
# Example input: band_600_3_0.008317,/home/zgh23/code/SparseWS/bench-hash
# Extract from input
# name=band_600_3_0.008317
# exec=/home/zgh23/code/SparseWS/bench-hash
name=$(echo $input | cut -d ',' -f 1)
exec=$(echo $input | cut -d ',' -f 2)

matrix_root=/home/zgh23/code/SparseWS/data/origin/
trans_root=/home/zgh23/code/SparseWS/data/shift/
# get the binary name from exec
# e.g. bench-coord-cf
exec_name=$(echo $exec | cut -d '/' -f 6)
# output=/home/zgh23/code/SparseWS/data/ablation/CSC_CSR_T/nofuse/$exec_name.csv
# output=/home/zgh23/code/SparseWS/data/ablation/CSC_CSR/cc/$exec_name.csv
output=/home/zgh23/code/SparseWS/data/ablation/CSR_CSR/shift/$exec_name.csv
# name is the name of the matrix, e.g. powerlaw_5000_10_0.1.mtx
# Store the "powerlaw" part in the variable "type"
type=$(echo $name | cut -d '_' -f 1)
$exec $matrix_root/$type/$name.mtx $trans_root/$type/$name-st.mtx 10 0 0 $output