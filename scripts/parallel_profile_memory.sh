#!/bin/bash
input=$1
# Example input: band_600_3_0.008317,/home/zgh23/code/SparseWS/bench-hash
# Extract from input
# name=band_600_3_0.008317
# exec=/home/zgh23/code/SparseWS/bench-hash
name=$(echo $input | cut -d ',' -f 1)
exec=$(echo $input | cut -d ',' -f 2)
algname=$(echo $exec | cut -d '-' -f 4)

matrix_root=/scratch/zgh23/sparse_mat/
trans_root=/scratch/zgh23/sparse_mat_t/
# get the binary name from exec
# e.g. bench-coord-cf
exec_name=$(echo $exec | cut -d '/' -f 6)

profile_root=/home/zgh23/code/SparseWS/data/profile/CSC_CSR/cc/memory
profile_name=$algname-$name
input=$profile_root/xtmem-$profile_name.txt
output=/home/zgh23/code/SparseWS/data/profile/CSC_CSR/cc/$exec_name.csv

valgrind --tool=massif --xtree-memory=full --xtree-memory-file=$profile_root/$profile_name.kcg --massif-out-file=$profile_root/$profile_name.out $exec $matrix_root/$name.mtx $trans_root/$name-st.mtx
callgrind_annotate --auto=yes --inclusive=yes --sort=curB:100,curBk:100,totB:100,totBk:100,totFdB:100,totFdBk:100 $profile_root/$profile_name.kcg > $input
python /home/zgh23/code/SparseWS/scripts/extract_memory.py --matrix $name --input $input --output $output