#!/bin/bash
matrix_root=/scratch/zgh23/sparse_mat
matrix_t_root=/scratch/zgh23/sparse_mat_t
matrix_name=$1
profile_root=/home/zgh23/code/SparseWS/data/profile/memory
profile_name=hash_$matrix_name
exec=/home/zgh23/code/SparseWS/memory-hash
matrixA=$matrix_root/$matrix_name.mtx
matrixB=$matrix_t_root/$matrix_name-st.mtx
valgrind --tool=massif --xtree-memory=full --xtree-memory-file=$profile_root/$profile_name.kcg --massif-out-file=$profile_root/$profile_name.out  $exec $matrixA $matrixB
# ms_print $profile_root/$profile_name.out > $profile_root/$profile_name.txt
callgrind_annotate --auto=yes --inclusive=yes --sort=curB:100,curBk:100,totB:100,totBk:100,totFdB:100,totFdBk:100 $profile_root/$profile_name.kcg > $profile_root/xtmem_$profile_name.txt
