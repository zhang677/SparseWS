name=$1
matrix_root=/home/eva_share/datasets/sparse_mat
trans_root=/home/eva_share/datasets/sparse_mat_t
output=/home/nfs_data/zhanggh/SparseWS/data/results/CSC_CSR_T/csc_csr_hash_flex_memory.csv
result_root=/home/nfs_data/zhanggh/SparseWS/data/hash
valgrind --tool=massif --max-snapshots=200 --time-unit=i --detailed-freq=1 --massif-out-file=${result_root}/massif_${name}.out /home/nfs_data/zhanggh/SparseWS/test-hash $matrix_root/$name/$name.mtx $trans_root/$name-st.mtx 1 0 0 $output
ms_print ${result_root}/massif_${name}.out > ${result_root}/massif_${name}.txt