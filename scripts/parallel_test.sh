name=$1
matrix_root=/home/eva_share/datasets/sparse_mat
trans_root=/home/eva_share/datasets/sparse_mat_t
output=/home/nfs_data/zhanggh/SparseWS/data/results/CSC_CSR_T/test-coord-bucket.csv
/home/nfs_data/zhanggh/SparseWS/test-coord $matrix_root/$name/$name.mtx $trans_root/$name-st.mtx 10 0 0 $output