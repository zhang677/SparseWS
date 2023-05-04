name=$1
matrix_root=/home/eva_share/datasets/sparse_mat
trans_root=/home/eva_share/datasets/sparse_mat_t
output=/home/nfs_data/zhanggh/SparseWS/data/test.csv
/home/nfs_data/zhanggh/SparseWS/test-hash $matrix_root/$name/$name.mtx $trans_root/$name-st.mtx 10 0 0 $output