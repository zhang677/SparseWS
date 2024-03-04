name=$1
exec=/home/zgh23/code/SparseWS/count-outer-cc-noT-FLOPS
matrix_root=/scratch/zgh23/sparse_mat
trans_root=/scratch/zgh23/sparse_mat_t
output=/home/zgh23/code/SparseWS/data/profile/CSC_CSR/cc/count-cc-noT-FLOPS.csv
$exec $matrix_root/$name.mtx $trans_root/$name-st.mtx $output