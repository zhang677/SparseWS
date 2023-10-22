name=$1
exec=/home/zgh23/code/SparseWS/bench-nofuse-hash
matrix_root=/scratch/zgh23/sparse_mat
trans_root=/scratch/zgh23/sparse_mat_t
output=/home/zgh23/code/SparseWS/data/ablation/CSC_CSR_T/nofuse/bench-nofuse-hash.csv
# $exec -matrixA $matrix_root/$name.mtx -matrixB $trans_root/$name-st.mtx -iter 20 -result $output
$exec $matrix_root/$name.mtx $trans_root/$name-st.mtx 20 0 0 $output

# /home/zgh23/code/SparseWS/test-profile /scratch/zgh23/sparse_mat/CAG_mat1916.mtx /scratch/zgh23/sparse_mat_t/CAG_mat1916-st.mtx 10 0 0 /home/zgh23/code/SparseWS/data/profile/CSC_CSR_T/profile-coord-bucket
# /home/zgh23/code/SparseWS/test-profile /scratch/zgh23/sparse_mat/M80PI_n.mtx /scratch/zgh23/sparse_mat_t/M80PI_n-st.mtx 5 0 0 /home/zgh23/code/SparseWS/data/profile/CSC_CSR_T/M80PI_n-coord-bucket