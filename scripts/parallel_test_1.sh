name=$1
exec=/home/zgh23/code/SparseWS/bench-gust-taco
matrix_root=/scratch/zgh23/sparse_mat
trans_root=/scratch/zgh23/sparse_mat_t
output=/home/zgh23/code/SparseWS/data/results/CSR_CSR/bench-gust-taco.csv
$exec $matrix_root/$name.mtx $trans_root/$name-st.mtx 20 0 0 $output