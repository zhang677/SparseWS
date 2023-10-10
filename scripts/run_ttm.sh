exec=/home/zgh23/code/SparseWS/test-ttm
tensor=/scratch/zgh23/sparse_ten/nell-2.tns
# tensor=/home/zgh23/code/SparseWS/data/origin/5-5-5-9.tns
Cname=/home/zgh23/code/ctf/C-nell-2-16.mtx
#Cname=/home/zgh23/code/SparseWS/data/origin/rec_2_5_2_6_0.1333.mtx
Aname=/home/zgh23/code/ctf/ttm-A-16.txt
# 5,5,5
$exec $tensor $Cname 1 0 0 $Aname 12092,9184,28818