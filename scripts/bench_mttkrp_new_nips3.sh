#!/bin/bash
Tname=nips3
Alg=$1
K=2482
L=2862
I=14036
Js=(4 8 16 32 64 128 256 512 1024 2048)
exec=/home/zgh23/code/SparseWS/bench-mttkrp-$Alg
tensor=/scratch/zgh23/sparse_ten/$Tname.tns
warmup=5
repeat=20
result=/home/zgh23/code/SparseWS/data/profile/mttkrp/ablation/001/mttkrp-$Alg-$Tname.csv
# result=/home/zgh23/code/SparseWS/data/result.csv

# for J in ${Js[@]}; do
#     Dname="/scratch/zgh23/sparse_ten/generated/D-$Tname-4-$J.mtx"
#     Cname="/scratch/zgh23/sparse_ten/generated/C-$Tname-4-$J.mtx"
#     $exec $tensor $Cname $Dname $repeat $warmup $result $K,$L,$I
# done

for J in ${Js[@]}; do
    Dname="/scratch/zgh23/sparse_ten/ctf_generated/001/D-$Tname-$J.mtx"
    Cname="/scratch/zgh23/sparse_ten/ctf_generated/001/C-$Tname-$J.mtx"
    $exec $tensor $Cname $Dname $repeat $warmup $result $K,$L,$I
done