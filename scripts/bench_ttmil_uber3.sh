#!/bin/bash
Tname=uber3
Alg=$1
K=183
L=1140
I=1717
Js=(4 8 16 32 64 128 256 512 1024 2048)
exec=/home/zgh23/code/SparseWS/bench-ttm-$Alg
tensor=/scratch/zgh23/sparse_ten/$Tname.tns
# result=/home/zgh23/code/SparseWS/data/profile/ttm/ttm_il_$Tname.csv
warmup=5
repeat=20
result=/home/zgh23/code/SparseWS/data/test.csv

for J in ${Js[@]}; do
    Cname="/scratch/zgh23/sparse_ten/ctf_generated/001/D-$Tname-$J.mtx"
    $exec $tensor $Cname $repeat $warmup $result $K,$L,$I
done