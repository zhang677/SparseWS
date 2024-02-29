#!/bin/bash
Tname=nips3
K=2482
L=2862
I=14036
Js=(4 8 16 32 64 128) # 256 512 1024 2048
exec=/home/zgh23/code/SparseWS/bench-mttkrp-HASH
tensor=/scratch/zgh23/sparse_ten/$Tname.tns
warmup=5
repeat=20
result=/home/zgh23/code/SparseWS/data/profile/mttkrp/mttkrp_new_$Tname.csv

for J in ${Js[@]}; do
    Dname="/scratch/zgh23/sparse_ten/generated/D-$Tname-4-$J.mtx"
    Cname="/scratch/zgh23/sparse_ten/generated/C-$Tname-4-$J.mtx"
    $exec $tensor $Cname $Dname $repeat $warmup $result $K,$L,$I
done