#!/bin/bash
Tname=freebase_sampled
K=38955429
L=38955429
I=532
Js=(4) # 8 16 32 64 128) # 256 512 1024 2048
exec=/home/zgh23/code/SparseWS/bench-mttkrp-HASH
tensor=/scratch/zgh23/sparse_ten/$Tname.tns
warmup=0
repeat=1
result=/home/zgh23/code/SparseWS/data/profile/mttkrp/mttkrp_new_fb_sampled.csv

for J in ${Js[@]}; do
    Dname="/scratch/zgh23/sparse_ten/generated/D-$Tname-4.mtx"
    Cname="/scratch/zgh23/sparse_ten/generated/C-$Tname-4.mtx"
    $exec $tensor $Cname $Dname $repeat $warmup $result $K,$L,$I
done