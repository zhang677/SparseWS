#!/bin/bash
Tname=facebook
K=1504
L=42390
I=39986
Js=(4 8 16 32 64 128) # 256 512 1024 2048
exec=/home/zgh23/code/SparseWS/bench-mttkrp
tensor=/scratch/zgh23/sparse_ten/$Tname.tns

for J in ${Js[@]}; do
    Dname="/scratch/zgh23/sparse_ten/generated/D-$Tname-4-$J.mtx"
    Cname="/scratch/zgh23/sparse_ten/generated/C-$Tname-4-$J.mtx"
    $exec $tensor $Cname 20 0 0 $Dname $K,$L,$I
done