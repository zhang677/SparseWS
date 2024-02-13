#!/bin/bash
Tname=uber3
K=183
L=1140
I=1717
Js=(4 8 16 32 64 128 256 512 1024 2048)
exec=/home/zgh23/code/SparseWS/bench-mttkrp
tensor=/scratch/zgh23/sparse_ten/$Tname.tns

for J in ${Js[@]}; do
    Dname="/home/zgh23/code/ctf/D-$Tname-$J.mtx"
    Cname="/home/zgh23/code/ctf/C-$Tname-$J.mtx"
    $exec $tensor $Cname 20 0 0 $Dname $K,$L,$I
done