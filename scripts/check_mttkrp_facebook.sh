#!/bin/bash
Tname=facebook
K=1504
L=42390
I=39986
Js=(4) # 8 16 32 64 128 256 512 1024 2048
exec=/home/zgh23/code/SparseWS/check-mttkrp-HASH
tensor=/scratch/zgh23/sparse_ten/$Tname.tns

for J in ${Js[@]}; do
    Dname="/home/zgh23/code/ctf/D-$Tname-$J.mtx"
    Cname="/home/zgh23/code/ctf/C-$Tname-$J.mtx"
    Aname="/home/zgh23/code/ctf/A-$Tname-$J.mtx"
    $exec $tensor $Cname $Dname $Aname 1 $K,$L,$I
done