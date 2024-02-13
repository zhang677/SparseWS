#!/bin/bash
# Tname=uber3
# K=183
# I=1140
# J=1717
# Tname=nell-2
# K=12092
# I=9184
# J=28818
Tname=facebook
K=1504
I=42390
J=39986
Ls=(4 8 16 32 64 128 256 512 1024 2048)
exec=/home/zgh23/code/SparseWS/bench-ttm
tensor=/scratch/zgh23/sparse_ten/$Tname.tns
for L in ${Ls[@]}; do
    Cname="/home/zgh23/code/ctf/C-$Tname-$L.mtx"
    $exec $tensor $Cname 20 0 0 none $K,$I,$J
done