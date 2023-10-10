#!/bin/bash
# Js=(16 64 128 256 512 1024 2048)
Js=(4 8 32)
exec=/home/zgh23/code/SparseWS/bench-mttkrp
tensor=/scratch/zgh23/sparse_ten/nell-2.tns

for J in ${Js[@]}; do
    Dname="/home/zgh23/code/ctf/D-nell-2-$J.mtx"
    Cname="/home/zgh23/code/ctf/C-nell-2-$J.mtx"
    $exec $tensor $Cname 20 0 0 $Dname 12092,9184,28818
done