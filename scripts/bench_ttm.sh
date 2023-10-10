#!/bin/bash
Js=(4 8 16 32 64 128 256 512 1024 2048)
exec=/home/zgh23/code/SparseWS/bench-ttm
tensor=/scratch/zgh23/sparse_ten/nell-2.tns
for J in ${Js[@]}; do
    Cname="/home/zgh23/code/ctf/C-nell-2-$J.mtx"
    $exec $tensor $Cname 20 0 0 none 12092,9184,28818
done