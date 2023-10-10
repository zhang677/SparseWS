#!/bin/bash
Js=(16 64 128 256 512 1024 2048)
exec=/home/zgh23/code/SparseWS/bench-mttkrp-splatt
tensor=/scratch/zgh23/sparse_ten/nell-2.tns

for J in ${Js[@]}; do
    $exec $tensor $tensor 1 $J 0 
done