#!/bin/bash
Tname=nell-2
K=12092
L=9184
I=28818
# Tname=uber3
# K=183
# L=1140
# I=1717
J=4
exec=/home/zgh23/code/SparseWS/check-ttm-COORDBUCKET
tensor=/scratch/zgh23/sparse_ten/$Tname.tns
Cname="/scratch/zgh23/sparse_ten/ctf_generated/001/D-$Tname-$J.mtx"
$exec $tensor $Cname $K,$L,$I