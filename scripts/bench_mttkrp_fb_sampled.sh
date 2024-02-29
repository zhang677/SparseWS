#!/bin/bash
Tname=freebase_sampled
J=$1
K=38955429
L=38955429
I=532
it=1
wp=0
exec=/home/zgh23/code/SparseWS/bench-mttkrp
tensor=/scratch/zgh23/sparse_ten/$Tname.tns
Dname="/scratch/zgh23/sparse_ten/generated/D-$Tname-$J.mtx"
Cname="/scratch/zgh23/sparse_ten/generated/C-$Tname-$J.mtx"
$exec $tensor $Cname 20 0 0 $Dname $K,$L,$I