#!/bin/bash
Tname=facebook
profile_root=/home/zgh23/code/SparseWS/data/profile/ttm/memory
K=1504
I=42390
J=39986
it=0
wp=1
L=4
profile_name=hash_${Tname}_$L
exec=/home/zgh23/code/SparseWS/bench-ttm-HASHIL
Bname="/scratch/zgh23/sparse_ten/$Tname.tns"
Cname="/scratch/zgh23/sparse_ten/generated/D-$Tname-4-$L.mtx"
result=/home/zgh23/code/SparseWS/data/test.csv
valgrind --tool=massif --xtree-memory=full --xtree-memory-file=$profile_root/$profile_name.kcg --massif-out-file=$profile_root/$profile_name.out  $exec $Bname $Cname $it $wp $result $K,$I,$J 
# ms_print $profile_root/$profile_name.out > $profile_root/$profile_name.txt
callgrind_annotate --auto=yes --inclusive=yes --sort=curB:100,curBk:100,totB:100,totBk:100,totFdB:100,totFdBk:100 $profile_root/$profile_name.kcg > $profile_root/xtmem_$profile_name.txt
