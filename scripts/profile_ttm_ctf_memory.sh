#!/bin/bash
Tname=facebook
profile_root=/home/zgh23/code/SparseWS/data/profile/ttm/memory
K=1504
I=42390
J=39986
it=0
wp=1
L=4
profile_name=ctf_${Tname}_$L
exec=/home/zgh23/code/ctf/bin/myttm_il
Bname="/scratch/zgh23/sparse_ten/$Tname-zero.tns"
Cname="/home/zgh23/code/ctf/D-$Tname-4.txt"
valgrind --tool=massif --xtree-memory=full --xtree-memory-file=$profile_root/$profile_name.kcg --massif-out-file=$profile_root/$profile_name.out mpirun -n 1 $exec -tensor $Bname -dims $K,$I,$J -iter $it -warmup $wp -ttmLDim $L -mode 1 -matrixC $Cname
# ms_print $profile_root/$profile_name.out > $profile_root/$profile_name.txt
callgrind_annotate --auto=yes --inclusive=yes --sort=curB:100,curBk:100,totB:100,totBk:100,totFdB:100,totFdBk:100 $profile_root/$profile_name.kcg > $profile_root/xtmem_$profile_name.txt
