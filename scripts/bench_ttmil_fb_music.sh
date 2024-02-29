#!/bin/bash
Tname=freebase_music
K=23343790
L=23344784
I=166
Js=(4) # 8 16 32 64 128) # 256 512 1024 2048)
exec=/home/zgh23/code/SparseWS/bench-ttm-HASHIL
tensor=/scratch/zgh23/sparse_ten/$Tname.tns
result=/home/zgh23/code/SparseWS/data/profile/ttm/ttm_il_fb_music.csv

for J in ${Js[@]}; do
    Cname="/scratch/zgh23/sparse_ten/spws_generated/D-freebase_music-4.mtx"
    $exec $tensor $Cname 1 0 $result $K,$L,$I
done