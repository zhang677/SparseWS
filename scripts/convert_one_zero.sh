#!/bin/bash
Root=/scratch/zgh23/sparse_ten
Names=(nell-1 flickr-3d delicious-3d freebase_music freebase_sampled)
for name in ${Names[@]}; do
    python convert_one_zero_tns.py --input $Root/$name.tns --output $Root/$name-zero.tns
done