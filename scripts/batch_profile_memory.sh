#!/bin/bash
# input=/home/zgh23/code/SparseWS/data/profile/memory/xtmem_hash_nofuse_$1.txt
# output=/home/zgh23/code/SparseWS/data/ablation/CSC_CSR_T/memory/hash_nofuse.csv
# python /home/zgh23/code/SparseWS/scripts/extract_memory.py --matrix $1 --input $input --output $output
input=/home/zgh23/code/SparseWS/data/profile/memory/xtmem_hash_$1.txt
output=/home/zgh23/code/SparseWS/data/ablation/CSC_CSR_T/memory/hash.csv
python /home/zgh23/code/SparseWS/scripts/extract_memory.py --matrix $1 --input $input --output $output