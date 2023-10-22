#!/bin/bash
# name1=/home/zgh23/code/SparseWS/data/ablation/CSC_CSR/cc/bench-cc-noT-bucket.csv
# name2=/home/zgh23/code/SparseWS/data/ablation/CSC_CSR/cc/bench-cc-noT-bucket-sorted.csv
# python sort_parallel_csv.py --filename /home/zgh23/code/SparseWS/info/temp_names.txt --table_name $name1 --new_table_name $name2
# name1=/home/zgh23/code/SparseWS/data/ablation/CSC_CSR/cc/bench-cc-noT-coord-c.csv
# name2=/home/zgh23/code/SparseWS/data/ablation/CSC_CSR/cc/bench-cc-noT-coord-c-sorted.csv
# python sort_parallel_csv.py --filename /home/zgh23/code/SparseWS/info/temp_names.txt --table_name $name1 --new_table_name $name2
# name1=/home/zgh23/code/SparseWS/data/ablation/CSC_CSR/cc/bench-cc-noT-coord-cf.csv
# name2=/home/zgh23/code/SparseWS/data/ablation/CSC_CSR/cc/bench-cc-noT-coord-cf-sorted.csv
# python sort_parallel_csv.py --filename /home/zgh23/code/SparseWS/info/temp_names.txt --table_name $name1 --new_table_name $name2
# name1=/home/zgh23/code/SparseWS/data/ablation/CSC_CSR/cc/bench-cc-noT-hash.csv
# name2=/home/zgh23/code/SparseWS/data/ablation/CSC_CSR/cc/bench-cc-noT-hash-sorted.csv
# python sort_parallel_csv.py --filename /home/zgh23/code/SparseWS/info/temp_names.txt --table_name $name1 --new_table_name $name2
# name1=/home/zgh23/code/SparseWS/data/ablation/CSC_CSR/cc/bench-cc-noT-hash-mt.csv
# name2=/home/zgh23/code/SparseWS/data/ablation/CSC_CSR/cc/bench-cc-noT-hash-mt-sorted.csv
# python sort_parallel_csv.py --filename /home/zgh23/code/SparseWS/info/temp_names.txt --table_name $name1 --new_table_name $name2
# name1=/home/zgh23/code/SparseWS/data/ablation/CSC_CSR/cc/bench-cc-noT-hash-flex.csv
# name2=/home/zgh23/code/SparseWS/data/ablation/CSC_CSR/cc/bench-cc-noT-hash-flex-sorted.csv
# python sort_parallel_csv.py --filename /home/zgh23/code/SparseWS/info/temp_names.txt --table_name $name1 --new_table_name $name2
# name1=/home/zgh23/code/SparseWS/data/ablation/CSC_CSR/cc/bench-cc-noT-hash-flex-mt.csv
# name2=/home/zgh23/code/SparseWS/data/ablation/CSC_CSR/cc/bench-cc-noT-hash-flex-mt-sorted.csv
# python sort_parallel_csv.py --filename /home/zgh23/code/SparseWS/info/temp_names.txt --table_name $name1 --new_table_name $name2
# name1=/home/zgh23/code/SparseWS/data/ablation/CSC_CSR/cc/bench-cc-noT-outer-taco.csv
# name2=/home/zgh23/code/SparseWS/data/ablation/CSC_CSR/cc/bench-cc-noT-outer-taco-sorted.csv
# python sort_parallel_csv.py --filename /home/zgh23/code/SparseWS/info/temp_names.txt --table_name $name1 --new_table_name $name2

Js=(bucket coord-c coord-cf hash hash-mt hash-flex hash-flex-mt taco)

for J in ${Js[@]}; do
    name1=/home/zgh23/code/SparseWS/data/ablation/CSR_CSR/shift/bench-noT-row-$J.csv
    name2=/home/zgh23/code/SparseWS/data/ablation/CSR_CSR/shift/bench-noT-row-$J-sorted.csv
    python sort_parallel_csv.py --filename /home/zgh23/code/SparseWS/info/temp_names_2.txt --table_name $name1 --new_table_name $name2
done