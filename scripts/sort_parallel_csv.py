import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Sort csv by names')
parser.add_argument('--filename', default="/home/zgh23/code/SparseWS/info/test_matrix_names_merged_sorted.txt", type=str, help='Read path')
parser.add_argument('--table_name', default="/home/zgh23/code/SparseWS/data/results/CSR_CSR/bench-taco-60.csv", type=str, help='Read path')
parser.add_argument('--new_table_name', default="/home/zgh23/code/SparseWS/data/results/CSR_CSR/bench-taco-60-sorted.csv", type=str, help='Read path')
args = parser.parse_args()
filename = args.filename # "/home/zgh23/code/SparseWS/info/test_matrix_names_merged_sorted.txt"
table_name = args.table_name # "/home/zgh23/code/SparseWS/data/results/CSR_CSR/bench-taco-60.csv"
new_table_name = args.new_table_name # "/home/zgh23/code/SparseWS/data/results/CSR_CSR/bench-taco-60-sorted.csv"
names = open(filename, "r").read().split("\n")
ori_table = pd.read_csv(table_name, header=None)
# Sort ori_table by names
ori_table = ori_table.set_index(0)
ori_table = ori_table.reindex(names)
ori_table = ori_table.reset_index()
ori_table.to_csv(new_table_name, header=False, index=False)
