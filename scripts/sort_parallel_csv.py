import pandas as pd

filename = "/home/nfs_data/zhanggh/SparseWS/info/sampled_matrix_names_sorted.txt"
table_name = "/home/nfs_data/zhanggh/SparseWS/data/results/CSC_CSR_T/test-hash-mt.csv"
new_table_name = "/home/nfs_data/zhanggh/SparseWS/data/results/CSC_CSR_T/test-hash-mt-sorted.csv"
names = open(filename, "r").read().split("\n")
ori_table = pd.read_csv(table_name, header=None)
# Sort ori_table by names
ori_table = ori_table.set_index(0)
ori_table = ori_table.reindex(names)
ori_table = ori_table.reset_index()
ori_table.to_csv(new_table_name, header=False, index=False)
