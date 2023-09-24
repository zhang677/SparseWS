from io import BytesIO
from scipy.io import mmread, mminfo
import random
import argparse
import os
import tqdm
import pandas as pd

def calc_sparsity(sym, row, col, nnz):
  if sym == "symmetric":
    return nnz / (row * col / 2)
  else:
    return nnz / (row * col)
    
def ablation_bin_split_old():
  parser = argparse.ArgumentParser(description='Filter out matrix coordiante real/integer general/symmetric')
  parser.add_argument('--input_path', default="/home/eva_share/datasets/sparse_mat", type=str, help='Read path')
  parser.add_argument('--output_file', default="/home/nfs_data/zhanggh/SparseWS/data/sampled_matrix_names.txt", type=str, help='Store path')
  parser.add_argument('--output_info', default="/home/nfs_data/zhanggh/SparseWS/data/sampled_matrix_info.csv", type=str, help='Store Info path')
  args = parser.parse_args()

  bins = [12000,180000]
  sample_num = 10
  bin_name = {"small":[],"medium":[],"large":[]}
  sampled_name = {"small":[],"medium":[],"large":[]}
  in_dirname = args.input_path
  mtxnames = [fname for fname in os.listdir(in_dirname)]
  mtx_files = [os.path.join(in_dirname, fname, fname + ".mtx") for fname in mtxnames]  
  for i, mtx in enumerate(mtx_files):
    rows, cols, entries, form, field, sym = mminfo(mtx)
    if form == "coordinate" and (field == "real" or field == "integer") and (sym == "general" or sym == "symmetric"):
      if entries < bins[0]:
        bin_name["small"].append(mtxnames[i])
      elif entries < bins[1]:
        bin_name["medium"].append(mtxnames[i])
      else:
        bin_name["large"].append(mtxnames[i]) 
  print(len(bin_name["small"]),len(bin_name["medium"]),len(bin_name["large"]))
  for key in bin_name.keys():
    sampled_name[key] = random.sample(bin_name[key],sample_num)
  print(sampled_name)
  out_file = open(args.output_file, 'w')
  out_info = open(args.output_info, 'w')
  print("name,rows,cols,entries,form,field,sym,sparsity", file = out_info)
  for key in sampled_name.keys():
    for name in sampled_name[key]:
      mtx_path = os.path.join(in_dirname, name, name + ".mtx")
      rows, cols, entries, form, field, sym = mminfo(mtx_path)
      print(",".join([name,str(rows),str(cols),str(entries),form,field,sym,str(calc_sparsity(sym, rows, cols, entries))]), file = out_info)
      print(name, file = out_file)
  out_file.close()
  out_info.close()

def ablation_bin_split_new():
  parser = argparse.ArgumentParser(description='Sample test matrices')
  parser.add_argument('--database', default="/home/zgh23/code/SparseWS/info/matrix_info.csv", type=str, help='All matrix info')
  parser.add_argument('--output_info', default="/home/zgh23/code/SparseWS/info/test_matrix_info", type=str, help='Store Info path')
  parser.add_argument('--sample_num', default=10, type=int, help='Number of matrices to sample')
  parser.add_argument('--seed', default=0, type=int, help='Random seed')
  args = parser.parse_args()

  bins = [12000,180000]
  bin_name = {"small":[],"medium":[],"large":[]}
  sampled_name = {"small":[],"medium":[],"large":[]}
  # Ablation_bin_split_new is different from Ablation_bin_split_old in that the new one reads all matrix info from a csv file in a format of: name,row,col,nnz,sparsity. 
  # The old one reads matrix info from .mtx files.
  # Enumerate all matrices and put them into different bins.
  df = pd.read_csv(args.database)
  for i in range(len(df)):
    row = df.iloc[i]
    if row["nnz"] < bins[0]:
      bin_name["small"].append(row["name"])
    elif row["nnz"] < bins[1]:
      bin_name["medium"].append(row["name"])
    else:
      bin_name["large"].append(row["name"])
  print(len(bin_name["small"]),len(bin_name["medium"]),len(bin_name["large"]))
  # Set random seed to make sure that the sampled matrices are the same for each run.
  random.seed(args.seed)
  # Sample matrices from each bin.
  for key in bin_name.keys():
    sampled_name[key] = random.sample(bin_name[key],args.sample_num)
  print(sampled_name)
  # Write the info of sampled matrices to a csv file.
  filename = args.output_info + str(args.seed) + "_" + str(3 * args.sample_num) + ".csv"
  out_info = open(filename, 'w')
  print("name,rows,cols,entries,sparsity", file = out_info)
  for key in sampled_name.keys():
    for name in sampled_name[key]:
      row = df.loc[df["name"] == name]
      print(",".join([name,str(row["row"].values[0]),str(row["col"].values[0]),str(row["nnz"].values[0]),str(row["sparsity"].values[0])]), file = out_info)
  out_info.close()


def matrix_collection():
  parser = argparse.ArgumentParser(description='Filter out matrix coordiante real/integer general/symmetric')
  parser.add_argument('--input_path', default="/home/nfs_data/zhanggh/SparseWS/data", type=str, help='Read path')
  parser.add_argument('--output_file', default="/home/nfs_data/zhanggh/SparseWS/data/matrix_names.txt", type=str, help='Store path')
  args = parser.parse_args()

  in_dirname = args.input_path
  out_file = open(args.output_file, 'w')
  print("name,row,col,nnz,sparsity", file = out_file)
  mtxnames = [fname for fname in os.listdir(in_dirname)]
  mtx_files = [os.path.join(in_dirname, fname, fname + ".mtx") for fname in mtxnames]
  for i, mtx in enumerate(mtx_files):
    rows, cols, entries, form, field, sym = mminfo(mtx)
    if form == "coordinate" and (field == "real" or field == "integer") and (sym == "general" or sym == "symmetric"):        
      print(",".join([mtxnames[i],str(rows),str(cols),str(entries),str(calc_sparsity(sym, rows, cols, entries))]), file = out_file)   
  out_file.close()

def find_largest_matrix():
  parser = argparse.ArgumentParser(description='Filter out matrix coordiante real/integer general/symmetric')
  parser.add_argument('--input_path', default="/home/eva_share/datasets/sparse_mat", type=str, help='Read path')
  args = parser.parse_args()

  in_dirname = args.input_path
  mtxnames = [fname for fname in os.listdir(in_dirname)]
  mtx_files = [os.path.join(in_dirname, fname, fname + ".mtx") for fname in mtxnames]  
  max_nnz = 0
  max_file = ""
  for mtx in mtx_files:
    _, _, entries, form, field, sym = mminfo(mtx)
    if form == "coordinate" and (field == "real" or field == "integer") and (sym == "general" or sym == "symmetric"):
      if entries > max_nnz:
        max_nnz = entries
        max_file = mtx
  print(max_file, max_nnz)

def print_matrix_info():
  parser = argparse.ArgumentParser(description='Filter out matrix coordiante real/integer general/symmetric')
  parser.add_argument('--input_path', default="/home/eva_share/datasets/sparse_mat/Freescale2/Freescale2.mtx", type=str, help='Read path')
  args = parser.parse_args()

  filename = args.input_path
  rows, cols, entries, form, field, sym = mminfo(filename)
  print(rows, cols, entries, form, field, sym)

def output_names_csv():
  parser = argparse.ArgumentParser(description='Output matrix names from csv to a txt file.')
  parser.add_argument('--input_path', default="/home/zgh23/code/SparseWS/info/test_matrix_info42_30_sorted.csv", type=str, help='Read path')
  parser.add_argument('--output_file', default="/home/zgh23/code/SparseWS/info/test_matrix_names42_30_sorted.txt", type=str, help='Store path')
  args = parser.parse_args()

  df = pd.read_csv(args.input_path)
  output_file = open(args.output_file, 'w')
  for i in range(len(df)):
    print(df.iloc[i]["name"], file = output_file)
  output_file.close()


if __name__ == '__main__':
  random.seed(0)
  # ablation_bin_split()
  # find_largest_matrix()
  # print_matrix_info()
  # ablation_bin_split_new()
  output_names_csv()