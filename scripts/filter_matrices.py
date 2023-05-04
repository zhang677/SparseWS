from io import BytesIO
from scipy.io import mmread, mminfo
import random
import argparse
import os
import tqdm

def calc_sparsity(sym, row, col, nnz):
  if sym == "symmetric":
    return nnz / (row * col / 2)
  else:
    return nnz / (row * col)
    
def ablation_bin_split():
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

if __name__ == '__main__':
  random.seed(0)
  ablation_bin_split()
