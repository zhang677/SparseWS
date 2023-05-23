import sparse
import numpy as np
import os
import argparse

def read_tns_file(filename, tensor_name="", info=False):
  """
  Read a sparse tensor from a file in the .tns format.
  """
  
  with open(filename, "r") as f:
    lines = f.readlines()
    nnz = len(lines)
    # Read nonzeros
    indices = []
    values = []
    for i in range(nnz):
      line = lines[i].split()
      indices.append([int(x) - 1 for x in line[:-1]]) # Convert to 0-based indexing
      values.append(float(line[-1]))
    # Convert to sparse tensor
    indices = np.array(indices)
    values = np.array(values)
    dims = np.max(indices, axis=0) + 1
    print("dims: ", dims)
    print("nnz: ", nnz)
    if info:
      return tensor_name + "," + ",".join([str(x) for x in dims])+","+str(nnz)
    else:
      return sparse.COO(indices.T, values, shape=dims)
  
def generate_random_tensor(dims, nnz):
  """
  Generate a random sparse tensor with the given dimensions and number of nonzeros.
  """
  indices = np.random.randint(1, dims, size=(nnz, len(dims))) 
  values = np.random.rand(nnz)
  return sparse.COO(indices.T, values, shape=dims)

def write_tns_file(tensor, filename):
  """
  Write a sparse tensor to a file in the .tns format.
  """
  with open(filename, "w") as f:
    #f.write(" ".join([str(x) for x in tensor.shape]) + "\n")
    #f.write(str(tensor.nnz) + "\n")
    for i in range(tensor.nnz):
      f.write(" ".join([str(x) for x in tensor.coords[:, i]]) + " " + str(tensor.data[i]) + "\n")

def test_random_ttm_t():
  """
  Test the ttm_t function on a random tensor.
  """
  subscripts = ["ijk,kl->jli", "ijk,kl->jil", "ijk,kl->lij", "ijk,kl->lji", "ijk,kl->ijl", "ijk,kl->ilj"]
  dims = [5, 5, 5]
  nnz = 10
  B = generate_random_tensor(dims, nnz)

  cdims = [5, 5]
  cnnz = 5
  C = generate_random_tensor(cdims, cnnz)
  for sub in subscripts:
    A = sparse.einsum(sub, B, C)
    A_true = np.einsum(sub, B.todense(), C.todense())
    assert np.allclose(A_true, A.todense())

def test_subscripts(subscripts, arrays):
  output = sparse.einsum(subscripts, *arrays)
  true_output = np.einsum(subscripts, *(s.todense() for s in arrays))
  assert np.allclose(output, true_output)

def test_tensor_generation():
  root = "/home/nfs_data/zhanggh/SparseWS/data/origin"
  dims = [5, 5, 5]
  nnz = 10
  B = generate_random_tensor(dims, nnz)
  write_tns_file(B, os.path.join(root, "5-5-5-10.tns"))

  cdims = [5, 5]
  cnnz = 5
  C = generate_random_tensor(cdims, cnnz)
  write_tns_file(C, os.path.join(root, "5-5-5.tns"))

def test_all_tensor_info():
  root = "/home/eva_share/datasets/datasets_zhr_nfs/datasets/sparse_ten"
  entry_name_file = "/home/nfs_data/zhanggh/SparseWS/data/tensor_names.txt"
  names = open(entry_name_file, "r").read().split("\n")
  output_file = "/home/nfs_data/zhanggh/SparseWS/data/tensor_info.csv"
  output = open(output_file, "w")
  for name in names:
    filename = os.path.join(root, name, name + ".tns")
    print(read_tns_file(filename, info=True), file=output)

def test_single_tensor_info(name, output):
  root = "/home/eva_share/datasets/datasets_zhr_nfs/datasets/sparse_ten"
  filename = os.path.join(root, name, name + ".tns")
  print(name)
  print(read_tns_file(filename, name, info=True), file=output)

if __name__ == "__main__":
  #test_random_ttm_t()
  #filename = "/home/nfs_data/zhanggh/SparseWS/data/origin/5-5-5-10.tns"
  #filename = "/home/eva_share/datasets/datasets_zhr_nfs/datasets/sparse_ten/nell-2/nell-2.tns"
  #tensor = read_tns_file(filename)
  #print(tensor.coords)
  #print(tensor.data)
  parser = argparse.ArgumentParser(description='Test eigen and coord')
  parser.add_argument("--name", type=str, default="nell-2")
  args = parser.parse_args()
  f = open("/home/nfs_data/zhanggh/SparseWS/data/tensor_info.csv", "a")
  test_single_tensor_info(args.name, f)
  f.close()
