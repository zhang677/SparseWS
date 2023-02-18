from io import BytesIO
import numpy as np
from scipy.sparse import coo_matrix
from scipy.io import mmwrite
import argparse
import os

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Generate sparse matrix in mtx format')
  parser.add_argument('--row', default=5, type=int, help='# of rows')
  parser.add_argument('--col', default=5, type=int, help='# of columns')
  parser.add_argument('--density', default=0.1, type=float, help='density')
  parser.add_argument('--seed', default=42, type=int, help='Random seed')
  parser.add_argument('--path', default="/home/nfs_data/zhanggh/SparseWS/data", type=str, help='Store path')

  config = parser.parse_args()

  rng = np.random.default_rng(config.seed)
  row = config.row
  col = config.col
  tot = row * col
  z = int((1 - config.density) * tot)
  val = np.random.rand(row, col)
  z_ids = rng.choice(np.arange(tot), size = z, replace=False)
  val = val.flatten()
  val[z_ids] = 0.000000
  val = val.reshape(row, col)
  target = BytesIO()
  mmwrite(target, coo_matrix(val), precision=4)
  output_file = open(os.path.join(config.path, f"{row}-{col}-{config.density}.mtx"), 'w')
  print(target.getvalue().decode('latin1'), file = output_file)

