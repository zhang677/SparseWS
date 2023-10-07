import numpy as np
import argparse
import os
from tqdm import tqdm

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Generate sparse matrix in mtx format')
  parser.add_argument('--row', default=5, type=int, help='# of rows')
  parser.add_argument('--col', default=5, type=int, help='# of cols')
  parser.add_argument('--nnz', default=5, type=int, help='# of non-zeros')
  parser.add_argument('--seed', default=42, type=int, help='Random seed')
  parser.add_argument('--path', default="/home/eva_share/user_file/zhanggh/SpWS/data", type=str, help='Store path')

  config = parser.parse_args()
  name = f"{config.row}-{config.col}-{config.nnz}-{config.seed}.mtx"
  file = open(os.path.join(config.path, name), 'w')
  print("%%MatrixMarket matrix coordinate real general", file = file)
  print("%", file = file)
  print(f"{config.row} {config.row} {config.nnz}", file = file)
  avg = config.nnz // config.row
  cur = 0
  col_candidates = list(range(config.row))
  rng = np.random.default_rng(config.seed)

  for i in tqdm(range(config.row - 1)):
    col_ids = rng.choice(col_candidates, size = avg, replace=False)
    for j in col_ids:
      print(f"{i+1} {j+1} 2", file = file)
      cur += 1

  col_ids = rng.choice(col_candidates, size = config.nnz - cur, replace=False)
  for j in col_ids:
    print(f"{config.row} {j+1} 2", file = file)
  file.close()

