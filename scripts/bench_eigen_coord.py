import os
import argparse
import subprocess

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Test eigen and coord')
  parser.add_argument('--name', metavar='ssname', type=str, action='store', help='tensor name to run format conversion')
  parser.add_argument('--input_path', default="/home/nfs_data/zhanggh/SparseWS/data/origin", type=str, help='Store path')
  parser.add_argument('--input_tran_path', default="/home/nfs_data/zhanggh/SparseWS/data/shift", type=str, help='Store path')
  parser.add_argument('--output_path', default="/home/nfs_data/zhanggh/SparseWS/data/results/bench_eigen_coord.csv", type=str, help='Store path')
  parser.add_argument('--tiles', action='store_true')
  args = parser.parse_args()

  cwd = os.getcwd()

  if args.name is None and args.tiles is False:
      print("Please enter a matrix name")
      exit()

  in_dirname = args.input_path
  in_tran_dirname = args.input_tran_path

  if args.output_path is None:
      out_dirname = os.getenv('RESULT_PATH', default=cwd)
  else:
      out_dirname = args.output_path

  if args.tiles:
      out_file = open(args.output_path, "w")
      mtxnames = [fname for fname in os.listdir(in_dirname) if fname.endswith(".mtx")]
      print(mtxnames)
      mtx_files = [os.path.join(in_dirname, fname) for fname in mtxnames]
      tran_mtx_files = [os.path.join(in_tran_dirname, fname.replace(".mtx", "-st.mtx")) for fname in mtxnames]
      subprocess.run(["echo", "name,eigen,coord"], stdout=out_file)
      for i in range(len(mtx_files)):
          subprocess.run(["/home/nfs_data/zhanggh/SparseWS/test", mtx_files[i], tran_mtx_files[i]], stdout=out_file)
