from io import BytesIO
from scipy.sparse import coo_matrix
from scipy.io import mmwrite
import argparse
import os
from pathlib import Path
from tqdm import tqdm

from util import SuiteSparseTensor, InputCacheSuiteSparse
from sam.util import SUITESPARSE_FORMATTED_PATH, ScipyTensorShifter

class ScipyMatrixMarketTensorStore:
    def __init__(self, outdir: str):
        self.outdir = outdir

    def store(self, coo: coo_matrix, name: str):
        target = BytesIO()
        mmwrite(target, coo)
        output_file = open(os.path.join(self.outdir, name), 'w')
        print(target.getvalue().decode('latin1'), file = output_file)
        output_file.close()


  
if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Shift and transpose sparse matrix in mtx format')
  parser.add_argument('--name', metavar='ssname', type=str, action='store', help='tensor name to run format conversion')
  parser.add_argument('--input_path', default="/scratch/zgh23/sparse_mat", type=str, help='Store path')
  parser.add_argument('--output_path', default="/scratch/zgh23/sparse_mat_t", type=str, help='Store path')
  parser.add_argument('--tiles', action='store_true')
  args = parser.parse_args()

  cwd = os.getcwd()

  if args.name is None and args.tiles is False:
      print("Please enter a matrix name")
      exit()

  in_dirname = args.input_path

  if args.output_path is None:
      out_dirname = SUITESPARSE_FORMATTED_PATH
  else:
      out_dirname = args.output_path

  out_path = Path(out_dirname)
  out_path.mkdir(parents=True, exist_ok=True, mode=0o777)
  if args.name is not None:
    out_name = args.name + "-st.mtx"

  outputStore = ScipyMatrixMarketTensorStore(out_path)
  inputCache = InputCacheSuiteSparse()

  tensor = None
  mtx_files = None

  if args.tiles:
      print(in_dirname)
      #mtxnames = [fname for fname in os.listdir(in_dirname) if fname.endswith(".mtx")]
      mtxnames = [fname for fname in os.listdir(in_dirname)]
      exist_mtxnames = [fname.replace("-st.mtx", ".mtx") for fname in os.listdir(out_dirname)]
      mtxnames = list(set(mtxnames) - set(exist_mtxnames))
      print(mtxnames)
      mtx_files = [os.path.join(in_dirname, fname) for fname in mtxnames]
      #mtx_files = [os.path.join(in_dirname, fname, fname + ".mtx") for fname in mtxnames]
      out_mtx_files = [os.path.join(out_dirname, fname.replace(".mtx", "-st.mtx")) for fname in mtxnames]
      #out_mtx_files = [os.path.join(out_dirname, fname + "-st.mtx") for fname in mtxnames]
      tensors = [SuiteSparseTensor(mtx_file) for mtx_file in mtx_files]
      for i, ten in tqdm(enumerate(tensors)):
           out_name = out_mtx_files[i]
           coo = inputCache.load(ten, False)
           shifted = ScipyTensorShifter().shiftLastMode(coo)
           trans_shifted = shifted.transpose()
           outputStore.store(trans_shifted, out_name)
  else:
      tensor = SuiteSparseTensor(os.path.join(in_dirname, args.name + ".mtx"))
      coo = inputCache.load(tensor, False)
      shifted = ScipyTensorShifter().shiftLastMode(coo)
      trans_shifted = shifted.transpose()
      outputStore.store(trans_shifted, out_name)




