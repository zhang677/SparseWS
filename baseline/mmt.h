#ifndef MMT_H
#define MMT_H

#include "../utils/dataloader.h"
#include "../utils/lib.h"
#include <mkl_spblas.h>

void csr2csr(int64_t rows, int64_t cols, int *rowptr, int *colidx,
             float *values, int *colptr, int *rowidx, float *values_t) {
  // Step 1. Construct sparse matrix
  sparse_matrix_t A;
  int rows_start[rows], rows_end[rows]; // arrays, could infer from rowptr
  for (int i = 0; i < rows; i++) {
    rows_start[i] = rowptr[i];
    rows_end[i] = rowptr[i + 1];
  }
  checkMKLError(mkl_sparse_s_create_csr(&A, SPARSE_INDEX_BASE_ZERO, (int)rows,
                                        (int)cols, rows_start, rows_end, colidx,
                                        values));

  // Step 2. Transpose
  sparse_matrix_t B;
  int rows_t = (int)cols, cols_t = (int)rows;
  int *cols_start, *cols_end; // prepare for colptr
  checkMKLError(mkl_sparse_convert_csr(A, SPARSE_OPERATION_TRANSPOSE, &B));

  sparse_index_base_t indextype;
  checkMKLError(mkl_sparse_s_export_csr(B, &indextype, &rows_t, &cols_t,
                                        &cols_start, &cols_end, &rowidx,
                                        &values_t));
  for (int i = 0; i < cols; i++) {
    colptr[i] = cols_start[i];
  }
  colptr[cols + 1] = cols_end[cols];
}

void mmt(const string& A_name, bool print = false) {
  vector<int> indptr;
  vector<int> indices;
  vector<int> id_buffer;
  vector<float> value;
  int nrow;
  int ncol;
  int nnz;
  int* colptr;
  int* rowidx;
  float* values_t;
  colptr = (int*)malloc(sizeof(int) * (ncol + 2));
  read_mtx_csr(A_name.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
  csr2csr(nrow, ncol, indptr.data(), indices.data(), value.data(), colptr, rowidx, values_t);
  // print colptr, rowidx, values_t
  if (print)
    print_csr(ncol, nrow, colptr, rowidx, values_t);
  free(colptr);
  return;
}

#endif