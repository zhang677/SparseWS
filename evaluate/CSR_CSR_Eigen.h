#include "../utils/dataloader.h"
#include "../utils/lib.h"

void CSR_CSR_Eigen(const string& A_name, const string& B_name) {
  vector<int> indptr;
  vector<int> indices;
  vector<int> id_buffer;
  vector<float> value;
  int nrow;
  int ncol;
  int nnz;
  read_mtx_csr(A_name.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
  EigenCSR A = to_eigen_csr(nrow, ncol, nnz, id_buffer, indices, value);
  read_mtx_csr(B_name.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
  EigenCSR B = to_eigen_csr(nrow, ncol, nnz, id_buffer, indices, value);
  EigenCSR C = A * B;

  print_eigen_csr(A);
  print_eigen_csr(B);
  print_eigen_csr(C);
}

void CSR_CSC_Eigen(const string& A_name, const string& B_name, EigenCSR& C, bool print = false) {
  vector<int> indptr;
  vector<int> indices;
  vector<int> id_buffer;
  vector<float> value;
  int nrow;
  int ncol;
  int nnz;
  read_mtx_csr(A_name.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
  EigenCSR A = to_eigen_csr(nrow, ncol, nnz, id_buffer, indices, value);
  read_mtx_csr(B_name.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
  EigenCSC B = to_eigen_csc(nrow, ncol, nnz, id_buffer, indices, value);
  C = A * B;
  if (print) {
    print_eigen_csr(A);
    print_eigen_csc(B);
    print_eigen_csr(C);
  }
}