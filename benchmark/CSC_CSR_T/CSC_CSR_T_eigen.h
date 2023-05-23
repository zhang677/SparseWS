#include "../../utils/dataloader.h"
#include "../../utils/lib.h"

double CSC_CSR_T_Eigen(EigenCSC& A, EigenCSR& B, EigenCSR& C, int32_t warmup, int32_t bench, bool print = false) {
  for (int i = 0; i < warmup; i++) {
    C = (A * B).transpose();
  }
  //auto start = std::chrono::high_resolution_clock::now();
  double start = clock();
  for (int i = 0; i < bench; i++) {
    C = (A * B).transpose();
  }
  //auto end = std::chrono::high_resolution_clock::now();
  double end = clock();
  double duration = (double)(end - start) / (CLOCKS_PER_SEC * bench);
  if (print) {
    print_eigen_csc(A);
    print_eigen_csr(B);
    print_eigen_csr(C);
  }
  return duration;
}