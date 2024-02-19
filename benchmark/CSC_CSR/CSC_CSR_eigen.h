#include "../../utils/dataloader.h"
#include "../../utils/lib.h"

double CSC_CSR_Eigen(EigenCSC& A, EigenCSR& B, EigenCSR& C, int32_t warmup, int32_t repeat, bool bench = false, bool print = false) {
  printf("Eigen threads: %d\n", Eigen::nbThreads());
  for (int i = 0; i < warmup; i++) {
    C = A * B;
    if (bench) {
      EigenCSR tmp;
      C = tmp;
    }
  }
  //auto start = std::chrono::high_resolution_clock::now();
  double start = clock();
  for (int i = 0; i < repeat; i++) {
    C = A * B;
    if (bench) {
      EigenCSR tmp;
      C = tmp;
    }
  }
  //auto end = std::chrono::high_resolution_clock::now();
  double end = clock();
  double duration = (double)(end - start) / (CLOCKS_PER_SEC * repeat);
  if (print) {
    print_eigen_csc(A);
    print_eigen_csr(B);
    print_eigen_csr(C);
  }
  return duration;
}