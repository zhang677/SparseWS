#include "../../utils/dataloader.h"
#include "../../utils/lib.h"

double CSR_CSR_Eigen(EigenCSR& A, EigenCSR& B, EigenCSR& C, int32_t warmup, int32_t bench, bool print = false) {
  for (int i = 0; i < warmup; i++) {
    C = A * B;
  }
  Timer timer;
  timer.reset();
  for (int i = 0; i < bench; i++) {
    C = A * B;
  }
  double duration = timer.elapsed() / bench;
  if (print) {
    print_eigen_csr(A);
    print_eigen_csr(B);
    print_eigen_csr(C);
  }
  return duration;
}