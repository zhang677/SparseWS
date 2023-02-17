//#include "seperate/CSR_CSC_0.h" // failed
//#include "evaluate/CSR_CSR_1.h"
//#include "seperate/DCSR_DCSR_4.h"
//#include "seperate/DCSR_DCSC_4.h"
//#include "evaluate/DCSR_DCSC_1.h"
#include "seperate/CSR_CSR_4.h"
//#include "seperate/CSR_CSR_5.h"
#include "evaluate/CSR_CSR_Eigen.h"
#include <iostream>
using namespace std;

int main(int argc, char* argv[]) {
    const string filename = (argc > 1) ? argv[1] : "./data/test5.mtx";
    taco_tensor_t C;
    CSR_CSR_4(filename, filename, &C); // C(i,k) = A(i,j) * B(j,k); C: CSR, A: CSR, B: CSC
    EigenCSR C_true;
    CSR_CSR_Eigen(filename, filename, C_true, true);
    check_csr_taco_eigen(C, C_true);
    return 0;    
}