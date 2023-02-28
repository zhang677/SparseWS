//#include "seperate/CSR_CSC_0.h" // failed
//#include "evaluate/CSR_CSR_1.h"
//#include "seperate/DCSR_DCSR_4.h"
//#include "seperate/DCSR_DCSC_4.h"
//#include "evaluate/DCSR_DCSC_1.h"
//#include "seperate/CSR_CSR_4.h"
//#include "seperate/CSR_CSR_5.h"
//#include "evaluate/CSR_CSR_Eigen.h"
//#include "benchmark/CSR_CSR_T_coord.h"
#include "benchmark/CSR_CSR_T_generated.h"
#include "benchmark/CSR_CSR_T_eigen.h"
#include <iostream>
using namespace std;

int main(int argc, char* argv[]) {
    const string filename = (argc > 1) ? argv[1] : "./data/test5.mtx";
    taco_tensor_t C;
    //CSR_CSR_T_coord(filename, filename, &C); 
    CSR_CSR_T_generated(filename, filename, &C); 
    EigenCSR C_true;
    CSR_CSR_T_Eigen(filename, filename, C_true);
    check_csr_taco_eigen(C, C_true);
    return 0;    
}