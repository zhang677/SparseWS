//#include "seperate/CSR_CSC_0.h" // failed
//#include "evaluate/CSR_CSR_1.h"
//#include "seperate/DCSR_DCSR_4.h"
//#include "seperate/DCSR_DCSC_4.h"
//#include "evaluate/DCSR_DCSC_1.h"
//#include "seperate/CSR_CSR_4.h"
//#include "seperate/CSR_CSR_5.h"
//#include "evaluate/CSR_CSR_Eigen.h"
#include "benchmark/CSR_CSR_T_coord_sort.h"
//#include "benchmark/CSR_CSR_T_generated.h"
//#include "baseline/mmt.h"
//#include "benchmark/CSR_CSR_T_hash.h"
#include "benchmark/CSR_CSR_T_eigen.h"
#include <iostream>
#include <time.h>
#include <string>
using namespace std;

int main(int argc, char* argv[]) {
    const string filename = (argc > 1) ? argv[1] : "./data/test5.mtx";
    const int w_cap = (argc > 2) ? stoi(argv[2]) : 16;
    //mmt(filename);
    taco_tensor_t C;
    clock_t start, finish;
    double duration;

    EigenCSR C_true;
    start = clock();
    CSR_CSR_T_Eigen(filename, filename, C_true);
    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    std::cout << duration << " seconds" << std::endl;

    start = clock();
    CSR_CSR_T_coord_sort(filename, filename, &C, w_cap);
    //CSR_CSR_T_hash(filename, filename, &C); 
    //CSR_CSR_5(filename, filename, &C, w_cap); 
    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    std::cout << duration << " seconds" << std::endl;
    
    check_csr_taco_eigen(C, C_true);
    return 0;    
}