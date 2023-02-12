//#include "seperate/CSR_CSC_0.h" // failed
//#include "evaluate/CSR_CSR_1.h"
//#include "seperate/DCSR_DCSR_4.h"
//#include "seperate/DCSR_DCSC_4.h"
//#include "evaluate/DCSR_DCSC_1.h"
#include "seperate/CSR_CSC_4.h"
//#include "seperate/CSR_CSR_5.h"
#include "evaluate/CSR_CSR_Eigen.h"
#include <iostream>
using namespace std;
int main(int argc, char* argv[]) {
    const string filename = (argc > 1) ? argv[1] : "./data/test5.mtx";
    taco_tensor_t C;
    CSR_CSC_4(filename, filename, &C);
    EigenCSR C_true;
    CSR_CSC_Eigen(filename, filename, C_true);
    int row = C_true.outerSize();
    int nnz = C_true.outerIndexPtr()[row];
    cout << (C.indices[0][0][0] == row) << endl;
    cout << (C.indices[1][0][C.indices[0][0][0]] == nnz) << endl;
    compare_array<int>(C.indices[1][0], C_true.outerIndexPtr(), row + 1);
    compare_array<int>(C.indices[1][1], C_true.innerIndexPtr(), nnz);
    compare_array<float>(C.vals, C_true.valuePtr(), nnz);
    return 0;    
}