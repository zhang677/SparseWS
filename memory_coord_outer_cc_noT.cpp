#if defined(SPWSCOORDC)
    #include "benchmark/DCSC_DCSR/DCSC_DCSR_coord_index_chase.h"
#elif defined(SPWSCOORDCF)
    #include "benchmark/DCSC_DCSR/DCSC_DCSR_coord_index_chase_flex.h"
#else
#endif
#include <iostream>
#include <time.h>
#include <string>
#include <fstream>
using namespace std;

void single_coord_outer_CC_noT(const string filename1, const string filename2) {
    taco_tensor_t C;
    vector<int> indptr;
    vector<int> indices;
    vector<int> id_buffer;
    vector<float> value;
    int nrow;
    int ncol;
    int nnz;
    read_mtx_csr(filename1.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
    EigenCSC A_true = to_eigen_csc(nrow, ncol, nnz, id_buffer, indices, value);
    indptr.clear();
    id_buffer.clear();
    indices.clear();
    value.clear();
    taco_tensor_t A = CC_to_taco_tensor(A_true.outerIndexPtr(),A_true.innerIndexPtr(),A_true.valuePtr(),nrow,ncol,nnz,{1,0});
    indptr.clear();
    id_buffer.clear();
    indices.clear();
    value.clear();
    read_mtx_csr(filename2.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
    taco_tensor_t B = CC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
    init_taco_tensor_DC(&C, ncol, ncol, {0,1});
    indptr.clear();
    id_buffer.clear();
    indices.clear();
    value.clear();

    int w_cap = pow(2,int(log2(nnz))); // heuristic
    DCSC_DCSR_coord(&A, &B, &C, w_cap);
}

int main(int argc, char* argv[]) {
    const string filename1 = (argc > 1) ? argv[1] : "./data/test5.mtx";
    const string filename2 = (argc > 2) ? argv[2] : "./data/test5.mtx";

    single_coord_outer_CC_noT(filename1, filename2);

    return 0;
}