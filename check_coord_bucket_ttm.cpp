#include "benchmark/TTM/COO_CSF_DCSR_bucket_il.h"
#include "benchmark/TTM/COO_CSF_DCSR_coord_il.h"
#include <iostream>
#include <time.h>
#include <string>
#include <fstream>
using namespace std;

void check_ttmil_coord_bucket(const string filename1, const string filename2, vector<int>& dims) {
    taco_tensor_t B;
    read_tns_csf(filename1, B, dims);
    vector<int> indptr;
    vector<int> indices;
    vector<int> id_buffer;
    vector<float> value;
    int nrow;
    int ncol;
    int nnz;
    read_mtx_csr_keep_value(filename2.data(), nrow, ncol, nnz, indptr, indices, value, false /*one_base*/);
    taco_tensor_t C = CC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
    //print_taco_tensor_CC(&C);
    indptr.clear();
    id_buffer.clear();
    indices.clear();
    value.clear();
    taco_tensor_t A;
    if (nrow == dims[0]) {
        init_taco_tensor_COO(&A, {dims[1],dims[2],ncol}, {0,1,2});
    } else if (nrow == dims[1]) {
        init_taco_tensor_COO(&A, {dims[0],dims[2],ncol}, {0,1,2});
    } else {
        init_taco_tensor_COO(&A, {dims[0],dims[1],ncol}, {0,1,2});
    }

    taco_tensor_t A_coord;
    if (nrow == dims[0]) {
        init_taco_tensor_COO(&A_coord, {dims[1],dims[2],ncol}, {0,1,2});
    } else if (nrow == dims[1]) {
        init_taco_tensor_COO(&A_coord, {dims[0],dims[2],ncol}, {0,1,2});
    } else {
        init_taco_tensor_COO(&A_coord, {dims[0],dims[1],ncol}, {0,1,2});
    }
    printf("nrow: %d, dims[0]: %d, dims[1]: %d, dims[2]: %d\n", nrow, dims[0], dims[1], dims[2]);
    int w_cap = pow(2,int(log2(1.0 * dims[1] * ncol * 0.05)));
    COO_CSF_DCSR_bucket(&A, &B, &C, w_cap, 1, 0, false/*in bench*/);
    COO_CSF_DCSR_coord(&A_coord, &B, &C, w_cap);

    check_coo_taco_taco_3d(A, A_coord);
    store_taco_tensor_coo("/home/zgh23/code/SparseWS/data/A_bucket.txt", &A);
    store_taco_tensor_coo("/home/zgh23/code/SparseWS/data/A_coord.txt", &A_coord);
}

int main(int argc, char* argv[]) {
    const string filename1 = (argc > 1) ? argv[1] : "./data/test5.mtx";
    const string filename2 = (argc > 2) ? argv[2] : "./data/test5.mtx";
    const string tensorDims = (argc > 3) ? argv[3] : "0,1";

    auto dimsStr = split(tensorDims, ",", false /* keepDelim */);
    std::vector<int> dims;
    for (auto it : dimsStr) {
        dims.push_back(atoi(it.c_str()));
    }

    check_ttmil_coord_bucket(filename1, filename2, dims);

    return 0;
}