#include "benchmark/MTTKRP/COO_CSF_DCSR_DCSR_hash.h"
#include <iostream>
#include <time.h>
#include <string>
#include <fstream>
using namespace std;

void check_mttkrp_dcsf_dcsr_dcsr_hash(const string filename1, const string filename2, const string filename3, const string filename4, const int repeat, bool verbose, vector<int>& dims) {
    taco_tensor_t B;
    read_tns_csf(filename1, B, dims);
    std::cout << "B: " << B.dimensions[0] << "," << B.dimensions[1] << "," << B.dimensions[2] << "," << B.vals_size << std::endl;
    vector<int> indptr;
    vector<int> indices;
    vector<int> id_buffer;
    vector<float> value;
    int nrow;
    int ncol;
    int nnz;
    read_mtx_csr_keep_value(filename2.data(), nrow, ncol, nnz, indptr, indices, value, false /*one_base*/);
    taco_tensor_t C = CC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
    std::cout << "C: " << nrow << "," << ncol << "," << nnz << std::endl; // 12092,32,
    indptr.clear();
    id_buffer.clear();
    indices.clear();
    value.clear();
    read_mtx_csr_keep_value(filename3.data(), nrow, ncol, nnz, indptr, indices, value, false /*one_base*/);
    taco_tensor_t D = CC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
    std::cout << "D: " << nrow << "," << ncol << "," << nnz << std::endl; // 9184,32,
    taco_tensor_t A;
    init_taco_tensor_DC(&A, dims[2], ncol, {0,1});
    std::cout << "A: " << nrow << "," << ncol << std::endl; // 28818,32

    int w_cap = pow(2,int(log2(1.0 * nnz * dims[2] / nrow))); // heuristic
    std::cout << "w_cap: " << w_cap << std::endl; 
    int warmup = 1;
    double duration_taco = COO_CSF_DCSR_DCSR_hash(&A, &B, &C, &D, w_cap, warmup, repeat, false /*in bench*/, verbose);

    id_buffer.clear();
    indices.clear();
    value.clear();
    read_mtx_coo(filename4.data(), nrow, ncol, nnz, id_buffer, indices, value, false/*one_base*/);
    cout << "nnz: " << A.vals_size << "," << nnz << endl;
    compare_array<int>(A.indices[1][0], id_buffer.data(), nnz);
    compare_array<int>(A.indices[1][1], indices.data(), nnz);
    compare_array<float>(A.vals, value.data(), nnz);
}

int main(int argc, char* argv[]) {
    if (argc < 6) {
        std::cerr << "Usage: " << argv[0] << " B.mtx C.mtx D.mtx A.mtx repeat" << std::endl;
        return 1;
    }
    const string Bname(argv[1]);
    const string Cname(argv[2]);
    const string Dname(argv[3]);
    const string Aname(argv[4]);
    int repeat = atoi(argv[5]);
    const string tensorDims(argv[6]); 
    bool verbose = false;
    auto dimsStr = split(tensorDims, ",", false /* keepDelim */);
    std::vector<int> dims;
    for (auto it : dimsStr) {
        dims.push_back(atoi(it.c_str()));
    }
    check_mttkrp_dcsf_dcsr_dcsr_hash(Bname, Cname, Dname, Aname, repeat, verbose, dims);
    return 0;
}