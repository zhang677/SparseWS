#include "benchmark/MTTKRP/COO_CSF_DCSR_DCSR_hash_free.h"
#include <iostream>
#include <time.h>
#include <string>
#include <fstream>
using namespace std;

void bench_mttkrp_dcsf_dcsr_dcsr_hash(const string filename1, const string filename2, const string filename3, const int repeat, const int warmup, vector<int>& dims, const string& result_name) {
    taco_tensor_t B;
    read_tns_csf(filename1, B, dims);
    vector<int> indptr;
    vector<int> indices;
    vector<int> id_buffer;
    vector<float> value;
    int nrow;
    int ncol;
    int nnz;
    read_mtx_csr(filename2.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value, false);
    taco_tensor_t C = CC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
    indptr.clear();
    id_buffer.clear();
    indices.clear();
    value.clear();
    read_mtx_csr(filename3.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value, false);
    taco_tensor_t D = CC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
    taco_tensor_t A;
    init_taco_tensor_COO(&A, {dims[2],ncol}, {0,1});

    vector<string> result;
    boost::split(result, filename1, boost::is_any_of("/"));
    vector<string> result2;
    boost::split(result2, result[result.size()-1], boost::is_any_of("."));

    // int w_cap = pow(2,int(log2(dims[2] * ncol * 0.05))); // heuristic
    int w_cap = pow(2,int(log2(dims[2] * ncol * 0.1)));
    printf("w_cap: %d\n", w_cap);
    double duration_taco = COO_CSF_DCSR_DCSR_hash(&A, &B, &C, &D, w_cap, warmup, repeat, true /* in bench*/, false/* verbose*/);
    
    std::cout << result2[0] << "," << ncol << "," << duration_taco * 1000 << std::endl; // ms
    std::ofstream outfile;
    outfile.open(result_name, std::ios_base::app);
    outfile << result2[0] << "," << ncol << "," << duration_taco * 1000 << std::endl;
    outfile.close();
}

int main(int argc, char* argv[]) {
    const string filename1 = (argc > 1) ? argv[1] : "./data/test5.mtx";
    const string filename2 = (argc > 2) ? argv[2] : "./data/test5.mtx";
    const string filename3 = (argc > 3) ? argv[3] : "./data/test5.mtx";
    const int repeat = (argc > 4) ? stoi(argv[4]) : 10;
    const int warmup = (argc > 5) ? stoi(argv[5]) : 5;
    const string result_name = (argc > 6) ? argv[6] : "./data/test.csv";
    const string tensorDims = (argc > 7) ? argv[7] : "0,1";

    auto dimsStr = split(tensorDims, ",", false /* keepDelim */);
    std::vector<int> dims;
    for (auto it : dimsStr) {
        dims.push_back(atoi(it.c_str()));
    }

    bench_mttkrp_dcsf_dcsr_dcsr_hash(filename1, filename2, filename3, repeat, warmup, dims, result_name);

    return 0;
}