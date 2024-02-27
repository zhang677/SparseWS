#if defined(SPWSHASH)
    #include "benchmark/TTM/COO_CSF_DCSR_hash.h"
#elif defined(SPWSHASHIL)
    #include "benchmark/TTM/COO_CSF_DCSR_hash_il.h"
#else
#endif
#include <iostream>
#include <time.h>
#include <string>
#include <fstream>
using namespace std;

void bench_ttm_dcsf_dcsr_hash(const string filename1, const string filename2, const int repeat, vector<int>& dims, const string& result_name) {
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
    

    vector<string> result;
    boost::split(result, filename1, boost::is_any_of("/"));
    vector<string> result2;
    boost::split(result2, result[result.size()-1], boost::is_any_of("."));

    int w_cap = pow(2,int(log2(1.0 * dims[1] * 4))); // heuristic
    // printf("w_cap: %d\n", w_cap);
    int warmup = 5;
    double duration_taco = COO_CSF_DCSR_hash(&A, &B, &C, w_cap, warmup, repeat, true/*in bench*/);
    std::cout << result2[0] << "," << ncol << "," << duration_taco * 1000 << std::endl; // ms
    std::ofstream outfile;
    outfile.open(result_name, std::ios_base::app);
    outfile << result2[0] << "," << ncol << "," << duration_taco * 1000 << std::endl;
    outfile.close();
}

int main(int argc, char* argv[]) {
    const string filename1 = (argc > 1) ? argv[1] : "./data/test5.mtx";
    const string filename2 = (argc > 2) ? argv[2] : "./data/test5.mtx";
    const int repeat = (argc > 3) ? stoi(argv[3]) : 10;
    const string result_name = (argc > 4) ? argv[4] : "./data/test.csv";
    const string tensorDims = (argc > 5) ? argv[5] : "0,1";

    auto dimsStr = split(tensorDims, ",", false /* keepDelim */);
    std::vector<int> dims;
    for (auto it : dimsStr) {
        dims.push_back(atoi(it.c_str()));
    }

    bench_ttm_dcsf_dcsr_hash(filename1, filename2, repeat, dims, result_name);

    return 0;
}