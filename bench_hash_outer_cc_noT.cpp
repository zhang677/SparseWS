// #include "benchmark/DCSC_DCSR/DCSC_DCSR_hash.h"
// #include "benchmark/DCSC_DCSR/DCSC_DCSR_hash_mt.h"
// #include "benchmark/DCSC_DCSR/DCSC_DCSR_hash_flex.h"
// #include "benchmark/DCSC_DCSR/DCSC_DCSR_hash_flex_mt.h"
#if defined(SPWSBUCKET)
    #include "benchmark/DCSC_DCSR/DCSC_DCSR_coord_bucket.h"
#elif defined(SPWSHASH)
    #include "benchmark/DCSC_DCSR/DCSC_DCSR_hash.h"
#elif defined(SPWSHASHMT)
    #include "benchmark/DCSC_DCSR/DCSC_DCSR_hash_mt.h"
#elif defined(SPWSHASHFLEX)
    #include "benchmark/DCSC_DCSR/DCSC_DCSR_hash_flex.h"
#elif defined(SPWSHASHFLEXMT)
    #include "benchmark/DCSC_DCSR/DCSC_DCSR_hash_flex_mt.h"
#else
#endif
#include <iostream>
#include <time.h>
#include <string>
#include <fstream>
using namespace std;

void bench_hash_outer_CC_noT(const string filename1, const string filename2, const int repeat, const string& result_name, pid_t pid) {
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

    clock_t start, finish;
    double duration_taco;
    const int warmup = 5;

    vector<string> result;
    boost::split(result, filename1, boost::is_any_of("/"));
    vector<string> result2;
    boost::split(result2, result[result.size()-1], boost::is_any_of("."));

    int w_cap = pow(2,int(log2(nnz))); // heuristic
    duration_taco = DCSC_DCSR_hash(&A, &B, &C, w_cap, warmup, repeat, true);
    std::cout << result2[0] << "," << duration_taco << "," << pid << std::endl;
    std::ofstream outfile;
    outfile.open(result_name, std::ios_base::app);
    outfile << result2[0] << "," << duration_taco << std::endl;
    outfile.close();
}

int main(int argc, char* argv[]) {
    const string filename1 = (argc > 1) ? argv[1] : "./data/test5.mtx";
    const string filename2 = (argc > 2) ? argv[2] : "./data/test5.mtx";
    const int repeat = (argc > 3) ? stoi(argv[3]) : 10;
    const string result_name = (argc > 4) ? argv[4] : "./data/test.csv";

    pid_t pid = getpid();


    bench_hash_outer_CC_noT(filename1, filename2, repeat, result_name, pid);

    return 0;
}