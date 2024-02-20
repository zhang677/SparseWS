#include "benchmark/CSR_CSR/CSR_CSR_eigen.h"
#include "benchmark/Parallel/CSR_CSR_hash_parallel.h"
#include <iostream>
#include <time.h>
#include <string>
#include <fstream>
using namespace std;

void bench_eigen_hash(const string filename1, const string filename2, const int repeat, const string& result_name, pid_t pid) {
    taco_tensor_t C;
    EigenCSR C_true;
    vector<int> indptr;
    vector<int> indices;
    vector<int> id_buffer;
    vector<float> value;
    int nrow;
    int ncol;
    int nnz;
    read_mtx_csr(filename1.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
    taco_tensor_t A = DC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
    EigenCSR A_true = to_eigen_csr(nrow, ncol, nnz, id_buffer, indices, value, false);
    indptr.clear();
    id_buffer.clear();
    indices.clear();
    value.clear();
    read_mtx_csr(filename2.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
    taco_tensor_t B = DC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
    EigenCSR B_true = to_eigen_csr(nrow, ncol, nnz, id_buffer, indices, value, false);
    init_taco_tensor_DC(&C, ncol, ncol, {0,1});
    indptr.clear();
    id_buffer.clear();
    indices.clear();
    value.clear();

    printf("TACO threads: %d\n", omp_get_max_threads());
    printf("Eigen threads = {%d}\n", Eigen::nbThreads());
    clock_t start, finish;
    double duration_taco, duration_eigen;
    const int warmup = 5;

    vector<string> result;
    boost::split(result, filename1, boost::is_any_of("/"));
    vector<string> result2;
    boost::split(result2, result[result.size()-1], boost::is_any_of("."));

    int w_cap = pow(2,int(log2(nnz))); // heuristic
    duration_taco = CSR_CSR_hash(&A, &B, &C, w_cap, warmup, repeat, true);
    duration_eigen = CSR_CSR_Eigen(A_true, B_true, C_true, warmup, repeat, true);
    std::cout << result2[0] << "," << duration_eigen << "," << duration_taco << "," << duration_eigen / duration_taco << "," << pid << std::endl;
    std::ofstream outfile;
    outfile.open(result_name, std::ios_base::app);
    outfile << result2[0] << "," << duration_eigen << "," << duration_taco << "," << duration_eigen / duration_taco << std::endl;
    outfile.close();
}

int main(int argc, char* argv[]) {
    const string filename1 = (argc > 1) ? argv[1] : "./data/test5.mtx";
    const string filename2 = (argc > 2) ? argv[2] : "./data/test5.mtx";
    const int repeat = (argc > 3) ? stoi(argv[3]) : 10;
    const string result_name = (argc > 4) ? argv[4] : "result.csv";
    

    bench_eigen_hash(filename1, filename2, repeat, result_name, getpid());
    return 0;
}