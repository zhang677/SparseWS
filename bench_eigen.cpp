#include "benchmark/CSR_CSR/CSR_CSR_eigen.h"
#include <iostream>
#include <time.h>
#include <string>
#include <fstream>
using namespace std;

void bench_eigen_noT(const string filename1, const string filename2, const int repeat, const string& result_name, pid_t pid) {
    EigenCSR C_true;
    vector<int> indptr;
    vector<int> indices;
    vector<int> id_buffer;
    vector<float> value;
    int nrow;
    int ncol;
    int nnz;
    read_mtx_csr(filename1.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
    EigenCSR A_use = to_eigen_csr(nrow, ncol, nnz, id_buffer, indices, value);
    indptr.clear();
    id_buffer.clear();
    indices.clear();
    value.clear();
    read_mtx_csr(filename2.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
    EigenCSR B_true = to_eigen_csr(nrow, ncol, nnz, id_buffer, indices, value, false);
    indptr.clear();
    id_buffer.clear();
    indices.clear();
    value.clear();

    clock_t start, finish;
    double duration_eigen;
    const int warmup = 5;

    vector<string> result;
    boost::split(result, filename1, boost::is_any_of("/"));
    vector<string> result2;
    boost::split(result2, result[result.size()-1], boost::is_any_of("."));

    duration_eigen = CSR_CSR_Eigen(A_use, B_true, C_true, warmup, repeat, true);
    std::cout << result2[0] << "," << duration_eigen << "," << pid << std::endl;
    std::ofstream outfile;
    outfile.open(result_name, std::ios_base::app);
    outfile << result2[0] << "," << duration_eigen << std::endl;
    outfile.close();
}

int main(int argc, char* argv[]) {
    const string filename1 = (argc > 1) ? argv[1] : "./data/test5.mtx";
    const string filename2 = (argc > 2) ? argv[2] : "./data/test5.mtx";
    const int repeat = (argc > 3) ? stoi(argv[3]) : 10;
    const string result_name = (argc > 4) ? argv[4] : "./data/test.csv";

    pid_t pid = getpid();


    bench_eigen_noT(filename1, filename2, repeat, result_name, pid);

    return 0;
}