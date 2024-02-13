# include "benchmark/DCSC_DCSR/DCSC_DCSR_taco.h"
# include "benchmark/CSR_CSR/CSR_CSR_eigen.h"

#include <iostream>
#include <time.h>
#include <string>
#include <fstream>
using namespace std;

void benchmark_eigen_taco(const string filename1, const string filename2, const int repeat, const string& result_name) {
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
    EigenCSC A_tmp = to_eigen_csc(nrow, ncol, nnz, id_buffer, indices, value);
    EigenCSR A_true = to_eigen_csr(nrow, ncol, nnz, id_buffer, indices, value, false);
    taco_tensor_t A = CC_to_taco_tensor(A_true.outerIndexPtr(),A_true.innerIndexPtr(),A_true.valuePtr(),nrow,ncol,nnz,{1,0});
    //print_taco_tensor_DC(&A);
    indptr.clear();
    id_buffer.clear();
    indices.clear();
    value.clear();
    read_mtx_csr(filename2.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
    taco_tensor_t B = CC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
    EigenCSR B_true = to_eigen_csr(nrow, ncol, nnz, id_buffer, indices, value, false);
    init_taco_tensor_DC(&C, ncol, ncol, {0,1});

    double duration_eigen, duration_taco;
    const int warmup = 5;

    vector<string> result;
    boost::split(result, filename1, boost::is_any_of("/"));
    vector<string> result2;
    boost::split(result2, result[result.size()-1], boost::is_any_of("."));

    duration_taco = DCSC_DCSR_taco(&A, &B, &C, warmup, repeat, true);
    duration_eigen = CSR_CSR_Eigen(A_true, B_true, C_true, warmup, repeat, true);
    std::cout << result2[0] << "," << duration_eigen << "," << duration_taco << "," << duration_eigen / duration_taco << std::endl;
    std::ofstream outfile;
    outfile.open(result_name, std::ios_base::app);
    outfile << result2[0] << "," << duration_eigen << "," << duration_taco << "," << duration_eigen / duration_taco << std::endl;
    outfile.close();
}

void check_eigen_taco(const string filename1, const string filename2) {
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
    EigenCSC A_tmp = to_eigen_csc(nrow, ncol, nnz, id_buffer, indices, value);
    EigenCSR A_true = to_eigen_csr(nrow, ncol, nnz, id_buffer, indices, value);
    taco_tensor_t A = CC_to_taco_tensor(A_tmp.outerIndexPtr(),A_tmp.innerIndexPtr(),A_tmp.valuePtr(),nrow,ncol,nnz,{1,0});
    //print_taco_tensor_DC(&A);
    indptr.clear();
    id_buffer.clear();
    indices.clear();
    value.clear();
    read_mtx_csr(filename2.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
    taco_tensor_t B = CC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
    EigenCSR B_true = to_eigen_csr(nrow, ncol, nnz, id_buffer, indices, value);
    init_taco_tensor_DC(&C, ncol, ncol, {0,1});

    DCSC_DCSR_taco(&A, &B, &C, 1, 0);
    CSR_CSR_Eigen(A_true, B_true, C_true, 1,0);
    check_csr_taco_eigen(C, C_true);
}

int main(int argc, char** argv) {
  const string filename1 = (argc > 1) ? argv[1] : "./data/test5.mtx";
  const string filename2 = (argc > 2) ? argv[2] : "./data/test5.mtx";
  const int repeat = (argc > 3) ? stoi(argv[3]) : 10;
  const string result_name = (argc > 4) ? argv[4] : "./data/test.csv";
  // benchmark_eigen_taco(filename1, filename2, repeat, result_name);
  check_eigen_taco(filename1, filename2);
  return 0;
}