#include "benchmark/DCSC_DCSR/DCSC_DCSR_noacc_map.h"
#include "benchmark/DCSC_DCSR/DCSC_DCSR_noacc_coord.h"
#include "benchmark/DCSC_DCSR/CSC_CSR_eigen.h"

void check_noacc(const string filename1, const string filename2) {
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
    EigenCSC A_true = to_eigen_csc(nrow, ncol, nnz, id_buffer, indices, value);
    taco_tensor_t A = CC_to_taco_tensor(A_true.outerIndexPtr(),A_true.innerIndexPtr(),A_true.valuePtr(),nrow,ncol,nnz,{1,0});
    //print_taco_tensor_DC(&A);
    indptr.clear();
    id_buffer.clear();
    indices.clear();
    value.clear();
    read_mtx_csr(filename2.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
    taco_tensor_t B = CC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
    EigenCSR B_true = to_eigen_csr(nrow, ncol, nnz, id_buffer, indices, value);
    init_taco_tensor_DC(&C, ncol, ncol, {0,1});

    int w_cap = pow(2,int(log2(nnz))); // heuristic
    CSC_CSR_Eigen(A_true, B_true, C_true, 1, 0);
    std::cout << "coord: " << std::endl;
    DCSC_DCSR_noacc_coord(&A, &B, &C, w_cap, 1, 0);
    check_csr_taco_eigen(C, C_true);   
    free(C.vals);
    free(C.indices[1][0]);
    free(C.indices[1][1]); 
    std::cout << "map: " << std::endl;
    DCSC_DCSR_noacc_map(&A, &B, &C, 1, 0);
    check_csr_taco_eigen(C, C_true);   

}


int main(int argc, char** argv) {
  const string filename1 = (argc > 1) ? argv[1] : "./data/test5.mtx";
  const string filename2 = (argc > 2) ? argv[2] : "./data/test5.mtx";
  const int repeat = (argc > 3) ? stoi(argv[3]) : 10;
  const int w_cap = (argc > 4) ? stoi(argv[4]) : 16;
  const int verbose = (argc > 5) ? stoi(argv[5]) : 0;
  const string result_name = (argc > 6) ? argv[6] : "./data/test.csv";
  // benchmark_eigen_taco(filename1, filename2, repeat, result_name);
  check_noacc(filename1, filename2);
  return 0;
}