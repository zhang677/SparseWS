#include "benchmark/CSR_CSR/CSR_CSR_eigen.h"
#include "benchmark/Parallel/DCSC_DCSR_hash_parallel.h"

void check_eigen_hash_outer_CC_noT(const string filename1, const string filename2, bool verbose) {
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
    EigenCSR A_use = to_eigen_csr(nrow, ncol, nnz, id_buffer, indices, value);
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
    EigenCSR B_true = to_eigen_csr(nrow, ncol, nnz, id_buffer, indices, value, false);
    init_taco_tensor_DC(&C, ncol, ncol, {0,1});
    indptr.clear();
    id_buffer.clear();
    indices.clear();
    value.clear();

    std::cout << "Hash" << std::endl;
    int w_cap = pow(2,int(log2(nnz))); // heuristic
    DCSC_DCSR_hash(&A, &B, &C, w_cap, 0, 1, false);
    std::cout << "Eigen" << std::endl;
    CSR_CSR_Eigen(A_use, B_true, C_true, 0, 1, false);
    int outSize = C_true.outerSize();
    std::cout << "Filename " << filename1 << std::endl;
    std::cout << "Output Nnz: " << C_true.outerIndexPtr()[outSize] << std::endl;
    check_csr_taco_eigen(C, C_true);
}
int main(int argc, char* argv[]) {
    const string filename1 = (argc > 1) ? argv[1] : "./data/test5.mtx";
    const string filename2 = (argc > 2) ? argv[2] : "./data/test5.mtx";
    const int verbose = (argc > 3) ? stoi(argv[3]) : 0;

    check_eigen_hash_outer_CC_noT(filename1, filename2, verbose);
    return 0;
}