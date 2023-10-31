#include "evaluate/CSC_CSR_T_hash_record_nofuse.h"

void record_hash_outer_transpose_nofuse(const string filename1, const string filename2, const string result_name) {
    taco_tensor_t C;
    EigenCSR C_true;
    vector<int> indptr;
    vector<int> indices;
    vector<int> id_buffer;
    vector<float> value;
    int nrow;
    int ncol;
    int nnz;
    read_mtx_csc(filename1.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
    EigenCSC A_true = to_eigen_csc(nrow, ncol, nnz, id_buffer, indices, value);
    indptr.clear();
    id_buffer.clear();
    indices.clear();
    value.clear();
    taco_tensor_t A = DC_to_taco_tensor(A_true.outerIndexPtr(),A_true.innerIndexPtr(),A_true.valuePtr(),nrow,ncol,nnz,{1,0});
    indptr.clear();
    id_buffer.clear();
    indices.clear();
    value.clear();
    read_mtx_csr(filename2.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
    taco_tensor_t B = DC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
    init_taco_tensor_DC(&C, ncol, ncol, {0,1});
    taco_tensor_t C_noT;
    init_taco_tensor_DC(&C_noT, ncol, ncol, {0,1});
    indptr.clear();
    id_buffer.clear();
    indices.clear();
    value.clear();
    //print_taco_tensor_DC(&A);
    //print_taco_tensor_DC(&B);
    // Check result

    int collisions = 0;
    std::cout << "Hash" << std::endl;
    int w_cap = pow(2,int(log2(nnz))); // heuristic
    //w_cap = 3;
    CSC_CSR_T_hash_nonfuse(&A, &B, &C_noT, &C, w_cap, collisions);
    int output_zeros = C.vals_size;
    // std::cout << collisions << "," << output_zeros << std::endl;
    vector<string> result;
    boost::split(result, filename1, boost::is_any_of("/"));
    vector<string> result2;
    boost::split(result2, result[result.size()-1], boost::is_any_of("."));
    std::ofstream outfile;
    outfile.open(result_name, std::ios_base::app);
    outfile << result2[0] << "," << ncol << "," << collisions << "," << output_zeros << std::endl;
    outfile.close();
}

int main(int argc, char* argv[]) {
    const string filename1 = (argc > 1) ? argv[1] : "./data/test5.mtx";
    const string filename2 = (argc > 2) ? argv[2] : "./data/test5.mtx";
    const string filename3 = (argc > 3) ? argv[3] : "./data/test.csv";
    record_hash_outer_transpose_nofuse(filename1, filename2, filename3);
    return 0;
}