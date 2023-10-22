//#include "seperate/CSR_CSC_0.h" // failed
//#include "evaluate/CSR_CSR_1.h"
//#include "seperate/DCSR_DCSR_4.h"
//#include "seperate/DCSR_DCSC_4.h"
//#include "evaluate/DCSR_DCSC_1.h"
//#include "seperate/CSR_CSR_4.h"
//#include "seperate/CSR_CSR_5.h"
//#include "evaluate/CSR_CSR_Eigen.h"
//#include "benchmark/CSR_CSR_T_coord.h"
//#include "benchmark/CSR_CSR_T_generated.h"
//#include "benchmark/CSR_CSR_T_coord_index_double_chase.h"
//#include "benchmark/CSR_CSR_T_coord_index_double.h"
//#include "benchmark/CSR_CSR_T_coord_index.h"
//#include "benchmark/CSR_CSR_T_coord_index_flex.h"
//#include "benchmark/CSR_CSR_T_hash.h"
//#include "benchmark/CSR_CSR_T_hash_flex.h"
//#include "benchmark/CSR_CSR_T_hash_index.h"
//#include "benchmark/CSR_CSR_T_map.h"
//#include "benchmark/CSR_CSR_T_eigen.h"
// #include "benchmark/CSC_CSR_T/CSC_CSR_T_coord_index_chase.h"
// #include "benchmark/CSC_CSR_T/CSC_CSR_T_coord_index_chase_flex.h"
// #include "benchmark/CSC_CSR_T/CSC_CSR_T_hash.h"
// #include "benchmark/CSC_CSR_T/CSC_CSR_T_hash_flex.h"
// #include "benchmark/CSC_CSR_T/CSC_CSR_T_hash_mt.h"
// #include "benchmark/CSC_CSR_T/CSC_CSR_T_hash_flex_mt.h"
// #include "benchmark/CSC_CSR_T/CSC_CSR_T_coord_bucket.h"
//#include "benchmark/CSC_CSR_T/CSC_CSR_T_map.h"
//#include "benchmark/CSC_CSR_T/CSC_CSR_T_map_bucket_rb.h"
// #include "benchmark/DCSC_DCSR_T/DCSC_DCSR_T_hash.h"
// #include "benchmark/DCSC_DCSR_T/DCSC_DCSR_T_hash_mt.h"
// #include "benchmark/DCSC_DCSR_T/DCSC_DCSR_T_hash_flex.h"
// #include "benchmark/DCSC_DCSR_T/DCSC_DCSR_T_hash_flex_mt.h"
// #include "benchmark/DCSC_DCSR_T/DCSC_DCSR_T_coord_index_chase.h"
// #include "benchmark/DCSC_DCSR_T/DCSC_DCSR_T_coord_index_chase_flex.h"
// #include "benchmark/DCSC_DCSR_T/DCSC_DCSR_T_coord_bucket.h"
// #include "benchmark/MTTKRP/COO_CSF_DCSR_DCSR_hash.h"
// #include "benchmark/MTTKRP/MTTKRP_splatt.h"
#include "benchmark/CSC_CSR_T/nonfuse/CSC_CSR_T_hash.h"
// #include "benchmark/TTM/COO_CSF_DCSR_hash.h"
#include "benchmark/CSC_CSR_T/CSC_CSR_T_eigen.h"
// #include "benchmark/DCSC_DCSR/DCSC_DCSR_hash.h"
// #include "benchmark/DCSC_DCSR/DCSC_DCSR_hash_mt.h"
// #include "benchmark/DCSC_DCSR/DCSC_DCSR_hash_flex.h"
// #include "benchmark/DCSC_DCSR/DCSC_DCSR_hash_flex_mt.h"
// #include "benchmark/DCSC_DCSR/DCSC_DCSR_coord_bucket.h"
// #include "benchmark/DCSC_DCSR/DCSC_DCSR_coord_index_chase.h"
// #include "benchmark/DCSC_DCSR/DCSC_DCSR_coord_index_chase_flex.h"
// #include "benchmark/CSR_CSR/CSR_CSR_eigen.h"
//#include "benchmark/CSR_CSR/CSR_CSR_taco.h"
// #include "benchmark/CSR_CSR/CSR_CSR_coord_index_chase_flex.h"
// #include "benchmark/CSR_CSR/CSR_CSR_coord_index_chase.h"
// #include "benchmark/CSR_CSR/CSR_CSR_coord_bucket.h"
// #include "benchmark/CSR_CSR/CSR_CSR_hash_mt.h"
// #include "benchmark/CSR_CSR/CSR_CSR_hash_flex.h"
// #include "benchmark/CSR_CSR/CSR_CSR_hash_flex_mt.h"
// #include "benchmark/CSR_CSR/CSR_CSR_hash.h"
// #include "benchmark/CSC_CSR_T/CSC_CSR_T_coord_bucket_time.h"
#include <iostream>
#include <time.h>
#include <string>
#include <fstream>
using namespace std;

// void assemble_hash(const string filename1, const string filename2, const int repeat, int w_cap, bool verbose) {
//     taco_tensor_t C;
//     EigenCSR C_true;
//     vector<int> indptr;
//     vector<int> indices;
//     vector<int> id_buffer;
//     vector<float> value;
//     int nrow;
//     int ncol;
//     int nnz;
//     read_mtx_csr(filename1.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
//     EigenCSC A_true = to_eigen_csc(nrow, ncol, nnz, id_buffer, indices, value);
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     taco_tensor_t A = DC_to_taco_tensor(A_true.outerIndexPtr(),A_true.innerIndexPtr(),A_true.valuePtr(),nrow,ncol,nnz,{1,0});
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     read_mtx_csr(filename2.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
//     taco_tensor_t B = DC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
//     EigenCSR B_true = to_eigen_csr(nrow, ncol, nnz, id_buffer, indices, value, false);
//     init_taco_tensor_DC(&C, ncol, ncol, {0,1});

//     w_cap = pow(2,int(log2(nnz))); // heuristic
//     std::cout << "Hash" << std::endl;
//     double duration = CSC_CSR_T_hash(&A, &B, &C, w_cap, 0, repeat, verbose);
//     std::cout << "Hash Duration: " << duration << "s" << std::endl;
//     duration = CSC_CSR_T_Eigen(A_true, B_true, C_true, 0, repeat);
//     std::cout << "Eigen Duration: " << duration << "s" << std::endl;
// }

// void get_output_nnz(const string filename1, const string filename2, const string result_name) {
//     EigenCSR C_true;
//     vector<int> indptr;
//     vector<int> indices;
//     vector<int> id_buffer;
//     vector<float> value;
//     int nrow;
//     int ncol;
//     int nnz;
//     read_mtx_csr(filename1.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
//     EigenCSC A_true = to_eigen_csc(nrow, ncol, nnz, id_buffer, indices, value);
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     read_mtx_csr(filename2.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
//     EigenCSR B_true = to_eigen_csr(nrow, ncol, nnz, id_buffer, indices, value, false);

//     vector<string> result;
//     boost::split(result, filename1, boost::is_any_of("/"));
//     vector<string> result2;
//     boost::split(result2, result[result.size()-1], boost::is_any_of("."));
//     CSC_CSR_T_Eigen(A_true, B_true, C_true, 0, 1);
//     std::ofstream outfile;
//     outfile.open(result_name, std::ios_base::app);
//     outfile << result2[0] << "," << ncol << "," << C_true.outerIndexPtr()[C_true.outerSize()] << std::endl;
//     outfile.close();
// }

// void check_eigen_taco(const string filename1, const string filename2, const int repeat, bool verbose) {
//     taco_tensor_t C;
//     EigenCSR C_true;
//     vector<int> indptr;
//     vector<int> indices;
//     vector<int> id_buffer;
//     vector<float> value;
//     int nrow;
//     int ncol;
//     int nnz;
//     read_mtx_csr(filename1.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
//     taco_tensor_t A = DC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
//     EigenCSR A_true = to_eigen_csr(nrow, ncol, nnz, id_buffer, indices, value, false);
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     read_mtx_csr(filename2.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
//     taco_tensor_t B = DC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
//     EigenCSR B_true = to_eigen_csr(nrow, ncol, nnz, id_buffer, indices, value, false);
//     init_taco_tensor_DC(&C, ncol, ncol, {0,1});

//     std::cout << "Eigen" << std::endl;
//     CSR_CSR_Eigen(A_true, B_true, C_true, 0, 1, verbose);
//     std::cout << "TACO" << std::endl;
//     CSR_CSR_taco(&A, &B, &C, 0, 1, verbose);
//     check_csr_taco_eigen(C, C_true);
// }

// void check_eigen_coord(const string filename1, const string filename2, const int repeat, int w_cap, bool verbose) {
//     taco_tensor_t C;
//     EigenCSR C_true;
//     vector<int> indptr;
//     vector<int> indices;
//     vector<int> id_buffer;
//     vector<float> value;
//     int nrow;
//     int ncol;
//     int nnz;
//     read_mtx_csr(filename1.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
//     taco_tensor_t A = DC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
//     EigenCSR A_true = to_eigen_csr(nrow, ncol, nnz, id_buffer, indices, value, false);
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     read_mtx_csr(filename2.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
//     taco_tensor_t B = DC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
//     EigenCSR B_true = to_eigen_csr(nrow, ncol, nnz, id_buffer, indices, value, false);
//     init_taco_tensor_DC(&C, ncol, ncol, {0,1});

//     std::cout << "Eigen" << std::endl;
//     w_cap = pow(2,int(log2(nnz))); // heuristic
//     CSR_CSR_Eigen(A_true, B_true, C_true, 0, 1, verbose);
//     std::cout << "Coord" << std::endl;
//     CSR_CSR_coord(&A, &B, &C, w_cap, verbose);
//     check_csr_taco_eigen(C, C_true);
// }

// void check_eigen_hash(const string filename1, const string filename2, const int repeat, int w_cap, bool verbose) {
//     taco_tensor_t C;
//     EigenCSR C_true;
//     vector<int> indptr;
//     vector<int> indices;
//     vector<int> id_buffer;
//     vector<float> value;
//     int nrow;
//     int ncol;
//     int nnz;
//     read_mtx_csr(filename1.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
//     taco_tensor_t A = DC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
//     EigenCSR A_true = to_eigen_csr(nrow, ncol, nnz, id_buffer, indices, value, false);
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     read_mtx_csr(filename2.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
//     taco_tensor_t B = DC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
//     EigenCSR B_true = to_eigen_csr(nrow, ncol, nnz, id_buffer, indices, value, false);
//     init_taco_tensor_DC(&C, ncol, ncol, {0,1});

//     std::cout << "Eigen" << std::endl;
//     w_cap = pow(2,int(log2(nnz))); // heuristic
//     CSR_CSR_Eigen(A_true, B_true, C_true, 0, 1, verbose);
//     std::cout << "Hash" << std::endl;
//     CSR_CSR_hash(&A, &B, &C, w_cap, 0, 1, verbose);
//     check_csr_taco_eigen(C, C_true);
// }

// void check_eigen_hash_outer(const string filename1, const string filename2, const int repeat, int w_cap, bool verbose) {
//     taco_tensor_t C;
//     EigenCSR C_true;
//     vector<int> indptr;
//     vector<int> indices;
//     vector<int> id_buffer;
//     vector<float> value;
//     int nrow;
//     int ncol;
//     int nnz;
//     read_mtx_csc(filename1.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
//     EigenCSC A_true = to_eigen_csc(nrow, ncol, nnz, id_buffer, indices, value);
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     taco_tensor_t A = DC_to_taco_tensor(A_true.outerIndexPtr(),A_true.innerIndexPtr(),A_true.valuePtr(),nrow,ncol,nnz,{1,0});
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     read_mtx_csr(filename2.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
//     taco_tensor_t B = DC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
//     EigenCSR B_true = to_eigen_csr(nrow, ncol, nnz, id_buffer, indices, value, false);
//     init_taco_tensor_DC(&C, ncol, ncol, {0,1});
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     //print_taco_tensor_DC(&A);
//     //print_taco_tensor_DC(&B);
//     // Check result


//     std::cout << "Hash" << std::endl;
//     w_cap = pow(2,int(log2(nnz))); // heuristic
//     //w_cap = 3;
//     CSC_CSR_T_hash(&A, &B, &C, w_cap, 0, 1, false, verbose);
//     std::cout << "Eigen" << std::endl;
//     CSC_CSR_T_Eigen(A_true, B_true, C_true, 0, 1, false, verbose);
//     int outSize = C_true.outerSize();
//     std::cout << "Filename " << filename1 << std::endl;
//     std::cout << "Output Nnz: " << C_true.outerIndexPtr()[outSize] << std::endl;
//     check_csr_taco_eigen(C, C_true);
// }



// void check_eigen_coord_outer(const string filename1, const string filename2, const int repeat, int w_cap, bool verbose) {
//     taco_tensor_t C;
//     EigenCSR C_true;
//     vector<int> indptr;
//     vector<int> indices;
//     vector<int> id_buffer;
//     vector<float> value;
//     int nrow;
//     int ncol;
//     int nnz;
//     read_mtx_csc(filename1.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
//     EigenCSC A_true = to_eigen_csc(nrow, ncol, nnz, id_buffer, indices, value);
//     EigenCSR A_true_csr = to_eigen_csr(nrow, ncol, nnz, id_buffer, indices, value);
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     taco_tensor_t A = DC_to_taco_tensor(A_true.outerIndexPtr(),A_true.innerIndexPtr(),A_true.valuePtr(),nrow,ncol,nnz,{1,0});
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     read_mtx_csr(filename2.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
//     taco_tensor_t B = DC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
//     EigenCSR B_true = to_eigen_csr(nrow, ncol, nnz, id_buffer, indices, value, false);
//     init_taco_tensor_DC(&C, ncol, ncol, {0,1});

//     std::cout << "Eigen" << std::endl;
//     CSC_CSR_T_Eigen(A_true, B_true, C_true, 0, 1, verbose);
//     int outSize = C_true.outerSize();
//     std::cout << "Filename " << filename1 << std::endl;
//     //std::cout << "Output sparsity " << double(C_true.outerIndexPtr()[outSize]) / double(outSize * outSize) << std::endl; // overflow
//     std::cout << "Output Nnz: " << C_true.outerIndexPtr()[outSize] << std::endl;

//     std::cout << "Coord" << std::endl;
//     w_cap = pow(2,int(log2(nnz))); // heuristic
//     CSC_CSR_T_coord(&A, &B, &C, w_cap, verbose);
//     check_csr_taco_eigen(C, C_true);
// }



// void check_eigen_map_outer(const string filename1, const string filename2, const int repeat, bool verbose) {
//     taco_tensor_t C;
//     EigenCSR C_true;
//     vector<int> indptr;
//     vector<int> indices;
//     vector<int> id_buffer;
//     vector<float> value;
//     int nrow;
//     int ncol;
//     int nnz;
//     read_mtx_csc(filename1.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
//     EigenCSC A_true = to_eigen_csc(nrow, ncol, nnz, id_buffer, indices, value);
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     taco_tensor_t A = DC_to_taco_tensor(A_true.outerIndexPtr(),A_true.innerIndexPtr(),A_true.valuePtr(),nrow,ncol,nnz,{1,0});
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     read_mtx_csr(filename2.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
//     taco_tensor_t B = DC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
//     EigenCSR B_true = to_eigen_csr(nrow, ncol, nnz, id_buffer, indices, value, false);
//     init_taco_tensor_DC(&C, ncol, ncol, {0,1});
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     //print_taco_tensor_DC(&A);
//     //print_taco_tensor_DC(&B);
//     // Check result


//     std::cout << "Map" << std::endl;
//     CSC_CSR_T_map(&A, &B, &C, 0, 1, false, verbose);
//     std::cout << "Eigen" << std::endl;
//     CSC_CSR_T_Eigen(A_true, B_true, C_true, 0, 1, false, verbose);
//     int outSize = C_true.outerSize();
//     std::cout << "Filename " << filename1 << std::endl;
//     std::cout << "Output Nnz: " << C_true.outerIndexPtr()[outSize] << std::endl;
//     check_csr_taco_eigen(C, C_true);
// }

// void check_eigen_map(const string filename1, const string filename2, const int repeat, const int w_cap, bool verbose) {
//     taco_tensor_t C;
//     EigenCSR C_true;
//     std::cout << "Eigen" << std::endl;
//     CSR_CSR_T_Eigen(filename1, filename2, C_true, verbose);
//     //CSR_CSR_5(filename1, filename2, &C, w_cap);
//     std::cout << "Hash" << std::endl;
//     CSR_CSR_T_map(filename1, filename2, &C, w_cap, verbose);
//     check_csr_taco_eigen(C, C_true);
// }

// TODO: Sparse TTM
// void check_splatt_coord(const string filename1, const string filename2, const string chkfilename,const int w_cap, bool verbose) {
//     taco_tensor_t C;
//     splatt_csf* A = read_tensor(filename1);
//     taco_tensor_t A_taco = to_taco_tensor(A);
//     splatt_csf* C_true = read_tensor(chkfilename);

// }

// void benchmark_eigen_taco(const string filename1, const string filename2, const int repeat, const string& result_name) {
//     taco_tensor_t C;
//     EigenCSR C_true;
//     vector<int> indptr;
//     vector<int> indices;
//     vector<int> id_buffer;
//     vector<float> value;
//     int nrow;
//     int ncol;
//     int nnz;
//     read_mtx_csr(filename1.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
//     taco_tensor_t A = DC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
//     EigenCSR A_true = to_eigen_csr(nrow, ncol, nnz, id_buffer, indices, value, false);
//     //print_taco_tensor_DC(&A);
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     read_mtx_csr(filename2.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
//     taco_tensor_t B = DC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
//     EigenCSR B_true = to_eigen_csr(nrow, ncol, nnz, id_buffer, indices, value, false);
//     init_taco_tensor_DC(&C, ncol, ncol, {0,1});

//     clock_t start, finish;
//     double duration_eigen, duration_taco;
//     const int warmup = 5;

//     vector<string> result;
//     boost::split(result, filename1, boost::is_any_of("/"));
//     vector<string> result2;
//     boost::split(result2, result[result.size()-1], boost::is_any_of("."));

//     duration_taco = CSR_CSR_taco(&A, &B, &C, warmup, repeat);
//     duration_eigen = CSR_CSR_Eigen(A_true, B_true, C_true, warmup, repeat);
//     std::cout << result2[0] << "," << duration_eigen << "," << duration_taco << "," << duration_eigen / duration_taco << std::endl;
//     std::ofstream outfile;
//     outfile.open(result_name, std::ios_base::app);
//     outfile << result2[0] << "," << duration_eigen << "," << duration_taco << "," << duration_eigen / duration_taco << std::endl;
//     outfile.close();
// }

// void benchmark_eigen_coord(const string filename1, const string filename2, const int repeat, int w_cap, const string& result_name) {
//     taco_tensor_t C;
//     EigenCSR C_true;
//     vector<int> indptr;
//     vector<int> indices;
//     vector<int> id_buffer;
//     vector<float> value;
//     int nrow;
//     int ncol;
//     int nnz;
//     read_mtx_csr(filename1.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
//     taco_tensor_t A = DC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
//     EigenCSR A_true = to_eigen_csr(nrow, ncol, nnz, id_buffer, indices, value, false);
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     read_mtx_csr(filename2.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
//     taco_tensor_t B = DC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
//     EigenCSR B_true = to_eigen_csr(nrow, ncol, nnz, id_buffer, indices, value, false);
//     init_taco_tensor_DC(&C, ncol, ncol, {0,1});
    
//     clock_t start, finish;
//     double duration_eigen, duration_taco;
//     const int warmup = 5;

//     vector<string> result;
//     boost::split(result, filename1, boost::is_any_of("/"));
//     vector<string> result2;
//     boost::split(result2, result[result.size()-1], boost::is_any_of("."));


//     duration_eigen = CSR_CSR_Eigen(A_true, B_true, C_true, warmup, repeat, true);

//     w_cap = pow(2,int(log2(nnz))); // heuristic
//     for (int i = 0; i < warmup; i++) {
//         CSR_CSR_coord(&A, &B, &C, w_cap);
//     }
//     start = clock();
//     for (int i = 0; i < repeat; i++) { 
//         CSR_CSR_coord(&A, &B, &C, w_cap);
//     }
//     finish = clock();
//     duration_taco = (double)(finish - start) / (CLOCKS_PER_SEC * repeat);
//     std::cout << result2[0] << "," << duration_eigen << "," << duration_taco << "," << duration_eigen / duration_taco << std::endl;
//     std::ofstream outfile;
//     outfile.open(result_name, std::ios_base::app);
//     outfile << result2[0] << "," << duration_eigen << "," << duration_taco << "," << duration_eigen / duration_taco << std::endl;
//     outfile.close();
// }

// void benchmark_eigen_hash(const string filename1, const string filename2, const int repeat, int w_cap, const string& result_name) {
//     taco_tensor_t C;
//     EigenCSR C_true;
//     vector<int> indptr;
//     vector<int> indices;
//     vector<int> id_buffer;
//     vector<float> value;
//     int nrow;
//     int ncol;
//     int nnz;
//     read_mtx_csr(filename1.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
//     taco_tensor_t A = DC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
//     EigenCSR A_true = to_eigen_csr(nrow, ncol, nnz, id_buffer, indices, value, false);
//     //print_taco_tensor_DC(&A);
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     read_mtx_csr(filename2.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
//     taco_tensor_t B = DC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
//     EigenCSR B_true = to_eigen_csr(nrow, ncol, nnz, id_buffer, indices, value, false);
//     init_taco_tensor_DC(&C, ncol, ncol, {0,1});

//     clock_t start, finish;
//     double duration_eigen, duration_taco;
//     const int warmup = 5;

//     vector<string> result;
//     boost::split(result, filename1, boost::is_any_of("/"));
//     vector<string> result2;
//     boost::split(result2, result[result.size()-1], boost::is_any_of("."));

//     w_cap = pow(2,int(log2(nnz))); // heuristic
//     duration_taco = CSR_CSR_hash(&A, &B, &C, w_cap, warmup, repeat);
//     duration_eigen = CSR_CSR_Eigen(A_true, B_true, C_true, warmup, repeat);
//     std::cout << result2[0] << "," << duration_eigen << "," << duration_taco << "," << duration_eigen / duration_taco << std::endl;
//     std::ofstream outfile;
//     outfile.open(result_name, std::ios_base::app);
//     outfile << result2[0] << "," << duration_eigen << "," << duration_taco << "," << duration_eigen / duration_taco << std::endl;
//     outfile.close();
// }

// void benchmark_eigen_map(const string filename1, const string filename2, const int repeat, const int w_cap) {
//     taco_tensor_t C;
//     clock_t start, finish;
//     double duration_eigen, duration_taco;
//     const int warmup = 5;

//     vector<string> result;
//     boost::split(result, filename1, boost::is_any_of("/"));
//     vector<string> result2;
//     boost::split(result2, result[result.size()-1], boost::is_any_of("."));

//     double duration_eigen = CSR_CSR_T_Eigen(filename1, filename2, C_true, warmup, repeat);
//     duration_taco = CSR_CSR_T_map(filename1, filename2, &C, w_cap, warmup, repeat);
//     std::cout << result2[0] << "," << duration_eigen << "," << duration_taco << "," << duration_eigen / duration_taco << std::endl;
// }

// void benchmark_eigen_hash_outer(const string filename1, const string filename2, const int repeat, int w_cap, const string& result_name, pid_t pid) {
//     taco_tensor_t C;
//     EigenCSR C_true;
//     vector<int> indptr;
//     vector<int> indices;
//     vector<int> id_buffer;
//     vector<float> value;
//     int nrow;
//     int ncol;
//     int nnz;
//     read_mtx_csr(filename1.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
//     EigenCSC A_true = to_eigen_csc(nrow, ncol, nnz, id_buffer, indices, value);
//     // print_eigen_csc(A_true);
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     taco_tensor_t A = DC_to_taco_tensor(A_true.outerIndexPtr(),A_true.innerIndexPtr(),A_true.valuePtr(),nrow,ncol,nnz,{1,0});
//     //print_taco_tensor_DC(&A);
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     read_mtx_csr(filename2.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
//     taco_tensor_t B = DC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
//     EigenCSR B_true = to_eigen_csr(nrow, ncol, nnz, id_buffer, indices, value, false);
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
    

//     // sleep 2s
//     std::this_thread::sleep_for(std::chrono::milliseconds(100));
//     init_taco_tensor_DC(&C, ncol, ncol, {0,1});

//     clock_t start, finish;
//     double duration_eigen, duration_taco;
//     const int warmup = 5;

//     vector<string> result;
//     boost::split(result, filename1, boost::is_any_of("/"));
//     vector<string> result2;
//     boost::split(result2, result[result.size()-1], boost::is_any_of("."));

//     w_cap = pow(2,int(log2(nnz))); // heuristic
//     duration_taco = CSC_CSR_T_hash(&A, &B, &C, w_cap, warmup, repeat, true);
//     // delete_taco_tensor_DC(&C);
//     // delete_taco_tensor_DC(&A);
//     // delete_taco_tensor_DC(&B);
//     // sleep 2s
//     std::this_thread::sleep_for(std::chrono::milliseconds(2000));

//     duration_eigen = CSC_CSR_T_Eigen(A_true, B_true, C_true, warmup, repeat, true);
//     std::cout << result2[0] << "," << duration_eigen << "," << duration_taco << "," << duration_eigen / duration_taco << "," << pid << std::endl;
//     std::ofstream outfile;
//     outfile.open(result_name, std::ios_base::app);
//     outfile << result2[0] << "," << duration_eigen << "," << duration_taco << "," << duration_eigen / duration_taco << std::endl;
//     outfile.close();
// }

// void benchmark_eigen_coord_outer(const string filename1, const string filename2, const int repeat, int w_cap, const string& result_name) {
//     taco_tensor_t C;
//     EigenCSR C_true;
//     vector<int> indptr;
//     vector<int> indices;
//     vector<int> id_buffer;
//     vector<float> value;
//     int nrow;
//     int ncol;
//     int nnz;
//     read_mtx_csr(filename1.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
//     EigenCSC A_true = to_eigen_csc(nrow, ncol, nnz, id_buffer, indices, value);
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     taco_tensor_t A = DC_to_taco_tensor(A_true.outerIndexPtr(),A_true.innerIndexPtr(),A_true.valuePtr(),nrow,ncol,nnz,{1,0});
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     read_mtx_csr(filename2.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
//     taco_tensor_t B = DC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
//     EigenCSR B_true = to_eigen_csr(nrow, ncol, nnz, id_buffer, indices, value, false);
//     init_taco_tensor_DC(&C, ncol, ncol, {0,1});
    
//     clock_t start, finish;
//     double duration_eigen, duration_taco;
//     const int warmup = 5;

//     vector<string> result;
//     boost::split(result, filename1, boost::is_any_of("/"));
//     vector<string> result2;
//     boost::split(result2, result[result.size()-1], boost::is_any_of("."));


//     duration_eigen = CSC_CSR_T_Eigen(A_true, B_true, C_true, warmup, repeat);

//     w_cap = pow(2,int(log2(nnz))); // heuristic
//     for (int i = 0; i < warmup; i++) {
//         CSC_CSR_T_coord(&A, &B, &C, w_cap);
//     }
//     start = clock();
//     for (int i = 0; i < repeat; i++) { 
//         CSC_CSR_T_coord(&A, &B, &C, w_cap);
//     }
//     finish = clock();
//     duration_taco = (double)(finish - start) / (CLOCKS_PER_SEC * repeat);
//     std::cout << result2[0] << "," << duration_eigen << "," << duration_taco << "," << duration_eigen / duration_taco << std::endl;
//     std::ofstream outfile;
//     outfile.open(result_name, std::ios_base::app);
//     outfile << result2[0] << "," << duration_eigen << "," << duration_taco << "," << duration_eigen / duration_taco << std::endl;
//     outfile.close();
// }

// void benchmark_eigen_map_outer(const string filename1, const string filename2, const int repeat, const string& result_name, pid_t pid) {
//     taco_tensor_t C;
//     EigenCSR C_true;
//     vector<int> indptr;
//     vector<int> indices;
//     vector<int> id_buffer;
//     vector<float> value;
//     int nrow;
//     int ncol;
//     int nnz;
//     read_mtx_csr(filename1.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
//     EigenCSC A_true = to_eigen_csc(nrow, ncol, nnz, id_buffer, indices, value);
//     // print_eigen_csc(A_true);
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     taco_tensor_t A = DC_to_taco_tensor(A_true.outerIndexPtr(),A_true.innerIndexPtr(),A_true.valuePtr(),nrow,ncol,nnz,{1,0});
//     //print_taco_tensor_DC(&A);
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     read_mtx_csr(filename2.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
//     taco_tensor_t B = DC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
//     EigenCSR B_true = to_eigen_csr(nrow, ncol, nnz, id_buffer, indices, value, false);
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
    

//     // sleep 2s
//     std::this_thread::sleep_for(std::chrono::milliseconds(100));
//     init_taco_tensor_DC(&C, ncol, ncol, {0,1});

//     double duration_eigen, duration_taco;
//     const int warmup = 0;

//     vector<string> result;
//     boost::split(result, filename1, boost::is_any_of("/"));
//     vector<string> result2;
//     boost::split(result2, result[result.size()-1], boost::is_any_of("."));

//     duration_taco = CSC_CSR_T_map(&A, &B, &C, warmup, repeat, true);
//     // delete_taco_tensor_DC(&C);
//     // delete_taco_tensor_DC(&A);
//     // delete_taco_tensor_DC(&B);
//     // sleep 2s
//     std::this_thread::sleep_for(std::chrono::milliseconds(2000));

//     duration_eigen = CSC_CSR_T_Eigen(A_true, B_true, C_true, warmup, repeat, true);
//     std::cout << result2[0] << "," << duration_eigen << "," << duration_taco << "," << duration_eigen / duration_taco << "," << pid << std::endl;
//     std::ofstream outfile;
//     outfile.open(result_name, std::ios_base::app);
//     outfile << result2[0] << "," << duration_eigen << "," << duration_taco << "," << duration_eigen / duration_taco << std::endl;
//     outfile.close();
// }

// void benchmark_eigen_hash_outer_CC(const string filename1, const string filename2, const int repeat, int w_cap, const string& result_name, pid_t pid) {
//     taco_tensor_t C;
//     EigenCSR C_true;
//     vector<int> indptr;
//     vector<int> indices;
//     vector<int> id_buffer;
//     vector<float> value;
//     int nrow;
//     int ncol;
//     int nnz;
//     read_mtx_csr(filename1.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
//     EigenCSC A_true = to_eigen_csc(nrow, ncol, nnz, id_buffer, indices, value);
//     // print_eigen_csc(A_true);
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     taco_tensor_t A = CC_to_taco_tensor(A_true.outerIndexPtr(),A_true.innerIndexPtr(),A_true.valuePtr(),nrow,ncol,nnz,{1,0});
//     //print_taco_tensor_DC(&A);
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     read_mtx_csr(filename2.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
//     taco_tensor_t B = CC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
//     EigenCSR B_true = to_eigen_csr(nrow, ncol, nnz, id_buffer, indices, value, false);
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
    

//     // sleep 2s
//     std::this_thread::sleep_for(std::chrono::milliseconds(100));
//     init_taco_tensor_DC(&C, ncol, ncol, {0,1});

//     clock_t start, finish;
//     double duration_eigen, duration_taco;
//     const int warmup = 5;

//     vector<string> result;
//     boost::split(result, filename1, boost::is_any_of("/"));
//     vector<string> result2;
//     boost::split(result2, result[result.size()-1], boost::is_any_of("."));

//     w_cap = pow(2,int(log2(nnz))); // heuristic
//     duration_taco = CSC_CSR_T_hash(&A, &B, &C, w_cap, warmup, repeat, true);
//     // delete_taco_tensor_DC(&C);
//     // delete_taco_tensor_DC(&A);
//     // delete_taco_tensor_DC(&B);
//     // sleep 2s
//     std::this_thread::sleep_for(std::chrono::milliseconds(2000));

//     duration_eigen = CSC_CSR_T_Eigen(A_true, B_true, C_true, warmup, repeat, true);
//     std::cout << result2[0] << "," << duration_eigen << "," << duration_taco << "," << duration_eigen / duration_taco << "," << pid << std::endl;
//     std::ofstream outfile;
//     outfile.open(result_name, std::ios_base::app);
//     outfile << result2[0] << "," << duration_eigen << "," << duration_taco << "," << duration_eigen / duration_taco << std::endl;
//     outfile.close();
// }

// void profile_coord_bucket(const string filename1, const string filename2, const int repeat, int w_cap, const string& result_name) {
//     taco_tensor_t C;
//     EigenCSR C_true;
//     vector<int> indptr;
//     vector<int> indices;
//     vector<int> id_buffer;
//     vector<float> value;
//     int nrow;
//     int ncol;
//     int nnz;
//     read_mtx_csc(filename1.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
//     EigenCSC A_true = to_eigen_csc(nrow, ncol, nnz, id_buffer, indices, value);
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     taco_tensor_t A = DC_to_taco_tensor(A_true.outerIndexPtr(),A_true.innerIndexPtr(),A_true.valuePtr(),nrow,ncol,nnz,{1,0});
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     read_mtx_csr(filename2.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
//     taco_tensor_t B = DC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
//     init_taco_tensor_DC(&C, ncol, ncol, {0,1});
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();

//     w_cap = pow(2,int(log2(nnz))); // heuristic
//     CSC_CSR_T_coord(&A, &B, &C, w_cap, repeat, result_name);
// }

// void check_mttkrp_dcsf_dcsr_dcsr_hash(const string filename1, const string filename2, const string filename3, const string filename4, const int repeat, bool verbose, vector<int>& dims) {
//     taco_tensor_t B;
//     read_tns_csf(filename1, B, dims);
//     std::cout << "B: " << B.dimensions[0] << "," << B.dimensions[1] << "," << B.dimensions[2] << "," << B.vals_size << std::endl;
//     vector<int> indptr;
//     vector<int> indices;
//     vector<int> id_buffer;
//     vector<float> value;
//     int nrow;
//     int ncol;
//     int nnz;
//     read_mtx_csr_keep_value(filename2.data(), nrow, ncol, nnz, indptr, indices, value, false /*one_base*/);
//     taco_tensor_t C = CC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
//     std::cout << "C: " << nrow << "," << ncol << "," << nnz << std::endl; // 12092,32,
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     read_mtx_csr_keep_value(filename3.data(), nrow, ncol, nnz, indptr, indices, value, false /*one_base*/);
//     taco_tensor_t D = CC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
//     std::cout << "D: " << nrow << "," << ncol << "," << nnz << std::endl; // 9184,32,
//     taco_tensor_t A;
//     init_taco_tensor_DC(&A, dims[2], ncol, {0,1});
//     std::cout << "A: " << nrow << "," << ncol << std::endl; // 28818,32

//     int w_cap = pow(2,int(log2(1.0 * nnz * dims[2] / nrow))); // heuristic
//     std::cout << "w_cap: " << w_cap << std::endl; 
//     int warmup = 1;
//     double duration_taco = COO_CSF_DCSR_DCSR_hash(&A, &B, &C, &D, w_cap, warmup, repeat, false /*in bench*/, verbose);

//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     read_mtx_coo(filename4.data(), nrow, ncol, nnz, id_buffer, indices, value, false/*one_base*/);
//     cout << "nnz: " << A.vals_size << "," << nnz << endl;
//     compare_array<int>(A.indices[1][0], id_buffer.data(), nnz);
//     compare_array<int>(A.indices[1][1], indices.data(), nnz);
//     compare_array<float>(A.vals, value.data(), nnz);
// }

// void bench_mttkrp_dcsf_dcsr_dcsr_hash(const string filename1, const string filename2, const string filename3, const int repeat, bool verbose, vector<int>& dims) {
//     taco_tensor_t B;
//     read_tns_csf(filename1, B, dims);
//     vector<int> indptr;
//     vector<int> indices;
//     vector<int> id_buffer;
//     vector<float> value;
//     int nrow;
//     int ncol;
//     int nnz;
//     read_mtx_csr(filename2.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value, false);
//     taco_tensor_t C = CC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     read_mtx_csr(filename3.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value, false);
//     taco_tensor_t D = CC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
//     taco_tensor_t A;
//     init_taco_tensor_DC(&A, dims[2], ncol, {0,1});

//     vector<string> result;
//     boost::split(result, filename1, boost::is_any_of("/"));
//     vector<string> result2;
//     boost::split(result2, result[result.size()-1], boost::is_any_of("."));

//     int w_cap = pow(2,int(log2(nnz * dims[2] / nrow))); // heuristic
//     int warmup = 5;
//     double duration_taco = COO_CSF_DCSR_DCSR_hash(&A, &B, &C, &D, w_cap, warmup, repeat, true /* in bench*/, verbose);
    
//     std::cout << result2[0] << "," << ncol << "," << duration_taco * 1000 << std::endl; // ms
// }

// void bench_mttkrp_splatt(const string& filename1, int ncol, const int repeat) {
//     // OS Error: Invalid instruction
//     vector<string> result;
//     boost::split(result, filename1, boost::is_any_of("/"));
//     vector<string> result2;
//     boost::split(result2, result[result.size()-1], boost::is_any_of("."));

//     int warmup = 1;
//     double duration = mttkrp_splatt(filename1, ncol, warmup, repeat);
//     std::cout << result2[0] << "," << ncol << "," << duration * 1000 << std::endl; // ms
// }

// void check_ttm_dcsf_dcsr_hash(const string filename1, const string filename2, const string filename3, const int repeat, bool verbose, vector<int>& dims) {
//     taco_tensor_t B;
//     read_tns_csf(filename1, B, dims);
//     std::cout << "B: " << B.dimensions[0] << "," << B.dimensions[1] << "," << B.dimensions[2] << "," << B.vals_size << std::endl;
//     vector<int> indptr;
//     vector<int> indices;
//     vector<int> id_buffer;
//     vector<float> value;
//     int nrow;
//     int ncol;
//     int nnz;
//     read_mtx_csr_keep_value(filename2.data(), nrow, ncol, nnz, indptr, indices, value, false /*one_base*/);
//     taco_tensor_t C = CC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
//     std::cout << "C: " << C.dimensions[0] << "," << C.dimensions[1] << "," << C.vals_size << std::endl; // 12092,32,
//     //print_taco_tensor_CC(&C);
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     taco_tensor_t A;
//     init_taco_tensor_COO(&A, {dims[1],dims[2],ncol}, {0,1,2});
//     std::cout << "A: " << A.dimensions[0] << "," << A.dimensions[1] << "," << A.dimensions[2] << std::endl; // 28818,32

//     int w_cap = pow(2,int(log2(1.0 * dims[1] * dims[2]))); // heuristic
//     std::cout << "w_cap: " << w_cap << std::endl; 
//     int warmup = 1;
//     double duration_taco = COO_CSF_DCSR_hash(&A, &B, &C, w_cap, warmup, repeat, false /*in bench*/, verbose);
//     std::cout << "duration_taco: " << duration_taco * 1000 << " ms" << std::endl;
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     read_txt_coo3(filename3.data(), indptr, id_buffer, indices, value, false/*one_base*/);
//     nnz = indptr.size();
//     cout << "nnz: " << A.vals_size << "," << nnz << endl;
//     compare_array<int>(A.indices[1][0], indptr.data(), nnz);
//     compare_array<int>(A.indices[2][0], id_buffer.data(), nnz);
//     compare_array<int>(A.indices[3][0], indices.data(), nnz);
//     compare_array<float>(A.vals, value.data(), nnz);
// }

// void bench_ttm_dcsf_dcsr_hash(const string filename1, const string filename2, const int repeat, vector<int>& dims) {
//     taco_tensor_t B;
//     read_tns_csf(filename1, B, dims);
//     vector<int> indptr;
//     vector<int> indices;
//     vector<int> id_buffer;
//     vector<float> value;
//     int nrow;
//     int ncol;
//     int nnz;
//     read_mtx_csr_keep_value(filename2.data(), nrow, ncol, nnz, indptr, indices, value, false /*one_base*/);
//     taco_tensor_t C = CC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
//     //print_taco_tensor_CC(&C);
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     taco_tensor_t A;
//     init_taco_tensor_COO(&A, {dims[1],dims[2],ncol}, {0,1,2});

//     vector<string> result;
//     boost::split(result, filename1, boost::is_any_of("/"));
//     vector<string> result2;
//     boost::split(result2, result[result.size()-1], boost::is_any_of("."));

//     int w_cap = pow(2,int(log2(1.0 * dims[1] * dims[2]))); // heuristic
//     int warmup = 5;
//     double duration_taco = COO_CSF_DCSR_hash(&A, &B, &C, w_cap, warmup, repeat, true/*in bench*/);
//     std::cout << result2[0] << "," << ncol << "," << duration_taco * 1000 << std::endl; // ms

// }

// void check_eigen_hash_outer_nonfuse(const string filename1, const string filename2, const int repeat, int w_cap, bool verbose) {
//     taco_tensor_t C;
//     EigenCSR C_true;
//     vector<int> indptr;
//     vector<int> indices;
//     vector<int> id_buffer;
//     vector<float> value;
//     int nrow;
//     int ncol;
//     int nnz;
//     read_mtx_csr(filename1.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
//     EigenCSC A_true = to_eigen_csc(nrow, ncol, nnz, id_buffer, indices, value);
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     taco_tensor_t A = DC_to_taco_tensor(A_true.outerIndexPtr(),A_true.innerIndexPtr(),A_true.valuePtr(),nrow,ncol,nnz,{1,0});
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     read_mtx_csr(filename2.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
//     taco_tensor_t B = DC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
//     EigenCSR B_true = to_eigen_csr(nrow, ncol, nnz, id_buffer, indices, value, false);
//     init_taco_tensor_DC(&C, ncol, ncol, {0,1});
//     taco_tensor_t C_noT;
//     init_taco_tensor_DC(&C_noT, ncol, ncol, {0,1});
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     //print_taco_tensor_DC(&A);
//     //print_taco_tensor_DC(&B);
//     // Check result


//     std::cout << "Hash" << std::endl;
//     w_cap = pow(2,int(log2(nnz))); // heuristic
//     //std::cout << "w_cap: " << w_cap << std::endl;  
//     //w_cap = 3;
//     CSC_CSR_T_hash_nonfuse(&A, &B, &C_noT, &C, w_cap, 0, 1, false, verbose);
//     std::cout << "Eigen" << std::endl;
//     CSC_CSR_T_Eigen(A_true, B_true, C_true, 0, 1, false, verbose);
//     int outSize = C_true.outerSize();
//     std::cout << "Filename " << filename1 << std::endl;
//     std::cout << "Output Nnz: " << C_true.outerIndexPtr()[outSize] << std::endl;
//     check_csr_taco_eigen(C, C_true);
// }

void bench_eigen_hash_outer_nonfuse(const string filename1, const string filename2, const int repeat, int w_cap, const string& result_name, pid_t pid) {
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
    EigenCSR B_true = to_eigen_csr(nrow, ncol, nnz, id_buffer, indices, value, false);
    init_taco_tensor_DC(&C, ncol, ncol, {0,1});
    taco_tensor_t C_noT;
    init_taco_tensor_DC(&C_noT, ncol, ncol, {0,1});
    indptr.clear();
    id_buffer.clear();
    indices.clear();
    value.clear();

    clock_t start, finish;
    double duration_eigen, duration_taco;
    const int warmup = 5;

    vector<string> result;
    boost::split(result, filename1, boost::is_any_of("/"));
    vector<string> result2;
    boost::split(result2, result[result.size()-1], boost::is_any_of("."));

    w_cap = pow(2,int(log2(nnz))); // heuristic
    duration_taco = CSC_CSR_T_hash_nonfuse(&A, &B, &C_noT, &C, w_cap, warmup, repeat, true);
    duration_eigen = CSC_CSR_T_Eigen(A_true, B_true, C_true, warmup, repeat, true);
    std::cout << result2[0] << "," << duration_eigen << "," << duration_taco << "," << duration_eigen / duration_taco << "," << pid << std::endl;
    std::ofstream outfile;
    outfile.open(result_name, std::ios_base::app);
    outfile << result2[0] << "," << duration_eigen << "," << duration_taco << "," << duration_eigen / duration_taco << std::endl;
    outfile.close();
}

// void check_eigen_coord_outer_CC(const string filename1, const string filename2, const int repeat, int w_cap, bool verbose) {
//     taco_tensor_t C;
//     EigenCSR C_true;
//     vector<int> indptr;
//     vector<int> indices;
//     vector<int> id_buffer;
//     vector<float> value;
//     int nrow;
//     int ncol;
//     int nnz;
//     read_mtx_csr(filename1.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
//     EigenCSC A_true = to_eigen_csc(nrow, ncol, nnz, id_buffer, indices, value);
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     taco_tensor_t A = CC_to_taco_tensor(A_true.outerIndexPtr(),A_true.innerIndexPtr(),A_true.valuePtr(),nrow,ncol,nnz,{1,0});
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     read_mtx_csr(filename2.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
//     taco_tensor_t B = CC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
//     EigenCSR B_true = to_eigen_csr(nrow, ncol, nnz, id_buffer, indices, value, false);
//     init_taco_tensor_DC(&C, ncol, ncol, {0,1});

//     std::cout << "Eigen" << std::endl;
//     CSC_CSR_T_Eigen(A_true, B_true, C_true, 0, 1, false, verbose);
//     int outSize = C_true.outerSize();
//     std::cout << "Filename " << filename1 << std::endl;
//     //std::cout << "Output sparsity " << double(C_true.outerIndexPtr()[outSize]) / double(outSize * outSize) << std::endl; // overflow
//     std::cout << "Output Nnz: " << C_true.outerIndexPtr()[outSize] << std::endl;

//     std::cout << "Coord" << std::endl;
//     w_cap = pow(2,int(log2(nnz))); // heuristic
//     CSC_CSR_T_coord(&A, &B, &C, w_cap, verbose);
//     check_csr_taco_eigen(C, C_true);
// }

// void bench_eigen_coord_outer_CC(const string filename1, const string filename2, const int repeat, int w_cap, const string& result_name) {
//     taco_tensor_t C;
//     EigenCSR C_true;
//     vector<int> indptr;
//     vector<int> indices;
//     vector<int> id_buffer;
//     vector<float> value;
//     int nrow;
//     int ncol;
//     int nnz;
//     read_mtx_csr(filename1.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
//     EigenCSC A_true = to_eigen_csc(nrow, ncol, nnz, id_buffer, indices, value);
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     taco_tensor_t A = CC_to_taco_tensor(A_true.outerIndexPtr(),A_true.innerIndexPtr(),A_true.valuePtr(),nrow,ncol,nnz,{1,0});
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     read_mtx_csr(filename2.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
//     taco_tensor_t B = CC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
//     EigenCSR B_true = to_eigen_csr(nrow, ncol, nnz, id_buffer, indices, value, false);
//     init_taco_tensor_DC(&C, ncol, ncol, {0,1});

//     clock_t start, finish;
//     double duration_eigen, duration_taco;
//     const int warmup = 5;

//     vector<string> result;
//     boost::split(result, filename1, boost::is_any_of("/"));
//     vector<string> result2;
//     boost::split(result2, result[result.size()-1], boost::is_any_of("."));

//     duration_eigen = CSC_CSR_T_Eigen(A_true, B_true, C_true, warmup, repeat);

//     w_cap = pow(2,int(log2(nnz))); // heuristic
//     for (int i = 0; i < warmup; i++) {
//         CSC_CSR_T_coord(&A, &B, &C, w_cap);
//     }
//     start = clock();
//     for (int i = 0; i < repeat; i++) { 
//         CSC_CSR_T_coord(&A, &B, &C, w_cap);
//     }
//     finish = clock();
//     duration_taco = (double)(finish - start) / (CLOCKS_PER_SEC * repeat);
//     std::cout << result2[0] << "," << duration_eigen << "," << duration_taco << "," << duration_eigen / duration_taco << std::endl;
//     std::ofstream outfile;
//     outfile.open(result_name, std::ios_base::app);
//     outfile << result2[0] << "," << duration_eigen << "," << duration_taco << "," << duration_eigen / duration_taco << std::endl;
//     outfile.close();
// }

// void check_eigen_hash_outer_CC(const string filename1, const string filename2, const int repeat, int w_cap, bool verbose) {
//     taco_tensor_t C;
//     EigenCSR C_true;
//     vector<int> indptr;
//     vector<int> indices;
//     vector<int> id_buffer;
//     vector<float> value;
//     int nrow;
//     int ncol;
//     int nnz;
//     read_mtx_csr(filename1.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
//     EigenCSC A_true = to_eigen_csc(nrow, ncol, nnz, id_buffer, indices, value);
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     taco_tensor_t A = CC_to_taco_tensor(A_true.outerIndexPtr(),A_true.innerIndexPtr(),A_true.valuePtr(),nrow,ncol,nnz,{1,0});
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     read_mtx_csr(filename2.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
//     taco_tensor_t B = CC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
//     EigenCSR B_true = to_eigen_csr(nrow, ncol, nnz, id_buffer, indices, value, false);
//     init_taco_tensor_DC(&C, ncol, ncol, {0,1});
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     //print_taco_tensor_DC(&AD);
//     //print_taco_tensor_CC(&A);
//     //print_taco_tensor_CC(&B);
//     // Check result


//     std::cout << "Hash" << std::endl;
//     w_cap = pow(2,int(log2(nnz))); // heuristic
//     //w_cap = 3;
//     CSC_CSR_T_hash(&A, &B, &C, w_cap, 0, 1, false, verbose);
//     std::cout << "Eigen" << std::endl;
//     CSC_CSR_T_Eigen(A_true, B_true, C_true, 0, 1, false, verbose);
//     int outSize = C_true.outerSize();
//     std::cout << "Filename " << filename1 << std::endl;
//     std::cout << "Output Nnz: " << C_true.outerIndexPtr()[outSize] << std::endl;
//     check_csr_taco_eigen(C, C_true);
// }

// void bench_eigen_hash_outer_CC(const string filename1, const string filename2, const int repeat, int w_cap, const string& result_name, pid_t pid) {
//     taco_tensor_t C;
//     EigenCSR C_true;
//     vector<int> indptr;
//     vector<int> indices;
//     vector<int> id_buffer;
//     vector<float> value;
//     int nrow;
//     int ncol;
//     int nnz;
//     read_mtx_csr(filename1.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
//     EigenCSC A_true = to_eigen_csc(nrow, ncol, nnz, id_buffer, indices, value);
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     taco_tensor_t A = CC_to_taco_tensor(A_true.outerIndexPtr(),A_true.innerIndexPtr(),A_true.valuePtr(),nrow,ncol,nnz,{1,0});
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     read_mtx_csr(filename2.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
//     taco_tensor_t B = CC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
//     EigenCSR B_true = to_eigen_csr(nrow, ncol, nnz, id_buffer, indices, value, false);
//     init_taco_tensor_DC(&C, ncol, ncol, {0,1});
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();

//     clock_t start, finish;
//     double duration_eigen, duration_taco;
//     const int warmup = 5;

//     vector<string> result;
//     boost::split(result, filename1, boost::is_any_of("/"));
//     vector<string> result2;
//     boost::split(result2, result[result.size()-1], boost::is_any_of("."));

//     w_cap = pow(2,int(log2(nnz))); // heuristic
//     duration_taco = CSC_CSR_T_hash(&A, &B, &C, w_cap, warmup, repeat, true);
//     duration_eigen = CSC_CSR_T_Eigen(A_true, B_true, C_true, warmup, repeat, true);
//     std::cout << result2[0] << "," << duration_eigen << "," << duration_taco << "," << duration_eigen / duration_taco << "," << pid << std::endl;
//     std::ofstream outfile;
//     outfile.open(result_name, std::ios_base::app);
//     outfile << result2[0] << "," << duration_eigen << "," << duration_taco << "," << duration_eigen / duration_taco << std::endl;
//     outfile.close();
// }

// void bench_eigen_hash_outer_CC_noT(const string filename1, const string filename2, const int repeat, int w_cap, const string& result_name, pid_t pid) {
//     taco_tensor_t C;
//     EigenCSR C_true;
//     vector<int> indptr;
//     vector<int> indices;
//     vector<int> id_buffer;
//     vector<float> value;
//     int nrow;
//     int ncol;
//     int nnz;
//     read_mtx_csr(filename1.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
//     EigenCSC A_true = to_eigen_csc(nrow, ncol, nnz, id_buffer, indices, value);
//     EigenCSR A_use = to_eigen_csr(nrow, ncol, nnz, id_buffer, indices, value);
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     taco_tensor_t A = CC_to_taco_tensor(A_true.outerIndexPtr(),A_true.innerIndexPtr(),A_true.valuePtr(),nrow,ncol,nnz,{1,0});
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     read_mtx_csr(filename2.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
//     taco_tensor_t B = CC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
//     EigenCSR B_true = to_eigen_csr(nrow, ncol, nnz, id_buffer, indices, value, false);
//     init_taco_tensor_DC(&C, ncol, ncol, {0,1});
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();

//     clock_t start, finish;
//     double duration_eigen, duration_taco;
//     const int warmup = 5;

//     vector<string> result;
//     boost::split(result, filename1, boost::is_any_of("/"));
//     vector<string> result2;
//     boost::split(result2, result[result.size()-1], boost::is_any_of("."));

//     w_cap = pow(2,int(log2(nnz))); // heuristic
//     duration_taco = DCSC_DCSR_hash(&A, &B, &C, w_cap, warmup, repeat, true);
//     duration_eigen = CSR_CSR_Eigen(A_use, B_true, C_true, warmup, repeat, true);
//     std::cout << result2[0] << "," << duration_eigen << "," << duration_taco << "," << duration_eigen / duration_taco << "," << pid << std::endl;
//     std::ofstream outfile;
//     outfile.open(result_name, std::ios_base::app);
//     outfile << result2[0] << "," << duration_eigen << "," << duration_taco << "," << duration_eigen / duration_taco << std::endl;
//     outfile.close();
// }

// void bench_eigen_coord_outer_CC_noT(const string filename1, const string filename2, const int repeat, int w_cap, const string& result_name, pid_t pid) {
//     taco_tensor_t C;
//     EigenCSR C_true;
//     vector<int> indptr;
//     vector<int> indices;
//     vector<int> id_buffer;
//     vector<float> value;
//     int nrow;
//     int ncol;
//     int nnz;
//     read_mtx_csr(filename1.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
//     EigenCSC A_true = to_eigen_csc(nrow, ncol, nnz, id_buffer, indices, value);
//     EigenCSR A_use = to_eigen_csr(nrow, ncol, nnz, id_buffer, indices, value);
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     taco_tensor_t A = CC_to_taco_tensor(A_true.outerIndexPtr(),A_true.innerIndexPtr(),A_true.valuePtr(),nrow,ncol,nnz,{1,0});
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     read_mtx_csr(filename2.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
//     taco_tensor_t B = CC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
//     EigenCSR B_true = to_eigen_csr(nrow, ncol, nnz, id_buffer, indices, value, false);
//     init_taco_tensor_DC(&C, ncol, ncol, {0,1});
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();

//     clock_t start, finish;
//     double duration_eigen, duration_taco;
//     const int warmup = 5;

//     vector<string> result;
//     boost::split(result, filename1, boost::is_any_of("/"));
//     vector<string> result2;
//     boost::split(result2, result[result.size()-1], boost::is_any_of("."));

//     w_cap = pow(2,int(log2(nnz))); // heuristic
//     for (int i = 0; i < warmup; i++) {
//         DCSC_DCSR_coord(&A, &B, &C, w_cap);
//     }
//     start = clock();
//     for (int i = 0; i < repeat; i++) { 
//         DCSC_DCSR_coord(&A, &B, &C, w_cap);
//     }
//     finish = clock();
//     duration_taco = (double)(finish - start) / (CLOCKS_PER_SEC * repeat);
//     duration_eigen = CSR_CSR_Eigen(A_use, B_true, C_true, warmup, repeat, true);
//     std::cout << result2[0] << "," << duration_eigen << "," << duration_taco << "," << duration_eigen / duration_taco << "," << pid << std::endl;
//     std::ofstream outfile;
//     outfile.open(result_name, std::ios_base::app);
//     outfile << result2[0] << "," << duration_eigen << "," << duration_taco << "," << duration_eigen / duration_taco << std::endl;
//     outfile.close();
// }

// void check_eigen_hash_outer_CC_noT(const string filename1, const string filename2, const int repeat, int w_cap, bool verbose) {
//     taco_tensor_t C;
//     EigenCSR C_true;
//     vector<int> indptr;
//     vector<int> indices;
//     vector<int> id_buffer;
//     vector<float> value;
//     int nrow;
//     int ncol;
//     int nnz;
//     read_mtx_csr(filename1.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
//     EigenCSC A_true = to_eigen_csc(nrow, ncol, nnz, id_buffer, indices, value);
//     EigenCSR A_use = to_eigen_csr(nrow, ncol, nnz, id_buffer, indices, value);
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     taco_tensor_t A = CC_to_taco_tensor(A_true.outerIndexPtr(),A_true.innerIndexPtr(),A_true.valuePtr(),nrow,ncol,nnz,{1,0});
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();
//     read_mtx_csr(filename2.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
//     taco_tensor_t B = CC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
//     EigenCSR B_true = to_eigen_csr(nrow, ncol, nnz, id_buffer, indices, value, false);
//     init_taco_tensor_DC(&C, ncol, ncol, {0,1});
//     indptr.clear();
//     id_buffer.clear();
//     indices.clear();
//     value.clear();

//     std::cout << "Hash" << std::endl;
//     w_cap = pow(2,int(log2(nnz))); // heuristic
//     DCSC_DCSR_hash(&A, &B, &C, w_cap, 0, 1, false);
//     std::cout << "Eigen" << std::endl;
//     CSR_CSR_Eigen(A_use, B_true, C_true, 0, 1, false);
//     int outSize = C_true.outerSize();
//     std::cout << "Filename " << filename1 << std::endl;
//     std::cout << "Output Nnz: " << C_true.outerIndexPtr()[outSize] << std::endl;
//     check_csr_taco_eigen(C, C_true);
// }

int main(int argc, char* argv[]) {
    const string filename1 = (argc > 1) ? argv[1] : "./data/test5.mtx";
    const string filename2 = (argc > 2) ? argv[2] : "./data/test5.mtx";
    const int repeat = (argc > 3) ? stoi(argv[3]) : 10;
    const int w_cap = (argc > 4) ? stoi(argv[4]) : 16;
    const int verbose = (argc > 5) ? stoi(argv[5]) : 0;
    const string result_name = (argc > 6) ? argv[6] : "./data/test.csv";
    const string tensorDims = (argc > 7) ? argv[7] : "0,1";
    const string gold_name = (argc > 8) ? argv[8] : "./data/test_gold.mtx";

    pid_t pid = getpid();

    auto dimsStr = split(tensorDims, ",", false /* keepDelim */);
    std::vector<int> dims;
    for (auto it : dimsStr) {
        dims.push_back(atoi(it.c_str()));
    }

    //get_output_nnz(filename1, filename2, result_name);
    //assemble_hash(filename1, filename2, repeat, w_cap, verbose);

    // benchmark_eigen_coord(filename1, filename2, repeat, w_cap, result_name);
    //benchmark_eigen_hash(filename1, filename2, repeat, w_cap, result_name);
    //benchmark_eigen_taco(filename1, filename2, repeat, result_name);
    //benchmark_eigen_map(filename1, filename2, repeat, w_cap);
    //benchmark_eigen_hash_outer(filename1, filename2, repeat, w_cap, result_name, pid);
    //benchmark_eigen_coord_outer(filename1, filename2, repeat, w_cap, result_name);
    //benchmark_eigen_map_outer(filename1, filename2, repeat, result_name, pid);
    //bench_mttkrp_dcsf_dcsr_dcsr_hash(filename1, filename2, result_name, repeat, verbose, dims); // A(i,j) = B(k,l,i) * C(k,j) * D(l,j);
    //bench_mttkrp_splatt(filename1, w_cap, repeat); // w_cap is ncol here
    //bench_eigen_coord_outer_CC(filename1, filename2, repeat, w_cap, result_name);
    //bench_eigen_hash_outer_CC(filename1, filename2, repeat, w_cap, result_name, pid);
    bench_eigen_hash_outer_nonfuse(filename1, filename2, repeat, w_cap, result_name, pid);
    //bench_ttm_dcsf_dcsr_hash(filename1, filename2, repeat, dims);
    //bench_eigen_hash_outer_CC_noT(filename1, filename2, repeat, w_cap, result_name, pid);
    //bench_eigen_coord_outer_CC_noT(filename1, filename2, repeat, w_cap, result_name, pid);



    //check_eigen_hash(filename1, filename2, repeat, w_cap, verbose);
    //check_eigen_coord(filename1, filename2, repeat, w_cap, verbose);
    //check_eigen_taco(filename1, filename2, repeat, verbose);
    //check_eigen_map(filename1, filename2, repeat, w_cap, verbose);
    //check_eigen_hash_outer(filename1, filename2, repeat, w_cap, verbose);
    //check_eigen_coord_outer(filename1, filename2, repeat, w_cap, verbose);
    //check_eigen_map_outer(filename1, filename2, repeat, verbose);
    //check_eigen_hash_outer_CC(filename1, filename2, repeat, w_cap, verbose);
    //check_eigen_coord_outer_CC(filename1, filename2, repeat, w_cap, verbose);
    //check_mttkrp_dcsf_dcsr_dcsr_hash(filename1, filename2, result_name, gold_name, repeat, verbose, dims); // A(i,j) = B(k,l,i) * C(k,j) * D(l,j);
    //check_eigen_hash_outer_nonfuse(filename1, filename2, repeat, w_cap, verbose);
    //check_ttm_dcsf_dcsr_hash(filename1, filename2, result_name, repeat, verbose, dims);
    // check_eigen_hash_outer_CC_noT(filename1, filename2, repeat, w_cap, verbose);


    //profile_coord_bucket(filename1, filename2, repeat, w_cap, result_name);

    return 0;
}