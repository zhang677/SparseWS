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
//#include "benchmark/CSR_CSR_T_hash_index.h"
#include "benchmark/CSR_CSR_T_map.h"
#include "benchmark/CSR_CSR_T_eigen.h"
#include <iostream>
#include <time.h>
#include <string>
using namespace std;

// void check_eigen_coord(const string filename1, const string filename2, const int repeat, const int w_cap) {
//     taco_tensor_t C;
//     EigenCSR C_true;
//     CSR_CSR_T_Eigen(filename1, filename2, C_true);
//     CSR_CSR_T_coord(filename1, filename2, &C, w_cap);
//     check_csr_taco_eigen(C, C_true);
// }

// void check_eigen_hash(const string filename1, const string filename2, const int repeat, const int w_cap, bool verbose) {
//     taco_tensor_t C;
//     EigenCSR C_true;
//     std::cout << "Eigen" << std::endl;
//     CSR_CSR_T_Eigen(filename1, filename2, C_true, verbose);
//     //CSR_CSR_5(filename1, filename2, &C, w_cap);
//     std::cout << "Hash" << std::endl;
//     CSR_CSR_T_hash(filename1, filename2, &C, w_cap, verbose);
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

// void benchmark_eigen_coord_verbose(const string filename1, const string filename2, const int repeat, const int w_cap) {
//     taco_tensor_t C;
//     clock_t start, finish;
//     double duration;
//     const int warmup = 5;

//     EigenCSR C_true;
//     std::cout << "Eigen" << std::endl;
//     for (int i = 0; i < warmup; i++) {
//         CSR_CSR_T_Eigen(filename1, filename2, C_true);
//     }
//     start = clock();
//     for (int i = 0; i < repeat; i++) { 
//         CSR_CSR_T_Eigen(filename1, filename2, C_true);
//     }
//     finish = clock();
//     duration = (double)(finish - start) / (CLOCKS_PER_SEC * repeat);
//     std::cout << duration << " seconds" << std::endl;

//     std::cout << "CSR_CSR_T_coord_index" << std::endl;
//     for (int i = 0; i < warmup; i++) {
//         CSR_CSR_T_coord_index_double(filename1, filename2, &C, w_cap);
//     }
//     start = clock();
//     for (int i = 0; i < repeat; i++) { 
//         CSR_CSR_T_coord_index_double(filename1, filename2, &C, w_cap);
//     }
//     finish = clock();
//     duration = (double)(finish - start) / (CLOCKS_PER_SEC * repeat);
//     std::cout << duration << " seconds" << std::endl;
// }

// void benchmark_eigen_coord(const string filename1, const string filename2, const int repeat, const int w_cap) {
//     taco_tensor_t C;
//     clock_t start, finish;
//     double duration_eigen, duration_taco;
//     const int warmup = 5;

//     vector<string> result;
//     boost::split(result, filename1, boost::is_any_of("/"));
//     vector<string> result2;
//     boost::split(result2, result[result.size()-1], boost::is_any_of("."));

//     EigenCSR C_true;
//     for (int i = 0; i < warmup; i++) {
//         CSR_CSR_T_Eigen(filename1, filename2, C_true);
//     }
//     start = clock();
//     for (int i = 0; i < repeat; i++) { 
//         CSR_CSR_T_Eigen(filename1, filename2, C_true);
//     }
//     finish = clock();
//     duration_eigen = (double)(finish - start) / (CLOCKS_PER_SEC * repeat);

//     for (int i = 0; i < warmup; i++) {
//         CSR_CSR_T_coord(filename1, filename2, &C, w_cap);
//     }
//     start = clock();
//     for (int i = 0; i < repeat; i++) { 
//         CSR_CSR_T_coord(filename1, filename2, &C, w_cap);
//     }
//     finish = clock();
//     duration_taco = (double)(finish - start) / (CLOCKS_PER_SEC * repeat);
//     std::cout << result2[0] << "," << duration_eigen << "," << duration_taco << "," << duration_eigen / duration_taco << std::endl;
// }

// void benchmark_eigen_hash(const string filename1, const string filename2, const int repeat, const int w_cap) {
//     taco_tensor_t C;
//     clock_t start, finish;
//     double duration_eigen, duration_taco;
//     const int warmup = 5;

//     vector<string> result;
//     boost::split(result, filename1, boost::is_any_of("/"));
//     vector<string> result2;
//     boost::split(result2, result[result.size()-1], boost::is_any_of("."));

//     EigenCSR C_true;
//     for (int i = 0; i < warmup; i++) {
//         CSR_CSR_T_Eigen(filename1, filename2, C_true);
//     }
//     start = clock();
//     for (int i = 0; i < repeat; i++) { 
//         CSR_CSR_T_Eigen(filename1, filename2, C_true);
//     }
//     finish = clock();
//     duration_eigen = (double)(finish - start) / (CLOCKS_PER_SEC * repeat);

//     for (int i = 0; i < warmup; i++) {
//         CSR_CSR_T_hash(filename1, filename2, &C, w_cap);
//     }
//     start = clock();
//     for (int i = 0; i < repeat; i++) { 
//         CSR_CSR_T_hash(filename1, filename2, &C, w_cap);
//     }
//     finish = clock();
//     duration_taco = (double)(finish - start) / (CLOCKS_PER_SEC * repeat);
//     std::cout << result2[0] << "," << duration_eigen << "," << duration_taco << "," << duration_eigen / duration_taco << std::endl;
// }

void benchmark_eigen_map(const string filename1, const string filename2, const int repeat, const int w_cap) {
    taco_tensor_t C;
    clock_t start, finish;
    double duration_eigen, duration_taco;
    const int warmup = 5;

    vector<string> result;
    boost::split(result, filename1, boost::is_any_of("/"));
    vector<string> result2;
    boost::split(result2, result[result.size()-1], boost::is_any_of("."));

    EigenCSR C_true;
    for (int i = 0; i < warmup; i++) {
        CSR_CSR_T_Eigen(filename1, filename2, C_true);
    }
    start = clock();
    for (int i = 0; i < repeat; i++) { 
        CSR_CSR_T_Eigen(filename1, filename2, C_true);
    }
    finish = clock();
    duration_eigen = (double)(finish - start) / (CLOCKS_PER_SEC * repeat);

    for (int i = 0; i < warmup; i++) {
        CSR_CSR_T_map(filename1, filename2, &C, w_cap);
    }
    start = clock();
    for (int i = 0; i < repeat; i++) { 
        CSR_CSR_T_map(filename1, filename2, &C, w_cap);
    }
    finish = clock();
    duration_taco = (double)(finish - start) / (CLOCKS_PER_SEC * repeat);
    std::cout << result2[0] << "," << duration_eigen << "," << duration_taco << "," << duration_eigen / duration_taco << std::endl;
}

int main(int argc, char* argv[]) {
    const string filename1 = (argc > 1) ? argv[1] : "./data/test5.mtx";
    const string filename2 = (argc > 2) ? argv[2] : "./data/test5.mtx";
    const int repeat = (argc > 3) ? stoi(argv[3]) : 10;
    const int w_cap = (argc > 4) ? stoi(argv[4]) : 16;
    const int verbose = (argc > 5) ? stoi(argv[5]) : 0;

    //benchmark_eigen_coord(filename1, filename2, repeat, w_cap);
    //benchmark_eigen_hash(filename1, filename2, repeat, w_cap);
    benchmark_eigen_map(filename1, filename2, repeat, w_cap);
    //check_eigen_hash(filename1, filename2, repeat, w_cap, verbose);
    //check_eigen_coord(filename1, filename2, repeat, w_cap);
    //check_eigen_map(filename1, filename2, repeat, w_cap, verbose);
    return 0;    
}