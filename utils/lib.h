#ifndef LIB_H
#define LIB_H
#define _OPENMP 0
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <cstdio>
#include <climits>
#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <utility>
#include <string>
#include <vector>
#include <cmath>
#include <limits>
#include <Eigen/Sparse>

#define TACO_MIN(_a,_b) ((_a) < (_b) ? (_a) : (_b))
#define TACO_MAX(_a,_b) ((_a) > (_b) ? (_a) : (_b))
#define TACO_DEREF(_a) (((___context___*)(*__ctx__))->_a)
#define restrict __restrict__


typedef enum { taco_mode_dense, taco_mode_sparse } taco_mode_t;

struct taco_tensor_t {
  int32_t      order;         // tensor order (number of modes)
  int32_t*     dimensions;    // tensor dimensions
  int32_t      csize;         // component size
  int32_t*     mode_ordering; // mode storage ordering
  taco_mode_t* mode_types;    // mode storage types
  int32_t***   indices;       // tensor index data (per mode)
  float*     vals;          // tensor values
  int32_t      vals_size;     // values array size
};


int omp_get_thread_num() { return 0; }
int omp_get_max_threads() { return 1; }

typedef Eigen::SparseMatrix<float,Eigen::RowMajor,int> EigenCSR;
typedef Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> EigenRowMajor;
typedef Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> EigenColMajor;
typedef Eigen::Matrix<float,Eigen::Dynamic,1> EigenVector;

struct COO {
  int     m, n, nnz;
  int*    rows = nullptr;
  int*    cols = nullptr;
  float* vals = nullptr;
};

void deinit_taco_tensor_t(taco_tensor_t* t) {
  for (int i = 0; i < t->order; i++) {
    free(t->indices[i]);
  }
  free(t->indices);
  free(t->dimensions);
  free(t->mode_ordering);
  free(t->mode_types);
  free(t);
}
int cmp(const void *a, const void *b) {
  return *((const int*)a) - *((const int*)b);
}
int taco_binarySearchAfter(int *array, int arrayStart, int arrayEnd, int target) {
  if (array[arrayStart] >= target) {
    return arrayStart;
  }
  int lowerBound = arrayStart; // always < target
  int upperBound = arrayEnd; // always >= target
  while (upperBound - lowerBound > 1) {
    int mid = (upperBound + lowerBound) / 2;
    int midValue = array[mid];
    if (midValue < target) {
      lowerBound = mid;
    }
    else if (midValue > target) {
      upperBound = mid;
    }
    else {
      return mid;
    }
  }
  return upperBound;
}
int taco_binarySearchBefore(int *array, int arrayStart, int arrayEnd, int target) {
  if (array[arrayEnd] <= target) {
    return arrayEnd;
  }
  int lowerBound = arrayStart; // always <= target
  int upperBound = arrayEnd; // always > target
  while (upperBound - lowerBound > 1) {
    int mid = (upperBound + lowerBound) / 2;
    int midValue = array[mid];
    if (midValue < target) {
      lowerBound = mid;
    }
    else if (midValue > target) {
      upperBound = mid;
    }
    else {
      return mid;
    }
  }
  return lowerBound;
}

EigenRowMajor gen_row_major_matrix(int rows, int cols) {
  return EigenRowMajor::Ones(rows, cols);
  //return EigenRowMajor::Random(rows, cols) + 200 * EigenRowMajor::Ones(rows, cols);
}
EigenColMajor gen_col_major_matrix(int rows, int cols) {
  return EigenColMajor::Random(rows, cols) + EigenColMajor::Ones(rows, cols);
}
EigenVector gen_vector(int size) {
  return EigenVector::Random(size) + 2.0 * EigenVector::Ones(size);
}
EigenVector gen_one_vector(int size) {
    return EigenVector::Ones(size);
}

EigenCSR to_eigen_csr(const COO matrix, int shift = 0) {
  EigenCSR dst(matrix.m, matrix.n);
  std::vector<Eigen::Triplet<float>> tripletList;
  tripletList.reserve(matrix.nnz);
  if (shift == 0) {
    for (size_t i = 0; i < matrix.nnz; i++) {
      tripletList.push_back({matrix.rows[i], matrix.cols[i], matrix.vals[i]});
    }
  } else {
    for (size_t i = 0; i < matrix.nnz; i++) {
      const int col = (matrix.cols[i] + shift) % matrix.n;
      tripletList.push_back({matrix.rows[i], col, matrix.vals[i]});
    }
  }
  dst.setFromTriplets(tripletList.begin(), tripletList.end());
  dst.makeCompressed();
  return dst;
}


taco_tensor_t to_taco_tensor(const EigenVector& vector) {
  taco_tensor_t vt;

  vt.dimensions = new int32_t[1];
  vt.dimensions[0] = vector.innerSize();
  //vt.vals = (int32_t*)vector.data();
  vt.vals = (float*)vector.data();
  return vt;
}
taco_tensor_t to_taco_tensor(const EigenRowMajor& matrix) {
  taco_tensor_t mt;

  mt.dimensions = new int32_t[2];
  mt.dimensions[0] = matrix.rows();
  mt.dimensions[1] = matrix.cols();
  //mt.vals = (int32_t*)matrix.data();
  mt.vals = (float*)matrix.data();
  int32_t order[2] = {0,1};
  mt.mode_ordering = order;

  return mt;
}

taco_tensor_t to_taco_tensor(const EigenColMajor& matrix) {
  taco_tensor_t vt;

  vt.dimensions = new int32_t[2];
  vt.dimensions[0] = matrix.rows();
  vt.dimensions[1] = matrix.cols();
  //vt.vals = (int32_t*)matrix.data();
  vt.vals = (float*)matrix.data();
  int32_t order[2] = {1,0};
  vt.mode_ordering = order;

  return vt;
}

taco_tensor_t to_taco_tensor(const COO matrix) {
  int32_t *pos = (int32_t*)malloc(sizeof(int32_t) * 2);
  pos[0] = 0; pos[1] = matrix.nnz;

  taco_tensor_t coot;

  coot.dimensions = new int32_t[2];
  coot.dimensions[0] = matrix.m;
  coot.dimensions[1] = matrix.n;
  coot.indices = new int32_t**[2];
  coot.indices[0] = new int32_t*[2];
  coot.indices[1] = new int32_t*[2];
  coot.indices[0][0] = (int32_t*)pos;
  coot.indices[0][1] = (int32_t*)matrix.rows;
  coot.indices[1][1] = (int32_t*)matrix.cols;
  coot.vals = (float*)matrix.vals;
  return coot;
}

taco_tensor_t get_csr_taco_tensor(int rows, int cols) {
  taco_tensor_t csrt;

  csrt.dimensions = new int32_t[2];
  csrt.dimensions[0] = rows;
  csrt.dimensions[1] = cols;
  csrt.indices = new int32_t**[2];
  csrt.indices[0] = new int32_t*[2];

  return csrt;
}

taco_tensor_t to_taco_tensor(const EigenCSR& matrix) {
  taco_tensor_t csrt = get_csr_taco_tensor(matrix.rows(), matrix.cols());

  csrt.indices[0][0] = (int32_t*)matrix.outerIndexPtr();
  csrt.indices[0][1] = (int32_t*)matrix.innerIndexPtr();
  csrt.vals = (float*)matrix.valuePtr();
  return csrt;
}

float compare_array_my(const float *x, const float *y, const size_t N){
  float thres = 0.001;
  float count = 0.0;
  float col = 0;
  for (int i=0;i<N;i++){
    col = i%128;
    if(abs(x[i]-y[i])>thres){
       printf("col: %.0f, x: %.3f, y: %.3f\n",col,x[i],y[i]);
       count ++ ;
    }
    if(count>10000)
      break;
  }
  return count;
}
float compare_vectors(const taco_tensor_t a, const taco_tensor_t b) {
  assert(a.dimensions[0] == b.dimensions[0]);
  return compare_array_my((float*)a.vals, (float*)b.vals, a.dimensions[0]);
}

double compare_matrices_float(const taco_tensor_t a, const taco_tensor_t b) {
  assert(a.dimensions[0] == b.dimensions[0]);
  assert(a.dimensions[1] == b.dimensions[1]);
  const size_t N = a.dimensions[0] * a.dimensions[1];
  return compare_array_my((float*)a.vals, (float*)b.vals, N);
}

void fill_random(float array[], int size) {
  for (int i = 0; i < size; i++) {
    array[i] = (float)(std::rand() % 3) / 10;
  }
}

void fill_one(float array[], int size) {
  for (int i = 0; i < size; i++) {
    array[i] = 1.0f;
  }  
}

COO gen_sparse_vector(int N) {
  COO coo;
  int nrow = 1;
  int ncol = N;
  std::vector<int> csr_indices_buffer;
  std::vector<int> coo_rowind_buffer;
  for (int i=0;i<N;i++){
    if((i % 2)==1){
      csr_indices_buffer.push_back(i);
      coo_rowind_buffer.push_back(0);
    }
  }
  int nnz = coo_rowind_buffer.size();

  coo.vals = (float *)malloc(sizeof(float) * nnz);
  fill_random(coo.vals, nnz);
  coo.m = nrow;
  coo.n = ncol;
  coo.nnz = nnz;
  coo.rows = (int *)malloc(sizeof(int) * nnz);
  coo.cols = (int *)malloc(sizeof(int) * nnz);
  
  for (int i=0 ; i<nnz ; i++){
    coo.cols[i] = csr_indices_buffer[i];
    coo.rows[i] = coo_rowind_buffer[i];
  }
  return coo;
}

void print_coo_matrix(const COO coo, const std::string log_path) {
  std::ofstream log_file;
  log_file.open(log_path, std::ofstream::out);
  log_file << "row: " << std::endl;
  int i=0;
  for (i=0;i<coo.nnz;i++) {
    if (i==coo.nnz-1) {
      log_file<<coo.rows[i]<<std::endl;
    }
    else {
      log_file<<coo.rows[i]<<",";
    }
  }
  log_file << "col: " << std::endl;
  for (i=0;i<coo.nnz;i++) {
    if (i==coo.nnz-1) {
      log_file<<coo.cols[i]<<std::endl;
    }
    else {
      log_file<<coo.cols[i]<<",";
    }
  }
  log_file << "vals: " << std::endl;
  for (i=0;i<coo.nnz;i++) {
    if (i==coo.nnz-1) {
      log_file<<coo.vals[i]<<std::endl;
    }
    else {
      log_file<<coo.vals[i]<<",";
    }
  }  
}

void print_coo_vector(taco_tensor_t* A){
  int* restrict A1_pos = (int*)(A->indices[0][0]);
  int* restrict A1_crd = (int*)(A->indices[1][1]);
  float* restrict A_vals = (float*)(A->vals);
  std::cout<<"["<<A1_pos[0]<<","<<A1_pos[1]<<"]"<<std::endl;
  for (int i=0; i<A1_pos[1]; i++){
    if(i==0) {
      std::cout<<"["<<A1_crd[i]<<",";
    }
    else if(i==(A1_pos[1]-1)) {
      std::cout<<A1_crd[i]<<"]"<<std::endl;
    } else {
      std::cout<<A1_crd[i]<<",";
    }
  }
  for (int i=0; i<A1_pos[1]; i++){
    if(i==0) {
      std::cout<<"["<<A_vals[i]<<",";
    }
    else if(i==(A1_pos[1]-1)) {
      std::cout<<A_vals[i]<<"]"<<std::endl;
    } else {
      std::cout<<A_vals[i]<<",";
    }
  }  
}

void print_csr_vector(taco_tensor_t* A) {
  int* restrict A1_pos = (int*)(A->indices[0][0]);
  int* restrict A1_crd = (int*)(A->indices[0][1]);
  float* restrict A_vals = (float*)(A->vals);
  std::cout<<"["<<A1_pos[0]<<","<<A1_pos[1]<<"]"<<std::endl;
  for (int i=0; i<A1_pos[1]; i++){
    if(i==0) {
      std::cout<<"["<<A1_crd[i]<<",";
    }
    else if(i==(A1_pos[1]-1)) {
      std::cout<<A1_crd[i]<<"]"<<std::endl;
    } else {
      std::cout<<A1_crd[i]<<",";
    }
  }
  for (int i=0; i<A1_pos[1]; i++){
    if(i==0) {
      std::cout<<"["<<A_vals[i]<<",";
    }
    else if(i==(A1_pos[1]-1)) {
      std::cout<<A_vals[i]<<"]"<<std::endl;
    } else {
      std::cout<<A_vals[i]<<",";
    }
  } 

}

#endif