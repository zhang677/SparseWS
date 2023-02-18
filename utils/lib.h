#ifndef LIB_H
#define LIB_H
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
#include <omp.h>

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



typedef Eigen::SparseMatrix<float,Eigen::RowMajor,int> EigenCSR;
typedef Eigen::SparseMatrix<float,Eigen::ColMajor,int> EigenCSC;
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
/*
 * The `pack` functions convert coordinate and value arrays in COO format,
 * with nonzeros sorted lexicographically by their coordinates, to the
 * specified input format.
 *
 * The `unpack` function converts the specified output format to coordinate
 * and value arrays in COO format.
 *
 * For both, the `_COO_pos` arrays contain two elements, where the first is 0
 * and the second is the number of nonzeros in the tensor.
 */

int pack_C(taco_tensor_t *C, int* C_COO1_pos, int32_t* C_COO1_crd, int32_t* C_COO2_crd, float* C_COO_vals) {
  int C1_dimension = (int)(C->dimensions[0]);
  int* restrict C2_pos = (int*)(C->indices[1][0]);
  int* restrict C2_crd = (int*)(C->indices[1][1]);
  float* restrict C_vals = (float*)(C->vals);

  C2_pos = (int32_t*)malloc(sizeof(int32_t) * (C1_dimension + 1));
  C2_pos[0] = 0;
  for (int32_t pC2 = 1; pC2 < (C1_dimension + 1); pC2++) {
    C2_pos[pC2] = 0;
  }
  int32_t C2_crd_size = 1048576;
  C2_crd = (int32_t*)malloc(sizeof(int32_t) * C2_crd_size);
  int32_t jC = 0;
  int32_t C_capacity = 1048576;
  C_vals = (float*)malloc(sizeof(float) * C_capacity);

  int32_t iC_COO = C_COO1_pos[0];
  int32_t pC_COO1_end = C_COO1_pos[1];

  while (iC_COO < pC_COO1_end) {
    int32_t i = C_COO1_crd[iC_COO];
    int32_t C_COO1_segend = iC_COO + 1;
    while (C_COO1_segend < pC_COO1_end && C_COO1_crd[C_COO1_segend] == i) {
      C_COO1_segend++;
    }
    int32_t pC2_begin = jC;

    int32_t jC_COO = iC_COO;

    while (jC_COO < C_COO1_segend) {
      int32_t j = C_COO2_crd[jC_COO];
      float C_COO_val = C_COO_vals[jC_COO];
      jC_COO++;
      while (jC_COO < C_COO1_segend && C_COO2_crd[jC_COO] == j) {
        C_COO_val += C_COO_vals[jC_COO];
        jC_COO++;
      }
      if (C_capacity <= jC) {
        C_vals = (float*)realloc(C_vals, sizeof(float) * (C_capacity * 2));
        C_capacity *= 2;
      }
      C_vals[jC] = C_COO_val;
      if (C2_crd_size <= jC) {
        C2_crd = (int32_t*)realloc(C2_crd, sizeof(int32_t) * (C2_crd_size * 2));
        C2_crd_size *= 2;
      }
      C2_crd[jC] = j;
      jC++;
    }

    C2_pos[i + 1] = jC - pC2_begin;
    iC_COO = C_COO1_segend;
  }

  int32_t csC2 = 0;
  for (int32_t pC20 = 1; pC20 < (C1_dimension + 1); pC20++) {
    csC2 += C2_pos[pC20];
    C2_pos[pC20] = csC2;
  }

  C->indices[1][0] = (int32_t*)(C2_pos);
  C->indices[1][1] = (int32_t*)(C2_crd);
  C->vals = (float*)C_vals;
  return 0;
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

EigenCSR to_eigen_csr(int nrow, int ncol, int nnz, std::vector<int> &outer_buffer, std::vector<int> &inner_buffer, std::vector<float> &vals) {
  EigenCSR dst(nrow, ncol);
  std::vector<Eigen::Triplet<float>> tripletList;
  tripletList.reserve(nnz);
  for (size_t i = 0; i < nnz; i++) {
    tripletList.push_back({outer_buffer[i], inner_buffer[i], vals[i]});
  }
  dst.setFromTriplets(tripletList.begin(), tripletList.end());
  dst.makeCompressed();
  return dst;
}

EigenCSC to_eigen_csc(int &nrow, int &ncol, int &nnz, std::vector<int> &outer_buffer, std::vector<int> &inner_buffer, std::vector<float> &vals) {
  EigenCSC dst(nrow, ncol);
  std::vector<Eigen::Triplet<float>> tripletList;
  tripletList.reserve(nnz);
  for (size_t i = 0; i < nnz; i++) {
    tripletList.push_back({outer_buffer[i], inner_buffer[i], vals[i]});
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

template <typename Type>
void compare_array(Type* A, Type* B, int len, float thres = 0.0001) {
  int count = 0;
  for (int i = 0; i < len; i++) {
    if (abs(A[i] - B[i]) > thres) {
      count++;
      std::cout<<"id: "<< i << "," << "A_val: " << A[i] << "B_val: " << B[i] << std::endl; 
    }
    if (count > 10) {
      std::cout << "Wrong!" <<std::endl;
      return;
    }
  }
  std::cout << "Correct!" <<std::endl;
}

template <typename Type>
void print_array(Type* array, int len) {
  std::cout<<"[";
  for (int i=0; i<len; i++) {
    std::cout<<array[i]<<",";
  }
  std::cout<<"]"<<std::endl;
}

#endif