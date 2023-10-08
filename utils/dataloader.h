#ifndef DATA_LOADER
#define DATA_LOADER

#include "lib.h"
#include "mmio.h"
#include <time.h>
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <tuple>
#include <typeinfo>
#include <vector>
#include <set>
#include <chrono>
#include "splatt.h"
// #include "sptensor.h"
// #include "csf.h"


using namespace std;
typedef int Index;
typedef float DType;

std::vector<std::string> split(const std::string &str, const std::string &delim, bool keepDelim) {
  std::vector<std::string> results;
  size_t prev = 0;
  size_t next = 0;

  while ((next = str.find(delim, prev)) != std::string::npos) {
    if (next - prev != 0) {
      std::string substr = ((keepDelim) ? delim : "")
                         + str.substr(prev, next-prev);
      results.push_back(substr);
    }
    prev = next + delim.size();
  }

  if (prev < str.size()) {
    string substr = ((keepDelim) ? delim : "") + str.substr(prev);
    results.push_back(substr);
  }

  return results;
}


// Load sparse matrix from an mtx file. Only non-zero positions are loaded,
// and values are dropped.
void read_mtx_csr(const char *filename, int &nrow, int &ncol, int &nnz,
                   std::vector<int> &csr_indptr_buffer,
                   std::vector<int> &csr_indices_buffer,
                   std::vector<int> &coo_rowind_buffer,
                   std::vector<float> &csr_value_buffer,
                   bool one_based = true) {
  FILE *f;

  if ((f = fopen(filename, "r")) == NULL) {
    printf("File %s not found", filename);
    exit(EXIT_FAILURE);
  }

  MM_typecode matcode;
  // Read MTX banner
  if (mm_read_banner(f, &matcode) != 0) {
    printf("Could not process this file.\n");
    exit(EXIT_FAILURE);
  }
  if (mm_read_mtx_crd_size(f, &nrow, &ncol, &nnz) != 0) {
    printf("Could not process this file.\n");
    exit(EXIT_FAILURE);
  }
  if (!mm_is_matrix(matcode) || !mm_is_sparse(matcode) ||
      !mm_is_coordinate(matcode)) {
    printf("Could not process this file.\n");
    exit(EXIT_FAILURE);
  }
  // printf("Reading matrix %d rows, %d columns, %d nnz.\n", nrow, ncol, nnz);

  /// read tuples
  std::vector<float> value_temp;
  std::vector<std::tuple<int, int>> coords;
  int row_id, col_id;
  float dummy;
  for (int64_t i = 0; i < nnz; i++) {
    if (fscanf(f, "%d", &row_id) == EOF) {
      std::cout << "Error: not enough rows in mtx file.\n";
      exit(EXIT_FAILURE);
    } else {
      fscanf(f, "%d", &col_id);
      fscanf(f, "%f", &dummy);
      // mtx format is 1-based
      if (one_based) {
        coords.push_back(std::make_tuple(row_id - 1, col_id - 1));
      } else {
        coords.push_back(std::make_tuple(row_id, col_id));
      }
      value_temp.push_back(dummy);
    }
  }

  /// make symmetric
  std::vector<int> index;
  std::vector<float> new_value;
  std::vector<std::tuple<int, int>> new_coords;
  if (mm_is_symmetric(matcode)) {
    int cur_nz = 0;
    int cur_ptr = 0;
    for (auto iter = coords.begin(); iter != coords.end(); iter++) {
      int i = std::get<0>(*iter);
      int j = std::get<1>(*iter);
      if (i != j) {
        new_coords.push_back(std::make_tuple(i, j));
        index.push_back(cur_nz);
        new_value.push_back(value_temp[cur_ptr]);
        cur_nz++;
        new_coords.push_back(std::make_tuple(j, i));
        index.push_back(cur_nz);
        new_value.push_back(value_temp[cur_ptr]);
        cur_nz++;
      } else {
        new_coords.push_back(std::make_tuple(i, j));
        index.push_back(cur_nz);
        new_value.push_back(value_temp[cur_ptr]);
        cur_nz++;
      }
      cur_ptr++;
    }
    //std::sort(new_coords.begin(), new_coords.end());
    std::sort(index.begin(), index.end(),
              [&new_coords](int i1, int i2) {
                return std::get<0>(new_coords[i1]) == std::get<0>(new_coords[i2]) ? std::get<1>(new_coords[i1]) < std::get<1>(new_coords[i2]) : std::get<0>(new_coords[i1]) < std::get<0>(new_coords[i2]);
              });
    nnz = cur_nz;
  } else {
    boost::range::push_back(index, boost::irange(0, nnz));
    std::sort(index.begin(), index.end(),
          [&coords](int i1, int i2) {
            return std::get<0>(coords[i1]) == std::get<0>(coords[i2]) ? std::get<1>(coords[i1]) < std::get<1>(coords[i2]) : std::get<0>(coords[i1]) < std::get<0>(coords[i2]);
          });
  }
  /// generate csr from coo

  csr_indptr_buffer.clear();
  csr_indices_buffer.clear();
  csr_value_buffer.clear();
  coo_rowind_buffer.clear();

  int curr_pos = 0;
  csr_indptr_buffer.push_back(0);
  if (mm_is_symmetric(matcode)) {
    for (int64_t row = 0; row < nrow; row++) {
      while ((curr_pos < nnz) && (std::get<0>(new_coords[index[curr_pos]]) == row)) {
        csr_indices_buffer.push_back(std::get<1>(new_coords[index[curr_pos]]));
        coo_rowind_buffer.push_back(std::get<0>(new_coords[index[curr_pos]]));
        //csr_value_buffer.push_back(new_value[index[curr_pos]]);
        csr_value_buffer.push_back((float)rand()/RAND_MAX);
        curr_pos++;
      }
      // assert((std::get<0>(coords[curr_pos]) > row || curr_pos == nnz));
      csr_indptr_buffer.push_back(curr_pos);
    }
  } else {
    for (int64_t row = 0; row < nrow; row++) {
      while ((curr_pos < nnz) && (std::get<0>(coords[index[curr_pos]]) == row)) {
        csr_indices_buffer.push_back(std::get<1>(coords[index[curr_pos]]));
        coo_rowind_buffer.push_back(std::get<0>(coords[index[curr_pos]]));
        //csr_value_buffer.push_back(value_temp[index[curr_pos]]);
        csr_value_buffer.push_back((float)rand()/RAND_MAX);
        curr_pos++;
      }
      // assert((std::get<0>(coords[curr_pos]) > row || curr_pos == nnz));
      csr_indptr_buffer.push_back(curr_pos);
    }
  }
  nnz = csr_indices_buffer.size();
  fclose(f);
}

template <typename Index>
void compressedRow(Index nrow, std::vector<Index> &rowind,
                   std::vector<Index> &rowptr) {
  int curr_pos = 0;
  rowptr.push_back(0);
  int nnz = rowind.size();
  for (int64_t row = 0; row < nrow; row++) {
    while ((curr_pos < nnz) && ((rowind[curr_pos]) == row)) {
      curr_pos++;
    }
    rowptr.push_back(curr_pos);
  }
}

template <typename Index>
void transpose(Index ncol, std::vector<Index> &row, std::vector<Index> &col,
               std::vector<Index> &row_t, std::vector<Index> &col_t) {
  int nnz = col.size();
  Index *hist = new Index[ncol];
  Index *col_tmp = new Index[ncol + 1];
  memset(hist, 0x0, sizeof(Index) * ncol);
  for (int t = 0; t < nnz; t++)
    hist[col[t]]++;
  col_tmp[0] = 1;
  for (int c = 1; c <= ncol; ++c)
    col_tmp[c] = col_tmp[c - 1] + hist[c - 1];
  for (int nid = 0; nid < nnz; nid++) {
    int col_ = col[nid];
    int q = col_tmp[col_];
    row_t[q - 1] = col[nid];
    col_t[q - 1] = row[nid];
    col_tmp[col_]++;
  }
  delete[] hist;
  delete[] col_tmp;
}

// void read_mtx_csc(const char *filename, int &nrow, int &ncol, int &nnz,
//                    std::vector<int> &csc_indptr_buffer,
//                    std::vector<int> &csc_indices_buffer,
//                    std::vector<int> &coo_colind_buffer,
//                    std::vector<float> &csc_value_buffer) {
//   read_mtx_csr(filename, nrow, ncol, nnz, csc_indptr_buffer, csc_indices_buffer, coo_colind_buffer, csc_value_buffer);

//   EigenCSC temp = to_eigen_csc(nrow, ncol, nnz, coo_colind_buffer, csc_indices_buffer, csc_value_buffer);
//   csc_indptr_buffer.clear();
//   csc_indices_buffer.clear();
//   csc_value_buffer.clear();
//   for (int k = 0; k <= ncol; k++) {
//     csc_indptr_buffer.push_back(temp.outerIndexPtr()[k]);
//   }
//   for (int k = 0; k < nnz; k++) {
//     csc_indices_buffer.push_back(temp.innerIndexPtr()[k]);
//     csc_value_buffer.push_back(temp.valuePtr()[k]);
//   }

// }

// Load sparse matrix from an mtx file. Only non-zero positions are loaded,
// and values are dropped.
void read_mtx_csc(const char *filename, int &nrow, int &ncol, int &nnz,
                   std::vector<int> &csc_indptr_buffer,
                   std::vector<int> &csc_indices_buffer,
                   std::vector<int> &coo_rowind_buffer,
                   std::vector<float> &csc_value_buffer) {
  FILE *f;

  if ((f = fopen(filename, "r")) == NULL) {
    printf("File %s not found", filename);
    exit(EXIT_FAILURE);
  }

  MM_typecode matcode;
  // Read MTX banner
  if (mm_read_banner(f, &matcode) != 0) {
    printf("Could not process this file.\n");
    exit(EXIT_FAILURE);
  }
  if (mm_read_mtx_crd_size(f, &nrow, &ncol, &nnz) != 0) {
    printf("Could not process this file.\n");
    exit(EXIT_FAILURE);
  }
  if (!mm_is_matrix(matcode) || !mm_is_sparse(matcode) ||
      !mm_is_coordinate(matcode)) {
    printf("Could not process this file.\n");
    exit(EXIT_FAILURE);
  }
  // printf("Reading matrix %d rows, %d columns, %d nnz.\n", nrow, ncol, nnz);

  /// read tuples
  std::vector<float> value_temp;
  std::vector<std::tuple<int, int>> coords;
  int row_id, col_id;
  float dummy;
  for (int64_t i = 0; i < nnz; i++) {
    if (fscanf(f, "%d", &row_id) == EOF) {
      std::cout << "Error: not enough rows in mtx file.\n";
      exit(EXIT_FAILURE);
    } else {
      fscanf(f, "%d", &col_id);
      fscanf(f, "%f", &dummy);
      // mtx format is 1-based
      coords.push_back(std::make_tuple(row_id - 1, col_id - 1));
      value_temp.push_back(dummy);
    }
  }

  /// make symmetric
  std::vector<int> index;
  std::vector<float> new_value;
  std::vector<std::tuple<int, int>> new_coords;
  if (mm_is_symmetric(matcode)) {
    int cur_nz = 0;
    int cur_ptr = 0;
    for (auto iter = coords.begin(); iter != coords.end(); iter++) {
      int i = std::get<0>(*iter);
      int j = std::get<1>(*iter);
      if (i != j) {
        new_coords.push_back(std::make_tuple(i, j));
        index.push_back(cur_nz);
        new_value.push_back(value_temp[cur_ptr]);
        cur_nz++;
        new_coords.push_back(std::make_tuple(j, i));
        index.push_back(cur_nz);
        new_value.push_back(value_temp[cur_ptr]);
        cur_nz++;
      } else {
        new_coords.push_back(std::make_tuple(i, j));
        index.push_back(cur_nz);
        new_value.push_back(value_temp[cur_ptr]);
        cur_nz++;
      }
      cur_ptr++;
    }
    //std::sort(new_coords.begin(), new_coords.end());
    std::sort(index.begin(), index.end(),
              [&new_coords](int i1, int i2) {
                return std::get<0>(new_coords[i1]) == std::get<0>(new_coords[i2]) ? std::get<1>(new_coords[i1]) < std::get<1>(new_coords[i2]) : std::get<0>(new_coords[i1]) < std::get<0>(new_coords[i2]);
              });
    nnz = cur_nz;
  } else {
    boost::range::push_back(index, boost::irange(0, nnz));
    std::sort(index.begin(), index.end(),
          [&coords](int i1, int i2) {
            return std::get<0>(coords[i1]) == std::get<0>(coords[i2]) ? std::get<1>(coords[i1]) < std::get<1>(coords[i2]) : std::get<0>(coords[i1]) < std::get<0>(coords[i2]);
          });
  }
  /// generate csr from coo

  csc_indptr_buffer.clear();
  csc_indices_buffer.clear();
  csc_value_buffer.clear();
  coo_rowind_buffer.clear();

  int curr_pos = 0;
  csc_indptr_buffer.push_back(0);
  if (mm_is_symmetric(matcode)) {
    for (int64_t col = 0; col < ncol; col++) {
      while ((curr_pos < nnz) && (std::get<0>(new_coords[index[curr_pos]]) == col)) {
        csc_indices_buffer.push_back(std::get<1>(new_coords[index[curr_pos]]));
        coo_rowind_buffer.push_back(std::get<0>(new_coords[index[curr_pos]]));
        //csc_value_buffer.push_back(new_value[index[curr_pos]]);
        csc_value_buffer.push_back((float)rand()/RAND_MAX);
        curr_pos++;
      }
      // assert((std::get<0>(coords[curr_pos]) > row || curr_pos == nnz));
      csc_indptr_buffer.push_back(curr_pos);
    }
  } else {
    for (int64_t col = 0; col < ncol; col++) {
      while ((curr_pos < nnz) && (std::get<0>(coords[index[curr_pos]]) == col)) {
        csc_indices_buffer.push_back(std::get<1>(coords[index[curr_pos]]));
        coo_rowind_buffer.push_back(std::get<0>(coords[index[curr_pos]]));
        //csc_value_buffer.push_back(value_temp[index[curr_pos]]);
        csc_value_buffer.push_back((float)rand()/RAND_MAX);
        curr_pos++;
      }
      // assert((std::get<0>(coords[curr_pos]) > row || curr_pos == nnz));
      csc_indptr_buffer.push_back(curr_pos);
    }
  }
  nnz = csc_indices_buffer.size();
  fclose(f);
}


taco_tensor_t DC_to_taco_tensor(std::vector<int> &indptr_buffer,std::vector<int> &indices_buffer,
                                std::vector<float> &value_buffer,int nrow,int ncol,int nnz, vector<int> mode_ordering) {
    taco_tensor_t t;
    t.order = 2;
    t.mode_ordering = new int32_t[2];
    t.mode_ordering[0] = mode_ordering[0];
    t.mode_ordering[1] = mode_ordering[1];
    t.dimensions = new int32_t[2];
    t.dimensions[0] = nrow;
    t.dimensions[1] = ncol;
    t.indices = new int32_t**[2];
    t.indices[0] = new int32_t*[1];
    t.indices[1] = new int32_t*[2];
    t.indices[0][0] = new int32_t[1];
    t.indices[1][0] = new int32_t[(t.dimensions[mode_ordering[0]]+1)];
    t.indices[1][1] = new int32_t[nnz];
    t.vals = new float[nnz];
    t.vals_size = nnz;

    t.indices[0][0][0] = t.dimensions[mode_ordering[0]];
    std::copy(indptr_buffer.begin(), indptr_buffer.end(), t.indices[1][0]);
    std::copy(indices_buffer.begin(), indices_buffer.end(), t.indices[1][1]);
    std::copy(value_buffer.begin(), value_buffer.end(), t.vals);
    return t;

}

taco_tensor_t DC_to_taco_tensor(int* indptr_buffer, int* indices_buffer,
                                float* value_buffer,int nrow,int ncol,int nnz, vector<int> mode_ordering) {
    taco_tensor_t t;
    t.order = 2;
    t.mode_ordering = new int32_t[2];
    t.mode_ordering[0] = mode_ordering[0];
    t.mode_ordering[1] = mode_ordering[1];
    t.dimensions = new int32_t[2];
    t.dimensions[0] = nrow;
    t.dimensions[1] = ncol;
    t.indices = new int32_t**[2];
    t.indices[0] = new int32_t*[1];
    t.indices[1] = new int32_t*[2];
    t.indices[0][0] = new int32_t[1];
    t.indices[1][0] = new int32_t[(t.dimensions[mode_ordering[0]]+1)];
    t.indices[1][1] = new int32_t[nnz];
    t.vals = new float[nnz];
    t.vals_size = nnz;

    t.indices[0][0][0] = t.dimensions[mode_ordering[0]];
    memcpy(t.indices[1][0], indptr_buffer, (t.dimensions[mode_ordering[0]]+1)*sizeof(int32_t));
    memcpy(t.indices[1][1], indices_buffer, nnz*sizeof(int32_t));
    memcpy(t.vals, value_buffer, nnz*sizeof(float));
    return t;
}


void init_taco_tensor_DC(taco_tensor_t* t, int nrow, int ncol, vector<int> mode_ordering){
    t->order = 2;
    t->mode_ordering = new int32_t[2];
    t->mode_ordering[0] = mode_ordering[0];
    t->mode_ordering[1] = mode_ordering[1];
    t->dimensions = new int32_t[2];
    t->dimensions[0] = nrow;
    t->dimensions[1] = ncol;
    t->indices = new int32_t**[2];
    t->indices[0] = new int32_t*[1];
    t->indices[1] = new int32_t*[2];
    t->indices[0][0] = new int32_t[1 ];
    t->indices[0][0][0] = t->dimensions[mode_ordering[0]];
    t->vals_size = 0;
}


void print_taco_tensor_DC(taco_tensor_t* t) {
    cout<<"("<<t->dimensions[0]<<","<<t->dimensions[1]<<")"<<endl;
    print_array<int>(t->indices[1][0], t->indices[0][0][0] + 1);
    print_array<int>(t->indices[1][1], t->indices[1][0][t->indices[0][0][0]]);
    print_array<float>(t->vals, t->indices[1][0][t->indices[0][0][0]]); 
}

void print_eigen_csr(EigenCSR& C) {
  cout << "(" << C.outerSize() << "," << C.innerSize() << ")" << endl;
  print_array<int>(C.outerIndexPtr(), C.outerSize() + 1);
  print_array<int>(C.innerIndexPtr(), C.outerIndexPtr()[C.outerSize()]);
  print_array<float>(C.valuePtr(), C.outerIndexPtr()[C.outerSize()]);
}

void print_eigen_csc(EigenCSC& C) {
  cout << "(" << C.outerSize() << "," << C.innerSize() << ")" << endl;
  print_array<int>(C.outerIndexPtr(), C.outerSize() + 1);
  print_array<int>(C.innerIndexPtr(), C.outerIndexPtr()[C.outerSize()]);
  print_array<float>(C.valuePtr(), C.outerIndexPtr()[C.outerSize()]);
}

void print_csr(int nrow, int ncol, int* indptr_buffer, int* indices_buffer, float* value_buffer) {
  cout << "(" << nrow << "," << ncol << ")" << endl;
  int nnz = indptr_buffer[nrow];
  print_array<int>(indptr_buffer, nrow + 1);
  print_array<int>(indices_buffer, nnz);
  print_array<float>(value_buffer, nnz);
}

void dense_to_compressed(vector<int>& hindptr, vector<int>& hindices, vector<int>& new_indptr, int* indptr, int nnz) {
    hindices.clear();
    hindptr.clear();
    new_indptr.clear();
    new_indptr.push_back(0);
    for(int i=0; i<nnz-1; i++){
        if (indptr[i+1]!=indptr[i]) {
            hindices.push_back(i);
            new_indptr.push_back(indptr[i+1]);
        }
    }
    hindptr.push_back(0);
    hindptr.push_back(hindices.size());
}

void dense_to_compressed(vector<int>& hindptr, vector<int>& hindices, vector<int>& new_indptr, vector<int>& indptr) {
    hindices.clear();
    hindptr.clear();
    new_indptr.clear();
    new_indptr.push_back(0);
    for(int i=0; i<indptr.size()-1; i++){
        if (indptr[i+1]!=indptr[i]) {
            hindices.push_back(i);
            new_indptr.push_back(indptr[i+1]);
        }
    }
    hindptr.push_back(0);
    hindptr.push_back(hindices.size());
}

taco_tensor_t CC_to_taco_tensor(int* old_indptr_buffer, int* indices_buffer, float* value_buffer,
int nrow, int ncol, int nnz, vector<int> mode_ordering) {
    vector<int> hindptr;
    vector<int> hindices;
    vector<int> indptr_buffer;
    
    taco_tensor_t t;
    t.order = 2;
    t.mode_ordering = new int32_t[2];
    t.mode_ordering[0] = mode_ordering[0];
    t.mode_ordering[1] = mode_ordering[1];
    t.dimensions = new int32_t[2];
    t.dimensions[0] = nrow;
    t.dimensions[1] = ncol;
    dense_to_compressed(hindptr, hindices, indptr_buffer, old_indptr_buffer, t.dimensions[mode_ordering[0]] + 1);
    t.indices = new int32_t**[2];
    t.indices[0] = new int32_t*[2];
    t.indices[1] = new int32_t*[2];
    t.indices[0][0] = new int32_t[2];
    t.indices[0][1] = new int32_t[hindices.size()];
    t.indices[1][0] = new int32_t[hindices.size() + 1];
    t.indices[1][1] = new int32_t[nnz];
    t.vals = new float[nnz];
    t.vals_size = nnz;

    std::copy(hindptr.begin(), hindptr.end(), t.indices[0][0]);
    std::copy(hindices.begin(), hindices.end(), t.indices[0][1]);
    std::copy(indptr_buffer.begin(), indptr_buffer.end(), t.indices[1][0]);
    //std::copy(indices_buffer.begin(), indices_buffer.end(), t.indices[1][1]);
    memcpy(t.indices[1][1], indices_buffer, nnz*sizeof(int32_t));
    //std::copy(value_buffer.begin(), value_buffer.end(), t.vals);
    memcpy(t.vals, value_buffer, nnz*sizeof(float));
    return t;    
}

taco_tensor_t CC_to_taco_tensor(std::vector<int> old_indptr_buffer, std::vector<int> indices_buffer, std::vector<float> value_buffer,
int nrow, int ncol, int nnz, vector<int> mode_ordering) {
    vector<int> hindptr;
    vector<int> hindices;
    vector<int> indptr_buffer;
    dense_to_compressed(hindptr, hindices, indptr_buffer, old_indptr_buffer);
    taco_tensor_t t;
    t.order = 2;
    t.mode_ordering = new int32_t[2];
    t.mode_ordering[0] = mode_ordering[0];
    t.mode_ordering[1] = mode_ordering[1];
    t.dimensions = new int32_t[2];
    t.dimensions[0] = nrow;
    t.dimensions[1] = ncol;
    t.indices = new int32_t**[2];
    t.indices[0] = new int32_t*[2];
    t.indices[1] = new int32_t*[2];
    t.indices[0][0] = new int32_t[2];
    t.indices[0][1] = new int32_t[hindices.size()];
    t.indices[1][0] = new int32_t[hindices.size() + 1];
    t.indices[1][1] = new int32_t[nnz];
    t.vals = new float[nnz];
    t.vals_size = nnz;

    std::copy(hindptr.begin(), hindptr.end(), t.indices[0][0]);
    std::copy(hindices.begin(), hindices.end(), t.indices[0][1]);
    std::copy(indptr_buffer.begin(), indptr_buffer.end(), t.indices[1][0]);
    std::copy(indices_buffer.begin(), indices_buffer.end(), t.indices[1][1]);
    std::copy(value_buffer.begin(), value_buffer.end(), t.vals);
    return t;    
}

void print_taco_tensor_CC(taco_tensor_t* t) {
    cout<<"("<<t->dimensions[0]<<","<<t->dimensions[1]<<")"<<endl;
    cout<<"["<<t->indices[0][0][0]<<","<<t->indices[0][0][1]<<"]"<<endl;
    cout<<"[";
    for (int i=0; i<t->indices[0][0][1]; i++){
      cout<<t->indices[0][1][i]<<",";
    }
    cout<<"]"<<endl;
    cout<<"[";
    for (int i=0; i<=t->indices[0][0][1]; i++){
      cout<<t->indices[1][0][i]<<",";
    }
    cout<<"]"<<endl;
    cout<<"[";
    for (int i=0; i<t->vals_size; i++) {
      cout<<t->indices[1][1][i]<<",";
    }
    cout<<"]"<<endl;
    cout<<"[";
    for (int i=0; i<t->vals_size; i++) {
      cout<<t->vals[i]<<",";
    }
    cout<<"]"<<endl;           

}

void check_csr_taco_eigen(taco_tensor_t& C, EigenCSR& C_true){
  int row = C_true.outerSize();
  int nnz = C_true.outerIndexPtr()[row];
  cout << "row: " << C.indices[0][0][0] << " , " << row << (C.indices[0][0][0] == row ? "True" : "False") << endl;
  cout << "nnz: " << (C.indices[1][0][C.indices[0][0][0]] == nnz ? "True": "False") << endl;
  compare_array<int>(C.indices[1][0], C_true.outerIndexPtr(), row + 1);
  compare_array<int>(C.indices[1][1], C_true.innerIndexPtr(), nnz);
  compare_array<float>(C.vals, C_true.valuePtr(), nnz);
}

void check_csc_taco_eigen(taco_tensor_t& C, EigenCSC& C_true){
  int row = C_true.outerSize();
  int nnz = C_true.outerIndexPtr()[row];
  cout << "row: " << C.indices[0][0][0] << " , " << row << (C.indices[0][0][0] == row ? "True" : "False") << endl;
  cout << "nnz: " << (C.indices[1][0][C.indices[0][0][0]] == nnz ? "True": "False") << endl;
  compare_array<int>(C.indices[1][0], C_true.outerIndexPtr(), row + 1);
  compare_array<int>(C.indices[1][1], C_true.innerIndexPtr(), nnz);
  compare_array<float>(C.vals, C_true.valuePtr(), nnz);
}

taco_tensor_t read_tns_csf(const std::string tensor_path, taco_tensor_t& csft, vector<int>& dims) {
  // 
  // sptensor_t * tt;
  // tt = tt_read(tensor_path.c_str());
  // splatt_csf tensor;
  //   for(idx_t m=0; m < tt->nmodes; ++m) {
  //   tensor.dim_perm[m] = m;
  // }
  // csf_alloc_mode(tt, CSF_MODE_CUSTOM, 0, &tensor, opts);
  
  double *cpd_opts = splatt_default_opts();
  cpd_opts[SPLATT_OPTION_NTHREADS] = omp_get_num_threads();
  cpd_opts[SPLATT_OPTION_NITER] = 0;
  cpd_opts[SPLATT_OPTION_CSF_ALLOC] = SPLATT_CSF_ONEMODE;
  cpd_opts[SPLATT_OPTION_TILE] = SPLATT_DENSETILE;
  cpd_opts[SPLATT_OPTION_VERBOSITY] = SPLATT_VERBOSITY_NONE;

  splatt_idx_t nmodes;
  splatt_csf *tensor;
  splatt_csf_load(tensor_path.c_str(), &nmodes, &tensor, cpd_opts);
  splatt_free_opts(cpd_opts);

  csft.dimensions = new int32_t[tensor->nmodes];
  // for (int i = 0; i < tensor->nmodes; ++i) {
  //   csft.dimensions[i] = tensor->dims[tensor->dim_perm[i]];
  //   std::cout << tensor->dim_perm[i] << std::endl;
  // }
  for (int i = 0; i < tensor->nmodes; ++i) {
    csft.dimensions[i] = dims[i];
    // std::cout << tensor->dim_perm[i] << std::endl;
  }
  csft.indices = new int32_t**[tensor->nmodes];

  // cout<<"[";
  // for (int k=0; k<tensor->pt->nfibs[0]; k++){
  //   cout<<tensor->pt->fids[0][k]<<",";
  // }
  // cout<<"]"<<endl;
  // cout<< (tensor->pt->fids[0] == NULL) << endl;


  for (int i = 1; i < tensor->nmodes; ++i) {
    csft.indices[i] = new int32_t*[2];
    csft.indices[i][0] = new int32_t[tensor->pt->nfibs[i-1] + 1];
    for (int k = 0; k <= tensor->pt->nfibs[i-1]; ++k) {
      csft.indices[i][0][k] = (int32_t)tensor->pt->fptr[i - 1][k];
    }
    csft.indices[i][1] = new int32_t[tensor->pt->nfibs[i]];
    for (int k = 0; k < tensor->pt->nfibs[i]; ++k) {
      csft.indices[i][1][k] = (int32_t)tensor->pt->fids[i][k];
    }
  }
  csft.vals_size = tensor->pt->nfibs[tensor->nmodes - 1];
  csft.vals = new float[csft.vals_size];
  for (int i = 0; i < csft.vals_size; ++i) {
    csft.vals[i] = (float)tensor->pt->vals[i];
  }
  csft.order = tensor->nmodes;

  std::set<int> coords;
  int row_id;
  FILE *f;
  if ((f = fopen(tensor_path.c_str(), "r")) == NULL) {
    printf("File %s not found", tensor_path.c_str());
    exit(EXIT_FAILURE);
  }
  for (int64_t i = 0; i < csft.vals_size; i++) {
    if (fscanf(f, "%d", &row_id) == EOF) {
      std::cout << "Error: not enough rows in tns file.\n";
      exit(EXIT_FAILURE);
    } else {
      // Dummp the remaining line
      fscanf(f, "%*[^\n]\n");
      coords.insert(row_id - 1);
    }
  }
  fclose(f);
  csft.indices[0] = new int32_t*[2];
  csft.indices[0][0] = new int32_t[2];
  csft.indices[0][1] = new int32_t[coords.size()];
  csft.indices[0][0][0] = 0;
  csft.indices[0][0][1] = coords.size();
  int i = 0;
  for (auto it = coords.begin(); it != coords.end(); ++it) {
    csft.indices[0][1][i++] = (*it);
  }

  return csft;
}

void print_taco_tensor_CSF(taco_tensor_t* t) {
    cout<<"("<<t->dimensions[0]<<","<<t->dimensions[1]<<","<<t->dimensions[2]<<")"<<endl;
    int nnz = 1;
    for (int j = 0; j < t->order; j++) {
      cout<<"[";
      for (int i=0; i<=nnz; i++){
        cout<<t->indices[j][0][i]<<",";
      }
      cout<<"]"<<endl;
      cout<<"[";
      for (int i=0; i<t->indices[j][0][nnz]; i++) {
        cout<<t->indices[j][1][i]<<",";
      }
      cout<<"]"<<endl;
      nnz = t->indices[j][0][nnz];
    }
    cout<<"[";
    for (int i=0; i<t->vals_size; i++) {
      cout<<t->vals[i]<<",";
    }
    cout<<"]"<<endl;           

}

void print_taco_tensor_COO(taco_tensor_t* t) {
  cout<<"[";
  for (int i = 0; i < t->vals_size; ++i){
    cout << t->indices[1][0][i] << ",";
  }
  cout<<"]"<<endl; 
  cout<<"[";
  for (int i = 0; i < t->vals_size; ++i){
    cout << t->indices[1][1][i] << ",";
  }
  cout<<"]"<<endl; 
  cout<<"[";
  for (int i = 0; i < t->vals_size; ++i){
    cout << t->vals[i] << ",";
  }
  cout<<"]"<<endl; 
}

class Timer {
public:
    Timer() : beg_(clock_::now()) {}
    void reset() { beg_ = clock_::now(); }
    double elapsed() const {
        return std::chrono::duration_cast<second_>
            (clock_::now() - beg_).count(); }
private:
    typedef std::chrono::high_resolution_clock clock_;
    typedef std::chrono::duration<double, std::ratio<1> > second_;
    std::chrono::time_point<clock_> beg_;
};

#endif