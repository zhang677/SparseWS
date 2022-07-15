#ifndef DATA_LOADER
#define DATA_LOADER

#include "lib.h"
#include "mmio.h"
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
using namespace std;
typedef int Index;
typedef float DType;
// Load sparse matrix from an mtx file. Only non-zero positions are loaded,
// and values are dropped.
void read_mtx_csr(const char *filename, int &nrow, int &ncol, int &nnz,
                   std::vector<int> &csr_indptr_buffer,
                   std::vector<int> &csr_indices_buffer,
                   std::vector<int> &coo_rowind_buffer,
                   std::vector<float> &csr_value_buffer) {
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
  // printf("Reading matrix %d rows, %d columns, %d nnz.\n", nrow, ncol, nnz);

  /// read tuples

  std::vector<std::tuple<int, int>> coords;
  int row_id, col_id;
  float dummy;
  for (int64_t i = 0; i < nnz; i++) {
    if (fscanf(f, "%d", &row_id) == EOF) {
      std::cout << "Error: not enough rows in mtx file.\n";
      exit(EXIT_FAILURE);
    } else {
      fscanf(f, "%d", &col_id);
      if (mm_is_integer(matcode) || mm_is_real(matcode)) {
        fscanf(f, "%f", &dummy);
      } else if (mm_is_complex(matcode)) {
        fscanf(f, "%f", &dummy);
        fscanf(f, "%f", &dummy);
      }
      // mtx format is 1-based
      coords.push_back(std::make_tuple(row_id - 1, col_id - 1));
    }
  }

  /// make symmetric

  if (mm_is_symmetric(matcode)) {
    std::vector<std::tuple<int, int>> new_coords;
    for (auto iter = coords.begin(); iter != coords.end(); iter++) {
      int i = std::get<0>(*iter);
      int j = std::get<1>(*iter);
      if (i != j) {
        new_coords.push_back(std::make_tuple(i, j));
        new_coords.push_back(std::make_tuple(j, i));
      } else
        new_coords.push_back(std::make_tuple(i, j));
    }
    std::sort(new_coords.begin(), new_coords.end());
    coords.clear();
    for (auto iter = new_coords.begin(); iter != new_coords.end(); iter++) {
      if ((iter + 1) == new_coords.end() || (*iter != *(iter + 1))) {
        coords.push_back(*iter);
      }
    }
  } else {
    std::sort(coords.begin(), coords.end());
  }
  /// generate csr from coo

  csr_indptr_buffer.clear();
  csr_indices_buffer.clear();
  csr_value_buffer.clear();

  int curr_pos = 0;
  csr_indptr_buffer.push_back(0);
  for (int64_t row = 0; row < nrow; row++) {
    while ((curr_pos < nnz) && (std::get<0>(coords[curr_pos]) == row)) {
      csr_indices_buffer.push_back(std::get<1>(coords[curr_pos]));
      coo_rowind_buffer.push_back(std::get<0>(coords[curr_pos]));
      csr_value_buffer.push_back((float)(std::rand() % 10) / 10);
      curr_pos++;
    }
    // assert((std::get<0>(coords[curr_pos]) > row || curr_pos == nnz));
    csr_indptr_buffer.push_back(curr_pos);
  }
  nnz = csr_indices_buffer.size();
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

void read_mtx_csc(const char *filename, int &nrow, int &ncol, int &nnz,
                   std::vector<int> &csc_indptr_buffer,
                   std::vector<int> &csc_indices_buffer,
                   std::vector<int> &coo_colind_buffer,
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
  // printf("Reading matrix %d rows, %d columns, %d nnz.\n", nrow, ncol, nnz);

  /// read tuples

  std::vector<std::tuple<int, int>> coords;
  int row_id, col_id;
  float dummy;
  for (int64_t i = 0; i < nnz; i++) {
    if (fscanf(f, "%d", &row_id) == EOF) {
      std::cout << "Error: not enough rows in mtx file.\n";
      exit(EXIT_FAILURE);
    } else {
      fscanf(f, "%d", &col_id);
      if (mm_is_integer(matcode) || mm_is_real(matcode)) {
        fscanf(f, "%f", &dummy);
      } else if (mm_is_complex(matcode)) {
        fscanf(f, "%f", &dummy);
        fscanf(f, "%f", &dummy);
      }
      // mtx format is 1-based
      coords.push_back(std::make_tuple(col_id - 1, row_id - 1));
    }
  }

  /// make symmetric

  if (mm_is_symmetric(matcode)) {
    std::vector<std::tuple<int, int>> new_coords;
    for (auto iter = coords.begin(); iter != coords.end(); iter++) {
      int i = std::get<0>(*iter);
      int j = std::get<1>(*iter);
      if (i != j) {
        new_coords.push_back(std::make_tuple(i, j));
        new_coords.push_back(std::make_tuple(j, i));
      } else
        new_coords.push_back(std::make_tuple(i, j));
    }
    std::sort(new_coords.begin(), new_coords.end());
    coords.clear();
    for (auto iter = new_coords.begin(); iter != new_coords.end(); iter++) {
      if ((iter + 1) == new_coords.end() || (*iter != *(iter + 1))) {
        coords.push_back(*iter);
      }
    }
  } else {
    std::sort(coords.begin(), coords.end());
  }
  /// generate csr from coo

  csc_indptr_buffer.clear();
  csc_indices_buffer.clear();
  csc_value_buffer.clear();

  int curr_pos = 0;
  csc_indptr_buffer.push_back(0);
  for (int64_t col = 0; col < ncol; col++) {
    while ((curr_pos < nnz) && (std::get<0>(coords[curr_pos]) == col)) {
      csc_indices_buffer.push_back(std::get<1>(coords[curr_pos]));
      coo_colind_buffer.push_back(std::get<0>(coords[curr_pos]));
      csc_value_buffer.push_back((float)(std::rand() % 10) / 10);
      curr_pos++;
    }
    // assert((std::get<0>(coords[curr_pos]) > row || curr_pos == nnz));
    csc_indptr_buffer.push_back(curr_pos);
  }
  nnz = csc_indices_buffer.size();
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
    cout<<"[";
    for (int i=0; i<=t->indices[0][0][0]; i++) {
      cout<<t->indices[1][0][i]<<",";
    }
    cout<<"]"<<endl;
    cout<<"[";
    for (int i=0; i<t->indices[1][0][t->indices[0][0][0]]; i++) {
      cout<<t->indices[1][1][i]<<",";
    }
    cout<<"]"<<endl;
    cout<<"[";
    for (int i=0; i<t->indices[1][0][t->indices[0][0][0]]; i++) {
      cout<<t->vals[i]<<",";
    }
    cout<<"]"<<endl;    
}

void dense_to_compressed(vector<int>& hindptr, vector<int>& hindices, vector<int>& new_indptr, const vector<int>& indptr) {
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

taco_tensor_t CC_to_taco_tensor(vector<int>& old_indptr_buffer, vector<int>& indices_buffer, vector<float>& value_buffer,
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

}

#endif