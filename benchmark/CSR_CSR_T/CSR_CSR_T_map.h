#include "../utils/dataloader.h"
#include <map>

auto comp = [](const pair<int,int>& a, const pair<int,int>& b) {
  if (a.first < b.first) {
    return true;
  } else if (a.first == b.first) {
    return a.second < b.second;
  } else {
    return false;
  }
};
typedef std::map<pair<int,int>,float,decltype(comp)> MapType;
typedef std::map<pair<int,int>,float,decltype(comp)>::const_iterator MapIter;
MapType w_acc(comp);

int w_cmp(const pair<int,int>& a, const pair<int,int>& b) {
  if (a.first < b.first) {
    return -1;
  } else if (a.first == b.first) {
    if (a.second < b.second) {
      return -1;
    } else if (a.second == b.second) {
      return 0;
    } else {
      return 1;
    }
  } else {
    return 1;
  }
}

int32_t TryInsert_map(bool* insertFail, int32_t w_cap, int32_t* crds, float val) {
  if (w_acc.size() == w_cap) {
    *insertFail = true;
    return w_acc.size();
  } else {
    const auto [it, success] = w_acc.insert({{crds[0], crds[1]}, val});
    if (!success) {
      w_acc[{crds[0], crds[1]}] += val;
    }
    *insertFail = false;
    return w_acc.size(); 
  }
}

int32_t Merge_map(int32_t* COO1_crd, int32_t* COO2_crd, float* COO_vals, int32_t COO_size) {
  if (COO_size == 0) {
    int i = 0;
    for (MapIter it = w_acc.begin(); it != w_acc.end(); ++it) {
      COO1_crd[i] = it->first.first;
      COO2_crd[i] = it->first.second;
      COO_vals[i] = it->second;
      i++;
    }
    return w_acc.size();
  }
  int32_t* tmp_COO_crd[2];
  float* tmp_COO_vals;
  tmp_COO_crd[0] = (int32_t*)malloc(sizeof(int32_t) * (w_acc.size() + COO_size));
  tmp_COO_crd[1] = (int32_t*)malloc(sizeof(int32_t) * (w_acc.size()+ COO_size));
  tmp_COO_vals = (float*)malloc(sizeof(float) * (w_acc.size() + COO_size));
  MapIter accumulator_pointer = w_acc.begin();
  int content_pointer = 0;
  int target_pointer = 0;
  while(accumulator_pointer != w_acc.end() && content_pointer < COO_size) {
    if (w_cmp(accumulator_pointer->first, {COO1_crd[content_pointer], COO2_crd[content_pointer]}) == 0) {
      tmp_COO_crd[0][target_pointer] = accumulator_pointer->first.first;
      tmp_COO_crd[1][target_pointer] = accumulator_pointer->first.second;
      tmp_COO_vals[target_pointer] = accumulator_pointer->second + COO_vals[content_pointer];
      ++accumulator_pointer;
      content_pointer ++;
      target_pointer ++;
    } else if (w_cmp(accumulator_pointer->first, {COO1_crd[content_pointer], COO2_crd[content_pointer]})  < 0) {
      tmp_COO_crd[0][target_pointer] = accumulator_pointer->first.first;
      tmp_COO_crd[1][target_pointer] = accumulator_pointer->first.second;
      tmp_COO_vals[target_pointer] = accumulator_pointer->second;
      ++accumulator_pointer;
      target_pointer ++;
    } else {
      tmp_COO_crd[0][target_pointer] = COO1_crd[content_pointer];
      tmp_COO_crd[1][target_pointer] = COO2_crd[content_pointer];
      tmp_COO_vals[target_pointer] = COO_vals[content_pointer];
      content_pointer ++;
      target_pointer ++;
    }
  }
  while(accumulator_pointer != w_acc.end()) {
    tmp_COO_crd[0][target_pointer] = accumulator_pointer->first.first;
    tmp_COO_crd[1][target_pointer] = accumulator_pointer->first.second;
    tmp_COO_vals[target_pointer] = accumulator_pointer->second;
    ++accumulator_pointer;
    target_pointer ++;
  }
  while(content_pointer < COO_size) {
    tmp_COO_crd[0][target_pointer] = COO1_crd[content_pointer];
    tmp_COO_crd[1][target_pointer] = COO2_crd[content_pointer];
    tmp_COO_vals[target_pointer] = COO_vals[content_pointer];
    content_pointer ++;
    target_pointer ++;
  }
  for (int i = 0; i < target_pointer; i++) {
    COO1_crd[i] = tmp_COO_crd[0][i];
    COO2_crd[i] = tmp_COO_crd[1][i];
    COO_vals[i] = tmp_COO_vals[i];
  }
  free(tmp_COO_crd[0]);
  free(tmp_COO_crd[1]);
  free(tmp_COO_vals);
  return target_pointer;
}


int compute(taco_tensor_t *C, taco_tensor_t *A, taco_tensor_t *B, int32_t w_cap) {
  int C1_dimension = (int)(C->dimensions[0]);
  int* restrict C2_pos = (int*)(C->indices[1][0]);
  int* restrict C2_crd = (int*)(C->indices[1][1]);
  float* restrict C_vals = (float*)(C->vals);
  int A1_dimension = (int)(A->dimensions[0]);
  int* restrict A2_pos = (int*)(A->indices[1][0]);
  int* restrict A2_crd = (int*)(A->indices[1][1]);
  float* restrict A_vals = (float*)(A->vals);
  int B2_dimension = (int)(B->dimensions[1]);
  int* restrict B2_pos = (int*)(B->indices[1][0]);
  int* restrict B2_crd = (int*)(B->indices[1][1]);
  float* restrict B_vals = (float*)(B->vals);

  C2_pos = (int32_t*)malloc(sizeof(int32_t) * (C1_dimension + 1));
  C2_pos[0] = 0;
  for (int32_t pC2 = 1; pC2 < (C1_dimension + 1); pC2++) {
    C2_pos[pC2] = 0;
  }
  int32_t C2_crd_size = 1048576;
  C2_crd = (int32_t*)malloc(sizeof(int32_t) * C2_crd_size);
  int32_t iC = 0;
  int32_t C_capacity = 1048576;
  C_vals = (float*)malloc(sizeof(float) * C_capacity);

  int32_t w_accumulator_size = 0;
  int32_t w_all_capacity = w_cap;
  int32_t w_all_size = 0;
  int32_t* restrict w1_pos = 0;
  w1_pos = (int32_t*)malloc(sizeof(int32_t) * 2);
  int32_t* restrict w1_crd = 0;
  int32_t* restrict w2_crd = 0;
  w1_crd = (int32_t*)malloc(sizeof(int32_t) * w_all_capacity);
  w2_crd = (int32_t*)malloc(sizeof(int32_t) * w_all_capacity);
  float* restrict w_vals = 0;
  w_vals = (float*)malloc(sizeof(float) * w_all_capacity);
  bool* restrict w_insertFail = 0;
  w_insertFail = (bool*)malloc(sizeof(bool) * 1);
  w_insertFail[0] = 0;
  int32_t* restrict w_point = 0;
  w_point = (int32_t*)malloc(sizeof(int32_t) * 2);
  for (int32_t i = 0; i < A1_dimension; i++) {
    for (int32_t jA = A2_pos[i]; jA < A2_pos[i+1]; jA++) {
      int32_t j = A2_crd[jA];
      for (int32_t kB = B2_pos[j]; kB < B2_pos[j+1]; kB++) {
        int32_t k = B2_crd[kB];
        w_point[0] = k;
        w_point[1] = i; // Transpose
        // Try to insert to the accumulator array
        w_accumulator_size = TryInsert_map(w_insertFail, w_cap, w_point, A_vals[jA] * B_vals[kB]);
        if (w_insertFail[0]) {
          //std::cout << "Seperate point: (" << w_point[0]<< "," << w_point[1] << ")" << std::endl;
          // Enlarge
          if(w_accumulator_size + w_all_size > w_all_capacity) {
            w_all_capacity = w_accumulator_size + w_all_size + 1;
            // std::cout << w_all_capacity << std::endl;
            w1_crd = (int32_t*)realloc(w1_crd, sizeof(int32_t) * w_all_capacity);
            w2_crd = (int32_t*)realloc(w2_crd, sizeof(int32_t) * w_all_capacity);
            w_vals = (float*)realloc(w_vals, sizeof(float) * w_all_capacity);
          }
          // Merge
          w_all_size = Merge_map(w1_crd, w2_crd, w_vals, w_all_size);
          w_cap = w_cap * 2;
          // Clear
          w_acc.clear();
          // Insert the wspace that conflicts
          w_accumulator_size = TryInsert_map(w_insertFail, w_cap, w_point, A_vals[jA] * B_vals[kB]);
        }
      }
    }
  }
  if (w_accumulator_size > 0) {
    // Sort
    // The map is alreadt sorted by key.
    // Enlarge
    if(w_accumulator_size + w_all_size > w_all_capacity) {
        w_all_capacity = w_all_capacity * 2;
        //std::cout << w_all_capacity << std::endl;
        w1_crd = (int32_t*)realloc(w1_crd, sizeof(int32_t) * w_all_capacity);
        //std::cout << w_all_capacity << std::endl;
        w2_crd = (int32_t*)realloc(w2_crd, sizeof(int32_t) * w_all_capacity);
        w_vals = (float*)realloc(w_vals, sizeof(float) * w_all_capacity);
    }
    // Merge
    w_all_size = Merge_map(w1_crd, w2_crd, w_vals, w_all_size);
    // Clear
    w_acc.clear();
  }

  /*
  int w1_pos[2] = {0,w_all_size};
  pack_C(C, w1_pos, w1_crd, w2_crd, w_vals);
  */
  //std::cout << "Start packing" << std::endl;
  w1_pos[0] = 0;
  w1_pos[1] = w_all_size;
  int32_t kw = w1_pos[0];
  int32_t pw1_end = w1_pos[1];


  while (kw < pw1_end) {
    int32_t k = w1_crd[kw];
    int32_t w1_segend = kw + 1;
    while (w1_segend < pw1_end && w1_crd[w1_segend] == k) {
      w1_segend++;
    }
    int32_t kC = k;
    int32_t pC2_begin = iC;

    for (int32_t iw = kw; iw < w1_segend; iw++) {
      int32_t i = w2_crd[iw];
      if (C_capacity <= iC) {
        C_vals = (float*)realloc(C_vals, sizeof(float) * (C_capacity * 2));
        C_capacity *= 2;
      }
      C_vals[iC] = w_vals[iw];
      if (C2_crd_size <= iC) {
        C2_crd = (int32_t*)realloc(C2_crd, sizeof(int32_t) * (C2_crd_size * 2));
        C2_crd_size *= 2;
      }
      C2_crd[iC] = i;
      iC++;
    }

    C2_pos[k + 1] = iC - pC2_begin;
    kw = w1_segend;
  }

  free(w_vals);
  free(w_insertFail);
  free(w_point);
  free(w1_crd);
  free(w2_crd);
  free(w1_pos);

  int32_t csC2 = 0;
  for (int32_t pC2 = 1; pC2 < (C1_dimension + 1); pC2++) {
    csC2 += C2_pos[pC2];
    C2_pos[pC2] = csC2;
  }

  C->indices[1][0] = (int32_t*)(C2_pos);
  C->indices[1][1] = (int32_t*)(C2_crd);
  C->vals = (float*)C_vals;
  return 0;
}

void CSR_CSR_T_map(const string& A_name, const string& B_name, taco_tensor_t* C, int32_t w_cap, bool print = false) {
  // C(k,i) = A(i,j) * B(j,k); C: CSR, A: CSR, B: CSR
  vector<int> indptr;
  vector<int> indices;
  vector<int> id_buffer;
  vector<float> value;
  int nrow;
  int ncol;
  int nnz;
  read_mtx_csr(A_name.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
  taco_tensor_t A = DC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
  read_mtx_csr(B_name.data(), nrow, ncol, nnz, indptr, indices, id_buffer, value);
  taco_tensor_t B = DC_to_taco_tensor(indptr,indices,value,nrow,ncol,nnz,{0,1});
  init_taco_tensor_DC(C, nrow, ncol, {0,1});
  //w_cap = pow(2,int(log2(nnz / 2))); // heuristic
  compute(C,&A,&B,w_cap);
  if (print) {
    print_taco_tensor_DC(&A);
    print_taco_tensor_DC(&B);
    print_taco_tensor_DC(C);
  }
}