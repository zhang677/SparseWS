#include "../../utils/dataloader.h"
#include <math.h> 

struct index_t {
    int32_t crd[2];

    bool operator < (const index_t& other) const {
        for (int i = 0; i < 2; ++i) {
            if (crd[i] != other.crd[i]) {
                return crd[i] < other.crd[i];
            }
        }
        return false;
    }
};

int compute_map(taco_tensor_t *C, taco_tensor_t *A, taco_tensor_t *B) {
  int C1_dimension = (int)(C->dimensions[0]);
  int B2_dimension = (int)(B->dimensions[1]);
  int* restrict C2_pos = (int*)(C->indices[1][0]);
  int* restrict C2_crd = (int*)(C->indices[1][1]);
  float* restrict C_vals = (float*)(C->vals);
  int* restrict A1_pos = (int*)(A->indices[0][0]);
  int* restrict A1_crd = (int*)(A->indices[0][1]);
  int* restrict A2_pos = (int*)(A->indices[1][0]);
  int* restrict A2_crd = (int*)(A->indices[1][1]);
  float* restrict A_vals = (float*)(A->vals);
  int* restrict B1_pos = (int*)(B->indices[0][0]);
  int* restrict B1_crd = (int*)(B->indices[0][1]);
  int* restrict B2_pos = (int*)(B->indices[1][0]);
  int* restrict B2_crd = (int*)(B->indices[1][1]);
  float* restrict B_vals = (float*)(B->vals);

  C2_pos = (int32_t*)malloc(sizeof(int32_t) * (C1_dimension + 1));
  C2_pos[0] = 0;
  for (int32_t pC2 = 1; pC2 < (C1_dimension + 1); pC2++) {
    C2_pos[pC2] = 0;
  }

  std::map<index_t, float> w_all;
  index_t w_point;

  int32_t jA = A1_pos[0];
  int32_t pA1_end = A1_pos[1];
  int32_t jB = B1_pos[0];
  int32_t pB1_end = B1_pos[1];

  while (jA < pA1_end && jB < pB1_end) {
    int32_t jA0 = A1_crd[jA];
    int32_t jB0 = B1_crd[jB];
    int32_t j = TACO_MIN(jA0, jB0);
    if (jA0 == j && jB0 == j) {
      for (int32_t iA = A2_pos[jA]; iA < A2_pos[jA + 1]; iA++) {
        int32_t i = A2_crd[iA];
        w_point.crd[0] = i;
        for (int32_t kB = B2_pos[jB]; kB < B2_pos[jB + 1]; kB++) {
          int32_t k = B2_crd[kB];
          w_point.crd[1] = k;
          auto handle = w_all.insert(std::make_pair(w_point, A_vals[iA] * B_vals[kB]));
          if (!handle.second) {
            handle.first->second += A_vals[iA] * B_vals[kB];
          }
        }
      }
    }
    jA += (int32_t)(jA0 == j);
    jB += (int32_t)(jB0 == j);
  }

  int w_all_size = w_all.size();

  C2_crd = (int32_t*)malloc(sizeof(int32_t) * w_all_size);
  float* w_vals = (float*)malloc(sizeof(float) * w_all_size);
  int iw_id = 0;
  for (auto iw = w_all.begin(); iw != w_all.end(); ++iw) {
    C2_crd[iw_id] = iw->first.crd[1];
    w_vals[iw_id] = iw->second;
    iw_id++;
  }

  auto iw = w_all.begin();
  auto pw1_end = w_all.end();

  while (iw != pw1_end) {
    int32_t i = iw->first.crd[0];
    auto curr = iw;
    auto w1_segend = iw++;
    while (w1_segend != pw1_end && w1_segend->first.crd[0] == i) {
      w1_segend++;
    }
    C2_pos[i + 1] = std::distance(curr, w1_segend);
    iw = w1_segend;
  }

  int32_t csC2 = 0;
  for (int32_t pC2 = 1; pC2 < (C1_dimension + 1); pC2++) {
    csC2 += C2_pos[pC2];
    C2_pos[pC2] = csC2;
  }

  C->indices[1][0] = (int32_t*)(C2_pos);
  C->indices[1][1] = (int32_t*)(C2_crd);
  C->vals = (float*)w_vals;
  return 0;
}

double DCSC_DCSR_noacc_map(taco_tensor_t *A, taco_tensor_t *B, taco_tensor_t* C, int32_t warmup, int32_t repeat, bool bench = false, bool print = false) {
  // std::cout << "Capacity: " << w_cap << std::endl;
  for (int i = 0; i < warmup; i++) {
    compute_map(C,A,B);
    if (bench) {
      free(C->vals);
      free(C->indices[1][0]);
      free(C->indices[1][1]);
    }
  }
  double start = clock();
  for (int i = 0; i < repeat; i++) {
    compute_map(C,A,B);
    if (bench && i != repeat - 1) {
      free(C->vals);
      free(C->indices[1][0]);
      free(C->indices[1][1]);
    }
  }
  double end = clock();
  double duration = (double)(end - start) / (CLOCKS_PER_SEC * repeat);
  if (print) {
    print_taco_tensor_DC(A);
    print_taco_tensor_DC(B);
    print_taco_tensor_DC(C);
  }
  return duration;
}