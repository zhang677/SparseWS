#include "../../utils/dataloader.h"
#include <math.h> 

int compute(taco_tensor_t *C, taco_tensor_t *A, taco_tensor_t *B, int32_t w_accumulator_capacity) {
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

  std::vector< std::tuple< std::tuple<int32_t, int32_t>, float> > w_all;
  w_all.reserve(w_accumulator_capacity);
  int w_point[2];

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
        w_point[0] = i;
        for (int32_t kB = B2_pos[jB]; kB < B2_pos[jB + 1]; kB++) {
          int32_t k = B2_crd[kB];
          w_point[1] = k;
          w_all.push_back(std::make_tuple(std::make_tuple(i, k), A_vals[iA] * B_vals[kB]));
        }
      }
    }
    jA += (int32_t)(jA0 == j);
    jB += (int32_t)(jB0 == j);
  }

  // Lexicographically sort w_all
  std::sort(w_all.begin(), w_all.end(), [](const std::tuple< std::tuple<int32_t, int32_t>, float> &a, const std::tuple< std::tuple< int32_t, int32_t>, float> &b) {
    if (std::get<0>(std::get<0>(a)) != std::get<0>(std::get<0>(b))) {
      return std::get<0>(std::get<0>(a)) < std::get<0>(std::get<0>(b));
    } else {
      return std::get<1>(std::get<0>(a)) < std::get<1>(std::get<0>(b));
    }
  });

  // In-place reduce-by-key w_all
  int p1 = 0;
  int p2 = 1;
  std::tuple< std::tuple<int32_t, int32_t>, float> tmp_p1;
  std::tuple< std::tuple<int32_t, int32_t>, float> tmp_p2;
  while (p2 < w_all.size()) {
    tmp_p1 = w_all[p1];
    tmp_p2 = w_all[p2];
    if (std::get<0>(std::get<0>(tmp_p1)) == std::get<0>(std::get<0>(tmp_p2)) && std::get<1>(std::get<0>(tmp_p1)) == std::get<1>(std::get<0>(tmp_p2))) {
      w_all[p1] = std::make_tuple(std::get<0>(tmp_p1), std::get<1>(tmp_p1) + std::get<1>(tmp_p2));
    } else {
      if (p2 - p1 > 1) {
        p1++;
        w_all[p1] = tmp_p2;
      } else {
        p1++;
      }
    }
    p2++;
  }
  int w_all_size = p1 + 1;

  C2_crd = (int32_t*)malloc(sizeof(int32_t) * w_all_size);
  float* w_vals = (float*)malloc(sizeof(float) * w_all_size);
  for (int32_t kw = 0; kw < w_all_size; kw++) {
    C2_crd[kw] = std::get<1>(std::get<0>(w_all[kw]));
    w_vals[kw] = std::get<1>(w_all[kw]);
  }

  int w1_pos[2];
  w1_pos[0] = 0;
  w1_pos[1] = w_all_size;
  int32_t iw = w1_pos[0];
  int32_t pw1_end = w1_pos[1];

  while (iw < pw1_end) {
    // int32_t k = w1_crd[kw];
    int32_t i = std::get<0>(std::get<0>(w_all[iw]));
    int32_t w1_segend = iw + 1;
    while (w1_segend < pw1_end && std::get<0>(std::get<0>(w_all[w1_segend])) == i) {
      w1_segend++;
    }
    C2_pos[i + 1] = w1_segend - iw;
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

double DCSC_DCSR_noacc_coord(taco_tensor_t *A, taco_tensor_t *B, taco_tensor_t* C, int32_t w_cap, int32_t warmup, int32_t repeat, bool bench = false, bool print = false) {
  // std::cout << "Capacity: " << w_cap << std::endl;
  for (int i = 0; i < warmup; i++) {
    compute(C,A,B,w_cap);
    if (bench) {
      free(C->vals);
      free(C->indices[1][0]);
      free(C->indices[1][1]);
    }
  }
  double start = clock();
  for (int i = 0; i < repeat; i++) {
    compute(C,A,B,w_cap);
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