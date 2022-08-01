#include "dataloader.h"
#include "lib.h"
using namespace std;
const int w_order = 2;
const int HASH_SCAL = 1;
struct wspace {
  int32_t crd[2];
  float val;
};
// Initialize wspace
void init_wspace(wspace* w) {
    w->crd[0] = -1;
    //for (int i=0; i<w_order; i++) {
    //    w->crd[i] = -1;
    //}
}
// cmp function is generated when knowing the dimension of pos in wspace
int esc_cmp(const void *b, const void *a) {
  const int psize = w_order;
  for (int i=0; i<psize; i++) {
    if (((wspace*)b)->crd[i] == ((wspace*)a)->crd[i]) continue;
    return (((wspace*)b)->crd[i] - ((wspace*)a)->crd[i]);
  }
  return (((wspace*)b)->crd[psize-1] - ((wspace*)a)->crd[psize-1]);
}
int esc_cmp_rev(const void *a, const void *b) {
  const int psize = w_order;
  for (int i=0; i<psize; i++) {
    if (((wspace*)b)->crd[i] == ((wspace*)a)->crd[i]) continue;
    return (((wspace*)b)->crd[i] - ((wspace*)a)->crd[i]);
  }
  return (((wspace*)b)->crd[psize-1] - ((wspace*)a)->crd[psize-1]);
}
void print_wspace(wspace* w, const int len) {
  for(int i=0; i<len; i++) {
    cout<<"("<<w[i].crd[0]<<","<<w[i].crd[1]<<","<<w[i].val<<"),";
  }
  cout<<endl;
}
void print_coo(int32_t** C_COO_crd, float* C_COO_vals,const int len) {
  for(int i=0; i<len; i++) {
    cout<<"("<<C_COO_crd[0][i]<<","<<C_COO_crd[1][i]<<","<<C_COO_vals[i]<<"),";
  }
  cout<<endl;
}
int Enlarge(int32_t** C_COO_crd, float** C_COO_vals, int C_COO_capacity) {
  C_COO_capacity = C_COO_capacity * 2;
  C_COO_crd[0] = (int32_t*)realloc(C_COO_crd[0], sizeof(int32_t) * C_COO_capacity);
  C_COO_crd[1] = (int32_t*)realloc(C_COO_crd[1], sizeof(int32_t) * C_COO_capacity);
  *C_COO_vals = (float*)realloc(*C_COO_vals, sizeof(float) * C_COO_capacity);
  return C_COO_capacity;
}

int Merge(int32_t** COO_crd, float* COO_vals, int32_t COO_size, wspace* accumulator, int32_t accumulator_size, bool rev=false) {
  if (rev) {
    if (COO_size == 0) {
      for (int i=0; i<accumulator_size; i++) {
        COO_crd[0][i] = accumulator[accumulator_size-1-i].crd[0];
        COO_crd[1][i] = accumulator[accumulator_size-1-i].crd[1];
        COO_vals[i] = accumulator[accumulator_size-1-i].val;
      }
      return accumulator_size;
    }
    int32_t* tmp_COO_crd[2];
    float* tmp_COO_vals;
    tmp_COO_crd[0] = (int32_t*)malloc(sizeof(int32_t) * (accumulator_size + COO_size));
    tmp_COO_crd[1] = (int32_t*)malloc(sizeof(int32_t) * (accumulator_size + COO_size));
    tmp_COO_vals = (float*)malloc(sizeof(float) * (accumulator_size + COO_size));
    int accumulator_pointer = accumulator_size - 1;
    int content_pointer = 0;
    int target_pointer = 0;
    wspace tmp_con;
    while(accumulator_pointer >= 0 && content_pointer < COO_size) {
      tmp_con.crd[0] = COO_crd[0][content_pointer];
      tmp_con.crd[1] = COO_crd[1][content_pointer];
      if (esc_cmp(&accumulator[accumulator_pointer], &tmp_con) == 0) {
        tmp_COO_crd[0][target_pointer] = accumulator[accumulator_pointer].crd[0];
        tmp_COO_crd[1][target_pointer] = accumulator[accumulator_pointer].crd[1];
        tmp_COO_vals[target_pointer] = accumulator[accumulator_pointer].val + COO_vals[content_pointer];
        accumulator_pointer --;
        content_pointer ++;
        target_pointer ++;
      } else if (esc_cmp(&accumulator[accumulator_pointer], &tmp_con) < 0) {
        tmp_COO_crd[0][target_pointer] = accumulator[accumulator_pointer].crd[0];
        tmp_COO_crd[1][target_pointer] = accumulator[accumulator_pointer].crd[1];
        tmp_COO_vals[target_pointer] = accumulator[accumulator_pointer].val;
        accumulator_pointer --;
        target_pointer ++;
      } else {
        tmp_COO_crd[0][target_pointer] = COO_crd[0][content_pointer];
        tmp_COO_crd[1][target_pointer] = COO_crd[1][content_pointer];
        tmp_COO_vals[target_pointer] = COO_vals[content_pointer];
        content_pointer ++;
        target_pointer ++;
      }
    }
    while(accumulator_pointer >= 0) {
      tmp_COO_crd[0][target_pointer] = accumulator[accumulator_pointer].crd[0];
      tmp_COO_crd[1][target_pointer] = accumulator[accumulator_pointer].crd[1];
      tmp_COO_vals[target_pointer] = accumulator[accumulator_pointer].val;
      accumulator_pointer --;
      target_pointer ++;
    }
    while(content_pointer < COO_size) {
      tmp_COO_crd[0][target_pointer] = COO_crd[0][content_pointer];
      tmp_COO_crd[1][target_pointer] = COO_crd[1][content_pointer];
      tmp_COO_vals[target_pointer] = COO_vals[content_pointer];
      content_pointer ++;
      target_pointer ++;
    }
    for (int i = 0; i < target_pointer; i++) {
      COO_crd[0][i] = tmp_COO_crd[0][i];
      COO_crd[1][i] = tmp_COO_crd[1][i];
      COO_vals[i] = tmp_COO_vals[i];
    }
    free(tmp_COO_crd[0]);
    free(tmp_COO_crd[1]);
    free(tmp_COO_vals);
    return target_pointer;
  }
  else {
    if (COO_size == 0) {
      for (int i=0; i<accumulator_size; i++) {
        COO_crd[0][i] = accumulator[i].crd[0];
        COO_crd[1][i] = accumulator[i].crd[1];
        COO_vals[i] = accumulator[i].val;
      }
      return accumulator_size;
    }
    int32_t* tmp_COO_crd[2];
    float* tmp_COO_vals;
    tmp_COO_crd[0] = (int32_t*)malloc(sizeof(int32_t) * (accumulator_size + COO_size));
    tmp_COO_crd[1] = (int32_t*)malloc(sizeof(int32_t) * (accumulator_size + COO_size));
    tmp_COO_vals = (float*)malloc(sizeof(float) * (accumulator_size + COO_size));
    int accumulator_pointer = 0;
    int content_pointer = 0;
    int target_pointer = 0;
    wspace tmp_con;
    while(accumulator_pointer < accumulator_size && content_pointer < COO_size) {
      tmp_con.crd[0] = COO_crd[0][content_pointer];
      tmp_con.crd[1] = COO_crd[1][content_pointer];
      if (esc_cmp(&accumulator[accumulator_pointer], &tmp_con) == 0) {
        tmp_COO_crd[0][target_pointer] = accumulator[accumulator_pointer].crd[0];
        tmp_COO_crd[1][target_pointer] = accumulator[accumulator_pointer].crd[1];
        tmp_COO_vals[target_pointer] = accumulator[accumulator_pointer].val + COO_vals[content_pointer];
        accumulator_pointer ++;
        content_pointer ++;
        target_pointer ++;
      } else if (esc_cmp(&accumulator[accumulator_pointer], &tmp_con) < 0) {
        tmp_COO_crd[0][target_pointer] = accumulator[accumulator_pointer].crd[0];
        tmp_COO_crd[1][target_pointer] = accumulator[accumulator_pointer].crd[1];
        tmp_COO_vals[target_pointer] = accumulator[accumulator_pointer].val;
        accumulator_pointer ++;
        target_pointer ++;
      } else {
        tmp_COO_crd[0][target_pointer] = COO_crd[0][content_pointer];
        tmp_COO_crd[1][target_pointer] = COO_crd[1][content_pointer];
        tmp_COO_vals[target_pointer] = COO_vals[content_pointer];
        content_pointer ++;
        target_pointer ++;
      }
    }
    while(accumulator_pointer<accumulator_size) {
      tmp_COO_crd[0][target_pointer] = accumulator[accumulator_pointer].crd[0];
      tmp_COO_crd[1][target_pointer] = accumulator[accumulator_pointer].crd[1];
      tmp_COO_vals[target_pointer] = accumulator[accumulator_pointer].val;
      accumulator_pointer ++;
      target_pointer ++;
    }
    while(content_pointer<COO_size) {
      tmp_COO_crd[0][target_pointer] = COO_crd[0][content_pointer];
      tmp_COO_crd[1][target_pointer] = COO_crd[1][content_pointer];
      tmp_COO_vals[target_pointer] = COO_vals[content_pointer];
      content_pointer ++;
      target_pointer ++;
    }
    for (int i = 0; i < target_pointer; i++) {
      COO_crd[0][i] = tmp_COO_crd[0][i];
      COO_crd[1][i] = tmp_COO_crd[1][i];
      COO_vals[i] = tmp_COO_vals[i];
    }
    free(tmp_COO_crd[0]);
    free(tmp_COO_crd[1]);
    free(tmp_COO_vals);
    return target_pointer;
  }
}

int Sort(void * array, size_t size, bool rev) {
  if (rev) {
    qsort(array, size, sizeof(wspace), esc_cmp_rev);
    for (int i = size - 1; i >= 0; i--) {
      if(((wspace*)array)[i].crd[0] != -1) {
        return i + 1;
      }
    }
  } 
  qsort(array, size, sizeof(wspace), esc_cmp);
  return size;
}