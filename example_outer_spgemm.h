#include "Coord.h"
#include "Insert.h"
#include "Sort.h"
#include "Merge.h"

// Allocate
Coord Acc(2,cap);
AllArray All(2,cap);
int point[2];

// Insert-Sort-Merge
for (int k = 0; k < K; k++) {
  for (int iB = B.pos[k]; iB < B.pos[k + 1]; iB++) {
    int i = B.crd[iB];
    point[0] = i;
    for (int jC = C.pos[k]; jC < C.pos[k + 1]; jC++) {
      int j = C.crd[jC];
      point[1] = j;
      Insert(w_point, (B.vals[iB] * C.vals[jC]), &Acc);
      if (Acc.full) {
        All.realloc(Acc.size);
        Sort(&Acc);
        Merge(&Acc, &All);
        Acc.refresh();
        Insert(w_point, (B.vals[iB] * C.vals[jC]), &Acc);
      }
    }
  }
}
if (Acc.size > 0) {
  All.realloc(Acc.size);
  Sort(&Acc);
  Merge(&Acc, &All);
}

// Compress
A.crd = All.crd[1];
A.vals = All.vals;
int* A.pos = (int*)calloc(I + 1, sizeof(int));
int iw = 0;
while (iw < All.size) {
  int i = All.crd[0][iw];
  int segend = iw + 1;
  while (segend < All.size && All.crd[0][segend] == i) {
    segend++;
  }
  A.pos[i + 1] = segend - iw;
  iw = segend;
}
int cnt = 0;
for (int pA = 1; pA < I + 1; pA++) {
  cnt += A.pos[pA];
  A.pos[pA] = cnt;
}



