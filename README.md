# SpGEMM

## IndexStmt

With workspace:

```
assemble(forall(qi, where(forall(qk, C2_nnz(qi) += cast<int32_t>(qw(qk))), forall(qj, forall(qk, qw(qk) += (A_struct(qi,qj) * B_struct(qj,qk))))), CPUThread, NoRaces), forall(i, where(forall(k, C(i,k) = w(k)), forall(j, forall(k, w(k) += A(i,j) * B(j,k)))), CPUThread, NoRaces))
```

Without workspace:

```
assemble(forall(qi, forall(qk, where(C2_nnz(qi) += cast<int32_t>(qtqjC), forall(qj, qtqjC += (A_struct(qi,qj) * B_struct(qj,qk))))), CPUThread, NoRaces), forall(i, forall(k, forall(j, C(i,k) += A(i,j) * B(j,k))), CPUThread, NoRaces))
```
