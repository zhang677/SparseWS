# How to run

```
make clean
make 
./test
```

Change the header file included in `main.cpp` to test different methods.

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

## Design Phyilosophy

1. Keep modularity in mind
2. Keep generalization to multi-thread in mind
3. Keep consistent with current TACO assumption (code is generated "as is", tensor level type, ...)
4. Reuse current TACO's generated code as much as possible

## Current Design Principles

1. No temporary reduction on the workspace. The workspace will be written back to the output only if all the data points in the output have been prepared.
2. No assumption of the order of input data. Method 4 fails in CSR_CSR, because it assumes the output workspace array is inserted in the lexicographical order of coordinate.

## Development State
- [ ] Method 4 (coordinate) needs merge sort
- [ ] Method 4 (hash) needs wspace struct in a more modular way