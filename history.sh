make clean
make CPPFILE=bench_hash_ttm EXEPREFIX=bench-ttm ALGNAME=HASHIL
make CPPFILE=bench_hash_ttm EXEPREFIX=bench-ttm ALGNAME=HASH
make CPPFILE=bench_hash_mttkrp EXEPREFIX=bench-mttkrp ALGNAME=HASH
make CPPFILE=bench_hash_ttm EXEPREFIX=iteration-ttm ALGNAME=HASH
make CPPFILE=count_iteration_outer_cc_noT EXEPREFIX=count-outer-cc-noT ALGNAME=FLOPS
make CPPFILE=bench_coord_mttkrp EXEPREFIX=bench-mttkrp ALGNAME=COORDC
make CPPFILE=bench_coord_mttkrp EXEPREFIX=bench-mttkrp ALGNAME=COORDCF
make CPPFILE=bench_bucket_mttkrp EXEPREFIX=bench-mttkrp ALGNAME=BUCKET