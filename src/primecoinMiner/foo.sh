#!/bin/sh
#clang -S -emit-llvm -o test.ll -x cl foo.cl -O3 -target NVPTX -Dcl_clang_storage_class_specifiers -Dcl_khr_fp64 -Dcles_khr_int64  
clang -S -emit-llvm -o test.ll -x cl foo.cl -O3 -target r600-- -Dcl_clang_storage_class_specifiers -Dcl_khr_fp64 -Dcles_khr_int64 -mcpu=pitcairn -Xclang -mlink-bitcode-file -Xclang pitcairn-r600--.bc -include "clc/clc.h"
llvm-link -o foo.ll -S test.ll 
llc foo.ll -march=r600 -mcpu=pitcairn -O3 -o -
#-debug
#llc foo.ll -march=nvptx64 -O3 -o -




