#!/bin/sh
scl enable devtoolset-7 bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/src/tbb/build/linux_intel64_gcc_cc7.3.1_libc2.17_kernel3.10.0_release/
make clean
make -f makefile
echo miao