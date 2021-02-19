# set environment vars for Intel® Fortran Compiler
source /opt/intel/oneapi/setvars.sh intel64

export CPATH="/opt/intel/oneapi/mkl/2021.1.1/include:/opt/intel/oneapi/mkl/2021.1.1/include/intel64/lp64:/Users/hans/F_code_mkl:$CPATH"
export CPATH="/usr/local/Cellar/mpich/3.3.2_1/include:$CPATH"
export CPATH="/opt/intel/oneapi/mkl/2021.1.1/lib:$C_PATH"
export DYLD_LIBRARY_PATH="/usr/local/Cellar/mpich/3.3.2_1/lib:$DYLD_LIBRARY_PATH"
export DYLD_LIBRARY_PATH="/opt/intel/oneapi/mkl/2021.1.1/include/intel64/lp64:/Users/hans/F_code_mkl:/opt/intel/oneapi/mkl/2021.1.1/lib:$DYLD_LIBRARY_PATH"
export LIBRARY_PATH="/usr/local/Cellar/mpich/3.3.2_1/lib:$LIBRARY_PATH"
export LIBRARY_PATH="/opt/intel/oneapi/mkl/2021.1.1/include/intel64/lp64:/Users/hans/F_code_mkl:$LIBRARY_PATH"
export LIBRARY_PATH="/opt/intel/oneapi/mkl/2021.1.1/lib:$LIBRARY_PATH"
