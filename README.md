# Cholesky

**Name**: Cholesky Factorization Kernel  
**Contact Person**: PM Group, pm-tools@bsc.es  
**License Agreement**: GPL  
**Platform**: OmpSs  

### Description
This application performs a cholesky decomposition/factorization over a square matrix.
The matrix is distributed by blocks of contiguous memory.

The task implementation requires support for one external library that implements the mathematical operations. Supported ones are:
 - [Intel MKL](https://software.intel.com/en-us/mkl)
 - [Open BLAS](http://www.openblas.net/)
 - [ARM PL](https://developer.arm.com/products/software-development-tools/hpc/arm-performance-libraries)

### Build instructions
Clone the repository:
```
git clone https://pm.bsc.es/gitlab/applications/ompss/cholesky.git
cd cholesky
```

Build the application binaries:
```
make
```
##### Build variables
You can change the build process defining or modifying some environment variables.
The supported ones are:
  - `CFLAGS`
  - `LDFLAGS`
  - `MCC`. If not defined, the default value is: `mcc`. However, for SMP machines we recommend the use of `smpcc`.
  - `CROSS_COMPILE`
  - `MKL_DIR`. Installation directory of MKL library. The default value is: `$MKLROOT`.
    - `MKL_INC_DIR`. Installation directory of includes for MKL library. The default value is: `$MKL_DIR/include`.
    - `MKL_LIB_DIR`. Installation directory of OS libraries for MKL library. The default value is: `$MKL_DIR/lib`.
  - `OPENBLAS_DIR`. Installation directory of OpenBLAS library. The default value is: `$OPENBLAS_HOME`.
    - `OPENBLAS_INC_DIR`. Installation directory of includes for OpenBLAS library. The default value is: `$OPENBLAS_DIR/include`.
    - `OPENBLAS_LIB_DIR`. Installation directory of OS libraries for OpenBLAS library. The default value is: `$OPENBLAS_DIR/lib`.
  - `ARMPL_DIR`. Installation directory of ARM PL library. The default value is: `$ARMPL_HOME`.
    - `ARMPL_INC_DIR`. Installation directory of includes for ARM PL library. The default value is: `$ARMPL_DIR/include`.
    - `ARMPL_LIB_DIR`. Installation directory of OS libraries for ARM PL library. The default value is: `$ARMPL_DIR/lib`.

Note that in order to compile the application either `MKL_DIR`, `OPENBLAS_DIR` or `ARMPL_DIR` (or the derivate variables) must point to a valid installation.

For example, the build step to cross-compile the application for ARM using the `smpcc` profile and OpenBLAS library may be:
```
export MCC=smpcc
export CROSS_COMPILE=arm-linux-gnueabihf-
export ARMPL_DIR=/opt/install-arm/openblas
make
```


### Run instructions
The name of each binary file created by build step ends with a suffix which determines the version:
 - program-p: performance version
 - program-i: instrumented version
 - program-d: debug version

All versions use the same arguments structure:
```
./cholesky <matrix size> <block size> [<check>]
```
where:
 - `matrix size` is the dimension of the matrices. (Mandatory)
 - `block size` is the dimension of the matrices sub-blocks. (Mandatory)
 - `check` defines if the result must be checked. Default is: TRUE. (Optional)
