#include "cholesky.h"
#include "cholesky.fpga.h"

#if defined(OPENBLAS_IMPL) || defined(POTRF_SMP)
#pragma oss task inout([ts*ts]A)
#else
#pragma oss task device(fpga) inout([ts*ts]A)
#endif
void omp_potrf(type_t *A)
{
#if defined(OPENBLAS_IMPL) || defined(POTRF_SMP)
#else
   #pragma HLS inline
   #pragma HLS array_partition variable=A cyclic factor=FPGA_PWIDTH/64
   for (int j = 0; j < ts; ++j) {
      type_t tmp = A[j*ts + j];
      for (int k = 0; k < j; ++k) {
         #pragma HLS pipeline II=1
         type_t Akj = A[k*ts + j];
         tmp -= Akj*Akj;
      }

      A[j*ts + j] = sqrtf(tmp);

      for (int i = j + 1; i < ts; ++i) {
         type_t tmp = A[j*ts + i];
         for (int k = 0; k < j; ++k) {
            #pragma HLS pipeline II=1
            tmp -= A[k*ts + i]*A[k*ts + j];
         }
         A[j*ts + i] = tmp/A[j*ts + j];
      }
   }
#endif
}

#pragma oss task device(fpga) num_instances(TRSM_NUMACCS) copy_deps in([ts*ts]A) inout([ts*ts]B)
void omp_trsm(const type_t *A, type_t *B)
{
   #pragma HLS inline
   #pragma HLS array_partition variable=A cyclic factor=FPGA_PWIDTH/64
   #pragma HLS array_partition variable=B cyclic factor=ts/FPGA_OTHER_II
   #pragma HLS array_partition variable=tmp_row cyclic factor=ts/(2*FPGA_OTHER_II)
   type_t tmp_row[ts];

   for (int k = 0; k < ts; ++k) {
      type_t temp = 1. / A[k*ts + k];
      for (int i = 0; i < ts; ++i) {
         #pragma HLS unroll factor=ts/FPGA_OTHER_II
         #pragma HLS pipeline II=1
         //Sometimes Vivado HLS doesn't achieve II=1 because it detects
         //some false dependence on B, this fixes the issue. Same for the other loop
         #pragma HLS DEPENDENCE variable=B inter false
         B[k*ts + i] = tmp_row[i] = temp * B[k*ts + i];
      }

      for (int j = k + 1; j < ts; ++j) {
         #pragma HLS pipeline II=FPGA_OTHER_II
         #pragma HLS DEPENDENCE variable=B inter false
         for (int i = 0; i < ts; ++i) {
            B[j*ts + i] -= A[k*ts + j] * tmp_row[i];
         }
      }
   }
}

#pragma oss task device(fpga) num_instances(SYRK_NUMACCS) copy_deps in([ts*ts]A) inout([ts*ts]B)
void omp_syrk(const type_t *A, type_t *B)
{
   #pragma HLS inline
   #pragma HLS array_partition variable=A cyclic factor=ts/FPGA_OTHER_II
   #pragma HLS array_partition variable=B cyclic factor=ts/FPGA_OTHER_II

   for (int k = 0; k < ts; ++k) {
      for (int i = 0; i < ts; ++i) {
         #pragma HLS pipeline II=FPGA_OTHER_II
         for (int j = 0; j < ts; ++j) {
            //NOTE: Instead of reduce the 'i' iterations, multiply by 0
            B[i*ts + j] += -A[k*ts + i] * (j < i ? 0 : A[k*ts + j]);
         }
      }
   }
}

#pragma oss task device(fpga) num_instances(GEMM_NUMACCS) copy_deps in([ts*ts]A, [ts*ts]B) inout([ts*ts]C)
void omp_gemm(const type_t *A, const type_t *B, type_t *C)
{
   #pragma HLS inline
   #pragma HLS array_partition variable=A cyclic factor=ts/(2*FPGA_GEMM_II)
   #pragma HLS array_partition variable=B cyclic factor=FPGA_PWIDTH/64
   #pragma HLS array_partition variable=C cyclic factor=ts/FPGA_GEMM_II
   #ifdef USE_URAM
   #pragma HLS resource variable=A core=XPM_MEMORY uram
   #pragma HLS resource variable=B core=XPM_MEMORY uram
   #endif

   for (int k = 0; k < ts; ++k) {
      for (int i = 0; i < ts; ++i) {
         #pragma HLS pipeline II=FPGA_GEMM_II
         for (int j = 0; j < ts; ++j) {
            C[i*ts + j] += A[k*ts + j] * -B[k*ts + i];
         }
      }
   }
}

#pragma oss task device(fpga) inout([nt*nt*ts*ts]A)
void cholesky_blocked(const int nt, type_t* A)
{
   for (int k = 0; k < nt; k++) {

      // Diagonal Block factorization
      omp_potrf( A + (k*nt + k)*ts*ts );

      // Triangular systems
      // Create in inverse order because Picos wakes up ready
      // chain tasks in that order
      for (int i = nt-1; i >= k+1; i--) {
         omp_trsm( A + (k*nt + k)*ts*ts,
                   A + (k*nt + i)*ts*ts );
      }

      // Update trailing matrix
      for (int i = k + 1; i < nt; i++) {
         for (int j = k + 1; j < i; j++) {
            omp_gemm( A + (k*nt + i)*ts*ts,
                      A + (k*nt + j)*ts*ts,
                      A + (j*nt + i)*ts*ts );
         }
         omp_syrk( A + (k*nt + i)*ts*ts,
                   A + (i*nt + i)*ts*ts );
      }
   }
   #pragma oss taskwait
}
