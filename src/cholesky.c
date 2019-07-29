/*
* Copyright (c) 2017, BSC (Barcelona Supercomputing Center)
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*     * Redistributions of source code must retain the above copyright
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above copyright
*       notice, this list of conditions and the following disclaimer in the
*       documentation and/or other materials provided with the distribution.
*     * Neither the name of the <organization> nor the
*       names of its contributors may be used to endorse or promote products
*       derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY BSC ''AS IS'' AND ANY
* EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL <copyright holder> BE LIABLE FOR ANY
* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <assert.h>

#include "cholesky.h"
#include "cholesky.fpga.h"

#pragma omp target device(fpga) copy_inout([ts*ts]A)
#pragma omp task
void omp_potrf(type_t *A)
{
   #pragma HLS INLINE // off
   #pragma HLS array_partition variable=A block factor=ts
   for (int j = 0; j < ts; ++j) {
      for (int k = 0; k < j; ++k) {
         A[j*ts + j] -= A[k*ts + j]*A[k*ts + j];
      }
      //A[j][j] = sqrt(A[j][j]);
      const type_t A_jj = A[j*ts + j];
      type_t l = 0, h = A_jj, m;
      for (int i = 0; i < 256/*max_iters_to_converge*/; ++i) {
         m = (l+h)/2;
         if (m*m == A_jj) break;
         else if (m*m > A_jj) h = m;
         else l = m;
      }
      A[j*ts + j] = m;

      for (int i = j + 1; i < ts; ++i) {
         #pragma HLS pipeline II=1
         type_t tmp = A[j*ts + i];
         for (int k = 0; k < j; ++k) {
            tmp -= A[k*ts + i]*A[k*ts + j];
         }
         A[j*ts + i] = tmp/A[j*ts + j];
      }
   }
}

#pragma omp target device(fpga) copy_in([ts*ts]A) copy_inout([ts*ts]B)
#pragma omp task
void omp_trsm(const type_t *A, type_t *B)
{
   type_t tmp_row[ts];
   #pragma HLS INLINE // off
   #pragma HLS array_partition variable=A cyclic factor=4
   #pragma HLS array_partition variable=B cyclic factor=ts
   #pragma HLS array_partition variable=tmp_row cyclic factor=ts/2

   for (int k = 0; k < ts; ++k) {
      type_t temp = 1. / A[k*ts + k];
      for (int i__ = 0; i__ < ts; ++i__) {
         B[k*ts + i__] = tmp_row[i__] = temp * B[k*ts + i__];
      }
      for (int j = k + 1 ; j < ts; ++j) {
         #pragma HLS pipeline II=1
         temp = A[k*ts + j];
         for (int i__ = 0; i__ < ts; ++i__) {
            B[j*ts + i__] -= temp * tmp_row[i__];
         }
      }
   }
}

#pragma omp target device(fpga) copy_in([ts*ts]A) copy_inout([ts*ts]B)
#pragma omp task
void omp_syrk(const type_t *A, type_t *B)
{
   #pragma HLS INLINE // off
   #pragma HLS array_partition variable=A block factor=ts
   #pragma HLS array_partition variable=B cyclic factor=4

   for (int j = 0; j < ts; ++j) {
      for (int i__ = j; i__ < ts; ++i__) {
         #pragma HLS pipeline II=1
         type_t temp = B[j*ts + i__];
         for (int l = 0; l < ts; ++l) {
            temp += -A[l*ts + j] * A[l*ts + i__];
         }
         B[j*ts + i__] = temp;
      }
   }
}

#pragma omp target device(fpga) num_instances(3) copy_in([ts*ts]A, [ts*ts]B) copy_inout([ts*ts]C)
#pragma omp task
void omp_gemm(const type_t *A, const type_t *B, type_t *C)
{
   #pragma HLS INLINE // off
   #pragma HLS array_partition variable=A cyclic factor=4
   #pragma HLS array_partition variable=B cyclic factor=ts/2
   #pragma HLS array_partition variable=C cyclic factor=ts

   for (int k = 0; k < ts; ++k) {
      for (int i = 0; i < ts; ++i) {
         #pragma HLS pipeline II=1
         for (int j = 0; j < ts; ++j) {
            C[i*ts + j] += A[i*ts + k] * -B[k*ts + j];
         }
      }
   }
}

#pragma omp target device(fpga) copy_inout([nt*nt*ts*ts]A)
#pragma omp task
void cholesky_blocked(const int nt, type_t* A)
{
   for (int k = 0; k < nt; k++) {

      // Diagonal Block factorization
      omp_potrf( A + (k*nt + k)*ts*ts );
      #pragma omp taskwait

      // Triangular systems
      for (int i = k + 1; i < nt; i++) {
         omp_trsm( A + (k*nt + k)*ts*ts,
                   A + (k*nt + i)*ts*ts );
      }
      #pragma omp taskwait

      // Update trailing matrix
      for (int i = k + 1; i < nt; i++) {
         for (int j = k + 1; j < i; j++) {
            omp_gemm( A + (k*nt + i)*ts*ts,
                      A + (k*nt + j)*ts*ts,
                      A + (j*nt + i)*ts*ts );
         }
         omp_syrk( A + (k*nt + i)*ts*ts,
                   A + (i*nt + i)*ts*ts );
         #pragma omp taskwait
      }
   }
}

// Robust Check the factorization of the matrix A2
static int check_factorization(int N, type_t *A1, type_t *A2, int LDA, char uplo)
{
#ifdef VERBOSE
   printf ("Checking result ...\n");
#endif

   char NORM = 'I', ALL = 'A', UP = 'U', LO = 'L', TR = 'T', NU = 'N', RI = 'R';
   type_t alpha = 1.0;
   type_t const b = 2.0;
#ifdef USE_DOUBLE
   const int t = 53;
#else
   const int t = 24;
#endif
   type_t const eps = pow_di( b, -t );

   type_t *Residual = (type_t *)malloc(N*N*sizeof(type_t));
   type_t *L1       = (type_t *)malloc(N*N*sizeof(type_t));
   type_t *L2       = (type_t *)malloc(N*N*sizeof(type_t));
   type_t *work     = (type_t *)malloc(N*sizeof(type_t));

   memset((void*)L1, 0, N*N*sizeof(type_t));
   memset((void*)L2, 0, N*N*sizeof(type_t));

   lacpy(&ALL, &N, &N, A1, &LDA, Residual, &N);

   /* Dealing with L'L or U'U  */
   if (uplo == 'U'){
      lacpy(&UP, &N, &N, A2, &LDA, L1, &N);
      lacpy(&UP, &N, &N, A2, &LDA, L2, &N);
      trmm(CBLAS_MAT_ORDER, CBLAS_LF, CBLAS_UP, CBLAS_T, CBLAS_NU,
         N, N, alpha, L1, N, L2, N);
   } else {
      lacpy(&LO, &N, &N, A2, &LDA, L1, &N);
      lacpy(&LO, &N, &N, A2, &LDA, L2, &N);
      trmm(CBLAS_MAT_ORDER, CBLAS_RI, CBLAS_LO, CBLAS_T, CBLAS_NU,
         N, N, alpha, L1, N, L2, N);
   }

   /* Compute the Residual || A -L'L|| */
   for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
         Residual[j*N+i] = L2[j*N+i] - Residual[j*N+i];
      }
   }

   type_t Rnorm = lange(&NORM, &N, &N, Residual, &N, work);
   type_t Anorm = lange(&NORM, &N, &N, A1, &N, work);

   printf("==================================================\n");
   printf("Checking the Cholesky Factorization \n");
#ifdef VERBOSE
   printf("-- Rnorm = %e \n", Rnorm);
   printf("-- Anorm = %e \n", Anorm);
   printf("-- Anorm*N*eps = %e \n", Anorm*N*eps);
   printf("-- ||L'L-A||_oo/(||A||_oo.N.eps) = %e \n",Rnorm/(Anorm*N*eps));
#endif

   const int info_factorization = isnan(Rnorm/(Anorm*N*eps)) ||
      isinf(Rnorm/(Anorm*N*eps)) || (Rnorm/(Anorm*N*eps) > 60.0);

   if ( info_factorization ){
      fprintf(stderr, "\n-- Factorization is suspicious ! \n\n");
   } else {
      printf("\n-- Factorization is CORRECT ! \n\n");
   }

   free(work);
   free(L2);
   free(L1);
   free(Residual);

   return info_factorization;
}

void initialize_matrix(const int n, type_t *matrix)
{
   int ISEED[4] = {0,0,0,1};
   int intONE=1;

#ifdef VERBOSE
   printf("Initializing matrix with random values ...\n");
#endif

   for (int i = 0; i < n*n; i+=n) {
      larnv(&intONE, &ISEED[0], &n, &matrix[i]);
   }

   type_t a = (type_t)n;
   for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
         matrix[j*n + i] = matrix[j*n + i] + matrix[i*n + j];
         matrix[i*n + j] = matrix[j*n + i];
      }
      //add_to_diag
      matrix[i*n + i] += a;
   }
}

int main(int argc, char* argv[])
{
   char *result[3] = {"n/a","sucessful","UNSUCCESSFUL"};

   if ( argc < 3 ) {
      fprintf( stderr, "USAGE:\t%s <matrix size> [<check>]\n", argv[0] );
      exit( -1 );
   }
   const int  n = atoi(argv[1]); // matrix size
   int check    = argc > 2 ? atoi(argv[2]) : 1; // check result?
   const int nt = n / ts; // number of tiles
   if ( n % ts != 0 ) {
      fprintf( stderr, "ERROR:\t<matrix size> is not multiple of <block size>\n" );
      exit( -1 );
   }

   // Allocate matrix
   type_t * const matrix = (type_t *) malloc(n * n * sizeof(type_t));
   assert(matrix != NULL);

   // Init matrix
   initialize_matrix(n, matrix);

   type_t * original_matrix = NULL;
   if ( check ) {
      // Allocate matrix
      original_matrix = (type_t *) malloc(n * n * sizeof(type_t));
      assert(original_matrix != NULL);
      memcpy(original_matrix, matrix, n * n * sizeof(type_t));
   }

   // Allocate blocked matrix
   type_t *Ah[nt][nt];
   type_t *Ab;
   const size_t s = ts * ts * sizeof(type_t);
#ifdef USE_DMA_MEM
   Ab = (type_t *)nanos_fpga_alloc_dma_mem(s*nt*nt);
#else
   Ab = malloc(s*nt*nt);
#endif
   assert(Ab != NULL);

   for (int i = 0; i < nt; i++) {
      for (int j = 0; j < nt; j++) {
         Ah[i][j] = Ab + (i*nt + j)*ts*ts;
      }
   }

#ifdef VERBOSE
   printf ("Executing ...\n");
#endif

   convert_to_blocks(nt, n, (type_t(*)[n]) matrix, Ah);

   const float secs1 = get_time();
   cholesky_blocked(nt, Ab);

   const float secs2 = get_time();
   convert_to_linear(nt, n, Ah, (type_t (*)[n]) matrix);

   if ( check ) {
      const char uplo = 'L';
      if ( check_factorization(n, original_matrix, matrix, n, uplo) ) check++;
      free(original_matrix);
   }

   float time = secs2 - secs1;
   float gflops = (((1.0 / 3.0) * n * n * n) / (time));

   // Print results
   printf( "==================== RESULTS ===================== \n" );
   printf( "  Benchmark: %s (%s)\n", "Cholesky", "OmpSs" );
   printf( "  Elements type: %s\n", ELEM_T_STR );
#ifdef VERBOSE
   printf( "  Matrix size:           %dx%d\n", n, n);
   printf( "  Block size:            %dx%d\n", ts, ts);
#endif
   printf( "  Performance (gflops):  %f\n", gflops);
   printf( "  Execution time (secs): %f\n", time );
   printf( "================================================== \n" );

   // Free blocked matrix
#ifdef USE_DMA_MEM
   nanos_fpga_free_dma_mem(Ab);
#else
   free(Ab);
#endif

   // Free matrix
   free(matrix);

   return 0;
}
