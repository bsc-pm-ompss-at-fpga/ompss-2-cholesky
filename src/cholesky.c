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

#pragma omp target device(smp) copy_deps
#pragma omp task inout([ts]A)
void omp_potrf(type_t (*A)[ts])
{
   static const char L = 'L';
   int info;
   potrf(&L, &ts, (type_t *)A, &ts, &info);
}

#pragma omp target device(fpga) copy_deps onto(2)
#pragma omp task in([ts]A) inout([ts]B)
void omp_trsm(type_t (*A)[ts], type_t (*B)[ts])
{
#if 0
   trsm(CBLAS_MAT_ORDER, CBLAS_RI, CBLAS_LO, CBLAS_T, CBLAS_NU,
      ts, ts, 1.0, (type_t *)A, ts, (type_t *)B, ts);
#else
   type_t tmp_row[ts];
   for (int k = 0; k < ts; ++k) {
      type_t temp = 1. / A[k][k];
      for (int i__ = 0; i__ < ts; ++i__) {
         B[k][i__] = tmp_row[i__] = temp * B[k][i__];
      }
      for (int j = k + 1 ; j < ts; ++j) {
         temp = A[k][j];
         for (int i__ = 0; i__ < ts; ++i__) {
            B[j][i__] -= temp * tmp_row[i__];
         }
      }
   }
#endif
}

#pragma omp target device(fpga) copy_deps onto(1)
#pragma omp task in([ts]A) inout([ts]B)
void omp_syrk(type_t (*A)[ts], type_t (*B)[ts])
{
#if 0
   syrk(CBLAS_MAT_ORDER, CBLAS_LO, CBLAS_NT,
      ts, ts, -1.0, (type_t *)A, ts, 1.0, (type_t *)B, ts);
#else
   for (int j = 0; j < ts; ++j) {
      for (int i__ = j; i__ < ts; ++i__) {
         type_t temp = B[j][i__];
         for (int l = 0; l < ts; ++l) {
            temp += -A[l][j] * A[l][i__];
         }
         B[j][i__] = temp;
      }
   }
#endif
}

#pragma omp target device(fpga) copy_deps onto(0)
#pragma omp task in([ts]A, [ts]B) inout([ts]C)
void omp_gemm(type_t (*A)[ts], type_t (*B)[ts], type_t (*C)[ts])
{
#if 0
   gemm(CBLAS_MAT_ORDER, CBLAS_NT, CBLAS_T,
      ts, ts, ts, -1.0, (type_t *)A, ts, (type_t *)B, ts, 1.0, (type_t *)C, ts);
#else
   for (int j = 0; j < ts; ++j) {
      for (int i__ = 0; i__ < ts; ++i__) {
         type_t temp = C[j][i__];
         for (int l = 0; l < ts; ++l) {
            temp += -B[l][j] * A[l][i__];
         }
         C[j][i__] = temp;
      }
   }
#endif
}

void cholesky_blocked(const int nt, type_t* Ah[nt][nt])
{
   for (int k = 0; k < nt; k++) {

      // Diagonal Block factorization
      omp_potrf( (type_t (*)[ts])Ah[k][k] );

      // Triangular systems
      for (int i = k + 1; i < nt; i++) {
         omp_trsm( (type_t (*)[ts])Ah[k][k],
                   (type_t (*)[ts])Ah[k][i] );
      }

      // Update trailing matrix
      for (int i = k + 1; i < nt; i++) {
         for (int j = k + 1; j < i; j++) {
            omp_gemm( (type_t (*)[ts])Ah[k][i],
                      (type_t (*)[ts])Ah[k][j],
                      (type_t (*)[ts])Ah[j][i] );
         }
         omp_syrk( (type_t (*)[ts])Ah[k][i],
                   (type_t (*)[ts])Ah[i][i] );
      }

   }
   #pragma omp taskwait
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
#ifdef USE_FLOAT
   const int t = 24;
#else
   const int t = 53;
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

   for (int i = 0; i < nt; i++) {
      for (int j = 0; j < nt; j++) {
         Ah[i][j] = malloc(ts * ts * sizeof(type_t));
         assert(Ah[i][j] != NULL);
      }
   }

#ifdef VERBOSE
   printf ("Executing ...\n");
#endif

   convert_to_blocks(nt, n, (type_t(*)[n]) matrix, Ah);

   const float secs1 = get_time();
   cholesky_blocked(nt, (type_t* (*)[nt]) Ah);

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
#ifdef VERBOSE
   printf( "  Matrix size:           %dx%d\n", n, n);
   printf( "  Block size:            %dx%d\n", ts, ts);
#endif
   printf( "  Performance (gflops):  %f\n", gflops);
   printf( "  Execution time (secs): %f\n", time );
   printf( "================================================== \n" );

   // Free blocked matrix
   for (int i = 0; i < nt; i++) {
      for (int j = 0; j < nt; j++) {
         assert(Ah[i][j] != NULL);
         free(Ah[i][j]);
      }
   }

   // Free matrix
   free(matrix);

   return 0;
}
