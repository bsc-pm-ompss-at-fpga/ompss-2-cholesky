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

#include <sys/time.h>
#include <sys/times.h>
#include <math.h>

#define VERBOSE

#if USE_MKL
# include <mkl_lapacke.h>
# include <mkl.h>
#elif USE_OPENBLAS
# include <lapacke.h>
# include <cblas.h>
#elif USE_ARMPL
# include <armpl.h>
#else
# error No backend library found. See README for more information
#endif

#if defined(USE_FLOAT)
#  define type_t float
#  define gemm  cblas_sgemm
#  define trsm  cblas_strsm
#  define trmm  cblas_strmm
#  define syrk  cblas_ssyrk
#  define potrf LAPACKE_spotrf
#  define lacpy LAPACKE_slacpy
#  define lange LAPACKE_slange
#  define larnv LAPACKE_slarnv
#else
#  define type_t double
#  define gemm  cblas_dgemm
#  define trsm  cblas_dtrsm
#  define trmm  cblas_dtrmm
#  define syrk  cblas_dsyrk
#  define potrf LAPACKE_dpotrf
#  define lacpy LAPACKE_dlacpy
#  define lange LAPACKE_dlange
#  define larnv LAPACKE_dlarnv
#endif
#define LAPACKE_MAT_ORDER LAPACK_COL_MAJOR
#define CBLAS_MAT_ORDER   CblasColMajor
#define CBLAS_T           CblasTrans
#define CBLAS_NT          CblasNoTrans
#define CBLAS_LO          CblasLower
#define CBLAS_RI          CblasRight
#define CBLAS_NU          CblasNonUnit

float get_time()
{
   static double gtod_ref_time_sec = 0.0;

   struct timeval tv;
   gettimeofday(&tv, NULL);

   // If this is the first invocation of through dclock(), then initialize the
   // "reference time" global variable to the seconds field of the tv struct.
   if (gtod_ref_time_sec == 0.0)
      gtod_ref_time_sec = (double) tv.tv_sec;

   // Normalize the seconds field of the tv struct so that it is relative to the
   // "reference time" that was recorded during the first invocation of dclock().
   const double norm_sec = (double) tv.tv_sec - gtod_ref_time_sec;

   // Compute the number of seconds since the reference time.
   const double t = norm_sec + tv.tv_usec * 1.0e-6;

   return (float) t;
}

static type_t pow_di(type_t x, int n)
{
   type_t rv = 1.0;

   if (n < 0) {
      n = -n;
      x = 1.0 / x;
   }

   for (; n; n >>= 1, x *= x) {
      if (n & 1) rv *= x;
   }

   return rv;
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

   memset((void*)L1, 0, N*N*sizeof(type_t));
   memset((void*)L2, 0, N*N*sizeof(type_t));

   lacpy(LAPACKE_MAT_ORDER, ALL, N, N, A1, LDA, Residual, N);

   /* Dealing with L'L or U'U  */
   if (uplo == 'U'){
      lacpy(LAPACKE_MAT_ORDER, UP, N, N, A2, LDA, L1, N);
      lacpy(LAPACKE_MAT_ORDER, UP, N, N, A2, LDA, L2, N);
      trmm(CBLAS_MAT_ORDER, CBLAS_LO, CBLAS_LO, CBLAS_T, CBLAS_NU,
         N, N, alpha, L1, N, L2, N);
   } else {
      lacpy(LAPACKE_MAT_ORDER, LO, N, N, A2, LDA, L1, N);
      lacpy(LAPACKE_MAT_ORDER, LO, N, N, A2, LDA, L2, N);
      trmm(CBLAS_MAT_ORDER, CBLAS_RI, CBLAS_LO, CBLAS_T, CBLAS_NU,
         N, N, alpha, L1, N, L2, N);
   }

   /* Compute the Residual || A -L'L|| */
   for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
         Residual[j*N+i] = L2[j*N+i] - Residual[j*N+i];
      }
   }

   type_t Rnorm = lange(LAPACKE_MAT_ORDER, NORM, N, N, Residual, N);
   type_t Anorm = lange(LAPACKE_MAT_ORDER, NORM, N, N, A1, N);

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

   free(Residual);
   free(L1);
   free(L2);

   return info_factorization;
}

void initialize_matrix(const int n, const int ts, type_t *matrix)
{
   int ISEED[4] = {0,0,0,1};
   int intONE=1;

#ifdef VERBOSE
   printf("Initializing matrix with random values ...\n");
#endif

   for (int i = 0; i < n*n; i+=n) {
      larnv(intONE, &ISEED[0], n, &matrix[i]);
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

static void gather_block(const int N, const int ts, type_t *Alin, type_t *A)
{
   for (int i = 0; i < ts; i++) {
      for (int j = 0; j < ts; j++) {
         A[i*ts + j] = Alin[i*N + j];
      }
   }
}

static void scatter_block(const int N, const int ts, type_t *A, type_t *Alin)
{
   for (int i = 0; i < ts; i++) {
      for (int j = 0; j < ts; j++) {
         Alin[i*N + j] = A[i*ts + j];
      }
   }
}

static void convert_to_blocks(const int ts, const int DIM, const int N, type_t Alin[N][N], type_t *A[DIM][DIM])
{
   for (int i = 0; i < DIM; i++) {
      for (int j = 0; j < DIM; j++) {
         gather_block ( N, ts, &Alin[i*ts][j*ts], A[i][j]);
      }
   }
}

static void convert_to_linear(const int ts, const int DIM, const int N, type_t *A[DIM][DIM], type_t Alin[N][N])
{
   for (int i = 0; i < DIM; i++) {
      for (int j = 0; j < DIM; j++) {
         scatter_block ( N, ts, A[i][j], (type_t *) &Alin[i*ts][j*ts]);
      }
   }
}
