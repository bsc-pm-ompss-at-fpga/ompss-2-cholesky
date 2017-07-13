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

#include <math.h>
#if USE_MKL
# include <mkl.h>
#elif USE_OPENBLAS
# include <lapacke.h>
#elif USE_ARMPL
//# include <armpl.h>
#else
# error No backend library found. See README for more information
#endif

#include <sys/time.h>
#include <sys/times.h>

#if defined(USE_FLOAT)
#  define type_t float
#  define gemm_  sgemm_
#  define trsm_  strsm_
#  define trmm_  strmm_
#  define syrk_  ssyrk_
#  define potrf_ spotrf_
#  define lacpy_ slacpy_
#  define lange_ slange_
#  define larnv_ slarnv_
#else
#  define type_t double
#  define gemm_  dgemm_
#  define trsm_  dtrsm_
#  define trmm_  dtrmm_
#  define syrk_  dsyrk_
#  define potrf_ dpotrf_
#  define lacpy_ dlacpy_
#  define lange_ dlange_
#  define larnv_ dlarnv_
#endif

void dgemm_ (const char *transa, const char *transb, int *l, int *n, int *m, double *alpha,
             const void *a, int *lda, void *b, int *ldb, double *beta, void *c, int *ldc);
void dtrsm_ (char *side, char *uplo, char *transa, char *diag, int *m, int *n, double *alpha,
             double *a, int *lda, double *b, int *ldb);
void dtrmm_ (char *side, char *uplo, char *transa, char *diag, int *m, int *n, double *alpha,
             double *a, int *lda, double *b, int *ldb);
void dsyrk_ (char *uplo, char *trans, int *n, int *k, double *alpha, double *a, int *lda,
             double *beta, double *c, int *ldc);
void sgemm_ (const char *transa, const char *transb, int *l, int *n, int *m, float *alpha,
             const void *a, int *lda, void *b, int *ldb, float *beta, void *c, int *ldc);
void strsm_ (char *side, char *uplo, char *transa, char *diag, int *m, int *n, float *alpha,
             float *a, int *lda, float *b, int *ldb);
void strmm_ (char *side, char *uplo, char *transa, char *diag, int *m, int *n, float *alpha,
             float *a, int *lda, float *b, int *ldb);
void ssyrk_ (char *uplo, char *trans, int *n, int *k, float *alpha, float *a, int *lda,
             float *beta, float *c, int *ldc);

enum blas_order_type {
   blas_rowmajor = 101,
   blas_colmajor = 102
};

enum blas_cmach_type {
   blas_base      = 151,
   blas_t         = 152,
   blas_rnd       = 153,
   blas_ieee      = 154,
   blas_emin      = 155,
   blas_emax      = 156,
   blas_eps       = 157,
   blas_prec      = 158,
   blas_underflow = 159,
   blas_overflow  = 160,
   blas_sfmin     = 161
};

enum blas_norm_type {
   blas_one_norm       = 171,
   blas_real_one_norm  = 172,
   blas_two_norm       = 173,
   blas_frobenius_norm = 174,
   blas_inf_norm       = 175,
   blas_real_inf_norm  = 176,
   blas_max_norm       = 177,
   blas_real_max_norm  = 178
};

static void BLAS_error(char *rname, int err, int val, int x)
{
   fprintf( stderr, "%s %d %d %d\n", rname, err, val, x );
   abort();
}

static void BLAS_ge_norm(enum blas_order_type order, enum blas_norm_type norm,
      const int m, const int n, const type_t *a, const int lda, type_t *res)
{
   char rname[] = "BLAS_ge_norm";

   if (order != blas_colmajor) BLAS_error( rname, -1, order, 0 );

   float anorm, v;
   if (norm == blas_frobenius_norm) {
      anorm = 0.0f;
      for (int j = n; j; --j) {
         for (int i = m; i; --i) {
            v = a[0];
            anorm += v * v;
            a++;
         }
         a += lda - m;
      }
      anorm = sqrt( anorm );
   } else if (norm == blas_inf_norm) {
      anorm = 0.0f;
      for (int i = 0; i < m; ++i) {
         v = 0.0f;
         for (int j = 0; j < n; ++j) {
            v += abs( a[i + j * lda] );
         }
         if (v > anorm)
            anorm = v;
      }
   } else {
      BLAS_error( rname, -2, norm, 0 );
      return;
   }

   if (res) *res = anorm;
}

static double BLAS_dpow_di(type_t x, int n)
{
   double rv = 1.0;

   if (n < 0) {
      n = -n;
      x = 1.0 / x;
   }

   for (; n; n >>= 1, x *= x) {
      if (n & 1)
         rv *= x;
   }

   return rv;
}

static float BLAS_spow_di(type_t x, int n)
{
   float rv = 1.0;

   if (n < 0) {
      n = -n;
      x = 1.0 / x;
   }

   for (; n; n >>= 1, x *= x) {
      if (n & 1)
         rv *= x;
   }

   return rv;
}

static double BLAS_dfpinfo(enum blas_cmach_type cmach)
{
   const double b = 2.0;
   const int t = 53, l = 1024, m = -1021;
   char rname[] = "BLAS_dfpinfo";

   // for (i = 0; i < t; ++i) eps *= half;
   const double eps = BLAS_dpow_di( b, -t );
   // for (i = 0; i >= m; --i) r *= half;
   const double r = BLAS_dpow_di( b, m-1 );

   double o = 1.0;
   o -= eps;
   // for (i = 0; i < l; ++i) o *= b;
   o = (o * BLAS_dpow_di( b, l-1 )) * b;

   switch (cmach) {
      case blas_eps: return eps;
      case blas_sfmin: return r;
      default:
         BLAS_error( rname, -1, cmach, 0 );
         break;
   }
   return 0.0;
}

static float BLAS_sfpinfo(enum blas_cmach_type cmach)
{
   const float b = 2.0;
   const int t = 24, l = 1024, m = -1021;
   char rname[] = "BLAS_dfpinfo";

   // for (i = 0; i < t; ++i) eps *= half;
   const float eps = BLAS_spow_di( b, -t );
   // for (i = 0; i >= m; --i) r *= half;
   const float r = BLAS_spow_di( b, m-1 );

   float o = 1.0;
   o -= eps;
   // for (i = 0; i < l; ++i) o *= b;
   o = (o * BLAS_spow_di( b, l-1 )) * b;

   switch (cmach) {
      case blas_eps: return eps;
      case blas_sfmin: return r;
      default:
         BLAS_error( rname, -1, cmach, 0 );
         break;
   }
   return 0.0;
}

void add_to_diag_hierarchical (type_t ** matrix, const int ts, const int nt, const float alpha)
{
   for (int i = 0; i < nt * ts; i++) {
      matrix[(i/ts) * nt + (i/ts)][(i%ts) * ts + (i%ts)] += alpha;
   }
}

void add_to_diag(type_t * matrix, const int n, const double alpha)
{
   for (int i = 0; i < n; i++) {
      matrix[ i + i * n ] += alpha;
   }
}

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

// Robust Check the factorization of the matrix A2
static int check_factorization(int N, type_t *A1, type_t *A2, int LDA, char uplo, type_t eps)
{
   char NORM = 'I', ALL = 'A', UP = 'U', LO = 'L', TR = 'T', NU = 'N', RI = 'R';

#ifdef VERBOSE
   printf ("Checking result ...\n");
#endif

   type_t *Residual = (type_t *)malloc(N*N*sizeof(type_t));
   type_t *L1       = (type_t *)malloc(N*N*sizeof(type_t));
   type_t *L2       = (type_t *)malloc(N*N*sizeof(type_t));
   type_t *work     = (type_t *)malloc(N*sizeof(type_t));

   memset((void*)L1, 0, N*N*sizeof(type_t));
   memset((void*)L2, 0, N*N*sizeof(type_t));

   type_t alpha= 1.0;

   lacpy_(&ALL, &N, &N, A1, &LDA, Residual, &N);

   /* Dealing with L'L or U'U  */
   if (uplo == 'U'){
      lacpy_(&UP, &N, &N, A2, &LDA, L1, &N);
      lacpy_(&UP, &N, &N, A2, &LDA, L2, &N);
      trmm_(&LO, &uplo, &TR, &NU, &N, &N, &alpha, L1, &N, L2, &N);
   }
   else{
      lacpy_(&LO, &N, &N, A2, &LDA, L1, &N);
      lacpy_(&LO, &N, &N, A2, &LDA, L2, &N);
      trmm_(&RI, &LO, &TR, &NU, &N, &N, &alpha, L1, &N, L2, &N);
   }

   /* Compute the Residual || A -L'L|| */
   for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
         Residual[j*N+i] = L2[j*N+i] - Residual[j*N+i];
      }
   }

   type_t Rnorm = lange_(&NORM, &N, &N, Residual, &N, work);
   type_t Anorm = lange_(&NORM, &N, &N, A1, &N, work);

#ifdef VERBOSE
   printf("============\n");
   printf("Checking the Cholesky Factorization \n");
   printf("-- Rnorm = %e \n", Rnorm);
   printf("-- Anorm = %e \n", Anorm);
   printf("-- Anorm*N*eps = %e \n", Anorm*N*eps);
   printf("-- ||L'L-A||_oo/(||A||_oo.N.eps) = %e \n",Rnorm/(Anorm*N*eps));
#endif

   const int info_factorization = isnan(Rnorm/(Anorm*N*eps)) ||
      isinf(Rnorm/(Anorm*N*eps)) || (Rnorm/(Anorm*N*eps) > 60.0);

#ifdef VERBOSE
   if ( info_factorization){
      printf("\n-- Factorization is suspicious ! \n\n");
   }
   else{
      printf("\n-- Factorization is CORRECT ! \n\n");
   }
#endif

   free(Residual); free(L1); free(L2); free(work);

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
      larnv_(&intONE, &ISEED[0], &n, &matrix[i]);
   }

   for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
         matrix[j*n + i] = matrix[j*n + i] + matrix[i*n + j];
         matrix[i*n + j] = matrix[j*n + i];
      }
   }

   add_to_diag(matrix, n, (double) n);
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
   for (int i = 0; i < DIM; i++)
      for (int j = 0; j < DIM; j++) {
         gather_block ( N, ts, &Alin[i*ts][j*ts], A[i][j]);
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

static type_t * malloc_block (const int ts)
{
   type_t * const block = (type_t *) malloc(ts * ts * sizeof(type_t));
   assert(block != NULL);
   return block;
}
