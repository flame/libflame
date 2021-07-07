#ifndef FLAME_GBENCH_MAIN_H
#define FLAME_GBENCH_MAIN_H

#include <stdio.h>
#include "complex.h"
#include <gtest/gtest.h>
#include <dlfcn.h>

/* macros */
 #define NUM_SUB_TESTS (4)
 #define FLAME_BIG_MATRIX_SIZE (2000)
 #define FLAME_GBENCH_FORCED_ITERATION_COUNT (20)
 
 #define MODE_RANGE (2)
 #define MODE_COMBINATIONAL (1)
 #define MODE_DISCRETE (0)
 

/* Global variables declaration */

#ifndef _DEFINED_DCOMPLEX
#define _DEFINED_DCOMPLEX
typedef struct dcomplex
{
  double real, imag;
} dcomplex;
#endif


typedef struct { float real, imag; } _lapack_complex_float;
typedef struct { double real, imag; } _lapack_complex_double;


#define lapack_complex_float    float _Complex
#define lapack_complex_double   double _Complex

typedef struct Lin_driver_paramlist_t
{
   int nparams; // algorithm parameters - applicable for LIN driver APIs only
   int nerrbnds; // number of error bounds 1, 2, 3 - applicable for LIN driver APIs only
   char fact;  // Must be 'F' or 'N'.
   
} Lin_driver_paramlist;
typedef struct Lin_solver_paramlist_t
{
   int mode; // Any one of these:- 0: discrete, 1: Combinational, 2:Range with steps of increment of Matrix sizes
   int num_ranges; // number of ranges to run
   int m_range_start;
   int m_range_end;
   int m_range_step_size;
   int n_range_start;
   int n_range_end;
   int n_range_step_size;
   int num_tests;
   int matrix_layout;
   char Uplo; 
   char transr; // Must be 'N' or 'T' or 'C'.
   int m;
   int n;  // The order of A; the number of rows in B
   int nrhs; // number of rhight hand sides
   int lda; //  leading dimension of the array a
   int ldb; //  leading dimension of the array b
   int ldab;  //  leading dimension of the array ab
   int kl; // number of subdiagonals 
   int ku; // number of superdiagonals
   int kd; // number of super or sub diagonals
   char diag; // flag to indicate unit diagonal
   
   // below params are used only by Lin solver driver APIs.
   char fact;  // Must be 'F', 'N', or 'E'.
   char equed; // Must be 'N', 'R'. 'C', 'B'              
   char symm; // if symmetric 'S' or Hermitian 'H'
   float solver_threhold;// threshold to verify PASS/FAIL criteria
   char equed_porfsx; // Must be 'N', 'Y'.
   int  n_err_bnds_porfsx;
   int  nparams_porfsx;
   char norm_gbcon; // norm param for gbcon API
   int kl_gbcon; // number of subdiagonals 
   int ku_gbcon; // number of superdiagonals
   int ldab_gbcon;  //  leading dimension of the array ab
} Lin_solver_paramlist;

extern Lin_solver_paramlist lin_solver_paramslist[30];
extern Lin_driver_paramlist lin_driver_paramslist[NUM_SUB_TESTS];


inline int Circular_Increment_Index( int idx) {
   idx = idx+1;
   if (idx >= NUM_SUB_TESTS){
      idx = 0;
   }
   return idx;
}
/* struct to hold eigen parameters */
typedef struct EIG_paramlist_t
{
   int matrix_layout;
   char trans; // Must be 'N' or 'T' or 'C'.
   char uplo; // Must be 'U' or 'L' 
   char job; // Must be 'N', 'P', 'S' or 'B'
   char vect; // Vector must be 'Q' or  'P'
   int m;   // 
   int n;  // The order of A; the number of rows in B
   int p; //
   int nrhs; // number of rhight hand sides
   int lda; //  leading dimension of the array a
   int ldb; //  leading dimension of the array b
   int nb;  //  leading dimension of the array ab
   int ldt; // number of subdiagonals
   int k;
   int isgn;
   char compz;
   int kb;
   int itype;
   char vect_rd;
   char side;
   char job_seqr; // Must be 'E', 'S'
   char eigsrc; // Must be 'Q' or  'N'.
   char initv; // Must be 'N' or 'U'.
   char norm;
   char diag;
   char storev;
   int tsize;
   int threshold_value; // threshold value for EIG
}EIG_paramlist;

extern EIG_paramlist eig_paramslist[NUM_SUB_TESTS];

/* struct to hold eigen parameters */
typedef struct EIG_Non_symmetric_paramlist_t
{
   char howmny; // Must be 'A' or 'B' or 'S'.
   char initv; // Must be 'N' or 'U'.
   char job_seqr; // Must be 'E', 'S'
   char eigsrc; // Must be 'Q' or  'N'.
   char side;  // Must be 'R' or 'L' or 'B'.
   char job; //  Must be 'E' or 'V' or 'B'.
   char howmny_trsna; // Must be 'A' or 'S'.

   /* used params for trsen API  */
   char job_trsen;// must bie 'N' or 'E' or 'V' or 'B'.
   char compq; // Must be 'V' or 'N' .

   /* used params for trsyl API  */
   char trana_real; //  Must be 'N' or 'T' or 'C'.
   char trana_complex; //  Must be 'N' or 'T' or 'C'.
   char tranb_real; //  Must be 'N' or 'T' or 'C'.
   char tranb_complex; //  Must be 'N' or 'T' or 'C'.
   int isgn; //  +1 or -1
   
   /* Thresholds for the APIs  */
   float gghrd_threshold; // threshold for the gghrd API
   float ggbal_threshold; // threshold for the ggbal API
   float GenNonSymEigProblem_threshold; // threshold for the ggbal API

   char compq_hgeqz; // Must be 'I' or 'V' or 'N'
   char compz_hgeqz; // Must be 'I' or 'V' or 'N'
   
	char side_tgevc; // Must be 'R', 'L', or 'B'.
	char jobvsl; // Must be 'N', or 'V'.
	char jobvsr; // Must be 'N', or 'V'.
	char sort; // Must be 'N', or 'S'.
	char sense_ggesx;// must be 'N' or 'E' or 'V' or 'B'.
	char balance_ggevx;// must be 'N' or 'P' or 'S' or 'B'.
	char sense_ggevx;// must be 'N' or 'E' or 'V' or 'B'.
	char sort_gees; // Must be 'N', or 'S'.
	int  wantz; // Must be 1 or 0
	int  wantq; // Must be 1 or 0
    int  tgsen_ijob; // Must be between 0 to 5
	char unmhr_trans; // Must be N or C
}EIG_Non_symmetric_paramlist;

extern EIG_Non_symmetric_paramlist eig_non_sym_paramslist[NUM_SUB_TESTS];

/* struct to hold SVD parameters */
typedef struct SVD_paramlist_t
{
	int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
	char jobu; // Must be 'U' or 'N'.
	char jobv; // Must be 'V' or 'N'.
	char jobq; // Must be 'Q' or 'N'.
	int m; // The number of rows of the matrix A
	int p; // The number of rows of the matrix B
	int n; // The number of columns of the matrices A and B
    float tola;
    float tolb;
    char jobu_gesvd; // Must be 'A', 'S', 'O', or 'N'. 
    char jobvt_gesvd; // Must be 'A', 'S', 'O', or 'N'. 

	char joba_gejsv; //  Must be 'C', 'E', 'F', 'G', 'A', or 'R'.
	char jobu_gejsv; // Must be 'U', 'F', 'W', or 'N'.
	char jobv_gejsv; // Must be 'V', 'J', 'W', or 'N'.
	char jobr_gejsv; // Must be 'N' or 'R'.
	char jobt_gejsv; // Must be 'T' or 'N'.
	char jobp_gejsv; //  Must be 'P' or 'N'.
	int m_gejsv; // The number of rows of the matrix A
	int n_gejsv; // The number of rows of the matrix B

	/* Parameters for 'gesvj' API  */
	char joba_gesvj; //  Must be 'L', 'U' or 'G'.
	char jobu_gesvj; // Must be 'U', 'C' or 'N'.
	char jobv_gesvj; // Must be 'V', 'A' or 'N'.
	int m_gesvj; // The number of rows of the matrix A
	int n_gesvj; // The number of rows of the matrix B
	int mv_gesvj;
	float ctol_gesvj; // convergence of threshold

	/* Parameters for 'gesvdx' API  */
	char jobu_gesvdx; //  Must be 'V', or 'N'.
	char jobvt_gesvdx; // Must be 'V', or 'N'.
	char range_gesvdx; // Must be 'A', 'V', 'I'. 
	int il, iu; // the indices of the smallest and largest singular values.
	float vl, vu; //  the lower and upper bounds of the interval.

	/* Parameters for 'gesvdq' API  */
        char joba_gesvdq; //  Must be 'A', 'H', 'M' , 'E'
        char jobu_gesvdq; // Must be 'A', 'S', 'R' , 'N' 
        char jobv_gesvdq; // Must be 'A', 'V', 'R' , 'N'.

	/* Thresholds for the APIs  */
	float svd_threshold; // threshold for the gghrd API

}SVD_paramlist;

extern SVD_paramlist svd_paramslist[NUM_SUB_TESTS];

/*Function to read EIG parametes from config file */
void Read_EIG_params( const char *file_name );

/*Function to read EIG non-symmetric parametes from config file */
void Read_EIG_non_sym_params( const char *file_name );

/*Function to read SVD parametes from config file */
void Read_SVD_param ( const char *file_name );

// --- nrm2 ---

void bl1_snrm2( int n, float*    x, int incx, float*  norm );
void bl1_dnrm2( int n, double*   x, int incx, double* norm );
//void bl1_cnrm2( int n, scomplex* x, int incx, float*  norm );
//void bl1_znrm2( int n, dcomplex* x, int incx, double* norm );

typedef int FLA_Error;
//FLA_Error FLA_Nrm2_external( FLA_Obj x, FLA_Obj nrm_x );

void Read_config_parameters ();
#endif // FLAME_GBENCH_MAIN_H
