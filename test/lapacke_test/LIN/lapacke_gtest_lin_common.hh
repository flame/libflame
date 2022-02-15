#ifndef LAPACKE_GTEST_LIN_H
#define LAPACKE_GTEST_LIN_H

//#include "lapacke_utils.h"
#include "lapacke.h"

class lin_solver_double_parameters{
	
   public:
   /* input params to the API **/
   int 	matrix_layout;
   char uplo;
   int 	n;
   int 	nrhs;
   double  *a, *aref;
   int 	lda;
   int     *ipiv, *ipivref;
   double  *b, *bref;
   int 	ldb;
   
   int info, inforef;
   
   /* Below variables are specific to 'gbcon' API */
   char norm;
   int kl;
   int ku;
   double anorm;
   double* ab;
   int ldab;
   double rcond; 
   
   public: 
      lin_solver_double_parameters (int matrix_layout, char uplo, int n,
                            int nrhs, int lda, int ldb);
      ~lin_solver_double_parameters ();

}; /* end of lin_solver_double_parameters  class definition */


class lin_solver_float_parameters{
	
   public:
   /* input params to the API **/
   int 	matrix_layout;
   char uplo;
   int 	n;
   int 	nrhs;
   float  *a, *aref;
   int 	lda;
   int     *ipiv, *ipivref;
   float  *b, *bref;
   int 	ldb;
   
   int info, inforef;

   public: 
      lin_solver_float_parameters (int matrix_layout, char uplo, int n,
                            int nrhs, int lda, int ldb);
      ~lin_solver_float_parameters ();

}; /* end of lin_solver_float_parameters  class definition */


/* Begin single precision complex_common_parameters  class definition */
class lin_solver_scomplex_parameters{
   public:
   /* input params to the API **/
   int 	matrix_layout;
   char uplo;
   int 	n;
   int 	nrhs;
   int 	lda;
   int  *ipiv, *ipivref;
   int 	ldb;
   
   int info, inforef;
   /* Local arrays */
   lapack_complex_float *a, *b, *aref, *bref;

   public: 
      lin_solver_scomplex_parameters (int matrix_layout, char uplo, int n,
                            int nrhs, int lda, int ldb);
      ~lin_solver_scomplex_parameters ();   
   
	
};  /* end of scomplex_common_parameters  class definition */


/* Begin single precision complex_common_parameters  class definition */
class lin_solver_dcomplex_parameters{
   public:
   /* input params to the API **/
   int 	matrix_layout;
   char uplo;
   int 	n;
   int 	nrhs;
   int 	lda;
   int  *ipiv, *ipivref;
   int 	ldb;
   
   int info, inforef;
   /* Local arrays */
   lapack_complex_double *a, *b, *aref, *bref;

   public: 
      lin_solver_dcomplex_parameters (int matrix_layout, char uplo, int n,
                            int nrhs, int lda, int ldb);
      ~lin_solver_dcomplex_parameters ();   
   
	
};  /* end of scomplex_common_parameters  class definition */


#endif