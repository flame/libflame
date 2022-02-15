#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"

#define gbrfs_free() \
  if (a != NULL)    free (a  ); \
  if (aref != NULL) free (aref); \
  if (b != NULL)    free (b   ); \
  if (bref != NULL) free (bref); \
  if (x != NULL)    free (x  ); \
  if (xref != NULL) free (xref); \
  if (af != NULL)    free (af  ); \
  if (afref != NULL) free (afref); \
  if (ferr != NULL)    free (ferr  ); \
  if (ferrref != NULL) free (ferrref); \
  if (ipiv != NULL)    free (ipiv  ); \
  if (ipivref != NULL) free (ipivref); \
  if (berr != NULL)    free (berr  ); \
  if (berrref != NULL) free (berrref); \
  if( hModule != NULL) dlclose(hModule); \
  if(dModule != NULL) dlclose(dModule)


// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin gbrfs_float_parameters  class definition */
class gbrfs_float_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      float diff; // capture difference between ref o/p & libflame lapacke o/p.
      float diff_berr, diff_ferr, diff_xerr;
      void *hModule, *dModule;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char trans; //  Must be 'N' , 'T' or 'C'.

      lapack_int n; // The order of A (Square matrix).
      lapack_int nrhs; // The number of right-hand sides
      lapack_int kl;// The number of subdiagonals within the band of A
      lapack_int ku; // The number of superdiagonals within the band of A
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_int ldaf;  //  leading dimension of 'af'
      lapack_int ldx; //  leading dimension of 'x'

      /* Input/ Output parameters */
      float *b, *bref; //right-hand sides for the systems of equations.
      float *a, *aref; //The array ab contains the matrix A
      float *af, *afref; //contains the factored form of the matrix A
      
      /* Output parameters */
      lapack_int *ipiv, *ipivref; // pivot buffer
      float* x, *xref;
      float* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      float* berr, *berrref; // component-wise relative backward error for each solution vector xj.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      gbrfs_float_parameters ( int matrix_layout_i, char trans_i,
                                 lapack_int kl_i, lapack_int ku_i,
                               lapack_int n_i, lapack_int nrhs_i);
              
      ~gbrfs_float_parameters (); 
};  /* end of gbrfs_float_parameters  class definition */


/* Constructor gbrfs_float_parameters definition */
gbrfs_float_parameters:: gbrfs_float_parameters ( int matrix_layout_i, 
                       char trans_i, lapack_int kl_i, lapack_int ku_i,
						lapack_int n_i, lapack_int nrhs_i) {
                                    
    matrix_layout = matrix_layout_i;
    trans = trans_i;
    kl = kl_i;
    ku = ku_i;
    // equed acts as input when " fact = 'F' " else output
    n = n_i;
    nrhs = nrhs_i;
    lda = n;
    ldaf = n;
    
#if LAPACKE_TEST_VERBOSE
   printf(" \n gbrfs float:  n: %d,  trans: %c  lda: %d  \
ldb: %d nrhs: %d  ldaf: %d  ldx: %d \n",  n, trans, lda, 
                                          ldb, nrhs);
#endif
    b_bufsize = n*nrhs;
    x_bufsize = b_bufsize;

    if(matrix_layout==LAPACK_COL_MAJOR){
        ldb = n;
        ldx = n;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        ldb = nrhs;
        ldx = nrhs;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_float_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_float_buffer_pair( &x, &xref, x_bufsize);
    lapacke_gtest_alloc_float_buffer_pair( &af, &afref, (n*n));
    lapacke_gtest_alloc_float_buffer_pair( &ferr, &ferrref, nrhs);
    lapacke_gtest_alloc_float_buffer_pair( &berr, &berrref, nrhs);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);
    

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (x==NULL) || (xref==NULL) || \
        (af==NULL) || (afref==NULL) || \
        (ferr==NULL) || (ferrref==NULL) || \
        (ipiv==NULL) || (ipivref==NULL) || \
        (berr==NULL) || (berrref==NULL) ){
       EXPECT_FALSE( true) << "pbtrs_float_parameters object: malloc error.";
       gbrfs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, n*n);
    memcpy(af, a, (n*n*sizeof(float)));
    memcpy(afref, a, (n*n*sizeof(float)));
    
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, b_bufsize);
     memcpy(x, b, ( b_bufsize*sizeof(float)));
    memcpy(xref, b, ( b_bufsize*sizeof(float)));
   
    lapacke_gtest_init_float_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);
    lapacke_gtest_init_int_buffer_pair_with_constant(ipiv, ipivref, nrhs, 0);

    
   } /* end of Constructor  */

gbrfs_float_parameters:: ~gbrfs_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gbrfs_float_parameters object: destructor invoked. \n");
#endif
   gbrfs_free();
}

//  Test fixture class definition
class sgbrfs_test  : public  ::testing::Test {
public:
   gbrfs_float_parameters  *sgbrfs_obj;
   void SetUp();  
   void TearDown () { delete sgbrfs_obj; }
};


void sgbrfs_test::SetUp(){

    /* LAPACKE SGBRFS prototype */
    typedef int (*Fptr_NL_LAPACKE_sgbrfs) (int matrix_layout, char trans,
		lapack_int n, lapack_int kl, lapack_int ku, lapack_int nrhs,
		const float* ab, lapack_int ldab, const float* afb,
		lapack_int ldafb, const lapack_int* ipiv, const float* b,
		lapack_int ldb, float* x, lapack_int ldx, float* ferr, float* berr);

    Fptr_NL_LAPACKE_sgbrfs SGBRFS;

     /* LAPACKE SGBTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_sgbtrf) ( int matrix_layout,lapack_int m,
                                lapack_int n, lapack_int kl, lapack_int ku,
                              float *ab,lapack_int ldab,lapack_int* ipiv );

    Fptr_NL_LAPACKE_sgbtrf SGBTRF;

    typedef int (*Fptr_NL_LAPACKE_sgbtrs) (int matrix_layout , char trans,
		lapack_int n , lapack_int kl , lapack_int ku , lapack_int nrhs ,
		const float * ab , lapack_int ldab , const lapack_int * ipiv ,
		float * b , lapack_int ldb );

    Fptr_NL_LAPACKE_sgbtrs SGBTRS;

    sgbrfs_obj = new gbrfs_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].kl,
                           lin_solver_paramslist[idx].ku,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    sgbrfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgbrfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgbrfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgbrfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SGBTRF = (Fptr_NL_LAPACKE_sgbtrf)dlsym(sgbrfs_obj->hModule,"LAPACKE_sgbtrf");
    ASSERT_TRUE(SGBTRF != NULL) << "failed to get the Netlib LAPACKE_sgbtrf symbol";

    SGBTRS = (Fptr_NL_LAPACKE_sgbtrs)dlsym(sgbrfs_obj->hModule,"LAPACKE_sgbtrs");
    ASSERT_TRUE(SGBTRS != NULL) << "failed to get the Netlib LAPACKE_sgbtrs symbol";

    SGBRFS = (Fptr_NL_LAPACKE_sgbrfs)dlsym(sgbrfs_obj->hModule, "LAPACKE_sgbrfs");
    ASSERT_TRUE(SGBRFS != NULL) << "failed to get the Netlib LAPACKE_sgbrfs symbol";
    
	sgbrfs_obj->inforef = SGBTRF( sgbrfs_obj->matrix_layout,
								  sgbrfs_obj->n, sgbrfs_obj->n,
								  sgbrfs_obj->kl, sgbrfs_obj->ku,
								  sgbrfs_obj->afref,
								  sgbrfs_obj->lda,
								  sgbrfs_obj->ipivref);
						   
	sgbrfs_obj->info = LAPACKE_sgbtrf( sgbrfs_obj->matrix_layout,
									sgbrfs_obj->n, sgbrfs_obj->n,
								  sgbrfs_obj->kl, sgbrfs_obj->ku,
									   sgbrfs_obj->af,
									   sgbrfs_obj->lda,
									   sgbrfs_obj->ipiv);

	sgbrfs_obj->inforef = SGBTRS( sgbrfs_obj->matrix_layout,
								  sgbrfs_obj->trans,
								  sgbrfs_obj->n,
								  sgbrfs_obj->kl, sgbrfs_obj->ku,
								  sgbrfs_obj->nrhs,
								  sgbrfs_obj->afref,
								  sgbrfs_obj->lda,
								  sgbrfs_obj->ipivref,
                                  sgbrfs_obj->xref, sgbrfs_obj->ldx
								  );

	sgbrfs_obj->inforef = SGBTRS( sgbrfs_obj->matrix_layout,
								  sgbrfs_obj->trans,
								  sgbrfs_obj->n,
								  sgbrfs_obj->kl, sgbrfs_obj->ku,
								  sgbrfs_obj->nrhs,
								  sgbrfs_obj->af,
								  sgbrfs_obj->lda,
								  sgbrfs_obj->ipiv,
                                  sgbrfs_obj->x, sgbrfs_obj->ldx
								  );
								  
    /* Compute the reference o/p by invoking Netlib-Lapack's API */
    sgbrfs_obj->inforef = SGBRFS( sgbrfs_obj->matrix_layout,
                                  sgbrfs_obj->trans, sgbrfs_obj->n,
                                  sgbrfs_obj->kl, sgbrfs_obj->ku,
                                  sgbrfs_obj->nrhs,
                                  sgbrfs_obj->aref, sgbrfs_obj->lda, 
                                  sgbrfs_obj->afref, sgbrfs_obj->ldaf,
                                  sgbrfs_obj->ipivref,
                                  sgbrfs_obj->bref, sgbrfs_obj->ldb,
                                  sgbrfs_obj->xref, sgbrfs_obj->ldx,
                                  sgbrfs_obj->ferrref,
                                  sgbrfs_obj->berrref
								  );
								  
    /* Compute libflame's Lapacke o/p  */								  
    sgbrfs_obj->info = LAPACKE_sgbrfs( sgbrfs_obj->matrix_layout,
                                  sgbrfs_obj->trans, sgbrfs_obj->n,
                                  sgbrfs_obj->kl, sgbrfs_obj->ku,
                                  sgbrfs_obj->nrhs,
                                  sgbrfs_obj->a, sgbrfs_obj->lda, 
                                  sgbrfs_obj->af, sgbrfs_obj->ldaf,
                                  sgbrfs_obj->ipiv,
                                  sgbrfs_obj->b, sgbrfs_obj->ldb,
                                  sgbrfs_obj->x, sgbrfs_obj->ldx,
                                  sgbrfs_obj->ferr,
                                  sgbrfs_obj->berr
								  );


    if( sgbrfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_sgbrfs is wrong\n", sgbrfs_obj->info );
    }
    if( sgbrfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sgbrfs is wrong\n", 
        sgbrfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    sgbrfs_obj->diff_xerr =  computeDiff_s( sgbrfs_obj->x_bufsize, 
                sgbrfs_obj->x, sgbrfs_obj->xref );

    sgbrfs_obj->diff_berr =  computeDiff_s( sgbrfs_obj->nrhs, 
                sgbrfs_obj->berr, sgbrfs_obj->berrref );
                
    sgbrfs_obj->diff_ferr =  computeDiff_s( sgbrfs_obj->nrhs, 
                sgbrfs_obj->ferr, sgbrfs_obj->ferrref );
}

TEST_F(sgbrfs_test, sgbrfs1) {
    EXPECT_NEAR(0.0, sgbrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, sgbrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, sgbrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(sgbrfs_test, sgbrfs2) {
    EXPECT_NEAR(0.0, sgbrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, sgbrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, sgbrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(sgbrfs_test, sgbrfs3) {
    EXPECT_NEAR(0.0, sgbrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, sgbrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, sgbrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(sgbrfs_test, sgbrfs4) {
    EXPECT_NEAR(0.0, sgbrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, sgbrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, sgbrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

/* Begin gbrfs_double_parameters  class definition */
class gbrfs_double_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      double diff; // capture difference between ref o/p & libflame lapacke o/p.
      double diff_berr, diff_ferr, diff_xerr;
      void *hModule, *dModule;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char trans; //  Must be 'N' , 'T' or 'C'.

      lapack_int n; // The order of A (Square matrix).
      lapack_int nrhs; // The number of right-hand sides
      lapack_int kl;// The number of subdiagonals within the band of A
      lapack_int ku; // The number of superdiagonals within the band of A
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_int ldaf;  //  leading dimension of 'af'
      lapack_int ldx; //  leading dimension of 'x'

      /* Input/ Output parameters */
      double *b, *bref; //right-hand sides for the systems of equations.
      double *a, *aref; //The array ab contains the matrix A
      double *af, *afref; //contains the factored form of the matrix A
      
      /* Output parameters */
      lapack_int *ipiv, *ipivref; // pivot buffer
      double* x, *xref;
      double* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      double* berr, *berrref; // component-wise relative backward error for each solution vector xj.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      gbrfs_double_parameters ( int matrix_layout_i, char trans_i,
                                 lapack_int kl_i, lapack_int ku_i,
                               lapack_int n_i, lapack_int nrhs_i);
              
      ~gbrfs_double_parameters (); 
};  /* end of gbrfs_double_parameters  class definition */


/* Constructor gbrfs_double_parameters definition */
gbrfs_double_parameters:: gbrfs_double_parameters ( int matrix_layout_i, 
                       char trans_i, lapack_int kl_i, lapack_int ku_i,
						lapack_int n_i, lapack_int nrhs_i) {
                                    
    matrix_layout = matrix_layout_i;
    trans = trans_i;
    kl = kl_i;
    ku = ku_i;
    // equed acts as input when " fact = 'F' " else output
    n = n_i;
    nrhs = nrhs_i;
    lda = n;
    ldaf = n;
    
#if LAPACKE_TEST_VERBOSE
   printf(" \n gbrfs double:  n: %d,  trans: %c  lda: %d  \
ldb: %d nrhs: %d  ldaf: %d  ldx: %d \n",  n, trans, lda, 
                                          ldb, nrhs);
#endif
    b_bufsize = n*nrhs;
    x_bufsize = b_bufsize;

    if(matrix_layout==LAPACK_COL_MAJOR){
        ldb = n;
        ldx = n;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        ldb = nrhs;
        ldx = nrhs;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_double_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_double_buffer_pair( &x, &xref, x_bufsize);
    lapacke_gtest_alloc_double_buffer_pair( &af, &afref, (n*n));
    lapacke_gtest_alloc_double_buffer_pair( &ferr, &ferrref, nrhs);
    lapacke_gtest_alloc_double_buffer_pair( &berr, &berrref, nrhs);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);
    

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (x==NULL) || (xref==NULL) || \
        (af==NULL) || (afref==NULL) || \
        (ferr==NULL) || (ferrref==NULL) || \
        (ipiv==NULL) || (ipivref==NULL) || \
        (berr==NULL) || (berrref==NULL) ){
       EXPECT_FALSE( true) << "pbtrs_double_parameters object: malloc error.";
       gbrfs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, n*n);
    memcpy(af, a, (n*n*sizeof(double)));
    memcpy(afref, a, (n*n*sizeof(double)));
    
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, b_bufsize);
     memcpy(x, b, ( b_bufsize*sizeof(double)));
    memcpy(xref, b, ( b_bufsize*sizeof(double)));
   
    lapacke_gtest_init_double_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);
    lapacke_gtest_init_int_buffer_pair_with_constant(ipiv, ipivref, nrhs, 0);

    
   } /* end of Constructor  */

gbrfs_double_parameters:: ~gbrfs_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gbrfs_double_parameters object: destructor invoked. \n");
#endif
   gbrfs_free();
}

//  Test fixture class definition
class dgbrfs_test  : public  ::testing::Test {
public:
   gbrfs_double_parameters  *dgbrfs_obj;
   void SetUp();  
   void TearDown () { delete dgbrfs_obj; }
};


void dgbrfs_test::SetUp(){

    /* LAPACKE DGBRFS prototype */
    typedef int (*Fptr_NL_LAPACKE_dgbrfs) (int matrix_layout, char trans,
		lapack_int n, lapack_int kl, lapack_int ku, lapack_int nrhs,
		const double* ab, lapack_int ldab, const double* afb,
		lapack_int ldafb, const lapack_int* ipiv, const double* b,
		lapack_int ldb, double* x, lapack_int ldx, double* ferr, double* berr);

    Fptr_NL_LAPACKE_dgbrfs DGBRFS;

     /* LAPACKE DGBTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_dgbtrf) ( int matrix_layout,lapack_int m,
                                lapack_int n, lapack_int kl, lapack_int ku,
                              double *ab,lapack_int ldab,lapack_int* ipiv );

    Fptr_NL_LAPACKE_dgbtrf DGBTRF;

    typedef int (*Fptr_NL_LAPACKE_dgbtrs) (int matrix_layout , char trans,
		lapack_int n , lapack_int kl , lapack_int ku , lapack_int nrhs ,
		const double * ab , lapack_int ldab , const lapack_int * ipiv ,
		double * b , lapack_int ldb );

    Fptr_NL_LAPACKE_dgbtrs DGBTRS;

    dgbrfs_obj = new gbrfs_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].kl,
                           lin_solver_paramslist[idx].ku,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    dgbrfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgbrfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgbrfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgbrfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DGBTRF = (Fptr_NL_LAPACKE_dgbtrf)dlsym(dgbrfs_obj->hModule,"LAPACKE_dgbtrf");
    ASSERT_TRUE(DGBTRF != NULL) << "failed to get the Netlib LAPACKE_dgbtrf symbol";

    DGBTRS = (Fptr_NL_LAPACKE_dgbtrs)dlsym(dgbrfs_obj->hModule,"LAPACKE_dgbtrs");
    ASSERT_TRUE(DGBTRS != NULL) << "failed to get the Netlib LAPACKE_dgbtrs symbol";

    DGBRFS = (Fptr_NL_LAPACKE_dgbrfs)dlsym(dgbrfs_obj->hModule, "LAPACKE_dgbrfs");
    ASSERT_TRUE(DGBRFS != NULL) << "failed to get the Netlib LAPACKE_dgbrfs symbol";
    
	dgbrfs_obj->inforef = DGBTRF( dgbrfs_obj->matrix_layout,
								  dgbrfs_obj->n, dgbrfs_obj->n,
								  dgbrfs_obj->kl, dgbrfs_obj->ku,
								  dgbrfs_obj->afref,
								  dgbrfs_obj->lda,
								  dgbrfs_obj->ipivref);
						   
	dgbrfs_obj->info = LAPACKE_dgbtrf( dgbrfs_obj->matrix_layout,
									dgbrfs_obj->n, dgbrfs_obj->n,
								  dgbrfs_obj->kl, dgbrfs_obj->ku,
									   dgbrfs_obj->af,
									   dgbrfs_obj->lda,
									   dgbrfs_obj->ipiv);

	dgbrfs_obj->inforef = DGBTRS( dgbrfs_obj->matrix_layout,
								  dgbrfs_obj->trans,
								  dgbrfs_obj->n,
								  dgbrfs_obj->kl, dgbrfs_obj->ku,
								  dgbrfs_obj->nrhs,
								  dgbrfs_obj->afref,
								  dgbrfs_obj->lda,
								  dgbrfs_obj->ipivref,
                                  dgbrfs_obj->xref, dgbrfs_obj->ldx
								  );

	dgbrfs_obj->inforef = DGBTRS( dgbrfs_obj->matrix_layout,
								  dgbrfs_obj->trans,
								  dgbrfs_obj->n,
								  dgbrfs_obj->kl, dgbrfs_obj->ku,
								  dgbrfs_obj->nrhs,
								  dgbrfs_obj->af,
								  dgbrfs_obj->lda,
								  dgbrfs_obj->ipiv,
                                  dgbrfs_obj->x, dgbrfs_obj->ldx
								  );
								  
    /* Compute the reference o/p by invoking Netlib-Lapack's API */
    dgbrfs_obj->inforef = DGBRFS( dgbrfs_obj->matrix_layout,
                                  dgbrfs_obj->trans, dgbrfs_obj->n,
                                  dgbrfs_obj->kl, dgbrfs_obj->ku,
                                  dgbrfs_obj->nrhs,
                                  dgbrfs_obj->aref, dgbrfs_obj->lda, 
                                  dgbrfs_obj->afref, dgbrfs_obj->ldaf,
                                  dgbrfs_obj->ipivref,
                                  dgbrfs_obj->bref, dgbrfs_obj->ldb,
                                  dgbrfs_obj->xref, dgbrfs_obj->ldx,
                                  dgbrfs_obj->ferrref,
                                  dgbrfs_obj->berrref
								  );
								  
    /* Compute libflame's Lapacke o/p  */								  
    dgbrfs_obj->info = LAPACKE_dgbrfs( dgbrfs_obj->matrix_layout,
                                  dgbrfs_obj->trans, dgbrfs_obj->n,
                                  dgbrfs_obj->kl, dgbrfs_obj->ku,
                                  dgbrfs_obj->nrhs,
                                  dgbrfs_obj->a, dgbrfs_obj->lda, 
                                  dgbrfs_obj->af, dgbrfs_obj->ldaf,
                                  dgbrfs_obj->ipiv,
                                  dgbrfs_obj->b, dgbrfs_obj->ldb,
                                  dgbrfs_obj->x, dgbrfs_obj->ldx,
                                  dgbrfs_obj->ferr,
                                  dgbrfs_obj->berr
								  );


    if( dgbrfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dgbrfs is wrong\n", dgbrfs_obj->info );
    }
    if( dgbrfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dgbrfs is wrong\n", 
        dgbrfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dgbrfs_obj->diff_xerr =  computeDiff_d( dgbrfs_obj->x_bufsize, 
                dgbrfs_obj->x, dgbrfs_obj->xref );

    dgbrfs_obj->diff_berr =  computeDiff_d( dgbrfs_obj->nrhs, 
                dgbrfs_obj->berr, dgbrfs_obj->berrref );
                
    dgbrfs_obj->diff_ferr =  computeDiff_d( dgbrfs_obj->nrhs, 
                dgbrfs_obj->ferr, dgbrfs_obj->ferrref );
}

TEST_F(dgbrfs_test, dgbrfs1) {
    EXPECT_NEAR(0.0, dgbrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dgbrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dgbrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(dgbrfs_test, dgbrfs2) {
    EXPECT_NEAR(0.0, dgbrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dgbrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dgbrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(dgbrfs_test, dgbrfs3) {
    EXPECT_NEAR(0.0, dgbrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dgbrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dgbrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(dgbrfs_test, dgbrfs4) {
    EXPECT_NEAR(0.0, dgbrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dgbrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dgbrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

/* Begin gbrfs_scomplex_parameters  class definition */
class gbrfs_scomplex_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      float diff_berr, diff_ferr, diff_xerr;
      void *hModule, *dModule;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char trans; //  Must be 'N' , 'T' or 'C'.

      lapack_int n; // The order of A (Square matrix).
      lapack_int nrhs; // The number of right-hand sides
      lapack_int kl;// The number of subdiagonals within the band of A
      lapack_int ku; // The number of superdiagonals within the band of A
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_int ldaf;  //  leading dimension of 'af'
      lapack_int ldx; //  leading dimension of 'x'

      /* Input/ Output parameters */
      lapack_complex_float *b, *bref; //right-hand sides for the systems of equations.
      lapack_complex_float *a, *aref; //The array ab contains the matrix A
      lapack_complex_float *af, *afref; //contains the factored form of the matrix A
      
      /* Output parameters */
      lapack_int *ipiv, *ipivref; // pivot buffer
      lapack_complex_float* x, *xref;
      float* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      float* berr, *berrref; // component-wise relative backward error for each solution vector xj.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      gbrfs_scomplex_parameters ( int matrix_layout_i, char trans_i,
                                 lapack_int kl_i, lapack_int ku_i,
                               lapack_int n_i, lapack_int nrhs_i);
              
      ~gbrfs_scomplex_parameters (); 
};  /* end of gbrfs_scomplex_parameters  class definition */


/* Constructor gbrfs_scomplex_parameters definition */
gbrfs_scomplex_parameters:: gbrfs_scomplex_parameters ( int matrix_layout_i, 
                       char trans_i, lapack_int kl_i, lapack_int ku_i,
						lapack_int n_i, lapack_int nrhs_i) {
                                    
    matrix_layout = matrix_layout_i;
    trans = trans_i;
    kl = kl_i;
    ku = ku_i;
    // equed acts as input when " fact = 'F' " else output
    n = n_i;
    nrhs = nrhs_i;
    lda = n;
    ldaf = n;
    
#if LAPACKE_TEST_VERBOSE
   printf(" \n gbrfs lapack_complex_float:  n: %d,  trans: %c  lda: %d  \
ldb: %d nrhs: %d  ldaf: %d  ldx: %d \n",  n, trans, lda, 
                                          ldb, nrhs);
#endif
    b_bufsize = n*nrhs;
    x_bufsize = b_bufsize;

    if(matrix_layout==LAPACK_COL_MAJOR){
        ldb = n;
        ldx = n;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        ldb = nrhs;
        ldx = nrhs;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &x, &xref, x_bufsize);
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &af, &afref, (n*n));
    lapacke_gtest_alloc_float_buffer_pair( &ferr, &ferrref, nrhs);
    lapacke_gtest_alloc_float_buffer_pair( &berr, &berrref, nrhs);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);
    

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (x==NULL) || (xref==NULL) || \
        (af==NULL) || (afref==NULL) || \
        (ferr==NULL) || (ferrref==NULL) || \
        (ipiv==NULL) || (ipivref==NULL) || \
        (berr==NULL) || (berrref==NULL) ){
       EXPECT_FALSE( true) << "pbtrs_scomplex_parameters object: malloc error.";
       gbrfs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, n*n);
    memcpy(af, a, (n*n*sizeof(lapack_complex_float)));
    memcpy(afref, a, (n*n*sizeof(lapack_complex_float)));
    
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
     memcpy(x, b, ( b_bufsize*sizeof(lapack_complex_float)));
    memcpy(xref, b, ( b_bufsize*sizeof(lapack_complex_float)));
   
    lapacke_gtest_init_float_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);
    lapacke_gtest_init_int_buffer_pair_with_constant(ipiv, ipivref, nrhs, 0);

    
   } /* end of Constructor  */

gbrfs_scomplex_parameters:: ~gbrfs_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gbrfs_scomplex_parameters object: destructor invoked. \n");
#endif
   gbrfs_free();
}

//  Test fixture class definition
class cgbrfs_test  : public  ::testing::Test {
public:
   gbrfs_scomplex_parameters  *cgbrfs_obj;
   void SetUp();  
   void TearDown () { delete cgbrfs_obj; }
};


void cgbrfs_test::SetUp(){

    /* LAPACKE CGBRFS prototype */
    typedef int (*Fptr_NL_LAPACKE_cgbrfs) (int matrix_layout, char trans,
		lapack_int n, lapack_int kl, lapack_int ku, lapack_int nrhs,
		const lapack_complex_float* ab, lapack_int ldab, const lapack_complex_float* afb,
		lapack_int ldafb, const lapack_int* ipiv, const lapack_complex_float* b,
		lapack_int ldb, lapack_complex_float* x, lapack_int ldx, float* ferr, float* berr);

    Fptr_NL_LAPACKE_cgbrfs CGBRFS;

     /* LAPACKE CGBTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_cgbtrf) ( int matrix_layout,lapack_int m,
                                lapack_int n, lapack_int kl, lapack_int ku,
                              lapack_complex_float *ab,lapack_int ldab,lapack_int* ipiv );

    Fptr_NL_LAPACKE_cgbtrf CGBTRF;

    typedef int (*Fptr_NL_LAPACKE_cgbtrs) (int matrix_layout , char trans,
		lapack_int n , lapack_int kl , lapack_int ku , lapack_int nrhs ,
		const lapack_complex_float * ab , lapack_int ldab , const lapack_int * ipiv ,
		lapack_complex_float * b , lapack_int ldb );

    Fptr_NL_LAPACKE_cgbtrs CGBTRS;

    cgbrfs_obj = new gbrfs_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].kl,
                           lin_solver_paramslist[idx].ku,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    cgbrfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgbrfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgbrfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgbrfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CGBTRF = (Fptr_NL_LAPACKE_cgbtrf)dlsym(cgbrfs_obj->hModule,"LAPACKE_cgbtrf");
    ASSERT_TRUE(CGBTRF != NULL) << "failed to get the Netlib LAPACKE_cgbtrf symbol";

    CGBTRS = (Fptr_NL_LAPACKE_cgbtrs)dlsym(cgbrfs_obj->hModule,"LAPACKE_cgbtrs");
    ASSERT_TRUE(CGBTRS != NULL) << "failed to get the Netlib LAPACKE_cgbtrs symbol";

    CGBRFS = (Fptr_NL_LAPACKE_cgbrfs)dlsym(cgbrfs_obj->hModule, "LAPACKE_cgbrfs");
    ASSERT_TRUE(CGBRFS != NULL) << "failed to get the Netlib LAPACKE_cgbrfs symbol";
    
	cgbrfs_obj->inforef = CGBTRF( cgbrfs_obj->matrix_layout,
								  cgbrfs_obj->n, cgbrfs_obj->n,
								  cgbrfs_obj->kl, cgbrfs_obj->ku,
								  cgbrfs_obj->afref,
								  cgbrfs_obj->lda,
								  cgbrfs_obj->ipivref);
						   
	cgbrfs_obj->info = LAPACKE_cgbtrf( cgbrfs_obj->matrix_layout,
									cgbrfs_obj->n, cgbrfs_obj->n,
								  cgbrfs_obj->kl, cgbrfs_obj->ku,
									   cgbrfs_obj->af,
									   cgbrfs_obj->lda,
									   cgbrfs_obj->ipiv);

	cgbrfs_obj->inforef = CGBTRS( cgbrfs_obj->matrix_layout,
								  cgbrfs_obj->trans,
								  cgbrfs_obj->n,
								  cgbrfs_obj->kl, cgbrfs_obj->ku,
								  cgbrfs_obj->nrhs,
								  cgbrfs_obj->afref,
								  cgbrfs_obj->lda,
								  cgbrfs_obj->ipivref,
                                  cgbrfs_obj->xref, cgbrfs_obj->ldx
								  );

	cgbrfs_obj->inforef = CGBTRS( cgbrfs_obj->matrix_layout,
								  cgbrfs_obj->trans,
								  cgbrfs_obj->n,
								  cgbrfs_obj->kl, cgbrfs_obj->ku,
								  cgbrfs_obj->nrhs,
								  cgbrfs_obj->af,
								  cgbrfs_obj->lda,
								  cgbrfs_obj->ipiv,
                                  cgbrfs_obj->x, cgbrfs_obj->ldx
								  );
								  
    /* Compute the reference o/p by invoking Netlib-Lapack's API */
    cgbrfs_obj->inforef = CGBRFS( cgbrfs_obj->matrix_layout,
                                  cgbrfs_obj->trans, cgbrfs_obj->n,
                                  cgbrfs_obj->kl, cgbrfs_obj->ku,
                                  cgbrfs_obj->nrhs,
                                  cgbrfs_obj->aref, cgbrfs_obj->lda, 
                                  cgbrfs_obj->afref, cgbrfs_obj->ldaf,
                                  cgbrfs_obj->ipivref,
                                  cgbrfs_obj->bref, cgbrfs_obj->ldb,
                                  cgbrfs_obj->xref, cgbrfs_obj->ldx,
                                  cgbrfs_obj->ferrref,
                                  cgbrfs_obj->berrref
								  );
								  
    /* Compute libflame's Lapacke o/p  */								  
    cgbrfs_obj->info = LAPACKE_cgbrfs( cgbrfs_obj->matrix_layout,
                                  cgbrfs_obj->trans, cgbrfs_obj->n,
                                  cgbrfs_obj->kl, cgbrfs_obj->ku,
                                  cgbrfs_obj->nrhs,
                                  cgbrfs_obj->a, cgbrfs_obj->lda, 
                                  cgbrfs_obj->af, cgbrfs_obj->ldaf,
                                  cgbrfs_obj->ipiv,
                                  cgbrfs_obj->b, cgbrfs_obj->ldb,
                                  cgbrfs_obj->x, cgbrfs_obj->ldx,
                                  cgbrfs_obj->ferr,
                                  cgbrfs_obj->berr
								  );


    if( cgbrfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cgbrfs is wrong\n", cgbrfs_obj->info );
    }
    if( cgbrfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cgbrfs is wrong\n", 
        cgbrfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cgbrfs_obj->diff_xerr =  computeDiff_c( cgbrfs_obj->x_bufsize, 
                cgbrfs_obj->x, cgbrfs_obj->xref );

    cgbrfs_obj->diff_berr =  computeDiff_s( cgbrfs_obj->nrhs, 
                cgbrfs_obj->berr, cgbrfs_obj->berrref );
                
    cgbrfs_obj->diff_ferr =  computeDiff_s( cgbrfs_obj->nrhs, 
                cgbrfs_obj->ferr, cgbrfs_obj->ferrref );
}

TEST_F(cgbrfs_test, cgbrfs1) {
    EXPECT_NEAR(0.0, cgbrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, cgbrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, cgbrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(cgbrfs_test, cgbrfs2) {
    EXPECT_NEAR(0.0, cgbrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, cgbrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, cgbrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(cgbrfs_test, cgbrfs3) {
    EXPECT_NEAR(0.0, cgbrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, cgbrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, cgbrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(cgbrfs_test, cgbrfs4) {
    EXPECT_NEAR(0.0, cgbrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, cgbrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, cgbrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

/* Begin gbrfs_dcomplex_parameters  class definition */
class gbrfs_dcomplex_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      double diff_berr, diff_ferr, diff_xerr;
      void *hModule, *dModule;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char trans; //  Must be 'N' , 'T' or 'C'.

      lapack_int n; // The order of A (Square matrix).
      lapack_int nrhs; // The number of right-hand sides
      lapack_int kl;// The number of subdiagonals within the band of A
      lapack_int ku; // The number of superdiagonals within the band of A
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_int ldaf;  //  leading dimension of 'af'
      lapack_int ldx; //  leading dimension of 'x'

      /* Input/ Output parameters */
      lapack_complex_double *b, *bref; //right-hand sides for the systems of equations.
      lapack_complex_double *a, *aref; //The array ab contains the matrix A
      lapack_complex_double *af, *afref; //contains the factored form of the matrix A
      
      /* Output parameters */
      lapack_int *ipiv, *ipivref; // pivot buffer
      lapack_complex_double* x, *xref;
      double* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      double* berr, *berrref; // component-wise relative backward error for each solution vector xj.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      gbrfs_dcomplex_parameters ( int matrix_layout_i, char trans_i,
                                 lapack_int kl_i, lapack_int ku_i,
                               lapack_int n_i, lapack_int nrhs_i);
              
      ~gbrfs_dcomplex_parameters (); 
};  /* end of gbrfs_dcomplex_parameters  class definition */


/* Constructor gbrfs_dcomplex_parameters definition */
gbrfs_dcomplex_parameters:: gbrfs_dcomplex_parameters ( int matrix_layout_i, 
                       char trans_i, lapack_int kl_i, lapack_int ku_i,
						lapack_int n_i, lapack_int nrhs_i) {
                                    
    matrix_layout = matrix_layout_i;
    trans = trans_i;
    kl = kl_i;
    ku = ku_i;
    // equed acts as input when " fact = 'F' " else output
    n = n_i;
    nrhs = nrhs_i;
    lda = n;
    ldaf = n;
    
#if LAPACKE_TEST_VERBOSE
   printf(" \n gbrfs lapack_complex_double:  n: %d,  trans: %c  lda: %d  \
ldb: %d nrhs: %d  ldaf: %d  ldx: %d \n",  n, trans, lda, 
                                          ldb, nrhs);
#endif
    b_bufsize = n*nrhs;
    x_bufsize = b_bufsize;

    if(matrix_layout==LAPACK_COL_MAJOR){
        ldb = n;
        ldx = n;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        ldb = nrhs;
        ldx = nrhs;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &x, &xref, x_bufsize);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &af, &afref, (n*n));
    lapacke_gtest_alloc_double_buffer_pair( &ferr, &ferrref, nrhs);
    lapacke_gtest_alloc_double_buffer_pair( &berr, &berrref, nrhs);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);
    

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (x==NULL) || (xref==NULL) || \
        (af==NULL) || (afref==NULL) || \
        (ferr==NULL) || (ferrref==NULL) || \
        (ipiv==NULL) || (ipivref==NULL) || \
        (berr==NULL) || (berrref==NULL) ){
       EXPECT_FALSE( true) << "pbtrs_dcomplex_parameters object: malloc error.";
       gbrfs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, n*n);
    memcpy(af, a, (n*n*sizeof(lapack_complex_double)));
    memcpy(afref, a, (n*n*sizeof(lapack_complex_double)));
    
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
     memcpy(x, b, ( b_bufsize*sizeof(lapack_complex_double)));
    memcpy(xref, b, ( b_bufsize*sizeof(lapack_complex_double)));
   
    lapacke_gtest_init_double_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);
    lapacke_gtest_init_int_buffer_pair_with_constant(ipiv, ipivref, nrhs, 0);

    
   } /* end of Constructor  */

gbrfs_dcomplex_parameters:: ~gbrfs_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gbrfs_dcomplex_parameters object: destructor invoked. \n");
#endif
   gbrfs_free();
}

//  Test fixture class definition
class zgbrfs_test  : public  ::testing::Test {
public:
   gbrfs_dcomplex_parameters  *zgbrfs_obj;
   void SetUp();  
   void TearDown () { delete zgbrfs_obj; }
};


void zgbrfs_test::SetUp(){

    /* LAPACKE ZGBRFS prototype */
    typedef int (*Fptr_NL_LAPACKE_zgbrfs) (int matrix_layout, char trans,
		lapack_int n, lapack_int kl, lapack_int ku, lapack_int nrhs,
		const lapack_complex_double* ab, lapack_int ldab, const lapack_complex_double* afb,
		lapack_int ldafb, const lapack_int* ipiv, const lapack_complex_double* b,
		lapack_int ldb, lapack_complex_double* x, lapack_int ldx, double* ferr, double* berr);

    Fptr_NL_LAPACKE_zgbrfs ZGBRFS;

     /* LAPACKE ZGBTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zgbtrf) ( int matrix_layout,lapack_int m,
                                lapack_int n, lapack_int kl, lapack_int ku,
                              lapack_complex_double *ab,lapack_int ldab,lapack_int* ipiv );

    Fptr_NL_LAPACKE_zgbtrf ZGBTRF;

    typedef int (*Fptr_NL_LAPACKE_zgbtrs) (int matrix_layout , char trans,
		lapack_int n , lapack_int kl , lapack_int ku , lapack_int nrhs ,
		const lapack_complex_double * ab , lapack_int ldab , const lapack_int * ipiv ,
		lapack_complex_double * b , lapack_int ldb );

    Fptr_NL_LAPACKE_zgbtrs ZGBTRS;

    zgbrfs_obj = new gbrfs_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].kl,
                           lin_solver_paramslist[idx].ku,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    zgbrfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgbrfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgbrfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgbrfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZGBTRF = (Fptr_NL_LAPACKE_zgbtrf)dlsym(zgbrfs_obj->hModule,"LAPACKE_zgbtrf");
    ASSERT_TRUE(ZGBTRF != NULL) << "failed to get the Netlib LAPACKE_zgbtrf symbol";

    ZGBTRS = (Fptr_NL_LAPACKE_zgbtrs)dlsym(zgbrfs_obj->hModule,"LAPACKE_zgbtrs");
    ASSERT_TRUE(ZGBTRS != NULL) << "failed to get the Netlib LAPACKE_zgbtrs symbol";

    ZGBRFS = (Fptr_NL_LAPACKE_zgbrfs)dlsym(zgbrfs_obj->hModule, "LAPACKE_zgbrfs");
    ASSERT_TRUE(ZGBRFS != NULL) << "failed to get the Netlib LAPACKE_zgbrfs symbol";
    
	zgbrfs_obj->inforef = ZGBTRF( zgbrfs_obj->matrix_layout,
								  zgbrfs_obj->n, zgbrfs_obj->n,
								  zgbrfs_obj->kl, zgbrfs_obj->ku,
								  zgbrfs_obj->afref,
								  zgbrfs_obj->lda,
								  zgbrfs_obj->ipivref);
						   
	zgbrfs_obj->info = LAPACKE_zgbtrf( zgbrfs_obj->matrix_layout,
									zgbrfs_obj->n, zgbrfs_obj->n,
								  zgbrfs_obj->kl, zgbrfs_obj->ku,
									   zgbrfs_obj->af,
									   zgbrfs_obj->lda,
									   zgbrfs_obj->ipiv);

	zgbrfs_obj->inforef = ZGBTRS( zgbrfs_obj->matrix_layout,
								  zgbrfs_obj->trans,
								  zgbrfs_obj->n,
								  zgbrfs_obj->kl, zgbrfs_obj->ku,
								  zgbrfs_obj->nrhs,
								  zgbrfs_obj->afref,
								  zgbrfs_obj->lda,
								  zgbrfs_obj->ipivref,
                                  zgbrfs_obj->xref, zgbrfs_obj->ldx
								  );

	zgbrfs_obj->inforef = ZGBTRS( zgbrfs_obj->matrix_layout,
								  zgbrfs_obj->trans,
								  zgbrfs_obj->n,
								  zgbrfs_obj->kl, zgbrfs_obj->ku,
								  zgbrfs_obj->nrhs,
								  zgbrfs_obj->af,
								  zgbrfs_obj->lda,
								  zgbrfs_obj->ipiv,
                                  zgbrfs_obj->x, zgbrfs_obj->ldx
								  );
								  
    /* Compute the reference o/p by invoking Netlib-Lapack's API */
    zgbrfs_obj->inforef = ZGBRFS( zgbrfs_obj->matrix_layout,
                                  zgbrfs_obj->trans, zgbrfs_obj->n,
                                  zgbrfs_obj->kl, zgbrfs_obj->ku,
                                  zgbrfs_obj->nrhs,
                                  zgbrfs_obj->aref, zgbrfs_obj->lda, 
                                  zgbrfs_obj->afref, zgbrfs_obj->ldaf,
                                  zgbrfs_obj->ipivref,
                                  zgbrfs_obj->bref, zgbrfs_obj->ldb,
                                  zgbrfs_obj->xref, zgbrfs_obj->ldx,
                                  zgbrfs_obj->ferrref,
                                  zgbrfs_obj->berrref
								  );
								  
    /* Compute libflame's Lapacke o/p  */								  
    zgbrfs_obj->info = LAPACKE_zgbrfs( zgbrfs_obj->matrix_layout,
                                  zgbrfs_obj->trans, zgbrfs_obj->n,
                                  zgbrfs_obj->kl, zgbrfs_obj->ku,
                                  zgbrfs_obj->nrhs,
                                  zgbrfs_obj->a, zgbrfs_obj->lda, 
                                  zgbrfs_obj->af, zgbrfs_obj->ldaf,
                                  zgbrfs_obj->ipiv,
                                  zgbrfs_obj->b, zgbrfs_obj->ldb,
                                  zgbrfs_obj->x, zgbrfs_obj->ldx,
                                  zgbrfs_obj->ferr,
                                  zgbrfs_obj->berr
								  );


    if( zgbrfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zgbrfs is wrong\n", zgbrfs_obj->info );
    }
    if( zgbrfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zgbrfs is wrong\n", 
        zgbrfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zgbrfs_obj->diff_xerr =  computeDiff_z( zgbrfs_obj->x_bufsize, 
                zgbrfs_obj->x, zgbrfs_obj->xref );

    zgbrfs_obj->diff_berr =  computeDiff_d( zgbrfs_obj->nrhs, 
                zgbrfs_obj->berr, zgbrfs_obj->berrref );
                
    zgbrfs_obj->diff_ferr =  computeDiff_d( zgbrfs_obj->nrhs, 
                zgbrfs_obj->ferr, zgbrfs_obj->ferrref );
}

TEST_F(zgbrfs_test, zgbrfs1) {
    EXPECT_NEAR(0.0, zgbrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zgbrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zgbrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(zgbrfs_test, zgbrfs2) {
    EXPECT_NEAR(0.0, zgbrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zgbrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zgbrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(zgbrfs_test, zgbrfs3) {
    EXPECT_NEAR(0.0, zgbrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zgbrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zgbrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(zgbrfs_test, zgbrfs4) {
    EXPECT_NEAR(0.0, zgbrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zgbrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zgbrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}
