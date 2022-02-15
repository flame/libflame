#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"

#define syrfs_free() \
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
  if (berr != NULL)    free (berr  ); \
  if (berrref != NULL) free (berrref); \
  if (ipiv != NULL)    free (ipiv  ); \
  if (ipivref != NULL) free (ipivref); \
  if( hModule != NULL) dlclose(hModule); \
  if(dModule != NULL) dlclose(dModule)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin syrfs_float_parameters  class definition */
class syrfs_float_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      float diff; // capture difference between ref o/p & libflame lapacke o/p.
      float diff_berr, diff_ferr, diff_xerr;
      void *hModule, *dModule;
      int ipiv_diff;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.

      lapack_int n; // The order of A (Square matrix).
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_int ldaf;  //  leading dimension of 'af'
      lapack_int ldx; //  leading dimension of 'x'

      /* Input/ Output parameters */
      float *b, *bref; //right-hand sides for the systems of equations.
      float *a, *aref; //The array ab contains the matrix A
      float *af, *afref; //contains the ored form of the matrix A
      
      /* Output parameters */
      float* x, *xref;
      float* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      float* berr, *berrref; // component-wise relative backward error for each solution vector xj.
      lapack_int *ipiv, *ipivref; // pivot buffer

      /* Return Values */
      lapack_int info, inforef;

   public: 
      syrfs_float_parameters ( int matrix_layout_i, char uplo_i,
                               lapack_int n_i, lapack_int nrhs_i);
              
      ~syrfs_float_parameters (); 
};  /* end of syrfs_float_parameters  class definition */


/* Constructor syrfs_float_parameters definition */
syrfs_float_parameters:: syrfs_float_parameters ( int matrix_layout_i, 
                 char uplo_i, lapack_int n_i,
                  lapack_int nrhs_i) {
                                    
    matrix_layout = matrix_layout_i;
    uplo = uplo_i;
    // equed acts as input when "  = 'F' " else output
    n = n_i;
    nrhs = nrhs_i;
    lda = n;
    ldaf = n;
    
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

#if LAPACKE_TEST_VERBOSE
   printf(" \n syrfs float:  n: %d, uplo: %c  lda: %d  \
ldb: %d nrhs: %d  ldaf: %d  ldx: %d \n",  n, , uplo, lda, 
ldb, nrhs, ldaf, ldx);
#endif

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
       syrfs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand_custom_matrix ( a, aref, n, n, uplo);
    memcpy(af, a, (n*n*sizeof(float)));
    memcpy(afref, a, (n*n*sizeof(float)));
    
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, b_bufsize);
    memcpy(x, b, ( b_bufsize*sizeof(float)));
    memcpy(xref, b, ( b_bufsize*sizeof(float)));
    
    lapacke_gtest_init_float_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);
    lapacke_gtest_init_int_buffer_pair_with_constant(ipiv, ipivref, nrhs, 0);

    
   } /* end of Constructor  */

syrfs_float_parameters:: ~syrfs_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" syrfs_float_parameters object: destructor invoked. \n");
#endif
   syrfs_free();
}

//  Test fixture class definition
class ssyrfs_test  : public  ::testing::Test {
public:
   syrfs_float_parameters  *ssyrfs_obj;
   void SetUp();  
   void TearDown () { delete ssyrfs_obj; }
};


void ssyrfs_test::SetUp(){

     /* LAPACKE SSYTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_ssytrf) ( int matrix_layout , char uplo,
            lapack_int n , float * a , lapack_int lda , lapack_int * ipiv );
    Fptr_NL_LAPACKE_ssytrf SSYTRF;

    /* LAPACKE SSYRFS prototype */
    typedef int (*Fptr_NL_LAPACKE_ssytrs) ( int matrix_layout , char uplo , lapack_int n , lapack_int nrhs , const float * a , lapack_int lda , const lapack_int * ipiv , float * b , lapack_int ldb );
    Fptr_NL_LAPACKE_ssytrs SSYTRS;

    typedef int (*Fptr_NL_LAPACKE_ssyrfs) ( int matrix_layout,
		char uplo, lapack_int n, lapack_int nrhs, const float* a,
		lapack_int lda, const float* af, lapack_int ldaf,
		const lapack_int* ipiv, const float* b, lapack_int ldb,
		float* x, lapack_int ldx, float* ferr, float* berr );
    Fptr_NL_LAPACKE_ssyrfs SSYRFS;

    ssyrfs_obj = new syrfs_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    ssyrfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ssyrfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ssyrfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ssyrfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	SSYTRF = (Fptr_NL_LAPACKE_ssytrf)dlsym(ssyrfs_obj->hModule,"LAPACKE_ssytrf");
	ASSERT_TRUE(SSYTRF != NULL) << "failed to syt the Netlib LAPACKE_ssytrf symbol";

	SSYTRS = (Fptr_NL_LAPACKE_ssytrs)dlsym(ssyrfs_obj->hModule,"LAPACKE_ssytrs");
	ASSERT_TRUE(SSYTRS != NULL) << "failed to syt the Netlib LAPACKE_ssytrs symbol";

    SSYRFS = (Fptr_NL_LAPACKE_ssyrfs)dlsym(ssyrfs_obj->hModule, "LAPACKE_ssyrfs");
    ASSERT_TRUE(SSYRFS != NULL) << "failed to syt the Netlib LAPACKE_ssyrfs symbol";

    /* Generate i/ps to syrfs API from 'sytrf' & 'sytrs' API sequence */
	ssyrfs_obj->inforef = SSYTRF( ssyrfs_obj->matrix_layout,
								  ssyrfs_obj->uplo, ssyrfs_obj->n,
								  ssyrfs_obj->afref,
								  ssyrfs_obj->lda,
								  ssyrfs_obj->ipivref);
						   
	ssyrfs_obj->info = LAPACKE_ssytrf( ssyrfs_obj->matrix_layout,
									   ssyrfs_obj->uplo, ssyrfs_obj->n,
									   ssyrfs_obj->af,
									   ssyrfs_obj->lda,
									   ssyrfs_obj->ipiv);

	ssyrfs_obj->inforef = SSYTRS( ssyrfs_obj->matrix_layout,
								  ssyrfs_obj->uplo,
								  ssyrfs_obj->n,
								  ssyrfs_obj->nrhs,
								  ssyrfs_obj->afref,
								  ssyrfs_obj->lda,
								  ssyrfs_obj->ipivref,
								  ssyrfs_obj->xref,
								  ssyrfs_obj->ldx
								  );
						   
	ssyrfs_obj->info = LAPACKE_ssytrs( ssyrfs_obj->matrix_layout,
								  ssyrfs_obj->uplo,
								  ssyrfs_obj->n,
								  ssyrfs_obj->nrhs,
								  ssyrfs_obj->af,
								  ssyrfs_obj->lda,
								  ssyrfs_obj->ipiv,
								  ssyrfs_obj->x,
								  ssyrfs_obj->ldx
								  );

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    ssyrfs_obj->inforef = SSYRFS( ssyrfs_obj->matrix_layout,
                                  ssyrfs_obj->uplo, ssyrfs_obj->n,
                                  ssyrfs_obj->nrhs,
                                  ssyrfs_obj->aref, ssyrfs_obj->lda, 
                                  ssyrfs_obj->afref, ssyrfs_obj->ldaf,
                                  ssyrfs_obj->ipivref,
                                  ssyrfs_obj->bref, ssyrfs_obj->ldb,
                                  ssyrfs_obj->xref, ssyrfs_obj->ldx,
                                  ssyrfs_obj->ferrref,
                                  ssyrfs_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    ssyrfs_obj->info = LAPACKE_ssyrfs( ssyrfs_obj->matrix_layout,
                                  ssyrfs_obj->uplo, ssyrfs_obj->n,
                                  ssyrfs_obj->nrhs,
                                  ssyrfs_obj->a, ssyrfs_obj->lda, 
                                  ssyrfs_obj->af, ssyrfs_obj->ldaf,
                                  ssyrfs_obj->ipiv,
                                  ssyrfs_obj->b, ssyrfs_obj->ldb,
                                  ssyrfs_obj->x, ssyrfs_obj->ldx,
                                  ssyrfs_obj->ferr,
                                  ssyrfs_obj->berr);

    if( ssyrfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ssyrfs is wrong\n", ssyrfs_obj->info );
    }
    if( ssyrfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ssyrfs is wrong\n", 
        ssyrfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    ssyrfs_obj->ipiv_diff = computeDiff_i( ssyrfs_obj->n, ssyrfs_obj->ipiv, ssyrfs_obj->ipivref );
    
    ssyrfs_obj->diff_xerr =  computeDiff_s( ssyrfs_obj->x_bufsize, 
                ssyrfs_obj->x, ssyrfs_obj->xref );

    ssyrfs_obj->diff_berr =  computeDiff_s( ssyrfs_obj->nrhs, 
                ssyrfs_obj->berr, ssyrfs_obj->berrref );
                
    ssyrfs_obj->diff_ferr =  computeDiff_s( ssyrfs_obj->nrhs, 
                ssyrfs_obj->ferr, ssyrfs_obj->ferrref );
}

TEST_F(ssyrfs_test, ssyrfs1) {
    EXPECT_NEAR(0.0, ssyrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, ssyrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, ssyrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, ssyrfs_obj->ipiv_diff );
    idx = Circular_Increment_Index(idx);
}

TEST_F(ssyrfs_test, ssyrfs2) {
    EXPECT_NEAR(0.0, ssyrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, ssyrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, ssyrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, ssyrfs_obj->ipiv_diff );
    idx = Circular_Increment_Index(idx);
}

TEST_F(ssyrfs_test, ssyrfs3) {
    EXPECT_NEAR(0.0, ssyrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, ssyrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, ssyrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, ssyrfs_obj->ipiv_diff );
    idx = Circular_Increment_Index(idx);
}

TEST_F(ssyrfs_test, ssyrfs4) {
    EXPECT_NEAR(0.0, ssyrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, ssyrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, ssyrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, ssyrfs_obj->ipiv_diff );
    idx = Circular_Increment_Index(idx);
}

/* Begin syrfs_double_parameters  class definition */
class syrfs_double_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      double diff; // capture difference between ref o/p & libflame lapacke o/p.
      double diff_berr, diff_ferr, diff_xerr;
      void *hModule, *dModule;
      int ipiv_diff;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.

      lapack_int n; // The order of A (Square matrix).
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_int ldaf;  //  leading dimension of 'af'
      lapack_int ldx; //  leading dimension of 'x'

      /* Input/ Output parameters */
      double *b, *bref; //right-hand sides for the systems of equations.
      double *a, *aref; //The array ab contains the matrix A
      double *af, *afref; //contains the ored form of the matrix A
      
      /* Output parameters */
      double* x, *xref;
      double* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      double* berr, *berrref; // component-wise relative backward error for each solution vector xj.
      lapack_int *ipiv, *ipivref; // pivot buffer

      /* Return Values */
      lapack_int info, inforef;

   public: 
      syrfs_double_parameters ( int matrix_layout_i, char uplo_i,
                               lapack_int n_i, lapack_int nrhs_i);
              
      ~syrfs_double_parameters (); 
};  /* end of syrfs_double_parameters  class definition */


/* Constructor syrfs_double_parameters definition */
syrfs_double_parameters:: syrfs_double_parameters ( int matrix_layout_i, 
                 char uplo_i, lapack_int n_i,
                  lapack_int nrhs_i) {
                                    
    matrix_layout = matrix_layout_i;
    uplo = uplo_i;
    // equed acts as input when "  = 'F' " else output
    n = n_i;
    nrhs = nrhs_i;
    lda = n;
    ldaf = n;
    
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

#if LAPACKE_TEST_VERBOSE
   printf(" \n syrfs double:  n: %d, uplo: %c  lda: %d  \
ldb: %d nrhs: %d  ldaf: %d  ldx: %d \n",  n, , uplo, lda, 
ldb, nrhs, ldaf, ldx);
#endif

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
       syrfs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand_custom_matrix ( a, aref, n, n, uplo);
    memcpy(af, a, (n*n*sizeof(double)));
    memcpy(afref, a, (n*n*sizeof(double)));
    
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, b_bufsize);
    memcpy(x, b, ( b_bufsize*sizeof(double)));
    memcpy(xref, b, ( b_bufsize*sizeof(double)));
    
    lapacke_gtest_init_double_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);
    lapacke_gtest_init_int_buffer_pair_with_constant(ipiv, ipivref, nrhs, 0);

    
   } /* end of Constructor  */

syrfs_double_parameters:: ~syrfs_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" syrfs_double_parameters object: destructor invoked. \n");
#endif
   syrfs_free();
}

//  Test fixture class definition
class dsyrfs_test  : public  ::testing::Test {
public:
   syrfs_double_parameters  *dsyrfs_obj;
   void SetUp();  
   void TearDown () { delete dsyrfs_obj; }
};


void dsyrfs_test::SetUp(){

     /* LAPACKE DSYTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_dsytrf) ( int matrix_layout , char uplo,
            lapack_int n , double * a , lapack_int lda , lapack_int * ipiv );
    Fptr_NL_LAPACKE_dsytrf DSYTRF;

    /* LAPACKE DSYRFS prototype */
    typedef int (*Fptr_NL_LAPACKE_dsytrs) ( int matrix_layout , char uplo , lapack_int n , lapack_int nrhs , const double * a , lapack_int lda , const lapack_int * ipiv , double * b , lapack_int ldb );
    Fptr_NL_LAPACKE_dsytrs DSYTRS;

    typedef int (*Fptr_NL_LAPACKE_dsyrfs) ( int matrix_layout,
		char uplo, lapack_int n, lapack_int nrhs, const double* a,
		lapack_int lda, const double* af, lapack_int ldaf,
		const lapack_int* ipiv, const double* b, lapack_int ldb,
		double* x, lapack_int ldx, double* ferr, double* berr );
    Fptr_NL_LAPACKE_dsyrfs DSYRFS;

    dsyrfs_obj = new syrfs_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    dsyrfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dsyrfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dsyrfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dsyrfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	DSYTRF = (Fptr_NL_LAPACKE_dsytrf)dlsym(dsyrfs_obj->hModule,"LAPACKE_dsytrf");
	ASSERT_TRUE(DSYTRF != NULL) << "failed to syt the Netlib LAPACKE_dsytrf symbol";

	DSYTRS = (Fptr_NL_LAPACKE_dsytrs)dlsym(dsyrfs_obj->hModule,"LAPACKE_dsytrs");
	ASSERT_TRUE(DSYTRS != NULL) << "failed to syt the Netlib LAPACKE_dsytrs symbol";

    DSYRFS = (Fptr_NL_LAPACKE_dsyrfs)dlsym(dsyrfs_obj->hModule, "LAPACKE_dsyrfs");
    ASSERT_TRUE(DSYRFS != NULL) << "failed to syt the Netlib LAPACKE_dsyrfs symbol";

    /* Generate i/ps to syrfs API from 'sytrf' & 'sytrs' API sequence */
	dsyrfs_obj->inforef = DSYTRF( dsyrfs_obj->matrix_layout,
								  dsyrfs_obj->uplo, dsyrfs_obj->n,
								  dsyrfs_obj->afref,
								  dsyrfs_obj->lda,
								  dsyrfs_obj->ipivref);
						   
	dsyrfs_obj->info = LAPACKE_dsytrf( dsyrfs_obj->matrix_layout,
									   dsyrfs_obj->uplo, dsyrfs_obj->n,
									   dsyrfs_obj->af,
									   dsyrfs_obj->lda,
									   dsyrfs_obj->ipiv);

	dsyrfs_obj->inforef = DSYTRS( dsyrfs_obj->matrix_layout,
								  dsyrfs_obj->uplo,
								  dsyrfs_obj->n,
								  dsyrfs_obj->nrhs,
								  dsyrfs_obj->afref,
								  dsyrfs_obj->lda,
								  dsyrfs_obj->ipivref,
								  dsyrfs_obj->xref,
								  dsyrfs_obj->ldx
								  );
						   
	dsyrfs_obj->info = LAPACKE_dsytrs( dsyrfs_obj->matrix_layout,
								  dsyrfs_obj->uplo,
								  dsyrfs_obj->n,
								  dsyrfs_obj->nrhs,
								  dsyrfs_obj->af,
								  dsyrfs_obj->lda,
								  dsyrfs_obj->ipiv,
								  dsyrfs_obj->x,
								  dsyrfs_obj->ldx
								  );

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    dsyrfs_obj->inforef = DSYRFS( dsyrfs_obj->matrix_layout,
                                  dsyrfs_obj->uplo, dsyrfs_obj->n,
                                  dsyrfs_obj->nrhs,
                                  dsyrfs_obj->aref, dsyrfs_obj->lda, 
                                  dsyrfs_obj->afref, dsyrfs_obj->ldaf,
                                  dsyrfs_obj->ipivref,
                                  dsyrfs_obj->bref, dsyrfs_obj->ldb,
                                  dsyrfs_obj->xref, dsyrfs_obj->ldx,
                                  dsyrfs_obj->ferrref,
                                  dsyrfs_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    dsyrfs_obj->info = LAPACKE_dsyrfs( dsyrfs_obj->matrix_layout,
                                  dsyrfs_obj->uplo, dsyrfs_obj->n,
                                  dsyrfs_obj->nrhs,
                                  dsyrfs_obj->a, dsyrfs_obj->lda, 
                                  dsyrfs_obj->af, dsyrfs_obj->ldaf,
                                  dsyrfs_obj->ipiv,
                                  dsyrfs_obj->b, dsyrfs_obj->ldb,
                                  dsyrfs_obj->x, dsyrfs_obj->ldx,
                                  dsyrfs_obj->ferr,
                                  dsyrfs_obj->berr);

    if( dsyrfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dsyrfs is wrong\n", dsyrfs_obj->info );
    }
    if( dsyrfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dsyrfs is wrong\n", 
        dsyrfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dsyrfs_obj->ipiv_diff = computeDiff_i( dsyrfs_obj->n, dsyrfs_obj->ipiv, dsyrfs_obj->ipivref );
    
    dsyrfs_obj->diff_xerr =  computeDiff_d( dsyrfs_obj->x_bufsize, 
                dsyrfs_obj->x, dsyrfs_obj->xref );

    dsyrfs_obj->diff_berr =  computeDiff_d( dsyrfs_obj->nrhs, 
                dsyrfs_obj->berr, dsyrfs_obj->berrref );
                
    dsyrfs_obj->diff_ferr =  computeDiff_d( dsyrfs_obj->nrhs, 
                dsyrfs_obj->ferr, dsyrfs_obj->ferrref );
}

TEST_F(dsyrfs_test, dsyrfs1) {
    EXPECT_NEAR(0.0, dsyrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dsyrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dsyrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, dsyrfs_obj->ipiv_diff );
    idx = Circular_Increment_Index(idx);
}

TEST_F(dsyrfs_test, dsyrfs2) {
    EXPECT_NEAR(0.0, dsyrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dsyrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dsyrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, dsyrfs_obj->ipiv_diff );
    idx = Circular_Increment_Index(idx);
}

TEST_F(dsyrfs_test, dsyrfs3) {
    EXPECT_NEAR(0.0, dsyrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dsyrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dsyrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, dsyrfs_obj->ipiv_diff );
    idx = Circular_Increment_Index(idx);
}

TEST_F(dsyrfs_test, dsyrfs4) {
    EXPECT_NEAR(0.0, dsyrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dsyrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dsyrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, dsyrfs_obj->ipiv_diff );
    idx = Circular_Increment_Index(idx);
}

/* Begin syrfs_scomplex_parameters  class definition */
class syrfs_scomplex_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      float diff_berr, diff_ferr, diff_xerr;
      void *hModule, *dModule;
      int ipiv_diff;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.

      lapack_int n; // The order of A (Square matrix).
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_int ldaf;  //  leading dimension of 'af'
      lapack_int ldx; //  leading dimension of 'x'

      /* Input/ Output parameters */
      lapack_complex_float *b, *bref; //right-hand sides for the systems of equations.
      lapack_complex_float *a, *aref; //The array ab contains the matrix A
      lapack_complex_float *af, *afref; //contains the ored form of the matrix A
      
      /* Output parameters */
      lapack_complex_float* x, *xref;
      float* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      float* berr, *berrref; // component-wise relative backward error for each solution vector xj.
      lapack_int *ipiv, *ipivref; // pivot buffer

      /* Return Values */
      lapack_int info, inforef;

   public: 
      syrfs_scomplex_parameters ( int matrix_layout_i, char uplo_i,
                               lapack_int n_i, lapack_int nrhs_i);
              
      ~syrfs_scomplex_parameters (); 
};  /* end of syrfs_scomplex_parameters  class definition */


/* Constructor syrfs_scomplex_parameters definition */
syrfs_scomplex_parameters:: syrfs_scomplex_parameters ( int matrix_layout_i, 
                 char uplo_i, lapack_int n_i,
                  lapack_int nrhs_i) {
                                    
    matrix_layout = matrix_layout_i;
    uplo = uplo_i;
    // equed acts as input when "  = 'F' " else output
    n = n_i;
    nrhs = nrhs_i;
    lda = n;
    ldaf = n;
    
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

#if LAPACKE_TEST_VERBOSE
   printf(" \n syrfs lapack_complex_float:  n: %d, uplo: %c  lda: %d  \
ldb: %d nrhs: %d  ldaf: %d  ldx: %d \n",  n, , uplo, lda, 
ldb, nrhs, ldaf, ldx);
#endif

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
       syrfs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix ( a, aref, n, n, uplo);
    memcpy(af, a, (n*n*sizeof(lapack_complex_float)));
    memcpy(afref, a, (n*n*sizeof(lapack_complex_float)));
    
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
    memcpy(x, b, ( b_bufsize*sizeof(lapack_complex_float)));
    memcpy(xref, b, ( b_bufsize*sizeof(lapack_complex_float)));
    
    lapacke_gtest_init_float_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);
    lapacke_gtest_init_int_buffer_pair_with_constant(ipiv, ipivref, nrhs, 0);

    
   } /* end of Constructor  */

syrfs_scomplex_parameters:: ~syrfs_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" syrfs_scomplex_parameters object: destructor invoked. \n");
#endif
   syrfs_free();
}

//  Test fixture class definition
class csyrfs_test  : public  ::testing::Test {
public:
   syrfs_scomplex_parameters  *csyrfs_obj;
   void SetUp();  
   void TearDown () { delete csyrfs_obj; }
};


void csyrfs_test::SetUp(){

     /* LAPACKE CSYTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_csytrf) ( int matrix_layout , char uplo,
            lapack_int n , lapack_complex_float * a , lapack_int lda , lapack_int * ipiv );
    Fptr_NL_LAPACKE_csytrf CSYTRF;

    /* LAPACKE CSYRFS prototype */
    typedef int (*Fptr_NL_LAPACKE_csytrs) ( int matrix_layout , char uplo , lapack_int n , lapack_int nrhs , const lapack_complex_float * a , lapack_int lda , const lapack_int * ipiv , lapack_complex_float * b , lapack_int ldb );
    Fptr_NL_LAPACKE_csytrs CSYTRS;

    typedef int (*Fptr_NL_LAPACKE_csyrfs) ( int matrix_layout,
		char uplo, lapack_int n, lapack_int nrhs, const lapack_complex_float* a,
		lapack_int lda, const lapack_complex_float* af, lapack_int ldaf,
		const lapack_int* ipiv, const lapack_complex_float* b, lapack_int ldb,
		lapack_complex_float* x, lapack_int ldx, float* ferr, float* berr );
    Fptr_NL_LAPACKE_csyrfs CSYRFS;

    csyrfs_obj = new syrfs_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    csyrfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    csyrfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(csyrfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(csyrfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	CSYTRF = (Fptr_NL_LAPACKE_csytrf)dlsym(csyrfs_obj->hModule,"LAPACKE_csytrf");
	ASSERT_TRUE(CSYTRF != NULL) << "failed to syt the Netlib LAPACKE_csytrf symbol";

	CSYTRS = (Fptr_NL_LAPACKE_csytrs)dlsym(csyrfs_obj->hModule,"LAPACKE_csytrs");
	ASSERT_TRUE(CSYTRS != NULL) << "failed to syt the Netlib LAPACKE_csytrs symbol";

    CSYRFS = (Fptr_NL_LAPACKE_csyrfs)dlsym(csyrfs_obj->hModule, "LAPACKE_csyrfs");
    ASSERT_TRUE(CSYRFS != NULL) << "failed to syt the Netlib LAPACKE_csyrfs symbol";

    /* Generate i/ps to syrfs API from 'sytrf' & 'sytrs' API sequence */
	csyrfs_obj->inforef = CSYTRF( csyrfs_obj->matrix_layout,
								  csyrfs_obj->uplo, csyrfs_obj->n,
								  csyrfs_obj->afref,
								  csyrfs_obj->lda,
								  csyrfs_obj->ipivref);
						   
	csyrfs_obj->info = LAPACKE_csytrf( csyrfs_obj->matrix_layout,
									   csyrfs_obj->uplo, csyrfs_obj->n,
									   csyrfs_obj->af,
									   csyrfs_obj->lda,
									   csyrfs_obj->ipiv);

	csyrfs_obj->inforef = CSYTRS( csyrfs_obj->matrix_layout,
								  csyrfs_obj->uplo,
								  csyrfs_obj->n,
								  csyrfs_obj->nrhs,
								  csyrfs_obj->afref,
								  csyrfs_obj->lda,
								  csyrfs_obj->ipivref,
								  csyrfs_obj->xref,
								  csyrfs_obj->ldx
								  );
						   
	csyrfs_obj->info = LAPACKE_csytrs( csyrfs_obj->matrix_layout,
								  csyrfs_obj->uplo,
								  csyrfs_obj->n,
								  csyrfs_obj->nrhs,
								  csyrfs_obj->af,
								  csyrfs_obj->lda,
								  csyrfs_obj->ipiv,
								  csyrfs_obj->x,
								  csyrfs_obj->ldx
								  );

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    csyrfs_obj->inforef = CSYRFS( csyrfs_obj->matrix_layout,
                                  csyrfs_obj->uplo, csyrfs_obj->n,
                                  csyrfs_obj->nrhs,
                                  csyrfs_obj->aref, csyrfs_obj->lda, 
                                  csyrfs_obj->afref, csyrfs_obj->ldaf,
                                  csyrfs_obj->ipivref,
                                  csyrfs_obj->bref, csyrfs_obj->ldb,
                                  csyrfs_obj->xref, csyrfs_obj->ldx,
                                  csyrfs_obj->ferrref,
                                  csyrfs_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    csyrfs_obj->info = LAPACKE_csyrfs( csyrfs_obj->matrix_layout,
                                  csyrfs_obj->uplo, csyrfs_obj->n,
                                  csyrfs_obj->nrhs,
                                  csyrfs_obj->a, csyrfs_obj->lda, 
                                  csyrfs_obj->af, csyrfs_obj->ldaf,
                                  csyrfs_obj->ipiv,
                                  csyrfs_obj->b, csyrfs_obj->ldb,
                                  csyrfs_obj->x, csyrfs_obj->ldx,
                                  csyrfs_obj->ferr,
                                  csyrfs_obj->berr);

    if( csyrfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_csyrfs is wrong\n", csyrfs_obj->info );
    }
    if( csyrfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_csyrfs is wrong\n", 
        csyrfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    csyrfs_obj->ipiv_diff = computeDiff_i( csyrfs_obj->n, csyrfs_obj->ipiv, csyrfs_obj->ipivref );
    
    csyrfs_obj->diff_xerr =  computeDiff_c( csyrfs_obj->x_bufsize, 
                csyrfs_obj->x, csyrfs_obj->xref );

    csyrfs_obj->diff_berr =  computeDiff_s( csyrfs_obj->nrhs, 
                csyrfs_obj->berr, csyrfs_obj->berrref );
                
    csyrfs_obj->diff_ferr =  computeDiff_s( csyrfs_obj->nrhs, 
                csyrfs_obj->ferr, csyrfs_obj->ferrref );
}

TEST_F(csyrfs_test, csyrfs1) {
    EXPECT_NEAR(0.0, csyrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, csyrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, csyrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, csyrfs_obj->ipiv_diff );
    idx = Circular_Increment_Index(idx);
}

TEST_F(csyrfs_test, csyrfs2) {
    EXPECT_NEAR(0.0, csyrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, csyrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, csyrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, csyrfs_obj->ipiv_diff );
    idx = Circular_Increment_Index(idx);
}

TEST_F(csyrfs_test, csyrfs3) {
    EXPECT_NEAR(0.0, csyrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, csyrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, csyrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, csyrfs_obj->ipiv_diff );
    idx = Circular_Increment_Index(idx);
}

TEST_F(csyrfs_test, csyrfs4) {
    EXPECT_NEAR(0.0, csyrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, csyrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, csyrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, csyrfs_obj->ipiv_diff );
    idx = Circular_Increment_Index(idx);
}

/* Begin syrfs_dcomplex_parameters  class definition */
class syrfs_dcomplex_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      double diff_berr, diff_ferr, diff_xerr;
      void *hModule, *dModule;
      int ipiv_diff;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.

      lapack_int n; // The order of A (Square matrix).
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_int ldaf;  //  leading dimension of 'af'
      lapack_int ldx; //  leading dimension of 'x'

      /* Input/ Output parameters */
      lapack_complex_double *b, *bref; //right-hand sides for the systems of equations.
      lapack_complex_double *a, *aref; //The array ab contains the matrix A
      lapack_complex_double *af, *afref; //contains the ored form of the matrix A
      
      /* Output parameters */
      lapack_complex_double* x, *xref;
      double* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      double* berr, *berrref; // component-wise relative backward error for each solution vector xj.
      lapack_int *ipiv, *ipivref; // pivot buffer

      /* Return Values */
      lapack_int info, inforef;

   public: 
      syrfs_dcomplex_parameters ( int matrix_layout_i, char uplo_i,
                               lapack_int n_i, lapack_int nrhs_i);
              
      ~syrfs_dcomplex_parameters (); 
};  /* end of syrfs_dcomplex_parameters  class definition */


/* Constructor syrfs_dcomplex_parameters definition */
syrfs_dcomplex_parameters:: syrfs_dcomplex_parameters ( int matrix_layout_i, 
                 char uplo_i, lapack_int n_i,
                  lapack_int nrhs_i) {
                                    
    matrix_layout = matrix_layout_i;
    uplo = uplo_i;
    // equed acts as input when "  = 'F' " else output
    n = n_i;
    nrhs = nrhs_i;
    lda = n;
    ldaf = n;
    
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

#if LAPACKE_TEST_VERBOSE
   printf(" \n syrfs lapack_complex_double:  n: %d, uplo: %c  lda: %d  \
ldb: %d nrhs: %d  ldaf: %d  ldx: %d \n",  n, , uplo, lda, 
ldb, nrhs, ldaf, ldx);
#endif

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
       syrfs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix ( a, aref, n, n, uplo);
    memcpy(af, a, (n*n*sizeof(lapack_complex_double)));
    memcpy(afref, a, (n*n*sizeof(lapack_complex_double)));
    
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
    memcpy(x, b, ( b_bufsize*sizeof(lapack_complex_double)));
    memcpy(xref, b, ( b_bufsize*sizeof(lapack_complex_double)));
    
    lapacke_gtest_init_double_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);
    lapacke_gtest_init_int_buffer_pair_with_constant(ipiv, ipivref, nrhs, 0);

    
   } /* end of Constructor  */

syrfs_dcomplex_parameters:: ~syrfs_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" syrfs_dcomplex_parameters object: destructor invoked. \n");
#endif
   syrfs_free();
}

//  Test fixture class definition
class zsyrfs_test  : public  ::testing::Test {
public:
   syrfs_dcomplex_parameters  *zsyrfs_obj;
   void SetUp();  
   void TearDown () { delete zsyrfs_obj; }
};


void zsyrfs_test::SetUp(){

     /* LAPACKE ZSYTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zsytrf) ( int matrix_layout , char uplo,
            lapack_int n , lapack_complex_double * a , lapack_int lda , lapack_int * ipiv );
    Fptr_NL_LAPACKE_zsytrf ZSYTRF;

    /* LAPACKE ZSYRFS prototype */
    typedef int (*Fptr_NL_LAPACKE_zsytrs) ( int matrix_layout , char uplo , lapack_int n , lapack_int nrhs , const lapack_complex_double * a , lapack_int lda , const lapack_int * ipiv , lapack_complex_double * b , lapack_int ldb );
    Fptr_NL_LAPACKE_zsytrs ZSYTRS;

    typedef int (*Fptr_NL_LAPACKE_zsyrfs) ( int matrix_layout,
		char uplo, lapack_int n, lapack_int nrhs, const lapack_complex_double* a,
		lapack_int lda, const lapack_complex_double* af, lapack_int ldaf,
		const lapack_int* ipiv, const lapack_complex_double* b, lapack_int ldb,
		lapack_complex_double* x, lapack_int ldx, double* ferr, double* berr );
    Fptr_NL_LAPACKE_zsyrfs ZSYRFS;

    zsyrfs_obj = new syrfs_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    zsyrfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zsyrfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zsyrfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zsyrfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	ZSYTRF = (Fptr_NL_LAPACKE_zsytrf)dlsym(zsyrfs_obj->hModule,"LAPACKE_zsytrf");
	ASSERT_TRUE(ZSYTRF != NULL) << "failed to syt the Netlib LAPACKE_zsytrf symbol";

	ZSYTRS = (Fptr_NL_LAPACKE_zsytrs)dlsym(zsyrfs_obj->hModule,"LAPACKE_zsytrs");
	ASSERT_TRUE(ZSYTRS != NULL) << "failed to syt the Netlib LAPACKE_zsytrs symbol";

    ZSYRFS = (Fptr_NL_LAPACKE_zsyrfs)dlsym(zsyrfs_obj->hModule, "LAPACKE_zsyrfs");
    ASSERT_TRUE(ZSYRFS != NULL) << "failed to syt the Netlib LAPACKE_zsyrfs symbol";

    /* Generate i/ps to syrfs API from 'sytrf' & 'sytrs' API sequence */
	zsyrfs_obj->inforef = ZSYTRF( zsyrfs_obj->matrix_layout,
								  zsyrfs_obj->uplo, zsyrfs_obj->n,
								  zsyrfs_obj->afref,
								  zsyrfs_obj->lda,
								  zsyrfs_obj->ipivref);
						   
	zsyrfs_obj->info = LAPACKE_zsytrf( zsyrfs_obj->matrix_layout,
									   zsyrfs_obj->uplo, zsyrfs_obj->n,
									   zsyrfs_obj->af,
									   zsyrfs_obj->lda,
									   zsyrfs_obj->ipiv);

	zsyrfs_obj->inforef = ZSYTRS( zsyrfs_obj->matrix_layout,
								  zsyrfs_obj->uplo,
								  zsyrfs_obj->n,
								  zsyrfs_obj->nrhs,
								  zsyrfs_obj->afref,
								  zsyrfs_obj->lda,
								  zsyrfs_obj->ipivref,
								  zsyrfs_obj->xref,
								  zsyrfs_obj->ldx
								  );
						   
	zsyrfs_obj->info = LAPACKE_zsytrs( zsyrfs_obj->matrix_layout,
								  zsyrfs_obj->uplo,
								  zsyrfs_obj->n,
								  zsyrfs_obj->nrhs,
								  zsyrfs_obj->af,
								  zsyrfs_obj->lda,
								  zsyrfs_obj->ipiv,
								  zsyrfs_obj->x,
								  zsyrfs_obj->ldx
								  );

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    zsyrfs_obj->inforef = ZSYRFS( zsyrfs_obj->matrix_layout,
                                  zsyrfs_obj->uplo, zsyrfs_obj->n,
                                  zsyrfs_obj->nrhs,
                                  zsyrfs_obj->aref, zsyrfs_obj->lda, 
                                  zsyrfs_obj->afref, zsyrfs_obj->ldaf,
                                  zsyrfs_obj->ipivref,
                                  zsyrfs_obj->bref, zsyrfs_obj->ldb,
                                  zsyrfs_obj->xref, zsyrfs_obj->ldx,
                                  zsyrfs_obj->ferrref,
                                  zsyrfs_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    zsyrfs_obj->info = LAPACKE_zsyrfs( zsyrfs_obj->matrix_layout,
                                  zsyrfs_obj->uplo, zsyrfs_obj->n,
                                  zsyrfs_obj->nrhs,
                                  zsyrfs_obj->a, zsyrfs_obj->lda, 
                                  zsyrfs_obj->af, zsyrfs_obj->ldaf,
                                  zsyrfs_obj->ipiv,
                                  zsyrfs_obj->b, zsyrfs_obj->ldb,
                                  zsyrfs_obj->x, zsyrfs_obj->ldx,
                                  zsyrfs_obj->ferr,
                                  zsyrfs_obj->berr);

    if( zsyrfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zsyrfs is wrong\n", zsyrfs_obj->info );
    }
    if( zsyrfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zsyrfs is wrong\n", 
        zsyrfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zsyrfs_obj->ipiv_diff = computeDiff_i( zsyrfs_obj->n, zsyrfs_obj->ipiv, zsyrfs_obj->ipivref );
    
    zsyrfs_obj->diff_xerr =  computeDiff_z( zsyrfs_obj->x_bufsize, 
                zsyrfs_obj->x, zsyrfs_obj->xref );

    zsyrfs_obj->diff_berr =  computeDiff_d( zsyrfs_obj->nrhs, 
                zsyrfs_obj->berr, zsyrfs_obj->berrref );
                
    zsyrfs_obj->diff_ferr =  computeDiff_d( zsyrfs_obj->nrhs, 
                zsyrfs_obj->ferr, zsyrfs_obj->ferrref );
}

TEST_F(zsyrfs_test, zsyrfs1) {
    EXPECT_NEAR(0.0, zsyrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zsyrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zsyrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, zsyrfs_obj->ipiv_diff );
    idx = Circular_Increment_Index(idx);
}

TEST_F(zsyrfs_test, zsyrfs2) {
    EXPECT_NEAR(0.0, zsyrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zsyrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zsyrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, zsyrfs_obj->ipiv_diff );
    idx = Circular_Increment_Index(idx);
}

TEST_F(zsyrfs_test, zsyrfs3) {
    EXPECT_NEAR(0.0, zsyrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zsyrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zsyrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, zsyrfs_obj->ipiv_diff );
    idx = Circular_Increment_Index(idx);
}

TEST_F(zsyrfs_test, zsyrfs4) {
    EXPECT_NEAR(0.0, zsyrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zsyrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zsyrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, zsyrfs_obj->ipiv_diff );
    idx = Circular_Increment_Index(idx);
}
