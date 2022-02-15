#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"

#define sprfs_free() \
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

/* Begin sprfs_float_parameters  class definition */
class sprfs_float_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      int a_bufsize;
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
      float *af, *afref; //contains the factored form of the matrix A
      
      /* Output parameters */
      float* x, *xref;
      float* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      float* berr, *berrref; // component-wise relative backward error for each solution vector xj.
      lapack_int *ipiv, *ipivref; // pivot buffer

      /* Return Values */
      lapack_int info, inforef;

   public: 
      sprfs_float_parameters ( int matrix_layout_i, char uplo_i,
                               lapack_int n_i, lapack_int nrhs_i);
              
      ~sprfs_float_parameters (); 
};  /* end of sprfs_float_parameters  class definition */


/* Constructor sprfs_float_parameters definition */
sprfs_float_parameters:: sprfs_float_parameters ( int matrix_layout_i, 
					char uplo_i, lapack_int n_i, lapack_int nrhs_i) {
                                    
    matrix_layout = matrix_layout_i;
    uplo = uplo_i;
    // equed acts as input when " fact = 'F' " else output
    n = n_i;
    nrhs = nrhs_i;
    lda = n;
    ldaf = n;
    
    a_bufsize = n*(n+1)/2;
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
   printf(" \n sprfs float:  n: %d, uplo: %c  lda: %d  \
ldb: %d nrhs: %d  ldaf: %d  ldx: %d \n",  n, uplo, lda, 
                                          ldb, nrhs, ldaf, ldx);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, a_bufsize);
    lapacke_gtest_alloc_float_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_float_buffer_pair( &x, &xref, x_bufsize);
    lapacke_gtest_alloc_float_buffer_pair( &af, &afref, a_bufsize);
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
       sprfs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand ( a, aref, a_bufsize);
    memcpy(af, a, (a_bufsize*sizeof(float)));
    memcpy(afref, a, (a_bufsize*sizeof(float)));
    
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, b_bufsize);
    memcpy(x, b, ( b_bufsize*sizeof(float)));
    memcpy(xref, b, ( b_bufsize*sizeof(float)));
    
    lapacke_gtest_init_float_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);

    
   } /* end of Constructor  */

sprfs_float_parameters:: ~sprfs_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sprfs_float_parameters object: destructor invoked. \n");
#endif
   sprfs_free();
}

//  Test fixture class definition
class ssprfs_test  : public  ::testing::Test {
public:
   sprfs_float_parameters  *ssprfs_obj;
   void SetUp();  
   void TearDown () { delete ssprfs_obj; }
};


void ssprfs_test::SetUp(){

     /* LAPACKE SSPTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_ssptrf) ( int matrix_layout , char uplo,
                            lapack_int n , float * ap , lapack_int * ipiv );

    Fptr_NL_LAPACKE_ssptrf SSPTRF;

    typedef int (*Fptr_NL_LAPACKE_ssptrs) ( int matrix_layout , char uplo, 
					   lapack_int n , lapack_int nrhs , const float * ap , 
					const lapack_int * ipiv , float * b , lapack_int ldb );

    Fptr_NL_LAPACKE_ssptrs SSPTRS;

    /* LAPACKE SSPRFS prototype */
    typedef int (*Fptr_NL_LAPACKE_ssprfs) (int matrix_layout,
		char uplo, lapack_int n, lapack_int nrhs, const float* ap,
		const float* afp, const lapack_int* ipiv, const float* b,
		lapack_int ldb, float* x, lapack_int ldx, float* ferr,
		float* berr);

    Fptr_NL_LAPACKE_ssprfs SSPRFS;

    ssprfs_obj = new sprfs_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    ssprfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ssprfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ssprfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ssprfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    SSPRFS = (Fptr_NL_LAPACKE_ssprfs)dlsym(ssprfs_obj->hModule, "LAPACKE_ssprfs");
    ASSERT_TRUE(SSPRFS != NULL) << "failed to syt the Netlib LAPACKE_ssprfs symbol";
    
    SSPTRS = (Fptr_NL_LAPACKE_ssptrs)dlsym(ssprfs_obj->hModule, "LAPACKE_ssptrs");
    ASSERT_TRUE(SSPTRS != NULL) << "failed to syt the Netlib LAPACKE_ssptrs symbol";

	SSPTRF = (Fptr_NL_LAPACKE_ssptrf)dlsym(ssprfs_obj->hModule,"LAPACKE_ssptrf");
	ASSERT_TRUE(SSPTRF != NULL) << "failed to syt the Netlib LAPACKE_ssptrf symbol";

    /* Generate i/ps to syrfs API from 'sytrf' & 'sytrs' API sequence */
	ssprfs_obj->inforef = SSPTRF( ssprfs_obj->matrix_layout,
								  ssprfs_obj->uplo, ssprfs_obj->n,
								  ssprfs_obj->afref,
								  ssprfs_obj->ipivref);
						   
	ssprfs_obj->info = LAPACKE_ssptrf( ssprfs_obj->matrix_layout,
									   ssprfs_obj->uplo, ssprfs_obj->n,
									   ssprfs_obj->af,
									   ssprfs_obj->ipiv);

	ssprfs_obj->inforef = SSPTRS( ssprfs_obj->matrix_layout,
								  ssprfs_obj->uplo,
								  ssprfs_obj->n,
								  ssprfs_obj->nrhs,
								  ssprfs_obj->afref,
								  ssprfs_obj->ipivref,
                                  ssprfs_obj->bref,
								  ssprfs_obj->ldb
								  );
						   
	ssprfs_obj->info = LAPACKE_ssptrs( ssprfs_obj->matrix_layout,
									   ssprfs_obj->uplo,
									   ssprfs_obj->n,
									   ssprfs_obj->nrhs,
									   ssprfs_obj->af,
									   ssprfs_obj->ipiv,
									   ssprfs_obj->b,
									   ssprfs_obj->ldb
								    );
    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    ssprfs_obj->inforef = SSPRFS( ssprfs_obj->matrix_layout,
                                  ssprfs_obj->uplo, ssprfs_obj->n,
                                  ssprfs_obj->nrhs,
                                  ssprfs_obj->aref,
                                  ssprfs_obj->afref,
                                  ssprfs_obj->ipivref,
                                  ssprfs_obj->bref, ssprfs_obj->ldb,
                                  ssprfs_obj->xref, ssprfs_obj->ldx,
                                  ssprfs_obj->ferrref,
                                  ssprfs_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    ssprfs_obj->info = LAPACKE_ssprfs( ssprfs_obj->matrix_layout,
                                  ssprfs_obj->uplo, ssprfs_obj->n,
                                  ssprfs_obj->nrhs,
                                  ssprfs_obj->a,
                                  ssprfs_obj->af,
                                  ssprfs_obj->ipiv,
                                  ssprfs_obj->b, ssprfs_obj->ldb,
                                  ssprfs_obj->x, ssprfs_obj->ldx,
                                  ssprfs_obj->ferr,
                                  ssprfs_obj->berr);

    if( ssprfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ssprfs is wrong\n", ssprfs_obj->info );
    }
    if( ssprfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ssprfs is wrong\n", 
        ssprfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    ssprfs_obj->ipiv_diff = computeDiff_i( ssprfs_obj->n, ssprfs_obj->ipiv, ssprfs_obj->ipivref );
    
    ssprfs_obj->diff_xerr =  computeDiff_s( ssprfs_obj->x_bufsize, 
                ssprfs_obj->x, ssprfs_obj->xref );

    ssprfs_obj->diff_berr =  computeDiff_s( ssprfs_obj->nrhs, 
                ssprfs_obj->berr, ssprfs_obj->berrref );
                
    ssprfs_obj->diff_ferr =  computeDiff_s( ssprfs_obj->nrhs, 
                ssprfs_obj->ferr, ssprfs_obj->ferrref );
}

TEST_F(ssprfs_test, ssprfs1) {
    EXPECT_NEAR(0.0, ssprfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, ssprfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, ssprfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, ssprfs_obj->ipiv_diff );
    idx = Circular_Increment_Index(idx);
}

TEST_F(ssprfs_test, ssprfs2) {
    EXPECT_NEAR(0.0, ssprfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, ssprfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, ssprfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, ssprfs_obj->ipiv_diff );
    idx = Circular_Increment_Index(idx);
}

TEST_F(ssprfs_test, ssprfs3) {
    EXPECT_NEAR(0.0, ssprfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, ssprfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, ssprfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, ssprfs_obj->ipiv_diff );
    idx = Circular_Increment_Index(idx);
}

TEST_F(ssprfs_test, ssprfs4) {
    EXPECT_NEAR(0.0, ssprfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, ssprfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, ssprfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, ssprfs_obj->ipiv_diff );
    idx = Circular_Increment_Index(idx);
}

/* Begin sprfs_double_parameters  class definition */
class sprfs_double_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      int a_bufsize;
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
      double *af, *afref; //contains the factored form of the matrix A
      
      /* Output parameters */
      double* x, *xref;
      double* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      double* berr, *berrref; // component-wise relative backward error for each solution vector xj.
      lapack_int *ipiv, *ipivref; // pivot buffer

      /* Return Values */
      lapack_int info, inforef;

   public: 
      sprfs_double_parameters ( int matrix_layout_i, char uplo_i,
                               lapack_int n_i, lapack_int nrhs_i);
              
      ~sprfs_double_parameters (); 
};  /* end of sprfs_double_parameters  class definition */


/* Constructor sprfs_double_parameters definition */
sprfs_double_parameters:: sprfs_double_parameters ( int matrix_layout_i, 
					char uplo_i, lapack_int n_i, lapack_int nrhs_i) {
                                    
    matrix_layout = matrix_layout_i;
    uplo = uplo_i;
    // equed acts as input when " fact = 'F' " else output
    n = n_i;
    nrhs = nrhs_i;
    lda = n;
    ldaf = n;
    
    a_bufsize = n*(n+1)/2;
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
   printf(" \n sprfs double:  n: %d, uplo: %c  lda: %d  \
ldb: %d nrhs: %d  ldaf: %d  ldx: %d \n",  n, uplo, lda, 
                                          ldb, nrhs, ldaf, ldx);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, a_bufsize);
    lapacke_gtest_alloc_double_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_double_buffer_pair( &x, &xref, x_bufsize);
    lapacke_gtest_alloc_double_buffer_pair( &af, &afref, a_bufsize);
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
       sprfs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand ( a, aref, a_bufsize);
    memcpy(af, a, (a_bufsize*sizeof(double)));
    memcpy(afref, a, (a_bufsize*sizeof(double)));
    
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, b_bufsize);
    memcpy(x, b, ( b_bufsize*sizeof(double)));
    memcpy(xref, b, ( b_bufsize*sizeof(double)));
    
    lapacke_gtest_init_double_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);

    
   } /* end of Constructor  */

sprfs_double_parameters:: ~sprfs_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sprfs_double_parameters object: destructor invoked. \n");
#endif
   sprfs_free();
}

//  Test fixture class definition
class dsprfs_test  : public  ::testing::Test {
public:
   sprfs_double_parameters  *dsprfs_obj;
   void SetUp();  
   void TearDown () { delete dsprfs_obj; }
};


void dsprfs_test::SetUp(){

     /* LAPACKE DSPTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_dsptrf) ( int matrix_layout , char uplo,
                            lapack_int n , double * ap , lapack_int * ipiv );

    Fptr_NL_LAPACKE_dsptrf DSPTRF;

    typedef int (*Fptr_NL_LAPACKE_dsptrs) ( int matrix_layout , char uplo, 
					   lapack_int n , lapack_int nrhs , const double * ap , 
					const lapack_int * ipiv , double * b , lapack_int ldb );

    Fptr_NL_LAPACKE_dsptrs DSPTRS;

    /* LAPACKE DSPRFS prototype */
    typedef int (*Fptr_NL_LAPACKE_dsprfs) (int matrix_layout,
		char uplo, lapack_int n, lapack_int nrhs, const double* ap,
		const double* afp, const lapack_int* ipiv, const double* b,
		lapack_int ldb, double* x, lapack_int ldx, double* ferr,
		double* berr);

    Fptr_NL_LAPACKE_dsprfs DSPRFS;

    dsprfs_obj = new sprfs_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    dsprfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dsprfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dsprfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dsprfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    DSPRFS = (Fptr_NL_LAPACKE_dsprfs)dlsym(dsprfs_obj->hModule, "LAPACKE_dsprfs");
    ASSERT_TRUE(DSPRFS != NULL) << "failed to syt the Netlib LAPACKE_dsprfs symbol";
    
    DSPTRS = (Fptr_NL_LAPACKE_dsptrs)dlsym(dsprfs_obj->hModule, "LAPACKE_dsptrs");
    ASSERT_TRUE(DSPTRS != NULL) << "failed to syt the Netlib LAPACKE_dsptrs symbol";

	DSPTRF = (Fptr_NL_LAPACKE_dsptrf)dlsym(dsprfs_obj->hModule,"LAPACKE_dsptrf");
	ASSERT_TRUE(DSPTRF != NULL) << "failed to syt the Netlib LAPACKE_dsptrf symbol";

    /* Generate i/ps to syrfs API from 'sytrf' & 'sytrs' API sequence */
	dsprfs_obj->inforef = DSPTRF( dsprfs_obj->matrix_layout,
								  dsprfs_obj->uplo, dsprfs_obj->n,
								  dsprfs_obj->afref,
								  dsprfs_obj->ipivref);
						   
	dsprfs_obj->info = LAPACKE_dsptrf( dsprfs_obj->matrix_layout,
									   dsprfs_obj->uplo, dsprfs_obj->n,
									   dsprfs_obj->af,
									   dsprfs_obj->ipiv);

	dsprfs_obj->inforef = DSPTRS( dsprfs_obj->matrix_layout,
								  dsprfs_obj->uplo,
								  dsprfs_obj->n,
								  dsprfs_obj->nrhs,
								  dsprfs_obj->afref,
								  dsprfs_obj->ipivref,
                                  dsprfs_obj->bref,
								  dsprfs_obj->ldb
								  );
						   
	dsprfs_obj->info = LAPACKE_dsptrs( dsprfs_obj->matrix_layout,
									   dsprfs_obj->uplo,
									   dsprfs_obj->n,
									   dsprfs_obj->nrhs,
									   dsprfs_obj->af,
									   dsprfs_obj->ipiv,
									   dsprfs_obj->b,
									   dsprfs_obj->ldb
								    );
    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    dsprfs_obj->inforef = DSPRFS( dsprfs_obj->matrix_layout,
                                  dsprfs_obj->uplo, dsprfs_obj->n,
                                  dsprfs_obj->nrhs,
                                  dsprfs_obj->aref,
                                  dsprfs_obj->afref,
                                  dsprfs_obj->ipivref,
                                  dsprfs_obj->bref, dsprfs_obj->ldb,
                                  dsprfs_obj->xref, dsprfs_obj->ldx,
                                  dsprfs_obj->ferrref,
                                  dsprfs_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    dsprfs_obj->info = LAPACKE_dsprfs( dsprfs_obj->matrix_layout,
                                  dsprfs_obj->uplo, dsprfs_obj->n,
                                  dsprfs_obj->nrhs,
                                  dsprfs_obj->a,
                                  dsprfs_obj->af,
                                  dsprfs_obj->ipiv,
                                  dsprfs_obj->b, dsprfs_obj->ldb,
                                  dsprfs_obj->x, dsprfs_obj->ldx,
                                  dsprfs_obj->ferr,
                                  dsprfs_obj->berr);

    if( dsprfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dsprfs is wrong\n", dsprfs_obj->info );
    }
    if( dsprfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dsprfs is wrong\n", 
        dsprfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dsprfs_obj->ipiv_diff = computeDiff_i( dsprfs_obj->n, dsprfs_obj->ipiv, dsprfs_obj->ipivref );
    
    dsprfs_obj->diff_xerr =  computeDiff_d( dsprfs_obj->x_bufsize, 
                dsprfs_obj->x, dsprfs_obj->xref );

    dsprfs_obj->diff_berr =  computeDiff_d( dsprfs_obj->nrhs, 
                dsprfs_obj->berr, dsprfs_obj->berrref );
                
    dsprfs_obj->diff_ferr =  computeDiff_d( dsprfs_obj->nrhs, 
                dsprfs_obj->ferr, dsprfs_obj->ferrref );
}

TEST_F(dsprfs_test, dsprfs1) {
    EXPECT_NEAR(0.0, dsprfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dsprfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dsprfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, dsprfs_obj->ipiv_diff );
    idx = Circular_Increment_Index(idx);
}

TEST_F(dsprfs_test, dsprfs2) {
    EXPECT_NEAR(0.0, dsprfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dsprfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dsprfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, dsprfs_obj->ipiv_diff );
    idx = Circular_Increment_Index(idx);
}

TEST_F(dsprfs_test, dsprfs3) {
    EXPECT_NEAR(0.0, dsprfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dsprfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dsprfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, dsprfs_obj->ipiv_diff );
    idx = Circular_Increment_Index(idx);
}

TEST_F(dsprfs_test, dsprfs4) {
    EXPECT_NEAR(0.0, dsprfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dsprfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dsprfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, dsprfs_obj->ipiv_diff );
    idx = Circular_Increment_Index(idx);
}

/* Begin sprfs_scomplex_parameters  class definition */
class sprfs_scomplex_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      int a_bufsize;
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
      lapack_complex_float *af, *afref; //contains the factored form of the matrix A
      
      /* Output parameters */
      lapack_complex_float* x, *xref;
      float* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      float* berr, *berrref; // component-wise relative backward error for each solution vector xj.
      lapack_int *ipiv, *ipivref; // pivot buffer

      /* Return Values */
      lapack_int info, inforef;

   public: 
      sprfs_scomplex_parameters ( int matrix_layout_i, char uplo_i,
                               lapack_int n_i, lapack_int nrhs_i);
              
      ~sprfs_scomplex_parameters (); 
};  /* end of sprfs_scomplex_parameters  class definition */


/* Constructor sprfs_scomplex_parameters definition */
sprfs_scomplex_parameters:: sprfs_scomplex_parameters ( int matrix_layout_i, 
					char uplo_i, lapack_int n_i, lapack_int nrhs_i) {
                                    
    matrix_layout = matrix_layout_i;
    uplo = uplo_i;
    // equed acts as input when " fact = 'F' " else output
    n = n_i;
    nrhs = nrhs_i;
    lda = n;
    ldaf = n;
    
    a_bufsize = n*(n+1)/2;
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
   printf(" \n sprfs lapack_complex_float:  n: %d, uplo: %c  lda: %d  \
ldb: %d nrhs: %d  ldaf: %d  ldx: %d \n",  n, uplo, lda, 
                                          ldb, nrhs, ldaf, ldx);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, a_bufsize);
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &x, &xref, x_bufsize);
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &af, &afref, a_bufsize);
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
       sprfs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand ( a, aref, a_bufsize);
    memcpy(af, a, (a_bufsize*sizeof(lapack_complex_float)));
    memcpy(afref, a, (a_bufsize*sizeof(lapack_complex_float)));
    
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
    memcpy(x, b, ( b_bufsize*sizeof(lapack_complex_float)));
    memcpy(xref, b, ( b_bufsize*sizeof(lapack_complex_float)));
    
    lapacke_gtest_init_float_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);

    
   } /* end of Constructor  */

sprfs_scomplex_parameters:: ~sprfs_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sprfs_scomplex_parameters object: destructor invoked. \n");
#endif
   sprfs_free();
}

//  Test fixture class definition
class csprfs_test  : public  ::testing::Test {
public:
   sprfs_scomplex_parameters  *csprfs_obj;
   void SetUp();  
   void TearDown () { delete csprfs_obj; }
};


void csprfs_test::SetUp(){

     /* LAPACKE CSPTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_csptrf) ( int matrix_layout , char uplo,
                            lapack_int n , lapack_complex_float * ap , lapack_int * ipiv );

    Fptr_NL_LAPACKE_csptrf CSPTRF;

    typedef int (*Fptr_NL_LAPACKE_csptrs) ( int matrix_layout , char uplo, 
					   lapack_int n , lapack_int nrhs , const lapack_complex_float * ap , 
					const lapack_int * ipiv , lapack_complex_float * b , lapack_int ldb );

    Fptr_NL_LAPACKE_csptrs CSPTRS;

    /* LAPACKE CSPRFS prototype */
    typedef int (*Fptr_NL_LAPACKE_csprfs) (int matrix_layout,
		char uplo, lapack_int n, lapack_int nrhs, const lapack_complex_float* ap,
		const lapack_complex_float* afp, const lapack_int* ipiv, const lapack_complex_float* b,
		lapack_int ldb, lapack_complex_float* x, lapack_int ldx, float* ferr,
		float* berr);

    Fptr_NL_LAPACKE_csprfs CSPRFS;

    csprfs_obj = new sprfs_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    csprfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    csprfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(csprfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(csprfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    CSPRFS = (Fptr_NL_LAPACKE_csprfs)dlsym(csprfs_obj->hModule, "LAPACKE_csprfs");
    ASSERT_TRUE(CSPRFS != NULL) << "failed to syt the Netlib LAPACKE_csprfs symbol";
    
    CSPTRS = (Fptr_NL_LAPACKE_csptrs)dlsym(csprfs_obj->hModule, "LAPACKE_csptrs");
    ASSERT_TRUE(CSPTRS != NULL) << "failed to syt the Netlib LAPACKE_csptrs symbol";

	CSPTRF = (Fptr_NL_LAPACKE_csptrf)dlsym(csprfs_obj->hModule,"LAPACKE_csptrf");
	ASSERT_TRUE(CSPTRF != NULL) << "failed to syt the Netlib LAPACKE_csptrf symbol";

    /* Generate i/ps to syrfs API from 'sytrf' & 'sytrs' API sequence */
	csprfs_obj->inforef = CSPTRF( csprfs_obj->matrix_layout,
								  csprfs_obj->uplo, csprfs_obj->n,
								  csprfs_obj->afref,
								  csprfs_obj->ipivref);
						   
	csprfs_obj->info = LAPACKE_csptrf( csprfs_obj->matrix_layout,
									   csprfs_obj->uplo, csprfs_obj->n,
									   csprfs_obj->af,
									   csprfs_obj->ipiv);

	csprfs_obj->inforef = CSPTRS( csprfs_obj->matrix_layout,
								  csprfs_obj->uplo,
								  csprfs_obj->n,
								  csprfs_obj->nrhs,
								  csprfs_obj->afref,
								  csprfs_obj->ipivref,
                                  csprfs_obj->bref,
								  csprfs_obj->ldb
								  );
						   
	csprfs_obj->info = LAPACKE_csptrs( csprfs_obj->matrix_layout,
									   csprfs_obj->uplo,
									   csprfs_obj->n,
									   csprfs_obj->nrhs,
									   csprfs_obj->af,
									   csprfs_obj->ipiv,
									   csprfs_obj->b,
									   csprfs_obj->ldb
								    );
    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    csprfs_obj->inforef = CSPRFS( csprfs_obj->matrix_layout,
                                  csprfs_obj->uplo, csprfs_obj->n,
                                  csprfs_obj->nrhs,
                                  csprfs_obj->aref,
                                  csprfs_obj->afref,
                                  csprfs_obj->ipivref,
                                  csprfs_obj->bref, csprfs_obj->ldb,
                                  csprfs_obj->xref, csprfs_obj->ldx,
                                  csprfs_obj->ferrref,
                                  csprfs_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    csprfs_obj->info = LAPACKE_csprfs( csprfs_obj->matrix_layout,
                                  csprfs_obj->uplo, csprfs_obj->n,
                                  csprfs_obj->nrhs,
                                  csprfs_obj->a,
                                  csprfs_obj->af,
                                  csprfs_obj->ipiv,
                                  csprfs_obj->b, csprfs_obj->ldb,
                                  csprfs_obj->x, csprfs_obj->ldx,
                                  csprfs_obj->ferr,
                                  csprfs_obj->berr);

    if( csprfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_csprfs is wrong\n", csprfs_obj->info );
    }
    if( csprfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_csprfs is wrong\n", 
        csprfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    csprfs_obj->ipiv_diff = computeDiff_i( csprfs_obj->n, csprfs_obj->ipiv, csprfs_obj->ipivref );
    
    csprfs_obj->diff_xerr =  computeDiff_c( csprfs_obj->x_bufsize, 
                csprfs_obj->x, csprfs_obj->xref );

    csprfs_obj->diff_berr =  computeDiff_s( csprfs_obj->nrhs, 
                csprfs_obj->berr, csprfs_obj->berrref );
                
    csprfs_obj->diff_ferr =  computeDiff_s( csprfs_obj->nrhs, 
                csprfs_obj->ferr, csprfs_obj->ferrref );
}

TEST_F(csprfs_test, csprfs1) {
    EXPECT_NEAR(0.0, csprfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, csprfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, csprfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, csprfs_obj->ipiv_diff );
    idx = Circular_Increment_Index(idx);
}

TEST_F(csprfs_test, csprfs2) {
    EXPECT_NEAR(0.0, csprfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, csprfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, csprfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, csprfs_obj->ipiv_diff );
    idx = Circular_Increment_Index(idx);
}

TEST_F(csprfs_test, csprfs3) {
    EXPECT_NEAR(0.0, csprfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, csprfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, csprfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, csprfs_obj->ipiv_diff );
    idx = Circular_Increment_Index(idx);
}

TEST_F(csprfs_test, csprfs4) {
    EXPECT_NEAR(0.0, csprfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, csprfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, csprfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, csprfs_obj->ipiv_diff );
    idx = Circular_Increment_Index(idx);
}

/* Begin sprfs_dcomplex_parameters  class definition */
class sprfs_dcomplex_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      int a_bufsize;
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
      lapack_complex_double *af, *afref; //contains the factored form of the matrix A
      
      /* Output parameters */
      lapack_complex_double* x, *xref;
      double* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      double* berr, *berrref; // component-wise relative backward error for each solution vector xj.
      lapack_int *ipiv, *ipivref; // pivot buffer

      /* Return Values */
      lapack_int info, inforef;

   public: 
      sprfs_dcomplex_parameters ( int matrix_layout_i, char uplo_i,
                               lapack_int n_i, lapack_int nrhs_i);
              
      ~sprfs_dcomplex_parameters (); 
};  /* end of sprfs_dcomplex_parameters  class definition */


/* Constructor sprfs_dcomplex_parameters definition */
sprfs_dcomplex_parameters:: sprfs_dcomplex_parameters ( int matrix_layout_i, 
					char uplo_i, lapack_int n_i, lapack_int nrhs_i) {
                                    
    matrix_layout = matrix_layout_i;
    uplo = uplo_i;
    // equed acts as input when " fact = 'F' " else output
    n = n_i;
    nrhs = nrhs_i;
    lda = n;
    ldaf = n;
    
    a_bufsize = n*(n+1)/2;
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
   printf(" \n sprfs lapack_complex_double:  n: %d, uplo: %c  lda: %d  \
ldb: %d nrhs: %d  ldaf: %d  ldx: %d \n",  n, uplo, lda, 
                                          ldb, nrhs, ldaf, ldx);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, a_bufsize);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &x, &xref, x_bufsize);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &af, &afref, a_bufsize);
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
       sprfs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand ( a, aref, a_bufsize);
    memcpy(af, a, (a_bufsize*sizeof(lapack_complex_double)));
    memcpy(afref, a, (a_bufsize*sizeof(lapack_complex_double)));
    
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
    memcpy(x, b, ( b_bufsize*sizeof(lapack_complex_double)));
    memcpy(xref, b, ( b_bufsize*sizeof(lapack_complex_double)));
    
    lapacke_gtest_init_double_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);

    
   } /* end of Constructor  */

sprfs_dcomplex_parameters:: ~sprfs_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sprfs_dcomplex_parameters object: destructor invoked. \n");
#endif
   sprfs_free();
}

//  Test fixture class definition
class zsprfs_test  : public  ::testing::Test {
public:
   sprfs_dcomplex_parameters  *zsprfs_obj;
   void SetUp();  
   void TearDown () { delete zsprfs_obj; }
};


void zsprfs_test::SetUp(){

     /* LAPACKE ZSPTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zsptrf) ( int matrix_layout , char uplo,
                            lapack_int n , lapack_complex_double * ap , lapack_int * ipiv );

    Fptr_NL_LAPACKE_zsptrf ZSPTRF;

    typedef int (*Fptr_NL_LAPACKE_zsptrs) ( int matrix_layout , char uplo, 
					   lapack_int n , lapack_int nrhs , const lapack_complex_double * ap , 
					const lapack_int * ipiv , lapack_complex_double * b , lapack_int ldb );

    Fptr_NL_LAPACKE_zsptrs ZSPTRS;

    /* LAPACKE ZSPRFS prototype */
    typedef int (*Fptr_NL_LAPACKE_zsprfs) (int matrix_layout,
		char uplo, lapack_int n, lapack_int nrhs, const lapack_complex_double* ap,
		const lapack_complex_double* afp, const lapack_int* ipiv, const lapack_complex_double* b,
		lapack_int ldb, lapack_complex_double* x, lapack_int ldx, double* ferr,
		double* berr);

    Fptr_NL_LAPACKE_zsprfs ZSPRFS;

    zsprfs_obj = new sprfs_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    zsprfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zsprfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zsprfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zsprfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    ZSPRFS = (Fptr_NL_LAPACKE_zsprfs)dlsym(zsprfs_obj->hModule, "LAPACKE_zsprfs");
    ASSERT_TRUE(ZSPRFS != NULL) << "failed to syt the Netlib LAPACKE_zsprfs symbol";
    
    ZSPTRS = (Fptr_NL_LAPACKE_zsptrs)dlsym(zsprfs_obj->hModule, "LAPACKE_zsptrs");
    ASSERT_TRUE(ZSPTRS != NULL) << "failed to syt the Netlib LAPACKE_zsptrs symbol";

	ZSPTRF = (Fptr_NL_LAPACKE_zsptrf)dlsym(zsprfs_obj->hModule,"LAPACKE_zsptrf");
	ASSERT_TRUE(ZSPTRF != NULL) << "failed to syt the Netlib LAPACKE_zsptrf symbol";

    /* Generate i/ps to syrfs API from 'sytrf' & 'sytrs' API sequence */
	zsprfs_obj->inforef = ZSPTRF( zsprfs_obj->matrix_layout,
								  zsprfs_obj->uplo, zsprfs_obj->n,
								  zsprfs_obj->afref,
								  zsprfs_obj->ipivref);
						   
	zsprfs_obj->info = LAPACKE_zsptrf( zsprfs_obj->matrix_layout,
									   zsprfs_obj->uplo, zsprfs_obj->n,
									   zsprfs_obj->af,
									   zsprfs_obj->ipiv);

	zsprfs_obj->inforef = ZSPTRS( zsprfs_obj->matrix_layout,
								  zsprfs_obj->uplo,
								  zsprfs_obj->n,
								  zsprfs_obj->nrhs,
								  zsprfs_obj->afref,
								  zsprfs_obj->ipivref,
                                  zsprfs_obj->bref,
								  zsprfs_obj->ldb
								  );
						   
	zsprfs_obj->info = LAPACKE_zsptrs( zsprfs_obj->matrix_layout,
									   zsprfs_obj->uplo,
									   zsprfs_obj->n,
									   zsprfs_obj->nrhs,
									   zsprfs_obj->af,
									   zsprfs_obj->ipiv,
									   zsprfs_obj->b,
									   zsprfs_obj->ldb
								    );
    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    zsprfs_obj->inforef = ZSPRFS( zsprfs_obj->matrix_layout,
                                  zsprfs_obj->uplo, zsprfs_obj->n,
                                  zsprfs_obj->nrhs,
                                  zsprfs_obj->aref,
                                  zsprfs_obj->afref,
                                  zsprfs_obj->ipivref,
                                  zsprfs_obj->bref, zsprfs_obj->ldb,
                                  zsprfs_obj->xref, zsprfs_obj->ldx,
                                  zsprfs_obj->ferrref,
                                  zsprfs_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    zsprfs_obj->info = LAPACKE_zsprfs( zsprfs_obj->matrix_layout,
                                  zsprfs_obj->uplo, zsprfs_obj->n,
                                  zsprfs_obj->nrhs,
                                  zsprfs_obj->a,
                                  zsprfs_obj->af,
                                  zsprfs_obj->ipiv,
                                  zsprfs_obj->b, zsprfs_obj->ldb,
                                  zsprfs_obj->x, zsprfs_obj->ldx,
                                  zsprfs_obj->ferr,
                                  zsprfs_obj->berr);

    if( zsprfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zsprfs is wrong\n", zsprfs_obj->info );
    }
    if( zsprfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zsprfs is wrong\n", 
        zsprfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zsprfs_obj->ipiv_diff = computeDiff_i( zsprfs_obj->n, zsprfs_obj->ipiv, zsprfs_obj->ipivref );
    
    zsprfs_obj->diff_xerr =  computeDiff_z( zsprfs_obj->x_bufsize, 
                zsprfs_obj->x, zsprfs_obj->xref );

    zsprfs_obj->diff_berr =  computeDiff_d( zsprfs_obj->nrhs, 
                zsprfs_obj->berr, zsprfs_obj->berrref );
                
    zsprfs_obj->diff_ferr =  computeDiff_d( zsprfs_obj->nrhs, 
                zsprfs_obj->ferr, zsprfs_obj->ferrref );
}

TEST_F(zsprfs_test, zsprfs1) {
    EXPECT_NEAR(0.0, zsprfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zsprfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zsprfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, zsprfs_obj->ipiv_diff );
    idx = Circular_Increment_Index(idx);
}

TEST_F(zsprfs_test, zsprfs2) {
    EXPECT_NEAR(0.0, zsprfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zsprfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zsprfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, zsprfs_obj->ipiv_diff );
    idx = Circular_Increment_Index(idx);
}

TEST_F(zsprfs_test, zsprfs3) {
    EXPECT_NEAR(0.0, zsprfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zsprfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zsprfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, zsprfs_obj->ipiv_diff );
    idx = Circular_Increment_Index(idx);
}

TEST_F(zsprfs_test, zsprfs4) {
    EXPECT_NEAR(0.0, zsprfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zsprfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zsprfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, zsprfs_obj->ipiv_diff );
    idx = Circular_Increment_Index(idx);
}
