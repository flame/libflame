#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define herfs_free() \
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

/* Begin herfs_scomplex_parameters  class definition */
class herfs_scomplex_parameters{
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
      lapack_complex_float *af, *afref; //contains the factored form of the matrix A
      
      /* Output parameters */
      lapack_complex_float* x, *xref;
      float* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      float* berr, *berrref; // component-wise relative backward error for each solution vector xj.
      lapack_int *ipiv, *ipivref; // pivot buffer

      /* Return Values */
      lapack_int info, inforef;

   public: 
      herfs_scomplex_parameters ( int matrix_layout_i, char uplo_i,
                               lapack_int n_i, lapack_int nrhs_i);
              
      ~herfs_scomplex_parameters (); 
};  /* end of herfs_scomplex_parameters  class definition */


/* Constructor herfs_scomplex_parameters definition */
herfs_scomplex_parameters:: herfs_scomplex_parameters ( int matrix_layout_i, 
						  char uplo_i, lapack_int n_i, lapack_int nrhs_i) {
                                    
    matrix_layout = matrix_layout_i;
    uplo = uplo_i;
    // equed acts as input when " fact = 'F' " else output
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
   printf(" \n herfs lapack_complex_float:  n: %d, fact: %c uplo: %c  lda: %d  \
ldb: %d nrhs: %d  ldaf: %d  ldx: %d \n",  n, fact, uplo, lda, 
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
       herfs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix ( a, aref, n, n, 'H');
    memcpy(af, a, (n*n*sizeof(lapack_complex_float)));
    memcpy(afref, a, (n*n*sizeof(lapack_complex_float)));
    
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
    memcpy(x, b, ( b_bufsize*sizeof(lapack_complex_float)));
    memcpy(xref, b, ( b_bufsize*sizeof(lapack_complex_float)));
    
    lapacke_gtest_init_float_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);

    
   } /* end of Constructor  */

herfs_scomplex_parameters:: ~herfs_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" herfs_scomplex_parameters object: destructor invoked. \n");
#endif
   herfs_free();
}

//  Test fixture class definition
class cherfs_test  : public  ::testing::Test {
public:
   herfs_scomplex_parameters  *cherfs_obj;
   void SetUp();  
   void TearDown () { delete cherfs_obj; }
};


void cherfs_test::SetUp(){

     /* LAPACKE CHETRF prototype */
    typedef int (*Fptr_NL_LAPACKE_chetrf) ( int matrix_layout,
			char uplo, lapack_int n , lapack_complex_float * a,
			lapack_int lda , lapack_int * ipiv );

    Fptr_NL_LAPACKE_chetrf CHETRF;

    typedef int (*Fptr_NL_LAPACKE_chetrs) ( int matrix_layout,
		char uplo , lapack_int n , lapack_int nrhs ,
		const lapack_complex_float * a , lapack_int lda ,
		const lapack_int * ipiv , lapack_complex_float * b ,
		lapack_int ldb );

    Fptr_NL_LAPACKE_chetrs CHETRS;

    /* LAPACKE CHERFS prototype */
    typedef int (*Fptr_NL_LAPACKE_cherfs) (int matrix_layout,
		char uplo, lapack_int n, lapack_int nrhs,
		const lapack_complex_float* a, lapack_int lda,
		const lapack_complex_float* af, lapack_int ldaf,
		const lapack_int* ipiv, const lapack_complex_float* b,
		lapack_int ldb, lapack_complex_float* x,
		lapack_int ldx, float* ferr, float* berr);

    Fptr_NL_LAPACKE_cherfs CHERFS;

    cherfs_obj = new herfs_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);


    cherfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cherfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cherfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cherfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    CHERFS = (Fptr_NL_LAPACKE_cherfs)dlsym(cherfs_obj->hModule, "LAPACKE_cherfs");
    ASSERT_TRUE(CHERFS != NULL) << "failed to syt the Netlib LAPACKE_cherfs symbol";
    
	CHETRF = (Fptr_NL_LAPACKE_chetrf)dlsym(cherfs_obj->hModule,"LAPACKE_chetrf");
	ASSERT_TRUE(CHETRF != NULL) << "failed to syt the Netlib LAPACKE_chetrf symbol";
		
	CHETRS = (Fptr_NL_LAPACKE_chetrs)dlsym(cherfs_obj->hModule,"LAPACKE_chetrs");
	ASSERT_TRUE(CHETRS != NULL) << "failed to syt the Netlib LAPACKE_chetrs symbol";

    /* Generate i/ps to herfs API from 'hetrf' & 'hetrs' API sequence */
	cherfs_obj->inforef = CHETRF( cherfs_obj->matrix_layout,
								  cherfs_obj->uplo, cherfs_obj->n,
								  cherfs_obj->afref,
								  cherfs_obj->lda,
								  cherfs_obj->ipivref);

	cherfs_obj->info = LAPACKE_chetrf( cherfs_obj->matrix_layout,
									   cherfs_obj->uplo, cherfs_obj->n,
									   cherfs_obj->af,
									   cherfs_obj->lda,
									   cherfs_obj->ipiv);
		
	cherfs_obj->inforef = CHETRS( cherfs_obj->matrix_layout,
								  cherfs_obj->uplo, cherfs_obj->n,
                                  cherfs_obj->nrhs,
								  cherfs_obj->afref,
								  cherfs_obj->lda,
								  cherfs_obj->ipivref,
								  cherfs_obj->xref,
								  cherfs_obj->ldx
								  );

	cherfs_obj->info = LAPACKE_chetrs( cherfs_obj->matrix_layout,
								  cherfs_obj->uplo, cherfs_obj->n,
                                  cherfs_obj->nrhs,
								  cherfs_obj->af,
								  cherfs_obj->lda,
								  cherfs_obj->ipiv,
								  cherfs_obj->x,
								  cherfs_obj->ldx
								  );

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    cherfs_obj->inforef = CHERFS( cherfs_obj->matrix_layout,
                                  cherfs_obj->uplo, cherfs_obj->n,
                                  cherfs_obj->nrhs,
                                  cherfs_obj->aref, cherfs_obj->lda, 
                                  cherfs_obj->afref, cherfs_obj->ldaf,
                                  cherfs_obj->ipivref,
                                  cherfs_obj->bref, cherfs_obj->ldb,
                                  cherfs_obj->xref, cherfs_obj->ldx,
                                  cherfs_obj->ferrref,
                                  cherfs_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    cherfs_obj->info = LAPACKE_cherfs( cherfs_obj->matrix_layout,
                                  cherfs_obj->uplo, cherfs_obj->n,
                                  cherfs_obj->nrhs,
                                  cherfs_obj->a, cherfs_obj->lda, 
                                  cherfs_obj->af, cherfs_obj->ldaf,
                                  cherfs_obj->ipiv,
                                  cherfs_obj->b, cherfs_obj->ldb,
                                  cherfs_obj->x, cherfs_obj->ldx,
                                  cherfs_obj->ferr,
                                  cherfs_obj->berr);

    if( cherfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cherfs is wrong\n", cherfs_obj->info );
    }
    if( cherfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cherfs is wrong\n", 
        cherfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cherfs_obj->ipiv_diff = computeDiff_i( cherfs_obj->n, cherfs_obj->ipiv, cherfs_obj->ipivref );
    
    cherfs_obj->diff_xerr =  computeDiff_c( cherfs_obj->x_bufsize, 
                cherfs_obj->x, cherfs_obj->xref );

    cherfs_obj->diff_berr =  computeDiff_s( cherfs_obj->nrhs, 
                cherfs_obj->berr, cherfs_obj->berrref );
                
    cherfs_obj->diff_ferr =  computeDiff_s( cherfs_obj->nrhs, 
                cherfs_obj->ferr, cherfs_obj->ferrref );
}

TEST_F(cherfs_test, cherfs1) {
    EXPECT_NEAR(0.0, cherfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, cherfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, cherfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, cherfs_obj->ipiv_diff );
    idx = Circular_Increment_Index(idx);
}

TEST_F(cherfs_test, cherfs2) {
    EXPECT_NEAR(0.0, cherfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, cherfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, cherfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, cherfs_obj->ipiv_diff );
    idx = Circular_Increment_Index(idx);
}

TEST_F(cherfs_test, cherfs3) {
    EXPECT_NEAR(0.0, cherfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, cherfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, cherfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, cherfs_obj->ipiv_diff );
    idx = Circular_Increment_Index(idx);
}

TEST_F(cherfs_test, cherfs4) {
    EXPECT_NEAR(0.0, cherfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, cherfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, cherfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, cherfs_obj->ipiv_diff );
    idx = Circular_Increment_Index(idx);
}

/* Begin herfs_dcomplex_parameters  class definition */
class herfs_dcomplex_parameters{
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
      lapack_complex_double *af, *afref; //contains the factored form of the matrix A
      
      /* Output parameters */
      lapack_complex_double* x, *xref;
      double* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      double* berr, *berrref; // component-wise relative backward error for each solution vector xj.
      lapack_int *ipiv, *ipivref; // pivot buffer

      /* Return Values */
      lapack_int info, inforef;

   public: 
      herfs_dcomplex_parameters ( int matrix_layout_i, char uplo_i,
                               lapack_int n_i, lapack_int nrhs_i);
              
      ~herfs_dcomplex_parameters (); 
};  /* end of herfs_dcomplex_parameters  class definition */


/* Constructor herfs_dcomplex_parameters definition */
herfs_dcomplex_parameters:: herfs_dcomplex_parameters ( int matrix_layout_i, 
						  char uplo_i, lapack_int n_i, lapack_int nrhs_i) {
                                    
    matrix_layout = matrix_layout_i;
    uplo = uplo_i;
    // equed acts as input when " fact = 'F' " else output
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
   printf(" \n herfs lapack_complex_double:  n: %d, fact: %c uplo: %c  lda: %d  \
ldb: %d nrhs: %d  ldaf: %d  ldx: %d \n",  n, fact, uplo, lda, 
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
       herfs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix ( a, aref, n, n, 'H');
    memcpy(af, a, (n*n*sizeof(lapack_complex_double)));
    memcpy(afref, a, (n*n*sizeof(lapack_complex_double)));
    
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
    memcpy(x, b, ( b_bufsize*sizeof(lapack_complex_double)));
    memcpy(xref, b, ( b_bufsize*sizeof(lapack_complex_double)));
    
    lapacke_gtest_init_double_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);

    
   } /* end of Constructor  */

herfs_dcomplex_parameters:: ~herfs_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" herfs_dcomplex_parameters object: destructor invoked. \n");
#endif
   herfs_free();
}

//  Test fixture class definition
class zherfs_test  : public  ::testing::Test {
public:
   herfs_dcomplex_parameters  *zherfs_obj;
   void SetUp();  
   void TearDown () { delete zherfs_obj; }
};


void zherfs_test::SetUp(){

     /* LAPACKE ZHETRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zhetrf) ( int matrix_layout,
			char uplo, lapack_int n , lapack_complex_double * a,
			lapack_int lda , lapack_int * ipiv );

    Fptr_NL_LAPACKE_zhetrf ZHETRF;

    typedef int (*Fptr_NL_LAPACKE_zhetrs) ( int matrix_layout,
		char uplo , lapack_int n , lapack_int nrhs ,
		const lapack_complex_double * a , lapack_int lda ,
		const lapack_int * ipiv , lapack_complex_double * b ,
		lapack_int ldb );

    Fptr_NL_LAPACKE_zhetrs ZHETRS;

    /* LAPACKE ZHERFS prototype */
    typedef int (*Fptr_NL_LAPACKE_zherfs) (int matrix_layout,
		char uplo, lapack_int n, lapack_int nrhs,
		const lapack_complex_double* a, lapack_int lda,
		const lapack_complex_double* af, lapack_int ldaf,
		const lapack_int* ipiv, const lapack_complex_double* b,
		lapack_int ldb, lapack_complex_double* x,
		lapack_int ldx, double* ferr, double* berr);

    Fptr_NL_LAPACKE_zherfs ZHERFS;

    zherfs_obj = new herfs_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);


    zherfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zherfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zherfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zherfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    ZHERFS = (Fptr_NL_LAPACKE_zherfs)dlsym(zherfs_obj->hModule, "LAPACKE_zherfs");
    ASSERT_TRUE(ZHERFS != NULL) << "failed to syt the Netlib LAPACKE_zherfs symbol";
    
	ZHETRF = (Fptr_NL_LAPACKE_zhetrf)dlsym(zherfs_obj->hModule,"LAPACKE_zhetrf");
	ASSERT_TRUE(ZHETRF != NULL) << "failed to syt the Netlib LAPACKE_zhetrf symbol";
		
	ZHETRS = (Fptr_NL_LAPACKE_zhetrs)dlsym(zherfs_obj->hModule,"LAPACKE_zhetrs");
	ASSERT_TRUE(ZHETRS != NULL) << "failed to syt the Netlib LAPACKE_zhetrs symbol";

    /* Generate i/ps to herfs API from 'hetrf' & 'hetrs' API sequence */
	zherfs_obj->inforef = ZHETRF( zherfs_obj->matrix_layout,
								  zherfs_obj->uplo, zherfs_obj->n,
								  zherfs_obj->afref,
								  zherfs_obj->lda,
								  zherfs_obj->ipivref);

	zherfs_obj->info = LAPACKE_zhetrf( zherfs_obj->matrix_layout,
									   zherfs_obj->uplo, zherfs_obj->n,
									   zherfs_obj->af,
									   zherfs_obj->lda,
									   zherfs_obj->ipiv);
		
	zherfs_obj->inforef = ZHETRS( zherfs_obj->matrix_layout,
								  zherfs_obj->uplo, zherfs_obj->n,
                                  zherfs_obj->nrhs,
								  zherfs_obj->afref,
								  zherfs_obj->lda,
								  zherfs_obj->ipivref,
								  zherfs_obj->xref,
								  zherfs_obj->ldx
								  );

	zherfs_obj->info = LAPACKE_zhetrs( zherfs_obj->matrix_layout,
								  zherfs_obj->uplo, zherfs_obj->n,
                                  zherfs_obj->nrhs,
								  zherfs_obj->af,
								  zherfs_obj->lda,
								  zherfs_obj->ipiv,
								  zherfs_obj->x,
								  zherfs_obj->ldx
								  );

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    zherfs_obj->inforef = ZHERFS( zherfs_obj->matrix_layout,
                                  zherfs_obj->uplo, zherfs_obj->n,
                                  zherfs_obj->nrhs,
                                  zherfs_obj->aref, zherfs_obj->lda, 
                                  zherfs_obj->afref, zherfs_obj->ldaf,
                                  zherfs_obj->ipivref,
                                  zherfs_obj->bref, zherfs_obj->ldb,
                                  zherfs_obj->xref, zherfs_obj->ldx,
                                  zherfs_obj->ferrref,
                                  zherfs_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    zherfs_obj->info = LAPACKE_zherfs( zherfs_obj->matrix_layout,
                                  zherfs_obj->uplo, zherfs_obj->n,
                                  zherfs_obj->nrhs,
                                  zherfs_obj->a, zherfs_obj->lda, 
                                  zherfs_obj->af, zherfs_obj->ldaf,
                                  zherfs_obj->ipiv,
                                  zherfs_obj->b, zherfs_obj->ldb,
                                  zherfs_obj->x, zherfs_obj->ldx,
                                  zherfs_obj->ferr,
                                  zherfs_obj->berr);

    if( zherfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zherfs is wrong\n", zherfs_obj->info );
    }
    if( zherfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zherfs is wrong\n", 
        zherfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zherfs_obj->ipiv_diff = computeDiff_i( zherfs_obj->n, zherfs_obj->ipiv, zherfs_obj->ipivref );
    
    zherfs_obj->diff_xerr =  computeDiff_z( zherfs_obj->x_bufsize, 
                zherfs_obj->x, zherfs_obj->xref );

    zherfs_obj->diff_berr =  computeDiff_d( zherfs_obj->nrhs, 
                zherfs_obj->berr, zherfs_obj->berrref );
                
    zherfs_obj->diff_ferr =  computeDiff_d( zherfs_obj->nrhs, 
                zherfs_obj->ferr, zherfs_obj->ferrref );
}

TEST_F(zherfs_test, zherfs1) {
    EXPECT_NEAR(0.0, zherfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zherfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zherfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, zherfs_obj->ipiv_diff );
    idx = Circular_Increment_Index(idx);
}

TEST_F(zherfs_test, zherfs2) {
    EXPECT_NEAR(0.0, zherfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zherfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zherfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, zherfs_obj->ipiv_diff );
    idx = Circular_Increment_Index(idx);
}

TEST_F(zherfs_test, zherfs3) {
    EXPECT_NEAR(0.0, zherfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zherfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zherfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, zherfs_obj->ipiv_diff );
    idx = Circular_Increment_Index(idx);
}

TEST_F(zherfs_test, zherfs4) {
    EXPECT_NEAR(0.0, zherfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zherfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zherfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, zherfs_obj->ipiv_diff );
    idx = Circular_Increment_Index(idx);
}


