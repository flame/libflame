#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define gerfs_free() \
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

/* Begin gerfs_float_parameters  class definition */
class gerfs_float_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      float diff; // capture difference between ref o/p & libflame lapacke o/p.
	  float diff_berr, diff_ferr, diff_xerr;
      void *hModule, *dModule;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char fact; //  Must be 'F', 'N', or 'E'.
      char trans; //  Must be 'N' , 'T' or 'C'.

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
      lapack_int *ipiv, *ipivref; // pivot buffer
      
      /* Output parameters */
      float rpivot, rpivotref; // reciprocal pivot growth factor.
      float rcond, rcondref; // estimate of the reciprocal condition number of the matrix A.
      float* x, *xref;
      float* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      float* berr, *berrref; // component-wise relative backward error for each solution vector xj.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      gerfs_float_parameters ( int matrix_layout_i, char trans_i,
                                lapack_int n_i, lapack_int nrhs_i);
              
      ~gerfs_float_parameters (); 
};  /* end of gerfs_float_parameters  class definition */


/* Constructor gerfs_float_parameters definition */
gerfs_float_parameters:: gerfs_float_parameters ( int matrix_layout_i, 
                 char trans_i,  lapack_int n_i, lapack_int nrhs_i) {
                                    
    matrix_layout = matrix_layout_i;
    trans = trans_i;
    n = n_i;
    nrhs = nrhs_i;
    lda = n;
	ldaf = n;

    if(matrix_layout==LAPACK_COL_MAJOR){
		ldb = n;
		ldx = n;
		
        b_bufsize = ldb*nrhs;
        x_bufsize = ldx*nrhs;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
		ldb = nrhs;
		ldx = nrhs;
		
        b_bufsize = ldb*n;
        x_bufsize = ldx*n;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }
	
#if LAPACKE_TEST_VERBOSE
   printf(" \n gerfs Double:  n: %d, fact: %c trans: %c  lda: %d  \
ldb: %d nrhs: %d  ldaf: %d  ldx: %d \n",  n, fact, trans, lda, 
                                          ldb, nrhs);
#endif
	
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, n*n);
    lapacke_gtest_alloc_float_buffer_pair( &b, &bref, n*nrhs);
    lapacke_gtest_alloc_float_buffer_pair( &x, &xref, n*nrhs);
    lapacke_gtest_alloc_float_buffer_pair( &af, &afref, (n*n));
    lapacke_gtest_alloc_float_buffer_pair( &ferr, &ferrref, nrhs);
    lapacke_gtest_alloc_float_buffer_pair( &berr, &berrref, nrhs);
    
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n+2);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (x==NULL) || (xref==NULL) || \
        (af==NULL) || (afref==NULL) || \
        (ferr==NULL) || (ferrref==NULL) || \
        (berr==NULL) || (berrref==NULL) || \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "pbtrs_float_parameters object: malloc error.";
       gerfs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, n*n);
	
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, n*nrhs);
	lapacke_gtest_init_float_buffer_pair_with_constant(x, xref, n*nrhs, 0.0);
	lapacke_gtest_init_float_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
	lapacke_gtest_init_float_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);
    memcpy(af, a, (n*n*sizeof(float)));
    memcpy(afref, a, (n*n*sizeof(float)));
    
   } /* end of Constructor  */

gerfs_float_parameters:: ~gerfs_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gerfs_float_parameters object: destructor invoked. \n");
#endif
   gerfs_free();
//  if (a != NULL)    free (a  );
//  if (aref != NULL) free (aref);
//  if (b != NULL)    free (b   );
//  if (bref != NULL) free (bref);
//  if (x != NULL)    free (x  );
//  if (xref != NULL) free (xref);
//  if (af != NULL)    free (af  );
//  if (afref != NULL) free (afref);
//  if (ferr != NULL)    free (ferr  );
//  if (ferrref != NULL) free (ferrref);
//  if (berr != NULL)    free (berr  );
//  if (berrref != NULL) free (berrref);
//  if (ipiv != NULL)    free (ipiv  );
//  if (ipivref != NULL) free (ipivref);
//  if( hModule != NULL) dlclose(hModule);
//  if(dModule != NULL) dlclose(dModule);

}

//  Test fixture class definition
class sgerfs_test  : public  ::testing::Test {
public:
   gerfs_float_parameters  *sgerfs_obj;
   void SetUp();  
   void TearDown () { delete sgerfs_obj; }
};


void sgerfs_test::SetUp(){

    /* LAPACKE SGERFS prototype */
    typedef int (*Fptr_NL_LAPACKE_sgerfs) (int matrix_layout,
											char trans,
											lapack_int n,
											lapack_int nrhs,
											const float* a,
											lapack_int lda,
											const float* af,
											lapack_int ldaf,
											const lapack_int* ipiv,
											const float* b,
											lapack_int ldb,
											float* x,
											lapack_int ldx,
											float* ferr,
											float* berr);

    Fptr_NL_LAPACKE_sgerfs SGERFS;


     /* LAPACKE SGETRF prototype */
    typedef int (*Fptr_NL_LAPACKE_sgetrf) ( int matrix_layout,lapack_int m,lapack_int n,
                                    float* a,lapack_int lda,lapack_int* ipiv );

    Fptr_NL_LAPACKE_sgetrf SGETRF;


    sgerfs_obj = new gerfs_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    sgerfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgerfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgerfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgerfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    SGERFS = (Fptr_NL_LAPACKE_sgerfs)dlsym(sgerfs_obj->hModule, "LAPACKE_sgerfs");
    ASSERT_TRUE(SGERFS != NULL) << "failed to get the Netlib LAPACKE_sgerfs symbol";
    
        SGETRF = (Fptr_NL_LAPACKE_sgetrf)dlsym(sgerfs_obj->hModule,"LAPACKE_sgetrf");
        ASSERT_TRUE(SGETRF != NULL) << "failed to get the Netlib LAPACKE_sgetrf symbol";
			

        sgerfs_obj->inforef = SGETRF( sgerfs_obj->matrix_layout,
                                      sgerfs_obj->n, sgerfs_obj->n,
                                      sgerfs_obj->afref,
                                      sgerfs_obj->lda,
									  sgerfs_obj->ipivref);
							   
        sgerfs_obj->info = LAPACKE_sgetrf( sgerfs_obj->matrix_layout,
                                           sgerfs_obj->n, sgerfs_obj->n,
                                           sgerfs_obj->af,
                                           sgerfs_obj->lda,
										   sgerfs_obj->ipiv);

    /* Compute the reference o/p by invoking Netlib-Lapack's API */	
    sgerfs_obj->inforef = SGERFS( sgerfs_obj->matrix_layout,
                                  sgerfs_obj->trans,
								  sgerfs_obj->n,
                                  sgerfs_obj->nrhs,
                                  sgerfs_obj->aref,
								  sgerfs_obj->lda, 
                                  sgerfs_obj->afref,
								  sgerfs_obj->ldaf,
                                  sgerfs_obj->ipivref,								  
								  sgerfs_obj->bref,
								  sgerfs_obj->ldb,
								  sgerfs_obj->xref,
								  sgerfs_obj->ldx,
								  sgerfs_obj->ferrref,
								  sgerfs_obj->berrref );

    /* Compute libflame's Lapacke o/p  */
    sgerfs_obj->info = LAPACKE_sgerfs( sgerfs_obj->matrix_layout,
                                  sgerfs_obj->trans,
								  sgerfs_obj->n,
                                  sgerfs_obj->nrhs,
                                  sgerfs_obj->a,
								  sgerfs_obj->lda, 
                                  sgerfs_obj->af,
								  sgerfs_obj->ldaf,
                                  sgerfs_obj->ipiv,
								  sgerfs_obj->b,
								  sgerfs_obj->ldb,
								  sgerfs_obj->x,
								  sgerfs_obj->ldx,
								  sgerfs_obj->ferr,
								  sgerfs_obj->berr);

    if( sgerfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_sgerfs is wrong\n", sgerfs_obj->info );
    }
    if( sgerfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sgerfs is wrong\n", 
        sgerfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    sgerfs_obj->diff_xerr =  computeDiff_s( sgerfs_obj->x_bufsize, 
                sgerfs_obj->x, sgerfs_obj->xref );

    sgerfs_obj->diff_berr =  computeDiff_s( sgerfs_obj->nrhs, 
                sgerfs_obj->berr, sgerfs_obj->berrref );
				
    sgerfs_obj->diff_ferr =  computeDiff_s( sgerfs_obj->nrhs, 
                sgerfs_obj->ferr, sgerfs_obj->ferrref );
}

TEST_F(sgerfs_test, sgerfs1) {
    EXPECT_NEAR(0.0, sgerfs_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgerfs_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgerfs_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgerfs_test, sgerfs2) {
    EXPECT_NEAR(0.0, sgerfs_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgerfs_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgerfs_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgerfs_test, sgerfs3) {
    EXPECT_NEAR(0.0, sgerfs_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgerfs_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgerfs_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgerfs_test, sgerfs4) {
    EXPECT_NEAR(0.0, sgerfs_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgerfs_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgerfs_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}


/* Begin gerfs_double_parameters  class definition */
class gerfs_double_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      double diff; // capture difference between ref o/p & libflame lapacke o/p.
	  double diff_berr, diff_ferr, diff_xerr;
      void *hModule, *dModule;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char fact; //  Must be 'F', 'N', or 'E'.
      char trans; //  Must be 'N' , 'T' or 'C'.

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
      lapack_int *ipiv, *ipivref; // pivot buffer
      
      /* Output parameters */
      double rpivot, rpivotref; // reciprocal pivot growth factor.
      double rcond, rcondref; // estimate of the reciprocal condition number of the matrix A.
      double* x, *xref;
      double* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      double* berr, *berrref; // component-wise relative backward error for each solution vector xj.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      gerfs_double_parameters ( int matrix_layout_i, char trans_i,
                                lapack_int n_i, lapack_int nrhs_i);
              
      ~gerfs_double_parameters (); 
};  /* end of gerfs_double_parameters  class definition */


/* Constructor gerfs_double_parameters definition */
gerfs_double_parameters:: gerfs_double_parameters ( int matrix_layout_i, 
                 char trans_i,  lapack_int n_i, lapack_int nrhs_i) {
                                    
    matrix_layout = matrix_layout_i;
    trans = trans_i;
    n = n_i;
    nrhs = nrhs_i;
    lda = n;
	ldaf = n;

    if(matrix_layout==LAPACK_COL_MAJOR){
		ldb = n;
		ldx = n;
		
        b_bufsize = ldb*nrhs;
        x_bufsize = ldx*nrhs;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
		ldb = nrhs;
		ldx = nrhs;
		
        b_bufsize = ldb*n;
        x_bufsize = ldx*n;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }
	
#if LAPACKE_TEST_VERBOSE
   printf(" \n gerfs Double:  n: %d, fact: %c trans: %c  lda: %d  \
ldb: %d nrhs: %d  ldaf: %d  ldx: %d \n",  n, fact, trans, lda, 
                                          ldb, nrhs);
#endif
	
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, n*n);
    lapacke_gtest_alloc_double_buffer_pair( &b, &bref, n*nrhs);
    lapacke_gtest_alloc_double_buffer_pair( &x, &xref, n*nrhs);
    lapacke_gtest_alloc_double_buffer_pair( &af, &afref, (n*n));
    lapacke_gtest_alloc_double_buffer_pair( &ferr, &ferrref, nrhs);
    lapacke_gtest_alloc_double_buffer_pair( &berr, &berrref, nrhs);
    
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n+2);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (x==NULL) || (xref==NULL) || \
        (af==NULL) || (afref==NULL) || \
        (ferr==NULL) || (ferrref==NULL) || \
        (berr==NULL) || (berrref==NULL) || \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "pbtrs_double_parameters object: malloc error.";
       gerfs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, n*n);
	
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, n*nrhs);
	lapacke_gtest_init_double_buffer_pair_with_constant(x, xref, n*nrhs, 0.0);
	lapacke_gtest_init_double_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
	lapacke_gtest_init_double_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);
    memcpy(af, a, (n*n*sizeof(double)));
    memcpy(afref, a, (n*n*sizeof(double)));
    
   } /* end of Constructor  */

gerfs_double_parameters:: ~gerfs_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gerfs_double_parameters object: destructor invoked. \n");
#endif
   gerfs_free();
}

//  Test fixture class definition
class dgerfs_test  : public  ::testing::Test {
public:
   gerfs_double_parameters  *dgerfs_obj;
   void SetUp();  
   void TearDown () { delete dgerfs_obj; }
};


void dgerfs_test::SetUp(){

    /* LAPACKE DGERFS prototype */
    typedef int (*Fptr_NL_LAPACKE_dgerfs) (int matrix_layout,
											char trans,
											lapack_int n,
											lapack_int nrhs,
											const double* a,
											lapack_int lda,
											const double* af,
											lapack_int ldaf,
											const lapack_int* ipiv,
											const double* b,
											lapack_int ldb,
											double* x,
											lapack_int ldx,
											double* ferr,
											double* berr);

    Fptr_NL_LAPACKE_dgerfs DGERFS;


     /* LAPACKE DGETRF prototype */
    typedef int (*Fptr_NL_LAPACKE_dgetrf) ( int matrix_layout,lapack_int m,lapack_int n,
                                    double* a,lapack_int lda,lapack_int* ipiv );

    Fptr_NL_LAPACKE_dgetrf DGETRF;


    dgerfs_obj = new gerfs_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    dgerfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgerfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgerfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgerfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    DGERFS = (Fptr_NL_LAPACKE_dgerfs)dlsym(dgerfs_obj->hModule, "LAPACKE_dgerfs");
    ASSERT_TRUE(DGERFS != NULL) << "failed to get the Netlib LAPACKE_dgerfs symbol";
    
        DGETRF = (Fptr_NL_LAPACKE_dgetrf)dlsym(dgerfs_obj->hModule,"LAPACKE_dgetrf");
        ASSERT_TRUE(DGETRF != NULL) << "failed to get the Netlib LAPACKE_dgetrf symbol";
			

        dgerfs_obj->inforef = DGETRF( dgerfs_obj->matrix_layout,
                                      dgerfs_obj->n, dgerfs_obj->n,
                                      dgerfs_obj->afref,
                                      dgerfs_obj->lda,
									  dgerfs_obj->ipivref);
							   
        dgerfs_obj->info = LAPACKE_dgetrf( dgerfs_obj->matrix_layout,
                                           dgerfs_obj->n, dgerfs_obj->n,
                                           dgerfs_obj->af,
                                           dgerfs_obj->lda,
										   dgerfs_obj->ipiv);

    /* Compute the reference o/p by invoking Netlib-Lapack's API */	
    dgerfs_obj->inforef = DGERFS( dgerfs_obj->matrix_layout,
                                  dgerfs_obj->trans,
								  dgerfs_obj->n,
                                  dgerfs_obj->nrhs,
                                  dgerfs_obj->aref,
								  dgerfs_obj->lda, 
                                  dgerfs_obj->afref,
								  dgerfs_obj->ldaf,
                                  dgerfs_obj->ipivref,								  
								  dgerfs_obj->bref,
								  dgerfs_obj->ldb,
								  dgerfs_obj->xref,
								  dgerfs_obj->ldx,
								  dgerfs_obj->ferrref,
								  dgerfs_obj->berrref );

    /* Compute libflame's Lapacke o/p  */
    dgerfs_obj->info = LAPACKE_dgerfs( dgerfs_obj->matrix_layout,
                                  dgerfs_obj->trans,
								  dgerfs_obj->n,
                                  dgerfs_obj->nrhs,
                                  dgerfs_obj->a,
								  dgerfs_obj->lda, 
                                  dgerfs_obj->af,
								  dgerfs_obj->ldaf,
                                  dgerfs_obj->ipiv,
								  dgerfs_obj->b,
								  dgerfs_obj->ldb,
								  dgerfs_obj->x,
								  dgerfs_obj->ldx,
								  dgerfs_obj->ferr,
								  dgerfs_obj->berr);

    if( dgerfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dgerfs is wrong\n", dgerfs_obj->info );
    }
    if( dgerfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dgerfs is wrong\n", 
        dgerfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dgerfs_obj->diff_xerr =  computeDiff_d( dgerfs_obj->x_bufsize, 
                dgerfs_obj->x, dgerfs_obj->xref );

    dgerfs_obj->diff_berr =  computeDiff_d( dgerfs_obj->nrhs, 
                dgerfs_obj->berr, dgerfs_obj->berrref );
				
    dgerfs_obj->diff_ferr =  computeDiff_d( dgerfs_obj->nrhs, 
                dgerfs_obj->ferr, dgerfs_obj->ferrref );
}

TEST_F(dgerfs_test, dgerfs1) {
    EXPECT_NEAR(0.0, dgerfs_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgerfs_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgerfs_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgerfs_test, dgerfs2) {
    EXPECT_NEAR(0.0, dgerfs_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgerfs_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgerfs_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgerfs_test, dgerfs3) {
    EXPECT_NEAR(0.0, dgerfs_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgerfs_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgerfs_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgerfs_test, dgerfs4) {
    EXPECT_NEAR(0.0, dgerfs_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgerfs_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgerfs_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

/* Begin gerfs_scomplex_parameters  class definition */
class gerfs_scomplex_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      float diff; // capture difference between ref o/p & libflame lapacke o/p.
	  float diff_berr, diff_ferr, diff_xerr;
      void *hModule, *dModule;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char fact; //  Must be 'F', 'N', or 'E'.
      char trans; //  Must be 'N' , 'T' or 'C'.

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
      lapack_int *ipiv, *ipivref; // pivot buffer
      
      /* Output parameters */
      lapack_complex_float* x, *xref;
      float* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      float* berr, *berrref; // component-wise relative backward error for each solution vector xj.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      gerfs_scomplex_parameters ( int matrix_layout_i, char trans_i,
                                lapack_int n_i, lapack_int nrhs_i);
              
      ~gerfs_scomplex_parameters (); 
};  /* end of gerfs_scomplex_parameters  class definition */


/* Constructor gerfs_scomplex_parameters definition */
gerfs_scomplex_parameters:: gerfs_scomplex_parameters ( int matrix_layout_i, 
                 char trans_i,  lapack_int n_i, lapack_int nrhs_i) {
                                    
    matrix_layout = matrix_layout_i;
    trans = trans_i;
    n = n_i;
    nrhs = nrhs_i;
    lda = n;
	ldaf = n;

    if(matrix_layout==LAPACK_COL_MAJOR){
		ldb = n;
		ldx = n;
		
        b_bufsize = ldb*nrhs;
        x_bufsize = ldx*nrhs;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
		ldb = nrhs;
		ldx = nrhs;
		
        b_bufsize = ldb*n;
        x_bufsize = ldx*n;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }
	
#if LAPACKE_TEST_VERBOSE
   printf(" \n gerfs Double:  n: %d, fact: %c trans: %c  lda: %d  \
ldb: %d nrhs: %d  ldaf: %d  ldx: %d \n",  n, fact, trans, lda, 
                                          ldb, nrhs);
#endif
	
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, n*n);
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &b, &bref, n*nrhs);
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &x, &xref, n*nrhs);
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &af, &afref, (n*n));
    lapacke_gtest_alloc_float_buffer_pair( &ferr, &ferrref, nrhs);
    lapacke_gtest_alloc_float_buffer_pair( &berr, &berrref, nrhs);
    
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n+2);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (x==NULL) || (xref==NULL) || \
        (af==NULL) || (afref==NULL) || \
        (ferr==NULL) || (ferrref==NULL) || \
        (berr==NULL) || (berrref==NULL) || \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "pbtrs_scomplex_parameters object: malloc error.";
       gerfs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, n*n);
	
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, n*nrhs);
	lapacke_gtest_init_scomplex_buffer_pair_with_constant(x, xref, n*nrhs, 0.0);
	lapacke_gtest_init_float_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
	lapacke_gtest_init_float_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);
    memcpy(af, a, (n*n*sizeof(lapack_complex_float)));
    memcpy(afref, a, (n*n*sizeof(lapack_complex_float)));
    
   } /* end of Constructor  */

gerfs_scomplex_parameters:: ~gerfs_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gerfs_scomplex_parameters object: destructor invoked. \n");
#endif
   gerfs_free();
}

//  Test fixture class definition
class cgerfs_test  : public  ::testing::Test {
public:
   gerfs_scomplex_parameters  *cgerfs_obj;
   void SetUp();  
   void TearDown () { delete cgerfs_obj; }
};


void cgerfs_test::SetUp(){

    /* LAPACKE CGERFS prototype */
    typedef int (*Fptr_NL_LAPACKE_cgerfs) (int matrix_layout,
											char trans,
											lapack_int n,
											lapack_int nrhs,
											const lapack_complex_float* a,
											lapack_int lda,
											const lapack_complex_float* af,
											lapack_int ldaf,
											const lapack_int* ipiv,
											const lapack_complex_float* b,
											lapack_int ldb,
											lapack_complex_float* x,
											lapack_int ldx,
											float* ferr,
											float* berr);

    Fptr_NL_LAPACKE_cgerfs CGERFS;


     /* LAPACKE CGETRF prototype */
    typedef int (*Fptr_NL_LAPACKE_cgetrf) ( int matrix_layout,lapack_int m,lapack_int n,
                                    lapack_complex_float* a,lapack_int lda,lapack_int* ipiv );

    Fptr_NL_LAPACKE_cgetrf CGETRF;


    cgerfs_obj = new gerfs_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    cgerfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgerfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgerfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgerfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    CGERFS = (Fptr_NL_LAPACKE_cgerfs)dlsym(cgerfs_obj->hModule, "LAPACKE_cgerfs");
    ASSERT_TRUE(CGERFS != NULL) << "failed to get the Netlib LAPACKE_cgerfs symbol";
    
        CGETRF = (Fptr_NL_LAPACKE_cgetrf)dlsym(cgerfs_obj->hModule,"LAPACKE_cgetrf");
        ASSERT_TRUE(CGETRF != NULL) << "failed to get the Netlib LAPACKE_cgetrf symbol";
			

        cgerfs_obj->inforef = CGETRF( cgerfs_obj->matrix_layout,
                                      cgerfs_obj->n, cgerfs_obj->n,
                                      cgerfs_obj->afref,
                                      cgerfs_obj->lda,
									  cgerfs_obj->ipivref);
							   
        cgerfs_obj->info = LAPACKE_cgetrf( cgerfs_obj->matrix_layout,
                                           cgerfs_obj->n, cgerfs_obj->n,
                                           cgerfs_obj->af,
                                           cgerfs_obj->lda,
										   cgerfs_obj->ipiv);

    /* Compute the reference o/p by invoking Netlib-Lapack's API */	
    cgerfs_obj->inforef = CGERFS( cgerfs_obj->matrix_layout,
                                  cgerfs_obj->trans,
								  cgerfs_obj->n,
                                  cgerfs_obj->nrhs,
                                  cgerfs_obj->aref,
								  cgerfs_obj->lda, 
                                  cgerfs_obj->afref,
								  cgerfs_obj->ldaf,
                                  cgerfs_obj->ipivref,								  
								  cgerfs_obj->bref,
								  cgerfs_obj->ldb,
								  cgerfs_obj->xref,
								  cgerfs_obj->ldx,
								  cgerfs_obj->ferrref,
								  cgerfs_obj->berrref );

    /* Compute libflame's Lapacke o/p  */
    cgerfs_obj->info = LAPACKE_cgerfs( cgerfs_obj->matrix_layout,
                                  cgerfs_obj->trans,
								  cgerfs_obj->n,
                                  cgerfs_obj->nrhs,
                                  cgerfs_obj->a,
								  cgerfs_obj->lda, 
                                  cgerfs_obj->af,
								  cgerfs_obj->ldaf,
                                  cgerfs_obj->ipiv,
								  cgerfs_obj->b,
								  cgerfs_obj->ldb,
								  cgerfs_obj->x,
								  cgerfs_obj->ldx,
								  cgerfs_obj->ferr,
								  cgerfs_obj->berr);

    if( cgerfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cgerfs is wrong\n", cgerfs_obj->info );
    }
    if( cgerfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cgerfs is wrong\n", 
        cgerfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cgerfs_obj->diff_xerr =  computeDiff_c( cgerfs_obj->x_bufsize, 
                cgerfs_obj->x, cgerfs_obj->xref );

    cgerfs_obj->diff_berr =  computeDiff_s( cgerfs_obj->nrhs, 
                cgerfs_obj->berr, cgerfs_obj->berrref );
				
    cgerfs_obj->diff_ferr =  computeDiff_s( cgerfs_obj->nrhs, 
                cgerfs_obj->ferr, cgerfs_obj->ferrref );
}

TEST_F(cgerfs_test, cgerfs1) {
    EXPECT_NEAR(0.0, cgerfs_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgerfs_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgerfs_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgerfs_test, cgerfs2) {
    EXPECT_NEAR(0.0, cgerfs_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgerfs_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgerfs_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgerfs_test, cgerfs3) {
    EXPECT_NEAR(0.0, cgerfs_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgerfs_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgerfs_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgerfs_test, cgerfs4) {
    EXPECT_NEAR(0.0, cgerfs_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgerfs_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgerfs_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

/* Begin gerfs_dcomplex_parameters  class definition */
class gerfs_dcomplex_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      double diff; // capture difference between ref o/p & libflame lapacke o/p.
	  double diff_berr, diff_ferr, diff_xerr;
      void *hModule, *dModule;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char fact; //  Must be 'F', 'N', or 'E'.
      char trans; //  Must be 'N' , 'T' or 'C'.

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
      lapack_int *ipiv, *ipivref; // pivot buffer
      
      /* Output parameters */
      lapack_complex_double* x, *xref;
      double* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      double* berr, *berrref; // component-wise relative backward error for each solution vector xj.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      gerfs_dcomplex_parameters ( int matrix_layout_i, char trans_i,
                                lapack_int n_i, lapack_int nrhs_i);
              
      ~gerfs_dcomplex_parameters (); 
};  /* end of gerfs_dcomplex_parameters  class definition */


/* Constructor gerfs_dcomplex_parameters definition */
gerfs_dcomplex_parameters:: gerfs_dcomplex_parameters ( int matrix_layout_i, 
                 char trans_i,  lapack_int n_i, lapack_int nrhs_i) {
                                    
    matrix_layout = matrix_layout_i;
    trans = trans_i;
    n = n_i;
    nrhs = nrhs_i;
    lda = n;
	ldaf = n;

    if(matrix_layout==LAPACK_COL_MAJOR){
		ldb = n;
		ldx = n;
		
        b_bufsize = ldb*nrhs;
        x_bufsize = ldx*nrhs;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
		ldb = nrhs;
		ldx = nrhs;
		
        b_bufsize = ldb*n;
        x_bufsize = ldx*n;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }
	
#if LAPACKE_TEST_VERBOSE
   printf(" \n gerfs Double:  n: %d, fact: %c trans: %c  lda: %d  \
ldb: %d nrhs: %d  ldaf: %d  ldx: %d \n",  n, fact, trans, lda, 
                                          ldb, nrhs);
#endif
	
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, n*n);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &b, &bref, n*nrhs);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &x, &xref, n*nrhs);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &af, &afref, (n*n));
    lapacke_gtest_alloc_double_buffer_pair( &ferr, &ferrref, nrhs);
    lapacke_gtest_alloc_double_buffer_pair( &berr, &berrref, nrhs);
    
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n+2);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (x==NULL) || (xref==NULL) || \
        (af==NULL) || (afref==NULL) || \
        (ferr==NULL) || (ferrref==NULL) || \
        (berr==NULL) || (berrref==NULL) || \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "pbtrs_dcomplex_parameters object: malloc error.";
       gerfs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, n*n);
	
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, n*nrhs);
	lapacke_gtest_init_dcomplex_buffer_pair_with_constant(x, xref, n*nrhs, 0.0);
	lapacke_gtest_init_double_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
	lapacke_gtest_init_double_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);
    memcpy(af, a, (n*n*sizeof(lapack_complex_double)));
    memcpy(afref, a, (n*n*sizeof(lapack_complex_double)));
    
   } /* end of Constructor  */

gerfs_dcomplex_parameters:: ~gerfs_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gerfs_dcomplex_parameters object: destructor invoked. \n");
#endif
   gerfs_free();
}

//  Test fixture class definition
class zgerfs_test  : public  ::testing::Test {
public:
   gerfs_dcomplex_parameters  *zgerfs_obj;
   void SetUp();  
   void TearDown () { delete zgerfs_obj; }
};


void zgerfs_test::SetUp(){

    /* LAPACKE ZGERFS prototype */
    typedef int (*Fptr_NL_LAPACKE_zgerfs) (int matrix_layout,
											char trans,
											lapack_int n,
											lapack_int nrhs,
											const lapack_complex_double* a,
											lapack_int lda,
											const lapack_complex_double* af,
											lapack_int ldaf,
											const lapack_int* ipiv,
											const lapack_complex_double* b,
											lapack_int ldb,
											lapack_complex_double* x,
											lapack_int ldx,
											double* ferr,
											double* berr);

    Fptr_NL_LAPACKE_zgerfs ZGERFS;


     /* LAPACKE ZGETRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zgetrf) ( int matrix_layout,lapack_int m,lapack_int n,
                                    lapack_complex_double* a,lapack_int lda,lapack_int* ipiv );

    Fptr_NL_LAPACKE_zgetrf ZGETRF;


    zgerfs_obj = new gerfs_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    zgerfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgerfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgerfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgerfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    ZGERFS = (Fptr_NL_LAPACKE_zgerfs)dlsym(zgerfs_obj->hModule, "LAPACKE_zgerfs");
    ASSERT_TRUE(ZGERFS != NULL) << "failed to get the Netlib LAPACKE_zgerfs symbol";
    
        ZGETRF = (Fptr_NL_LAPACKE_zgetrf)dlsym(zgerfs_obj->hModule,"LAPACKE_zgetrf");
        ASSERT_TRUE(ZGETRF != NULL) << "failed to get the Netlib LAPACKE_zgetrf symbol";
			

        zgerfs_obj->inforef = ZGETRF( zgerfs_obj->matrix_layout,
                                      zgerfs_obj->n, zgerfs_obj->n,
                                      zgerfs_obj->afref,
                                      zgerfs_obj->lda,
									  zgerfs_obj->ipivref);
							   
        zgerfs_obj->info = LAPACKE_zgetrf( zgerfs_obj->matrix_layout,
                                           zgerfs_obj->n, zgerfs_obj->n,
                                           zgerfs_obj->af,
                                           zgerfs_obj->lda,
										   zgerfs_obj->ipiv);

    /* Compute the reference o/p by invoking Netlib-Lapack's API */	
    zgerfs_obj->inforef = ZGERFS( zgerfs_obj->matrix_layout,
                                  zgerfs_obj->trans,
								  zgerfs_obj->n,
                                  zgerfs_obj->nrhs,
                                  zgerfs_obj->aref,
								  zgerfs_obj->lda, 
                                  zgerfs_obj->afref,
								  zgerfs_obj->ldaf,
                                  zgerfs_obj->ipivref,								  
								  zgerfs_obj->bref,
								  zgerfs_obj->ldb,
								  zgerfs_obj->xref,
								  zgerfs_obj->ldx,
								  zgerfs_obj->ferrref,
								  zgerfs_obj->berrref );

    /* Compute libflame's Lapacke o/p  */
    zgerfs_obj->info = LAPACKE_zgerfs( zgerfs_obj->matrix_layout,
                                  zgerfs_obj->trans,
								  zgerfs_obj->n,
                                  zgerfs_obj->nrhs,
                                  zgerfs_obj->a,
								  zgerfs_obj->lda, 
                                  zgerfs_obj->af,
								  zgerfs_obj->ldaf,
                                  zgerfs_obj->ipiv,
								  zgerfs_obj->b,
								  zgerfs_obj->ldb,
								  zgerfs_obj->x,
								  zgerfs_obj->ldx,
								  zgerfs_obj->ferr,
								  zgerfs_obj->berr);

    if( zgerfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zgerfs is wrong\n", zgerfs_obj->info );
    }
    if( zgerfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zgerfs is wrong\n", 
        zgerfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zgerfs_obj->diff_xerr =  computeDiff_z( zgerfs_obj->x_bufsize, 
                zgerfs_obj->x, zgerfs_obj->xref );

    zgerfs_obj->diff_berr =  computeDiff_d( zgerfs_obj->nrhs, 
                zgerfs_obj->berr, zgerfs_obj->berrref );
				
    zgerfs_obj->diff_ferr =  computeDiff_d( zgerfs_obj->nrhs, 
                zgerfs_obj->ferr, zgerfs_obj->ferrref );
}

TEST_F(zgerfs_test, zgerfs1) {
    EXPECT_NEAR(0.0, zgerfs_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgerfs_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgerfs_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgerfs_test, zgerfs2) {
    EXPECT_NEAR(0.0, zgerfs_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgerfs_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgerfs_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgerfs_test, zgerfs3) {
    EXPECT_NEAR(0.0, zgerfs_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgerfs_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgerfs_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgerfs_test, zgerfs4) {
    EXPECT_NEAR(0.0, zgerfs_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgerfs_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgerfs_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}


