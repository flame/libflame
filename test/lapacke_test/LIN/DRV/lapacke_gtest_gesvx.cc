#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"

#define gesvx_free() \
  if (a != NULL)    free (a  ); \
  if (aref != NULL) free (aref); \
  if (b != NULL)    free (b   ); \
  if (bref != NULL) free (bref); \
  if (x != NULL)    free (x  ); \
  if (xref != NULL) free (xref); \
  if (af != NULL)    free (af  ); \
  if (afref != NULL) free (afref); \
  if (r != NULL)    free (r  ); \
  if (rref != NULL) free (rref); \
  if (c != NULL)    free (c  ); \
  if (cref != NULL) free (cref); \
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

/* Begin gesvx_double_parameters  class definition */
class gesvx_double_parameters{
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
      double* r, *rref; // the row scale factors for A.
      double* c, *cref; // column scale factors for A.
      lapack_int *ipiv, *ipivref; // pivot buffer
      
      /* Output parameters */
      char   equed, equedref; //  Must be 'N', 'R', 'C', or 'B'.
      double rpivot, rpivotref; // reciprocal pivot growth factor.
      double rcond, rcondref; // estimate of the reciprocal condition number of the matrix A.
      double* x, *xref;
      double* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      double* berr, *berrref; // component-wise relative backward error for each solution vector xj.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      gesvx_double_parameters ( int matrix_layout_i, char fact_i, char trans_i,
                               char equed_i, lapack_int n_i, lapack_int nrhs_i);
              
      ~gesvx_double_parameters (); 
};  /* end of gesvx_double_parameters  class definition */


/* Constructor gesvx_double_parameters definition */
gesvx_double_parameters:: gesvx_double_parameters ( int matrix_layout_i, 
                char fact_i, char trans_i, char equed_i, lapack_int n_i,
                  lapack_int nrhs_i) {
                                    
    matrix_layout = matrix_layout_i;
    fact= fact_i;
    trans = trans_i;
    // equed acts as input when " fact = 'F' " else output
    equed = equed_i;
    equedref = equed_i;
    n = n_i;
    nrhs = nrhs_i;
    lda = n;
    ldaf = n;
    
#if LAPACKE_TEST_VERBOSE
   printf(" \n gesvx Double:  n: %d, fact: %c trans: %c  lda: %d  \
ldb: %d nrhs: %d  ldaf: %d  ldx: %d \n",  n, fact, trans, lda, 
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
    lapacke_gtest_alloc_double_buffer_pair( &r, &rref, n);
    lapacke_gtest_alloc_double_buffer_pair( &c, &cref, n);
    lapacke_gtest_alloc_double_buffer_pair( &ferr, &ferrref, nrhs);
    lapacke_gtest_alloc_double_buffer_pair( &berr, &berrref, nrhs);
    
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (x==NULL) || (xref==NULL) || \
        (af==NULL) || (afref==NULL) || \
        (r==NULL) || (rref==NULL) || \
        (c==NULL) || (cref==NULL) || \
        (ferr==NULL) || (ferrref==NULL) || \
        (berr==NULL) || (berrref==NULL) || \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "pbtrs_double_parameters object: malloc error.";
       gesvx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, n*n);
    memcpy(af, a, (n*n*sizeof(double)));
    memcpy(afref, a, (n*n*sizeof(double)));
    
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, b_bufsize);
    
    lapacke_gtest_init_double_buffer_pair_with_constant(x, xref, x_bufsize, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(r, rref, n, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(c, cref, n, 0.0);

    
   } /* end of Constructor  */

gesvx_double_parameters:: ~gesvx_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gesvx_double_parameters object: destructor invoked. \n");
#endif
   gesvx_free();
}

//  Test fixture class definition
class dgesvx_test  : public  ::testing::Test {
public:
   gesvx_double_parameters  *dgesvx_obj;
   void SetUp();  
   void TearDown () { delete dgesvx_obj; }
};


void dgesvx_test::SetUp(){

    /* LAPACKE DGESVX prototype */
    typedef int (*Fptr_NL_LAPACKE_dgesvx) (int matrix_layout, char fact,
                   char trans, lapack_int n, lapack_int nrhs, double *a,
                            lapack_int lda, double *af, lapack_int ldaf,
                            lapack_int *ipiv, char *equed, double *r,
                            double *c, double *b, lapack_int ldb,
                            double *x, lapack_int ldx, double *rcond,
                            double *ferr, double *berr, double *rpivot);

    Fptr_NL_LAPACKE_dgesvx DGESVX;


     /* LAPACKE DGETRF prototype */
    typedef int (*Fptr_NL_LAPACKE_dgetrf) ( int matrix_layout,lapack_int m,lapack_int n,
                                    double* a,lapack_int lda,lapack_int* ipiv );

    Fptr_NL_LAPACKE_dgetrf DGETRF;


    dgesvx_obj = new gesvx_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].fact,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].equed,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    dgesvx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgesvx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgesvx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgesvx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    DGESVX = (Fptr_NL_LAPACKE_dgesvx)dlsym(dgesvx_obj->hModule, "LAPACKE_dgesvx");
    ASSERT_TRUE(DGESVX != NULL) << "failed to get the Netlib LAPACKE_dgesvx symbol";
    
    /* if  'fact' == 'F', then the i/p matrix A should be in 
       factored form.
       so, invoke the getrf API to compute the factorized A.  */
    if(dgesvx_obj->fact == 'F') {
        DGETRF = (Fptr_NL_LAPACKE_dgetrf)dlsym(dgesvx_obj->hModule,"LAPACKE_dgetrf");
        ASSERT_TRUE(DGETRF != NULL) << "failed to get the Netlib LAPACKE_dgetrf symbol";
            

        dgesvx_obj->inforef = DGETRF( dgesvx_obj->matrix_layout,
                                      dgesvx_obj->n, dgesvx_obj->n,
                                      dgesvx_obj->afref,
                                      dgesvx_obj->lda,
                                      dgesvx_obj->ipivref);
                               
        dgesvx_obj->info = LAPACKE_dgetrf( dgesvx_obj->matrix_layout,
                                           dgesvx_obj->n, dgesvx_obj->n,
                                           dgesvx_obj->af,
                                           dgesvx_obj->lda,
                                           dgesvx_obj->ipiv);
    }

    typedef int (*Fptr_NL_LAPACKE_dgesvx) (int matrix_layout, char fact,
                   char trans, lapack_int n, lapack_int nrhs, double *a,
                            lapack_int lda, double *af, lapack_int ldaf,
                            lapack_int *ipiv, char *equed, double *r,
                            double *c, double *b, lapack_int ldb,
                            double *x, lapack_int ldx, double *rcond,
                            double *ferr, double *berr, double *rpivot);

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    dgesvx_obj->inforef = DGESVX( dgesvx_obj->matrix_layout, dgesvx_obj->fact,
                                  dgesvx_obj->trans, dgesvx_obj->n,
                                  dgesvx_obj->nrhs,
                                  dgesvx_obj->aref, dgesvx_obj->lda, 
                                  dgesvx_obj->afref, dgesvx_obj->ldaf,
                                  dgesvx_obj->ipivref, &dgesvx_obj->equedref,
                                  dgesvx_obj->rref, dgesvx_obj->cref,                                 
                                  dgesvx_obj->bref, dgesvx_obj->ldb,
                                  dgesvx_obj->xref, dgesvx_obj->ldx,
                                  &dgesvx_obj->rcondref, 
                                  dgesvx_obj->ferrref,
                                  dgesvx_obj->berrref,
                                  &dgesvx_obj->rpivotref);

    /* Compute libflame's Lapacke o/p  */
    dgesvx_obj->info = LAPACKE_dgesvx( dgesvx_obj->matrix_layout, dgesvx_obj->fact,
                                  dgesvx_obj->trans, dgesvx_obj->n,
                                  dgesvx_obj->nrhs,
                                  dgesvx_obj->a, dgesvx_obj->lda, 
                                  dgesvx_obj->af, dgesvx_obj->ldaf,
                                  dgesvx_obj->ipiv, &dgesvx_obj->equed,
                                  dgesvx_obj->r, dgesvx_obj->c,                               
                                  dgesvx_obj->b, dgesvx_obj->ldb,
                                  dgesvx_obj->x, dgesvx_obj->ldx,
                                  &dgesvx_obj->rcond, 
                                  dgesvx_obj->ferr,
                                  dgesvx_obj->berr,
                                  &dgesvx_obj->rpivot);

    if( dgesvx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dgesvx is wrong\n", dgesvx_obj->info );
    }
    if( dgesvx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dgesvx is wrong\n", 
        dgesvx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dgesvx_obj->diff =  computeDiff_d( dgesvx_obj->b_bufsize, 
                dgesvx_obj->b, dgesvx_obj->bref );

    dgesvx_obj->diff_xerr =  computeDiff_d( dgesvx_obj->x_bufsize, 
                dgesvx_obj->x, dgesvx_obj->xref );

    dgesvx_obj->diff_berr =  computeDiff_d( dgesvx_obj->nrhs, 
                dgesvx_obj->berr, dgesvx_obj->berrref );
                
    dgesvx_obj->diff_ferr =  computeDiff_d( dgesvx_obj->nrhs, 
                dgesvx_obj->ferr, dgesvx_obj->ferrref );
}

TEST_F(dgesvx_test, dgesvx1) {
    EXPECT_NEAR(0.0, dgesvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgesvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgesvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgesvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (dgesvx_obj->rpivot - dgesvx_obj->rpivotref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgesvx_test, dgesvx2) {
    EXPECT_NEAR(0.0, dgesvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgesvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgesvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgesvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (dgesvx_obj->rpivot - dgesvx_obj->rpivotref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgesvx_test, dgesvx3) {
    EXPECT_NEAR(0.0, dgesvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgesvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgesvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgesvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (dgesvx_obj->rpivot - dgesvx_obj->rpivotref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgesvx_test, dgesvx4) {
    EXPECT_NEAR(0.0, dgesvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgesvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgesvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgesvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (dgesvx_obj->rpivot - dgesvx_obj->rpivotref),
                LAPACKE_GTEST_THRESHOLD);
}

/* Begin gesvx_float_parameters  class definition */
class gesvx_float_parameters{
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
      float* r, *rref; // the row scale factors for A.
      float* c, *cref; // column scale factors for A.
      lapack_int *ipiv, *ipivref; // pivot buffer
      
      /* Output parameters */
      char   equed, equedref; //  Must be 'N', 'R', 'C', or 'B'.
      float rpivot, rpivotref; // reciprocal pivot growth factor.
      float rcond, rcondref; // estimate of the reciprocal condition number of the matrix A.
      float* x, *xref;
      float* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      float* berr, *berrref; // component-wise relative backward error for each solution vector xj.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      gesvx_float_parameters ( int matrix_layout_i, char fact_i, char trans_i,
                               char equed_i, lapack_int n_i, lapack_int nrhs_i);
              
      ~gesvx_float_parameters (); 
};  /* end of gesvx_float_parameters  class definition */


/* Constructor gesvx_float_parameters definition */
gesvx_float_parameters:: gesvx_float_parameters ( int matrix_layout_i, 
                char fact_i, char trans_i, char equed_i, lapack_int n_i,
                  lapack_int nrhs_i) {
                                    
    matrix_layout = matrix_layout_i;
    fact= fact_i;
    trans = trans_i;
    // equed acts as input when " fact = 'F' " else output
    equed = equed_i;
    equedref = equed_i;
    n = n_i;
    nrhs = nrhs_i;
    lda = n;
    ldaf = n;
    
#if LAPACKE_TEST_VERBOSE
   printf(" \n gesvx Double:  n: %d, fact: %c trans: %c  lda: %d  \
ldb: %d nrhs: %d  ldaf: %d  ldx: %d \n",  n, fact, trans, lda, 
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
    lapacke_gtest_alloc_float_buffer_pair( &r, &rref, n);
    lapacke_gtest_alloc_float_buffer_pair( &c, &cref, n);
    lapacke_gtest_alloc_float_buffer_pair( &ferr, &ferrref, nrhs);
    lapacke_gtest_alloc_float_buffer_pair( &berr, &berrref, nrhs);
    
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (x==NULL) || (xref==NULL) || \
        (af==NULL) || (afref==NULL) || \
        (r==NULL) || (rref==NULL) || \
        (c==NULL) || (cref==NULL) || \
        (ferr==NULL) || (ferrref==NULL) || \
        (berr==NULL) || (berrref==NULL) || \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "pbtrs_float_parameters object: malloc error.";
       gesvx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, n*n);
    memcpy(af, a, (n*n*sizeof(float)));
    memcpy(afref, a, (n*n*sizeof(float)));
    
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, b_bufsize);
    
    lapacke_gtest_init_float_buffer_pair_with_constant(x, xref, x_bufsize, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(r, rref, n, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(c, cref, n, 0.0);

    
   } /* end of Constructor  */

gesvx_float_parameters:: ~gesvx_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gesvx_float_parameters object: destructor invoked. \n");
#endif
   gesvx_free();
}

//  Test fixture class definition
class sgesvx_test  : public  ::testing::Test {
public:
   gesvx_float_parameters  *sgesvx_obj;
   void SetUp();  
   void TearDown () { delete sgesvx_obj; }
};


void sgesvx_test::SetUp(){

    /* LAPACKE SGESVX prototype */
    typedef int (*Fptr_NL_LAPACKE_sgesvx) (int matrix_layout, char fact,
                   char trans, lapack_int n, lapack_int nrhs, float *a,
                            lapack_int lda, float *af, lapack_int ldaf,
                            lapack_int *ipiv, char *equed, float *r,
                            float *c, float *b, lapack_int ldb,
                            float *x, lapack_int ldx, float *rcond,
                            float *ferr, float *berr, float *rpivot);

    Fptr_NL_LAPACKE_sgesvx SGESVX;


     /* LAPACKE SGETRF prototype */
    typedef int (*Fptr_NL_LAPACKE_sgetrf) ( int matrix_layout,lapack_int m,lapack_int n,
                                    float* a,lapack_int lda,lapack_int* ipiv );

    Fptr_NL_LAPACKE_sgetrf SGETRF;


    sgesvx_obj = new gesvx_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].fact,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].equed,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    sgesvx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgesvx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgesvx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgesvx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    SGESVX = (Fptr_NL_LAPACKE_sgesvx)dlsym(sgesvx_obj->hModule, "LAPACKE_sgesvx");
    ASSERT_TRUE(SGESVX != NULL) << "failed to get the Netlib LAPACKE_sgesvx symbol";
    
    /* if  'fact' == 'F', then the i/p matrix A should be in 
       factored form.
       so, invoke the getrf API to compute the factorized A.  */
    if(sgesvx_obj->fact == 'F') {
        SGETRF = (Fptr_NL_LAPACKE_sgetrf)dlsym(sgesvx_obj->hModule,"LAPACKE_sgetrf");
        ASSERT_TRUE(SGETRF != NULL) << "failed to get the Netlib LAPACKE_sgetrf symbol";
            

        sgesvx_obj->inforef = SGETRF( sgesvx_obj->matrix_layout,
                                      sgesvx_obj->n, sgesvx_obj->n,
                                      sgesvx_obj->afref,
                                      sgesvx_obj->lda,
                                      sgesvx_obj->ipivref);
                               
        sgesvx_obj->info = LAPACKE_sgetrf( sgesvx_obj->matrix_layout,
                                           sgesvx_obj->n, sgesvx_obj->n,
                                           sgesvx_obj->af,
                                           sgesvx_obj->lda,
                                           sgesvx_obj->ipiv);
    }

    typedef int (*Fptr_NL_LAPACKE_sgesvx) (int matrix_layout, char fact,
                   char trans, lapack_int n, lapack_int nrhs, float *a,
                            lapack_int lda, float *af, lapack_int ldaf,
                            lapack_int *ipiv, char *equed, float *r,
                            float *c, float *b, lapack_int ldb,
                            float *x, lapack_int ldx, float *rcond,
                            float *ferr, float *berr, float *rpivot);

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    sgesvx_obj->inforef = SGESVX( sgesvx_obj->matrix_layout, sgesvx_obj->fact,
                                  sgesvx_obj->trans, sgesvx_obj->n,
                                  sgesvx_obj->nrhs,
                                  sgesvx_obj->aref, sgesvx_obj->lda, 
                                  sgesvx_obj->afref, sgesvx_obj->ldaf,
                                  sgesvx_obj->ipivref, &sgesvx_obj->equedref,
                                  sgesvx_obj->rref, sgesvx_obj->cref,                                 
                                  sgesvx_obj->bref, sgesvx_obj->ldb,
                                  sgesvx_obj->xref, sgesvx_obj->ldx,
                                  &sgesvx_obj->rcondref, 
                                  sgesvx_obj->ferrref,
                                  sgesvx_obj->berrref,
                                  &sgesvx_obj->rpivotref);

    /* Compute libflame's Lapacke o/p  */
    sgesvx_obj->info = LAPACKE_sgesvx( sgesvx_obj->matrix_layout, sgesvx_obj->fact,
                                  sgesvx_obj->trans, sgesvx_obj->n,
                                  sgesvx_obj->nrhs,
                                  sgesvx_obj->a, sgesvx_obj->lda, 
                                  sgesvx_obj->af, sgesvx_obj->ldaf,
                                  sgesvx_obj->ipiv, &sgesvx_obj->equed,
                                  sgesvx_obj->r, sgesvx_obj->c,                               
                                  sgesvx_obj->b, sgesvx_obj->ldb,
                                  sgesvx_obj->x, sgesvx_obj->ldx,
                                  &sgesvx_obj->rcond, 
                                  sgesvx_obj->ferr,
                                  sgesvx_obj->berr,
                                  &sgesvx_obj->rpivot);

    if( sgesvx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_sgesvx is wrong\n", sgesvx_obj->info );
    }
    if( sgesvx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sgesvx is wrong\n", 
        sgesvx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    sgesvx_obj->diff =  computeDiff_s( sgesvx_obj->b_bufsize, 
                sgesvx_obj->b, sgesvx_obj->bref );

    sgesvx_obj->diff_xerr =  computeDiff_s( sgesvx_obj->x_bufsize, 
                sgesvx_obj->x, sgesvx_obj->xref );

    sgesvx_obj->diff_berr =  computeDiff_s( sgesvx_obj->nrhs, 
                sgesvx_obj->berr, sgesvx_obj->berrref );
                
    sgesvx_obj->diff_ferr =  computeDiff_s( sgesvx_obj->nrhs, 
                sgesvx_obj->ferr, sgesvx_obj->ferrref );
}

TEST_F(sgesvx_test, sgesvx1) {
    EXPECT_NEAR(0.0, sgesvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgesvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgesvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgesvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (sgesvx_obj->rpivot - sgesvx_obj->rpivotref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgesvx_test, sgesvx2) {
    EXPECT_NEAR(0.0, sgesvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgesvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgesvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgesvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (sgesvx_obj->rpivot - sgesvx_obj->rpivotref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgesvx_test, sgesvx3) {
    EXPECT_NEAR(0.0, sgesvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgesvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgesvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgesvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (sgesvx_obj->rpivot - sgesvx_obj->rpivotref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgesvx_test, sgesvx4) {
    EXPECT_NEAR(0.0, sgesvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgesvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgesvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgesvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (sgesvx_obj->rpivot - sgesvx_obj->rpivotref),
                LAPACKE_GTEST_THRESHOLD);
}

/* Begin gesvx_scomplex_parameters  class definition */
class gesvx_scomplex_parameters{
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
      float* r, *rref; // the row scale factors for A.
      float* c, *cref; // column scale factors for A.
      lapack_int *ipiv, *ipivref; // pivot buffer
      
      /* Output parameters */
      char   equed, equedref; //  Must be 'N', 'R', 'C', or 'B'.
      float rpivot, rpivotref; // reciprocal pivot growth factor.
      float rcond, rcondref; // estimate of the reciprocal condition number of the matrix A.
      lapack_complex_float* x, *xref; // solution matrix
      float* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      float* berr, *berrref; // component-wise relative backward error for each solution vector xj.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      gesvx_scomplex_parameters ( int matrix_layout_i, char fact_i, char trans_i,
                               char equed_i, lapack_int n_i, lapack_int nrhs_i);
              
      ~gesvx_scomplex_parameters (); 
};  /* end of gesvx_scomplex_parameters  class definition */


/* Constructor gesvx_scomplex_parameters definition */
gesvx_scomplex_parameters:: gesvx_scomplex_parameters ( int matrix_layout_i, 
                char fact_i, char trans_i, char equed_i, lapack_int n_i,
                  lapack_int nrhs_i) {
                                    
    matrix_layout = matrix_layout_i;
    fact= fact_i;
    trans = trans_i;
    // equed acts as input when " fact = 'F' " else output
    equed = equed_i;
    equedref = equed_i;
    n = n_i;
    nrhs = nrhs_i;
    lda = n;
    ldaf = n;
    
#if LAPACKE_TEST_VERBOSE
   printf(" \n gesvx Double:  n: %d, fact: %c trans: %c  lda: %d  \
ldb: %d nrhs: %d  ldaf: %d  ldx: %d \n",  n, fact, trans, lda, 
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
    lapacke_gtest_alloc_float_buffer_pair( &r, &rref, n);
    lapacke_gtest_alloc_float_buffer_pair( &c, &cref, n);
    lapacke_gtest_alloc_float_buffer_pair( &ferr, &ferrref, nrhs);
    lapacke_gtest_alloc_float_buffer_pair( &berr, &berrref, nrhs);
    
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (x==NULL) || (xref==NULL) || \
        (af==NULL) || (afref==NULL) || \
        (r==NULL) || (rref==NULL) || \
        (c==NULL) || (cref==NULL) || \
        (ferr==NULL) || (ferrref==NULL) || \
        (berr==NULL) || (berrref==NULL) || \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "pbtrs_scomplex_parameters object: malloc error.";
       gesvx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, n*n);
    memcpy(af, a, (n*n*sizeof(lapack_complex_float)));
    memcpy(afref, a, (n*n*sizeof(lapack_complex_float)));
    
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
    
    lapacke_gtest_init_scomplex_buffer_pair_with_constant(x, xref, x_bufsize, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(r, rref, n, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(c, cref, n, 0.0);

    
   } /* end of Constructor  */

gesvx_scomplex_parameters:: ~gesvx_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gesvx_scomplex_parameters object: destructor invoked. \n");
#endif
   gesvx_free();
}

//  Test fixture class definition
class cgesvx_test  : public  ::testing::Test {
public:
   gesvx_scomplex_parameters  *cgesvx_obj;
   void SetUp();  
   void TearDown () { delete cgesvx_obj; }
};


void cgesvx_test::SetUp(){

    /* LAPACKE CGESVX prototype */
    typedef int (*Fptr_NL_LAPACKE_cgesvx) (int matrix_layout,
                                            char fact,
                                            char trans,
                                            lapack_int n,
                                            lapack_int nrhs,
                                            lapack_complex_float* a,
                                            lapack_int lda,
                                            lapack_complex_float* af,
                                            lapack_int ldaf,
                                            lapack_int* ipiv,
                                            char* equed,
                                            float* r,
                                            float* c,
                                            lapack_complex_float* b,
                                            lapack_int ldb,
                                            lapack_complex_float* x,
                                            lapack_int ldx,
                                            float* rcond,
                                            float* ferr,
                                            float* berr,
                                            float* rpivot);

    Fptr_NL_LAPACKE_cgesvx CGESVX;


     /* LAPACKE CGETRF prototype */
    typedef int (*Fptr_NL_LAPACKE_cgetrf) ( int matrix_layout,lapack_int m,lapack_int n,
                                    lapack_complex_float* a,lapack_int lda,lapack_int* ipiv );

    Fptr_NL_LAPACKE_cgetrf CGETRF;


    cgesvx_obj = new gesvx_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].fact,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].equed,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    cgesvx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgesvx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgesvx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgesvx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    CGESVX = (Fptr_NL_LAPACKE_cgesvx)dlsym(cgesvx_obj->hModule, "LAPACKE_cgesvx");
    ASSERT_TRUE(CGESVX != NULL) << "failed to get the Netlib LAPACKE_cgesvx symbol";
    
    /* if  'fact' == 'F', then the i/p matrix A should be in 
       factored form.
       so, invoke the getrf API to compute the factorized A.  */
    if(cgesvx_obj->fact == 'F') {
        CGETRF = (Fptr_NL_LAPACKE_cgetrf)dlsym(cgesvx_obj->hModule,"LAPACKE_cgetrf");
        ASSERT_TRUE(CGETRF != NULL) << "failed to get the Netlib LAPACKE_cgetrf symbol";
            

        cgesvx_obj->inforef = CGETRF( cgesvx_obj->matrix_layout,
                                      cgesvx_obj->n, cgesvx_obj->n,
                                      cgesvx_obj->afref,
                                      cgesvx_obj->lda,
                                      cgesvx_obj->ipivref);
                               
        cgesvx_obj->info = LAPACKE_cgetrf( cgesvx_obj->matrix_layout,
                                           cgesvx_obj->n, cgesvx_obj->n,
                                           cgesvx_obj->af,
                                           cgesvx_obj->lda,
                                           cgesvx_obj->ipiv);
    }


    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    cgesvx_obj->inforef = CGESVX( cgesvx_obj->matrix_layout, cgesvx_obj->fact,
                                  cgesvx_obj->trans, cgesvx_obj->n,
                                  cgesvx_obj->nrhs,
                                  cgesvx_obj->aref, cgesvx_obj->lda, 
                                  cgesvx_obj->afref, cgesvx_obj->ldaf,
                                  cgesvx_obj->ipivref, &cgesvx_obj->equedref,
                                  cgesvx_obj->rref, cgesvx_obj->cref,                                 
                                  cgesvx_obj->bref, cgesvx_obj->ldb,
                                  cgesvx_obj->xref, cgesvx_obj->ldx,
                                  &cgesvx_obj->rcondref, 
                                  cgesvx_obj->ferrref,
                                  cgesvx_obj->berrref,
                                  &cgesvx_obj->rpivotref);

    /* Compute libflame's Lapacke o/p  */
    cgesvx_obj->info = LAPACKE_cgesvx( cgesvx_obj->matrix_layout, cgesvx_obj->fact,
                                  cgesvx_obj->trans, cgesvx_obj->n,
                                  cgesvx_obj->nrhs,
                                  cgesvx_obj->a, cgesvx_obj->lda, 
                                  cgesvx_obj->af, cgesvx_obj->ldaf,
                                  cgesvx_obj->ipiv, &cgesvx_obj->equed,
                                  cgesvx_obj->r, cgesvx_obj->c,                               
                                  cgesvx_obj->b, cgesvx_obj->ldb,
                                  cgesvx_obj->x, cgesvx_obj->ldx,
                                  &cgesvx_obj->rcond, 
                                  cgesvx_obj->ferr,
                                  cgesvx_obj->berr,
                                  &cgesvx_obj->rpivot);

    if( cgesvx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cgesvx is wrong\n", cgesvx_obj->info );
    }
    if( cgesvx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cgesvx is wrong\n", 
        cgesvx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cgesvx_obj->diff =  computeDiff_c( cgesvx_obj->b_bufsize, 
                cgesvx_obj->b, cgesvx_obj->bref );

    cgesvx_obj->diff_xerr =  computeDiff_c( cgesvx_obj->x_bufsize, 
                cgesvx_obj->x, cgesvx_obj->xref );

    cgesvx_obj->diff_berr =  computeDiff_s( cgesvx_obj->nrhs, 
                cgesvx_obj->berr, cgesvx_obj->berrref );
                
    cgesvx_obj->diff_ferr =  computeDiff_s( cgesvx_obj->nrhs, 
                cgesvx_obj->ferr, cgesvx_obj->ferrref );
}

TEST_F(cgesvx_test, cgesvx1) {
    EXPECT_NEAR(0.0, cgesvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgesvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgesvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgesvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (cgesvx_obj->rpivot - cgesvx_obj->rpivotref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgesvx_test, cgesvx2) {
    EXPECT_NEAR(0.0, cgesvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgesvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgesvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgesvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (cgesvx_obj->rpivot - cgesvx_obj->rpivotref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgesvx_test, cgesvx3) {
    EXPECT_NEAR(0.0, cgesvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgesvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgesvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgesvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (cgesvx_obj->rpivot - cgesvx_obj->rpivotref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgesvx_test, cgesvx4) {
    EXPECT_NEAR(0.0, cgesvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgesvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgesvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgesvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (cgesvx_obj->rpivot - cgesvx_obj->rpivotref),
                LAPACKE_GTEST_THRESHOLD);
}

/* Begin gesvx_dcomplex_parameters  class definition */
class gesvx_dcomplex_parameters{
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
      double* r, *rref; // the row scale factors for A.
      double* c, *cref; // column scale factors for A.
      lapack_int *ipiv, *ipivref; // pivot buffer
      
      /* Output parameters */
      char   equed, equedref; //  Must be 'N', 'R', 'C', or 'B'.
      double rpivot, rpivotref; // reciprocal pivot growth factor.
      double rcond, rcondref; // estimate of the reciprocal condition number of the matrix A.
      lapack_complex_double* x, *xref; // solution matrix
      double* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      double* berr, *berrref; // component-wise relative backward error for each solution vector xj.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      gesvx_dcomplex_parameters ( int matrix_layout_i, char fact_i, char trans_i,
                               char equed_i, lapack_int n_i, lapack_int nrhs_i);
              
      ~gesvx_dcomplex_parameters (); 
};  /* end of gesvx_dcomplex_parameters  class definition */


/* Constructor gesvx_dcomplex_parameters definition */
gesvx_dcomplex_parameters:: gesvx_dcomplex_parameters ( int matrix_layout_i, 
                char fact_i, char trans_i, char equed_i, lapack_int n_i,
                  lapack_int nrhs_i) {
                                    
    matrix_layout = matrix_layout_i;
    fact= fact_i;
    trans = trans_i;
    // equed acts as input when " fact = 'F' " else output
    equed = equed_i;
    equedref = equed_i;
    n = n_i;
    nrhs = nrhs_i;
    lda = n;
    ldaf = n;
    
#if LAPACKE_TEST_VERBOSE
   printf(" \n gesvx Double:  n: %d, fact: %c trans: %c  lda: %d  \
ldb: %d nrhs: %d  ldaf: %d  ldx: %d \n",  n, fact, trans, lda, 
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
    lapacke_gtest_alloc_double_buffer_pair( &r, &rref, n);
    lapacke_gtest_alloc_double_buffer_pair( &c, &cref, n);
    lapacke_gtest_alloc_double_buffer_pair( &ferr, &ferrref, nrhs);
    lapacke_gtest_alloc_double_buffer_pair( &berr, &berrref, nrhs);
    
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (x==NULL) || (xref==NULL) || \
        (af==NULL) || (afref==NULL) || \
        (r==NULL) || (rref==NULL) || \
        (c==NULL) || (cref==NULL) || \
        (ferr==NULL) || (ferrref==NULL) || \
        (berr==NULL) || (berrref==NULL) || \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "pbtrs_dcomplex_parameters object: malloc error.";
       gesvx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, n*n);
    memcpy(af, a, (n*n*sizeof(lapack_complex_double)));
    memcpy(afref, a, (n*n*sizeof(lapack_complex_double)));
    
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
    
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant(x, xref, x_bufsize, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(r, rref, n, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(c, cref, n, 0.0);

    
   } /* end of Constructor  */

gesvx_dcomplex_parameters:: ~gesvx_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gesvx_dcomplex_parameters object: destructor invoked. \n");
#endif
   gesvx_free();
}

//  Test fixture class definition
class zgesvx_test  : public  ::testing::Test {
public:
   gesvx_dcomplex_parameters  *zgesvx_obj;
   void SetUp();  
   void TearDown () { delete zgesvx_obj; }
};


void zgesvx_test::SetUp(){

    /* LAPACKE ZGESVX prototype */
    typedef int (*Fptr_NL_LAPACKE_zgesvx) (int matrix_layout,
                                            char fact,
                                            char trans,
                                            lapack_int n,
                                            lapack_int nrhs,
                                            lapack_complex_double* a,
                                            lapack_int lda,
                                            lapack_complex_double* af,
                                            lapack_int ldaf,
                                            lapack_int* ipiv,
                                            char* equed,
                                            double* r,
                                            double* c,
                                            lapack_complex_double* b,
                                            lapack_int ldb,
                                            lapack_complex_double* x,
                                            lapack_int ldx,
                                            double* rcond,
                                            double* ferr,
                                            double* berr,
                                            double* rpivot);

    Fptr_NL_LAPACKE_zgesvx ZGESVX;


     /* LAPACKE ZGETRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zgetrf) ( int matrix_layout,lapack_int m,lapack_int n,
                                    lapack_complex_double* a,lapack_int lda,lapack_int* ipiv );

    Fptr_NL_LAPACKE_zgetrf ZGETRF;


    zgesvx_obj = new gesvx_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].fact,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].equed,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    zgesvx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgesvx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgesvx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgesvx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    ZGESVX = (Fptr_NL_LAPACKE_zgesvx)dlsym(zgesvx_obj->hModule, "LAPACKE_zgesvx");
    ASSERT_TRUE(ZGESVX != NULL) << "failed to get the Netlib LAPACKE_zgesvx symbol";
    
    /* if  'fact' == 'F', then the i/p matrix A should be in 
       factored form.
       so, invoke the getrf API to compute the factorized A.  */
    if(zgesvx_obj->fact == 'F') {
        ZGETRF = (Fptr_NL_LAPACKE_zgetrf)dlsym(zgesvx_obj->hModule,"LAPACKE_zgetrf");
        ASSERT_TRUE(ZGETRF != NULL) << "failed to get the Netlib LAPACKE_zgetrf symbol";
            

        zgesvx_obj->inforef = ZGETRF( zgesvx_obj->matrix_layout,
                                      zgesvx_obj->n, zgesvx_obj->n,
                                      zgesvx_obj->afref,
                                      zgesvx_obj->lda,
                                      zgesvx_obj->ipivref);
                               
        zgesvx_obj->info = LAPACKE_zgetrf( zgesvx_obj->matrix_layout,
                                           zgesvx_obj->n, zgesvx_obj->n,
                                           zgesvx_obj->af,
                                           zgesvx_obj->lda,
                                           zgesvx_obj->ipiv);
    }


    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    zgesvx_obj->inforef = ZGESVX( zgesvx_obj->matrix_layout, zgesvx_obj->fact,
                                  zgesvx_obj->trans, zgesvx_obj->n,
                                  zgesvx_obj->nrhs,
                                  zgesvx_obj->aref, zgesvx_obj->lda, 
                                  zgesvx_obj->afref, zgesvx_obj->ldaf,
                                  zgesvx_obj->ipivref, &zgesvx_obj->equedref,
                                  zgesvx_obj->rref, zgesvx_obj->cref,                                 
                                  zgesvx_obj->bref, zgesvx_obj->ldb,
                                  zgesvx_obj->xref, zgesvx_obj->ldx,
                                  &zgesvx_obj->rcondref, 
                                  zgesvx_obj->ferrref,
                                  zgesvx_obj->berrref,
                                  &zgesvx_obj->rpivotref);

    /* Compute libflame's Lapacke o/p  */
    zgesvx_obj->info = LAPACKE_zgesvx( zgesvx_obj->matrix_layout, zgesvx_obj->fact,
                                  zgesvx_obj->trans, zgesvx_obj->n,
                                  zgesvx_obj->nrhs,
                                  zgesvx_obj->a, zgesvx_obj->lda, 
                                  zgesvx_obj->af, zgesvx_obj->ldaf,
                                  zgesvx_obj->ipiv, &zgesvx_obj->equed,
                                  zgesvx_obj->r, zgesvx_obj->c,                               
                                  zgesvx_obj->b, zgesvx_obj->ldb,
                                  zgesvx_obj->x, zgesvx_obj->ldx,
                                  &zgesvx_obj->rcond, 
                                  zgesvx_obj->ferr,
                                  zgesvx_obj->berr,
                                  &zgesvx_obj->rpivot);

    if( zgesvx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zgesvx is wrong\n", zgesvx_obj->info );
    }
    if( zgesvx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zgesvx is wrong\n", 
        zgesvx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zgesvx_obj->diff =  computeDiff_z( zgesvx_obj->b_bufsize, 
                zgesvx_obj->b, zgesvx_obj->bref );

    zgesvx_obj->diff_xerr =  computeDiff_z( zgesvx_obj->x_bufsize, 
                zgesvx_obj->x, zgesvx_obj->xref );

    zgesvx_obj->diff_berr =  computeDiff_d( zgesvx_obj->nrhs, 
                zgesvx_obj->berr, zgesvx_obj->berrref );
                
    zgesvx_obj->diff_ferr =  computeDiff_d( zgesvx_obj->nrhs, 
                zgesvx_obj->ferr, zgesvx_obj->ferrref );
}

TEST_F(zgesvx_test, zgesvx1) {
    EXPECT_NEAR(0.0, zgesvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgesvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgesvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgesvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (zgesvx_obj->rpivot - zgesvx_obj->rpivotref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgesvx_test, zgesvx2) {
    EXPECT_NEAR(0.0, zgesvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgesvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgesvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgesvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (zgesvx_obj->rpivot - zgesvx_obj->rpivotref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgesvx_test, zgesvx3) {
    EXPECT_NEAR(0.0, zgesvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgesvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgesvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgesvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (zgesvx_obj->rpivot - zgesvx_obj->rpivotref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgesvx_test, zgesvx4) {
    EXPECT_NEAR(0.0, zgesvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgesvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgesvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgesvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (zgesvx_obj->rpivot - zgesvx_obj->rpivotref),
                LAPACKE_GTEST_THRESHOLD);
}
