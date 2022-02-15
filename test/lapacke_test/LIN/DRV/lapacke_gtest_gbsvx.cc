#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"

#define gbsvx_free() \
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

/* Begin gbsvx_float_parameters  class definition */
class gbsvx_float_parameters{
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
      gbsvx_float_parameters ( int matrix_layout_i, char fact_i, char trans_i,
                                 lapack_int kl_i, lapack_int ku_i,
                               char equed_i, lapack_int n_i, lapack_int nrhs_i);
              
      ~gbsvx_float_parameters (); 
};  /* end of gbsvx_float_parameters  class definition */


/* Constructor gbsvx_float_parameters definition */
gbsvx_float_parameters:: gbsvx_float_parameters ( int matrix_layout_i, 
                                            char fact_i, char trans_i,
                                     lapack_int kl_i, lapack_int ku_i,
                    char equed_i, lapack_int n_i, lapack_int nrhs_i) {
                                    
    matrix_layout = matrix_layout_i;
    fact= fact_i;
    trans = trans_i;
    kl = kl_i;
    ku = ku_i;
    // equed acts as input when " fact = 'F' " else output
    equed = equed_i;
    equedref = equed_i;
    n = n_i;
    nrhs = nrhs_i;
    lda = n;
    ldaf = n;
    
#if LAPACKE_TEST_VERBOSE
   printf(" \n gbsvx float:  n: %d, fact: %c trans: %c  lda: %d  \
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
       gbsvx_free();
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

gbsvx_float_parameters:: ~gbsvx_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gbsvx_float_parameters object: destructor invoked. \n");
#endif
   gbsvx_free();
}

//  Test fixture class definition
class sgbsvx_test  : public  ::testing::Test {
public:
   gbsvx_float_parameters  *sgbsvx_obj;
   void SetUp();  
   void TearDown () { delete sgbsvx_obj; }
};


void sgbsvx_test::SetUp(){

    /* LAPACKE SGBSVX prototype */
    typedef int (*Fptr_NL_LAPACKE_sgbsvx) (int matrix_layout, char fact,
                            char trans, lapack_int n, lapack_int kl, 
                            lapack_int ku,lapack_int nrhs, float *a,
                            lapack_int lda, float *af, lapack_int ldaf,
                            lapack_int *ipiv, char *equed, float *r,
                            float *c, float *b, lapack_int ldb,
                            float *x, lapack_int ldx, float *rcond,
                            float *ferr, float *berr, float *rpivot);

    Fptr_NL_LAPACKE_sgbsvx SGBSVX;

     /* LAPACKE SGBTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_sgbtrf) ( int matrix_layout,lapack_int m,
                                lapack_int n, lapack_int kl, lapack_int ku,
                              float *ab,lapack_int ldab,lapack_int* ipiv );

    Fptr_NL_LAPACKE_sgbtrf SGBTRF;

    sgbsvx_obj = new gbsvx_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].fact,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].kl,
                           lin_solver_paramslist[idx].ku,
                           lin_solver_paramslist[idx].equed,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    sgbsvx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgbsvx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgbsvx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgbsvx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    SGBSVX = (Fptr_NL_LAPACKE_sgbsvx)dlsym(sgbsvx_obj->hModule, "LAPACKE_sgbsvx");
    ASSERT_TRUE(SGBSVX != NULL) << "failed to get the Netlib LAPACKE_sgbsvx symbol";
    
    SGBTRF = (Fptr_NL_LAPACKE_sgbtrf)dlsym(sgbsvx_obj->hModule,"LAPACKE_sgbtrf");
    ASSERT_TRUE(SGBTRF != NULL) << "failed to get the Netlib LAPACKE_sgbtrf symbol";
    
    /* if  'fact' == 'F', then the i/p matrix A should be in 
       factored form.
       so, invoke the gbtrf API to compute the factorized A.  */
    if(sgbsvx_obj->fact == 'F') {
        sgbsvx_obj->inforef = SGBTRF( sgbsvx_obj->matrix_layout,
                                      sgbsvx_obj->n, sgbsvx_obj->n,
                                      sgbsvx_obj->kl, sgbsvx_obj->ku,
                                      sgbsvx_obj->afref,
                                      sgbsvx_obj->lda,
                                      sgbsvx_obj->ipivref);
                               
        sgbsvx_obj->info = LAPACKE_sgbtrf( sgbsvx_obj->matrix_layout,
                                        sgbsvx_obj->n, sgbsvx_obj->n,
                                      sgbsvx_obj->kl, sgbsvx_obj->ku,
                                           sgbsvx_obj->af,
                                           sgbsvx_obj->lda,
                                           sgbsvx_obj->ipiv);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    sgbsvx_obj->inforef = SGBSVX( sgbsvx_obj->matrix_layout, sgbsvx_obj->fact,
                                  sgbsvx_obj->trans, sgbsvx_obj->n,
                                  sgbsvx_obj->kl, sgbsvx_obj->ku,
                                  sgbsvx_obj->nrhs,
                                  sgbsvx_obj->aref, sgbsvx_obj->lda, 
                                  sgbsvx_obj->afref, sgbsvx_obj->ldaf,
                                  sgbsvx_obj->ipivref, &sgbsvx_obj->equedref,
                                  sgbsvx_obj->rref, sgbsvx_obj->cref,                                 
                                  sgbsvx_obj->bref, sgbsvx_obj->ldb,
                                  sgbsvx_obj->xref, sgbsvx_obj->ldx,
                                  &sgbsvx_obj->rcondref, 
                                  sgbsvx_obj->ferrref,
                                  sgbsvx_obj->berrref,
                                  &sgbsvx_obj->rpivotref);

    /* Compute libflame's Lapacke o/p  */
    sgbsvx_obj->info = LAPACKE_sgbsvx( sgbsvx_obj->matrix_layout, sgbsvx_obj->fact,
                                  sgbsvx_obj->trans, sgbsvx_obj->n,
                                  sgbsvx_obj->kl, sgbsvx_obj->ku,
                                  sgbsvx_obj->nrhs,
                                  sgbsvx_obj->a, sgbsvx_obj->lda, 
                                  sgbsvx_obj->af, sgbsvx_obj->ldaf,
                                  sgbsvx_obj->ipiv, &sgbsvx_obj->equed,
                                  sgbsvx_obj->r, sgbsvx_obj->c,                               
                                  sgbsvx_obj->b, sgbsvx_obj->ldb,
                                  sgbsvx_obj->x, sgbsvx_obj->ldx,
                                  &sgbsvx_obj->rcond, 
                                  sgbsvx_obj->ferr,
                                  sgbsvx_obj->berr,
                                  &sgbsvx_obj->rpivot);

    if( sgbsvx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_sgbsvx is wrong\n", sgbsvx_obj->info );
    }
    if( sgbsvx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sgbsvx is wrong\n", 
        sgbsvx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    sgbsvx_obj->diff =  computeDiff_s( sgbsvx_obj->b_bufsize, 
                sgbsvx_obj->b, sgbsvx_obj->bref );

    sgbsvx_obj->diff_xerr =  computeDiff_s( sgbsvx_obj->x_bufsize, 
                sgbsvx_obj->x, sgbsvx_obj->xref );

    sgbsvx_obj->diff_berr =  computeDiff_s( sgbsvx_obj->nrhs, 
                sgbsvx_obj->berr, sgbsvx_obj->berrref );
                
    sgbsvx_obj->diff_ferr =  computeDiff_s( sgbsvx_obj->nrhs, 
                sgbsvx_obj->ferr, sgbsvx_obj->ferrref );
}

TEST_F(sgbsvx_test, sgbsvx1) {
    EXPECT_NEAR(0.0, sgbsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgbsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgbsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgbsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (sgbsvx_obj->rpivot - sgbsvx_obj->rpivotref),
                     LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgbsvx_test, sgbsvx2) {
    EXPECT_NEAR(0.0, sgbsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgbsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgbsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgbsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (sgbsvx_obj->rpivot - sgbsvx_obj->rpivotref),
                     LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgbsvx_test, sgbsvx3) {
    EXPECT_NEAR(0.0, sgbsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgbsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgbsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgbsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (sgbsvx_obj->rpivot - sgbsvx_obj->rpivotref),
                     LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgbsvx_test, sgbsvx4) {
    EXPECT_NEAR(0.0, sgbsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgbsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgbsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgbsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (sgbsvx_obj->rpivot - sgbsvx_obj->rpivotref),
                     LAPACKE_GTEST_THRESHOLD);
}

/* Begin gbsvx_double_parameters  class definition */
class gbsvx_double_parameters{
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
      gbsvx_double_parameters ( int matrix_layout_i, char fact_i, char trans_i,
                                 lapack_int kl_i, lapack_int ku_i,
                               char equed_i, lapack_int n_i, lapack_int nrhs_i);
              
      ~gbsvx_double_parameters (); 
};  /* end of gbsvx_double_parameters  class definition */


/* Constructor gbsvx_double_parameters definition */
gbsvx_double_parameters:: gbsvx_double_parameters ( int matrix_layout_i, 
                                            char fact_i, char trans_i,
                                     lapack_int kl_i, lapack_int ku_i,
                    char equed_i, lapack_int n_i, lapack_int nrhs_i) {
                                    
    matrix_layout = matrix_layout_i;
    fact= fact_i;
    trans = trans_i;
    kl = kl_i;
    ku = ku_i;
    // equed acts as input when " fact = 'F' " else output
    equed = equed_i;
    equedref = equed_i;
    n = n_i;
    nrhs = nrhs_i;
    lda = n;
    ldaf = n;
    
#if LAPACKE_TEST_VERBOSE
   printf(" \n gbsvx double:  n: %d, fact: %c trans: %c  lda: %d  \
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
       gbsvx_free();
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

gbsvx_double_parameters:: ~gbsvx_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gbsvx_double_parameters object: destructor invoked. \n");
#endif
   gbsvx_free();
}

//  Test fixture class definition
class dgbsvx_test  : public  ::testing::Test {
public:
   gbsvx_double_parameters  *dgbsvx_obj;
   void SetUp();  
   void TearDown () { delete dgbsvx_obj; }
};


void dgbsvx_test::SetUp(){

    /* LAPACKE DGBSVX prototype */
    typedef int (*Fptr_NL_LAPACKE_dgbsvx) (int matrix_layout, char fact,
                            char trans, lapack_int n, lapack_int kl, 
                            lapack_int ku,lapack_int nrhs, double *a,
                            lapack_int lda, double *af, lapack_int ldaf,
                            lapack_int *ipiv, char *equed, double *r,
                            double *c, double *b, lapack_int ldb,
                            double *x, lapack_int ldx, double *rcond,
                            double *ferr, double *berr, double *rpivot);

    Fptr_NL_LAPACKE_dgbsvx DGBSVX;


     /* LAPACKE DGBTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_dgbtrf) ( int matrix_layout,lapack_int m,
                                lapack_int n, lapack_int kl, lapack_int ku,
                                    double* a,lapack_int lda,lapack_int* ipiv );

    Fptr_NL_LAPACKE_dgbtrf DGBTRF;


    dgbsvx_obj = new gbsvx_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].fact,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].kl,
                           lin_solver_paramslist[idx].ku,
                           lin_solver_paramslist[idx].equed,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    dgbsvx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgbsvx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgbsvx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgbsvx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    DGBSVX = (Fptr_NL_LAPACKE_dgbsvx)dlsym(dgbsvx_obj->hModule, "LAPACKE_dgbsvx");
    ASSERT_TRUE(DGBSVX != NULL) << "failed to get the Netlib LAPACKE_dgbsvx symbol";
    DGBTRF = (Fptr_NL_LAPACKE_dgbtrf)dlsym(dgbsvx_obj->hModule,"LAPACKE_dgbtrf");
    ASSERT_TRUE(DGBTRF != NULL) << "failed to get the Netlib LAPACKE_dgbtrf symbol";
    
    /* if  'fact' == 'F', then the i/p matrix A should be in 
       factored form.
       so, invoke the gbtrf API to compute the factorized A.  */
    if(dgbsvx_obj->fact == 'F') {
        dgbsvx_obj->inforef = DGBTRF( dgbsvx_obj->matrix_layout,
                                      dgbsvx_obj->n, dgbsvx_obj->n,
                                      dgbsvx_obj->kl, dgbsvx_obj->ku,
                                      dgbsvx_obj->afref,
                                      dgbsvx_obj->lda,
                                      dgbsvx_obj->ipivref);
                               
        dgbsvx_obj->info = LAPACKE_dgbtrf( dgbsvx_obj->matrix_layout,
                                           dgbsvx_obj->n, dgbsvx_obj->n,
                                      dgbsvx_obj->kl, dgbsvx_obj->ku,
                                           dgbsvx_obj->af,
                                           dgbsvx_obj->lda,
                                           dgbsvx_obj->ipiv);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    dgbsvx_obj->inforef = DGBSVX( dgbsvx_obj->matrix_layout, dgbsvx_obj->fact,
                                  dgbsvx_obj->trans, dgbsvx_obj->n,
                                  dgbsvx_obj->kl, dgbsvx_obj->ku,
                                  dgbsvx_obj->nrhs,
                                  dgbsvx_obj->aref, dgbsvx_obj->lda, 
                                  dgbsvx_obj->afref, dgbsvx_obj->ldaf,
                                  dgbsvx_obj->ipivref, &dgbsvx_obj->equedref,
                                  dgbsvx_obj->rref, dgbsvx_obj->cref,                                 
                                  dgbsvx_obj->bref, dgbsvx_obj->ldb,
                                  dgbsvx_obj->xref, dgbsvx_obj->ldx,
                                  &dgbsvx_obj->rcondref, 
                                  dgbsvx_obj->ferrref,
                                  dgbsvx_obj->berrref,
                                  &dgbsvx_obj->rpivotref);

    /* Compute libflame's Lapacke o/p  */
    dgbsvx_obj->info = LAPACKE_dgbsvx( dgbsvx_obj->matrix_layout, dgbsvx_obj->fact,
                                  dgbsvx_obj->trans, dgbsvx_obj->n,
                                  dgbsvx_obj->kl, dgbsvx_obj->ku,
                                  dgbsvx_obj->nrhs,
                                  dgbsvx_obj->a, dgbsvx_obj->lda, 
                                  dgbsvx_obj->af, dgbsvx_obj->ldaf,
                                  dgbsvx_obj->ipiv, &dgbsvx_obj->equed,
                                  dgbsvx_obj->r, dgbsvx_obj->c,                               
                                  dgbsvx_obj->b, dgbsvx_obj->ldb,
                                  dgbsvx_obj->x, dgbsvx_obj->ldx,
                                  &dgbsvx_obj->rcond, 
                                  dgbsvx_obj->ferr,
                                  dgbsvx_obj->berr,
                                  &dgbsvx_obj->rpivot);

    if( dgbsvx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dgbsvx is wrong\n", dgbsvx_obj->info );
    }
    if( dgbsvx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dgbsvx is wrong\n", 
        dgbsvx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dgbsvx_obj->diff =  computeDiff_d( dgbsvx_obj->b_bufsize, 
                dgbsvx_obj->b, dgbsvx_obj->bref );

    dgbsvx_obj->diff_xerr =  computeDiff_d( dgbsvx_obj->x_bufsize, 
                dgbsvx_obj->x, dgbsvx_obj->xref );

    dgbsvx_obj->diff_berr =  computeDiff_d( dgbsvx_obj->nrhs, 
                dgbsvx_obj->berr, dgbsvx_obj->berrref );
                
    dgbsvx_obj->diff_ferr =  computeDiff_d( dgbsvx_obj->nrhs, 
                dgbsvx_obj->ferr, dgbsvx_obj->ferrref );
}

TEST_F(dgbsvx_test, dgbsvx1) {
    EXPECT_NEAR(0.0, dgbsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgbsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgbsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgbsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (dgbsvx_obj->rpivot - dgbsvx_obj->rpivotref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgbsvx_test, dgbsvx2) {
    EXPECT_NEAR(0.0, dgbsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgbsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgbsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgbsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (dgbsvx_obj->rpivot - dgbsvx_obj->rpivotref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgbsvx_test, dgbsvx3) {
    EXPECT_NEAR(0.0, dgbsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgbsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgbsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgbsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (dgbsvx_obj->rpivot - dgbsvx_obj->rpivotref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgbsvx_test, dgbsvx4) {
    EXPECT_NEAR(0.0, dgbsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgbsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgbsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgbsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (dgbsvx_obj->rpivot - dgbsvx_obj->rpivotref),
                LAPACKE_GTEST_THRESHOLD);
}

/* Begin gbsvx_scomplex_parameters  class definition */
class gbsvx_scomplex_parameters{
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
      gbsvx_scomplex_parameters ( int matrix_layout_i, char fact_i, char trans_i,
                                 lapack_int kl_i, lapack_int ku_i,
                               char equed_i, lapack_int n_i, lapack_int nrhs_i);
              
      ~gbsvx_scomplex_parameters (); 
};  /* end of gbsvx_scomplex_parameters  class definition */


/* Constructor gbsvx_scomplex_parameters definition */
gbsvx_scomplex_parameters:: gbsvx_scomplex_parameters ( int matrix_layout_i, 
                                            char fact_i, char trans_i,
                                     lapack_int kl_i, lapack_int ku_i,
                    char equed_i, lapack_int n_i, lapack_int nrhs_i) {
                                    
    matrix_layout = matrix_layout_i;
    fact= fact_i;
    trans = trans_i;
    kl = kl_i;
    ku = ku_i;
    // equed acts as input when " fact = 'F' " else output
    equed = equed_i;
    equedref = equed_i;
    n = n_i;
    nrhs = nrhs_i;
    lda = n;
    ldaf = n;
    
#if LAPACKE_TEST_VERBOSE
   printf(" \n gbsvx scomplex:  n: %d, fact: %c trans: %c  lda: %d  \
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
       gbsvx_free();
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

gbsvx_scomplex_parameters:: ~gbsvx_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gbsvx_scomplex_parameters object: destructor invoked. \n");
#endif
   gbsvx_free();
}

//  Test fixture class definition
class cgbsvx_test  : public  ::testing::Test {
public:
   gbsvx_scomplex_parameters  *cgbsvx_obj;
   void SetUp();  
   void TearDown () { delete cgbsvx_obj; }
};


void cgbsvx_test::SetUp(){

    /* LAPACKE CGBSVX prototype */
    typedef int (*Fptr_NL_LAPACKE_cgbsvx) (int matrix_layout,
                                            char fact,
                                            char trans,
                                            lapack_int n,
                                            lapack_int kl,
                                            lapack_int ku,
                                            lapack_int nrhs,
                                            lapack_complex_float *a,
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

    Fptr_NL_LAPACKE_cgbsvx CGBSVX;


     /* LAPACKE CGBTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_cgbtrf) ( int matrix_layout,lapack_int m,
                                lapack_int n, lapack_int kl, lapack_int ku,
                                lapack_complex_float *ab, lapack_int ldab,
                                lapack_int* ipiv );

    Fptr_NL_LAPACKE_cgbtrf CGBTRF;


    cgbsvx_obj = new gbsvx_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].fact,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].kl,
                           lin_solver_paramslist[idx].ku,
                           lin_solver_paramslist[idx].equed,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    cgbsvx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgbsvx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgbsvx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgbsvx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    CGBSVX = (Fptr_NL_LAPACKE_cgbsvx)dlsym(cgbsvx_obj->hModule, "LAPACKE_cgbsvx");
    ASSERT_TRUE(CGBSVX != NULL) << "failed to get the Netlib LAPACKE_cgbsvx symbol";
    CGBTRF = (Fptr_NL_LAPACKE_cgbtrf)dlsym(cgbsvx_obj->hModule,"LAPACKE_cgbtrf");
    ASSERT_TRUE(CGBTRF != NULL) << "failed to get the Netlib LAPACKE_cgbtrf symbol";
    
    /* if  'fact' == 'F', then the i/p matrix A should be in 
       factored form.
       so, invoke the getrf API to compute the factorized A.  */
    if(cgbsvx_obj->fact == 'F') {
        cgbsvx_obj->inforef = CGBTRF( cgbsvx_obj->matrix_layout,
                                      cgbsvx_obj->n, cgbsvx_obj->n,
                                      cgbsvx_obj->kl, cgbsvx_obj->ku,
                                      cgbsvx_obj->afref,
                                      cgbsvx_obj->lda,
                                      cgbsvx_obj->ipivref);
                               
        cgbsvx_obj->info = LAPACKE_cgbtrf( cgbsvx_obj->matrix_layout,
                                           cgbsvx_obj->n, cgbsvx_obj->n,
                                      cgbsvx_obj->kl, cgbsvx_obj->ku,
                                           cgbsvx_obj->af,
                                           cgbsvx_obj->lda,
                                           cgbsvx_obj->ipiv);
    }


    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    cgbsvx_obj->inforef = CGBSVX( cgbsvx_obj->matrix_layout, cgbsvx_obj->fact,
                                  cgbsvx_obj->trans, cgbsvx_obj->n,
                                      cgbsvx_obj->kl, cgbsvx_obj->ku,
                                  cgbsvx_obj->nrhs,
                                  cgbsvx_obj->aref, cgbsvx_obj->lda, 
                                  cgbsvx_obj->afref, cgbsvx_obj->ldaf,
                                  cgbsvx_obj->ipivref, &cgbsvx_obj->equedref,
                                  cgbsvx_obj->rref, cgbsvx_obj->cref,                                 
                                  cgbsvx_obj->bref, cgbsvx_obj->ldb,
                                  cgbsvx_obj->xref, cgbsvx_obj->ldx,
                                  &cgbsvx_obj->rcondref, 
                                  cgbsvx_obj->ferrref,
                                  cgbsvx_obj->berrref,
                                  &cgbsvx_obj->rpivotref);

    /* Compute libflame's Lapacke o/p  */
    cgbsvx_obj->info = LAPACKE_cgbsvx( cgbsvx_obj->matrix_layout, cgbsvx_obj->fact,
                                  cgbsvx_obj->trans, cgbsvx_obj->n,
                                      cgbsvx_obj->kl, cgbsvx_obj->ku,
                                  cgbsvx_obj->nrhs,
                                  cgbsvx_obj->a, cgbsvx_obj->lda, 
                                  cgbsvx_obj->af, cgbsvx_obj->ldaf,
                                  cgbsvx_obj->ipiv, &cgbsvx_obj->equed,
                                  cgbsvx_obj->r, cgbsvx_obj->c,                               
                                  cgbsvx_obj->b, cgbsvx_obj->ldb,
                                  cgbsvx_obj->x, cgbsvx_obj->ldx,
                                  &cgbsvx_obj->rcond, 
                                  cgbsvx_obj->ferr,
                                  cgbsvx_obj->berr,
                                  &cgbsvx_obj->rpivot);

    if( cgbsvx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cgbsvx is wrong\n", cgbsvx_obj->info );
    }
    if( cgbsvx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cgbsvx is wrong\n", 
        cgbsvx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cgbsvx_obj->diff =  computeDiff_c( cgbsvx_obj->b_bufsize, 
                cgbsvx_obj->b, cgbsvx_obj->bref );

    cgbsvx_obj->diff_xerr =  computeDiff_c( cgbsvx_obj->x_bufsize, 
                cgbsvx_obj->x, cgbsvx_obj->xref );

    cgbsvx_obj->diff_berr =  computeDiff_s( cgbsvx_obj->nrhs, 
                cgbsvx_obj->berr, cgbsvx_obj->berrref );
                
    cgbsvx_obj->diff_ferr =  computeDiff_s( cgbsvx_obj->nrhs, 
                cgbsvx_obj->ferr, cgbsvx_obj->ferrref );
}

TEST_F(cgbsvx_test, cgbsvx1) {
    EXPECT_NEAR(0.0, cgbsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgbsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgbsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgbsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (cgbsvx_obj->rpivot - cgbsvx_obj->rpivotref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgbsvx_test, cgbsvx2) {
    EXPECT_NEAR(0.0, cgbsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgbsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgbsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgbsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (cgbsvx_obj->rpivot - cgbsvx_obj->rpivotref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgbsvx_test, cgbsvx3) {
    EXPECT_NEAR(0.0, cgbsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgbsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgbsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgbsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (cgbsvx_obj->rpivot - cgbsvx_obj->rpivotref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgbsvx_test, cgbsvx4) {
    EXPECT_NEAR(0.0, cgbsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgbsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgbsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgbsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (cgbsvx_obj->rpivot - cgbsvx_obj->rpivotref),
                LAPACKE_GTEST_THRESHOLD);
}

/* Begin gbsvx_dcomplex_parameters  class definition */
class gbsvx_dcomplex_parameters{
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
      gbsvx_dcomplex_parameters ( int matrix_layout_i, char fact_i, char trans_i,
                                 lapack_int kl_i, lapack_int ku_i,
                               char equed_i, lapack_int n_i, lapack_int nrhs_i);
              
      ~gbsvx_dcomplex_parameters (); 
};  /* end of gbsvx_dcomplex_parameters  class definition */


/* Constructor gbsvx_dcomplex_parameters definition */
gbsvx_dcomplex_parameters:: gbsvx_dcomplex_parameters ( int matrix_layout_i, 
                                            char fact_i, char trans_i,
                                     lapack_int kl_i, lapack_int ku_i,
                    char equed_i, lapack_int n_i, lapack_int nrhs_i) {
                                    
    matrix_layout = matrix_layout_i;
    fact= fact_i;
    trans = trans_i;
    kl = kl_i;
    ku = ku_i;
    // equed acts as input when " fact = 'F' " else output
    equed = equed_i;
    equedref = equed_i;
    n = n_i;
    nrhs = nrhs_i;
    lda = n;
    ldaf = n;
    
#if LAPACKE_TEST_VERBOSE
   printf(" \n gbsvx dcomplex:  n: %d, fact: %c trans: %c  lda: %d  \
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
       gbsvx_free();
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

gbsvx_dcomplex_parameters:: ~gbsvx_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gbsvx_dcomplex_parameters object: destructor invoked. \n");
#endif
   gbsvx_free();
}

//  Test fixture class definition
class zgbsvx_test  : public  ::testing::Test {
public:
   gbsvx_dcomplex_parameters  *zgbsvx_obj;
   void SetUp();  
   void TearDown () { delete zgbsvx_obj; }
};


void zgbsvx_test::SetUp(){

    /* LAPACKE ZGBSVX prototype */
    typedef int (*Fptr_NL_LAPACKE_zgbsvx) (int matrix_layout,
                                            char fact,
                                            char trans,
                                            lapack_int n,
                                            lapack_int kl,
                                            lapack_int ku,
                                            lapack_int nrhs,
                                            lapack_complex_double *a,
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

    Fptr_NL_LAPACKE_zgbsvx ZGBSVX;


     /* LAPACKE ZGBTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zgbtrf) ( int matrix_layout,lapack_int m,
                                lapack_int n, lapack_int kl, lapack_int ku,
                                lapack_complex_double *ab, lapack_int ldab,
                                lapack_int* ipiv );

    Fptr_NL_LAPACKE_zgbtrf ZGBTRF;


    zgbsvx_obj = new gbsvx_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].fact,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].kl,
                           lin_solver_paramslist[idx].ku,
                           lin_solver_paramslist[idx].equed,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    zgbsvx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgbsvx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgbsvx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgbsvx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    ZGBSVX = (Fptr_NL_LAPACKE_zgbsvx)dlsym(zgbsvx_obj->hModule, "LAPACKE_zgbsvx");
    ASSERT_TRUE(ZGBSVX != NULL) << "failed to get the Netlib LAPACKE_zgbsvx symbol";
    ZGBTRF = (Fptr_NL_LAPACKE_zgbtrf)dlsym(zgbsvx_obj->hModule,"LAPACKE_zgbtrf");
    ASSERT_TRUE(ZGBTRF != NULL) << "failed to get the Netlib LAPACKE_zgbtrf symbol";
    
    /* if  'fact' == 'F', then the i/p matrix A should be in 
       factored form.
       so, invoke the getrf API to compute the factorized A.  */
    if(zgbsvx_obj->fact == 'F') {
        zgbsvx_obj->inforef = ZGBTRF( zgbsvx_obj->matrix_layout,
                                      zgbsvx_obj->n, zgbsvx_obj->n,
                                      zgbsvx_obj->kl, zgbsvx_obj->ku,
                                      zgbsvx_obj->afref,
                                      zgbsvx_obj->lda,
                                      zgbsvx_obj->ipivref);
                               
        zgbsvx_obj->info = LAPACKE_zgbtrf( zgbsvx_obj->matrix_layout,
                                           zgbsvx_obj->n, zgbsvx_obj->n,
                                      zgbsvx_obj->kl, zgbsvx_obj->ku,
                                           zgbsvx_obj->af,
                                           zgbsvx_obj->lda,
                                           zgbsvx_obj->ipiv);
    }


    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    zgbsvx_obj->inforef = ZGBSVX( zgbsvx_obj->matrix_layout, zgbsvx_obj->fact,
                                  zgbsvx_obj->trans, zgbsvx_obj->n,
                                      zgbsvx_obj->kl, zgbsvx_obj->ku,
                                  zgbsvx_obj->nrhs,
                                  zgbsvx_obj->aref, zgbsvx_obj->lda, 
                                  zgbsvx_obj->afref, zgbsvx_obj->ldaf,
                                  zgbsvx_obj->ipivref, &zgbsvx_obj->equedref,
                                  zgbsvx_obj->rref, zgbsvx_obj->cref,                                 
                                  zgbsvx_obj->bref, zgbsvx_obj->ldb,
                                  zgbsvx_obj->xref, zgbsvx_obj->ldx,
                                  &zgbsvx_obj->rcondref, 
                                  zgbsvx_obj->ferrref,
                                  zgbsvx_obj->berrref,
                                  &zgbsvx_obj->rpivotref);

    /* Compute libflame's Lapacke o/p  */
    zgbsvx_obj->info = LAPACKE_zgbsvx( zgbsvx_obj->matrix_layout, zgbsvx_obj->fact,
                                  zgbsvx_obj->trans, zgbsvx_obj->n,
                                      zgbsvx_obj->kl, zgbsvx_obj->ku,
                                  zgbsvx_obj->nrhs,
                                  zgbsvx_obj->a, zgbsvx_obj->lda, 
                                  zgbsvx_obj->af, zgbsvx_obj->ldaf,
                                  zgbsvx_obj->ipiv, &zgbsvx_obj->equed,
                                  zgbsvx_obj->r, zgbsvx_obj->c,                               
                                  zgbsvx_obj->b, zgbsvx_obj->ldb,
                                  zgbsvx_obj->x, zgbsvx_obj->ldx,
                                  &zgbsvx_obj->rcond, 
                                  zgbsvx_obj->ferr,
                                  zgbsvx_obj->berr,
                                  &zgbsvx_obj->rpivot);

    if( zgbsvx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zgbsvx is wrong\n", zgbsvx_obj->info );
    }
    if( zgbsvx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zgbsvx is wrong\n", 
        zgbsvx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zgbsvx_obj->diff =  computeDiff_z( zgbsvx_obj->b_bufsize, 
                zgbsvx_obj->b, zgbsvx_obj->bref );

    zgbsvx_obj->diff_xerr =  computeDiff_z( zgbsvx_obj->x_bufsize, 
                zgbsvx_obj->x, zgbsvx_obj->xref );

    zgbsvx_obj->diff_berr =  computeDiff_d( zgbsvx_obj->nrhs, 
                zgbsvx_obj->berr, zgbsvx_obj->berrref );
                
    zgbsvx_obj->diff_ferr =  computeDiff_d( zgbsvx_obj->nrhs, 
                zgbsvx_obj->ferr, zgbsvx_obj->ferrref );
}

TEST_F(zgbsvx_test, zgbsvx1) {
    EXPECT_NEAR(0.0, zgbsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgbsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgbsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgbsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (zgbsvx_obj->rpivot - zgbsvx_obj->rpivotref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgbsvx_test, zgbsvx2) {
    EXPECT_NEAR(0.0, zgbsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgbsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgbsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgbsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (zgbsvx_obj->rpivot - zgbsvx_obj->rpivotref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgbsvx_test, zgbsvx3) {
    EXPECT_NEAR(0.0, zgbsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgbsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgbsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgbsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (zgbsvx_obj->rpivot - zgbsvx_obj->rpivotref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgbsvx_test, zgbsvx4) {
    EXPECT_NEAR(0.0, zgbsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgbsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgbsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgbsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (zgbsvx_obj->rpivot - zgbsvx_obj->rpivotref),
                LAPACKE_GTEST_THRESHOLD);
}
