#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"

#define LAPACKE_TEST_VERBOSE (1)
#define gbtrs_free() \
       if (ipiv != NULL) free (ipiv); \
       if (ipivref != NULL)free (ipivref); \
       if (b != NULL)    free (b   ); \
       if (bref != NULL) free (bref); \
       if (ab != NULL)    free (ab  ); \
       if (abref != NULL) free (abref); \
       if( hModule != NULL) dlclose(hModule); \
       if(dModule != NULL) dlclose(dModule)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;


/* Begin gbtrs_double_parameters  class definition */
class gbtrs_double_parameters{
   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char trans; // 'N' or 'T' or 'C'. Indicates the form of the equations:
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int kl;// The number of subdiagonals within the band of A
      lapack_int ku; // The number of superdiagonals within the band of A
      lapack_int ldab;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      double *ab, *abref; //The array 'a' contains the matrix A
      lapack_int *ipiv, *ipivref; // The pivot indices
      void *hModule, *dModule;
      int b_bufsize;

      /* Input/ Output parameters */
      double *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      gbtrs_double_parameters ( int matrix_layout_i, char trans_i,
                                lapack_int n_i, lapack_int ldab_i,
								 lapack_int kl_i, lapack_int ku_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~gbtrs_double_parameters ();
};  /* end of gbtrs_double_parameters  class definition */


/* Constructor gbtrs_double_parameters definition */
gbtrs_double_parameters:: gbtrs_double_parameters ( int matrix_layout_i,
                        char trans_i, lapack_int n_i, lapack_int ldab_i,
								       lapack_int kl_i, lapack_int ku_i,
                                 lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    trans = trans_i;
	kl = kl_i;
	ku = ku_i;
    ldab = ldab_i;
    ldb = ldb_i;
    nrhs = nrhs_i;

    hModule = NULL;
    dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n gbtrs Double:  n: %d, trans: %c lda: %d ldb: %d nrhs: %d \n",
             n, trans, ldab, ldb, nrhs);
#endif


    if(matrix_layout==LAPACK_COL_MAJOR){
        b_bufsize = ldb*nrhs;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        b_bufsize = ldb*n;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &ab, &abref, (n*ldab));
    lapacke_gtest_alloc_double_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (ab==NULL) || (abref==NULL) ||  \
        (b==NULL) || (bref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "gbtrs_double_parameters object: malloc error.";
       gbtrs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( ab, abref, n*ldab);
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

gbtrs_double_parameters:: ~gbtrs_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gbtrs_double_parameters object: destructor invoked. \n");
#endif
   gbtrs_free();
}


//  Test fixture class definition
class dgbtrs_test  : public  ::testing::Test {
public:
   gbtrs_double_parameters  *dgbtrs_obj;
   void SetUp();  
   void TearDown () { delete dgbtrs_obj; }
};


void dgbtrs_test::SetUp(){

    /* LAPACKE DGBTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_dgbtrs) ( int matrix_layout, char trans,
                               lapack_int n, lapack_int kl, lapack_int ku,
						                lapack_int nrhs, const double * ab,
                                  lapack_int ldab, const lapack_int * ipiv,
                                            double * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_dgbtrs DGBTRS;

     /* LAPACKE DGBTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_dgbtrf) ( int matrix_layout,lapack_int m,
	                          lapack_int n, lapack_int kl , lapack_int ku ,
                             double* ab,lapack_int ldab,lapack_int* ipiv );

    Fptr_NL_LAPACKE_dgbtrf DGBTRF;

    dgbtrs_obj = new  gbtrs_double_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].transr,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].ldab,
                                         lin_solver_paramslist[idx].kl,
                                         lin_solver_paramslist[idx].ku,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );

    idx = Circular_Increment_Index(idx);

    dgbtrs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgbtrs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgbtrs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgbtrs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DGBTRS = (Fptr_NL_LAPACKE_dgbtrs)dlsym(dgbtrs_obj->hModule, "LAPACKE_dgbtrs");
    ASSERT_TRUE(DGBTRS != NULL) << "failed to get the Netlib LAPACKE_dgbtrs symbol";

    DGBTRF = (Fptr_NL_LAPACKE_dgbtrf)dlsym(dgbtrs_obj->hModule,"LAPACKE_dgbtrf");
    ASSERT_TRUE(DGBTRF != NULL) << "failed to get the Netlib LAPACKE_dgbtrf symbol";

    /* Pre condition: need to call gbtrf - before calling gbtrs function */

    /* Compute the Netlib-Lapacke's reference o/p */
    dgbtrs_obj->inforef = DGBTRF( dgbtrs_obj->matrix_layout,
                                    dgbtrs_obj->n, dgbtrs_obj->n,
                                    dgbtrs_obj->kl, dgbtrs_obj->ku,
                                     dgbtrs_obj->abref,
                               dgbtrs_obj->ldab, dgbtrs_obj->ipivref);

    dgbtrs_obj->inforef = DGBTRS( dgbtrs_obj->matrix_layout,
                                  dgbtrs_obj->trans, dgbtrs_obj->n,
                                    dgbtrs_obj->kl, dgbtrs_obj->ku,
                                  dgbtrs_obj->nrhs,
                                  (const double *)dgbtrs_obj->abref,
                                  dgbtrs_obj->ldab, dgbtrs_obj->ipivref,
                                  dgbtrs_obj->bref, dgbtrs_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    dgbtrs_obj->info = LAPACKE_dgbtrf( dgbtrs_obj->matrix_layout,
                                    dgbtrs_obj->n, dgbtrs_obj->n,
                                    dgbtrs_obj->kl, dgbtrs_obj->ku,
                                     dgbtrs_obj->ab,
                               dgbtrs_obj->ldab, dgbtrs_obj->ipiv);

    dgbtrs_obj->info = LAPACKE_dgbtrs( dgbtrs_obj->matrix_layout,
                                dgbtrs_obj->trans, dgbtrs_obj->n, 
                dgbtrs_obj->kl, dgbtrs_obj->ku,	dgbtrs_obj->nrhs,
                                  (const double *)dgbtrs_obj->ab,
                               dgbtrs_obj->ldab, dgbtrs_obj->ipiv,
                                 dgbtrs_obj->b, dgbtrs_obj->ldb );


    if( dgbtrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_dgbtrs is wrong\n",
                    dgbtrs_obj->info );
    }
    if( dgbtrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dgbtrs is wrong\n",
        dgbtrs_obj->inforef );
    }
}

TEST_F(dgbtrs_test, dgbtrs1) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dgbtrs_obj->b_bufsize,
                           dgbtrs_obj->b, dgbtrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(dgbtrs_test, dgbtrs2) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dgbtrs_obj->b_bufsize,
                           dgbtrs_obj->b, dgbtrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(dgbtrs_test, dgbtrs3) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dgbtrs_obj->b_bufsize,
                           dgbtrs_obj->b, dgbtrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(dgbtrs_test, dgbtrs4) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dgbtrs_obj->b_bufsize,
                           dgbtrs_obj->b, dgbtrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin gbtrs_float_parameters  class definition */
class gbtrs_float_parameters{

   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char trans; // 'N' or 'T' or 'C'. Indicates the form of the equations:
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
	  lapack_int kl;// The number of subdiagonals within the band of A
	  lapack_int ku; // The number of superdiagonals within the band of A
      lapack_int ldab;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      float *ab, *abref; //The array 'a' contains the matrix A
      lapack_int *ipiv, *ipivref; // The pivot indices
      void *hModule, *dModule;
      int b_bufsize;

      /* Input/ Output parameters */
      float *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      gbtrs_float_parameters (int matrix_layout_i, char trans_i,
                               lapack_int n_i, lapack_int ldab_i,
								lapack_int kl_i, lapack_int ku_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
      ~gbtrs_float_parameters ();
};  /* end of gbtrs_float_parameters  class definition */


/* Constructor gbtrs_float_parameters definition */
gbtrs_float_parameters:: gbtrs_float_parameters ( int matrix_layout_i,
                      char trans_i, lapack_int n_i, lapack_int ldab_i,
								     lapack_int kl_i, lapack_int ku_i,
                                lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
	kl = kl_i;
	ku = ku_i;
    trans = trans_i;
    ldab = ldab_i;
    ldb = ldb_i;
    nrhs = nrhs_i;

    hModule = NULL;
    dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n gbtrs float:  n: %d, trans: %c lda: %d ldb: %d nrhs: %d \n",
             n, trans, ldab, ldb, nrhs);
#endif

    if(matrix_layout==LAPACK_COL_MAJOR){
        b_bufsize = ldb*nrhs;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        b_bufsize = ldb*n;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &ab, &abref, (ldab*n));
    lapacke_gtest_alloc_float_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (ab==NULL) || (abref==NULL) ||  \
        (b==NULL) || (bref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "gbtrs_double_parameters object: malloc error.";
       gbtrs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( ab, abref, ldab*n);
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, b_bufsize);

   } /* end of Constructor  */

gbtrs_float_parameters:: ~gbtrs_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gbtrs_float_parameters object: destructor invoked. \n");
#endif
   gbtrs_free();
}


//  Test fixture class definition
class sgbtrs_test  : public  ::testing::Test {
public:
   gbtrs_float_parameters  *sgbtrs_obj;
   void SetUp();  
   void TearDown () { delete sgbtrs_obj; }
};


void sgbtrs_test::SetUp(){

    /* LAPACKE SGBTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_sgbtrs) ( int matrix_layout, char trans,
                               lapack_int n, lapack_int kl, lapack_int ku, 
						                 lapack_int nrhs, const float *ab,
                                 lapack_int ldab, const lapack_int * ipiv,
                                            float *b, lapack_int ldb  );

    Fptr_NL_LAPACKE_sgbtrs SGBTRS;

     /* LAPACKE SGBTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_sgbtrf) ( int matrix_layout,lapack_int m,
	                            lapack_int n, lapack_int kl, lapack_int ku,
                              float *ab,lapack_int ldab,lapack_int* ipiv );

    Fptr_NL_LAPACKE_sgbtrf SGBTRF;

    sgbtrs_obj = new  gbtrs_float_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].transr,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].ldab,
                                         lin_solver_paramslist[idx].kl,
                                         lin_solver_paramslist[idx].ku,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    sgbtrs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgbtrs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgbtrs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgbtrs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SGBTRS = (Fptr_NL_LAPACKE_sgbtrs)dlsym(sgbtrs_obj->hModule, "LAPACKE_sgbtrs");
    ASSERT_TRUE(SGBTRS != NULL) << "failed to get the Netlib LAPACKE_sgbtrs symbol";

    SGBTRF = (Fptr_NL_LAPACKE_sgbtrf)dlsym(sgbtrs_obj->hModule,"LAPACKE_sgbtrf");
    ASSERT_TRUE(SGBTRF != NULL) << "failed to get the Netlib LAPACKE_sgbtrf symbol";

    /* Pre condition: need to call gbtrf - before calling gbtrs function */

    /* Compute the Netlib-Lapacke's reference o/p */
    sgbtrs_obj->inforef = SGBTRF( sgbtrs_obj->matrix_layout,
                                    sgbtrs_obj->n, sgbtrs_obj->n,
                                    sgbtrs_obj->kl, sgbtrs_obj->ku,
                                     sgbtrs_obj->abref,
                               sgbtrs_obj->ldab, sgbtrs_obj->ipivref);

    sgbtrs_obj->inforef = SGBTRS( sgbtrs_obj->matrix_layout,
                                  sgbtrs_obj->trans, sgbtrs_obj->n,
                                    sgbtrs_obj->kl, sgbtrs_obj->ku,
                                  sgbtrs_obj->nrhs,
                                  (const float *)sgbtrs_obj->abref,
                                  sgbtrs_obj->ldab, sgbtrs_obj->ipivref,
                                  sgbtrs_obj->bref, sgbtrs_obj->ldb);

    /* Compute the Libflme lapacke o/p by invoking Netlib-Lapack's API */
    sgbtrs_obj->info = LAPACKE_sgbtrf( sgbtrs_obj->matrix_layout,
                                    sgbtrs_obj->n, sgbtrs_obj->n,
                                    sgbtrs_obj->kl, sgbtrs_obj->ku,
                                     sgbtrs_obj->ab,
                               sgbtrs_obj->ldab, sgbtrs_obj->ipiv);

    sgbtrs_obj->info = LAPACKE_sgbtrs( sgbtrs_obj->matrix_layout,
                                sgbtrs_obj->trans, sgbtrs_obj->n, 
                sgbtrs_obj->kl, sgbtrs_obj->ku,	sgbtrs_obj->nrhs,
                                  (const float *)sgbtrs_obj->ab,
                               sgbtrs_obj->ldab, sgbtrs_obj->ipiv,
                                 sgbtrs_obj->b, sgbtrs_obj->ldb );
    if( sgbtrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_sgbtrs is wrong\n",
                    sgbtrs_obj->info );
    }
    if( sgbtrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sgbtrs is wrong\n",
        sgbtrs_obj->inforef );
    }
}

TEST_F(sgbtrs_test, sgbtrs1) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( sgbtrs_obj->b_bufsize,
                           sgbtrs_obj->b, sgbtrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin gbtrs_scomplex_parameters  class definition */
class gbtrs_scomplex_parameters{
   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char trans; // 'N' or 'T' or 'C'. Indicates the form of the equations:
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
	  lapack_int kl;// The number of subdiagonals within the band of A
	  lapack_int ku; // The number of superdiagonals within the band of A
      lapack_int ldab;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_complex_float *ab, *abref; //The array 'a' contains the matrix A
      lapack_int *ipiv, *ipivref; // The pivot indices
      void *hModule, *dModule;
      int b_bufsize;

      /* Input/ Output parameters */
      lapack_complex_float *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      gbtrs_scomplex_parameters ( int matrix_layout_i, char trans_i,
                                lapack_int n_i, lapack_int ldab_i,
								 lapack_int kl_i, lapack_int ku_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~gbtrs_scomplex_parameters ();
};  /* end of gbtrs_scomplex_parameters  class definition */


/* Constructor gbtrs_scomplex_parameters definition */
gbtrs_scomplex_parameters:: gbtrs_scomplex_parameters ( int matrix_layout_i,
                         char trans_i, lapack_int n_i, lapack_int ldab_i,
								       lapack_int kl_i, lapack_int ku_i,
                                lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
	kl = kl_i;
	ku = ku_i;
    trans = trans_i;
    ldab = ldab_i;
    ldb = ldb_i;
    nrhs = nrhs_i;

    hModule = NULL;
    dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n gbtrs scomplex:  n: %d, trans: %c lda: %d ldb: %d nrhs: %d \n",
             n, trans, ldab, ldb, nrhs);
#endif

    if(matrix_layout==LAPACK_COL_MAJOR){
        b_bufsize = ldb*nrhs;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        b_bufsize = ldb*n;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &ab, &abref, (n*ldab));
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (ab==NULL) || (abref==NULL) ||  \
        (b==NULL) || (bref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "gbtrs_scomplex_parameters object: malloc error.";
       gbtrs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( ab, abref, n*ldab);
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

gbtrs_scomplex_parameters:: ~gbtrs_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gbtrs_scomplex_parameters object: destructor invoked. \n");
#endif
   gbtrs_free();
}


//  Test fixture class definition
class cgbtrs_test  : public  ::testing::Test {
public:
   gbtrs_scomplex_parameters  *cgbtrs_obj;
   void SetUp();  
   void TearDown () { delete cgbtrs_obj; }
};


void cgbtrs_test::SetUp(){

    /* LAPACKE CGBTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_cgbtrs) ( int matrix_layout, char trans,
                               lapack_int n, lapack_int kl, lapack_int ku, 
						 lapack_int nrhs, const lapack_complex_float * ab,
                                 lapack_int ldab, const lapack_int * ipiv,
                              lapack_complex_float * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_cgbtrs CGBTRS;

     /* LAPACKE CGBTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_cgbtrf) ( int matrix_layout,lapack_int m,
	                            lapack_int n, lapack_int kl, lapack_int ku,
               lapack_complex_float* ab,lapack_int ldab,lapack_int* ipiv );

    Fptr_NL_LAPACKE_cgbtrf CGBTRF;


    cgbtrs_obj = new  gbtrs_scomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].transr,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].ldab,
                                         lin_solver_paramslist[idx].kl,
                                         lin_solver_paramslist[idx].ku,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    cgbtrs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgbtrs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgbtrs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgbtrs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CGBTRS = (Fptr_NL_LAPACKE_cgbtrs)dlsym(cgbtrs_obj->hModule, "LAPACKE_cgbtrs");
    ASSERT_TRUE(CGBTRS != NULL) << "failed to get the Netlib LAPACKE_cgbtrs symbol";

    CGBTRF = (Fptr_NL_LAPACKE_cgbtrf)dlsym(cgbtrs_obj->hModule,"LAPACKE_cgbtrf");
    ASSERT_TRUE(CGBTRF != NULL) << "failed to get the Netlib LAPACKE_cgbtrf symbol";

    /* Pre condition: need to call gbtrf - before calling gbtrs function */

    /* Compute the Netlib-Lapacke's reference o/p */
    cgbtrs_obj->inforef = CGBTRF( cgbtrs_obj->matrix_layout,
                                    cgbtrs_obj->n, cgbtrs_obj->n,
                                    cgbtrs_obj->kl, cgbtrs_obj->ku,
                                     cgbtrs_obj->abref,
                               cgbtrs_obj->ldab, cgbtrs_obj->ipivref);

    cgbtrs_obj->inforef = CGBTRS( cgbtrs_obj->matrix_layout,
                                  cgbtrs_obj->trans, cgbtrs_obj->n,
                                  cgbtrs_obj->kl, cgbtrs_obj->ku,
                                  cgbtrs_obj->nrhs,
                                  (const lapack_complex_float *)cgbtrs_obj->abref,
                                  cgbtrs_obj->ldab, cgbtrs_obj->ipivref,
                                  cgbtrs_obj->bref, cgbtrs_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    cgbtrs_obj->info = LAPACKE_cgbtrf( cgbtrs_obj->matrix_layout,
                                    cgbtrs_obj->n, cgbtrs_obj->n,
                                    cgbtrs_obj->kl, cgbtrs_obj->ku,
                                     cgbtrs_obj->ab,
                               cgbtrs_obj->ldab, cgbtrs_obj->ipiv);

    cgbtrs_obj->info = LAPACKE_cgbtrs( cgbtrs_obj->matrix_layout,
                cgbtrs_obj->trans, cgbtrs_obj->n, cgbtrs_obj->kl,
                                     cgbtrs_obj->ku,cgbtrs_obj->nrhs,
                       (const lapack_complex_float *)cgbtrs_obj->ab,
                               cgbtrs_obj->ldab, cgbtrs_obj->ipiv,
                                 cgbtrs_obj->b, cgbtrs_obj->ldb );

    if( cgbtrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_cgbtrs is wrong\n",
                    cgbtrs_obj->info );
    }
    if( cgbtrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cgbtrs is wrong\n",
        cgbtrs_obj->inforef );
    }
}

TEST_F(cgbtrs_test, cgbtrs1) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( cgbtrs_obj->b_bufsize,
                           cgbtrs_obj->b, cgbtrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(cgbtrs_test, cgbtrs2) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( cgbtrs_obj->b_bufsize,
                           cgbtrs_obj->b, cgbtrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(cgbtrs_test, cgbtrs3) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( cgbtrs_obj->b_bufsize,
                           cgbtrs_obj->b, cgbtrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(cgbtrs_test, cgbtrs4) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( cgbtrs_obj->b_bufsize,
                           cgbtrs_obj->b, cgbtrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin gbtrs_dcomplex_parameters  class definition */
class gbtrs_dcomplex_parameters{
   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char trans; // 'N' or 'T' or 'C'. Indicates the form of the equations:
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
	  lapack_int kl;// The number of subdiagonals within the band of A
	  lapack_int ku; // The number of superdiagonals within the band of A
      lapack_int ldab;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_complex_double *ab, *abref; //The array 'a' contains the matrix A
      lapack_int *ipiv, *ipivref; // The pivot indices
      void *hModule, *dModule;
      int b_bufsize;

      /* Input/ Output parameters */
      lapack_complex_double *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      gbtrs_dcomplex_parameters ( int matrix_layout_i, char trans_i,
                                lapack_int n_i, lapack_int ldab_i,
								 lapack_int kl_i, lapack_int ku_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~gbtrs_dcomplex_parameters ();
};  /* end of gbtrs_dcomplex_parameters  class definition */


/* Constructor gbtrs_dcomplex_parameters definition */
gbtrs_dcomplex_parameters:: gbtrs_dcomplex_parameters ( int matrix_layout_i,
                         char trans_i, lapack_int n_i, lapack_int ldab_i,
								       lapack_int kl_i, lapack_int ku_i,
                                lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
	kl = kl_i;
	ku = ku_i;
    trans = trans_i;
    ldab = ldab_i;
    ldb = ldb_i;
    nrhs = nrhs_i;

    hModule = NULL;
    dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n gbtrs DComplex:  n: %d, trans: %c lda: %d ldb: %d nrhs: %d \n",
             n, trans, ldab, ldb, nrhs);
#endif

    if(matrix_layout==LAPACK_COL_MAJOR){
        b_bufsize = ldb*nrhs;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        b_bufsize = ldb*n;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &ab, &abref, (n*ldab));
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (ab==NULL) || (abref==NULL) ||  \
        (b==NULL) || (bref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "gbtrs_dcomplex_parameters object: malloc error.";
       gbtrs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( ab, abref, n*ldab);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

gbtrs_dcomplex_parameters:: ~gbtrs_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gbtrs_dcomplex_parameters object: destructor invoked. \n");
#endif
   gbtrs_free();
}


//  Test fixture class definition
class zgbtrs_test  : public  ::testing::Test {
public:
   gbtrs_dcomplex_parameters  *zgbtrs_obj;
   void SetUp();  
   void TearDown () { delete zgbtrs_obj; }
};


void zgbtrs_test::SetUp(){

    /* LAPACKE ZGBTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_zgbtrs) ( int matrix_layout, char trans,
                               lapack_int n, lapack_int kl, lapack_int ku,
						lapack_int nrhs, const lapack_complex_double * ab,
                                 lapack_int ldab, const lapack_int * ipiv,
                              lapack_complex_double * b, lapack_int ldb );

    Fptr_NL_LAPACKE_zgbtrs ZGBTRS;

     /* LAPACKE ZGBTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zgbtrf) ( int matrix_layout,lapack_int m,
	                            lapack_int n, lapack_int kl, lapack_int ku,
              lapack_complex_double* ab,lapack_int ldab,lapack_int* ipiv );

    Fptr_NL_LAPACKE_zgbtrf ZGBTRF;


    zgbtrs_obj = new  gbtrs_dcomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].transr,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].ldab,
                                         lin_solver_paramslist[idx].kl,
                                         lin_solver_paramslist[idx].ku,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    zgbtrs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgbtrs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgbtrs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgbtrs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZGBTRS = (Fptr_NL_LAPACKE_zgbtrs)dlsym(zgbtrs_obj->hModule, "LAPACKE_zgbtrs");
    ASSERT_TRUE(ZGBTRS != NULL) << "failed to get the Netlib LAPACKE_zgbtrs symbol";

    ZGBTRF = (Fptr_NL_LAPACKE_zgbtrf)dlsym(zgbtrs_obj->hModule,"LAPACKE_zgbtrf");
    ASSERT_TRUE(ZGBTRF != NULL) << "failed to get the Netlib LAPACKE_zgbtrf symbol";

    /* Pre condition: need to call gbtrf - before calling gbtrs function */

    /* Compute the Netlib-Lapacke's reference o/p */
    zgbtrs_obj->inforef = ZGBTRF( zgbtrs_obj->matrix_layout,
                                    zgbtrs_obj->n, zgbtrs_obj->n,
                                    zgbtrs_obj->kl, zgbtrs_obj->ku,
                                     zgbtrs_obj->abref,
                               zgbtrs_obj->ldab, zgbtrs_obj->ipivref);

    zgbtrs_obj->inforef = ZGBTRS( zgbtrs_obj->matrix_layout,
                                  zgbtrs_obj->trans, zgbtrs_obj->n,
                                  zgbtrs_obj->kl, zgbtrs_obj->ku,
                                  zgbtrs_obj->nrhs,
                       (const lapack_complex_double *)zgbtrs_obj->abref,
                                  zgbtrs_obj->ldab, zgbtrs_obj->ipivref,
                                  zgbtrs_obj->bref, zgbtrs_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    zgbtrs_obj->info = LAPACKE_zgbtrf( zgbtrs_obj->matrix_layout,
                                    zgbtrs_obj->n, zgbtrs_obj->n,
                                  zgbtrs_obj->kl, zgbtrs_obj->ku,
                                     zgbtrs_obj->ab,
                               zgbtrs_obj->ldab, zgbtrs_obj->ipiv);

    zgbtrs_obj->info = LAPACKE_zgbtrs( zgbtrs_obj->matrix_layout,
                zgbtrs_obj->trans, zgbtrs_obj->n, zgbtrs_obj->kl,
                                     zgbtrs_obj->ku,zgbtrs_obj->nrhs,
                                  (const lapack_complex_double *)zgbtrs_obj->ab,
                               zgbtrs_obj->ldab, zgbtrs_obj->ipiv,
                                 zgbtrs_obj->b, zgbtrs_obj->ldb );


    if( zgbtrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zgbtrs is wrong\n",
                    zgbtrs_obj->info );
    }
    if( zgbtrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zgbtrs is wrong\n",
        zgbtrs_obj->inforef );
    }
}

TEST_F(zgbtrs_test, zgbtrs1) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zgbtrs_obj->b_bufsize,
                           zgbtrs_obj->b, zgbtrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zgbtrs_test, zgbtrs2) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zgbtrs_obj->b_bufsize,
                           zgbtrs_obj->b, zgbtrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zgbtrs_test, zgbtrs3) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zgbtrs_obj->b_bufsize,
                           zgbtrs_obj->b, zgbtrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zgbtrs_test, zgbtrs4) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zgbtrs_obj->b_bufsize,
                           zgbtrs_obj->b, zgbtrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
