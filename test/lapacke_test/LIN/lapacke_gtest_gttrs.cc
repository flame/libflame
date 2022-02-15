#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"

#define gttrs_free() \
       if (ipiv != NULL) free (ipiv); \
       if (ipivref != NULL)free (ipivref); \
       if (b != NULL)    free (b   ); \
       if (bref != NULL) free (bref); \
       if (d != NULL)    free (d   ); \
       if (dref != NULL) free (dref); \
       if (du != NULL)    free (du   ); \
       if (duref != NULL) free (duref); \
       if (du2 != NULL)    free (du2   ); \
       if (du2ref != NULL) free (du2ref); \
       if (dl != NULL)    free (dl   ); \
       if (dlref != NULL) free (dlref); \
       if( hModule != NULL) dlclose(hModule); \
       if(dModule != NULL) dlclose(dModule)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;


/* Begin gttrs_double_parameters  class definition */
class gttrs_double_parameters{
   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char trans; // 'N' or 'T' or 'C'. Indicates the form of the equations:
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      double * dl, *dlref; // (n - 1) multipliers that define the matrix L
	  double * d, *dref; //   n diagonal elements of the upper triangular matrix U
	  double * du, *duref; //  (n - 1) elements of the first superdiagonal of U.
	  double * du2, *du2ref; //	(n - 2) elements of the second superdiagonal of U
      lapack_int *ipiv, *ipivref; // The pivot indices
      void *hModule, *dModule;
      int b_bufsize;

      /* Input/ Output parameters */
      double *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      gttrs_double_parameters ( int matrix_layout_i, char trans_i,
              lapack_int n_i, lapack_int nrhs_i, lapack_int ldb_i);
             
      ~gttrs_double_parameters ();
};  /* end of gttrs_double_parameters  class definition */


/* Constructor gttrs_double_parameters definition */
gttrs_double_parameters:: gttrs_double_parameters ( int matrix_layout_i, 
     char trans_i, lapack_int n_i, lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    trans = trans_i;
    ldb = ldb_i;
    nrhs = nrhs_i;

    hModule = NULL;
    dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n gttrs Double:  n: %d, trans: %c ldb: %d nrhs: %d \n",
             n, trans, ldb, nrhs);
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
    lapacke_gtest_alloc_double_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_double_buffer_pair( &d, &dref, n);
    lapacke_gtest_alloc_double_buffer_pair( &dl, &dlref, (n-1));
    lapacke_gtest_alloc_double_buffer_pair( &du, &duref, (n-1));
    lapacke_gtest_alloc_double_buffer_pair( &du2, &du2ref, (n-2));
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (b==NULL) || (bref==NULL) ||  \
        (dl==NULL) || (dlref==NULL) ||  \
        (d==NULL) || (dref==NULL) ||  \
        (du==NULL) || (duref==NULL) ||  \
        (du2==NULL) || (du2ref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       gttrs_free();
       EXPECT_FALSE( true) << "gttrs_double_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, b_bufsize);
    lapacke_gtest_init_double_buffer_pair_rand( d, dref, n);
    lapacke_gtest_init_double_buffer_pair_rand( dl, dlref, (n-1));
    lapacke_gtest_init_double_buffer_pair_rand( du, duref, (n-1));
    lapacke_gtest_init_double_buffer_pair_rand( du2, du2ref, (n-2));
   

   } /* end of Constructor  */

gttrs_double_parameters:: ~gttrs_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gttrs_double_parameters object: destructor invoked. \n");
#endif
   gttrs_free();
}


//  Test fixture class definition
class dgttrs_test  : public  ::testing::Test {
public:
   gttrs_double_parameters  *dgttrs_obj;
   void SetUp();  
   void TearDown () { delete dgttrs_obj; }
};


void dgttrs_test::SetUp(){

    /* LAPACKE DGTTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_dgttrs) ( int matrix_layout, char trans, 
	lapack_int n, lapack_int nrhs, const double * dl, const double * d, 
	const double * du, const double * du2, const lapack_int * ipiv, 
	double * b, lapack_int ldb );

    Fptr_NL_LAPACKE_dgttrs DGTTRS;

     /* LAPACKE DGTTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_dgttrf) ( lapack_int n, double * dl,
	double * d, double * du, double * du2, lapack_int * ipiv);

    Fptr_NL_LAPACKE_dgttrf DGTTRF;

    dgttrs_obj = new  gttrs_double_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].transr,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);

    dgttrs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgttrs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgttrs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgttrs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DGTTRS = (Fptr_NL_LAPACKE_dgttrs)dlsym(dgttrs_obj->hModule, "LAPACKE_dgttrs");
    ASSERT_TRUE(DGTTRS != NULL) << "failed to get the Netlib LAPACKE_dgttrs symbol";

    DGTTRF = (Fptr_NL_LAPACKE_dgttrf)dlsym(dgttrs_obj->hModule,"LAPACKE_dgttrf");
    ASSERT_TRUE(DGTTRF != NULL) << "failed to get the Netlib LAPACKE_dgttrf symbol";

    /* Pre condition: need to call gttrf - before calling gttrs function */

    /* Compute the Netlib-Lapacke's reference o/p */
    dgttrs_obj->inforef = DGTTRF( dgttrs_obj->n,
                                    dgttrs_obj->dlref, dgttrs_obj->dref,
                                    dgttrs_obj->duref, dgttrs_obj->du2ref,
                           dgttrs_obj->ipivref);

    dgttrs_obj->inforef = DGTTRS( dgttrs_obj->matrix_layout,
                                  dgttrs_obj->trans, dgttrs_obj->n,
                                  dgttrs_obj->nrhs,
                                  (const double *)dgttrs_obj->dlref,
                                  (const double *)dgttrs_obj->dref,
                                  (const double *)dgttrs_obj->duref,
                                  (const double *)dgttrs_obj->du2ref,
                                  dgttrs_obj->ipivref, 
                                  dgttrs_obj->bref, dgttrs_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    dgttrs_obj->info = LAPACKE_dgttrf( dgttrs_obj->n,
                                    dgttrs_obj->dl, dgttrs_obj->d,
                                     dgttrs_obj->du,
                               dgttrs_obj->du2, dgttrs_obj->ipiv);

    dgttrs_obj->info = LAPACKE_dgttrs( dgttrs_obj->matrix_layout,
                dgttrs_obj->trans, dgttrs_obj->n, dgttrs_obj->nrhs,
                                  (const double *)dgttrs_obj->dl,
                                  (const double *)dgttrs_obj->d,
                                  (const double *)dgttrs_obj->du,
                                  (const double *)dgttrs_obj->du2,
                                   dgttrs_obj->ipiv,
                                 dgttrs_obj->b, dgttrs_obj->ldb );

    if( dgttrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
		           LAPACKE_dgttrs is wrong\n", dgttrs_obj->info );
    }
    if( dgttrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dgttrs \
		                    is wrong\n",  dgttrs_obj->inforef );
    }
}

TEST_F(dgttrs_test, dgttrs1) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dgttrs_obj->b_bufsize,
                           dgttrs_obj->b, dgttrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(dgttrs_test, dgttrs2) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dgttrs_obj->b_bufsize,
                           dgttrs_obj->b, dgttrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(dgttrs_test, dgttrs3) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dgttrs_obj->b_bufsize,
                           dgttrs_obj->b, dgttrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(dgttrs_test, dgttrs4) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dgttrs_obj->b_bufsize,
                           dgttrs_obj->b, dgttrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin gttrs_float_parameters  class definition */
class gttrs_float_parameters{
   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char trans; // 'N' or 'T' or 'C'. Indicates the form of the equations:
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      float * dl, *dlref; // (n - 1) multipliers that define the matrix L
	  float * d, *dref; //   n diagonal elements of the upper triangular matrix U
	  float * du, *duref; //  (n - 1) elements of the first superdiagonal of U.
	  float * du2, *du2ref; //	(n - 2) elements of the second superdiagonal of U
      lapack_int *ipiv, *ipivref; // The pivot indices
      void *hModule, *dModule;
      int b_bufsize;

      /* Input/ Output parameters */
      float *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      gttrs_float_parameters ( int matrix_layout_i, char trans_i,
              lapack_int n_i, lapack_int nrhs_i, lapack_int ldb_i);
             
      ~gttrs_float_parameters ();
};  /* end of gttrs_float_parameters  class definition */


/* Constructor gttrs_float_parameters definition */
gttrs_float_parameters:: gttrs_float_parameters ( int matrix_layout_i, 
     char trans_i, lapack_int n_i, lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    trans = trans_i;
    ldb = ldb_i;
    nrhs = nrhs_i;

    hModule = NULL;
    dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n gttrs Double:  n: %d, trans: %c ldb: %d nrhs: %d \n",
             n, trans, ldb, nrhs);
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
    lapacke_gtest_alloc_float_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_float_buffer_pair( &d, &dref, n);
    lapacke_gtest_alloc_float_buffer_pair( &dl, &dlref, (n-1));
    lapacke_gtest_alloc_float_buffer_pair( &du, &duref, (n-1));
    lapacke_gtest_alloc_float_buffer_pair( &du2, &du2ref, (n-2));
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (b==NULL) || (bref==NULL) ||  \
        (dl==NULL) || (dlref==NULL) ||  \
        (d==NULL) || (dref==NULL) ||  \
        (du==NULL) || (duref==NULL) ||  \
        (du2==NULL) || (du2ref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       gttrs_free();
       EXPECT_FALSE( true) << "gttrs_float_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, b_bufsize);
    lapacke_gtest_init_float_buffer_pair_rand( d, dref, n);
    lapacke_gtest_init_float_buffer_pair_rand( dl, dlref, (n-1));
    lapacke_gtest_init_float_buffer_pair_rand( du, duref, (n-1));
    lapacke_gtest_init_float_buffer_pair_rand( du2, du2ref, (n-2));
   

   } /* end of Constructor  */

gttrs_float_parameters:: ~gttrs_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gttrs_float_parameters object: destructor invoked. \n");
#endif
   gttrs_free();
}


//  Test fixture class definition
class sgttrs_test  : public  ::testing::Test {
public:
   gttrs_float_parameters  *sgttrs_obj;
   void SetUp();  
   void TearDown () { delete sgttrs_obj; }
};


void sgttrs_test::SetUp(){

    /* LAPACKE SGTTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_sgttrs) ( int matrix_layout, char trans, 
	lapack_int n, lapack_int nrhs, const float * dl, const float * d, 
	const float * du, const float * du2, const lapack_int * ipiv, 
	float * b, lapack_int ldb );

    Fptr_NL_LAPACKE_sgttrs SGTTRS;

     /* LAPACKE SGTTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_sgttrf) ( lapack_int n, float * dl,
	float * d, float * du, float * du2, lapack_int * ipiv);

    Fptr_NL_LAPACKE_sgttrf SGTTRF;

    sgttrs_obj = new  gttrs_float_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].transr,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);

    sgttrs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgttrs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgttrs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgttrs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SGTTRS = (Fptr_NL_LAPACKE_sgttrs)dlsym(sgttrs_obj->hModule, "LAPACKE_sgttrs");
    ASSERT_TRUE(SGTTRS != NULL) << "failed to get the Netlib LAPACKE_sgttrs symbol";

    SGTTRF = (Fptr_NL_LAPACKE_sgttrf)dlsym(sgttrs_obj->hModule,"LAPACKE_sgttrf");
    ASSERT_TRUE(SGTTRF != NULL) << "failed to get the Netlib LAPACKE_sgttrf symbol";

    /* Pre condition: need to call gttrf - before calling gttrs function */

    /* Compute the Netlib-Lapacke's reference o/p */
    sgttrs_obj->inforef = SGTTRF( sgttrs_obj->n,
                                    sgttrs_obj->dlref, sgttrs_obj->dref,
                                    sgttrs_obj->duref, sgttrs_obj->du2ref,
                           sgttrs_obj->ipivref);

    sgttrs_obj->inforef = SGTTRS( sgttrs_obj->matrix_layout,
                                  sgttrs_obj->trans, sgttrs_obj->n,
                                  sgttrs_obj->nrhs,
                                  (const float *)sgttrs_obj->dlref,
                                  (const float *)sgttrs_obj->dref,
                                  (const float *)sgttrs_obj->duref,
                                  (const float *)sgttrs_obj->du2ref,
                                  sgttrs_obj->ipivref, 
                                  sgttrs_obj->bref, sgttrs_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    sgttrs_obj->info = LAPACKE_sgttrf( sgttrs_obj->n,
                                    sgttrs_obj->dl, sgttrs_obj->d,
                                     sgttrs_obj->du,
                               sgttrs_obj->du2, sgttrs_obj->ipiv);

    sgttrs_obj->info = LAPACKE_sgttrs( sgttrs_obj->matrix_layout,
                sgttrs_obj->trans, sgttrs_obj->n, sgttrs_obj->nrhs,
                                  (const float *)sgttrs_obj->dl,
                                  (const float *)sgttrs_obj->d,
                                  (const float *)sgttrs_obj->du,
                                  (const float *)sgttrs_obj->du2,
                                   sgttrs_obj->ipiv,
                                 sgttrs_obj->b, sgttrs_obj->ldb );


    if( sgttrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_sgttrs \
		is wrong\n", sgttrs_obj->info );
    }
    if( sgttrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sgttrs is wrong\n",
        sgttrs_obj->inforef );
    }
}

TEST_F(sgttrs_test, sgttrs1) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( sgttrs_obj->b_bufsize,
                           sgttrs_obj->b, sgttrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(sgttrs_test, sgttrs2) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( sgttrs_obj->b_bufsize,
                           sgttrs_obj->b, sgttrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(sgttrs_test, sgttrs3) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( sgttrs_obj->b_bufsize,
                           sgttrs_obj->b, sgttrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(sgttrs_test, sgttrs4) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( sgttrs_obj->b_bufsize,
                           sgttrs_obj->b, sgttrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin gttrs_scomplex_parameters  class definition */
class gttrs_scomplex_parameters{
   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char trans; // 'N' or 'T' or 'C'. Indicates the form of the equations:
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_complex_float *dl, *dlref; // (n - 1) multipliers that define the matrix L
	  lapack_complex_float *d, *dref; //   n diagonal elements of the upper triangular matrix U
	  lapack_complex_float *du, *duref; //  (n - 1) elements of the first superdiagonal of U.
	  lapack_complex_float * du2, *du2ref; //	(n - 2) elements of the second superdiagonal of U
      lapack_int *ipiv, *ipivref; // The pivot indices
      void *hModule, *dModule;
      int b_bufsize;

      /* Input/ Output parameters */
      lapack_complex_float *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      gttrs_scomplex_parameters ( int matrix_layout_i, char trans_i,
              lapack_int n_i, lapack_int nrhs_i, lapack_int ldb_i);
             
      ~gttrs_scomplex_parameters ();
};  /* end of gttrs_scomplex_parameters  class definition */


/* Constructor gttrs_scomplex_parameters definition */
gttrs_scomplex_parameters:: gttrs_scomplex_parameters ( int matrix_layout_i, 
     char trans_i, lapack_int n_i, lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    trans = trans_i;
    ldb = ldb_i;
    nrhs = nrhs_i;

    hModule = NULL;
    dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n gttrs Double:  n: %d, trans: %c ldb: %d nrhs: %d \n",
             n, trans, ldb, nrhs);
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
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &d, &dref, n);
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &dl, &dlref, (n-1));
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &du, &duref, (n-1));
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &du2, &du2ref, (n-2));
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (b==NULL) || (bref==NULL) ||  \
        (dl==NULL) || (dlref==NULL) ||  \
        (d==NULL) || (dref==NULL) ||  \
        (du==NULL) || (duref==NULL) ||  \
        (du2==NULL) || (du2ref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       gttrs_free();
       EXPECT_FALSE( true) << "gttrs_scomplex_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
    lapacke_gtest_init_scomplex_buffer_pair_rand( d, dref, n);
    lapacke_gtest_init_scomplex_buffer_pair_rand( dl, dlref, (n-1));
    lapacke_gtest_init_scomplex_buffer_pair_rand( du, duref, (n-1));
    lapacke_gtest_init_scomplex_buffer_pair_rand( du2, du2ref, (n-2));
   

   } /* end of Constructor  */

gttrs_scomplex_parameters:: ~gttrs_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gttrs_scomplex_parameters object: destructor invoked. \n");
#endif
   gttrs_free();
}


//  Test fixture class definition
class cgttrs_test  : public  ::testing::Test {
public:
   gttrs_scomplex_parameters  *cgttrs_obj;
   void SetUp();  
   void TearDown () { delete cgttrs_obj; }
};


void cgttrs_test::SetUp(){

    /* LAPACKE CGTTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_cgttrs) ( int matrix_layout, char trans, 
	lapack_int n, lapack_int nrhs, const lapack_complex_float * dl, 
	const lapack_complex_float * d, 
	const lapack_complex_float * du, const lapack_complex_float * du2, 
	const lapack_int * ipiv, 
	lapack_complex_float * b, lapack_int ldb );

    Fptr_NL_LAPACKE_cgttrs CGTTRS;

     /* LAPACKE CGTTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_cgttrf) ( lapack_int n,
	lapack_complex_float * dl, lapack_complex_float * d,
	lapack_complex_float * du, lapack_complex_float * du2,
	lapack_int * ipiv);

    Fptr_NL_LAPACKE_cgttrf CGTTRF;

    cgttrs_obj = new  gttrs_scomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].transr,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);

    cgttrs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgttrs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgttrs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgttrs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CGTTRS = (Fptr_NL_LAPACKE_cgttrs)dlsym(cgttrs_obj->hModule, "LAPACKE_cgttrs");
    ASSERT_TRUE(CGTTRS != NULL) << "failed to get the Netlib LAPACKE_cgttrs symbol";

    CGTTRF = (Fptr_NL_LAPACKE_cgttrf)dlsym(cgttrs_obj->hModule,"LAPACKE_cgttrf");
    ASSERT_TRUE(CGTTRF != NULL) << "failed to get the Netlib LAPACKE_cgttrf symbol";

    /* Pre condition: need to call gttrf - before calling gttrs function */

    /* Compute the Netlib-Lapacke's reference o/p */
    cgttrs_obj->inforef = CGTTRF( cgttrs_obj->n,
                                    cgttrs_obj->dlref, cgttrs_obj->dref,
                                    cgttrs_obj->duref, cgttrs_obj->du2ref,
                           cgttrs_obj->ipivref);

    cgttrs_obj->inforef = CGTTRS( cgttrs_obj->matrix_layout,
                                  cgttrs_obj->trans, cgttrs_obj->n,
                                  cgttrs_obj->nrhs,
                                  (const lapack_complex_float *)cgttrs_obj->dlref,
                                  (const lapack_complex_float *)cgttrs_obj->dref,
                                  (const lapack_complex_float *)cgttrs_obj->duref,
                                  (const lapack_complex_float *)cgttrs_obj->du2ref,
                                  cgttrs_obj->ipivref, 
                                  cgttrs_obj->bref, cgttrs_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    cgttrs_obj->info = LAPACKE_cgttrf( cgttrs_obj->n,
                                    cgttrs_obj->dl, cgttrs_obj->d,
                                     cgttrs_obj->du,
                               cgttrs_obj->du2, cgttrs_obj->ipiv);

    cgttrs_obj->info = LAPACKE_cgttrs( cgttrs_obj->matrix_layout,
                cgttrs_obj->trans, cgttrs_obj->n, cgttrs_obj->nrhs,
                                  (const lapack_complex_float *)cgttrs_obj->dl,
                                  (const lapack_complex_float *)cgttrs_obj->d,
                                  (const lapack_complex_float *)cgttrs_obj->du,
                                  (const lapack_complex_float *)cgttrs_obj->du2,
                                   cgttrs_obj->ipiv,
                                 cgttrs_obj->b, cgttrs_obj->ldb );

    if( cgttrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_cgttrs is wrong\n",
                    cgttrs_obj->info );
    }
    if( cgttrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cgttrs is wrong\n",
        cgttrs_obj->inforef );
    }
}

TEST_F(cgttrs_test, cgttrs1) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( cgttrs_obj->b_bufsize,
                           cgttrs_obj->b, cgttrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(cgttrs_test, cgttrs2) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( cgttrs_obj->b_bufsize,
                           cgttrs_obj->b, cgttrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(cgttrs_test, cgttrs3) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( cgttrs_obj->b_bufsize,
                           cgttrs_obj->b, cgttrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(cgttrs_test, cgttrs4) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( cgttrs_obj->b_bufsize,
                           cgttrs_obj->b, cgttrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin gttrs_dcomplex_parameters  class definition */
class gttrs_dcomplex_parameters{
   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char trans; // 'N' or 'T' or 'C'. Indicates the form of the equations:
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_complex_double *dl, *dlref; // (n - 1) multipliers that define the matrix L
	  lapack_complex_double *d, *dref; //   n diagonal elements of the upper triangular matrix U
	  lapack_complex_double *du, *duref; //  (n - 1) elements of the first superdiagonal of U.
	  lapack_complex_double * du2, *du2ref; //	(n - 2) elements of the second superdiagonal of U
      lapack_int *ipiv, *ipivref; // The pivot indices
      void *hModule, *dModule;
      int b_bufsize;

      /* Input/ Output parameters */
      lapack_complex_double *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      gttrs_dcomplex_parameters ( int matrix_layout_i, char trans_i,
              lapack_int n_i, lapack_int nrhs_i, lapack_int ldb_i);
             
      ~gttrs_dcomplex_parameters ();
};  /* end of gttrs_dcomplex_parameters  class definition */


/* Constructor gttrs_dcomplex_parameters definition */
gttrs_dcomplex_parameters:: gttrs_dcomplex_parameters ( int matrix_layout_i, 
     char trans_i, lapack_int n_i, lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    trans = trans_i;
    ldb = ldb_i;
    nrhs = nrhs_i;

    hModule = NULL;
    dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n gttrs Double:  n: %d, trans: %c ldb: %d nrhs: %d \n",
             n, trans, ldb, nrhs);
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
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &d, &dref, n);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &dl, &dlref, (n-1));
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &du, &duref, (n-1));
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &du2, &du2ref, (n-2));
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (b==NULL) || (bref==NULL) ||  \
        (dl==NULL) || (dlref==NULL) ||  \
        (d==NULL) || (dref==NULL) ||  \
        (du==NULL) || (duref==NULL) ||  \
        (du2==NULL) || (du2ref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       gttrs_free();
       EXPECT_FALSE( true) << "gttrs_dcomplex_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( d, dref, n);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( dl, dlref, (n-1));
    lapacke_gtest_init_dcomplex_buffer_pair_rand( du, duref, (n-1));
    lapacke_gtest_init_dcomplex_buffer_pair_rand( du2, du2ref, (n-2));
   

   } /* end of Constructor  */

gttrs_dcomplex_parameters:: ~gttrs_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gttrs_dcomplex_parameters object: destructor invoked. \n");
#endif
   gttrs_free();
}


//  Test fixture class definition
class zgttrs_test  : public  ::testing::Test {
public:
   gttrs_dcomplex_parameters  *zgttrs_obj;
   void SetUp();  
   void TearDown () { delete zgttrs_obj; }
};


void zgttrs_test::SetUp(){

    /* LAPACKE ZGTTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_zgttrs) ( int matrix_layout, char trans, 
	lapack_int n, lapack_int nrhs, const lapack_complex_double * dl,
	const lapack_complex_double * d, const lapack_complex_double * du,
	const lapack_complex_double * du2, const lapack_int * ipiv, 
	lapack_complex_double * b, lapack_int ldb );

    Fptr_NL_LAPACKE_zgttrs ZGTTRS;

     /* LAPACKE ZGTTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zgttrf) ( lapack_int n,
	lapack_complex_double * dl,	lapack_complex_double * d,
	lapack_complex_double * du, lapack_complex_double * du2,
	lapack_int * ipiv);

    Fptr_NL_LAPACKE_zgttrf ZGTTRF;

    zgttrs_obj = new  gttrs_dcomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].transr,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);

    zgttrs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgttrs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgttrs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgttrs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZGTTRS = (Fptr_NL_LAPACKE_zgttrs)dlsym(zgttrs_obj->hModule, "LAPACKE_zgttrs");
    ASSERT_TRUE(ZGTTRS != NULL) << "failed to get the Netlib LAPACKE_zgttrs symbol";

    ZGTTRF = (Fptr_NL_LAPACKE_zgttrf)dlsym(zgttrs_obj->hModule,"LAPACKE_zgttrf");
    ASSERT_TRUE(ZGTTRF != NULL) << "failed to get the Netlib LAPACKE_zgttrf symbol";

    /* Pre condition: need to call gttrf - before calling gttrs function */

    /* Compute the Netlib-Lapacke's reference o/p */
    zgttrs_obj->inforef = ZGTTRF( zgttrs_obj->n,
                                    zgttrs_obj->dlref, zgttrs_obj->dref,
                                    zgttrs_obj->duref, zgttrs_obj->du2ref,
                           zgttrs_obj->ipivref);

    zgttrs_obj->inforef = ZGTTRS( zgttrs_obj->matrix_layout,
                                  zgttrs_obj->trans, zgttrs_obj->n,
                                  zgttrs_obj->nrhs,
                                  (const lapack_complex_double *)zgttrs_obj->dlref,
                                  (const lapack_complex_double *)zgttrs_obj->dref,
                                  (const lapack_complex_double *)zgttrs_obj->duref,
                                  (const lapack_complex_double *)zgttrs_obj->du2ref,
                                  zgttrs_obj->ipivref, 
                                  zgttrs_obj->bref, zgttrs_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    zgttrs_obj->info = LAPACKE_zgttrf( zgttrs_obj->n,
                                    zgttrs_obj->dl, zgttrs_obj->d,
                                     zgttrs_obj->du,
                               zgttrs_obj->du2, zgttrs_obj->ipiv);

    zgttrs_obj->info = LAPACKE_zgttrs( zgttrs_obj->matrix_layout,
                zgttrs_obj->trans, zgttrs_obj->n, zgttrs_obj->nrhs,
                                  (const lapack_complex_double *)zgttrs_obj->dl,
                                  (const lapack_complex_double *)zgttrs_obj->d,
                                  (const lapack_complex_double *)zgttrs_obj->du,
                                  (const lapack_complex_double *)zgttrs_obj->du2,
                                   zgttrs_obj->ipiv,
                                 zgttrs_obj->b, zgttrs_obj->ldb );

    if( zgttrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zgttrs \
		is wrong\n", zgttrs_obj->info );
    }
    if( zgttrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zgttrs is wrong\n",
        zgttrs_obj->inforef );
    }
}

TEST_F(zgttrs_test, zgttrs1) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zgttrs_obj->b_bufsize,
                           zgttrs_obj->b, zgttrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zgttrs_test, zgttrs2) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zgttrs_obj->b_bufsize,
                           zgttrs_obj->b, zgttrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zgttrs_test, zgttrs3) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zgttrs_obj->b_bufsize,
                           zgttrs_obj->b, zgttrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zgttrs_test, zgttrs4) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zgttrs_obj->b_bufsize,
                           zgttrs_obj->b, zgttrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
