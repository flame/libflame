#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"

#define gghrd_free() \
    if (a!=NULL)        free(a); \
    if (aref!=NULL)     free(aref); \
    if (b!=NULL)        free(b); \
    if (bref!=NULL)     free(bref); \
    if (z!=NULL)        free(z); \
    if (zref!=NULL)     free(zref); \
    if (q!=NULL)        free(q); \
    if (qref!=NULL)     free(qref)


// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin gghrd_float_common_parameters  class definition */
class gghrd_float_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_a, diff_b, diff_q, diff_z;
	float threshold;
    void *hModule, *dModule;
    
    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    
    char compq; // Must be 'N', 'I', or 'V'.
    char compz; // Must be 'N', 'I', or 'V'.
    lapack_int n; // The order of the matrices A and B (n≥ 0).
    lapack_int ilo; // ilo and ihi mark the rows and columns of A which are to be reduced.
    lapack_int ihi;
    
    /* Input / Output parameters */
    
    float* a, *aref; // contains the n-by-n general matrix A.
    lapack_int lda; //  The leading dimension of a
    float* b, *bref; // contains the n-by-n upper triangular matrix B.
    lapack_int ldb; //  The leading dimension of b
    float* q, *qref; // If compq = 'V', then q is the orthogonal/unitary matrix Q1
              // , typically from the QR factorization of B.
    
    lapack_int ldq; // The leading dimension of q;
    float* z, *zref; //
    lapack_int ldz; // The leading dimension of z;
    
    /*Return Values */
    int info, inforef;

   public:
      gghrd_float_parameters (int matrix_layout_i, char compq, char compz,
                              lapack_int n);

      ~gghrd_float_parameters ();
};

/* Constructor definition  float_common_parameters */
gghrd_float_parameters:: gghrd_float_parameters (int matrix_layout_i,
                            char compq_i, char compz_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    compq = compq_i;
    compz = compz_i;
    n = n_i;

    lda = n;
    ldb = n;
    ldq = n;
    ldz = n;
    
    hModule = NULL;
    dModule = NULL;
    diff_a = 0;
    diff_b = 0;
    diff_q = 0;
    diff_z = 0;
	// setting 'ilo', 'ihi' to standard values.
	ilo = 1;
	ihi = n;

#if LAPACKE_TEST_VERBOSE
   printf(" \n gghrd float: matrix_layout: %d n: %d  compz: %c \
       compq: %c  \n", matrix_layout, n, compz, compq);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, lda*n );
    lapacke_gtest_alloc_float_buffer_pair( &b, &bref, ldb*n );
    lapacke_gtest_alloc_float_buffer_pair( &q, &qref, ldq*n );
    lapacke_gtest_alloc_float_buffer_pair( &z, &zref, ldz*n );

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (z==NULL) || (zref==NULL) || \
        (q==NULL) || (qref==NULL) ){
       EXPECT_FALSE( true) << "gghrd_float_parameters object: malloc error. Exiting ";
       gghrd_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, lda*n );
    lapacke_gtest_init_float_buffer_pair_rand_custom_matrix( b, bref, n,ldb, 'U');
    //lapacke_gtest_init_float_buffer_pair_rand( b, bref, ldb*n );
    lapacke_gtest_init_float_buffer_pair_rand( q, qref, ldq*n );
    lapacke_gtest_init_float_buffer_pair_rand( z, zref, ldz*n );
    
   } /* end of Constructor  */
    

/* Destructor definition  'gghrd_float_common_parameters' */
gghrd_float_parameters :: ~gghrd_float_parameters ()
{
   /* De-Allocate memory for the input matrices */
   gghrd_free();
} 

//  Test fixture class definition
class sgghrd_test  : public  ::testing::Test {
public:
   gghrd_float_parameters  *sgghrd_obj;
   void SetUp();  
   void TearDown () { delete sgghrd_obj; }
};

void sgghrd_test::SetUp()
{
    /* LAPACKE STRSNA prototype */
    typedef int (*Fptr_NL_LAPACKE_sgghrd) ( int matrix_layout, char compq,
                 char compz, lapack_int n, lapack_int ilo, lapack_int ihi,
                 float* a, lapack_int lda, float* b, lapack_int ldb, 
                 float* q, lapack_int ldq, float* z, lapack_int ldz);
				 
    Fptr_NL_LAPACKE_sgghrd SGGHRD;

    sgghrd_obj = new  gghrd_float_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_paramslist[idx].compz,
                                         eig_paramslist[idx].compz,
                                         eig_paramslist[idx].n );
										 
    idx = Circular_Increment_Index(idx);
    sgghrd_obj->threshold = eig_non_sym_paramslist[idx].gghrd_threshold;
    sgghrd_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgghrd_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgghrd_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgghrd_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SGGHRD = (Fptr_NL_LAPACKE_sgghrd)dlsym(sgghrd_obj->hModule, "LAPACKE_sgghrd");
    ASSERT_TRUE(SGGHRD != NULL) << "failed to ppt the Netlib LAPACKE_sgghrd symbol";
 

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    sgghrd_obj->inforef = SGGHRD(   sgghrd_obj->matrix_layout,
                                    sgghrd_obj->compq,
                                    sgghrd_obj->compz, 
                                    sgghrd_obj->n,
                                    sgghrd_obj->ilo,
                                    sgghrd_obj->ihi, 
                                    sgghrd_obj->aref,
                                    sgghrd_obj->lda,
                                    sgghrd_obj->bref,
                                    sgghrd_obj->ldb,
                                    sgghrd_obj->qref,
                                    sgghrd_obj->ldq,
                                    sgghrd_obj->zref,
                                    sgghrd_obj->ldz
                                    );

    /* Compute libflame's Lapacke o/p  */
    sgghrd_obj->info = LAPACKE_sgghrd(  sgghrd_obj->matrix_layout,
                                    sgghrd_obj->compq,
                                    sgghrd_obj->compz, 
                                    sgghrd_obj->n,
                                    sgghrd_obj->ilo,
                                    sgghrd_obj->ihi, 
                                    sgghrd_obj->a,
                                    sgghrd_obj->lda,
                                    sgghrd_obj->b,
                                    sgghrd_obj->ldb,
                                    sgghrd_obj->q,
                                    sgghrd_obj->ldq,
                                    sgghrd_obj->z,
                                    sgghrd_obj->ldz
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    /* If howmny = 'A' or 'S', then vl need not be set. */
    sgghrd_obj->diff_a =  computeDiff_s( (sgghrd_obj->lda)*(sgghrd_obj->n), 
                sgghrd_obj->a, sgghrd_obj->aref );

    sgghrd_obj->diff_b =  computeDiff_s( (sgghrd_obj->ldb)*(sgghrd_obj->n), 
                sgghrd_obj->b, sgghrd_obj->bref );

    sgghrd_obj->diff_q =  computeDiff_s( (sgghrd_obj->ldq)*(sgghrd_obj->n), 
                sgghrd_obj->q, sgghrd_obj->qref );

    sgghrd_obj->diff_z =  computeDiff_s( (sgghrd_obj->ldz)*(sgghrd_obj->n), 
                sgghrd_obj->z, sgghrd_obj->zref );
}

TEST_F(sgghrd_test, sgghrd1) {
    EXPECT_NEAR(0.0, sgghrd_obj->diff_a, sgghrd_obj->threshold);
    EXPECT_NEAR(0.0, sgghrd_obj->diff_b, sgghrd_obj->threshold);
    EXPECT_NEAR(0.0, sgghrd_obj->diff_q, sgghrd_obj->threshold);
    EXPECT_NEAR(0.0, sgghrd_obj->diff_z, sgghrd_obj->threshold);
}

TEST_F(sgghrd_test, sgghrd2) {
    EXPECT_NEAR(0.0, sgghrd_obj->diff_a, sgghrd_obj->threshold);
    EXPECT_NEAR(0.0, sgghrd_obj->diff_b, sgghrd_obj->threshold);
    EXPECT_NEAR(0.0, sgghrd_obj->diff_q, sgghrd_obj->threshold);
    EXPECT_NEAR(0.0, sgghrd_obj->diff_z, sgghrd_obj->threshold);
}

TEST_F(sgghrd_test, sgghrd3) {
    EXPECT_NEAR(0.0, sgghrd_obj->diff_a, sgghrd_obj->threshold);
    EXPECT_NEAR(0.0, sgghrd_obj->diff_b, sgghrd_obj->threshold);
    EXPECT_NEAR(0.0, sgghrd_obj->diff_q, sgghrd_obj->threshold);
    EXPECT_NEAR(0.0, sgghrd_obj->diff_z, sgghrd_obj->threshold);
}

TEST_F(sgghrd_test, sgghrd4) {
    EXPECT_NEAR(0.0, sgghrd_obj->diff_a, sgghrd_obj->threshold);
    EXPECT_NEAR(0.0, sgghrd_obj->diff_b, sgghrd_obj->threshold);
    EXPECT_NEAR(0.0, sgghrd_obj->diff_q, sgghrd_obj->threshold);
    EXPECT_NEAR(0.0, sgghrd_obj->diff_z, sgghrd_obj->threshold);
}

/* Begin gghrd_double__common_parameters  class definition */
class gghrd_double__parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    double diff_a, diff_b, diff_q, diff_z;
	float threshold;
    void *hModule, *dModule;
    
    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    
    char compq; // Must be 'N', 'I', or 'V'.
    char compz; // Must be 'N', 'I', or 'V'.
    lapack_int n; // The order of the matrices A and B (n≥ 0).
    lapack_int ilo; // ilo and ihi mark the rows and columns of A which are to be reduced.
    lapack_int ihi;
    
    /* Input / Output parameters */
    
    double* a, *aref; // contains the n-by-n general matrix A.
    lapack_int lda; //  The leading dimension of a
    double* b, *bref; // contains the n-by-n upper triangular matrix B.
    lapack_int ldb; //  The leading dimension of b
    double* q, *qref; // If compq = 'V', then q is the orthogonal/unitary matrix Q1
              // , typically from the QR factorization of B.
    
    lapack_int ldq; // The leading dimension of q;
    double* z, *zref; //
    lapack_int ldz; // The leading dimension of z;
    
    /*Return Values */
    int info, inforef;

   public:
      gghrd_double__parameters (int matrix_layout_i, char compq, char compz,
                              lapack_int n);

      ~gghrd_double__parameters ();
};

/* Constructor definition  double_common_parameters */
gghrd_double__parameters:: gghrd_double__parameters (int matrix_layout_i,
                            char compq_i, char compz_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    compq = compq_i;
    compz = compz_i;
    n = n_i;

    lda = n;
    ldb = n;
    ldq = n;
    ldz = n;
    
    hModule = NULL;
    dModule = NULL;
    diff_a = 0;
    diff_b = 0;
    diff_q = 0;
    diff_z = 0;
	// setting 'ilo', 'ihi' to standard values.
	ilo = 1;
	ihi = n;

#if LAPACKE_TEST_VERBOSE
   printf(" \n gghrd double: matrix_layout: %d n: %d  compz: %c \
       compq: %c  \n", matrix_layout, n, compz, compq);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, lda*n );
    lapacke_gtest_alloc_double_buffer_pair( &b, &bref, ldb*n );
    lapacke_gtest_alloc_double_buffer_pair( &q, &qref, ldq*n );
    lapacke_gtest_alloc_double_buffer_pair( &z, &zref, ldz*n );

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (z==NULL) || (zref==NULL) || \
        (q==NULL) || (qref==NULL) ){
       EXPECT_FALSE( true) << "gghrd_double__parameters object: malloc error. Exiting ";
       gghrd_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, lda*n );
    lapacke_gtest_init_double_buffer_pair_rand_custom_matrix( b, bref, n,ldb, 'U');
    //lapacke_gtest_init_double_buffer_pair_rand( b, bref, ldb*n );
    lapacke_gtest_init_double_buffer_pair_rand( q, qref, ldq*n );
    lapacke_gtest_init_double_buffer_pair_rand( z, zref, ldz*n );
    
   } /* end of Constructor  */
    

/* Destructor definition  'gghrd_double__common_parameters' */
gghrd_double__parameters :: ~gghrd_double__parameters ()
{
   /* De-Allocate memory for the input matrices */
   gghrd_free();
} 

//  Test fixture class definition
class dgghrd_test  : public  ::testing::Test {
public:
   gghrd_double__parameters  *dgghrd_obj;
   void SetUp();  
   void TearDown () { delete dgghrd_obj; }
};

void dgghrd_test::SetUp()
{
    /* LAPACKE STRSNA prototype */
    typedef int (*Fptr_NL_LAPACKE_dgghrd) ( int matrix_layout, char compq,
                 char compz, lapack_int n, lapack_int ilo, lapack_int ihi,
                 double* a, lapack_int lda, double* b, lapack_int ldb, 
                 double* q, lapack_int ldq, double* z, lapack_int ldz);
				 
    Fptr_NL_LAPACKE_dgghrd DGGHRD;

    dgghrd_obj = new  gghrd_double__parameters( eig_paramslist[idx].matrix_layout,
                                         eig_paramslist[idx].compz,
                                         eig_paramslist[idx].compz,
                                         eig_paramslist[idx].n );
										 
    idx = Circular_Increment_Index(idx);
    dgghrd_obj->threshold = eig_non_sym_paramslist[idx].gghrd_threshold;
    dgghrd_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgghrd_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgghrd_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgghrd_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DGGHRD = (Fptr_NL_LAPACKE_dgghrd)dlsym(dgghrd_obj->hModule, "LAPACKE_dgghrd");
    ASSERT_TRUE(DGGHRD != NULL) << "failed to ppt the Netlib LAPACKE_dgghrd symbol";
 

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    dgghrd_obj->inforef = DGGHRD(   dgghrd_obj->matrix_layout,
                                    dgghrd_obj->compq,
                                    dgghrd_obj->compz, 
                                    dgghrd_obj->n,
                                    dgghrd_obj->ilo,
                                    dgghrd_obj->ihi, 
                                    dgghrd_obj->aref,
                                    dgghrd_obj->lda,
                                    dgghrd_obj->bref,
                                    dgghrd_obj->ldb,
                                    dgghrd_obj->qref,
                                    dgghrd_obj->ldq,
                                    dgghrd_obj->zref,
                                    dgghrd_obj->ldz
                                    );

    /* Compute libflame's Lapacke o/p  */
    dgghrd_obj->info = LAPACKE_dgghrd(  dgghrd_obj->matrix_layout,
                                    dgghrd_obj->compq,
                                    dgghrd_obj->compz, 
                                    dgghrd_obj->n,
                                    dgghrd_obj->ilo,
                                    dgghrd_obj->ihi, 
                                    dgghrd_obj->a,
                                    dgghrd_obj->lda,
                                    dgghrd_obj->b,
                                    dgghrd_obj->ldb,
                                    dgghrd_obj->q,
                                    dgghrd_obj->ldq,
                                    dgghrd_obj->z,
                                    dgghrd_obj->ldz
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    /* If howmny = 'A' or 'S', then vl need not be set. */
    dgghrd_obj->diff_a =  computeDiff_d( (dgghrd_obj->lda)*(dgghrd_obj->n), 
                dgghrd_obj->a, dgghrd_obj->aref );

    dgghrd_obj->diff_b =  computeDiff_d( (dgghrd_obj->ldb)*(dgghrd_obj->n), 
                dgghrd_obj->b, dgghrd_obj->bref );

    dgghrd_obj->diff_q =  computeDiff_d( (dgghrd_obj->ldq)*(dgghrd_obj->n), 
                dgghrd_obj->q, dgghrd_obj->qref );

    dgghrd_obj->diff_z =  computeDiff_d( (dgghrd_obj->ldz)*(dgghrd_obj->n), 
                dgghrd_obj->z, dgghrd_obj->zref );
}

TEST_F(dgghrd_test, dgghrd1) {
    EXPECT_NEAR(0.0, dgghrd_obj->diff_a, dgghrd_obj->threshold);
    EXPECT_NEAR(0.0, dgghrd_obj->diff_b, dgghrd_obj->threshold);
    EXPECT_NEAR(0.0, dgghrd_obj->diff_q, dgghrd_obj->threshold);
    EXPECT_NEAR(0.0, dgghrd_obj->diff_z, dgghrd_obj->threshold);
}

TEST_F(dgghrd_test, dgghrd2) {
    EXPECT_NEAR(0.0, dgghrd_obj->diff_a, dgghrd_obj->threshold);
    EXPECT_NEAR(0.0, dgghrd_obj->diff_b, dgghrd_obj->threshold);
    EXPECT_NEAR(0.0, dgghrd_obj->diff_q, dgghrd_obj->threshold);
    EXPECT_NEAR(0.0, dgghrd_obj->diff_z, dgghrd_obj->threshold);
}

TEST_F(dgghrd_test, dgghrd3) {
    EXPECT_NEAR(0.0, dgghrd_obj->diff_a, dgghrd_obj->threshold);
    EXPECT_NEAR(0.0, dgghrd_obj->diff_b, dgghrd_obj->threshold);
    EXPECT_NEAR(0.0, dgghrd_obj->diff_q, dgghrd_obj->threshold);
    EXPECT_NEAR(0.0, dgghrd_obj->diff_z, dgghrd_obj->threshold);
}

TEST_F(dgghrd_test, dgghrd4) {
    EXPECT_NEAR(0.0, dgghrd_obj->diff_a, dgghrd_obj->threshold);
    EXPECT_NEAR(0.0, dgghrd_obj->diff_b, dgghrd_obj->threshold);
    EXPECT_NEAR(0.0, dgghrd_obj->diff_q, dgghrd_obj->threshold);
    EXPECT_NEAR(0.0, dgghrd_obj->diff_z, dgghrd_obj->threshold);
}

/* Begin gghrd_lapack_complex_float_common_parameters  class definition */
class gghrd_scomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_a, diff_b, diff_q, diff_z;
	float threshold;
    void *hModule, *dModule;
    
    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    
    char compq; // Must be 'N', 'I', or 'V'.
    char compz; // Must be 'N', 'I', or 'V'.
    lapack_int n; // The order of the matrices A and B (n≥ 0).
    lapack_int ilo; // ilo and ihi mark the rows and columns of A which are to be reduced.
    lapack_int ihi;
    
    /* Input / Output parameters */
    
    lapack_complex_float* a, *aref; // contains the n-by-n general matrix A.
    lapack_int lda; //  The leading dimension of a
    lapack_complex_float* b, *bref; // contains the n-by-n upper triangular matrix B.
    lapack_int ldb; //  The leading dimension of b
    lapack_complex_float* q, *qref; // If compq = 'V', then q is the orthogonal/unitary matrix Q1
              // , typically from the QR factorization of B.
    
    lapack_int ldq; // The leading dimension of q;
    lapack_complex_float* z, *zref; //
    lapack_int ldz; // The leading dimension of z;
    
    /*Return Values */
    int info, inforef;

   public:
      gghrd_scomplex_parameters (int matrix_layout_i, char compq, char compz,
                              lapack_int n);

      ~gghrd_scomplex_parameters ();
};

/* Constructor definition  lapack_complex_float_common_parameters */
gghrd_scomplex_parameters:: gghrd_scomplex_parameters (int matrix_layout_i,
                            char compq_i, char compz_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    compq = compq_i;
    compz = compz_i;
    n = n_i;

    lda = n;
    ldb = n;
    ldq = n;
    ldz = n;
    
    hModule = NULL;
    dModule = NULL;
    diff_a = 0;
    diff_b = 0;
    diff_q = 0;
    diff_z = 0;
	// setting 'ilo', 'ihi' to standard values.
	ilo = 1;
	ihi = n;

#if LAPACKE_TEST_VERBOSE
   printf(" \n gghrd lapack_complex_float: matrix_layout: %d n: %d  compz: %c \
       compq: %c  \n", matrix_layout, n, compz, compq);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, lda*n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &b, &bref, ldb*n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &q, &qref, ldq*n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &z, &zref, ldz*n );

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (z==NULL) || (zref==NULL) || \
        (q==NULL) || (qref==NULL) ){
       EXPECT_FALSE( true) << "gghrd_scomplex_parameters object: malloc error. Exiting ";
       gghrd_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, lda*n );
    lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix( b, bref, n,ldb, 'U');
    //lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, ldb*n );
    lapacke_gtest_init_scomplex_buffer_pair_rand( q, qref, ldq*n );
    lapacke_gtest_init_scomplex_buffer_pair_rand( z, zref, ldz*n );
    
   } /* end of Constructor  */
    

/* Destructor definition  'gghrd_lapack_complex_float_common_parameters' */
gghrd_scomplex_parameters :: ~gghrd_scomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   gghrd_free();
} 

//  Test fixture class definition
class cgghrd_test  : public  ::testing::Test {
public:
   gghrd_scomplex_parameters  *cgghrd_obj;
   void SetUp();  
   void TearDown () { delete cgghrd_obj; }
};

void cgghrd_test::SetUp()
{
    /* LAPACKE STRSNA prototype */
    typedef int (*Fptr_NL_LAPACKE_cgghrd) ( int matrix_layout, char compq,
                 char compz, lapack_int n, lapack_int ilo, lapack_int ihi,
                 lapack_complex_float* a, lapack_int lda, lapack_complex_float* b, lapack_int ldb, 
                 lapack_complex_float* q, lapack_int ldq, lapack_complex_float* z, lapack_int ldz);
				 
    Fptr_NL_LAPACKE_cgghrd CGGHRD;

    cgghrd_obj = new  gghrd_scomplex_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_paramslist[idx].compz,
                                         eig_paramslist[idx].compz,
                                         eig_paramslist[idx].n );
										 
    idx = Circular_Increment_Index(idx);
    cgghrd_obj->threshold = eig_non_sym_paramslist[idx].gghrd_threshold;
    cgghrd_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgghrd_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgghrd_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgghrd_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CGGHRD = (Fptr_NL_LAPACKE_cgghrd)dlsym(cgghrd_obj->hModule, "LAPACKE_cgghrd");
    ASSERT_TRUE(CGGHRD != NULL) << "failed to ppt the Netlib LAPACKE_cgghrd symbol";
 

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    cgghrd_obj->inforef = CGGHRD(   cgghrd_obj->matrix_layout,
                                    cgghrd_obj->compq,
                                    cgghrd_obj->compz, 
                                    cgghrd_obj->n,
                                    cgghrd_obj->ilo,
                                    cgghrd_obj->ihi, 
                                    cgghrd_obj->aref,
                                    cgghrd_obj->lda,
                                    cgghrd_obj->bref,
                                    cgghrd_obj->ldb,
                                    cgghrd_obj->qref,
                                    cgghrd_obj->ldq,
                                    cgghrd_obj->zref,
                                    cgghrd_obj->ldz
                                    );

    /* Compute libflame's Lapacke o/p  */
    cgghrd_obj->info = LAPACKE_cgghrd(  cgghrd_obj->matrix_layout,
                                    cgghrd_obj->compq,
                                    cgghrd_obj->compz, 
                                    cgghrd_obj->n,
                                    cgghrd_obj->ilo,
                                    cgghrd_obj->ihi, 
                                    cgghrd_obj->a,
                                    cgghrd_obj->lda,
                                    cgghrd_obj->b,
                                    cgghrd_obj->ldb,
                                    cgghrd_obj->q,
                                    cgghrd_obj->ldq,
                                    cgghrd_obj->z,
                                    cgghrd_obj->ldz
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    /* If howmny = 'A' or 'S', then vl need not be set. */
    cgghrd_obj->diff_a =  computeDiff_c( (cgghrd_obj->lda)*(cgghrd_obj->n), 
                cgghrd_obj->a, cgghrd_obj->aref );

    cgghrd_obj->diff_b =  computeDiff_c( (cgghrd_obj->ldb)*(cgghrd_obj->n), 
                cgghrd_obj->b, cgghrd_obj->bref );

    cgghrd_obj->diff_q =  computeDiff_c( (cgghrd_obj->ldq)*(cgghrd_obj->n), 
                cgghrd_obj->q, cgghrd_obj->qref );

    cgghrd_obj->diff_z =  computeDiff_c( (cgghrd_obj->ldz)*(cgghrd_obj->n), 
                cgghrd_obj->z, cgghrd_obj->zref );
}

TEST_F(cgghrd_test, cgghrd1) {
    EXPECT_NEAR(0.0, cgghrd_obj->diff_a, cgghrd_obj->threshold);
    EXPECT_NEAR(0.0, cgghrd_obj->diff_b, cgghrd_obj->threshold);
    EXPECT_NEAR(0.0, cgghrd_obj->diff_q, cgghrd_obj->threshold);
    EXPECT_NEAR(0.0, cgghrd_obj->diff_z, cgghrd_obj->threshold);
}

TEST_F(cgghrd_test, cgghrd2) {
    EXPECT_NEAR(0.0, cgghrd_obj->diff_a, cgghrd_obj->threshold);
    EXPECT_NEAR(0.0, cgghrd_obj->diff_b, cgghrd_obj->threshold);
    EXPECT_NEAR(0.0, cgghrd_obj->diff_q, cgghrd_obj->threshold);
    EXPECT_NEAR(0.0, cgghrd_obj->diff_z, cgghrd_obj->threshold);
}

TEST_F(cgghrd_test, cgghrd3) {
    EXPECT_NEAR(0.0, cgghrd_obj->diff_a, cgghrd_obj->threshold);
    EXPECT_NEAR(0.0, cgghrd_obj->diff_b, cgghrd_obj->threshold);
    EXPECT_NEAR(0.0, cgghrd_obj->diff_q, cgghrd_obj->threshold);
    EXPECT_NEAR(0.0, cgghrd_obj->diff_z, cgghrd_obj->threshold);
}

TEST_F(cgghrd_test, cgghrd4) {
    EXPECT_NEAR(0.0, cgghrd_obj->diff_a, cgghrd_obj->threshold);
    EXPECT_NEAR(0.0, cgghrd_obj->diff_b, cgghrd_obj->threshold);
    EXPECT_NEAR(0.0, cgghrd_obj->diff_q, cgghrd_obj->threshold);
    EXPECT_NEAR(0.0, cgghrd_obj->diff_z, cgghrd_obj->threshold);
}

/* Begin gghrd_lapack_complex_double_common_parameters  class definition */
class gghrd_dcomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    double diff_a, diff_b, diff_q, diff_z;
	double threshold;
    void *hModule, *dModule;
    
    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    
    char compq; // Must be 'N', 'I', or 'V'.
    char compz; // Must be 'N', 'I', or 'V'.
    lapack_int n; // The order of the matrices A and B (n≥ 0).
    lapack_int ilo; // ilo and ihi mark the rows and columns of A which are to be reduced.
    lapack_int ihi;
    
    /* Input / Output parameters */
    
    lapack_complex_double* a, *aref; // contains the n-by-n general matrix A.
    lapack_int lda; //  The leading dimension of a
    lapack_complex_double* b, *bref; // contains the n-by-n upper triangular matrix B.
    lapack_int ldb; //  The leading dimension of b
    lapack_complex_double* q, *qref; // If compq = 'V', then q is the orthogonal/unitary matrix Q1
              // , typically from the QR factorization of B.
    
    lapack_int ldq; // The leading dimension of q;
    lapack_complex_double* z, *zref; //
    lapack_int ldz; // The leading dimension of z;
    
    /*Return Values */
    int info, inforef;

   public:
      gghrd_dcomplex_parameters (int matrix_layout_i, char compq, char compz,
                              lapack_int n);

      ~gghrd_dcomplex_parameters ();
};

/* Constructor definition  lapack_complex_double_common_parameters */
gghrd_dcomplex_parameters:: gghrd_dcomplex_parameters (int matrix_layout_i,
                            char compq_i, char compz_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    compq = compq_i;
    compz = compz_i;
    n = n_i;

    lda = n;
    ldb = n;
    ldq = n;
    ldz = n;
    
    hModule = NULL;
    dModule = NULL;
    diff_a = 0;
    diff_b = 0;
    diff_q = 0;
    diff_z = 0;
	// setting 'ilo', 'ihi' to standard values.
	ilo = 1;
	ihi = n;

#if LAPACKE_TEST_VERBOSE
   printf(" \n gghrd lapack_complex_double: matrix_layout: %d n: %d  compz: %c \
       compq: %c  \n", matrix_layout, n, compz, compq);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, lda*n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &b, &bref, ldb*n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &q, &qref, ldq*n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &z, &zref, ldz*n );

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (z==NULL) || (zref==NULL) || \
        (q==NULL) || (qref==NULL) ){
       EXPECT_FALSE( true) << "gghrd_dcomplex_parameters object: malloc error. Exiting ";
       gghrd_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, lda*n );
    lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix( b, bref, n,ldb, 'U');
    //lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, ldb*n );
    lapacke_gtest_init_dcomplex_buffer_pair_rand( q, qref, ldq*n );
    lapacke_gtest_init_dcomplex_buffer_pair_rand( z, zref, ldz*n );
    
   } /* end of Constructor  */
    

/* Destructor definition  'gghrd_lapack_complex_double_common_parameters' */
gghrd_dcomplex_parameters :: ~gghrd_dcomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   gghrd_free();
} 

//  Test fixture class definition
class zgghrd_test  : public  ::testing::Test {
public:
   gghrd_dcomplex_parameters  *zgghrd_obj;
   void SetUp();  
   void TearDown () { delete zgghrd_obj; }
};

void zgghrd_test::SetUp()
{
    /* LAPACKE STRSNA prototype */
    typedef int (*Fptr_NL_LAPACKE_zgghrd) ( int matrix_layout, char compq,
                 char compz, lapack_int n, lapack_int ilo, lapack_int ihi,
                 lapack_complex_double* a, lapack_int lda, lapack_complex_double* b, lapack_int ldb, 
                 lapack_complex_double* q, lapack_int ldq, lapack_complex_double* z, lapack_int ldz);
				 
    Fptr_NL_LAPACKE_zgghrd ZGGHRD;

    zgghrd_obj = new  gghrd_dcomplex_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_paramslist[idx].compz,
                                         eig_paramslist[idx].compz,
                                         eig_paramslist[idx].n );
										 
    idx = Circular_Increment_Index(idx);
    zgghrd_obj->threshold = eig_non_sym_paramslist[idx].gghrd_threshold;
    zgghrd_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgghrd_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgghrd_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgghrd_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZGGHRD = (Fptr_NL_LAPACKE_zgghrd)dlsym(zgghrd_obj->hModule, "LAPACKE_zgghrd");
    ASSERT_TRUE(ZGGHRD != NULL) << "failed to ppt the Netlib LAPACKE_zgghrd symbol";
 

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    zgghrd_obj->inforef = ZGGHRD(   zgghrd_obj->matrix_layout,
                                    zgghrd_obj->compq,
                                    zgghrd_obj->compz, 
                                    zgghrd_obj->n,
                                    zgghrd_obj->ilo,
                                    zgghrd_obj->ihi, 
                                    zgghrd_obj->aref,
                                    zgghrd_obj->lda,
                                    zgghrd_obj->bref,
                                    zgghrd_obj->ldb,
                                    zgghrd_obj->qref,
                                    zgghrd_obj->ldq,
                                    zgghrd_obj->zref,
                                    zgghrd_obj->ldz
                                    );

    /* Compute libflame's Lapacke o/p  */
    zgghrd_obj->info = LAPACKE_zgghrd(  zgghrd_obj->matrix_layout,
                                    zgghrd_obj->compq,
                                    zgghrd_obj->compz, 
                                    zgghrd_obj->n,
                                    zgghrd_obj->ilo,
                                    zgghrd_obj->ihi, 
                                    zgghrd_obj->a,
                                    zgghrd_obj->lda,
                                    zgghrd_obj->b,
                                    zgghrd_obj->ldb,
                                    zgghrd_obj->q,
                                    zgghrd_obj->ldq,
                                    zgghrd_obj->z,
                                    zgghrd_obj->ldz
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    /* If howmny = 'A' or 'S', then vl need not be set. */
    zgghrd_obj->diff_a =  computeDiff_z( (zgghrd_obj->lda)*(zgghrd_obj->n), 
                zgghrd_obj->a, zgghrd_obj->aref );

    zgghrd_obj->diff_b =  computeDiff_z( (zgghrd_obj->ldb)*(zgghrd_obj->n), 
                zgghrd_obj->b, zgghrd_obj->bref );

    zgghrd_obj->diff_q =  computeDiff_z( (zgghrd_obj->ldq)*(zgghrd_obj->n), 
                zgghrd_obj->q, zgghrd_obj->qref );

    zgghrd_obj->diff_z =  computeDiff_z( (zgghrd_obj->ldz)*(zgghrd_obj->n), 
                zgghrd_obj->z, zgghrd_obj->zref );
}

TEST_F(zgghrd_test, zgghrd1) {
    EXPECT_NEAR(0.0, zgghrd_obj->diff_a, zgghrd_obj->threshold);
    EXPECT_NEAR(0.0, zgghrd_obj->diff_b, zgghrd_obj->threshold);
    EXPECT_NEAR(0.0, zgghrd_obj->diff_q, zgghrd_obj->threshold);
    EXPECT_NEAR(0.0, zgghrd_obj->diff_z, zgghrd_obj->threshold);
}

TEST_F(zgghrd_test, zgghrd2) {
    EXPECT_NEAR(0.0, zgghrd_obj->diff_a, zgghrd_obj->threshold);
    EXPECT_NEAR(0.0, zgghrd_obj->diff_b, zgghrd_obj->threshold);
    EXPECT_NEAR(0.0, zgghrd_obj->diff_q, zgghrd_obj->threshold);
    EXPECT_NEAR(0.0, zgghrd_obj->diff_z, zgghrd_obj->threshold);
}

TEST_F(zgghrd_test, zgghrd3) {
    EXPECT_NEAR(0.0, zgghrd_obj->diff_a, zgghrd_obj->threshold);
    EXPECT_NEAR(0.0, zgghrd_obj->diff_b, zgghrd_obj->threshold);
    EXPECT_NEAR(0.0, zgghrd_obj->diff_q, zgghrd_obj->threshold);
    EXPECT_NEAR(0.0, zgghrd_obj->diff_z, zgghrd_obj->threshold);
}

TEST_F(zgghrd_test, zgghrd4) {
    EXPECT_NEAR(0.0, zgghrd_obj->diff_a, zgghrd_obj->threshold);
    EXPECT_NEAR(0.0, zgghrd_obj->diff_b, zgghrd_obj->threshold);
    EXPECT_NEAR(0.0, zgghrd_obj->diff_q, zgghrd_obj->threshold);
    EXPECT_NEAR(0.0, zgghrd_obj->diff_z, zgghrd_obj->threshold);
}