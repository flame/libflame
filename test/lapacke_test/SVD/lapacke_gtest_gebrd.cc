#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"
#define LAPACKE_TEST_VERBOSE  (1)
#define gebrd_free() \
    if (a!=NULL)        free(a); \
    if (aref!=NULL)     free(aref); \
    if (d!=NULL)        free(d); \
    if (dref!=NULL)     free(dref); \
    if (e!=NULL)        free(e); \
    if (eref!=NULL)     free(eref); \
    if (tauq!=NULL)        free(tauq); \
    if (tauqref!=NULL)     free(tauqref) ; \
    if (taup!=NULL)        free(taup); \
    if (taupref!=NULL)     free(taupref)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin gebrd_float_common_parameters  class definition */
class gebrd_float_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_a, diff_d, diff_e, diff_tauq, diff_taup;
    float threshold;
    void *hModule, *dModule;
	int min_mn;
    
    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char jobz; // Must be 'A', 'S', 'O', or 'N'. 

    lapack_int m; // The number of rows of the matrix A
    lapack_int n; // The number of columns of the matrix A 
    lapack_int lda; //  The leading dimension of a

    /* Input / Output parameters */
    float* a, *aref; // contains m-by-n matrix A.

    /* Output parameters */
    // Below buffers contain the o/p orthogonal/unitary matrices
    float* d, *dref;
    float* e, *eref;
    float* tauq, *tauqref;
    float* taup, *taupref;
   
    /*Return Values */
    int info, inforef;

   public:
      gebrd_float_parameters (int matrix_layout_i, lapack_int m, lapack_int n);
      ~gebrd_float_parameters ();
};

/* Constructor definition  gebrd float_common_parameters */
gebrd_float_parameters:: gebrd_float_parameters (int matrix_layout_i,
                    lapack_int m_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    n = n_i;
    m = m_i;
	min_mn = (m<n)?m:n;
    

    if (matrix_layout == LAPACK_ROW_MAJOR)
    {
        lda = n;
    }   else
    {
        lda = m;
    }

    // Initialize with default values
    hModule = NULL;
    dModule = NULL;
    diff_a = 0;
    diff_d = 0;
    diff_e = 0;
    diff_tauq = 0;
    diff_taup = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n gebrd float: matrix_layout: %d   \t \
m: %d \t n: %d  \t lda: %d  \n",
matrix_layout, m_i, n_i, lda);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, m*n );
    lapacke_gtest_alloc_float_buffer_pair( &d, &dref, min_mn );
    lapacke_gtest_alloc_float_buffer_pair( &e, &eref, min_mn );
    lapacke_gtest_alloc_float_buffer_pair( &tauq, &tauqref, min_mn );
    lapacke_gtest_alloc_float_buffer_pair( &taup, &taupref, min_mn );

    if( (a==NULL) || (aref==NULL) ||  \
        (d==NULL) || (dref==NULL) || \
        (e==NULL) || (eref==NULL) || \
        (tauq==NULL) || (tauqref==NULL) || \
        (taup==NULL) || (taupref==NULL) ){
       EXPECT_FALSE( true) << "gebrd_float_parameters object: malloc error. Exiting ";
       gebrd_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, m*n );
    lapacke_gtest_init_float_buffer_pair_with_constant( d, dref, min_mn, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant( e, eref, min_mn, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant( tauq, tauqref, min_mn, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant( taup, taupref, min_mn, 0.0);

   } /* end of Constructor  */

/* Destructor definition  'gebrd_float_common_parameters' */
gebrd_float_parameters :: ~gebrd_float_parameters ()
{
   /* De-Allocate memory for the input matrices */
   gebrd_free();
} 

//  Test fixture class definition
class sgebrd_test  : public  ::testing::Test {
public:
   gebrd_float_parameters  *sgebrd_obj;
   void SetUp();  
   void TearDown () { delete sgebrd_obj; }
};

void sgebrd_test::SetUp()
{
    /* LAPACKE sgebrd prototype */
    typedef int (*Fptr_NL_LAPACKE_sgebrd) ( int matrix_layout,
		lapack_int m, lapack_int n, float* a, lapack_int lda,
		float* d, float* e, float* tauq, float* taup);
            
    Fptr_NL_LAPACKE_sgebrd SGEBRD;

    sgebrd_obj = new  gebrd_float_parameters( svd_paramslist[idx].matrix_layout,
                                         svd_paramslist[idx].m,
                                         svd_paramslist[idx].n );

    sgebrd_obj->threshold = svd_paramslist[idx].svd_threshold;
    idx = Circular_Increment_Index(idx);
    sgebrd_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgebrd_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgebrd_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgebrd_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SGEBRD = (Fptr_NL_LAPACKE_sgebrd)dlsym(sgebrd_obj->hModule, "LAPACKE_sgebrd");
    ASSERT_TRUE(SGEBRD != NULL) << "failed to ppt the Netlib LAPACKE_sgebrd symbol";

    /* Compute libflame's Lapacke o/p  */
    sgebrd_obj->info = LAPACKE_sgebrd(  sgebrd_obj->matrix_layout,
                                    sgebrd_obj->m,
                                    sgebrd_obj->n,
                                    sgebrd_obj->a,
                                    sgebrd_obj->lda,
                                    sgebrd_obj->d,
                                    sgebrd_obj->e,
                                    sgebrd_obj->tauq,
                                    sgebrd_obj->taup
                                    );

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    sgebrd_obj->inforef = SGEBRD( sgebrd_obj->matrix_layout,
                                    sgebrd_obj->m,
                                    sgebrd_obj->n,
                                    sgebrd_obj->aref,
                                    sgebrd_obj->lda,
                                    sgebrd_obj->dref,
                                    sgebrd_obj->eref,
                                    sgebrd_obj->tauqref,
                                    sgebrd_obj->taupref
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    sgebrd_obj->diff_d =  computeDiff_s( sgebrd_obj->min_mn, 
                sgebrd_obj->d, sgebrd_obj->dref );

    sgebrd_obj->diff_a =  computeDiff_s( (sgebrd_obj->m)*(sgebrd_obj->n), 
                sgebrd_obj->a, sgebrd_obj->aref );

    sgebrd_obj->diff_e =  computeDiff_s( sgebrd_obj->min_mn-1, 
                sgebrd_obj->e, sgebrd_obj->eref );

    sgebrd_obj->diff_tauq =  computeDiff_s( sgebrd_obj->min_mn, 
                sgebrd_obj->tauq, sgebrd_obj->tauqref );

    sgebrd_obj->diff_taup =  computeDiff_s( sgebrd_obj->min_mn, 
                sgebrd_obj->taup, sgebrd_obj->taupref );


#if LAPACKE_TEST_VERBOSE
   printf(" \n gebrd float: \n diff_a: %f \n diff_d: %f \n \
diff_e: %f \n diff_tauq: %f \n diff_taup: %f \n info: %d \t inforef: %d\n \
gebrd_threshold: %f  \t  \n",
       sgebrd_obj->diff_a, sgebrd_obj->diff_d, sgebrd_obj->diff_e,
       sgebrd_obj->diff_tauq, sgebrd_obj->diff_taup, sgebrd_obj->info, 
	   sgebrd_obj->inforef, sgebrd_obj->threshold);
#endif
}

TEST_F(sgebrd_test, sgebrd_1) {
    EXPECT_NEAR(0.0, sgebrd_obj->diff_e, sgebrd_obj->threshold);
    EXPECT_NEAR(0.0, sgebrd_obj->diff_d, sgebrd_obj->threshold);
    EXPECT_NEAR(0.0, sgebrd_obj->diff_taup, sgebrd_obj->threshold);
    EXPECT_NEAR(0.0, sgebrd_obj->diff_tauq, sgebrd_obj->threshold);
}

TEST_F(sgebrd_test, sgebrd_2) {
    EXPECT_NEAR(0.0, sgebrd_obj->diff_e, sgebrd_obj->threshold);
    EXPECT_NEAR(0.0, sgebrd_obj->diff_d, sgebrd_obj->threshold);
    EXPECT_NEAR(0.0, sgebrd_obj->diff_taup, sgebrd_obj->threshold);
    EXPECT_NEAR(0.0, sgebrd_obj->diff_tauq, sgebrd_obj->threshold);
}
TEST_F(sgebrd_test, sgebrd_3) {
    EXPECT_NEAR(0.0, sgebrd_obj->diff_e, sgebrd_obj->threshold);
    EXPECT_NEAR(0.0, sgebrd_obj->diff_d, sgebrd_obj->threshold);
    EXPECT_NEAR(0.0, sgebrd_obj->diff_taup, sgebrd_obj->threshold);
    EXPECT_NEAR(0.0, sgebrd_obj->diff_tauq, sgebrd_obj->threshold);
}
TEST_F(sgebrd_test, sgebrd_4) {
    EXPECT_NEAR(0.0, sgebrd_obj->diff_e, sgebrd_obj->threshold);
    EXPECT_NEAR(0.0, sgebrd_obj->diff_d, sgebrd_obj->threshold);
    EXPECT_NEAR(0.0, sgebrd_obj->diff_taup, sgebrd_obj->threshold);
    EXPECT_NEAR(0.0, sgebrd_obj->diff_tauq, sgebrd_obj->threshold);
}

/* Begin gebrd_double_common_parameters  class definition */
class gebrd_double_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    double diff_a, diff_d, diff_e, diff_tauq, diff_taup;
    float threshold;
    void *hModule, *dModule;
	int min_mn;
    
    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char jobz; // Must be 'A', 'S', 'O', or 'N'. 

    lapack_int m; // The number of rows of the matrix A
    lapack_int n; // The number of columns of the matrix A 
    lapack_int lda; //  The leading dimension of a

    /* Input / Output parameters */
    double* a, *aref; // contains m-by-n matrix A.

    /* Output parameters */
    // Below buffers contain the o/p orthogonal/unitary matrices
    double* d, *dref;
    double* e, *eref;
    double* tauq, *tauqref;
    double* taup, *taupref;
   
    /*Return Values */
    int info, inforef;

   public:
      gebrd_double_parameters (int matrix_layout_i, lapack_int m, lapack_int n);
      ~gebrd_double_parameters ();
};

/* Constructor definition  gebrd double_common_parameters */
gebrd_double_parameters:: gebrd_double_parameters (int matrix_layout_i,
                    lapack_int m_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    n = n_i;
    m = m_i;
	min_mn = (m<n)?m:n;
    

    if (matrix_layout == LAPACK_ROW_MAJOR)
    {
        lda = n;
    }   else
    {
        lda = m;
    }

    // Initialize with default values
    hModule = NULL;
    dModule = NULL;
    diff_a = 0;
    diff_d = 0;
    diff_e = 0;
    diff_tauq = 0;
    diff_taup = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n gebrd double: matrix_layout: %d   \t \
m: %d \t n: %d  \t lda: %d  \n",
matrix_layout, m_i, n_i, lda);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, m*n );
    lapacke_gtest_alloc_double_buffer_pair( &d, &dref, min_mn );
    lapacke_gtest_alloc_double_buffer_pair( &e, &eref, min_mn );
    lapacke_gtest_alloc_double_buffer_pair( &tauq, &tauqref, min_mn );
    lapacke_gtest_alloc_double_buffer_pair( &taup, &taupref, min_mn );

    if( (a==NULL) || (aref==NULL) ||  \
        (d==NULL) || (dref==NULL) || \
        (e==NULL) || (eref==NULL) || \
        (tauq==NULL) || (tauqref==NULL) || \
        (taup==NULL) || (taupref==NULL) ){
       EXPECT_FALSE( true) << "gebrd_double_parameters object: malloc error. Exiting ";
       gebrd_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, m*n );
    lapacke_gtest_init_double_buffer_pair_with_constant( d, dref, min_mn, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant( e, eref, min_mn, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant( tauq, tauqref, min_mn, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant( taup, taupref, min_mn, 0.0);

   } /* end of Constructor  */

/* Destructor definition  'gebrd_double_common_parameters' */
gebrd_double_parameters :: ~gebrd_double_parameters ()
{
   /* De-Allocate memory for the input matrices */
   gebrd_free();
} 

//  Test fixture class definition
class dgebrd_test  : public  ::testing::Test {
public:
   gebrd_double_parameters  *dgebrd_obj;
   void SetUp();  
   void TearDown () { delete dgebrd_obj; }
};

void dgebrd_test::SetUp()
{
    /* LAPACKE dgebrd prototype */
    typedef int (*Fptr_NL_LAPACKE_dgebrd) ( int matrix_layout,
		lapack_int m, lapack_int n, double* a, lapack_int lda,
		double* d, double* e, double* tauq, double* taup);
            
    Fptr_NL_LAPACKE_dgebrd DGEBRD;

    dgebrd_obj = new  gebrd_double_parameters( svd_paramslist[idx].matrix_layout,
                                         svd_paramslist[idx].m,
                                         svd_paramslist[idx].n );

    dgebrd_obj->threshold = svd_paramslist[idx].svd_threshold;
    idx = Circular_Increment_Index(idx);
    dgebrd_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgebrd_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgebrd_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgebrd_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DGEBRD = (Fptr_NL_LAPACKE_dgebrd)dlsym(dgebrd_obj->hModule, "LAPACKE_dgebrd");
    ASSERT_TRUE(DGEBRD != NULL) << "failed to ppt the Netlib LAPACKE_dgebrd symbol";

    /* Compute libflame's Lapacke o/p  */
    dgebrd_obj->info = LAPACKE_dgebrd(  dgebrd_obj->matrix_layout,
                                    dgebrd_obj->m,
                                    dgebrd_obj->n,
                                    dgebrd_obj->a,
                                    dgebrd_obj->lda,
                                    dgebrd_obj->d,
                                    dgebrd_obj->e,
                                    dgebrd_obj->tauq,
                                    dgebrd_obj->taup
                                    );

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    dgebrd_obj->inforef = DGEBRD( dgebrd_obj->matrix_layout,
                                    dgebrd_obj->m,
                                    dgebrd_obj->n,
                                    dgebrd_obj->aref,
                                    dgebrd_obj->lda,
                                    dgebrd_obj->dref,
                                    dgebrd_obj->eref,
                                    dgebrd_obj->tauqref,
                                    dgebrd_obj->taupref
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    dgebrd_obj->diff_d =  computeDiff_d( dgebrd_obj->min_mn, 
                dgebrd_obj->d, dgebrd_obj->dref );

    dgebrd_obj->diff_a =  computeDiff_d( (dgebrd_obj->m)*(dgebrd_obj->n), 
                dgebrd_obj->a, dgebrd_obj->aref );

    dgebrd_obj->diff_e =  computeDiff_d( dgebrd_obj->min_mn-1, 
                dgebrd_obj->e, dgebrd_obj->eref );

    dgebrd_obj->diff_tauq =  computeDiff_d( dgebrd_obj->min_mn, 
                dgebrd_obj->tauq, dgebrd_obj->tauqref );

    dgebrd_obj->diff_taup =  computeDiff_d( dgebrd_obj->min_mn, 
                dgebrd_obj->taup, dgebrd_obj->taupref );


#if LAPACKE_TEST_VERBOSE
   printf(" \n gebrd double: \n diff_a: %f \n diff_d: %f \n \
diff_e: %f \n diff_tauq: %f \n diff_taup: %f \n info: %d \t inforef: %d\n \
gebrd_threshold: %f  \t  \n",
       dgebrd_obj->diff_a, dgebrd_obj->diff_d, dgebrd_obj->diff_e,
       dgebrd_obj->diff_tauq, dgebrd_obj->diff_taup, dgebrd_obj->info, 
	   dgebrd_obj->inforef, dgebrd_obj->threshold);
#endif
}

TEST_F(dgebrd_test, dgebrd_1) {
    EXPECT_NEAR(0.0, dgebrd_obj->diff_e, dgebrd_obj->threshold);
    EXPECT_NEAR(0.0, dgebrd_obj->diff_d, dgebrd_obj->threshold);
    EXPECT_NEAR(0.0, dgebrd_obj->diff_taup, dgebrd_obj->threshold);
    EXPECT_NEAR(0.0, dgebrd_obj->diff_tauq, dgebrd_obj->threshold);
}

TEST_F(dgebrd_test, dgebrd_2) {
    EXPECT_NEAR(0.0, dgebrd_obj->diff_e, dgebrd_obj->threshold);
    EXPECT_NEAR(0.0, dgebrd_obj->diff_d, dgebrd_obj->threshold);
    EXPECT_NEAR(0.0, dgebrd_obj->diff_taup, dgebrd_obj->threshold);
    EXPECT_NEAR(0.0, dgebrd_obj->diff_tauq, dgebrd_obj->threshold);
}
TEST_F(dgebrd_test, dgebrd_3) {
    EXPECT_NEAR(0.0, dgebrd_obj->diff_e, dgebrd_obj->threshold);
    EXPECT_NEAR(0.0, dgebrd_obj->diff_d, dgebrd_obj->threshold);
    EXPECT_NEAR(0.0, dgebrd_obj->diff_taup, dgebrd_obj->threshold);
    EXPECT_NEAR(0.0, dgebrd_obj->diff_tauq, dgebrd_obj->threshold);
}
TEST_F(dgebrd_test, dgebrd_4) {
    EXPECT_NEAR(0.0, dgebrd_obj->diff_e, dgebrd_obj->threshold);
    EXPECT_NEAR(0.0, dgebrd_obj->diff_d, dgebrd_obj->threshold);
    EXPECT_NEAR(0.0, dgebrd_obj->diff_taup, dgebrd_obj->threshold);
    EXPECT_NEAR(0.0, dgebrd_obj->diff_tauq, dgebrd_obj->threshold);
}

/* Begin gebrd_scomplex_common_parameters  class definition */
class gebrd_scomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_a, diff_d, diff_e, diff_tauq, diff_taup;
    float threshold;
    void *hModule, *dModule;
	int min_mn;
    
    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char jobz; // Must be 'A', 'S', 'O', or 'N'. 

    lapack_int m; // The number of rows of the matrix A
    lapack_int n; // The number of columns of the matrix A 
    lapack_int lda; //  The leading dimension of a

    /* Input / Output parameters */
    lapack_complex_float* a, *aref; // contains m-by-n matrix A.

    /* Output parameters */
    // Below buffers contain the o/p orthogonal/unitary matrices
    float* d, *dref;
    float* e, *eref;
    lapack_complex_float* tauq, *tauqref;
    lapack_complex_float* taup, *taupref;
   
    /*Return Values */
    int info, inforef;

   public:
      gebrd_scomplex_parameters (int matrix_layout_i, lapack_int m, lapack_int n);
      ~gebrd_scomplex_parameters ();
};

/* Constructor definition  gebrd lapack_complex_float_common_parameters */
gebrd_scomplex_parameters:: gebrd_scomplex_parameters (int matrix_layout_i,
                    lapack_int m_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    n = n_i;
    m = m_i;
	min_mn = (m<n)?m:n;
    

    if (matrix_layout == LAPACK_ROW_MAJOR)
    {
        lda = n;
    }   else
    {
        lda = m;
    }

    // Initialize with default values
    hModule = NULL;
    dModule = NULL;
    diff_a = 0;
    diff_d = 0;
    diff_e = 0;
    diff_tauq = 0;
    diff_taup = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n gebrd lapack_complex_float: matrix_layout: %d   \t \
m: %d \t n: %d  \t lda: %d  \n",
matrix_layout, m_i, n_i, lda);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, m*n );
    lapacke_gtest_alloc_float_buffer_pair( &d, &dref, min_mn );
    lapacke_gtest_alloc_float_buffer_pair( &e, &eref, min_mn );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &tauq, &tauqref, min_mn );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &taup, &taupref, min_mn );

    if( (a==NULL) || (aref==NULL) ||  \
        (d==NULL) || (dref==NULL) || \
        (e==NULL) || (eref==NULL) || \
        (tauq==NULL) || (tauqref==NULL) || \
        (taup==NULL) || (taupref==NULL) ){
       EXPECT_FALSE( true) << "gebrd_scomplex_parameters object: malloc error. Exiting ";
       gebrd_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, m*n );
    lapacke_gtest_init_float_buffer_pair_with_constant( d, dref, min_mn, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant( e, eref, min_mn, 0.0);
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( tauq, tauqref, min_mn, 0.0);
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( taup, taupref, min_mn, 0.0);

   } /* end of Constructor  */

/* Destructor definition  'gebrd_scomplex_common_parameters' */
gebrd_scomplex_parameters :: ~gebrd_scomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   gebrd_free();
} 

//  Test fixture class definition
class cgebrd_test  : public  ::testing::Test {
public:
   gebrd_scomplex_parameters  *cgebrd_obj;
   void SetUp();  
   void TearDown () { delete cgebrd_obj; }
};

void cgebrd_test::SetUp()
{
    /* LAPACKE cgebrd prototype */
    typedef int (*Fptr_NL_LAPACKE_cgebrd) ( int matrix_layout,
		lapack_int m, lapack_int n, lapack_complex_float* a, lapack_int lda,
		float* d, float* e, lapack_complex_float* tauq, lapack_complex_float* taup);
            
    Fptr_NL_LAPACKE_cgebrd CGEBRD;

    cgebrd_obj = new  gebrd_scomplex_parameters( svd_paramslist[idx].matrix_layout,
                                         svd_paramslist[idx].m,
                                         svd_paramslist[idx].n );

    cgebrd_obj->threshold = svd_paramslist[idx].svd_threshold;
    idx = Circular_Increment_Index(idx);
    cgebrd_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgebrd_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgebrd_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgebrd_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CGEBRD = (Fptr_NL_LAPACKE_cgebrd)dlsym(cgebrd_obj->hModule, "LAPACKE_cgebrd");
    ASSERT_TRUE(CGEBRD != NULL) << "failed to ppt the Netlib LAPACKE_cgebrd symbol";

    /* Compute libflame's Lapacke o/p  */
    cgebrd_obj->info = LAPACKE_cgebrd(  cgebrd_obj->matrix_layout,
                                    cgebrd_obj->m,
                                    cgebrd_obj->n,
                                    cgebrd_obj->a,
                                    cgebrd_obj->lda,
                                    cgebrd_obj->d,
                                    cgebrd_obj->e,
                                    cgebrd_obj->tauq,
                                    cgebrd_obj->taup
                                    );

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    cgebrd_obj->inforef = CGEBRD( cgebrd_obj->matrix_layout,
                                    cgebrd_obj->m,
                                    cgebrd_obj->n,
                                    cgebrd_obj->aref,
                                    cgebrd_obj->lda,
                                    cgebrd_obj->dref,
                                    cgebrd_obj->eref,
                                    cgebrd_obj->tauqref,
                                    cgebrd_obj->taupref
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    cgebrd_obj->diff_d =  computeDiff_s( cgebrd_obj->min_mn, 
                cgebrd_obj->d, cgebrd_obj->dref );

    cgebrd_obj->diff_a =  computeDiff_c( (cgebrd_obj->m)*(cgebrd_obj->n), 
                cgebrd_obj->a, cgebrd_obj->aref );

    cgebrd_obj->diff_e =  computeDiff_s( cgebrd_obj->min_mn-1, 
                cgebrd_obj->e, cgebrd_obj->eref );

    cgebrd_obj->diff_tauq =  computeDiff_c( cgebrd_obj->min_mn, 
                cgebrd_obj->tauq, cgebrd_obj->tauqref );

    cgebrd_obj->diff_taup =  computeDiff_c( cgebrd_obj->min_mn, 
                cgebrd_obj->taup, cgebrd_obj->taupref );


#if LAPACKE_TEST_VERBOSE
   printf(" \n gebrd lapack_complex_float: \n diff_a: %f \n diff_d: %f \n \
diff_e: %f \n diff_tauq: %f \n diff_taup: %f \n info: %d \t inforef: %d\n \
gebrd_threshold: %f  \t  \n",
       cgebrd_obj->diff_a, cgebrd_obj->diff_d, cgebrd_obj->diff_e,
       cgebrd_obj->diff_tauq, cgebrd_obj->diff_taup, cgebrd_obj->info, 
	   cgebrd_obj->inforef, cgebrd_obj->threshold);
#endif
}

TEST_F(cgebrd_test, cgebrd_1) {
    EXPECT_NEAR(0.0, cgebrd_obj->diff_e, cgebrd_obj->threshold);
    EXPECT_NEAR(0.0, cgebrd_obj->diff_d, cgebrd_obj->threshold);
    EXPECT_NEAR(0.0, cgebrd_obj->diff_taup, cgebrd_obj->threshold);
    EXPECT_NEAR(0.0, cgebrd_obj->diff_tauq, cgebrd_obj->threshold);
}

TEST_F(cgebrd_test, cgebrd_2) {
    EXPECT_NEAR(0.0, cgebrd_obj->diff_e, cgebrd_obj->threshold);
    EXPECT_NEAR(0.0, cgebrd_obj->diff_d, cgebrd_obj->threshold);
    EXPECT_NEAR(0.0, cgebrd_obj->diff_taup, cgebrd_obj->threshold);
    EXPECT_NEAR(0.0, cgebrd_obj->diff_tauq, cgebrd_obj->threshold);
}
TEST_F(cgebrd_test, cgebrd_3) {
    EXPECT_NEAR(0.0, cgebrd_obj->diff_e, cgebrd_obj->threshold);
    EXPECT_NEAR(0.0, cgebrd_obj->diff_d, cgebrd_obj->threshold);
    EXPECT_NEAR(0.0, cgebrd_obj->diff_taup, cgebrd_obj->threshold);
    EXPECT_NEAR(0.0, cgebrd_obj->diff_tauq, cgebrd_obj->threshold);
}
TEST_F(cgebrd_test, cgebrd_4) {
    EXPECT_NEAR(0.0, cgebrd_obj->diff_e, cgebrd_obj->threshold);
    EXPECT_NEAR(0.0, cgebrd_obj->diff_d, cgebrd_obj->threshold);
    EXPECT_NEAR(0.0, cgebrd_obj->diff_taup, cgebrd_obj->threshold);
    EXPECT_NEAR(0.0, cgebrd_obj->diff_tauq, cgebrd_obj->threshold);
}

/* Begin gebrd_dcomplex_common_parameters  class definition */
class gebrd_dcomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    double diff_a, diff_d, diff_e, diff_tauq, diff_taup;
    double threshold;
    void *hModule, *dModule;
	int min_mn;
    
    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char jobz; // Must be 'A', 'S', 'O', or 'N'. 

    lapack_int m; // The number of rows of the matrix A
    lapack_int n; // The number of columns of the matrix A 
    lapack_int lda; //  The leading dimension of a

    /* Input / Output parameters */
    lapack_complex_double* a, *aref; // contains m-by-n matrix A.

    /* Output parameters */
    // Below buffers contain the o/p orthogonal/unitary matrices
    double* d, *dref;
    double* e, *eref;
    lapack_complex_double* tauq, *tauqref;
    lapack_complex_double* taup, *taupref;
   
    /*Return Values */
    int info, inforef;

   public:
      gebrd_dcomplex_parameters (int matrix_layout_i, lapack_int m, lapack_int n);
      ~gebrd_dcomplex_parameters ();
};

/* Constructor definition  gebrd lapack_complex_double_common_parameters */
gebrd_dcomplex_parameters:: gebrd_dcomplex_parameters (int matrix_layout_i,
                    lapack_int m_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    n = n_i;
    m = m_i;
	min_mn = (m<n)?m:n;
    

    if (matrix_layout == LAPACK_ROW_MAJOR)
    {
        lda = n;
    }   else
    {
        lda = m;
    }

    // Initialize with default values
    hModule = NULL;
    dModule = NULL;
    diff_a = 0;
    diff_d = 0;
    diff_e = 0;
    diff_tauq = 0;
    diff_taup = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n gebrd lapack_complex_double: matrix_layout: %d   \t \
m: %d \t n: %d  \t lda: %d  \n",
matrix_layout, m_i, n_i, lda);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, m*n );
    lapacke_gtest_alloc_double_buffer_pair( &d, &dref, min_mn );
    lapacke_gtest_alloc_double_buffer_pair( &e, &eref, min_mn );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &tauq, &tauqref, min_mn );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &taup, &taupref, min_mn );

    if( (a==NULL) || (aref==NULL) ||  \
        (d==NULL) || (dref==NULL) || \
        (e==NULL) || (eref==NULL) || \
        (tauq==NULL) || (tauqref==NULL) || \
        (taup==NULL) || (taupref==NULL) ){
       EXPECT_FALSE( true) << "gebrd_dcomplex_parameters object: malloc error. Exiting ";
       gebrd_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, m*n );
    lapacke_gtest_init_double_buffer_pair_with_constant( d, dref, min_mn, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant( e, eref, min_mn, 0.0);
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( tauq, tauqref, min_mn, 0.0);
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( taup, taupref, min_mn, 0.0);

   } /* end of Constructor  */

/* Destructor definition  'gebrd_dcomplex_common_parameters' */
gebrd_dcomplex_parameters :: ~gebrd_dcomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   gebrd_free();
} 

//  Test fixture class definition
class zgebrd_test  : public  ::testing::Test {
public:
   gebrd_dcomplex_parameters  *zgebrd_obj;
   void SetUp();  
   void TearDown () { delete zgebrd_obj; }
};

void zgebrd_test::SetUp()
{
    /* LAPACKE zgebrd prototype */
    typedef int (*Fptr_NL_LAPACKE_zgebrd) ( int matrix_layout,
		lapack_int m, lapack_int n, lapack_complex_double* a, lapack_int lda,
		double* d, double* e, lapack_complex_double* tauq, lapack_complex_double* taup);
            
    Fptr_NL_LAPACKE_zgebrd ZGEBRD;

    zgebrd_obj = new  gebrd_dcomplex_parameters( svd_paramslist[idx].matrix_layout,
                                         svd_paramslist[idx].m,
                                         svd_paramslist[idx].n );

    zgebrd_obj->threshold = svd_paramslist[idx].svd_threshold;
    idx = Circular_Increment_Index(idx);
    zgebrd_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgebrd_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgebrd_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgebrd_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZGEBRD = (Fptr_NL_LAPACKE_zgebrd)dlsym(zgebrd_obj->hModule, "LAPACKE_zgebrd");
    ASSERT_TRUE(ZGEBRD != NULL) << "failed to ppt the Netlib LAPACKE_zgebrd symbol";

    /* Compute libflame's Lapacke o/p  */
    zgebrd_obj->info = LAPACKE_zgebrd(  zgebrd_obj->matrix_layout,
                                    zgebrd_obj->m,
                                    zgebrd_obj->n,
                                    zgebrd_obj->a,
                                    zgebrd_obj->lda,
                                    zgebrd_obj->d,
                                    zgebrd_obj->e,
                                    zgebrd_obj->tauq,
                                    zgebrd_obj->taup
                                    );

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    zgebrd_obj->inforef = ZGEBRD( zgebrd_obj->matrix_layout,
                                    zgebrd_obj->m,
                                    zgebrd_obj->n,
                                    zgebrd_obj->aref,
                                    zgebrd_obj->lda,
                                    zgebrd_obj->dref,
                                    zgebrd_obj->eref,
                                    zgebrd_obj->tauqref,
                                    zgebrd_obj->taupref
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    zgebrd_obj->diff_d =  computeDiff_d( zgebrd_obj->min_mn, 
                zgebrd_obj->d, zgebrd_obj->dref );

    zgebrd_obj->diff_a =  computeDiff_z( (zgebrd_obj->m)*(zgebrd_obj->n), 
                zgebrd_obj->a, zgebrd_obj->aref );

    zgebrd_obj->diff_e =  computeDiff_d( zgebrd_obj->min_mn-1, 
                zgebrd_obj->e, zgebrd_obj->eref );

    zgebrd_obj->diff_tauq =  computeDiff_z( zgebrd_obj->min_mn, 
                zgebrd_obj->tauq, zgebrd_obj->tauqref );

    zgebrd_obj->diff_taup =  computeDiff_z( zgebrd_obj->min_mn, 
                zgebrd_obj->taup, zgebrd_obj->taupref );


#if LAPACKE_TEST_VERBOSE
   printf(" \n gebrd lapack_complex_double: \n diff_a: %f \n diff_d: %f \n \
diff_e: %f \n diff_tauq: %f \n diff_taup: %f \n info: %d \t inforef: %d\n \
gebrd_threshold: %f  \t  \n",
       zgebrd_obj->diff_a, zgebrd_obj->diff_d, zgebrd_obj->diff_e,
       zgebrd_obj->diff_tauq, zgebrd_obj->diff_taup, zgebrd_obj->info, 
	   zgebrd_obj->inforef, zgebrd_obj->threshold);
#endif
}

TEST_F(zgebrd_test, zgebrd_1) {
    EXPECT_NEAR(0.0, zgebrd_obj->diff_e, zgebrd_obj->threshold);
    EXPECT_NEAR(0.0, zgebrd_obj->diff_d, zgebrd_obj->threshold);
    EXPECT_NEAR(0.0, zgebrd_obj->diff_taup, zgebrd_obj->threshold);
    EXPECT_NEAR(0.0, zgebrd_obj->diff_tauq, zgebrd_obj->threshold);
}

TEST_F(zgebrd_test, zgebrd_2) {
    EXPECT_NEAR(0.0, zgebrd_obj->diff_e, zgebrd_obj->threshold);
    EXPECT_NEAR(0.0, zgebrd_obj->diff_d, zgebrd_obj->threshold);
    EXPECT_NEAR(0.0, zgebrd_obj->diff_taup, zgebrd_obj->threshold);
    EXPECT_NEAR(0.0, zgebrd_obj->diff_tauq, zgebrd_obj->threshold);
}
TEST_F(zgebrd_test, zgebrd_3) {
    EXPECT_NEAR(0.0, zgebrd_obj->diff_e, zgebrd_obj->threshold);
    EXPECT_NEAR(0.0, zgebrd_obj->diff_d, zgebrd_obj->threshold);
    EXPECT_NEAR(0.0, zgebrd_obj->diff_taup, zgebrd_obj->threshold);
    EXPECT_NEAR(0.0, zgebrd_obj->diff_tauq, zgebrd_obj->threshold);
}
TEST_F(zgebrd_test, zgebrd_4) {
    EXPECT_NEAR(0.0, zgebrd_obj->diff_e, zgebrd_obj->threshold);
    EXPECT_NEAR(0.0, zgebrd_obj->diff_d, zgebrd_obj->threshold);
    EXPECT_NEAR(0.0, zgebrd_obj->diff_taup, zgebrd_obj->threshold);
    EXPECT_NEAR(0.0, zgebrd_obj->diff_tauq, zgebrd_obj->threshold);
}