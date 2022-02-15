#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"
#define LAPACKE_TEST_VERBOSE  (1)
#define gesvd_free() \
    if (a!=NULL)        free(a); \
    if (aref!=NULL)     free(aref); \
    if (s!=NULL)        free(s); \
    if (sref!=NULL)     free(sref); \
    if (u!=NULL)        free(u); \
    if (uref!=NULL)     free(uref); \
    if (vt!=NULL)        free(vt); \
    if (vtref!=NULL)     free(vtref); \
    if (superb!=NULL)        free(superb); \
    if (superbref!=NULL)     free(superbref)


// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin gesvd_float_common_parameters  class definition */
class gesvd_float_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_a, diff_s, diff_u, diff_vt, diff_superb;
    int min_mn;
    float threshold;
    void *hModule, *dModule;
    
    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char jobu; // Must be 'A', 'S', 'O', or 'N'. 
    char jobvt; // Must be 'A', 'S', 'O', or 'N'. 

    lapack_int m; // The number of rows of the matrix A
    lapack_int n; // The number of columns of the matrix A 
    lapack_int lda; //  The leading dimension of a

    lapack_int ldu; // leading dimension of the output array u . ldu≥ max(1, m)
    lapack_int ldvt; // The leading dimension of the output array vt . ldvt≥ max(1, p) 
    /* Input / Output parameters */
    float* a, *aref; // contains m-by-n matrix A.

    /* Output parameters */
    float *superb, *superbref; // the dimension of subblocks.
    // Below buffers contain the o/p orthogonal/unitary matrices
    float* u, *uref;
    float* vt, *vtref;
    float* s, *sref;
    
    /*Return Values */
    int info, inforef;

   public:
      gesvd_float_parameters (int matrix_layout_i, char jobu, char jobvt,
                                  lapack_int m, lapack_int n);
      ~gesvd_float_parameters ();
};

/* Constructor definition  gesvd float_common_parameters */
gesvd_float_parameters:: gesvd_float_parameters (int matrix_layout_i,
                    char jobu_i, char jobvt_i, lapack_int m_i,  
                    lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    jobu = jobu_i;
    jobvt = jobvt_i;
    n = n_i;
    m = m_i;
    min_mn = (m<n)?m:n;
    
    ldu = m;
    ldvt = n;

    if (matrix_layout == LAPACK_ROW_MAJOR)
    {
        lda = n;
        if( jobu== 'S')  ldu = min_mn; 

    }   else
    {
        lda = m;
        if( jobu== 'S')  ldvt = min_mn; 
    }

    // Initialize with default values
    hModule = NULL;
    dModule = NULL;
    diff_a = 0;
    diff_s = 0;
    diff_u = 0;
    diff_vt = 0;
    diff_superb = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n gesvd float: matrix_layout: %d   jobu: %c \t \
jobvt: %c \t m: %d \t n: %d \n",
matrix_layout, jobu_i, jobvt_i, m_i, n_i);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, m*n );
    lapacke_gtest_alloc_float_buffer_pair( &s, &sref, min_mn );
    lapacke_gtest_alloc_float_buffer_pair( &u, &uref, m*m );
    lapacke_gtest_alloc_float_buffer_pair( &vt, &vtref, n*n );
    lapacke_gtest_alloc_float_buffer_pair( &superb, &superbref, min_mn );

    if( (a==NULL) || (aref==NULL) ||  \
        (s==NULL) || (sref==NULL) || \
        (u==NULL) || (uref==NULL) || \
        (vt==NULL) || (vtref==NULL) || \
        (superb==NULL) || (superbref==NULL)  ){
       EXPECT_FALSE( true) << "gesvd_float_parameters object: malloc error. Exiting ";
       gesvd_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, m*n );
    lapacke_gtest_init_float_buffer_pair_with_constant( u, uref, m*m, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant( vt, vtref, n*n, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant( s, sref, min_mn, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant( superb, superbref, min_mn, 0.0);

   } /* end of Constructor  */

/* Destructor definition  'gesvd_float_common_parameters' */
gesvd_float_parameters :: ~gesvd_float_parameters ()
{
   /* De-Allocate memory for the input matrices */
   gesvd_free();
} 

//  Test fixture class definition
class sgesvd_test  : public  ::testing::Test {
public:
   gesvd_float_parameters  *sgesvd_obj;
   void SetUp();  
   void TearDown () { delete sgesvd_obj; }
};

void sgesvd_test::SetUp()
{
    /* LAPACKE sgesvd prototype */
    typedef int (*Fptr_NL_LAPACKE_sgesvd) ( int matrix_layout, char jobu,
        char jobvtt, lapack_int m, lapack_int n, float* a, lapack_int lda,
        float* s, float* u, lapack_int ldu, float* vt, lapack_int ldvtt,
        float* superb);
            
    Fptr_NL_LAPACKE_sgesvd SGESVD;

    sgesvd_obj = new  gesvd_float_parameters( svd_paramslist[idx].matrix_layout,
                                         svd_paramslist[idx].jobu_gesvd,
                                         svd_paramslist[idx].jobvt_gesvd,
                                         svd_paramslist[idx].m,
                                         svd_paramslist[idx].n );

    idx = Circular_Increment_Index(idx);
    sgesvd_obj->threshold = svd_paramslist[idx].svd_threshold;
    sgesvd_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgesvd_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgesvd_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgesvd_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SGESVD = (Fptr_NL_LAPACKE_sgesvd)dlsym(sgesvd_obj->hModule, "LAPACKE_sgesvd");
    ASSERT_TRUE(SGESVD != NULL) << "failed to ppt the Netlib LAPACKE_sgesvd symbol";


    /* Compute libflame's Lapacke o/p  */
    sgesvd_obj->info = LAPACKE_sgesvd(  sgesvd_obj->matrix_layout,
                                    sgesvd_obj->jobu,
                                    sgesvd_obj->jobvt,
                                    sgesvd_obj->m,
                                    sgesvd_obj->n,
                                    sgesvd_obj->a,
                                    sgesvd_obj->lda,
                                    sgesvd_obj->s,
                                    sgesvd_obj->u,
                                    sgesvd_obj->ldu,
                                    sgesvd_obj->vt,
                                    sgesvd_obj->ldvt,
                                    sgesvd_obj->superb
                                    );

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    sgesvd_obj->inforef = SGESVD( sgesvd_obj->matrix_layout,
                                    sgesvd_obj->jobu,
                                    sgesvd_obj->jobvt,
                                    sgesvd_obj->m,
                                    sgesvd_obj->n,
                                    sgesvd_obj->aref,
                                    sgesvd_obj->lda,
                                    sgesvd_obj->sref,
                                    sgesvd_obj->uref,
                                    sgesvd_obj->ldu,
                                    sgesvd_obj->vtref,
                                    sgesvd_obj->ldvt,
                                    sgesvd_obj->superbref
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    sgesvd_obj->diff_a =  computeDiff_s( (sgesvd_obj->m)*(sgesvd_obj->n), 
                sgesvd_obj->a, sgesvd_obj->aref );

    sgesvd_obj->diff_u =  computeDiff_s( (sgesvd_obj->m)*(sgesvd_obj->m), 
                sgesvd_obj->u, sgesvd_obj->uref );

    sgesvd_obj->diff_vt =  computeDiff_s( (sgesvd_obj->n)*(sgesvd_obj->n), 
                sgesvd_obj->vt, sgesvd_obj->vtref );

    sgesvd_obj->diff_s =  computeDiff_s( sgesvd_obj->min_mn, 
                sgesvd_obj->s, sgesvd_obj->sref );

    sgesvd_obj->diff_superb =  computeDiff_s( sgesvd_obj->min_mn, 
                sgesvd_obj->superb, sgesvd_obj->superbref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n gesvd float: \n diff_a: %f \n diff_s: %f \n \
diff_u: %f \n diff_vt: %f \n diff_superb: %f info: %d \t inforef: %d\n ",
       sgesvd_obj->diff_a, sgesvd_obj->diff_s, sgesvd_obj->diff_u,
       sgesvd_obj->diff_vt, sgesvd_obj->diff_superb,
	   sgesvd_obj->info, sgesvd_obj->inforef);
#endif
}

TEST_F(sgesvd_test, sgesvd_1) {
    EXPECT_NEAR(0.0, sgesvd_obj->diff_a, sgesvd_obj->threshold);
    EXPECT_NEAR(0.0, sgesvd_obj->diff_s, sgesvd_obj->threshold);
    EXPECT_NEAR(0.0, sgesvd_obj->diff_u, sgesvd_obj->threshold);
    EXPECT_NEAR(0.0, sgesvd_obj->diff_vt, sgesvd_obj->threshold);
    EXPECT_NEAR(0.0, sgesvd_obj->diff_superb, sgesvd_obj->threshold);
}

TEST_F(sgesvd_test, sgesvd_2) {
    EXPECT_NEAR(0.0, sgesvd_obj->diff_a, sgesvd_obj->threshold);
    EXPECT_NEAR(0.0, sgesvd_obj->diff_s, sgesvd_obj->threshold);
    EXPECT_NEAR(0.0, sgesvd_obj->diff_u, sgesvd_obj->threshold);
    EXPECT_NEAR(0.0, sgesvd_obj->diff_vt, sgesvd_obj->threshold);
    EXPECT_NEAR(0.0, sgesvd_obj->diff_superb, sgesvd_obj->threshold);}

TEST_F(sgesvd_test, sgesvd_3) {
    EXPECT_NEAR(0.0, sgesvd_obj->diff_a, sgesvd_obj->threshold);
    EXPECT_NEAR(0.0, sgesvd_obj->diff_s, sgesvd_obj->threshold);
    EXPECT_NEAR(0.0, sgesvd_obj->diff_u, sgesvd_obj->threshold);
    EXPECT_NEAR(0.0, sgesvd_obj->diff_vt, sgesvd_obj->threshold);
    EXPECT_NEAR(0.0, sgesvd_obj->diff_superb, sgesvd_obj->threshold);}

TEST_F(sgesvd_test, sgesvd_4) {
    EXPECT_NEAR(0.0, sgesvd_obj->diff_a, sgesvd_obj->threshold);
    EXPECT_NEAR(0.0, sgesvd_obj->diff_s, sgesvd_obj->threshold);
    EXPECT_NEAR(0.0, sgesvd_obj->diff_u, sgesvd_obj->threshold);
    EXPECT_NEAR(0.0, sgesvd_obj->diff_vt, sgesvd_obj->threshold);
    EXPECT_NEAR(0.0, sgesvd_obj->diff_superb, sgesvd_obj->threshold);}


/* Begin gesvd_double_common_parameters  class definition */
class gesvd_double_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    double diff_a, diff_s, diff_u, diff_vt, diff_superb;
    int min_mn;
    double threshold;
    void *hModule, *dModule;
    
    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char jobu; // Must be 'A', 'S', 'O', or 'N'. 
    char jobvt; // Must be 'A', 'S', 'O', or 'N'. 

    lapack_int m; // The number of rows of the matrix A
    lapack_int n; // The number of columns of the matrix A 
    lapack_int lda; //  The leading dimension of a

    lapack_int ldu; // leading dimension of the output array u . ldu≥ max(1, m)
    lapack_int ldvt; // The leading dimension of the output array vt . ldvt≥ max(1, p) 
    /* Input / Output parameters */
    double* a, *aref; // contains m-by-n matrix A.

    /* Output parameters */
    double *superb, *superbref; // the dimension of subblocks.
    // Below buffers contain the o/p orthogonal/unitary matrices
    double* u, *uref;
    double* vt, *vtref;
    double* s, *sref;
    
    /*Return Values */
    int info, inforef;

   public:
      gesvd_double_parameters (int matrix_layout_i, char jobu, char jobvt,
                                  lapack_int m, lapack_int n);
      ~gesvd_double_parameters ();
};

/* Constructor definition  gesvd double_common_parameters */
gesvd_double_parameters:: gesvd_double_parameters (int matrix_layout_i,
                    char jobu_i, char jobvt_i, lapack_int m_i,  
                    lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    jobu = jobu_i;
    jobvt = jobvt_i;
    n = n_i;
    m = m_i;
    min_mn = (m<n)?m:n;
    
    ldu = m;
    ldvt = n;

    if (matrix_layout == LAPACK_ROW_MAJOR)
    {
        lda = n;
        if( jobu== 'S')  ldu = min_mn; 

    }   else
    {
        lda = m;
        if( jobu== 'S')  ldvt = min_mn; 
    }

    // Initialize with default values
    hModule = NULL;
    dModule = NULL;
    diff_a = 0;
    diff_s = 0;
    diff_u = 0;
    diff_vt = 0;
    diff_superb = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n gesvd double: matrix_layout: %d   jobu: %c \t \
jobvt: %c \t m: %d \t n: %d \n",
matrix_layout, jobu_i, jobvt_i, m_i, n_i);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, m*n );
    lapacke_gtest_alloc_double_buffer_pair( &s, &sref, min_mn );
    lapacke_gtest_alloc_double_buffer_pair( &u, &uref, m*m );
    lapacke_gtest_alloc_double_buffer_pair( &vt, &vtref, n*n );
    lapacke_gtest_alloc_double_buffer_pair( &superb, &superbref, min_mn );

    if( (a==NULL) || (aref==NULL) ||  \
        (s==NULL) || (sref==NULL) || \
        (u==NULL) || (uref==NULL) || \
        (vt==NULL) || (vtref==NULL) || \
        (superb==NULL) || (superbref==NULL)  ){
       EXPECT_FALSE( true) << "gesvd_double_parameters object: malloc error. Exiting ";
       gesvd_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, m*n );
    lapacke_gtest_init_double_buffer_pair_with_constant( u, uref, m*m, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant( vt, vtref, n*n, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant( s, sref, min_mn, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant( superb, superbref, min_mn, 0.0);

   } /* end of Constructor  */

/* Destructor definition  'gesvd_double_common_parameters' */
gesvd_double_parameters :: ~gesvd_double_parameters ()
{
   /* De-Allocate memory for the input matrices */
   gesvd_free();
} 

//  Test fixture class definition
class dgesvd_test  : public  ::testing::Test {
public:
   gesvd_double_parameters  *dgesvd_obj;
   void SetUp();  
   void TearDown () { delete dgesvd_obj; }
};

void dgesvd_test::SetUp()
{
    /* LAPACKE dgesvd prototype */
    typedef int (*Fptr_NL_LAPACKE_dgesvd) ( int matrix_layout, char jobu,
        char jobvtt, lapack_int m, lapack_int n, double* a, lapack_int lda,
        double* s, double* u, lapack_int ldu, double* vt, lapack_int ldvtt,
        double* superb);
            
    Fptr_NL_LAPACKE_dgesvd DGESVD;

    dgesvd_obj = new  gesvd_double_parameters( svd_paramslist[idx].matrix_layout,
                                         svd_paramslist[idx].jobu_gesvd,
                                         svd_paramslist[idx].jobvt_gesvd,
                                         svd_paramslist[idx].m,
                                         svd_paramslist[idx].n );

    idx = Circular_Increment_Index(idx);
    dgesvd_obj->threshold = svd_paramslist[idx].svd_threshold;
    dgesvd_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgesvd_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgesvd_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgesvd_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DGESVD = (Fptr_NL_LAPACKE_dgesvd)dlsym(dgesvd_obj->hModule, "LAPACKE_dgesvd");
    ASSERT_TRUE(DGESVD != NULL) << "failed to ppt the Netlib LAPACKE_dgesvd symbol";


    /* Compute libflame's Lapacke o/p  */
    dgesvd_obj->info = LAPACKE_dgesvd(  dgesvd_obj->matrix_layout,
                                    dgesvd_obj->jobu,
                                    dgesvd_obj->jobvt,
                                    dgesvd_obj->m,
                                    dgesvd_obj->n,
                                    dgesvd_obj->a,
                                    dgesvd_obj->lda,
                                    dgesvd_obj->s,
                                    dgesvd_obj->u,
                                    dgesvd_obj->ldu,
                                    dgesvd_obj->vt,
                                    dgesvd_obj->ldvt,
                                    dgesvd_obj->superb
                                    );

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    dgesvd_obj->inforef = DGESVD( dgesvd_obj->matrix_layout,
                                    dgesvd_obj->jobu,
                                    dgesvd_obj->jobvt,
                                    dgesvd_obj->m,
                                    dgesvd_obj->n,
                                    dgesvd_obj->aref,
                                    dgesvd_obj->lda,
                                    dgesvd_obj->sref,
                                    dgesvd_obj->uref,
                                    dgesvd_obj->ldu,
                                    dgesvd_obj->vtref,
                                    dgesvd_obj->ldvt,
                                    dgesvd_obj->superbref
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    dgesvd_obj->diff_a =  computeDiff_d( (dgesvd_obj->m)*(dgesvd_obj->n), 
                dgesvd_obj->a, dgesvd_obj->aref );

    dgesvd_obj->diff_u =  computeDiff_d( (dgesvd_obj->m)*(dgesvd_obj->m), 
                dgesvd_obj->u, dgesvd_obj->uref );

    dgesvd_obj->diff_vt =  computeDiff_d( (dgesvd_obj->n)*(dgesvd_obj->n), 
                dgesvd_obj->vt, dgesvd_obj->vtref );

    dgesvd_obj->diff_s =  computeDiff_d( dgesvd_obj->min_mn, 
                dgesvd_obj->s, dgesvd_obj->sref );

    dgesvd_obj->diff_superb =  computeDiff_d( dgesvd_obj->min_mn, 
                dgesvd_obj->superb, dgesvd_obj->superbref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n gesvd double: \n diff_a: %f \n diff_s: %f \n \
diff_u: %f \n diff_vt: %f \n diff_superb: %f info: %d \t inforef: %d\n ",
       dgesvd_obj->diff_a, dgesvd_obj->diff_s, dgesvd_obj->diff_u,
       dgesvd_obj->diff_vt, dgesvd_obj->diff_superb,
	   dgesvd_obj->info, dgesvd_obj->inforef);
#endif
}

TEST_F(dgesvd_test, dgesvd_1) {
    EXPECT_NEAR(0.0, dgesvd_obj->diff_a, dgesvd_obj->threshold);
    EXPECT_NEAR(0.0, dgesvd_obj->diff_s, dgesvd_obj->threshold);
    EXPECT_NEAR(0.0, dgesvd_obj->diff_u, dgesvd_obj->threshold);
    EXPECT_NEAR(0.0, dgesvd_obj->diff_vt, dgesvd_obj->threshold);
    EXPECT_NEAR(0.0, dgesvd_obj->diff_superb, dgesvd_obj->threshold);
}

TEST_F(dgesvd_test, dgesvd_2) {
    EXPECT_NEAR(0.0, dgesvd_obj->diff_a, dgesvd_obj->threshold);
    EXPECT_NEAR(0.0, dgesvd_obj->diff_s, dgesvd_obj->threshold);
    EXPECT_NEAR(0.0, dgesvd_obj->diff_u, dgesvd_obj->threshold);
    EXPECT_NEAR(0.0, dgesvd_obj->diff_vt, dgesvd_obj->threshold);
    EXPECT_NEAR(0.0, dgesvd_obj->diff_superb, dgesvd_obj->threshold);}

TEST_F(dgesvd_test, dgesvd_3) {
    EXPECT_NEAR(0.0, dgesvd_obj->diff_a, dgesvd_obj->threshold);
    EXPECT_NEAR(0.0, dgesvd_obj->diff_s, dgesvd_obj->threshold);
    EXPECT_NEAR(0.0, dgesvd_obj->diff_u, dgesvd_obj->threshold);
    EXPECT_NEAR(0.0, dgesvd_obj->diff_vt, dgesvd_obj->threshold);
    EXPECT_NEAR(0.0, dgesvd_obj->diff_superb, dgesvd_obj->threshold);}

TEST_F(dgesvd_test, dgesvd_4) {
    EXPECT_NEAR(0.0, dgesvd_obj->diff_a, dgesvd_obj->threshold);
    EXPECT_NEAR(0.0, dgesvd_obj->diff_s, dgesvd_obj->threshold);
    EXPECT_NEAR(0.0, dgesvd_obj->diff_u, dgesvd_obj->threshold);
    EXPECT_NEAR(0.0, dgesvd_obj->diff_vt, dgesvd_obj->threshold);
    EXPECT_NEAR(0.0, dgesvd_obj->diff_superb, dgesvd_obj->threshold);}


/* Begin gesvd_scomplex_common_parameters  class definition */
class gesvd_scomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_a, diff_s, diff_u, diff_vt, diff_superb;
    int min_mn;
    float threshold;
    void *hModule, *dModule;
    
    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char jobu; // Must be 'A', 'S', 'O', or 'N'. 
    char jobvt; // Must be 'A', 'S', 'O', or 'N'. 

    lapack_int m; // The number of rows of the matrix A
    lapack_int n; // The number of columns of the matrix A 
    lapack_int lda; //  The leading dimension of a

    lapack_int ldu; // leading dimension of the output array u . ldu≥ max(1, m)
    lapack_int ldvt; // The leading dimension of the output array vt . ldvt≥ max(1, p) 
    /* Input / Output parameters */
    lapack_complex_float* a, *aref; // contains m-by-n matrix A.

    /* Output parameters */
    float *superb, *superbref; // the dimension of subblocks.
    // Below buffers contain the o/p orthogonal/unitary matrices
    lapack_complex_float* u, *uref;
    lapack_complex_float* vt, *vtref;
    float* s, *sref;
    
    /*Return Values */
    int info, inforef;

   public:
      gesvd_scomplex_parameters (int matrix_layout_i, char jobu, char jobvt,
                                  lapack_int m, lapack_int n);
      ~gesvd_scomplex_parameters ();
};

/* Constructor definition  gesvd lapack_complex_float_common_parameters */
gesvd_scomplex_parameters:: gesvd_scomplex_parameters (int matrix_layout_i,
                    char jobu_i, char jobvt_i, lapack_int m_i,  
                    lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    jobu = jobu_i;
    jobvt = jobvt_i;
    n = n_i;
    m = m_i;
    min_mn = (m<n)?m:n;
    
    ldu = m;
    ldvt = n;

    if (matrix_layout == LAPACK_ROW_MAJOR)
    {
        lda = n;
        if( jobu== 'S')  ldu = min_mn; 

    }   else
    {
        lda = m;
        if( jobu== 'S')  ldvt = min_mn; 
    }

    // Initialize with default values
    hModule = NULL;
    dModule = NULL;
    diff_a = 0;
    diff_s = 0;
    diff_u = 0;
    diff_vt = 0;
    diff_superb = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n gesvd lapack_complex_float: matrix_layout: %d   jobu: %c \t \
jobvt: %c \t m: %d \t n: %d \n",
matrix_layout, jobu_i, jobvt_i, m_i, n_i);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, m*n );
    lapacke_gtest_alloc_float_buffer_pair( &s, &sref, min_mn );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &u, &uref, m*m );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &vt, &vtref, n*n );
    lapacke_gtest_alloc_float_buffer_pair( &superb, &superbref, min_mn );

    if( (a==NULL) || (aref==NULL) ||  \
        (s==NULL) || (sref==NULL) || \
        (u==NULL) || (uref==NULL) || \
        (vt==NULL) || (vtref==NULL) || \
        (superb==NULL) || (superbref==NULL)  ){
       EXPECT_FALSE( true) << "gesvd_scomplex_parameters object: malloc error. Exiting ";
       gesvd_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, m*n );
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( u, uref, m*m, 0.0);
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( vt, vtref, n*n, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant( s, sref, min_mn, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant( superb, superbref, min_mn, 0.0);

   } /* end of Constructor  */

/* Destructor definition  'gesvd_scomplex_common_parameters' */
gesvd_scomplex_parameters :: ~gesvd_scomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   gesvd_free();
} 

//  Test fixture class definition
class cgesvd_test  : public  ::testing::Test {
public:
   gesvd_scomplex_parameters  *cgesvd_obj;
   void SetUp();  
   void TearDown () { delete cgesvd_obj; }
};

void cgesvd_test::SetUp()
{
    /* LAPACKE cgesvd prototype */
    typedef int (*Fptr_NL_LAPACKE_cgesvd) ( int matrix_layout, char jobu,
		char jobvt, lapack_int m, lapack_int n, lapack_complex_float* a,
		lapack_int lda, float* s, lapack_complex_float* u, lapack_int ldu,
		lapack_complex_float* vt, lapack_int ldvt, float* superb);
            
    Fptr_NL_LAPACKE_cgesvd CGESVD;

    cgesvd_obj = new  gesvd_scomplex_parameters( svd_paramslist[idx].matrix_layout,
                                         svd_paramslist[idx].jobu_gesvd,
                                         svd_paramslist[idx].jobvt_gesvd,
                                         svd_paramslist[idx].m,
                                         svd_paramslist[idx].n );

    idx = Circular_Increment_Index(idx);
    cgesvd_obj->threshold = svd_paramslist[idx].svd_threshold;
    cgesvd_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgesvd_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgesvd_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgesvd_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CGESVD = (Fptr_NL_LAPACKE_cgesvd)dlsym(cgesvd_obj->hModule, "LAPACKE_cgesvd");
    ASSERT_TRUE(CGESVD != NULL) << "failed to ppt the Netlib LAPACKE_cgesvd symbol";


    /* Compute libflame's Lapacke o/p  */
    cgesvd_obj->info = LAPACKE_cgesvd(  cgesvd_obj->matrix_layout,
                                    cgesvd_obj->jobu,
                                    cgesvd_obj->jobvt,
                                    cgesvd_obj->m,
                                    cgesvd_obj->n,
                                    cgesvd_obj->a,
                                    cgesvd_obj->lda,
                                    cgesvd_obj->s,
                                    cgesvd_obj->u,
                                    cgesvd_obj->ldu,
                                    cgesvd_obj->vt,
                                    cgesvd_obj->ldvt,
                                    cgesvd_obj->superb
                                    );

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    cgesvd_obj->inforef = CGESVD( cgesvd_obj->matrix_layout,
                                    cgesvd_obj->jobu,
                                    cgesvd_obj->jobvt,
                                    cgesvd_obj->m,
                                    cgesvd_obj->n,
                                    cgesvd_obj->aref,
                                    cgesvd_obj->lda,
                                    cgesvd_obj->sref,
                                    cgesvd_obj->uref,
                                    cgesvd_obj->ldu,
                                    cgesvd_obj->vtref,
                                    cgesvd_obj->ldvt,
                                    cgesvd_obj->superbref
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    cgesvd_obj->diff_a =  computeDiff_c( (cgesvd_obj->m)*(cgesvd_obj->n), 
                cgesvd_obj->a, cgesvd_obj->aref );

    cgesvd_obj->diff_u =  computeDiff_c( (cgesvd_obj->m)*(cgesvd_obj->m), 
                cgesvd_obj->u, cgesvd_obj->uref );

    cgesvd_obj->diff_vt =  computeDiff_c( (cgesvd_obj->n)*(cgesvd_obj->n), 
                cgesvd_obj->vt, cgesvd_obj->vtref );

    cgesvd_obj->diff_s =  computeDiff_s( cgesvd_obj->min_mn, 
                cgesvd_obj->s, cgesvd_obj->sref );

    cgesvd_obj->diff_superb =  computeDiff_s( cgesvd_obj->min_mn, 
                cgesvd_obj->superb, cgesvd_obj->superbref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n gesvd lapack_complex_float: \n diff_a: %f \n diff_s: %f \n \
diff_u: %f \n diff_vt: %f \n diff_superb: %f info: %d \t inforef: %d\n ",
       cgesvd_obj->diff_a, cgesvd_obj->diff_s, cgesvd_obj->diff_u,
       cgesvd_obj->diff_vt, cgesvd_obj->diff_superb,
	   cgesvd_obj->info, cgesvd_obj->inforef);
#endif
}

TEST_F(cgesvd_test, cgesvd_1) {
    EXPECT_NEAR(0.0, cgesvd_obj->diff_a, cgesvd_obj->threshold);
    EXPECT_NEAR(0.0, cgesvd_obj->diff_s, cgesvd_obj->threshold);
    EXPECT_NEAR(0.0, cgesvd_obj->diff_u, cgesvd_obj->threshold);
    EXPECT_NEAR(0.0, cgesvd_obj->diff_vt, cgesvd_obj->threshold);
    EXPECT_NEAR(0.0, cgesvd_obj->diff_superb, cgesvd_obj->threshold);
}

TEST_F(cgesvd_test, cgesvd_2) {
    EXPECT_NEAR(0.0, cgesvd_obj->diff_a, cgesvd_obj->threshold);
    EXPECT_NEAR(0.0, cgesvd_obj->diff_s, cgesvd_obj->threshold);
    EXPECT_NEAR(0.0, cgesvd_obj->diff_u, cgesvd_obj->threshold);
    EXPECT_NEAR(0.0, cgesvd_obj->diff_vt, cgesvd_obj->threshold);
    EXPECT_NEAR(0.0, cgesvd_obj->diff_superb, cgesvd_obj->threshold);}

TEST_F(cgesvd_test, cgesvd_3) {
    EXPECT_NEAR(0.0, cgesvd_obj->diff_a, cgesvd_obj->threshold);
    EXPECT_NEAR(0.0, cgesvd_obj->diff_s, cgesvd_obj->threshold);
    EXPECT_NEAR(0.0, cgesvd_obj->diff_u, cgesvd_obj->threshold);
    EXPECT_NEAR(0.0, cgesvd_obj->diff_vt, cgesvd_obj->threshold);
    EXPECT_NEAR(0.0, cgesvd_obj->diff_superb, cgesvd_obj->threshold);}

TEST_F(cgesvd_test, cgesvd_4) {
    EXPECT_NEAR(0.0, cgesvd_obj->diff_a, cgesvd_obj->threshold);
    EXPECT_NEAR(0.0, cgesvd_obj->diff_s, cgesvd_obj->threshold);
    EXPECT_NEAR(0.0, cgesvd_obj->diff_u, cgesvd_obj->threshold);
    EXPECT_NEAR(0.0, cgesvd_obj->diff_vt, cgesvd_obj->threshold);
    EXPECT_NEAR(0.0, cgesvd_obj->diff_superb, cgesvd_obj->threshold);
}

/* Begin gesvd_dcomplex_common_parameters  class definition */
class gesvd_dcomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    double diff_a, diff_s, diff_u, diff_vt, diff_superb;
    int min_mn;
    float threshold;
    void *hModule, *dModule;
    
    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char jobu; // Must be 'A', 'S', 'O', or 'N'. 
    char jobvt; // Must be 'A', 'S', 'O', or 'N'. 

    lapack_int m; // The number of rows of the matrix A
    lapack_int n; // The number of columns of the matrix A 
    lapack_int lda; //  The leading dimension of a

    lapack_int ldu; // leading dimension of the output array u . ldu≥ max(1, m)
    lapack_int ldvt; // The leading dimension of the output array vt . ldvt≥ max(1, p) 
    /* Input / Output parameters */
    lapack_complex_double* a, *aref; // contains m-by-n matrix A.

    /* Output parameters */
    double *superb, *superbref; // the dimension of subblocks.
    // Below buffers contain the o/p orthogonal/unitary matrices
    lapack_complex_double* u, *uref;
    lapack_complex_double* vt, *vtref;
    double* s, *sref;
    
    /*Return Values */
    int info, inforef;

   public:
      gesvd_dcomplex_parameters (int matrix_layout_i, char jobu, char jobvt,
                                  lapack_int m, lapack_int n);
      ~gesvd_dcomplex_parameters ();
};

/* Constructor definition  gesvd lapack_complex_double_common_parameters */
gesvd_dcomplex_parameters:: gesvd_dcomplex_parameters (int matrix_layout_i,
                    char jobu_i, char jobvt_i, lapack_int m_i,  
                    lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    jobu = jobu_i;
    jobvt = jobvt_i;
    n = n_i;
    m = m_i;
    min_mn = (m<n)?m:n;
    
    ldu = m;
    ldvt = n;

    if (matrix_layout == LAPACK_ROW_MAJOR)
    {
        lda = n;
        if( jobu== 'S')  ldu = min_mn; 

    }   else
    {
        lda = m;
        if( jobu== 'S')  ldvt = min_mn; 
    }

    // Initialize with default values
    hModule = NULL;
    dModule = NULL;
    diff_a = 0;
    diff_s = 0;
    diff_u = 0;
    diff_vt = 0;
    diff_superb = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n gesvd lapack_complex_double: matrix_layout: %d   jobu: %c \t \
jobvt: %c \t m: %d \t n: %d \n",
matrix_layout, jobu_i, jobvt_i, m_i, n_i);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, m*n );
    lapacke_gtest_alloc_double_buffer_pair( &s, &sref, min_mn );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &u, &uref, m*m );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &vt, &vtref, n*n );
    lapacke_gtest_alloc_double_buffer_pair( &superb, &superbref, min_mn );

    if( (a==NULL) || (aref==NULL) ||  \
        (s==NULL) || (sref==NULL) || \
        (u==NULL) || (uref==NULL) || \
        (vt==NULL) || (vtref==NULL) || \
        (superb==NULL) || (superbref==NULL)  ){
       EXPECT_FALSE( true) << "gesvd_dcomplex_parameters object: malloc error. Exiting ";
       gesvd_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, m*n );
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( u, uref, m*m, 0.0);
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( vt, vtref, n*n, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant( s, sref, min_mn, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant( superb, superbref, min_mn, 0.0);

   } /* end of Constructor  */

/* Destructor definition  'gesvd_dcomplex_common_parameters' */
gesvd_dcomplex_parameters :: ~gesvd_dcomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   gesvd_free();
} 

//  Test fixture class definition
class zgesvd_test  : public  ::testing::Test {
public:
   gesvd_dcomplex_parameters  *zgesvd_obj;
   void SetUp();  
   void TearDown () { delete zgesvd_obj; }
};

void zgesvd_test::SetUp()
{
    /* LAPACKE zgesvd prototype */
    typedef int (*Fptr_NL_LAPACKE_zgesvd) ( int matrix_layout, char jobu,
		char jobvt, lapack_int m, lapack_int n, lapack_complex_double* a,
		lapack_int lda, double* s, lapack_complex_double* u, lapack_int ldu,
		lapack_complex_double* vt, lapack_int ldvt, double* superb);
            
    Fptr_NL_LAPACKE_zgesvd ZGESVD;

    zgesvd_obj = new  gesvd_dcomplex_parameters( svd_paramslist[idx].matrix_layout,
                                         svd_paramslist[idx].jobu_gesvd,
                                         svd_paramslist[idx].jobvt_gesvd,
                                         svd_paramslist[idx].m,
                                         svd_paramslist[idx].n );

    idx = Circular_Increment_Index(idx);
    zgesvd_obj->threshold = svd_paramslist[idx].svd_threshold;
    zgesvd_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgesvd_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgesvd_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgesvd_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZGESVD = (Fptr_NL_LAPACKE_zgesvd)dlsym(zgesvd_obj->hModule, "LAPACKE_zgesvd");
    ASSERT_TRUE(ZGESVD != NULL) << "failed to ppt the Netlib LAPACKE_zgesvd symbol";


    /* Compute libflame's Lapacke o/p  */
    zgesvd_obj->info = LAPACKE_zgesvd(  zgesvd_obj->matrix_layout,
                                    zgesvd_obj->jobu,
                                    zgesvd_obj->jobvt,
                                    zgesvd_obj->m,
                                    zgesvd_obj->n,
                                    zgesvd_obj->a,
                                    zgesvd_obj->lda,
                                    zgesvd_obj->s,
                                    zgesvd_obj->u,
                                    zgesvd_obj->ldu,
                                    zgesvd_obj->vt,
                                    zgesvd_obj->ldvt,
                                    zgesvd_obj->superb
                                    );

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    zgesvd_obj->inforef = ZGESVD( zgesvd_obj->matrix_layout,
                                    zgesvd_obj->jobu,
                                    zgesvd_obj->jobvt,
                                    zgesvd_obj->m,
                                    zgesvd_obj->n,
                                    zgesvd_obj->aref,
                                    zgesvd_obj->lda,
                                    zgesvd_obj->sref,
                                    zgesvd_obj->uref,
                                    zgesvd_obj->ldu,
                                    zgesvd_obj->vtref,
                                    zgesvd_obj->ldvt,
                                    zgesvd_obj->superbref
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    zgesvd_obj->diff_a =  computeDiff_z( (zgesvd_obj->m)*(zgesvd_obj->n), 
                zgesvd_obj->a, zgesvd_obj->aref );

    zgesvd_obj->diff_u =  computeDiff_z( (zgesvd_obj->m)*(zgesvd_obj->m), 
                zgesvd_obj->u, zgesvd_obj->uref );

    zgesvd_obj->diff_vt =  computeDiff_z( (zgesvd_obj->n)*(zgesvd_obj->n), 
                zgesvd_obj->vt, zgesvd_obj->vtref );

    zgesvd_obj->diff_s =  computeDiff_d( zgesvd_obj->min_mn, 
                zgesvd_obj->s, zgesvd_obj->sref );

    zgesvd_obj->diff_superb =  computeDiff_d( zgesvd_obj->min_mn, 
                zgesvd_obj->superb, zgesvd_obj->superbref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n gesvd lapack_complex_double: \n diff_a: %f \n diff_s: %f \n \
diff_u: %f \n diff_vt: %f \n diff_superb: %f info: %d \t inforef: %d\n ",
       zgesvd_obj->diff_a, zgesvd_obj->diff_s, zgesvd_obj->diff_u,
       zgesvd_obj->diff_vt, zgesvd_obj->diff_superb,
	   zgesvd_obj->info, zgesvd_obj->inforef);
#endif
}

TEST_F(zgesvd_test, zgesvd_1) {
    EXPECT_NEAR(0.0, zgesvd_obj->diff_a, zgesvd_obj->threshold);
    EXPECT_NEAR(0.0, zgesvd_obj->diff_s, zgesvd_obj->threshold);
    EXPECT_NEAR(0.0, zgesvd_obj->diff_u, zgesvd_obj->threshold);
    EXPECT_NEAR(0.0, zgesvd_obj->diff_vt, zgesvd_obj->threshold);
    EXPECT_NEAR(0.0, zgesvd_obj->diff_superb, zgesvd_obj->threshold);
}

TEST_F(zgesvd_test, zgesvd_2) {
    EXPECT_NEAR(0.0, zgesvd_obj->diff_a, zgesvd_obj->threshold);
    EXPECT_NEAR(0.0, zgesvd_obj->diff_s, zgesvd_obj->threshold);
    EXPECT_NEAR(0.0, zgesvd_obj->diff_u, zgesvd_obj->threshold);
    EXPECT_NEAR(0.0, zgesvd_obj->diff_vt, zgesvd_obj->threshold);
    EXPECT_NEAR(0.0, zgesvd_obj->diff_superb, zgesvd_obj->threshold);}

TEST_F(zgesvd_test, zgesvd_3) {
    EXPECT_NEAR(0.0, zgesvd_obj->diff_a, zgesvd_obj->threshold);
    EXPECT_NEAR(0.0, zgesvd_obj->diff_s, zgesvd_obj->threshold);
    EXPECT_NEAR(0.0, zgesvd_obj->diff_u, zgesvd_obj->threshold);
    EXPECT_NEAR(0.0, zgesvd_obj->diff_vt, zgesvd_obj->threshold);
    EXPECT_NEAR(0.0, zgesvd_obj->diff_superb, zgesvd_obj->threshold);}

TEST_F(zgesvd_test, zgesvd_4) {
    EXPECT_NEAR(0.0, zgesvd_obj->diff_a, zgesvd_obj->threshold);
    EXPECT_NEAR(0.0, zgesvd_obj->diff_s, zgesvd_obj->threshold);
    EXPECT_NEAR(0.0, zgesvd_obj->diff_u, zgesvd_obj->threshold);
    EXPECT_NEAR(0.0, zgesvd_obj->diff_vt, zgesvd_obj->threshold);
    EXPECT_NEAR(0.0, zgesvd_obj->diff_superb, zgesvd_obj->threshold);
}

