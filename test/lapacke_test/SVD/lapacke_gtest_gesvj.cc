#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"
#define LAPACKE_TEST_VERBOSE  (1)
#define gesvj_free() \
    if (a!=NULL)        free(a); \
    if (aref!=NULL)     free(aref); \
    if (sva!=NULL)        free(sva); \
    if (svaref!=NULL)     free(svaref); \
    if (v!=NULL)        free(v); \
    if (vref!=NULL)     free(vref)


// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin gesvj_float_common_parameters  class definition */
class gesvj_float_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_a, diff_sva, diff_v, diff_stat;
    float threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char joba; //  Must be 'L', 'U' or 'G'.
    char jobu; // Must be 'U', 'C' or 'N'.
    char jobv; // Must be 'V', 'A' or 'N'.

    lapack_int m; // The number of rows of the matrix A
    lapack_int n; // The number of columns of the matrix A 
    lapack_int lda; //  The leading dimension of a

    lapack_int mv; // Jacobi rotations applied to the first mv rows of v.
    lapack_int ldv; // The leading dimension of the output array v . ldv≥ max(1, p)
    
    /* Input / Output parameters */
    float* a, *aref; // contains m-by-n matrix A.
    
    /* Output parameters */
    float* sva, *svaref; //   the singular values of A.
    float* v, *vref; //  contains the n-by-n matrix of the right singular vectors.
    float stat[6], statref[6]; //  scaling factor
    
    /*Return Values */
    int info, inforef;

   public:
      gesvj_float_parameters (int matrix_layout_i, char joba, char jobu, 
                    char jobv, lapack_int m, lapack_int n, lapack_int mv );
      ~gesvj_float_parameters ();
};

/* Constructor definition  gesvj float_common_parameters */
gesvj_float_parameters:: gesvj_float_parameters (int matrix_layout_i, char joba_i, 
                          char jobu_i, char jobv_i, lapack_int m_i, lapack_int n_i,
                          lapack_int mv_i)
{
    matrix_layout = matrix_layout_i;
    jobu = jobu_i;
    jobv = jobv_i;
    joba = joba_i;

    n = n_i;
    m = m_i;
    
    mv = mv_i;
    ldv = n;
    if (matrix_layout == LAPACK_ROW_MAJOR)
    {
        lda = n;

    }   else
    {
        lda = m;
		if( jobv == 'A')
           ldv = mv;
    }

    // Initialize with default values
    hModule = NULL;
    dModule = NULL;
    diff_a = 0;
    diff_sva = 0;
    diff_v = 0;
    diff_stat = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n gesvj float: matrix_layout: %d   jobu: %c \t \
jobv: %c \t joba: %c \n  m: %d \t n: %d \t mv: %d \n",
matrix_layout, jobu_i, jobv_i, joba_i,  m_i, n_i, mv);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, m*n );
    lapacke_gtest_alloc_float_buffer_pair( &sva, &svaref, n );
    lapacke_gtest_alloc_float_buffer_pair( &v, &vref, n*n );

    if( (a==NULL) || (aref==NULL) ||  \
        (sva==NULL) || (svaref==NULL) || \
        (v==NULL) || (vref==NULL) ){
       EXPECT_FALSE( true) << "gesvj_float_parameters object: malloc error. Exiting ";
       gesvj_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    if (joba=='G')
    {
       lapacke_gtest_init_float_buffer_pair_rand( a, aref, m*n );
    } else {
       lapacke_gtest_init_float_buffer_pair_rand_custom_matrix( a, aref, m, n, joba);
    }
    lapacke_gtest_init_float_buffer_pair_with_constant( v, vref, n*n, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant( sva, svaref, n, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant( stat, statref, 6, 0.0);
    stat[0] = statref[0] = svd_paramslist[idx].ctol_gesvj;
   } /* end of Constructor  */

/* Destructor definition  'gesvj_float_common_parameters' */
gesvj_float_parameters :: ~gesvj_float_parameters ()
{
   /* De-Allocate memory for the input matrices */
   gesvj_free();
} 

//  Test fixture class definition
class sgesvj_test  : public  ::testing::Test {
public:
   gesvj_float_parameters  *sgesvj_obj;
   void SetUp();  
   void TearDown () { delete sgesvj_obj; }
};

void sgesvj_test::SetUp()
{
    /* LAPACKE sgesvj prototype */
    typedef int (*Fptr_NL_LAPACKE_sgesvj) ( int matrix_layout, char joba,
        char jobu, char jobv, lapack_int m, lapack_int n, float * a,
        lapack_int lda, float * sva, lapack_int mv, float * v,
        lapack_int ldv, float * stat);
            
    Fptr_NL_LAPACKE_sgesvj SGESVJ;
    sgesvj_obj = new  gesvj_float_parameters( svd_paramslist[idx].matrix_layout,
                                         svd_paramslist[idx].joba_gesvj,
                                         svd_paramslist[idx].jobu_gesvj,
                                         svd_paramslist[idx].jobv_gesvj,
                                         svd_paramslist[idx].m_gesvj,
                                         svd_paramslist[idx].n_gesvj,
                                         svd_paramslist[idx].mv_gesvj);

    idx = Circular_Increment_Index(idx);
    sgesvj_obj->threshold = svd_paramslist[idx].svd_threshold;
    sgesvj_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgesvj_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgesvj_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgesvj_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SGESVJ = (Fptr_NL_LAPACKE_sgesvj)dlsym(sgesvj_obj->hModule, "LAPACKE_sgesvj");
    ASSERT_TRUE(SGESVJ != NULL) << "failed to ppt the Netlib LAPACKE_sgesvj symbol";

    /* Compute libflame's Lapacke o/p  */
    sgesvj_obj->info = LAPACKE_sgesvj(  sgesvj_obj->matrix_layout,
                                    sgesvj_obj->joba,
                                    sgesvj_obj->jobu,
                                    sgesvj_obj->jobv,
                                    sgesvj_obj->m,
                                    sgesvj_obj->n,
                                    sgesvj_obj->a,
                                    sgesvj_obj->lda,
                                    sgesvj_obj->sva,
                                    sgesvj_obj->mv,
                                    sgesvj_obj->v,
                                    sgesvj_obj->ldv,
                                    sgesvj_obj->stat
                                    );

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    sgesvj_obj->inforef = SGESVJ( sgesvj_obj->matrix_layout,
                                    sgesvj_obj->joba,
                                    sgesvj_obj->jobu,
                                    sgesvj_obj->jobv,
                                    sgesvj_obj->m,
                                    sgesvj_obj->n,
                                    sgesvj_obj->aref,
                                    sgesvj_obj->lda,
                                    sgesvj_obj->svaref,
                                    sgesvj_obj->mv,
                                    sgesvj_obj->vref,
                                    sgesvj_obj->ldv,
                                    sgesvj_obj->statref
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    if( sgesvj_obj->jobu != 'N'){
       sgesvj_obj->diff_a =  computeDiff_s( (sgesvj_obj->m)*(sgesvj_obj->n), 
                sgesvj_obj->a, sgesvj_obj->aref );
    }
    sgesvj_obj->diff_v =  computeDiff_s( (sgesvj_obj->n)*(sgesvj_obj->n), 
                sgesvj_obj->v, sgesvj_obj->vref );

    sgesvj_obj->diff_sva =  computeDiff_s( sgesvj_obj->n, 
                sgesvj_obj->sva, sgesvj_obj->svaref );

    sgesvj_obj->diff_stat =  computeDiff_s( 6 , 
                sgesvj_obj->stat, sgesvj_obj->statref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n gesvj float: \n  diff_a: %f \n \
diff_v: %f \n diff_sva: %f \n diff_stat: %f info: %d \t inforef: %d\n ",
       sgesvj_obj->diff_a, sgesvj_obj->diff_v,
       sgesvj_obj->diff_sva, sgesvj_obj->diff_stat,
       sgesvj_obj->info, sgesvj_obj->inforef);
#endif
}

TEST_F(sgesvj_test, sgesvj_1) {
    EXPECT_NEAR(0.0, sgesvj_obj->diff_sva, sgesvj_obj->threshold);
    EXPECT_NEAR(0.0, sgesvj_obj->diff_a, sgesvj_obj->threshold);
    EXPECT_NEAR(0.0, sgesvj_obj->diff_v, sgesvj_obj->threshold);
    EXPECT_NEAR(0.0, sgesvj_obj->diff_stat, sgesvj_obj->threshold);
}

TEST_F(sgesvj_test, sgesvj_2) {
    EXPECT_NEAR(0.0, sgesvj_obj->diff_sva, sgesvj_obj->threshold);
    EXPECT_NEAR(0.0, sgesvj_obj->diff_a, sgesvj_obj->threshold);
    EXPECT_NEAR(0.0, sgesvj_obj->diff_v, sgesvj_obj->threshold);
    EXPECT_NEAR(0.0, sgesvj_obj->diff_stat, sgesvj_obj->threshold);
}

TEST_F(sgesvj_test, sgesvj_3) {
    EXPECT_NEAR(0.0, sgesvj_obj->diff_sva, sgesvj_obj->threshold);
    EXPECT_NEAR(0.0, sgesvj_obj->diff_a, sgesvj_obj->threshold);
    EXPECT_NEAR(0.0, sgesvj_obj->diff_v, sgesvj_obj->threshold);
    EXPECT_NEAR(0.0, sgesvj_obj->diff_stat, sgesvj_obj->threshold);
}

TEST_F(sgesvj_test, sgesvj_4) {
    EXPECT_NEAR(0.0, sgesvj_obj->diff_sva, sgesvj_obj->threshold);
    EXPECT_NEAR(0.0, sgesvj_obj->diff_a, sgesvj_obj->threshold);
    EXPECT_NEAR(0.0, sgesvj_obj->diff_v, sgesvj_obj->threshold);
    EXPECT_NEAR(0.0, sgesvj_obj->diff_stat, sgesvj_obj->threshold);
}

/* Begin gesvj_double_common_parameters  class definition */
class gesvj_double_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    double diff_a, diff_sva, diff_v, diff_stat;
    float threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char joba; //  Must be 'L', 'U' or 'G'.
    char jobu; // Must be 'U', 'C' or 'N'.
    char jobv; // Must be 'V', 'A' or 'N'.

    lapack_int m; // The number of rows of the matrix A
    lapack_int n; // The number of columns of the matrix A 
    lapack_int lda; //  The leading dimension of a

    lapack_int mv; // Jacobi rotations applied to the first mv rows of v.
    lapack_int ldv; // The leading dimension of the output array v . ldv≥ max(1, p)
    
    /* Input / Output parameters */
    double* a, *aref; // contains m-by-n matrix A.
    
    /* Output parameters */
    double* sva, *svaref; //   the singular values of A.
    double* v, *vref; //  contains the n-by-n matrix of the right singular vectors.
    double stat[6], statref[6]; //  scaling factor
    
    /*Return Values */
    int info, inforef;

   public:
      gesvj_double_parameters (int matrix_layout_i, char joba, char jobu, 
                    char jobv, lapack_int m, lapack_int n, lapack_int mv );
      ~gesvj_double_parameters ();
};

/* Constructor definition  gesvj double_common_parameters */
gesvj_double_parameters:: gesvj_double_parameters (int matrix_layout_i, char joba_i, 
                          char jobu_i, char jobv_i, lapack_int m_i, lapack_int n_i,
                          lapack_int mv_i)
{
    matrix_layout = matrix_layout_i;
    jobu = jobu_i;
    jobv = jobv_i;
    joba = joba_i;

    n = n_i;
    m = m_i;
    
    mv = mv_i;
    ldv = n;
    if (matrix_layout == LAPACK_ROW_MAJOR)
    {
        lda = n;

    }   else
    {
        lda = m;
		if( jobv == 'A')
           ldv = mv;
    }

    // Initialize with default values
    hModule = NULL;
    dModule = NULL;
    diff_a = 0;
    diff_sva = 0;
    diff_v = 0;
    diff_stat = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n gesvj double: matrix_layout: %d   jobu: %c \t \
jobv: %c \t joba: %c \n  m: %d \t n: %d \t mv: %d \n",
matrix_layout, jobu_i, jobv_i, joba_i,  m_i, n_i, mv);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, m*n );
    lapacke_gtest_alloc_double_buffer_pair( &sva, &svaref, n );
    lapacke_gtest_alloc_double_buffer_pair( &v, &vref, n*n );

    if( (a==NULL) || (aref==NULL) ||  \
        (sva==NULL) || (svaref==NULL) || \
        (v==NULL) || (vref==NULL) ){
       EXPECT_FALSE( true) << "gesvj_double_parameters object: malloc error. Exiting ";
       gesvj_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    if (joba=='G')
    {
       lapacke_gtest_init_double_buffer_pair_rand( a, aref, m*n );
    } else {
       lapacke_gtest_init_double_buffer_pair_rand_custom_matrix( a, aref, m, n, joba);
    }
    lapacke_gtest_init_double_buffer_pair_with_constant( v, vref, n*n, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant( sva, svaref, n, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant( stat, statref, 6, 0.0);
    stat[0] = statref[0] = svd_paramslist[idx].ctol_gesvj;
   } /* end of Constructor  */

/* Destructor definition  'gesvj_double_common_parameters' */
gesvj_double_parameters :: ~gesvj_double_parameters ()
{
   /* De-Allocate memory for the input matrices */
   gesvj_free();
} 

//  Test fixture class definition
class dgesvj_test  : public  ::testing::Test {
public:
   gesvj_double_parameters  *dgesvj_obj;
   void SetUp();  
   void TearDown () { delete dgesvj_obj; }
};

void dgesvj_test::SetUp()
{
    /* LAPACKE dgesvj prototype */
    typedef int (*Fptr_NL_LAPACKE_dgesvj) ( int matrix_layout, char joba,
        char jobu, char jobv, lapack_int m, lapack_int n, double * a,
        lapack_int lda, double * sva, lapack_int mv, double * v,
        lapack_int ldv, double * stat);
            
    Fptr_NL_LAPACKE_dgesvj DGESVJ;
    dgesvj_obj = new  gesvj_double_parameters( svd_paramslist[idx].matrix_layout,
                                         svd_paramslist[idx].joba_gesvj,
                                         svd_paramslist[idx].jobu_gesvj,
                                         svd_paramslist[idx].jobv_gesvj,
                                         svd_paramslist[idx].m_gesvj,
                                         svd_paramslist[idx].n_gesvj,
                                         svd_paramslist[idx].mv_gesvj);

    idx = Circular_Increment_Index(idx);
    dgesvj_obj->threshold = svd_paramslist[idx].svd_threshold;
    dgesvj_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgesvj_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgesvj_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgesvj_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DGESVJ = (Fptr_NL_LAPACKE_dgesvj)dlsym(dgesvj_obj->hModule, "LAPACKE_dgesvj");
    ASSERT_TRUE(DGESVJ != NULL) << "failed to ppt the Netlib LAPACKE_dgesvj symbol";

    /* Compute libflame's Lapacke o/p  */
    dgesvj_obj->info = LAPACKE_dgesvj(  dgesvj_obj->matrix_layout,
                                    dgesvj_obj->joba,
                                    dgesvj_obj->jobu,
                                    dgesvj_obj->jobv,
                                    dgesvj_obj->m,
                                    dgesvj_obj->n,
                                    dgesvj_obj->a,
                                    dgesvj_obj->lda,
                                    dgesvj_obj->sva,
                                    dgesvj_obj->mv,
                                    dgesvj_obj->v,
                                    dgesvj_obj->ldv,
                                    dgesvj_obj->stat
                                    );

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    dgesvj_obj->inforef = DGESVJ( dgesvj_obj->matrix_layout,
                                    dgesvj_obj->joba,
                                    dgesvj_obj->jobu,
                                    dgesvj_obj->jobv,
                                    dgesvj_obj->m,
                                    dgesvj_obj->n,
                                    dgesvj_obj->aref,
                                    dgesvj_obj->lda,
                                    dgesvj_obj->svaref,
                                    dgesvj_obj->mv,
                                    dgesvj_obj->vref,
                                    dgesvj_obj->ldv,
                                    dgesvj_obj->statref
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    if( dgesvj_obj->jobu != 'N'){
       dgesvj_obj->diff_a =  computeDiff_d( (dgesvj_obj->m)*(dgesvj_obj->n), 
                dgesvj_obj->a, dgesvj_obj->aref );
    }
    dgesvj_obj->diff_v =  computeDiff_d( (dgesvj_obj->n)*(dgesvj_obj->n), 
                dgesvj_obj->v, dgesvj_obj->vref );

    dgesvj_obj->diff_sva =  computeDiff_d( dgesvj_obj->n, 
                dgesvj_obj->sva, dgesvj_obj->svaref );

    dgesvj_obj->diff_stat =  computeDiff_d( 6 , 
                dgesvj_obj->stat, dgesvj_obj->statref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n gesvj double: \n  diff_a: %f \n \
diff_v: %f \n diff_sva: %f \n diff_stat: %f info: %d \t inforef: %d\n ",
       dgesvj_obj->diff_a, dgesvj_obj->diff_v,
       dgesvj_obj->diff_sva, dgesvj_obj->diff_stat,
       dgesvj_obj->info, dgesvj_obj->inforef);
#endif
}

TEST_F(dgesvj_test, dgesvj_1) {
    EXPECT_NEAR(0.0, dgesvj_obj->diff_sva, dgesvj_obj->threshold);
    EXPECT_NEAR(0.0, dgesvj_obj->diff_a, dgesvj_obj->threshold);
    EXPECT_NEAR(0.0, dgesvj_obj->diff_v, dgesvj_obj->threshold);
    EXPECT_NEAR(0.0, dgesvj_obj->diff_stat, dgesvj_obj->threshold);
}

TEST_F(dgesvj_test, dgesvj_2) {
    EXPECT_NEAR(0.0, dgesvj_obj->diff_sva, dgesvj_obj->threshold);
    EXPECT_NEAR(0.0, dgesvj_obj->diff_a, dgesvj_obj->threshold);
    EXPECT_NEAR(0.0, dgesvj_obj->diff_v, dgesvj_obj->threshold);
    EXPECT_NEAR(0.0, dgesvj_obj->diff_stat, dgesvj_obj->threshold);
}

TEST_F(dgesvj_test, dgesvj_3) {
    EXPECT_NEAR(0.0, dgesvj_obj->diff_sva, dgesvj_obj->threshold);
    EXPECT_NEAR(0.0, dgesvj_obj->diff_a, dgesvj_obj->threshold);
    EXPECT_NEAR(0.0, dgesvj_obj->diff_v, dgesvj_obj->threshold);
    EXPECT_NEAR(0.0, dgesvj_obj->diff_stat, dgesvj_obj->threshold);
}

TEST_F(dgesvj_test, dgesvj_4) {
    EXPECT_NEAR(0.0, dgesvj_obj->diff_sva, dgesvj_obj->threshold);
    EXPECT_NEAR(0.0, dgesvj_obj->diff_a, dgesvj_obj->threshold);
    EXPECT_NEAR(0.0, dgesvj_obj->diff_v, dgesvj_obj->threshold);
    EXPECT_NEAR(0.0, dgesvj_obj->diff_stat, dgesvj_obj->threshold);
}


/* Begin gesvj_scomplex_common_parameters  class definition */
class gesvj_scomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_a, diff_sva, diff_v, diff_stat;
    float threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char joba; //  Must be 'L', 'U' or 'G'.
    char jobu; // Must be 'U', 'C' or 'N'.
    char jobv; // Must be 'V', 'A' or 'N'.

    lapack_int m; // The number of rows of the matrix A
    lapack_int n; // The number of columns of the matrix A 
    lapack_int lda; //  The leading dimension of a

    lapack_int mv; // Jacobi rotations applied to the first mv rows of v.
    lapack_int ldv; // The leading dimension of the output array v . ldv≥ max(1, p)
    
    /* Input / Output parameters */
    lapack_complex_float* a, *aref; // contains m-by-n matrix A.
    
    /* Output parameters */
    float* sva, *svaref; //   the singular values of A.
    lapack_complex_float* v, *vref; //  contains the n-by-n matrix of the right singular vectors.
    float stat[6], statref[6]; //  scaling factor
    
    /*Return Values */
    int info, inforef;

   public:
      gesvj_scomplex_parameters (int matrix_layout_i, char joba, char jobu, 
                    char jobv, lapack_int m, lapack_int n, lapack_int mv );
      ~gesvj_scomplex_parameters ();
};

/* Constructor definition  gesvj lapack_complex_float_common_parameters */
gesvj_scomplex_parameters:: gesvj_scomplex_parameters (int matrix_layout_i, char joba_i, 
                          char jobu_i, char jobv_i, lapack_int m_i, lapack_int n_i,
                          lapack_int mv_i)
{
    matrix_layout = matrix_layout_i;
    jobu = jobu_i;
    jobv = jobv_i;
    joba = joba_i;

    n = n_i;
    m = m_i;
    
    mv = mv_i;
    ldv = n;
    if (matrix_layout == LAPACK_ROW_MAJOR)
    {
        lda = n;

    }   else
    {
        lda = m;
		if( jobv == 'A')
           ldv = mv;
    }

    // Initialize with default values
    hModule = NULL;
    dModule = NULL;
    diff_a = 0;
    diff_sva = 0;
    diff_v = 0;
    diff_stat = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n gesvj lapack_complex_float: matrix_layout: %d   jobu: %c \t \
jobv: %c \t joba: %c \n  m: %d \t n: %d \t mv: %d \n",
matrix_layout, jobu_i, jobv_i, joba_i,  m_i, n_i, mv);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, m*n );
    lapacke_gtest_alloc_float_buffer_pair( &sva, &svaref, n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &v, &vref, n*n );

    if( (a==NULL) || (aref==NULL) ||  \
        (sva==NULL) || (svaref==NULL) || \
        (v==NULL) || (vref==NULL) ){
       EXPECT_FALSE( true) << "gesvj_scomplex_parameters object: malloc error. Exiting ";
       gesvj_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    if (joba=='G')
    {
       lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, m*n );
    } else {
       lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix( a, aref, m, n, joba);
    }
	if (jobv=='V') {
		lapacke_gtest_init_scomplex_buffer_pair_with_constant( v, vref, n*n, 0.0);
	} else {
		lapacke_gtest_init_scomplex_buffer_pair_with_constant( v, vref, mv*n, 0.0);
	}
    lapacke_gtest_init_float_buffer_pair_with_constant( sva, svaref, n, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant( stat, statref, 6, 0.0);
    stat[0] = statref[0] = svd_paramslist[idx].ctol_gesvj;
   } /* end of Constructor  */

/* Destructor definition  'gesvj_scomplex_common_parameters' */
gesvj_scomplex_parameters :: ~gesvj_scomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   gesvj_free();
} 

//  Test fixture class definition
class cgesvj_test  : public  ::testing::Test {
public:
   gesvj_scomplex_parameters  *cgesvj_obj;
   void SetUp();  
   void TearDown () { delete cgesvj_obj; }
};

void cgesvj_test::SetUp()
{
    /* LAPACKE cgesvj prototype */
    typedef int (*Fptr_NL_LAPACKE_cgesvj) ( int matrix_layout, char joba,
            char jobu, char jobv, lapack_int m, lapack_int n,
            lapack_complex_float * a, lapack_int lda, float * sva,
            lapack_int mv, lapack_complex_float * v,
            lapack_int ldv, float * stat);
            
    Fptr_NL_LAPACKE_cgesvj CGESVJ;
    cgesvj_obj = new  gesvj_scomplex_parameters( svd_paramslist[idx].matrix_layout,
                                         svd_paramslist[idx].joba_gesvj,
                                         svd_paramslist[idx].jobu_gesvj,
                                         svd_paramslist[idx].jobv_gesvj,
                                         svd_paramslist[idx].m_gesvj,
                                         svd_paramslist[idx].n_gesvj,
                                         svd_paramslist[idx].mv_gesvj);

    idx = Circular_Increment_Index(idx);
    cgesvj_obj->threshold = svd_paramslist[idx].svd_threshold;
    cgesvj_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgesvj_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgesvj_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgesvj_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CGESVJ = (Fptr_NL_LAPACKE_cgesvj)dlsym(cgesvj_obj->hModule, "LAPACKE_cgesvj");
    ASSERT_TRUE(CGESVJ != NULL) << "failed to ppt the Netlib LAPACKE_cgesvj symbol";

    /* Compute libflame's Lapacke o/p  */
    cgesvj_obj->info = LAPACKE_cgesvj(  cgesvj_obj->matrix_layout,
                                    cgesvj_obj->joba,
                                    cgesvj_obj->jobu,
                                    cgesvj_obj->jobv,
                                    cgesvj_obj->m,
                                    cgesvj_obj->n,
                                    cgesvj_obj->a,
                                    cgesvj_obj->lda,
                                    cgesvj_obj->sva,
                                    cgesvj_obj->mv,
                                    cgesvj_obj->v,
                                    cgesvj_obj->ldv,
                                    cgesvj_obj->stat
                                    );

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    cgesvj_obj->inforef = CGESVJ( cgesvj_obj->matrix_layout,
                                    cgesvj_obj->joba,
                                    cgesvj_obj->jobu,
                                    cgesvj_obj->jobv,
                                    cgesvj_obj->m,
                                    cgesvj_obj->n,
                                    cgesvj_obj->aref,
                                    cgesvj_obj->lda,
                                    cgesvj_obj->svaref,
                                    cgesvj_obj->mv,
                                    cgesvj_obj->vref,
                                    cgesvj_obj->ldv,
                                    cgesvj_obj->statref
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    if( cgesvj_obj->jobu != 'N'){
       cgesvj_obj->diff_a =  computeDiff_c( (cgesvj_obj->m)*(cgesvj_obj->n), 
                cgesvj_obj->a, cgesvj_obj->aref );
    }
	if (cgesvj_obj->jobv=='V') {
		cgesvj_obj->diff_v =  computeDiff_c( (cgesvj_obj->n)*(cgesvj_obj->n), 
                cgesvj_obj->v, cgesvj_obj->vref );
	} else {
		cgesvj_obj->diff_v =  computeDiff_c( (cgesvj_obj->mv)*(cgesvj_obj->n), 
                cgesvj_obj->v, cgesvj_obj->vref );
	}

    cgesvj_obj->diff_sva =  computeDiff_s( cgesvj_obj->n, 
                cgesvj_obj->sva, cgesvj_obj->svaref );

    cgesvj_obj->diff_stat =  computeDiff_s( 6 , 
                cgesvj_obj->stat, cgesvj_obj->statref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n gesvj lapack_complex_float: \n  diff_a: %f \n \
diff_v: %f \n diff_sva: %f \n diff_stat: %f info: %d \t inforef: %d\n ",
       cgesvj_obj->diff_a, cgesvj_obj->diff_v,
       cgesvj_obj->diff_sva, cgesvj_obj->diff_stat,
       cgesvj_obj->info, cgesvj_obj->inforef);
#endif
}

TEST_F(cgesvj_test, cgesvj_1) {
    EXPECT_NEAR(0.0, cgesvj_obj->diff_sva, cgesvj_obj->threshold);
    EXPECT_NEAR(0.0, cgesvj_obj->diff_a, cgesvj_obj->threshold);
    EXPECT_NEAR(0.0, cgesvj_obj->diff_v, cgesvj_obj->threshold);
    EXPECT_NEAR(0.0, cgesvj_obj->diff_stat, cgesvj_obj->threshold);
}

TEST_F(cgesvj_test, cgesvj_2) {
    EXPECT_NEAR(0.0, cgesvj_obj->diff_sva, cgesvj_obj->threshold);
    EXPECT_NEAR(0.0, cgesvj_obj->diff_a, cgesvj_obj->threshold);
    EXPECT_NEAR(0.0, cgesvj_obj->diff_v, cgesvj_obj->threshold);
    EXPECT_NEAR(0.0, cgesvj_obj->diff_stat, cgesvj_obj->threshold);
}

TEST_F(cgesvj_test, cgesvj_3) {
    EXPECT_NEAR(0.0, cgesvj_obj->diff_sva, cgesvj_obj->threshold);
    EXPECT_NEAR(0.0, cgesvj_obj->diff_a, cgesvj_obj->threshold);
    EXPECT_NEAR(0.0, cgesvj_obj->diff_v, cgesvj_obj->threshold);
    EXPECT_NEAR(0.0, cgesvj_obj->diff_stat, cgesvj_obj->threshold);
}

TEST_F(cgesvj_test, cgesvj_4) {
    EXPECT_NEAR(0.0, cgesvj_obj->diff_sva, cgesvj_obj->threshold);
    EXPECT_NEAR(0.0, cgesvj_obj->diff_a, cgesvj_obj->threshold);
    EXPECT_NEAR(0.0, cgesvj_obj->diff_v, cgesvj_obj->threshold);
    EXPECT_NEAR(0.0, cgesvj_obj->diff_stat, cgesvj_obj->threshold);
}

/* Begin gesvj_dcomplex_common_parameters  class definition */
class gesvj_dcomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    double diff_a, diff_sva, diff_v, diff_stat;
    float threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char joba; //  Must be 'L', 'U' or 'G'.
    char jobu; // Must be 'U', 'C' or 'N'.
    char jobv; // Must be 'V', 'A' or 'N'.

    lapack_int m; // The number of rows of the matrix A
    lapack_int n; // The number of columns of the matrix A 
    lapack_int lda; //  The leading dimension of a

    lapack_int mv; // Jacobi rotations applied to the first mv rows of v.
    lapack_int ldv; // The leading dimension of the output array v . ldv≥ max(1, p)
    
    /* Input / Output parameters */
    lapack_complex_double* a, *aref; // contains m-by-n matrix A.
    
    /* Output parameters */
    double* sva, *svaref; //   the singular values of A.
    lapack_complex_double* v, *vref; //  contains the n-by-n matrix of the right singular vectors.
    double stat[6], statref[6]; //  scaling factor
    
    /*Return Values */
    int info, inforef;

   public:
      gesvj_dcomplex_parameters (int matrix_layout_i, char joba, char jobu, 
                    char jobv, lapack_int m, lapack_int n, lapack_int mv );
      ~gesvj_dcomplex_parameters ();
};

/* Constructor definition  gesvj lapack_complex_double_common_parameters */
gesvj_dcomplex_parameters:: gesvj_dcomplex_parameters (int matrix_layout_i, char joba_i, 
                          char jobu_i, char jobv_i, lapack_int m_i, lapack_int n_i,
                          lapack_int mv_i)
{
    matrix_layout = matrix_layout_i;
    jobu = jobu_i;
    jobv = jobv_i;
    joba = joba_i;

    n = n_i;
    m = m_i;
    
    mv = mv_i;
    ldv = n;
    if (matrix_layout == LAPACK_ROW_MAJOR)
    {
        lda = n;

    }   else
    {
        lda = m;
		if( jobv == 'A')
           ldv = mv;
    }

    // Initialize with default values
    hModule = NULL;
    dModule = NULL;
    diff_a = 0;
    diff_sva = 0;
    diff_v = 0;
    diff_stat = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n gesvj lapack_complex_double: matrix_layout: %d   jobu: %c \t \
jobv: %c \t joba: %c \n  m: %d \t n: %d \t mv: %d \n",
matrix_layout, jobu_i, jobv_i, joba_i,  m_i, n_i, mv);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, m*n );
    lapacke_gtest_alloc_double_buffer_pair( &sva, &svaref, n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &v, &vref, n*n );

    if( (a==NULL) || (aref==NULL) ||  \
        (sva==NULL) || (svaref==NULL) || \
        (v==NULL) || (vref==NULL) ){
       EXPECT_FALSE( true) << "gesvj_dcomplex_parameters object: malloc error. Exiting ";
       gesvj_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    if (joba=='G')
    {
       lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, m*n );
    } else {
       lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix( a, aref, m, n, joba);
    }
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( v, vref, n*n, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant( sva, svaref, n, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant( stat, statref, 6, 0.0);
    stat[0] = statref[0] = svd_paramslist[idx].ctol_gesvj;
   } /* end of Constructor  */

/* Destructor definition  'gesvj_dcomplex_common_parameters' */
gesvj_dcomplex_parameters :: ~gesvj_dcomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   gesvj_free();
} 

//  Test fixture class definition
class zgesvj_test  : public  ::testing::Test {
public:
   gesvj_dcomplex_parameters  *zgesvj_obj;
   void SetUp();  
   void TearDown () { delete zgesvj_obj; }
};

void zgesvj_test::SetUp()
{
    /* LAPACKE zgesvj prototype */
    typedef int (*Fptr_NL_LAPACKE_zgesvj) ( int matrix_layout, char joba,
            char jobu, char jobv, lapack_int m, lapack_int n,
            lapack_complex_double * a, lapack_int lda, double * sva,
            lapack_int mv, lapack_complex_double * v,
            lapack_int ldv, double * stat);
            
    Fptr_NL_LAPACKE_zgesvj ZGESVJ;
    zgesvj_obj = new  gesvj_dcomplex_parameters( svd_paramslist[idx].matrix_layout,
                                         svd_paramslist[idx].joba_gesvj,
                                         svd_paramslist[idx].jobu_gesvj,
                                         svd_paramslist[idx].jobv_gesvj,
                                         svd_paramslist[idx].m_gesvj,
                                         svd_paramslist[idx].n_gesvj,
                                         svd_paramslist[idx].mv_gesvj);

    idx = Circular_Increment_Index(idx);
    zgesvj_obj->threshold = svd_paramslist[idx].svd_threshold;
    zgesvj_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgesvj_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgesvj_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgesvj_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZGESVJ = (Fptr_NL_LAPACKE_zgesvj)dlsym(zgesvj_obj->hModule, "LAPACKE_zgesvj");
    ASSERT_TRUE(ZGESVJ != NULL) << "failed to ppt the Netlib LAPACKE_zgesvj symbol";

    /* Compute libflame's Lapacke o/p  */
    zgesvj_obj->info = LAPACKE_zgesvj(  zgesvj_obj->matrix_layout,
                                    zgesvj_obj->joba,
                                    zgesvj_obj->jobu,
                                    zgesvj_obj->jobv,
                                    zgesvj_obj->m,
                                    zgesvj_obj->n,
                                    zgesvj_obj->a,
                                    zgesvj_obj->lda,
                                    zgesvj_obj->sva,
                                    zgesvj_obj->mv,
                                    zgesvj_obj->v,
                                    zgesvj_obj->ldv,
                                    zgesvj_obj->stat
                                    );

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    zgesvj_obj->inforef = ZGESVJ( zgesvj_obj->matrix_layout,
                                    zgesvj_obj->joba,
                                    zgesvj_obj->jobu,
                                    zgesvj_obj->jobv,
                                    zgesvj_obj->m,
                                    zgesvj_obj->n,
                                    zgesvj_obj->aref,
                                    zgesvj_obj->lda,
                                    zgesvj_obj->svaref,
                                    zgesvj_obj->mv,
                                    zgesvj_obj->vref,
                                    zgesvj_obj->ldv,
                                    zgesvj_obj->statref
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    if( zgesvj_obj->jobu != 'N'){
       zgesvj_obj->diff_a =  computeDiff_z( (zgesvj_obj->m)*(zgesvj_obj->n), 
                zgesvj_obj->a, zgesvj_obj->aref );
    }
	
	if (zgesvj_obj->jobv=='V') {
		zgesvj_obj->diff_v =  computeDiff_z( (zgesvj_obj->n)*(zgesvj_obj->n), 
                zgesvj_obj->v, zgesvj_obj->vref );
	} else {
		zgesvj_obj->diff_v =  computeDiff_z( (zgesvj_obj->mv)*(zgesvj_obj->n), 
                zgesvj_obj->v, zgesvj_obj->vref );
    }
    zgesvj_obj->diff_sva =  computeDiff_d( zgesvj_obj->n, 
                zgesvj_obj->sva, zgesvj_obj->svaref );

    zgesvj_obj->diff_stat =  computeDiff_d( 6 , 
                zgesvj_obj->stat, zgesvj_obj->statref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n gesvj lapack_complex_double: \n  diff_a: %f \n \
diff_v: %f \n diff_sva: %f \n diff_stat: %f info: %d \t inforef: %d\n ",
       zgesvj_obj->diff_a, zgesvj_obj->diff_v,
       zgesvj_obj->diff_sva, zgesvj_obj->diff_stat,
       zgesvj_obj->info, zgesvj_obj->inforef);
#endif
}

TEST_F(zgesvj_test, zgesvj_1) {
    EXPECT_NEAR(0.0, zgesvj_obj->diff_sva, zgesvj_obj->threshold);
    EXPECT_NEAR(0.0, zgesvj_obj->diff_a, zgesvj_obj->threshold);
    EXPECT_NEAR(0.0, zgesvj_obj->diff_v, zgesvj_obj->threshold);
    EXPECT_NEAR(0.0, zgesvj_obj->diff_stat, zgesvj_obj->threshold);
}

TEST_F(zgesvj_test, zgesvj_2) {
    EXPECT_NEAR(0.0, zgesvj_obj->diff_sva, zgesvj_obj->threshold);
    EXPECT_NEAR(0.0, zgesvj_obj->diff_a, zgesvj_obj->threshold);
    EXPECT_NEAR(0.0, zgesvj_obj->diff_v, zgesvj_obj->threshold);
    EXPECT_NEAR(0.0, zgesvj_obj->diff_stat, zgesvj_obj->threshold);
}

TEST_F(zgesvj_test, zgesvj_3) {
    EXPECT_NEAR(0.0, zgesvj_obj->diff_sva, zgesvj_obj->threshold);
    EXPECT_NEAR(0.0, zgesvj_obj->diff_a, zgesvj_obj->threshold);
    EXPECT_NEAR(0.0, zgesvj_obj->diff_v, zgesvj_obj->threshold);
    EXPECT_NEAR(0.0, zgesvj_obj->diff_stat, zgesvj_obj->threshold);
}

TEST_F(zgesvj_test, zgesvj_4) {
    EXPECT_NEAR(0.0, zgesvj_obj->diff_sva, zgesvj_obj->threshold);
    EXPECT_NEAR(0.0, zgesvj_obj->diff_a, zgesvj_obj->threshold);
    EXPECT_NEAR(0.0, zgesvj_obj->diff_v, zgesvj_obj->threshold);
    EXPECT_NEAR(0.0, zgesvj_obj->diff_stat, zgesvj_obj->threshold);
}
