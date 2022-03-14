#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"
#define LAPACKE_TEST_VERBOSE  (1)
#define gejsv_free() \
    if (a!=NULL)        free(a); \
    if (aref!=NULL)     free(aref); \
    if (all!=NULL)        free(all); \
    if (allref!=NULL)     free(allref); \
    if (sva!=NULL)        free(sva); \
    if (svaref!=NULL)     free(svaref); \
    if (u!=NULL)        free(u); \
    if (uref!=NULL)     free(uref); \
    if (v!=NULL)        free(v); \
    if (vref!=NULL)     free(vref); \
    if (istat!=NULL)        free(istat); \
    if (istatref!=NULL)        free(istatref); \
    if (stat!=NULL)        free(stat); \
    if (statref!=NULL)     free(statref)


// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin gejsv_float_common_parameters  class definition */
class gejsv_float_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_a, diff_sva, diff_u, diff_v, diff_stat;
    int min_mn;
    float threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
	char joba; //  Must be 'C', 'E', 'F', 'G', 'A', or 'R'.
	char jobu; // Must be 'U', 'F', 'W', or 'N'.
	char jobv; // Must be 'V', 'J', 'W', or 'N'.
	char jobr; // Must be 'N' or 'R'.
	char jobt; // Must be 'T' or 'N'.
	char jobp; //  Must be 'P' or 'N'.

    lapack_int m; // The number of rows of the matrix A
    lapack_int n; // The number of columns of the matrix A 
    lapack_int lda; //  The leading dimension of a

    lapack_int ldu; // leading dimension of the output array u . ldu≥ max(1, m)
    lapack_int ldv; // The leading dimension of the output array v . ldv≥ max(1, p)
	
    /* Input / Output parameters */
    float* a, *aref; // contains m-by-n matrix A.
    float* u, *uref; //  contains the m-by-n matrix of the left singular vectors.
    float* v, *vref; //  contains the n-by-n matrix of the right singular vectors.
    float* all, *allref; // contains work array.
	
    /* Output parameters */
    float* sva, *svaref; //   the singular values of A.
    float* stat, *statref; //  scaling factor
    lapack_int * istat, *istatref; // the numerical rank determined after the initial QR factorization
    
    /*Return Values */
    int info, inforef;

   public:
      gejsv_float_parameters (int matrix_layout_i, char joba, char jobu, 
							  char jobv,char jobr, char jobt, char jobp,
                                             lapack_int m, lapack_int n);
      ~gejsv_float_parameters ();
};

/* Constructor definition  gejsv float_common_parameters */
gejsv_float_parameters:: gejsv_float_parameters (int matrix_layout_i, char joba_i, 
							  char jobu_i, char jobv_i,char jobr_i, char jobt_i,
							  char jobp_i, lapack_int m_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    jobu = jobu_i;
    jobv = jobv_i;
    joba = joba_i;
    jobp = jobp_i;
    jobt = jobt_i;
    jobr = jobr_i;
    n = n_i;
    m = m_i;
    min_mn = (m<n)?m:n;
    
    ldu = m;
    ldv = n;

    if (matrix_layout == LAPACK_ROW_MAJOR)
    {
        lda = n;
        if( jobu== 'F')  ldu = m; 

    }   else
    {
        lda = m;
    }

    // Initialize with default values
    hModule = NULL;
    dModule = NULL;
    diff_a = 0;
    diff_sva = 0;
    diff_u = 0;
    diff_v = 0;
    diff_stat = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n gejsv float: matrix_layout: %d   jobu: %c \t \
jobv: %c \t joba: %c \n jobt: %c \t jobp: %c \t jobr: %c \t m: %d \t n: %d \n",
matrix_layout, jobu_i, jobv_i, joba_i, jobt_i, jobp_i, jobr_i, m_i, n_i);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, m*n );
    lapacke_gtest_alloc_float_buffer_pair( &all, &allref, n );
    lapacke_gtest_alloc_float_buffer_pair( &sva, &svaref, n );
    lapacke_gtest_alloc_float_buffer_pair( &u, &uref, m*m ); // allocating largest size
    lapacke_gtest_alloc_float_buffer_pair( &v, &vref, ldv*n );
    lapacke_gtest_alloc_float_buffer_pair( &stat, &statref, n );
    lapacke_gtest_alloc_int_buffer_pair( &istat, &istatref, min_mn );

    if( (a==NULL) || (aref==NULL) ||  \
        (sva==NULL) || (svaref==NULL) || \
        (u==NULL) || (uref==NULL) || \
        (all==NULL) || (allref==NULL) || \
        (v==NULL) || (vref==NULL) || \
        (istat==NULL) || (istatref==NULL) || \
        (stat==NULL) || (statref==NULL)  ){
       EXPECT_FALSE( true) << "gejsv_float_parameters object: malloc error. Exiting ";
       gejsv_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, m*n );
    lapacke_gtest_init_float_buffer_pair_with_constant( u, uref, m*ldu, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant( v, vref, ldv*n, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant( sva, svaref, n, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant( stat, statref, n, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant( all, allref, n, 0.0);
    lapacke_gtest_init_int_buffer_pair_with_constant( istat, istatref, n, 0);

   } /* end of Constructor  */

/* Destructor definition  'gejsv_float_common_parameters' */
gejsv_float_parameters :: ~gejsv_float_parameters ()
{
   /* De-Allocate memory for the input matrices */
   gejsv_free();
} 

//  Test fixture class definition
class sgejsv_test  : public  ::testing::Test {
public:
   gejsv_float_parameters  *sgejsv_obj;
   void SetUp();  
   void TearDown () { delete sgejsv_obj; }
};

void sgejsv_test::SetUp()
{
    /* LAPACKE sgejsv prototype */
    typedef int (*Fptr_NL_LAPACKE_sgejsv) ( int matrix_layout, char joba,
			char jobu, char jobv, char jobr, char jobt, char jobp,
			lapack_int m, lapack_int n, float * a, lapack_int lda,
			float * all, float * u, lapack_int ldu, float * v,
			lapack_int ldv, float * stat, lapack_int * state);
            
    Fptr_NL_LAPACKE_sgejsv SGEJSV;
    sgejsv_obj = new  gejsv_float_parameters( svd_paramslist[idx].matrix_layout,
                                         svd_paramslist[idx].joba_gejsv,
                                         svd_paramslist[idx].jobu_gejsv,
                                         svd_paramslist[idx].jobv_gejsv,
                                         svd_paramslist[idx].jobr_gejsv,
                                         svd_paramslist[idx].jobt_gejsv,
                                         svd_paramslist[idx].jobp_gejsv,
                                         svd_paramslist[idx].m_gejsv,
                                         svd_paramslist[idx].n_gejsv );

    sgejsv_obj->threshold = svd_paramslist[idx].svd_threshold;
    idx = Circular_Increment_Index(idx);
    sgejsv_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgejsv_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgejsv_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgejsv_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SGEJSV = (Fptr_NL_LAPACKE_sgejsv)dlsym(sgejsv_obj->hModule, "LAPACKE_sgejsv");
    ASSERT_TRUE(SGEJSV != NULL) << "failed to ppt the Netlib LAPACKE_sgejsv symbol";

    /* Compute libflame's Lapacke o/p  */
    sgejsv_obj->info = LAPACKE_sgejsv(  sgejsv_obj->matrix_layout,
                                    sgejsv_obj->joba,
                                    sgejsv_obj->jobu,
                                    sgejsv_obj->jobv,
                                    sgejsv_obj->jobr,
                                    sgejsv_obj->jobt,
                                    sgejsv_obj->jobp,
                                    sgejsv_obj->m,
                                    sgejsv_obj->n,
                                    sgejsv_obj->a,
                                    sgejsv_obj->lda,
                                    sgejsv_obj->all,
                                    sgejsv_obj->u,
                                    sgejsv_obj->ldu,
                                    sgejsv_obj->v,
                                    sgejsv_obj->ldv,
                                    sgejsv_obj->stat,
									sgejsv_obj->istat
                                    );

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    sgejsv_obj->inforef = SGEJSV( sgejsv_obj->matrix_layout,
                                    sgejsv_obj->joba,
                                    sgejsv_obj->jobu,
                                    sgejsv_obj->jobv,
                                    sgejsv_obj->jobr,
                                    sgejsv_obj->jobt,
                                    sgejsv_obj->jobp,
                                    sgejsv_obj->m,
                                    sgejsv_obj->n,
                                    sgejsv_obj->aref,
                                    sgejsv_obj->lda,
                                    sgejsv_obj->allref,
                                    sgejsv_obj->uref,
                                    sgejsv_obj->ldu,
                                    sgejsv_obj->vref,
                                    sgejsv_obj->ldv,
                                    sgejsv_obj->statref,
									sgejsv_obj->istatref
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    sgejsv_obj->diff_u =  computeDiff_s( (sgejsv_obj->m)*(sgejsv_obj->m), 
                sgejsv_obj->u, sgejsv_obj->uref );

    sgejsv_obj->diff_v =  computeDiff_s( (sgejsv_obj->ldv)*(sgejsv_obj->n), 
                sgejsv_obj->v, sgejsv_obj->vref );

    sgejsv_obj->diff_sva =  computeDiff_s( sgejsv_obj->n, 
                sgejsv_obj->all, sgejsv_obj->allref );

    sgejsv_obj->diff_stat =  computeDiff_s( sgejsv_obj->min_mn, 
                sgejsv_obj->stat, sgejsv_obj->statref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n gejsv float: \n  diff_sva: %f \n \
diff_u: %f \n diff_v: %f \n diff_stat: %f info: %d \t inforef: %d\n ",
       sgejsv_obj->diff_sva, sgejsv_obj->diff_u,
       sgejsv_obj->diff_v, sgejsv_obj->diff_stat,
	   sgejsv_obj->info, sgejsv_obj->inforef);
#endif
}

TEST_F(sgejsv_test, sgejsv_1) {
    EXPECT_NEAR(0.0, sgejsv_obj->diff_sva, sgejsv_obj->threshold);
    EXPECT_NEAR(0.0, sgejsv_obj->diff_u, sgejsv_obj->threshold);
    EXPECT_NEAR(0.0, sgejsv_obj->diff_v, sgejsv_obj->threshold);
    EXPECT_NEAR(0.0, sgejsv_obj->diff_stat, sgejsv_obj->threshold);
}

TEST_F(sgejsv_test, sgejsv_2) {
    EXPECT_NEAR(0.0, sgejsv_obj->diff_sva, sgejsv_obj->threshold);
    EXPECT_NEAR(0.0, sgejsv_obj->diff_u, sgejsv_obj->threshold);
    EXPECT_NEAR(0.0, sgejsv_obj->diff_v, sgejsv_obj->threshold);
    EXPECT_NEAR(0.0, sgejsv_obj->diff_stat, sgejsv_obj->threshold);}

TEST_F(sgejsv_test, sgejsv_3) {
    EXPECT_NEAR(0.0, sgejsv_obj->diff_sva, sgejsv_obj->threshold);
    EXPECT_NEAR(0.0, sgejsv_obj->diff_u, sgejsv_obj->threshold);
    EXPECT_NEAR(0.0, sgejsv_obj->diff_v, sgejsv_obj->threshold);
    EXPECT_NEAR(0.0, sgejsv_obj->diff_stat, sgejsv_obj->threshold);}

TEST_F(sgejsv_test, sgejsv_4) {
    EXPECT_NEAR(0.0, sgejsv_obj->diff_sva, sgejsv_obj->threshold);
    EXPECT_NEAR(0.0, sgejsv_obj->diff_u, sgejsv_obj->threshold);
    EXPECT_NEAR(0.0, sgejsv_obj->diff_v, sgejsv_obj->threshold);
    EXPECT_NEAR(0.0, sgejsv_obj->diff_stat, sgejsv_obj->threshold);}


/* Begin gejsv_double_common_parameters  class definition */
class gejsv_double_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    double diff_a, diff_sva, diff_u, diff_v, diff_stat;
    int min_mn;
    float threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
	char joba; //  Must be 'C', 'E', 'F', 'G', 'A', or 'R'.
	char jobu; // Must be 'U', 'F', 'W', or 'N'.
	char jobv; // Must be 'V', 'J', 'W', or 'N'.
	char jobr; // Must be 'N' or 'R'.
	char jobt; // Must be 'T' or 'N'.
	char jobp; //  Must be 'P' or 'N'.

    lapack_int m; // The number of rows of the matrix A
    lapack_int n; // The number of columns of the matrix A 
    lapack_int lda; //  The leading dimension of a

    lapack_int ldu; // leading dimension of the output array u . ldu≥ max(1, m)
    lapack_int ldv; // The leading dimension of the output array v . ldv≥ max(1, p)
	
    /* Input / Output parameters */
    double* a, *aref; // contains m-by-n matrix A.
    double* u, *uref; //  contains the m-by-n matrix of the left singular vectors.
    double* v, *vref; //  contains the n-by-n matrix of the right singular vectors.
    double* all, *allref; // contains work array.
	
    /* Output parameters */
    double* sva, *svaref; //   the singular values of A.
    double* stat, *statref; //  scaling factor
    lapack_int * istat, *istatref; // the numerical rank determined after the initial QR factorization
    
    /*Return Values */
    int info, inforef;

   public:
      gejsv_double_parameters (int matrix_layout_i, char joba, char jobu, 
							  char jobv,char jobr, char jobt, char jobp,
                                             lapack_int m, lapack_int n);
      ~gejsv_double_parameters ();
};

/* Constructor definition  gejsv double_common_parameters */
gejsv_double_parameters:: gejsv_double_parameters (int matrix_layout_i, char joba_i, 
							  char jobu_i, char jobv_i,char jobr_i, char jobt_i,
							  char jobp_i, lapack_int m_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    jobu = jobu_i;
    jobv = jobv_i;
    joba = joba_i;
    jobp = jobp_i;
    jobt = jobt_i;
    jobr = jobr_i;
    n = n_i;
    m = m_i;
    min_mn = (m<n)?m:n;
    
    ldu = m;
    ldv = n;

    if (matrix_layout == LAPACK_ROW_MAJOR)
    {
        lda = n;
        if( jobu== 'F')  ldu = m; 

    }   else
    {
        lda = m;
    }

    // Initialize with default values
    hModule = NULL;
    dModule = NULL;
    diff_a = 0;
    diff_sva = 0;
    diff_u = 0;
    diff_v = 0;
    diff_stat = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n gejsv double: matrix_layout: %d   jobu: %c \t \
jobv: %c \t joba: %c \n jobt: %c \t jobp: %c \t jobr: %c \t m: %d \t n: %d \n",
matrix_layout, jobu_i, jobv_i, joba_i, jobt_i, jobp_i, jobr_i, m_i, n_i);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, m*n );
    lapacke_gtest_alloc_double_buffer_pair( &all, &allref, n );
    lapacke_gtest_alloc_double_buffer_pair( &sva, &svaref, n );
    lapacke_gtest_alloc_double_buffer_pair( &u, &uref, m*m ); // allocating largest size
    lapacke_gtest_alloc_double_buffer_pair( &v, &vref, ldv*n );
    lapacke_gtest_alloc_double_buffer_pair( &stat, &statref, n );
    lapacke_gtest_alloc_int_buffer_pair( &istat, &istatref, min_mn );

    if( (a==NULL) || (aref==NULL) ||  \
        (sva==NULL) || (svaref==NULL) || \
        (u==NULL) || (uref==NULL) || \
        (all==NULL) || (allref==NULL) || \
        (v==NULL) || (vref==NULL) || \
        (istat==NULL) || (istatref==NULL) || \
        (stat==NULL) || (statref==NULL)  ){
       EXPECT_FALSE( true) << "gejsv_double_parameters object: malloc error. Exiting ";
       gejsv_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, m*n );
    lapacke_gtest_init_double_buffer_pair_with_constant( u, uref, m*ldu, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant( v, vref, ldv*n, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant( sva, svaref, n, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant( stat, statref, n, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant( all, allref, n, 0.0);
    lapacke_gtest_init_int_buffer_pair_with_constant( istat, istatref, n, 0);

   } /* end of Constructor  */

/* Destructor definition  'gejsv_double_common_parameters' */
gejsv_double_parameters :: ~gejsv_double_parameters ()
{
   /* De-Allocate memory for the input matrices */
   gejsv_free();
} 

//  Test fixture class definition
class dgejsv_test  : public  ::testing::Test {
public:
   gejsv_double_parameters  *dgejsv_obj;
   void SetUp();  
   void TearDown () { delete dgejsv_obj; }
};

void dgejsv_test::SetUp()
{
    /* LAPACKE dgejsv prototype */
    typedef int (*Fptr_NL_LAPACKE_dgejsv) ( int matrix_layout, char joba,
			char jobu, char jobv, char jobr, char jobt, char jobp,
			lapack_int m, lapack_int n, double * a, lapack_int lda,
			double * all, double * u, lapack_int ldu, double * v,
			lapack_int ldv, double * stat, lapack_int * state);
            
    Fptr_NL_LAPACKE_dgejsv DGEJSV;
    dgejsv_obj = new  gejsv_double_parameters( svd_paramslist[idx].matrix_layout,
                                         svd_paramslist[idx].joba_gejsv,
                                         svd_paramslist[idx].jobu_gejsv,
                                         svd_paramslist[idx].jobv_gejsv,
                                         svd_paramslist[idx].jobr_gejsv,
                                         svd_paramslist[idx].jobt_gejsv,
                                         svd_paramslist[idx].jobp_gejsv,
                                         svd_paramslist[idx].m_gejsv,
                                         svd_paramslist[idx].n_gejsv );

    dgejsv_obj->threshold = svd_paramslist[idx].svd_threshold;
    idx = Circular_Increment_Index(idx);
    dgejsv_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgejsv_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgejsv_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgejsv_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DGEJSV = (Fptr_NL_LAPACKE_dgejsv)dlsym(dgejsv_obj->hModule, "LAPACKE_dgejsv");
    ASSERT_TRUE(DGEJSV != NULL) << "failed to ppt the Netlib LAPACKE_dgejsv symbol";

    /* Compute libflame's Lapacke o/p  */
    dgejsv_obj->info = LAPACKE_dgejsv(  dgejsv_obj->matrix_layout,
                                    dgejsv_obj->joba,
                                    dgejsv_obj->jobu,
                                    dgejsv_obj->jobv,
                                    dgejsv_obj->jobr,
                                    dgejsv_obj->jobt,
                                    dgejsv_obj->jobp,
                                    dgejsv_obj->m,
                                    dgejsv_obj->n,
                                    dgejsv_obj->a,
                                    dgejsv_obj->lda,
                                    dgejsv_obj->all,
                                    dgejsv_obj->u,
                                    dgejsv_obj->ldu,
                                    dgejsv_obj->v,
                                    dgejsv_obj->ldv,
                                    dgejsv_obj->stat,
									dgejsv_obj->istat
                                    );

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    dgejsv_obj->inforef = DGEJSV( dgejsv_obj->matrix_layout,
                                    dgejsv_obj->joba,
                                    dgejsv_obj->jobu,
                                    dgejsv_obj->jobv,
                                    dgejsv_obj->jobr,
                                    dgejsv_obj->jobt,
                                    dgejsv_obj->jobp,
                                    dgejsv_obj->m,
                                    dgejsv_obj->n,
                                    dgejsv_obj->aref,
                                    dgejsv_obj->lda,
                                    dgejsv_obj->allref,
                                    dgejsv_obj->uref,
                                    dgejsv_obj->ldu,
                                    dgejsv_obj->vref,
                                    dgejsv_obj->ldv,
                                    dgejsv_obj->statref,
									dgejsv_obj->istatref
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    dgejsv_obj->diff_u =  computeDiff_d( (dgejsv_obj->m)*(dgejsv_obj->m), 
                dgejsv_obj->u, dgejsv_obj->uref );

    dgejsv_obj->diff_v =  computeDiff_d( (dgejsv_obj->ldv)*(dgejsv_obj->n), 
                dgejsv_obj->v, dgejsv_obj->vref );

    dgejsv_obj->diff_sva =  computeDiff_d( dgejsv_obj->n, 
                dgejsv_obj->all, dgejsv_obj->allref );

    dgejsv_obj->diff_stat =  computeDiff_d( dgejsv_obj->min_mn, 
                dgejsv_obj->stat, dgejsv_obj->statref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n gejsv double: \n  diff_sva: %f \n \
diff_u: %f \n diff_v: %f \n diff_stat: %f info: %d \t inforef: %d\n ",
       dgejsv_obj->diff_sva, dgejsv_obj->diff_u,
       dgejsv_obj->diff_v, dgejsv_obj->diff_stat,
	   dgejsv_obj->info, dgejsv_obj->inforef);
#endif
}

TEST_F(dgejsv_test, dgejsv_1) {
    EXPECT_NEAR(0.0, dgejsv_obj->diff_sva, dgejsv_obj->threshold);
    EXPECT_NEAR(0.0, dgejsv_obj->diff_u, dgejsv_obj->threshold);
    EXPECT_NEAR(0.0, dgejsv_obj->diff_v, dgejsv_obj->threshold);
    EXPECT_NEAR(0.0, dgejsv_obj->diff_stat, dgejsv_obj->threshold);
}

TEST_F(dgejsv_test, dgejsv_2) {
    EXPECT_NEAR(0.0, dgejsv_obj->diff_sva, dgejsv_obj->threshold);
    EXPECT_NEAR(0.0, dgejsv_obj->diff_u, dgejsv_obj->threshold);
    EXPECT_NEAR(0.0, dgejsv_obj->diff_v, dgejsv_obj->threshold);
    EXPECT_NEAR(0.0, dgejsv_obj->diff_stat, dgejsv_obj->threshold);}

TEST_F(dgejsv_test, dgejsv_3) {
    EXPECT_NEAR(0.0, dgejsv_obj->diff_sva, dgejsv_obj->threshold);
    EXPECT_NEAR(0.0, dgejsv_obj->diff_u, dgejsv_obj->threshold);
    EXPECT_NEAR(0.0, dgejsv_obj->diff_v, dgejsv_obj->threshold);
    EXPECT_NEAR(0.0, dgejsv_obj->diff_stat, dgejsv_obj->threshold);}

TEST_F(dgejsv_test, dgejsv_4) {
    EXPECT_NEAR(0.0, dgejsv_obj->diff_sva, dgejsv_obj->threshold);
    EXPECT_NEAR(0.0, dgejsv_obj->diff_u, dgejsv_obj->threshold);
    EXPECT_NEAR(0.0, dgejsv_obj->diff_v, dgejsv_obj->threshold);
    EXPECT_NEAR(0.0, dgejsv_obj->diff_stat, dgejsv_obj->threshold);}

/* Begin gejsv_scomplex_common_parameters  class definition */
class gejsv_scomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_a, diff_sva, diff_u, diff_v, diff_stat;
    int min_mn;
    float threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
	char joba; //  Must be 'C', 'E', 'F', 'G', 'A', or 'R'.
	char jobu; // Must be 'U', 'F', 'W', or 'N'.
	char jobv; // Must be 'V', 'J', 'W', or 'N'.
	char jobr; // Must be 'N' or 'R'.
	char jobt; // Must be 'T' or 'N'.
	char jobp; //  Must be 'P' or 'N'.

    lapack_int m; // The number of rows of the matrix A
    lapack_int n; // The number of columns of the matrix A 
    lapack_int lda; //  The leading dimension of a

    lapack_int ldu; // leading dimension of the output array u . ldu≥ max(1, m)
    lapack_int ldv; // The leading dimension of the output array v . ldv≥ max(1, p)
	
    /* Input / Output parameters */
    lapack_complex_float* a, *aref; // contains m-by-n matrix A.
    lapack_complex_float* u, *uref; //  contains the m-by-n matrix of the left singular vectors.
    lapack_complex_float* v, *vref; //  contains the n-by-n matrix of the right singular vectors.
    float* all, *allref; // contains work array.
	
    /* Output parameters */
    float* sva, *svaref; //   the singular values of A.
    float* stat, *statref; //  scaling factor
    lapack_int * istat, *istatref; // the numerical rank determined after the initial QR factorization
    
    /*Return Values */
    int info, inforef;

   public:
      gejsv_scomplex_parameters (int matrix_layout_i, char joba, char jobu, 
							  char jobv,char jobr, char jobt, char jobp,
                                             lapack_int m, lapack_int n);
      ~gejsv_scomplex_parameters ();
};

/* Constructor definition  gejsv lapack_complex_float_common_parameters */
gejsv_scomplex_parameters:: gejsv_scomplex_parameters (int matrix_layout_i, char joba_i, 
							  char jobu_i, char jobv_i,char jobr_i, char jobt_i,
							  char jobp_i, lapack_int m_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    jobu = jobu_i;
    jobv = jobv_i;
    joba = joba_i;
    jobp = jobp_i;
    jobt = jobt_i;
    jobr = jobr_i;
    n = n_i;
    m = m_i;
    min_mn = (m<n)?m:n;
    
    ldu = m;
    ldv = n;

    if (matrix_layout == LAPACK_ROW_MAJOR)
    {
        lda = n;
        if( jobu== 'F')  ldu = m; 

    }   else
    {
        lda = m;
    }

    // Initialize with default values
    hModule = NULL;
    dModule = NULL;
    diff_a = 0;
    diff_sva = 0;
    diff_u = 0;
    diff_v = 0;
    diff_stat = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n gejsv lapack_complex_float: matrix_layout: %d   jobu: %c \t \
jobv: %c \t joba: %c \n jobt: %c \t jobp: %c \t jobr: %c \t m: %d \t n: %d \n",
matrix_layout, jobu_i, jobv_i, joba_i, jobt_i, jobp_i, jobr_i, m_i, n_i);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, m*n );
    lapacke_gtest_alloc_float_buffer_pair( &all, &allref, n );
    lapacke_gtest_alloc_float_buffer_pair( &sva, &svaref, n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &u, &uref, m*m ); // allocating largest size
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &v, &vref, ldv*n );
    lapacke_gtest_alloc_float_buffer_pair( &stat, &statref, n );
    lapacke_gtest_alloc_int_buffer_pair( &istat, &istatref, min_mn );

    if( (a==NULL) || (aref==NULL) ||  \
        (sva==NULL) || (svaref==NULL) || \
        (u==NULL) || (uref==NULL) || \
        (all==NULL) || (allref==NULL) || \
        (v==NULL) || (vref==NULL) || \
        (istat==NULL) || (istatref==NULL) || \
        (stat==NULL) || (statref==NULL)  ){
       EXPECT_FALSE( true) << "gejsv_scomplex_parameters object: malloc error. Exiting ";
       gejsv_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, m*n );
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( u, uref, m*ldu, 0.0);
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( v, vref, ldv*n, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant( sva, svaref, n, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant( stat, statref, n, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant( all, allref, n, 0.0);
    lapacke_gtest_init_int_buffer_pair_with_constant( istat, istatref, n, 0);

   } /* end of Constructor  */

/* Destructor definition  'gejsv_scomplex_common_parameters' */
gejsv_scomplex_parameters :: ~gejsv_scomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   gejsv_free();
} 

//  Test fixture class definition
class cgejsv_test  : public  ::testing::Test {
public:
   gejsv_scomplex_parameters  *cgejsv_obj;
   void SetUp();  
   void TearDown () { delete cgejsv_obj; }
};

void cgejsv_test::SetUp()
{
    /* LAPACKE cgejsv prototype */
    typedef int (*Fptr_NL_LAPACKE_cgejsv) ( int matrix_layout, char joba,
			char jobu, char jobv, char jobr, char jobt, char jobp,
			lapack_int m, lapack_int n, lapack_complex_float * a, lapack_int lda,
			float * all, lapack_complex_float * u, lapack_int ldu, lapack_complex_float * v,
			lapack_int ldv, float * stat, lapack_int * state);
            
    Fptr_NL_LAPACKE_cgejsv CGEJSV;
    cgejsv_obj = new  gejsv_scomplex_parameters( svd_paramslist[idx].matrix_layout,
                                         svd_paramslist[idx].joba_gejsv,
                                         svd_paramslist[idx].jobu_gejsv,
                                         svd_paramslist[idx].jobv_gejsv,
                                         svd_paramslist[idx].jobr_gejsv,
                                         svd_paramslist[idx].jobt_gejsv,
                                         svd_paramslist[idx].jobp_gejsv,
                                         svd_paramslist[idx].m_gejsv,
                                         svd_paramslist[idx].n_gejsv );

    //idx = Circular_Increment_Index(idx);
    cgejsv_obj->threshold = svd_paramslist[idx].svd_threshold;
    cgejsv_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgejsv_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgejsv_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgejsv_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CGEJSV = (Fptr_NL_LAPACKE_cgejsv)dlsym(cgejsv_obj->hModule, "LAPACKE_cgejsv");
    ASSERT_TRUE(CGEJSV != NULL) << "failed to ppt the Netlib LAPACKE_cgejsv symbol";

    /* Compute libflame's Lapacke o/p  */
    cgejsv_obj->info = LAPACKE_cgejsv(  cgejsv_obj->matrix_layout,
                                    cgejsv_obj->joba,
                                    cgejsv_obj->jobu,
                                    cgejsv_obj->jobv,
                                    cgejsv_obj->jobr,
                                    cgejsv_obj->jobt,
                                    cgejsv_obj->jobp,
                                    cgejsv_obj->m,
                                    cgejsv_obj->n,
                                    cgejsv_obj->a,
                                    cgejsv_obj->lda,
                                    cgejsv_obj->all,
                                    cgejsv_obj->u,
                                    cgejsv_obj->ldu,
                                    cgejsv_obj->v,
                                    cgejsv_obj->ldv,
                                    cgejsv_obj->stat,
									cgejsv_obj->istat
                                    );

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    cgejsv_obj->inforef = CGEJSV( cgejsv_obj->matrix_layout,
                                    cgejsv_obj->joba,
                                    cgejsv_obj->jobu,
                                    cgejsv_obj->jobv,
                                    cgejsv_obj->jobr,
                                    cgejsv_obj->jobt,
                                    cgejsv_obj->jobp,
                                    cgejsv_obj->m,
                                    cgejsv_obj->n,
                                    cgejsv_obj->aref,
                                    cgejsv_obj->lda,
                                    cgejsv_obj->allref,
                                    cgejsv_obj->uref,
                                    cgejsv_obj->ldu,
                                    cgejsv_obj->vref,
                                    cgejsv_obj->ldv,
                                    cgejsv_obj->statref,
									cgejsv_obj->istatref
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    cgejsv_obj->diff_u =  computeDiff_c( (cgejsv_obj->m)*(cgejsv_obj->m), 
                cgejsv_obj->u, cgejsv_obj->uref );

    cgejsv_obj->diff_v =  computeDiff_c( (cgejsv_obj->ldv)*(cgejsv_obj->n), 
                cgejsv_obj->v, cgejsv_obj->vref );

    cgejsv_obj->diff_sva =  computeDiff_s( cgejsv_obj->n, 
                cgejsv_obj->all, cgejsv_obj->allref );

    cgejsv_obj->diff_stat =  computeDiff_s( cgejsv_obj->min_mn, 
                cgejsv_obj->stat, cgejsv_obj->statref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n gejsv lapack_complex_float: \n  diff_sva: %f \n \
diff_u: %f \n diff_v: %f \n diff_stat: %f info: %d \t inforef: %d\n ",
       cgejsv_obj->diff_sva, cgejsv_obj->diff_u,
       cgejsv_obj->diff_v, cgejsv_obj->diff_stat,
	   cgejsv_obj->info, cgejsv_obj->inforef);
#endif
}

TEST_F(cgejsv_test, cgejsv_1) {
    EXPECT_NEAR(0.0, cgejsv_obj->diff_sva, cgejsv_obj->threshold);
    EXPECT_NEAR(0.0, cgejsv_obj->diff_u, cgejsv_obj->threshold);
    EXPECT_NEAR(0.0, cgejsv_obj->diff_v, cgejsv_obj->threshold);
    EXPECT_NEAR(0.0, cgejsv_obj->diff_stat, cgejsv_obj->threshold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(cgejsv_test, cgejsv_2) {
    EXPECT_NEAR(0.0, cgejsv_obj->diff_sva, cgejsv_obj->threshold);
    EXPECT_NEAR(0.0, cgejsv_obj->diff_u, cgejsv_obj->threshold);
    EXPECT_NEAR(0.0, cgejsv_obj->diff_v, cgejsv_obj->threshold);
    EXPECT_NEAR(0.0, cgejsv_obj->diff_stat, cgejsv_obj->threshold);
    idx = Circular_Increment_Index(idx);
	
	// The below additional increment to avoid the option 'jobu =W' 
	// which is not suppoted in complex, double complex types
    idx = Circular_Increment_Index(idx);
}

TEST_F(cgejsv_test, cgejsv_3) {

    EXPECT_NEAR(0.0, cgejsv_obj->diff_sva, cgejsv_obj->threshold);
    EXPECT_NEAR(0.0, cgejsv_obj->diff_u, cgejsv_obj->threshold);
    EXPECT_NEAR(0.0, cgejsv_obj->diff_v, cgejsv_obj->threshold);
    EXPECT_NEAR(0.0, cgejsv_obj->diff_stat, cgejsv_obj->threshold);
    idx = Circular_Increment_Index(idx);

	}

TEST_F(cgejsv_test, cgejsv_4) {
    EXPECT_NEAR(0.0, cgejsv_obj->diff_sva, cgejsv_obj->threshold);
    EXPECT_NEAR(0.0, cgejsv_obj->diff_u, cgejsv_obj->threshold);
    EXPECT_NEAR(0.0, cgejsv_obj->diff_v, cgejsv_obj->threshold);
    EXPECT_NEAR(0.0, cgejsv_obj->diff_stat, cgejsv_obj->threshold);
    idx = 0;
	}

/* Begin gejsv_dcomplex_common_parameters  class definition */
class gejsv_dcomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    double diff_a, diff_sva, diff_u, diff_v, diff_stat;
    int min_mn;
    float threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
	char joba; //  Must be 'C', 'E', 'F', 'G', 'A', or 'R'.
	char jobu; // Must be 'U', 'F', 'W', or 'N'.
	char jobv; // Must be 'V', 'J', 'W', or 'N'.
	char jobr; // Must be 'N' or 'R'.
	char jobt; // Must be 'T' or 'N'.
	char jobp; //  Must be 'P' or 'N'.

    lapack_int m; // The number of rows of the matrix A
    lapack_int n; // The number of columns of the matrix A 
    lapack_int lda; //  The leading dimension of a

    lapack_int ldu; // leading dimension of the output array u . ldu≥ max(1, m)
    lapack_int ldv; // The leading dimension of the output array v . ldv≥ max(1, p)
	
    /* Input / Output parameters */
    lapack_complex_double* a, *aref; // contains m-by-n matrix A.
    lapack_complex_double* u, *uref; //  contains the m-by-n matrix of the left singular vectors.
    lapack_complex_double* v, *vref; //  contains the n-by-n matrix of the right singular vectors.
    double* all, *allref; // contains work array.
	
    /* Output parameters */
    double* sva, *svaref; //   the singular values of A.
    double* stat, *statref; //  scaling factor
    lapack_int * istat, *istatref; // the numerical rank determined after the initial QR factorization
    
    /*Return Values */
    int info, inforef;

   public:
      gejsv_dcomplex_parameters (int matrix_layout_i, char joba, char jobu, 
							  char jobv,char jobr, char jobt, char jobp,
                                             lapack_int m, lapack_int n);
      ~gejsv_dcomplex_parameters ();
};

/* Constructor definition  gejsv lapack_complex_double_common_parameters */
gejsv_dcomplex_parameters:: gejsv_dcomplex_parameters (int matrix_layout_i, char joba_i, 
							  char jobu_i, char jobv_i,char jobr_i, char jobt_i,
							  char jobp_i, lapack_int m_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    jobu = jobu_i;
    jobv = jobv_i;
    joba = joba_i;
    jobp = jobp_i;
    jobt = jobt_i;
    jobr = jobr_i;
    n = n_i;
    m = m_i;
    min_mn = (m<n)?m:n;
    
    ldu = m;
    ldv = n;

    if (matrix_layout == LAPACK_ROW_MAJOR)
    {
        lda = n;
        if( jobu== 'F')  ldu = m; 

    }   else
    {
        lda = m;
    }

    // Initialize with default values
    hModule = NULL;
    dModule = NULL;
    diff_a = 0;
    diff_sva = 0;
    diff_u = 0;
    diff_v = 0;
    diff_stat = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n gejsv lapack_complex_double: matrix_layout: %d   jobu: %c \t \
jobv: %c \t joba: %c \n jobt: %c \t jobp: %c \t jobr: %c \t m: %d \t n: %d \n",
matrix_layout, jobu_i, jobv_i, joba_i, jobt_i, jobp_i, jobr_i, m_i, n_i);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, m*n );
    lapacke_gtest_alloc_double_buffer_pair( &all, &allref, n );
    lapacke_gtest_alloc_double_buffer_pair( &sva, &svaref, n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &u, &uref, m*m ); // allocating largest size
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &v, &vref, ldv*n );
    lapacke_gtest_alloc_double_buffer_pair( &stat, &statref, n );
    lapacke_gtest_alloc_int_buffer_pair( &istat, &istatref, min_mn );

    if( (a==NULL) || (aref==NULL) ||  \
        (sva==NULL) || (svaref==NULL) || \
        (u==NULL) || (uref==NULL) || \
        (all==NULL) || (allref==NULL) || \
        (v==NULL) || (vref==NULL) || \
        (istat==NULL) || (istatref==NULL) || \
        (stat==NULL) || (statref==NULL)  ){
       EXPECT_FALSE( true) << "gejsv_dcomplex_parameters object: malloc error. Exiting ";
       gejsv_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, m*n );
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( u, uref, m*ldu, 0.0);
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( v, vref, ldv*n, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant( sva, svaref, n, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant( stat, statref, n, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant( all, allref, n, 0.0);
    lapacke_gtest_init_int_buffer_pair_with_constant( istat, istatref, n, 0);

   } /* end of Constructor  */

/* Destructor definition  'gejsv_dcomplex_common_parameters' */
gejsv_dcomplex_parameters :: ~gejsv_dcomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   gejsv_free();
} 

//  Test fixture class definition
class zgejsv_test  : public  ::testing::Test {
public:
   gejsv_dcomplex_parameters  *zgejsv_obj;
   void SetUp();  
   void TearDown () { delete zgejsv_obj; }
};

void zgejsv_test::SetUp()
{
    /* LAPACKE zgejsv prototype */
    typedef int (*Fptr_NL_LAPACKE_zgejsv) ( int matrix_layout, char joba,
			char jobu, char jobv, char jobr, char jobt, char jobp,
			lapack_int m, lapack_int n, lapack_complex_double * a, lapack_int lda,
			double * all, lapack_complex_double * u, lapack_int ldu, lapack_complex_double * v,
			lapack_int ldv, double * stat, lapack_int * state);
            
    Fptr_NL_LAPACKE_zgejsv ZGEJSV;
    zgejsv_obj = new  gejsv_dcomplex_parameters( svd_paramslist[idx].matrix_layout,
                                         svd_paramslist[idx].joba_gejsv,
                                         svd_paramslist[idx].jobu_gejsv,
                                         svd_paramslist[idx].jobv_gejsv,
                                         svd_paramslist[idx].jobr_gejsv,
                                         svd_paramslist[idx].jobt_gejsv,
                                         svd_paramslist[idx].jobp_gejsv,
                                         svd_paramslist[idx].m_gejsv,
                                         svd_paramslist[idx].n_gejsv );

    //idx = Circular_Increment_Index(idx);
    zgejsv_obj->threshold = svd_paramslist[idx].svd_threshold;
    zgejsv_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgejsv_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgejsv_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgejsv_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZGEJSV = (Fptr_NL_LAPACKE_zgejsv)dlsym(zgejsv_obj->hModule, "LAPACKE_zgejsv");
    ASSERT_TRUE(ZGEJSV != NULL) << "failed to ppt the Netlib LAPACKE_zgejsv symbol";

    /* Compute libflame's Lapacke o/p  */
    zgejsv_obj->info = LAPACKE_zgejsv(  zgejsv_obj->matrix_layout,
                                    zgejsv_obj->joba,
                                    zgejsv_obj->jobu,
                                    zgejsv_obj->jobv,
                                    zgejsv_obj->jobr,
                                    zgejsv_obj->jobt,
                                    zgejsv_obj->jobp,
                                    zgejsv_obj->m,
                                    zgejsv_obj->n,
                                    zgejsv_obj->a,
                                    zgejsv_obj->lda,
                                    zgejsv_obj->all,
                                    zgejsv_obj->u,
                                    zgejsv_obj->ldu,
                                    zgejsv_obj->v,
                                    zgejsv_obj->ldv,
                                    zgejsv_obj->stat,
									zgejsv_obj->istat
                                    );

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    zgejsv_obj->inforef = ZGEJSV( zgejsv_obj->matrix_layout,
                                    zgejsv_obj->joba,
                                    zgejsv_obj->jobu,
                                    zgejsv_obj->jobv,
                                    zgejsv_obj->jobr,
                                    zgejsv_obj->jobt,
                                    zgejsv_obj->jobp,
                                    zgejsv_obj->m,
                                    zgejsv_obj->n,
                                    zgejsv_obj->aref,
                                    zgejsv_obj->lda,
                                    zgejsv_obj->allref,
                                    zgejsv_obj->uref,
                                    zgejsv_obj->ldu,
                                    zgejsv_obj->vref,
                                    zgejsv_obj->ldv,
                                    zgejsv_obj->statref,
									zgejsv_obj->istatref
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    zgejsv_obj->diff_u =  computeDiff_z( (zgejsv_obj->m)*(zgejsv_obj->m), 
                zgejsv_obj->u, zgejsv_obj->uref );

    zgejsv_obj->diff_v =  computeDiff_z( (zgejsv_obj->ldv)*(zgejsv_obj->n), 
                zgejsv_obj->v, zgejsv_obj->vref );

    zgejsv_obj->diff_sva =  computeDiff_d( zgejsv_obj->n, 
                zgejsv_obj->all, zgejsv_obj->allref );

    zgejsv_obj->diff_stat =  computeDiff_d( zgejsv_obj->min_mn, 
                zgejsv_obj->stat, zgejsv_obj->statref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n gejsv lapack_complex_double: \n  diff_sva: %f \n \
diff_u: %f \n diff_v: %f \n diff_stat: %f info: %d \t inforef: %d\n ",
       zgejsv_obj->diff_sva, zgejsv_obj->diff_u,
       zgejsv_obj->diff_v, zgejsv_obj->diff_stat,
	   zgejsv_obj->info, zgejsv_obj->inforef);
#endif
}

TEST_F(zgejsv_test, zgejsv_1) {
    EXPECT_NEAR(0.0, zgejsv_obj->diff_sva, zgejsv_obj->threshold);
    EXPECT_NEAR(0.0, zgejsv_obj->diff_u, zgejsv_obj->threshold);
    EXPECT_NEAR(0.0, zgejsv_obj->diff_v, zgejsv_obj->threshold);
    EXPECT_NEAR(0.0, zgejsv_obj->diff_stat, zgejsv_obj->threshold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(zgejsv_test, zgejsv_2) {
    EXPECT_NEAR(0.0, zgejsv_obj->diff_sva, zgejsv_obj->threshold);
    EXPECT_NEAR(0.0, zgejsv_obj->diff_u, zgejsv_obj->threshold);
    EXPECT_NEAR(0.0, zgejsv_obj->diff_v, zgejsv_obj->threshold);
    EXPECT_NEAR(0.0, zgejsv_obj->diff_stat, zgejsv_obj->threshold);
    idx = Circular_Increment_Index(idx);
	
	// The below additional increment to avoid the option 'jobu =W' 
	// which is not suppoted in complex, double complex types
    idx = Circular_Increment_Index(idx);
}

TEST_F(zgejsv_test, zgejsv_3) {

    EXPECT_NEAR(0.0, zgejsv_obj->diff_sva, zgejsv_obj->threshold);
    EXPECT_NEAR(0.0, zgejsv_obj->diff_u, zgejsv_obj->threshold);
    EXPECT_NEAR(0.0, zgejsv_obj->diff_v, zgejsv_obj->threshold);
    EXPECT_NEAR(0.0, zgejsv_obj->diff_stat, zgejsv_obj->threshold);
    idx = Circular_Increment_Index(idx);

	}

TEST_F(zgejsv_test, zgejsv_4) {
    EXPECT_NEAR(0.0, zgejsv_obj->diff_sva, zgejsv_obj->threshold);
    EXPECT_NEAR(0.0, zgejsv_obj->diff_u, zgejsv_obj->threshold);
    EXPECT_NEAR(0.0, zgejsv_obj->diff_v, zgejsv_obj->threshold);
    EXPECT_NEAR(0.0, zgejsv_obj->diff_stat, zgejsv_obj->threshold);
    idx = 0;
	}
