#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"
#define LAPACKE_TEST_VERBOSE  (1)
#define gesvdx_free() \
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

/* Begin gesvdx_float_common_parameters  class definition */
class gesvdx_float_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_a, diff_s, diff_u, diff_vt, diff_superb;
    int min_mn;
    float threshold;
    void *hModule, *dModule;
    
    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char jobu; // Must be 'V', or 'N'. 
    char jobvt; // Must be 'V', 'N'. 
    char range; // Must be 'A', 'V', 'I'. 

    lapack_int m; // The number of rows of the matrix A
    lapack_int n; // The number of columns of the matrix A 
    lapack_int lda; //  The leading dimension of a

    lapack_int ldu; // leading dimension of the output array u . ldu≥ max(1, m)
    lapack_int ldvt; // The leading dimension of the output array vt . ldvt≥ max(1, p) 
    float vl; // the lower and upper bounds of the interval to be searched for singular values
	float vu;
   /* Input / Output parameters */
    float* a, *aref; // contains m-by-n matrix A.
	lapack_int il; // the indices of the smallest and largest singular values to be returned.
	lapack_int iu;
    /* Output parameters */
	lapack_int ns, nsref;
    lapack_int *superb, *superbref; // the dimension of subblocks.
    // Below buffers contain the o/p orthogonal/unitary matrices
    float* u, *uref;
    float* vt, *vtref;
    float* s, *sref;
    
    /*Return Values */
    int info, inforef;

   public:
      gesvdx_float_parameters (int matrix_layout_i, char jobu, char jobvt,
                                  char range, lapack_int m, lapack_int n,
								  lapack_int il, lapack_int iu,
								  float vl, float vu);
      ~gesvdx_float_parameters ();
};

/* Constructor definition  gesvdx float_common_parameters */
gesvdx_float_parameters:: gesvdx_float_parameters (int matrix_layout_i,
                    char jobu_i, char jobvt_i, char range_i, lapack_int m_i,  
                    lapack_int n_i, lapack_int il_i, lapack_int iu_i,
								  float vl_i, float vu_i)
{
    matrix_layout = matrix_layout_i;
    jobu = jobu_i;
    jobvt = jobvt_i;
	range = range_i;
    n = n_i;
    m = m_i;
    min_mn = (m<n)?m:n;
	il = il_i;
	iu = iu_i;
	vl = vl_i;
	vu = vu_i;
    
    ldu = m;
    ldvt = n;
	ns = 0;
	nsref = 0;

    if (matrix_layout == LAPACK_ROW_MAJOR)
    {
        lda = n;
    }   else   {
        lda = m;
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
   printf(" \n gesvdx float: matrix_layout: %d   jobu: %c \t \
jobvt: %c \t range: %c \n m: %d \t n: %d \t il:%d \t iu:%d \t \
vl: %f \t vu: %f \n",
matrix_layout, jobu_i, jobvt_i,range, m_i, n_i, il, iu, vl, vu);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, m*n );
    lapacke_gtest_alloc_float_buffer_pair( &s, &sref, min_mn );
    lapacke_gtest_alloc_float_buffer_pair( &u, &uref, m*m );
    lapacke_gtest_alloc_float_buffer_pair( &vt, &vtref, n*n );
    lapacke_gtest_alloc_int_buffer_pair( &superb, &superbref, 12*min_mn );

    if( (a==NULL) || (aref==NULL) ||  \
        (s==NULL) || (sref==NULL) || \
        (u==NULL) || (uref==NULL) || \
        (vt==NULL) || (vtref==NULL) || \
        (superb==NULL) || (superbref==NULL)  ){
       EXPECT_FALSE( true) << "gesvdx_float_parameters object: malloc error. Exiting ";
       gesvdx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, m*n );
    lapacke_gtest_init_float_buffer_pair_with_constant( u, uref, m*m, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant( vt, vtref, n*n, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant( s, sref, min_mn, 0.0);
    lapacke_gtest_init_int_buffer_pair_with_constant( superb, superbref, 12*min_mn, 0);

   } /* end of Constructor  */

/* Destructor definition  'gesvdx_float_common_parameters' */
gesvdx_float_parameters :: ~gesvdx_float_parameters ()
{
   /* De-Allocate memory for the input matrices */
   gesvdx_free();
} 

//  Test fixture class definition
class sgesvdx_test  : public  ::testing::Test {
public:
   gesvdx_float_parameters  *sgesvdx_obj;
   void SetUp();  
   void TearDown () { delete sgesvdx_obj; }
};

void sgesvdx_test::SetUp()
{
    /* LAPACKE sgesvdx prototype */

	typedef int (*Fptr_NL_LAPACKE_sgesvdx) (int matrix_layout, char jobu,
		char jobvt, char range, lapack_int m, lapack_int n, float * a, 
		lapack_int lda, float vl, float vu, lapack_int il, lapack_int iu,
		lapack_int * ns, float * s, float * u, lapack_int ldu, float * vt,
		lapack_int ldvt, lapack_int * superb);

    Fptr_NL_LAPACKE_sgesvdx SGESVDX;

    sgesvdx_obj = new  gesvdx_float_parameters( svd_paramslist[idx].matrix_layout,
                                         svd_paramslist[idx].jobu_gesvdx,
                                         svd_paramslist[idx].jobvt_gesvdx,
                                         svd_paramslist[idx].range_gesvdx,
                                         svd_paramslist[idx].m_gesvj,
                                         svd_paramslist[idx].n_gesvj,
                                         svd_paramslist[idx].il,
                                         svd_paramslist[idx].iu,
                                         svd_paramslist[idx].vl,
                                         svd_paramslist[idx].vu
										 );

    idx = Circular_Increment_Index(idx);
    sgesvdx_obj->threshold = svd_paramslist[idx].svd_threshold;
    sgesvdx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgesvdx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgesvdx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgesvdx_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SGESVDX = (Fptr_NL_LAPACKE_sgesvdx)dlsym(sgesvdx_obj->hModule, "LAPACKE_sgesvdx");
    ASSERT_TRUE(SGESVDX != NULL) << "failed to ppt the Netlib LAPACKE_sgesvdx symbol";
 
 /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    sgesvdx_obj->inforef = SGESVDX( sgesvdx_obj->matrix_layout,
                                    sgesvdx_obj->jobu,
                                    sgesvdx_obj->jobvt,
                                    sgesvdx_obj->range,
                                    sgesvdx_obj->m,
                                    sgesvdx_obj->n,
                                    sgesvdx_obj->aref,
                                    sgesvdx_obj->lda,
                                    sgesvdx_obj->vl,
                                    sgesvdx_obj->vu,
                                    sgesvdx_obj->il,
                                    sgesvdx_obj->iu,
                                    &sgesvdx_obj->nsref,
                                    sgesvdx_obj->sref,
                                    sgesvdx_obj->uref,
                                    sgesvdx_obj->ldu,
                                    sgesvdx_obj->vtref,
                                    sgesvdx_obj->ldvt,
                                    sgesvdx_obj->superbref
                                    );

    /* Compute libflame's Lapacke o/p  */
    sgesvdx_obj->info = LAPACKE_sgesvdx(  sgesvdx_obj->matrix_layout,
                                    sgesvdx_obj->jobu,
                                    sgesvdx_obj->jobvt,
                                    sgesvdx_obj->range,
                                    sgesvdx_obj->m,
                                    sgesvdx_obj->n,
                                    sgesvdx_obj->a,
                                    sgesvdx_obj->lda,
                                    sgesvdx_obj->vl,
                                    sgesvdx_obj->vu,
                                    sgesvdx_obj->il,
                                    sgesvdx_obj->iu,
                                    &sgesvdx_obj->ns,
                                    sgesvdx_obj->s,
                                    sgesvdx_obj->u,
                                    sgesvdx_obj->ldu,
                                    sgesvdx_obj->vt,
                                    sgesvdx_obj->ldvt,
                                    sgesvdx_obj->superb
                                    );


    /* Capture the Netlib, libflame o/p buffers' differences */
    sgesvdx_obj->diff_u =  computeDiff_s( (sgesvdx_obj->iu)-(sgesvdx_obj->il), 
                &sgesvdx_obj->u[sgesvdx_obj->il], &sgesvdx_obj->uref[sgesvdx_obj->il] );

    sgesvdx_obj->diff_vt =  computeDiff_s( (sgesvdx_obj->iu)-(sgesvdx_obj->il), 
                &sgesvdx_obj->vt[sgesvdx_obj->il], &sgesvdx_obj->vtref[sgesvdx_obj->il] );

    sgesvdx_obj->diff_s =  computeDiff_s( sgesvdx_obj->min_mn, 
                sgesvdx_obj->s, sgesvdx_obj->sref );

    sgesvdx_obj->diff_superb =  computeDiff_i( sgesvdx_obj->min_mn, 
                sgesvdx_obj->superb, sgesvdx_obj->superbref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n gesvdx float: \n diff_a: %f \n diff_s: %f \n \
diff_u: %f \n diff_vt: %f \n diff_superb: %f info: %d \t inforef: %d\n ",
       sgesvdx_obj->diff_a, sgesvdx_obj->diff_s, sgesvdx_obj->diff_u,
       sgesvdx_obj->diff_vt, sgesvdx_obj->diff_superb,
	   sgesvdx_obj->info, sgesvdx_obj->inforef);
#endif
}

TEST_F(sgesvdx_test, sgesvdx_1) {
    EXPECT_NEAR(0.0, sgesvdx_obj->diff_s, sgesvdx_obj->threshold);
    EXPECT_NEAR(0.0, sgesvdx_obj->diff_u, sgesvdx_obj->threshold);
    EXPECT_NEAR(0.0, sgesvdx_obj->diff_vt, sgesvdx_obj->threshold);
    EXPECT_EQ(0, sgesvdx_obj->diff_superb);}

TEST_F(sgesvdx_test, sgesvdx_2) {
    EXPECT_NEAR(0.0, sgesvdx_obj->diff_s, sgesvdx_obj->threshold);
    EXPECT_NEAR(0.0, sgesvdx_obj->diff_u, sgesvdx_obj->threshold);
    EXPECT_NEAR(0.0, sgesvdx_obj->diff_vt, sgesvdx_obj->threshold);
    EXPECT_EQ(0, sgesvdx_obj->diff_superb);}

TEST_F(sgesvdx_test, sgesvdx_3) {
    EXPECT_NEAR(0.0, sgesvdx_obj->diff_s, sgesvdx_obj->threshold);
    EXPECT_NEAR(0.0, sgesvdx_obj->diff_u, sgesvdx_obj->threshold);
    EXPECT_NEAR(0.0, sgesvdx_obj->diff_vt, sgesvdx_obj->threshold);
    EXPECT_EQ(0, sgesvdx_obj->diff_superb);}

TEST_F(sgesvdx_test, sgesvdx_4) {
    EXPECT_NEAR(0.0, sgesvdx_obj->diff_s, sgesvdx_obj->threshold);
    EXPECT_NEAR(0.0, sgesvdx_obj->diff_u, sgesvdx_obj->threshold);
    EXPECT_NEAR(0.0, sgesvdx_obj->diff_vt, sgesvdx_obj->threshold);
    EXPECT_EQ(0, sgesvdx_obj->diff_superb);
}

/* Begin gesvdx_double_common_parameters  class definition */
class gesvdx_double_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    double diff_a, diff_s, diff_u, diff_vt, diff_superb;
    int min_mn;
    double threshold;
    void *hModule, *dModule;
    
    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char jobu; // Must be 'V', or 'N'. 
    char jobvt; // Must be 'V', 'N'. 
    char range; // Must be 'A', 'V', 'I'. 

    lapack_int m; // The number of rows of the matrix A
    lapack_int n; // The number of columns of the matrix A 
    lapack_int lda; //  The leading dimension of a

    lapack_int ldu; // leading dimension of the output array u . ldu≥ max(1, m)
    lapack_int ldvt; // The leading dimension of the output array vt . ldvt≥ max(1, p) 
    double vl; // the lower and upper bounds of the interval to be searched for singular values
	double vu;
   /* Input / Output parameters */
    double* a, *aref; // contains m-by-n matrix A.
	lapack_int il; // the indices of the smallest and largest singular values to be returned.
	lapack_int iu;
    /* Output parameters */
	lapack_int ns, nsref;
    lapack_int *superb, *superbref; // the dimension of subblocks.
    // Below buffers contain the o/p orthogonal/unitary matrices
    double* u, *uref;
    double* vt, *vtref;
    double* s, *sref;
    
    /*Return Values */
    int info, inforef;

   public:
      gesvdx_double_parameters (int matrix_layout_i, char jobu, char jobvt,
                                  char range, lapack_int m, lapack_int n,
								  lapack_int il, lapack_int iu,
								  double vl, double vu);
      ~gesvdx_double_parameters ();
};

/* Constructor definition  gesvdx double_common_parameters */
gesvdx_double_parameters:: gesvdx_double_parameters (int matrix_layout_i,
                    char jobu_i, char jobvt_i, char range_i, lapack_int m_i,  
                    lapack_int n_i, lapack_int il_i, lapack_int iu_i,
								  double vl_i, double vu_i)
{
    matrix_layout = matrix_layout_i;
    jobu = jobu_i;
    jobvt = jobvt_i;
	range = range_i;
    n = n_i;
    m = m_i;
    min_mn = (m<n)?m:n;
	il = il_i;
	iu = iu_i;
	vl = vl_i;
	vu = vu_i;
    
    ldu = m;
    ldvt = n;
	ns = 0;
	nsref = 0;

    if (matrix_layout == LAPACK_ROW_MAJOR)
    {
        lda = n;
    }   else   {
        lda = m;
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
   printf(" \n gesvdx double: matrix_layout: %d   jobu: %c \t \
jobvt: %c \t range: %c \n m: %d \t n: %d \t il:%d \t iu:%d \t \
vl: %f \t vu: %f \n",
matrix_layout, jobu_i, jobvt_i,range, m_i, n_i, il, iu, vl, vu);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, m*n );
    lapacke_gtest_alloc_double_buffer_pair( &s, &sref, min_mn );
    lapacke_gtest_alloc_double_buffer_pair( &u, &uref, m*m );
    lapacke_gtest_alloc_double_buffer_pair( &vt, &vtref, n*n );
    lapacke_gtest_alloc_int_buffer_pair( &superb, &superbref, 12*min_mn );

    if( (a==NULL) || (aref==NULL) ||  \
        (s==NULL) || (sref==NULL) || \
        (u==NULL) || (uref==NULL) || \
        (vt==NULL) || (vtref==NULL) || \
        (superb==NULL) || (superbref==NULL)  ){
       EXPECT_FALSE( true) << "gesvdx_double_parameters object: malloc error. Exiting ";
       gesvdx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, m*n );
    lapacke_gtest_init_double_buffer_pair_with_constant( u, uref, m*m, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant( vt, vtref, n*n, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant( s, sref, min_mn, 0.0);
    lapacke_gtest_init_int_buffer_pair_with_constant( superb, superbref, 12*min_mn, 0);

   } /* end of Constructor  */

/* Destructor definition  'gesvdx_double_common_parameters' */
gesvdx_double_parameters :: ~gesvdx_double_parameters ()
{
   /* De-Allocate memory for the input matrices */
   gesvdx_free();
} 

//  Test fixture class definition
class dgesvdx_test  : public  ::testing::Test {
public:
   gesvdx_double_parameters  *dgesvdx_obj;
   void SetUp();  
   void TearDown () { delete dgesvdx_obj; }
};

void dgesvdx_test::SetUp()
{
    /* LAPACKE dgesvdx prototype */

	typedef int (*Fptr_NL_LAPACKE_dgesvdx) (int matrix_layout, char jobu,
		char jobvt, char range, lapack_int m, lapack_int n, double * a, 
		lapack_int lda, double vl, double vu, lapack_int il, lapack_int iu,
		lapack_int * ns, double * s, double * u, lapack_int ldu, double * vt,
		lapack_int ldvt, lapack_int * superb);

    Fptr_NL_LAPACKE_dgesvdx DGESVDX;

    dgesvdx_obj = new  gesvdx_double_parameters( svd_paramslist[idx].matrix_layout,
                                         svd_paramslist[idx].jobu_gesvdx,
                                         svd_paramslist[idx].jobvt_gesvdx,
                                         svd_paramslist[idx].range_gesvdx,
                                         svd_paramslist[idx].m_gesvj,
                                         svd_paramslist[idx].n_gesvj,
                                         svd_paramslist[idx].il,
                                         svd_paramslist[idx].iu,
                                         svd_paramslist[idx].vl,
                                         svd_paramslist[idx].vu
										 );

    idx = Circular_Increment_Index(idx);
    dgesvdx_obj->threshold = svd_paramslist[idx].svd_threshold;
    dgesvdx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgesvdx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgesvdx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgesvdx_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DGESVDX = (Fptr_NL_LAPACKE_dgesvdx)dlsym(dgesvdx_obj->hModule, "LAPACKE_dgesvdx");
    ASSERT_TRUE(DGESVDX != NULL) << "failed to ppt the Netlib LAPACKE_dgesvdx symbol";
 
 /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    dgesvdx_obj->inforef = DGESVDX( dgesvdx_obj->matrix_layout,
                                    dgesvdx_obj->jobu,
                                    dgesvdx_obj->jobvt,
                                    dgesvdx_obj->range,
                                    dgesvdx_obj->m,
                                    dgesvdx_obj->n,
                                    dgesvdx_obj->aref,
                                    dgesvdx_obj->lda,
                                    dgesvdx_obj->vl,
                                    dgesvdx_obj->vu,
                                    dgesvdx_obj->il,
                                    dgesvdx_obj->iu,
                                    &dgesvdx_obj->nsref,
                                    dgesvdx_obj->sref,
                                    dgesvdx_obj->uref,
                                    dgesvdx_obj->ldu,
                                    dgesvdx_obj->vtref,
                                    dgesvdx_obj->ldvt,
                                    dgesvdx_obj->superbref
                                    );

    /* Compute libflame's Lapacke o/p  */
    dgesvdx_obj->info = LAPACKE_dgesvdx(  dgesvdx_obj->matrix_layout,
                                    dgesvdx_obj->jobu,
                                    dgesvdx_obj->jobvt,
                                    dgesvdx_obj->range,
                                    dgesvdx_obj->m,
                                    dgesvdx_obj->n,
                                    dgesvdx_obj->a,
                                    dgesvdx_obj->lda,
                                    dgesvdx_obj->vl,
                                    dgesvdx_obj->vu,
                                    dgesvdx_obj->il,
                                    dgesvdx_obj->iu,
                                    &dgesvdx_obj->ns,
                                    dgesvdx_obj->s,
                                    dgesvdx_obj->u,
                                    dgesvdx_obj->ldu,
                                    dgesvdx_obj->vt,
                                    dgesvdx_obj->ldvt,
                                    dgesvdx_obj->superb
                                    );


    /* Capture the Netlib, libflame o/p buffers' differences */
    dgesvdx_obj->diff_u =  computeDiff_d( (dgesvdx_obj->iu)-(dgesvdx_obj->il), 
                &dgesvdx_obj->u[dgesvdx_obj->il], &dgesvdx_obj->uref[dgesvdx_obj->il] );

    dgesvdx_obj->diff_vt =  computeDiff_d( (dgesvdx_obj->iu)-(dgesvdx_obj->il), 
                &dgesvdx_obj->vt[dgesvdx_obj->il], &dgesvdx_obj->vtref[dgesvdx_obj->il] );

    dgesvdx_obj->diff_s =  computeDiff_d( dgesvdx_obj->min_mn, 
                dgesvdx_obj->s, dgesvdx_obj->sref );

    dgesvdx_obj->diff_superb =  computeDiff_i( dgesvdx_obj->min_mn, 
                dgesvdx_obj->superb, dgesvdx_obj->superbref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n gesvdx double: \n diff_a: %f \n diff_s: %f \n \
diff_u: %f \n diff_vt: %f \n diff_superb: %f info: %d \t inforef: %d\n ",
       dgesvdx_obj->diff_a, dgesvdx_obj->diff_s, dgesvdx_obj->diff_u,
       dgesvdx_obj->diff_vt, dgesvdx_obj->diff_superb,
	   dgesvdx_obj->info, dgesvdx_obj->inforef);
#endif
}

TEST_F(dgesvdx_test, dgesvdx_1) {
    EXPECT_NEAR(0.0, dgesvdx_obj->diff_s, dgesvdx_obj->threshold);
    EXPECT_NEAR(0.0, dgesvdx_obj->diff_u, dgesvdx_obj->threshold);
    EXPECT_NEAR(0.0, dgesvdx_obj->diff_vt, dgesvdx_obj->threshold);
    EXPECT_EQ(0, dgesvdx_obj->diff_superb);}

TEST_F(dgesvdx_test, dgesvdx_2) {
    EXPECT_NEAR(0.0, dgesvdx_obj->diff_s, dgesvdx_obj->threshold);
    EXPECT_NEAR(0.0, dgesvdx_obj->diff_u, dgesvdx_obj->threshold);
    EXPECT_NEAR(0.0, dgesvdx_obj->diff_vt, dgesvdx_obj->threshold);
    EXPECT_EQ(0, dgesvdx_obj->diff_superb);}

TEST_F(dgesvdx_test, dgesvdx_3) {
    EXPECT_NEAR(0.0, dgesvdx_obj->diff_s, dgesvdx_obj->threshold);
    EXPECT_NEAR(0.0, dgesvdx_obj->diff_u, dgesvdx_obj->threshold);
    EXPECT_NEAR(0.0, dgesvdx_obj->diff_vt, dgesvdx_obj->threshold);
    EXPECT_EQ(0, dgesvdx_obj->diff_superb);}

TEST_F(dgesvdx_test, dgesvdx_4) {
    EXPECT_NEAR(0.0, dgesvdx_obj->diff_s, dgesvdx_obj->threshold);
    EXPECT_NEAR(0.0, dgesvdx_obj->diff_u, dgesvdx_obj->threshold);
    EXPECT_NEAR(0.0, dgesvdx_obj->diff_vt, dgesvdx_obj->threshold);
    EXPECT_EQ(0, dgesvdx_obj->diff_superb);
}


/* Begin gesvdx_scomplex_common_parameters  class definition */
class gesvdx_scomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_a, diff_s, diff_u, diff_vt, diff_superb;
    int min_mn;
    float threshold;
    void *hModule, *dModule;
    
    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char jobu; // Must be 'V', or 'N'. 
    char jobvt; // Must be 'V', 'N'. 
    char range; // Must be 'A', 'V', 'I'. 

    lapack_int m; // The number of rows of the matrix A
    lapack_int n; // The number of columns of the matrix A 
    lapack_int lda; //  The leading dimension of a

    lapack_int ldu; // leading dimension of the output array u . ldu≥ max(1, m)
    lapack_int ldvt; // The leading dimension of the output array vt . ldvt≥ max(1, p) 
    float vl; // the lower and upper bounds of the interval to be searched for singular values
	float vu;
    /* Input / Output parameters */
    lapack_complex_float* a, *aref; // contains m-by-n matrix A.
	lapack_int il; // the indices of the smallest and largest singular values to be returned.
	lapack_int iu;
    /* Output parameters */
	lapack_int ns, nsref;
    lapack_int *superb, *superbref; // the dimension of subblocks.
    // Below buffers contain the o/p orthogonal/unitary matrices
    lapack_complex_float* u, *uref;
    lapack_complex_float* vt, *vtref;
    float* s, *sref;
    
    /*Return Values */
    int info, inforef;

   public:
      gesvdx_scomplex_parameters (int matrix_layout_i, char jobu, char jobvt,
                                  char range, lapack_int m, lapack_int n,
								  lapack_int il, lapack_int iu,
								  float vl, float vu);
      ~gesvdx_scomplex_parameters ();
};

/* Constructor definition  gesvdx lapack_complex_float_common_parameters */
gesvdx_scomplex_parameters:: gesvdx_scomplex_parameters (int matrix_layout_i,
                    char jobu_i, char jobvt_i, char range_i, lapack_int m_i,  
                    lapack_int n_i, lapack_int il_i, lapack_int iu_i,
					float vl_i, float vu_i)
{
    matrix_layout = matrix_layout_i;
    jobu = jobu_i;
    jobvt = jobvt_i;
	range = range_i;
    n = n_i;
    m = m_i;
    min_mn = (m<n)?m:n;
	il = il_i;
	iu = iu_i;
	vl = vl_i;
	vu = vu_i;
    
    ldu = m;
    ldvt = n;
	ns = 0;
	nsref = 0;

    if (matrix_layout == LAPACK_ROW_MAJOR)
    {
        lda = n;
    }   else   {
        lda = m;
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
   printf(" \n gesvdx lapack_complex_float: matrix_layout: %d   jobu: %c \t \
jobvt: %c \t range: %c \n m: %d \t n: %d \t il:%d \t iu:%d \t \
vl: %f \t vu: %f \n",
matrix_layout, jobu_i, jobvt_i,range, m_i, n_i, il, iu, vl, vu);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, m*n );
    lapacke_gtest_alloc_float_buffer_pair( &s, &sref, min_mn );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &u, &uref, m*m );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &vt, &vtref, n*n );
    lapacke_gtest_alloc_int_buffer_pair( &superb, &superbref, 12*min_mn );

    if( (a==NULL) || (aref==NULL) ||  \
        (s==NULL) || (sref==NULL) || \
        (u==NULL) || (uref==NULL) || \
        (vt==NULL) || (vtref==NULL) || \
        (superb==NULL) || (superbref==NULL)  ){
       EXPECT_FALSE( true) << "gesvdx_scomplex_parameters object: malloc error. Exiting ";
       gesvdx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, m*n );
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( u, uref, m*m, 0.0);
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( vt, vtref, n*n, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant( s, sref, min_mn, 0.0);
    lapacke_gtest_init_int_buffer_pair_with_constant( superb, superbref, 12*min_mn, 0);

   } /* end of Constructor  */

/* Destructor definition  'gesvdx_scomplex_common_parameters' */
gesvdx_scomplex_parameters :: ~gesvdx_scomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   gesvdx_free();
} 

//  Test fixture class definition
class cgesvdx_test  : public  ::testing::Test {
public:
   gesvdx_scomplex_parameters  *cgesvdx_obj;
   void SetUp();  
   void TearDown () { delete cgesvdx_obj; }
};

void cgesvdx_test::SetUp()
{
    /* LAPACKE cgesvdx prototype */

	typedef int (*Fptr_NL_LAPACKE_cgesvdx) (int matrix_layout, char jobu,
		char jobvt, char range, lapack_int m, lapack_int n,
		lapack_complex_float * a, lapack_int lda, float vl, float vu,
		lapack_int il, lapack_int iu, lapack_int * ns, float * s,
		lapack_complex_float * u, lapack_int ldu, lapack_complex_float * vt,
		lapack_int ldvt, lapack_int * superb);

    Fptr_NL_LAPACKE_cgesvdx CGESVDX;

    cgesvdx_obj = new  gesvdx_scomplex_parameters( svd_paramslist[idx].matrix_layout,
                                         svd_paramslist[idx].jobu_gesvdx,
                                         svd_paramslist[idx].jobvt_gesvdx,
                                         svd_paramslist[idx].range_gesvdx,
                                         svd_paramslist[idx].m_gesvj,
                                         svd_paramslist[idx].n_gesvj,
                                         svd_paramslist[idx].il,
                                         svd_paramslist[idx].iu,
                                         svd_paramslist[idx].vl,
                                         svd_paramslist[idx].vu
										 );

    idx = Circular_Increment_Index(idx);
    cgesvdx_obj->threshold = svd_paramslist[idx].svd_threshold;
    cgesvdx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgesvdx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgesvdx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgesvdx_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CGESVDX = (Fptr_NL_LAPACKE_cgesvdx)dlsym(cgesvdx_obj->hModule, "LAPACKE_cgesvdx");
    ASSERT_TRUE(CGESVDX != NULL) << "failed to ppt the Netlib LAPACKE_cgesvdx symbol";
 
 /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    cgesvdx_obj->inforef = CGESVDX( cgesvdx_obj->matrix_layout,
                                    cgesvdx_obj->jobu,
                                    cgesvdx_obj->jobvt,
                                    cgesvdx_obj->range,
                                    cgesvdx_obj->m,
                                    cgesvdx_obj->n,
                                    cgesvdx_obj->aref,
                                    cgesvdx_obj->lda,
                                    cgesvdx_obj->vl,
                                    cgesvdx_obj->vu,
                                    cgesvdx_obj->il,
                                    cgesvdx_obj->iu,
                                    &cgesvdx_obj->nsref,
                                    cgesvdx_obj->sref,
                                    cgesvdx_obj->uref,
                                    cgesvdx_obj->ldu,
                                    cgesvdx_obj->vtref,
                                    cgesvdx_obj->ldvt,
                                    cgesvdx_obj->superbref
                                    );

    /* Compute libflame's Lapacke o/p  */
    cgesvdx_obj->info = LAPACKE_cgesvdx(  cgesvdx_obj->matrix_layout,
                                    cgesvdx_obj->jobu,
                                    cgesvdx_obj->jobvt,
                                    cgesvdx_obj->range,
                                    cgesvdx_obj->m,
                                    cgesvdx_obj->n,
                                    cgesvdx_obj->a,
                                    cgesvdx_obj->lda,
                                    cgesvdx_obj->vl,
                                    cgesvdx_obj->vu,
                                    cgesvdx_obj->il,
                                    cgesvdx_obj->iu,
                                    &cgesvdx_obj->ns,
                                    cgesvdx_obj->s,
                                    cgesvdx_obj->u,
                                    cgesvdx_obj->ldu,
                                    cgesvdx_obj->vt,
                                    cgesvdx_obj->ldvt,
                                    cgesvdx_obj->superb
                                    );


    /* Capture the Netlib, libflame o/p buffers' differences */
    cgesvdx_obj->diff_u =  computeDiff_c( (cgesvdx_obj->iu)-(cgesvdx_obj->il), 
                &cgesvdx_obj->u[cgesvdx_obj->il], &cgesvdx_obj->uref[cgesvdx_obj->il] );

    cgesvdx_obj->diff_vt =  computeDiff_c( (cgesvdx_obj->iu)-(cgesvdx_obj->il), 
                &cgesvdx_obj->vt[cgesvdx_obj->il], &cgesvdx_obj->vtref[cgesvdx_obj->il] );

    cgesvdx_obj->diff_s =  computeDiff_s( cgesvdx_obj->min_mn, 
                cgesvdx_obj->s, cgesvdx_obj->sref );

    cgesvdx_obj->diff_superb =  computeDiff_i( cgesvdx_obj->min_mn, 
                cgesvdx_obj->superb, cgesvdx_obj->superbref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n gesvdx lapack_complex_float: \n diff_a: %f \n diff_s: %f \n \
diff_u: %f \n diff_vt: %f \n diff_superb: %f info: %d \t inforef: %d\n ",
       cgesvdx_obj->diff_a, cgesvdx_obj->diff_s, cgesvdx_obj->diff_u,
       cgesvdx_obj->diff_vt, cgesvdx_obj->diff_superb,
	   cgesvdx_obj->info, cgesvdx_obj->inforef);
#endif
}

TEST_F(cgesvdx_test, cgesvdx_1) {
    EXPECT_NEAR(0.0, cgesvdx_obj->diff_s, cgesvdx_obj->threshold);
    EXPECT_NEAR(0.0, cgesvdx_obj->diff_u, cgesvdx_obj->threshold);
    EXPECT_NEAR(0.0, cgesvdx_obj->diff_vt, cgesvdx_obj->threshold);
    EXPECT_EQ(0, cgesvdx_obj->diff_superb);}

TEST_F(cgesvdx_test, cgesvdx_2) {
    EXPECT_NEAR(0.0, cgesvdx_obj->diff_s, cgesvdx_obj->threshold);
    EXPECT_NEAR(0.0, cgesvdx_obj->diff_u, cgesvdx_obj->threshold);
    EXPECT_NEAR(0.0, cgesvdx_obj->diff_vt, cgesvdx_obj->threshold);
    EXPECT_EQ(0, cgesvdx_obj->diff_superb);}

TEST_F(cgesvdx_test, cgesvdx_3) {
    EXPECT_NEAR(0.0, cgesvdx_obj->diff_s, cgesvdx_obj->threshold);
    EXPECT_NEAR(0.0, cgesvdx_obj->diff_u, cgesvdx_obj->threshold);
    EXPECT_NEAR(0.0, cgesvdx_obj->diff_vt, cgesvdx_obj->threshold);
    EXPECT_EQ(0, cgesvdx_obj->diff_superb);}

TEST_F(cgesvdx_test, cgesvdx_4) {
    EXPECT_NEAR(0.0, cgesvdx_obj->diff_s, cgesvdx_obj->threshold);
    EXPECT_NEAR(0.0, cgesvdx_obj->diff_u, cgesvdx_obj->threshold);
    EXPECT_NEAR(0.0, cgesvdx_obj->diff_vt, cgesvdx_obj->threshold);
    EXPECT_EQ(0, cgesvdx_obj->diff_superb);
}

/* Begin gesvdx_dcomplex_common_parameters  class definition */
class gesvdx_dcomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    double diff_a, diff_s, diff_u, diff_vt, diff_superb;
    int min_mn;
    double threshold;
    void *hModule, *dModule;
    
    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char jobu; // Must be 'V', or 'N'. 
    char jobvt; // Must be 'V', 'N'. 
    char range; // Must be 'A', 'V', 'I'. 

    lapack_int m; // The number of rows of the matrix A
    lapack_int n; // The number of columns of the matrix A 
    lapack_int lda; //  The leading dimension of a

    lapack_int ldu; // leading dimension of the output array u . ldu≥ max(1, m)
    lapack_int ldvt; // The leading dimension of the output array vt . ldvt≥ max(1, p) 
    double vl; // the lower and upper bounds of the interval to be searched for singular values
	double vu;
    /* Input / Output parameters */
    lapack_complex_double* a, *aref; // contains m-by-n matrix A.
	lapack_int il; // the indices of the smallest and largest singular values to be returned.
	lapack_int iu;
    /* Output parameters */
	lapack_int ns, nsref;
    lapack_int *superb, *superbref; // the dimension of subblocks.
    // Below buffers contain the o/p orthogonal/unitary matrices
    lapack_complex_double* u, *uref;
    lapack_complex_double* vt, *vtref;
    double* s, *sref;
    
    /*Return Values */
    int info, inforef;

   public:
      gesvdx_dcomplex_parameters (int matrix_layout_i, char jobu, char jobvt,
                                  char range, lapack_int m, lapack_int n,
								  lapack_int il, lapack_int iu,
								  double vl, double vu);
      ~gesvdx_dcomplex_parameters ();
};

/* Constructor definition  gesvdx lapack_complex_double_common_parameters */
gesvdx_dcomplex_parameters:: gesvdx_dcomplex_parameters (int matrix_layout_i,
                    char jobu_i, char jobvt_i, char range_i, lapack_int m_i,  
                    lapack_int n_i, lapack_int il_i, lapack_int iu_i,
					double vl_i, double vu_i)
{
    matrix_layout = matrix_layout_i;
    jobu = jobu_i;
    jobvt = jobvt_i;
	range = range_i;
    n = n_i;
    m = m_i;
    min_mn = (m<n)?m:n;
	il = il_i;
	iu = iu_i;
	vl = vl_i;
	vu = vu_i;
    
    ldu = m;
    ldvt = n;
	ns = 0;
	nsref = 0;

    if (matrix_layout == LAPACK_ROW_MAJOR)
    {
        lda = n;
    }   else   {
        lda = m;
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
   printf(" \n gesvdx lapack_complex_double: matrix_layout: %d   jobu: %c \t \
jobvt: %c \t range: %c \n m: %d \t n: %d \t il:%d \t iu:%d \t \
vl: %f \t vu: %f \n",
matrix_layout, jobu_i, jobvt_i,range, m_i, n_i, il, iu, vl, vu);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, m*n );
    lapacke_gtest_alloc_double_buffer_pair( &s, &sref, min_mn );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &u, &uref, m*m );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &vt, &vtref, n*n );
    lapacke_gtest_alloc_int_buffer_pair( &superb, &superbref, 12*min_mn );

    if( (a==NULL) || (aref==NULL) ||  \
        (s==NULL) || (sref==NULL) || \
        (u==NULL) || (uref==NULL) || \
        (vt==NULL) || (vtref==NULL) || \
        (superb==NULL) || (superbref==NULL)  ){
       EXPECT_FALSE( true) << "gesvdx_dcomplex_parameters object: malloc error. Exiting ";
       gesvdx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, m*n );
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( u, uref, m*m, 0.0);
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( vt, vtref, n*n, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant( s, sref, min_mn, 0.0);
    lapacke_gtest_init_int_buffer_pair_with_constant( superb, superbref, 12*min_mn, 0);

   } /* end of Constructor  */

/* Destructor definition  'gesvdx_dcomplex_common_parameters' */
gesvdx_dcomplex_parameters :: ~gesvdx_dcomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   gesvdx_free();
} 

//  Test fixture class definition
class zgesvdx_test  : public  ::testing::Test {
public:
   gesvdx_dcomplex_parameters  *zgesvdx_obj;
   void SetUp();  
   void TearDown () { delete zgesvdx_obj; }
};

void zgesvdx_test::SetUp()
{
    /* LAPACKE zgesvdx prototype */

	typedef int (*Fptr_NL_LAPACKE_zgesvdx) (int matrix_layout, char jobu,
		char jobvt, char range, lapack_int m, lapack_int n,
		lapack_complex_double * a, lapack_int lda, double vl, double vu,
		lapack_int il, lapack_int iu, lapack_int * ns, double * s,
		lapack_complex_double * u, lapack_int ldu, lapack_complex_double * vt,
		lapack_int ldvt, lapack_int * superb);

    Fptr_NL_LAPACKE_zgesvdx ZGESVDX;

    zgesvdx_obj = new  gesvdx_dcomplex_parameters( svd_paramslist[idx].matrix_layout,
                                         svd_paramslist[idx].jobu_gesvdx,
                                         svd_paramslist[idx].jobvt_gesvdx,
                                         svd_paramslist[idx].range_gesvdx,
                                         svd_paramslist[idx].m_gesvj,
                                         svd_paramslist[idx].n_gesvj,
                                         svd_paramslist[idx].il,
                                         svd_paramslist[idx].iu,
                                         svd_paramslist[idx].vl,
                                         svd_paramslist[idx].vu
										 );

    idx = Circular_Increment_Index(idx);
    zgesvdx_obj->threshold = svd_paramslist[idx].svd_threshold;
    zgesvdx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgesvdx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgesvdx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgesvdx_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZGESVDX = (Fptr_NL_LAPACKE_zgesvdx)dlsym(zgesvdx_obj->hModule, "LAPACKE_zgesvdx");
    ASSERT_TRUE(ZGESVDX != NULL) << "failed to ppt the Netlib LAPACKE_zgesvdx symbol";
 
 /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    zgesvdx_obj->inforef = ZGESVDX( zgesvdx_obj->matrix_layout,
                                    zgesvdx_obj->jobu,
                                    zgesvdx_obj->jobvt,
                                    zgesvdx_obj->range,
                                    zgesvdx_obj->m,
                                    zgesvdx_obj->n,
                                    zgesvdx_obj->aref,
                                    zgesvdx_obj->lda,
                                    zgesvdx_obj->vl,
                                    zgesvdx_obj->vu,
                                    zgesvdx_obj->il,
                                    zgesvdx_obj->iu,
                                    &zgesvdx_obj->nsref,
                                    zgesvdx_obj->sref,
                                    zgesvdx_obj->uref,
                                    zgesvdx_obj->ldu,
                                    zgesvdx_obj->vtref,
                                    zgesvdx_obj->ldvt,
                                    zgesvdx_obj->superbref
                                    );

    /* Compute libflame's Lapacke o/p  */
    zgesvdx_obj->info = LAPACKE_zgesvdx(  zgesvdx_obj->matrix_layout,
                                    zgesvdx_obj->jobu,
                                    zgesvdx_obj->jobvt,
                                    zgesvdx_obj->range,
                                    zgesvdx_obj->m,
                                    zgesvdx_obj->n,
                                    zgesvdx_obj->a,
                                    zgesvdx_obj->lda,
                                    zgesvdx_obj->vl,
                                    zgesvdx_obj->vu,
                                    zgesvdx_obj->il,
                                    zgesvdx_obj->iu,
                                    &zgesvdx_obj->ns,
                                    zgesvdx_obj->s,
                                    zgesvdx_obj->u,
                                    zgesvdx_obj->ldu,
                                    zgesvdx_obj->vt,
                                    zgesvdx_obj->ldvt,
                                    zgesvdx_obj->superb
                                    );


    /* Capture the Netlib, libflame o/p buffers' differences */
    zgesvdx_obj->diff_u =  computeDiff_z( (zgesvdx_obj->iu)-(zgesvdx_obj->il), 
                &zgesvdx_obj->u[zgesvdx_obj->il], &zgesvdx_obj->uref[zgesvdx_obj->il] );

    zgesvdx_obj->diff_vt =  computeDiff_z( (zgesvdx_obj->iu)-(zgesvdx_obj->il), 
                &zgesvdx_obj->vt[zgesvdx_obj->il], &zgesvdx_obj->vtref[zgesvdx_obj->il] );

    zgesvdx_obj->diff_s =  computeDiff_d( zgesvdx_obj->min_mn, 
                zgesvdx_obj->s, zgesvdx_obj->sref );

    zgesvdx_obj->diff_superb =  computeDiff_i( zgesvdx_obj->min_mn, 
                zgesvdx_obj->superb, zgesvdx_obj->superbref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n gesvdx lapack_complex_double: \n diff_a: %f \n diff_s: %f \n \
diff_u: %f \n diff_vt: %f \n diff_superb: %f info: %d \t inforef: %d\n ",
       zgesvdx_obj->diff_a, zgesvdx_obj->diff_s, zgesvdx_obj->diff_u,
       zgesvdx_obj->diff_vt, zgesvdx_obj->diff_superb,
	   zgesvdx_obj->info, zgesvdx_obj->inforef);
#endif
}

TEST_F(zgesvdx_test, zgesvdx_1) {
    EXPECT_NEAR(0.0, zgesvdx_obj->diff_s, zgesvdx_obj->threshold);
    EXPECT_NEAR(0.0, zgesvdx_obj->diff_u, zgesvdx_obj->threshold);
    EXPECT_NEAR(0.0, zgesvdx_obj->diff_vt, zgesvdx_obj->threshold);
    EXPECT_EQ(0, zgesvdx_obj->diff_superb);}

TEST_F(zgesvdx_test, zgesvdx_2) {
    EXPECT_NEAR(0.0, zgesvdx_obj->diff_s, zgesvdx_obj->threshold);
    EXPECT_NEAR(0.0, zgesvdx_obj->diff_u, zgesvdx_obj->threshold);
    EXPECT_NEAR(0.0, zgesvdx_obj->diff_vt, zgesvdx_obj->threshold);
    EXPECT_EQ(0, zgesvdx_obj->diff_superb);}

TEST_F(zgesvdx_test, zgesvdx_3) {
    EXPECT_NEAR(0.0, zgesvdx_obj->diff_s, zgesvdx_obj->threshold);
    EXPECT_NEAR(0.0, zgesvdx_obj->diff_u, zgesvdx_obj->threshold);
    EXPECT_NEAR(0.0, zgesvdx_obj->diff_vt, zgesvdx_obj->threshold);
    EXPECT_EQ(0, zgesvdx_obj->diff_superb);}

TEST_F(zgesvdx_test, zgesvdx_4) {
    EXPECT_NEAR(0.0, zgesvdx_obj->diff_s, zgesvdx_obj->threshold);
    EXPECT_NEAR(0.0, zgesvdx_obj->diff_u, zgesvdx_obj->threshold);
    EXPECT_NEAR(0.0, zgesvdx_obj->diff_vt, zgesvdx_obj->threshold);
    EXPECT_EQ(0, zgesvdx_obj->diff_superb);
}
