#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"
#define LAPACKE_TEST_VERBOSE  (1)
#define bdsvdx_free() \
    if (d!=NULL)        free(d); \
    if (dref!=NULL)     free(dref); \
    if (s!=NULL)        free(s); \
    if (sref!=NULL)     free(sref); \
    if (e!=NULL)        free(e); \
    if (eref!=NULL)     free(eref); \
    if (z!=NULL)        free(z); \
    if (zref!=NULL)     free(zref); \
    if (superb!=NULL)        free(superb); \
    if (superbref!=NULL)     free(superbref)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin bdsvdx_float_common_parameters  class definition */
class bdsvdx_float_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_s, diff_z;
	int diff_superb;
    float threshold;
    void *hModule, *dModule;
    
    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char jobz; // Must be 'V', or 'N'. 
    char uplo; // Must be 'L', 'U'. 
    char range; // Must be 'A', 'V', 'I'.
    lapack_int n; // The order of bidiagonal matrix A 
    lapack_int ldz; //  The leading dimension of Z	
	lapack_int il; // the indices of the smallest and largest singular values to be returned.
	lapack_int iu;
    float vl; // the lower and upper bounds of the interval to be searched for singular values
	float vu;
    float* d, *dref; // The n diagonal elements of the bidiagonal matrix B.
    float* e, *eref; // The (n - 1) superdiagonal elements of the bidiagonal matrix B
	
    /* Output parameters */
	lapack_int ns, nsref;
    lapack_int *superb, *superbref; // contains the indices of the eigenvectors that failed to converge.
    float* s, *sref; // The first ns elements contain the selected singular values in ascending order.
    float* z, *zref; // contain the singular vectors of the matrix B corresponding to the selected singular values.
    
    /*Return Values */
    int info, inforef;

   public:
      bdsvdx_float_parameters (int matrix_layout_i, char jobz, char uplo,
                                  char range, lapack_int n,
								  lapack_int il, lapack_int iu,
								  float vl, float vu);
      ~bdsvdx_float_parameters ();
};

/* Constructor definition  bdsvdx float_common_parameters */
bdsvdx_float_parameters:: bdsvdx_float_parameters (int matrix_layout_i,
                    char jobz_i, char uplo_i, char range_i, lapack_int n_i, 
					lapack_int il_i, lapack_int iu_i, float vl_i, float vu_i)
{
    matrix_layout = matrix_layout_i;
    jobz = jobz_i;
    uplo = uplo_i;
	range = range_i;
    n = n_i;
	il = il_i;
	iu = iu_i;
	vl = vl_i;
	vu = vu_i;
    
    ldz = 2*n;
	ns = 0;
	nsref = 0;

    // Initialize with default values
    hModule = NULL;
    dModule = NULL;
    diff_s = 0.0;
    diff_z = 0.0;
    diff_superb = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n bdsvdx float: matrix_layout: %d   jobz: %c \t \
uplo: %c \t range: %c \n n: %d \t il:%d \t iu:%d \t \
vl: %f \t vu: %f \n",
matrix_layout, jobz_i, uplo_i,range, n_i, il, iu, vl, vu);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &d, &dref, n );
    lapacke_gtest_alloc_float_buffer_pair( &s, &sref, n );
    lapacke_gtest_alloc_float_buffer_pair( &e, &eref, n-1);
    lapacke_gtest_alloc_float_buffer_pair( &z, &zref, (n*(n+2)*2) );
    lapacke_gtest_alloc_int_buffer_pair( &superb, &superbref, 12*n );

    if( (d==NULL) || (dref==NULL) ||  \
        (s==NULL) || (sref==NULL) || \
        (e==NULL) || (eref==NULL) || \
        (z==NULL) || (zref==NULL) || \
        (superb==NULL) || (superbref==NULL)  ){
       EXPECT_FALSE( true) << "bdsvdx_float_parameters object: malloc error. Exiting ";
       bdsvdx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( d, dref, n );
    lapacke_gtest_init_float_buffer_pair_rand( e, eref, n-1);
    lapacke_gtest_init_float_buffer_pair_with_constant( z, zref, (n*(n+2)*2), 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant( s, sref, n, 0.0);
	//superb = superbref = NULL;
    lapacke_gtest_init_int_buffer_pair_with_constant( superb, superbref, 12*n, 0);

   } /* end of Constructor  */

/* Destructor definition  'bdsvdx_float_common_parameters' */
bdsvdx_float_parameters :: ~bdsvdx_float_parameters ()
{
   /* De-Allocate memory for the input matrices */
   //bdsvdx_free();
    if (d!=NULL)        free(d);
    if (dref!=NULL)     free(dref);
    if (s!=NULL)        free(s);
    if (sref!=NULL)     free(sref);
    if (e!=NULL)        free(e);
    if (eref!=NULL)     free(eref);
    if (z!=NULL)        free(z);
    if (zref!=NULL)     free(zref);
    if (superb!=NULL)        free(superb);
    if (superbref!=NULL)     free(superbref);

   
} 

//  Test fixture class definition
class sbdsvdx_test  : public  ::testing::Test {
public:
   bdsvdx_float_parameters  *sbdsvdx_obj;
   void SetUp();  
   void TearDown () { delete sbdsvdx_obj; }
};

void sbdsvdx_test::SetUp()
{
    /* LAPACKE sbdsvdx prototype */

	typedef int (*Fptr_NL_LAPACKE_sbdsvdx) (int matrix_layout,
		char uplo, char jobz, char range, lapack_int n,
		float * d, float * e, float vl, float vu, lapack_int il,
		lapack_int iu, lapack_int * ns, float * s, float * z,
		lapack_int ldz, lapack_int * superb);

    Fptr_NL_LAPACKE_sbdsvdx SBDSVDX;

    sbdsvdx_obj = new  bdsvdx_float_parameters( svd_paramslist[idx].matrix_layout,
                                         svd_paramslist[idx].jobu_gesvdx,
                                         eig_paramslist[idx].uplo,
                                         svd_paramslist[idx].range_gesvdx,
                                         svd_paramslist[idx].n_gesvj,
                                         svd_paramslist[idx].il,
                                         svd_paramslist[idx].iu,
                                         svd_paramslist[idx].vl,
                                         svd_paramslist[idx].vu
										 );

    idx = Circular_Increment_Index(idx);
    sbdsvdx_obj->threshold = svd_paramslist[idx].svd_threshold;
    sbdsvdx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sbdsvdx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sbdsvdx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sbdsvdx_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SBDSVDX = (Fptr_NL_LAPACKE_sbdsvdx)dlsym(sbdsvdx_obj->hModule, "LAPACKE_sbdsvdx");
    ASSERT_TRUE(SBDSVDX != NULL) << "failed to ppt the Netlib LAPACKE_sbdsvdx symbol";
 

    /* Compute libflame's Lapacke o/p  */
    sbdsvdx_obj->info = LAPACKE_sbdsvdx(  sbdsvdx_obj->matrix_layout,
                                    sbdsvdx_obj->uplo,
                                    sbdsvdx_obj->jobz,
                                    sbdsvdx_obj->range,
                                    sbdsvdx_obj->n,
                                    sbdsvdx_obj->d,
                                    sbdsvdx_obj->e,
                                    sbdsvdx_obj->vl,
                                    sbdsvdx_obj->vu,
                                    sbdsvdx_obj->il,
                                    sbdsvdx_obj->iu,
                                    &sbdsvdx_obj->ns,
                                    sbdsvdx_obj->s,
                                    sbdsvdx_obj->z,
                                    sbdsvdx_obj->ldz,
                                    sbdsvdx_obj->superb
                                    );

 /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    sbdsvdx_obj->inforef = SBDSVDX( sbdsvdx_obj->matrix_layout,
                                    sbdsvdx_obj->uplo,
                                    sbdsvdx_obj->jobz,
                                    sbdsvdx_obj->range,
                                    sbdsvdx_obj->n,
                                    sbdsvdx_obj->dref,
                                    sbdsvdx_obj->eref,
                                    sbdsvdx_obj->vl,
                                    sbdsvdx_obj->vu,
                                    sbdsvdx_obj->il,
                                    sbdsvdx_obj->iu,
                                    &sbdsvdx_obj->nsref,
                                    sbdsvdx_obj->sref,
                                    sbdsvdx_obj->zref,
                                    sbdsvdx_obj->ldz,
                                    sbdsvdx_obj->superbref
                                    );

    EXPECT_LE(sbdsvdx_obj->ns, sbdsvdx_obj->n);
    EXPECT_LE(sbdsvdx_obj->nsref, sbdsvdx_obj->n);
    if ( sbdsvdx_obj->ns <= sbdsvdx_obj->n) {
    /* Capture the Netlib, libflame o/p buffers' differences */
       sbdsvdx_obj->diff_s =  computeDiff_s( sbdsvdx_obj->ns, 
                          sbdsvdx_obj->s, sbdsvdx_obj->sref );
	   
       sbdsvdx_obj->diff_z =  computeDiff_s( sbdsvdx_obj->ns, 
                          sbdsvdx_obj->z, sbdsvdx_obj->zref );
	   
	}
#if LAPACKE_TEST_VERBOSE
   printf(" \n bdsvdx float: \n diff_s: %f \n \
diff_z: %f \n info: %d \t inforef: %d\n ",
       sbdsvdx_obj->diff_s, sbdsvdx_obj->diff_z,
	   sbdsvdx_obj->info, sbdsvdx_obj->inforef);
#endif
}

TEST_F(sbdsvdx_test, sbdsvdx_1) {
    EXPECT_NEAR(0.0, sbdsvdx_obj->diff_s, sbdsvdx_obj->threshold);
    EXPECT_NEAR(0.0, sbdsvdx_obj->diff_z, sbdsvdx_obj->threshold);
    EXPECT_EQ(sbdsvdx_obj->ns, sbdsvdx_obj->nsref);
}

TEST_F(sbdsvdx_test, sbdsvdx_2) {
    EXPECT_NEAR(0.0, sbdsvdx_obj->diff_s, sbdsvdx_obj->threshold);
    EXPECT_NEAR(0.0, sbdsvdx_obj->diff_z, sbdsvdx_obj->threshold);
    EXPECT_EQ(sbdsvdx_obj->ns, sbdsvdx_obj->nsref);
}

TEST_F(sbdsvdx_test, sbdsvdx_3) {
    EXPECT_NEAR(0.0, sbdsvdx_obj->diff_s, sbdsvdx_obj->threshold);
    EXPECT_NEAR(0.0, sbdsvdx_obj->diff_z, sbdsvdx_obj->threshold);
    EXPECT_EQ(sbdsvdx_obj->ns, sbdsvdx_obj->nsref);
}

TEST_F(sbdsvdx_test, sbdsvdx_4) {
    EXPECT_NEAR(0.0, sbdsvdx_obj->diff_s, sbdsvdx_obj->threshold);
    EXPECT_NEAR(0.0, sbdsvdx_obj->diff_z, sbdsvdx_obj->threshold);
    EXPECT_EQ(sbdsvdx_obj->ns, sbdsvdx_obj->nsref);
}

