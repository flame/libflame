#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"




#define hpgst_free() \
if (ap!=NULL)    free(ap); \
if (apref!=NULL) free(apref);\
if (bp!=NULL)  free(bp);\
if (bpref!=NULL) free(bpref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin scomplex_common_parameters  class definition */
class hpgst_scomplex_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;	
	lapack_complex_float* ap;
	lapack_int itype;
	char uplo;
	/*Output Parameter*/
	lapack_complex_float* bp;	
	lapack_complex_float *apref, *bpref;	
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      hpgst_scomplex_parameters (int matrix_layout, lapack_int itype, char uplo, lapack_int n);
      ~hpgst_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
hpgst_scomplex_parameters:: hpgst_scomplex_parameters (int matrix_layout_i, lapack_int itype_i, char uplo_i, lapack_int n_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;	
	itype = itype_i;
	uplo = uplo_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n hpgst scomplex: matrix_layout = %d, uplo:%c, n: %d lda: %d \n", matrix_layout, uplo, n, lda);
	#endif

	bufsize = (n*(n+1)/2);
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&ap, &apref, bufsize);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&bp, &bpref, bufsize);
	
	if ((ap==NULL) || (apref==NULL) || \
		(bp==NULL) || (bpref==NULL)) {
		EXPECT_FALSE( true) << "hpgst_float_parameters object: malloc error.";
		hpgst_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand(ap, apref, bufsize);
	lapacke_gtest_init_scomplex_buffer_pair_rand(bp, bpref, bufsize);	

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
hpgst_scomplex_parameters :: ~hpgst_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" hpgst_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   hpgst_free();

}
/*  Test fixture class definition */
class chpgst_test  : public  ::testing::Test {
public:
   hpgst_scomplex_parameters  *chpgst_obj;
   void SetUp();
   void TearDown () { delete chpgst_obj; }
};

void chpgst_test::SetUp(){

    /* LAPACKE chpgst prototype */
    typedef int (*Fptr_NL_LAPACKE_chpgst) (int matrix_layout, lapack_int itype, char uplo, \
										lapack_int n,lapack_complex_float *ap, lapack_complex_float* bp);

    Fptr_NL_LAPACKE_chpgst chpgst;

    chpgst_obj = new hpgst_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].itype,
                           eig_paramslist[idx].uplo,
                           eig_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    chpgst_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    chpgst_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(chpgst_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(chpgst_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    chpgst = (Fptr_NL_LAPACKE_chpgst)dlsym(chpgst_obj->hModule, "LAPACKE_chpgst");
    ASSERT_TRUE(chpgst != NULL) << "failed to get the Netlib LAPACKE_chpgst symbol";
    

    chpgst_obj->inforef = chpgst( chpgst_obj->matrix_layout, chpgst_obj->itype, chpgst_obj->uplo,
								chpgst_obj->n, chpgst_obj->apref, chpgst_obj->bpref);

    /* Compute libflame's Lapacke o/p  */
    chpgst_obj->info = LAPACKE_chpgst( chpgst_obj->matrix_layout, chpgst_obj->itype, chpgst_obj->uplo,
										chpgst_obj->n,chpgst_obj->ap,chpgst_obj->bp);

    if( chpgst_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_chpgst is wrong\n", chpgst_obj->info );
    }
    if( chpgst_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chpgst is wrong\n", 
        chpgst_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    chpgst_obj->diff =  computeDiff_c( chpgst_obj->bufsize, 
                chpgst_obj->ap, chpgst_obj->apref );

}

TEST_F(chpgst_test, chpgst1) {
    EXPECT_NEAR(0.0, chpgst_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chpgst_test, chpgst2) {
    EXPECT_NEAR(0.0, chpgst_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chpgst_test, chpgst3) {
    EXPECT_NEAR(0.0, chpgst_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chpgst_test, chpgst4) {
    EXPECT_NEAR(0.0, chpgst_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class hpgst_dcomplex_parameters{
	
	 public:
	int bufsize;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;	
	lapack_complex_double* ap;
	lapack_int itype;
	char uplo;
	/*Output Parameter*/
	lapack_complex_double* bp;	
	lapack_complex_double *apref, *bpref;	
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      hpgst_dcomplex_parameters (int matrix_layout, lapack_int itype, char uplo, lapack_int n);
      ~hpgst_dcomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
hpgst_dcomplex_parameters:: hpgst_dcomplex_parameters (int matrix_layout_i, lapack_int itype_i, char uplo_i, lapack_int n_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;	
	itype = itype_i;
	uplo = uplo_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n hpgst scomplex:  n: %d lda: %d \n", n, lda);
	#endif

	bufsize = (n*(n+1)/2);
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&ap, &apref, bufsize);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&bp, &bpref, bufsize);	
	
	if ((ap==NULL) || (apref==NULL) || \
		(bp==NULL) || (bpref==NULL)) {
		EXPECT_FALSE( true) << "hpgst_float_parameters object: malloc error.";
		hpgst_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( ap, apref, bufsize);
	lapacke_gtest_init_dcomplex_buffer_pair_rand( bp, bpref, bufsize);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
hpgst_dcomplex_parameters :: ~hpgst_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" hpgst_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   hpgst_free();

}
/*  Test fixture class definition */
class zhpgst_test  : public  ::testing::Test {
public:
   hpgst_dcomplex_parameters  *zhpgst_obj;
   void SetUp();
   void TearDown () { delete zhpgst_obj; }
};

void zhpgst_test::SetUp(){

    /* LAPACKE zhpgst prototype */
    typedef int (*Fptr_NL_LAPACKE_zhpgst) (int matrix_layout, lapack_int itype, char uplo, lapack_int n, lapack_complex_double *ap,\
										 lapack_complex_double* bp);

    Fptr_NL_LAPACKE_zhpgst zhpgst;

    zhpgst_obj = new hpgst_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].itype,
                           eig_paramslist[idx].uplo,
                           eig_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    zhpgst_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zhpgst_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zhpgst_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zhpgst_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zhpgst = (Fptr_NL_LAPACKE_zhpgst)dlsym(zhpgst_obj->hModule, "LAPACKE_zhpgst");
    ASSERT_TRUE(zhpgst != NULL) << "failed to get the Netlib LAPACKE_zhpgst symbol";
    

    zhpgst_obj->inforef = zhpgst( zhpgst_obj->matrix_layout, zhpgst_obj->itype,  zhpgst_obj->uplo,
								zhpgst_obj->n, zhpgst_obj->apref, zhpgst_obj->bpref);

    /* Compute libflame's Lapacke o/p  */
    zhpgst_obj->info = LAPACKE_zhpgst( zhpgst_obj->matrix_layout, zhpgst_obj->itype, zhpgst_obj->uplo,
										zhpgst_obj->n,zhpgst_obj->ap, zhpgst_obj->bp);

    if( zhpgst_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zhpgst is wrong\n", zhpgst_obj->info );
    }
    if( zhpgst_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhpgst is wrong\n", 
        zhpgst_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zhpgst_obj->diff =  computeDiff_z( zhpgst_obj->bufsize, 
                zhpgst_obj->ap, zhpgst_obj->apref );

}

TEST_F(zhpgst_test, zhpgst1) {
    EXPECT_NEAR(0.0, zhpgst_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhpgst_test, zhpgst2) {
    EXPECT_NEAR(0.0, zhpgst_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhpgst_test, zhpgst3) {
    EXPECT_NEAR(0.0, zhpgst_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhpgst_test, zhpgst4) {
    EXPECT_NEAR(0.0, zhpgst_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

	