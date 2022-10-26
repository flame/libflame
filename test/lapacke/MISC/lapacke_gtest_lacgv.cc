#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"


#define lacgv_free() \
if (x!=NULL)    free(x); \
if (xref!=NULL) free(xref);\
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin scomplex_common_parameters  class definition */
class lacgv_scomplex_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	lapack_int n;	
	lapack_complex_float* x;
	lapack_int incx;
	/*Output Parameter*/	
	lapack_complex_float *xref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      lacgv_scomplex_parameters (lapack_int n, lapack_int incx);
      ~lacgv_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
lacgv_scomplex_parameters:: lacgv_scomplex_parameters (lapack_int n_i, lapack_int incx_i)
{
	n = n_i;
	incx = incx_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n lacgv scomplex:  n: %d incx: %d \n", n, incx);
	#endif
	
	bufsize =  (1 +(n-1)*(incx));
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&x, &xref, bufsize);
	if ((x==NULL) || (xref==NULL)){
		EXPECT_FALSE( true) << "lacgv_float_parameters object: malloc error.";
		lacgv_free();
		exit(0);
	}

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
lacgv_scomplex_parameters :: ~lacgv_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" lacgv_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lacgv_free();

}
/*  Test fixture class definition */
class clacgv_test  : public  ::testing::Test {
public:
   lacgv_scomplex_parameters  *clacgv_obj;
   void SetUp();
   void TearDown () { delete clacgv_obj; }
};

void clacgv_test::SetUp(){

    /* LAPACKE clacgv prototype */
    typedef int (*Fptr_NL_LAPACKE_clacgv) (lapack_int n, lapack_complex_float *x, lapack_int incx);

    Fptr_NL_LAPACKE_clacgv clacgv;

    clacgv_obj = new lacgv_scomplex_parameters ( lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    clacgv_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    clacgv_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(clacgv_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(clacgv_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    clacgv = (Fptr_NL_LAPACKE_clacgv)dlsym(clacgv_obj->hModule, "LAPACKE_clacgv");
    ASSERT_TRUE(clacgv != NULL) << "failed to get the Netlib LAPACKE_clacgv symbol";
    

    clacgv_obj->inforef = clacgv( clacgv_obj->n, clacgv_obj->xref,
								clacgv_obj->incx);

    /* Compute libflame's Lapacke o/p  */
    clacgv_obj->info = LAPACKE_clacgv(clacgv_obj->n, clacgv_obj->x, 
										clacgv_obj->incx);

    if( clacgv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_clacgv is wrong\n", clacgv_obj->info );
    }
    if( clacgv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_clacgv is wrong\n", 
        clacgv_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    clacgv_obj->diff =  computeDiff_c( clacgv_obj->bufsize, 
                clacgv_obj->x, clacgv_obj->xref );

}

TEST_F(clacgv_test, clacgv1) {
    EXPECT_NEAR(0.0, clacgv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(clacgv_test, clacgv2) {
    EXPECT_NEAR(0.0, clacgv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(clacgv_test, clacgv3) {
    EXPECT_NEAR(0.0, clacgv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(clacgv_test, clacgv4) {
    EXPECT_NEAR(0.0, clacgv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}
/* Begin dcomplex_common_parameters  class definition */
class lacgv_dcomplex_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	lapack_int n;	
	lapack_complex_double* x;
	lapack_int incx;
	/*Output Parameter*/	
	lapack_complex_double *xref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      lacgv_dcomplex_parameters (lapack_int n, lapack_int incx);
      ~lacgv_dcomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
lacgv_dcomplex_parameters:: lacgv_dcomplex_parameters (lapack_int n_i, lapack_int incx_i)
{
	n = n_i;
	incx = incx_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n lacgv dcomplex:  n: %d incx: %d \n", n, incx);
	#endif
	
	bufsize =  (1 +(n-1)*(incx));
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&x, &xref, bufsize);
	if ((x==NULL) || (xref==NULL)){
		EXPECT_FALSE( true) << "lacgv_float_parameters object: malloc error.";
		lacgv_free();
		exit(0);
	}

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
lacgv_dcomplex_parameters :: ~lacgv_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" lacgv_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lacgv_free();

}
/*  Test fixture class definition */
class zlacgv_test  : public  ::testing::Test {
public:
   lacgv_dcomplex_parameters  *zlacgv_obj;
   void SetUp();
   void TearDown () { delete zlacgv_obj; }
};

void zlacgv_test::SetUp(){

    /* LAPACKE zlacgv prototype */
    typedef int (*Fptr_NL_LAPACKE_zlacgv) (lapack_int n, lapack_complex_double *x, lapack_int incx);

    Fptr_NL_LAPACKE_zlacgv zlacgv;

    zlacgv_obj = new lacgv_dcomplex_parameters ( lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    zlacgv_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zlacgv_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zlacgv_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zlacgv_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zlacgv = (Fptr_NL_LAPACKE_zlacgv)dlsym(zlacgv_obj->hModule, "LAPACKE_zlacgv");
    ASSERT_TRUE(zlacgv != NULL) << "failed to get the Netlib LAPACKE_zlacgv symbol";
    

    zlacgv_obj->inforef = zlacgv( zlacgv_obj->n, zlacgv_obj->xref,
								zlacgv_obj->incx);

    /* Compute libflame's Lapacke o/p  */
    zlacgv_obj->info = LAPACKE_zlacgv(zlacgv_obj->n, zlacgv_obj->x, 
										zlacgv_obj->incx);

    if( zlacgv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zlacgv is wrong\n", zlacgv_obj->info );
    }
    if( zlacgv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zlacgv is wrong\n", 
        zlacgv_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zlacgv_obj->diff =  computeDiff_z( zlacgv_obj->bufsize, 
                zlacgv_obj->x, zlacgv_obj->xref );

}

TEST_F(zlacgv_test, zlacgv1) {
    EXPECT_NEAR(0.0, zlacgv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlacgv_test, zlacgv2) {
    EXPECT_NEAR(0.0, zlacgv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlacgv_test, zlacgv3) {
    EXPECT_NEAR(0.0, zlacgv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlacgv_test, zlacgv4) {
    EXPECT_NEAR(0.0, zlacgv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}
