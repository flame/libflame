#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"


#define larfg_free() \
if (x!=NULL)    free(x); \
if (xref!=NULL) free(xref);\
if (tau!=NULL)    free(tau); \
if (tauref!=NULL) free(tauref);\
if (alpha!=NULL)  free(alpha); \
if (alpharef!=NULL)    free(alpharef); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin float_common_parameters  class definition */
class larfg_float_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	lapack_int n;	
	float* x, *alpha;
	lapack_int incx;
	float* tau;
	/*Output Parameter*/	
	float *xref, *tauref, *alpharef;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      larfg_float_parameters (lapack_int n, lapack_int incx);
      ~larfg_float_parameters ();

};

/* Constructor definition  float_common_parameters */
larfg_float_parameters:: larfg_float_parameters (lapack_int n_i, lapack_int incx_i)
{
	n = n_i;
	incx = incx_i;
	//alpha = alpha_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n larfg float:  n: %d incx: %d \n", n, incx);
	#endif
	
	bufsize =  (1 +(n-2)*abs(incx));
	
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&x, &xref, bufsize);
	lapacke_gtest_alloc_float_buffer_pair(&alpha, &alpharef, bufsize);
	lapacke_gtest_alloc_float_buffer_pair(&tau, &tauref, bufsize);
	if ((x==NULL) || (xref==NULL) || \
		(alpha == NULL) || (alpharef == NULL) || \
		(tau == NULL) || (tauref == NULL)){
		EXPECT_FALSE( true) << "larfg_float_parameters object: malloc error.";
		larfg_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand( x, xref, bufsize);
	lapacke_gtest_init_float_buffer_pair_with_constant(alpha , alpharef, bufsize, 1);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
larfg_float_parameters :: ~larfg_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" larfg_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   larfg_free();

}
/*  Test fixture class definition */
class slarfg_test  : public  ::testing::Test {
public:
   larfg_float_parameters  *slarfg_obj;
   void SetUp();
   void TearDown () { delete slarfg_obj; }
};

void slarfg_test::SetUp(){

    /* LAPACKE slarfg prototype */
    typedef int (*Fptr_NL_LAPACKE_slarfg) (lapack_int n , float * alpha , float * x , lapack_int incx , float * tau);

    Fptr_NL_LAPACKE_slarfg slarfg;

    slarfg_obj = new larfg_float_parameters ( lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);
						 //  #lin_solver_paramslist[idx].kl);

    idx = Circular_Increment_Index(idx);

    slarfg_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    slarfg_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(slarfg_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(slarfg_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    slarfg = (Fptr_NL_LAPACKE_slarfg)dlsym(slarfg_obj->hModule, "LAPACKE_slarfg");
    ASSERT_TRUE(slarfg != NULL) << "failed to get the Netlib LAPACKE_slarfg symbol";
    

    slarfg_obj->inforef = slarfg( slarfg_obj->n, slarfg_obj->alpharef, slarfg_obj->xref,
								slarfg_obj->incx, slarfg_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    slarfg_obj->info = LAPACKE_slarfg(slarfg_obj->n, slarfg_obj->alpha, slarfg_obj->x, 
										slarfg_obj->incx, slarfg_obj->tau);

    if( slarfg_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_slarfg is wrong\n", slarfg_obj->info );
    }
    if( slarfg_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_slarfg is wrong\n", 
        slarfg_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    slarfg_obj->diff =  computeDiff_s( slarfg_obj->bufsize, 
                slarfg_obj->x, slarfg_obj->xref );

}

TEST_F(slarfg_test, slarfg1) {
    EXPECT_NEAR(0.0, slarfg_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(slarfg_test, slarfg2) {
    EXPECT_NEAR(0.0, slarfg_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(slarfg_test, slarfg3) {
    EXPECT_NEAR(0.0, slarfg_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(slarfg_test, slarfg4) {
    EXPECT_NEAR(0.0, slarfg_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin double_common_parameters  class definition */
class larfg_double_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	lapack_int n;	
	double* x, *alpha;
	lapack_int incx;
	double* tau;
	/*Output Parameter*/	
	double *xref, *tauref, *alpharef;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      larfg_double_parameters (lapack_int n, lapack_int incx);
      ~larfg_double_parameters ();

};
/* Constructor definition  float_common_parameters */
larfg_double_parameters:: larfg_double_parameters (lapack_int n_i, lapack_int incx_i)
{
	n = n_i;
	incx = incx_i;
	//alpha = alpha_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n larfg double:  n: %d incx: %d \n", n, incx);
	#endif
	
	bufsize =  (1 +(n-2)*abs(incx));
	
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&x, &xref, bufsize);
	lapacke_gtest_alloc_double_buffer_pair(&alpha, &alpharef, bufsize);
	lapacke_gtest_alloc_double_buffer_pair(&tau, &tauref, bufsize);
	if ((x==NULL) || (xref==NULL) || \
		(alpha == NULL) || (alpharef == NULL) || \
		(tau == NULL) || (tauref == NULL)){
		EXPECT_FALSE( true) << "larfg_double_parameters object: malloc error.";
		larfg_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand( x, xref, bufsize);
	lapacke_gtest_init_double_buffer_pair_with_constant(alpha , alpharef, bufsize, 1);

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
larfg_double_parameters :: ~larfg_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" larfg_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   larfg_free();

}
/*  Test fixture class definition */
class dlarfg_test  : public  ::testing::Test {
public:
   larfg_double_parameters  *dlarfg_obj;
   void SetUp();
   void TearDown () { delete dlarfg_obj; }
};

void dlarfg_test::SetUp(){

    /* LAPACKE dlarfg prototype */
    typedef int (*Fptr_NL_LAPACKE_dlarfg) (lapack_int n , double * alpha , double * x , lapack_int incx , double * tau);

    Fptr_NL_LAPACKE_dlarfg dlarfg;

    dlarfg_obj = new larfg_double_parameters ( lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);
						 //  #lin_solver_paramslist[idx].kl);

    idx = Circular_Increment_Index(idx);

    dlarfg_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dlarfg_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dlarfg_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dlarfg_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    dlarfg = (Fptr_NL_LAPACKE_dlarfg)dlsym(dlarfg_obj->hModule, "LAPACKE_dlarfg");
    ASSERT_TRUE(dlarfg != NULL) << "failed to get the Netlib LAPACKE_dlarfg symbol";
    

    dlarfg_obj->inforef = dlarfg( dlarfg_obj->n, dlarfg_obj->alpharef, dlarfg_obj->xref,
								dlarfg_obj->incx, dlarfg_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    dlarfg_obj->info = LAPACKE_dlarfg(dlarfg_obj->n, dlarfg_obj->alpha, dlarfg_obj->x, 
										dlarfg_obj->incx, dlarfg_obj->tau);

    if( dlarfg_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dlarfg is wrong\n", dlarfg_obj->info );
    }
    if( dlarfg_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dlarfg is wrong\n", 
        dlarfg_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dlarfg_obj->diff =  computeDiff_d( dlarfg_obj->bufsize, 
                dlarfg_obj->x, dlarfg_obj->xref );

}

TEST_F(dlarfg_test, dlarfg1) {
    EXPECT_NEAR(0.0, dlarfg_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dlarfg_test, dlarfg2) {
    EXPECT_NEAR(0.0, dlarfg_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dlarfg_test, dlarfg3) {
    EXPECT_NEAR(0.0, dlarfg_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dlarfg_test, dlarfg4) {
    EXPECT_NEAR(0.0, dlarfg_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin scomplex_common_parameters  class definition */
class larfg_scomplex_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	lapack_int n;	
	lapack_complex_float* x, *alpha;
	lapack_int incx;
	lapack_complex_float* tau;
	/*Output Parameter*/	
	lapack_complex_float *xref, *tauref, *alpharef;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      larfg_scomplex_parameters (lapack_int n, lapack_int incx);
      ~larfg_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
larfg_scomplex_parameters:: larfg_scomplex_parameters (lapack_int n_i, lapack_int incx_i)
{
	n = n_i;
	incx = incx_i;
	//alpha = alpha_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n larfg scomplex:  n: %d incx: %d \n", n, incx);
	#endif
	
	bufsize =  (1 +(n-2)*abs(incx));
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&x, &xref, bufsize);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&alpha, &alpharef, bufsize);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&tau, &tauref, bufsize);
	if ((x==NULL) || (xref==NULL) || \
		(alpha == NULL) || (alpharef == NULL) || \
		(tau == NULL) || (tauref == NULL)){
		EXPECT_FALSE( true) << "larfg_float_parameters object: malloc error.";
		larfg_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( x, xref, bufsize);
	lapacke_gtest_init_scomplex_buffer_pair_with_constant(alpha , alpharef, bufsize, 1);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
larfg_scomplex_parameters :: ~larfg_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" larfg_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   larfg_free();

}
/*  Test fixture class definition */
class clarfg_test  : public  ::testing::Test {
public:
   larfg_scomplex_parameters  *clarfg_obj;
   void SetUp();
   void TearDown () { delete clarfg_obj; }
};

void clarfg_test::SetUp(){

    /* LAPACKE clarfg prototype */
    typedef int (*Fptr_NL_LAPACKE_clarfg) (lapack_int n , lapack_complex_float * alpha , lapack_complex_float * x , lapack_int incx , lapack_complex_float * tau);

    Fptr_NL_LAPACKE_clarfg clarfg;

    clarfg_obj = new larfg_scomplex_parameters ( lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);
						 //  #lin_solver_paramslist[idx].kl);

    idx = Circular_Increment_Index(idx);

    clarfg_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    clarfg_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(clarfg_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(clarfg_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    clarfg = (Fptr_NL_LAPACKE_clarfg)dlsym(clarfg_obj->hModule, "LAPACKE_clarfg");
    ASSERT_TRUE(clarfg != NULL) << "failed to get the Netlib LAPACKE_clarfg symbol";
    

    clarfg_obj->inforef = clarfg( clarfg_obj->n, clarfg_obj->alpharef, clarfg_obj->xref,
								clarfg_obj->incx, clarfg_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    clarfg_obj->info = LAPACKE_clarfg(clarfg_obj->n, clarfg_obj->alpha, clarfg_obj->x, 
										clarfg_obj->incx, clarfg_obj->tau);

    if( clarfg_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_clarfg is wrong\n", clarfg_obj->info );
    }
    if( clarfg_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_clarfg is wrong\n", 
        clarfg_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    clarfg_obj->diff =  computeDiff_c( clarfg_obj->bufsize, 
                clarfg_obj->x, clarfg_obj->xref );

}

TEST_F(clarfg_test, clarfg1) {
    EXPECT_NEAR(0.0, clarfg_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(clarfg_test, clarfg2) {
    EXPECT_NEAR(0.0, clarfg_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(clarfg_test, clarfg3) {
    EXPECT_NEAR(0.0, clarfg_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(clarfg_test, clarfg4) {
    EXPECT_NEAR(0.0, clarfg_obj->diff, LAPACKE_GTEST_THRESHOLD);
}
/* Begin dcomplex_common_parameters  class definition */
class larfg_dcomplex_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	lapack_int n;	
	lapack_complex_double* x, *alpha;
	lapack_int incx;
	lapack_complex_double* tau;
	/*Output Parameter*/	
	lapack_complex_double *xref, *tauref, *alpharef;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      larfg_dcomplex_parameters (lapack_int n, lapack_int incx);
      ~larfg_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
larfg_dcomplex_parameters:: larfg_dcomplex_parameters (lapack_int n_i, lapack_int incx_i)
{
	n = n_i;
	incx = incx_i;
	//alpha = alpha_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n larfg dcomplex:  n: %d incx: %d \n", n, incx);
	#endif
	
	bufsize =  (1 +(n-2)*abs(incx));
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&x, &xref, bufsize);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&alpha, &alpharef, bufsize);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&tau, &tauref, bufsize);
	if ((x==NULL) || (xref==NULL) || \
		(alpha == NULL) || (alpharef == NULL) || \
		(tau == NULL) || (tauref == NULL)){
		EXPECT_FALSE( true) << "larfg_double_parameters object: malloc error.";
		larfg_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( x, xref, bufsize);
	lapacke_gtest_init_dcomplex_buffer_pair_with_constant(alpha , alpharef, bufsize, 1);

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
larfg_dcomplex_parameters :: ~larfg_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" larfg_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   larfg_free();

}
/*  Test fixture class definition */
class zlarfg_test  : public  ::testing::Test {
public:
   larfg_dcomplex_parameters  *zlarfg_obj;
   void SetUp();
   void TearDown () { delete zlarfg_obj; }
};

void zlarfg_test::SetUp(){

    /* LAPACKE zlarfg prototype */
    typedef int (*Fptr_NL_LAPACKE_zlarfg) (lapack_int n , lapack_complex_double * alpha , lapack_complex_double * x , lapack_int incx , lapack_complex_double * tau);

    Fptr_NL_LAPACKE_zlarfg zlarfg;

    zlarfg_obj = new larfg_dcomplex_parameters ( lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);
						 //  #lin_solver_paramslist[idx].kl);

    idx = Circular_Increment_Index(idx);

    zlarfg_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zlarfg_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zlarfg_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zlarfg_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zlarfg = (Fptr_NL_LAPACKE_zlarfg)dlsym(zlarfg_obj->hModule, "LAPACKE_zlarfg");
    ASSERT_TRUE(zlarfg != NULL) << "failed to get the Netlib LAPACKE_zlarfg symbol";
    

    zlarfg_obj->inforef = zlarfg( zlarfg_obj->n, zlarfg_obj->alpharef, zlarfg_obj->xref,
								zlarfg_obj->incx, zlarfg_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    zlarfg_obj->info = LAPACKE_zlarfg(zlarfg_obj->n, zlarfg_obj->alpha, zlarfg_obj->x, 
										zlarfg_obj->incx, zlarfg_obj->tau);

    if( zlarfg_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zlarfg is wrong\n", zlarfg_obj->info );
    }
    if( zlarfg_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zlarfg is wrong\n", 
        zlarfg_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zlarfg_obj->diff =  computeDiff_z( zlarfg_obj->bufsize, 
                zlarfg_obj->x, zlarfg_obj->xref );

}

TEST_F(zlarfg_test, zlarfg1) {
    EXPECT_NEAR(0.0, zlarfg_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlarfg_test, zlarfg2) {
    EXPECT_NEAR(0.0, zlarfg_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlarfg_test, zlarfg3) {
    EXPECT_NEAR(0.0, zlarfg_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlarfg_test, zlarfg4) {
    EXPECT_NEAR(0.0, zlarfg_obj->diff, LAPACKE_GTEST_THRESHOLD);
}
