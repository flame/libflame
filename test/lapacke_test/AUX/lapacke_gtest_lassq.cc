#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"

#define lassq_free() \
if (x!=NULL)    free(x); \
if (xref!=NULL)    free(xref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
/* Begin float_common_parameters  class definition */
class lassq_float_parameters{

   public:
	int bufsize_x;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n, incx;
	float scale, sumsq;
	float* x;
	/*Output Parameter*/
	float  *xref;
	/*Return Values*/	
	lapack_int info, inforef;
   public:
      lassq_float_parameters (lapack_int n,  lapack_int incx);
      ~lassq_float_parameters ();

};

/* Constructor definition  float_common_parameters */
lassq_float_parameters:: lassq_float_parameters (lapack_int n_i, lapack_int incx_i)
{	
	n = n_i;
	incx = incx_i;
	scale = sumsq = 1;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n lassq float: n = %d,  incx: %d \n", n, incx);
	#endif
	
	/*Matrix sizes*/
	bufsize_x = 1+(n-1)*incx;
	
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&x, &xref, bufsize_x);	
	if ((x == NULL) || (xref == NULL)){
		EXPECT_FALSE( true) << "lassq_float_parameters object: malloc error.";
		lassq_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand(x, xref, bufsize_x);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
lassq_float_parameters :: ~lassq_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" lassq_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lassq_free();

}
/*  Test fixture class definition */
class slassq_test  : public  ::testing::Test {
public:
   lassq_float_parameters  *slassq_obj;
   void SetUp();
   void TearDown () { delete slassq_obj;}
};

void slassq_test::SetUp(){

	 /* LAPACKE slassq prototype */
    typedef int (*Fptr_NL_LAPACKE_slassq) (lapack_int n,  float* x, lapack_int incx,  float* scale,  float* sumsq );

    Fptr_NL_LAPACKE_slassq slassq;

    slassq_obj = new lassq_float_parameters (lin_solver_paramslist[idx].n,
						   lin_solver_paramslist[idx].kl);
						   

    idx = Circular_Increment_Index(idx);

    slassq_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    slassq_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(slassq_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(slassq_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*slassq library call */
    slassq = (Fptr_NL_LAPACKE_slassq)dlsym(slassq_obj->hModule, "LAPACKE_slassq");
    ASSERT_TRUE(slassq != NULL) << "failed to get the Netlib LAPACKE_slassq symbol";
    
/*Compute slassq's  o/p */
    slassq_obj->inforef = slassq(slassq_obj->n, slassq_obj->xref, slassq_obj->incx, &slassq_obj->scale, &slassq_obj->sumsq);

    /* Compute libflame's Lapacke o/p  */
    slassq_obj->info = LAPACKE_slassq( slassq_obj->n, slassq_obj->x, slassq_obj->incx, &slassq_obj->scale, &slassq_obj->sumsq);

    if( slassq_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_slassq is wrong\n", slassq_obj->info );
    }
    if( slassq_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_slassq is wrong\n", 
        slassq_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    slassq_obj->diff =  computeDiff_s( slassq_obj->bufsize_x, 
                slassq_obj->x, slassq_obj->xref );

}

TEST_F(slassq_test, slassq1) {
    EXPECT_NEAR(0.0, slassq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(slassq_test, slassq2) {
    EXPECT_NEAR(0.0, slassq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(slassq_test, slassq3) {
    EXPECT_NEAR(0.0, slassq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(slassq_test, slassq4) {
    EXPECT_NEAR(0.0, slassq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin double_common_parameters  class definition */
class lassq_double_parameters{

   public:
	int bufsize_x;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n, incx;
	double scale, sumsq;
	double* x;
	/*Output Parameter*/
	double  *xref;
	/*Return Values*/	
	lapack_int info, inforef;
   public:
      lassq_double_parameters (lapack_int n,  lapack_int incx);
      ~lassq_double_parameters ();

};

/* Constructor definition  double_common_parameters */
lassq_double_parameters:: lassq_double_parameters (lapack_int n_i, lapack_int incx_i)
{	
	n = n_i;
	incx = incx_i;
	scale = sumsq = 1;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n lassq double: n = %d,  incx: %d \n", n, incx);
	#endif
	
	/*Matrix sizes*/
	bufsize_x = 1+(n-1)*incx;
	
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&x, &xref, bufsize_x);	
	if ((x == NULL) || (xref == NULL)){
		EXPECT_FALSE( true) << "lassq_double_parameters object: malloc error.";
		lassq_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand(x, xref, bufsize_x);

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
lassq_double_parameters :: ~lassq_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" lassq_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lassq_free();

}
/*  Test fixture class definition */
class dlassq_test  : public  ::testing::Test {
public:
   lassq_double_parameters  *dlassq_obj;
   void SetUp();
   void TearDown () { delete dlassq_obj;}
};

void dlassq_test::SetUp(){

	 /* LAPACKE dlassq prototype */
    typedef int (*Fptr_NL_LAPACKE_dlassq) (lapack_int n,  double* x, lapack_int incx,  double* scale,  double* sumsq );

    Fptr_NL_LAPACKE_dlassq dlassq;

    dlassq_obj = new lassq_double_parameters (lin_solver_paramslist[idx].n,
						   lin_solver_paramslist[idx].kl);
						   

    idx = Circular_Increment_Index(idx);

    dlassq_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dlassq_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dlassq_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dlassq_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*dlassq library call */
    dlassq = (Fptr_NL_LAPACKE_dlassq)dlsym(dlassq_obj->hModule, "LAPACKE_dlassq");
    ASSERT_TRUE(dlassq != NULL) << "failed to get the Netlib LAPACKE_dlassq symbol";
    
/*Compute dlassq's  o/p */
    dlassq_obj->inforef = dlassq(dlassq_obj->n, dlassq_obj->xref, dlassq_obj->incx, &dlassq_obj->scale, &dlassq_obj->sumsq);

    /* Compute libflame's Lapacke o/p  */
    dlassq_obj->info = LAPACKE_dlassq( dlassq_obj->n, dlassq_obj->x, dlassq_obj->incx, &dlassq_obj->scale, &dlassq_obj->sumsq);

    if( dlassq_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dlassq is wrong\n", dlassq_obj->info );
    }
    if( dlassq_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dlassq is wrong\n", 
        dlassq_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dlassq_obj->diff =  computeDiff_d( dlassq_obj->bufsize_x, 
                dlassq_obj->x, dlassq_obj->xref );

}

TEST_F(dlassq_test, dlassq1) {
    EXPECT_NEAR(0.0, dlassq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dlassq_test, dlassq2) {
    EXPECT_NEAR(0.0, dlassq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dlassq_test, dlassq3) {
    EXPECT_NEAR(0.0, dlassq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dlassq_test, dlassq4) {
    EXPECT_NEAR(0.0, dlassq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin scomplex_common_parameters  class definition */
class lassq_scomplex_parameters{

   public:
	int bufsize_x;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n, incx;
	float scale, sumsq;
	lapack_complex_float* x;
	/*Output Parameter*/
	lapack_complex_float  *xref;
	/*Return Values*/	
	lapack_int info, inforef;
   public:
      lassq_scomplex_parameters (lapack_int n,  lapack_int incx);
      ~lassq_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
lassq_scomplex_parameters:: lassq_scomplex_parameters (lapack_int n_i, lapack_int incx_i)
{	
	n = n_i;
	incx = incx_i;
	scale = sumsq = 1;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n lassq scomplex: n = %d,  incx: %d \n", n, incx);
	#endif
	
	/*Matrix sizes*/
	bufsize_x = 1+(n-1)*incx;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&x, &xref, bufsize_x);	
	if ((x == NULL) || (xref == NULL)){
		EXPECT_FALSE( true) << "lassq_float_parameters object: malloc error.";
		lassq_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand(x, xref, bufsize_x);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
lassq_scomplex_parameters :: ~lassq_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" lassq_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lassq_free();

}
/*  Test fixture class definition */
class classq_test  : public  ::testing::Test {
public:
   lassq_scomplex_parameters  *classq_obj;
   void SetUp();
   void TearDown () { delete classq_obj;}
};

void classq_test::SetUp(){

	 /* LAPACKE classq prototype */
    typedef int (*Fptr_NL_LAPACKE_classq) (lapack_int n,  lapack_complex_float* x, lapack_int incx,  float* scale,  float* sumsq );

    Fptr_NL_LAPACKE_classq classq;

    classq_obj = new lassq_scomplex_parameters (lin_solver_paramslist[idx].n,
						   lin_solver_paramslist[idx].kl);
						   

    idx = Circular_Increment_Index(idx);

    classq_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    classq_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(classq_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(classq_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*classq library call */
    classq = (Fptr_NL_LAPACKE_classq)dlsym(classq_obj->hModule, "LAPACKE_classq");
    ASSERT_TRUE(classq != NULL) << "failed to get the Netlib LAPACKE_classq symbol";
    
/*Compute classq's  o/p */
    classq_obj->inforef = classq(classq_obj->n, classq_obj->xref, classq_obj->incx, &classq_obj->scale, &classq_obj->sumsq);

    /* Compute libflame's Lapacke o/p  */
    classq_obj->info = LAPACKE_classq( classq_obj->n, classq_obj->x, classq_obj->incx, &classq_obj->scale, &classq_obj->sumsq);

    if( classq_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_classq is wrong\n", classq_obj->info );
    }
    if( classq_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_classq is wrong\n", 
        classq_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    classq_obj->diff =  computeDiff_c( classq_obj->bufsize_x, 
                classq_obj->x, classq_obj->xref );

}

TEST_F(classq_test, classq1) {
    EXPECT_NEAR(0.0, classq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(classq_test, classq2) {
    EXPECT_NEAR(0.0, classq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(classq_test, classq3) {
    EXPECT_NEAR(0.0, classq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(classq_test, classq4) {
    EXPECT_NEAR(0.0, classq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class lassq_dcomplex_parameters{

   public:
	int bufsize_x;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n, incx;
	double scale, sumsq;
	lapack_complex_double* x;
	/*Output Parameter*/
	lapack_complex_double  *xref;
	/*Return Values*/	
	lapack_int info, inforef;
   public:
      lassq_dcomplex_parameters (lapack_int n,  lapack_int incx);
      ~lassq_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
lassq_dcomplex_parameters:: lassq_dcomplex_parameters (lapack_int n_i, lapack_int incx_i)
{	
	n = n_i;
	incx = incx_i;
	scale = sumsq = 1;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n lassq dcomplex: n = %d,  incx: %d \n", n, incx);
	#endif
	
	/*Matrix sizes*/
	bufsize_x = 1+(n-1)*incx;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&x, &xref, bufsize_x);	
	if ((x == NULL) || (xref == NULL)){
		EXPECT_FALSE( true) << "lassq_double_parameters object: malloc error.";
		lassq_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand(x, xref, bufsize_x);

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
lassq_dcomplex_parameters :: ~lassq_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" lassq_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lassq_free();

}
/*  Test fixture class definition */
class zlassq_test  : public  ::testing::Test {
public:
   lassq_dcomplex_parameters  *zlassq_obj;
   void SetUp();
   void TearDown () { delete zlassq_obj;}
};

void zlassq_test::SetUp(){

	 /* LAPACKE zlassq prototype */
    typedef int (*Fptr_NL_LAPACKE_zlassq) (lapack_int n,  lapack_complex_double* x, lapack_int incx,  double* scale,  double* sumsq );

    Fptr_NL_LAPACKE_zlassq zlassq;

    zlassq_obj = new lassq_dcomplex_parameters (lin_solver_paramslist[idx].n,
						   lin_solver_paramslist[idx].kl);
						   

    idx = Circular_Increment_Index(idx);

    zlassq_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zlassq_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zlassq_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zlassq_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*zlassq library call */
    zlassq = (Fptr_NL_LAPACKE_zlassq)dlsym(zlassq_obj->hModule, "LAPACKE_zlassq");
    ASSERT_TRUE(zlassq != NULL) << "failed to get the Netlib LAPACKE_zlassq symbol";
    
/*Compute zlassq's  o/p */
    zlassq_obj->inforef = zlassq(zlassq_obj->n, zlassq_obj->xref, zlassq_obj->incx, &zlassq_obj->scale, &zlassq_obj->sumsq);

    /* Compute libflame's Lapacke o/p  */
    zlassq_obj->info = LAPACKE_zlassq( zlassq_obj->n, zlassq_obj->x, zlassq_obj->incx, &zlassq_obj->scale, &zlassq_obj->sumsq);

    if( zlassq_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zlassq is wrong\n", zlassq_obj->info );
    }
    if( zlassq_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zlassq is wrong\n", 
        zlassq_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zlassq_obj->diff =  computeDiff_z( zlassq_obj->bufsize_x, 
                zlassq_obj->x, zlassq_obj->xref );

}

TEST_F(zlassq_test, zlassq1) {
    EXPECT_NEAR(0.0, zlassq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlassq_test, zlassq2) {
    EXPECT_NEAR(0.0, zlassq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlassq_test, zlassq3) {
    EXPECT_NEAR(0.0, zlassq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlassq_test, zlassq4) {
    EXPECT_NEAR(0.0, zlassq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}