#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"


#define larnv_free() \
if (x!=NULL)    free(x); \
if (xref!=NULL) free(xref);\
if (iseed!=NULL)  free(iseed);\
if (iseedref!=NULL) free(iseedref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
/* Begin float_common_parameters  class definition */
class larnv_float_parameters{

   public:
	int bufsize;
	int bufsize_iseed;
	void *hModule, *dModule;
	float diff;
   /*input parameters */	
	lapack_int n;
	lapack_int idist; 
	float * x;
	lapack_int larnv_threshold;
	/*Output Parameter*/
	lapack_int* iseed, *iseedref;
	float *xref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      larnv_float_parameters (lapack_int n, lapack_int idist, lapack_int larnv_threshold);
      ~larnv_float_parameters ();

};

/* Constructor definition  float_common_parameters */
larnv_float_parameters:: larnv_float_parameters ( lapack_int n_i, lapack_int idist_i, lapack_int larnv_threshold_i)
{
	n = n_i;
	idist = idist_i;
	larnv_threshold = larnv_threshold_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n larnv float: n: %d , idist: %d \n", n, idist);
	#endif
	
	bufsize = n;
	bufsize_iseed = 4;
	
	
	if (idist == 3)
		idist = 1; 
	
	
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&x, &xref, bufsize);
	lapacke_gtest_alloc_int_buffer_pair(&iseed, &iseedref, bufsize_iseed);	
	if ((x==NULL) || (xref==NULL) || \
		(iseed==NULL) || (iseedref==NULL)){
		EXPECT_FALSE( true) << "larnv_float_parameters object: malloc error.";
		larnv_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_int_buffer_pair_rand(iseed, iseedref, bufsize_iseed);
	if (!(iseed[3] %2))
		iseed[3] += 1;

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
larnv_float_parameters :: ~larnv_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" larnv_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   larnv_free();

}
/*  Test fixture class definition */
class slarnv_test  : public  ::testing::Test {
public:
   larnv_float_parameters  *slarnv_obj;
   void SetUp();
   void TearDown () { delete slarnv_obj; }
};

void slarnv_test::SetUp(){

    /* LAPACKE slarnv prototype */
    typedef int (*Fptr_NL_LAPACKE_slarnv) (lapack_int idist , lapack_int * iseed , lapack_int n , float * x );

    Fptr_NL_LAPACKE_slarnv slarnv;

    slarnv_obj = new larnv_float_parameters ( eig_paramslist[idx].n,
                           eig_paramslist[idx].itype,
						   eig_paramslist[idx].threshold_value);

    idx = Circular_Increment_Index(idx);
	

    slarnv_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    slarnv_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(slarnv_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(slarnv_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    slarnv = (Fptr_NL_LAPACKE_slarnv)dlsym(slarnv_obj->hModule, "LAPACKE_slarnv");
    ASSERT_TRUE(slarnv != NULL) << "failed to get the Netlib LAPACKE_slarnv symbol";
    

    slarnv_obj->inforef = slarnv(slarnv_obj->idist, slarnv_obj->iseedref,
								slarnv_obj->n, slarnv_obj->xref);

    /* Compute libflame's Lapacke o/p  */
    slarnv_obj->info = LAPACKE_slarnv( slarnv_obj->idist, slarnv_obj->iseed, slarnv_obj->n, slarnv_obj->x);

    if( slarnv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_slarnv is wrong\n", slarnv_obj->info );
    }
    if( slarnv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_slarnv is wrong\n", 
        slarnv_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    slarnv_obj->diff =  computeDiff_s( slarnv_obj->bufsize, 
                slarnv_obj->x, slarnv_obj->xref );

}

TEST_F(slarnv_test, slarnv1) {
    EXPECT_NEAR(0.0, slarnv_obj->diff, slarnv_obj->larnv_threshold);
}

TEST_F(slarnv_test, slarnv2) {
    EXPECT_NEAR(0.0, slarnv_obj->diff, slarnv_obj->larnv_threshold);
}

TEST_F(slarnv_test, slarnv3) {
    EXPECT_NEAR(0.0, slarnv_obj->diff, slarnv_obj->larnv_threshold);
}

TEST_F(slarnv_test, slarnv4) {
    EXPECT_NEAR(0.0, slarnv_obj->diff, slarnv_obj->larnv_threshold);
}


/* Begin double_common_parameters  class definition */
class larnv_double_parameters{

   public:
	int bufsize;
	int bufsize_iseed;
	void *hModule, *dModule;
	double diff;
   /*input parameters */	
	lapack_int n;
	lapack_int idist; 
	lapack_int larnv_threshold;
	double * x;
	/*Output Parameter*/
	lapack_int* iseed, *iseedref;
	double *xref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      larnv_double_parameters (lapack_int n, lapack_int idist, lapack_int larnv_threshold);
      ~larnv_double_parameters ();

};

/* Constructor definition  double_common_parameters */
larnv_double_parameters:: larnv_double_parameters ( lapack_int n_i, lapack_int idist_i, lapack_int larnv_threshold_i)
{
	n = n_i;
	idist = idist_i;
	larnv_threshold = larnv_threshold_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n larnv double: n: %d , idist: %d \n", n, idist);
	#endif
	
	bufsize = n;
	bufsize_iseed = 4;
	
	if (idist == 3)
		idist = 1;
	
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&x, &xref, bufsize);
	lapacke_gtest_alloc_int_buffer_pair(&iseed, &iseedref, bufsize_iseed);	
	if ((x==NULL) || (xref==NULL) || \
		(iseed==NULL) || (iseedref==NULL)){
		EXPECT_FALSE( true) << "larnv_double_parameters object: malloc error.";
		larnv_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_int_buffer_pair_rand(iseed, iseedref, bufsize_iseed);
	if (!(iseed[3] %2))
		iseed[3] += 1;
} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
larnv_double_parameters :: ~larnv_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" larnv_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   larnv_free();

}
/*  Test fixture class definition */
class dlarnv_test  : public  ::testing::Test {
public:
   larnv_double_parameters  *dlarnv_obj;
   void SetUp();
   void TearDown () { delete dlarnv_obj; }
};

void dlarnv_test::SetUp(){

    /* LAPACKE dlarnv prototype */
    typedef int (*Fptr_NL_LAPACKE_dlarnv) (lapack_int idist , lapack_int * iseed , lapack_int n , double * x );

    Fptr_NL_LAPACKE_dlarnv dlarnv;

    dlarnv_obj = new larnv_double_parameters ( eig_paramslist[idx].n,
                           eig_paramslist[idx].itype,
						   eig_paramslist[idx].threshold_value);

    idx = Circular_Increment_Index(idx);
	

    dlarnv_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dlarnv_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dlarnv_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dlarnv_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    dlarnv = (Fptr_NL_LAPACKE_dlarnv)dlsym(dlarnv_obj->hModule, "LAPACKE_dlarnv");
    ASSERT_TRUE(dlarnv != NULL) << "failed to get the Netlib LAPACKE_dlarnv symbol";
    

    dlarnv_obj->inforef = dlarnv(dlarnv_obj->idist, dlarnv_obj->iseedref,
								dlarnv_obj->n, dlarnv_obj->xref);

    /* Compute libflame's Lapacke o/p  */
    dlarnv_obj->info = LAPACKE_dlarnv( dlarnv_obj->idist, dlarnv_obj->iseed, dlarnv_obj->n, dlarnv_obj->x);

    if( dlarnv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dlarnv is wrong\n", dlarnv_obj->info );
    }
    if( dlarnv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dlarnv is wrong\n", 
        dlarnv_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dlarnv_obj->diff =  computeDiff_d( dlarnv_obj->bufsize, 
                dlarnv_obj->x, dlarnv_obj->xref );

}

TEST_F(dlarnv_test, dlarnv1) {
    EXPECT_NEAR(0.0, dlarnv_obj->diff, dlarnv_obj->larnv_threshold);
}

TEST_F(dlarnv_test, dlarnv2) {
    EXPECT_NEAR(0.0, dlarnv_obj->diff, dlarnv_obj->larnv_threshold);
}

TEST_F(dlarnv_test, dlarnv3) {
    EXPECT_NEAR(0.0, dlarnv_obj->diff, dlarnv_obj->larnv_threshold);
}

TEST_F(dlarnv_test, dlarnv4) {
    EXPECT_NEAR(0.0, dlarnv_obj->diff, dlarnv_obj->larnv_threshold);
}


/* Begin scomplex_common_parameters  class definition */
class larnv_scomplex_parameters{

   public:
	int bufsize;
	int bufsize_iseed;
	void *hModule, *dModule;
	float diff;
   /*input parameters */	
	lapack_int n;
	lapack_int idist; 
	lapack_int larnv_threshold;
	lapack_complex_float * x;
	/*Output Parameter*/
	lapack_int* iseed, *iseedref;
	lapack_complex_float *xref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      larnv_scomplex_parameters (lapack_int n, lapack_int idist, lapack_int larnv_threshold);
      ~larnv_scomplex_parameters ();

};

/* Constructor definition  scomplex_common_parameters */
larnv_scomplex_parameters:: larnv_scomplex_parameters ( lapack_int n_i, lapack_int idist_i, lapack_int larnv_threshold_i)
{
	n = n_i;
	idist = idist_i;
	larnv_threshold = larnv_threshold_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n larnv scomplex: n: %d , idist: %d \n", n, idist);
	#endif
	
	bufsize = n;
	bufsize_iseed = 4;
	
	if (idist == 3)
		idist = 1;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&x, &xref, bufsize);
	lapacke_gtest_alloc_int_buffer_pair(&iseed, &iseedref, bufsize_iseed);	
	if ((x==NULL) || (xref==NULL) || \
		(iseed==NULL) || (iseedref==NULL)){
		EXPECT_FALSE( true) << "larnv_scomplex_parameters object: malloc error.";
		larnv_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_int_buffer_pair_rand(iseed, iseedref, bufsize_iseed);
	if (!(iseed[3] %2))
		iseed[3] += 1;

} /* end of Constructor  */

/* Destructor definition  'scomplex_common_parameters' class  */
larnv_scomplex_parameters :: ~larnv_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" larnv_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   larnv_free();

}
/*  Test fixture class definition */
class clarnv_test  : public  ::testing::Test {
public:
   larnv_scomplex_parameters  *clarnv_obj;
   void SetUp();
   void TearDown () { delete clarnv_obj; }
};

void clarnv_test::SetUp(){

    /* LAPACKE clarnv prototype */
    typedef int (*Fptr_NL_LAPACKE_clarnv) (lapack_int idist , lapack_int * iseed , lapack_int n , lapack_complex_float * x );

    Fptr_NL_LAPACKE_clarnv clarnv;

    clarnv_obj = new larnv_scomplex_parameters ( eig_paramslist[idx].n,
                           eig_paramslist[idx].itype, 
						   eig_paramslist[idx].threshold_value);

    idx = Circular_Increment_Index(idx);
	

    clarnv_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    clarnv_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(clarnv_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(clarnv_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    clarnv = (Fptr_NL_LAPACKE_clarnv)dlsym(clarnv_obj->hModule, "LAPACKE_clarnv");
    ASSERT_TRUE(clarnv != NULL) << "failed to get the Netlib LAPACKE_clarnv symbol";
    

    clarnv_obj->inforef = clarnv(clarnv_obj->idist, clarnv_obj->iseedref,
								clarnv_obj->n, clarnv_obj->xref);

    /* Compute libflame's Lapacke o/p  */
    clarnv_obj->info = LAPACKE_clarnv( clarnv_obj->idist, clarnv_obj->iseed, clarnv_obj->n, clarnv_obj->x);

    if( clarnv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_clarnv is wrong\n", clarnv_obj->info );
    }
    if( clarnv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_clarnv is wrong\n", 
        clarnv_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    clarnv_obj->diff =  computeDiff_c( clarnv_obj->bufsize, 
                clarnv_obj->x, clarnv_obj->xref );

}

TEST_F(clarnv_test, clarnv1) {
    EXPECT_NEAR(0.0, clarnv_obj->diff, clarnv_obj->larnv_threshold);
}

TEST_F(clarnv_test, clarnv2) {
    EXPECT_NEAR(0.0, clarnv_obj->diff, clarnv_obj->larnv_threshold);
}

TEST_F(clarnv_test, clarnv3) {
    EXPECT_NEAR(0.0, clarnv_obj->diff, clarnv_obj->larnv_threshold);
}

TEST_F(clarnv_test, clarnv4) {
    EXPECT_NEAR(0.0, clarnv_obj->diff, clarnv_obj->larnv_threshold);
}


/* Begin lapack_complex_dcomplex_common_parameters  class definition */
class larnv_dcomplex_parameters{

   public:
	int bufsize;
	int bufsize_iseed;
	void *hModule, *dModule;
	double diff;
   /*input parameters */	
	lapack_int n;
	lapack_int idist;
	lapack_int larnv_threshold;
	lapack_complex_double * x;
	/*Output Parameter*/
	lapack_int* iseed, *iseedref;
	lapack_complex_double *xref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      larnv_dcomplex_parameters (lapack_int n, lapack_int idist, lapack_int larnv_threshold);
      ~larnv_dcomplex_parameters ();

};

/* Constructor definition  dcomplex_common_parameters */
larnv_dcomplex_parameters:: larnv_dcomplex_parameters ( lapack_int n_i, lapack_int idist_i, lapack_int larnv_threshold_i)
{
	n = n_i;
	idist = idist_i;
	larnv_threshold = larnv_threshold_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n larnv dcomplex: n: %d , idist: %d \n", n, idist);
	#endif
	
	bufsize = n;
	bufsize_iseed = 4;
	
	if (idist == 3)
		idist = 1;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&x, &xref, bufsize);
	lapacke_gtest_alloc_int_buffer_pair(&iseed, &iseedref, bufsize_iseed);	
	if ((x==NULL) || (xref==NULL) || \
		(iseed==NULL) || (iseedref==NULL)){
		EXPECT_FALSE( true) << "larnv_dcomplex_parameters object: malloc error.";
		larnv_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_int_buffer_pair_rand(iseed, iseedref, bufsize_iseed);
	if (!(iseed[3] %2))
		iseed[3] += 1;

} /* end of Constructor  */

/* Destructor definition  'dcomplex_common_parameters' class  */
larnv_dcomplex_parameters :: ~larnv_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" larnv_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   larnv_free();

}
/*  Test fixture class definition */
class zlarnv_test  : public  ::testing::Test {
public:
   larnv_dcomplex_parameters  *zlarnv_obj;
   void SetUp();
   void TearDown () { delete zlarnv_obj; }
};

void zlarnv_test::SetUp(){

    /* LAPACKE zlarnv prototype */
    typedef int (*Fptr_NL_LAPACKE_zlarnv) (lapack_int idist , lapack_int * iseed , lapack_int n , lapack_complex_double * x );

    Fptr_NL_LAPACKE_zlarnv zlarnv;

    zlarnv_obj = new larnv_dcomplex_parameters ( eig_paramslist[idx].n,
                            eig_paramslist[idx].itype,
							eig_paramslist[idx].threshold_value);

    idx = Circular_Increment_Index(idx);
	

    zlarnv_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zlarnv_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zlarnv_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zlarnv_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zlarnv = (Fptr_NL_LAPACKE_zlarnv)dlsym(zlarnv_obj->hModule, "LAPACKE_zlarnv");
    ASSERT_TRUE(zlarnv != NULL) << "failed to get the Netlib LAPACKE_zlarnv symbol";
    

    zlarnv_obj->inforef = zlarnv(zlarnv_obj->idist, zlarnv_obj->iseedref,
								zlarnv_obj->n, zlarnv_obj->xref);

    /* Compute libflame's Lapacke o/p  */
    zlarnv_obj->info = LAPACKE_zlarnv( zlarnv_obj->idist, zlarnv_obj->iseed, zlarnv_obj->n, zlarnv_obj->x);

    if( zlarnv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zlarnv is wrong\n", zlarnv_obj->info );
    }
    if( zlarnv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zlarnv is wrong\n", 
        zlarnv_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zlarnv_obj->diff =  computeDiff_z( zlarnv_obj->bufsize, 
                zlarnv_obj->x, zlarnv_obj->xref );

}

TEST_F(zlarnv_test, zlarnv1) {
    EXPECT_NEAR(0.0, zlarnv_obj->diff, zlarnv_obj->larnv_threshold);
}

TEST_F(zlarnv_test, zlarnv2) {
    EXPECT_NEAR(0.0, zlarnv_obj->diff, zlarnv_obj->larnv_threshold);
}

TEST_F(zlarnv_test, zlarnv3) {
    EXPECT_NEAR(0.0, zlarnv_obj->diff, zlarnv_obj->larnv_threshold);
}

TEST_F(zlarnv_test, zlarnv4) {
    EXPECT_NEAR(0.0, zlarnv_obj->diff, zlarnv_obj->larnv_threshold);
}
