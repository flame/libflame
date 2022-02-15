#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"


#define pteqr_free() \
if (z!=NULL)    free(z); \
if (zref!=NULL) free(zref);\
if (e!=NULL)  free(e);\
if (eref!=NULL) free(eref); \
if (d!=NULL)  free(d);\
if (dref!=NULL) free(dref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
/* Begin float_common_parameters  class definition */
class pteqr_float_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	float* z;
	float* e;
	float* d;
	lapack_int ldz;
	char compz;
	/*Output Parameter*/	
	float *zref, *dref;
	float *eref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      pteqr_float_parameters (int matrix_layout,char compz, lapack_int n, lapack_int ldz);
      ~pteqr_float_parameters ();

};

/* Constructor definition  float_common_parameters */
pteqr_float_parameters:: pteqr_float_parameters (int matrix_layout_i, char compz_i,lapack_int n_i, lapack_int ldz_i)
{
	matrix_layout = matrix_layout_i;
	n = n_i;	
	ldz = ldz_i;
	compz = compz_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n pteqr float:  n: %d ldz: %d \n", n, ldz);
	#endif

	bufsize = ldz*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&z, &zref, bufsize);
	lapacke_gtest_alloc_float_buffer_pair(&d, &dref, (n));
	lapacke_gtest_alloc_float_buffer_pair(&e, &eref, ((n-1)));
	if ((z==NULL) || (zref==NULL) || \
		(d==NULL) || (dref==NULL) || \
		(e == NULL) ||(eref == NULL)){
		EXPECT_FALSE( true) << "pteqr_float_parameters object: malloc error.";
		pteqr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand( z, zref, bufsize);
	lapacke_gtest_init_float_buffer_pair_rand( d, dref, (n));
	lapacke_gtest_init_float_buffer_pair_rand( e, eref, ((n-1)));	

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
pteqr_float_parameters :: ~pteqr_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" pteqr_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   pteqr_free();

}
/*  Test fixture class definition */
class spteqr_test  : public  ::testing::Test {
public:
   pteqr_float_parameters  *spteqr_obj;
   void SetUp();
   void TearDown () { delete spteqr_obj; }
};

void spteqr_test::SetUp(){

    /* LAPACKE spteqr prototype */
    typedef int (*Fptr_NL_LAPACKE_spteqr) (int matrix_layout, char compz, lapack_int n, \
											float* d, float* e, float *z, lapack_int ldz);

    Fptr_NL_LAPACKE_spteqr spteqr;

    spteqr_obj = new pteqr_float_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].compz,
						   eig_paramslist[idx].n,
                           eig_paramslist[idx].ldb);

    idx = Circular_Increment_Index(idx);
	

    spteqr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    spteqr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(spteqr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(spteqr_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    spteqr = (Fptr_NL_LAPACKE_spteqr)dlsym(spteqr_obj->hModule, "LAPACKE_spteqr");
    ASSERT_TRUE(spteqr != NULL) << "failed to get the Netlib LAPACKE_spteqr symbol";
    

    spteqr_obj->inforef = spteqr( spteqr_obj->matrix_layout, spteqr_obj->compz, spteqr_obj->n,\
								spteqr_obj->dref, spteqr_obj->eref, spteqr_obj->zref,\
								spteqr_obj->ldz);

    /* Compute libflame's Lapacke o/p  */
    spteqr_obj->info = LAPACKE_spteqr( spteqr_obj->matrix_layout, spteqr_obj->compz, spteqr_obj->n,\
										spteqr_obj->d, spteqr_obj->e, spteqr_obj->z,\
										spteqr_obj->ldz);

    if( spteqr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_spteqr is wrong\n", spteqr_obj->info );
    }
    if( spteqr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_spteqr is wrong\n", 
        spteqr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    spteqr_obj->diff =  computeDiff_s( spteqr_obj->bufsize, 
                spteqr_obj->z, spteqr_obj->zref );

}

TEST_F(spteqr_test, spteqr1) {
    EXPECT_NEAR(0.0, spteqr_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(spteqr_test, spteqr2) {
    EXPECT_NEAR(0.0, spteqr_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(spteqr_test, spteqr3) {
    EXPECT_NEAR(0.0, spteqr_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(spteqr_test, spteqr4) {
    EXPECT_NEAR(0.0, spteqr_obj->diff, LAPACKE_EIG_THRESHOLD);
}


/* Begin double_common_parameters  class definition */
class pteqr_double_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	double* z;
	double* e;
	double* d;
	lapack_int ldz;
	char compz;
	/*Output Parameter*/	
	double *zref, *dref;
	double *eref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      pteqr_double_parameters (int matrix_layout,char compz, lapack_int n, lapack_int ldz);
      ~pteqr_double_parameters ();

};

/* Constructor definition  double_common_parameters */
pteqr_double_parameters:: pteqr_double_parameters (int matrix_layout_i, char compz_i,lapack_int n_i, lapack_int ldz_i)
{
	matrix_layout = matrix_layout_i;
	n = n_i;	
	ldz = ldz_i;
	compz = compz_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n pteqr double:  n: %d ldz: %d \n", n, ldz);
	#endif

	bufsize = ldz*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&z, &zref, bufsize);
	lapacke_gtest_alloc_double_buffer_pair(&d, &dref, (n));
	lapacke_gtest_alloc_double_buffer_pair(&e, &eref, ((n-1)));
	if ((z==NULL) || (zref==NULL) || \
		(d==NULL) || (dref==NULL) || \
		(e == NULL) ||(eref == NULL)){
		EXPECT_FALSE( true) << "pteqr_double_parameters object: malloc error.";
		pteqr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand( z, zref, bufsize);
	lapacke_gtest_init_double_buffer_pair_rand( d, dref, (n));
	lapacke_gtest_init_double_buffer_pair_rand( e, eref, ((n-1)));	

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
pteqr_double_parameters :: ~pteqr_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" pteqr_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   pteqr_free();

}
/*  Test fixture class definition */
class dpteqr_test  : public  ::testing::Test {
public:
   pteqr_double_parameters  *dpteqr_obj;
   void SetUp();
   void TearDown () { delete dpteqr_obj; }
};

void dpteqr_test::SetUp(){

    /* LAPACKE dpteqr prototype */
    typedef int (*Fptr_NL_LAPACKE_dpteqr) (int matrix_layout, char compz, lapack_int n, \
											double* d, double* e, double *z, lapack_int ldz);

    Fptr_NL_LAPACKE_dpteqr dpteqr;

    dpteqr_obj = new pteqr_double_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].compz,
						   eig_paramslist[idx].n,
                           eig_paramslist[idx].ldb);

    idx = Circular_Increment_Index(idx);
	

    dpteqr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dpteqr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dpteqr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dpteqr_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    dpteqr = (Fptr_NL_LAPACKE_dpteqr)dlsym(dpteqr_obj->hModule, "LAPACKE_dpteqr");
    ASSERT_TRUE(dpteqr != NULL) << "failed to get the Netlib LAPACKE_dpteqr symbol";
    

    dpteqr_obj->inforef = dpteqr( dpteqr_obj->matrix_layout, dpteqr_obj->compz, dpteqr_obj->n,\
								dpteqr_obj->dref, dpteqr_obj->eref, dpteqr_obj->zref,\
								dpteqr_obj->ldz);

    /* Compute libflame's Lapacke o/p  */
    dpteqr_obj->info = LAPACKE_dpteqr( dpteqr_obj->matrix_layout, dpteqr_obj->compz, dpteqr_obj->n,\
										dpteqr_obj->d, dpteqr_obj->e, dpteqr_obj->z,\
										dpteqr_obj->ldz);

    if( dpteqr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dpteqr is wrong\n", dpteqr_obj->info );
    }
    if( dpteqr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dpteqr is wrong\n", 
        dpteqr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dpteqr_obj->diff =  computeDiff_d( dpteqr_obj->bufsize, 
                dpteqr_obj->z, dpteqr_obj->zref );

}

TEST_F(dpteqr_test, dpteqr1) {
    EXPECT_NEAR(0.0, dpteqr_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(dpteqr_test, dpteqr2) {
    EXPECT_NEAR(0.0, dpteqr_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(dpteqr_test, dpteqr3) {
    EXPECT_NEAR(0.0, dpteqr_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(dpteqr_test, dpteqr4) {
    EXPECT_NEAR(0.0, dpteqr_obj->diff, LAPACKE_EIG_THRESHOLD);
}


/* Begin scomplex_common_parameters  class definition */
class pteqr_scomplex_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_complex_float* z;
	float* e;
	float* d;
	lapack_int ldz;
	char compz;
	/*Output Parameter*/	
	lapack_complex_float *zref;
	float *eref, *dref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      pteqr_scomplex_parameters (int matrix_layout,char compz, lapack_int n, lapack_int ldz);
      ~pteqr_scomplex_parameters ();

};

/* Constructor definition  scomplex_common_parameters */
pteqr_scomplex_parameters:: pteqr_scomplex_parameters (int matrix_layout_i, char compz_i,lapack_int n_i, lapack_int ldz_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;	
	ldz = ldz_i;
	compz = compz_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n pteqr scomplex:  n: %d ldz: %d \n", n, ldz);
	#endif

	bufsize = ldz*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&z, &zref, bufsize);
	lapacke_gtest_alloc_float_buffer_pair(&d, &dref, (n));
	lapacke_gtest_alloc_float_buffer_pair(&e, &eref, ((n-1)));
	if ((z==NULL) || (zref==NULL) || \
		(d==NULL) || (dref==NULL) || \
		(e == NULL) ||(eref == NULL)){
		EXPECT_FALSE( true) << "pteqr_scomplex_parameters object: malloc error.";
		pteqr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( z, zref, bufsize);
	lapacke_gtest_init_float_buffer_pair_rand( d, dref, (n));
	lapacke_gtest_init_float_buffer_pair_rand( e, eref, ((n-1)));	

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
pteqr_scomplex_parameters :: ~pteqr_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" pteqr_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   pteqr_free();

}
/*  Test fixture class definition */
class cpteqr_test  : public  ::testing::Test {
public:
   pteqr_scomplex_parameters  *cpteqr_obj;
   void SetUp();
   void TearDown () { delete cpteqr_obj; }
};

void cpteqr_test::SetUp(){

    /* LAPACKE cpteqr prototype */
    typedef int (*Fptr_NL_LAPACKE_cpteqr) (int matrix_layout, char compz, lapack_int n, \
											float* d, float* e, lapack_complex_float *z, lapack_int ldz);

    Fptr_NL_LAPACKE_cpteqr cpteqr;

    cpteqr_obj = new pteqr_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].compz,
						   eig_paramslist[idx].n,
                           eig_paramslist[idx].ldb);

    idx = Circular_Increment_Index(idx);
	

    cpteqr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cpteqr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cpteqr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cpteqr_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    cpteqr = (Fptr_NL_LAPACKE_cpteqr)dlsym(cpteqr_obj->hModule, "LAPACKE_cpteqr");
    ASSERT_TRUE(cpteqr != NULL) << "failed to get the Netlib LAPACKE_cpteqr symbol";
    

    cpteqr_obj->inforef = cpteqr( cpteqr_obj->matrix_layout, cpteqr_obj->compz, cpteqr_obj->n,\
								cpteqr_obj->dref, cpteqr_obj->eref, cpteqr_obj->zref,\
								cpteqr_obj->ldz);

    /* Compute libflame's Lapacke o/p  */
    cpteqr_obj->info = LAPACKE_cpteqr( cpteqr_obj->matrix_layout, cpteqr_obj->compz, cpteqr_obj->n,\
										cpteqr_obj->d, cpteqr_obj->e, cpteqr_obj->z,\
										cpteqr_obj->ldz);

    if( cpteqr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cpteqr is wrong\n", cpteqr_obj->info );
    }
    if( cpteqr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cpteqr is wrong\n", 
        cpteqr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cpteqr_obj->diff =  computeDiff_c( cpteqr_obj->bufsize, 
                cpteqr_obj->z, cpteqr_obj->zref );

}

TEST_F(cpteqr_test, cpteqr1) {
    EXPECT_NEAR(0.0, cpteqr_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(cpteqr_test, cpteqr2) {
    EXPECT_NEAR(0.0, cpteqr_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(cpteqr_test, cpteqr3) {
    EXPECT_NEAR(0.0, cpteqr_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(cpteqr_test, cpteqr4) {
    EXPECT_NEAR(0.0, cpteqr_obj->diff, LAPACKE_EIG_THRESHOLD);
}



/* Begin dcomplex_common_parameters  class definition */
class pteqr_dcomplex_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_complex_double* z;
	double* e;
	double* d;
	lapack_int ldz;
	char compz;
	/*Output Parameter*/	
	lapack_complex_double *zref;
	double *eref, *dref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      pteqr_dcomplex_parameters (int matrix_layout,char compz, lapack_int n, lapack_int ldz);
      ~pteqr_dcomplex_parameters ();

};

/* Constructor definition  dcomplex_common_parameters */
pteqr_dcomplex_parameters:: pteqr_dcomplex_parameters (int matrix_layout_i, char compz_i,lapack_int n_i, lapack_int ldz_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;	
	ldz = ldz_i;
	compz = compz_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n pteqr dcomplex:  n: %d ldz: %d \n", n, ldz);
	#endif

	bufsize = ldz*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&z, &zref, bufsize);
	lapacke_gtest_alloc_double_buffer_pair(&d, &dref, (n));
	lapacke_gtest_alloc_double_buffer_pair(&e, &eref, ((n-1)));
	if ((z==NULL) || (zref==NULL) || \
		(d==NULL) || (dref==NULL) || \
		(e == NULL) ||(eref == NULL)){
		EXPECT_FALSE( true) << "pteqr_dcomplex_parameters object: malloc error.";
		pteqr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( z, zref, bufsize);
	lapacke_gtest_init_double_buffer_pair_rand( d, dref, (n));
	lapacke_gtest_init_double_buffer_pair_rand( e, eref, ((n-1)));	

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
pteqr_dcomplex_parameters :: ~pteqr_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" pteqr_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   pteqr_free();

}
/*  Test fixture class definition */
class zpteqr_test  : public  ::testing::Test {
public:
   pteqr_dcomplex_parameters  *zpteqr_obj;
   void SetUp();
   void TearDown () { delete zpteqr_obj; }
};

void zpteqr_test::SetUp(){

    /* LAPACKE zpteqr prototype */
    typedef int (*Fptr_NL_LAPACKE_zpteqr) (int matrix_layout, char compz, lapack_int n, \
											double* d, double* e, lapack_complex_double *z, lapack_int ldz);

    Fptr_NL_LAPACKE_zpteqr zpteqr;

    zpteqr_obj = new pteqr_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].compz,
						   eig_paramslist[idx].n,
                           eig_paramslist[idx].ldb);

    idx = Circular_Increment_Index(idx);
	

    zpteqr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zpteqr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zpteqr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zpteqr_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zpteqr = (Fptr_NL_LAPACKE_zpteqr)dlsym(zpteqr_obj->hModule, "LAPACKE_zpteqr");
    ASSERT_TRUE(zpteqr != NULL) << "failed to get the Netlib LAPACKE_zpteqr symbol";
    

    zpteqr_obj->inforef = zpteqr( zpteqr_obj->matrix_layout, zpteqr_obj->compz, zpteqr_obj->n,\
								zpteqr_obj->dref, zpteqr_obj->eref, zpteqr_obj->zref,\
								zpteqr_obj->ldz);

    /* Compute libflame's Lapacke o/p  */
    zpteqr_obj->info = LAPACKE_zpteqr( zpteqr_obj->matrix_layout, zpteqr_obj->compz, zpteqr_obj->n,\
										zpteqr_obj->d, zpteqr_obj->e, zpteqr_obj->z,\
										zpteqr_obj->ldz);

    if( zpteqr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zpteqr is wrong\n", zpteqr_obj->info );
    }
    if( zpteqr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zpteqr is wrong\n", 
        zpteqr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zpteqr_obj->diff =  computeDiff_z( zpteqr_obj->bufsize, 
                zpteqr_obj->z, zpteqr_obj->zref );

}

TEST_F(zpteqr_test, zpteqr1) {
    EXPECT_NEAR(0.0, zpteqr_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zpteqr_test, zpteqr2) {
    EXPECT_NEAR(0.0, zpteqr_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zpteqr_test, zpteqr3) {
    EXPECT_NEAR(0.0, zpteqr_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zpteqr_test, zpteqr4) {
    EXPECT_NEAR(0.0, zpteqr_obj->diff, LAPACKE_EIG_THRESHOLD);
}
