#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"


#define steqr_free() \
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
class steqr_float_parameters{

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
      steqr_float_parameters (int matrix_layout,char compz, lapack_int n, lapack_int ldz);
      ~steqr_float_parameters ();

};

/* Constructor definition  float_common_parameters */
steqr_float_parameters:: steqr_float_parameters (int matrix_layout_i, char compz_i,lapack_int n_i, lapack_int ldz_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;	
	ldz = ldz_i;
	compz = compz_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n steqr float:  n: %d ldz: %d \n", n, ldz);
	#endif

	bufsize = ldz*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&z, &zref, bufsize);
	lapacke_gtest_alloc_float_buffer_pair(&d, &dref, (n));
	lapacke_gtest_alloc_float_buffer_pair(&e, &eref, ((n-1)));
	if ((z==NULL) || (zref==NULL) || \
		(d==NULL) || (dref==NULL) || \
		(e == NULL) ||(eref == NULL)){
		EXPECT_FALSE( true) << "steqr_float_parameters object: malloc error.";
		steqr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand( z, zref, bufsize);
	lapacke_gtest_init_float_buffer_pair_rand( d, dref, (n));
	lapacke_gtest_init_float_buffer_pair_rand( e, eref, ((n-1)));	

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
steqr_float_parameters :: ~steqr_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" steqr_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   steqr_free();

}
/*  Test fixture class definition */
class ssteqr_test  : public  ::testing::Test {
public:
   steqr_float_parameters  *ssteqr_obj;
   void SetUp();
   void TearDown () { delete ssteqr_obj; }
};

void ssteqr_test::SetUp(){

    /* LAPACKE ssteqr prototype */
    typedef int (*Fptr_NL_LAPACKE_ssteqr) (int matrix_layout, char compz, lapack_int n, \
											float* d, float* e, float *z, lapack_int ldz);

    Fptr_NL_LAPACKE_ssteqr ssteqr;

    ssteqr_obj = new steqr_float_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].compz,
						   eig_paramslist[idx].n,
                           eig_paramslist[idx].ldb);

    idx = Circular_Increment_Index(idx);
	

    ssteqr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ssteqr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ssteqr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ssteqr_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    ssteqr = (Fptr_NL_LAPACKE_ssteqr)dlsym(ssteqr_obj->hModule, "LAPACKE_ssteqr");
    ASSERT_TRUE(ssteqr != NULL) << "failed to get the Netlib LAPACKE_ssteqr symbol";
    

    ssteqr_obj->inforef = ssteqr( ssteqr_obj->matrix_layout, ssteqr_obj->compz, ssteqr_obj->n,\
								ssteqr_obj->dref, ssteqr_obj->eref, ssteqr_obj->zref,\
								ssteqr_obj->ldz);

    /* Compute libflame's Lapacke o/p  */
    ssteqr_obj->info = LAPACKE_ssteqr( ssteqr_obj->matrix_layout, ssteqr_obj->compz, ssteqr_obj->n,\
										ssteqr_obj->d, ssteqr_obj->e, ssteqr_obj->z,\
										ssteqr_obj->ldz);

    if( ssteqr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ssteqr is wrong\n", ssteqr_obj->info );
    }
    if( ssteqr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ssteqr is wrong\n", 
        ssteqr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    ssteqr_obj->diff =  computeDiff_s( ssteqr_obj->bufsize, 
                ssteqr_obj->z, ssteqr_obj->zref );

}

TEST_F(ssteqr_test, ssteqr1) {
    EXPECT_NEAR(0.0, ssteqr_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(ssteqr_test, ssteqr2) {
    EXPECT_NEAR(0.0, ssteqr_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(ssteqr_test, ssteqr3) {
    EXPECT_NEAR(0.0, ssteqr_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(ssteqr_test, ssteqr4) {
    EXPECT_NEAR(0.0, ssteqr_obj->diff, LAPACKE_EIG_THRESHOLD);
}


/* Begin double_common_parameters  class definition */
class steqr_double_parameters{

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
      steqr_double_parameters (int matrix_layout,char compz, lapack_int n, lapack_int ldz);
      ~steqr_double_parameters ();

};

/* Constructor definition  double_common_parameters */
steqr_double_parameters:: steqr_double_parameters (int matrix_layout_i, char compz_i,lapack_int n_i, lapack_int ldz_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;	
	ldz = ldz_i;
	compz = compz_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n steqr double:  n: %d ldz: %d \n", n, ldz);
	#endif

	bufsize = ldz*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&z, &zref, bufsize);
	lapacke_gtest_alloc_double_buffer_pair(&d, &dref, (n));
	lapacke_gtest_alloc_double_buffer_pair(&e, &eref, ((n-1)));
	if ((z==NULL) || (zref==NULL) || \
		(d==NULL) || (dref==NULL) || \
		(e == NULL) ||(eref == NULL)){
		EXPECT_FALSE( true) << "steqr_double_parameters object: malloc error.";
		steqr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand( z, zref, bufsize);
	lapacke_gtest_init_double_buffer_pair_rand( d, dref, (n));
	lapacke_gtest_init_double_buffer_pair_rand( e, eref, ((n-1)));	

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
steqr_double_parameters :: ~steqr_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" steqr_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   steqr_free();

}
/*  Test fixture class definition */
class dsteqr_test  : public  ::testing::Test {
public:
   steqr_double_parameters  *dsteqr_obj;
   void SetUp();
   void TearDown () { delete dsteqr_obj; }
};

void dsteqr_test::SetUp(){

    /* LAPACKE dsteqr prototype */
    typedef int (*Fptr_NL_LAPACKE_dsteqr) (int matrix_layout, char compz, lapack_int n, \
											double* d, double* e, double *z, lapack_int ldz);

    Fptr_NL_LAPACKE_dsteqr dsteqr;

    dsteqr_obj = new steqr_double_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].compz,
						   eig_paramslist[idx].n,
                           eig_paramslist[idx].ldb);

    idx = Circular_Increment_Index(idx);
	

    dsteqr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dsteqr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dsteqr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dsteqr_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    dsteqr = (Fptr_NL_LAPACKE_dsteqr)dlsym(dsteqr_obj->hModule, "LAPACKE_dsteqr");
    ASSERT_TRUE(dsteqr != NULL) << "failed to get the Netlib LAPACKE_dsteqr symbol";
    

    dsteqr_obj->inforef = dsteqr( dsteqr_obj->matrix_layout, dsteqr_obj->compz, dsteqr_obj->n,\
								dsteqr_obj->dref, dsteqr_obj->eref, dsteqr_obj->zref,\
								dsteqr_obj->ldz);

    /* Compute libflame's Lapacke o/p  */
    dsteqr_obj->info = LAPACKE_dsteqr( dsteqr_obj->matrix_layout, dsteqr_obj->compz, dsteqr_obj->n,\
										dsteqr_obj->d, dsteqr_obj->e, dsteqr_obj->z,\
										dsteqr_obj->ldz);

    if( dsteqr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dsteqr is wrong\n", dsteqr_obj->info );
    }
    if( dsteqr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dsteqr is wrong\n", 
        dsteqr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dsteqr_obj->diff =  computeDiff_d( dsteqr_obj->bufsize, 
                dsteqr_obj->z, dsteqr_obj->zref );

}

TEST_F(dsteqr_test, dsteqr1) {
    EXPECT_NEAR(0.0, dsteqr_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(dsteqr_test, dsteqr2) {
    EXPECT_NEAR(0.0, dsteqr_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(dsteqr_test, dsteqr3) {
    EXPECT_NEAR(0.0, dsteqr_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(dsteqr_test, dsteqr4) {
    EXPECT_NEAR(0.0, dsteqr_obj->diff, LAPACKE_EIG_THRESHOLD);
}


/* Begin scomplex_common_parameters  class definition */
class steqr_scomplex_parameters{

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
      steqr_scomplex_parameters (int matrix_layout,char compz, lapack_int n, lapack_int ldz);
      ~steqr_scomplex_parameters ();

};

/* Constructor definition  scomplex_common_parameters */
steqr_scomplex_parameters:: steqr_scomplex_parameters (int matrix_layout_i, char compz_i,lapack_int n_i, lapack_int ldz_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;	
	ldz = ldz_i;
	compz = compz_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n steqr scomplex:  n: %d ldz: %d \n", n, ldz);
	#endif

	bufsize = ldz*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&z, &zref, bufsize);
	lapacke_gtest_alloc_float_buffer_pair(&d, &dref, (n));
	lapacke_gtest_alloc_float_buffer_pair(&e, &eref, ((n-1)));
	if ((z==NULL) || (zref==NULL) || \
		(d==NULL) || (dref==NULL) || \
		(e == NULL) ||(eref == NULL)){
		EXPECT_FALSE( true) << "steqr_scomplex_parameters object: malloc error.";
		steqr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( z, zref, bufsize);
	lapacke_gtest_init_float_buffer_pair_rand( d, dref, (n));
	lapacke_gtest_init_float_buffer_pair_rand( e, eref, ((n-1)));	

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
steqr_scomplex_parameters :: ~steqr_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" steqr_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   steqr_free();

}
/*  Test fixture class definition */
class csteqr_test  : public  ::testing::Test {
public:
   steqr_scomplex_parameters  *csteqr_obj;
   void SetUp();
   void TearDown () { delete csteqr_obj; }
};

void csteqr_test::SetUp(){

    /* LAPACKE csteqr prototype */
    typedef int (*Fptr_NL_LAPACKE_csteqr) (int matrix_layout, char compz, lapack_int n, \
											float* d, float* e, lapack_complex_float *z, lapack_int ldz);

    Fptr_NL_LAPACKE_csteqr csteqr;

    csteqr_obj = new steqr_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].compz,
						   eig_paramslist[idx].n,
                           eig_paramslist[idx].ldb);

    idx = Circular_Increment_Index(idx);
	

    csteqr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    csteqr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(csteqr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(csteqr_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    csteqr = (Fptr_NL_LAPACKE_csteqr)dlsym(csteqr_obj->hModule, "LAPACKE_csteqr");
    ASSERT_TRUE(csteqr != NULL) << "failed to get the Netlib LAPACKE_csteqr symbol";
    

    csteqr_obj->inforef = csteqr( csteqr_obj->matrix_layout, csteqr_obj->compz, csteqr_obj->n,\
								csteqr_obj->dref, csteqr_obj->eref, csteqr_obj->zref,\
								csteqr_obj->ldz);

    /* Compute libflame's Lapacke o/p  */
    csteqr_obj->info = LAPACKE_csteqr( csteqr_obj->matrix_layout, csteqr_obj->compz, csteqr_obj->n,\
										csteqr_obj->d, csteqr_obj->e, csteqr_obj->z,\
										csteqr_obj->ldz);

    if( csteqr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_csteqr is wrong\n", csteqr_obj->info );
    }
    if( csteqr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_csteqr is wrong\n", 
        csteqr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    csteqr_obj->diff =  computeDiff_c( csteqr_obj->bufsize, 
                csteqr_obj->z, csteqr_obj->zref );

}

TEST_F(csteqr_test, csteqr1) {
    EXPECT_NEAR(0.0, csteqr_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(csteqr_test, csteqr2) {
    EXPECT_NEAR(0.0, csteqr_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(csteqr_test, csteqr3) {
    EXPECT_NEAR(0.0, csteqr_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(csteqr_test, csteqr4) {
    EXPECT_NEAR(0.0, csteqr_obj->diff, LAPACKE_EIG_THRESHOLD);
}



/* Begin dcomplex_common_parameters  class definition */
class steqr_dcomplex_parameters{

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
      steqr_dcomplex_parameters (int matrix_layout,char compz, lapack_int n, lapack_int ldz);
      ~steqr_dcomplex_parameters ();

};

/* Constructor definition  dcomplex_common_parameters */
steqr_dcomplex_parameters:: steqr_dcomplex_parameters (int matrix_layout_i, char compz_i,lapack_int n_i, lapack_int ldz_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;	
	ldz = ldz_i;
	compz = compz_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n steqr dcomplex:  n: %d ldz: %d \n", n, ldz);
	#endif

	bufsize = ldz*n;	
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&z, &zref, bufsize);
	lapacke_gtest_alloc_double_buffer_pair(&d, &dref, (n));
	lapacke_gtest_alloc_double_buffer_pair(&e, &eref, ((n-1)));
	if ((z==NULL) || (zref==NULL) || \
		(d==NULL) || (dref==NULL) || \
		(e == NULL) ||(eref == NULL)){
		EXPECT_FALSE( true) << "steqr_dcomplex_parameters object: malloc error.";
		steqr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( z, zref, bufsize);
	lapacke_gtest_init_double_buffer_pair_rand( d, dref, (n));
	lapacke_gtest_init_double_buffer_pair_rand( e, eref, ((n-1)));	

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
steqr_dcomplex_parameters :: ~steqr_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" steqr_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   steqr_free();

}
/*  Test fixture class definition */
class zsteqr_test  : public  ::testing::Test {
public:
   steqr_dcomplex_parameters  *zsteqr_obj;
   void SetUp();
   void TearDown () { delete zsteqr_obj; }
};

void zsteqr_test::SetUp(){

    /* LAPACKE zsteqr prototype */
    typedef int (*Fptr_NL_LAPACKE_zsteqr) (int matrix_layout, char compz, lapack_int n, \
											double* d, double* e, lapack_complex_double *z, lapack_int ldz);

    Fptr_NL_LAPACKE_zsteqr zsteqr;

    zsteqr_obj = new steqr_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].compz,
						   eig_paramslist[idx].n,
                           eig_paramslist[idx].ldb);

    idx = Circular_Increment_Index(idx);
	

    zsteqr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zsteqr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zsteqr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zsteqr_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zsteqr = (Fptr_NL_LAPACKE_zsteqr)dlsym(zsteqr_obj->hModule, "LAPACKE_zsteqr");
    ASSERT_TRUE(zsteqr != NULL) << "failed to get the Netlib LAPACKE_zsteqr symbol";
    

    zsteqr_obj->inforef = zsteqr( zsteqr_obj->matrix_layout, zsteqr_obj->compz, zsteqr_obj->n,\
								zsteqr_obj->dref, zsteqr_obj->eref, zsteqr_obj->zref,\
								zsteqr_obj->ldz);

    /* Compute libflame's Lapacke o/p  */
    zsteqr_obj->info = LAPACKE_zsteqr( zsteqr_obj->matrix_layout, zsteqr_obj->compz, zsteqr_obj->n,\
										zsteqr_obj->d, zsteqr_obj->e, zsteqr_obj->z,\
										zsteqr_obj->ldz);

    if( zsteqr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zsteqr is wrong\n", zsteqr_obj->info );
    }
    if( zsteqr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zsteqr is wrong\n", 
        zsteqr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zsteqr_obj->diff =  computeDiff_z( zsteqr_obj->bufsize, 
                zsteqr_obj->z, zsteqr_obj->zref );

}

TEST_F(zsteqr_test, zsteqr1) {
    EXPECT_NEAR(0.0, zsteqr_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zsteqr_test, zsteqr2) {
    EXPECT_NEAR(0.0, zsteqr_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zsteqr_test, zsteqr3) {
    EXPECT_NEAR(0.0, zsteqr_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zsteqr_test, zsteqr4) {
    EXPECT_NEAR(0.0, zsteqr_obj->diff, LAPACKE_EIG_THRESHOLD);
}
