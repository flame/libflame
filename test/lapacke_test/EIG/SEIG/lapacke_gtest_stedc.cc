#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"


#define stedc_free() \
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
class stedc_float_parameters{

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
      stedc_float_parameters (int matrix_layout,char compz, lapack_int n, lapack_int ldz);
      ~stedc_float_parameters ();

};

/* Constructor definition  float_common_parameters */
stedc_float_parameters:: stedc_float_parameters (int matrix_layout_i, char compz_i,lapack_int n_i, lapack_int ldz_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;	
	ldz = ldz_i;
	compz = compz_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n stedc float:  n: %d ldz: %d \n", n, ldz);
	#endif

	bufsize = ldz*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&z, &zref, bufsize);
	lapacke_gtest_alloc_float_buffer_pair(&d, &dref, (n));
	lapacke_gtest_alloc_float_buffer_pair(&e, &eref, ((n-1)));
	if ((z==NULL) || (zref==NULL) || \
		(d==NULL) || (dref==NULL) || \
		(e == NULL) ||(eref == NULL)){
		EXPECT_FALSE( true) << "stedc_float_parameters object: malloc error.";
		stedc_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand( z, zref, bufsize);
	lapacke_gtest_init_float_buffer_pair_rand( d, dref, (n));
	lapacke_gtest_init_float_buffer_pair_rand( e, eref, ((n-1)));	

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
stedc_float_parameters :: ~stedc_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" stedc_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   stedc_free();

}
/*  Test fixture class definition */
class sstedc_test  : public  ::testing::Test {
public:
   stedc_float_parameters  *sstedc_obj;
   void SetUp();
   void TearDown () { delete sstedc_obj; }
};

void sstedc_test::SetUp(){

    /* LAPACKE sstedc prototype */
    typedef int (*Fptr_NL_LAPACKE_sstedc) (int matrix_layout, char compz, lapack_int n, \
											float* d, float* e, float *z, lapack_int ldz);

    Fptr_NL_LAPACKE_sstedc sstedc;

    sstedc_obj = new stedc_float_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].compz,
						   eig_paramslist[idx].n,
                           eig_paramslist[idx].ldb);

    idx = Circular_Increment_Index(idx);
	

    sstedc_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sstedc_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sstedc_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sstedc_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    sstedc = (Fptr_NL_LAPACKE_sstedc)dlsym(sstedc_obj->hModule, "LAPACKE_sstedc");
    ASSERT_TRUE(sstedc != NULL) << "failed to get the Netlib LAPACKE_sstedc symbol";
    

    sstedc_obj->inforef = sstedc( sstedc_obj->matrix_layout, sstedc_obj->compz, sstedc_obj->n,\
								sstedc_obj->dref, sstedc_obj->eref, sstedc_obj->zref,\
								sstedc_obj->ldz);

    /* Compute libflame's Lapacke o/p  */
    sstedc_obj->info = LAPACKE_sstedc( sstedc_obj->matrix_layout, sstedc_obj->compz, sstedc_obj->n,\
										sstedc_obj->d, sstedc_obj->e, sstedc_obj->z,\
										sstedc_obj->ldz);

    if( sstedc_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_sstedc is wrong\n", sstedc_obj->info );
    }
    if( sstedc_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sstedc is wrong\n", 
        sstedc_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    sstedc_obj->diff =  computeDiff_s( sstedc_obj->bufsize, 
                sstedc_obj->z, sstedc_obj->zref );

}

TEST_F(sstedc_test, sstedc1) {
    EXPECT_NEAR(0.0, sstedc_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sstedc_test, sstedc2) {
    EXPECT_NEAR(0.0, sstedc_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sstedc_test, sstedc3) {
    EXPECT_NEAR(0.0, sstedc_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sstedc_test, sstedc4) {
    EXPECT_NEAR(0.0, sstedc_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin double_common_parameters  class definition */
class stedc_double_parameters{

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
      stedc_double_parameters (int matrix_layout,char compz, lapack_int n, lapack_int ldz);
      ~stedc_double_parameters ();

};

/* Constructor definition  double_common_parameters */
stedc_double_parameters:: stedc_double_parameters (int matrix_layout_i, char compz_i,lapack_int n_i, lapack_int ldz_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;	
	ldz = ldz_i;
	compz = compz_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n stedc double:  n: %d ldz: %d \n", n, ldz);
	#endif

	bufsize = ldz*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&z, &zref, bufsize);
	lapacke_gtest_alloc_double_buffer_pair(&d, &dref, (n));
	lapacke_gtest_alloc_double_buffer_pair(&e, &eref, ((n-1)));
	if ((z==NULL) || (zref==NULL) || \
		(d==NULL) || (dref==NULL) || \
		(e == NULL) ||(eref == NULL)){
		EXPECT_FALSE( true) << "stedc_double_parameters object: malloc error.";
		stedc_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand( z, zref, bufsize);
	lapacke_gtest_init_double_buffer_pair_rand( d, dref, (n));
	lapacke_gtest_init_double_buffer_pair_rand( e, eref, ((n-1)));	

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
stedc_double_parameters :: ~stedc_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" stedc_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   stedc_free();

}
/*  Test fixture class definition */
class dstedc_test  : public  ::testing::Test {
public:
   stedc_double_parameters  *dstedc_obj;
   void SetUp();
   void TearDown () { delete dstedc_obj; }
};

void dstedc_test::SetUp(){

    /* LAPACKE dstedc prototype */
    typedef int (*Fptr_NL_LAPACKE_dstedc) (int matrix_layout, char compz, lapack_int n, \
											double* d, double* e, double *z, lapack_int ldz);

    Fptr_NL_LAPACKE_dstedc dstedc;

    dstedc_obj = new stedc_double_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].compz,
						   eig_paramslist[idx].n,
                           eig_paramslist[idx].ldb);

    idx = Circular_Increment_Index(idx);
	

    dstedc_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dstedc_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dstedc_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dstedc_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    dstedc = (Fptr_NL_LAPACKE_dstedc)dlsym(dstedc_obj->hModule, "LAPACKE_dstedc");
    ASSERT_TRUE(dstedc != NULL) << "failed to get the Netlib LAPACKE_dstedc symbol";
    

    dstedc_obj->inforef = dstedc( dstedc_obj->matrix_layout, dstedc_obj->compz, dstedc_obj->n,\
								dstedc_obj->dref, dstedc_obj->eref, dstedc_obj->zref,\
								dstedc_obj->ldz);

    /* Compute libflame's Lapacke o/p  */
    dstedc_obj->info = LAPACKE_dstedc( dstedc_obj->matrix_layout, dstedc_obj->compz, dstedc_obj->n,\
										dstedc_obj->d, dstedc_obj->e, dstedc_obj->z,\
										dstedc_obj->ldz);

    if( dstedc_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dstedc is wrong\n", dstedc_obj->info );
    }
    if( dstedc_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dstedc is wrong\n", 
        dstedc_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dstedc_obj->diff =  computeDiff_d( dstedc_obj->bufsize, 
                dstedc_obj->z, dstedc_obj->zref );

}

TEST_F(dstedc_test, dstedc1) {
    EXPECT_NEAR(0.0, dstedc_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dstedc_test, dstedc2) {
    EXPECT_NEAR(0.0, dstedc_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dstedc_test, dstedc3) {
    EXPECT_NEAR(0.0, dstedc_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dstedc_test, dstedc4) {
    EXPECT_NEAR(0.0, dstedc_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin scomplex_common_parameters  class definition */
class stedc_scomplex_parameters{

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
      stedc_scomplex_parameters (int matrix_layout,char compz, lapack_int n, lapack_int ldz);
      ~stedc_scomplex_parameters ();

};

/* Constructor definition  scomplex_common_parameters */
stedc_scomplex_parameters:: stedc_scomplex_parameters (int matrix_layout_i, char compz_i,lapack_int n_i, lapack_int ldz_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;	
	ldz = ldz_i;
	compz = compz_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n stedc scomplex:  n: %d ldz: %d \n", n, ldz);
	#endif

	bufsize = ldz*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&z, &zref, bufsize);
	lapacke_gtest_alloc_float_buffer_pair(&d, &dref, (n));
	lapacke_gtest_alloc_float_buffer_pair(&e, &eref, ((n-1)));
	if ((z==NULL) || (zref==NULL) || \
		(d==NULL) || (dref==NULL) || \
		(e == NULL) ||(eref == NULL)){
		EXPECT_FALSE( true) << "stedc_scomplex_parameters object: malloc error.";
		stedc_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( z, zref, bufsize);
	lapacke_gtest_init_float_buffer_pair_rand( d, dref, (n));
	lapacke_gtest_init_float_buffer_pair_rand( e, eref, ((n-1)));	

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
stedc_scomplex_parameters :: ~stedc_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" stedc_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   stedc_free();

}
/*  Test fixture class definition */
class cstedc_test  : public  ::testing::Test {
public:
   stedc_scomplex_parameters  *cstedc_obj;
   void SetUp();
   void TearDown () { delete cstedc_obj; }
};

void cstedc_test::SetUp(){

    /* LAPACKE cstedc prototype */
    typedef int (*Fptr_NL_LAPACKE_cstedc) (int matrix_layout, char compz, lapack_int n, \
											float* d, float* e, lapack_complex_float *z, lapack_int ldz);

    Fptr_NL_LAPACKE_cstedc cstedc;

    cstedc_obj = new stedc_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].compz,
						   eig_paramslist[idx].n,
                           eig_paramslist[idx].ldb);

    idx = Circular_Increment_Index(idx);
	

    cstedc_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cstedc_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cstedc_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cstedc_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    cstedc = (Fptr_NL_LAPACKE_cstedc)dlsym(cstedc_obj->hModule, "LAPACKE_cstedc");
    ASSERT_TRUE(cstedc != NULL) << "failed to get the Netlib LAPACKE_cstedc symbol";
    

    cstedc_obj->inforef = cstedc( cstedc_obj->matrix_layout, cstedc_obj->compz, cstedc_obj->n,\
								cstedc_obj->dref, cstedc_obj->eref, cstedc_obj->zref,\
								cstedc_obj->ldz);

    /* Compute libflame's Lapacke o/p  */
    cstedc_obj->info = LAPACKE_cstedc( cstedc_obj->matrix_layout, cstedc_obj->compz, cstedc_obj->n,\
										cstedc_obj->d, cstedc_obj->e, cstedc_obj->z,\
										cstedc_obj->ldz);

    if( cstedc_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cstedc is wrong\n", cstedc_obj->info );
    }
    if( cstedc_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cstedc is wrong\n", 
        cstedc_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cstedc_obj->diff =  computeDiff_c( cstedc_obj->bufsize, 
                cstedc_obj->z, cstedc_obj->zref );

}

TEST_F(cstedc_test, cstedc1) {
    EXPECT_NEAR(0.0, cstedc_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cstedc_test, cstedc2) {
    EXPECT_NEAR(0.0, cstedc_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cstedc_test, cstedc3) {
    EXPECT_NEAR(0.0, cstedc_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cstedc_test, cstedc4) {
    EXPECT_NEAR(0.0, cstedc_obj->diff, LAPACKE_GTEST_THRESHOLD);
}



/* Begin dcomplex_common_parameters  class definition */
class stedc_dcomplex_parameters{

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
      stedc_dcomplex_parameters (int matrix_layout,char compz, lapack_int n, lapack_int ldz);
      ~stedc_dcomplex_parameters ();

};

/* Constructor definition  dcomplex_common_parameters */
stedc_dcomplex_parameters:: stedc_dcomplex_parameters (int matrix_layout_i, char compz_i,lapack_int n_i, lapack_int ldz_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;	
	ldz = ldz_i;
	compz = compz_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n stedc dcomplex:  n: %d ldz: %d \n", n, ldz);
	#endif

	bufsize = ldz*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&z, &zref, bufsize);
	lapacke_gtest_alloc_double_buffer_pair(&d, &dref, (n));
	lapacke_gtest_alloc_double_buffer_pair(&e, &eref, ((n-1)));
	if ((z==NULL) || (zref==NULL) || \
		(d==NULL) || (dref==NULL) || \
		(e == NULL) ||(eref == NULL)){
		EXPECT_FALSE( true) << "stedc_dcomplex_parameters object: malloc error.";
		stedc_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( z, zref, bufsize);
	lapacke_gtest_init_double_buffer_pair_rand( d, dref, (n));
	lapacke_gtest_init_double_buffer_pair_rand( e, eref, ((n-1)));	

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
stedc_dcomplex_parameters :: ~stedc_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" stedc_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   stedc_free();

}
/*  Test fixture class definition */
class zstedc_test  : public  ::testing::Test {
public:
   stedc_dcomplex_parameters  *zstedc_obj;
   void SetUp();
   void TearDown () { delete zstedc_obj; }
};

void zstedc_test::SetUp(){

    /* LAPACKE zstedc prototype */
    typedef int (*Fptr_NL_LAPACKE_zstedc) (int matrix_layout, char compz, lapack_int n, \
											double* d, double* e, lapack_complex_double *z, lapack_int ldz);

    Fptr_NL_LAPACKE_zstedc zstedc;

    zstedc_obj = new stedc_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].compz,
						   eig_paramslist[idx].n,
                           eig_paramslist[idx].ldb);

    idx = Circular_Increment_Index(idx);
	

    zstedc_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zstedc_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zstedc_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zstedc_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zstedc = (Fptr_NL_LAPACKE_zstedc)dlsym(zstedc_obj->hModule, "LAPACKE_zstedc");
    ASSERT_TRUE(zstedc != NULL) << "failed to get the Netlib LAPACKE_zstedc symbol";
    

    zstedc_obj->inforef = zstedc( zstedc_obj->matrix_layout, zstedc_obj->compz, zstedc_obj->n,\
								zstedc_obj->dref, zstedc_obj->eref, zstedc_obj->zref,\
								zstedc_obj->ldz);

    /* Compute libflame's Lapacke o/p  */
    zstedc_obj->info = LAPACKE_zstedc( zstedc_obj->matrix_layout, zstedc_obj->compz, zstedc_obj->n,\
										zstedc_obj->d, zstedc_obj->e, zstedc_obj->z,\
										zstedc_obj->ldz);

    if( zstedc_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zstedc is wrong\n", zstedc_obj->info );
    }
    if( zstedc_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zstedc is wrong\n", 
        zstedc_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zstedc_obj->diff =  computeDiff_z( zstedc_obj->bufsize, 
                zstedc_obj->z, zstedc_obj->zref );

}

TEST_F(zstedc_test, zstedc1) {
    EXPECT_NEAR(0.0, zstedc_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zstedc_test, zstedc2) {
    EXPECT_NEAR(0.0, zstedc_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zstedc_test, zstedc3) {
    EXPECT_NEAR(0.0, zstedc_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zstedc_test, zstedc4) {
    EXPECT_NEAR(0.0, zstedc_obj->diff, LAPACKE_GTEST_THRESHOLD);
}
