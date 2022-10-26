#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"


#define lansy_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
class lansy_float_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	char norm;
	char uplo;
	float* A;
	lapack_int lda;
	/*Output Parameter*/
	float *Aref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      lansy_float_parameters (int matrix_layout, char uplo, char norm, lapack_int n, lapack_int lda);
      ~lansy_float_parameters ();

};

/* Constructor definition  float_common_parameters */
lansy_float_parameters::lansy_float_parameters (int matrix_layout_i, char uplo_i, char norm_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	norm = norm_i;
	lda = lda_i;
	uplo = uplo_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n lansy float:   n: %d lda: %d, uplo:%c, ldb:%d\n",  n, lda, norm);
	#endif
	bufsize = lda*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&A, &Aref, bufsize);
	if ((A==NULL) || (Aref==NULL)){
		EXPECT_FALSE( true) << "lansy_float_parameters object: malloc error.";
		lansy_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand( A, Aref, bufsize);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
lansy_float_parameters :: ~lansy_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" lansy_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lansy_free();

}
/*  Test fixture class definition */
class slansy_test  : public  ::testing::Test {
public:
   lansy_float_parameters  *slansy_obj;
   void SetUp();
   void TearDown () { delete slansy_obj; }
};

void slansy_test::SetUp(){

    /* LAPACKE slansy prototype */
    typedef int (*Fptr_NL_LAPACKE_slansy) (int matrix_layout, char uplo, char norm, lapack_int n, const float *A, lapack_int lda);

    Fptr_NL_LAPACKE_slansy slansy;

    slansy_obj = new lansy_float_parameters ( eig_paramslist[idx].matrix_layout,
						eig_paramslist[idx].uplo,
						eig_paramslist[idx].norm,
						eig_paramslist[idx].n,
						eig_paramslist[idx].ldb);

    idx = Circular_Increment_Index(idx);

    slansy_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    slansy_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(slansy_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(slansy_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    slansy = (Fptr_NL_LAPACKE_slansy)dlsym(slansy_obj->hModule, "LAPACKE_slansy");
    ASSERT_TRUE(slansy != NULL) << "failed to get the Netlib LAPACKE_slansy symbol";
    

    slansy_obj->inforef = slansy( slansy_obj->matrix_layout, slansy_obj->uplo, slansy_obj->norm,
								slansy_obj->n, (const float *)slansy_obj->Aref,
								slansy_obj->lda);

    /* Compute libflame's Lapacke o/p  */
    slansy_obj->info = LAPACKE_slansy( slansy_obj->matrix_layout, slansy_obj->uplo, slansy_obj->norm,
										slansy_obj->n, (const float*)slansy_obj->A, 
										slansy_obj->lda);

    if( slansy_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_slansy is wrong\n", slansy_obj->info );
    }
    if( slansy_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_slansy is wrong\n", 
        slansy_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    slansy_obj->diff =  computeDiff_s( slansy_obj->bufsize, 
                slansy_obj->A, slansy_obj->Aref );

}

TEST_F(slansy_test, slansy1) {
    EXPECT_NEAR(0.0, slansy_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(slansy_test, slansy2) {
    EXPECT_NEAR(0.0, slansy_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(slansy_test, slansy3) {
    EXPECT_NEAR(0.0, slansy_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(slansy_test, slansy4) {
    EXPECT_NEAR(0.0, slansy_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin double_common_parameters  class definition */
class lansy_double_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	char norm;
	char uplo;
	double* A;
	lapack_int lda;
	/*Output Parameter*/
	double *Aref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      lansy_double_parameters (int matrix_layout, char uplo, char norm, lapack_int n, lapack_int lda);
      ~lansy_double_parameters ();

};

/* Constructor definition  double_common_parameters */
lansy_double_parameters::lansy_double_parameters (int matrix_layout_i,char uplo_i, char norm_i,  lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	norm = norm_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n lansy double:  m: %d, n: %d lda: %d, uplo:%c, ldb:%d\n",  m, n, lda, norm);
	#endif
	
	bufsize = lda*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&A, &Aref, bufsize);
	if ((A==NULL) || (Aref==NULL)){
		EXPECT_FALSE( true) << "lansy_double_parameters object: malloc error.";
		lansy_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand( A, Aref, bufsize);

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
lansy_double_parameters :: ~lansy_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" lansy_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lansy_free();

}
/*  Test fixture class definition */
class dlansy_test  : public  ::testing::Test {
public:
   lansy_double_parameters  *dlansy_obj;
   void SetUp();
   void TearDown () { delete dlansy_obj; }
};

void dlansy_test::SetUp(){

    /* LAPACKE dlansy prototype */
    typedef int (*Fptr_NL_LAPACKE_dlansy) (int matrix_layout, char uplo, char norm, lapack_int n, const double *A, lapack_int lda);

    Fptr_NL_LAPACKE_dlansy dlansy;

    dlansy_obj = new lansy_double_parameters ( eig_paramslist[idx].matrix_layout,
						eig_paramslist[idx].uplo,
						eig_paramslist[idx].norm,
						eig_paramslist[idx].n,
						eig_paramslist[idx].ldb);

    idx = Circular_Increment_Index(idx);

    dlansy_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dlansy_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dlansy_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dlansy_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    dlansy = (Fptr_NL_LAPACKE_dlansy)dlsym(dlansy_obj->hModule, "LAPACKE_dlansy");
    ASSERT_TRUE(dlansy != NULL) << "failed to get the Netlib LAPACKE_dlansy symbol";
    

    dlansy_obj->inforef = dlansy( dlansy_obj->matrix_layout, dlansy_obj->uplo,dlansy_obj->norm,
								dlansy_obj->n, (const double *)dlansy_obj->Aref,
								dlansy_obj->lda);

    /* Compute libflame's Lapacke o/p  */
    dlansy_obj->info = LAPACKE_dlansy( dlansy_obj->matrix_layout, dlansy_obj->uplo, dlansy_obj->norm,
										 dlansy_obj->n, (const double*)dlansy_obj->A, 
										dlansy_obj->lda);

    if( dlansy_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dlansy is wrong\n", dlansy_obj->info );
    }
    if( dlansy_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dlansy is wrong\n", 
        dlansy_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dlansy_obj->diff =  computeDiff_d( dlansy_obj->bufsize, 
                dlansy_obj->A, dlansy_obj->Aref );

}

TEST_F(dlansy_test, dlansy1) {
    EXPECT_NEAR(0.0, dlansy_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dlansy_test, dlansy2) {
    EXPECT_NEAR(0.0, dlansy_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dlansy_test, dlansy3) {
    EXPECT_NEAR(0.0, dlansy_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dlansy_test, dlansy4) {
    EXPECT_NEAR(0.0, dlansy_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin scomplex_common_parameters  class definition */
class lansy_scomplex_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	char uplo;
	char norm;
	lapack_complex_float* A;
	lapack_int lda;
	/*Output Parameter*/
	lapack_complex_float *Aref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      lansy_scomplex_parameters (int matrix_layout, char uplo, char norm, lapack_int n, lapack_int lda);
      ~lansy_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
lansy_scomplex_parameters::lansy_scomplex_parameters (int matrix_layout_i,char uplo_i, char norm_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	norm = norm_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n lansy scomplex:  m: %d, n: %d lda: %d, uplo:%c, ldb:%d\n", n, lda, norm );
	#endif
	
	bufsize = lda *n;
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize);
	if ((A==NULL) || (Aref==NULL)){
		EXPECT_FALSE( true) << "lansy_float_parameters object: malloc error.";
		lansy_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( A, Aref, bufsize);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
lansy_scomplex_parameters :: ~lansy_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" lansy_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lansy_free();

}
/*  Test fixture class definition */
class clansy_test  : public  ::testing::Test {
public:
   lansy_scomplex_parameters  *clansy_obj;
   void SetUp();
   void TearDown () { delete clansy_obj; }
};

void clansy_test::SetUp(){

    /* LAPACKE clansy prototype */
    typedef int (*Fptr_NL_LAPACKE_clansy) (int matrix_layout, char uplo, char norm, lapack_int n, const lapack_complex_float *A, lapack_int lda);

    Fptr_NL_LAPACKE_clansy clansy;

    clansy_obj = new lansy_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
						eig_paramslist[idx].uplo,
						eig_paramslist[idx].norm,
						eig_paramslist[idx].n,
						eig_paramslist[idx].ldb);

    idx = Circular_Increment_Index(idx);

    clansy_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    clansy_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(clansy_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(clansy_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    clansy = (Fptr_NL_LAPACKE_clansy)dlsym(clansy_obj->hModule, "LAPACKE_clansy");
    ASSERT_TRUE(clansy != NULL) << "failed to get the Netlib LAPACKE_clansy symbol";
    

    clansy_obj->inforef = clansy( clansy_obj->matrix_layout,clansy_obj->uplo, clansy_obj->norm,
								clansy_obj->n, (const lapack_complex_float *)clansy_obj->Aref,
								clansy_obj->lda);

    /* Compute libflame's Lapacke o/p  */
    clansy_obj->info = LAPACKE_clansy( clansy_obj->matrix_layout, clansy_obj->uplo,clansy_obj->norm,
										clansy_obj->n, (const lapack_complex_float*)clansy_obj->A, 
										clansy_obj->lda);

    if( clansy_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_clansy is wrong\n", clansy_obj->info );
    }
    if( clansy_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_clansy is wrong\n", 
        clansy_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    clansy_obj->diff =  computeDiff_c( clansy_obj->bufsize, 
                clansy_obj->A, clansy_obj->Aref );

}

TEST_F(clansy_test, clansy1) {
    EXPECT_NEAR(0.0, clansy_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(clansy_test, clansy2) {
    EXPECT_NEAR(0.0, clansy_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(clansy_test, clansy3) {
    EXPECT_NEAR(0.0, clansy_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(clansy_test, clansy4) {
    EXPECT_NEAR(0.0, clansy_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class lansy_dcomplex_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	char norm;
	char uplo;
	lapack_complex_double* A;
	lapack_int lda;
	/*Output Parameter*/
	lapack_complex_double *Aref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      lansy_dcomplex_parameters (int matrix_layout, char uplo, char norm, lapack_int n, lapack_int lda);
      ~lansy_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
lansy_dcomplex_parameters::lansy_dcomplex_parameters (int matrix_layout_i, char uplo_i, char norm_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	norm = norm_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n lansy dcomplex:  m: %d, n: %d lda: %d, uplo:%c, ldb:%d\n",  m, n, lda, norm );
	#endif
	
	bufsize = lda *n;
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize);
	if ((A==NULL) || (Aref==NULL)){
		EXPECT_FALSE( true) << "lansy_double_parameters object: malloc error.";
		lansy_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( A, Aref, bufsize);

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
lansy_dcomplex_parameters :: ~lansy_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" lansy_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lansy_free();

}
/*  Test fixture class definition */
class zlansy_test  : public  ::testing::Test {
public:
   lansy_dcomplex_parameters  *zlansy_obj;
   void SetUp();
   void TearDown () { delete zlansy_obj; }
};

void zlansy_test::SetUp(){

    /* LAPACKE zlansy prototype */
    typedef int (*Fptr_NL_LAPACKE_zlansy) (int matrix_layout, char uplo, char norm, lapack_int n, const lapack_complex_double *A, lapack_int lda);

    Fptr_NL_LAPACKE_zlansy zlansy;

    zlansy_obj = new lansy_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
						eig_paramslist[idx].uplo,
						eig_paramslist[idx].norm,
						eig_paramslist[idx].n,
						eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);

    zlansy_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zlansy_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zlansy_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zlansy_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zlansy = (Fptr_NL_LAPACKE_zlansy)dlsym(zlansy_obj->hModule, "LAPACKE_zlansy");
    ASSERT_TRUE(zlansy != NULL) << "failed to get the Netlib LAPACKE_zlansy symbol";
    

    zlansy_obj->inforef = zlansy( zlansy_obj->matrix_layout, zlansy_obj->uplo, zlansy_obj->norm,
								zlansy_obj->n, (const lapack_complex_double *)zlansy_obj->Aref,
								zlansy_obj->lda);

    /* Compute libflame's Lapacke o/p  */
    zlansy_obj->info = LAPACKE_zlansy( zlansy_obj->matrix_layout, zlansy_obj->uplo, zlansy_obj->norm,
										zlansy_obj->n, (const lapack_complex_double*)zlansy_obj->A, 
										zlansy_obj->lda);

    if( zlansy_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zlansy is wrong\n", zlansy_obj->info );
    }
    if( zlansy_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zlansy is wrong\n", 
        zlansy_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zlansy_obj->diff =  computeDiff_z( zlansy_obj->bufsize, 
                zlansy_obj->A, zlansy_obj->Aref );

}

TEST_F(zlansy_test, zlansy1) {
    EXPECT_NEAR(0.0, zlansy_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlansy_test, zlansy2) {
    EXPECT_NEAR(0.0, zlansy_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlansy_test, zlansy3) {
    EXPECT_NEAR(0.0, zlansy_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlansy_test, zlansy4) {
    EXPECT_NEAR(0.0, zlansy_obj->diff, LAPACKE_GTEST_THRESHOLD);
}
