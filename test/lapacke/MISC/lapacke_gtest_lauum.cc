#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"


#define lauum_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
/* Begin float_common_parameters  class definition */
class lauum_float_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	char uplo;
	float* A;
	lapack_int lda;
	/*Output Parameter*/
	float *Aref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      lauum_float_parameters (int matrix_layout, char uplo, lapack_int n, lapack_int lda);
      ~lauum_float_parameters ();

};

/* Constructor definition  float_common_parameters */
lauum_float_parameters::lauum_float_parameters (int matrix_layout_i,char uplo_i, lapack_int n_i, lapack_int lda_i)
{
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n lauum float:  m: %d, n: %d lda: %d, uplo:%c\n",  m, n, lda, uplo);
	#endif

	bufsize = lda*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&A, &Aref, bufsize);	
	if ((A==NULL) || (Aref==NULL)){
		EXPECT_FALSE( true) << "lauum_float_parameters object: malloc error.";
		lauum_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand( A, Aref, bufsize);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
lauum_float_parameters :: ~lauum_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" lauum_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lauum_free();

}
/*  Test fixture class definition */
class slauum_test  : public  ::testing::Test {
public:
   lauum_float_parameters  *slauum_obj;
   void SetUp();
   void TearDown () { delete slauum_obj; }
};

void slauum_test::SetUp(){

    /* LAPACKE slauum prototype */
    typedef int (*Fptr_NL_LAPACKE_slauum) (int matrix_layout, char uplo, lapack_int n, float *A, lapack_int lda);

    Fptr_NL_LAPACKE_slauum slauum;

    slauum_obj = new lauum_float_parameters ( eig_paramslist[idx].matrix_layout,
						eig_paramslist[idx].uplo,
						eig_paramslist[idx].n,
						eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);

    slauum_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    slauum_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(slauum_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(slauum_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    slauum = (Fptr_NL_LAPACKE_slauum)dlsym(slauum_obj->hModule, "LAPACKE_slauum");
    ASSERT_TRUE(slauum != NULL) << "failed to get the Netlib LAPACKE_slauum symbol";
    

    slauum_obj->inforef = slauum( slauum_obj->matrix_layout, slauum_obj->uplo,
								slauum_obj->n, slauum_obj->Aref,
								slauum_obj->lda);

    /* Compute libflame's Lapacke o/p  */
    slauum_obj->info = LAPACKE_slauum( slauum_obj->matrix_layout, slauum_obj->uplo,
										slauum_obj->n, slauum_obj->A, 
										slauum_obj->lda);

    if( slauum_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_slauum is wrong\n", slauum_obj->info );
    }
    if( slauum_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_slauum is wrong\n", 
        slauum_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    slauum_obj->diff =  computeDiff_s( slauum_obj->bufsize, 
                slauum_obj->A, slauum_obj->Aref );

}

TEST_F(slauum_test, slauum1) {
    EXPECT_NEAR(0.0, slauum_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(slauum_test, slauum2) {
    EXPECT_NEAR(0.0, slauum_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(slauum_test, slauum3) {
    EXPECT_NEAR(0.0, slauum_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(slauum_test, slauum4) {
    EXPECT_NEAR(0.0, slauum_obj->diff, LAPACKE_EIG_THRESHOLD);
}

/* Begin double_common_parameters  class definition */
class lauum_double_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	char uplo;
	double* A;
	lapack_int lda;
	/*Output Parameter*/
	double *Aref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      lauum_double_parameters (int matrix_layout, char uplo, lapack_int n, lapack_int lda);
      ~lauum_double_parameters ();

};

/* Constructor definition  double_common_parameters */
lauum_double_parameters::lauum_double_parameters (int matrix_layout_i,char uplo_i, lapack_int n_i, lapack_int lda_i)
{
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n lauum double:  m: %d, n: %d lda: %d, uplo:%c\n",  m, n, lda, uplo);
	#endif

	bufsize = lda*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&A, &Aref, bufsize);	
	if ((A==NULL) || (Aref==NULL)){
		EXPECT_FALSE( true) << "lauum_double_parameters object: malloc error.";
		lauum_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand( A, Aref, bufsize);

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
lauum_double_parameters :: ~lauum_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" lauum_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lauum_free();

}
/*  Test fixture class definition */
class dlauum_test  : public  ::testing::Test {
public:
   lauum_double_parameters  *dlauum_obj;
   void SetUp();
   void TearDown () { delete dlauum_obj; }
};

void dlauum_test::SetUp(){

    /* LAPACKE dlauum prototype */
    typedef int (*Fptr_NL_LAPACKE_dlauum) (int matrix_layout, char uplo, lapack_int n, double *A, lapack_int lda);

    Fptr_NL_LAPACKE_dlauum dlauum;

    dlauum_obj = new lauum_double_parameters ( eig_paramslist[idx].matrix_layout,
						eig_paramslist[idx].uplo,
						eig_paramslist[idx].n,
						eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);

    dlauum_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dlauum_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dlauum_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dlauum_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    dlauum = (Fptr_NL_LAPACKE_dlauum)dlsym(dlauum_obj->hModule, "LAPACKE_dlauum");
    ASSERT_TRUE(dlauum != NULL) << "failed to get the Netlib LAPACKE_dlauum symbol";
    

    dlauum_obj->inforef = dlauum( dlauum_obj->matrix_layout, dlauum_obj->uplo,
								dlauum_obj->n, dlauum_obj->Aref,
								dlauum_obj->lda);

    /* Compute libflame's Lapacke o/p  */
    dlauum_obj->info = LAPACKE_dlauum( dlauum_obj->matrix_layout, dlauum_obj->uplo,
										dlauum_obj->n, dlauum_obj->A, 
										dlauum_obj->lda);

    if( dlauum_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dlauum is wrong\n", dlauum_obj->info );
    }
    if( dlauum_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dlauum is wrong\n", 
        dlauum_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dlauum_obj->diff =  computeDiff_d( dlauum_obj->bufsize, 
                dlauum_obj->A, dlauum_obj->Aref );

}

TEST_F(dlauum_test, dlauum1) {
    EXPECT_NEAR(0.0, dlauum_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(dlauum_test, dlauum2) {
    EXPECT_NEAR(0.0, dlauum_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(dlauum_test, dlauum3) {
    EXPECT_NEAR(0.0, dlauum_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(dlauum_test, dlauum4) {
    EXPECT_NEAR(0.0, dlauum_obj->diff, LAPACKE_EIG_THRESHOLD);
}



/* Begin scomplex_common_parameters  class definition */
class lauum_scomplex_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	char uplo;
	lapack_complex_float* A;
	lapack_int lda;
	/*Output Parameter*/
	lapack_complex_float *Aref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      lauum_scomplex_parameters (int matrix_layout, char uplo, lapack_int n, lapack_int lda);
      ~lauum_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
lauum_scomplex_parameters::lauum_scomplex_parameters (int matrix_layout_i,char uplo_i, lapack_int n_i, lapack_int lda_i)
{
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n lauum scomplex:  m: %d, n: %d lda: %d, uplo:%c\n",  m, n, lda, uplo);
	#endif

	bufsize = lda*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize);	
	if ((A==NULL) || (Aref==NULL)){
		EXPECT_FALSE( true) << "lauum_float_parameters object: malloc error.";
		lauum_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( A, Aref, bufsize);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
lauum_scomplex_parameters :: ~lauum_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" lauum_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lauum_free();

}
/*  Test fixture class definition */
class clauum_test  : public  ::testing::Test {
public:
   lauum_scomplex_parameters  *clauum_obj;
   void SetUp();
   void TearDown () { delete clauum_obj; }
};

void clauum_test::SetUp(){

    /* LAPACKE clauum prototype */
    typedef int (*Fptr_NL_LAPACKE_clauum) (int matrix_layout, char uplo, lapack_int n, lapack_complex_float *A, lapack_int lda);

    Fptr_NL_LAPACKE_clauum clauum;

    clauum_obj = new lauum_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
						eig_paramslist[idx].uplo,
						eig_paramslist[idx].n,
						eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);

    clauum_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    clauum_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(clauum_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(clauum_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    clauum = (Fptr_NL_LAPACKE_clauum)dlsym(clauum_obj->hModule, "LAPACKE_clauum");
    ASSERT_TRUE(clauum != NULL) << "failed to get the Netlib LAPACKE_clauum symbol";
		

    clauum_obj->inforef = clauum( clauum_obj->matrix_layout, clauum_obj->uplo,
								clauum_obj->n, clauum_obj->Aref,
								clauum_obj->lda);

    /* Compute libflame's Lapacke o/p  */
    clauum_obj->info = LAPACKE_clauum( clauum_obj->matrix_layout, clauum_obj->uplo,
										clauum_obj->n, clauum_obj->A, 
										clauum_obj->lda);

    if( clauum_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_clauum is wrong\n", clauum_obj->info );
    }
    if( clauum_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_clauum is wrong\n", 
        clauum_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    clauum_obj->diff =  computeDiff_c( clauum_obj->bufsize, 
                clauum_obj->A, clauum_obj->Aref );

}

TEST_F(clauum_test, clauum1) {
    EXPECT_NEAR(0.0, clauum_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(clauum_test, clauum2) {
    EXPECT_NEAR(0.0, clauum_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(clauum_test, clauum3) {
    EXPECT_NEAR(0.0, clauum_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(clauum_test, clauum4) {
    EXPECT_NEAR(0.0, clauum_obj->diff, LAPACKE_EIG_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class lauum_dcomplex_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	char uplo;
	lapack_complex_double* A;
	lapack_int lda;
	/*Output Parameter*/
	lapack_complex_double *Aref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      lauum_dcomplex_parameters (int matrix_layout, char uplo, lapack_int n, lapack_int lda);
      ~lauum_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
lauum_dcomplex_parameters::lauum_dcomplex_parameters (int matrix_layout_i,char uplo_i, lapack_int n_i, lapack_int lda_i)
{
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n lauum dcomplex:   n: %d lda: %d, uplo:%c\n",   n, lda, uplo);
	#endif

	bufsize = lda*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize);	
	if ((A==NULL) || (Aref==NULL)){
		EXPECT_FALSE( true) << "lauum_double_parameters object: malloc error.";
		lauum_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( A, Aref, bufsize);

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
lauum_dcomplex_parameters :: ~lauum_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" lauum_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lauum_free();

}
/*  Test fixture class definition */
class zlauum_test  : public  ::testing::Test {
public:
   lauum_dcomplex_parameters  *zlauum_obj;
   void SetUp();
   void TearDown () { delete zlauum_obj; }
};

void zlauum_test::SetUp(){

    /* LAPACKE zlauum prototype */
    typedef int (*Fptr_NL_LAPACKE_zlauum) (int matrix_layout, char uplo, lapack_int n, lapack_complex_double *A, lapack_int lda);

    Fptr_NL_LAPACKE_zlauum zlauum;

    zlauum_obj = new lauum_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
						eig_paramslist[idx].uplo,
						eig_paramslist[idx].n,
						eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);

    zlauum_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zlauum_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zlauum_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zlauum_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zlauum = (Fptr_NL_LAPACKE_zlauum)dlsym(zlauum_obj->hModule, "LAPACKE_zlauum");
    ASSERT_TRUE(zlauum != NULL) << "failed to get the Netlib LAPACKE_zlauum symbol";
    

    zlauum_obj->inforef = zlauum( zlauum_obj->matrix_layout, zlauum_obj->uplo,
								zlauum_obj->n, zlauum_obj->Aref,
								zlauum_obj->lda);

    /* Compute libflame's Lapacke o/p  */
    zlauum_obj->info = LAPACKE_zlauum( zlauum_obj->matrix_layout, zlauum_obj->uplo,
										zlauum_obj->n, zlauum_obj->A, 
										zlauum_obj->lda);

    if( zlauum_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zlauum is wrong\n", zlauum_obj->info );
    }
    if( zlauum_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zlauum is wrong\n", 
        zlauum_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zlauum_obj->diff =  computeDiff_z( zlauum_obj->bufsize, 
                zlauum_obj->A, zlauum_obj->Aref );

}

TEST_F(zlauum_test, zlauum1) {
    EXPECT_NEAR(0.0, zlauum_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zlauum_test, zlauum2) {
    EXPECT_NEAR(0.0, zlauum_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zlauum_test, zlauum3) {
    EXPECT_NEAR(0.0, zlauum_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zlauum_test, zlauum4) {
    EXPECT_NEAR(0.0, zlauum_obj->diff, LAPACKE_EIG_THRESHOLD);
}
