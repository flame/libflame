#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"


#define lantr_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
class lantr_float_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	char norm;
	char uplo;
	char diag;
	float* A;
	lapack_int lda;
	/*Output Parameter*/
	float *Aref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      lantr_float_parameters (int matrix_layout, char norm, char uplo,  char diag, lapack_int m ,lapack_int n);
      ~lantr_float_parameters ();

};

/* Constructor definition  float_common_parameters */
lantr_float_parameters::lantr_float_parameters (int matrix_layout_i, char norm_i, char uplo_i, char diag_i, lapack_int m_i , lapack_int n_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	norm = norm_i;
	uplo = uplo_i;
	diag = diag_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n lantr float:   n: %d lda: %d, uplo:%c, ldb:%d\n",  n, lda, norm);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{	lda = m;
		bufsize = lda*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{	lda = n;
		bufsize = lda*m;
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";

	
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&A, &Aref, bufsize);
	if ((A==NULL) || (Aref==NULL)){
		EXPECT_FALSE( true) << "lantr_float_parameters object: malloc error.";
		lantr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand( A, Aref, bufsize);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
lantr_float_parameters :: ~lantr_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" lantr_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lantr_free();

}
/*  Test fixture class definition */
class slantr_test  : public  ::testing::Test {
public:
   lantr_float_parameters  *slantr_obj;
   void SetUp();
   void TearDown () { delete slantr_obj; }
};

void slantr_test::SetUp(){

    /* LAPACKE slantr prototype */
    typedef int (*Fptr_NL_LAPACKE_slantr) (int matrix_layout, char norm, char uplo,  char diag, lapack_int m , lapack_int n, const float *A, lapack_int lda);

    Fptr_NL_LAPACKE_slantr slantr;

    slantr_obj = new lantr_float_parameters ( eig_paramslist[idx].matrix_layout,
						eig_paramslist[idx].norm,
						eig_paramslist[idx].uplo,
						eig_paramslist[idx].diag,
						eig_paramslist[idx].m,
						eig_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    slantr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    slantr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(slantr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(slantr_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    slantr = (Fptr_NL_LAPACKE_slantr)dlsym(slantr_obj->hModule, "LAPACKE_slantr");
    ASSERT_TRUE(slantr != NULL) << "failed to get the Netlib LAPACKE_slantr symbol";
    

    slantr_obj->inforef = slantr( slantr_obj->matrix_layout, slantr_obj->norm, slantr_obj->uplo, slantr_obj->diag, 
								slantr_obj->m, slantr_obj->n, (const float *)slantr_obj->Aref,
								slantr_obj->lda);

    /* Compute libflame's Lapacke o/p  */
    slantr_obj->info = LAPACKE_slantr( slantr_obj->matrix_layout, slantr_obj->norm, slantr_obj->uplo,slantr_obj->diag,
										slantr_obj->m, slantr_obj->n, (const float*)slantr_obj->A, 
										slantr_obj->lda);

    if( slantr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_slantr is wrong\n", slantr_obj->info );
    }
    if( slantr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_slantr is wrong\n", 
        slantr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    slantr_obj->diff =  computeDiff_s( slantr_obj->bufsize, 
                slantr_obj->A, slantr_obj->Aref );

}

TEST_F(slantr_test, slantr1) {
    EXPECT_NEAR(0.0, slantr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(slantr_test, slantr2) {
    EXPECT_NEAR(0.0, slantr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(slantr_test, slantr3) {
    EXPECT_NEAR(0.0, slantr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(slantr_test, slantr4) {
    EXPECT_NEAR(0.0, slantr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin double_common_parameters  class definition */
class lantr_double_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	char norm;
	char uplo;
	char diag;
	double* A;
	lapack_int lda;
	/*Output Parameter*/
	double *Aref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      lantr_double_parameters (int matrix_layout, char norm, char uplo,  char diag, lapack_int m ,lapack_int n);
      ~lantr_double_parameters ();

};

/* Constructor definition  double_common_parameters */
lantr_double_parameters::lantr_double_parameters (int matrix_layout_i, char norm_i, char uplo_i, char diag_i, lapack_int m_i , lapack_int n_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	norm = norm_i;
	uplo = uplo_i;
	diag = diag_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n lantr double:   n: %d lda: %d, uplo:%c, ldb:%d\n",  n, lda, norm);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{	lda = m;
		bufsize = lda*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{	lda = n;
		bufsize = lda*m;
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";

	
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&A, &Aref, bufsize);
	if ((A==NULL) || (Aref==NULL)){
		EXPECT_FALSE( true) << "lantr_double_parameters object: malloc error.";
		lantr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand( A, Aref, bufsize);

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
lantr_double_parameters :: ~lantr_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" lantr_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lantr_free();

}
/*  Test fixture class definition */
class dlantr_test  : public  ::testing::Test {
public:
   lantr_double_parameters  *dlantr_obj;
   void SetUp();
   void TearDown () { delete dlantr_obj; }
};

void dlantr_test::SetUp(){

    /* LAPACKE dlantr prototype */
    typedef int (*Fptr_NL_LAPACKE_dlantr) (int matrix_layout, char norm, char uplo,  char diag, lapack_int m , lapack_int n, const double *A, lapack_int lda);

    Fptr_NL_LAPACKE_dlantr dlantr;

    dlantr_obj = new lantr_double_parameters ( eig_paramslist[idx].matrix_layout,
						eig_paramslist[idx].norm,
						eig_paramslist[idx].uplo,
						eig_paramslist[idx].diag,
						eig_paramslist[idx].m,
						eig_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    dlantr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dlantr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dlantr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dlantr_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    dlantr = (Fptr_NL_LAPACKE_dlantr)dlsym(dlantr_obj->hModule, "LAPACKE_dlantr");
    ASSERT_TRUE(dlantr != NULL) << "failed to get the Netlib LAPACKE_dlantr symbol";
    

    dlantr_obj->inforef = dlantr( dlantr_obj->matrix_layout, dlantr_obj->norm, dlantr_obj->uplo, dlantr_obj->diag, 
								dlantr_obj->m, dlantr_obj->n, (const double *)dlantr_obj->Aref,
								dlantr_obj->lda);

    /* Compute libflame's Lapacke o/p  */
    dlantr_obj->info = LAPACKE_dlantr( dlantr_obj->matrix_layout, dlantr_obj->norm, dlantr_obj->uplo,dlantr_obj->diag,
										dlantr_obj->m, dlantr_obj->n, (const double*)dlantr_obj->A, 
										dlantr_obj->lda);

    if( dlantr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dlantr is wrong\n", dlantr_obj->info );
    }
    if( dlantr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dlantr is wrong\n", 
        dlantr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dlantr_obj->diff =  computeDiff_d( dlantr_obj->bufsize, 
                dlantr_obj->A, dlantr_obj->Aref );

}

TEST_F(dlantr_test, dlantr1) {
    EXPECT_NEAR(0.0, dlantr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dlantr_test, dlantr2) {
    EXPECT_NEAR(0.0, dlantr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dlantr_test, dlantr3) {
    EXPECT_NEAR(0.0, dlantr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dlantr_test, dlantr4) {
    EXPECT_NEAR(0.0, dlantr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}



/* Begin scomplex_common_parameters  class definition */
class lantr_scomplex_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	char uplo;
	char norm;
	char diag;
	lapack_complex_float* A;
	lapack_int lda;
	/*Output Parameter*/
	lapack_complex_float *Aref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      lantr_scomplex_parameters (int matrix_layout, char norm, char uplo, char diag, lapack_int m , lapack_int n);
      ~lantr_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
lantr_scomplex_parameters::lantr_scomplex_parameters (int matrix_layout_i,char norm_i, char uplo_i, char diag_i, lapack_int m_i ,lapack_int n_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	norm = norm_i;
	uplo = uplo_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n lantr scomplex:  m: %d, n: %d, uplo:%c, ldb:%d\n", n,  norm );
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{	lda = m;
		bufsize = lda*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{	lda = n;
		bufsize = lda*m;
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";

	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize);
	if ((A==NULL) || (Aref==NULL)){
		EXPECT_FALSE( true) << "lantr_float_parameters object: malloc error.";
		lantr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( A, Aref, bufsize);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
lantr_scomplex_parameters :: ~lantr_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" lantr_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lantr_free();

}
/*  Test fixture class definition */
class clantr_test  : public  ::testing::Test {
public:
   lantr_scomplex_parameters  *clantr_obj;
   void SetUp();
   void TearDown () { delete clantr_obj; }
};

void clantr_test::SetUp(){

    /* LAPACKE clantr prototype */
    typedef int (*Fptr_NL_LAPACKE_clantr) (int matrix_layout, char norm, char uplo, char diag, lapack_int m , lapack_int n, const lapack_complex_float *A, lapack_int lda);

    Fptr_NL_LAPACKE_clantr clantr;

    clantr_obj = new lantr_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
						eig_paramslist[idx].norm,
						eig_paramslist[idx].uplo,
						eig_paramslist[idx].diag,
						eig_paramslist[idx].m,
						eig_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    clantr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    clantr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(clantr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(clantr_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    clantr = (Fptr_NL_LAPACKE_clantr)dlsym(clantr_obj->hModule, "LAPACKE_clantr");
    ASSERT_TRUE(clantr != NULL) << "failed to get the Netlib LAPACKE_clantr symbol";
    

    clantr_obj->inforef = clantr( clantr_obj->matrix_layout,clantr_obj->norm, clantr_obj->uplo, clantr_obj->diag,
								clantr_obj->m, clantr_obj->n, (const lapack_complex_float *)clantr_obj->Aref,
								clantr_obj->lda);

    /* Compute libflame's Lapacke o/p  */
    clantr_obj->info = LAPACKE_clantr( clantr_obj->matrix_layout, clantr_obj->norm, clantr_obj->uplo, clantr_obj->diag,
										clantr_obj->m, clantr_obj->n, (const lapack_complex_float*)clantr_obj->A, 
										clantr_obj->lda);

    if( clantr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_clantr is wrong\n", clantr_obj->info );
    }
    if( clantr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_clantr is wrong\n", 
        clantr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    clantr_obj->diff =  computeDiff_c( clantr_obj->bufsize, 
                clantr_obj->A, clantr_obj->Aref );

}

TEST_F(clantr_test, clantr1) {
    EXPECT_NEAR(0.0, clantr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(clantr_test, clantr2) {
    EXPECT_NEAR(0.0, clantr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(clantr_test, clantr3) {
    EXPECT_NEAR(0.0, clantr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(clantr_test, clantr4) {
    EXPECT_NEAR(0.0, clantr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class lantr_dcomplex_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	char uplo;
	char norm;
	char diag;
	lapack_complex_double* A;
	lapack_int lda;
	/*Output Parameter*/
	lapack_complex_double *Aref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      lantr_dcomplex_parameters (int matrix_layout, char norm, char uplo, char diag, lapack_int m , lapack_int n);
      ~lantr_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
lantr_dcomplex_parameters::lantr_dcomplex_parameters (int matrix_layout_i,char norm_i, char uplo_i, char diag_i, lapack_int m_i ,lapack_int n_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	norm = norm_i;
	uplo = uplo_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n lantr scomplex:  m: %d, n: %d, uplo:%c, ldb:%d\n", n,  norm );
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{	lda = m;
		bufsize = lda*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{	lda = n;
		bufsize = lda*m;
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";

	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize);
	if ((A==NULL) || (Aref==NULL)){
		EXPECT_FALSE( true) << "lantr_double_parameters object: malloc error.";
		lantr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( A, Aref, bufsize);

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
lantr_dcomplex_parameters :: ~lantr_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" lantr_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lantr_free();

}
/*  Test fixture class definition */
class zlantr_test  : public  ::testing::Test {
public:
   lantr_dcomplex_parameters  *zlantr_obj;
   void SetUp();
   void TearDown () { delete zlantr_obj; }
};

void zlantr_test::SetUp(){

    /* LAPACKE zlantr prototype */
    typedef int (*Fptr_NL_LAPACKE_zlantr) (int matrix_layout, char norm, char uplo, char diag, lapack_int m , lapack_int n, const lapack_complex_double *A, lapack_int lda);

    Fptr_NL_LAPACKE_zlantr zlantr;

    zlantr_obj = new lantr_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
						eig_paramslist[idx].norm,
						eig_paramslist[idx].uplo,
						eig_paramslist[idx].diag,
						eig_paramslist[idx].m,
						eig_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    zlantr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zlantr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zlantr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zlantr_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zlantr = (Fptr_NL_LAPACKE_zlantr)dlsym(zlantr_obj->hModule, "LAPACKE_zlantr");
    ASSERT_TRUE(zlantr != NULL) << "failed to get the Netlib LAPACKE_zlantr symbol";
    

    zlantr_obj->inforef = zlantr( zlantr_obj->matrix_layout,zlantr_obj->norm, zlantr_obj->uplo, zlantr_obj->diag,
								zlantr_obj->m, zlantr_obj->n, (const lapack_complex_double *)zlantr_obj->Aref,
								zlantr_obj->lda);

    /* Compute libflame's Lapacke o/p  */
    zlantr_obj->info = LAPACKE_zlantr( zlantr_obj->matrix_layout, zlantr_obj->norm, zlantr_obj->uplo, zlantr_obj->diag,
										zlantr_obj->m, zlantr_obj->n, (const lapack_complex_double*)zlantr_obj->A, 
										zlantr_obj->lda);

    if( zlantr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zlantr is wrong\n", zlantr_obj->info );
    }
    if( zlantr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zlantr is wrong\n", 
        zlantr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zlantr_obj->diff =  computeDiff_z( zlantr_obj->bufsize, 
                zlantr_obj->A, zlantr_obj->Aref );

}

TEST_F(zlantr_test, zlantr1) {
    EXPECT_NEAR(0.0, zlantr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlantr_test, zlantr2) {
    EXPECT_NEAR(0.0, zlantr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlantr_test, zlantr3) {
    EXPECT_NEAR(0.0, zlantr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlantr_test, zlantr4) {
    EXPECT_NEAR(0.0, zlantr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}
