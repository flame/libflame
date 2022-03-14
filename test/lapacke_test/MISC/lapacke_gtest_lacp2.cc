#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"


#define lacp2_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (b!=NULL)  free(b);\
if (bref!=NULL) free(bref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin scomplex_common_parameters  class definition */
class lacp2_scomplex_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	char uplo;
	float* A;
	lapack_int lda, ldb;
	/*Output Parameter*/
	lapack_complex_float* b;
	float *Aref;
	lapack_complex_float *bref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      lacp2_scomplex_parameters (int matrix_layout, char uplo, lapack_int m , lapack_int n, lapack_int lda, lapack_int ldb);
      ~lacp2_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
lacp2_scomplex_parameters::lacp2_scomplex_parameters (int matrix_layout_i,char uplo_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i, lapack_int ldb_i)
{
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	uplo = uplo_i;
	lda = lda_i;
	ldb = ldb_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n lacp2 scomplex:  m: %d, n: %d lda: %d, uplo:%c, ldb:%d\n",  m, n, lda, uplo, ldb);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{	lda = m;
		ldb = m;
		bufsize = lda*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{	lda = n;
		ldb = n;
		bufsize = lda*m;
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&b, &bref, bufsize);	
	if ((A==NULL) || (Aref==NULL) || \
		(b==NULL) || (bref==NULL)){
		EXPECT_FALSE( true) << "lacp2_float_parameters object: malloc error.";
		lacp2_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand( A, Aref, bufsize);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
lacp2_scomplex_parameters :: ~lacp2_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" lacp2_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lacp2_free();

}
/*  Test fixture class definition */
class clacp2_test  : public  ::testing::Test {
public:
   lacp2_scomplex_parameters  *clacp2_obj;
   void SetUp();
   void TearDown () { delete clacp2_obj; }
};

void clacp2_test::SetUp(){

    /* LAPACKE clacp2 prototype */
    typedef int (*Fptr_NL_LAPACKE_clacp2) (int matrix_layout, char uplo, lapack_int m,lapack_int n, const float *A, lapack_int lda,  lapack_complex_float* b, lapack_int ldb);

    Fptr_NL_LAPACKE_clacp2 clacp2;

    clacp2_obj = new lacp2_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
						eig_paramslist[idx].uplo,
						eig_paramslist[idx].m,
						eig_paramslist[idx].n,
						eig_paramslist[idx].lda,
						eig_paramslist[idx].ldb);

    idx = Circular_Increment_Index(idx);

    clacp2_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    clacp2_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(clacp2_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(clacp2_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    clacp2 = (Fptr_NL_LAPACKE_clacp2)dlsym(clacp2_obj->hModule, "LAPACKE_clacp2");
    ASSERT_TRUE(clacp2 != NULL) << "failed to get the Netlib LAPACKE_clacp2 symbol";
    

    clacp2_obj->inforef = clacp2( clacp2_obj->matrix_layout, clacp2_obj->uplo,
								clacp2_obj->m, clacp2_obj->n, (const float *)clacp2_obj->Aref,
								clacp2_obj->lda, clacp2_obj->bref, clacp2_obj->ldb);

    /* Compute libflame's Lapacke o/p  */
    clacp2_obj->info = LAPACKE_clacp2( clacp2_obj->matrix_layout, clacp2_obj->uplo,
										clacp2_obj->m, clacp2_obj->n, (const float*)clacp2_obj->A, 
										clacp2_obj->lda, clacp2_obj->b, clacp2_obj->ldb);

    if( clacp2_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_clacp2 is wrong\n", clacp2_obj->info );
    }
    if( clacp2_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_clacp2 is wrong\n", 
        clacp2_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    clacp2_obj->diff =  computeDiff_c( clacp2_obj->bufsize, 
                clacp2_obj->b, clacp2_obj->bref );

}

TEST_F(clacp2_test, clacp21) {
    EXPECT_NEAR(0.0, clacp2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(clacp2_test, clacp22) {
    EXPECT_NEAR(0.0, clacp2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(clacp2_test, clacp23) {
    EXPECT_NEAR(0.0, clacp2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(clacp2_test, clacp24) {
    EXPECT_NEAR(0.0, clacp2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class lacp2_dcomplex_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	char uplo;
	double* A;
	lapack_int lda, ldb;
	/*Output Parameter*/
	lapack_complex_double* b;
	double *Aref;
	lapack_complex_double *bref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      lacp2_dcomplex_parameters (int matrix_layout, char uplo, lapack_int m , lapack_int n, lapack_int lda, lapack_int ldb);
      ~lacp2_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
lacp2_dcomplex_parameters::lacp2_dcomplex_parameters (int matrix_layout_i,char uplo_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i, lapack_int ldb_i)
{
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	uplo = uplo_i;
	lda = lda_i;
	ldb = ldb_i;
//	#if LAPACKE_TEST_VERBOSE
	printf(" \n lacp2 dcomplex:  m: %d, n: %d lda: %d, uplo:%c, ldb:%d\n",  m, n, lda, uplo, ldb);
	//#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{	lda = m;
		ldb = m;
		bufsize = lda*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{	lda = n;
		ldb = n;
		bufsize = lda*m;
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&b, &bref, bufsize);	
	if ((A==NULL) || (Aref==NULL) || \
		(b==NULL) || (bref==NULL)){
		EXPECT_FALSE( true) << "lacp2_double_parameters object: malloc error.";
		lacp2_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand( A, Aref, bufsize);

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
lacp2_dcomplex_parameters :: ~lacp2_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" lacp2_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lacp2_free();

}
/*  Test fixture class definition */
class zlacp2_test  : public  ::testing::Test {
public:
   lacp2_dcomplex_parameters  *zlacp2_obj;
   void SetUp();
   void TearDown () { delete zlacp2_obj; }
};

void zlacp2_test::SetUp(){

    /* LAPACKE zlacp2 prototype */
    typedef int (*Fptr_NL_LAPACKE_zlacp2) (int matrix_layout, char uplo, lapack_int m,lapack_int n, const double *A, lapack_int lda,  lapack_complex_double* b, lapack_int ldb);

    Fptr_NL_LAPACKE_zlacp2 zlacp2;

    zlacp2_obj = new lacp2_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
						eig_paramslist[idx].uplo,
						eig_paramslist[idx].m,
						eig_paramslist[idx].n,
						eig_paramslist[idx].lda,
						eig_paramslist[idx].ldb);

    idx = Circular_Increment_Index(idx);

    zlacp2_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zlacp2_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zlacp2_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zlacp2_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zlacp2 = (Fptr_NL_LAPACKE_zlacp2)dlsym(zlacp2_obj->hModule, "LAPACKE_zlacp2");
    ASSERT_TRUE(zlacp2 != NULL) << "failed to get the Netlib LAPACKE_zlacp2 symbol";
    

    zlacp2_obj->inforef = zlacp2( zlacp2_obj->matrix_layout, zlacp2_obj->uplo,
								zlacp2_obj->m, zlacp2_obj->n, (const double *)zlacp2_obj->Aref,
								zlacp2_obj->lda, zlacp2_obj->bref, zlacp2_obj->ldb);

    /* Compute libflame's Lapacke o/p  */
    zlacp2_obj->info = LAPACKE_zlacp2( zlacp2_obj->matrix_layout, zlacp2_obj->uplo,
										zlacp2_obj->m, zlacp2_obj->n, (const double*)zlacp2_obj->A, 
										zlacp2_obj->lda, zlacp2_obj->b, zlacp2_obj->ldb);

    if( zlacp2_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zlacp2 is wrong\n", zlacp2_obj->info );
    }
    if( zlacp2_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zlacp2 is wrong\n", 
        zlacp2_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zlacp2_obj->diff =  computeDiff_z( zlacp2_obj->bufsize, 
                zlacp2_obj->b, zlacp2_obj->bref );

}

TEST_F(zlacp2_test, zlacp21) {
    EXPECT_NEAR(0.0, zlacp2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlacp2_test, zlacp22) {
    EXPECT_NEAR(0.0, zlacp2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlacp2_test, zlacp23) {
    EXPECT_NEAR(0.0, zlacp2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlacp2_test, zlacp24) {
    EXPECT_NEAR(0.0, zlacp2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}
