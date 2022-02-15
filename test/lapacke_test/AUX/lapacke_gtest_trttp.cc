#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"

#define trttp_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (ap!=NULL)    free(ap); \
if (apref!=NULL) free(apref);\
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
class trttp_float_parameters{

   public:
	int bufsize_a;
	int bufsize_ap;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	float* A;
	float* ap;
	lapack_int lda;
	char uplo;
	/*Output Parameter*/
	float *Aref, *apref;
	/*Return Values*/	
	lapack_int info, inforef;
   public:
      trttp_float_parameters (int matrix_layout, char uplo, lapack_int n, lapack_int lda);
      ~trttp_float_parameters ();

};

/* Constructor definition  float_common_parameters */
trttp_float_parameters:: trttp_float_parameters (int matrix_layout_i, char uplo_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n trttp float: matrix_layout = %d, uplo:%c, n: %d \n", matrix_layout, uplo, n);
	#endif	
	
	bufsize_a = lda *n;
	bufsize_ap = (n*(n+1)/2);
	
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_float_buffer_pair(&ap, &apref, bufsize_ap);
	if ((A==NULL) || (Aref==NULL) || \
		(ap == NULL) || (apref == NULL)){
		EXPECT_FALSE( true) << "trttp_float_parameters object: malloc error.";
		trttp_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand_custom_matrix( A, Aref, lda, n, uplo);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
trttp_float_parameters :: ~trttp_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" trttp_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   trttp_free();

}
/*  Test fixture class definition */
class strttp_test  : public  ::testing::Test {
public:
   trttp_float_parameters  *strttp_obj;
   void SetUp();
   void TearDown () { delete strttp_obj;}
};

void strttp_test::SetUp(){

	 /* LAPACKE strttp prototype */
    typedef int (*Fptr_NL_LAPACKE_strttp) (int matrix_layout , char uplo , lapack_int n , const float * A , lapack_int lda , float * ap);

    Fptr_NL_LAPACKE_strttp strttp;

    strttp_obj = new trttp_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].n, 
						   lin_solver_paramslist[idx].ldb);
						   

    idx = Circular_Increment_Index(idx);

    strttp_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    strttp_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(strttp_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(strttp_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*strttp library call */
    strttp = (Fptr_NL_LAPACKE_strttp)dlsym(strttp_obj->hModule, "LAPACKE_strttp");
    ASSERT_TRUE(strttp != NULL) << "failed to get the Netlib LAPACKE_strttp symbol";
    
/*Compute strttp's  o/p */
    strttp_obj->inforef = strttp( strttp_obj->matrix_layout,  strttp_obj->uplo,  strttp_obj->n,
								 (const float*)strttp_obj->Aref, strttp_obj->lda, strttp_obj->apref);

    /* Compute libflame's Lapacke o/p  */
    strttp_obj->info = LAPACKE_strttp( strttp_obj->matrix_layout,strttp_obj->uplo, strttp_obj->n, 
										(const float*)strttp_obj->A,  strttp_obj->lda, strttp_obj->ap);

    if( strttp_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_strttp is wrong\n", strttp_obj->info );
    }
    if( strttp_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_strttp is wrong\n", 
        strttp_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    strttp_obj->diff =  computeDiff_s( strttp_obj->bufsize_ap, 
                strttp_obj->ap, strttp_obj->apref );

}

TEST_F(strttp_test, strttp1) {
    EXPECT_NEAR(0.0, strttp_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(strttp_test, strttp2) {
    EXPECT_NEAR(0.0, strttp_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(strttp_test, strttp3) {
    EXPECT_NEAR(0.0, strttp_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(strttp_test, strttp4) {
    EXPECT_NEAR(0.0, strttp_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin double_common_parameters  class definition */

class trttp_double_parameters{

   public:
	int bufsize_a;
	int bufsize_ap;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	double* A;
	double* ap;
	lapack_int lda;
	char uplo;	
	/*Output Parameter*/
	double *Aref, *apref;
	/*Return Values*/	
	lapack_int info, inforef;
   public:
      trttp_double_parameters (int matrix_layout, char uplo, lapack_int n, lapack_int lda);
      ~trttp_double_parameters ();

};

/* Constructor definition  double_common_parameters */
trttp_double_parameters:: trttp_double_parameters (int matrix_layout_i, char uplo_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n trttp double: matrix_layout = %d, uplo:%c, n: %d \n", matrix_layout, uplo, n);
	#endif	
	
	bufsize_a = lda *n;
	bufsize_ap = (n*(n+1)/2);
	
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_double_buffer_pair(&ap, &apref, bufsize_ap);
	if ((A==NULL) || (Aref==NULL) || \
		(ap == NULL) || (apref == NULL)){
		EXPECT_FALSE( true) << "trttp_double_parameters object: malloc error.";
		trttp_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand_custom_matrix( A, Aref, lda, n  , uplo);
	

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
trttp_double_parameters :: ~trttp_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" trttp_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   trttp_free();

}
/*  Test fixture class definition */
class dtrttp_test  : public  ::testing::Test {
public:
   trttp_double_parameters  *dtrttp_obj;
   void SetUp();
   void TearDown () { delete dtrttp_obj;}
};

void dtrttp_test::SetUp(){

	 /* LAPACKE dtrttp prototype */
    typedef int (*Fptr_NL_LAPACKE_dtrttp) (int matrix_layout , char uplo , lapack_int n , const double * A , lapack_int lda , double * ap );

    Fptr_NL_LAPACKE_dtrttp dtrttp;

    dtrttp_obj = new trttp_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].n,
						   lin_solver_paramslist[idx].ldb);
						   

    idx = Circular_Increment_Index(idx);

    dtrttp_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dtrttp_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dtrttp_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dtrttp_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*dtrttp library call */
    dtrttp = (Fptr_NL_LAPACKE_dtrttp)dlsym(dtrttp_obj->hModule, "LAPACKE_dtrttp");
    ASSERT_TRUE(dtrttp != NULL) << "failed to get the Netlib LAPACKE_dtrttp symbol";
    
/*Compute dtrttp's  o/p */
    dtrttp_obj->inforef = dtrttp( dtrttp_obj->matrix_layout,  dtrttp_obj->uplo,  dtrttp_obj->n,
								 (const double*)dtrttp_obj->Aref, dtrttp_obj->lda, dtrttp_obj->apref);

    /* Compute libflame's Lapacke o/p  */
    dtrttp_obj->info = LAPACKE_dtrttp( dtrttp_obj->matrix_layout,dtrttp_obj->uplo, dtrttp_obj->n, 
										(const double*)dtrttp_obj->A,  dtrttp_obj->lda, dtrttp_obj->ap);

    if( dtrttp_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dtrttp is wrong\n", dtrttp_obj->info );
    }
    if( dtrttp_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dtrttp is wrong\n", 
        dtrttp_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dtrttp_obj->diff =  computeDiff_d( dtrttp_obj->bufsize_ap, 
                dtrttp_obj->ap, dtrttp_obj->apref );

}

TEST_F(dtrttp_test, dtrttp1) {
    EXPECT_NEAR(0.0, dtrttp_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dtrttp_test, dtrttp2) {
    EXPECT_NEAR(0.0, dtrttp_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dtrttp_test, dtrttp3) {
    EXPECT_NEAR(0.0, dtrttp_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dtrttp_test, dtrttp4) {
    EXPECT_NEAR(0.0, dtrttp_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin scomplex_common_parameters  class definition */

class trttp_scomplex_parameters{

   public:
	int bufsize_a;
	int bufsize_ap;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_complex_float* A;
	lapack_complex_float* ap;
	lapack_int lda;
	char uplo;	
	/*Output Parameter*/
	lapack_complex_float *Aref, *apref;
	/*Return Values*/	
	lapack_int info, inforef;
   public:
      trttp_scomplex_parameters (int matrix_layout, char uplo, lapack_int n, lapack_int lda);
      ~trttp_scomplex_parameters ();

};

/* Constructor definition  scomplex_common_parameters */
trttp_scomplex_parameters:: trttp_scomplex_parameters (int matrix_layout_i, char uplo_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n trttp scomplex: matrix_layout = %d, uplo:%c, n: %d \n", matrix_layout, uplo, n);
	#endif	
	
	bufsize_a = lda *n;
	bufsize_ap = (n*(n+1)/2);
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&ap, &apref, bufsize_ap);
	if ((A==NULL) || (Aref==NULL) || \
		(ap == NULL) || (apref == NULL)){
		EXPECT_FALSE( true) << "trttp_scomplex_parameters object: malloc error.";
		trttp_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix( A, Aref, lda, n, uplo);


} /* end of Constructor  */

/* Destructor definition  'scomplex_common_parameters' class  */
trttp_scomplex_parameters :: ~trttp_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" trttp_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   trttp_free();

}
/*  Test fixture class definition */
class ctrttp_test  : public  ::testing::Test {
public:
   trttp_scomplex_parameters  *ctrttp_obj;
   void SetUp();
   void TearDown () { delete ctrttp_obj;}
};

void ctrttp_test::SetUp(){

	 /* LAPACKE ctrttp prototype */
    typedef int (*Fptr_NL_LAPACKE_ctrttp) (int matrix_layout , char uplo , lapack_int n , const lapack_complex_float * A , lapack_int lda , lapack_complex_float * ap);

    Fptr_NL_LAPACKE_ctrttp ctrttp;

    ctrttp_obj = new trttp_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].n,
						   lin_solver_paramslist[idx].ldb);
						   

    idx = Circular_Increment_Index(idx);

    ctrttp_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ctrttp_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ctrttp_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ctrttp_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*ctrttp library call */
    ctrttp = (Fptr_NL_LAPACKE_ctrttp)dlsym(ctrttp_obj->hModule, "LAPACKE_ctrttp");
    ASSERT_TRUE(ctrttp != NULL) << "failed to get the Netlib LAPACKE_ctrttp symbol";
    
/*Compute ctrttp's  o/p */
    ctrttp_obj->inforef = ctrttp( ctrttp_obj->matrix_layout,  ctrttp_obj->uplo,  ctrttp_obj->n,
								 (const lapack_complex_float*)ctrttp_obj->Aref, ctrttp_obj->lda, ctrttp_obj->apref);

    /* Compute libflame's Lapacke o/p  */
    ctrttp_obj->info = LAPACKE_ctrttp( ctrttp_obj->matrix_layout,ctrttp_obj->uplo, ctrttp_obj->n, 
										 (const lapack_complex_float*)ctrttp_obj->A,  ctrttp_obj->lda, ctrttp_obj->ap);

    if( ctrttp_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ctrttp is wrong\n", ctrttp_obj->info );
    }
    if( ctrttp_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ctrttp is wrong\n", 
        ctrttp_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    ctrttp_obj->diff =  computeDiff_c( ctrttp_obj->bufsize_ap, 
                ctrttp_obj->ap, ctrttp_obj->apref );

}

TEST_F(ctrttp_test, ctrttp1) {
    EXPECT_NEAR(0.0, ctrttp_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ctrttp_test, ctrttp2) {
    EXPECT_NEAR(0.0, ctrttp_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ctrttp_test, ctrttp3) {
    EXPECT_NEAR(0.0, ctrttp_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ctrttp_test, ctrttp4) {
    EXPECT_NEAR(0.0, ctrttp_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */

class trttp_dcomplex_parameters{

   public:
	int bufsize_a;
	int bufsize_ap;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_complex_double* A;
	lapack_complex_double* ap;
	lapack_int lda;
	char uplo;	
	/*Output Parameter*/
	lapack_complex_double *Aref, *apref;
	/*Return Values*/	
	lapack_int info, inforef;
   public:
      trttp_dcomplex_parameters (int matrix_layout, char uplo, lapack_int n, lapack_int lda);
      ~trttp_dcomplex_parameters ();

};

/* Constructor definition  dcomplex_common_parameters */
trttp_dcomplex_parameters:: trttp_dcomplex_parameters (int matrix_layout_i, char uplo_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n trttp dcomplex: matrix_layout = %d, uplo:%c, n: %d \n", matrix_layout, uplo, n);
	#endif	
	
	bufsize_a = lda *n;
	bufsize_ap = (n*(n+1)/2);
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&ap, &apref, bufsize_ap);
	if ((A==NULL) || (Aref==NULL) || \
		(ap == NULL) || (apref == NULL)){
		EXPECT_FALSE( true) << "trttp_dcomplex_parameters object: malloc error.";
		trttp_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix( A, Aref, lda, n, uplo);

} /* end of Constructor  */

/* Destructor definition  'dcomplex_common_parameters' class  */
trttp_dcomplex_parameters :: ~trttp_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" trttp_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   trttp_free();

}
/*  Test fixture class definition */
class ztrttp_test  : public  ::testing::Test {
public:
   trttp_dcomplex_parameters  *ztrttp_obj;
   void SetUp();
   void TearDown () { delete ztrttp_obj;}
};

void ztrttp_test::SetUp(){

	 /* LAPACKE ztrttp prototype */
    typedef int (*Fptr_NL_LAPACKE_ztrttp) (int matrix_layout , char uplo , lapack_int n , const lapack_complex_double * A , lapack_int lda , lapack_complex_double * ap );

    Fptr_NL_LAPACKE_ztrttp ztrttp;

    ztrttp_obj = new trttp_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].n, 
						   lin_solver_paramslist[idx].ldb);
						   

    idx = Circular_Increment_Index(idx);

    ztrttp_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ztrttp_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ztrttp_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ztrttp_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*ztrttp library call */
    ztrttp = (Fptr_NL_LAPACKE_ztrttp)dlsym(ztrttp_obj->hModule, "LAPACKE_ztrttp");
    ASSERT_TRUE(ztrttp != NULL) << "failed to get the Netlib LAPACKE_ztrttp symbol";
    
/*Compute ztrttp's  o/p */
    ztrttp_obj->inforef = ztrttp( ztrttp_obj->matrix_layout,  ztrttp_obj->uplo,  ztrttp_obj->n,
								 (const lapack_complex_double*)ztrttp_obj->Aref, ztrttp_obj->lda, ztrttp_obj->apref);

    /* Compute libflame's Lapacke o/p  */
    ztrttp_obj->info = LAPACKE_ztrttp( ztrttp_obj->matrix_layout,ztrttp_obj->uplo, ztrttp_obj->n, 
										(const lapack_complex_double*)ztrttp_obj->A,  ztrttp_obj->lda, ztrttp_obj->ap);

    if( ztrttp_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ztrttp is wrong\n", ztrttp_obj->info );
    }
    if( ztrttp_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ztrttp is wrong\n", 
        ztrttp_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    ztrttp_obj->diff =  computeDiff_z( ztrttp_obj->bufsize_ap, 
                ztrttp_obj->ap, ztrttp_obj->apref );

}

TEST_F(ztrttp_test, ztrttp1) {
    EXPECT_NEAR(0.0, ztrttp_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ztrttp_test, ztrttp2) {
    EXPECT_NEAR(0.0, ztrttp_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ztrttp_test, ztrttp3) {
    EXPECT_NEAR(0.0, ztrttp_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ztrttp_test, ztrttp4) {
    EXPECT_NEAR(0.0, ztrttp_obj->diff, LAPACKE_GTEST_THRESHOLD);
}