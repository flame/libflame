#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"

#define trttf_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (arf != NULL) free(arf); \
if (arfref != NULL) free(arfref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin float_common_parameters  class definition */
class trttf_float_parameters{

   public:
	int bufsize_a;
	int bufsize_arf;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	float* A;
	lapack_int lda;
	char uplo;
	char transr;
	/*Output Parameter*/
	float *arf;
	float *Aref, *arfref;
	/*Return Values*/	
	lapack_int info, inforef;
   public:
      trttf_float_parameters (int matrix_layout,  char transr, char uplo,  lapack_int n, lapack_int lda);
      ~trttf_float_parameters ();

};

/* Constructor definition  float_common_parameters */
trttf_float_parameters:: trttf_float_parameters (int matrix_layout_i, char transr_i, char uplo_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;	
	transr = transr_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n trttf float: matrix_layout = %d, uplo:%c, n: %d, lda:%d, transr:%c \n", matrix_layout, uplo, n, lda, transr);
	#endif	

	bufsize_arf = (n*(n +1)/2);
	bufsize_a = lda*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&A, &Aref, bufsize_a);	
	lapacke_gtest_alloc_float_buffer_pair(&arf, &arfref, bufsize_arf);
	if ((A==NULL) || (Aref==NULL) || \
		(arf == NULL) || (arfref == NULL)){
		EXPECT_FALSE( true) << "trttf_float_parameters object: malloc error.";
		trttf_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand_custom_matrix(A, Aref, lda, n, uplo);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
trttf_float_parameters :: ~trttf_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" trttf_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   trttf_free();

}
/*  Test fixture class definition */
class strttf_test  : public  ::testing::Test {
public:
   trttf_float_parameters  *strttf_obj;
   void SetUp();
   void TearDown () { delete strttf_obj;}
};

void strttf_test::SetUp(){

	 /* LAPACKE strttf prototype */
    typedef int (*Fptr_NL_LAPACKE_strttf) (int matrix_layout , char transr , char uplo , lapack_int n , const float * A , lapack_int lda , float * arf );

    Fptr_NL_LAPACKE_strttf strttf;

    strttf_obj = new trttf_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].n,						   
						   lin_solver_paramslist[idx].n);
						   

    idx = Circular_Increment_Index(idx);

    strttf_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    strttf_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(strttf_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(strttf_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*strttf library call */
    strttf = (Fptr_NL_LAPACKE_strttf)dlsym(strttf_obj->hModule, "LAPACKE_strttf");
    ASSERT_TRUE(strttf != NULL) << "failed to get the Netlib LAPACKE_strttf symbol";
    
/*Compute strttf's  o/p */
    strttf_obj->inforef = strttf( strttf_obj->matrix_layout, strttf_obj->transr, strttf_obj->uplo,  strttf_obj->n,
								(const float *)strttf_obj->Aref, strttf_obj->lda, strttf_obj->arfref);

    /* Compute libflame's Lapacke o/p  */
    strttf_obj->info = LAPACKE_strttf( strttf_obj->matrix_layout, strttf_obj->transr, strttf_obj->uplo, strttf_obj->n, 
										(const float *)strttf_obj->A,  strttf_obj->lda, strttf_obj->arf);

    if( strttf_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_strttf is wrong\n", strttf_obj->info );
    }
    if( strttf_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_strttf is wrong\n", 
        strttf_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    strttf_obj->diff =  computeDiff_s( strttf_obj->bufsize_arf, 
                strttf_obj->arf, strttf_obj->arfref );

}

TEST_F(strttf_test, strttf1) {
    EXPECT_NEAR(0.0, strttf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(strttf_test, strttf2) {
    EXPECT_NEAR(0.0, strttf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(strttf_test, strttf3) {
    EXPECT_NEAR(0.0, strttf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(strttf_test, strttf4) {
    EXPECT_NEAR(0.0, strttf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin double_common_parameters  class definition */
class trttf_double_parameters{

   public:
	int bufsize_a;
	int bufsize_arf;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	double* A;
	lapack_int lda;
	char uplo;
	char transr;
	/*Output Parameter*/
	double *arf;
	double *Aref, *arfref;
	/*Return Values*/	
	lapack_int info, inforef;
   public:
      trttf_double_parameters (int matrix_layout,  char transr, char uplo,  lapack_int n, lapack_int lda);
      ~trttf_double_parameters ();

};

/* Constructor definition  float_common_parameters */
trttf_double_parameters:: trttf_double_parameters (int matrix_layout_i, char transr_i, char uplo_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;	
	transr = transr_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n trttf float: matrix_layout = %d, uplo:%c, n: %d, lda:%d, transr:%c \n", matrix_layout, uplo, n, lda, transr);
	#endif	

	bufsize_arf = (n*(n +1)/2);
	bufsize_a = lda*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&A, &Aref, bufsize_a);	
	lapacke_gtest_alloc_double_buffer_pair(&arf, &arfref, bufsize_arf);
	if ((A==NULL) || (Aref==NULL) || \
		(arf == NULL) || (arfref == NULL)){
		EXPECT_FALSE( true) << "trttf_float_parameters object: malloc error.";
		trttf_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand_custom_matrix(A, Aref, lda, n, uplo);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
trttf_double_parameters :: ~trttf_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" trttf_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   trttf_free();

}
/*  Test fixture class definition */
class dtrttf_test  : public  ::testing::Test {
public:
   trttf_double_parameters  *dtrttf_obj;
   void SetUp();
   void TearDown () { delete dtrttf_obj;}
};

void dtrttf_test::SetUp(){

	 /* LAPACKE dtrttf prototype */
    typedef int (*Fptr_NL_LAPACKE_dtrttf) (int matrix_layout, char transr, char uplo,  lapack_int n,  
											const double *A, lapack_int lda, double* arf);

    Fptr_NL_LAPACKE_dtrttf dtrttf;

    dtrttf_obj = new trttf_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].n,						   
						   lin_solver_paramslist[idx].n);
						   

    idx = Circular_Increment_Index(idx);

    dtrttf_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dtrttf_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dtrttf_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dtrttf_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*dtrttf library call */
    dtrttf = (Fptr_NL_LAPACKE_dtrttf)dlsym(dtrttf_obj->hModule, "LAPACKE_dtrttf");
    ASSERT_TRUE(dtrttf != NULL) << "failed to get the Netlib LAPACKE_dtrttf symbol";
    
/*Compute dtrttf's  o/p */
    dtrttf_obj->inforef = dtrttf( dtrttf_obj->matrix_layout, dtrttf_obj->transr, dtrttf_obj->uplo,  dtrttf_obj->n,
								(const double *)dtrttf_obj->Aref, dtrttf_obj->lda, dtrttf_obj->arfref);

    /* Compute libflame's Lapacke o/p  */
    dtrttf_obj->info = LAPACKE_dtrttf( dtrttf_obj->matrix_layout, dtrttf_obj->transr, dtrttf_obj->uplo, dtrttf_obj->n, 
										(const double *)dtrttf_obj->A,  dtrttf_obj->lda, dtrttf_obj->arf);

    if( dtrttf_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dtrttf is wrong\n", dtrttf_obj->info );
    }
    if( dtrttf_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dtrttf is wrong\n", 
        dtrttf_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dtrttf_obj->diff =  computeDiff_d( dtrttf_obj->bufsize_arf, 
                dtrttf_obj->arf, dtrttf_obj->arfref );

}

TEST_F(dtrttf_test, dtrttf1) {
    EXPECT_NEAR(0.0, dtrttf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dtrttf_test, dtrttf2) {
    EXPECT_NEAR(0.0, dtrttf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dtrttf_test, dtrttf3) {
    EXPECT_NEAR(0.0, dtrttf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dtrttf_test, dtrttf4) {
    EXPECT_NEAR(0.0, dtrttf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}
/* Begin scomplex_common_parameters  class definition */
class trttf_scomplex_parameters{

   public:
	int bufsize_a;
	int bufsize_arf;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_complex_float* A;
	lapack_int lda;
	char uplo;
	char transr;
	/*Output Parameter*/
	lapack_complex_float *arf;
	lapack_complex_float *Aref, *arfref;
	/*Return Values*/	
	lapack_int info, inforef;
   public:
      trttf_scomplex_parameters (int matrix_layout,  char transr, char uplo,  lapack_int n, lapack_int lda);
      ~trttf_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
trttf_scomplex_parameters:: trttf_scomplex_parameters (int matrix_layout_i, char transr_i, char uplo_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;	
	transr = transr_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n trttf float: matrix_layout = %d, uplo:%c, n: %d, lda:%d, transr:%c \n", matrix_layout, uplo, n, lda, transr);
	#endif	

	bufsize_arf = (n*(n +1)/2);
	bufsize_a = lda*n;
	
	if (transr == 'T')
		transr = 'C';
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize_a);	
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&arf, &arfref, bufsize_arf);
	if ((A==NULL) || (Aref==NULL) || \
		(arf == NULL) || (arfref == NULL)){
		EXPECT_FALSE( true) << "trttf_float_parameters object: malloc error.";
		trttf_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix(A, Aref, lda, n, uplo);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
trttf_scomplex_parameters :: ~trttf_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" trttf_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   trttf_free();

}
/*  Test fixture class definition */
class ctrttf_test  : public  ::testing::Test {
public:
   trttf_scomplex_parameters  *ctrttf_obj;
   void SetUp();
   void TearDown () { delete ctrttf_obj;}
};

void ctrttf_test::SetUp(){

	 /* LAPACKE ctrttf prototype */
    typedef int (*Fptr_NL_LAPACKE_ctrttf) (int matrix_layout, char transr, char uplo,  lapack_int n,  
											const lapack_complex_float *A, lapack_int lda,
											 lapack_complex_float* arf);

    Fptr_NL_LAPACKE_ctrttf ctrttf;

    ctrttf_obj = new trttf_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].n,						   
						   lin_solver_paramslist[idx].n);
						   

    idx = Circular_Increment_Index(idx);

    ctrttf_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ctrttf_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ctrttf_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ctrttf_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*ctrttf library call */
    ctrttf = (Fptr_NL_LAPACKE_ctrttf)dlsym(ctrttf_obj->hModule, "LAPACKE_ctrttf");
    ASSERT_TRUE(ctrttf != NULL) << "failed to get the Netlib LAPACKE_ctrttf symbol";
    
/*Compute ctrttf's  o/p */
    ctrttf_obj->inforef = ctrttf( ctrttf_obj->matrix_layout, ctrttf_obj->transr, ctrttf_obj->uplo,  ctrttf_obj->n,
								(const lapack_complex_float *)ctrttf_obj->Aref, ctrttf_obj->lda, ctrttf_obj->arfref);

    /* Compute libflame's Lapacke o/p  */
    ctrttf_obj->info = LAPACKE_ctrttf( ctrttf_obj->matrix_layout, ctrttf_obj->transr, ctrttf_obj->uplo, ctrttf_obj->n, 
										(const lapack_complex_float *)ctrttf_obj->A,  ctrttf_obj->lda, ctrttf_obj->arf);

    if( ctrttf_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ctrttf is wrong\n", ctrttf_obj->info );
    }
    if( ctrttf_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ctrttf is wrong\n", 
        ctrttf_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    ctrttf_obj->diff =  computeDiff_c( ctrttf_obj->bufsize_arf, 
                ctrttf_obj->arf, ctrttf_obj->arfref );

}

TEST_F(ctrttf_test, ctrttf1) {
    EXPECT_NEAR(0.0, ctrttf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ctrttf_test, ctrttf2) {
    EXPECT_NEAR(0.0, ctrttf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ctrttf_test, ctrttf3) {
    EXPECT_NEAR(0.0, ctrttf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ctrttf_test, ctrttf4) {
    EXPECT_NEAR(0.0, ctrttf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class trttf_dcomplex_parameters{

   public:
	int bufsize_a;
	int bufsize_arf;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_complex_double* A;
	lapack_int lda;
	char uplo;
	char transr;
	/*Output Parameter*/
	lapack_complex_double *arf;
	lapack_complex_double *Aref, *arfref;
	/*Return Values*/	
	lapack_int info, inforef;
   public:
      trttf_dcomplex_parameters (int matrix_layout,  char transr, char uplo,  lapack_int n, lapack_int lda);
      ~trttf_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
trttf_dcomplex_parameters:: trttf_dcomplex_parameters (int matrix_layout_i, char transr_i, char uplo_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;	
	transr = transr_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n trttf float: matrix_layout = %d, uplo:%c, n: %d, lda:%d, transr:%c \n", matrix_layout, uplo, n, lda, transr);
	#endif	

	bufsize_arf = (n*(n +1)/2);
	bufsize_a = lda*n;
	
	if (transr == 'T')
		transr = 'C';
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize_a);	
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&arf, &arfref, bufsize_arf);
	if ((A==NULL) || (Aref==NULL) || \
		(arf == NULL) || (arfref == NULL)){
		EXPECT_FALSE( true) << "trttf_double_parameters object: malloc error.";
		trttf_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix(A, Aref, lda, n, uplo);

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
trttf_dcomplex_parameters :: ~trttf_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" trttf_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   trttf_free();

}
/*  Test fixture class definition */
class ztrttf_test  : public  ::testing::Test {
public:
   trttf_dcomplex_parameters  *ztrttf_obj;
   void SetUp();
   void TearDown () { delete ztrttf_obj;}
};

void ztrttf_test::SetUp(){

	 /* LAPACKE ztrttf prototype */
    typedef int (*Fptr_NL_LAPACKE_ztrttf) (int matrix_layout, char transr, char uplo,  lapack_int n,  
											const lapack_complex_double *A, lapack_int lda,
											 lapack_complex_double* arf);

    Fptr_NL_LAPACKE_ztrttf ztrttf;

    ztrttf_obj = new trttf_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].n,						   
						   lin_solver_paramslist[idx].n);
						   

    idx = Circular_Increment_Index(idx);

    ztrttf_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ztrttf_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ztrttf_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ztrttf_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*ztrttf library call */
    ztrttf = (Fptr_NL_LAPACKE_ztrttf)dlsym(ztrttf_obj->hModule, "LAPACKE_ztrttf");
    ASSERT_TRUE(ztrttf != NULL) << "failed to get the Netlib LAPACKE_ztrttf symbol";
    
/*Compute ztrttf's  o/p */
    ztrttf_obj->inforef = ztrttf( ztrttf_obj->matrix_layout, ztrttf_obj->transr, ztrttf_obj->uplo,  ztrttf_obj->n,
								(const lapack_complex_double *)ztrttf_obj->Aref, ztrttf_obj->lda, 
								ztrttf_obj->arfref);

    /* Compute libflame's Lapacke o/p  */
    ztrttf_obj->info = LAPACKE_ztrttf( ztrttf_obj->matrix_layout, ztrttf_obj->transr, ztrttf_obj->uplo, ztrttf_obj->n, 
										(const lapack_complex_double *)ztrttf_obj->A,  ztrttf_obj->lda, ztrttf_obj->arf);

    if( ztrttf_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ztrttf is wrong\n", ztrttf_obj->info );
    }
    if( ztrttf_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ztrttf is wrong\n", 
        ztrttf_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    ztrttf_obj->diff =  computeDiff_z( ztrttf_obj->bufsize_arf, 
                ztrttf_obj->arf, ztrttf_obj->arfref );

}

TEST_F(ztrttf_test, ztrttf1) {
    EXPECT_NEAR(0.0, ztrttf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ztrttf_test, ztrttf2) {
    EXPECT_NEAR(0.0, ztrttf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ztrttf_test, ztrttf3) {
    EXPECT_NEAR(0.0, ztrttf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ztrttf_test, ztrttf4) {
    EXPECT_NEAR(0.0, ztrttf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}