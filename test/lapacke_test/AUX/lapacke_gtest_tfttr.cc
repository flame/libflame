#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"

#define tfttr_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (ap!=NULL)    free(ap); \
if (apref!=NULL) free(apref);\
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
class tfttr_float_parameters{

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
	char uplo, transr;
	/*Output Parameter*/
	float *Aref, *apref;
	/*Return Values*/	
	lapack_int info, inforef;
   public:
      tfttr_float_parameters (int matrix_layout, char transr, char uplo, lapack_int n, lapack_int lda);
      ~tfttr_float_parameters ();

};

/* Constructor definition  float_common_parameters */
tfttr_float_parameters:: tfttr_float_parameters (int matrix_layout_i, char transr_i, char uplo_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	lda = lda_i;
	transr = transr_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n tfttr float: matrix_layout = %d, uplo:%c, n: %d \n", matrix_layout, uplo, n);
	#endif	
	
	bufsize_a = lda *n;
	bufsize_ap = (n*(n+1)/2);
	
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_float_buffer_pair(&ap, &apref, bufsize_ap);
	if ((A==NULL) || (Aref==NULL) || \
		(ap == NULL) || (apref == NULL)){
		EXPECT_FALSE( true) << "tfttr_float_parameters object: malloc error.";
		tfttr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand_custom_matrix( A, Aref, lda, n, uplo);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
tfttr_float_parameters :: ~tfttr_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" tfttr_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   tfttr_free();

}
/*  Test fixture class definition */
class stfttr_test  : public  ::testing::Test {
public:
   tfttr_float_parameters  *stfttr_obj;
   void SetUp();
   void TearDown () { delete stfttr_obj;}
};

void stfttr_test::SetUp(){

	 /* LAPACKE stfttr prototype */
    typedef int (*Fptr_NL_LAPACKE_stfttr) (int matrix_layout , char transr , char uplo , lapack_int n , const float * arf , float * a , lapack_int lda );

    Fptr_NL_LAPACKE_stfttr stfttr;

    stfttr_obj = new tfttr_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
							lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].n, 
						   lin_solver_paramslist[idx].ldb);
						   

    idx = Circular_Increment_Index(idx);

    stfttr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    stfttr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(stfttr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(stfttr_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*stfttr library call */
    stfttr = (Fptr_NL_LAPACKE_stfttr)dlsym(stfttr_obj->hModule, "LAPACKE_stfttr");
    ASSERT_TRUE(stfttr != NULL) << "failed to get the Netlib LAPACKE_stfttr symbol";
    
/*Compute stfttr's  o/p */
    stfttr_obj->inforef = stfttr( stfttr_obj->matrix_layout,  stfttr_obj->transr, stfttr_obj->uplo,  stfttr_obj->n,
								 (const float*)stfttr_obj->apref, stfttr_obj->Aref, stfttr_obj->lda);

    /* Compute libflame's Lapacke o/p  */
    stfttr_obj->info = LAPACKE_stfttr( stfttr_obj->matrix_layout, stfttr_obj->transr, stfttr_obj->uplo, stfttr_obj->n, 
										(const float*)stfttr_obj->ap,  stfttr_obj->A, stfttr_obj->lda);

    if( stfttr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_stfttr is wrong\n", stfttr_obj->info );
    }
    if( stfttr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_stfttr is wrong\n", 
        stfttr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    stfttr_obj->diff =  computeDiff_s( stfttr_obj->bufsize_a, 
                stfttr_obj->A, stfttr_obj->Aref );

}

TEST_F(stfttr_test, stfttr1) {
    EXPECT_NEAR(0.0, stfttr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(stfttr_test, stfttr2) {
    EXPECT_NEAR(0.0, stfttr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(stfttr_test, stfttr3) {
    EXPECT_NEAR(0.0, stfttr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(stfttr_test, stfttr4) {
    EXPECT_NEAR(0.0, stfttr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin double_common_parameters  class definition */

class tfttr_double_parameters{

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
	char uplo, transr;	
	/*Output Parameter*/
	double *Aref, *apref;
	/*Return Values*/	
	lapack_int info, inforef;
   public:
      tfttr_double_parameters (int matrix_layout, char transr, char uplo, lapack_int n, lapack_int lda);
      ~tfttr_double_parameters ();

};

/* Constructor definition  double_common_parameters */
tfttr_double_parameters:: tfttr_double_parameters (int matrix_layout_i, char transr_i, char uplo_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	lda = lda_i;
	transr = transr_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n tfttr double: matrix_layout = %d, uplo:%c, n: %d \n", matrix_layout, uplo, n);
	#endif	
	
	bufsize_a = lda *n;
	bufsize_ap = (n*(n+1)/2);
	
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_double_buffer_pair(&ap, &apref, bufsize_ap);
	if ((A==NULL) || (Aref==NULL) || \
		(ap == NULL) || (apref == NULL)){
		EXPECT_FALSE( true) << "tfttr_double_parameters object: malloc error.";
		tfttr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand_custom_matrix( A, Aref, lda, n  , uplo);
	

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
tfttr_double_parameters :: ~tfttr_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" tfttr_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   tfttr_free();

}
/*  Test fixture class definition */
class dtfttr_test  : public  ::testing::Test {
public:
   tfttr_double_parameters  *dtfttr_obj;
   void SetUp();
   void TearDown () { delete dtfttr_obj;}
};

void dtfttr_test::SetUp(){

	 /* LAPACKE dtfttr prototype */
    typedef int (*Fptr_NL_LAPACKE_dtfttr) (int matrix_layout , char transr , char uplo , lapack_int n , const double * arf , double * a , lapack_int lda );

    Fptr_NL_LAPACKE_dtfttr dtfttr;

    dtfttr_obj = new tfttr_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
							lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].n,
						   lin_solver_paramslist[idx].ldb);
						   

    idx = Circular_Increment_Index(idx);

    dtfttr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dtfttr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dtfttr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dtfttr_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*dtfttr library call */
    dtfttr = (Fptr_NL_LAPACKE_dtfttr)dlsym(dtfttr_obj->hModule, "LAPACKE_dtfttr");
    ASSERT_TRUE(dtfttr != NULL) << "failed to get the Netlib LAPACKE_dtfttr symbol";
    
/*Compute dtfttr's  o/p */
    dtfttr_obj->inforef = dtfttr( dtfttr_obj->matrix_layout, dtfttr_obj->transr,  dtfttr_obj->uplo,  dtfttr_obj->n,
								 (const double*)dtfttr_obj->apref, dtfttr_obj->Aref, dtfttr_obj->lda);

    /* Compute libflame's Lapacke o/p  */
    dtfttr_obj->info = LAPACKE_dtfttr( dtfttr_obj->matrix_layout, dtfttr_obj->transr, dtfttr_obj->uplo, dtfttr_obj->n, 
										(const double*)dtfttr_obj->ap,  dtfttr_obj->A, dtfttr_obj->lda);

    if( dtfttr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dtfttr is wrong\n", dtfttr_obj->info );
    }
    if( dtfttr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dtfttr is wrong\n", 
        dtfttr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dtfttr_obj->diff =  computeDiff_d( dtfttr_obj->bufsize_a, 
                dtfttr_obj->A, dtfttr_obj->Aref );

}

TEST_F(dtfttr_test, dtfttr1) {
    EXPECT_NEAR(0.0, dtfttr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dtfttr_test, dtfttr2) {
    EXPECT_NEAR(0.0, dtfttr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dtfttr_test, dtfttr3) {
    EXPECT_NEAR(0.0, dtfttr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dtfttr_test, dtfttr4) {
    EXPECT_NEAR(0.0, dtfttr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin scomplex_common_parameters  class definition */

class tfttr_scomplex_parameters{

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
	char uplo, transr;	
	/*Output Parameter*/
	lapack_complex_float *Aref, *apref;
	/*Return Values*/	
	lapack_int info, inforef;
   public:
      tfttr_scomplex_parameters (int matrix_layout, char transr, char uplo, lapack_int n, lapack_int lda);
      ~tfttr_scomplex_parameters ();

};

/* Constructor definition  scomplex_common_parameters */
tfttr_scomplex_parameters:: tfttr_scomplex_parameters (int matrix_layout_i, char transr_i, char uplo_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	lda = lda_i;
	transr = transr_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n tfttr scomplex: matrix_layout = %d, uplo:%c, n: %d \n", matrix_layout, uplo, n);
	#endif	
	
	if (transr == 'T')
		transr = 'C';
	
	bufsize_a = lda *n;
	bufsize_ap = (n*(n+1)/2);
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&ap, &apref, bufsize_ap);
	if ((A==NULL) || (Aref==NULL) || \
		(ap == NULL) || (apref == NULL)){
		EXPECT_FALSE( true) << "tfttr_scomplex_parameters object: malloc error.";
		tfttr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix( A, Aref, lda, n, uplo);


} /* end of Constructor  */

/* Destructor definition  'scomplex_common_parameters' class  */
tfttr_scomplex_parameters :: ~tfttr_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" tfttr_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   tfttr_free();

}
/*  Test fixture class definition */
class ctfttr_test  : public  ::testing::Test {
public:
   tfttr_scomplex_parameters  *ctfttr_obj;
   void SetUp();
   void TearDown () { delete ctfttr_obj;}
};

void ctfttr_test::SetUp(){

	 /* LAPACKE ctfttr prototype */
    typedef int (*Fptr_NL_LAPACKE_ctfttr) (int matrix_layout , char transr , char uplo , lapack_int n , const lapack_complex_float * arf , lapack_complex_float * a , lapack_int lda );

    Fptr_NL_LAPACKE_ctfttr ctfttr;

    ctfttr_obj = new tfttr_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].transr,
						   lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].n,
						   lin_solver_paramslist[idx].ldb);
						   

    idx = Circular_Increment_Index(idx);

    ctfttr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ctfttr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ctfttr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ctfttr_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*ctfttr library call */
    ctfttr = (Fptr_NL_LAPACKE_ctfttr)dlsym(ctfttr_obj->hModule, "LAPACKE_ctfttr");
    ASSERT_TRUE(ctfttr != NULL) << "failed to get the Netlib LAPACKE_ctfttr symbol";
    
/*Compute ctfttr's  o/p */
    ctfttr_obj->inforef = ctfttr( ctfttr_obj->matrix_layout, ctfttr_obj->transr,  ctfttr_obj->uplo,  ctfttr_obj->n,
								 (const lapack_complex_float*)ctfttr_obj->ap, ctfttr_obj->Aref, ctfttr_obj->lda);

    /* Compute libflame's Lapacke o/p  */
    ctfttr_obj->info = LAPACKE_ctfttr( ctfttr_obj->matrix_layout,ctfttr_obj->transr, ctfttr_obj->uplo, ctfttr_obj->n, 
										 (const lapack_complex_float*)ctfttr_obj->ap,  ctfttr_obj->A, ctfttr_obj->lda);

    if( ctfttr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ctfttr is wrong\n", ctfttr_obj->info );
    }
    if( ctfttr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ctfttr is wrong\n", 
        ctfttr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    ctfttr_obj->diff =  computeDiff_c( ctfttr_obj->bufsize_a, 
                ctfttr_obj->A, ctfttr_obj->Aref );

}

TEST_F(ctfttr_test, ctfttr1) {
    EXPECT_NEAR(0.0, ctfttr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ctfttr_test, ctfttr2) {
    EXPECT_NEAR(0.0, ctfttr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ctfttr_test, ctfttr3) {
    EXPECT_NEAR(0.0, ctfttr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ctfttr_test, ctfttr4) {
    EXPECT_NEAR(0.0, ctfttr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */

class tfttr_dcomplex_parameters{

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
	char uplo, transr;	
	/*Output Parameter*/
	lapack_complex_double *Aref, *apref;
	/*Return Values*/	
	lapack_int info, inforef;
   public:
      tfttr_dcomplex_parameters (int matrix_layout, char transr, char uplo, lapack_int n, lapack_int lda);
      ~tfttr_dcomplex_parameters ();

};

/* Constructor definition  dcomplex_common_parameters */
tfttr_dcomplex_parameters:: tfttr_dcomplex_parameters (int matrix_layout_i, char transr_i, char uplo_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	lda = lda_i;
	transr = transr_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n tfttr dcomplex: matrix_layout = %d, uplo:%c, n: %d \n", matrix_layout, uplo, n);
	#endif	
	
	if (transr == 'T')
		transr = 'C';
	
	bufsize_a = lda *n;
	bufsize_ap = (n*(n+1)/2);
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&ap, &apref, bufsize_ap);
	if ((A==NULL) || (Aref==NULL) || \
		(ap == NULL) || (apref == NULL)){
		EXPECT_FALSE( true) << "tfttr_dcomplex_parameters object: malloc error.";
		tfttr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix( A, Aref, lda, n, uplo);

} /* end of Constructor  */

/* Destructor definition  'dcomplex_common_parameters' class  */
tfttr_dcomplex_parameters :: ~tfttr_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" tfttr_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   tfttr_free();

}
/*  Test fixture class definition */
class ztfttr_test  : public  ::testing::Test {
public:
   tfttr_dcomplex_parameters  *ztfttr_obj;
   void SetUp();
   void TearDown () { delete ztfttr_obj;}
};

void ztfttr_test::SetUp(){

	 /* LAPACKE ztfttr prototype */
    typedef int (*Fptr_NL_LAPACKE_ztfttr) (int matrix_layout , char transr , char uplo , lapack_int n , const lapack_complex_double * arf , lapack_complex_double * a , lapack_int lda);

    Fptr_NL_LAPACKE_ztfttr ztfttr;

    ztfttr_obj = new tfttr_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
							lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].n, 
						   lin_solver_paramslist[idx].ldb);
						   

    idx = Circular_Increment_Index(idx);

    ztfttr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ztfttr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ztfttr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ztfttr_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*ztfttr library call */
    ztfttr = (Fptr_NL_LAPACKE_ztfttr)dlsym(ztfttr_obj->hModule, "LAPACKE_ztfttr");
    ASSERT_TRUE(ztfttr != NULL) << "failed to get the Netlib LAPACKE_ztfttr symbol";
    
/*Compute ztfttr's  o/p */
    ztfttr_obj->inforef = ztfttr( ztfttr_obj->matrix_layout,  ztfttr_obj->transr, ztfttr_obj->uplo,  ztfttr_obj->n,
								 (const lapack_complex_double*)ztfttr_obj->apref, ztfttr_obj->Aref, ztfttr_obj->lda);

    /* Compute libflame's Lapacke o/p  */
    ztfttr_obj->info = LAPACKE_ztfttr( ztfttr_obj->matrix_layout,ztfttr_obj->transr, ztfttr_obj->uplo, ztfttr_obj->n, 
										(const lapack_complex_double*)ztfttr_obj->ap,  ztfttr_obj->A, ztfttr_obj->lda);

    if( ztfttr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ztfttr is wrong\n", ztfttr_obj->info );
    }
    if( ztfttr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ztfttr is wrong\n", 
        ztfttr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    ztfttr_obj->diff =  computeDiff_z( ztfttr_obj->bufsize_a, 
                ztfttr_obj->A, ztfttr_obj->Aref );

}

TEST_F(ztfttr_test, ztfttr1) {
    EXPECT_NEAR(0.0, ztfttr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ztfttr_test, ztfttr2) {
    EXPECT_NEAR(0.0, ztfttr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ztfttr_test, ztfttr3) {
    EXPECT_NEAR(0.0, ztfttr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ztfttr_test, ztfttr4) {
    EXPECT_NEAR(0.0, ztfttr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}