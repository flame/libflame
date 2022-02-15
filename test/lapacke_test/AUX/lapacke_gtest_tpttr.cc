#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"

#define tpttr_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (ap!=NULL)    free(ap); \
if (apref!=NULL) free(apref);\
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
class tpttr_float_parameters{

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
      tpttr_float_parameters (int matrix_layout, char uplo, lapack_int n, lapack_int lda);
      ~tpttr_float_parameters ();

};

/* Constructor definition  float_common_parameters */
tpttr_float_parameters:: tpttr_float_parameters (int matrix_layout_i, char uplo_i, lapack_int n_i, lapack_int lda_i)
{
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n tpttr float: matrix_layout = %d, uplo:%c, n: %d \n", matrix_layout, uplo, n);
	#endif	
	
	bufsize_a = lda *n;
	bufsize_ap = (n*(n+1)/2);
	
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_float_buffer_pair(&ap, &apref, bufsize_ap);
	if ((A==NULL) || (Aref==NULL) || \
		(ap == NULL) || (apref == NULL)){
		EXPECT_FALSE( true) << "tpttr_float_parameters object: malloc error.";
		tpttr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand_custom_matrix( A, Aref, lda, n, uplo);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
tpttr_float_parameters :: ~tpttr_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" tpttr_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   tpttr_free();

}
/*  Test fixture class definition */
class stpttr_test  : public  ::testing::Test {
public:
   tpttr_float_parameters  *stpttr_obj;
   void SetUp();
   void TearDown () { delete stpttr_obj;}
};

void stpttr_test::SetUp(){

	 /* LAPACKE stpttr prototype */
    typedef int (*Fptr_NL_LAPACKE_stpttr) (int matrix_layout , char uplo , lapack_int n , const float * ap , float * A , lapack_int lda);

    Fptr_NL_LAPACKE_stpttr stpttr;

    stpttr_obj = new tpttr_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].n, 
						   lin_solver_paramslist[idx].ldb);
						   

    idx = Circular_Increment_Index(idx);

    stpttr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    stpttr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(stpttr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(stpttr_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*stpttr library call */
    stpttr = (Fptr_NL_LAPACKE_stpttr)dlsym(stpttr_obj->hModule, "LAPACKE_stpttr");
    ASSERT_TRUE(stpttr != NULL) << "failed to get the Netlib LAPACKE_stpttr symbol";
    
/*Compute stpttr's  o/p */
    stpttr_obj->inforef = stpttr( stpttr_obj->matrix_layout,  stpttr_obj->uplo,  stpttr_obj->n,
								 (const float*)stpttr_obj->apref, stpttr_obj->Aref, stpttr_obj->lda);

    /* Compute libflame's Lapacke o/p  */
    stpttr_obj->info = LAPACKE_stpttr( stpttr_obj->matrix_layout,stpttr_obj->uplo, stpttr_obj->n, 
										(const float*)stpttr_obj->ap,  stpttr_obj->A, stpttr_obj->lda);

    if( stpttr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_stpttr is wrong\n", stpttr_obj->info );
    }
    if( stpttr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_stpttr is wrong\n", 
        stpttr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    stpttr_obj->diff =  computeDiff_s( stpttr_obj->bufsize_a, 
                stpttr_obj->A, stpttr_obj->Aref );

}

TEST_F(stpttr_test, stpttr1) {
    EXPECT_NEAR(0.0, stpttr_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(stpttr_test, stpttr2) {
    EXPECT_NEAR(0.0, stpttr_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(stpttr_test, stpttr3) {
    EXPECT_NEAR(0.0, stpttr_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(stpttr_test, stpttr4) {
    EXPECT_NEAR(0.0, stpttr_obj->diff, LAPACKE_EIG_THRESHOLD);
}

/* Begin double_common_parameters  class definition */

class tpttr_double_parameters{

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
      tpttr_double_parameters (int matrix_layout, char uplo, lapack_int n, lapack_int lda);
      ~tpttr_double_parameters ();

};

/* Constructor definition  double_common_parameters */
tpttr_double_parameters:: tpttr_double_parameters (int matrix_layout_i, char uplo_i, lapack_int n_i, lapack_int lda_i)
{
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n tpttr double: matrix_layout = %d, uplo:%c, n: %d \n", matrix_layout, uplo, n);
	#endif	
	
	bufsize_a = lda *n;
	bufsize_ap = (n*(n+1)/2);
	
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_double_buffer_pair(&ap, &apref, bufsize_ap);
	if ((A==NULL) || (Aref==NULL) || \
		(ap == NULL) || (apref == NULL)){
		EXPECT_FALSE( true) << "tpttr_double_parameters object: malloc error.";
		tpttr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand_custom_matrix( A, Aref, lda, n  , uplo);
	

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
tpttr_double_parameters :: ~tpttr_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" tpttr_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   tpttr_free();

}
/*  Test fixture class definition */
class dtpttr_test  : public  ::testing::Test {
public:
   tpttr_double_parameters  *dtpttr_obj;
   void SetUp();
   void TearDown () { delete dtpttr_obj;}
};

void dtpttr_test::SetUp(){

	 /* LAPACKE dtpttr prototype */
    typedef int (*Fptr_NL_LAPACKE_dtpttr) (int matrix_layout , char uplo , lapack_int n , const double * ap , double * A , lapack_int lda );

    Fptr_NL_LAPACKE_dtpttr dtpttr;

    dtpttr_obj = new tpttr_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].n,
						   lin_solver_paramslist[idx].ldb);
						   

    idx = Circular_Increment_Index(idx);

    dtpttr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dtpttr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dtpttr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dtpttr_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*dtpttr library call */
    dtpttr = (Fptr_NL_LAPACKE_dtpttr)dlsym(dtpttr_obj->hModule, "LAPACKE_dtpttr");
    ASSERT_TRUE(dtpttr != NULL) << "failed to get the Netlib LAPACKE_dtpttr symbol";
    
/*Compute dtpttr's  o/p */
    dtpttr_obj->inforef = dtpttr( dtpttr_obj->matrix_layout,  dtpttr_obj->uplo,  dtpttr_obj->n,
								 (const double*)dtpttr_obj->apref, dtpttr_obj->Aref, dtpttr_obj->lda);

    /* Compute libflame's Lapacke o/p  */
    dtpttr_obj->info = LAPACKE_dtpttr( dtpttr_obj->matrix_layout,dtpttr_obj->uplo, dtpttr_obj->n, 
										(const double*)dtpttr_obj->ap,  dtpttr_obj->A, dtpttr_obj->lda);

    if( dtpttr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dtpttr is wrong\n", dtpttr_obj->info );
    }
    if( dtpttr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dtpttr is wrong\n", 
        dtpttr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dtpttr_obj->diff =  computeDiff_d( dtpttr_obj->bufsize_a, 
                dtpttr_obj->A, dtpttr_obj->Aref );

}

TEST_F(dtpttr_test, dtpttr1) {
    EXPECT_NEAR(0.0, dtpttr_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(dtpttr_test, dtpttr2) {
    EXPECT_NEAR(0.0, dtpttr_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(dtpttr_test, dtpttr3) {
    EXPECT_NEAR(0.0, dtpttr_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(dtpttr_test, dtpttr4) {
    EXPECT_NEAR(0.0, dtpttr_obj->diff, LAPACKE_EIG_THRESHOLD);
}

/* Begin scomplex_common_parameters  class definition */

class tpttr_scomplex_parameters{

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
      tpttr_scomplex_parameters (int matrix_layout, char uplo, lapack_int n, lapack_int lda);
      ~tpttr_scomplex_parameters ();

};

/* Constructor definition  scomplex_common_parameters */
tpttr_scomplex_parameters:: tpttr_scomplex_parameters (int matrix_layout_i, char uplo_i, lapack_int n_i, lapack_int lda_i)
{
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n tpttr scomplex: matrix_layout = %d, uplo:%c, n: %d \n", matrix_layout, uplo, n);
	#endif	
	
	bufsize_a = lda *n;
	bufsize_ap = (n*(n+1)/2);
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&ap, &apref, bufsize_ap);
	if ((A==NULL) || (Aref==NULL) || \
		(ap == NULL) || (apref == NULL)){
		EXPECT_FALSE( true) << "tpttr_scomplex_parameters object: malloc error.";
		tpttr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix( A, Aref, lda, n, uplo);


} /* end of Constructor  */

/* Destructor definition  'scomplex_common_parameters' class  */
tpttr_scomplex_parameters :: ~tpttr_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" tpttr_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   tpttr_free();

}
/*  Test fixture class definition */
class ctpttr_test  : public  ::testing::Test {
public:
   tpttr_scomplex_parameters  *ctpttr_obj;
   void SetUp();
   void TearDown () { delete ctpttr_obj;}
};

void ctpttr_test::SetUp(){

	 /* LAPACKE ctpttr prototype */
    typedef int (*Fptr_NL_LAPACKE_ctpttr) (int matrix_layout , char uplo , lapack_int n , const lapack_complex_float * ap , lapack_complex_float * A , lapack_int lda);

    Fptr_NL_LAPACKE_ctpttr ctpttr;

    ctpttr_obj = new tpttr_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].n,
						   lin_solver_paramslist[idx].ldb);
						   

    idx = Circular_Increment_Index(idx);

    ctpttr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ctpttr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ctpttr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ctpttr_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*ctpttr library call */
    ctpttr = (Fptr_NL_LAPACKE_ctpttr)dlsym(ctpttr_obj->hModule, "LAPACKE_ctpttr");
    ASSERT_TRUE(ctpttr != NULL) << "failed to get the Netlib LAPACKE_ctpttr symbol";
    
/*Compute ctpttr's  o/p */
    ctpttr_obj->inforef = ctpttr( ctpttr_obj->matrix_layout,  ctpttr_obj->uplo,  ctpttr_obj->n,
								 (const lapack_complex_float*)ctpttr_obj->ap, ctpttr_obj->Aref, ctpttr_obj->lda);

    /* Compute libflame's Lapacke o/p  */
    ctpttr_obj->info = LAPACKE_ctpttr( ctpttr_obj->matrix_layout,ctpttr_obj->uplo, ctpttr_obj->n, 
										 (const lapack_complex_float*)ctpttr_obj->ap,  ctpttr_obj->A, ctpttr_obj->lda);

    if( ctpttr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ctpttr is wrong\n", ctpttr_obj->info );
    }
    if( ctpttr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ctpttr is wrong\n", 
        ctpttr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    ctpttr_obj->diff =  computeDiff_c( ctpttr_obj->bufsize_a, 
                ctpttr_obj->A, ctpttr_obj->Aref );

}

TEST_F(ctpttr_test, ctpttr1) {
    EXPECT_NEAR(0.0, ctpttr_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(ctpttr_test, ctpttr2) {
    EXPECT_NEAR(0.0, ctpttr_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(ctpttr_test, ctpttr3) {
    EXPECT_NEAR(0.0, ctpttr_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(ctpttr_test, ctpttr4) {
    EXPECT_NEAR(0.0, ctpttr_obj->diff, LAPACKE_EIG_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */

class tpttr_dcomplex_parameters{

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
      tpttr_dcomplex_parameters (int matrix_layout, char uplo, lapack_int n, lapack_int lda);
      ~tpttr_dcomplex_parameters ();

};

/* Constructor definition  dcomplex_common_parameters */
tpttr_dcomplex_parameters:: tpttr_dcomplex_parameters (int matrix_layout_i, char uplo_i, lapack_int n_i, lapack_int lda_i)
{
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n tpttr dcomplex: matrix_layout = %d, uplo:%c, n: %d \n", matrix_layout, uplo, n);
	#endif	
	
	bufsize_a = lda *n;
	bufsize_ap = (n*(n+1)/2);
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&ap, &apref, bufsize_ap);
	if ((A==NULL) || (Aref==NULL) || \
		(ap == NULL) || (apref == NULL)){
		EXPECT_FALSE( true) << "tpttr_dcomplex_parameters object: malloc error.";
		tpttr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix( A, Aref, lda, n, uplo);

} /* end of Constructor  */

/* Destructor definition  'dcomplex_common_parameters' class  */
tpttr_dcomplex_parameters :: ~tpttr_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" tpttr_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   tpttr_free();

}
/*  Test fixture class definition */
class ztpttr_test  : public  ::testing::Test {
public:
   tpttr_dcomplex_parameters  *ztpttr_obj;
   void SetUp();
   void TearDown () { delete ztpttr_obj;}
};

void ztpttr_test::SetUp(){

	 /* LAPACKE ztpttr prototype */
    typedef int (*Fptr_NL_LAPACKE_ztpttr) (int matrix_layout , char uplo , lapack_int n , const lapack_complex_double * ap , lapack_complex_double * A , lapack_int lda);

    Fptr_NL_LAPACKE_ztpttr ztpttr;

    ztpttr_obj = new tpttr_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].n, 
						   lin_solver_paramslist[idx].ldb);
						   

    idx = Circular_Increment_Index(idx);

    ztpttr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ztpttr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ztpttr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ztpttr_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*ztpttr library call */
    ztpttr = (Fptr_NL_LAPACKE_ztpttr)dlsym(ztpttr_obj->hModule, "LAPACKE_ztpttr");
    ASSERT_TRUE(ztpttr != NULL) << "failed to get the Netlib LAPACKE_ztpttr symbol";
    
/*Compute ztpttr's  o/p */
    ztpttr_obj->inforef = ztpttr( ztpttr_obj->matrix_layout,  ztpttr_obj->uplo,  ztpttr_obj->n,
								 (const lapack_complex_double*)ztpttr_obj->apref, ztpttr_obj->Aref, ztpttr_obj->lda);

    /* Compute libflame's Lapacke o/p  */
    ztpttr_obj->info = LAPACKE_ztpttr( ztpttr_obj->matrix_layout,ztpttr_obj->uplo, ztpttr_obj->n, 
										(const lapack_complex_double*)ztpttr_obj->ap,  ztpttr_obj->A, ztpttr_obj->lda);

    if( ztpttr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ztpttr is wrong\n", ztpttr_obj->info );
    }
    if( ztpttr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ztpttr is wrong\n", 
        ztpttr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    ztpttr_obj->diff =  computeDiff_z( ztpttr_obj->bufsize_a, 
                ztpttr_obj->A, ztpttr_obj->Aref );

}

TEST_F(ztpttr_test, ztpttr1) {
    EXPECT_NEAR(0.0, ztpttr_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(ztpttr_test, ztpttr2) {
    EXPECT_NEAR(0.0, ztpttr_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(ztpttr_test, ztpttr3) {
    EXPECT_NEAR(0.0, ztpttr_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(ztpttr_test, ztpttr4) {
    EXPECT_NEAR(0.0, ztpttr_obj->diff, LAPACKE_EIG_THRESHOLD);
}
