#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"

#define tpttf_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (arf != NULL) free(arf); \
if (arfref != NULL) free(arfref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin float_common_parameters  class definition */
class tpttf_float_parameters{

   public:
	int bufsize_a;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	float* A;	
	char uplo;
	char transr;
	/*Output Parameter*/
	float *arf;
	float *Aref, *arfref;
	/*Return Values*/	
	lapack_int info, inforef;
   public:
      tpttf_float_parameters (int matrix_layout,  char transr, char uplo,  lapack_int n);
      ~tpttf_float_parameters ();

};

/* Constructor definition  float_common_parameters */
tpttf_float_parameters:: tpttf_float_parameters (int matrix_layout_i, char transr_i, char uplo_i, lapack_int n_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;	
	transr = transr_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n tpttf float: matrix_layout = %d, uplo:%c, n: %d,  transr:%c \n", matrix_layout, uplo, n,transr);
	#endif


	bufsize_a = (n*(n +1)/2);
	
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&A, &Aref, bufsize_a);	
	lapacke_gtest_alloc_float_buffer_pair(&arf, &arfref, bufsize_a);
	if ((A==NULL) || (Aref==NULL) || \
		(arf == NULL) || (arfref == NULL)){
		EXPECT_FALSE( true) << "tpttf_float_parameters object: malloc error.";
		tpttf_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand(A, Aref, bufsize_a);
} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
tpttf_float_parameters :: ~tpttf_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" tpttf_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   tpttf_free();

}
/*  Test fixture class definition */
class stpttf_test  : public  ::testing::Test {
public:
   tpttf_float_parameters  *stpttf_obj;
   void SetUp();
   void TearDown () { delete stpttf_obj;}
};

void stpttf_test::SetUp(){

	 /* LAPACKE stpttf prototype */
    typedef int (*Fptr_NL_LAPACKE_stpttf) (int matrix_layout , char transr , char uplo , lapack_int n , const float * ap , float * arf);

    Fptr_NL_LAPACKE_stpttf stpttf;

    stpttf_obj = new tpttf_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].n);
						   

    idx = Circular_Increment_Index(idx);

    stpttf_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    stpttf_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(stpttf_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(stpttf_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*stpttf library call */
    stpttf = (Fptr_NL_LAPACKE_stpttf)dlsym(stpttf_obj->hModule, "LAPACKE_stpttf");
    ASSERT_TRUE(stpttf != NULL) << "failed to get the Netlib LAPACKE_stpttf symbol";
    
/*Compute stpttf's  o/p */
    stpttf_obj->inforef = stpttf( stpttf_obj->matrix_layout, stpttf_obj->transr, stpttf_obj->uplo,  stpttf_obj->n,
								(const float *)stpttf_obj->Aref, stpttf_obj->arfref);

    /* Compute libflame's Lapacke o/p  */
    stpttf_obj->info = LAPACKE_stpttf( stpttf_obj->matrix_layout, stpttf_obj->transr, stpttf_obj->uplo, stpttf_obj->n, 
										(const float *)stpttf_obj->A,  stpttf_obj->arf);

    if( stpttf_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_stpttf is wrong\n", stpttf_obj->info );
    }
    if( stpttf_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_stpttf is wrong\n", 
        stpttf_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    stpttf_obj->diff =  computeDiff_s( stpttf_obj->bufsize_a, 
                stpttf_obj->arf, stpttf_obj->arfref );

}

TEST_F(stpttf_test, stpttf1) {
    EXPECT_NEAR(0.0, stpttf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(stpttf_test, stpttf2) {
    EXPECT_NEAR(0.0, stpttf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(stpttf_test, stpttf3) {
    EXPECT_NEAR(0.0, stpttf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(stpttf_test, stpttf4) {
    EXPECT_NEAR(0.0, stpttf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin double_common_parameters  class definition */
class tpttf_double_parameters{

   public:
	int bufsize_a;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	double* A;	
	char uplo;
	char transr;
	/*Output Parameter*/
	double *arf;
	double *Aref, *arfref;
	/*Return Values*/	
	lapack_int info, inforef;
   public:
      tpttf_double_parameters (int matrix_layout,  char transr, char uplo,  lapack_int n);
      ~tpttf_double_parameters ();

};

/* Constructor definition  float_common_parameters */
tpttf_double_parameters:: tpttf_double_parameters (int matrix_layout_i, char transr_i, char uplo_i, lapack_int n_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;	
	transr = transr_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n tpttf float: matrix_layout = %d, uplo:%c, n: %d,  transr:%c \n", matrix_layout, uplo, n, transr);
	#endif	

	bufsize_a = (n*(n +1)/2);
	
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&A, &Aref, bufsize_a);	
	lapacke_gtest_alloc_double_buffer_pair(&arf, &arfref, bufsize_a);
	if ((A==NULL) || (Aref==NULL) || \
		(arf == NULL) || (arfref == NULL)){
		EXPECT_FALSE( true) << "tpttf_float_parameters object: malloc error.";
		tpttf_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand(A, Aref, bufsize_a);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
tpttf_double_parameters :: ~tpttf_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" tpttf_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   tpttf_free();

}
/*  Test fixture class definition */
class dtpttf_test  : public  ::testing::Test {
public:
   tpttf_double_parameters  *dtpttf_obj;
   void SetUp();
   void TearDown () { delete dtpttf_obj;}
};

void dtpttf_test::SetUp(){

	 /* LAPACKE dtpttf prototype */
    typedef int (*Fptr_NL_LAPACKE_dtpttf) (int matrix_layout , char transr , char uplo , lapack_int n , const double * ap , double * arf);

    Fptr_NL_LAPACKE_dtpttf dtpttf;

    dtpttf_obj = new tpttf_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].n);
						   

    idx = Circular_Increment_Index(idx);

    dtpttf_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dtpttf_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dtpttf_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dtpttf_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*dtpttf library call */
    dtpttf = (Fptr_NL_LAPACKE_dtpttf)dlsym(dtpttf_obj->hModule, "LAPACKE_dtpttf");
    ASSERT_TRUE(dtpttf != NULL) << "failed to get the Netlib LAPACKE_dtpttf symbol";
    
/*Compute dtpttf's  o/p */
    dtpttf_obj->inforef = dtpttf( dtpttf_obj->matrix_layout, dtpttf_obj->transr, dtpttf_obj->uplo,  dtpttf_obj->n,
								(const double *)dtpttf_obj->Aref, dtpttf_obj->arfref);

    /* Compute libflame's Lapacke o/p  */
    dtpttf_obj->info = LAPACKE_dtpttf( dtpttf_obj->matrix_layout, dtpttf_obj->transr, dtpttf_obj->uplo, dtpttf_obj->n, 
										(const double *)dtpttf_obj->A,  dtpttf_obj->arf);

    if( dtpttf_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dtpttf is wrong\n", dtpttf_obj->info );
    }
    if( dtpttf_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dtpttf is wrong\n", 
        dtpttf_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dtpttf_obj->diff =  computeDiff_d( dtpttf_obj->bufsize_a, 
                dtpttf_obj->arf, dtpttf_obj->arfref );

}

TEST_F(dtpttf_test, dtpttf1) {
    EXPECT_NEAR(0.0, dtpttf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dtpttf_test, dtpttf2) {
    EXPECT_NEAR(0.0, dtpttf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dtpttf_test, dtpttf3) {
    EXPECT_NEAR(0.0, dtpttf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dtpttf_test, dtpttf4) {
    EXPECT_NEAR(0.0, dtpttf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}
/* Begin scomplex_common_parameters  class definition */
class tpttf_scomplex_parameters{

   public:
	int bufsize_a;
	int bufsize_arf;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_complex_float* A;
	char uplo;
	char transr;
	/*Output Parameter*/
	lapack_complex_float *arf;
	lapack_complex_float *Aref, *arfref;
	/*Return Values*/	
	lapack_int info, inforef;
   public:
      tpttf_scomplex_parameters (int matrix_layout,  char transr, char uplo,  lapack_int n);
      ~tpttf_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
tpttf_scomplex_parameters:: tpttf_scomplex_parameters (int matrix_layout_i, char transr_i, char uplo_i, lapack_int n_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;	
	transr = transr_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n tpttf float: matrix_layout = %d, uplo:%c, n: %d,  transr:%c \n", matrix_layout, uplo, n, transr);
	#endif	

	bufsize_a = (n*(n +1)/2);
	
	if (transr == 'T')
		transr = 'C';
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize_a);	
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&arf, &arfref, bufsize_a);
	if ((A==NULL) || (Aref==NULL) || \
		(arf == NULL) || (arfref == NULL)){
		EXPECT_FALSE( true) << "tpttf_float_parameters object: malloc error.";
		tpttf_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand(A, Aref, bufsize_a);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
tpttf_scomplex_parameters :: ~tpttf_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" tpttf_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   tpttf_free();

}
/*  Test fixture class definition */
class ctpttf_test  : public  ::testing::Test {
public:
   tpttf_scomplex_parameters  *ctpttf_obj;
   void SetUp();
   void TearDown () { delete ctpttf_obj;}
};

void ctpttf_test::SetUp(){

	 /* LAPACKE ctpttf prototype */
    typedef int (*Fptr_NL_LAPACKE_ctpttf) (int matrix_layout , char transr , char uplo , lapack_int n , const lapack_complex_float * ap , lapack_complex_float * arf);

    Fptr_NL_LAPACKE_ctpttf ctpttf;

    ctpttf_obj = new tpttf_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].n);
						   

    idx = Circular_Increment_Index(idx);

    ctpttf_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ctpttf_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ctpttf_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ctpttf_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*ctpttf library call */
    ctpttf = (Fptr_NL_LAPACKE_ctpttf)dlsym(ctpttf_obj->hModule, "LAPACKE_ctpttf");
    ASSERT_TRUE(ctpttf != NULL) << "failed to get the Netlib LAPACKE_ctpttf symbol";
    
/*Compute ctpttf's  o/p */
    ctpttf_obj->inforef = ctpttf( ctpttf_obj->matrix_layout, ctpttf_obj->transr, ctpttf_obj->uplo,  ctpttf_obj->n,
								(const lapack_complex_float *)ctpttf_obj->Aref, ctpttf_obj->arfref);

    /* Compute libflame's Lapacke o/p  */
    ctpttf_obj->info = LAPACKE_ctpttf( ctpttf_obj->matrix_layout, ctpttf_obj->transr, ctpttf_obj->uplo, ctpttf_obj->n, 
										(const lapack_complex_float *)ctpttf_obj->A,  ctpttf_obj->arf);

    if( ctpttf_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ctpttf is wrong\n", ctpttf_obj->info );
    }
    if( ctpttf_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ctpttf is wrong\n", 
        ctpttf_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    ctpttf_obj->diff =  computeDiff_c( ctpttf_obj->bufsize_a, 
                ctpttf_obj->arf, ctpttf_obj->arfref );

}

TEST_F(ctpttf_test, ctpttf1) {
    EXPECT_NEAR(0.0, ctpttf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ctpttf_test, ctpttf2) {
    EXPECT_NEAR(0.0, ctpttf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ctpttf_test, ctpttf3) {
    EXPECT_NEAR(0.0, ctpttf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ctpttf_test, ctpttf4) {
    EXPECT_NEAR(0.0, ctpttf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class tpttf_dcomplex_parameters{

   public:
	int bufsize_a;
	int bufsize_arf;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_complex_double* A;	
	char uplo;
	char transr;
	/*Output Parameter*/
	lapack_complex_double *arf;
	lapack_complex_double *Aref, *arfref;
	/*Return Values*/	
	lapack_int info, inforef;
   public:
      tpttf_dcomplex_parameters (int matrix_layout,  char transr, char uplo,  lapack_int n);
      ~tpttf_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
tpttf_dcomplex_parameters:: tpttf_dcomplex_parameters (int matrix_layout_i, char transr_i, char uplo_i, lapack_int n_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;	
	transr = transr_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n tpttf float: matrix_layout = %d, uplo:%c, n: %d, lda:%d, transr:%c \n", matrix_layout, uplo, n, transr);
	#endif	

	bufsize_a = (n*(n +1)/2);
	
	if (transr == 'T')
		transr = 'C';
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize_a);	
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&arf, &arfref, bufsize_a);
	if ((A==NULL) || (Aref==NULL) || \
		(arf == NULL) || (arfref == NULL)){
		EXPECT_FALSE( true) << "tpttf_double_parameters object: malloc error.";
		tpttf_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand(A, Aref, bufsize_a);

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
tpttf_dcomplex_parameters :: ~tpttf_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" tpttf_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   tpttf_free();

}
/*  Test fixture class definition */
class ztpttf_test  : public  ::testing::Test {
public:
   tpttf_dcomplex_parameters  *ztpttf_obj;
   void SetUp();
   void TearDown () { delete ztpttf_obj;}
};

void ztpttf_test::SetUp(){

	 /* LAPACKE ztpttf prototype */
    typedef int (*Fptr_NL_LAPACKE_ztpttf) (int matrix_layout , char transr , char uplo , lapack_int n , const lapack_complex_double * ap , lapack_complex_double * arf );

    Fptr_NL_LAPACKE_ztpttf ztpttf;

    ztpttf_obj = new tpttf_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].n);
						   

    idx = Circular_Increment_Index(idx);

    ztpttf_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ztpttf_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ztpttf_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ztpttf_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*ztpttf library call */
    ztpttf = (Fptr_NL_LAPACKE_ztpttf)dlsym(ztpttf_obj->hModule, "LAPACKE_ztpttf");
    ASSERT_TRUE(ztpttf != NULL) << "failed to get the Netlib LAPACKE_ztpttf symbol";
    
/*Compute ztpttf's  o/p */
    ztpttf_obj->inforef = ztpttf( ztpttf_obj->matrix_layout, ztpttf_obj->transr, ztpttf_obj->uplo,  ztpttf_obj->n,
								(const lapack_complex_double *)ztpttf_obj->Aref, ztpttf_obj->arfref);

    /* Compute libflame's Lapacke o/p  */
    ztpttf_obj->info = LAPACKE_ztpttf( ztpttf_obj->matrix_layout, ztpttf_obj->transr, ztpttf_obj->uplo, ztpttf_obj->n, 
										(const lapack_complex_double *)ztpttf_obj->A, ztpttf_obj->arf);

    if( ztpttf_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ztpttf is wrong\n", ztpttf_obj->info );
    }
    if( ztpttf_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ztpttf is wrong\n", 
        ztpttf_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    ztpttf_obj->diff =  computeDiff_z( ztpttf_obj->bufsize_a, 
                ztpttf_obj->arf, ztpttf_obj->arfref );

}

TEST_F(ztpttf_test, ztpttf1) {
    EXPECT_NEAR(0.0, ztpttf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ztpttf_test, ztpttf2) {
    EXPECT_NEAR(0.0, ztpttf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ztpttf_test, ztpttf3) {
    EXPECT_NEAR(0.0, ztpttf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ztpttf_test, ztpttf4) {
    EXPECT_NEAR(0.0, ztpttf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}