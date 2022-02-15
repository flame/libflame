#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"

#define tfttp_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (ap!=NULL)    free(ap); \
if (apref!=NULL) free(apref);\
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
class tfttp_float_parameters{

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
	char uplo, transr;
	/*Output Parameter*/
	float *Aref, *apref;
	/*Return Values*/	
	lapack_int info, inforef;
   public:
      tfttp_float_parameters (int matrix_layout, char transr, char uplo, lapack_int n);
      ~tfttp_float_parameters ();

};

/* Constructor definition  float_common_parameters */
tfttp_float_parameters:: tfttp_float_parameters (int matrix_layout_i, char transr_i, char uplo_i, lapack_int n_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	transr = transr_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n tfttp float: matrix_layout = %d, uplo:%c, n: %d \n", matrix_layout, uplo, n);
	#endif	
	
	bufsize_a = (n*(n+1)/2);
	bufsize_ap = (n*(n+1)/2);
	
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_float_buffer_pair(&ap, &apref, bufsize_ap);
	if ((A==NULL) || (Aref==NULL) || \
		(ap == NULL) || (apref == NULL)){
		EXPECT_FALSE( true) << "tfttp_float_parameters object: malloc error.";
		tfttp_free();
		exit(0);
	}
	/* Initialization of input matrices */
	//lapacke_gtest_init_float_buffer_pair_rand_custom_matrix( ap, apref, bufsize_ap, n, uplo);
	lapacke_gtest_init_float_buffer_pair_rand( ap, apref, bufsize_ap);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
tfttp_float_parameters :: ~tfttp_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" tfttp_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   tfttp_free();

}
/*  Test fixture class definition */
class stfttp_test  : public  ::testing::Test {
public:
   tfttp_float_parameters  *stfttp_obj;
   void SetUp();
   void TearDown () { delete stfttp_obj;}
};

void stfttp_test::SetUp(){

	 /* LAPACKE stfttp prototype */
    typedef int (*Fptr_NL_LAPACKE_stfttp) (int matrix_layout , char transr , char uplo , lapack_int n , const float * arf , float * ap);

    Fptr_NL_LAPACKE_stfttp stfttp;

    stfttp_obj = new tfttp_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
							lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].n);
						   

    idx = Circular_Increment_Index(idx);

    stfttp_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    stfttp_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(stfttp_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(stfttp_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*stfttp library call */
    stfttp = (Fptr_NL_LAPACKE_stfttp)dlsym(stfttp_obj->hModule, "LAPACKE_stfttp");
    ASSERT_TRUE(stfttp != NULL) << "failed to get the Netlib LAPACKE_stfttp symbol";
    
/*Compute stfttp's  o/p */
    stfttp_obj->inforef = stfttp( stfttp_obj->matrix_layout, stfttp_obj->transr, stfttp_obj->uplo,  stfttp_obj->n,
								 (const float*)stfttp_obj->Aref, stfttp_obj->apref);

    /* Compute libflame's Lapacke o/p  */
    stfttp_obj->info = LAPACKE_stfttp( stfttp_obj->matrix_layout, stfttp_obj->transr, stfttp_obj->uplo, stfttp_obj->n, 
										(const float*)stfttp_obj->A,  stfttp_obj->ap);

    if( stfttp_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_stfttp is wrong\n", stfttp_obj->info );
    }
    if( stfttp_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_stfttp is wrong\n", 
        stfttp_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    stfttp_obj->diff =  computeDiff_s( stfttp_obj->bufsize_a, 
                stfttp_obj->A, stfttp_obj->Aref );

}

TEST_F(stfttp_test, stfttp1) {
    EXPECT_NEAR(0.0, stfttp_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(stfttp_test, stfttp2) {
    EXPECT_NEAR(0.0, stfttp_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(stfttp_test, stfttp3) {
    EXPECT_NEAR(0.0, stfttp_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(stfttp_test, stfttp4) {
    EXPECT_NEAR(0.0, stfttp_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin double_common_parameters  class definition */

class tfttp_double_parameters{

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
	char uplo, transr;
	/*Output Parameter*/
	double *Aref, *apref;
	/*Return Values*/	
	lapack_int info, inforef;
   public:
      tfttp_double_parameters (int matrix_layout, char transr, char uplo, lapack_int n);
      ~tfttp_double_parameters ();

};

/* Constructor definition  double_common_parameters */
tfttp_double_parameters:: tfttp_double_parameters (int matrix_layout_i, char transr_i, char uplo_i, lapack_int n_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	transr = transr_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n tfttp double: matrix_layout = %d, uplo:%c, n: %d \n", matrix_layout, uplo, n);
	#endif	
	
	bufsize_a = (n*(n+1)/2);
	bufsize_ap = (n*(n+1)/2);
	
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_double_buffer_pair(&ap, &apref, bufsize_ap);
	if ((A==NULL) || (Aref==NULL) || \
		(ap == NULL) || (apref == NULL)){
		EXPECT_FALSE( true) << "tfttp_double_parameters object: malloc error.";
		tfttp_free();
		exit(0);
	}
	/* Initialization of input matrices */
	//lapacke_gtest_init_double_buffer_pair_rand_custom_matrix( ap, apref, bufsize_ap, n, uplo);
	lapacke_gtest_init_double_buffer_pair_rand( ap, apref, bufsize_ap);

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
tfttp_double_parameters :: ~tfttp_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" tfttp_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   tfttp_free();

}
/*  Test fixture class definition */
class dtfttp_test  : public  ::testing::Test {
public:
   tfttp_double_parameters  *dtfttp_obj;
   void SetUp();
   void TearDown () { delete dtfttp_obj;}
};

void dtfttp_test::SetUp(){

	 /* LAPACKE dtfttp prototype */
    typedef int (*Fptr_NL_LAPACKE_dtfttp) (int matrix_layout , char transr , char uplo , lapack_int n , const double * arf , double * ap);

    Fptr_NL_LAPACKE_dtfttp dtfttp;

    dtfttp_obj = new tfttp_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
							lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].n);
						   

    idx = Circular_Increment_Index(idx);

    dtfttp_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dtfttp_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dtfttp_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dtfttp_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*dtfttp library call */
    dtfttp = (Fptr_NL_LAPACKE_dtfttp)dlsym(dtfttp_obj->hModule, "LAPACKE_dtfttp");
    ASSERT_TRUE(dtfttp != NULL) << "failed to get the Netlib LAPACKE_dtfttp symbol";
    
/*Compute dtfttp's  o/p */
    dtfttp_obj->inforef = dtfttp( dtfttp_obj->matrix_layout, dtfttp_obj->transr, dtfttp_obj->uplo,  dtfttp_obj->n,
								 (const double*)dtfttp_obj->Aref, dtfttp_obj->apref);

    /* Compute libflame's Lapacke o/p  */
    dtfttp_obj->info = LAPACKE_dtfttp( dtfttp_obj->matrix_layout, dtfttp_obj->transr, dtfttp_obj->uplo, dtfttp_obj->n, 
										(const double*)dtfttp_obj->A,  dtfttp_obj->ap);

    if( dtfttp_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dtfttp is wrong\n", dtfttp_obj->info );
    }
    if( dtfttp_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dtfttp is wrong\n", 
        dtfttp_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dtfttp_obj->diff =  computeDiff_d( dtfttp_obj->bufsize_a, 
                dtfttp_obj->A, dtfttp_obj->Aref );

}

TEST_F(dtfttp_test, dtfttp1) {
    EXPECT_NEAR(0.0, dtfttp_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dtfttp_test, dtfttp2) {
    EXPECT_NEAR(0.0, dtfttp_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dtfttp_test, dtfttp3) {
    EXPECT_NEAR(0.0, dtfttp_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dtfttp_test, dtfttp4) {
    EXPECT_NEAR(0.0, dtfttp_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin scomplex_common_parameters  class definition */
class tfttp_scomplex_parameters{

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
	char uplo, transr;
	/*Output Parameter*/
	lapack_complex_float *Aref, *apref;
	/*Return Values*/	
	lapack_int info, inforef;
   public:
      tfttp_scomplex_parameters (int matrix_layout, char transr, char uplo, lapack_int n);
      ~tfttp_scomplex_parameters ();

};

/* Constructor definition  scomplex_common_parameters */
tfttp_scomplex_parameters:: tfttp_scomplex_parameters (int matrix_layout_i, char transr_i, char uplo_i, lapack_int n_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	transr = transr_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n tfttp lapack_complex_float: matrix_layout = %d, uplo:%c, n: %d \n", matrix_layout, uplo, n);
	#endif	
	
	if (transr == 'T')
		transr = 'C';
	
	bufsize_a = (n*(n+1)/2);
	bufsize_ap = (n*(n+1)/2);
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&ap, &apref, bufsize_ap);
	if ((A==NULL) || (Aref==NULL) || \
		(ap == NULL) || (apref == NULL)){
		EXPECT_FALSE( true) << "tfttp_lapack_complex_float_parameters object: malloc error.";
		tfttp_free();
		exit(0);
	}
	/* Initialization of input matrices */
	//lapacke_gtest_init_lapack_complex_float_buffer_pair_rand_custom_matrix( ap, apref, bufsize_ap, n, uplo);
	lapacke_gtest_init_scomplex_buffer_pair_rand( ap, apref, bufsize_ap);

} /* end of Constructor  */

/* Destructor definition  'scomplex_common_parameters' class  */
tfttp_scomplex_parameters :: ~tfttp_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" tfttp_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   tfttp_free();

}
/*  Test fixture class definition */
class ctfttp_test  : public  ::testing::Test {
public:
   tfttp_scomplex_parameters  *ctfttp_obj;
   void SetUp();
   void TearDown () { delete ctfttp_obj;}
};

void ctfttp_test::SetUp(){

	 /* LAPACKE ctfttp prototype */
    typedef int (*Fptr_NL_LAPACKE_ctfttp) (int matrix_layout , char transr , char uplo , lapack_int n , const lapack_complex_float * arf , lapack_complex_float * ap);

    Fptr_NL_LAPACKE_ctfttp ctfttp;

    ctfttp_obj = new tfttp_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
							lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].n);
						   

    idx = Circular_Increment_Index(idx);

    ctfttp_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ctfttp_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ctfttp_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ctfttp_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*ctfttp library call */
    ctfttp = (Fptr_NL_LAPACKE_ctfttp)dlsym(ctfttp_obj->hModule, "LAPACKE_ctfttp");
    ASSERT_TRUE(ctfttp != NULL) << "failed to get the Netlib LAPACKE_ctfttp symbol";
    
/*Compute ctfttp's  o/p */
    ctfttp_obj->inforef = ctfttp( ctfttp_obj->matrix_layout, ctfttp_obj->transr, ctfttp_obj->uplo,  ctfttp_obj->n,
								 (const lapack_complex_float*)ctfttp_obj->Aref, ctfttp_obj->apref);

    /* Compute libflame's Lapacke o/p  */
    ctfttp_obj->info = LAPACKE_ctfttp( ctfttp_obj->matrix_layout, ctfttp_obj->transr, ctfttp_obj->uplo, ctfttp_obj->n, 
										(const lapack_complex_float*)ctfttp_obj->A,  ctfttp_obj->ap);

    if( ctfttp_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ctfttp is wrong\n", ctfttp_obj->info );
    }
    if( ctfttp_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ctfttp is wrong\n", 
        ctfttp_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    ctfttp_obj->diff =  computeDiff_c( ctfttp_obj->bufsize_a, 
                ctfttp_obj->A, ctfttp_obj->Aref );

}

TEST_F(ctfttp_test, ctfttp1) {
    EXPECT_NEAR(0.0, ctfttp_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ctfttp_test, ctfttp2) {
    EXPECT_NEAR(0.0, ctfttp_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ctfttp_test, ctfttp3) {
    EXPECT_NEAR(0.0, ctfttp_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ctfttp_test, ctfttp4) {
    EXPECT_NEAR(0.0, ctfttp_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */

class tfttp_dcomplex_parameters{

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
	char uplo, transr;
	/*Output Parameter*/
	lapack_complex_double *Aref, *apref;
	/*Return Values*/	
	lapack_int info, inforef;
   public:
      tfttp_dcomplex_parameters (int matrix_layout, char transr, char uplo, lapack_int n);
      ~tfttp_dcomplex_parameters ();

};

/* Constructor definition  dcomplex_common_parameters */
tfttp_dcomplex_parameters:: tfttp_dcomplex_parameters (int matrix_layout_i, char transr_i, char uplo_i, lapack_int n_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	transr = transr_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n tfttp lapack_complex_double: matrix_layout = %d, uplo:%c, n: %d \n", matrix_layout, uplo, n);
	#endif	
	
	if (transr == 'T')
		transr = 'C';
	
	bufsize_a = (n*(n+1)/2);
	bufsize_ap = (n*(n+1)/2);
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&ap, &apref, bufsize_ap);
	if ((A==NULL) || (Aref==NULL) || \
		(ap == NULL) || (apref == NULL)){
		EXPECT_FALSE( true) << "tfttp_lapack_complex_double_parameters object: malloc error.";
		tfttp_free();
		exit(0);
	}
	/* Initialization of input matrices */
	//lapacke_gtest_init_lapack_complex_double_buffer_pair_rand_custom_matrix( ap, apref, bufsize_ap, n, uplo);
	lapacke_gtest_init_dcomplex_buffer_pair_rand( ap, apref, bufsize_ap);

} /* end of Constructor  */

/* Destructor definition  'dcomplex_common_parameters' class  */
tfttp_dcomplex_parameters :: ~tfttp_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" tfttp_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   tfttp_free();

}
/*  Test fixture class definition */
class ztfttp_test  : public  ::testing::Test {
public:
   tfttp_dcomplex_parameters  *ztfttp_obj;
   void SetUp();
   void TearDown () { delete ztfttp_obj;}
};

void ztfttp_test::SetUp(){

	 /* LAPACKE ztfttp prototype */
    typedef int (*Fptr_NL_LAPACKE_ztfttp) (int matrix_layout , char transr , char uplo , lapack_int n , const lapack_complex_double * arf , lapack_complex_double * ap);

    Fptr_NL_LAPACKE_ztfttp ztfttp;

    ztfttp_obj = new tfttp_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
							lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].n);
						   

    idx = Circular_Increment_Index(idx);

    ztfttp_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ztfttp_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ztfttp_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ztfttp_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*ztfttp library call */
    ztfttp = (Fptr_NL_LAPACKE_ztfttp)dlsym(ztfttp_obj->hModule, "LAPACKE_ztfttp");
    ASSERT_TRUE(ztfttp != NULL) << "failed to get the Netlib LAPACKE_ztfttp symbol";
    
/*Compute ztfttp's  o/p */
    ztfttp_obj->inforef = ztfttp( ztfttp_obj->matrix_layout, ztfttp_obj->transr, ztfttp_obj->uplo,  ztfttp_obj->n,
								 (const lapack_complex_double*)ztfttp_obj->Aref, ztfttp_obj->apref);

    /* Compute libflame's Lapacke o/p  */
    ztfttp_obj->info = LAPACKE_ztfttp( ztfttp_obj->matrix_layout, ztfttp_obj->transr, ztfttp_obj->uplo, ztfttp_obj->n, 
										(const lapack_complex_double*)ztfttp_obj->A,  ztfttp_obj->ap);

    if( ztfttp_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ztfttp is wrong\n", ztfttp_obj->info );
    }
    if( ztfttp_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ztfttp is wrong\n", 
        ztfttp_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    ztfttp_obj->diff =  computeDiff_z( ztfttp_obj->bufsize_a, 
                ztfttp_obj->A, ztfttp_obj->Aref );

}

TEST_F(ztfttp_test, ztfttp1) {
    EXPECT_NEAR(0.0, ztfttp_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ztfttp_test, ztfttp2) {
    EXPECT_NEAR(0.0, ztfttp_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ztfttp_test, ztfttp3) {
    EXPECT_NEAR(0.0, ztfttp_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ztfttp_test, ztfttp4) {
    EXPECT_NEAR(0.0, ztfttp_obj->diff, LAPACKE_GTEST_THRESHOLD);
}
