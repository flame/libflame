#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"


#define pbstf_free() \
if (bb!=NULL)    free(bb); \
if (bbref!=NULL) free(bbref);\
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
/* Begin float_common_parameters  class definition */
class pbstf_float_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int kb;
	char uplo;
	float* bb;
	lapack_int ldbb;
	/*Output Parameter*/
	float *bbref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      pbstf_float_parameters (int matrix_layout, char uplo, lapack_int n, lapack_int kb, lapack_int ldbb);
      ~pbstf_float_parameters ();

};

/* Constructor definition  float_common_parameters */
pbstf_float_parameters:: pbstf_float_parameters (int matrix_layout_i, char uplo_i, lapack_int n_i, lapack_int kb_i, lapack_int ldbb_i)
{
	matrix_layout = matrix_layout_i;
	n = n_i;
	ldbb = ldbb_i;
	uplo = uplo_i;
	kb = kb_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n pbstf float: n: %d lda: %d \n",n, lda);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
		bufsize = ldbb*n;
	else if (matrix_layout == LAPACK_ROW_MAJOR)
		bufsize = ldbb*(kb+1);
	else 
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&bb, &bbref, bufsize);
	if ((bb==NULL) || (bbref==NULL)){
		EXPECT_FALSE( true) << "pbstf_float_parameters object: malloc error.";
		pbstf_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand( bb, bbref, bufsize);
} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
pbstf_float_parameters :: ~pbstf_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" pbstf_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   pbstf_free();

}
/*  Test fixture class definition */
class spbstf_test  : public  ::testing::Test {
public:
   pbstf_float_parameters  *spbstf_obj;
   void SetUp();
   void TearDown () { delete spbstf_obj; }
};

void spbstf_test::SetUp(){

    /* LAPACKE spbstf prototype */
    typedef int (*Fptr_NL_LAPACKE_spbstf) (int matrix_layout, char uplo,lapack_int n,\
											lapack_int kb, float *bb, lapack_int ldbb);

    Fptr_NL_LAPACKE_spbstf spbstf;

    spbstf_obj = new pbstf_float_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].uplo,
                           eig_paramslist[idx].n,
						   eig_paramslist[idx].kb,
                           eig_paramslist[idx].ldb);

    idx = Circular_Increment_Index(idx);
	

    spbstf_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    spbstf_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(spbstf_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(spbstf_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    spbstf = (Fptr_NL_LAPACKE_spbstf)dlsym(spbstf_obj->hModule, "LAPACKE_spbstf");
    ASSERT_TRUE(spbstf != NULL) << "failed to get the Netlib LAPACKE_spbstf symbol";
    

    spbstf_obj->inforef = spbstf( spbstf_obj->matrix_layout, spbstf_obj->uplo,
								spbstf_obj->n, spbstf_obj->kb, spbstf_obj->bbref,
								spbstf_obj->ldbb);

    /* Compute libflame's Lapacke o/p  */
    spbstf_obj->info = LAPACKE_spbstf( spbstf_obj->matrix_layout, spbstf_obj->uplo, 
										spbstf_obj->n,spbstf_obj->kb,spbstf_obj->bb, 
										spbstf_obj->ldbb);

    if( spbstf_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_spbstf is wrong\n", spbstf_obj->info );
    }
    if( spbstf_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_spbstf is wrong\n", 
        spbstf_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    spbstf_obj->diff =  computeDiff_s( spbstf_obj->bufsize, 
                spbstf_obj->bb, spbstf_obj->bbref );

}

TEST_F(spbstf_test, spbstf1) {
    EXPECT_NEAR(0.0, spbstf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(spbstf_test, spbstf2) {
    EXPECT_NEAR(0.0, spbstf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(spbstf_test, spbstf3) {
    EXPECT_NEAR(0.0, spbstf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(spbstf_test, spbstf4) {
    EXPECT_NEAR(0.0, spbstf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin double_common_parameters  class definition */
class pbstf_double_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int kb;
	char uplo;
	double* bb;
	lapack_int ldbb;
	/*Output Parameter*/
	double *bbref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      pbstf_double_parameters (int matrix_layout, char uplo, lapack_int n, lapack_int kb, lapack_int ldbb);
      ~pbstf_double_parameters ();

};

/* Constructor definition  double_common_parameters */
pbstf_double_parameters:: pbstf_double_parameters (int matrix_layout_i, char uplo_i, lapack_int n_i, lapack_int kb_i, lapack_int ldbb_i)
{
	matrix_layout = matrix_layout_i;
	n = n_i;
	ldbb = ldbb_i;
	uplo = uplo_i;
	kb = kb_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n pbstf double: n: %d lda: %d \n",n, lda);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
		bufsize = ldbb*n;
	else if (matrix_layout == LAPACK_ROW_MAJOR)
		bufsize = ldbb*(kb+1);
	else 
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&bb, &bbref, bufsize);
	if ((bb==NULL) || (bbref==NULL)){
		EXPECT_FALSE( true) << "pbstf_double_parameters object: malloc error.";
		pbstf_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand( bb, bbref, bufsize);
} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
pbstf_double_parameters :: ~pbstf_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" pbstf_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   pbstf_free();

}
/*  Test fixture class definition */
class dpbstf_test  : public  ::testing::Test {
public:
   pbstf_double_parameters  *dpbstf_obj;
   void SetUp();
   void TearDown () { delete dpbstf_obj; }
};

void dpbstf_test::SetUp(){

    /* LAPACKE dpbstf prototype */
    typedef int (*Fptr_NL_LAPACKE_dpbstf) (int matrix_layout, char uplo,lapack_int n,\
											lapack_int kb, double *bb, lapack_int ldbb);

    Fptr_NL_LAPACKE_dpbstf dpbstf;

    dpbstf_obj = new pbstf_double_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].uplo,
                           eig_paramslist[idx].n,
						   eig_paramslist[idx].kb,
                           eig_paramslist[idx].ldb);

    idx = Circular_Increment_Index(idx);
	

    dpbstf_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dpbstf_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dpbstf_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dpbstf_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    dpbstf = (Fptr_NL_LAPACKE_dpbstf)dlsym(dpbstf_obj->hModule, "LAPACKE_dpbstf");
    ASSERT_TRUE(dpbstf != NULL) << "failed to get the Netlib LAPACKE_dpbstf symbol";
    

    dpbstf_obj->inforef = dpbstf( dpbstf_obj->matrix_layout, dpbstf_obj->uplo,
								dpbstf_obj->n, dpbstf_obj->kb, dpbstf_obj->bbref,
								dpbstf_obj->ldbb);

    /* Compute libflame's Lapacke o/p  */
    dpbstf_obj->info = LAPACKE_dpbstf( dpbstf_obj->matrix_layout, dpbstf_obj->uplo, 
										dpbstf_obj->n,dpbstf_obj->kb,dpbstf_obj->bb, 
										dpbstf_obj->ldbb);

    if( dpbstf_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dpbstf is wrong\n", dpbstf_obj->info );
    }
    if( dpbstf_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dpbstf is wrong\n", 
        dpbstf_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dpbstf_obj->diff =  computeDiff_d( dpbstf_obj->bufsize, 
                dpbstf_obj->bb, dpbstf_obj->bbref );

}

TEST_F(dpbstf_test, dpbstf1) {
    EXPECT_NEAR(0.0, dpbstf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dpbstf_test, dpbstf2) {
    EXPECT_NEAR(0.0, dpbstf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dpbstf_test, dpbstf3) {
    EXPECT_NEAR(0.0, dpbstf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dpbstf_test, dpbstf4) {
    EXPECT_NEAR(0.0, dpbstf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}
/* Begin scomplex_common_parameters  class definition */
class pbstf_scomplex_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int kb;
	char uplo;
	lapack_complex_float* bb;
	lapack_int ldbb;
	/*Output Parameter*/
	lapack_complex_float* bbref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      pbstf_scomplex_parameters (int matrix_layout, char uplo, lapack_int n, lapack_int kb, lapack_int ldbb);
      ~pbstf_scomplex_parameters ();

};

/* Constructor definition  scomplex_common_parameters */
pbstf_scomplex_parameters:: pbstf_scomplex_parameters (int matrix_layout_i, char uplo_i, lapack_int n_i, lapack_int kb_i, lapack_int ldbb_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;
	ldbb = ldbb_i;
	uplo = uplo_i;
	kb = kb_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n pbstf scomplex: n: %d lda: %d \n",n, lda);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
		bufsize = ldbb*n;
	else if (matrix_layout == LAPACK_ROW_MAJOR)
		bufsize = ldbb*(kb+1);
	else 
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&bb, &bbref, bufsize);
	if ((bb==NULL) || (bbref==NULL)){
		EXPECT_FALSE( true) << "pbstf_double_parameters object: malloc error.";
		pbstf_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( bb, bbref, bufsize);
} /* end of Constructor  */

/* Destructor definition  'scomplex_common_parameters' class  */
pbstf_scomplex_parameters :: ~pbstf_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" pbstf_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   pbstf_free();

}
/*  Test fixture class definition */
class cpbstf_test  : public  ::testing::Test {
public:
   pbstf_scomplex_parameters  *cpbstf_obj;
   void SetUp();
   void TearDown () { delete cpbstf_obj; }
};

void cpbstf_test::SetUp(){

    /* LAPACKE cpbstf prototype */
    typedef int (*Fptr_NL_LAPACKE_cpbstf) (int matrix_layout, char uplo,lapack_int n,\
											lapack_int kb, lapack_complex_float *bb, lapack_int ldbb);

    Fptr_NL_LAPACKE_cpbstf cpbstf;

    cpbstf_obj = new pbstf_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].uplo,
                           eig_paramslist[idx].n,
						   eig_paramslist[idx].kb,
                           eig_paramslist[idx].ldb);

    idx = Circular_Increment_Index(idx);
	

    cpbstf_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cpbstf_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cpbstf_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cpbstf_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    cpbstf = (Fptr_NL_LAPACKE_cpbstf)dlsym(cpbstf_obj->hModule, "LAPACKE_cpbstf");
    ASSERT_TRUE(cpbstf != NULL) << "failed to get the Netlib LAPACKE_cpbstf symbol";
    

    cpbstf_obj->inforef = cpbstf( cpbstf_obj->matrix_layout, cpbstf_obj->uplo,
								cpbstf_obj->n, cpbstf_obj->kb, cpbstf_obj->bbref,
								cpbstf_obj->ldbb);

    /* Compute libflame's Lapacke o/p  */
    cpbstf_obj->info = LAPACKE_cpbstf( cpbstf_obj->matrix_layout, cpbstf_obj->uplo, 
										cpbstf_obj->n,cpbstf_obj->kb,cpbstf_obj->bb, 
										cpbstf_obj->ldbb);

    if( cpbstf_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cpbstf is wrong\n", cpbstf_obj->info );
    }
    if( cpbstf_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cpbstf is wrong\n", 
        cpbstf_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cpbstf_obj->diff =  computeDiff_c( cpbstf_obj->bufsize, 
                cpbstf_obj->bb, cpbstf_obj->bbref );

}

TEST_F(cpbstf_test, cpbstf1) {
    EXPECT_NEAR(0.0, cpbstf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cpbstf_test, cpbstf2) {
    EXPECT_NEAR(0.0, cpbstf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cpbstf_test, cpbstf3) {
    EXPECT_NEAR(0.0, cpbstf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cpbstf_test, cpbstf4) {
    EXPECT_NEAR(0.0, cpbstf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin scomplex_common_parameters  class definition */
class pbstf_dcomplex_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int kb;
	char uplo;
	lapack_complex_double* bb;
	lapack_int ldbb;
	/*Output Parameter*/
	lapack_complex_double* bbref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      pbstf_dcomplex_parameters (int matrix_layout, char uplo, lapack_int n, lapack_int kb, lapack_int ldbb);
      ~pbstf_dcomplex_parameters ();

};

/* Constructor definition  dcomplex_common_parameters */
pbstf_dcomplex_parameters:: pbstf_dcomplex_parameters (int matrix_layout_i, char uplo_i, lapack_int n_i, lapack_int kb_i, lapack_int ldbb_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;
	ldbb = ldbb_i;
	uplo = uplo_i;
	kb = kb_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n pbstf scomplex: n: %d lda: %d \n",n, lda);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
		bufsize = ldbb*n;
	else if (matrix_layout == LAPACK_ROW_MAJOR)
		bufsize = ldbb*(kb+1);
	else 
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&bb, &bbref, bufsize);
	if ((bb==NULL) || (bbref==NULL)){
		EXPECT_FALSE( true) << "pbstf_double_parameters object: malloc error.";
		pbstf_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( bb, bbref, bufsize);
} /* end of Constructor  */

/* Destructor definition  'dcomplex_common_parameters' class  */
pbstf_dcomplex_parameters :: ~pbstf_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" pbstf_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   pbstf_free();

}
/*  Test fixture class definition */
class zpbstf_test  : public  ::testing::Test {
public:
   pbstf_dcomplex_parameters  *zpbstf_obj;
   void SetUp();
   void TearDown () { delete zpbstf_obj; }
};

void zpbstf_test::SetUp(){

    /* LAPACKE zpbstf prototype */
    typedef int (*Fptr_NL_LAPACKE_zpbstf) (int matrix_layout, char uplo,lapack_int n,\
											lapack_int kb, lapack_complex_double *bb, lapack_int ldbb);

    Fptr_NL_LAPACKE_zpbstf zpbstf;

    zpbstf_obj = new pbstf_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].uplo,
                           eig_paramslist[idx].n,
						   eig_paramslist[idx].kb,
                           eig_paramslist[idx].ldb);

    idx = Circular_Increment_Index(idx);
	

    zpbstf_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zpbstf_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zpbstf_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zpbstf_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zpbstf = (Fptr_NL_LAPACKE_zpbstf)dlsym(zpbstf_obj->hModule, "LAPACKE_zpbstf");
    ASSERT_TRUE(zpbstf != NULL) << "failed to get the Netlib LAPACKE_zpbstf symbol";
    

    zpbstf_obj->inforef = zpbstf( zpbstf_obj->matrix_layout, zpbstf_obj->uplo,
								zpbstf_obj->n, zpbstf_obj->kb, zpbstf_obj->bbref,
								zpbstf_obj->ldbb);

    /* Compute libflame's Lapacke o/p  */
    zpbstf_obj->info = LAPACKE_zpbstf( zpbstf_obj->matrix_layout, zpbstf_obj->uplo, 
										zpbstf_obj->n,zpbstf_obj->kb,zpbstf_obj->bb, 
										zpbstf_obj->ldbb);

    if( zpbstf_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zpbstf is wrong\n", zpbstf_obj->info );
    }
    if( zpbstf_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zpbstf is wrong\n", 
        zpbstf_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zpbstf_obj->diff =  computeDiff_z( zpbstf_obj->bufsize, 
                zpbstf_obj->bb, zpbstf_obj->bbref );

}

TEST_F(zpbstf_test, zpbstf1) {
    EXPECT_NEAR(0.0, zpbstf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zpbstf_test, zpbstf2) {
    EXPECT_NEAR(0.0, zpbstf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zpbstf_test, zpbstf3) {
    EXPECT_NEAR(0.0, zpbstf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zpbstf_test, zpbstf4) {
    EXPECT_NEAR(0.0, zpbstf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

