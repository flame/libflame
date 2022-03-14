#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"

#define syswapr_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (ipiv!=NULL)    free(ipiv); \
if (ipivref!=NULL)    free(ipivref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule); \
if(bModule != NULL) dlclose(bModule); \
if(lModule != NULL) dlclose(lModule)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
/* Begin float_common_parameters  class definition */
class syswapr_float_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	void *bModule, *lModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;	
	float* A;
	char uplo;
	lapack_int i1, i2;
	lapack_int lda, *ipiv;
	float* Aref;
	lapack_int *ipivref;
	/*Return Values*/
	lapack_int info, inforef;
	lapack_int info_sytrf, inforef_sytrf;

   public:
      syswapr_float_parameters (int matrix_layout, char uplo, lapack_int n, lapack_int lda, lapack_int i1, lapack_int i2);
      ~syswapr_float_parameters ();

};

/* Constructor definition  float_common_parameters */
syswapr_float_parameters:: syswapr_float_parameters (int matrix_layout_i, char uplo_i, lapack_int n_i, lapack_int lda_i, lapack_int i1_i, lapack_int i2_i)
{

	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	lda = lda_i;
	i1 = i1_i;
	i2 = i2_i;
	#if LAPACKE_TEST_VERBOSE
		printf(" \n syswapr float: matrix_layout = %d, uplo:%c, n: %d \n", matrix_layout, uplo, n);
	#endif

	bufsize = lda*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_int_buffer_pair(&ipiv, &ipivref, n);
	if ((A==NULL) || (Aref==NULL) || \
		(ipiv == NULL) || (ipivref == NULL)) {	
		EXPECT_FALSE( true) << "syswapr_float_parameters object: malloc error.";
		syswapr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand_custom_matrix( A, Aref, bufsize, n, uplo);
	

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
syswapr_float_parameters :: ~syswapr_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" syswapr_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   syswapr_free();

}
/*  Test fixture class definition */
class ssyswapr_test  : public  ::testing::Test {
public:
   syswapr_float_parameters  *ssyswapr_obj;
   void SetUp();
   void TearDown () { delete ssyswapr_obj;}
};

void ssyswapr_test::SetUp(){

	/* LAPACKE ssytrf prototype */
    typedef int (*Fptr_NL_LAPACKE_ssytrf) (int matrix_layout , char uplo , lapack_int n , float * a , lapack_int lda , lapack_int * ipiv);
	 Fptr_NL_LAPACKE_ssytrf ssytrf;
	 /* LAPACKE ssyswapr prototype */
    typedef int (*Fptr_NL_LAPACKE_ssyswapr) (int matrix_layout, char uplo, lapack_int n, float* a, lapack_int lda,\
                             lapack_int i1, lapack_int i2 );

    Fptr_NL_LAPACKE_ssyswapr ssyswapr;

    ssyswapr_obj = new syswapr_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
						   lin_solver_paramslist[idx].ldb,
						   lin_solver_paramslist[idx].kl,
						   lin_solver_paramslist[idx].ku);

    idx = Circular_Increment_Index(idx);

    ssyswapr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ssyswapr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ssyswapr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ssyswapr_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*ssyswapr library call */
    ssyswapr = (Fptr_NL_LAPACKE_ssyswapr)dlsym(ssyswapr_obj->hModule, "LAPACKE_ssyswapr");
    ASSERT_TRUE(ssyswapr != NULL) << "failed to get the Netlib LAPACKE_ssyswapr symbol";

	/*ssytrf library call*/
	ssyswapr_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ssyswapr_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ssyswapr_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ssyswapr_obj->lModule != NULL) << "Netlib lapacke handle NULL";

    ssytrf = (Fptr_NL_LAPACKE_ssytrf)dlsym(ssyswapr_obj->lModule, "LAPACKE_ssytrf");
    ASSERT_TRUE(ssytrf != NULL) << "failed to get the Netlib LAPACKE_ssytrf symbol";    
    

    ssyswapr_obj->inforef_sytrf = ssytrf( ssyswapr_obj->matrix_layout, ssyswapr_obj->uplo,
								ssyswapr_obj->n, ssyswapr_obj->Aref,  ssyswapr_obj->lda, ssyswapr_obj->ipivref);

    /* Compute libflame's Lapacke o/p  */
    ssyswapr_obj->info_sytrf = LAPACKE_ssytrf( ssyswapr_obj->matrix_layout, ssyswapr_obj->uplo,
										ssyswapr_obj->n,ssyswapr_obj->A, ssyswapr_obj->lda, ssyswapr_obj->ipiv);

    if( ssyswapr_obj->info_sytrf < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ssytrf is wrong\n", ssyswapr_obj->info_sytrf );
    }
    if( ssyswapr_obj->inforef_sytrf < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ssytrf is wrong\n", 
        ssyswapr_obj->inforef_sytrf );
    }  
/*Compute ssyswapr's  o/p */
    ssyswapr_obj->inforef = ssyswapr( ssyswapr_obj->matrix_layout, ssyswapr_obj->uplo,\
							ssyswapr_obj->n, ssyswapr_obj->Aref, ssyswapr_obj->lda, ssyswapr_obj->i1, ssyswapr_obj->i2);

    /* Compute libflame's Lapacke o/p  */
	
    ssyswapr_obj->info = LAPACKE_ssyswapr( ssyswapr_obj->matrix_layout, ssyswapr_obj->uplo,\
											ssyswapr_obj->n,ssyswapr_obj->A, ssyswapr_obj->lda, ssyswapr_obj->i1, ssyswapr_obj->i2);

    if( ssyswapr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ssyswapr is wrong\n", ssyswapr_obj->info );
    }
    if( ssyswapr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ssyswapr is wrong\n", 
        ssyswapr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    ssyswapr_obj->diff =  computeDiff_s( ssyswapr_obj->bufsize, 
                ssyswapr_obj->A, ssyswapr_obj->Aref );

}

TEST_F(ssyswapr_test, ssyswapr1) {
    EXPECT_NEAR(0.0, ssyswapr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssyswapr_test, ssyswapr2) {
    EXPECT_NEAR(0.0, ssyswapr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssyswapr_test, ssyswapr3) {
    EXPECT_NEAR(0.0, ssyswapr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssyswapr_test, ssyswapr4) {
    EXPECT_NEAR(0.0, ssyswapr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin double_common_parameters  class definition */
class syswapr_double_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	void *bModule, *lModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;	
	double* A;
	char uplo;
	lapack_int i1, i2;
	lapack_int lda, *ipiv;
	double* Aref;
	lapack_int *ipivref;
	/*Return Values*/
	lapack_int info, inforef;
	lapack_int info_sytrf, inforef_sytrf;

   public:
      syswapr_double_parameters (int matrix_layout, char uplo, lapack_int n, lapack_int lda, lapack_int i1, lapack_int i2);
      ~syswapr_double_parameters ();

};

/* Constructor definition  double_common_parameters */
syswapr_double_parameters:: syswapr_double_parameters (int matrix_layout_i, char uplo_i, lapack_int n_i, lapack_int lda_i, lapack_int i1_i, lapack_int i2_i)
{

	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	lda = lda_i;
	i1 = i1_i;
	i2 = i2_i;
	#if LAPACKE_TEST_VERBOSE
		printf(" \n syswapr double: matrix_layout = %d, uplo:%c, n: %d \n", matrix_layout, uplo, n);
	#endif

	bufsize = lda*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_int_buffer_pair(&ipiv, &ipivref, n);
	if ((A==NULL) || (Aref==NULL) || \
		(ipiv == NULL) || (ipivref == NULL)) {	
		EXPECT_FALSE( true) << "syswapr_double_parameters object: malloc error.";
		syswapr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand_custom_matrix( A, Aref, bufsize, n, uplo);
	

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
syswapr_double_parameters :: ~syswapr_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" syswapr_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   syswapr_free();

}
/*  Test fixture class definition */
class dsyswapr_test  : public  ::testing::Test {
public:
   syswapr_double_parameters  *dsyswapr_obj;
   void SetUp();
   void TearDown () { delete dsyswapr_obj;}
};

void dsyswapr_test::SetUp(){

	/* LAPACKE dsytrf prototype */
    typedef int (*Fptr_NL_LAPACKE_dsytrf) (int matrix_layout , char uplo , lapack_int n , double * a , lapack_int lda , lapack_int * ipiv);
	 Fptr_NL_LAPACKE_dsytrf dsytrf;
	 /* LAPACKE dsyswapr prototype */
    typedef int (*Fptr_NL_LAPACKE_dsyswapr) (int matrix_layout, char uplo, lapack_int n, double* a, lapack_int lda,\
                             lapack_int i1, lapack_int i2 );

    Fptr_NL_LAPACKE_dsyswapr dsyswapr;

    dsyswapr_obj = new syswapr_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
						   lin_solver_paramslist[idx].ldb,
						   lin_solver_paramslist[idx].kl,
						   lin_solver_paramslist[idx].ku);

    idx = Circular_Increment_Index(idx);

    dsyswapr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dsyswapr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dsyswapr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dsyswapr_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*dsyswapr library call */
    dsyswapr = (Fptr_NL_LAPACKE_dsyswapr)dlsym(dsyswapr_obj->hModule, "LAPACKE_dsyswapr");
    ASSERT_TRUE(dsyswapr != NULL) << "failed to get the Netlib LAPACKE_dsyswapr symbol";

	/*dsytrf library call*/
	dsyswapr_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dsyswapr_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dsyswapr_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dsyswapr_obj->lModule != NULL) << "Netlib lapacke handle NULL";

    dsytrf = (Fptr_NL_LAPACKE_dsytrf)dlsym(dsyswapr_obj->lModule, "LAPACKE_dsytrf");
    ASSERT_TRUE(dsytrf != NULL) << "failed to get the Netlib LAPACKE_dsytrf symbol";    
    

    dsyswapr_obj->inforef_sytrf = dsytrf( dsyswapr_obj->matrix_layout, dsyswapr_obj->uplo,
								dsyswapr_obj->n, dsyswapr_obj->Aref,  dsyswapr_obj->lda, dsyswapr_obj->ipivref);

    /* Compute libflame's Lapacke o/p  */
    dsyswapr_obj->info_sytrf = LAPACKE_dsytrf( dsyswapr_obj->matrix_layout, dsyswapr_obj->uplo,
										dsyswapr_obj->n,dsyswapr_obj->A, dsyswapr_obj->lda, dsyswapr_obj->ipiv);

    if( dsyswapr_obj->info_sytrf < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dsytrf is wrong\n", dsyswapr_obj->info_sytrf );
    }
    if( dsyswapr_obj->inforef_sytrf < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dsytrf is wrong\n", 
        dsyswapr_obj->inforef_sytrf );
    }  
/*Compute dsyswapr's  o/p */
    dsyswapr_obj->inforef = dsyswapr( dsyswapr_obj->matrix_layout, dsyswapr_obj->uplo,\
							dsyswapr_obj->n, dsyswapr_obj->Aref, dsyswapr_obj->lda, dsyswapr_obj->i1, dsyswapr_obj->i2);

    /* Compute libflame's Lapacke o/p  */
	
    dsyswapr_obj->info = LAPACKE_dsyswapr( dsyswapr_obj->matrix_layout, dsyswapr_obj->uplo,\
											dsyswapr_obj->n,dsyswapr_obj->A, dsyswapr_obj->lda, dsyswapr_obj->i1, dsyswapr_obj->i2);

    if( dsyswapr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dsyswapr is wrong\n", dsyswapr_obj->info );
    }
    if( dsyswapr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dsyswapr is wrong\n", 
        dsyswapr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dsyswapr_obj->diff =  computeDiff_d( dsyswapr_obj->bufsize, 
                dsyswapr_obj->A, dsyswapr_obj->Aref );

}

TEST_F(dsyswapr_test, dsyswapr1) {
    EXPECT_NEAR(0.0, dsyswapr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dsyswapr_test, dsyswapr2) {
    EXPECT_NEAR(0.0, dsyswapr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dsyswapr_test, dsyswapr3) {
    EXPECT_NEAR(0.0, dsyswapr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dsyswapr_test, dsyswapr4) {
    EXPECT_NEAR(0.0, dsyswapr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}





/* Begin scomplex_common_parameters  class definition */
class syswapr_scomplex_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	void *bModule, *lModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;	
	lapack_complex_float* A;
	char uplo;
	lapack_int i1, i2;
	lapack_int lda, *ipiv;
	lapack_complex_float* Aref;
	lapack_int *ipivref;
	/*Return Values*/
	lapack_int info, inforef;
	lapack_int info_sytrf, inforef_sytrf;

   public:
      syswapr_scomplex_parameters (int matrix_layout, char uplo, lapack_int n, lapack_int lda, lapack_int i1, lapack_int i2);
      ~syswapr_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
syswapr_scomplex_parameters:: syswapr_scomplex_parameters (int matrix_layout_i, char uplo_i, lapack_int n_i, lapack_int lda_i, lapack_int i1_i, lapack_int i2_i)
{

	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	lda = lda_i;
	i1 = i1_i;
	i2 = i2_i;
	#if LAPACKE_TEST_VERBOSE
		printf(" \n syswapr scomplex: matrix_layout = %d, uplo:%c, n: %d \n", matrix_layout, uplo, n);
	#endif

	bufsize = lda*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_int_buffer_pair(&ipiv, &ipivref, n);
	if ((A==NULL) || (Aref==NULL) || \
		(ipiv == NULL) || (ipivref == NULL)) {	
		EXPECT_FALSE( true) << "syswapr_float_parameters object: malloc error.";
		syswapr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix( A, Aref, bufsize, n, uplo);
	

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
syswapr_scomplex_parameters :: ~syswapr_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" syswapr_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   syswapr_free();

}
/*  Test fixture class definition */
class csyswapr_test  : public  ::testing::Test {
public:
   syswapr_scomplex_parameters  *csyswapr_obj;
   void SetUp();
   void TearDown () { delete csyswapr_obj;}
};

void csyswapr_test::SetUp(){

	/* LAPACKE csytrf prototype */
    typedef int (*Fptr_NL_LAPACKE_csytrf) (int matrix_layout , char uplo , lapack_int n , lapack_complex_float * a , lapack_int lda , lapack_int * ipiv);
	 Fptr_NL_LAPACKE_csytrf csytrf;
	 /* LAPACKE csyswapr prototype */
    typedef int (*Fptr_NL_LAPACKE_csyswapr) (int matrix_layout, char uplo, lapack_int n, lapack_complex_float* a, lapack_int lda,\
                             lapack_int i1, lapack_int i2 );

    Fptr_NL_LAPACKE_csyswapr csyswapr;

    csyswapr_obj = new syswapr_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
						   lin_solver_paramslist[idx].ldb,
						   lin_solver_paramslist[idx].kl,
						   lin_solver_paramslist[idx].ku);

    idx = Circular_Increment_Index(idx);

    csyswapr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    csyswapr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(csyswapr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(csyswapr_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*csyswapr library call */
    csyswapr = (Fptr_NL_LAPACKE_csyswapr)dlsym(csyswapr_obj->hModule, "LAPACKE_csyswapr");
    ASSERT_TRUE(csyswapr != NULL) << "failed to get the Netlib LAPACKE_csyswapr symbol";

	/*csytrf library call*/
	csyswapr_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    csyswapr_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(csyswapr_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(csyswapr_obj->lModule != NULL) << "Netlib lapacke handle NULL";

    csytrf = (Fptr_NL_LAPACKE_csytrf)dlsym(csyswapr_obj->lModule, "LAPACKE_csytrf");
    ASSERT_TRUE(csytrf != NULL) << "failed to get the Netlib LAPACKE_csytrf symbol";    
    

    csyswapr_obj->inforef_sytrf = csytrf( csyswapr_obj->matrix_layout, csyswapr_obj->uplo,
								csyswapr_obj->n, csyswapr_obj->Aref,  csyswapr_obj->lda, csyswapr_obj->ipivref);

    /* Compute libflame's Lapacke o/p  */
    csyswapr_obj->info_sytrf = LAPACKE_csytrf( csyswapr_obj->matrix_layout, csyswapr_obj->uplo,
										csyswapr_obj->n,csyswapr_obj->A, csyswapr_obj->lda, csyswapr_obj->ipiv);

    if( csyswapr_obj->info_sytrf < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_csytrf is wrong\n", csyswapr_obj->info_sytrf );
    }
    if( csyswapr_obj->inforef_sytrf < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_csytrf is wrong\n", 
        csyswapr_obj->inforef_sytrf );
    }  
/*Compute csyswapr's  o/p */
    csyswapr_obj->inforef = csyswapr( csyswapr_obj->matrix_layout, csyswapr_obj->uplo,\
							csyswapr_obj->n, csyswapr_obj->Aref, csyswapr_obj->lda, csyswapr_obj->i1, csyswapr_obj->i2);

    /* Compute libflame's Lapacke o/p  */
	
    csyswapr_obj->info = LAPACKE_csyswapr( csyswapr_obj->matrix_layout, csyswapr_obj->uplo,\
											csyswapr_obj->n,csyswapr_obj->A, csyswapr_obj->lda, csyswapr_obj->i1, csyswapr_obj->i2);

    if( csyswapr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_csyswapr is wrong\n", csyswapr_obj->info );
    }
    if( csyswapr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_csyswapr is wrong\n", 
        csyswapr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    csyswapr_obj->diff =  computeDiff_c( csyswapr_obj->bufsize, 
                csyswapr_obj->A, csyswapr_obj->Aref );

}

TEST_F(csyswapr_test, csyswapr1) {
    EXPECT_NEAR(0.0, csyswapr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(csyswapr_test, csyswapr2) {
    EXPECT_NEAR(0.0, csyswapr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(csyswapr_test, csyswapr3) {
    EXPECT_NEAR(0.0, csyswapr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(csyswapr_test, csyswapr4) {
    EXPECT_NEAR(0.0, csyswapr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class syswapr_dcomplex_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	void *bModule, *lModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;	
	lapack_complex_double* A;
	char uplo;
	lapack_int i1, i2;
	lapack_int lda, *ipiv;
	lapack_complex_double* Aref;
	lapack_int *ipivref;
	/*Return Values*/
	lapack_int info, inforef;
	lapack_int info_sytrf, inforef_sytrf;

   public:
      syswapr_dcomplex_parameters (int matrix_layout, char uplo, lapack_int n, lapack_int lda, lapack_int i1, lapack_int i2);
      ~syswapr_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
syswapr_dcomplex_parameters:: syswapr_dcomplex_parameters (int matrix_layout_i, char uplo_i, lapack_int n_i, lapack_int lda_i, lapack_int i1_i, lapack_int i2_i)
{

	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	lda = lda_i;
	i1 = i1_i;
	i2 = i2_i;
	#if LAPACKE_TEST_VERBOSE
		printf(" \n syswapr dcomplex: matrix_layout = %d, uplo:%c, n: %d \n", matrix_layout, uplo, n);
	#endif

	bufsize = lda*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_int_buffer_pair(&ipiv, &ipivref, n);
	if ((A==NULL) || (Aref==NULL) || \
		(ipiv == NULL) || (ipivref == NULL)) {	
		EXPECT_FALSE( true) << "syswapr_double_parameters object: malloc error.";
		syswapr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix( A, Aref, bufsize, n, uplo);
	

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
syswapr_dcomplex_parameters :: ~syswapr_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" syswapr_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   syswapr_free();

}
/*  Test fixture class definition */
class zsyswapr_test  : public  ::testing::Test {
public:
   syswapr_dcomplex_parameters  *zsyswapr_obj;
   void SetUp();
   void TearDown () { delete zsyswapr_obj;}
};

void zsyswapr_test::SetUp(){

	/* LAPACKE zsytrf prototype */
    typedef int (*Fptr_NL_LAPACKE_zsytrf) (int matrix_layout , char uplo , lapack_int n , lapack_complex_double * a , lapack_int lda , lapack_int * ipiv);
	 Fptr_NL_LAPACKE_zsytrf zsytrf;
	 /* LAPACKE zsyswapr prototype */
    typedef int (*Fptr_NL_LAPACKE_zsyswapr) (int matrix_layout, char uplo, lapack_int n, lapack_complex_double* a, lapack_int lda,\
                             lapack_int i1, lapack_int i2 );

    Fptr_NL_LAPACKE_zsyswapr zsyswapr;

    zsyswapr_obj = new syswapr_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
						   lin_solver_paramslist[idx].ldb,
						   lin_solver_paramslist[idx].kl,
						   lin_solver_paramslist[idx].ku);

    idx = Circular_Increment_Index(idx);

    zsyswapr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zsyswapr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zsyswapr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zsyswapr_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*zsyswapr library call */
    zsyswapr = (Fptr_NL_LAPACKE_zsyswapr)dlsym(zsyswapr_obj->hModule, "LAPACKE_zsyswapr");
    ASSERT_TRUE(zsyswapr != NULL) << "failed to get the Netlib LAPACKE_zsyswapr symbol";

	/*zsytrf library call*/
	zsyswapr_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zsyswapr_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zsyswapr_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zsyswapr_obj->lModule != NULL) << "Netlib lapacke handle NULL";

    zsytrf = (Fptr_NL_LAPACKE_zsytrf)dlsym(zsyswapr_obj->lModule, "LAPACKE_zsytrf");
    ASSERT_TRUE(zsytrf != NULL) << "failed to get the Netlib LAPACKE_zsytrf symbol";    
    

    zsyswapr_obj->inforef_sytrf = zsytrf( zsyswapr_obj->matrix_layout, zsyswapr_obj->uplo,
								zsyswapr_obj->n, zsyswapr_obj->Aref,  zsyswapr_obj->lda, zsyswapr_obj->ipivref);

    /* Compute libflame's Lapacke o/p  */
    zsyswapr_obj->info_sytrf = LAPACKE_zsytrf( zsyswapr_obj->matrix_layout, zsyswapr_obj->uplo,
										zsyswapr_obj->n,zsyswapr_obj->A, zsyswapr_obj->lda, zsyswapr_obj->ipiv);

    if( zsyswapr_obj->info_sytrf < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zsytrf is wrong\n", zsyswapr_obj->info_sytrf );
    }
    if( zsyswapr_obj->inforef_sytrf < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zsytrf is wrong\n", 
        zsyswapr_obj->inforef_sytrf );
    }  
/*Compute zsyswapr's  o/p */
    zsyswapr_obj->inforef = zsyswapr( zsyswapr_obj->matrix_layout, zsyswapr_obj->uplo,\
							zsyswapr_obj->n, zsyswapr_obj->Aref, zsyswapr_obj->lda, zsyswapr_obj->i1, zsyswapr_obj->i2);

    /* Compute libflame's Lapacke o/p  */
	
    zsyswapr_obj->info = LAPACKE_zsyswapr( zsyswapr_obj->matrix_layout, zsyswapr_obj->uplo,\
											zsyswapr_obj->n,zsyswapr_obj->A, zsyswapr_obj->lda, zsyswapr_obj->i1, zsyswapr_obj->i2);

    if( zsyswapr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zsyswapr is wrong\n", zsyswapr_obj->info );
    }
    if( zsyswapr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zsyswapr is wrong\n", 
        zsyswapr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zsyswapr_obj->diff =  computeDiff_z( zsyswapr_obj->bufsize, 
                zsyswapr_obj->A, zsyswapr_obj->Aref );

}

TEST_F(zsyswapr_test, zsyswapr1) {
    EXPECT_NEAR(0.0, zsyswapr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zsyswapr_test, zsyswapr2) {
    EXPECT_NEAR(0.0, zsyswapr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zsyswapr_test, zsyswapr3) {
    EXPECT_NEAR(0.0, zsyswapr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zsyswapr_test, zsyswapr4) {
    EXPECT_NEAR(0.0, zsyswapr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}