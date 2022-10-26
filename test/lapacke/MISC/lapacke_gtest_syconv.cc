#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"

#define syconv_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (e!=NULL)    free(e); \
if (eref!=NULL)    free(eref); \
if (ipiv!=NULL)    free(ipiv); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule); \
if(bModule != NULL) dlclose(bModule); \
if(lModule != NULL) dlclose(lModule); \
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
/* Begin float_common_parameters  class definition */
class syconv_float_parameters{

   public:
	int bufsize;
	int bufsize_e;
	void *hModule, *dModule;
	void *bModule, *lModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;	
	float* A;
	float* e;
	lapack_int *ipiv;
	char uplo, way;
	lapack_int lda;
	float* Aref, *eref;
	lapack_int *ipivref;
	/*Return Values*/
	lapack_int info, inforef;
	lapack_int info_sytrf, inforef_sytrf;

   public:
      syconv_float_parameters (int matrix_layout, char uplo, char way, lapack_int n,  lapack_int lda);
      ~syconv_float_parameters ();

};

/* Constructor definition  float_common_parameters */
syconv_float_parameters:: syconv_float_parameters (int matrix_layout_i, char uplo_i, char way_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	way = way_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n syconv float: matrix_layout = %d, uplo:%c, n: %d, way:%c, lda:%d \n", matrix_layout, uplo, n, way, lda);
	#endif

	bufsize = lda*n;
	bufsize_e = n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_float_buffer_pair(&e, &eref, bufsize_e);
	lapacke_gtest_alloc_int_buffer_pair(&ipiv, &ipivref, n);
	if ((A==NULL) || (Aref==NULL)||
		(e == NULL) ||(eref == NULL)||
		(ipiv == NULL) || (ipivref == NULL)) {	
		EXPECT_FALSE( true) << "syconv_float_parameters object: malloc error.";
		syconv_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand_custom_matrix( A, Aref, bufsize, n, uplo);
	

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
syconv_float_parameters :: ~syconv_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" syconv_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   syconv_free();

}
/*  Test fixture class definition */
class ssyconv_test  : public  ::testing::Test {
public:
   syconv_float_parameters  *ssyconv_obj;
   void SetUp();
   void TearDown () { delete ssyconv_obj;}
};

void ssyconv_test::SetUp(){

	/* LAPACKE ssytrf prototype */
    typedef int (*Fptr_NL_LAPACKE_ssytrf) (int matrix_layout , char uplo , lapack_int n , float * A , lapack_int lda , lapack_int * ipiv);
	 Fptr_NL_LAPACKE_ssytrf ssytrf;
	 /* LAPACKE ssyconv prototype */
    typedef int (*Fptr_NL_LAPACKE_ssyconv) (int matrix_layout, char uplo, char way, lapack_int n, float * A,\
	lapack_int lda, const lapack_int * ipiv, float * e);

    Fptr_NL_LAPACKE_ssyconv ssyconv;

    ssyconv_obj = new syconv_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           eig_paramslist[idx].storev,
						   lin_solver_paramslist[idx].n,
						   lin_solver_paramslist[idx].ldb);

    idx = Circular_Increment_Index(idx);

    ssyconv_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ssyconv_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ssyconv_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ssyconv_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*ssyconv library call */
    ssyconv = (Fptr_NL_LAPACKE_ssyconv)dlsym(ssyconv_obj->hModule, "LAPACKE_ssyconv");
    ASSERT_TRUE(ssyconv != NULL) << "failed to get the Netlib LAPACKE_ssyconv symbol";

	/*ssytrf library call*/
	ssyconv_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ssyconv_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ssyconv_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ssyconv_obj->lModule != NULL) << "Netlib lapacke handle NULL";

    ssytrf = (Fptr_NL_LAPACKE_ssytrf)dlsym(ssyconv_obj->lModule, "LAPACKE_ssytrf");
    ASSERT_TRUE(ssytrf != NULL) << "failed to get the Netlib LAPACKE_ssytrf symbol";    
    

    ssyconv_obj->inforef_sytrf = ssytrf( ssyconv_obj->matrix_layout, ssyconv_obj->uplo,
								ssyconv_obj->n, ssyconv_obj->Aref, ssyconv_obj->lda, ssyconv_obj->ipivref);

    /* Compute libflame's Lapacke o/p  */
    ssyconv_obj->info_sytrf = LAPACKE_ssytrf( ssyconv_obj->matrix_layout, ssyconv_obj->uplo,
										ssyconv_obj->n,ssyconv_obj->A, ssyconv_obj->lda, ssyconv_obj->ipiv);

    if( ssyconv_obj->info_sytrf < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ssytrf is wrong\n", ssyconv_obj->info_sytrf );
    }
    if( ssyconv_obj->inforef_sytrf < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ssytrf is wrong\n", 
        ssyconv_obj->inforef_sytrf );
    }  
/*Compute ssyconv's  o/p */
    ssyconv_obj->inforef = ssyconv( ssyconv_obj->matrix_layout, ssyconv_obj->uplo, ssyconv_obj->way,
								ssyconv_obj->n, ssyconv_obj->Aref, ssyconv_obj->lda, (const int*)ssyconv_obj->ipivref, ssyconv_obj->eref);

    /* Compute libflame's Lapacke o/p  */
	
    ssyconv_obj->info = LAPACKE_ssyconv( ssyconv_obj->matrix_layout, ssyconv_obj->uplo,ssyconv_obj->way,
										ssyconv_obj->n,ssyconv_obj->A, ssyconv_obj->lda,  (const int*)ssyconv_obj->ipiv, ssyconv_obj->e);

    if( ssyconv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ssyconv is wrong\n", ssyconv_obj->info );
    }
    if( ssyconv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ssyconv is wrong\n", 
        ssyconv_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    ssyconv_obj->diff =  computeDiff_s( ssyconv_obj->bufsize_e, 
                ssyconv_obj->e, ssyconv_obj->eref );

}

TEST_F(ssyconv_test, ssyconv1) {
    EXPECT_NEAR(0.0, ssyconv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssyconv_test, ssyconv2) {
    EXPECT_NEAR(0.0, ssyconv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssyconv_test, ssyconv3) {
    EXPECT_NEAR(0.0, ssyconv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssyconv_test, ssyconv4) {
    EXPECT_NEAR(0.0, ssyconv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}
/* Begin double_common_parameters  class definition */
class syconv_double_parameters{

   public:
	int bufsize;
	int bufsize_e;
	void *hModule, *dModule;
	void *bModule, *lModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;	
	double* A;
	double* e;
	lapack_int *ipiv;
	char uplo, way;
	lapack_int lda;
	double* Aref, *eref;
	lapack_int *ipivref;
	/*Return Values*/
	lapack_int info, inforef;
	lapack_int info_sytrf, inforef_sytrf;

   public:
      syconv_double_parameters (int matrix_layout, char uplo, char way, lapack_int n,  lapack_int lda);
      ~syconv_double_parameters ();

};

/* Constructor definition  double_common_parameters */
syconv_double_parameters:: syconv_double_parameters (int matrix_layout_i, char uplo_i, char way_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	way = way_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n syconv double: matrix_layout = %d, uplo:%c, n: %d, way:%c, lda:%d \n", matrix_layout, uplo, n, way, lda);
	#endif

	bufsize = lda*n;
	bufsize_e = n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_double_buffer_pair(&e, &eref, bufsize_e);
	lapacke_gtest_alloc_int_buffer_pair(&ipiv, &ipivref, n);
	if ((A==NULL) || (Aref==NULL)||
		(e == NULL) ||(eref == NULL)||
		(ipiv == NULL) || (ipivref == NULL)) {	
		EXPECT_FALSE( true) << "syconv_double_parameters object: malloc error.";
		syconv_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand_custom_matrix( A, Aref, bufsize, n, uplo);
	

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
syconv_double_parameters :: ~syconv_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" syconv_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   syconv_free();

}
/*  Test fixture class definition */
class dsyconv_test  : public  ::testing::Test {
public:
   syconv_double_parameters  *dsyconv_obj;
   void SetUp();
   void TearDown () { delete dsyconv_obj;}
};

void dsyconv_test::SetUp(){

	/* LAPACKE dsytrf prototype */
    typedef int (*Fptr_NL_LAPACKE_dsytrf) (int matrix_layout , char uplo , lapack_int n , double * A , lapack_int lda , lapack_int * ipiv);
	 Fptr_NL_LAPACKE_dsytrf dsytrf;
	 /* LAPACKE dsyconv prototype */
    typedef int (*Fptr_NL_LAPACKE_dsyconv) (int matrix_layout, char uplo, char way, lapack_int n, double * A,\
	lapack_int lda, const lapack_int * ipiv, double * e);

    Fptr_NL_LAPACKE_dsyconv dsyconv;

    dsyconv_obj = new syconv_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           eig_paramslist[idx].storev,
						   lin_solver_paramslist[idx].n,
						   lin_solver_paramslist[idx].ldb);

    idx = Circular_Increment_Index(idx);

    dsyconv_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dsyconv_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dsyconv_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dsyconv_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*dsyconv library call */
    dsyconv = (Fptr_NL_LAPACKE_dsyconv)dlsym(dsyconv_obj->hModule, "LAPACKE_dsyconv");
    ASSERT_TRUE(dsyconv != NULL) << "failed to get the Netlib LAPACKE_dsyconv symbol";

	/*dsytrf library call*/
	dsyconv_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dsyconv_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dsyconv_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dsyconv_obj->lModule != NULL) << "Netlib lapacke handle NULL";

    dsytrf = (Fptr_NL_LAPACKE_dsytrf)dlsym(dsyconv_obj->lModule, "LAPACKE_dsytrf");
    ASSERT_TRUE(dsytrf != NULL) << "failed to get the Netlib LAPACKE_dsytrf symbol";    
    

    dsyconv_obj->inforef_sytrf = dsytrf( dsyconv_obj->matrix_layout, dsyconv_obj->uplo,
								dsyconv_obj->n, dsyconv_obj->Aref, dsyconv_obj->lda, dsyconv_obj->ipivref);

    /* Compute libflame's Lapacke o/p  */
    dsyconv_obj->info_sytrf = LAPACKE_dsytrf( dsyconv_obj->matrix_layout, dsyconv_obj->uplo,
										dsyconv_obj->n,dsyconv_obj->A, dsyconv_obj->lda, dsyconv_obj->ipiv);

    if( dsyconv_obj->info_sytrf < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dsytrf is wrong\n", dsyconv_obj->info_sytrf );
    }
    if( dsyconv_obj->inforef_sytrf < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dsytrf is wrong\n", 
        dsyconv_obj->inforef_sytrf );
    }  
/*Compute dsyconv's  o/p */
    dsyconv_obj->inforef = dsyconv( dsyconv_obj->matrix_layout, dsyconv_obj->uplo, dsyconv_obj->way,
								dsyconv_obj->n, dsyconv_obj->Aref, dsyconv_obj->lda, (const int*)dsyconv_obj->ipivref, dsyconv_obj->eref);

    /* Compute libflame's Lapacke o/p  */
	
    dsyconv_obj->info = LAPACKE_dsyconv( dsyconv_obj->matrix_layout, dsyconv_obj->uplo,dsyconv_obj->way,
										dsyconv_obj->n,dsyconv_obj->A, dsyconv_obj->lda,  (const int*)dsyconv_obj->ipiv, dsyconv_obj->e);

    if( dsyconv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dsyconv is wrong\n", dsyconv_obj->info );
    }
    if( dsyconv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dsyconv is wrong\n", 
        dsyconv_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dsyconv_obj->diff =  computeDiff_d( dsyconv_obj->bufsize_e, 
                dsyconv_obj->e, dsyconv_obj->eref );

}

TEST_F(dsyconv_test, dsyconv1) {
    EXPECT_NEAR(0.0, dsyconv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dsyconv_test, dsyconv2) {
    EXPECT_NEAR(0.0, dsyconv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dsyconv_test, dsyconv3) {
    EXPECT_NEAR(0.0, dsyconv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dsyconv_test, dsyconv4) {
    EXPECT_NEAR(0.0, dsyconv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}



/* Begin scomplex_common_parameters  class definition */
class syconv_scomplex_parameters{

   public:
	int bufsize;
	int bufsize_e;
	void *hModule, *dModule;
	void *bModule, *lModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;	
	lapack_complex_float* A;
	lapack_complex_float* e;
	lapack_int *ipiv;
	char uplo, way;
	lapack_int lda;
	lapack_complex_float* Aref, *eref;
	lapack_int *ipivref;
	/*Return Values*/
	lapack_int info, inforef;
	lapack_int info_sytrf, inforef_sytrf;

   public:
      syconv_scomplex_parameters (int matrix_layout, char uplo, char way, lapack_int n,  lapack_int lda);
      ~syconv_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
syconv_scomplex_parameters:: syconv_scomplex_parameters (int matrix_layout_i, char uplo_i, char way_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	way = way_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n syconv scomplex: matrix_layout = %d, uplo:%c, n: %d, way:%c, lda:%d \n", matrix_layout, uplo, n, way, lda);
	#endif

	bufsize = lda*n;
	bufsize_e = n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&e, &eref, bufsize_e);
	lapacke_gtest_alloc_int_buffer_pair(&ipiv, &ipivref, n);
	if ((A==NULL) || (Aref==NULL)||
		(e == NULL) ||(eref == NULL)||
		(ipiv == NULL) || (ipivref == NULL)) {	
		EXPECT_FALSE( true) << "syconv_float_parameters object: malloc error.";
		syconv_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix( A, Aref, bufsize, n, uplo);
	

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
syconv_scomplex_parameters :: ~syconv_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" syconv_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   syconv_free();

}
/*  Test fixture class definition */
class csyconv_test  : public  ::testing::Test {
public:
   syconv_scomplex_parameters  *csyconv_obj;
   void SetUp();
   void TearDown () { delete csyconv_obj;}
};

void csyconv_test::SetUp(){

	/* LAPACKE csytrf prototype */
    typedef int (*Fptr_NL_LAPACKE_csytrf) (int matrix_layout , char uplo , lapack_int n , lapack_complex_float * A , lapack_int lda , lapack_int * ipiv);
	 Fptr_NL_LAPACKE_csytrf csytrf;
	 /* LAPACKE csyconv prototype */
    typedef int (*Fptr_NL_LAPACKE_csyconv) (int matrix_layout, char uplo, char way, lapack_int n, lapack_complex_float * A,\
	lapack_int lda, const lapack_int * ipiv, lapack_complex_float * e);

    Fptr_NL_LAPACKE_csyconv csyconv;

    csyconv_obj = new syconv_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           eig_paramslist[idx].storev,
						   lin_solver_paramslist[idx].n,
						   lin_solver_paramslist[idx].ldb);

    idx = Circular_Increment_Index(idx);

    csyconv_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    csyconv_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(csyconv_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(csyconv_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*csyconv library call */
    csyconv = (Fptr_NL_LAPACKE_csyconv)dlsym(csyconv_obj->hModule, "LAPACKE_csyconv");
    ASSERT_TRUE(csyconv != NULL) << "failed to get the Netlib LAPACKE_csyconv symbol";

	/*csytrf library call*/
	csyconv_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    csyconv_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(csyconv_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(csyconv_obj->lModule != NULL) << "Netlib lapacke handle NULL";

    csytrf = (Fptr_NL_LAPACKE_csytrf)dlsym(csyconv_obj->lModule, "LAPACKE_csytrf");
    ASSERT_TRUE(csytrf != NULL) << "failed to get the Netlib LAPACKE_csytrf symbol";    
    

    csyconv_obj->inforef_sytrf = csytrf( csyconv_obj->matrix_layout, csyconv_obj->uplo,
								csyconv_obj->n, csyconv_obj->Aref, csyconv_obj->lda, csyconv_obj->ipivref);

    /* Compute libflame's Lapacke o/p  */
    csyconv_obj->info_sytrf = LAPACKE_csytrf( csyconv_obj->matrix_layout, csyconv_obj->uplo,
										csyconv_obj->n,csyconv_obj->A, csyconv_obj->lda, csyconv_obj->ipiv);

    if( csyconv_obj->info_sytrf < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_csytrf is wrong\n", csyconv_obj->info_sytrf );
    }
    if( csyconv_obj->inforef_sytrf < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_csytrf is wrong\n", 
        csyconv_obj->inforef_sytrf );
    }  
/*Compute csyconv's  o/p */
    csyconv_obj->inforef = csyconv( csyconv_obj->matrix_layout, csyconv_obj->uplo, csyconv_obj->way,
								csyconv_obj->n, csyconv_obj->Aref, csyconv_obj->lda, (const int*)csyconv_obj->ipivref, csyconv_obj->eref);

    /* Compute libflame's Lapacke o/p  */
	
    csyconv_obj->info = LAPACKE_csyconv( csyconv_obj->matrix_layout, csyconv_obj->uplo,csyconv_obj->way,
										csyconv_obj->n,csyconv_obj->A, csyconv_obj->lda,  (const int*)csyconv_obj->ipiv, csyconv_obj->e);

    if( csyconv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_csyconv is wrong\n", csyconv_obj->info );
    }
    if( csyconv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_csyconv is wrong\n", 
        csyconv_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    csyconv_obj->diff =  computeDiff_c( csyconv_obj->bufsize_e, 
                csyconv_obj->e, csyconv_obj->eref );

}

TEST_F(csyconv_test, csyconv1) {
    EXPECT_NEAR(0.0, csyconv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(csyconv_test, csyconv2) {
    EXPECT_NEAR(0.0, csyconv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(csyconv_test, csyconv3) {
    EXPECT_NEAR(0.0, csyconv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(csyconv_test, csyconv4) {
    EXPECT_NEAR(0.0, csyconv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class syconv_dcomplex_parameters{

   public:
	int bufsize;
	int bufsize_e;
	void *hModule, *dModule;
	void *bModule, *lModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;	
	lapack_complex_double* A;
	lapack_complex_double* e;
	lapack_int *ipiv;
	char uplo, way;
	lapack_int lda;
	lapack_complex_double* Aref, *eref;
	lapack_int *ipivref;
	/*Return Values*/
	lapack_int info, inforef;
	lapack_int info_sytrf, inforef_sytrf;

   public:
      syconv_dcomplex_parameters (int matrix_layout, char uplo, char way, lapack_int n,  lapack_int lda);
      ~syconv_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
syconv_dcomplex_parameters:: syconv_dcomplex_parameters (int matrix_layout_i, char uplo_i, char way_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	way = way_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n syconv dcomplex: matrix_layout = %d, uplo:%c, n: %d, way:%c, lda:%d \n", matrix_layout, uplo, n, way, lda);
	#endif

	bufsize = lda*n;
	bufsize_e = n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&e, &eref, bufsize_e);
	lapacke_gtest_alloc_int_buffer_pair(&ipiv, &ipivref, n);
	if ((A==NULL) || (Aref==NULL)||
		(e == NULL) ||(eref == NULL)||
		(ipiv == NULL) || (ipivref == NULL)) {	
		EXPECT_FALSE( true) << "syconv_double_parameters object: malloc error.";
		syconv_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix( A, Aref, bufsize, n, uplo);
	

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
syconv_dcomplex_parameters :: ~syconv_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" syconv_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   syconv_free();

}
/*  Test fixture class definition */
class zsyconv_test  : public  ::testing::Test {
public:
   syconv_dcomplex_parameters  *zsyconv_obj;
   void SetUp();
   void TearDown () { delete zsyconv_obj;}
};

void zsyconv_test::SetUp(){

	/* LAPACKE zsytrf prototype */
    typedef int (*Fptr_NL_LAPACKE_zsytrf) (int matrix_layout , char uplo , lapack_int n , lapack_complex_double * A , lapack_int lda , lapack_int * ipiv);
	 Fptr_NL_LAPACKE_zsytrf zsytrf;
	 /* LAPACKE zsyconv prototype */
    typedef int (*Fptr_NL_LAPACKE_zsyconv) (int matrix_layout, char uplo, char way, lapack_int n, lapack_complex_double * A,\
	lapack_int lda, const lapack_int * ipiv, lapack_complex_double * e);

    Fptr_NL_LAPACKE_zsyconv zsyconv;

    zsyconv_obj = new syconv_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           eig_paramslist[idx].storev,
						   lin_solver_paramslist[idx].n,
						   lin_solver_paramslist[idx].ldb);

    idx = Circular_Increment_Index(idx);

    zsyconv_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zsyconv_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zsyconv_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zsyconv_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*zsyconv library call */
    zsyconv = (Fptr_NL_LAPACKE_zsyconv)dlsym(zsyconv_obj->hModule, "LAPACKE_zsyconv");
    ASSERT_TRUE(zsyconv != NULL) << "failed to get the Netlib LAPACKE_zsyconv symbol";

	/*zsytrf library call*/
	zsyconv_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zsyconv_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zsyconv_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zsyconv_obj->lModule != NULL) << "Netlib lapacke handle NULL";

    zsytrf = (Fptr_NL_LAPACKE_zsytrf)dlsym(zsyconv_obj->lModule, "LAPACKE_zsytrf");
    ASSERT_TRUE(zsytrf != NULL) << "failed to get the Netlib LAPACKE_zsytrf symbol";    
    

    zsyconv_obj->inforef_sytrf = zsytrf( zsyconv_obj->matrix_layout, zsyconv_obj->uplo,
								zsyconv_obj->n, zsyconv_obj->Aref, zsyconv_obj->lda, zsyconv_obj->ipivref);

    /* Compute libflame's Lapacke o/p  */
    zsyconv_obj->info_sytrf = LAPACKE_zsytrf( zsyconv_obj->matrix_layout, zsyconv_obj->uplo,
										zsyconv_obj->n,zsyconv_obj->A, zsyconv_obj->lda, zsyconv_obj->ipiv);

    if( zsyconv_obj->info_sytrf < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zsytrf is wrong\n", zsyconv_obj->info_sytrf );
    }
    if( zsyconv_obj->inforef_sytrf < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zsytrf is wrong\n", 
        zsyconv_obj->inforef_sytrf );
    }  
/*Compute zsyconv's  o/p */
    zsyconv_obj->inforef = zsyconv( zsyconv_obj->matrix_layout, zsyconv_obj->uplo, zsyconv_obj->way,
								zsyconv_obj->n, zsyconv_obj->Aref, zsyconv_obj->lda, (const int*)zsyconv_obj->ipivref, zsyconv_obj->eref);

    /* Compute libflame's Lapacke o/p  */
	
    zsyconv_obj->info = LAPACKE_zsyconv( zsyconv_obj->matrix_layout, zsyconv_obj->uplo,zsyconv_obj->way,
										zsyconv_obj->n,zsyconv_obj->A, zsyconv_obj->lda,  (const int*)zsyconv_obj->ipiv, zsyconv_obj->e);

    if( zsyconv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zsyconv is wrong\n", zsyconv_obj->info );
    }
    if( zsyconv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zsyconv is wrong\n", 
        zsyconv_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zsyconv_obj->diff =  computeDiff_z( zsyconv_obj->bufsize_e, 
                zsyconv_obj->e, zsyconv_obj->eref );

}

TEST_F(zsyconv_test, zsyconv1) {
    EXPECT_NEAR(0.0, zsyconv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zsyconv_test, zsyconv2) {
    EXPECT_NEAR(0.0, zsyconv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zsyconv_test, zsyconv3) {
    EXPECT_NEAR(0.0, zsyconv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zsyconv_test, zsyconv4) {
    EXPECT_NEAR(0.0, zsyconv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}