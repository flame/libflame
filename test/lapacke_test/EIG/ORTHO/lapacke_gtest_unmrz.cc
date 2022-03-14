#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define unmrz_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (tau!=NULL)  free(tau);\
if (tauref!=NULL) free(tauref); \
if (c != NULL) free(c); \
if (cref != NULL) free(cref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule);\
if( bModule != NULL) dlclose(bModule); \
if(lModule != NULL) dlclose(lModule);
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
/* Begin scomplex_common_parameters  class definition */
class unmrz_scomplex_parameters{

   public:
	int bufsize_a;
	int bufsize_c;
	void *hModule, *dModule;
	void *bModule, *lModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int k;
	lapack_int l;
	lapack_complex_float* A;
	char side;
	char trans;
	lapack_int ldc, lda, lda_tzrzf;
	/*Output Parameter*/
	lapack_complex_float* tau, *c;	
	lapack_complex_float *Aref, *tauref, *cref;
	/*Return Values*/	
	lapack_int info, inforef;
	lapack_int info_tzrzf, inforef_tzrzf;

   public:
      unmrz_scomplex_parameters (int matrix_layout, char side , char trans, lapack_int m, lapack_int n);
      ~unmrz_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
unmrz_scomplex_parameters:: unmrz_scomplex_parameters (int matrix_layout_i, char side_i , char trans_i,lapack_int m_i, lapack_int n_i)
{
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	side = side_i;
	trans = trans_i;
	
	if (trans == 'T')
		trans = 'N';
	#if LAPACKE_TEST_VERBOSE
		printf(" \n unmrz scomplex: matrix_layout = %d, side :%c, trans:%c, m:%d, n: %d \n", matrix_layout, side, trans, m, n);
	#endif
	if (side == 'L')
	{
		k = m;
		l = m;
		if (matrix_layout == LAPACK_COL_MAJOR)
		{		
			ldc = m;
			//lda = k;
			lda_tzrzf = m;
			lda = m;
			bufsize_a = (lda_tzrzf*n);
			bufsize_c = (ldc*n);			
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{		
			ldc = n;
			//lda = m;
			lda_tzrzf =n;
			lda = n;
			bufsize_a = (lda_tzrzf*k);
			bufsize_c = (ldc*m);
		}else
			
		{
			EXPECT_TRUE(false) << "matrix_layout invalid";
		}
	}else if (side == 'R') {
		k = n;
		l = n;
		if (matrix_layout == LAPACK_COL_MAJOR)
		{	ldc = m;
			lda = k;			
			bufsize_c = (ldc*n);
			lda_tzrzf = m;

			bufsize_a = (lda_tzrzf*n);
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{	ldc = n;
			lda = n;
			lda_tzrzf =n;			
			bufsize_a = (lda_tzrzf*m);
			bufsize_c = (ldc*m);
		}else
		{
			EXPECT_TRUE(false) << "matrix_layout invalid";
		}
	}else 
		EXPECT_TRUE(false) << "side is invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&tau, &tauref, k);	
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&c, &cref, bufsize_c);
	
	if ((A==NULL) || (Aref==NULL) ||
		(tau==NULL) || (tauref==NULL) ||
		(c==NULL) || (cref==NULL))
	{
		EXPECT_FALSE( true) << "unmrz_float_parameters object: malloc error.";
		unmrz_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( A, Aref, bufsize_a);
	lapacke_gtest_init_scomplex_buffer_pair_rand( tau, tauref, k);


} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
unmrz_scomplex_parameters :: ~unmrz_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" unmrz_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   unmrz_free();

}
/*  Test fixture class definition */
class cunmrz_test  : public  ::testing::Test {
public:
   unmrz_scomplex_parameters  *cunmrz_obj;
   void SetUp();
   void TearDown () { delete cunmrz_obj;}
};

void cunmrz_test::SetUp(){

	/* LAPACKE ctzrzf prototype */
    typedef int (*Fptr_NL_LAPACKE_ctzrzf) (int matrix_layout, lapack_int m, lapack_int n, lapack_complex_float* a, lapack_int lda, lapack_complex_float* tau);
	 Fptr_NL_LAPACKE_ctzrzf ctzrzf;
	 /* LAPACKE cunmrz prototype */
    typedef int (*Fptr_NL_LAPACKE_cunmrz) (int matrix_layout, char side, char trans, lapack_int m, lapack_int n, lapack_int k, lapack_int l,\
	const lapack_complex_float* a, lapack_int lda, const lapack_complex_float* tau, lapack_complex_float* c, lapack_int ldc);

    Fptr_NL_LAPACKE_cunmrz cunmrz;
	
    cunmrz_obj = new unmrz_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].side,
						   eig_paramslist[idx].trans,
						   eig_paramslist[idx].n,
                           eig_paramslist[idx].m);
						   

    idx = Circular_Increment_Index(idx);

    cunmrz_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cunmrz_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cunmrz_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cunmrz_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*cunmrz library call */
    cunmrz = (Fptr_NL_LAPACKE_cunmrz)dlsym(cunmrz_obj->hModule, "LAPACKE_cunmrz");
    ASSERT_TRUE(cunmrz != NULL) << "failed to get the Netlib LAPACKE_cunmrz symbol";

	/*ctzrzf library call*/
	cunmrz_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cunmrz_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cunmrz_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cunmrz_obj->lModule != NULL) << "Netlib lapacke handle NULL";

    ctzrzf = (Fptr_NL_LAPACKE_ctzrzf)dlsym(cunmrz_obj->lModule, "LAPACKE_ctzrzf");
    ASSERT_TRUE(ctzrzf != NULL) << "failed to get the Netlib LAPACKE_ctzrzf symbol";    
    

    cunmrz_obj->inforef_tzrzf = ctzrzf( cunmrz_obj->matrix_layout,cunmrz_obj->m, cunmrz_obj->n,
								cunmrz_obj->Aref, cunmrz_obj->lda_tzrzf, cunmrz_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    cunmrz_obj->info_tzrzf = LAPACKE_ctzrzf( cunmrz_obj->matrix_layout,cunmrz_obj->m, cunmrz_obj->n,
											cunmrz_obj->A, cunmrz_obj->lda_tzrzf, cunmrz_obj->tau);
										
										

    if( cunmrz_obj->info_tzrzf < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ctzrzf is wrong\n", cunmrz_obj->info_tzrzf );
    }
    if( cunmrz_obj->inforef_tzrzf < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ctzrzf is wrong\n", 
        cunmrz_obj->inforef_tzrzf );
    }  
/*Compute cunmrz's  o/p */
    cunmrz_obj->inforef = cunmrz( cunmrz_obj->matrix_layout, cunmrz_obj->side, cunmrz_obj->trans,
								cunmrz_obj->m, cunmrz_obj->n, cunmrz_obj->k, cunmrz_obj->l, (const lapack_complex_float*)cunmrz_obj->Aref,
								cunmrz_obj->lda, (const lapack_complex_float*)cunmrz_obj->tauref, cunmrz_obj->cref, cunmrz_obj->ldc);

    /* Compute libflame's Lapacke o/p  */
    cunmrz_obj->info = LAPACKE_cunmrz( cunmrz_obj->matrix_layout, cunmrz_obj->side, cunmrz_obj->trans,
									cunmrz_obj->m, cunmrz_obj->n, cunmrz_obj->k, cunmrz_obj->l, (const lapack_complex_float*)cunmrz_obj->A,
									cunmrz_obj->lda, (const lapack_complex_float*)cunmrz_obj->tau, cunmrz_obj->c, cunmrz_obj->ldc);
    if( cunmrz_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame\
        LAPACKE_cunmrz is wrong\n", cunmrz_obj->info );
    }
    if( cunmrz_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cunmrz is wrong\n", 
        cunmrz_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cunmrz_obj->diff =  computeDiff_c( cunmrz_obj->bufsize_c, 
                cunmrz_obj->c, cunmrz_obj->cref );

}

TEST_F(cunmrz_test, cunmrz1) {
    EXPECT_NEAR(0.0, cunmrz_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(cunmrz_test, cunmrz2) {
    EXPECT_NEAR(0.0, cunmrz_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(cunmrz_test, cunmrz3) {
    EXPECT_NEAR(0.0, cunmrz_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(cunmrz_test, cunmrz4) {
    EXPECT_NEAR(0.0, cunmrz_obj->diff, LAPACKE_EIG_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class unmrz_dcomplex_parameters{

   public:
	int bufsize_a;
	int bufsize_c;
	void *hModule, *dModule;
	void *bModule, *lModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int k;
	lapack_int l;
	lapack_complex_double* A;
	char side;
	char trans;
	lapack_int ldc, lda, lda_tzrzf;
	/*Output Parameter*/
	lapack_complex_double* tau, *c;	
	lapack_complex_double *Aref, *tauref, *cref;
	/*Return Values*/	
	lapack_int info, inforef;
	lapack_int info_tzrzf, inforef_tzrzf;

   public:
      unmrz_dcomplex_parameters (int matrix_layout, char side , char trans, lapack_int m, lapack_int n);
      ~unmrz_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
unmrz_dcomplex_parameters:: unmrz_dcomplex_parameters (int matrix_layout_i, char side_i , char trans_i,lapack_int m_i, lapack_int n_i)
{
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	side = side_i;
	trans = trans_i;
	
	if (trans == 'T')
		trans = 'N';
	#if LAPACKE_TEST_VERBOSE
		printf(" \n unmrz dcomplex: matrix_layout = %d, side :%c, trans:%c, m:%d, n: %d \n", matrix_layout, side, trans, m, n);
	#endif
	if (side == 'L')
	{
		k = m;
		l = m;
		if (matrix_layout == LAPACK_COL_MAJOR)
		{		
			ldc = m;
			lda = k;
			lda_tzrzf = m;
			bufsize_a = (lda_tzrzf*n);
			bufsize_c = (ldc*n);			
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{		
			ldc = n;
			lda = m;
			lda_tzrzf =n;
			bufsize_a = (lda_tzrzf*k);
			bufsize_c = (ldc*m);
		}else
			
		{
			EXPECT_TRUE(false) << "matrix_layout invalid";
		}
	}else if (side == 'R') {
		k = n;
		l = n;
		if (matrix_layout == LAPACK_COL_MAJOR)
		{	ldc = m;
			lda = k;			
			bufsize_c = (ldc*n);
			lda_tzrzf = m;
			bufsize_a = (lda_tzrzf*n);
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{	ldc = n;
			lda = n;
			lda_tzrzf =n;
			bufsize_a = (lda_tzrzf*m);
			bufsize_c = (ldc*m);
		}else
		{
			EXPECT_TRUE(false) << "matrix_layout invalid";
		}
	}else 
		EXPECT_TRUE(false) << "side is invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&tau, &tauref, k);	
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&c, &cref, bufsize_c);
	
	if ((A==NULL) || (Aref==NULL) ||
		(tau==NULL) || (tauref==NULL) ||
		(c==NULL) || (cref==NULL))
	{
		EXPECT_FALSE( true) << "unmrz_double_parameters object: malloc error.";
		unmrz_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( A, Aref, bufsize_a);
	lapacke_gtest_init_dcomplex_buffer_pair_rand( tau, tauref, k);


} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
unmrz_dcomplex_parameters :: ~unmrz_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" unmrz_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   unmrz_free();

}
/*  Test fixture class definition */
class zunmrz_test  : public  ::testing::Test {
public:
   unmrz_dcomplex_parameters  *zunmrz_obj;
   void SetUp();
   void TearDown () { delete zunmrz_obj;}
};

void zunmrz_test::SetUp(){

	/* LAPACKE ztzrzf prototype */
    typedef int (*Fptr_NL_LAPACKE_ztzrzf) (int matrix_layout, lapack_int m, lapack_int n, lapack_complex_double* a, lapack_int lda, lapack_complex_double* tau);
	 Fptr_NL_LAPACKE_ztzrzf ztzrzf;
	 /* LAPACKE zunmrz prototype */
    typedef int (*Fptr_NL_LAPACKE_zunmrz) (int matrix_layout, char side, char trans, lapack_int m, lapack_int n, lapack_int k, lapack_int l,\
	const lapack_complex_double* a, lapack_int lda, const lapack_complex_double* tau, lapack_complex_double* c, lapack_int ldc);

    Fptr_NL_LAPACKE_zunmrz zunmrz;
	
    zunmrz_obj = new unmrz_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].side,
						   eig_paramslist[idx].trans,
						   eig_paramslist[idx].n,
                           eig_paramslist[idx].m);
						   

    idx = Circular_Increment_Index(idx);

    zunmrz_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zunmrz_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zunmrz_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zunmrz_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*zunmrz library call */
    zunmrz = (Fptr_NL_LAPACKE_zunmrz)dlsym(zunmrz_obj->hModule, "LAPACKE_zunmrz");
    ASSERT_TRUE(zunmrz != NULL) << "failed to get the Netlib LAPACKE_zunmrz symbol";

	/*ztzrzf library call*/
	zunmrz_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zunmrz_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zunmrz_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zunmrz_obj->lModule != NULL) << "Netlib lapacke handle NULL";

    ztzrzf = (Fptr_NL_LAPACKE_ztzrzf)dlsym(zunmrz_obj->lModule, "LAPACKE_ztzrzf");
    ASSERT_TRUE(ztzrzf != NULL) << "failed to get the Netlib LAPACKE_ztzrzf symbol";    
    

    zunmrz_obj->inforef_tzrzf = ztzrzf( zunmrz_obj->matrix_layout,zunmrz_obj->m, zunmrz_obj->n,
								zunmrz_obj->Aref, zunmrz_obj->lda_tzrzf, zunmrz_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    zunmrz_obj->info_tzrzf = LAPACKE_ztzrzf( zunmrz_obj->matrix_layout,zunmrz_obj->m, zunmrz_obj->n,
											zunmrz_obj->A, zunmrz_obj->lda_tzrzf, zunmrz_obj->tau);
										
										

    if( zunmrz_obj->info_tzrzf < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ztzrzf is wrong\n", zunmrz_obj->info_tzrzf );
    }
    if( zunmrz_obj->inforef_tzrzf < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ztzrzf is wrong\n", 
        zunmrz_obj->inforef_tzrzf );
    }  
/*Compute zunmrz's  o/p */
    zunmrz_obj->inforef = zunmrz( zunmrz_obj->matrix_layout, zunmrz_obj->side, zunmrz_obj->trans,
								zunmrz_obj->m, zunmrz_obj->n, zunmrz_obj->k, zunmrz_obj->l, (const lapack_complex_double*)zunmrz_obj->Aref,
								zunmrz_obj->lda, (const lapack_complex_double*)zunmrz_obj->tauref, zunmrz_obj->cref, zunmrz_obj->ldc);

    /* Compute libflame's Lapacke o/p  */
    zunmrz_obj->info = LAPACKE_zunmrz( zunmrz_obj->matrix_layout, zunmrz_obj->side, zunmrz_obj->trans,
									zunmrz_obj->m, zunmrz_obj->n, zunmrz_obj->k, zunmrz_obj->l, (const lapack_complex_double*)zunmrz_obj->A,
									zunmrz_obj->lda, (const lapack_complex_double*)zunmrz_obj->tau, zunmrz_obj->c, zunmrz_obj->ldc);
    if( zunmrz_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame\
        LAPACKE_zunmrz is wrong\n", zunmrz_obj->info );
    }
    if( zunmrz_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zunmrz is wrong\n", 
        zunmrz_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zunmrz_obj->diff =  computeDiff_z( zunmrz_obj->bufsize_c, 
                zunmrz_obj->c, zunmrz_obj->cref );

}

TEST_F(zunmrz_test, zunmrz1) {
    EXPECT_NEAR(0.0, zunmrz_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zunmrz_test, zunmrz2) {
    EXPECT_NEAR(0.0, zunmrz_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zunmrz_test, zunmrz3) {
    EXPECT_NEAR(0.0, zunmrz_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zunmrz_test, zunmrz4) {
    EXPECT_NEAR(0.0, zunmrz_obj->diff, LAPACKE_EIG_THRESHOLD);
}