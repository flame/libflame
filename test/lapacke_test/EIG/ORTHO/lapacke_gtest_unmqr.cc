#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"

#define unmqr_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (tau!=NULL)  free(tau);\
if (tauref!=NULL) free(tauref); \
if (c != NULL) free(c); \
if (cref != NULL) free(cref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule);\
if( bModule != NULL) dlclose(bModule); \
if(lModule != NULL) dlclose(lModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
/* Begin scomplex_common_parameters  class definition */
class unmqr_scomplex_parameters{

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
	lapack_complex_float* A;
	char side;
	char trans;
	lapack_int ldc, lda, lda_geqrf;
	/*Output Parameter*/
	lapack_complex_float* tau, *c;	
	lapack_complex_float *Aref, *tauref, *cref;
	/*Return Values*/	
	lapack_int info, inforef;
	lapack_int info_geqrf, inforef_geqrf;

   public:
      unmqr_scomplex_parameters (int matrix_layout, char side , char trans, lapack_int m, lapack_int n);
      ~unmqr_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
unmqr_scomplex_parameters:: unmqr_scomplex_parameters (int matrix_layout_i, char side_i , char trans_i,lapack_int m_i, lapack_int n_i)
{
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	side = side_i;
	trans = trans_i;
	
	if (trans == 'T')
		trans  = 'N' ;
	
	#if LAPACKE_TEST_VERBOSE
		printf(" \n unmqr scomplex: matrix_layout = %d, side :%c, trans:%c, m:%d, n: %d \n", matrix_layout, side, trans, m, n);
	#endif
	if (side == 'L')
	{
		k = m;
		if (matrix_layout == LAPACK_COL_MAJOR)
		{		
			ldc = m;
			lda = k;
			lda_geqrf = m;
			bufsize_a = (lda_geqrf*n);
			bufsize_c = (ldc*n);			
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{		
			ldc = n;
			lda = m;
			lda_geqrf =n;
			bufsize_a = (lda_geqrf*k);
			bufsize_c = (ldc*m);
		}else
		{
			EXPECT_TRUE(false) << "matrix_layout invalid";
		}
	}else if (side == 'R') {
		k = n;
		if (matrix_layout == LAPACK_COL_MAJOR)
		{	ldc = m;
			lda = k;
			lda_geqrf = m;
			bufsize_a = (lda_geqrf*n);
			bufsize_c = (ldc*n);			
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{	ldc = n;
			lda = n;
			lda_geqrf =n;
			bufsize_a = (lda_geqrf*m);
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
		EXPECT_FALSE( true) << "unmqr_float_parameters object: malloc error.";
		unmqr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( A, Aref, bufsize_a);
	lapacke_gtest_init_scomplex_buffer_pair_rand( c, cref, bufsize_c);


} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
unmqr_scomplex_parameters :: ~unmqr_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" unmqr_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   unmqr_free();

}
/*  Test fixture class definition */
class cunmqr_test  : public  ::testing::Test {
public:
   unmqr_scomplex_parameters  *cunmqr_obj;
   void SetUp();
   void TearDown () { delete cunmqr_obj;}
};

void cunmqr_test::SetUp(){

	/* LAPACKE cgeqrf prototype */
    typedef int (*Fptr_NL_LAPACKE_cgeqrf) (int matrix_layout, lapack_int m, lapack_int n,\
										   lapack_complex_float* a, lapack_int lda, lapack_complex_float* tau);
	 Fptr_NL_LAPACKE_cgeqrf cgeqrf;
	 /* LAPACKE cunmqr prototype */
    typedef int (*Fptr_NL_LAPACKE_cunmqr) (int matrix_layout, char side, char trans, lapack_int m, lapack_int n, \
											lapack_int k, const lapack_complex_float* a, lapack_int lda, \
											const lapack_complex_float* tau, lapack_complex_float* c, lapack_int ldc);

    Fptr_NL_LAPACKE_cunmqr cunmqr;

    cunmqr_obj = new unmqr_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].side,
						   eig_paramslist[idx].trans,
						   eig_paramslist[idx].m,
                           eig_paramslist[idx].n);
						   

    idx = Circular_Increment_Index(idx);

    cunmqr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cunmqr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cunmqr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cunmqr_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*cunmqr library call */
    cunmqr = (Fptr_NL_LAPACKE_cunmqr)dlsym(cunmqr_obj->hModule, "LAPACKE_cunmqr");
    ASSERT_TRUE(cunmqr != NULL) << "failed to get the Netlib LAPACKE_cunmqr symbol";

	/*cgeqrf library call*/
	cunmqr_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cunmqr_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cunmqr_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cunmqr_obj->lModule != NULL) << "Netlib lapacke handle NULL";

    cgeqrf = (Fptr_NL_LAPACKE_cgeqrf)dlsym(cunmqr_obj->lModule, "LAPACKE_cgeqrf");
    ASSERT_TRUE(cgeqrf != NULL) << "failed to get the Netlib LAPACKE_cgeqrf symbol";    
    

    cunmqr_obj->inforef_geqrf = cgeqrf( cunmqr_obj->matrix_layout,cunmqr_obj->m, cunmqr_obj->n,
								cunmqr_obj->Aref, cunmqr_obj->lda_geqrf, cunmqr_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    cunmqr_obj->info_geqrf = LAPACKE_cgeqrf( cunmqr_obj->matrix_layout,cunmqr_obj->m, cunmqr_obj->n,
											cunmqr_obj->A, cunmqr_obj->lda_geqrf, cunmqr_obj->tau);
										
										

    if( cunmqr_obj->info_geqrf < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cgeqrf is wrong\n", cunmqr_obj->info_geqrf );
    }
    if( cunmqr_obj->inforef_geqrf < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cgeqrf is wrong\n", 
        cunmqr_obj->inforef_geqrf );
    }  
/*Compute cunmqr's  o/p */
    cunmqr_obj->inforef = cunmqr( cunmqr_obj->matrix_layout, cunmqr_obj->side, cunmqr_obj->trans,
								cunmqr_obj->m, cunmqr_obj->n, cunmqr_obj->k, (const lapack_complex_float*)cunmqr_obj->Aref,
								cunmqr_obj->lda, (const lapack_complex_float*)cunmqr_obj->tauref, cunmqr_obj->cref, cunmqr_obj->ldc);

    /* Compute libflame's Lapacke o/p  */
    cunmqr_obj->info = LAPACKE_cunmqr( cunmqr_obj->matrix_layout, cunmqr_obj->side, cunmqr_obj->trans,
									cunmqr_obj->m, cunmqr_obj->n, cunmqr_obj->k, (const lapack_complex_float*)cunmqr_obj->A,
									cunmqr_obj->lda, (const lapack_complex_float*)cunmqr_obj->tau, cunmqr_obj->c, cunmqr_obj->ldc);
    if( cunmqr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame\
        LAPACKE_cunmqr is wrong\n", cunmqr_obj->info );
    }
    if( cunmqr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cunmqr is wrong\n", 
        cunmqr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cunmqr_obj->diff =  computeDiff_c( cunmqr_obj->bufsize_c, 
                cunmqr_obj->c, cunmqr_obj->cref );

}

TEST_F(cunmqr_test, cunmqr1) {
    EXPECT_NEAR(0.0, cunmqr_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(cunmqr_test, cunmqr2) {
    EXPECT_NEAR(0.0, cunmqr_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(cunmqr_test, cunmqr3) {
    EXPECT_NEAR(0.0, cunmqr_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(cunmqr_test, cunmqr4) {
    EXPECT_NEAR(0.0, cunmqr_obj->diff, LAPACKE_EIG_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class unmqr_dcomplex_parameters{

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
	lapack_complex_double* A;
	char side;
	char trans;
	lapack_int ldc, lda, lda_geqrf;
	/*Output Parameter*/
	lapack_complex_double* tau, *c;	
	lapack_complex_double *Aref, *tauref, *cref;
	/*Return Values*/	
	lapack_int info, inforef;
	lapack_int info_geqrf, inforef_geqrf;

   public:
      unmqr_dcomplex_parameters (int matrix_layout, char side , char trans, lapack_int m, lapack_int n);
      ~unmqr_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
unmqr_dcomplex_parameters:: unmqr_dcomplex_parameters (int matrix_layout_i, char side_i , char trans_i,lapack_int m_i, lapack_int n_i)
{
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	side = side_i;
	trans = trans_i;
	
	if (trans == 'T')
		trans  = 'N' ;
	
	#if LAPACKE_TEST_VERBOSE
		printf(" \n unmqr dcomplex: matrix_layout = %d, side :%c, trans:%c, m:%d, n: %d \n", matrix_layout, side, trans, m, n);
	#endif
	if (side == 'L')
	{
		k = m;
		if (matrix_layout == LAPACK_COL_MAJOR)
		{		
			ldc = m;
			lda = k;
			lda_geqrf = m;
			bufsize_a = (lda_geqrf*n);
			bufsize_c = (ldc*n);			
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{		
			ldc = n;
			lda = m;
			lda_geqrf =n;
			bufsize_a = (lda_geqrf*k);
			bufsize_c = (ldc*m);
		}else
		{
			EXPECT_TRUE(false) << "matrix_layout invalid";
		}
	}else if (side == 'R') {
		k = n;
		if (matrix_layout == LAPACK_COL_MAJOR)
		{	ldc = m;
			lda = k;
			lda_geqrf = m;
			bufsize_a = (lda_geqrf*n);
			bufsize_c = (ldc*n);			
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{	ldc = n;
			lda = n;
			lda_geqrf =n;
			bufsize_a = (lda_geqrf*m);
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
		EXPECT_FALSE( true) << "unmqr_double_parameters object: malloc error.";
		unmqr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( A, Aref, bufsize_a);
	lapacke_gtest_init_dcomplex_buffer_pair_rand( c, cref, bufsize_c);


} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
unmqr_dcomplex_parameters :: ~unmqr_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" unmqr_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   unmqr_free();

}
/*  Test fixture class definition */
class zunmqr_test  : public  ::testing::Test {
public:
   unmqr_dcomplex_parameters  *zunmqr_obj;
   void SetUp();
   void TearDown () { delete zunmqr_obj;}
};

void zunmqr_test::SetUp(){

	/* LAPACKE zgeqrf prototype */
    typedef int (*Fptr_NL_LAPACKE_zgeqrf) (int matrix_layout, lapack_int m, lapack_int n,\
										   lapack_complex_double* a, lapack_int lda, lapack_complex_double* tau);
	 Fptr_NL_LAPACKE_zgeqrf zgeqrf;
	 /* LAPACKE zunmqr prototype */
    typedef int (*Fptr_NL_LAPACKE_zunmqr) (int matrix_layout, char side, char trans, lapack_int m, lapack_int n, \
											lapack_int k, const lapack_complex_double* a, lapack_int lda, \
											const lapack_complex_double* tau, lapack_complex_double* c, lapack_int ldc);

    Fptr_NL_LAPACKE_zunmqr zunmqr;

    zunmqr_obj = new unmqr_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].side,
						   eig_paramslist[idx].trans,
						   eig_paramslist[idx].m,
                           eig_paramslist[idx].n);
						   

    idx = Circular_Increment_Index(idx);

    zunmqr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zunmqr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zunmqr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zunmqr_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*zunmqr library call */
    zunmqr = (Fptr_NL_LAPACKE_zunmqr)dlsym(zunmqr_obj->hModule, "LAPACKE_zunmqr");
    ASSERT_TRUE(zunmqr != NULL) << "failed to get the Netlib LAPACKE_zunmqr symbol";

	/*zgeqrf library call*/
	zunmqr_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zunmqr_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zunmqr_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zunmqr_obj->lModule != NULL) << "Netlib lapacke handle NULL";

    zgeqrf = (Fptr_NL_LAPACKE_zgeqrf)dlsym(zunmqr_obj->lModule, "LAPACKE_zgeqrf");
    ASSERT_TRUE(zgeqrf != NULL) << "failed to get the Netlib LAPACKE_zgeqrf symbol";    
    

    zunmqr_obj->inforef_geqrf = zgeqrf( zunmqr_obj->matrix_layout,zunmqr_obj->m, zunmqr_obj->n,
								zunmqr_obj->Aref, zunmqr_obj->lda_geqrf, zunmqr_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    zunmqr_obj->info_geqrf = LAPACKE_zgeqrf( zunmqr_obj->matrix_layout,zunmqr_obj->m, zunmqr_obj->n,
											zunmqr_obj->A, zunmqr_obj->lda_geqrf, zunmqr_obj->tau);
										
										

    if( zunmqr_obj->info_geqrf < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zgeqrf is wrong\n", zunmqr_obj->info_geqrf );
    }
    if( zunmqr_obj->inforef_geqrf < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zgeqrf is wrong\n", 
        zunmqr_obj->inforef_geqrf );
    }  
/*Compute zunmqr's  o/p */
    zunmqr_obj->inforef = zunmqr( zunmqr_obj->matrix_layout, zunmqr_obj->side, zunmqr_obj->trans,
								zunmqr_obj->m, zunmqr_obj->n, zunmqr_obj->k, (const lapack_complex_double*)zunmqr_obj->Aref,
								zunmqr_obj->lda, (const lapack_complex_double*)zunmqr_obj->tauref, zunmqr_obj->cref, zunmqr_obj->ldc);

    /* Compute libflame's Lapacke o/p  */
    zunmqr_obj->info = LAPACKE_zunmqr( zunmqr_obj->matrix_layout, zunmqr_obj->side, zunmqr_obj->trans,
									zunmqr_obj->m, zunmqr_obj->n, zunmqr_obj->k, (const lapack_complex_double*)zunmqr_obj->A,
									zunmqr_obj->lda, (const lapack_complex_double*)zunmqr_obj->tau, zunmqr_obj->c, zunmqr_obj->ldc);
    if( zunmqr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame\
        LAPACKE_zunmqr is wrong\n", zunmqr_obj->info );
    }
    if( zunmqr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zunmqr is wrong\n", 
        zunmqr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zunmqr_obj->diff =  computeDiff_z( zunmqr_obj->bufsize_c, 
                zunmqr_obj->c, zunmqr_obj->cref );

}

TEST_F(zunmqr_test, zunmqr1) {
    EXPECT_NEAR(0.0, zunmqr_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zunmqr_test, zunmqr2) {
    EXPECT_NEAR(0.0, zunmqr_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zunmqr_test, zunmqr3) {
    EXPECT_NEAR(0.0, zunmqr_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zunmqr_test, zunmqr4) {
    EXPECT_NEAR(0.0, zunmqr_obj->diff, LAPACKE_EIG_THRESHOLD);
}
