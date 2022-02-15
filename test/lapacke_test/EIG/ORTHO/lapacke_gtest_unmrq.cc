#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"

#define unmrq_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (tau!=NULL)  free(tau);\
if (tauref!=NULL) free(tauref); \
if (c != NULL) free(c); \
if (cref != NULL) free(cref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule);\
if( bModule != NULL) dlclose(bModule); \
if( lModule != NULL) dlclose(lModule); \
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
/* Begin scomplex_common_parameters  class definition */
class unmrq_scomplex_parameters{

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
	lapack_int ldc, lda, lda_gerqf;
	/*Output Parameter*/
	lapack_complex_float* tau, *c;
	lapack_complex_float *Aref, *tauref, *cref;
	/*Return Values*/
	lapack_int info, inforef;
	lapack_int info_gerqf, inforef_gerqf;

   public:
      unmrq_scomplex_parameters (int matrix_layout, char side , char trans, lapack_int m, lapack_int n);
      ~unmrq_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
unmrq_scomplex_parameters:: unmrq_scomplex_parameters (int matrix_layout_i, char side_i , char trans_i,lapack_int m_i, lapack_int n_i)
{
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	side = side_i;
	trans = trans_i;
	
	if (trans == 'T')
		trans = 'N';
	#if LAPACKE_TEST_VERBOSE
		printf(" \n unmrq scomplex: matrix_layout = %d, side :%c, trans:%c, m:%d, n: %d \n", matrix_layout, side, trans, m, n);
	#endif
	if (side == 'L')
	{
		k = m;
		if (matrix_layout == LAPACK_COL_MAJOR)
		{		
			ldc = m;
			lda = k;
			lda_gerqf = m;
			bufsize_a = (lda_gerqf*n);
			bufsize_c = (ldc*n);			
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{		
			ldc = n;
			lda = m;
			lda_gerqf =n;
			bufsize_a = (lda_gerqf*k);
			bufsize_c = (ldc*m);
		}else
		{
			EXPECT_TRUE(false) << "matrix_layout invalid";
		}
	}else if (side == 'R') {
		k = n;
		if (matrix_layout == LAPACK_COL_MAJOR)
		{	ldc = m;
			lda = m;			
			bufsize_c = (ldc*n);
			lda_gerqf = m;
			bufsize_a = (lda_gerqf*n);
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{	ldc = n;
			lda = m;
			lda_gerqf =n;
			bufsize_a = (lda_gerqf*m);
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
		EXPECT_FALSE( true) << "unmrq_float_parameters object: malloc error.";
		unmrq_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( A, Aref, bufsize_a);
	lapacke_gtest_init_scomplex_buffer_pair_rand( c, cref, bufsize_c);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
unmrq_scomplex_parameters :: ~unmrq_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" unmrq_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   unmrq_free();

}
/*  Test fixture class definition */
class cunmrq_test  : public  ::testing::Test {
public:
   unmrq_scomplex_parameters  *cunmrq_obj;
   void SetUp();
   void TearDown () { delete cunmrq_obj;}
};

void cunmrq_test::SetUp(){

	/* LAPACKE cgerqf prototype */
    typedef int (*Fptr_NL_LAPACKE_cgerqf) (int matrix_layout, lapack_int m, lapack_int n,\
										   lapack_complex_float* a, lapack_int lda, lapack_complex_float* tau);
	 Fptr_NL_LAPACKE_cgerqf cgerqf;
	 /* LAPACKE cunmrq prototype */
    typedef int (*Fptr_NL_LAPACKE_cunmrq) (int matrix_layout, char side, char trans, lapack_int m, lapack_int n, \
											lapack_int k, const lapack_complex_float* a, lapack_int lda, \
											const lapack_complex_float* tau, lapack_complex_float* c, lapack_int ldc);

    Fptr_NL_LAPACKE_cunmrq cunmrq;
	
    cunmrq_obj = new unmrq_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].side,
						   eig_paramslist[idx].trans,
						   eig_paramslist[idx].m,
                           eig_paramslist[idx].n);
						   

    idx = Circular_Increment_Index(idx);

    cunmrq_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cunmrq_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cunmrq_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cunmrq_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*cunmrq library call */
    cunmrq = (Fptr_NL_LAPACKE_cunmrq)dlsym(cunmrq_obj->hModule, "LAPACKE_cunmrq");
    ASSERT_TRUE(cunmrq != NULL) << "failed to get the Netlib LAPACKE_cunmrq symbol";

	/*cgerqf library call*/
	cunmrq_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cunmrq_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cunmrq_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cunmrq_obj->lModule != NULL) << "Netlib lapacke handle NULL";

    cgerqf = (Fptr_NL_LAPACKE_cgerqf)dlsym(cunmrq_obj->lModule, "LAPACKE_cgerqf");
    ASSERT_TRUE(cgerqf != NULL) << "failed to get the Netlib LAPACKE_cgerqf symbol";    
    

    cunmrq_obj->inforef_gerqf = cgerqf( cunmrq_obj->matrix_layout,cunmrq_obj->m, cunmrq_obj->n,
								cunmrq_obj->Aref, cunmrq_obj->lda_gerqf, cunmrq_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    cunmrq_obj->info_gerqf = LAPACKE_cgerqf( cunmrq_obj->matrix_layout,cunmrq_obj->m, cunmrq_obj->n,
											cunmrq_obj->A, cunmrq_obj->lda_gerqf, cunmrq_obj->tau);
										
										

    if( cunmrq_obj->info_gerqf < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cgerqf is wrong\n", cunmrq_obj->info_gerqf );
    }
    if( cunmrq_obj->inforef_gerqf < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cgerqf is wrong\n", 
        cunmrq_obj->inforef_gerqf );
    }  
/*Compute cunmrq's  o/p */
    cunmrq_obj->inforef = cunmrq( cunmrq_obj->matrix_layout, cunmrq_obj->side, cunmrq_obj->trans,
								cunmrq_obj->m, cunmrq_obj->n, cunmrq_obj->k, (const lapack_complex_float*)cunmrq_obj->Aref,
								cunmrq_obj->lda, (const lapack_complex_float*)cunmrq_obj->tauref, cunmrq_obj->cref, cunmrq_obj->ldc);

    /* Compute libflame's Lapacke o/p  */
    cunmrq_obj->info = LAPACKE_cunmrq( cunmrq_obj->matrix_layout, cunmrq_obj->side, cunmrq_obj->trans,
									cunmrq_obj->m, cunmrq_obj->n, cunmrq_obj->k, (const lapack_complex_float*)cunmrq_obj->A,
									cunmrq_obj->lda, (const lapack_complex_float*)cunmrq_obj->tau, cunmrq_obj->c, cunmrq_obj->ldc);
    if( cunmrq_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame\
        LAPACKE_cunmrq is wrong\n", cunmrq_obj->info );
    }
    if( cunmrq_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cunmrq is wrong\n", 
        cunmrq_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cunmrq_obj->diff =  computeDiff_c( cunmrq_obj->bufsize_c, 
                cunmrq_obj->c, cunmrq_obj->cref );

}

TEST_F(cunmrq_test, cunmrq1) {
    EXPECT_NEAR(0.0, cunmrq_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(cunmrq_test, cunmrq2) {
    EXPECT_NEAR(0.0, cunmrq_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(cunmrq_test, cunmrq3) {
    EXPECT_NEAR(0.0, cunmrq_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(cunmrq_test, cunmrq4) {
    EXPECT_NEAR(0.0, cunmrq_obj->diff, LAPACKE_EIG_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class unmrq_dcomplex_parameters{

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
	lapack_int ldc, lda, lda_gerqf;
	/*Output Parameter*/
	lapack_complex_double* tau, *c;
	lapack_complex_double *Aref, *tauref, *cref;
	/*Return Values*/
	lapack_int info, inforef;
	lapack_int info_gerqf, inforef_gerqf;

   public:
      unmrq_dcomplex_parameters (int matrix_layout, char side , char trans, lapack_int m, lapack_int n);
      ~unmrq_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
unmrq_dcomplex_parameters:: unmrq_dcomplex_parameters (int matrix_layout_i, char side_i , char trans_i,lapack_int m_i, lapack_int n_i)
{
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	side = side_i;
	trans = trans_i;
	
	if (trans == 'T')
		trans = 'N';
	#if LAPACKE_TEST_VERBOSE
		printf(" \n unmrq dcomplex: matrix_layout = %d, side :%c, trans:%c, m:%d, n: %d \n", matrix_layout, side, trans, m, n);
	#endif
	if (side == 'L')
	{
		k = m;
		if (matrix_layout == LAPACK_COL_MAJOR)
		{		
			ldc = m;
			lda = k;
			lda_gerqf = m;
			bufsize_a = (lda_gerqf*n);
			bufsize_c = (ldc*n);			
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{		
			ldc = n;
			lda = m;
			lda_gerqf =n;
			bufsize_a = (lda_gerqf*k);
			bufsize_c = (ldc*m);
		}else
		{
			EXPECT_TRUE(false) << "matrix_layout invalid";
		}
	}else if (side == 'R') {
		k = n;
		if (matrix_layout == LAPACK_COL_MAJOR)
		{	ldc = m;
			lda = m;			
			bufsize_c = (ldc*n);
			lda_gerqf = m;
			bufsize_a = (lda_gerqf*n);
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{	ldc = n;
			lda = m;
			lda_gerqf =n;
			bufsize_a = (lda_gerqf*m);
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
		EXPECT_FALSE( true) << "unmrq_double_parameters object: malloc error.";
		unmrq_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( A, Aref, bufsize_a);
	lapacke_gtest_init_dcomplex_buffer_pair_rand( c, cref, bufsize_c);

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
unmrq_dcomplex_parameters :: ~unmrq_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" unmrq_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   unmrq_free();

}
/*  Test fixture class definition */
class zunmrq_test  : public  ::testing::Test {
public:
   unmrq_dcomplex_parameters  *zunmrq_obj;
   void SetUp();
   void TearDown () { delete zunmrq_obj;}
};

void zunmrq_test::SetUp(){

	/* LAPACKE zgerqf prototype */
    typedef int (*Fptr_NL_LAPACKE_zgerqf) (int matrix_layout, lapack_int m, lapack_int n,\
										   lapack_complex_double* a, lapack_int lda, lapack_complex_double* tau);
	 Fptr_NL_LAPACKE_zgerqf zgerqf;
	 /* LAPACKE zunmrq prototype */
    typedef int (*Fptr_NL_LAPACKE_zunmrq) (int matrix_layout, char side, char trans, lapack_int m, lapack_int n, \
											lapack_int k, const lapack_complex_double* a, lapack_int lda, \
											const lapack_complex_double* tau, lapack_complex_double* c, lapack_int ldc);

    Fptr_NL_LAPACKE_zunmrq zunmrq;
	
    zunmrq_obj = new unmrq_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].side,
						   eig_paramslist[idx].trans,
						   eig_paramslist[idx].m,
                           eig_paramslist[idx].n);
						   

    idx = Circular_Increment_Index(idx);

    zunmrq_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zunmrq_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zunmrq_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zunmrq_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*zunmrq library call */
    zunmrq = (Fptr_NL_LAPACKE_zunmrq)dlsym(zunmrq_obj->hModule, "LAPACKE_zunmrq");
    ASSERT_TRUE(zunmrq != NULL) << "failed to get the Netlib LAPACKE_zunmrq symbol";

	/*zgerqf library call*/
	zunmrq_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zunmrq_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zunmrq_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zunmrq_obj->lModule != NULL) << "Netlib lapacke handle NULL";

    zgerqf = (Fptr_NL_LAPACKE_zgerqf)dlsym(zunmrq_obj->lModule, "LAPACKE_zgerqf");
    ASSERT_TRUE(zgerqf != NULL) << "failed to get the Netlib LAPACKE_zgerqf symbol";    
    

    zunmrq_obj->inforef_gerqf = zgerqf( zunmrq_obj->matrix_layout,zunmrq_obj->m, zunmrq_obj->n,
								zunmrq_obj->Aref, zunmrq_obj->lda_gerqf, zunmrq_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    zunmrq_obj->info_gerqf = LAPACKE_zgerqf( zunmrq_obj->matrix_layout,zunmrq_obj->m, zunmrq_obj->n,
											zunmrq_obj->A, zunmrq_obj->lda_gerqf, zunmrq_obj->tau);
										
										

    if( zunmrq_obj->info_gerqf < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zgerqf is wrong\n", zunmrq_obj->info_gerqf );
    }
    if( zunmrq_obj->inforef_gerqf < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zgerqf is wrong\n", 
        zunmrq_obj->inforef_gerqf );
    }  
/*Compute zunmrq's  o/p */
    zunmrq_obj->inforef = zunmrq( zunmrq_obj->matrix_layout, zunmrq_obj->side, zunmrq_obj->trans,
								zunmrq_obj->m, zunmrq_obj->n, zunmrq_obj->k, (const lapack_complex_double*)zunmrq_obj->Aref,
								zunmrq_obj->lda, (const lapack_complex_double*)zunmrq_obj->tauref, zunmrq_obj->cref, zunmrq_obj->ldc);

    /* Compute libflame's Lapacke o/p  */
    zunmrq_obj->info = LAPACKE_zunmrq( zunmrq_obj->matrix_layout, zunmrq_obj->side, zunmrq_obj->trans,
									zunmrq_obj->m, zunmrq_obj->n, zunmrq_obj->k, (const lapack_complex_double*)zunmrq_obj->A,
									zunmrq_obj->lda, (const lapack_complex_double*)zunmrq_obj->tau, zunmrq_obj->c, zunmrq_obj->ldc);
    if( zunmrq_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame\
        LAPACKE_zunmrq is wrong\n", zunmrq_obj->info );
    }
    if( zunmrq_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zunmrq is wrong\n", 
        zunmrq_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zunmrq_obj->diff =  computeDiff_z( zunmrq_obj->bufsize_c, 
                zunmrq_obj->c, zunmrq_obj->cref );

}

TEST_F(zunmrq_test, zunmrq1) {
    EXPECT_NEAR(0.0, zunmrq_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zunmrq_test, zunmrq2) {
    EXPECT_NEAR(0.0, zunmrq_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zunmrq_test, zunmrq3) {
    EXPECT_NEAR(0.0, zunmrq_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zunmrq_test, zunmrq4) {
    EXPECT_NEAR(0.0, zunmrq_obj->diff, LAPACKE_EIG_THRESHOLD);
}