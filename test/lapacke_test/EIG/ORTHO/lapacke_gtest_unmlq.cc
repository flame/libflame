#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"

#define unmlq_free() \
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
class unmlq_scomplex_parameters{

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
	lapack_int ldc, lda, lda_gelqf;
	/*Output Parameter*/
	lapack_complex_float* tau, *c;	
	lapack_complex_float *Aref, *tauref, *cref;
	/*Return Values*/	
	lapack_int info, inforef;
	lapack_int info_gelqf, inforef_gelqf;

   public:
      unmlq_scomplex_parameters (int matrix_layout, char side , char trans, lapack_int m, lapack_int n);
      ~unmlq_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
unmlq_scomplex_parameters:: unmlq_scomplex_parameters (int matrix_layout_i, char side_i , char trans_i,lapack_int m_i, lapack_int n_i)
{
	
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	side = side_i;
	trans = trans_i;
	
	if (trans == 'T')
		trans  = 'N' ;
	
	#if LAPACKE_TEST_VERBOSE
		printf(" \n unmlq scomplex: matrix_layout = %d, side :%c, trans:%c, m:%d, n: %d \n", matrix_layout, side, trans, m, n);
	#endif
	if (side == 'L')
	{
		k = m;
		if (matrix_layout == LAPACK_COL_MAJOR)
		{		
			ldc = m;
			lda = k;
			lda_gelqf = m;
			bufsize_a = (lda_gelqf*n);
			bufsize_c = (ldc*n);			
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{		
			ldc = n;
			lda = k;
			lda_gelqf =n;
			bufsize_a = (lda_gelqf*m);
			bufsize_c = (ldc*m);
		}else
		{
			EXPECT_TRUE(false) << "matrix_layout invalid";
		}
	}else if (side == 'R') {
		k = n;
		if (matrix_layout == LAPACK_COL_MAJOR)
		{	ldc = m;
			lda = n;
			lda_gelqf = m;
			bufsize_a = (lda_gelqf*k);
			bufsize_c = (ldc*n);			
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{	ldc = n;
			lda = k;
			lda_gelqf = n;
			bufsize_a = (lda_gelqf*m);
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
		EXPECT_FALSE( true) << "unmlq_float_parameters object: malloc error.";
		unmlq_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( A, Aref, bufsize_a);
	lapacke_gtest_init_scomplex_buffer_pair_rand( c, cref, bufsize_c);
	

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
unmlq_scomplex_parameters :: ~unmlq_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" unmlq_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   unmlq_free();

}
/*  Test fixture class definition */
class cunmlq_test  : public  ::testing::Test {
public:
   unmlq_scomplex_parameters  *cunmlq_obj;
   void SetUp();
   void TearDown () { delete cunmlq_obj;}
};

void cunmlq_test::SetUp(){

	/* LAPACKE cgelqf prototype */
    typedef int (*Fptr_NL_LAPACKE_cgelqf) (int matrix_layout, lapack_int m, lapack_int n, lapack_complex_float* a, lapack_int lda, lapack_complex_float* tau);
	 Fptr_NL_LAPACKE_cgelqf cgelqf;
	 /* LAPACKE cunmlq prototype */
    typedef int (*Fptr_NL_LAPACKE_cunmlq) (int matrix_layout, char side, char trans, lapack_int m, lapack_int n, lapack_int k,\
											const lapack_complex_float* a, lapack_int lda, const lapack_complex_float* tau, lapack_complex_float* c, lapack_int ldc);

    Fptr_NL_LAPACKE_cunmlq cunmlq;

    cunmlq_obj = new unmlq_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].side,
						   eig_paramslist[idx].trans,
						   eig_paramslist[idx].m,
                           eig_paramslist[idx].n);
						   

    idx = Circular_Increment_Index(idx);

    cunmlq_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cunmlq_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cunmlq_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cunmlq_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*cunmlq library call */
    cunmlq = (Fptr_NL_LAPACKE_cunmlq)dlsym(cunmlq_obj->hModule, "LAPACKE_cunmlq");
    ASSERT_TRUE(cunmlq != NULL) << "failed to get the Netlib LAPACKE_cunmlq symbol";

	/*cgelqf library call*/
	cunmlq_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cunmlq_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cunmlq_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cunmlq_obj->lModule != NULL) << "Netlib lapacke handle NULL";

    cgelqf = (Fptr_NL_LAPACKE_cgelqf)dlsym(cunmlq_obj->lModule, "LAPACKE_cgelqf");
    ASSERT_TRUE(cgelqf != NULL) << "failed to get the Netlib LAPACKE_cgelqf symbol";    
    

    cunmlq_obj->inforef_gelqf = cgelqf( cunmlq_obj->matrix_layout,cunmlq_obj->m, cunmlq_obj->n,
								cunmlq_obj->Aref, cunmlq_obj->lda_gelqf, cunmlq_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    cunmlq_obj->info_gelqf = LAPACKE_cgelqf( cunmlq_obj->matrix_layout,cunmlq_obj->m, cunmlq_obj->n,
											cunmlq_obj->A, cunmlq_obj->lda_gelqf, cunmlq_obj->tau);
										
										

    if( cunmlq_obj->info_gelqf < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cgelqf is wrong\n", cunmlq_obj->info_gelqf );
    }
    if( cunmlq_obj->inforef_gelqf < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cgelqf is wrong\n", 
        cunmlq_obj->inforef_gelqf );
    }  
/*Compute cunmlq's  o/p */
    cunmlq_obj->inforef = cunmlq( cunmlq_obj->matrix_layout, cunmlq_obj->side, cunmlq_obj->trans,
								cunmlq_obj->m, cunmlq_obj->n, cunmlq_obj->k, (const lapack_complex_float*)cunmlq_obj->Aref,
								cunmlq_obj->lda, (const lapack_complex_float*)cunmlq_obj->tauref, cunmlq_obj->cref, cunmlq_obj->ldc);

    /* Compute libflame's Lapacke o/p  */
    cunmlq_obj->info = LAPACKE_cunmlq( cunmlq_obj->matrix_layout, cunmlq_obj->side, cunmlq_obj->trans,
									cunmlq_obj->m, cunmlq_obj->n, cunmlq_obj->k, (const lapack_complex_float*)cunmlq_obj->A,
									cunmlq_obj->lda, (const lapack_complex_float*)cunmlq_obj->tau, cunmlq_obj->c, cunmlq_obj->ldc);
    if( cunmlq_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame\
        LAPACKE_cunmlq is wrong\n", cunmlq_obj->info );
    }
    if( cunmlq_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cunmlq is wrong\n", 
        cunmlq_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cunmlq_obj->diff =  computeDiff_c( cunmlq_obj->bufsize_c, 
                cunmlq_obj->c, cunmlq_obj->cref );

}

TEST_F(cunmlq_test, cunmlq1) {
    EXPECT_NEAR(0.0, cunmlq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cunmlq_test, cunmlq2) {
    EXPECT_NEAR(0.0, cunmlq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cunmlq_test, cunmlq3) {
    EXPECT_NEAR(0.0, cunmlq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cunmlq_test, cunmlq4) {
    EXPECT_NEAR(0.0, cunmlq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class unmlq_dcomplex_parameters{

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
	lapack_int ldc, lda, lda_gelqf;
	/*Output Parameter*/
	lapack_complex_double* tau, *c;	
	lapack_complex_double *Aref, *tauref, *cref;
	/*Return Values*/	
	lapack_int info, inforef;
	lapack_int info_gelqf, inforef_gelqf;

   public:
      unmlq_dcomplex_parameters (int matrix_layout, char side , char trans, lapack_int m, lapack_int n);
      ~unmlq_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
unmlq_dcomplex_parameters:: unmlq_dcomplex_parameters (int matrix_layout_i, char side_i , char trans_i,lapack_int m_i, lapack_int n_i)
{
	
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	side = side_i;
	trans = trans_i;
	
	if (trans == 'T')
		trans  = 'N' ;
	
	#if LAPACKE_TEST_VERBOSE
		printf(" \n unmlq dcomplex: matrix_layout = %d, side :%c, trans:%c, m:%d, n: %d \n", matrix_layout, side, trans, m, n);
	#endif
	if (side == 'L')
	{
		k = m;
		if (matrix_layout == LAPACK_COL_MAJOR)
		{		
			ldc = m;
			lda = k;
			lda_gelqf = m;
			bufsize_a = (lda_gelqf*n);
			bufsize_c = (ldc*n);			
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{		
			ldc = n;
			lda = k;
			lda_gelqf =n;
			bufsize_a = (lda_gelqf*m);
			bufsize_c = (ldc*m);
		}else
		{
			EXPECT_TRUE(false) << "matrix_layout invalid";
		}
	}else if (side == 'R') {
		k = n;
		if (matrix_layout == LAPACK_COL_MAJOR)
		{	ldc = m;
			lda = n;
			lda_gelqf = m;
			bufsize_a = (lda_gelqf*k);
			bufsize_c = (ldc*n);			
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{	ldc = n;
			lda = k;
			lda_gelqf = n;
			bufsize_a = (lda_gelqf*m);
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
		EXPECT_FALSE( true) << "unmlq_double_parameters object: malloc error.";
		unmlq_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( A, Aref, bufsize_a);
	lapacke_gtest_init_dcomplex_buffer_pair_rand( c, cref, bufsize_c);
	

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
unmlq_dcomplex_parameters :: ~unmlq_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" unmlq_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   unmlq_free();

}
/*  Test fixture class definition */
class zunmlq_test  : public  ::testing::Test {
public:
   unmlq_dcomplex_parameters  *zunmlq_obj;
   void SetUp();
   void TearDown () { delete zunmlq_obj;}
};

void zunmlq_test::SetUp(){

	/* LAPACKE zgelqf prototype */
    typedef int (*Fptr_NL_LAPACKE_zgelqf) (int matrix_layout, lapack_int m, lapack_int n, lapack_complex_double* a, lapack_int lda, lapack_complex_double* tau);
	 Fptr_NL_LAPACKE_zgelqf zgelqf;
	 /* LAPACKE zunmlq prototype */
    typedef int (*Fptr_NL_LAPACKE_zunmlq) (int matrix_layout, char side, char trans, lapack_int m, lapack_int n, lapack_int k,\
											const lapack_complex_double* a, lapack_int lda, const lapack_complex_double* tau, lapack_complex_double* c, lapack_int ldc);

    Fptr_NL_LAPACKE_zunmlq zunmlq;

    zunmlq_obj = new unmlq_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].side,
						   eig_paramslist[idx].trans,
						   eig_paramslist[idx].m,
                           eig_paramslist[idx].n);
						   

    idx = Circular_Increment_Index(idx);

    zunmlq_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zunmlq_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zunmlq_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zunmlq_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*zunmlq library call */
    zunmlq = (Fptr_NL_LAPACKE_zunmlq)dlsym(zunmlq_obj->hModule, "LAPACKE_zunmlq");
    ASSERT_TRUE(zunmlq != NULL) << "failed to get the Netlib LAPACKE_zunmlq symbol";

	/*zgelqf library call*/
	zunmlq_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zunmlq_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zunmlq_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zunmlq_obj->lModule != NULL) << "Netlib lapacke handle NULL";

    zgelqf = (Fptr_NL_LAPACKE_zgelqf)dlsym(zunmlq_obj->lModule, "LAPACKE_zgelqf");
    ASSERT_TRUE(zgelqf != NULL) << "failed to get the Netlib LAPACKE_zgelqf symbol";    
    

    zunmlq_obj->inforef_gelqf = zgelqf( zunmlq_obj->matrix_layout,zunmlq_obj->m, zunmlq_obj->n,
								zunmlq_obj->Aref, zunmlq_obj->lda_gelqf, zunmlq_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    zunmlq_obj->info_gelqf = LAPACKE_zgelqf( zunmlq_obj->matrix_layout,zunmlq_obj->m, zunmlq_obj->n,
											zunmlq_obj->A, zunmlq_obj->lda_gelqf, zunmlq_obj->tau);
										
										

    if( zunmlq_obj->info_gelqf < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zgelqf is wrong\n", zunmlq_obj->info_gelqf );
    }
    if( zunmlq_obj->inforef_gelqf < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zgelqf is wrong\n", 
        zunmlq_obj->inforef_gelqf );
    }  
/*Compute zunmlq's  o/p */
    zunmlq_obj->inforef = zunmlq( zunmlq_obj->matrix_layout, zunmlq_obj->side, zunmlq_obj->trans,
								zunmlq_obj->m, zunmlq_obj->n, zunmlq_obj->k, (const lapack_complex_double*)zunmlq_obj->Aref,
								zunmlq_obj->lda, (const lapack_complex_double*)zunmlq_obj->tauref, zunmlq_obj->cref, zunmlq_obj->ldc);

    /* Compute libflame's Lapacke o/p  */
    zunmlq_obj->info = LAPACKE_zunmlq( zunmlq_obj->matrix_layout, zunmlq_obj->side, zunmlq_obj->trans,
									zunmlq_obj->m, zunmlq_obj->n, zunmlq_obj->k, (const lapack_complex_double*)zunmlq_obj->A,
									zunmlq_obj->lda, (const lapack_complex_double*)zunmlq_obj->tau, zunmlq_obj->c, zunmlq_obj->ldc);
    if( zunmlq_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame\
        LAPACKE_zunmlq is wrong\n", zunmlq_obj->info );
    }
    if( zunmlq_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zunmlq is wrong\n", 
        zunmlq_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zunmlq_obj->diff =  computeDiff_z( zunmlq_obj->bufsize_c, 
                zunmlq_obj->c, zunmlq_obj->cref );

}

TEST_F(zunmlq_test, zunmlq1) {
    EXPECT_NEAR(0.0, zunmlq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zunmlq_test, zunmlq2) {
    EXPECT_NEAR(0.0, zunmlq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zunmlq_test, zunmlq3) {
    EXPECT_NEAR(0.0, zunmlq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zunmlq_test, zunmlq4) {
    EXPECT_NEAR(0.0, zunmlq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}