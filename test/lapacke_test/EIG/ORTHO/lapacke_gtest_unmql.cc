#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define unmql_free() \
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
class unmql_scomplex_parameters{

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
	lapack_int ldc, lda, lda_geqlf;
	/*Output Parameter*/
	lapack_complex_float* tau, *c;	
	lapack_complex_float *Aref, *tauref, *cref;
	/*Return Values*/	
	lapack_int info, inforef;
	lapack_int info_geqlf, inforef_geqlf;

   public:
      unmql_scomplex_parameters (int matrix_layout, char side , char trans, lapack_int m, lapack_int n);
      ~unmql_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
unmql_scomplex_parameters:: unmql_scomplex_parameters (int matrix_layout_i, char side_i , char trans_i,lapack_int m_i, lapack_int n_i)
{
	
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	side = side_i;
	trans = trans_i;
	
	if (trans == 'T')
		trans  = 'N' ;
	
	#if LAPACKE_TEST_VERBOSE
		printf(" \n unmql scomplex: matrix_layout = %d, side :%c, trans:%c, m:%d, n: %d \n", matrix_layout, side, trans, m, n);
	#endif
	if (side == 'L')
	{
		k = m;
		if (matrix_layout == LAPACK_COL_MAJOR)
		{		
			ldc = m;
			lda = k;
			lda_geqlf = m;
			bufsize_a = (lda_geqlf*n);
			bufsize_c = (ldc*n);			
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{		
			ldc = n;
			lda = k;
			lda_geqlf =n;
			bufsize_a = (lda_geqlf*m);
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
			lda_geqlf = m;
			bufsize_a = (lda_geqlf*k);
			bufsize_c = (ldc*n);			
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{	ldc = n;
			lda = k;
			lda_geqlf = n;
			bufsize_a = (lda_geqlf*m);
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
		EXPECT_FALSE( true) << "unmql_float_parameters object: malloc error.";
		unmql_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( A, Aref, bufsize_a);
	lapacke_gtest_init_scomplex_buffer_pair_rand( c, cref, bufsize_c);
	

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
unmql_scomplex_parameters :: ~unmql_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" unmql_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   unmql_free();

}
/*  Test fixture class definition */
class cunmql_test  : public  ::testing::Test {
public:
   unmql_scomplex_parameters  *cunmql_obj;
   void SetUp();
   void TearDown () { delete cunmql_obj;}
};

void cunmql_test::SetUp(){

	/* LAPACKE cgeqlf prototype */
    typedef int (*Fptr_NL_LAPACKE_cgeqlf) (int matrix_layout, lapack_int m, lapack_int n,\
										   lapack_complex_float* a, lapack_int lda, lapack_complex_float* tau);
	 Fptr_NL_LAPACKE_cgeqlf cgeqlf;
	 /* LAPACKE cunmql prototype */
    typedef int (*Fptr_NL_LAPACKE_cunmql) (int matrix_layout, char side, char trans, lapack_int m, lapack_int n, \
											lapack_int k, const lapack_complex_float* a, lapack_int lda, \
											const lapack_complex_float* tau, lapack_complex_float* c, lapack_int ldc);

    Fptr_NL_LAPACKE_cunmql cunmql;

    cunmql_obj = new unmql_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].side,
						   eig_paramslist[idx].trans,
						   eig_paramslist[idx].m,
                           eig_paramslist[idx].n);
						   

    idx = Circular_Increment_Index(idx);

    cunmql_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cunmql_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cunmql_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cunmql_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*cunmql library call */
    cunmql = (Fptr_NL_LAPACKE_cunmql)dlsym(cunmql_obj->hModule, "LAPACKE_cunmql");
    ASSERT_TRUE(cunmql != NULL) << "failed to get the Netlib LAPACKE_cunmql symbol";

	/*cgeqlf library call*/
	cunmql_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cunmql_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cunmql_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cunmql_obj->lModule != NULL) << "Netlib lapacke handle NULL";

    cgeqlf = (Fptr_NL_LAPACKE_cgeqlf)dlsym(cunmql_obj->lModule, "LAPACKE_cgeqlf");
    ASSERT_TRUE(cgeqlf != NULL) << "failed to get the Netlib LAPACKE_cgeqlf symbol";    
    

    cunmql_obj->inforef_geqlf = cgeqlf( cunmql_obj->matrix_layout,cunmql_obj->m, cunmql_obj->n,
								cunmql_obj->Aref, cunmql_obj->lda_geqlf, cunmql_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    cunmql_obj->info_geqlf = LAPACKE_cgeqlf( cunmql_obj->matrix_layout,cunmql_obj->m, cunmql_obj->n,
											cunmql_obj->A, cunmql_obj->lda_geqlf, cunmql_obj->tau);
										
										

    if( cunmql_obj->info_geqlf < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cgeqlf is wrong\n", cunmql_obj->info_geqlf );
    }
    if( cunmql_obj->inforef_geqlf < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cgeqlf is wrong\n", 
        cunmql_obj->inforef_geqlf );
    }  
/*Compute cunmql's  o/p */
    cunmql_obj->inforef = cunmql( cunmql_obj->matrix_layout, cunmql_obj->side, cunmql_obj->trans,
								cunmql_obj->m, cunmql_obj->n, cunmql_obj->k, (const lapack_complex_float*)cunmql_obj->Aref,
								cunmql_obj->lda, (const lapack_complex_float*)cunmql_obj->tauref, cunmql_obj->cref, cunmql_obj->ldc);

    /* Compute libflame's Lapacke o/p  */
    cunmql_obj->info = LAPACKE_cunmql( cunmql_obj->matrix_layout, cunmql_obj->side, cunmql_obj->trans,
									cunmql_obj->m, cunmql_obj->n, cunmql_obj->k, (const lapack_complex_float*)cunmql_obj->A,
									cunmql_obj->lda, (const lapack_complex_float*)cunmql_obj->tau, cunmql_obj->c, cunmql_obj->ldc);
    if( cunmql_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame\
        LAPACKE_cunmql is wrong\n", cunmql_obj->info );
    }
    if( cunmql_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cunmql is wrong\n", 
        cunmql_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cunmql_obj->diff =  computeDiff_c( cunmql_obj->bufsize_c, 
                cunmql_obj->c, cunmql_obj->cref );

}

TEST_F(cunmql_test, cunmql1) {
    EXPECT_NEAR(0.0, cunmql_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(cunmql_test, cunmql2) {
    EXPECT_NEAR(0.0, cunmql_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(cunmql_test, cunmql3) {
    EXPECT_NEAR(0.0, cunmql_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(cunmql_test, cunmql4) {
    EXPECT_NEAR(0.0, cunmql_obj->diff, LAPACKE_EIG_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class unmql_dcomplex_parameters{

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
	lapack_int ldc, lda, lda_geqlf;
	/*Output Parameter*/
	lapack_complex_double* tau, *c;	
	lapack_complex_double *Aref, *tauref, *cref;
	/*Return Values*/	
	lapack_int info, inforef;
	lapack_int info_geqlf, inforef_geqlf;

   public:
      unmql_dcomplex_parameters (int matrix_layout, char side , char trans, lapack_int m, lapack_int n);
      ~unmql_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
unmql_dcomplex_parameters:: unmql_dcomplex_parameters (int matrix_layout_i, char side_i , char trans_i,lapack_int m_i, lapack_int n_i)
{
	
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	side = side_i;
	trans = trans_i;
	
	if (trans == 'T')
		trans  = 'N' ;
	
	#if LAPACKE_TEST_VERBOSE
		printf(" \n unmql dcomplex: matrix_layout = %d, side :%c, trans:%c, m:%d, n: %d \n", matrix_layout, side, trans, m, n);
	#endif
	if (side == 'L')
	{
		k = m;
		if (matrix_layout == LAPACK_COL_MAJOR)
		{		
			ldc = m;
			lda = k;
			lda_geqlf = m;
			bufsize_a = (lda_geqlf*n);
			bufsize_c = (ldc*n);			
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{		
			ldc = n;
			lda = k;
			lda_geqlf =n;
			bufsize_a = (lda_geqlf*m);
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
			lda_geqlf = m;
			bufsize_a = (lda_geqlf*k);
			bufsize_c = (ldc*n);			
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{	ldc = n;
			lda = k;
			lda_geqlf = n;
			bufsize_a = (lda_geqlf*m);
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
		EXPECT_FALSE( true) << "unmql_double_parameters object: malloc error.";
		unmql_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( A, Aref, bufsize_a);
	lapacke_gtest_init_dcomplex_buffer_pair_rand( c, cref, bufsize_c);
	

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
unmql_dcomplex_parameters :: ~unmql_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" unmql_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   unmql_free();

}
/*  Test fixture class definition */
class zunmql_test  : public  ::testing::Test {
public:
   unmql_dcomplex_parameters  *zunmql_obj;
   void SetUp();
   void TearDown () { delete zunmql_obj;}
};

void zunmql_test::SetUp(){

	/* LAPACKE zgeqlf prototype */
    typedef int (*Fptr_NL_LAPACKE_zgeqlf) (int matrix_layout, lapack_int m, lapack_int n,\
										   lapack_complex_double* a, lapack_int lda, lapack_complex_double* tau);
	 Fptr_NL_LAPACKE_zgeqlf zgeqlf;
	 /* LAPACKE zunmql prototype */
    typedef int (*Fptr_NL_LAPACKE_zunmql) (int matrix_layout, char side, char trans, lapack_int m, lapack_int n, \
											lapack_int k, const lapack_complex_double* a, lapack_int lda, \
											const lapack_complex_double* tau, lapack_complex_double* c, lapack_int ldc);

    Fptr_NL_LAPACKE_zunmql zunmql;

    zunmql_obj = new unmql_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].side,
						   eig_paramslist[idx].trans,
						   eig_paramslist[idx].m,
                           eig_paramslist[idx].n);
						   

    idx = Circular_Increment_Index(idx);

    zunmql_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zunmql_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zunmql_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zunmql_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*zunmql library call */
    zunmql = (Fptr_NL_LAPACKE_zunmql)dlsym(zunmql_obj->hModule, "LAPACKE_zunmql");
    ASSERT_TRUE(zunmql != NULL) << "failed to get the Netlib LAPACKE_zunmql symbol";

	/*zgeqlf library call*/
	zunmql_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zunmql_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zunmql_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zunmql_obj->lModule != NULL) << "Netlib lapacke handle NULL";

    zgeqlf = (Fptr_NL_LAPACKE_zgeqlf)dlsym(zunmql_obj->lModule, "LAPACKE_zgeqlf");
    ASSERT_TRUE(zgeqlf != NULL) << "failed to get the Netlib LAPACKE_zgeqlf symbol";    
    

    zunmql_obj->inforef_geqlf = zgeqlf( zunmql_obj->matrix_layout,zunmql_obj->m, zunmql_obj->n,
								zunmql_obj->Aref, zunmql_obj->lda_geqlf, zunmql_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    zunmql_obj->info_geqlf = LAPACKE_zgeqlf( zunmql_obj->matrix_layout,zunmql_obj->m, zunmql_obj->n,
											zunmql_obj->A, zunmql_obj->lda_geqlf, zunmql_obj->tau);
										
										

    if( zunmql_obj->info_geqlf < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zgeqlf is wrong\n", zunmql_obj->info_geqlf );
    }
    if( zunmql_obj->inforef_geqlf < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zgeqlf is wrong\n", 
        zunmql_obj->inforef_geqlf );
    }  
/*Compute zunmql's  o/p */
    zunmql_obj->inforef = zunmql( zunmql_obj->matrix_layout, zunmql_obj->side, zunmql_obj->trans,
								zunmql_obj->m, zunmql_obj->n, zunmql_obj->k, (const lapack_complex_double*)zunmql_obj->Aref,
								zunmql_obj->lda, (const lapack_complex_double*)zunmql_obj->tauref, zunmql_obj->cref, zunmql_obj->ldc);

    /* Compute libflame's Lapacke o/p  */
    zunmql_obj->info = LAPACKE_zunmql( zunmql_obj->matrix_layout, zunmql_obj->side, zunmql_obj->trans,
									zunmql_obj->m, zunmql_obj->n, zunmql_obj->k, (const lapack_complex_double*)zunmql_obj->A,
									zunmql_obj->lda, (const lapack_complex_double*)zunmql_obj->tau, zunmql_obj->c, zunmql_obj->ldc);
    if( zunmql_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame\
        LAPACKE_zunmql is wrong\n", zunmql_obj->info );
    }
    if( zunmql_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zunmql is wrong\n", 
        zunmql_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zunmql_obj->diff =  computeDiff_z( zunmql_obj->bufsize_c, 
                zunmql_obj->c, zunmql_obj->cref );

}

TEST_F(zunmql_test, zunmql1) {
    EXPECT_NEAR(0.0, zunmql_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zunmql_test, zunmql2) {
    EXPECT_NEAR(0.0, zunmql_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zunmql_test, zunmql3) {
    EXPECT_NEAR(0.0, zunmql_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zunmql_test, zunmql4) {
    EXPECT_NEAR(0.0, zunmql_obj->diff, LAPACKE_EIG_THRESHOLD);
}