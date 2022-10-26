#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"




#define unmtr_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (tau!=NULL)  free(tau);\
if (tauref!=NULL) free(tauref); \
if (d != NULL) free(d); \
if (e != NULL) free(e); \
if (dref!=NULL) free(dref);\
if (eref!=NULL) free(eref);\
if (c != NULL) free(c); \
if (cref != NULL) free(cref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin scomplex_common_parameters  class definition */
class unmtr_scomplex_parameters{

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
	lapack_int r;
	lapack_complex_float* A;
	char uplo;
	char side;
	char trans;
	lapack_int ldc;
	lapack_int lda;
	/*Output Parameter*/
	lapack_complex_float* tau, *c;
	float* d;
	float* e;
	lapack_complex_float *Aref, *tauref, *cref;
	float* dref;
	float* eref;
	/*Return Values*/	
	lapack_int info, inforef;
	lapack_int info_hetrd, inforef_hetrd;

   public:
      unmtr_scomplex_parameters (int matrix_layout,  char uplo, lapack_int m, lapack_int n);
      ~unmtr_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
unmtr_scomplex_parameters:: unmtr_scomplex_parameters (int matrix_layout_i, char uplo_i, lapack_int m_i, lapack_int n_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	uplo = uplo_i;
	trans = 'N';
	#if LAPACKE_TEST_VERBOSE
	printf(" \n unmtr scomplex: matrix_layout = %d, uplo:%c,m:%d, n: %d \n", matrix_layout, uplo, m, n);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{	ldc = m;
		r = n;		
		side = 'R';		
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{	ldc = n;
		r = m;
		side = 'L';
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	lda = r;
	bufsize_c = ldc*r;
	bufsize_a =  lda *r;
	
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&tau, &tauref, (r-1));
	lapacke_gtest_alloc_float_buffer_pair(&d, &dref, r);
	lapacke_gtest_alloc_float_buffer_pair(&e, &eref, (r-1));
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&c, &cref, bufsize_c);
	if ((A==NULL) || (Aref==NULL) || \
		(tau==NULL) || (tauref==NULL) ||
		(d == NULL) || (dref == NULL) ||
		(e == NULL) || (eref == NULL) ||
		(c == NULL) || (cref == NULL)){
		EXPECT_FALSE( true) << "unmtr_float_parameters object: malloc error.";
		unmtr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( A, Aref, bufsize_a);
	lapacke_gtest_init_scomplex_buffer_pair_rand( c, cref, bufsize_c);
	lapacke_gtest_init_scomplex_buffer_pair_rand( tau, tauref, (r-1));
	/*initialize output matrix by 0 *
	for(i=0;i<(r-1);i++) {
		tau[i] = 0;
		tauref[i] = tau[i];
	} */
	/*initialize output matrix by 0 */
	for(i=0;i<(r-1);i++) {
		e[i] = 0;
		eref[i] = e[i];
	}
	/*initialize output matrix by 0 */
	for(i=0;i<r;i++) {
		d[i] = 0;
		dref[i] = d[i];
	}
	/*for(i=0;i<bufsize_c;i++) {
		c[i] = 0.0;
		cref[i] = c[i];
	}*/

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
unmtr_scomplex_parameters :: ~unmtr_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" unmtr_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   unmtr_free();

}
/*  Test fixture class definition */
class cunmtr_test  : public  ::testing::Test {
public:
   unmtr_scomplex_parameters  *cunmtr_obj;
   void SetUp();
   void TearDown () { delete cunmtr_obj;}
};

void cunmtr_test::SetUp(){

	/* LAPACKE chetrd prototype */
    typedef int (*Fptr_NL_LAPACKE_chetrd) (int matrix_layout, char uplo, lapack_int n,lapack_complex_float *A,
											lapack_int lda, float* d, float* e, lapack_complex_float* tau);
	 Fptr_NL_LAPACKE_chetrd chetrd;
	 /* LAPACKE cunmtr prototype */
    typedef int (*Fptr_NL_LAPACKE_cunmtr) (int matrix_layout, char side, char uplo, char trans, lapack_int m,  
											lapack_int n, lapack_complex_float *A, lapack_int lda, lapack_complex_float* tau,
											lapack_complex_float* c, lapack_int ldc);

    Fptr_NL_LAPACKE_cunmtr cunmtr;

    cunmtr_obj = new unmtr_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].uplo,
                           eig_paramslist[idx].m,
						   eig_paramslist[idx].n);
						   

    idx = Circular_Increment_Index(idx);

    cunmtr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cunmtr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cunmtr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cunmtr_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*cunmtr library call */
    cunmtr = (Fptr_NL_LAPACKE_cunmtr)dlsym(cunmtr_obj->hModule, "LAPACKE_cunmtr");
    ASSERT_TRUE(cunmtr != NULL) << "failed to get the Netlib LAPACKE_cunmtr symbol";

	/*chetrd library call*/
	cunmtr_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cunmtr_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cunmtr_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cunmtr_obj->lModule != NULL) << "Netlib lapacke handle NULL";

    chetrd = (Fptr_NL_LAPACKE_chetrd)dlsym(cunmtr_obj->lModule, "LAPACKE_chetrd");
    ASSERT_TRUE(chetrd != NULL) << "failed to get the Netlib LAPACKE_chetrd symbol";    
    

    cunmtr_obj->inforef_hetrd = chetrd( cunmtr_obj->matrix_layout, cunmtr_obj->uplo,
								cunmtr_obj->n, cunmtr_obj->Aref, cunmtr_obj->lda, cunmtr_obj->dref, 
								cunmtr_obj->eref, cunmtr_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    cunmtr_obj->info_hetrd = LAPACKE_chetrd( cunmtr_obj->matrix_layout, cunmtr_obj->uplo,
										cunmtr_obj->n,cunmtr_obj->A, cunmtr_obj->lda, cunmtr_obj->d,
										cunmtr_obj->e,cunmtr_obj->tau);

    if( cunmtr_obj->info_hetrd < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_chetrd is wrong\n", cunmtr_obj->info_hetrd );
    }
    if( cunmtr_obj->inforef_hetrd < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chetrd is wrong\n", 
        cunmtr_obj->inforef_hetrd );
    }  
/*Compute cunmtr's  o/p */
    cunmtr_obj->inforef = cunmtr( cunmtr_obj->matrix_layout, cunmtr_obj->side, cunmtr_obj->uplo, cunmtr_obj->trans, cunmtr_obj->m,
								cunmtr_obj->n, cunmtr_obj->Aref, cunmtr_obj->lda, cunmtr_obj->tauref, 
								cunmtr_obj->cref, cunmtr_obj->ldc);

    /* Compute libflame's Lapacke o/p  */
    cunmtr_obj->info = LAPACKE_cunmtr( cunmtr_obj->matrix_layout, cunmtr_obj->side, cunmtr_obj->uplo,cunmtr_obj->trans, cunmtr_obj->m, 
										cunmtr_obj->n, cunmtr_obj->A, cunmtr_obj->lda, cunmtr_obj->tau,
										cunmtr_obj->c,cunmtr_obj->ldc);

    if( cunmtr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cunmtr is wrong\n", cunmtr_obj->info );
    }
    if( cunmtr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cunmtr is wrong\n", 
        cunmtr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cunmtr_obj->diff =  computeDiff_c( cunmtr_obj->bufsize_c, 
                cunmtr_obj->c, cunmtr_obj->cref );

}

TEST_F(cunmtr_test, cunmtr1) {
    EXPECT_NEAR(0.0, cunmtr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cunmtr_test, cunmtr2) {
    EXPECT_NEAR(0.0, cunmtr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cunmtr_test, cunmtr3) {
    EXPECT_NEAR(0.0, cunmtr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cunmtr_test, cunmtr4) {
    EXPECT_NEAR(0.0, cunmtr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class unmtr_dcomplex_parameters{

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
	lapack_int r;
	lapack_complex_double* A;
	char uplo;
	char side;
	char trans;
	lapack_int ldc;
	lapack_int lda;
	/*Output Parameter*/
	lapack_complex_double* tau, *c;
	double* d;
	double* e;
	lapack_complex_double *Aref, *tauref, *cref;
	double* dref;
	double* eref;
	/*Return Values*/	
	lapack_int info, inforef;
	lapack_int info_hetrd, inforef_hetrd;

   public:
      unmtr_dcomplex_parameters (int matrix_layout,  char uplo, lapack_int m, lapack_int n);
      ~unmtr_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
unmtr_dcomplex_parameters:: unmtr_dcomplex_parameters (int matrix_layout_i, char uplo_i, lapack_int m_i, lapack_int n_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	uplo = uplo_i;
	trans = 'N';
	#if LAPACKE_TEST_VERBOSE
	printf(" \n unmtr scomplex: matrix_layout = %d, uplo:%c,m:%d, n: %d \n", matrix_layout, uplo, m, n);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{	ldc = m;
		r = n;		
		side = 'R';		
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{	ldc = n;
		r = m;
		side = 'L';
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	lda = r;
	bufsize_c = ldc*r;
	bufsize_a = lda*r;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&tau, &tauref, (r-1));
	lapacke_gtest_alloc_double_buffer_pair(&d, &dref, r);
	lapacke_gtest_alloc_double_buffer_pair(&e, &eref, (r-1));
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&c, &cref, bufsize_c);
	if ((A==NULL) || (Aref==NULL) || 
		(tau==NULL) || (tauref==NULL) ||
		(d == NULL) || (dref == NULL) ||
		(e == NULL) || (eref == NULL) ||
		(c == NULL) || (cref == NULL)){
		EXPECT_FALSE( true) << "unmtr_float_parameters object: malloc error.";
		unmtr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( A, Aref, bufsize_a);
	lapacke_gtest_init_dcomplex_buffer_pair_rand( c, cref, bufsize_c);
	lapacke_gtest_init_dcomplex_buffer_pair_rand( tau, tauref, (r-1));
	/*initialize output matrix by 0 *
	for(i=0;i<(r-1);i++) {
		tau[i] = 0;
		tauref[i] = tau[i];
	} */
	/*initialize output matrix by 0 */
	for(i=0;i<(r-1);i++) {
		e[i] = 0;
		eref[i] = e[i];
	}
	/*initialize output matrix by 0 */
	for(i=0;i<r;i++) {
		d[i] = 0;
		dref[i] = d[i];
	}
	/*for(i=0;i<bufsize_c;i++) {
		c[i] = 0.0;
		cref[i] = c[i];
	}*/

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
unmtr_dcomplex_parameters :: ~unmtr_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" unmtr_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   unmtr_free();

}
/*  Test fixture class definition */
class zunmtr_test  : public  ::testing::Test {
public:
   unmtr_dcomplex_parameters  *zunmtr_obj;
   void SetUp();
   void TearDown () { delete zunmtr_obj;}
};

void zunmtr_test::SetUp(){

	/* LAPACKE chetrd prototype */
    typedef int (*Fptr_NL_LAPACKE_zhetrd) (int matrix_layout, char uplo, lapack_int n,lapack_complex_double *A, lapack_int lda,
											double* d, double* e, lapack_complex_double* tau);
	 Fptr_NL_LAPACKE_zhetrd zhetrd;
	 /* LAPACKE zunmtr prototype */
    typedef int (*Fptr_NL_LAPACKE_zunmtr) (int matrix_layout, char side, char uplo, char trans, lapack_int m, 
											lapack_int n, lapack_complex_double *A, lapack_int lda,lapack_complex_double* tau,
											lapack_complex_double* c, lapack_int ldc);

    Fptr_NL_LAPACKE_zunmtr zunmtr;

    zunmtr_obj = new unmtr_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].uplo,
                           eig_paramslist[idx].m,
						   eig_paramslist[idx].n);
						   

    idx = Circular_Increment_Index(idx);

    zunmtr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zunmtr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zunmtr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zunmtr_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*zunmtr library call */
    zunmtr = (Fptr_NL_LAPACKE_zunmtr)dlsym(zunmtr_obj->hModule, "LAPACKE_zunmtr");
    ASSERT_TRUE(zunmtr != NULL) << "failed to get the Netlib LAPACKE_zunmtr symbol";

	/*chetrd library call*/
	zunmtr_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zunmtr_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zunmtr_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zunmtr_obj->lModule != NULL) << "Netlib lapacke handle NULL";

    zhetrd = (Fptr_NL_LAPACKE_zhetrd)dlsym(zunmtr_obj->lModule, "LAPACKE_zhetrd");
    ASSERT_TRUE(zhetrd != NULL) << "failed to get the Netlib LAPACKE_chetrd symbol";    
    

    zunmtr_obj->inforef_hetrd = zhetrd( zunmtr_obj->matrix_layout, zunmtr_obj->uplo,
								zunmtr_obj->n, zunmtr_obj->Aref, zunmtr_obj->lda, zunmtr_obj->dref, 
								zunmtr_obj->eref, zunmtr_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    zunmtr_obj->info_hetrd = LAPACKE_zhetrd( zunmtr_obj->matrix_layout, zunmtr_obj->uplo,
										zunmtr_obj->n,zunmtr_obj->A, zunmtr_obj->lda, zunmtr_obj->d,
										zunmtr_obj->e,zunmtr_obj->tau);

    if( zunmtr_obj->info_hetrd < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_chetrd is wrong\n", zunmtr_obj->info_hetrd );
    }
    if( zunmtr_obj->inforef_hetrd < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chetrd is wrong\n", 
        zunmtr_obj->inforef_hetrd );
    }  
/*Compute zunmtr's  o/p */
    zunmtr_obj->inforef = zunmtr( zunmtr_obj->matrix_layout, zunmtr_obj->side, zunmtr_obj->uplo, zunmtr_obj->trans, zunmtr_obj->m,
								zunmtr_obj->n, zunmtr_obj->Aref, zunmtr_obj->lda, zunmtr_obj->tauref, 
								zunmtr_obj->cref, zunmtr_obj->ldc);

    /* Compute libflame's Lapacke o/p  */
    zunmtr_obj->info = LAPACKE_zunmtr( zunmtr_obj->matrix_layout, zunmtr_obj->side, zunmtr_obj->uplo,zunmtr_obj->trans, zunmtr_obj->m, 
										zunmtr_obj->n,zunmtr_obj->A, zunmtr_obj->lda, zunmtr_obj->tau,
										zunmtr_obj->c,zunmtr_obj->ldc);

    if( zunmtr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zunmtr is wrong\n", zunmtr_obj->info );
    }
    if( zunmtr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zunmtr is wrong\n", 
        zunmtr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zunmtr_obj->diff =  computeDiff_z( zunmtr_obj->bufsize_c, 
                zunmtr_obj->c, zunmtr_obj->cref );

}

TEST_F(zunmtr_test, zunmtr1) {
    EXPECT_NEAR(0.0, zunmtr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zunmtr_test, zunmtr2) {
    EXPECT_NEAR(0.0, zunmtr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zunmtr_test, zunmtr3) {
    EXPECT_NEAR(0.0, zunmtr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zunmtr_test, zunmtr4) {
    EXPECT_NEAR(0.0, zunmtr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


	
