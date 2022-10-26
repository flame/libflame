#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"




#define upmtr_free() \
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
class upmtr_scomplex_parameters{

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
	/*Output Parameter*/
	lapack_complex_float* tau, *c;
	float* d;
	float* e;
	lapack_complex_float *Aref, *tauref, *cref;
	float* dref;
	float* eref;
	/*Return Values*/	
	lapack_int info, inforef;
	lapack_int info_hptrd, inforef_hptrd;

   public:
      upmtr_scomplex_parameters (int matrix_layout,  char uplo, lapack_int m, lapack_int n);
      ~upmtr_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
upmtr_scomplex_parameters:: upmtr_scomplex_parameters (int matrix_layout_i, char uplo_i, lapack_int m_i, lapack_int n_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	uplo = uplo_i;
	trans = 'N';
	#if LAPACKE_TEST_VERBOSE
	printf(" \n upmtr scomplex: matrix_layout = %d, uplo:%c,m:%d, n: %d \n", matrix_layout, uplo, m, n);
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
	
	bufsize_c = ldc*r;
	bufsize_a = (r*(r+1)/2);
	
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
		EXPECT_FALSE( true) << "upmtr_float_parameters object: malloc error.";
		upmtr_free();
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
upmtr_scomplex_parameters :: ~upmtr_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" upmtr_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   upmtr_free();

}
/*  Test fixture class definition */
class cupmtr_test  : public  ::testing::Test {
public:
   upmtr_scomplex_parameters  *cupmtr_obj;
   void SetUp();
   void TearDown () { delete cupmtr_obj;}
};

void cupmtr_test::SetUp(){

	/* LAPACKE chptrd prototype */
    typedef int (*Fptr_NL_LAPACKE_chptrd) (int matrix_layout, char uplo, lapack_int n,lapack_complex_float *A, 
											float* d, float* e, lapack_complex_float* tau);
	 Fptr_NL_LAPACKE_chptrd chptrd;
	 /* LAPACKE cupmtr prototype */
    typedef int (*Fptr_NL_LAPACKE_cupmtr) (int matrix_layout, char side, char uplo, char trans, lapack_int m, 
											lapack_int n,  const lapack_complex_float *A, const lapack_complex_float* tau,
											lapack_complex_float* c, lapack_int ldc);

    Fptr_NL_LAPACKE_cupmtr cupmtr;

    cupmtr_obj = new upmtr_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].uplo,
                           eig_paramslist[idx].m,
						   eig_paramslist[idx].n);
						   

    idx = Circular_Increment_Index(idx);

    cupmtr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cupmtr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cupmtr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cupmtr_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*cupmtr library call */
    cupmtr = (Fptr_NL_LAPACKE_cupmtr)dlsym(cupmtr_obj->hModule, "LAPACKE_cupmtr");
    ASSERT_TRUE(cupmtr != NULL) << "failed to get the Netlib LAPACKE_cupmtr symbol";

	/*chptrd library call*/
	cupmtr_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cupmtr_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cupmtr_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cupmtr_obj->lModule != NULL) << "Netlib lapacke handle NULL";

    chptrd = (Fptr_NL_LAPACKE_chptrd)dlsym(cupmtr_obj->lModule, "LAPACKE_chptrd");
    ASSERT_TRUE(chptrd != NULL) << "failed to get the Netlib LAPACKE_chptrd symbol";    
    

    cupmtr_obj->inforef_hptrd = chptrd( cupmtr_obj->matrix_layout, cupmtr_obj->uplo,
								cupmtr_obj->n, cupmtr_obj->Aref,cupmtr_obj->dref, 
								cupmtr_obj->eref, cupmtr_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    cupmtr_obj->info_hptrd = LAPACKE_chptrd( cupmtr_obj->matrix_layout, cupmtr_obj->uplo,
										cupmtr_obj->n,cupmtr_obj->A, cupmtr_obj->d,
										cupmtr_obj->e,cupmtr_obj->tau);

    if( cupmtr_obj->info_hptrd < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_chptrd is wrong\n", cupmtr_obj->info_hptrd );
    }
    if( cupmtr_obj->inforef_hptrd < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chptrd is wrong\n", 
        cupmtr_obj->inforef_hptrd );
    }  
/*Compute cupmtr's  o/p */
    cupmtr_obj->inforef = cupmtr( cupmtr_obj->matrix_layout, cupmtr_obj->side, cupmtr_obj->uplo, cupmtr_obj->trans, cupmtr_obj->m,
								cupmtr_obj->n, (const lapack_complex_float *)cupmtr_obj->Aref,(const lapack_complex_float *)cupmtr_obj->tauref, 
								cupmtr_obj->cref, cupmtr_obj->ldc);

    /* Compute libflame's Lapacke o/p  */
    cupmtr_obj->info = LAPACKE_cupmtr( cupmtr_obj->matrix_layout, cupmtr_obj->side, cupmtr_obj->uplo,cupmtr_obj->trans, cupmtr_obj->m, 
										cupmtr_obj->n,(const lapack_complex_float *)cupmtr_obj->A, (const lapack_complex_float*)cupmtr_obj->tau,
										cupmtr_obj->c,cupmtr_obj->ldc);

    if( cupmtr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cupmtr is wrong\n", cupmtr_obj->info );
    }
    if( cupmtr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cupmtr is wrong\n", 
        cupmtr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cupmtr_obj->diff =  computeDiff_c( cupmtr_obj->bufsize_c, 
                cupmtr_obj->c, cupmtr_obj->cref );

}

TEST_F(cupmtr_test, cupmtr1) {
    EXPECT_NEAR(0.0, cupmtr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cupmtr_test, cupmtr2) {
    EXPECT_NEAR(0.0, cupmtr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cupmtr_test, cupmtr3) {
    EXPECT_NEAR(0.0, cupmtr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cupmtr_test, cupmtr4) {
    EXPECT_NEAR(0.0, cupmtr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class upmtr_dcomplex_parameters{

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
	/*Output Parameter*/
	lapack_complex_double* tau, *c;
	double* d;
	double* e;
	lapack_complex_double *Aref, *tauref, *cref;
	double* dref;
	double* eref;
	/*Return Values*/	
	lapack_int info, inforef;
	lapack_int info_hptrd, inforef_hptrd;

   public:
      upmtr_dcomplex_parameters (int matrix_layout,  char uplo, lapack_int m, lapack_int n);
      ~upmtr_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
upmtr_dcomplex_parameters:: upmtr_dcomplex_parameters (int matrix_layout_i, char uplo_i, lapack_int m_i, lapack_int n_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	uplo = uplo_i;
	trans = 'N';
	#if LAPACKE_TEST_VERBOSE
	printf(" \n upmtr scomplex: matrix_layout = %d, uplo:%c,m:%d, n: %d \n", matrix_layout, uplo, m, n);
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
	
	bufsize_c = ldc*r;
	bufsize_a = (r*(r+1)/2);
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&tau, &tauref, (r-1));
	lapacke_gtest_alloc_double_buffer_pair(&d, &dref, r);
	lapacke_gtest_alloc_double_buffer_pair(&e, &eref, (r-1));
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&c, &cref, bufsize_c);
	if ((A==NULL) || (Aref==NULL) || \
		(tau==NULL) || (tauref==NULL) ||
		(d == NULL) || (dref == NULL) ||
		(e == NULL) || (eref == NULL) ||
		(c == NULL) || (cref == NULL)){
		EXPECT_FALSE( true) << "upmtr_float_parameters object: malloc error.";
		upmtr_free();
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
upmtr_dcomplex_parameters :: ~upmtr_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" upmtr_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   upmtr_free();

}
/*  Test fixture class definition */
class zupmtr_test  : public  ::testing::Test {
public:
   upmtr_dcomplex_parameters  *zupmtr_obj;
   void SetUp();
   void TearDown () { delete zupmtr_obj;}
};

void zupmtr_test::SetUp(){

	/* LAPACKE chptrd prototype */
    typedef int (*Fptr_NL_LAPACKE_zhptrd) (int matrix_layout, char uplo, lapack_int n,lapack_complex_double *A, 
											double* d, double* e, lapack_complex_double* tau);
	 Fptr_NL_LAPACKE_zhptrd zhptrd;
	 /* LAPACKE zupmtr prototype */
    typedef int (*Fptr_NL_LAPACKE_zupmtr) (int matrix_layout, char side, char uplo, char trans, lapack_int m, 
										lapack_int n,  const lapack_complex_double *A, const lapack_complex_double* tau,	
											lapack_complex_double* c, lapack_int ldc);

    Fptr_NL_LAPACKE_zupmtr zupmtr;

    zupmtr_obj = new upmtr_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].uplo,
                           eig_paramslist[idx].m,
						   eig_paramslist[idx].n);
						   

    idx = Circular_Increment_Index(idx);

    zupmtr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zupmtr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zupmtr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zupmtr_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*zupmtr library call */
    zupmtr = (Fptr_NL_LAPACKE_zupmtr)dlsym(zupmtr_obj->hModule, "LAPACKE_zupmtr");
    ASSERT_TRUE(zupmtr != NULL) << "failed to get the Netlib LAPACKE_zupmtr symbol";

	/*chptrd library call*/
	zupmtr_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zupmtr_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zupmtr_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zupmtr_obj->lModule != NULL) << "Netlib lapacke handle NULL";

    zhptrd = (Fptr_NL_LAPACKE_zhptrd)dlsym(zupmtr_obj->lModule, "LAPACKE_zhptrd");
    ASSERT_TRUE(zhptrd != NULL) << "failed to get the Netlib LAPACKE_chptrd symbol";    
    

    zupmtr_obj->inforef_hptrd = zhptrd( zupmtr_obj->matrix_layout, zupmtr_obj->uplo,
								zupmtr_obj->n, zupmtr_obj->Aref,zupmtr_obj->dref, 
								zupmtr_obj->eref, zupmtr_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    zupmtr_obj->info_hptrd = LAPACKE_zhptrd( zupmtr_obj->matrix_layout, zupmtr_obj->uplo,
										zupmtr_obj->n,zupmtr_obj->A, zupmtr_obj->d,
										zupmtr_obj->e,zupmtr_obj->tau);

    if( zupmtr_obj->info_hptrd < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_chptrd is wrong\n", zupmtr_obj->info_hptrd );
    }
    if( zupmtr_obj->inforef_hptrd < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chptrd is wrong\n", 
        zupmtr_obj->inforef_hptrd );
    }  
/*Compute zupmtr's  o/p */
    zupmtr_obj->inforef = zupmtr( zupmtr_obj->matrix_layout, zupmtr_obj->side, zupmtr_obj->uplo, zupmtr_obj->trans, zupmtr_obj->m,
								zupmtr_obj->n, (const lapack_complex_double *)zupmtr_obj->Aref,(const lapack_complex_double *)zupmtr_obj->tauref, 
								zupmtr_obj->cref, zupmtr_obj->ldc);

    /* Compute libflame's Lapacke o/p  */
    zupmtr_obj->info = LAPACKE_zupmtr( zupmtr_obj->matrix_layout, zupmtr_obj->side, zupmtr_obj->uplo,zupmtr_obj->trans, zupmtr_obj->m, 
										zupmtr_obj->n,(const lapack_complex_double *)zupmtr_obj->A, (const lapack_complex_double*)zupmtr_obj->tau,
										zupmtr_obj->c,zupmtr_obj->ldc);

    if( zupmtr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zupmtr is wrong\n", zupmtr_obj->info );
    }
    if( zupmtr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zupmtr is wrong\n", 
        zupmtr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zupmtr_obj->diff =  computeDiff_z( zupmtr_obj->bufsize_c, 
                zupmtr_obj->c, zupmtr_obj->cref );

}

TEST_F(zupmtr_test, zupmtr1) {
    EXPECT_NEAR(0.0, zupmtr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zupmtr_test, zupmtr2) {
    EXPECT_NEAR(0.0, zupmtr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zupmtr_test, zupmtr3) {
    EXPECT_NEAR(0.0, zupmtr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zupmtr_test, zupmtr4) {
    EXPECT_NEAR(0.0, zupmtr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


	
