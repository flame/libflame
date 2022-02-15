#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"

#define unmbr_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (tau!=NULL)  free(tau);\
if (tauref!=NULL) free(tauref); \
if (taup!=NULL)  free(taup);\
if (taupref!=NULL) free(taupref); \
if (tauq!=NULL)  free(tauq);\
if (tauqref!=NULL) free(tauqref); \
if (c != NULL) free(c); \
if (cref != NULL) free(cref); \
if ( d!= NULL) free(d); \
if ( dref!= NULL) free(dref); \
if ( e!= NULL) free(e); \
if ( eref!= NULL) free(eref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin scomplex_common_parameters  class definition */
class unmbr_scomplex_parameters{

   public:
	int bufsize_a;
	int bufsize_c;
	void *hModule, *dModule;
	void *bModule, *lModule;
	float diff;
	float* d, *e;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int k;
	lapack_int r;
	lapack_complex_float* A;
	char vect;
	char side;
	char trans;
	lapack_int ldc;
	lapack_int lda;
	/*Output Parameter*/
	lapack_complex_float* tau, *c;
	lapack_complex_float* taup, *tauq;
	lapack_complex_float* taupref, *tauqref;
	lapack_complex_float *Aref, *tauref, *cref;
	float* dref, *eref;
	/*Return Values*/	
	lapack_int info, inforef;
	lapack_int info_gebrd, inforef_gebrd;

   public:
      unmbr_scomplex_parameters (int matrix_layout,  char vect, char side , char trans, lapack_int m, lapack_int n);
      ~unmbr_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
unmbr_scomplex_parameters:: unmbr_scomplex_parameters (int matrix_layout_i, char vect_i, char side_i, char trans_i, lapack_int m_i, lapack_int n_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	vect = vect_i;
	side = side_i;
	trans = trans_i;
	
	//trans = 'N';
	#if LAPACKE_TEST_VERBOSE
	printf(" \n unmbr scomplex: matrix_layout = %d, vect:%c, side:%c, trans:%c, m:%d, n: %d \n", matrix_layout, vect, side, trans, m, n);
	#endif
	
	if (trans == 'T')
		trans = 'N';
	
	ldc = m;
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{
		bufsize_c = ldc *n;
		if (vect == 'Q')
		{	k = n;
			r = m;
			lda = r;
			bufsize_a  = lda *k;		
		}
		else if (vect == 'P')
		{
			k = m;
			
			if (side == 'L')
			{	
				r = m;
				lda  = min(r, k);
				bufsize_a  = lda *m;
			}else if (side == 'R')
			{	
				r = n;
				lda  = min(r, k);
				bufsize_a  = lda *n;
			}
			
		}			
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{
		bufsize_c = ldc*m;
		if (vect == 'P')
		{	k = m;
			r = n;
			lda = r;
			bufsize_a  = lda *k;
					
		}else if (vect == 'Q')
		{
			k = n;
			lda = k;
			if (side == 'L')
			{	bufsize_a  = lda *m;
				r = m; 
			}else if (side == 'R')
			{	bufsize_a  = lda *n;
				r = n;
			}
		}
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";

	printf(" \n unmbr scomplex: matrix_layout = %d, vect:%c, side:%c, trans:%c, m:%d, n: %d \n", matrix_layout, vect, side, trans, m, n);
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&tau, &tauref, min(r,k));
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&taup, &taupref, min(m,n));
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&tauq, &tauqref, min(m,n));
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&c, &cref, bufsize_c);
	lapacke_gtest_alloc_float_buffer_pair(&d, &dref, min(m,n));
	lapacke_gtest_alloc_float_buffer_pair(&e, &eref, (min(m,n) -1));
	if ((A==NULL) || (Aref==NULL) || \
		(tau==NULL) || (tauref==NULL) ||
		(taup==NULL) || (taupref==NULL) ||
		(tauq==NULL) || (tauqref==NULL) ||
		(c == NULL) || (cref == NULL)||
		(e == NULL) || (eref == NULL)||
		(d == NULL) || (dref == NULL)){
		EXPECT_FALSE( true) << "unmbr_float_parameters object: malloc error.";
		unmbr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( A, Aref, bufsize_a);
} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
unmbr_scomplex_parameters :: ~unmbr_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" unmbr_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   unmbr_free();

}
/*  Test fixture class definition */
class cunmbr_test  : public  ::testing::Test {
public:
   unmbr_scomplex_parameters  *cunmbr_obj;
   void SetUp();
   void TearDown () { delete cunmbr_obj;}
};

void cunmbr_test::SetUp(){

	/* LAPACKE cgebrd prototype */
    typedef int (*Fptr_NL_LAPACKE_cgebrd) (int matrix_layout, lapack_int m, lapack_int n, lapack_complex_float* a, lapack_int lda, float* d, float* e, 
	lapack_complex_float* tauq, lapack_complex_float* taup);
	
	 Fptr_NL_LAPACKE_cgebrd cgebrd;
	 /* LAPACKE cunmbr prototype */
    typedef int (*Fptr_NL_LAPACKE_cunmbr) (int matrix_layout, char vect, char side, char trans, lapack_int m, lapack_int n, lapack_int k, 
	const lapack_complex_float* a, lapack_int lda, const lapack_complex_float* tau, lapack_complex_float* c, lapack_int ldc);

    Fptr_NL_LAPACKE_cunmbr cunmbr;

    cunmbr_obj = new unmbr_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].vect,
						   eig_paramslist[idx].side,
						   eig_paramslist[idx].trans,
                           eig_paramslist[idx].m,
						   eig_paramslist[idx].n);						   

    idx = Circular_Increment_Index(idx);

    cunmbr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cunmbr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cunmbr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cunmbr_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*cunmbr library call */
    cunmbr = (Fptr_NL_LAPACKE_cunmbr)dlsym(cunmbr_obj->hModule, "LAPACKE_cunmbr");
    ASSERT_TRUE(cunmbr != NULL) << "failed to get the Netlib LAPACKE_cunmbr symbol";

	/*cgebrd library call*/
	cunmbr_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cunmbr_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cunmbr_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cunmbr_obj->lModule != NULL) << "Netlib lapacke handle NULL";

    cgebrd = (Fptr_NL_LAPACKE_cgebrd)dlsym(cunmbr_obj->lModule, "LAPACKE_cgebrd");
    ASSERT_TRUE(cgebrd != NULL) << "failed to get the Netlib LAPACKE_cgebrd symbol";    
    

    cunmbr_obj->inforef_gebrd = cgebrd( cunmbr_obj->matrix_layout, cunmbr_obj->m, cunmbr_obj->n, cunmbr_obj->Aref, cunmbr_obj->lda, cunmbr_obj->dref, 
								cunmbr_obj->eref, cunmbr_obj->tauqref, cunmbr_obj->taupref);

    /* Compute libflame's Lapacke o/p  */
    cunmbr_obj->info_gebrd = LAPACKE_cgebrd(cunmbr_obj->matrix_layout, cunmbr_obj->m, cunmbr_obj->n, cunmbr_obj->A, cunmbr_obj->lda, cunmbr_obj->d,
	cunmbr_obj->e, cunmbr_obj->tauq, cunmbr_obj->taup);

    if( cunmbr_obj->info_gebrd < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cgebrd is wrong\n", cunmbr_obj->info_gebrd );
    }
    if( cunmbr_obj->inforef_gebrd < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cgebrd is wrong\n", 
        cunmbr_obj->inforef_gebrd );
    }  
/*Compute cunmbr's  o/p */
    cunmbr_obj->inforef = cunmbr( cunmbr_obj->matrix_layout, cunmbr_obj->vect, cunmbr_obj->side , cunmbr_obj->trans, cunmbr_obj->m,
								cunmbr_obj->n, cunmbr_obj->k, (const lapack_complex_float*)cunmbr_obj->Aref, cunmbr_obj->lda, (const lapack_complex_float*)cunmbr_obj->tauref, 
								cunmbr_obj->cref, cunmbr_obj->ldc);

    /* Compute libflame's Lapacke o/p  */	
    cunmbr_obj->info = LAPACKE_cunmbr( cunmbr_obj->matrix_layout, cunmbr_obj->vect, cunmbr_obj->side,cunmbr_obj->trans, cunmbr_obj->m, 
										cunmbr_obj->n, cunmbr_obj->k, (const lapack_complex_float*) cunmbr_obj->A, cunmbr_obj->lda, (const lapack_complex_float*)cunmbr_obj->tau,
										cunmbr_obj->c,cunmbr_obj->ldc);

    if( cunmbr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cunmbr is wrong\n", cunmbr_obj->info );
    }
    if( cunmbr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cunmbr is wrong\n", 
        cunmbr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cunmbr_obj->diff =  computeDiff_c( cunmbr_obj->bufsize_c, 
                cunmbr_obj->c, cunmbr_obj->cref );

}

TEST_F(cunmbr_test, cunmbr1) {
    EXPECT_NEAR(0.0, cunmbr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cunmbr_test, cunmbr2) {
    EXPECT_NEAR(0.0, cunmbr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cunmbr_test, cunmbr3) {
    EXPECT_NEAR(0.0, cunmbr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cunmbr_test, cunmbr4) {
    EXPECT_NEAR(0.0, cunmbr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class unmbr_dcomplex_parameters{

   public:
	int bufsize_a;
	int bufsize_c;
	void *hModule, *dModule;
	void *bModule, *lModule;
	double diff;
	double* d, *e;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int k;
	lapack_int r;
	lapack_complex_double* A;
	char vect;
	char side;
	char trans;
	lapack_int ldc;
	lapack_int lda;
	/*Output Parameter*/
	lapack_complex_double* tau, *c;
	lapack_complex_double* taup, *tauq;
	lapack_complex_double* taupref, *tauqref;
	lapack_complex_double *Aref, *tauref, *cref;
	double* dref, *eref;
	/*Return Values*/	
	lapack_int info, inforef;
	lapack_int info_gebrd, inforef_gebrd;

   public:
      unmbr_dcomplex_parameters (int matrix_layout,  char vect, char side , char trans, lapack_int m, lapack_int n);
      ~unmbr_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
unmbr_dcomplex_parameters:: unmbr_dcomplex_parameters (int matrix_layout_i, char vect_i, char side_i, char trans_i, lapack_int m_i, lapack_int n_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	vect = vect_i;
	side = side_i;
	trans = trans_i;
	
	//trans = 'N';
	#if LAPACKE_TEST_VERBOSE
	printf(" \n unmbr dcomplex: matrix_layout = %d, vect:%c, side:%c, trans:%c, m:%d, n: %d \n", matrix_layout, vect, side, trans, m, n);
	#endif
	
	if (trans == 'T')
		trans = 'N';
	
	ldc = m;
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{
		bufsize_c = ldc *n;
		if (vect == 'Q')
		{	k = n;
			r = m;
			lda = r;
			bufsize_a  = lda *k;		
		}
		else if (vect == 'P')
		{
			k = m;
			
			if (side == 'L')
			{	
				r = m;
				lda  = min(r, k);
				bufsize_a  = lda *m;
			}else if (side == 'R')
			{	
				r = n;
				lda  = min(r, k);
				bufsize_a  = lda *n;
			}
			
		}			
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{
		bufsize_c = ldc*m;
		if (vect == 'P')
		{	k = m;
			r = n;
			lda = r;
			bufsize_a  = lda *k;
					
		}else if (vect == 'Q')
		{
			k = n;
			lda = k;
			if (side == 'L')
			{	bufsize_a  = lda *m;
				r = m; 
			}else if (side == 'R')
			{	bufsize_a  = lda *n;
				r = n;
			}
		}
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";

	printf(" \n unmbr dcomplex: matrix_layout = %d, vect:%c, side:%c, trans:%c, m:%d, n: %d \n", matrix_layout, vect, side, trans, m, n);
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&tau, &tauref, min(r,k));
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&taup, &taupref, min(m,n));
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&tauq, &tauqref, min(m,n));
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&c, &cref, bufsize_c);
	lapacke_gtest_alloc_double_buffer_pair(&d, &dref, min(m,n));
	lapacke_gtest_alloc_double_buffer_pair(&e, &eref, (min(m,n) -1));
	if ((A==NULL) || (Aref==NULL) || \
		(tau==NULL) || (tauref==NULL) ||
		(taup==NULL) || (taupref==NULL) ||
		(tauq==NULL) || (tauqref==NULL) ||
		(c == NULL) || (cref == NULL)||
		(e == NULL) || (eref == NULL)||
		(d == NULL) || (dref == NULL)){
		EXPECT_FALSE( true) << "unmbr_double_parameters object: malloc error.";
		unmbr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( A, Aref, bufsize_a);
} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
unmbr_dcomplex_parameters :: ~unmbr_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" unmbr_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   unmbr_free();

}
/*  Test fixture class definition */
class zunmbr_test  : public  ::testing::Test {
public:
   unmbr_dcomplex_parameters  *zunmbr_obj;
   void SetUp();
   void TearDown () { delete zunmbr_obj;}
};

void zunmbr_test::SetUp(){

	/* LAPACKE zgebrd prototype */
    typedef int (*Fptr_NL_LAPACKE_zgebrd) (int matrix_layout, lapack_int m, lapack_int n, lapack_complex_double* a, lapack_int lda, double* d, double* e, 
	lapack_complex_double* tauq, lapack_complex_double* taup);
	
	 Fptr_NL_LAPACKE_zgebrd zgebrd;
	 /* LAPACKE zunmbr prototype */
    typedef int (*Fptr_NL_LAPACKE_zunmbr) (int matrix_layout, char vect, char side, char trans, lapack_int m, lapack_int n, lapack_int k, 
	const lapack_complex_double* a, lapack_int lda, const lapack_complex_double* tau, lapack_complex_double* c, lapack_int ldc);

    Fptr_NL_LAPACKE_zunmbr zunmbr;

    zunmbr_obj = new unmbr_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].vect,
						   eig_paramslist[idx].side,
						   eig_paramslist[idx].trans,
                           eig_paramslist[idx].m,
						   eig_paramslist[idx].n);						   

    idx = Circular_Increment_Index(idx);

    zunmbr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zunmbr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zunmbr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zunmbr_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*zunmbr library call */
    zunmbr = (Fptr_NL_LAPACKE_zunmbr)dlsym(zunmbr_obj->hModule, "LAPACKE_zunmbr");
    ASSERT_TRUE(zunmbr != NULL) << "failed to get the Netlib LAPACKE_zunmbr symbol";

	/*zgebrd library call*/
	zunmbr_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zunmbr_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zunmbr_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zunmbr_obj->lModule != NULL) << "Netlib lapacke handle NULL";

    zgebrd = (Fptr_NL_LAPACKE_zgebrd)dlsym(zunmbr_obj->lModule, "LAPACKE_zgebrd");
    ASSERT_TRUE(zgebrd != NULL) << "failed to get the Netlib LAPACKE_zgebrd symbol";    
    

    zunmbr_obj->inforef_gebrd = zgebrd( zunmbr_obj->matrix_layout, zunmbr_obj->m, zunmbr_obj->n, zunmbr_obj->Aref, zunmbr_obj->lda, zunmbr_obj->dref, 
								zunmbr_obj->eref, zunmbr_obj->tauqref, zunmbr_obj->taupref);

    /* Compute libflame's Lapacke o/p  */
    zunmbr_obj->info_gebrd = LAPACKE_zgebrd(zunmbr_obj->matrix_layout, zunmbr_obj->m, zunmbr_obj->n, zunmbr_obj->A, zunmbr_obj->lda, zunmbr_obj->d,
	zunmbr_obj->e, zunmbr_obj->tauq, zunmbr_obj->taup);

    if( zunmbr_obj->info_gebrd < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zgebrd is wrong\n", zunmbr_obj->info_gebrd );
    }
    if( zunmbr_obj->inforef_gebrd < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zgebrd is wrong\n", 
        zunmbr_obj->inforef_gebrd );
    }  
/*Compute zunmbr's  o/p */
    zunmbr_obj->inforef = zunmbr( zunmbr_obj->matrix_layout, zunmbr_obj->vect, zunmbr_obj->side , zunmbr_obj->trans, zunmbr_obj->m,
								zunmbr_obj->n, zunmbr_obj->k, (const lapack_complex_double*)zunmbr_obj->Aref, zunmbr_obj->lda, (const lapack_complex_double*)zunmbr_obj->tauref, 
								zunmbr_obj->cref, zunmbr_obj->ldc);

    /* Compute libflame's Lapacke o/p  */	
    zunmbr_obj->info = LAPACKE_zunmbr( zunmbr_obj->matrix_layout, zunmbr_obj->vect, zunmbr_obj->side,zunmbr_obj->trans, zunmbr_obj->m, 
										zunmbr_obj->n, zunmbr_obj->k, (const lapack_complex_double*) zunmbr_obj->A, zunmbr_obj->lda, (const lapack_complex_double*)zunmbr_obj->tau,
										zunmbr_obj->c,zunmbr_obj->ldc);

    if( zunmbr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zunmbr is wrong\n", zunmbr_obj->info );
    }
    if( zunmbr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zunmbr is wrong\n", 
        zunmbr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zunmbr_obj->diff =  computeDiff_z( zunmbr_obj->bufsize_c, 
                zunmbr_obj->c, zunmbr_obj->cref );

}

TEST_F(zunmbr_test, zunmbr1) {
    EXPECT_NEAR(0.0, zunmbr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zunmbr_test, zunmbr2) {
    EXPECT_NEAR(0.0, zunmbr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zunmbr_test, zunmbr3) {
    EXPECT_NEAR(0.0, zunmbr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zunmbr_test, zunmbr4) {
    EXPECT_NEAR(0.0, zunmbr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}