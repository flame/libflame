#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"

#define gemlq_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (A_gelq!=NULL)    free(A_gelq); \
if (A_gelqref!=NULL) free(A_gelqref);\
if (t!=NULL)  free(t);\
if (tref!=NULL) free(tref); \
if (c != NULL) free(c); \
if (cref != NULL) free(cref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule);\
if( bModule != NULL) dlclose(bModule); \
if( lModule != NULL) dlclose(lModule); \


// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin float_common_parameters  class definition */
class gemlq_float_parameters{

   public:
	int bufsize_a, bufsize_agelq;
	int bufsize_c;
	void *hModule, *dModule;
	void *bModule, *lModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int k;
	lapack_int tsize;
	float* A, *A_gelq;
	char side;
	char trans;
	lapack_int ldc, lda, lda_gelq;
	/*Output Parameter*/
	float* t, *c;	
	float *Aref, *tref, *cref, *A_gelqref;
	/*Return Values*/	
	lapack_int info, inforef;
	lapack_int info_gelq, inforef_gelq;

   public:
      gemlq_float_parameters (int matrix_layout, char side , char trans, lapack_int m, lapack_int n, lapack_int tsize);
      ~gemlq_float_parameters ();

};

/* Constructor definition  float_common_parameters */
gemlq_float_parameters:: gemlq_float_parameters (int matrix_layout_i, char side_i , char trans_i,lapack_int m_i, lapack_int n_i, lapack_int tsize_i)
{
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	tsize = tsize_i;
	side = side_i;
	trans = trans_i;
	
	if (trans == 'C')
		trans  = 'N' ;
	
	#if LAPACKE_TEST_VERBOSE
		printf(" \n gemlq float: matrix_layout = %d, side :%c, trans:%c, m:%d, n: %d \n", matrix_layout, side, trans, m, n);
	#endif
	if (side == 'L')
	{
		k = m;
		if (matrix_layout == LAPACK_COL_MAJOR)
		{		
			ldc = m;
			lda = k;
			lda_gelq = m;
			bufsize_a = (lda*m);
			bufsize_agelq = lda_gelq*n;
			bufsize_c = (ldc*n);			
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{		
			ldc = n;
			lda = m;
			lda_gelq =n;
			bufsize_a = (lda*k);
			bufsize_agelq = lda_gelq*m;
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
			bufsize_a = (lda*n);
			bufsize_c = (ldc*n);
			lda_gelq = m;
			bufsize_agelq = lda_gelq*n;
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{	ldc = n;
			lda = n;
			lda_gelq =n;
			bufsize_a = (lda*k);
			bufsize_c = (ldc*m);
			bufsize_agelq = lda_gelq*m;
		}else
		{
			EXPECT_TRUE(false) << "matrix_layout invalid";
		}
	}else 
		EXPECT_TRUE(false) << "side is invalid";
	tsize = bufsize_a;

	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_float_buffer_pair(&A_gelq, &A_gelqref, bufsize_agelq);
	lapacke_gtest_alloc_float_buffer_pair(&t, &tref, bufsize_a);	
	lapacke_gtest_alloc_float_buffer_pair(&c, &cref, bufsize_c);
	
	if ((A==NULL) || (Aref==NULL) ||
		(t==NULL) || (tref==NULL) ||
		(c==NULL) || (cref==NULL) ||
		(A_gelq == NULL) || (A_gelqref==NULL))
	{
		EXPECT_FALSE( true) << "gemlq_float_parameters object: malloc error.";
		gemlq_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand( A, Aref, bufsize_a);
	lapacke_gtest_init_float_buffer_pair_rand( A_gelq, A_gelqref, bufsize_agelq);
	lapacke_gtest_init_float_buffer_pair_rand( c, cref, bufsize_c);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
gemlq_float_parameters :: ~gemlq_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" gemlq_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   gemlq_free();

}
/*  Test fixture class definition */
class sgemlq_test  : public  ::testing::Test {
public:
   gemlq_float_parameters  *sgemlq_obj;
   void SetUp();
   void TearDown () { delete sgemlq_obj;}
};

void sgemlq_test::SetUp(){

	/* LAPACKE sgelq prototype */
    typedef int (*Fptr_NL_LAPACKE_sgelq) (int matrix_layout, lapack_int m,lapack_int n,\ 
											float *A, lapack_int lda, float* t, lapack_int tsize);
	 Fptr_NL_LAPACKE_sgelq sgelq;
	 /* LAPACKE sgemlq prototype */
    typedef int (*Fptr_NL_LAPACKE_sgemlq) (int matrix_layout, char side, char trans,\
                           lapack_int m, lapack_int n, lapack_int k,\
                           const float* a, lapack_int lda,\
                           const float* t, lapack_int tsize,\
                           float* c, lapack_int ldc);

    Fptr_NL_LAPACKE_sgemlq sgemlq;

    sgemlq_obj = new gemlq_float_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].side,
						   eig_paramslist[idx].trans,
						   eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
						   eig_paramslist[idx].m);
						   

    idx = Circular_Increment_Index(idx);

    sgemlq_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgemlq_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgemlq_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgemlq_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*sgemlq library call */
    sgemlq = (Fptr_NL_LAPACKE_sgemlq)dlsym(sgemlq_obj->hModule, "LAPACKE_sgemlq");
    ASSERT_TRUE(sgemlq != NULL) << "failed to get the Netlib LAPACKE_sgemlq symbol";

	/*sgelq library call*/
	sgemlq_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgemlq_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgemlq_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgemlq_obj->lModule != NULL) << "Netlib lapacke handle NULL";

    sgelq = (Fptr_NL_LAPACKE_sgelq)dlsym(sgemlq_obj->lModule, "LAPACKE_sgelq");
    ASSERT_TRUE(sgelq != NULL) << "failed to get the Netlib LAPACKE_sgelq symbol";    
    

    sgemlq_obj->inforef_gelq = sgelq( sgemlq_obj->matrix_layout,sgemlq_obj->m, sgemlq_obj->n,
								sgemlq_obj->A_gelqref, sgemlq_obj->lda_gelq, sgemlq_obj->tref, sgemlq_obj->tsize);

    /* Compute libflame's Lapacke o/p  */
    sgemlq_obj->info_gelq = LAPACKE_sgelq( sgemlq_obj->matrix_layout,sgemlq_obj->m, sgemlq_obj->n,
											sgemlq_obj->A_gelq, sgemlq_obj->lda_gelq, sgemlq_obj->t, sgemlq_obj->tsize);
										
										

    if( sgemlq_obj->info_gelq < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_sgelq is wrong\n", sgemlq_obj->info_gelq );
    }
    if( sgemlq_obj->inforef_gelq < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sgelq is wrong\n", 
        sgemlq_obj->inforef_gelq );
    } 
/* copying input array from gelq to gemlq*/
	/*   Array A*/
	memcpy(sgemlq_obj->A, sgemlq_obj->A_gelq, sgemlq_obj->bufsize_a);
	memcpy(sgemlq_obj->Aref, sgemlq_obj->A_gelqref, sgemlq_obj->bufsize_a);	
	
/*Compute sgemlq's  o/p */
    sgemlq_obj->inforef = sgemlq( sgemlq_obj->matrix_layout, sgemlq_obj->side, sgemlq_obj->trans,
								sgemlq_obj->m, sgemlq_obj->n, sgemlq_obj->k, (const float*)sgemlq_obj->Aref,
								sgemlq_obj->lda, (const float*)sgemlq_obj->tref, sgemlq_obj->tsize, sgemlq_obj->cref, sgemlq_obj->ldc);

    /* Compute libflame's Lapacke o/p  */
    sgemlq_obj->info = LAPACKE_sgemlq( sgemlq_obj->matrix_layout, sgemlq_obj->side, sgemlq_obj->trans,
									sgemlq_obj->m, sgemlq_obj->n, sgemlq_obj->k, (const float*)sgemlq_obj->A,
									sgemlq_obj->lda, (const float*)sgemlq_obj->t, sgemlq_obj->tsize, sgemlq_obj->c, sgemlq_obj->ldc);
    if( sgemlq_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame\
        LAPACKE_sgemlq is wrong\n", sgemlq_obj->info );
    }
    if( sgemlq_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sgemlq is wrong\n", 
        sgemlq_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    sgemlq_obj->diff =  computeDiff_s( sgemlq_obj->bufsize_c, 
                sgemlq_obj->c, sgemlq_obj->cref );

}

TEST_F(sgemlq_test, sgemlq1) {
    EXPECT_NEAR(0.0, sgemlq_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(sgemlq_test, sgemlq2) {
    EXPECT_NEAR(0.0, sgemlq_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(sgemlq_test, sgemlq3) {
    EXPECT_NEAR(0.0, sgemlq_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(sgemlq_test, sgemlq4) {
    EXPECT_NEAR(0.0, sgemlq_obj->diff, LAPACKE_EIG_THRESHOLD);
}

/* Begin double_common_parameters  class definition */
class gemlq_double_parameters{

   public:
	int bufsize_a, bufsize_agelq;
	int bufsize_c;
	void *hModule, *dModule;
	void *bModule, *lModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int k;
	lapack_int tsize;
	double* A, *A_gelq;
	char side;
	char trans;
	lapack_int ldc, lda, lda_gelq;
	/*Output Parameter*/
	double* t, *c;	
	double *Aref, *tref, *cref, *A_gelqref;
	/*Return Values*/	
	lapack_int info, inforef;
	lapack_int info_gelq, inforef_gelq;

   public:
      gemlq_double_parameters (int matrix_layout, char side , char trans, lapack_int m, lapack_int n, lapack_int tsize);
      ~gemlq_double_parameters ();

};

/* Constructor definition  double_common_parameters */
gemlq_double_parameters:: gemlq_double_parameters (int matrix_layout_i, char side_i , char trans_i,lapack_int m_i, lapack_int n_i, lapack_int tsize_i)
{
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	tsize = tsize_i;
	side = side_i;
	trans = trans_i;
	
	if (trans == 'C')
		trans  = 'N' ;
	
	#if LAPACKE_TEST_VERBOSE
		printf(" \n gemlq double: matrix_layout = %d, side :%c, trans:%c, m:%d, n: %d \n", matrix_layout, side, trans, m, n);
	#endif
	if (side == 'L')
	{
		k = m;
		if (matrix_layout == LAPACK_COL_MAJOR)
		{		
			ldc = m;
			lda = k;
			lda_gelq = m;
			bufsize_a = (lda*m);
			bufsize_agelq = lda_gelq*n;
			bufsize_c = (ldc*n);			
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{		
			ldc = n;
			lda = m;
			lda_gelq =n;
			bufsize_a = (lda*k);
			bufsize_agelq = lda_gelq*m;
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
			bufsize_a = (lda*n);
			bufsize_c = (ldc*n);
			lda_gelq = m;
			bufsize_agelq = lda_gelq*n;
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{	ldc = n;
			lda = n;
			lda_gelq =n;
			bufsize_a = (lda*k);
			bufsize_c = (ldc*m);
			bufsize_agelq = lda_gelq*m;
		}else
		{
			EXPECT_TRUE(false) << "matrix_layout invalid";
		}
	}else 
		EXPECT_TRUE(false) << "side is invalid";
	tsize = bufsize_a;

	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_double_buffer_pair(&A_gelq, &A_gelqref, bufsize_agelq);
	lapacke_gtest_alloc_double_buffer_pair(&t, &tref, bufsize_a);	
	lapacke_gtest_alloc_double_buffer_pair(&c, &cref, bufsize_c);
	
	if ((A==NULL) || (Aref==NULL) ||
		(t==NULL) || (tref==NULL) ||
		(c==NULL) || (cref==NULL) ||
		(A_gelq == NULL) || (A_gelqref==NULL))
	{
		EXPECT_FALSE( true) << "gemlq_double_parameters object: malloc error.";
		gemlq_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand( A, Aref, bufsize_a);
	lapacke_gtest_init_double_buffer_pair_rand( A_gelq, A_gelqref, bufsize_agelq);
	lapacke_gtest_init_double_buffer_pair_rand( c, cref, bufsize_c);

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
gemlq_double_parameters :: ~gemlq_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" gemlq_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   gemlq_free();

}
/*  Test fixture class definition */
class dgemlq_test  : public  ::testing::Test {
public:
   gemlq_double_parameters  *dgemlq_obj;
   void SetUp();
   void TearDown () { delete dgemlq_obj;}
};

void dgemlq_test::SetUp(){

	/* LAPACKE dgelq prototype */
    typedef int (*Fptr_NL_LAPACKE_dgelq) (int matrix_layout, lapack_int m,lapack_int n,\ 
											double *A, lapack_int lda, double* t, lapack_int tsize);
	 Fptr_NL_LAPACKE_dgelq dgelq;
	 /* LAPACKE dgemlq prototype */
    typedef int (*Fptr_NL_LAPACKE_dgemlq) (int matrix_layout, char side, char trans,\
                           lapack_int m, lapack_int n, lapack_int k,\
                           const double* a, lapack_int lda,\
                           const double* t, lapack_int tsize,\
                           double* c, lapack_int ldc);

    Fptr_NL_LAPACKE_dgemlq dgemlq;

    dgemlq_obj = new gemlq_double_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].side,
						   eig_paramslist[idx].trans,
						   eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
						   eig_paramslist[idx].m);
						   

    idx = Circular_Increment_Index(idx);

    dgemlq_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgemlq_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgemlq_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgemlq_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*dgemlq library call */
    dgemlq = (Fptr_NL_LAPACKE_dgemlq)dlsym(dgemlq_obj->hModule, "LAPACKE_dgemlq");
    ASSERT_TRUE(dgemlq != NULL) << "failed to get the Netlib LAPACKE_dgemlq symbol";

	/*dgelq library call*/
	dgemlq_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgemlq_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgemlq_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgemlq_obj->lModule != NULL) << "Netlib lapacke handle NULL";

    dgelq = (Fptr_NL_LAPACKE_dgelq)dlsym(dgemlq_obj->lModule, "LAPACKE_dgelq");
    ASSERT_TRUE(dgelq != NULL) << "failed to get the Netlib LAPACKE_dgelq symbol";    
    

    dgemlq_obj->inforef_gelq = dgelq( dgemlq_obj->matrix_layout,dgemlq_obj->m, dgemlq_obj->n,
								dgemlq_obj->A_gelqref, dgemlq_obj->lda_gelq, dgemlq_obj->tref, dgemlq_obj->tsize);

    /* Compute libflame's Lapacke o/p  */
    dgemlq_obj->info_gelq = LAPACKE_dgelq( dgemlq_obj->matrix_layout,dgemlq_obj->m, dgemlq_obj->n,
											dgemlq_obj->A_gelq, dgemlq_obj->lda_gelq, dgemlq_obj->t, dgemlq_obj->tsize);
										
										

    if( dgemlq_obj->info_gelq < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dgelq is wrong\n", dgemlq_obj->info_gelq );
    }
    if( dgemlq_obj->inforef_gelq < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dgelq is wrong\n", 
        dgemlq_obj->inforef_gelq );
    } 
/* copying input array from gelq to gemlq*/
	/*   Array A*/
	memcpy(dgemlq_obj->A, dgemlq_obj->A_gelq, dgemlq_obj->bufsize_a);
	memcpy(dgemlq_obj->Aref, dgemlq_obj->A_gelqref, dgemlq_obj->bufsize_a);	
	
/*Compute dgemlq's  o/p */
    dgemlq_obj->inforef = dgemlq( dgemlq_obj->matrix_layout, dgemlq_obj->side, dgemlq_obj->trans,
								dgemlq_obj->m, dgemlq_obj->n, dgemlq_obj->k, (const double*)dgemlq_obj->Aref,
								dgemlq_obj->lda, (const double*)dgemlq_obj->tref, dgemlq_obj->tsize, dgemlq_obj->cref, dgemlq_obj->ldc);

    /* Compute libflame's Lapacke o/p  */
    dgemlq_obj->info = LAPACKE_dgemlq( dgemlq_obj->matrix_layout, dgemlq_obj->side, dgemlq_obj->trans,
									dgemlq_obj->m, dgemlq_obj->n, dgemlq_obj->k, (const double*)dgemlq_obj->A,
									dgemlq_obj->lda, (const double*)dgemlq_obj->t, dgemlq_obj->tsize, dgemlq_obj->c, dgemlq_obj->ldc);
    if( dgemlq_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame\
        LAPACKE_dgemlq is wrong\n", dgemlq_obj->info );
    }
    if( dgemlq_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dgemlq is wrong\n", 
        dgemlq_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dgemlq_obj->diff =  computeDiff_d( dgemlq_obj->bufsize_c, 
                dgemlq_obj->c, dgemlq_obj->cref );

}

TEST_F(dgemlq_test, dgemlq1) {
    EXPECT_NEAR(0.0, dgemlq_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(dgemlq_test, dgemlq2) {
    EXPECT_NEAR(0.0, dgemlq_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(dgemlq_test, dgemlq3) {
    EXPECT_NEAR(0.0, dgemlq_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(dgemlq_test, dgemlq4) {
    EXPECT_NEAR(0.0, dgemlq_obj->diff, LAPACKE_EIG_THRESHOLD);
}

/* Begin scomplex_common_parameters  class definition */
class gemlq_scomplex_parameters{

   public:
	int bufsize_a, bufsize_agelq;
	int bufsize_c;
	void *hModule, *dModule;
	void *bModule, *lModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int k;
	lapack_int tsize;
	lapack_complex_float* A, *A_gelq;
	char side;
	char trans;
	lapack_int ldc, lda, lda_gelq;
	/*Output Parameter*/
	lapack_complex_float* t, *c;	
	lapack_complex_float *Aref, *tref, *cref, *A_gelqref;
	/*Return Values*/	
	lapack_int info, inforef;
	lapack_int info_gelq, inforef_gelq;

   public:
      gemlq_scomplex_parameters (int matrix_layout, char side , char trans, lapack_int m, lapack_int n, lapack_int tsize);
      ~gemlq_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
gemlq_scomplex_parameters:: gemlq_scomplex_parameters (int matrix_layout_i, char side_i , char trans_i,lapack_int m_i, lapack_int n_i, lapack_int tsize_i)
{
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	tsize = tsize_i;
	side = side_i;
	trans = trans_i;
	
	if (trans == 'T')
		trans  = 'N' ;
	
	#if LAPACKE_TEST_VERBOSE
		printf(" \n gemlq scomplex: matrix_layout = %d, side :%c, trans:%c, m:%d, n: %d \n", matrix_layout, side, trans, m, n);
	#endif
	if (side == 'L')
	{
		k = m;
		if (matrix_layout == LAPACK_COL_MAJOR)
		{		
			ldc = m;
			lda = k;
			lda_gelq = m;
			bufsize_a = (lda*m);
			bufsize_agelq = lda_gelq*n;
			bufsize_c = (ldc*n);			
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{		
			ldc = n;
			lda = m;
			lda_gelq =n;
			bufsize_a = (lda*k);
			bufsize_agelq = lda_gelq*m;
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
			bufsize_a = (lda*n);
			bufsize_c = (ldc*n);
			lda_gelq = m;
			bufsize_agelq = lda_gelq*n;
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{	ldc = n;
			lda = n;
			lda_gelq =n;
			bufsize_a = (lda*k);
			bufsize_c = (ldc*m);
			bufsize_agelq = lda_gelq*m;
		}else
		{
			EXPECT_TRUE(false) << "matrix_layout invalid";
		}
	}else 
		EXPECT_TRUE(false) << "side is invalid";
	tsize = bufsize_a;

	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A_gelq, &A_gelqref, bufsize_agelq);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&t, &tref, bufsize_a);	
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&c, &cref, bufsize_c);
	
	if ((A==NULL) || (Aref==NULL) ||
		(t==NULL) || (tref==NULL) ||
		(c==NULL) || (cref==NULL) ||
		(A_gelq == NULL) || (A_gelqref==NULL))
	{
		EXPECT_FALSE( true) << "gemlq_float_parameters object: malloc error.";
		gemlq_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( A, Aref, bufsize_a);
	lapacke_gtest_init_scomplex_buffer_pair_rand( A_gelq, A_gelqref, bufsize_agelq);
	lapacke_gtest_init_scomplex_buffer_pair_rand( c, cref, bufsize_c);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
gemlq_scomplex_parameters :: ~gemlq_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" gemlq_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   gemlq_free();

}
/*  Test fixture class definition */
class cgemlq_test  : public  ::testing::Test {
public:
   gemlq_scomplex_parameters  *cgemlq_obj;
   void SetUp();
   void TearDown () { delete cgemlq_obj;}
};

void cgemlq_test::SetUp(){

	/* LAPACKE cgelq prototype */
    typedef int (*Fptr_NL_LAPACKE_cgelq) (int matrix_layout, lapack_int m,lapack_int n,\ 
											lapack_complex_float *A, lapack_int lda, lapack_complex_float* t, lapack_int tsize);
	 Fptr_NL_LAPACKE_cgelq cgelq;
	 /* LAPACKE cgemlq prototype */
    typedef int (*Fptr_NL_LAPACKE_cgemlq) (int matrix_layout, char side, char trans,\
                           lapack_int m, lapack_int n, lapack_int k,\
                           const lapack_complex_float* a, lapack_int lda,\
                           const lapack_complex_float* t, lapack_int tsize,\
                           lapack_complex_float* c, lapack_int ldc);

    Fptr_NL_LAPACKE_cgemlq cgemlq;

    cgemlq_obj = new gemlq_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].side,
						   eig_paramslist[idx].trans,
						   eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
						   eig_paramslist[idx].m);
						   

    idx = Circular_Increment_Index(idx);

    cgemlq_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgemlq_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgemlq_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgemlq_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*cgemlq library call */
    cgemlq = (Fptr_NL_LAPACKE_cgemlq)dlsym(cgemlq_obj->hModule, "LAPACKE_cgemlq");
    ASSERT_TRUE(cgemlq != NULL) << "failed to get the Netlib LAPACKE_cgemlq symbol";

	/*cgelq library call*/
	cgemlq_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgemlq_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgemlq_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgemlq_obj->lModule != NULL) << "Netlib lapacke handle NULL";

    cgelq = (Fptr_NL_LAPACKE_cgelq)dlsym(cgemlq_obj->lModule, "LAPACKE_cgelq");
    ASSERT_TRUE(cgelq != NULL) << "failed to get the Netlib LAPACKE_cgelq symbol";    
    

    cgemlq_obj->inforef_gelq = cgelq( cgemlq_obj->matrix_layout,cgemlq_obj->m, cgemlq_obj->n,
								cgemlq_obj->A_gelqref, cgemlq_obj->lda_gelq, cgemlq_obj->tref, cgemlq_obj->tsize);

    /* Compute libflame's Lapacke o/p  */
    cgemlq_obj->info_gelq = LAPACKE_cgelq( cgemlq_obj->matrix_layout,cgemlq_obj->m, cgemlq_obj->n,
											cgemlq_obj->A_gelq, cgemlq_obj->lda_gelq, cgemlq_obj->t, cgemlq_obj->tsize);
										
										

    if( cgemlq_obj->info_gelq < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cgelq is wrong\n", cgemlq_obj->info_gelq );
    }
    if( cgemlq_obj->inforef_gelq < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cgelq is wrong\n", 
        cgemlq_obj->inforef_gelq );
    } 
/* copying input array from gelq to gemlq*/
	/*   Array A*/
	memcpy(cgemlq_obj->A, cgemlq_obj->A_gelq, cgemlq_obj->bufsize_a);
	memcpy(cgemlq_obj->Aref, cgemlq_obj->A_gelqref, cgemlq_obj->bufsize_a);	
	
/*Compute cgemlq's  o/p */
    cgemlq_obj->inforef = cgemlq( cgemlq_obj->matrix_layout, cgemlq_obj->side, cgemlq_obj->trans,
								cgemlq_obj->m, cgemlq_obj->n, cgemlq_obj->k, (const lapack_complex_float*)cgemlq_obj->Aref,
								cgemlq_obj->lda, (const lapack_complex_float*)cgemlq_obj->tref, cgemlq_obj->tsize, cgemlq_obj->cref, cgemlq_obj->ldc);

    /* Compute libflame's Lapacke o/p  */
    cgemlq_obj->info = LAPACKE_cgemlq( cgemlq_obj->matrix_layout, cgemlq_obj->side, cgemlq_obj->trans,
									cgemlq_obj->m, cgemlq_obj->n, cgemlq_obj->k, (const lapack_complex_float*)cgemlq_obj->A,
									cgemlq_obj->lda, (const lapack_complex_float*)cgemlq_obj->t, cgemlq_obj->tsize, cgemlq_obj->c, cgemlq_obj->ldc);
    if( cgemlq_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame\
        LAPACKE_cgemlq is wrong\n", cgemlq_obj->info );
    }
    if( cgemlq_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cgemlq is wrong\n", 
        cgemlq_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cgemlq_obj->diff =  computeDiff_c( cgemlq_obj->bufsize_c, 
                cgemlq_obj->c, cgemlq_obj->cref );

}

TEST_F(cgemlq_test, cgemlq1) {
    EXPECT_NEAR(0.0, cgemlq_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(cgemlq_test, cgemlq2) {
    EXPECT_NEAR(0.0, cgemlq_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(cgemlq_test, cgemlq3) {
    EXPECT_NEAR(0.0, cgemlq_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(cgemlq_test, cgemlq4) {
    EXPECT_NEAR(0.0, cgemlq_obj->diff, LAPACKE_EIG_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class gemlq_dcomplex_parameters{

   public:
	int bufsize_a, bufsize_agelq;
	int bufsize_c;
	void *hModule, *dModule;
	void *bModule, *lModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int k;
	lapack_int tsize;
	lapack_complex_double* A, *A_gelq;
	char side;
	char trans;
	lapack_int ldc, lda, lda_gelq;
	/*Output Parameter*/
	lapack_complex_double* t, *c;	
	lapack_complex_double *Aref, *tref, *cref, *A_gelqref;
	/*Return Values*/	
	lapack_int info, inforef;
	lapack_int info_gelq, inforef_gelq;

   public:
      gemlq_dcomplex_parameters (int matrix_layout, char side , char trans, lapack_int m, lapack_int n, lapack_int tsize);
      ~gemlq_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
gemlq_dcomplex_parameters:: gemlq_dcomplex_parameters (int matrix_layout_i, char side_i , char trans_i,lapack_int m_i, lapack_int n_i, lapack_int tsize_i)
{
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	tsize = tsize_i;
	side = side_i;
	trans = trans_i;
	
	if (trans == 'T')
		trans  = 'N' ;
	
	#if LAPACKE_TEST_VERBOSE
		printf(" \n gemlq dcomplex: matrix_layout = %d, side :%c, trans:%c, m:%d, n: %d \n", matrix_layout, side, trans, m, n);
	#endif
	if (side == 'L')
	{
		k = m;
		if (matrix_layout == LAPACK_COL_MAJOR)
		{		
			ldc = m;
			lda = k;
			lda_gelq = m;
			bufsize_a = (lda*m);
			bufsize_agelq = lda_gelq*n;
			bufsize_c = (ldc*n);			
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{		
			ldc = n;
			lda = m;
			lda_gelq =n;
			bufsize_a = (lda*k);
			bufsize_agelq = lda_gelq*m;
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
			bufsize_a = (lda*n);
			bufsize_c = (ldc*n);
			lda_gelq = m;
			bufsize_agelq = lda_gelq*n;
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{	ldc = n;
			lda = n;
			lda_gelq =n;
			bufsize_a = (lda*k);
			bufsize_c = (ldc*m);
			bufsize_agelq = lda_gelq*m;
		}else
		{
			EXPECT_TRUE(false) << "matrix_layout invalid";
		}
	}else 
		EXPECT_TRUE(false) << "side is invalid";
	tsize = bufsize_a;

	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A_gelq, &A_gelqref, bufsize_agelq);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&t, &tref, bufsize_a);	
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&c, &cref, bufsize_c);
	
	if ((A==NULL) || (Aref==NULL) ||
		(t==NULL) || (tref==NULL) ||
		(c==NULL) || (cref==NULL) ||
		(A_gelq == NULL) || (A_gelqref==NULL))
	{
		EXPECT_FALSE( true) << "gemlq_double_parameters object: malloc error.";
		gemlq_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( A, Aref, bufsize_a);
	lapacke_gtest_init_dcomplex_buffer_pair_rand( A_gelq, A_gelqref, bufsize_agelq);
	lapacke_gtest_init_dcomplex_buffer_pair_rand( c, cref, bufsize_c);

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
gemlq_dcomplex_parameters :: ~gemlq_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" gemlq_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   gemlq_free();

}
/*  Test fixture class definition */
class zgemlq_test  : public  ::testing::Test {
public:
   gemlq_dcomplex_parameters  *zgemlq_obj;
   void SetUp();
   void TearDown () { delete zgemlq_obj;}
};

void zgemlq_test::SetUp(){

	/* LAPACKE zgelq prototype */
    typedef int (*Fptr_NL_LAPACKE_zgelq) (int matrix_layout, lapack_int m,lapack_int n,\ 
											lapack_complex_double *A, lapack_int lda, lapack_complex_double* t, lapack_int tsize);
	 Fptr_NL_LAPACKE_zgelq zgelq;
	 /* LAPACKE zgemlq prototype */
    typedef int (*Fptr_NL_LAPACKE_zgemlq) (int matrix_layout, char side, char trans,\
                           lapack_int m, lapack_int n, lapack_int k,\
                           const lapack_complex_double* a, lapack_int lda,\
                           const lapack_complex_double* t, lapack_int tsize,\
                           lapack_complex_double* c, lapack_int ldc);

    Fptr_NL_LAPACKE_zgemlq zgemlq;

    zgemlq_obj = new gemlq_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].side,
						   eig_paramslist[idx].trans,
						   eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
						   eig_paramslist[idx].m);
						   

    idx = Circular_Increment_Index(idx);

    zgemlq_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgemlq_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgemlq_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgemlq_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*zgemlq library call */
    zgemlq = (Fptr_NL_LAPACKE_zgemlq)dlsym(zgemlq_obj->hModule, "LAPACKE_zgemlq");
    ASSERT_TRUE(zgemlq != NULL) << "failed to get the Netlib LAPACKE_zgemlq symbol";

	/*zgelq library call*/
	zgemlq_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgemlq_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgemlq_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgemlq_obj->lModule != NULL) << "Netlib lapacke handle NULL";

    zgelq = (Fptr_NL_LAPACKE_zgelq)dlsym(zgemlq_obj->lModule, "LAPACKE_zgelq");
    ASSERT_TRUE(zgelq != NULL) << "failed to get the Netlib LAPACKE_zgelq symbol";    
    

    zgemlq_obj->inforef_gelq = zgelq( zgemlq_obj->matrix_layout,zgemlq_obj->m, zgemlq_obj->n,
								zgemlq_obj->A_gelqref, zgemlq_obj->lda_gelq, zgemlq_obj->tref, zgemlq_obj->tsize);

    /* Compute libflame's Lapacke o/p  */
    zgemlq_obj->info_gelq = LAPACKE_zgelq( zgemlq_obj->matrix_layout,zgemlq_obj->m, zgemlq_obj->n,
											zgemlq_obj->A_gelq, zgemlq_obj->lda_gelq, zgemlq_obj->t, zgemlq_obj->tsize);
										
										

    if( zgemlq_obj->info_gelq < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zgelq is wrong\n", zgemlq_obj->info_gelq );
    }
    if( zgemlq_obj->inforef_gelq < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zgelq is wrong\n", 
        zgemlq_obj->inforef_gelq );
    } 
/* copying input array from gelq to gemlq*/
	/*   Array A*/
	memcpy(zgemlq_obj->A, zgemlq_obj->A_gelq, zgemlq_obj->bufsize_a);
	memcpy(zgemlq_obj->Aref, zgemlq_obj->A_gelqref, zgemlq_obj->bufsize_a);	
	
/*Compute zgemlq's  o/p */
    zgemlq_obj->inforef = zgemlq( zgemlq_obj->matrix_layout, zgemlq_obj->side, zgemlq_obj->trans,
								zgemlq_obj->m, zgemlq_obj->n, zgemlq_obj->k, (const lapack_complex_double*)zgemlq_obj->Aref,
								zgemlq_obj->lda, (const lapack_complex_double*)zgemlq_obj->tref, zgemlq_obj->tsize, zgemlq_obj->cref, zgemlq_obj->ldc);

    /* Compute libflame's Lapacke o/p  */
    zgemlq_obj->info = LAPACKE_zgemlq( zgemlq_obj->matrix_layout, zgemlq_obj->side, zgemlq_obj->trans,
									zgemlq_obj->m, zgemlq_obj->n, zgemlq_obj->k, (const lapack_complex_double*)zgemlq_obj->A,
									zgemlq_obj->lda, (const lapack_complex_double*)zgemlq_obj->t, zgemlq_obj->tsize, zgemlq_obj->c, zgemlq_obj->ldc);
    if( zgemlq_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame\
        LAPACKE_zgemlq is wrong\n", zgemlq_obj->info );
    }
    if( zgemlq_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zgemlq is wrong\n", 
        zgemlq_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zgemlq_obj->diff =  computeDiff_z( zgemlq_obj->bufsize_c, 
                zgemlq_obj->c, zgemlq_obj->cref );

}

TEST_F(zgemlq_test, zgemlq1) {
    EXPECT_NEAR(0.0, zgemlq_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zgemlq_test, zgemlq2) {
    EXPECT_NEAR(0.0, zgemlq_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zgemlq_test, zgemlq3) {
    EXPECT_NEAR(0.0, zgemlq_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zgemlq_test, zgemlq4) {
    EXPECT_NEAR(0.0, zgemlq_obj->diff, LAPACKE_EIG_THRESHOLD);
}