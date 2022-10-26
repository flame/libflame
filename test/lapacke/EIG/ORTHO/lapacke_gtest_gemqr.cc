#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define gemqr_free() \
if (t!=NULL)  free(t);\
if (tref!=NULL) free(tref); \
if (c != NULL) free(c); \
if (cref != NULL) free(cref); \
if (A != NULL) free(A); \
if (Aref != NULL) free(Aref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule);\
if( bModule != NULL) dlclose(bModule); \
if( lModule != NULL) dlclose(lModule); 
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin float_common_parameters  class definition */
class gemqr_float_parameters{

   public:
	int bufsize_c;
	int bufsize_a;
	int bufsize_t;
	void *hModule, *dModule;
	void *bModule, *lModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int k;
	float *A;
	char side;
	char trans;
	lapack_int ldc, lda, lda_geqr, tsize;
	/*Output Parameter*/
	float *t, *c;	
	float *tref, *cref, *Aref;
	/*Return Values*/	
	lapack_int info, inforef;
	lapack_int info_geqr, inforef_geqr;

   public:
      gemqr_float_parameters (int matrix_layout, char side , char trans, lapack_int m, lapack_int n);
      ~gemqr_float_parameters ();

};

/* Constructor definition  float_common_parameters */
gemqr_float_parameters:: gemqr_float_parameters (int matrix_layout_i, char side_i , char trans_i,lapack_int m_i, lapack_int n_i)
{
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	side = side_i;
	trans = trans_i;

	if (trans == 'C')
		trans  = 'N' ;
	
	#if LAPACKE_TEST_VERBOSE
	printf(" \n gemqr float: matrix_layout = %d, side :%c, trans:%c, m:%d, n: %d \n", matrix_layout, side, trans, m, n);
	#endif
	if (side == 'L')
	{
		k = m;
		if (matrix_layout == LAPACK_COL_MAJOR)
		{		
			lda = k;
			ldc = m;
			tsize = m;
			lda_geqr = m;
			bufsize_t =  tsize*min(m, n);
			bufsize_a = lda_geqr*n;
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{		
			ldc = n;
			lda = k;
			tsize = m;
			lda_geqr =n;
			bufsize_t = tsize*min(m, n);
			bufsize_a = lda_geqr*m;
		}else
		{
			EXPECT_TRUE(false) << "matrix_layout invalid";
		}
	}else if (side == 'R') {
		k = n;
		if (matrix_layout == LAPACK_COL_MAJOR)
		{	ldc = m;
			lda = k;
			tsize = m;
			bufsize_t =  tsize*min(m, n);
			lda_geqr = m;
			bufsize_a = lda_geqr*n;
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{	ldc = n;
			lda = k;
			tsize = m;
			lda_geqr =n;
			bufsize_t = tsize*min(m, n);
			bufsize_a = lda_geqr*m;
		}else
		{
			EXPECT_TRUE(false) << "matrix_layout invalid";
		}
	}else 
		EXPECT_TRUE(false) << "side is invalid";
	
	bufsize_c = m*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&t, &tref, bufsize_t);	
	lapacke_gtest_alloc_float_buffer_pair(&c, &cref, bufsize_c);
	lapacke_gtest_alloc_float_buffer_pair(&A, &Aref, bufsize_a);
	
	if ((t==NULL) || (tref==NULL) ||
		(c==NULL) || (cref==NULL) ||
		(A==NULL) || (Aref==NULL))
	{
		EXPECT_FALSE( true) << "gemqr_float_parameters object: malloc error.";
		gemqr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand( A, Aref, bufsize_a);
	lapacke_gtest_init_float_buffer_pair_rand( c, cref, bufsize_c);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
gemqr_float_parameters :: ~gemqr_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" gemqr_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   gemqr_free();

}
/*  Test fixture class definition */
class sgemqr_test  : public  ::testing::Test {
public:
   gemqr_float_parameters  *sgemqr_obj;
   void SetUp();
   void TearDown () { delete sgemqr_obj;}
};

void sgemqr_test::SetUp(){

	/* LAPACKE sgeqr prototype */
    typedef int (*Fptr_NL_LAPACKE_sgeqr) (int matrix_layout, lapack_int m, lapack_int n,\
                          float* a, lapack_int lda,\
                          float* t, lapack_int tsize);
	
	 Fptr_NL_LAPACKE_sgeqr sgeqr;
	 /* LAPACKE sgemqr prototype */
    typedef int (*Fptr_NL_LAPACKE_sgemqr) ( int matrix_layout, char side, char trans,\
                           lapack_int m, lapack_int n, lapack_int k,\
                           const float* a, lapack_int lda,\
                           const float* t, lapack_int tsize,\
                           float* c, lapack_int ldc );

    Fptr_NL_LAPACKE_sgemqr sgemqr;

    sgemqr_obj = new gemqr_float_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].side,
						   eig_paramslist[idx].trans,
						   eig_paramslist[idx].m,
                           eig_paramslist[idx].n);
						   

    idx = Circular_Increment_Index(idx);

    sgemqr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgemqr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgemqr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgemqr_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*sgemqr library call */
    sgemqr = (Fptr_NL_LAPACKE_sgemqr)dlsym(sgemqr_obj->hModule, "LAPACKE_sgemqr");
    ASSERT_TRUE(sgemqr != NULL) << "failed to get the Netlib LAPACKE_sgemqr symbol";

	/*sgeqr library call*/
	sgemqr_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgemqr_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgemqr_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgemqr_obj->lModule != NULL) << "Netlib lapacke handle NULL";

    sgeqr = (Fptr_NL_LAPACKE_sgeqr)dlsym(sgemqr_obj->lModule, "LAPACKE_sgeqr");
    ASSERT_TRUE(sgeqr != NULL) << "failed to get the Netlib LAPACKE_sgeqr symbol";    
    

    sgemqr_obj->inforef_geqr = sgeqr( sgemqr_obj->matrix_layout,sgemqr_obj->m, sgemqr_obj->n,\
	sgemqr_obj->Aref, sgemqr_obj->lda_geqr, sgemqr_obj->tref, sgemqr_obj->tsize);

    /* Compute libflame's Lapacke o/p  */
    sgemqr_obj->info_geqr = LAPACKE_sgeqr( sgemqr_obj->matrix_layout,sgemqr_obj->m, sgemqr_obj->n,\
	sgemqr_obj->A, sgemqr_obj->lda_geqr, sgemqr_obj->t, sgemqr_obj->tsize);
										
										

    if( sgemqr_obj->info_geqr < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_sgeqr is wrong\n", sgemqr_obj->info_geqr );
    }
    if( sgemqr_obj->inforef_geqr < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sgeqr is wrong\n", 
        sgemqr_obj->inforef_geqr );
    }  
/*Compute sgemqr's  o/p */
    sgemqr_obj->inforef = sgemqr( sgemqr_obj->matrix_layout, sgemqr_obj->side, sgemqr_obj->trans,
								sgemqr_obj->m, sgemqr_obj->n, sgemqr_obj->k, (const float*)sgemqr_obj->Aref,
								sgemqr_obj->lda, (const float*)sgemqr_obj->tref, sgemqr_obj->tsize, sgemqr_obj->cref, sgemqr_obj->ldc);

    /* Compute libflame's Lapacke o/p  */
    sgemqr_obj->info = LAPACKE_sgemqr( sgemqr_obj->matrix_layout, sgemqr_obj->side, sgemqr_obj->trans,
									sgemqr_obj->m, sgemqr_obj->n, sgemqr_obj->k, (const float*)sgemqr_obj->A,
									sgemqr_obj->lda, (const float*)sgemqr_obj->t, sgemqr_obj->tsize, sgemqr_obj->c, sgemqr_obj->ldc);
    if( sgemqr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame\
        LAPACKE_sgemqr is wrong\n", sgemqr_obj->info );
    }
    if( sgemqr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sgemqr is wrong\n", 
        sgemqr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    sgemqr_obj->diff =  computeDiff_s( sgemqr_obj->bufsize_c, 
                sgemqr_obj->c, sgemqr_obj->cref );

}

TEST_F(sgemqr_test, sgemqr1) {
    EXPECT_NEAR(0.0, sgemqr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgemqr_test, sgemqr2) {
    EXPECT_NEAR(0.0, sgemqr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgemqr_test, sgemqr3) {
    EXPECT_NEAR(0.0, sgemqr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgemqr_test, sgemqr4) {
    EXPECT_NEAR(0.0, sgemqr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin double_common_parameters  class definition */
class gemqr_double_parameters{

   public:
	int bufsize_c;
	int bufsize_a;
	int bufsize_t;
	void *hModule, *dModule;
	void *bModule, *lModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int k;
	double *A;
	char side;
	char trans;
	lapack_int ldc, lda, lda_geqr, tsize;
	/*Output Parameter*/
	double *t, *c;	
	double *tref, *cref, *Aref;
	/*Return Values*/	
	lapack_int info, inforef;
	lapack_int info_geqr, inforef_geqr;

   public:
      gemqr_double_parameters (int matrix_layout, char side , char trans, lapack_int m, lapack_int n);
      ~gemqr_double_parameters ();

};

/* Constructor definition  double_common_parameters */
gemqr_double_parameters:: gemqr_double_parameters (int matrix_layout_i, char side_i , char trans_i,lapack_int m_i, lapack_int n_i)
{
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	side = side_i;
	trans = trans_i;

	if (trans == 'C')
		trans  = 'N' ;
	
	#if LAPACKE_TEST_VERBOSE
	printf(" \n gemqr double: matrix_layout = %d, side :%c, trans:%c, m:%d, n: %d \n", matrix_layout, side, trans, m, n);
	#endif
	if (side == 'L')
	{
		k = m;
		if (matrix_layout == LAPACK_COL_MAJOR)
		{		
			lda = k;
			ldc = m;
			tsize = m;
			lda_geqr = m;
			bufsize_t =  tsize*min(m, n);
			bufsize_a = lda_geqr*n;
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{		
			ldc = n;
			lda = k;
			tsize = m;
			lda_geqr =n;
			bufsize_t = tsize*min(m, n);
			bufsize_a = lda_geqr*m;
		}else
		{
			EXPECT_TRUE(false) << "matrix_layout invalid";
		}
	}else if (side == 'R') {
		k = n;
		if (matrix_layout == LAPACK_COL_MAJOR)
		{	ldc = m;
			lda = k;
			tsize = m;
			bufsize_t =  tsize*min(m, n);
			lda_geqr = m;
			bufsize_a = lda_geqr*n;
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{	ldc = n;
			lda = k;
			tsize = m;
			lda_geqr =n;
			bufsize_t = tsize*min(m, n);
			bufsize_a = lda_geqr*m;
		}else
		{
			EXPECT_TRUE(false) << "matrix_layout invalid";
		}
	}else 
		EXPECT_TRUE(false) << "side is invalid";
	
	bufsize_c = m*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&t, &tref, bufsize_t);	
	lapacke_gtest_alloc_double_buffer_pair(&c, &cref, bufsize_c);
	lapacke_gtest_alloc_double_buffer_pair(&A, &Aref, bufsize_a);
	
	if ((t==NULL) || (tref==NULL) ||
		(c==NULL) || (cref==NULL) ||
		(A==NULL) || (Aref==NULL))
	{
		EXPECT_FALSE( true) << "gemqr_double_parameters object: malloc error.";
		gemqr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand( A, Aref, bufsize_a);
	lapacke_gtest_init_double_buffer_pair_rand( c, cref, bufsize_c);

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
gemqr_double_parameters :: ~gemqr_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" gemqr_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   gemqr_free();

}
/*  Test fixture class definition */
class dgemqr_test  : public  ::testing::Test {
public:
   gemqr_double_parameters  *dgemqr_obj;
   void SetUp();
   void TearDown () { delete dgemqr_obj;}
};

void dgemqr_test::SetUp(){

	/* LAPACKE dgeqr prototype */
    typedef int (*Fptr_NL_LAPACKE_dgeqr) (int matrix_layout, lapack_int m, lapack_int n,\
                          double* a, lapack_int lda,\
                          double* t, lapack_int tsize);
	
	 Fptr_NL_LAPACKE_dgeqr dgeqr;
	 /* LAPACKE dgemqr prototype */
    typedef int (*Fptr_NL_LAPACKE_dgemqr) ( int matrix_layout, char side, char trans,\
                           lapack_int m, lapack_int n, lapack_int k,\
                           const double* a, lapack_int lda,\
                           const double* t, lapack_int tsize,\
                           double* c, lapack_int ldc );

    Fptr_NL_LAPACKE_dgemqr dgemqr;

    dgemqr_obj = new gemqr_double_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].side,
						   eig_paramslist[idx].trans,
						   eig_paramslist[idx].m,
                           eig_paramslist[idx].n);
						   

    idx = Circular_Increment_Index(idx);

    dgemqr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgemqr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgemqr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgemqr_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*dgemqr library call */
    dgemqr = (Fptr_NL_LAPACKE_dgemqr)dlsym(dgemqr_obj->hModule, "LAPACKE_dgemqr");
    ASSERT_TRUE(dgemqr != NULL) << "failed to get the Netlib LAPACKE_dgemqr symbol";

	/*dgeqr library call*/
	dgemqr_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgemqr_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgemqr_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgemqr_obj->lModule != NULL) << "Netlib lapacke handle NULL";

    dgeqr = (Fptr_NL_LAPACKE_dgeqr)dlsym(dgemqr_obj->lModule, "LAPACKE_dgeqr");
    ASSERT_TRUE(dgeqr != NULL) << "failed to get the Netlib LAPACKE_dgeqr symbol";    
    

    dgemqr_obj->inforef_geqr = dgeqr( dgemqr_obj->matrix_layout,dgemqr_obj->m, dgemqr_obj->n,\
	dgemqr_obj->Aref, dgemqr_obj->lda_geqr, dgemqr_obj->tref, dgemqr_obj->tsize);

    /* Compute libflame's Lapacke o/p  */
    dgemqr_obj->info_geqr = LAPACKE_dgeqr( dgemqr_obj->matrix_layout,dgemqr_obj->m, dgemqr_obj->n,\
	dgemqr_obj->A, dgemqr_obj->lda_geqr, dgemqr_obj->t, dgemqr_obj->tsize);
										
										

    if( dgemqr_obj->info_geqr < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dgeqr is wrong\n", dgemqr_obj->info_geqr );
    }
    if( dgemqr_obj->inforef_geqr < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dgeqr is wrong\n", 
        dgemqr_obj->inforef_geqr );
    }  
/*Compute dgemqr's  o/p */
    dgemqr_obj->inforef = dgemqr( dgemqr_obj->matrix_layout, dgemqr_obj->side, dgemqr_obj->trans,
								dgemqr_obj->m, dgemqr_obj->n, dgemqr_obj->k, (const double*)dgemqr_obj->Aref,
								dgemqr_obj->lda, (const double*)dgemqr_obj->tref, dgemqr_obj->tsize, dgemqr_obj->cref, dgemqr_obj->ldc);

    /* Compute libflame's Lapacke o/p  */
    dgemqr_obj->info = LAPACKE_dgemqr( dgemqr_obj->matrix_layout, dgemqr_obj->side, dgemqr_obj->trans,
									dgemqr_obj->m, dgemqr_obj->n, dgemqr_obj->k, (const double*)dgemqr_obj->A,
									dgemqr_obj->lda, (const double*)dgemqr_obj->t, dgemqr_obj->tsize, dgemqr_obj->c, dgemqr_obj->ldc);
    if( dgemqr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame\
        LAPACKE_dgemqr is wrong\n", dgemqr_obj->info );
    }
    if( dgemqr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dgemqr is wrong\n", 
        dgemqr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dgemqr_obj->diff =  computeDiff_d( dgemqr_obj->bufsize_c, 
                dgemqr_obj->c, dgemqr_obj->cref );

}

TEST_F(dgemqr_test, dgemqr1) {
    EXPECT_NEAR(0.0, dgemqr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgemqr_test, dgemqr2) {
    EXPECT_NEAR(0.0, dgemqr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgemqr_test, dgemqr3) {
    EXPECT_NEAR(0.0, dgemqr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgemqr_test, dgemqr4) {
    EXPECT_NEAR(0.0, dgemqr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin scomplex_common_parameters  class definition */
class gemqr_scomplex_parameters{

   public:
	int bufsize_c;
	int bufsize_a;
	int bufsize_t;
	void *hModule, *dModule;
	void *bModule, *lModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int k;
	lapack_complex_float *A;
	char side;
	char trans;
	lapack_int ldc, lda, lda_geqr, tsize;
	/*Output Parameter*/
	lapack_complex_float *t, *c;	
	lapack_complex_float *tref, *cref, *Aref;
	/*Return Values*/	
	lapack_int info, inforef;
	lapack_int info_geqr, inforef_geqr;

   public:
      gemqr_scomplex_parameters (int matrix_layout, char side , char trans, lapack_int m, lapack_int n);
      ~gemqr_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
gemqr_scomplex_parameters:: gemqr_scomplex_parameters (int matrix_layout_i, char side_i , char trans_i,lapack_int m_i, lapack_int n_i)
{
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	side = side_i;
	trans = trans_i;

	if (trans == 'T')
		trans  = 'N' ;
	
	#if LAPACKE_TEST_VERBOSE
	printf(" \n gemqr scomplex: matrix_layout = %d, side :%c, trans:%c, m:%d, n: %d \n", matrix_layout, side, trans, m, n);
	#endif
	if (side == 'L')
	{
		k = m;
		if (matrix_layout == LAPACK_COL_MAJOR)
		{		
			lda = k;
			ldc = m;
			tsize = m;
			lda_geqr = m;
			bufsize_t =  tsize*min(m, n);
			bufsize_a = lda_geqr*n;
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{		
			ldc = n;
			lda = k;
			tsize = m;
			lda_geqr =n;
			bufsize_t = tsize*min(m, n);
			bufsize_a = lda_geqr*m;
		}else
		{
			EXPECT_TRUE(false) << "matrix_layout invalid";
		}
	}else if (side == 'R') {
		k = n;
		if (matrix_layout == LAPACK_COL_MAJOR)
		{	ldc = m;
			lda = k;
			tsize = m;
			bufsize_t =  tsize*min(m, n);
			lda_geqr = m;
			bufsize_a = lda_geqr*n;
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{	ldc = n;
			lda = k;
			tsize = m;
			lda_geqr =n;
			bufsize_t = tsize*min(m, n);
			bufsize_a = lda_geqr*m;
		}else
		{
			EXPECT_TRUE(false) << "matrix_layout invalid";
		}
	}else 
		EXPECT_TRUE(false) << "side is invalid";
	
	bufsize_c = m*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&t, &tref, bufsize_t);	
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&c, &cref, bufsize_c);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize_a);
	
	if ((t==NULL) || (tref==NULL) ||
		(c==NULL) || (cref==NULL) ||
		(A==NULL) || (Aref==NULL))
	{
		EXPECT_FALSE( true) << "gemqr_float_parameters object: malloc error.";
		gemqr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( A, Aref, bufsize_a);
	lapacke_gtest_init_scomplex_buffer_pair_rand( c, cref, bufsize_c);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
gemqr_scomplex_parameters :: ~gemqr_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" gemqr_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   gemqr_free();

}
/*  Test fixture class definition */
class cgemqr_test  : public  ::testing::Test {
public:
   gemqr_scomplex_parameters  *cgemqr_obj;
   void SetUp();
   void TearDown () { delete cgemqr_obj;}
};

void cgemqr_test::SetUp(){

	/* LAPACKE cgeqr prototype */
    typedef int (*Fptr_NL_LAPACKE_cgeqr) (int matrix_layout, lapack_int m, lapack_int n,\
                          lapack_complex_float* a, lapack_int lda,\
                          lapack_complex_float* t, lapack_int tsize);
	
	 Fptr_NL_LAPACKE_cgeqr cgeqr;
	 /* LAPACKE cgemqr prototype */
    typedef int (*Fptr_NL_LAPACKE_cgemqr) ( int matrix_layout, char side, char trans,\
                           lapack_int m, lapack_int n, lapack_int k,\
                           const lapack_complex_float* a, lapack_int lda,\
                           const lapack_complex_float* t, lapack_int tsize,\
                           lapack_complex_float* c, lapack_int ldc );

    Fptr_NL_LAPACKE_cgemqr cgemqr;

    cgemqr_obj = new gemqr_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].side,
						   eig_paramslist[idx].trans,
						   eig_paramslist[idx].m,
                           eig_paramslist[idx].n);
						   

    idx = Circular_Increment_Index(idx);

    cgemqr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgemqr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgemqr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgemqr_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*cgemqr library call */
    cgemqr = (Fptr_NL_LAPACKE_cgemqr)dlsym(cgemqr_obj->hModule, "LAPACKE_cgemqr");
    ASSERT_TRUE(cgemqr != NULL) << "failed to get the Netlib LAPACKE_cgemqr symbol";

	/*cgeqr library call*/
	cgemqr_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgemqr_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgemqr_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgemqr_obj->lModule != NULL) << "Netlib lapacke handle NULL";

    cgeqr = (Fptr_NL_LAPACKE_cgeqr)dlsym(cgemqr_obj->lModule, "LAPACKE_cgeqr");
    ASSERT_TRUE(cgeqr != NULL) << "failed to get the Netlib LAPACKE_cgeqr symbol";    
    

    cgemqr_obj->inforef_geqr = cgeqr( cgemqr_obj->matrix_layout,cgemqr_obj->m, cgemqr_obj->n,\
	cgemqr_obj->Aref, cgemqr_obj->lda_geqr, cgemqr_obj->tref, cgemqr_obj->tsize);

    /* Compute libflame's Lapacke o/p  */
    cgemqr_obj->info_geqr = LAPACKE_cgeqr( cgemqr_obj->matrix_layout,cgemqr_obj->m, cgemqr_obj->n,\
	cgemqr_obj->A, cgemqr_obj->lda_geqr, cgemqr_obj->t, cgemqr_obj->tsize);
										
										

    if( cgemqr_obj->info_geqr < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cgeqr is wrong\n", cgemqr_obj->info_geqr );
    }
    if( cgemqr_obj->inforef_geqr < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cgeqr is wrong\n", 
        cgemqr_obj->inforef_geqr );
    }  
/*Compute cgemqr's  o/p */
    cgemqr_obj->inforef = cgemqr( cgemqr_obj->matrix_layout, cgemqr_obj->side, cgemqr_obj->trans,
								cgemqr_obj->m, cgemqr_obj->n, cgemqr_obj->k, (const lapack_complex_float*)cgemqr_obj->Aref,
								cgemqr_obj->lda, (const lapack_complex_float*)cgemqr_obj->tref, cgemqr_obj->tsize, cgemqr_obj->cref, cgemqr_obj->ldc);

    /* Compute libflame's Lapacke o/p  */
    cgemqr_obj->info = LAPACKE_cgemqr( cgemqr_obj->matrix_layout, cgemqr_obj->side, cgemqr_obj->trans,
									cgemqr_obj->m, cgemqr_obj->n, cgemqr_obj->k, (const lapack_complex_float*)cgemqr_obj->A,
									cgemqr_obj->lda, (const lapack_complex_float*)cgemqr_obj->t, cgemqr_obj->tsize, cgemqr_obj->c, cgemqr_obj->ldc);
    if( cgemqr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame\
        LAPACKE_cgemqr is wrong\n", cgemqr_obj->info );
    }
    if( cgemqr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cgemqr is wrong\n", 
        cgemqr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cgemqr_obj->diff =  computeDiff_c( cgemqr_obj->bufsize_c, 
                cgemqr_obj->c, cgemqr_obj->cref );

}

TEST_F(cgemqr_test, cgemqr1) {
    EXPECT_NEAR(0.0, cgemqr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgemqr_test, cgemqr2) {
    EXPECT_NEAR(0.0, cgemqr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgemqr_test, cgemqr3) {
    EXPECT_NEAR(0.0, cgemqr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgemqr_test, cgemqr4) {
    EXPECT_NEAR(0.0, cgemqr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class gemqr_dcomplex_parameters{

   public:
	int bufsize_c;
	int bufsize_a;
	int bufsize_t;
	void *hModule, *dModule;
	void *bModule, *lModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int k;
	lapack_complex_double *A;
	char side;
	char trans;
	lapack_int ldc, lda, lda_geqr, tsize;
	/*Output Parameter*/
	lapack_complex_double *t, *c;	
	lapack_complex_double *tref, *cref, *Aref;
	/*Return Values*/	
	lapack_int info, inforef;
	lapack_int info_geqr, inforef_geqr;

   public:
      gemqr_dcomplex_parameters (int matrix_layout, char side , char trans, lapack_int m, lapack_int n);
      ~gemqr_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
gemqr_dcomplex_parameters:: gemqr_dcomplex_parameters (int matrix_layout_i, char side_i , char trans_i,lapack_int m_i, lapack_int n_i)
{
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	side = side_i;
	trans = trans_i;

	if (trans == 'T')
		trans  = 'N' ;
	
	#if LAPACKE_TEST_VERBOSE
	printf(" \n gemqr dcomplex: matrix_layout = %d, side :%c, trans:%c, m:%d, n: %d \n", matrix_layout, side, trans, m, n);
	#endif
	if (side == 'L')
	{
		k = m;
		if (matrix_layout == LAPACK_COL_MAJOR)
		{		
			lda = k;
			ldc = m;
			tsize = m;
			lda_geqr = m;
			bufsize_t =  tsize*min(m, n);
			bufsize_a = lda_geqr*n;
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{		
			ldc = n;
			lda = k;
			tsize = m;
			lda_geqr =n;
			bufsize_t = tsize*min(m, n);
			bufsize_a = lda_geqr*m;
		}else
		{
			EXPECT_TRUE(false) << "matrix_layout invalid";
		}
	}else if (side == 'R') {
		k = n;
		if (matrix_layout == LAPACK_COL_MAJOR)
		{	ldc = m;
			lda = k;
			tsize = m;
			bufsize_t =  tsize*min(m, n);
			lda_geqr = m;
			bufsize_a = lda_geqr*n;
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{	ldc = n;
			lda = k;
			tsize = m;
			lda_geqr =n;
			bufsize_t = tsize*min(m, n);
			bufsize_a = lda_geqr*m;
		}else
		{
			EXPECT_TRUE(false) << "matrix_layout invalid";
		}
	}else 
		EXPECT_TRUE(false) << "side is invalid";
	
	bufsize_c = m*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&t, &tref, bufsize_t);	
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&c, &cref, bufsize_c);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize_a);
	
	if ((t==NULL) || (tref==NULL) ||
		(c==NULL) || (cref==NULL) ||
		(A==NULL) || (Aref==NULL))
	{
		EXPECT_FALSE( true) << "gemqr_double_parameters object: malloc error.";
		gemqr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( A, Aref, bufsize_a);
	lapacke_gtest_init_dcomplex_buffer_pair_rand( c, cref, bufsize_c);

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
gemqr_dcomplex_parameters :: ~gemqr_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" gemqr_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   gemqr_free();

}
/*  Test fixture class definition */
class zgemqr_test  : public  ::testing::Test {
public:
   gemqr_dcomplex_parameters  *zgemqr_obj;
   void SetUp();
   void TearDown () { delete zgemqr_obj;}
};

void zgemqr_test::SetUp(){

	/* LAPACKE zgeqr prototype */
    typedef int (*Fptr_NL_LAPACKE_zgeqr) (int matrix_layout, lapack_int m, lapack_int n,\
                          lapack_complex_double* a, lapack_int lda,\
                          lapack_complex_double* t, lapack_int tsize);
	
	 Fptr_NL_LAPACKE_zgeqr zgeqr;
	 /* LAPACKE zgemqr prototype */
    typedef int (*Fptr_NL_LAPACKE_zgemqr) ( int matrix_layout, char side, char trans,\
                           lapack_int m, lapack_int n, lapack_int k,\
                           const lapack_complex_double* a, lapack_int lda,\
                           const lapack_complex_double* t, lapack_int tsize,\
                           lapack_complex_double* c, lapack_int ldc );

    Fptr_NL_LAPACKE_zgemqr zgemqr;

    zgemqr_obj = new gemqr_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].side,
						   eig_paramslist[idx].trans,
						   eig_paramslist[idx].m,
                           eig_paramslist[idx].n);
						   

    idx = Circular_Increment_Index(idx);

    zgemqr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgemqr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgemqr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgemqr_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*zgemqr library call */
    zgemqr = (Fptr_NL_LAPACKE_zgemqr)dlsym(zgemqr_obj->hModule, "LAPACKE_zgemqr");
    ASSERT_TRUE(zgemqr != NULL) << "failed to get the Netlib LAPACKE_zgemqr symbol";

	/*zgeqr library call*/
	zgemqr_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgemqr_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgemqr_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgemqr_obj->lModule != NULL) << "Netlib lapacke handle NULL";

    zgeqr = (Fptr_NL_LAPACKE_zgeqr)dlsym(zgemqr_obj->lModule, "LAPACKE_zgeqr");
    ASSERT_TRUE(zgeqr != NULL) << "failed to get the Netlib LAPACKE_zgeqr symbol";    
    

    zgemqr_obj->inforef_geqr = zgeqr( zgemqr_obj->matrix_layout,zgemqr_obj->m, zgemqr_obj->n,\
	zgemqr_obj->Aref, zgemqr_obj->lda_geqr, zgemqr_obj->tref, zgemqr_obj->tsize);

    /* Compute libflame's Lapacke o/p  */
    zgemqr_obj->info_geqr = LAPACKE_zgeqr( zgemqr_obj->matrix_layout,zgemqr_obj->m, zgemqr_obj->n,\
	zgemqr_obj->A, zgemqr_obj->lda_geqr, zgemqr_obj->t, zgemqr_obj->tsize);
										
										

    if( zgemqr_obj->info_geqr < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zgeqr is wrong\n", zgemqr_obj->info_geqr );
    }
    if( zgemqr_obj->inforef_geqr < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zgeqr is wrong\n", 
        zgemqr_obj->inforef_geqr );
    }  
/*Compute zgemqr's  o/p */
    zgemqr_obj->inforef = zgemqr( zgemqr_obj->matrix_layout, zgemqr_obj->side, zgemqr_obj->trans,
								zgemqr_obj->m, zgemqr_obj->n, zgemqr_obj->k, (const lapack_complex_double*)zgemqr_obj->Aref,
								zgemqr_obj->lda, (const lapack_complex_double*)zgemqr_obj->tref, zgemqr_obj->tsize, zgemqr_obj->cref, zgemqr_obj->ldc);

    /* Compute libflame's Lapacke o/p  */
    zgemqr_obj->info = LAPACKE_zgemqr( zgemqr_obj->matrix_layout, zgemqr_obj->side, zgemqr_obj->trans,
									zgemqr_obj->m, zgemqr_obj->n, zgemqr_obj->k, (const lapack_complex_double*)zgemqr_obj->A,
									zgemqr_obj->lda, (const lapack_complex_double*)zgemqr_obj->t, zgemqr_obj->tsize, zgemqr_obj->c, zgemqr_obj->ldc);
    if( zgemqr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame\
        LAPACKE_zgemqr is wrong\n", zgemqr_obj->info );
    }
    if( zgemqr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zgemqr is wrong\n", 
        zgemqr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zgemqr_obj->diff =  computeDiff_z( zgemqr_obj->bufsize_c, 
                zgemqr_obj->c, zgemqr_obj->cref );

}

TEST_F(zgemqr_test, zgemqr1) {
    EXPECT_NEAR(0.0, zgemqr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgemqr_test, zgemqr2) {
    EXPECT_NEAR(0.0, zgemqr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgemqr_test, zgemqr3) {
    EXPECT_NEAR(0.0, zgemqr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgemqr_test, zgemqr4) {
    EXPECT_NEAR(0.0, zgemqr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}