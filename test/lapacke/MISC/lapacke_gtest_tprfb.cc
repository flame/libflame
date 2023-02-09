#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"

#define tprfb_free() \
if (v!=NULL)    free(v); \
if (vref!=NULL) free(vref);\
if (t!=NULL)  free(t);\
if (tref!=NULL) free(tref); \
if (a != NULL) free(a); \
if (aref != NULL) free(aref); \
if (b != NULL) free(b); \
if (bref != NULL) free(bref); \
if (A_tpqrt != NULL) free(A_tpqrt);\
if (A_tpqrtref != NULL) free(A_tpqrtref);\
if (b_tpqrt != NULL) free(b_tpqrt); \
if (b_tpqrtref != NULL) free(b_tpqrtref); \
if (t_tpqrt!=NULL)  free(t_tpqrt);\
if (t_tpqrtref!=NULL) free(t_tpqrtref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule);\
if( bModule != NULL) dlclose(bModule); \
if( lModule != NULL) dlclose(lModule); \

	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
/* Begin float_common_parameters  class definition */
class tprfb_float_parameters{

   public:
	int bufsize_t;
	int bufsize_v;
	int bufsize_a;
	int bufsize_b;
	int bufsize_t_tpqrt;
	void *hModule, *dModule;
	void *bModule, *lModule;
	float diff;
	float diff_b;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int k;
	lapack_int l;
	lapack_int lda, ldv, ldt, ldb;
	lapack_int ldt_tpqrt;
	float* v, *A_tpqrt, *b_tpqrt, *t_tpqrt;
	char side, trans, direct, storev;
	/*Output Parameter*/
	float* t, *a, *b;
	float *A_tpqrtref, *b_tpqrtref, *t_tpqrtref;
	float *vref, *tref, *aref, *bref;
	/*Return Values*/	
	lapack_int info, inforef;
	lapack_int info_tpqrt, inforef_tpqrt;

   public:
      tprfb_float_parameters (int matrix_layout, char side, char trans,  lapack_int m, lapack_int n, lapack_int k, lapack_int l, char storev);
      ~tprfb_float_parameters ();

};

/* Constructor definition  float_common_parameters */
tprfb_float_parameters:: tprfb_float_parameters (int matrix_layout_i, char side_i, char trans_i , lapack_int m_i, lapack_int n_i, lapack_int k_i, lapack_int l_i, char storev_i)
{
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	k = k_i;
	l = l_i;
	ldt = k;
	side = side_i;
	trans = trans_i;
	storev = storev_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n tprfb float: matrix_layout = %d,  m:%d, n: %d, trans:%c, l:%d \n", matrix_layout, side, m, n, trans, l);
	#endif

	if (matrix_layout == LAPACK_COL_MAJOR)
	{	/* checking sizes with storev*/	
		if (storev == 'C')
		{
			if ( side == 'L')
			{	
				ldv = m;
				bufsize_v = ldv*k;
				lda = k;
				bufsize_a = lda*n;
			} else if (side == 'R')
			{
				ldv = n;
				lda = m;
				bufsize_v = ldv*k;
				bufsize_a = lda*k;
			}
			/* checking sizes with storev*/	
		} else if (storev == 'R')
		{
			if ( side == 'L')
			{	
				ldv = k;
				bufsize_v = ldv*m;
				lda = k;
				bufsize_a = lda*n;
			} else if (side == 'R')
			{
				ldv = k;
				bufsize_v = ldv*n;
				lda = m;
				bufsize_a = lda*k;
			}
		}
		direct = 'F';
		ldb = m;
		bufsize_b = ldb*n;
		
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{		
		if (storev == 'R')
		{
			if ( side == 'L')
			{	
				ldv = m;
				lda =  n;
				bufsize_v = ldv*k;
				bufsize_a = lda*k;
			} else if (side == 'R')
			{
				ldv = n;
				lda = k;
				bufsize_v = ldv*k;
				bufsize_a = lda*m;
			}
			/* checking sizes with storev*/	
		}else if (storev == 'C')
		{	if ( side == 'L')
			{	
				ldv = k;
				lda =  n;
				bufsize_v = ldv*m;
				bufsize_a = lda*k;	
			} else if (side == 'R')
			{
				ldv = k;
				lda = k;
				bufsize_v = ldv*n;
				bufsize_a = lda*m;
			}
		}
		direct = 'B';
		ldb = n;
		bufsize_b = ldb*m;
	}else
	{
		EXPECT_TRUE(false) << "matrix_layout invalid";
	}
	
	bufsize_t = ldt*k;
	bufsize_t_tpqrt = n*n; 
	
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&v, &vref, bufsize_v);
	lapacke_gtest_alloc_float_buffer_pair(&t, &tref, bufsize_t);	
	lapacke_gtest_alloc_float_buffer_pair(&b, &bref, bufsize_b);
	lapacke_gtest_alloc_float_buffer_pair(&a, &aref, bufsize_a);
	lapacke_gtest_alloc_float_buffer_pair(&A_tpqrt, &A_tpqrtref, (n*n));
	lapacke_gtest_alloc_float_buffer_pair(&b_tpqrt, &b_tpqrtref, bufsize_b);
	lapacke_gtest_alloc_float_buffer_pair(&t_tpqrt, &t_tpqrtref, bufsize_t_tpqrt);
	
	if ((v==NULL) || (vref==NULL) ||
		(t==NULL) || (tref==NULL) ||
		(a==NULL) || (aref==NULL) ||
		(b==NULL) || (bref==NULL) ||
		(A_tpqrt ==NULL) || (A_tpqrtref==NULL) ||
		(b_tpqrt == NULL) || (b_tpqrtref == NULL) ||
		(t_tpqrt == NULL) || (t_tpqrtref == NULL))
	{
		EXPECT_FALSE( true) << "tprfb_float_parameters object: malloc error.";
		tprfb_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand( v, vref, bufsize_v);
	lapacke_gtest_init_float_buffer_pair_rand( b, bref, bufsize_b);
	lapacke_gtest_init_float_buffer_pair_rand( a, aref, bufsize_a);
	lapacke_gtest_init_float_buffer_pair_rand_custom_matrix( t, tref, ldt, k, 'S');
	lapacke_gtest_init_float_buffer_pair_rand( A_tpqrt, A_tpqrtref, (n*n));
	lapacke_gtest_init_float_buffer_pair_rand( b_tpqrt, b_tpqrtref, bufsize_b);


} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
tprfb_float_parameters :: ~tprfb_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" tprfb_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   tprfb_free();

}
/*  Test fixture class definition */
class stprfb_test  : public  ::testing::Test {
public:
   tprfb_float_parameters  *stprfb_obj;
   void SetUp();
   void TearDown () { delete stprfb_obj;}
};

void stprfb_test::SetUp(){
	/* LAPACKE stpqrt prototype */
	typedef int (*Fptr_NL_LAPACKE_stpqrt) (int matrix_layout, lapack_int m, lapack_int n,\
	lapack_int l, lapack_int nb, float* a, lapack_int lda, float* b, lapack_int ldb, float* t, lapack_int ldt);
	
	Fptr_NL_LAPACKE_stpqrt stpqrt;
	
	 /* LAPACKE stprfb prototype */
    typedef int (*Fptr_NL_LAPACKE_stprfb) (int matrix_layout, char side, char trans, char direct, char storev, \
	lapack_int m, lapack_int n, lapack_int k, lapack_int l, const float * v, lapack_int ldv, const float * t, \
	lapack_int ldt, float * a, lapack_int lda, float * b, lapack_int ldb);

    Fptr_NL_LAPACKE_stprfb stprfb;

    stprfb_obj = new tprfb_float_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].side,
						   eig_paramslist[idx].trans,
						   eig_paramslist[idx].m,
						   eig_paramslist[idx].n,
                           eig_paramslist[idx].m,
						   eig_paramslist[idx].k,
						   eig_paramslist[idx].storev);
						   

    idx = Circular_Increment_Index(idx);

    stprfb_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    stprfb_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(stprfb_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(stprfb_obj->hModule != NULL) << "Netlib lapacke handle NULL";
	
	/*stprfb library call */
    stprfb = (Fptr_NL_LAPACKE_stprfb)dlsym(stprfb_obj->hModule, "LAPACKE_stprfb");
    ASSERT_TRUE(stprfb != NULL) << "failed to get the Netlib LAPACKE_stprfb symbol";
	
	/*stpqrt library call*/
	stprfb_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    stprfb_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(stprfb_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(stprfb_obj->lModule != NULL) << "Netlib lapacke handle NULL";
	
	stpqrt = (Fptr_NL_LAPACKE_stpqrt)dlsym(stprfb_obj->hModule, "LAPACKE_stpqrt");
    ASSERT_TRUE(stpqrt != NULL) << "failed to get the Netlib LAPACKE_stpqrt symbol";
    

    stprfb_obj->inforef_tpqrt = stpqrt( stprfb_obj->matrix_layout, stprfb_obj->m,
								stprfb_obj->n, min(stprfb_obj->m, stprfb_obj->n), stprfb_obj->n, stprfb_obj->A_tpqrtref,
								stprfb_obj->n, stprfb_obj->b_tpqrtref, stprfb_obj->ldb, stprfb_obj->t_tpqrtref, stprfb_obj->n);

    /* Compute libflame's Lapacke o/p  */
    stprfb_obj->info_tpqrt = LAPACKE_stpqrt( stprfb_obj->matrix_layout, stprfb_obj->m, stprfb_obj->n,
											min(stprfb_obj->m, stprfb_obj->n), stprfb_obj->n, stprfb_obj->A_tpqrt, 
											stprfb_obj->n, stprfb_obj->b_tpqrt,	stprfb_obj->ldb, stprfb_obj->t_tpqrt, stprfb_obj->n);

    if( stprfb_obj->info_tpqrt < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_stpqrt is wrong\n", stprfb_obj->info );
    }
    if( stprfb_obj->inforef_tpqrt < 0 ) {
        printf( "The i:%d th argument with Netlib stpqrt is wrong\n", 
        stprfb_obj->inforef );
    }
	/*copy */
	
  
/*Compute stprfb's  o/p */
    stprfb_obj->inforef = stprfb( stprfb_obj->matrix_layout, stprfb_obj->side, stprfb_obj->trans, stprfb_obj->direct, stprfb_obj->storev,\
	stprfb_obj->m, stprfb_obj->n,stprfb_obj->k, stprfb_obj->l, (const float*)stprfb_obj->b_tpqrtref,stprfb_obj->ldv, (const float*)stprfb_obj->t_tpqrtref, stprfb_obj->ldt,
	stprfb_obj->aref, stprfb_obj->lda, stprfb_obj->b, stprfb_obj->ldb);

    /* Compute libflame's Lapacke o/p  */
    stprfb_obj->info = LAPACKE_stprfb( stprfb_obj->matrix_layout, stprfb_obj->side, stprfb_obj->trans, stprfb_obj->direct, stprfb_obj->storev,\
	stprfb_obj->m, stprfb_obj->n, stprfb_obj->k, stprfb_obj->l, (const float*)stprfb_obj->b_tpqrt, stprfb_obj->ldv, (const float*)stprfb_obj->t_tpqrt, stprfb_obj->ldt, \
	stprfb_obj->a, stprfb_obj->lda, stprfb_obj->b, stprfb_obj->ldb);
	
    if( stprfb_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame\
        LAPACKE_stprfb is wrong\n", stprfb_obj->info );
    }
    if( stprfb_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_stprfb is wrong\n", 
        stprfb_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    stprfb_obj->diff =  computeDiff_s( stprfb_obj->bufsize_a, 
                stprfb_obj->a, stprfb_obj->aref );
				
	stprfb_obj->diff_b =  computeDiff_s( stprfb_obj->bufsize_b,
						stprfb_obj->b, stprfb_obj->bref );
	

}

TEST_F(stprfb_test, stprfb1) {
    EXPECT_NEAR(0.0, stprfb_obj->diff, LAPACKE_EIG_THRESHOLD);
//	EXPECT_NEAR(0.0, stprfb_obj->diff_b, LAPACKE_EIG_THRESHOLD);
}

TEST_F(stprfb_test, stprfb2) {
    EXPECT_NEAR(0.0, stprfb_obj->diff, LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, stprfb_obj->diff_b, LAPACKE_EIG_THRESHOLD);
}

TEST_F(stprfb_test, stprfb3) {
    EXPECT_NEAR(0.0, stprfb_obj->diff, LAPACKE_EIG_THRESHOLD);
//	EXPECT_NEAR(0.0, stprfb_obj->diff_b, LAPACKE_EIG_THRESHOLD);
}

TEST_F(stprfb_test, stprfb4) {
    EXPECT_NEAR(0.0, stprfb_obj->diff, LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, stprfb_obj->diff_b, LAPACKE_EIG_THRESHOLD);
}

/* Begin double_common_parameters  class definition */
class tprfb_double_parameters{

   public:
	int bufsize_t;
	int bufsize_v;
	int bufsize_a;
	int bufsize_b;
	int bufsize_t_tpqrt;
	void *hModule, *dModule;
	void *bModule, *lModule;
	double diff;
	double diff_b;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int k;
	lapack_int l;
	lapack_int lda, ldv, ldt, ldb;
	lapack_int ldt_tpqrt;
	double* v, *A_tpqrt, *b_tpqrt, *t_tpqrt;
	char side, trans, direct, storev;
	/*Output Parameter*/
	double* t, *a, *b;
	double *A_tpqrtref, *b_tpqrtref, *t_tpqrtref;
	double *vref, *tref, *aref, *bref;
	/*Return Values*/	
	lapack_int info, inforef;
	lapack_int info_tpqrt, inforef_tpqrt;

   public:
      tprfb_double_parameters (int matrix_layout, char side, char trans,  lapack_int m, lapack_int n, lapack_int k, lapack_int l, char storev);
      ~tprfb_double_parameters ();

};

/* Constructor definition  double_common_parameters */
tprfb_double_parameters:: tprfb_double_parameters (int matrix_layout_i, char side_i, char trans_i , lapack_int m_i, lapack_int n_i, lapack_int k_i, lapack_int l_i, char storev_i)
{
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	k = k_i;
	l = l_i;
	ldt = k;
	side = side_i;
	trans = trans_i;
	storev = storev_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n tprfb double: matrix_layout = %d,  m:%d, n: %d, trans:%c, l:%d \n", matrix_layout, side, m, n, trans, l);
	#endif

	if (matrix_layout == LAPACK_COL_MAJOR)
	{	/* checking sizes with storev*/	
		if (storev == 'C')
		{
			if ( side == 'L')
			{	
				ldv = m;
				bufsize_v = ldv*k;
				lda = k;
				bufsize_a = lda*n;
			} else if (side == 'R')
			{
				ldv = n;
				lda = m;
				bufsize_v = ldv*k;
				bufsize_a = lda*k;
			}
			/* checking sizes with storev*/	
		} else if (storev == 'R')
		{
			if ( side == 'L')
			{	
				ldv = k;
				bufsize_v = ldv*m;
				lda = k;
				bufsize_a = lda*n;
			} else if (side == 'R')
			{
				ldv = k;
				bufsize_v = ldv*n;
				lda = m;
				bufsize_a = lda*k;
			}
		}
		direct = 'F';
		ldb = m;
		bufsize_b = ldb*n;
		
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{		
		if (storev == 'R')
		{
			if ( side == 'L')
			{	
				ldv = m;
				lda =  n;
				bufsize_v = ldv*k;
				bufsize_a = lda*k;
			} else if (side == 'R')
			{
				ldv = n;
				lda = k;
				bufsize_v = ldv*k;
				bufsize_a = lda*m;
			}
			/* checking sizes with storev*/	
		}else if (storev == 'C')
		{	if ( side == 'L')
			{	
				ldv = k;
				lda =  n;
				bufsize_v = ldv*m;
				bufsize_a = lda*k;	
			} else if (side == 'R')
			{
				ldv = k;
				lda = k;
				bufsize_v = ldv*n;
				bufsize_a = lda*m;
			}
		}
		direct = 'B';
		ldb = n;
		bufsize_b = ldb*m;
	}else
	{
		EXPECT_TRUE(false) << "matrix_layout invalid";
	}
	
	bufsize_t = ldt*k;
	bufsize_t_tpqrt = n*n; 
	
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&v, &vref, bufsize_v);
	lapacke_gtest_alloc_double_buffer_pair(&t, &tref, bufsize_t);	
	lapacke_gtest_alloc_double_buffer_pair(&b, &bref, bufsize_b);
	lapacke_gtest_alloc_double_buffer_pair(&a, &aref, bufsize_a);
	lapacke_gtest_alloc_double_buffer_pair(&A_tpqrt, &A_tpqrtref, (n*n));
	lapacke_gtest_alloc_double_buffer_pair(&b_tpqrt, &b_tpqrtref, bufsize_b);
	lapacke_gtest_alloc_double_buffer_pair(&t_tpqrt, &t_tpqrtref, bufsize_t_tpqrt);
	
	if ((v==NULL) || (vref==NULL) ||
		(t==NULL) || (tref==NULL) ||
		(a==NULL) || (aref==NULL) ||
		(b==NULL) || (bref==NULL) ||
		(A_tpqrt ==NULL) || (A_tpqrtref==NULL) ||
		(b_tpqrt == NULL) || (b_tpqrtref == NULL) ||
		(t_tpqrt == NULL) || (t_tpqrtref == NULL))
	{
		EXPECT_FALSE( true) << "tprfb_double_parameters object: malloc error.";
		tprfb_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand( v, vref, bufsize_v);
	lapacke_gtest_init_double_buffer_pair_rand( b, bref, bufsize_b);
	lapacke_gtest_init_double_buffer_pair_rand( a, aref, bufsize_a);
	lapacke_gtest_init_double_buffer_pair_rand_custom_matrix( t, tref, ldt, k, 'S');
	lapacke_gtest_init_double_buffer_pair_rand( A_tpqrt, A_tpqrtref, (n*n));
	lapacke_gtest_init_double_buffer_pair_rand( b_tpqrt, b_tpqrtref, bufsize_b);


} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
tprfb_double_parameters :: ~tprfb_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" tprfb_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   tprfb_free();

}
/*  Test fixture class definition */
class dtprfb_test  : public  ::testing::Test {
public:
   tprfb_double_parameters  *dtprfb_obj;
   void SetUp();
   void TearDown () { delete dtprfb_obj;}
};

void dtprfb_test::SetUp(){
	/* LAPACKE dtpqrt prototype */
	typedef int (*Fptr_NL_LAPACKE_dtpqrt) (int matrix_layout, lapack_int m, lapack_int n,\
	lapack_int l, lapack_int nb, double* a, lapack_int lda, double* b, lapack_int ldb, double* t, lapack_int ldt);
	
	Fptr_NL_LAPACKE_dtpqrt dtpqrt;
	
	 /* LAPACKE dtprfb prototype */
    typedef int (*Fptr_NL_LAPACKE_dtprfb) (int matrix_layout, char side, char trans, char direct, char storev, \
	lapack_int m, lapack_int n, lapack_int k, lapack_int l, const double * v, lapack_int ldv, const double * t, \
	lapack_int ldt, double * a, lapack_int lda, double * b, lapack_int ldb);

    Fptr_NL_LAPACKE_dtprfb dtprfb;

    dtprfb_obj = new tprfb_double_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].side,
						   eig_paramslist[idx].trans,
						   eig_paramslist[idx].m,
						   eig_paramslist[idx].n,
                           eig_paramslist[idx].m,
						   eig_paramslist[idx].k,
						   eig_paramslist[idx].storev);
						   

    idx = Circular_Increment_Index(idx);

    dtprfb_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dtprfb_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dtprfb_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dtprfb_obj->hModule != NULL) << "Netlib lapacke handle NULL";
	
	/*dtprfb library call */
    dtprfb = (Fptr_NL_LAPACKE_dtprfb)dlsym(dtprfb_obj->hModule, "LAPACKE_dtprfb");
    ASSERT_TRUE(dtprfb != NULL) << "failed to get the Netlib LAPACKE_dtprfb symbol";
	
	/*dtpqrt library call*/
	dtprfb_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dtprfb_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dtprfb_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dtprfb_obj->lModule != NULL) << "Netlib lapacke handle NULL";
	
	dtpqrt = (Fptr_NL_LAPACKE_dtpqrt)dlsym(dtprfb_obj->hModule, "LAPACKE_dtpqrt");
    ASSERT_TRUE(dtpqrt != NULL) << "failed to get the Netlib LAPACKE_dtpqrt symbol";
    

    dtprfb_obj->inforef_tpqrt = dtpqrt( dtprfb_obj->matrix_layout, dtprfb_obj->m,
								dtprfb_obj->n, min(dtprfb_obj->m, dtprfb_obj->n), dtprfb_obj->n, dtprfb_obj->A_tpqrtref,
								dtprfb_obj->n, dtprfb_obj->b_tpqrtref, dtprfb_obj->ldb, dtprfb_obj->t_tpqrtref, dtprfb_obj->n);

    /* Compute libflame's Lapacke o/p  */
    dtprfb_obj->info_tpqrt = LAPACKE_dtpqrt( dtprfb_obj->matrix_layout, dtprfb_obj->m, dtprfb_obj->n,
											min(dtprfb_obj->m, dtprfb_obj->n), dtprfb_obj->n, dtprfb_obj->A_tpqrt, 
											dtprfb_obj->n, dtprfb_obj->b_tpqrt,	dtprfb_obj->ldb, dtprfb_obj->t_tpqrt, dtprfb_obj->n);

    if( dtprfb_obj->info_tpqrt < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dtpqrt is wrong\n", dtprfb_obj->info );
    }
    if( dtprfb_obj->inforef_tpqrt < 0 ) {
        printf( "The i:%d th argument with Netlib dtpqrt is wrong\n", 
        dtprfb_obj->inforef );
    }
	/*copy */
	
  
/*Compute dtprfb's  o/p */
    dtprfb_obj->inforef = dtprfb( dtprfb_obj->matrix_layout, dtprfb_obj->side, dtprfb_obj->trans, dtprfb_obj->direct, dtprfb_obj->storev,\
	dtprfb_obj->m, dtprfb_obj->n,dtprfb_obj->k, dtprfb_obj->l, (const double*)dtprfb_obj->b_tpqrtref,dtprfb_obj->ldv, (const double*)dtprfb_obj->t_tpqrtref, dtprfb_obj->ldt,
	dtprfb_obj->aref, dtprfb_obj->lda, dtprfb_obj->b, dtprfb_obj->ldb);

    /* Compute libflame's Lapacke o/p  */
    dtprfb_obj->info = LAPACKE_dtprfb( dtprfb_obj->matrix_layout, dtprfb_obj->side, dtprfb_obj->trans, dtprfb_obj->direct, dtprfb_obj->storev,\
	dtprfb_obj->m, dtprfb_obj->n, dtprfb_obj->k, dtprfb_obj->l, (const double*)dtprfb_obj->b_tpqrt, dtprfb_obj->ldv, (const double*)dtprfb_obj->t_tpqrt, dtprfb_obj->ldt, \
	dtprfb_obj->a, dtprfb_obj->lda, dtprfb_obj->b, dtprfb_obj->ldb);
	
    if( dtprfb_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame\
        LAPACKE_dtprfb is wrong\n", dtprfb_obj->info );
    }
    if( dtprfb_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dtprfb is wrong\n", 
        dtprfb_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dtprfb_obj->diff =  computeDiff_d( dtprfb_obj->bufsize_a, 
                dtprfb_obj->a, dtprfb_obj->aref );
				
	dtprfb_obj->diff_b =  computeDiff_d( dtprfb_obj->bufsize_b,
						dtprfb_obj->b, dtprfb_obj->bref );
	

}

TEST_F(dtprfb_test, dtprfb1) {
    EXPECT_NEAR(0.0, dtprfb_obj->diff, LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, dtprfb_obj->diff_b, LAPACKE_EIG_THRESHOLD);
}

TEST_F(dtprfb_test, dtprfb2) {
    EXPECT_NEAR(0.0, dtprfb_obj->diff, LAPACKE_EIG_THRESHOLD);
//	EXPECT_NEAR(0.0, dtprfb_obj->diff_b, LAPACKE_EIG_THRESHOLD);
}

TEST_F(dtprfb_test, dtprfb3) {
    EXPECT_NEAR(0.0, dtprfb_obj->diff, LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, dtprfb_obj->diff_b, LAPACKE_EIG_THRESHOLD);
}

TEST_F(dtprfb_test, dtprfb4) {
    EXPECT_NEAR(0.0, dtprfb_obj->diff, LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, dtprfb_obj->diff_b, LAPACKE_EIG_THRESHOLD);
}


/* Begin scomplex_common_parameters  class definition */
class tprfb_scomplex_parameters{

   public:
	int bufsize_t;
	int bufsize_v;
	int bufsize_a;
	int bufsize_b;
	int bufsize_t_tpqrt;
	void *hModule, *dModule;
	void *bModule, *lModule;
	float diff;
	float diff_b;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int k;
	lapack_int l;
	lapack_int lda, ldv, ldt, ldb;
	lapack_int ldt_tpqrt;
	lapack_complex_float* v, *A_tpqrt, *b_tpqrt, *t_tpqrt;
	char side, trans, direct, storev;
	/*Output Parameter*/
	lapack_complex_float* t, *a, *b;
	lapack_complex_float *A_tpqrtref, *b_tpqrtref, *t_tpqrtref;
	lapack_complex_float *vref, *tref, *aref, *bref;
	/*Return Values*/	
	lapack_int info, inforef;
	lapack_int info_tpqrt, inforef_tpqrt;

   public:
      tprfb_scomplex_parameters (int matrix_layout, char side, char trans,  lapack_int m, lapack_int n, lapack_int k, lapack_int l, char storev);
      ~tprfb_scomplex_parameters ();

};

/* Constructor definition  scomplex_common_parameters */
tprfb_scomplex_parameters:: tprfb_scomplex_parameters (int matrix_layout_i, char side_i, char trans_i , lapack_int m_i, lapack_int n_i, lapack_int k_i, lapack_int l_i, char storev_i)
{
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	k = k_i;
	l = l_i;
	ldt = k;
	side = side_i;
	trans = trans_i;
	storev = storev_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n tprfb scomplex: matrix_layout = %d,  m:%d, n: %d, trans:%c, l:%d \n", matrix_layout, side, m, n, trans, l);
	#endif

	if (matrix_layout == LAPACK_COL_MAJOR)
	{	/* checking sizes with storev*/	
		if (storev == 'C')
		{
			if ( side == 'L')
			{	
				ldv = m;
				bufsize_v = ldv*k;
				lda = k;
				bufsize_a = lda*n;
			} else if (side == 'R')
			{
				ldv = n;
				lda = m;
				bufsize_v = ldv*k;
				bufsize_a = lda*k;
			}
			/* checking sizes with storev*/	
		} else if (storev == 'R')
		{
			if ( side == 'L')
			{	
				ldv = k;
				bufsize_v = ldv*m;
				lda = k;
				bufsize_a = lda*n;
			} else if (side == 'R')
			{
				ldv = k;
				bufsize_v = ldv*n;
				lda = m;
				bufsize_a = lda*k;
			}
		}
		direct = 'F';
		ldb = m;
		bufsize_b = ldb*n;
		
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{		
		if (storev == 'R')
		{
			if ( side == 'L')
			{	
				ldv = m;
				lda =  n;
				bufsize_v = ldv*k;
				bufsize_a = lda*k;
			} else if (side == 'R')
			{
				ldv = n;
				lda = k;
				bufsize_v = ldv*k;
				bufsize_a = lda*m;
			}
			/* checking sizes with storev*/	
		}else if (storev == 'C')
		{	if ( side == 'L')
			{	
				ldv = k;
				lda =  n;
				bufsize_v = ldv*m;
				bufsize_a = lda*k;	
			} else if (side == 'R')
			{
				ldv = k;
				lda = k;
				bufsize_v = ldv*n;
				bufsize_a = lda*m;
			}
		}
		direct = 'B';
		ldb = n;
		bufsize_b = ldb*m;
	}else
	{
		EXPECT_TRUE(false) << "matrix_layout invalid";
	}
	
	bufsize_t = ldt*k;
	bufsize_t_tpqrt = n*n; 
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&v, &vref, bufsize_v);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&t, &tref, bufsize_t);	
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&b, &bref, bufsize_b);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&a, &aref, bufsize_a);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A_tpqrt, &A_tpqrtref, (n*n));
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&b_tpqrt, &b_tpqrtref, bufsize_b);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&t_tpqrt, &t_tpqrtref, bufsize_t_tpqrt);
	
	if ((v==NULL) || (vref==NULL) ||
		(t==NULL) || (tref==NULL) ||
		(a==NULL) || (aref==NULL) ||
		(b==NULL) || (bref==NULL) ||
		(A_tpqrt ==NULL) || (A_tpqrtref==NULL) ||
		(b_tpqrt == NULL) || (b_tpqrtref == NULL) ||
		(t_tpqrt == NULL) || (t_tpqrtref == NULL))
	{
		EXPECT_FALSE( true) << "tprfb_scomplex_parameters object: malloc error.";
		tprfb_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( v, vref, bufsize_v);
	lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, bufsize_b);
	lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, bufsize_a);
	lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix( t, tref, ldt, k, 'S');
	lapacke_gtest_init_scomplex_buffer_pair_rand( A_tpqrt, A_tpqrtref, (n*n));
	lapacke_gtest_init_scomplex_buffer_pair_rand( b_tpqrt, b_tpqrtref, bufsize_b);


} /* end of Constructor  */

/* Destructor definition  'scomplex_common_parameters' class  */
tprfb_scomplex_parameters :: ~tprfb_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" tprfb_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   tprfb_free();

}
/*  Test fixture class definition */
class ctprfb_test  : public  ::testing::Test {
public:
   tprfb_scomplex_parameters  *ctprfb_obj;
   void SetUp();
   void TearDown () { delete ctprfb_obj;}
};

void ctprfb_test::SetUp(){
	/* LAPACKE stpqrt prototype */
	typedef int (*Fptr_NL_LAPACKE_ctpqrt) (int matrix_layout, lapack_int m, lapack_int n,\
	lapack_int l, lapack_int nb, lapack_complex_float *a, lapack_int lda, lapack_complex_float *b, lapack_int ldb, lapack_complex_float *t, lapack_int ldt);
	
	Fptr_NL_LAPACKE_ctpqrt ctpqrt;
	
	 /* LAPACKE ctprfb prototype */
    typedef int (*Fptr_NL_LAPACKE_ctprfb) (int matrix_layout, char side, char trans, char direct, char storev, \
	lapack_int m, lapack_int n, lapack_int k, lapack_int l, const lapack_complex_float *v, lapack_int ldv, const lapack_complex_float * t, \
	lapack_int ldt, lapack_complex_float *a, lapack_int lda, lapack_complex_float * b, lapack_int ldb);

    Fptr_NL_LAPACKE_ctprfb ctprfb;

    ctprfb_obj = new tprfb_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].side,
						   eig_paramslist[idx].trans,
						   eig_paramslist[idx].m,
						   eig_paramslist[idx].n,
                           eig_paramslist[idx].m,
						   eig_paramslist[idx].k,
						   eig_paramslist[idx].storev);
						   

    idx = Circular_Increment_Index(idx);

    ctprfb_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ctprfb_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ctprfb_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ctprfb_obj->hModule != NULL) << "Netlib lapacke handle NULL";
	
	/*ctprfb library call */
    ctprfb = (Fptr_NL_LAPACKE_ctprfb)dlsym(ctprfb_obj->hModule, "LAPACKE_ctprfb");
    ASSERT_TRUE(ctprfb != NULL) << "failed to get the Netlib LAPACKE_ctprfb symbol";
	
	/*stpqrt library call*/
	ctprfb_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ctprfb_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ctprfb_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ctprfb_obj->lModule != NULL) << "Netlib lapacke handle NULL";
	
	ctpqrt = (Fptr_NL_LAPACKE_ctpqrt)dlsym(ctprfb_obj->hModule, "LAPACKE_ctpqrt");
    ASSERT_TRUE(ctpqrt != NULL) << "failed to get the Netlib LAPACKE_ctpqrt symbol";
    

    ctprfb_obj->inforef_tpqrt = ctpqrt( ctprfb_obj->matrix_layout, ctprfb_obj->m,
								ctprfb_obj->n, min(ctprfb_obj->m, ctprfb_obj->n), ctprfb_obj->n, ctprfb_obj->A_tpqrtref,
								ctprfb_obj->n, ctprfb_obj->b_tpqrtref, ctprfb_obj->ldb, ctprfb_obj->t_tpqrtref, ctprfb_obj->n);

    /* Compute libflame's Lapacke o/p  */
    ctprfb_obj->info_tpqrt = LAPACKE_ctpqrt( ctprfb_obj->matrix_layout, ctprfb_obj->m, ctprfb_obj->n,
											min(ctprfb_obj->m, ctprfb_obj->n), ctprfb_obj->n, ctprfb_obj->A_tpqrt, 
											ctprfb_obj->n, ctprfb_obj->b_tpqrt,	ctprfb_obj->ldb, ctprfb_obj->t_tpqrt, ctprfb_obj->n);

    if( ctprfb_obj->info_tpqrt < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_stpqrt is wrong\n", ctprfb_obj->info );
    }
    if( ctprfb_obj->inforef_tpqrt < 0 ) {
        printf( "The i:%d th argument with Netlib stpqrt is wrong\n", 
        ctprfb_obj->inforef );
    }
	/*copy */
	
  
/*Compute ctprfb's  o/p */
    ctprfb_obj->inforef = ctprfb( ctprfb_obj->matrix_layout, ctprfb_obj->side, ctprfb_obj->trans, ctprfb_obj->direct, ctprfb_obj->storev,\
	ctprfb_obj->m, ctprfb_obj->n,ctprfb_obj->k, ctprfb_obj->l, (const lapack_complex_float*)ctprfb_obj->b_tpqrtref,ctprfb_obj->ldv, (const lapack_complex_float*)ctprfb_obj->t_tpqrtref, ctprfb_obj->ldt,
	ctprfb_obj->aref, ctprfb_obj->lda, ctprfb_obj->b, ctprfb_obj->ldb);

    /* Compute libflame's Lapacke o/p  */
    ctprfb_obj->info = LAPACKE_ctprfb( ctprfb_obj->matrix_layout, ctprfb_obj->side, ctprfb_obj->trans, ctprfb_obj->direct, ctprfb_obj->storev,\
	ctprfb_obj->m, ctprfb_obj->n, ctprfb_obj->k, ctprfb_obj->l, (const lapack_complex_float*)ctprfb_obj->b_tpqrt, ctprfb_obj->ldv, (const lapack_complex_float*)ctprfb_obj->t_tpqrt, ctprfb_obj->ldt, \
	ctprfb_obj->a, ctprfb_obj->lda, ctprfb_obj->b, ctprfb_obj->ldb);
	
    if( ctprfb_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame\
        LAPACKE_ctprfb is wrong\n", ctprfb_obj->info );
    }
    if( ctprfb_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ctprfb is wrong\n", 
        ctprfb_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    ctprfb_obj->diff =  computeDiff_c( ctprfb_obj->bufsize_a, 
                ctprfb_obj->a, ctprfb_obj->aref );
				
	ctprfb_obj->diff_b =  computeDiff_c( ctprfb_obj->bufsize_b,
						ctprfb_obj->b, ctprfb_obj->bref );
	

}

TEST_F(ctprfb_test, ctprfb1) {
    EXPECT_NEAR(0.0, ctprfb_obj->diff, LAPACKE_EIG_THRESHOLD);
//	EXPECT_NEAR(0.0, ctprfb_obj->diff_b, LAPACKE_EIG_THRESHOLD);
}

TEST_F(ctprfb_test, ctprfb2) {
    EXPECT_NEAR(0.0, ctprfb_obj->diff, LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, ctprfb_obj->diff_b, LAPACKE_EIG_THRESHOLD);
}

TEST_F(ctprfb_test, ctprfb3) {
    EXPECT_NEAR(0.0, ctprfb_obj->diff, LAPACKE_EIG_THRESHOLD);
//	EXPECT_NEAR(0.0, ctprfb_obj->diff_b, LAPACKE_EIG_THRESHOLD);
}

TEST_F(ctprfb_test, ctprfb4) {
    EXPECT_NEAR(0.0, ctprfb_obj->diff, LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, ctprfb_obj->diff_b, LAPACKE_EIG_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class tprfb_dcomplex_parameters{

   public:
	int bufsize_t;
	int bufsize_v;
	int bufsize_a;
	int bufsize_b;
	int bufsize_t_tpqrt;
	void *hModule, *dModule;
	void *bModule, *lModule;
	double diff;
	double diff_b;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int k;
	lapack_int l;
	lapack_int lda, ldv, ldt, ldb;
	lapack_int ldt_tpqrt;
	lapack_complex_double* v, *A_tpqrt, *b_tpqrt, *t_tpqrt;
	char side, trans, direct, storev;
	/*Output Parameter*/
	lapack_complex_double* t, *a, *b;
	lapack_complex_double *A_tpqrtref, *b_tpqrtref, *t_tpqrtref;
	lapack_complex_double *vref, *tref, *aref, *bref;
	/*Return Values*/	
	lapack_int info, inforef;
	lapack_int info_tpqrt, inforef_tpqrt;

   public:
      tprfb_dcomplex_parameters (int matrix_layout, char side, char trans,  lapack_int m, lapack_int n, lapack_int k, lapack_int l, char storev);
      ~tprfb_dcomplex_parameters ();

};

/* Constructor definition  dcomplex_common_parameters */
tprfb_dcomplex_parameters:: tprfb_dcomplex_parameters (int matrix_layout_i, char side_i, char trans_i , lapack_int m_i, lapack_int n_i, lapack_int k_i, lapack_int l_i, char storev_i)
{
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	k = k_i;
	l = l_i;
	ldt = k;
	side = side_i;
	trans = trans_i;
	storev = storev_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n tprfb dcomplex: matrix_layout = %d,  m:%d, n: %d, trans:%c, l:%d \n", matrix_layout, side, m, n, trans, l);
	#endif

	if (matrix_layout == LAPACK_COL_MAJOR)
	{	/* checking sizes with storev*/	
		if (storev == 'C')
		{
			if ( side == 'L')
			{	
				ldv = m;
				bufsize_v = ldv*k;
				lda = k;
				bufsize_a = lda*n;
			} else if (side == 'R')
			{
				ldv = n;
				lda = m;
				bufsize_v = ldv*k;
				bufsize_a = lda*k;
			}
			/* checking sizes with storev*/	
		} else if (storev == 'R')
		{
			if ( side == 'L')
			{	
				ldv = k;
				bufsize_v = ldv*m;
				lda = k;
				bufsize_a = lda*n;
			} else if (side == 'R')
			{
				ldv = k;
				bufsize_v = ldv*n;
				lda = m;
				bufsize_a = lda*k;
			}
		}
		direct = 'F';
		ldb = m;
		bufsize_b = ldb*n;
		
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{		
		if (storev == 'R')
		{
			if ( side == 'L')
			{	
				ldv = m;
				lda =  n;
				bufsize_v = ldv*k;
				bufsize_a = lda*k;
			} else if (side == 'R')
			{
				ldv = n;
				lda = k;
				bufsize_v = ldv*k;
				bufsize_a = lda*m;
			}
			/* checking sizes with storev*/	
		}else if (storev == 'C')
		{	if ( side == 'L')
			{	
				ldv = k;
				lda =  n;
				bufsize_v = ldv*m;
				bufsize_a = lda*k;	
			} else if (side == 'R')
			{
				ldv = k;
				lda = k;
				bufsize_v = ldv*n;
				bufsize_a = lda*m;
			}
		}
		direct = 'B';
		ldb = n;
		bufsize_b = ldb*m;
	}else
	{
		EXPECT_TRUE(false) << "matrix_layout invalid";
	}
	
	bufsize_t = ldt*k;
	bufsize_t_tpqrt = n*n; 
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&v, &vref, bufsize_v);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&t, &tref, bufsize_t);	
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&b, &bref, bufsize_b);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&a, &aref, bufsize_a);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A_tpqrt, &A_tpqrtref, (n*n));
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&b_tpqrt, &b_tpqrtref, bufsize_b);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&t_tpqrt, &t_tpqrtref, bufsize_t_tpqrt);
	
	if ((v==NULL) || (vref==NULL) ||
		(t==NULL) || (tref==NULL) ||
		(a==NULL) || (aref==NULL) ||
		(b==NULL) || (bref==NULL) ||
		(A_tpqrt ==NULL) || (A_tpqrtref==NULL) ||
		(b_tpqrt == NULL) || (b_tpqrtref == NULL) ||
		(t_tpqrt == NULL) || (t_tpqrtref == NULL))
	{
		EXPECT_FALSE( true) << "tprfb_dcomplex_parameters object: malloc error.";
		tprfb_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( v, vref, bufsize_v);
	lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, bufsize_b);
	lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, bufsize_a);
	lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix( t, tref, ldt, k, 'S');
	lapacke_gtest_init_dcomplex_buffer_pair_rand( A_tpqrt, A_tpqrtref, (n*n));
	lapacke_gtest_init_dcomplex_buffer_pair_rand( b_tpqrt, b_tpqrtref, bufsize_b);


} /* end of Constructor  */

/* Destructor definition  'dcomplex_common_parameters' class  */
tprfb_dcomplex_parameters :: ~tprfb_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" tprfb_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   tprfb_free();

}
/*  Test fixture class definition */
class ztprfb_test  : public  ::testing::Test {
public:
   tprfb_dcomplex_parameters  *ztprfb_obj;
   void SetUp();
   void TearDown () { delete ztprfb_obj;}
};

void ztprfb_test::SetUp(){
	/* LAPACKE stpqrt prototype */
	typedef int (*Fptr_NL_LAPACKE_ztpqrt) (int matrix_layout, lapack_int m, lapack_int n,\
	lapack_int l, lapack_int nb, lapack_complex_double *a, lapack_int lda, lapack_complex_double *b, lapack_int ldb, lapack_complex_double *t, lapack_int ldt);
	
	Fptr_NL_LAPACKE_ztpqrt ztpqrt;
	
	 /* LAPACKE ztprfb prototype */
    typedef int (*Fptr_NL_LAPACKE_ztprfb) (int matrix_layout, char side, char trans, char direct, char storev, \
	lapack_int m, lapack_int n, lapack_int k, lapack_int l, const lapack_complex_double *v, lapack_int ldv, const lapack_complex_double * t, \
	lapack_int ldt, lapack_complex_double *a, lapack_int lda, lapack_complex_double * b, lapack_int ldb);

    Fptr_NL_LAPACKE_ztprfb ztprfb;

    ztprfb_obj = new tprfb_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].side,
						   eig_paramslist[idx].trans,
						   eig_paramslist[idx].m,
						   eig_paramslist[idx].n,
                           eig_paramslist[idx].m,
						   eig_paramslist[idx].k,
						   eig_paramslist[idx].storev);
						   

    idx = Circular_Increment_Index(idx);

    ztprfb_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ztprfb_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ztprfb_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ztprfb_obj->hModule != NULL) << "Netlib lapacke handle NULL";
	
	/*ztprfb library call */
    ztprfb = (Fptr_NL_LAPACKE_ztprfb)dlsym(ztprfb_obj->hModule, "LAPACKE_ztprfb");
    ASSERT_TRUE(ztprfb != NULL) << "failed to get the Netlib LAPACKE_ztprfb symbol";
	
	/*stpqrt library call*/
	ztprfb_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ztprfb_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ztprfb_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ztprfb_obj->lModule != NULL) << "Netlib lapacke handle NULL";
	
	ztpqrt = (Fptr_NL_LAPACKE_ztpqrt)dlsym(ztprfb_obj->hModule, "LAPACKE_ztpqrt");
    ASSERT_TRUE(ztpqrt != NULL) << "failed to get the Netlib LAPACKE_ztpqrt symbol";
    

    ztprfb_obj->inforef_tpqrt = ztpqrt( ztprfb_obj->matrix_layout, ztprfb_obj->m,
								ztprfb_obj->n, min(ztprfb_obj->m, ztprfb_obj->n), ztprfb_obj->n, ztprfb_obj->A_tpqrtref,
								ztprfb_obj->n, ztprfb_obj->b_tpqrtref, ztprfb_obj->ldb, ztprfb_obj->t_tpqrtref, ztprfb_obj->n);

    /* Compute libflame's Lapacke o/p  */
    ztprfb_obj->info_tpqrt = LAPACKE_ztpqrt( ztprfb_obj->matrix_layout, ztprfb_obj->m, ztprfb_obj->n,
											min(ztprfb_obj->m, ztprfb_obj->n), ztprfb_obj->n, ztprfb_obj->A_tpqrt, 
											ztprfb_obj->n, ztprfb_obj->b_tpqrt,	ztprfb_obj->ldb, ztprfb_obj->t_tpqrt, ztprfb_obj->n);

    if( ztprfb_obj->info_tpqrt < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_stpqrt is wrong\n", ztprfb_obj->info );
    }
    if( ztprfb_obj->inforef_tpqrt < 0 ) {
        printf( "The i:%d th argument with Netlib stpqrt is wrong\n", 
        ztprfb_obj->inforef );
    }
	/*copy */
	
  
/*Compute ztprfb's  o/p */
    ztprfb_obj->inforef = ztprfb( ztprfb_obj->matrix_layout, ztprfb_obj->side, ztprfb_obj->trans, ztprfb_obj->direct, ztprfb_obj->storev,\
	ztprfb_obj->m, ztprfb_obj->n,ztprfb_obj->k, ztprfb_obj->l, (const lapack_complex_double*)ztprfb_obj->b_tpqrtref,ztprfb_obj->ldv, (const lapack_complex_double*)ztprfb_obj->t_tpqrtref, ztprfb_obj->ldt,
	ztprfb_obj->aref, ztprfb_obj->lda, ztprfb_obj->b, ztprfb_obj->ldb);

    /* Compute libflame's Lapacke o/p  */
    ztprfb_obj->info = LAPACKE_ztprfb( ztprfb_obj->matrix_layout, ztprfb_obj->side, ztprfb_obj->trans, ztprfb_obj->direct, ztprfb_obj->storev,\
	ztprfb_obj->m, ztprfb_obj->n, ztprfb_obj->k, ztprfb_obj->l, (const lapack_complex_double*)ztprfb_obj->b_tpqrt, ztprfb_obj->ldv, (const lapack_complex_double*)ztprfb_obj->t_tpqrt, ztprfb_obj->ldt, \
	ztprfb_obj->a, ztprfb_obj->lda, ztprfb_obj->b, ztprfb_obj->ldb);
	
    if( ztprfb_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame\
        LAPACKE_ztprfb is wrong\n", ztprfb_obj->info );
    }
    if( ztprfb_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ztprfb is wrong\n", 
        ztprfb_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    ztprfb_obj->diff =  computeDiff_z( ztprfb_obj->bufsize_a, 
                ztprfb_obj->a, ztprfb_obj->aref );
				
	ztprfb_obj->diff_b =  computeDiff_z( ztprfb_obj->bufsize_b,
						ztprfb_obj->b, ztprfb_obj->bref );
	

}

TEST_F(ztprfb_test, ztprfb1) {
    EXPECT_NEAR(0.0, ztprfb_obj->diff, LAPACKE_EIG_THRESHOLD);
//	EXPECT_NEAR(0.0, ztprfb_obj->diff_b, LAPACKE_EIG_THRESHOLD);
}

TEST_F(ztprfb_test, ztprfb2) {
    EXPECT_NEAR(0.0, ztprfb_obj->diff, LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, ztprfb_obj->diff_b, LAPACKE_EIG_THRESHOLD);
}

TEST_F(ztprfb_test, ztprfb3) {
    EXPECT_NEAR(0.0, ztprfb_obj->diff, LAPACKE_EIG_THRESHOLD);
//	EXPECT_NEAR(0.0, ztprfb_obj->diff_b, LAPACKE_EIG_THRESHOLD);
}

TEST_F(ztprfb_test, ztprfb4) {
    EXPECT_NEAR(0.0, ztprfb_obj->diff, LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, ztprfb_obj->diff_b, LAPACKE_EIG_THRESHOLD);
}