#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define gemqrt_free() \
if (v!=NULL)    free(v); \
if (vref!=NULL) free(vref);\
if (t!=NULL)  free(t);\
if (tref!=NULL) free(tref); \
if (c != NULL) free(c); \
if (cref != NULL) free(cref); \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule);\
if( bModule != NULL) dlclose(bModule); \
if(lModule != NULL) dlclose(lModule);\

	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin float_common_parameters  class definition */
class gemqrt_float_parameters{

   public:
	int bufsize_v;
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
	lapack_int nb;
	float* v, *A;
	char side;
	char trans;
	lapack_int ldc, ldv, lda_geqrt, ldt;
	/*Output Parameter*/
	float* t, *c;	
	float *vref, *tref, *cref, *Aref;
	/*Return Values*/	
	lapack_int info, inforef;
	lapack_int info_geqrt, inforef_geqrt;

   public:
      gemqrt_float_parameters (int matrix_layout, char side , char trans, lapack_int m, lapack_int n, lapack_int nb);
      ~gemqrt_float_parameters ();

};

/* Constructor definition  float_common_parameters */
gemqrt_float_parameters:: gemqrt_float_parameters (int matrix_layout_i, char side_i , char trans_i,lapack_int m_i, lapack_int n_i, lapack_int nb_i)
{
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	side = side_i;
	trans = trans_i;
	nb = min(m,n);

	if (trans == 'C')
		trans  = 'N' ;
	
	#if LAPACKE_TEST_VERBOSE
	printf(" \n gemqrt float: matrix_layout = %d, side :%c, trans:%c, m:%d, n: %d \n", matrix_layout, side, trans, m, n);
	#endif	
	
	/*Sizes based on side*/
	if (side  == 'L')
	{
		k = m;
		nb = k;
		if (matrix_layout == LAPACK_COL_MAJOR)
		{
			ldv = m;
			ldt = nb;
			ldc = m;
			lda_geqrt = m;
			bufsize_a = lda_geqrt*n;
			bufsize_v = ldv *k;
			bufsize_t = ldt*min(m,n);
			bufsize_c = ldc*n;
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{	
			ldv = k;
			ldt = k;
			ldc = n;
			lda_geqrt = n;
			bufsize_a = lda_geqrt*m;
			bufsize_v = ldv *m;
			bufsize_t = ldt*nb;
			bufsize_c = ldc*m;
		} else
		{
			EXPECT_TRUE(false) << "matrix_layout invalid";
		}
	}else if (side == 'R')
	{
		k = n;
		nb = k;
		if (matrix_layout == LAPACK_COL_MAJOR)
		{
			ldv = n;
			ldt = nb;
			ldc = m;
			lda_geqrt = m;
			bufsize_a = lda_geqrt*n;
			bufsize_v = ldv*n;
			bufsize_t = ldt*min(m,n);
			bufsize_c = ldc*n;
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{
			ldv = k;
			ldt = k;
			ldc = n;
			lda_geqrt = n;
			bufsize_a = lda_geqrt*m;
			bufsize_v = ldv*n;
			bufsize_t = ldt*nb;
			bufsize_c = ldc*m;
		} else
		{
			EXPECT_TRUE(false) << "matrix_layout invalid";
		}
	}else 
		EXPECT_TRUE(false) << "side is invalid";	
	
	bufsize_c = m*n;
	
	printf(" \n gemqrt float: matrix_layout = %d, side :%c, trans:%c, m:%d, n: %d , ldv:%d, ldc:%d, k:%d\n", matrix_layout, side, trans, m, n, ldv, ldc, k);
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_float_buffer_pair(&v, &vref, bufsize_v);
	lapacke_gtest_alloc_float_buffer_pair(&t, &tref, bufsize_t);	
	lapacke_gtest_alloc_float_buffer_pair(&c, &cref, bufsize_c);
	
	if ((v==NULL) || (vref==NULL) ||
		(t==NULL) || (tref==NULL) ||
		(c==NULL) || (cref==NULL) ||
		(A ==NULL) || (Aref ==NULL))
	{
		EXPECT_FALSE( true) << "gemqrt_float_parameters object: malloc error.";
		gemqrt_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand( v, vref, bufsize_v);
	lapacke_gtest_init_float_buffer_pair_rand( c, cref, bufsize_c);
	lapacke_gtest_init_float_buffer_pair_rand(A, Aref, bufsize_a);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
gemqrt_float_parameters :: ~gemqrt_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" gemqrt_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   gemqrt_free();

}
/*  Test fixture class definition */
class sgemqrt_test  : public  ::testing::Test {
public:
   gemqrt_float_parameters  *sgemqrt_obj;
   void SetUp();
   void TearDown () { delete sgemqrt_obj;}
};

void sgemqrt_test::SetUp(){

	/* LAPACKE sgeqrt prototype */
    typedef int (*Fptr_NL_LAPACKE_sgeqrt) (int matrix_layout, lapack_int m, lapack_int n, lapack_int nb,\
	float* a, lapack_int lda, float* t, lapack_int ldt);
	
	 Fptr_NL_LAPACKE_sgeqrt sgeqrt; 
	 /* LAPACKE sgemqrt prototype */
    typedef int (*Fptr_NL_LAPACKE_sgemqrt) (int matrix_layout, char side, char trans, lapack_int m, lapack_int n,\
	lapack_int k, lapack_int nb, const float* v, lapack_int ldv, const float* t,\
	lapack_int ldt, float* c, lapack_int ldc);

    Fptr_NL_LAPACKE_sgemqrt sgemqrt;

    sgemqrt_obj = new gemqrt_float_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].side,
						   eig_paramslist[idx].trans,
						   eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
						   eig_paramslist[idx].nb);
						   

    idx = Circular_Increment_Index(idx);

    sgemqrt_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgemqrt_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgemqrt_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgemqrt_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*sgemqrt library call */
    sgemqrt = (Fptr_NL_LAPACKE_sgemqrt)dlsym(sgemqrt_obj->hModule, "LAPACKE_sgemqrt");
    ASSERT_TRUE(sgemqrt != NULL) << "failed to get the Netlib LAPACKE_sgemqrt symbol";

	/*sgeqrt library call*/
	sgemqrt_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgemqrt_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgemqrt_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgemqrt_obj->lModule != NULL) << "Netlib lapacke handle NULL";

    sgeqrt = (Fptr_NL_LAPACKE_sgeqrt)dlsym(sgemqrt_obj->lModule, "LAPACKE_sgeqrt");
    ASSERT_TRUE(sgeqrt != NULL) << "failed to get the Netlib LAPACKE_sgeqrt symbol";    
    

    sgemqrt_obj->inforef_geqrt = sgeqrt( sgemqrt_obj->matrix_layout,sgemqrt_obj->m, sgemqrt_obj->n, min(sgemqrt_obj->m,sgemqrt_obj->n),
								sgemqrt_obj->Aref, sgemqrt_obj->lda_geqrt, sgemqrt_obj->tref, min(sgemqrt_obj->m,sgemqrt_obj->n));

    /* Compute libflame's Lapacke o/p  */
    sgemqrt_obj->info_geqrt = LAPACKE_sgeqrt( sgemqrt_obj->matrix_layout,sgemqrt_obj->m, sgemqrt_obj->n, min(sgemqrt_obj->m,sgemqrt_obj->n),
											sgemqrt_obj->A, sgemqrt_obj->lda_geqrt, sgemqrt_obj->t, min(sgemqrt_obj->m,sgemqrt_obj->n));
										
										

    if( sgemqrt_obj->info_geqrt < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_sgeqrt is wrong\n", sgemqrt_obj->info_geqrt );
    }
    if( sgemqrt_obj->inforef_geqrt < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sgeqrt is wrong\n",
        sgemqrt_obj->inforef_geqrt ); 
    }  
/*Compute sgemqrt's  o/p */
    sgemqrt_obj->inforef = sgemqrt( sgemqrt_obj->matrix_layout, sgemqrt_obj->side, sgemqrt_obj->trans,
								sgemqrt_obj->m, sgemqrt_obj->n, sgemqrt_obj->k, sgemqrt_obj->nb, (const float*)sgemqrt_obj->vref,
								sgemqrt_obj->ldv, (const float*)sgemqrt_obj->tref, sgemqrt_obj->ldt, sgemqrt_obj->cref, sgemqrt_obj->ldc);

    /* Compute libflame's Lapacke o/p  */
    sgemqrt_obj->info = LAPACKE_sgemqrt( sgemqrt_obj->matrix_layout, sgemqrt_obj->side, sgemqrt_obj->trans,
									sgemqrt_obj->m, sgemqrt_obj->n, sgemqrt_obj->k, sgemqrt_obj->nb, (const float*)sgemqrt_obj->v,
									sgemqrt_obj->ldv, (const float*)sgemqrt_obj->t, sgemqrt_obj->ldt, sgemqrt_obj->c, sgemqrt_obj->ldc);
    if( sgemqrt_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame\
        LAPACKE_sgemqrt is wrong\n", sgemqrt_obj->info );
    }
    if( sgemqrt_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sgemqrt is wrong\n", 
        sgemqrt_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps */
    sgemqrt_obj->diff =  computeDiff_s( sgemqrt_obj->bufsize_c, 
                sgemqrt_obj->c, sgemqrt_obj->cref );

}

TEST_F(sgemqrt_test, sgemqrt1) {
    EXPECT_NEAR(0.0, sgemqrt_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(sgemqrt_test, sgemqrt2) {
    EXPECT_NEAR(0.0, sgemqrt_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(sgemqrt_test, sgemqrt3) {
    EXPECT_NEAR(0.0, sgemqrt_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(sgemqrt_test, sgemqrt4) {
    EXPECT_NEAR(0.0, sgemqrt_obj->diff, LAPACKE_EIG_THRESHOLD);
}

/* Begin double_common_parameters  class definition */
class gemqrt_double_parameters{

   public:
	int bufsize_v;
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
	lapack_int nb;
	double* v, *A;
	char side;
	char trans;
	lapack_int ldc, ldv, lda_geqrt, ldt;
	/*Output Parameter*/
	double* t, *c;	
	double *vref, *tref, *cref, *Aref;
	/*Return Values*/	
	lapack_int info, inforef;
	lapack_int info_geqrt, inforef_geqrt;

   public:
      gemqrt_double_parameters (int matrix_layout, char side , char trans, lapack_int m, lapack_int n, lapack_int nb);
      ~gemqrt_double_parameters ();

};

/* Constructor definition  double_common_parameters */
gemqrt_double_parameters:: gemqrt_double_parameters (int matrix_layout_i, char side_i , char trans_i,lapack_int m_i, lapack_int n_i, lapack_int nb_i)
{
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	side = side_i;
	trans = trans_i;
	nb = min(m,n);

	if (trans == 'C')
		trans  = 'N' ;
	
	#if LAPACKE_TEST_VERBOSE
	printf(" \n gemqrt double: matrix_layout = %d, side :%c, trans:%c, m:%d, n: %d \n", matrix_layout, side, trans, m, n);
	#endif	
	
	/*Sizes based on side*/
	if (side  == 'L')
	{
		k = m;
		nb = k;
		if (matrix_layout == LAPACK_COL_MAJOR)
		{
			ldv = m;
			ldt = nb;
			ldc = m;
			lda_geqrt = m;
			bufsize_a = lda_geqrt*n;
			bufsize_v = ldv *k;
			bufsize_t = ldt*min(m,n);
			bufsize_c = ldc*n;
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{	
			ldv = k;
			ldt = k;
			ldc = n;
			lda_geqrt = n;
			bufsize_a = lda_geqrt*m;
			bufsize_v = ldv *m;
			bufsize_t = ldt*nb;
			bufsize_c = ldc*m;
		} else
		{
			EXPECT_TRUE(false) << "matrix_layout invalid";
		}
	}else if (side == 'R')
	{
		k = n;
		nb = k;
		if (matrix_layout == LAPACK_COL_MAJOR)
		{
			ldv = n;
			ldt = nb;
			ldc = m;
			lda_geqrt = m;
			bufsize_a = lda_geqrt*n;
			bufsize_v = ldv*n;
			bufsize_t = ldt*min(m,n);
			bufsize_c = ldc*n;
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{
			ldv = k;
			ldt = k;
			ldc = n;
			lda_geqrt = n;
			bufsize_a = lda_geqrt*m;
			bufsize_v = ldv*n;
			bufsize_t = ldt*nb;
			bufsize_c = ldc*m;
		} else
		{
			EXPECT_TRUE(false) << "matrix_layout invalid";
		}
	}else 
		EXPECT_TRUE(false) << "side is invalid";	
	
	bufsize_c = m*n;
	
	
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_double_buffer_pair(&v, &vref, bufsize_v);
	lapacke_gtest_alloc_double_buffer_pair(&t, &tref, bufsize_t);	
	lapacke_gtest_alloc_double_buffer_pair(&c, &cref, bufsize_c);
	
	if ((v==NULL) || (vref==NULL) ||
		(t==NULL) || (tref==NULL) ||
		(c==NULL) || (cref==NULL) ||
		(A ==NULL) || (Aref ==NULL))
	{
		EXPECT_FALSE( true) << "gemqrt_double_parameters object: malloc error.";
		gemqrt_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand( v, vref, bufsize_v);
	lapacke_gtest_init_double_buffer_pair_rand( c, cref, bufsize_c);
	lapacke_gtest_init_double_buffer_pair_rand(A, Aref, bufsize_a);

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
gemqrt_double_parameters :: ~gemqrt_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" gemqrt_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   gemqrt_free();

}
/*  Test fixture class definition */
class dgemqrt_test  : public  ::testing::Test {
public:
   gemqrt_double_parameters  *dgemqrt_obj;
   void SetUp();
   void TearDown () { delete dgemqrt_obj;}
};

void dgemqrt_test::SetUp(){

	/* LAPACKE dgeqrt prototype */
    typedef int (*Fptr_NL_LAPACKE_dgeqrt) (int matrix_layout, lapack_int m, lapack_int n, lapack_int nb,\
	double* a, lapack_int lda, double* t, lapack_int ldt);
	
	 Fptr_NL_LAPACKE_dgeqrt dgeqrt; 
	 /* LAPACKE dgemqrt prototype */
    typedef int (*Fptr_NL_LAPACKE_dgemqrt) (int matrix_layout, char side, char trans, lapack_int m, lapack_int n,\
	lapack_int k, lapack_int nb, const double* v, lapack_int ldv, const double* t,\
	lapack_int ldt, double* c, lapack_int ldc);

    Fptr_NL_LAPACKE_dgemqrt dgemqrt;

    dgemqrt_obj = new gemqrt_double_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].side,
						   eig_paramslist[idx].trans,
						   eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
						   eig_paramslist[idx].nb);
						   

    idx = Circular_Increment_Index(idx);

    dgemqrt_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgemqrt_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgemqrt_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgemqrt_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*dgemqrt library call */
    dgemqrt = (Fptr_NL_LAPACKE_dgemqrt)dlsym(dgemqrt_obj->hModule, "LAPACKE_dgemqrt");
    ASSERT_TRUE(dgemqrt != NULL) << "failed to get the Netlib LAPACKE_dgemqrt symbol";

	/*dgeqrt library call*/
	dgemqrt_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgemqrt_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgemqrt_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgemqrt_obj->lModule != NULL) << "Netlib lapacke handle NULL";

    dgeqrt = (Fptr_NL_LAPACKE_dgeqrt)dlsym(dgemqrt_obj->lModule, "LAPACKE_dgeqrt");
    ASSERT_TRUE(dgeqrt != NULL) << "failed to get the Netlib LAPACKE_dgeqrt symbol";    
    

    dgemqrt_obj->inforef_geqrt = dgeqrt( dgemqrt_obj->matrix_layout,dgemqrt_obj->m, dgemqrt_obj->n, min(dgemqrt_obj->m,dgemqrt_obj->n),
								dgemqrt_obj->Aref, dgemqrt_obj->lda_geqrt, dgemqrt_obj->tref, min(dgemqrt_obj->m,dgemqrt_obj->n));

    /* Compute libflame's Lapacke o/p  */
    dgemqrt_obj->info_geqrt = LAPACKE_dgeqrt( dgemqrt_obj->matrix_layout,dgemqrt_obj->m, dgemqrt_obj->n, min(dgemqrt_obj->m,dgemqrt_obj->n),
											dgemqrt_obj->A, dgemqrt_obj->lda_geqrt, dgemqrt_obj->t, min(dgemqrt_obj->m,dgemqrt_obj->n));
										
										

    if( dgemqrt_obj->info_geqrt < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dgeqrt is wrong\n", dgemqrt_obj->info_geqrt );
    }
    if( dgemqrt_obj->inforef_geqrt < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dgeqrt is wrong\n", 
        dgemqrt_obj->inforef_geqrt ); 
    }  
/*Compute dgemqrt's  o/p */
    dgemqrt_obj->inforef = dgemqrt( dgemqrt_obj->matrix_layout, dgemqrt_obj->side, dgemqrt_obj->trans,
								dgemqrt_obj->m, dgemqrt_obj->n, dgemqrt_obj->k, dgemqrt_obj->nb, (const double*)dgemqrt_obj->vref,
								dgemqrt_obj->ldv, (const double*)dgemqrt_obj->tref, dgemqrt_obj->ldt, dgemqrt_obj->cref, dgemqrt_obj->ldc);

    /* Compute libflame's Lapacke o/p  */
    dgemqrt_obj->info = LAPACKE_dgemqrt( dgemqrt_obj->matrix_layout, dgemqrt_obj->side, dgemqrt_obj->trans,
									dgemqrt_obj->m, dgemqrt_obj->n, dgemqrt_obj->k, dgemqrt_obj->nb, (const double*)dgemqrt_obj->v,
									dgemqrt_obj->ldv, (const double*)dgemqrt_obj->t, dgemqrt_obj->ldt, dgemqrt_obj->c, dgemqrt_obj->ldc);
    if( dgemqrt_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame\
        LAPACKE_dgemqrt is wrong\n", dgemqrt_obj->info );
    }
    if( dgemqrt_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dgemqrt is wrong\n", 
        dgemqrt_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps */
    dgemqrt_obj->diff =  computeDiff_d( dgemqrt_obj->bufsize_c, 
                dgemqrt_obj->c, dgemqrt_obj->cref );

}

TEST_F(dgemqrt_test, dgemqrt1) {
    EXPECT_NEAR(0.0, dgemqrt_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(dgemqrt_test, dgemqrt2) {
    EXPECT_NEAR(0.0, dgemqrt_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(dgemqrt_test, dgemqrt3) {
    EXPECT_NEAR(0.0, dgemqrt_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(dgemqrt_test, dgemqrt4) {
    EXPECT_NEAR(0.0, dgemqrt_obj->diff, LAPACKE_EIG_THRESHOLD);
}


/* Begin scomplex_common_parameters  class definition */
class gemqrt_scomplex_parameters{

   public:
	int bufsize_v;
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
	lapack_int nb;
	lapack_complex_float* v, *A;
	char side;
	char trans;
	lapack_int ldc, ldv, lda_geqrt, ldt;
	/*Output Parameter*/
	lapack_complex_float* t, *c;	
	lapack_complex_float *vref, *tref, *cref, *Aref;
	/*Return Values*/	
	lapack_int info, inforef;
	lapack_int info_geqrt, inforef_geqrt;

   public:
      gemqrt_scomplex_parameters (int matrix_layout, char side , char trans, lapack_int m, lapack_int n, lapack_int nb);
      ~gemqrt_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
gemqrt_scomplex_parameters:: gemqrt_scomplex_parameters (int matrix_layout_i, char side_i , char trans_i,lapack_int m_i, lapack_int n_i, lapack_int nb_i)
{
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	side = side_i;
	trans = trans_i;
	nb = min(m,n);

	if (trans == 'T')
		trans  = 'N' ;
	
	#if LAPACKE_TEST_VERBOSE
	printf(" \n gemqrt scomplex: matrix_layout = %d, side :%c, trans:%c, m:%d, n: %d \n", matrix_layout, side, trans, m, n);
	#endif	
	
	/*Sizes based on side*/
	if (side  == 'L')
	{
		k = m;
		nb = k;
		if (matrix_layout == LAPACK_COL_MAJOR)
		{
			ldv = m;
			ldt = nb;
			ldc = m;
			lda_geqrt = m;
			bufsize_a = lda_geqrt*n;
			bufsize_v = ldv *k;
			bufsize_t = ldt*min(m,n);
			bufsize_c = ldc*n;
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{	
			ldv = k;
			ldt = k;
			ldc = n;
			lda_geqrt = n;
			bufsize_a = lda_geqrt*m;
			bufsize_v = ldv *m;
			bufsize_t = ldt*nb;
			bufsize_c = ldc*m;
		} else
		{
			EXPECT_TRUE(false) << "matrix_layout invalid";
		}
	}else if (side == 'R')
	{
		k = n;
		nb = k;
		if (matrix_layout == LAPACK_COL_MAJOR)
		{
			ldv = n;
			ldt = nb;
			ldc = m;
			lda_geqrt = m;
			bufsize_a = lda_geqrt*n;
			bufsize_v = ldv*n;
			bufsize_t = ldt*min(m,n);
			bufsize_c = ldc*n;
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{
			ldv = k;
			ldt = k;
			ldc = n;
			lda_geqrt = n;
			bufsize_a = lda_geqrt*m;
			bufsize_v = ldv*n;
			bufsize_t = ldt*nb;
			bufsize_c = ldc*m;
		} else
		{
			EXPECT_TRUE(false) << "matrix_layout invalid";
		}
	}else 
		EXPECT_TRUE(false) << "side is invalid";	
	
	bufsize_c = m*n;
	
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&v, &vref, bufsize_v);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&t, &tref, bufsize_t);	
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&c, &cref, bufsize_c);
	
	if ((v==NULL) || (vref==NULL) ||
		(t==NULL) || (tref==NULL) ||
		(c==NULL) || (cref==NULL) ||
		(A ==NULL) || (Aref ==NULL))
	{
		EXPECT_FALSE( true) << "gemqrt_float_parameters object: malloc error.";
		gemqrt_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( v, vref, bufsize_v);
	lapacke_gtest_init_scomplex_buffer_pair_rand( c, cref, bufsize_c);
	lapacke_gtest_init_scomplex_buffer_pair_rand(A, Aref, bufsize_a);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
gemqrt_scomplex_parameters :: ~gemqrt_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" gemqrt_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   gemqrt_free();

}
/*  Test fixture class definition */
class cgemqrt_test  : public  ::testing::Test {
public:
   gemqrt_scomplex_parameters  *cgemqrt_obj;
   void SetUp();
   void TearDown () { delete cgemqrt_obj;}
};

void cgemqrt_test::SetUp(){

	/* LAPACKE cgeqrt prototype */
    typedef int (*Fptr_NL_LAPACKE_cgeqrt) (int matrix_layout, lapack_int m, lapack_int n, lapack_int nb,\
	lapack_complex_float* a, lapack_int lda, lapack_complex_float* t, lapack_int ldt);
	
	 Fptr_NL_LAPACKE_cgeqrt cgeqrt; 
	 /* LAPACKE cgemqrt prototype */
    typedef int (*Fptr_NL_LAPACKE_cgemqrt) (int matrix_layout, char side, char trans, lapack_int m, lapack_int n,\
	lapack_int k, lapack_int nb, const lapack_complex_float* v, lapack_int ldv, const lapack_complex_float* t,\
	lapack_int ldt, lapack_complex_float* c, lapack_int ldc);

    Fptr_NL_LAPACKE_cgemqrt cgemqrt;

    cgemqrt_obj = new gemqrt_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].side,
						   eig_paramslist[idx].trans,
						   eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
						   eig_paramslist[idx].nb);
						   

    idx = Circular_Increment_Index(idx);

    cgemqrt_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgemqrt_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgemqrt_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgemqrt_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*cgemqrt library call */
    cgemqrt = (Fptr_NL_LAPACKE_cgemqrt)dlsym(cgemqrt_obj->hModule, "LAPACKE_cgemqrt");
    ASSERT_TRUE(cgemqrt != NULL) << "failed to get the Netlib LAPACKE_cgemqrt symbol";

	/*cgeqrt library call*/
	cgemqrt_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgemqrt_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgemqrt_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgemqrt_obj->lModule != NULL) << "Netlib lapacke handle NULL";

    cgeqrt = (Fptr_NL_LAPACKE_cgeqrt)dlsym(cgemqrt_obj->lModule, "LAPACKE_cgeqrt");
    ASSERT_TRUE(cgeqrt != NULL) << "failed to get the Netlib LAPACKE_cgeqrt symbol";    
    

    cgemqrt_obj->inforef_geqrt = cgeqrt( cgemqrt_obj->matrix_layout,cgemqrt_obj->m, cgemqrt_obj->n, min(cgemqrt_obj->m,cgemqrt_obj->n),
								cgemqrt_obj->Aref, cgemqrt_obj->lda_geqrt, cgemqrt_obj->tref, min(cgemqrt_obj->m,cgemqrt_obj->n));

    /* Compute libflame's Lapacke o/p  */
    cgemqrt_obj->info_geqrt = LAPACKE_cgeqrt( cgemqrt_obj->matrix_layout,cgemqrt_obj->m, cgemqrt_obj->n, min(cgemqrt_obj->m,cgemqrt_obj->n),
											cgemqrt_obj->A, cgemqrt_obj->lda_geqrt, cgemqrt_obj->t, min(cgemqrt_obj->m,cgemqrt_obj->n));
										
										

    if( cgemqrt_obj->info_geqrt < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cgeqrt is wrong\n", cgemqrt_obj->info_geqrt );
    }
    if( cgemqrt_obj->inforef_geqrt < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cgeqrt is wrong\n", 
        cgemqrt_obj->inforef_geqrt ); 
    }  
/*Compute cgemqrt's  o/p */
    cgemqrt_obj->inforef = cgemqrt( cgemqrt_obj->matrix_layout, cgemqrt_obj->side, cgemqrt_obj->trans,
								cgemqrt_obj->m, cgemqrt_obj->n, cgemqrt_obj->k, cgemqrt_obj->nb, (const lapack_complex_float*)cgemqrt_obj->vref,
								cgemqrt_obj->ldv, (const lapack_complex_float*)cgemqrt_obj->tref, cgemqrt_obj->ldt, cgemqrt_obj->cref, cgemqrt_obj->ldc);

    /* Compute libflame's Lapacke o/p  */
    cgemqrt_obj->info = LAPACKE_cgemqrt( cgemqrt_obj->matrix_layout, cgemqrt_obj->side, cgemqrt_obj->trans,
									cgemqrt_obj->m, cgemqrt_obj->n, cgemqrt_obj->k, cgemqrt_obj->nb, (const lapack_complex_float*)cgemqrt_obj->v,
									cgemqrt_obj->ldv, (const lapack_complex_float*)cgemqrt_obj->t, cgemqrt_obj->ldt, cgemqrt_obj->c, cgemqrt_obj->ldc);
    if( cgemqrt_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame\
        LAPACKE_cgemqrt is wrong\n", cgemqrt_obj->info );
    }
    if( cgemqrt_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cgemqrt is wrong\n", 
        cgemqrt_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps */
    cgemqrt_obj->diff =  computeDiff_c( cgemqrt_obj->bufsize_c, 
                cgemqrt_obj->c, cgemqrt_obj->cref );

}

TEST_F(cgemqrt_test, cgemqrt1) {
    EXPECT_NEAR(0.0, cgemqrt_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(cgemqrt_test, cgemqrt2) {
    EXPECT_NEAR(0.0, cgemqrt_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(cgemqrt_test, cgemqrt3) {
    EXPECT_NEAR(0.0, cgemqrt_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(cgemqrt_test, cgemqrt4) {
    EXPECT_NEAR(0.0, cgemqrt_obj->diff, LAPACKE_EIG_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class gemqrt_dcomplex_parameters{

   public:
	int bufsize_v;
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
	lapack_int nb;
	lapack_complex_double* v, *A;
	char side;
	char trans;
	lapack_int ldc, ldv, lda_geqrt, ldt;
	/*Output Parameter*/
	lapack_complex_double* t, *c;	
	lapack_complex_double *vref, *tref, *cref, *Aref;
	/*Return Values*/	
	lapack_int info, inforef;
	lapack_int info_geqrt, inforef_geqrt;

   public:
      gemqrt_dcomplex_parameters (int matrix_layout, char side , char trans, lapack_int m, lapack_int n, lapack_int nb);
      ~gemqrt_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
gemqrt_dcomplex_parameters:: gemqrt_dcomplex_parameters (int matrix_layout_i, char side_i , char trans_i,lapack_int m_i, lapack_int n_i, lapack_int nb_i)
{
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	side = side_i;
	trans = trans_i;
	nb = min(m,n);

	if (trans == 'T')
		trans  = 'N' ;
	
	#if LAPACKE_TEST_VERBOSE
	printf(" \n gemqrt dcomplex: matrix_layout = %d, side :%c, trans:%c, m:%d, n: %d \n", matrix_layout, side, trans, m, n);
	#endif	
	
	/*Sizes based on side*/
	if (side  == 'L')
	{
		k = m;
		nb = k;
		if (matrix_layout == LAPACK_COL_MAJOR)
		{
			ldv = m;
			ldt = nb;
			ldc = m;
			lda_geqrt = m;
			bufsize_a = lda_geqrt*n;
			bufsize_v = ldv *k;
			bufsize_t = ldt*min(m,n);
			bufsize_c = ldc*n;
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{	
			ldv = k;
			ldt = k;
			ldc = n;
			lda_geqrt = n;
			bufsize_a = lda_geqrt*m;
			bufsize_v = ldv *m;
			bufsize_t = ldt*nb;
			bufsize_c = ldc*m;
		} else
		{
			EXPECT_TRUE(false) << "matrix_layout invalid";
		}
	}else if (side == 'R')
	{
		k = n;
		nb = k;
		if (matrix_layout == LAPACK_COL_MAJOR)
		{
			ldv = n;
			ldt = nb;
			ldc = m;
			lda_geqrt = m;
			bufsize_a = lda_geqrt*n;
			bufsize_v = ldv*n;
			bufsize_t = ldt*min(m,n);
			bufsize_c = ldc*n;
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{
			ldv = k;
			ldt = k;
			ldc = n;
			lda_geqrt = n;
			bufsize_a = lda_geqrt*m;
			bufsize_v = ldv*n;
			bufsize_t = ldt*nb;
			bufsize_c = ldc*m;
		} else
		{
			EXPECT_TRUE(false) << "matrix_layout invalid";
		}
	}else 
		EXPECT_TRUE(false) << "side is invalid";	
	
	bufsize_c = m*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&v, &vref, bufsize_v);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&t, &tref, bufsize_t);	
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&c, &cref, bufsize_c);
	
	if ((v==NULL) || (vref==NULL) ||
		(t==NULL) || (tref==NULL) ||
		(c==NULL) || (cref==NULL) ||
		(A ==NULL) || (Aref ==NULL))
	{
		EXPECT_FALSE( true) << "gemqrt_double_parameters object: malloc error.";
		gemqrt_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( v, vref, bufsize_v);
	lapacke_gtest_init_dcomplex_buffer_pair_rand( c, cref, bufsize_c);
	lapacke_gtest_init_dcomplex_buffer_pair_rand(A, Aref, bufsize_a);

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
gemqrt_dcomplex_parameters :: ~gemqrt_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" gemqrt_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   gemqrt_free();

}
/*  Test fixture class definition */
class zgemqrt_test  : public  ::testing::Test {
public:
   gemqrt_dcomplex_parameters  *zgemqrt_obj;
   void SetUp();
   void TearDown () { delete zgemqrt_obj;}
};

void zgemqrt_test::SetUp(){

	/* LAPACKE zgeqrt prototype */
    typedef int (*Fptr_NL_LAPACKE_zgeqrt) (int matrix_layout, lapack_int m, lapack_int n, lapack_int nb,\
	lapack_complex_double* a, lapack_int lda, lapack_complex_double* t, lapack_int ldt);
	
	 Fptr_NL_LAPACKE_zgeqrt zgeqrt; 
	 /* LAPACKE zgemqrt prototype */
    typedef int (*Fptr_NL_LAPACKE_zgemqrt) (int matrix_layout, char side, char trans, lapack_int m, lapack_int n,\
	lapack_int k, lapack_int nb, const lapack_complex_double* v, lapack_int ldv, const lapack_complex_double* t,\
	lapack_int ldt, lapack_complex_double* c, lapack_int ldc);

    Fptr_NL_LAPACKE_zgemqrt zgemqrt;

    zgemqrt_obj = new gemqrt_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].side,
						   eig_paramslist[idx].trans,
						   eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
						   eig_paramslist[idx].nb);
						   

    idx = Circular_Increment_Index(idx);

    zgemqrt_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgemqrt_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgemqrt_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgemqrt_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*zgemqrt library call */
    zgemqrt = (Fptr_NL_LAPACKE_zgemqrt)dlsym(zgemqrt_obj->hModule, "LAPACKE_zgemqrt");
    ASSERT_TRUE(zgemqrt != NULL) << "failed to get the Netlib LAPACKE_zgemqrt symbol";

	/*zgeqrt library call*/
	zgemqrt_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgemqrt_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgemqrt_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgemqrt_obj->lModule != NULL) << "Netlib lapacke handle NULL";

    zgeqrt = (Fptr_NL_LAPACKE_zgeqrt)dlsym(zgemqrt_obj->lModule, "LAPACKE_zgeqrt");
    ASSERT_TRUE(zgeqrt != NULL) << "failed to get the Netlib LAPACKE_zgeqrt symbol";    
    

    zgemqrt_obj->inforef_geqrt = zgeqrt( zgemqrt_obj->matrix_layout,zgemqrt_obj->m, zgemqrt_obj->n, min(zgemqrt_obj->m,zgemqrt_obj->n),
								zgemqrt_obj->Aref, zgemqrt_obj->lda_geqrt, zgemqrt_obj->tref, min(zgemqrt_obj->m,zgemqrt_obj->n));

    /* Compute libflame's Lapacke o/p  */
    zgemqrt_obj->info_geqrt = LAPACKE_zgeqrt( zgemqrt_obj->matrix_layout,zgemqrt_obj->m, zgemqrt_obj->n, min(zgemqrt_obj->m,zgemqrt_obj->n),
											zgemqrt_obj->A, zgemqrt_obj->lda_geqrt, zgemqrt_obj->t, min(zgemqrt_obj->m,zgemqrt_obj->n));
										
										

    if( zgemqrt_obj->info_geqrt < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zgeqrt is wrong\n", zgemqrt_obj->info_geqrt );
    }
    if( zgemqrt_obj->inforef_geqrt < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zgeqrt is wrong\n", 
        zgemqrt_obj->inforef_geqrt ); 
    }  
/*Compute zgemqrt's  o/p */
    zgemqrt_obj->inforef = zgemqrt( zgemqrt_obj->matrix_layout, zgemqrt_obj->side, zgemqrt_obj->trans,
								zgemqrt_obj->m, zgemqrt_obj->n, zgemqrt_obj->k, zgemqrt_obj->nb, (const lapack_complex_double*)zgemqrt_obj->vref,
								zgemqrt_obj->ldv, (const lapack_complex_double*)zgemqrt_obj->tref, zgemqrt_obj->ldt, zgemqrt_obj->cref, zgemqrt_obj->ldc);

    /* Compute libflame's Lapacke o/p  */
    zgemqrt_obj->info = LAPACKE_zgemqrt( zgemqrt_obj->matrix_layout, zgemqrt_obj->side, zgemqrt_obj->trans,
									zgemqrt_obj->m, zgemqrt_obj->n, zgemqrt_obj->k, zgemqrt_obj->nb, (const lapack_complex_double*)zgemqrt_obj->v,
									zgemqrt_obj->ldv, (const lapack_complex_double*)zgemqrt_obj->t, zgemqrt_obj->ldt, zgemqrt_obj->c, zgemqrt_obj->ldc);
    if( zgemqrt_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame\
        LAPACKE_zgemqrt is wrong\n", zgemqrt_obj->info );
    }
    if( zgemqrt_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zgemqrt is wrong\n", 
        zgemqrt_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps */
    zgemqrt_obj->diff =  computeDiff_z( zgemqrt_obj->bufsize_c, 
                zgemqrt_obj->c, zgemqrt_obj->cref );

}

TEST_F(zgemqrt_test, zgemqrt1) {
    EXPECT_NEAR(0.0, zgemqrt_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zgemqrt_test, zgemqrt2) {
    EXPECT_NEAR(0.0, zgemqrt_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zgemqrt_test, zgemqrt3) {
    EXPECT_NEAR(0.0, zgemqrt_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zgemqrt_test, zgemqrt4) {
    EXPECT_NEAR(0.0, zgemqrt_obj->diff, LAPACKE_EIG_THRESHOLD);
}