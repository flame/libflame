#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define tpmqrt_free() \
if (v!=NULL)    free(v); \
if (vref!=NULL) free(vref);\
if (t!=NULL)  free(t);\
if (tref!=NULL) free(tref); \
if (t_tpqrt!=NULL)  free(t_tpqrt);\
if (t_tpqrtref!=NULL) free(t_tpqrtref); \
if (b != NULL) free(b); \
if (bref != NULL) free(bref); \
if (A != NULL) free(A); \
if (Aref != NULL) free(Aref); \
if (A_tpqrt != NULL) free(A_tpqrt); \
if (A_tpqrtref != NULL) free(A_tpqrtref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule);\
if( bModule != NULL) dlclose(bModule); \
if( lModule != NULL) dlclose(lModule); 
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
/* Begin float_common_parameters  class definition */
class tpmqrt_float_parameters{

   public:
	int bufsize_v;
	int bufsize_a, bufsize_atpqrt;
	int bufsize_b;
	int bufsize_t, bufsize_tpqrt;
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
	lapack_int nb, nb_tpqrt;
	lapack_int lda;
	float* v, *A, *A_tpqrt, *t_tpqrt;
	char side;
	char trans;
	lapack_int ldb, ldv, lda_tpqrt, ldt, ldt_tpqrt;
	/*Output Parameter*/
	float* t, *b;	
	float *vref, *tref, *bref, *Aref, *A_tpqrtref, *t_tpqrtref;
	/*Return Values*/	
	lapack_int info, inforef;
	lapack_int info_tpqrt, inforef_tpqrt;

   public:
      tpmqrt_float_parameters (int matrix_layout, char side , char trans, lapack_int m, lapack_int n, lapack_int nb);
      ~tpmqrt_float_parameters ();

};

/* Constructor definition  float_common_parameters */
tpmqrt_float_parameters:: tpmqrt_float_parameters (int matrix_layout_i, char side_i , char trans_i,lapack_int m_i, lapack_int n_i, lapack_int nb_i)
{
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	side = side_i;
	trans = trans_i;
	nb = nb_i;
	l  = fla_min(m,n);
	
	if (trans == 'C')
		trans  = 'N' ;
	
	#if LAPACKE_TEST_VERBOSE
		printf(" \n tpmqrt float: matrix_layout = %d, side :%b, trans:%b, m:%d, n: %d \n", matrix_layout, side, trans, m, n);
	#endif
	if (side == 'L')
	{
		k = m;
		//nb = k;
		if (matrix_layout == LAPACK_COL_MAJOR)
		{
			lda = k;
			ldv = m;
			ldt = nb;
			ldt_tpqrt = nb;
			ldb = m;
			bufsize_v = ldv*k;
			bufsize_t =  ldt*k;
			bufsize_a = lda*m;
			bufsize_b = ldb*n;
			bufsize_tpqrt = ldt_tpqrt*n;
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{		 
			lda= m;
			ldv = k;
			ldt = k;
			ldt_tpqrt = n;
			ldb = n;
			bufsize_v = ldv*m;
			bufsize_t = ldt*nb;
			bufsize_a = lda*k;
			bufsize_b = ldb*m;
			bufsize_tpqrt = ldt_tpqrt*nb;
		}else
		{
			EXPECT_TRUE(false) << "matrix_layout invalid";
		}
	}else if (side == 'R') {
		k = m;
		//nb = k;
		if (matrix_layout == LAPACK_COL_MAJOR)
		{	lda= m;
			ldv = n;
			ldt = nb;
			ldb = m;
			ldt_tpqrt = nb;
			bufsize_v = ldv*n;
			bufsize_t =  ldt*k;
			bufsize_a = lda*k;;
			bufsize_b = ldb*n;
			bufsize_tpqrt = ldt_tpqrt*n;
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{	lda= m;
			ldv = k;
			ldt = k;
			ldb =n;
			ldt_tpqrt = n;
			bufsize_v = ldv*n;
			bufsize_t = ldt*nb;
			bufsize_a = lda*m;
			bufsize_b = ldb*m;
			bufsize_tpqrt = ldt_tpqrt*nb;
		}else
		{
			EXPECT_TRUE(false) << "matrix_layout invalid";
		}
	}else 
		EXPECT_TRUE(false) << "side is invalid";
	
	lda_tpqrt =n;
	bufsize_atpqrt =  lda_tpqrt*n;
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&v, &vref, bufsize_v);
	lapacke_gtest_alloc_float_buffer_pair(&t, &tref, bufsize_t);	
	lapacke_gtest_alloc_float_buffer_pair(&b, &bref, bufsize_b);
	lapacke_gtest_alloc_float_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_float_buffer_pair(&A_tpqrt, &A_tpqrtref, bufsize_atpqrt);
	lapacke_gtest_alloc_float_buffer_pair(&t_tpqrt, &t_tpqrtref, bufsize_tpqrt);
	
	if ((v==NULL) || (vref==NULL) ||
		(t==NULL) || (tref==NULL) ||
		(b==NULL) || (bref==NULL) ||
		(A==NULL) || (Aref==NULL) ||
		(A_tpqrt == NULL) || (A_tpqrtref == NULL) ||
		(t_tpqrt == NULL) || (t_tpqrtref == NULL))
	{
		EXPECT_FALSE( true) << "tpmqrt_float_parameters object: malloc error.";
		tpmqrt_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand( v, vref, bufsize_v);
	lapacke_gtest_init_float_buffer_pair_rand( A, Aref, bufsize_a);
	lapacke_gtest_init_float_buffer_pair_rand( t, tref, bufsize_t);
	lapacke_gtest_init_float_buffer_pair_rand( b, bref, bufsize_b);
	lapacke_gtest_init_float_buffer_pair_rand( A_tpqrt, A_tpqrtref, bufsize_atpqrt);
	

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
tpmqrt_float_parameters :: ~tpmqrt_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" tpmqrt_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   tpmqrt_free();

}
/*  Test fixture class definition */
class stpmqrt_test  : public  ::testing::Test {
public:
   tpmqrt_float_parameters  *stpmqrt_obj;
   void SetUp();
   void TearDown () { delete stpmqrt_obj;}
};

void stpmqrt_test::SetUp(){

	/* LAPACKE stpqrt prototype */
    typedef int (*Fptr_NL_LAPACKE_stpqrt) (int matrix_layout, lapack_int m, lapack_int n,\
                           lapack_int l, lapack_int nb,\
                           float* a, lapack_int lda,\
                           float* b, lapack_int ldb,\
                           float* t, lapack_int ldt);
	
	 Fptr_NL_LAPACKE_stpqrt stpqrt;
	 /* LAPACKE stpmqrt prototype */
    typedef int (*Fptr_NL_LAPACKE_stpmqrt) (int matrix_layout, char side, char trans,\
                            lapack_int m, lapack_int n, lapack_int k,\
                            lapack_int l, lapack_int nb,\
                            const float* v, lapack_int ldv,\
                            const float* t, lapack_int ldt,\
                            float* a, lapack_int lda,\
                            float* b, lapack_int ldb );

    Fptr_NL_LAPACKE_stpmqrt stpmqrt;

    stpmqrt_obj = new tpmqrt_float_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].side,
						   eig_paramslist[idx].trans,
						   eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
						   eig_paramslist[idx].nb);
						   

    idx = Circular_Increment_Index(idx);

    stpmqrt_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    stpmqrt_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(stpmqrt_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(stpmqrt_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*stpmqrt library call */
    stpmqrt = (Fptr_NL_LAPACKE_stpmqrt)dlsym(stpmqrt_obj->hModule, "LAPACKE_stpmqrt");
    ASSERT_TRUE(stpmqrt != NULL) << "failed to get the Netlib LAPACKE_stpmqrt symbol";

	/*stpqrt library call*/
	stpmqrt_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    stpmqrt_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(stpmqrt_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(stpmqrt_obj->lModule != NULL) << "Netlib lapacke handle NULL";

    stpqrt = (Fptr_NL_LAPACKE_stpqrt)dlsym(stpmqrt_obj->lModule, "LAPACKE_stpqrt");
    ASSERT_TRUE(stpqrt != NULL) << "failed to get the Netlib LAPACKE_stpqrt symbol";    
    

    stpmqrt_obj->inforef_tpqrt = stpqrt( stpmqrt_obj->matrix_layout,stpmqrt_obj->m, stpmqrt_obj->n, stpmqrt_obj->l, stpmqrt_obj->nb,
								stpmqrt_obj->A_tpqrtref, stpmqrt_obj->lda_tpqrt, stpmqrt_obj->bref, stpmqrt_obj->ldb, stpmqrt_obj->t_tpqrtref, stpmqrt_obj->ldt_tpqrt);

    /* Compute libflame's Lapacke o/p  */
    stpmqrt_obj->info_tpqrt = LAPACKE_stpqrt( stpmqrt_obj->matrix_layout,stpmqrt_obj->m, stpmqrt_obj->n, stpmqrt_obj->l, stpmqrt_obj->nb,
											stpmqrt_obj->A_tpqrt, stpmqrt_obj->lda_tpqrt, stpmqrt_obj->b, stpmqrt_obj->ldb, stpmqrt_obj->t_tpqrt, stpmqrt_obj->ldt_tpqrt);
										
										

    if( stpmqrt_obj->info_tpqrt < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_stpqrt is wrong\n", stpmqrt_obj->info_tpqrt );
    }
    if( stpmqrt_obj->inforef_tpqrt < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_stpqrt is wrong\n", 
        stpmqrt_obj->inforef_tpqrt );
    }  

/*Copy matrices from geqrt to gemqrt*
	memcpy((void*)stpmqrt_obj->tref, stpmqrt_obj->t_tpqrtref, stpmqrt_obj->bufsize_t);
	memcpy((void*)stpmqrt_obj->t, stpmqrt_obj->t_tpqrt, stpmqrt_obj->bufsize_t);
	
	memcpy((void*)stpmqrt_obj->vref, stpmqrt_obj->A_tpqrtref, stpmqrt_obj->bufsize_v);
	memcpy((void*)stpmqrt_obj->v, stpmqrt_obj->A_tpqrt, stpmqrt_obj->bufsize_v);*/
	
	
	
	
	/*Compute stpmqrt's  o/p */
	
    stpmqrt_obj->inforef = stpmqrt( stpmqrt_obj->matrix_layout, stpmqrt_obj->side, stpmqrt_obj->trans,\
								stpmqrt_obj->m, stpmqrt_obj->n, stpmqrt_obj->k, stpmqrt_obj->l,\
								stpmqrt_obj->nb, (const float*)stpmqrt_obj->vref,\
								stpmqrt_obj->ldv, (const float*)stpmqrt_obj->tref,\
								stpmqrt_obj->ldt, stpmqrt_obj->Aref, stpmqrt_obj->lda, stpmqrt_obj->bref, stpmqrt_obj->ldb);

    /* Compute libflame's Lapacke o/p  */
    stpmqrt_obj->info = LAPACKE_stpmqrt( stpmqrt_obj->matrix_layout, stpmqrt_obj->side, stpmqrt_obj->trans,\
									stpmqrt_obj->m, stpmqrt_obj->n, stpmqrt_obj->k, stpmqrt_obj->l,\
									stpmqrt_obj->nb, (const float*)stpmqrt_obj->v,\
									stpmqrt_obj->ldv, (const float*)stpmqrt_obj->t,\ 
									stpmqrt_obj->ldt, stpmqrt_obj->A, stpmqrt_obj->lda, stpmqrt_obj->b, stpmqrt_obj->ldb);
    if( stpmqrt_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame\
        LAPACKE_stpmqrt is wrong\n", stpmqrt_obj->info );
    }
    if( stpmqrt_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_stpmqrt is wrong\n", 
        stpmqrt_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    stpmqrt_obj->diff =  computeDiff_s( stpmqrt_obj->bufsize_a, 
                stpmqrt_obj->A, stpmqrt_obj->Aref );

}

TEST_F(stpmqrt_test, stpmqrt1) {
    EXPECT_NEAR(0.0, stpmqrt_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(stpmqrt_test, stpmqrt2) {
    EXPECT_NEAR(0.0, stpmqrt_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(stpmqrt_test, stpmqrt3) {
    EXPECT_NEAR(0.0, stpmqrt_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(stpmqrt_test, stpmqrt4) {
    EXPECT_NEAR(0.0, stpmqrt_obj->diff, LAPACKE_EIG_THRESHOLD);
}

/* Begin double_common_parameters  class definition */
class tpmqrt_double_parameters{

   public:
	int bufsize_v;
	int bufsize_a, bufsize_atpqrt;
	int bufsize_b;
	int bufsize_t, bufsize_tpqrt;
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
	lapack_int nb, nb_tpqrt;
	lapack_int lda;
	double* v, *A, *A_tpqrt, *t_tpqrt;
	char side;
	char trans;
	lapack_int ldb, ldv, lda_tpqrt, ldt, ldt_tpqrt;
	/*Output Parameter*/
	double* t, *b;	
	double *vref, *tref, *bref, *Aref, *A_tpqrtref, *t_tpqrtref;
	/*Return Values*/	
	lapack_int info, inforef;
	lapack_int info_tpqrt, inforef_tpqrt;

   public:
      tpmqrt_double_parameters (int matrix_layout, char side , char trans, lapack_int m, lapack_int n, lapack_int nb);
      ~tpmqrt_double_parameters ();

};

/* Constructor definition  double_common_parameters */
tpmqrt_double_parameters:: tpmqrt_double_parameters (int matrix_layout_i, char side_i , char trans_i,lapack_int m_i, lapack_int n_i, lapack_int nb_i)
{
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	side = side_i;
	trans = trans_i;
	nb = nb_i;
	l  = fla_min(m,n);
	k = m;
	
	if (trans == 'C')
		trans  = 'N' ;
	
	#if LAPACKE_TEST_VERBOSE
		printf(" \n tpmqrt double: matrix_layout = %d, side :%b, trans:%b, m:%d, n: %d \n", matrix_layout, side, trans, m, n);
	#endif
	if (side == 'L')
	{
		//k = m;
		//nb = k;
		if (matrix_layout == LAPACK_COL_MAJOR)
		{
			lda = k;
			ldv = m;
			ldt = nb;
			ldt_tpqrt = nb;
			ldb = m;
			bufsize_v = ldv*k;
			bufsize_t =  ldt*k;
			bufsize_a = lda*n;
			bufsize_b = ldb*n;
			bufsize_tpqrt = ldt_tpqrt*n;
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{		 
			lda= m;
			ldv = k;
			ldt = k;
			ldt_tpqrt = n;
			ldb = n;
			bufsize_v = ldv*m;
			bufsize_t = ldt*nb;
			bufsize_a = lda*k;
			bufsize_b = ldb*m;
			bufsize_tpqrt = ldt_tpqrt*nb;
		}else
		{
			EXPECT_TRUE(false) << "matrix_layout invalid";
		}
	}else if (side == 'R') {
		//k = n;
		//nb = k;
		if (matrix_layout == LAPACK_COL_MAJOR)
		{	lda= m;
			ldv = n;
			ldt = nb;
			ldb = m;
			ldt_tpqrt = nb;
			bufsize_v = ldv*n;
			bufsize_t =  ldt*k;
			bufsize_a = lda*k;;
			bufsize_b = ldb*n;
			bufsize_tpqrt = ldt_tpqrt*n;
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{	lda= k;
			ldv = k;
			ldt = k;
			ldb =n;
			ldt_tpqrt = n;
			bufsize_v = ldv*n;
			bufsize_t = ldt*nb;
			bufsize_a = lda*m;
			bufsize_b = ldb*m;
			bufsize_tpqrt = ldt_tpqrt*nb;
		}else
		{
			EXPECT_TRUE(false) << "matrix_layout invalid";
		}
	}else 
		EXPECT_TRUE(false) << "side is invalid";
	
	lda_tpqrt =n;
	bufsize_atpqrt =  lda_tpqrt*n;
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&v, &vref, bufsize_v);
	lapacke_gtest_alloc_double_buffer_pair(&t, &tref, bufsize_t);	
	lapacke_gtest_alloc_double_buffer_pair(&b, &bref, bufsize_b);
	lapacke_gtest_alloc_double_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_double_buffer_pair(&A_tpqrt, &A_tpqrtref, bufsize_atpqrt);
	lapacke_gtest_alloc_double_buffer_pair(&t_tpqrt, &t_tpqrtref, bufsize_tpqrt);
	
	if ((v==NULL) || (vref==NULL) ||
		(t==NULL) || (tref==NULL) ||
		(b==NULL) || (bref==NULL) ||
		(A==NULL) || (Aref==NULL) ||
		(A_tpqrt == NULL) || (A_tpqrtref == NULL) ||
		(t_tpqrt == NULL) || (t_tpqrtref == NULL))
	{
		EXPECT_FALSE( true) << "tpmqrt_double_parameters object: malloc error.";
		tpmqrt_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand( v, vref, bufsize_v);
	lapacke_gtest_init_double_buffer_pair_rand( A, Aref, bufsize_a);
	lapacke_gtest_init_double_buffer_pair_rand( t, tref, bufsize_t);
	lapacke_gtest_init_double_buffer_pair_rand( b, bref, bufsize_b);
	lapacke_gtest_init_double_buffer_pair_rand( A_tpqrt, A_tpqrtref, bufsize_atpqrt);
	

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
tpmqrt_double_parameters :: ~tpmqrt_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" tpmqrt_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   tpmqrt_free();

}
/*  Test fixture class definition */
class dtpmqrt_test  : public  ::testing::Test {
public:
   tpmqrt_double_parameters  *dtpmqrt_obj;
   void SetUp();
   void TearDown () { delete dtpmqrt_obj;}
};

void dtpmqrt_test::SetUp(){

	/* LAPACKE dtpqrt prototype */
    typedef int (*Fptr_NL_LAPACKE_dtpqrt) (int matrix_layout, lapack_int m, lapack_int n,\
                           lapack_int l, lapack_int nb,\
                           double* a, lapack_int lda,\
                           double* b, lapack_int ldb,\
                           double* t, lapack_int ldt);
	
	 Fptr_NL_LAPACKE_dtpqrt dtpqrt;
	 /* LAPACKE dtpmqrt prototype */
    typedef int (*Fptr_NL_LAPACKE_dtpmqrt) (int matrix_layout, char side, char trans,\
                            lapack_int m, lapack_int n, lapack_int k,\
                            lapack_int l, lapack_int nb,\
                            const double* v, lapack_int ldv,\
                            const double* t, lapack_int ldt,\
                            double* a, lapack_int lda,\
                            double* b, lapack_int ldb );

    Fptr_NL_LAPACKE_dtpmqrt dtpmqrt;

    dtpmqrt_obj = new tpmqrt_double_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].side,
						   eig_paramslist[idx].trans,
						   eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
						   eig_paramslist[idx].nb);
						   

    idx = Circular_Increment_Index(idx);

    dtpmqrt_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dtpmqrt_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dtpmqrt_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dtpmqrt_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*dtpmqrt library call */
    dtpmqrt = (Fptr_NL_LAPACKE_dtpmqrt)dlsym(dtpmqrt_obj->hModule, "LAPACKE_dtpmqrt");
    ASSERT_TRUE(dtpmqrt != NULL) << "failed to get the Netlib LAPACKE_dtpmqrt symbol";

	/*dtpqrt library call*/
	dtpmqrt_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dtpmqrt_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dtpmqrt_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dtpmqrt_obj->lModule != NULL) << "Netlib lapacke handle NULL";

    dtpqrt = (Fptr_NL_LAPACKE_dtpqrt)dlsym(dtpmqrt_obj->lModule, "LAPACKE_dtpqrt");
    ASSERT_TRUE(dtpqrt != NULL) << "failed to get the Netlib LAPACKE_dtpqrt symbol";    
    

    dtpmqrt_obj->inforef_tpqrt = dtpqrt( dtpmqrt_obj->matrix_layout,dtpmqrt_obj->m, dtpmqrt_obj->n, dtpmqrt_obj->l, dtpmqrt_obj->nb,
								dtpmqrt_obj->A_tpqrtref, dtpmqrt_obj->lda_tpqrt, dtpmqrt_obj->bref, dtpmqrt_obj->ldb, dtpmqrt_obj->t_tpqrtref, dtpmqrt_obj->ldt_tpqrt);

    /* Compute libflame's Lapacke o/p  */
    dtpmqrt_obj->info_tpqrt = LAPACKE_dtpqrt( dtpmqrt_obj->matrix_layout,dtpmqrt_obj->m, dtpmqrt_obj->n, dtpmqrt_obj->l, dtpmqrt_obj->nb,
											dtpmqrt_obj->A_tpqrt, dtpmqrt_obj->lda_tpqrt, dtpmqrt_obj->b, dtpmqrt_obj->ldb, dtpmqrt_obj->t_tpqrt, dtpmqrt_obj->ldt_tpqrt);
										
										

    if( dtpmqrt_obj->info_tpqrt < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dtpqrt is wrong\n", dtpmqrt_obj->info_tpqrt );
    }
    if( dtpmqrt_obj->inforef_tpqrt < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dtpqrt is wrong\n", 
        dtpmqrt_obj->inforef_tpqrt );
    }

	/*Copy matrices from geqrt to gemqrt*
	memcpy((void*)dtpmqrt_obj->tref, dtpmqrt_obj->t_tpqrtref, dtpmqrt_obj->bufsize_t);
	memcpy((void*)dtpmqrt_obj->t, dtpmqrt_obj->t_tpqrt, dtpmqrt_obj->bufsize_t);
	
	memcpy((void*)dtpmqrt_obj->vref, dtpmqrt_obj->A_tpqrtref, dtpmqrt_obj->bufsize_v);
	memcpy((void*)dtpmqrt_obj->v, dtpmqrt_obj->A_tpqrt, dtpmqrt_obj->bufsize_v);*/
	
/*Compute dtpmqrt's  o/p */
    dtpmqrt_obj->inforef = dtpmqrt( dtpmqrt_obj->matrix_layout, dtpmqrt_obj->side, dtpmqrt_obj->trans,\
								dtpmqrt_obj->m, dtpmqrt_obj->n, dtpmqrt_obj->k, dtpmqrt_obj->l,\
								dtpmqrt_obj->nb, (const double*)dtpmqrt_obj->vref,\
								dtpmqrt_obj->ldv, (const double*)dtpmqrt_obj->tref,\
								dtpmqrt_obj->ldt, dtpmqrt_obj->Aref, dtpmqrt_obj->lda, dtpmqrt_obj->bref, dtpmqrt_obj->ldb);

    /* Compute libflame's Lapacke o/p  */
    dtpmqrt_obj->info = LAPACKE_dtpmqrt( dtpmqrt_obj->matrix_layout, dtpmqrt_obj->side, dtpmqrt_obj->trans,\
									dtpmqrt_obj->m, dtpmqrt_obj->n, dtpmqrt_obj->k, dtpmqrt_obj->l,\
									dtpmqrt_obj->nb, (const double*)dtpmqrt_obj->v,\
									dtpmqrt_obj->ldv, (const double*)dtpmqrt_obj->t,\ 
									dtpmqrt_obj->ldt, dtpmqrt_obj->A, dtpmqrt_obj->lda, dtpmqrt_obj->b, dtpmqrt_obj->ldb);
    if( dtpmqrt_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame\
        LAPACKE_dtpmqrt is wrong\n", dtpmqrt_obj->info );
    }
    if( dtpmqrt_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dtpmqrt is wrong\n", 
        dtpmqrt_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dtpmqrt_obj->diff =  computeDiff_d( dtpmqrt_obj->bufsize_a, 
                dtpmqrt_obj->A, dtpmqrt_obj->Aref );

}

TEST_F(dtpmqrt_test, dtpmqrt1) {
    EXPECT_NEAR(0.0, dtpmqrt_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(dtpmqrt_test, dtpmqrt2) {
    EXPECT_NEAR(0.0, dtpmqrt_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(dtpmqrt_test, dtpmqrt3) {
    EXPECT_NEAR(0.0, dtpmqrt_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(dtpmqrt_test, dtpmqrt4) {
    EXPECT_NEAR(0.0, dtpmqrt_obj->diff, LAPACKE_EIG_THRESHOLD);
}

/* Begin scomplex_common_parameters  class definition */
class tpmqrt_scomplex_parameters{

   public:
	int bufsize_v;
	int bufsize_a, bufsize_atpqrt;
	int bufsize_b;
	int bufsize_t, bufsize_tpqrt;
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
	lapack_int nb, nb_tpqrt;
	lapack_int lda;
	lapack_complex_float* v, *A, *A_tpqrt, *t_tpqrt;
	char side;
	char trans;
	lapack_int ldb, ldv, lda_tpqrt, ldt, ldt_tpqrt;
	/*Output Parameter*/
	lapack_complex_float* t, *b;	
	lapack_complex_float *vref, *tref, *bref, *Aref, *A_tpqrtref, *t_tpqrtref;
	/*Return Values*/	
	lapack_int info, inforef;
	lapack_int info_tpqrt, inforef_tpqrt;

   public:
      tpmqrt_scomplex_parameters (int matrix_layout, char side , char trans, lapack_int m, lapack_int n, lapack_int nb);
      ~tpmqrt_scomplex_parameters ();

};

/* Constructor definition  scomplex_common_parameters */
tpmqrt_scomplex_parameters:: tpmqrt_scomplex_parameters (int matrix_layout_i, char side_i , char trans_i,lapack_int m_i, lapack_int n_i, lapack_int nb_i)
{
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	side = side_i;
	trans = trans_i;
	nb = nb_i;
	l  = fla_min(m,n);
	k = m;
	
	if (trans == 'T')
		trans  = 'N' ;
	
	#if LAPACKE_TEST_VERBOSE
		printf(" \n tpmqrt scomplex: matrix_layout = %d, side :%b, trans:%b, m:%d, n: %d \n", matrix_layout, side, trans, m, n);
	#endif
	if (side == 'L')
	{
		//k = m;
		//nb = k;
		if (matrix_layout == LAPACK_COL_MAJOR)
		{
			lda = k;
			ldv = m;
			ldt = nb;
			ldt_tpqrt = nb;
			ldb = m;
			bufsize_v = ldv*k;
			bufsize_t =  ldt*k;
			bufsize_a = lda*n;
			bufsize_b = ldb*n;
			bufsize_tpqrt = ldt_tpqrt*n;
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{		 
			lda= m;
			ldv = k;
			ldt = k;
			ldt_tpqrt = n;
			ldb = n;
			bufsize_v = ldv*m;
			bufsize_t = ldt*nb;
			bufsize_a = lda*k;
			bufsize_b = ldb*m;
			bufsize_tpqrt = ldt_tpqrt*nb;
		}else
		{
			EXPECT_TRUE(false) << "matrix_layout invalid";
		}
	}else if (side == 'R') {
		//k = n;
		//nb = k;
		if (matrix_layout == LAPACK_COL_MAJOR)
		{	lda= m;
			ldv = n;
			ldt = nb;
			ldb = m;
			ldt_tpqrt = nb;
			bufsize_v = ldv*n;
			bufsize_t =  ldt*k;
			bufsize_a = lda*k;;
			bufsize_b = ldb*n;
			bufsize_tpqrt = ldt_tpqrt*n;
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{	lda= k;
			ldv = k;
			ldt = k;
			ldb =n;
			ldt_tpqrt = n;
			bufsize_v = ldv*n;
			bufsize_t = ldt*nb;
			bufsize_a = lda*m;
			bufsize_b = ldb*m;
			bufsize_tpqrt = ldt_tpqrt*nb;
		}else
		{
			EXPECT_TRUE(false) << "matrix_layout invalid";
		}
	}else 
		EXPECT_TRUE(false) << "side is invalid";
	
	lda_tpqrt =n;
	bufsize_atpqrt =  lda_tpqrt*n;
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&v, &vref, bufsize_v);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&t, &tref, bufsize_t);	
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&b, &bref, bufsize_b);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A_tpqrt, &A_tpqrtref, bufsize_atpqrt);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&t_tpqrt, &t_tpqrtref, bufsize_tpqrt);
	
	if ((v==NULL) || (vref==NULL) ||
		(t==NULL) || (tref==NULL) ||
		(b==NULL) || (bref==NULL) ||
		(A==NULL) || (Aref==NULL) ||
		(A_tpqrt == NULL) || (A_tpqrtref == NULL) ||
		(t_tpqrt == NULL) || (t_tpqrtref == NULL))
	{
		EXPECT_FALSE( true) << "tpmqrt_scomplex_parameters object: malloc error.";
		tpmqrt_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( v, vref, bufsize_v);
	lapacke_gtest_init_scomplex_buffer_pair_rand( A, Aref, bufsize_a);
	lapacke_gtest_init_scomplex_buffer_pair_rand( t, tref, bufsize_t);
	lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, bufsize_b);
	lapacke_gtest_init_scomplex_buffer_pair_rand( A_tpqrt, A_tpqrtref, bufsize_atpqrt);
	

} /* end of Constructor  */

/* Destructor definition  'scomplex_common_parameters' class  */
tpmqrt_scomplex_parameters :: ~tpmqrt_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" tpmqrt_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   tpmqrt_free();

}
/*  Test fixture class definition */
class ctpmqrt_test  : public  ::testing::Test {
public:
   tpmqrt_scomplex_parameters  *ctpmqrt_obj;
   void SetUp();
   void TearDown () { delete ctpmqrt_obj;}
};

void ctpmqrt_test::SetUp(){

	/* LAPACKE ctpqrt prototype */
    typedef int (*Fptr_NL_LAPACKE_ctpqrt) (int matrix_layout, lapack_int m, lapack_int n,\
                           lapack_int l, lapack_int nb,\
                           lapack_complex_float* a, lapack_int lda,\
                           lapack_complex_float* b, lapack_int ldb,\
                           lapack_complex_float* t, lapack_int ldt);
	
	 Fptr_NL_LAPACKE_ctpqrt ctpqrt;
	 /* LAPACKE ctpmqrt prototype */
    typedef int (*Fptr_NL_LAPACKE_ctpmqrt) (int matrix_layout, char side, char trans,\
                            lapack_int m, lapack_int n, lapack_int k,\
                            lapack_int l, lapack_int nb,\
                            const lapack_complex_float* v, lapack_int ldv,\
                            const lapack_complex_float* t, lapack_int ldt,\
                            lapack_complex_float* a, lapack_int lda,\
                            lapack_complex_float* b, lapack_int ldb );

    Fptr_NL_LAPACKE_ctpmqrt ctpmqrt;

    ctpmqrt_obj = new tpmqrt_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].side,
						   eig_paramslist[idx].trans,
						   eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
						   eig_paramslist[idx].nb);
						   

    idx = Circular_Increment_Index(idx);

    ctpmqrt_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ctpmqrt_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ctpmqrt_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ctpmqrt_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*ctpmqrt library call */
    ctpmqrt = (Fptr_NL_LAPACKE_ctpmqrt)dlsym(ctpmqrt_obj->hModule, "LAPACKE_ctpmqrt");
    ASSERT_TRUE(ctpmqrt != NULL) << "failed to get the Netlib LAPACKE_ctpmqrt symbol";

	/*ctpqrt library call*/
	ctpmqrt_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ctpmqrt_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ctpmqrt_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ctpmqrt_obj->lModule != NULL) << "Netlib lapacke handle NULL";

    ctpqrt = (Fptr_NL_LAPACKE_ctpqrt)dlsym(ctpmqrt_obj->lModule, "LAPACKE_ctpqrt");
    ASSERT_TRUE(ctpqrt != NULL) << "failed to get the Netlib LAPACKE_ctpqrt symbol";    
    

    ctpmqrt_obj->inforef_tpqrt = ctpqrt( ctpmqrt_obj->matrix_layout,ctpmqrt_obj->m, ctpmqrt_obj->n, ctpmqrt_obj->l, ctpmqrt_obj->nb,
								ctpmqrt_obj->A_tpqrtref, ctpmqrt_obj->lda_tpqrt, ctpmqrt_obj->bref, ctpmqrt_obj->ldb, ctpmqrt_obj->t_tpqrtref, ctpmqrt_obj->ldt_tpqrt);

    /* Compute libflame's Lapacke o/p  */
    ctpmqrt_obj->info_tpqrt = LAPACKE_ctpqrt( ctpmqrt_obj->matrix_layout,ctpmqrt_obj->m, ctpmqrt_obj->n, ctpmqrt_obj->l, ctpmqrt_obj->nb,
											ctpmqrt_obj->A_tpqrt, ctpmqrt_obj->lda_tpqrt, ctpmqrt_obj->b, ctpmqrt_obj->ldb, ctpmqrt_obj->t_tpqrt, ctpmqrt_obj->ldt_tpqrt);
										
										

    if( ctpmqrt_obj->info_tpqrt < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ctpqrt is wrong\n", ctpmqrt_obj->info_tpqrt );
    }
    if( ctpmqrt_obj->inforef_tpqrt < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ctpqrt is wrong\n", 
        ctpmqrt_obj->inforef_tpqrt );
    }
	/*Copy matrices from geqrt to gemqrt*
	memcpy((void*)ctpmqrt_obj->tref, ctpmqrt_obj->t_tpqrtref, ctpmqrt_obj->bufsize_t);
	memcpy((void*)ctpmqrt_obj->t, ctpmqrt_obj->t_tpqrt, ctpmqrt_obj->bufsize_t);
	
	memcpy((void*)ctpmqrt_obj->vref, ctpmqrt_obj->A_tpqrtref, ctpmqrt_obj->bufsize_v);
	memcpy((void*)ctpmqrt_obj->v, ctpmqrt_obj->A_tpqrt, ctpmqrt_obj->bufsize_v);*/
	
	
/*Compute ctpmqrt's  o/p */
    ctpmqrt_obj->inforef = ctpmqrt( ctpmqrt_obj->matrix_layout, ctpmqrt_obj->side, ctpmqrt_obj->trans,\
								ctpmqrt_obj->m, ctpmqrt_obj->n, ctpmqrt_obj->k, ctpmqrt_obj->l,\
								ctpmqrt_obj->nb, (const lapack_complex_float*)ctpmqrt_obj->vref,\
								ctpmqrt_obj->ldv, (const lapack_complex_float*)ctpmqrt_obj->tref,\
								ctpmqrt_obj->ldt, ctpmqrt_obj->Aref, ctpmqrt_obj->lda, ctpmqrt_obj->bref, ctpmqrt_obj->ldb);

    /* Compute libflame's Lapacke o/p  */
    ctpmqrt_obj->info = LAPACKE_ctpmqrt( ctpmqrt_obj->matrix_layout, ctpmqrt_obj->side, ctpmqrt_obj->trans,\
									ctpmqrt_obj->m, ctpmqrt_obj->n, ctpmqrt_obj->k, ctpmqrt_obj->l,\
									ctpmqrt_obj->nb, (const lapack_complex_float*)ctpmqrt_obj->v,\
									ctpmqrt_obj->ldv, (const lapack_complex_float*)ctpmqrt_obj->t,\ 
									ctpmqrt_obj->ldt, ctpmqrt_obj->A, ctpmqrt_obj->lda, ctpmqrt_obj->b, ctpmqrt_obj->ldb);
    if( ctpmqrt_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame\
        LAPACKE_ctpmqrt is wrong\n", ctpmqrt_obj->info );
    }
    if( ctpmqrt_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ctpmqrt is wrong\n", 
        ctpmqrt_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    ctpmqrt_obj->diff =  computeDiff_c( ctpmqrt_obj->bufsize_a, 
                ctpmqrt_obj->A, ctpmqrt_obj->Aref );

}

TEST_F(ctpmqrt_test, ctpmqrt1) {
    EXPECT_NEAR(0.0, ctpmqrt_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(ctpmqrt_test, ctpmqrt2) {
    EXPECT_NEAR(0.0, ctpmqrt_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(ctpmqrt_test, ctpmqrt3) {
    EXPECT_NEAR(0.0, ctpmqrt_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(ctpmqrt_test, ctpmqrt4) {
    EXPECT_NEAR(0.0, ctpmqrt_obj->diff, LAPACKE_EIG_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class tpmqrt_dcomplex_parameters{

   public:
	int bufsize_v;
	int bufsize_a, bufsize_atpqrt;
	int bufsize_b;
	int bufsize_t, bufsize_tpqrt;
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
	lapack_int nb, nb_tpqrt;
	lapack_int lda;
	lapack_complex_double* v, *A, *A_tpqrt, *t_tpqrt;
	char side;
	char trans;
	lapack_int ldb, ldv, lda_tpqrt, ldt, ldt_tpqrt;
	/*Output Parameter*/
	lapack_complex_double* t, *b;	
	lapack_complex_double *vref, *tref, *bref, *Aref, *A_tpqrtref, *t_tpqrtref;
	/*Return Values*/	
	lapack_int info, inforef;
	lapack_int info_tpqrt, inforef_tpqrt;

   public:
      tpmqrt_dcomplex_parameters (int matrix_layout, char side , char trans, lapack_int m, lapack_int n, lapack_int nb);
      ~tpmqrt_dcomplex_parameters ();

};

/* Constructor definition  dcomplex_common_parameters */
tpmqrt_dcomplex_parameters:: tpmqrt_dcomplex_parameters (int matrix_layout_i, char side_i , char trans_i,lapack_int m_i, lapack_int n_i, lapack_int nb_i)
{
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	side = side_i;
	trans = trans_i;
	nb = nb_i;
	l  = fla_min(m,n);
	k = m;
	
	if (trans == 'T')
		trans  = 'N' ;
	
	#if LAPACKE_TEST_VERBOSE
		printf(" \n tpmqrt dcomplex: matrix_layout = %d, side :%b, trans:%b, m:%d, n: %d \n", matrix_layout, side, trans, m, n);
	#endif
	if (side == 'L')
	{
		//k = m;
		//nb = k;
		if (matrix_layout == LAPACK_COL_MAJOR)
		{
			lda = k;
			ldv = m;
			ldt = nb;
			ldt_tpqrt = nb;
			ldb = m;
			bufsize_v = ldv*k;
			bufsize_t =  ldt*k;
			bufsize_a = lda*n;
			bufsize_b = ldb*n;
			bufsize_tpqrt = ldt_tpqrt*n;
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{		 
			lda= m;
			ldv = k;
			ldt = k;
			ldt_tpqrt = n;
			ldb = n;
			bufsize_v = ldv*m;
			bufsize_t = ldt*nb;
			bufsize_a = lda*k;
			bufsize_b = ldb*m;
			bufsize_tpqrt = ldt_tpqrt*nb;
		}else
		{
			EXPECT_TRUE(false) << "matrix_layout invalid";
		}
	}else if (side == 'R') {
		//k = n;
		//nb = k;
		if (matrix_layout == LAPACK_COL_MAJOR)
		{	lda= m;
			ldv = n;
			ldt = nb;
			ldb = m;
			ldt_tpqrt = nb;
			bufsize_v = ldv*n;
			bufsize_t =  ldt*k;
			bufsize_a = lda*k;;
			bufsize_b = ldb*n;
			bufsize_tpqrt = ldt_tpqrt*n;
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{	lda= k;
			ldv = k;
			ldt = k;
			ldb =n;
			ldt_tpqrt = n;
			bufsize_v = ldv*n;
			bufsize_t = ldt*nb;
			bufsize_a = lda*m;
			bufsize_b = ldb*m;
			bufsize_tpqrt = ldt_tpqrt*nb;
		}else
		{
			EXPECT_TRUE(false) << "matrix_layout invalid";
		}
	}else 
		EXPECT_TRUE(false) << "side is invalid";
	
	lda_tpqrt =n;
	bufsize_atpqrt =  lda_tpqrt*n;
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&v, &vref, bufsize_v);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&t, &tref, bufsize_t);	
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&b, &bref, bufsize_b);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A_tpqrt, &A_tpqrtref, bufsize_atpqrt);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&t_tpqrt, &t_tpqrtref, bufsize_tpqrt);

	
	if ((v==NULL) || (vref==NULL) ||
		(t==NULL) || (tref==NULL) ||
		(b==NULL) || (bref==NULL) ||
		(A==NULL) || (Aref==NULL) ||
		(A_tpqrt == NULL) || (A_tpqrtref == NULL) ||
		(t_tpqrt == NULL) || (t_tpqrtref == NULL))
	{
		EXPECT_FALSE( true) << "tpmqrt_dcomplex_parameters object: malloc error.";
		tpmqrt_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( v, vref, bufsize_v);
	lapacke_gtest_init_dcomplex_buffer_pair_rand( A, Aref, bufsize_a);
	lapacke_gtest_init_dcomplex_buffer_pair_rand( t, tref, bufsize_t);
	lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, bufsize_b);
	lapacke_gtest_init_dcomplex_buffer_pair_rand( A_tpqrt, A_tpqrtref, bufsize_atpqrt);


} /* end of Constructor  */

/* Destructor definition  'dcomplex_common_parameters' class  */
tpmqrt_dcomplex_parameters :: ~tpmqrt_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" tpmqrt_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   tpmqrt_free();

}
/*  Test fixture class definition */
class ztpmqrt_test  : public  ::testing::Test {
public:
   tpmqrt_dcomplex_parameters  *ztpmqrt_obj;
   void SetUp();
   void TearDown () { delete ztpmqrt_obj;}
};

void ztpmqrt_test::SetUp(){

	/* LAPACKE ztpqrt prototype */
    typedef int (*Fptr_NL_LAPACKE_ztpqrt) (int matrix_layout, lapack_int m, lapack_int n,\
                           lapack_int l, lapack_int nb,\
                           lapack_complex_double* a, lapack_int lda,\
                           lapack_complex_double* b, lapack_int ldb,\
                           lapack_complex_double* t, lapack_int ldt);
	
	 Fptr_NL_LAPACKE_ztpqrt ztpqrt;
	 /* LAPACKE ztpmqrt prototype */
    typedef int (*Fptr_NL_LAPACKE_ztpmqrt) (int matrix_layout, char side, char trans,\
                            lapack_int m, lapack_int n, lapack_int k,\
                            lapack_int l, lapack_int nb,\
                            const lapack_complex_double* v, lapack_int ldv,\
                            const lapack_complex_double* t, lapack_int ldt,\
                            lapack_complex_double* a, lapack_int lda,\
                            lapack_complex_double* b, lapack_int ldb );

    Fptr_NL_LAPACKE_ztpmqrt ztpmqrt;

    ztpmqrt_obj = new tpmqrt_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].side,
						   eig_paramslist[idx].trans,
						   eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
						   eig_paramslist[idx].nb);
						   

    idx = Circular_Increment_Index(idx);

    ztpmqrt_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ztpmqrt_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ztpmqrt_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ztpmqrt_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*ztpmqrt library call */
    ztpmqrt = (Fptr_NL_LAPACKE_ztpmqrt)dlsym(ztpmqrt_obj->hModule, "LAPACKE_ztpmqrt");
    ASSERT_TRUE(ztpmqrt != NULL) << "failed to get the Netlib LAPACKE_ztpmqrt symbol";

	/*ztpqrt library call*/
	ztpmqrt_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ztpmqrt_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ztpmqrt_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ztpmqrt_obj->lModule != NULL) << "Netlib lapacke handle NULL";

    ztpqrt = (Fptr_NL_LAPACKE_ztpqrt)dlsym(ztpmqrt_obj->lModule, "LAPACKE_ztpqrt");
    ASSERT_TRUE(ztpqrt != NULL) << "failed to get the Netlib LAPACKE_ztpqrt symbol";    
    

    ztpmqrt_obj->inforef_tpqrt = ztpqrt( ztpmqrt_obj->matrix_layout,ztpmqrt_obj->m, ztpmqrt_obj->n, ztpmqrt_obj->l, ztpmqrt_obj->nb,
								ztpmqrt_obj->A_tpqrtref, ztpmqrt_obj->lda_tpqrt, ztpmqrt_obj->bref, ztpmqrt_obj->ldb, ztpmqrt_obj->t_tpqrtref, ztpmqrt_obj->ldt_tpqrt);

    /* Compute libflame's Lapacke o/p  */
    ztpmqrt_obj->info_tpqrt = LAPACKE_ztpqrt( ztpmqrt_obj->matrix_layout,ztpmqrt_obj->m, ztpmqrt_obj->n, ztpmqrt_obj->l, ztpmqrt_obj->nb,
											ztpmqrt_obj->A_tpqrt, ztpmqrt_obj->lda_tpqrt, ztpmqrt_obj->b, ztpmqrt_obj->ldb, ztpmqrt_obj->t_tpqrt, ztpmqrt_obj->ldt_tpqrt);
										
										

    if( ztpmqrt_obj->info_tpqrt < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ztpqrt is wrong\n", ztpmqrt_obj->info_tpqrt );
    }
    if( ztpmqrt_obj->inforef_tpqrt < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ztpqrt is wrong\n", 
        ztpmqrt_obj->inforef_tpqrt );
    }
	
	/*Copy matrices from geqrt to gemqrt*
	memcpy((void*)ztpmqrt_obj->tref, ztpmqrt_obj->t_tpqrtref, ztpmqrt_obj->bufsize_t);
	memcpy((void*)ztpmqrt_obj->t, ztpmqrt_obj->t_tpqrt, ztpmqrt_obj->bufsize_t);
	
	memcpy((void*)ztpmqrt_obj->vref, ztpmqrt_obj->A_tpqrtref, ztpmqrt_obj->bufsize_v);
	memcpy((void*)ztpmqrt_obj->v, ztpmqrt_obj->A_tpqrt, ztpmqrt_obj->bufsize_v);*/
	
/*Compute ztpmqrt's  o/p */
    ztpmqrt_obj->inforef = ztpmqrt( ztpmqrt_obj->matrix_layout, ztpmqrt_obj->side, ztpmqrt_obj->trans,\
								ztpmqrt_obj->m, ztpmqrt_obj->n, ztpmqrt_obj->k, ztpmqrt_obj->l,\
								ztpmqrt_obj->nb, (const lapack_complex_double*)ztpmqrt_obj->vref,\
								ztpmqrt_obj->ldv, (const lapack_complex_double*)ztpmqrt_obj->tref,\
								ztpmqrt_obj->ldt, ztpmqrt_obj->Aref, ztpmqrt_obj->lda, ztpmqrt_obj->bref, ztpmqrt_obj->ldb);

    /* Compute libflame's Lapacke o/p  */
    ztpmqrt_obj->info = LAPACKE_ztpmqrt( ztpmqrt_obj->matrix_layout, ztpmqrt_obj->side, ztpmqrt_obj->trans,\
									ztpmqrt_obj->m, ztpmqrt_obj->n, ztpmqrt_obj->k, ztpmqrt_obj->l,\
									ztpmqrt_obj->nb, (const lapack_complex_double*)ztpmqrt_obj->v,\
									ztpmqrt_obj->ldv, (const lapack_complex_double*)ztpmqrt_obj->t,\ 
									ztpmqrt_obj->ldt, ztpmqrt_obj->A, ztpmqrt_obj->lda, ztpmqrt_obj->b, ztpmqrt_obj->ldb);
    if( ztpmqrt_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame\
        LAPACKE_ztpmqrt is wrong\n", ztpmqrt_obj->info );
    }
    if( ztpmqrt_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ztpmqrt is wrong\n", 
        ztpmqrt_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    ztpmqrt_obj->diff =  computeDiff_z( ztpmqrt_obj->bufsize_a, 
                ztpmqrt_obj->A, ztpmqrt_obj->Aref );

}

TEST_F(ztpmqrt_test, ztpmqrt1) {
    EXPECT_NEAR(0.0, ztpmqrt_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(ztpmqrt_test, ztpmqrt2) {
    EXPECT_NEAR(0.0, ztpmqrt_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(ztpmqrt_test, ztpmqrt3) {
    EXPECT_NEAR(0.0, ztpmqrt_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(ztpmqrt_test, ztpmqrt4) {
    EXPECT_NEAR(0.0, ztpmqrt_obj->diff, LAPACKE_EIG_THRESHOLD);
}
