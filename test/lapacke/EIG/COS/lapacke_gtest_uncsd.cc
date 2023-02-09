#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define uncsd_free() \
if (x11!=NULL)    free(x11); \
if (x11ref!=NULL) free(x11ref);\
if (x12 != NULL) free(x12); \
if (x12ref != NULL) free(x12ref); \
if (x21 != NULL) free(x21); \
if (x21ref!=NULL) free(x21ref);\
if (u1!=NULL)  free(u1);\
if (u1ref!=NULL) free(u1ref); \
if (u2!=NULL)  free(u2);\
if (u2ref!=NULL) free(u2ref); \
if (v1t!=NULL)  free(v1t);\
if (v1tref!=NULL) free(v1tref); \
if (v2t!=NULL)  free(v2t);\
if (v2tref!=NULL) free(v2tref); \
if (theta != NULL) free(theta); \
if (thetaref != NULL) free(thetaref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin scomplex_common_parameters  class definition */
class uncsd_scomplex_parameters{

   public:
	int bufsize_x11, bufsize_x12, bufsize_x21, bufsize_x22, bufsize_theta;
	int bufsize_u1, bufsize_u2, bufsize_v1t, bufsize_v2t;
	void *hModule, *dModule;
	float diff_x11, diff_x12, diff_x21, diff_x22, diff_theta;
	float diff_u1, diff_u2, diff_v1t, diff_v2t;
   /*input parameters */
	int matrix_layout;
	lapack_int p;
	lapack_int m;
	lapack_int q;
	lapack_int uncsd_threshold;
	lapack_complex_float* x11, *x12, *x21, *x22;
	char signs;
	char trans;
	char jobu1, jobu2, jobv1t, jobv2t;
	lapack_int ldx11, ldx12, ldx21, ldx22;
	lapack_int ldu1, ldu2, ldv1t, ldv2t;
	/*Output Parameter*/
	float* theta;
	float* thetaref;
	lapack_complex_float *u1, *u2, *v1t, *v2t;
	lapack_complex_float *u1ref, *u2ref, *v1tref, *v2tref;
	lapack_complex_float *x11ref, *x12ref, *x21ref, *x22ref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      uncsd_scomplex_parameters (int matrix_layout, char trans, char signs, lapack_int m, lapack_int p, lapack_int q, lapack_int uncsd_threshold);
      ~uncsd_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
uncsd_scomplex_parameters:: uncsd_scomplex_parameters (int matrix_layout_i, char trans_i, char signs_i, lapack_int m_i, lapack_int p_i, lapack_int q_i, lapack_int uncsd_threshold_i)
{
	matrix_layout = matrix_layout_i;
	m = m_i;
	p = p_i;
	q = q_i;
	m = m+p;
	signs = signs_i;
	trans = trans_i;
	uncsd_threshold = uncsd_threshold_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n uncsd scomplex: matrix_layout = %d, trans:%c, signs:%c, m: %d, p:%d, q:%d \n", matrix_layout, trans, signs, m, p, q);
	#endif
	
	/*check for signs*/
	if (signs == 'L')
		signs = 'O';
	
	/* assigning 'Y' all jobz variable*/
	jobu1 = jobu2 = jobv1t = jobv2t = 'Y';
	
	if (trans == 'C') // Trans is either 'N' or 'T'
		trans = 'N';
	
	/*sizes for x11*/
	ldx11 = p;
	bufsize_x11 = ldx11*q;
	/*sizes for x12*/
	ldx12 = p;
	bufsize_x12 = ldx12*(m-q);
	/*sizes for x21*/
	ldx21 = (m-p);
	bufsize_x21 = ldx21*q;
	/*sizes for x22*/
	ldx22 = m-p;
	bufsize_x22 = ldx22*(m-q);
	/*sizes for u1*/
	ldu1 = p;
	bufsize_u1 = ldu1*p;
	/*sizes for u2*/
	ldu2 = (m-p);
	bufsize_u2 = ldu2*(m-p);
	/*sizes for v1t*/
	ldv1t = q;
	bufsize_v1t = ldv1t*q;
	/*sizes for v2t*/
	ldv2t = (m-q);
	bufsize_v2t = ldu2*(m-q);
	
	/*Bufsize for theta */
	p < min(m-p,m-q)? bufsize_theta = p: bufsize_theta = min(m-p,m-q);
	bufsize_theta = min(bufsize_theta, q);	
		
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&x11, &x11ref, bufsize_x11);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&x12, &x12ref, bufsize_x12);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&x21, &x21ref, bufsize_x21);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&x22, &x22ref, bufsize_x22);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&u1, &u1ref, bufsize_u1);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&u2, &u2ref, bufsize_u2);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&v1t, &v1tref, bufsize_v1t);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&v2t, &v2tref, bufsize_v2t);
	lapacke_gtest_alloc_float_buffer_pair(&theta, &thetaref, bufsize_theta);

	if ((x11==NULL) || (x11ref==NULL) || \
		(x12==NULL) || (x12ref==NULL) || \
		(x21==NULL) || (x21ref==NULL) || \
		(x22==NULL) || (x22ref==NULL) || \
		(u1==NULL) || (u1ref==NULL) ||\
		(u2==NULL) || (u2ref==NULL) ||\
		(v1t==NULL) || (v1tref==NULL) ||\
		(v2t==NULL) || (v2tref==NULL) ||\
		(theta==NULL) || (thetaref==NULL)){
		EXPECT_FALSE( true) << "uncsd_float_parameters object: malloc error.";
		uncsd_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( x11, x11ref, bufsize_x11);
	lapacke_gtest_init_scomplex_buffer_pair_rand( x12, x12ref, bufsize_x12);
	lapacke_gtest_init_scomplex_buffer_pair_rand( x21, x21ref, bufsize_x21);
	lapacke_gtest_init_scomplex_buffer_pair_rand( x22, x22ref, bufsize_x22);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
uncsd_scomplex_parameters :: ~uncsd_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" uncsd_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   uncsd_free();

}
/*  Test fixture class definition */
class cuncsd_test  : public  ::testing::Test {
public:
   uncsd_scomplex_parameters  *cuncsd_obj;
   void SetUp();
   void TearDown () { delete cuncsd_obj;}
};

void cuncsd_test::SetUp(){

	/* LAPACKE cuncsd prototype */
    typedef int (*Fptr_NL_LAPACKE_cuncsd) (int matrix_layout, char jobu1, char jobu2,\
                           char jobv1t, char jobv2t, char trans, char signs,\
                           lapack_int m, lapack_int p, lapack_int q,\
                           lapack_complex_float* x11, lapack_int ldx11,\
                           lapack_complex_float* x12, lapack_int ldx12,\
                           lapack_complex_float* x21, lapack_int ldx21,\
                           lapack_complex_float* x22, lapack_int ldx22,\
                           float* theta, lapack_complex_float* u1,\
                           lapack_int ldu1, lapack_complex_float* u2,\
                           lapack_int ldu2, lapack_complex_float* v1t,\
                           lapack_int ldv1t, lapack_complex_float* v2t,\
                           lapack_int ldv2t );

    Fptr_NL_LAPACKE_cuncsd cuncsd;

    cuncsd_obj = new uncsd_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].trans,
                           eig_paramslist[idx].uplo,
						   eig_paramslist[idx].m,
						   eig_paramslist[idx].m,
						   eig_paramslist[idx].m,
						   eig_paramslist[idx].threshold_value);
						   

    idx = Circular_Increment_Index(idx);

    cuncsd_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cuncsd_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cuncsd_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cuncsd_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*cuncsd library call */
    cuncsd = (Fptr_NL_LAPACKE_cuncsd)dlsym(cuncsd_obj->hModule, "LAPACKE_cuncsd");
    ASSERT_TRUE(cuncsd != NULL) << "failed to get the Netlib LAPACKE_cuncsd symbol";
	
/*Compute cuncsd's  o/p */
    cuncsd_obj->inforef = cuncsd( cuncsd_obj->matrix_layout, cuncsd_obj->jobu1, cuncsd_obj->jobu2, cuncsd_obj->jobv1t, cuncsd_obj->jobv2t,
								cuncsd_obj->trans, cuncsd_obj->signs, cuncsd_obj->m, cuncsd_obj->p,cuncsd_obj->q, cuncsd_obj->x11ref,
								cuncsd_obj->ldx11, cuncsd_obj->x12ref, cuncsd_obj->ldx12, cuncsd_obj->x21ref, cuncsd_obj->ldx21, 
								cuncsd_obj->x22ref, cuncsd_obj->ldx22, cuncsd_obj->thetaref, cuncsd_obj->u1ref, cuncsd_obj->ldu1,
								cuncsd_obj->u2ref, cuncsd_obj->ldu2, cuncsd_obj->v1tref, cuncsd_obj->ldv1t, cuncsd_obj->v2tref, cuncsd_obj->ldv2t); 

    /* Compute libflame's Lapacke o/p  */
    cuncsd_obj->info = LAPACKE_cuncsd(cuncsd_obj->matrix_layout, cuncsd_obj->jobu1, cuncsd_obj->jobu2, cuncsd_obj->jobv1t, cuncsd_obj->jobv2t,
								cuncsd_obj->trans, cuncsd_obj->signs, cuncsd_obj->m, cuncsd_obj->p,cuncsd_obj->q, cuncsd_obj->x11,
								cuncsd_obj->ldx11, cuncsd_obj->x12, cuncsd_obj->ldx12, cuncsd_obj->x21, cuncsd_obj->ldx21, 
								cuncsd_obj->x22, cuncsd_obj->ldx22, cuncsd_obj->theta, cuncsd_obj->u1, cuncsd_obj->ldu1, cuncsd_obj->u2,
								cuncsd_obj->ldu2, cuncsd_obj->v1t, cuncsd_obj->ldv1t,  cuncsd_obj->v2t, cuncsd_obj->ldv2t);

    if( cuncsd_obj->info < 0 ) {
        printf( "\n warning: The i:%x22 th argument with libflame \
        LAPACKE_cuncsd is wrong\n", cuncsd_obj->info );
    }
    if( cuncsd_obj->inforef < 0 ) {
        printf( "The i:%x22 th argument with Netlib LAPACKE_cuncsd is wrong\n", 
        cuncsd_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cuncsd_obj->diff_x11 =  computeDiff_c( cuncsd_obj->bufsize_x11, cuncsd_obj->x11, cuncsd_obj->x11ref );
	cuncsd_obj->diff_x12 =  computeDiff_c( cuncsd_obj->bufsize_x12, cuncsd_obj->x12, cuncsd_obj->x12ref );
	cuncsd_obj->diff_x21 =  computeDiff_c( cuncsd_obj->bufsize_x21, cuncsd_obj->x21, cuncsd_obj->x21ref );
	cuncsd_obj->diff_x22 =  computeDiff_c( cuncsd_obj->bufsize_x22, cuncsd_obj->x21, cuncsd_obj->x22ref );
	cuncsd_obj->diff_u1 =  computeDiff_c( cuncsd_obj->bufsize_u1, cuncsd_obj->u1, cuncsd_obj->u1ref );
	cuncsd_obj->diff_u2 =  computeDiff_c( cuncsd_obj->bufsize_u2, cuncsd_obj->u2, cuncsd_obj->u2ref );
	cuncsd_obj->diff_v1t =  computeDiff_c( cuncsd_obj->bufsize_v1t, cuncsd_obj->v1t, cuncsd_obj->v1tref );
	cuncsd_obj->diff_v2t =  computeDiff_c( cuncsd_obj->bufsize_v2t, cuncsd_obj->v2t, cuncsd_obj->v2tref );
	cuncsd_obj->diff_theta =  computeDiff_s( cuncsd_obj->bufsize_theta, cuncsd_obj->theta, cuncsd_obj->thetaref );

}

TEST_F(cuncsd_test, cuncsd1) {
    EXPECT_NEAR(0.0, cuncsd_obj->diff_x11, cuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, cuncsd_obj->diff_x12, cuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, cuncsd_obj->diff_x21, cuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, cuncsd_obj->diff_x22, cuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, cuncsd_obj->diff_u1, cuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, cuncsd_obj->diff_u2, cuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, cuncsd_obj->diff_v1t, cuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, cuncsd_obj->diff_v2t, cuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, cuncsd_obj->diff_theta, cuncsd_obj->uncsd_threshold);
}

TEST_F(cuncsd_test, cuncsd2) {
    EXPECT_NEAR(0.0, cuncsd_obj->diff_x11, cuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, cuncsd_obj->diff_x12, cuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, cuncsd_obj->diff_x21, cuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, cuncsd_obj->diff_x22, cuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, cuncsd_obj->diff_u1, cuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, cuncsd_obj->diff_u2, cuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, cuncsd_obj->diff_v1t, cuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, cuncsd_obj->diff_v2t, cuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, cuncsd_obj->diff_theta, cuncsd_obj->uncsd_threshold);
}

TEST_F(cuncsd_test, cuncsd3) {
    EXPECT_NEAR(0.0, cuncsd_obj->diff_x11, cuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, cuncsd_obj->diff_x12, cuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, cuncsd_obj->diff_x21, cuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, cuncsd_obj->diff_x22, cuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, cuncsd_obj->diff_u1, cuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, cuncsd_obj->diff_u2, cuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, cuncsd_obj->diff_v1t, cuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, cuncsd_obj->diff_v2t, cuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, cuncsd_obj->diff_theta, cuncsd_obj->uncsd_threshold);

}

TEST_F(cuncsd_test, cuncsd4) {
    EXPECT_NEAR(0.0, cuncsd_obj->diff_x11, cuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, cuncsd_obj->diff_x12, cuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, cuncsd_obj->diff_x21, cuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, cuncsd_obj->diff_x22, cuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, cuncsd_obj->diff_u1, cuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, cuncsd_obj->diff_u2, cuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, cuncsd_obj->diff_v1t, cuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, cuncsd_obj->diff_v2t, cuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, cuncsd_obj->diff_theta, cuncsd_obj->uncsd_threshold);
}

/* Begin dcomplex_common_parameters  class definition */
class uncsd_dcomplex_parameters{

   public:
	int bufsize_x11, bufsize_x12, bufsize_x21, bufsize_x22, bufsize_theta;
	int bufsize_u1, bufsize_u2, bufsize_v1t, bufsize_v2t;
	void *hModule, *dModule;
	double diff_x11, diff_x12, diff_x21, diff_x22, diff_theta;
	double diff_u1, diff_u2, diff_v1t, diff_v2t;
   /*input parameters */
	int matrix_layout;
	lapack_int p;
	lapack_int m;
	lapack_int q;
	lapack_int uncsd_threshold;
	lapack_complex_double* x11, *x12, *x21, *x22;
	char signs;
	char trans;
	char jobu1, jobu2, jobv1t, jobv2t;
	lapack_int ldx11, ldx12, ldx21, ldx22;
	lapack_int ldu1, ldu2, ldv1t, ldv2t;
	/*Output Parameter*/
	double* theta;
	double* thetaref;
	lapack_complex_double *u1, *u2, *v1t, *v2t;
	lapack_complex_double *u1ref, *u2ref, *v1tref, *v2tref;
	lapack_complex_double *x11ref, *x12ref, *x21ref, *x22ref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      uncsd_dcomplex_parameters (int matrix_layout, char trans, char signs, lapack_int m, lapack_int p, lapack_int q, lapack_int uncsd_threshold);
      ~uncsd_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
uncsd_dcomplex_parameters:: uncsd_dcomplex_parameters (int matrix_layout_i, char trans_i, char signs_i, lapack_int m_i, lapack_int p_i, lapack_int q_i, lapack_int uncsd_threshold_i)
{
	matrix_layout = matrix_layout_i;
	m = m_i;
	p = p_i;
	q = q_i;
	m = m+p;
	signs = signs_i;
	trans = trans_i;
	uncsd_threshold = uncsd_threshold_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n uncsd dcomplex: matrix_layout = %d, trans:%c, signs:%c, m: %d, p:%d, q:%d \n", matrix_layout, trans, signs, m, p, q);
	#endif
	
	/*check for signs*/
	if (signs == 'L')
		signs = 'O';
	
	/* assigning 'Y' all jobz variable*/
	jobu1 = jobu2 = jobv1t = jobv2t = 'Y';
	
	if (trans == 'C') // Trans is either 'N' or 'T'
		trans = 'N';
	
	
	/*sizes for x11*/
	ldx11 = p;
	bufsize_x11 = ldx11*q;
	/*sizes for x12*/
	ldx12 = p;
	bufsize_x12 = ldx12*(m-q);
	/*sizes for x21*/
	ldx21 = (m-p);
	bufsize_x21 = ldx21*q;
	/*sizes for x22*/
	ldx22 = m-p;
	bufsize_x22 = ldx22*(m-q);
	/*sizes for u1*/
	ldu1 = p;
	bufsize_u1 = ldu1*p;
	/*sizes for u2*/
	ldu2 = (m-p);
	bufsize_u2 = ldu2*(m-p);
	/*sizes for v1t*/
	ldv1t = q;
	bufsize_v1t = ldv1t*q;
	/*sizes for v2t*/
	ldv2t = (m-q);
	bufsize_v2t = ldu2*(m-q);
	
	/*Bufsize for theta */
	p < min(m-p,m-q)? bufsize_theta = p: bufsize_theta = min(m-p,m-q);
	bufsize_theta = min(bufsize_theta, q);	
	
		
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&x11, &x11ref, bufsize_x11);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&x12, &x12ref, bufsize_x12);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&x21, &x21ref, bufsize_x21);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&x22, &x22ref, bufsize_x22);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&u1, &u1ref, bufsize_u1);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&u2, &u2ref, bufsize_u2);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&v1t, &v1tref, bufsize_v1t);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&v2t, &v2tref, bufsize_v2t);
	lapacke_gtest_alloc_double_buffer_pair(&theta, &thetaref, bufsize_theta);

	if ((x11==NULL) || (x11ref==NULL) || \
		(x12==NULL) || (x12ref==NULL) || \
		(x21==NULL) || (x21ref==NULL) || \
		(x22==NULL) || (x22ref==NULL) || \
		(u1==NULL) || (u1ref==NULL) ||\
		(u2==NULL) || (u2ref==NULL) ||\
		(v1t==NULL) || (v1tref==NULL) ||\
		(v2t==NULL) || (v2tref==NULL) ||\
		(theta==NULL) || (thetaref==NULL)){
		EXPECT_FALSE( true) << "uncsd_double_parameters object: malloc error.";
		uncsd_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( x11, x11ref, bufsize_x11);
	lapacke_gtest_init_dcomplex_buffer_pair_rand( x12, x12ref, bufsize_x12);
	lapacke_gtest_init_dcomplex_buffer_pair_rand( x21, x21ref, bufsize_x21);
	lapacke_gtest_init_dcomplex_buffer_pair_rand( x22, x22ref, bufsize_x22);

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
uncsd_dcomplex_parameters :: ~uncsd_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" uncsd_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   uncsd_free();

}
/*  Test fixture class definition */
class zuncsd_test  : public  ::testing::Test {
public:
   uncsd_dcomplex_parameters  *zuncsd_obj;
   void SetUp();
   void TearDown () { delete zuncsd_obj;}
};

void zuncsd_test::SetUp(){

	/* LAPACKE zuncsd prototype */
    typedef int (*Fptr_NL_LAPACKE_zuncsd) (int matrix_layout, char jobu1, char jobu2,\
                           char jobv1t, char jobv2t, char trans, char signs,\
                           lapack_int m, lapack_int p, lapack_int q,\
                           lapack_complex_double* x11, lapack_int ldx11,\
                           lapack_complex_double* x12, lapack_int ldx12,\
                           lapack_complex_double* x21, lapack_int ldx21,\
                           lapack_complex_double* x22, lapack_int ldx22,\
                           double* theta, lapack_complex_double* u1,\
                           lapack_int ldu1, lapack_complex_double* u2,\
                           lapack_int ldu2, lapack_complex_double* v1t,\
                           lapack_int ldv1t, lapack_complex_double* v2t,\
                           lapack_int ldv2t );

    Fptr_NL_LAPACKE_zuncsd zuncsd;

    zuncsd_obj = new uncsd_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].trans,
                           eig_paramslist[idx].uplo,
						   eig_paramslist[idx].m,
						   eig_paramslist[idx].m,
						   eig_paramslist[idx].m,
						   eig_paramslist[idx].threshold_value);
						   

    idx = Circular_Increment_Index(idx);

    zuncsd_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zuncsd_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zuncsd_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zuncsd_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*zuncsd library call */
    zuncsd = (Fptr_NL_LAPACKE_zuncsd)dlsym(zuncsd_obj->hModule, "LAPACKE_zuncsd");
    ASSERT_TRUE(zuncsd != NULL) << "failed to get the Netlib LAPACKE_zuncsd symbol";
	
/*Compute zuncsd's  o/p */
    zuncsd_obj->inforef = zuncsd( zuncsd_obj->matrix_layout, zuncsd_obj->jobu1, zuncsd_obj->jobu2, zuncsd_obj->jobv1t, zuncsd_obj->jobv2t,
								zuncsd_obj->trans, zuncsd_obj->signs, zuncsd_obj->m, zuncsd_obj->p,zuncsd_obj->q, zuncsd_obj->x11ref,
								zuncsd_obj->ldx11, zuncsd_obj->x12ref, zuncsd_obj->ldx12, zuncsd_obj->x21ref, zuncsd_obj->ldx21, 
								zuncsd_obj->x22ref, zuncsd_obj->ldx22, zuncsd_obj->thetaref, zuncsd_obj->u1ref, zuncsd_obj->ldu1,
								zuncsd_obj->u2ref, zuncsd_obj->ldu2, zuncsd_obj->v1tref, zuncsd_obj->ldv1t, zuncsd_obj->v2tref, zuncsd_obj->ldv2t); 

    /* Compute libflame's Lapacke o/p  */
    zuncsd_obj->info = LAPACKE_zuncsd(zuncsd_obj->matrix_layout, zuncsd_obj->jobu1, zuncsd_obj->jobu2, zuncsd_obj->jobv1t, zuncsd_obj->jobv2t,
								zuncsd_obj->trans, zuncsd_obj->signs, zuncsd_obj->m, zuncsd_obj->p,zuncsd_obj->q, zuncsd_obj->x11,
								zuncsd_obj->ldx11, zuncsd_obj->x12, zuncsd_obj->ldx12, zuncsd_obj->x21, zuncsd_obj->ldx21, 
								zuncsd_obj->x22, zuncsd_obj->ldx22, zuncsd_obj->theta, zuncsd_obj->u1, zuncsd_obj->ldu1, zuncsd_obj->u2,
								zuncsd_obj->ldu2, zuncsd_obj->v1t, zuncsd_obj->ldv1t,  zuncsd_obj->v2t, zuncsd_obj->ldv2t);

    if( zuncsd_obj->info < 0 ) {
        printf( "\n warning: The i:%x22 th argument with libflame \
        LAPACKE_zuncsd is wrong\n", zuncsd_obj->info );
    }
    if( zuncsd_obj->inforef < 0 ) {
        printf( "The i:%x22 th argument with Netlib LAPACKE_zuncsd is wrong\n", 
        zuncsd_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zuncsd_obj->diff_x11 =  computeDiff_z( zuncsd_obj->bufsize_x11, zuncsd_obj->x11, zuncsd_obj->x11ref );
	zuncsd_obj->diff_x12 =  computeDiff_z( zuncsd_obj->bufsize_x12, zuncsd_obj->x12, zuncsd_obj->x12ref );
	zuncsd_obj->diff_x21 =  computeDiff_z( zuncsd_obj->bufsize_x21, zuncsd_obj->x21, zuncsd_obj->x21ref );
	zuncsd_obj->diff_x22 =  computeDiff_z( zuncsd_obj->bufsize_x22, zuncsd_obj->x21, zuncsd_obj->x22ref );
	zuncsd_obj->diff_u1 =  computeDiff_z( zuncsd_obj->bufsize_u1, zuncsd_obj->u1, zuncsd_obj->u1ref );
	zuncsd_obj->diff_u2 =  computeDiff_z( zuncsd_obj->bufsize_u2, zuncsd_obj->u2, zuncsd_obj->u2ref );
	zuncsd_obj->diff_v1t =  computeDiff_z( zuncsd_obj->bufsize_v1t, zuncsd_obj->v1t, zuncsd_obj->v1tref );
	zuncsd_obj->diff_v2t =  computeDiff_z( zuncsd_obj->bufsize_v2t, zuncsd_obj->v2t, zuncsd_obj->v2tref );
	zuncsd_obj->diff_theta =  computeDiff_d( zuncsd_obj->bufsize_theta, zuncsd_obj->theta, zuncsd_obj->thetaref );

}

TEST_F(zuncsd_test, zuncsd1) {
    EXPECT_NEAR(0.0, zuncsd_obj->diff_x11, zuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, zuncsd_obj->diff_x12, zuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, zuncsd_obj->diff_x21, zuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, zuncsd_obj->diff_x22, zuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, zuncsd_obj->diff_u1, zuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, zuncsd_obj->diff_u2, zuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, zuncsd_obj->diff_v1t, zuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, zuncsd_obj->diff_v2t, zuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, zuncsd_obj->diff_theta, zuncsd_obj->uncsd_threshold);
}

TEST_F(zuncsd_test, zuncsd2) {
    EXPECT_NEAR(0.0, zuncsd_obj->diff_x11, zuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, zuncsd_obj->diff_x12, zuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, zuncsd_obj->diff_x21, zuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, zuncsd_obj->diff_x22, zuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, zuncsd_obj->diff_u1, zuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, zuncsd_obj->diff_u2, zuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, zuncsd_obj->diff_v1t, zuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, zuncsd_obj->diff_v2t, zuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, zuncsd_obj->diff_theta, zuncsd_obj->uncsd_threshold);
}

TEST_F(zuncsd_test, zuncsd3) {
    EXPECT_NEAR(0.0, zuncsd_obj->diff_x11, zuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, zuncsd_obj->diff_x12, zuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, zuncsd_obj->diff_x21, zuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, zuncsd_obj->diff_x22, zuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, zuncsd_obj->diff_u1, zuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, zuncsd_obj->diff_u2, zuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, zuncsd_obj->diff_v1t, zuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, zuncsd_obj->diff_v2t, zuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, zuncsd_obj->diff_theta, zuncsd_obj->uncsd_threshold);

}

TEST_F(zuncsd_test, zuncsd4) {
    EXPECT_NEAR(0.0, zuncsd_obj->diff_x11, zuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, zuncsd_obj->diff_x12, zuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, zuncsd_obj->diff_x21, zuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, zuncsd_obj->diff_x22, zuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, zuncsd_obj->diff_u1, zuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, zuncsd_obj->diff_u2, zuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, zuncsd_obj->diff_v1t, zuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, zuncsd_obj->diff_v2t, zuncsd_obj->uncsd_threshold);
	EXPECT_NEAR(0.0, zuncsd_obj->diff_theta, zuncsd_obj->uncsd_threshold);
}