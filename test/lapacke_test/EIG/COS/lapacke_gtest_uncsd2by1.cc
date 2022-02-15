#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"

#define uncsd2by1_free() \
if (x11!=NULL)    free(x11); \
if (x11ref!=NULL) free(x11ref);\
if (x21 != NULL) free(x21); \
if (x21ref!=NULL) free(x21ref);\
if (u1!=NULL)  free(u1);\
if (u1ref!=NULL) free(u1ref); \
if (u2!=NULL)  free(u2);\
if (u2ref!=NULL) free(u2ref); \
if (v1t!=NULL)  free(v1t);\
if (v1tref!=NULL) free(v1tref); \
if (theta != NULL) free(theta); \
if (thetaref != NULL) free(thetaref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin scomplex_common_parameters  class definition */
class uncsd2by1_scomplex_parameters{

   public:
	int bufsize_x11, bufsize_x21, bufsize_theta;
	int bufsize_u1, bufsize_u2, bufsize_v1t;
	void *hModule, *dModule;
	float diff_x11, diff_x21, diff_theta;
	float diff_u1, diff_u2, diff_v1t;
   /*input parameters */
	int matrix_layout;
	lapack_int p;
	lapack_int m;
	lapack_int q;
	lapack_int uncsd2by1_threshold;
	lapack_complex_float* x11, *x21;
	char jobu1, jobu2, jobv1t;
	lapack_int ldx11,ldx21;
	lapack_int ldu1, ldu2, ldv1t;
	/*Output Parameter*/
	float* theta;
	float* thetaref;
	lapack_complex_float *u1, *u2, *v1t;
	lapack_complex_float *u1ref, *u2ref, *v1tref;
	lapack_complex_float *x11ref, *x21ref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      uncsd2by1_scomplex_parameters (int matrix_layout, lapack_int m, lapack_int p, lapack_int q, lapack_int uncsd2by1_threshold);
      ~uncsd2by1_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
uncsd2by1_scomplex_parameters:: uncsd2by1_scomplex_parameters (int matrix_layout_i, lapack_int m_i, lapack_int p_i, lapack_int q_i, lapack_int uncsd2by1_threshold_i)
{
	matrix_layout = matrix_layout_i;
	m = m_i;
	p = p_i;
	q = q_i;
	m = m+p;
	uncsd2by1_threshold = uncsd2by1_threshold_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n uncsd2by1 scomplex: matrix_layout = %d, m: %d, p:%d, q:%d \n", matrix_layout, m, p, q);
	#endif

	/* assigning 'Y' all jobz variable*/
	jobu1 = jobu2 = jobv1t  = 'Y';

	/*sizes for x11*/
	ldx11 = p;
	bufsize_x11 = ldx11*q;
	/*sizes for x21*/
	ldx21 = (m-p);
	bufsize_x21 = ldx21*q;
	/*sizes for u1*/
	ldu1 = p;
	bufsize_u1 = ldu1*p;
	/*sizes for u2*/
	ldu2 = (m-p);
	bufsize_u2 = ldu2*(m-p);
	/*sizes for v1t*/
	ldv1t = q;
	bufsize_v1t = ldv1t*q;
	
	/*Bufsize for theta */
	p < min(m-p,m-q)? bufsize_theta = p: bufsize_theta = min(m-p,m-q);
	bufsize_theta = min(bufsize_theta, q);	

		
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&x11, &x11ref, bufsize_x11);	
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&x21, &x21ref, bufsize_x21);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&u1, &u1ref, bufsize_u1);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&u2, &u2ref, bufsize_u2);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&v1t, &v1tref, bufsize_v1t);
	lapacke_gtest_alloc_float_buffer_pair(&theta, &thetaref, bufsize_theta);

	if ((x11==NULL) || (x11ref==NULL) || \
		(x21==NULL) || (x21ref==NULL) || \
		(u1==NULL) || (u1ref==NULL) ||\
		(u2==NULL) || (u2ref==NULL) ||\
		(v1t==NULL) || (v1tref==NULL) ||\
		(theta==NULL) || (thetaref==NULL)){
		EXPECT_FALSE( true) << "uncsd2by1_float_parameters object: malloc error.";
		uncsd2by1_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( x11, x11ref, bufsize_x11);	
	lapacke_gtest_init_scomplex_buffer_pair_rand( x21, x21ref, bufsize_x21);
	

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
uncsd2by1_scomplex_parameters :: ~uncsd2by1_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" uncsd2by1_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   uncsd2by1_free();

}
/*  Test fixture class definition */
class cuncsd2by1_test  : public  ::testing::Test {
public:
   uncsd2by1_scomplex_parameters  *cuncsd2by1_obj;
   void SetUp();
   void TearDown () { delete cuncsd2by1_obj;}
};

void cuncsd2by1_test::SetUp(){

	/* LAPACKE cuncsd2by1 prototype */
    typedef int (*Fptr_NL_LAPACKE_cuncsd2by1) (int matrix_layout, char jobu1, char jobu2,\
                           char jobv1t, lapack_int m, lapack_int p, lapack_int q,\
                           lapack_complex_float* x11, lapack_int ldx11,\
                           lapack_complex_float* x21, lapack_int ldx21,\
                           float* theta, lapack_complex_float* u1,\
                           lapack_int ldu1, lapack_complex_float* u2,\
                           lapack_int ldu2, lapack_complex_float* v1t, lapack_int ldv1t  );

    Fptr_NL_LAPACKE_cuncsd2by1 cuncsd2by1;

    cuncsd2by1_obj = new uncsd2by1_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
						   eig_paramslist[idx].m,
						   eig_paramslist[idx].m,
						   eig_paramslist[idx].m, 
						   eig_paramslist[idx].threshold_value);
						   

    idx = Circular_Increment_Index(idx);

    cuncsd2by1_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cuncsd2by1_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cuncsd2by1_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cuncsd2by1_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*cuncsd2by1 library call */
    cuncsd2by1 = (Fptr_NL_LAPACKE_cuncsd2by1)dlsym(cuncsd2by1_obj->hModule, "LAPACKE_cuncsd2by1");
    ASSERT_TRUE(cuncsd2by1 != NULL) << "failed to get the Netlib LAPACKE_cuncsd2by1 symbol";
	
/*Compute cuncsd2by1's  o/p */
    cuncsd2by1_obj->inforef = cuncsd2by1( cuncsd2by1_obj->matrix_layout, cuncsd2by1_obj->jobu1, cuncsd2by1_obj->jobu2, cuncsd2by1_obj->jobv1t,
								cuncsd2by1_obj->m, cuncsd2by1_obj->p,cuncsd2by1_obj->q, cuncsd2by1_obj->x11ref,
								cuncsd2by1_obj->ldx11, cuncsd2by1_obj->x21ref, cuncsd2by1_obj->ldx21, cuncsd2by1_obj->thetaref, cuncsd2by1_obj->u1ref, cuncsd2by1_obj->ldu1,
								cuncsd2by1_obj->u2ref, cuncsd2by1_obj->ldu2, cuncsd2by1_obj->v1tref, cuncsd2by1_obj->ldv1t);

    /* Compute libflame's Lapacke o/p  */
    cuncsd2by1_obj->info = LAPACKE_cuncsd2by1(cuncsd2by1_obj->matrix_layout, cuncsd2by1_obj->jobu1, cuncsd2by1_obj->jobu2, cuncsd2by1_obj->jobv1t,
								cuncsd2by1_obj->m, cuncsd2by1_obj->p,cuncsd2by1_obj->q, cuncsd2by1_obj->x11,
								cuncsd2by1_obj->ldx11, cuncsd2by1_obj->x21, cuncsd2by1_obj->ldx21, cuncsd2by1_obj->theta, cuncsd2by1_obj->u1, 
								cuncsd2by1_obj->ldu1, cuncsd2by1_obj->u2, cuncsd2by1_obj->ldu2, cuncsd2by1_obj->v1t, cuncsd2by1_obj->ldv1t);

    if( cuncsd2by1_obj->info < 0 ) {
        printf( "\n warning: The i:%x22 th argument with libflame \
        LAPACKE_cuncsd2by1 is wrong\n", cuncsd2by1_obj->info );
    }
    if( cuncsd2by1_obj->inforef < 0 ) {
        printf( "The i:%x22 th argument with Netlib LAPACKE_cuncsd2by1 is wrong\n", 
        cuncsd2by1_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cuncsd2by1_obj->diff_x11 =  computeDiff_c( cuncsd2by1_obj->bufsize_x11, cuncsd2by1_obj->x11, cuncsd2by1_obj->x11ref );
	cuncsd2by1_obj->diff_x21 =  computeDiff_c( cuncsd2by1_obj->bufsize_x21, cuncsd2by1_obj->x21, cuncsd2by1_obj->x21ref );
	cuncsd2by1_obj->diff_u1 =  computeDiff_c( cuncsd2by1_obj->bufsize_u1, cuncsd2by1_obj->u1, cuncsd2by1_obj->u1ref );
	cuncsd2by1_obj->diff_u2 =  computeDiff_c( cuncsd2by1_obj->bufsize_u2, cuncsd2by1_obj->u2, cuncsd2by1_obj->u2ref );
	cuncsd2by1_obj->diff_v1t =  computeDiff_c( cuncsd2by1_obj->bufsize_v1t, cuncsd2by1_obj->v1t, cuncsd2by1_obj->v1tref );
	cuncsd2by1_obj->diff_theta =  computeDiff_s( cuncsd2by1_obj->bufsize_theta, cuncsd2by1_obj->theta, cuncsd2by1_obj->thetaref );

}

TEST_F(cuncsd2by1_test, cuncsd2by11) {
    EXPECT_NEAR(0.0, cuncsd2by1_obj->diff_x11, cuncsd2by1_obj->uncsd2by1_threshold);
	EXPECT_NEAR(0.0, cuncsd2by1_obj->diff_x21, cuncsd2by1_obj->uncsd2by1_threshold);
	EXPECT_NEAR(0.0, cuncsd2by1_obj->diff_u1, cuncsd2by1_obj->uncsd2by1_threshold);
	EXPECT_NEAR(0.0, cuncsd2by1_obj->diff_u2, cuncsd2by1_obj->uncsd2by1_threshold);
	EXPECT_NEAR(0.0, cuncsd2by1_obj->diff_v1t, cuncsd2by1_obj->uncsd2by1_threshold);
	EXPECT_NEAR(0.0, cuncsd2by1_obj->diff_theta, cuncsd2by1_obj->uncsd2by1_threshold);
}

TEST_F(cuncsd2by1_test, cuncsd2by12) {
    EXPECT_NEAR(0.0, cuncsd2by1_obj->diff_x11, cuncsd2by1_obj->uncsd2by1_threshold);
	EXPECT_NEAR(0.0, cuncsd2by1_obj->diff_x21, cuncsd2by1_obj->uncsd2by1_threshold);
	EXPECT_NEAR(0.0, cuncsd2by1_obj->diff_u1, cuncsd2by1_obj->uncsd2by1_threshold);
	EXPECT_NEAR(0.0, cuncsd2by1_obj->diff_u2, cuncsd2by1_obj->uncsd2by1_threshold);
	EXPECT_NEAR(0.0, cuncsd2by1_obj->diff_v1t, cuncsd2by1_obj->uncsd2by1_threshold);
	EXPECT_NEAR(0.0, cuncsd2by1_obj->diff_theta, cuncsd2by1_obj->uncsd2by1_threshold);
}

TEST_F(cuncsd2by1_test, cuncsd2by13) {
    EXPECT_NEAR(0.0, cuncsd2by1_obj->diff_x11, cuncsd2by1_obj->uncsd2by1_threshold);
	EXPECT_NEAR(0.0, cuncsd2by1_obj->diff_x21, cuncsd2by1_obj->uncsd2by1_threshold);
	EXPECT_NEAR(0.0, cuncsd2by1_obj->diff_u1, cuncsd2by1_obj->uncsd2by1_threshold);
	EXPECT_NEAR(0.0, cuncsd2by1_obj->diff_u2, cuncsd2by1_obj->uncsd2by1_threshold);
	EXPECT_NEAR(0.0, cuncsd2by1_obj->diff_v1t, cuncsd2by1_obj->uncsd2by1_threshold);
	EXPECT_NEAR(0.0, cuncsd2by1_obj->diff_theta, cuncsd2by1_obj->uncsd2by1_threshold);

}

TEST_F(cuncsd2by1_test, cuncsd2by14) {
    EXPECT_NEAR(0.0, cuncsd2by1_obj->diff_x11, cuncsd2by1_obj->uncsd2by1_threshold);
	EXPECT_NEAR(0.0, cuncsd2by1_obj->diff_x21, cuncsd2by1_obj->uncsd2by1_threshold);
	EXPECT_NEAR(0.0, cuncsd2by1_obj->diff_u1, cuncsd2by1_obj->uncsd2by1_threshold);
	EXPECT_NEAR(0.0, cuncsd2by1_obj->diff_u2, cuncsd2by1_obj->uncsd2by1_threshold);
	EXPECT_NEAR(0.0, cuncsd2by1_obj->diff_v1t, cuncsd2by1_obj->uncsd2by1_threshold);
	EXPECT_NEAR(0.0, cuncsd2by1_obj->diff_theta, cuncsd2by1_obj->uncsd2by1_threshold);
}

/* Begin dcomplex_common_parameters  class definition */
class uncsd2by1_dcomplex_parameters{

   public:
	int bufsize_x11, bufsize_x21, bufsize_theta;
	int bufsize_u1, bufsize_u2, bufsize_v1t;
	void *hModule, *dModule;
	double diff_x11, diff_x21, diff_theta;
	double diff_u1, diff_u2, diff_v1t;
   /*input parameters */
	int matrix_layout;
	lapack_int p;
	lapack_int m;
	lapack_int q;
	lapack_int uncsd2by1_threshold;
	lapack_complex_double* x11, *x21;
	char jobu1, jobu2, jobv1t;
	lapack_int ldx11,ldx21;
	lapack_int ldu1, ldu2, ldv1t;
	/*Output Parameter*/
	double* theta;
	double* thetaref;
	lapack_complex_double *u1, *u2, *v1t;
	lapack_complex_double *u1ref, *u2ref, *v1tref;
	lapack_complex_double *x11ref, *x21ref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      uncsd2by1_dcomplex_parameters (int matrix_layout, lapack_int m, lapack_int p, lapack_int q, lapack_int uncsd2by1_threshold);
      ~uncsd2by1_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
uncsd2by1_dcomplex_parameters:: uncsd2by1_dcomplex_parameters (int matrix_layout_i, lapack_int m_i, lapack_int p_i, lapack_int q_i, lapack_int uncsd2by1_threshold_i)
{
	matrix_layout = matrix_layout_i;
	m = m_i;
	p = p_i;
	q = q_i;
	m = m+p;
	uncsd2by1_threshold = uncsd2by1_threshold_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n uncsd2by1 dcomplex: matrix_layout = %d, m: %d, p:%d, q:%d \n", matrix_layout, m, p, q);
	#endif

	/* assigning 'Y' all jobz variable*/
	jobu1 = jobu2 = jobv1t  = 'Y';
	
	/*sizes for x11*/
	ldx11 = p;
	bufsize_x11 = ldx11*q;
	/*sizes for x21*/
	ldx21 = (m-p);
	bufsize_x21 = ldx21*q;
	/*sizes for u1*/
	ldu1 = p;
	bufsize_u1 = ldu1*p;
	/*sizes for u2*/
	ldu2 = (m-p);
	bufsize_u2 = ldu2*(m-p);
	/*sizes for v1t*/
	ldv1t = q;
	bufsize_v1t = ldv1t*q;
	
	/*Bufsize for theta */
	p < min(m-p,m-q)? bufsize_theta = p: bufsize_theta = min(m-p,m-q);
	bufsize_theta = min(bufsize_theta, q);	
	
		
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&x11, &x11ref, bufsize_x11);	
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&x21, &x21ref, bufsize_x21);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&u1, &u1ref, bufsize_u1);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&u2, &u2ref, bufsize_u2);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&v1t, &v1tref, bufsize_v1t);
	lapacke_gtest_alloc_double_buffer_pair(&theta, &thetaref, bufsize_theta);

	if ((x11==NULL) || (x11ref==NULL) || \
		(x21==NULL) || (x21ref==NULL) || \
		(u1==NULL) || (u1ref==NULL) ||\
		(u2==NULL) || (u2ref==NULL) ||\
		(v1t==NULL) || (v1tref==NULL) ||\
		(theta==NULL) || (thetaref==NULL)){
		EXPECT_FALSE( true) << "uncsd2by1_double_parameters object: malloc error.";
		uncsd2by1_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( x11, x11ref, bufsize_x11);	
	lapacke_gtest_init_dcomplex_buffer_pair_rand( x21, x21ref, bufsize_x21);
	

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
uncsd2by1_dcomplex_parameters :: ~uncsd2by1_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" uncsd2by1_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   uncsd2by1_free();

}
/*  Test fixture class definition */
class zuncsd2by1_test  : public  ::testing::Test {
public:
   uncsd2by1_dcomplex_parameters  *zuncsd2by1_obj;
   void SetUp();
   void TearDown () { delete zuncsd2by1_obj;}
};

void zuncsd2by1_test::SetUp(){

	/* LAPACKE zuncsd2by1 prototype */
    typedef int (*Fptr_NL_LAPACKE_zuncsd2by1) (int matrix_layout, char jobu1, char jobu2,\
                           char jobv1t, lapack_int m, lapack_int p, lapack_int q,\
                           lapack_complex_double* x11, lapack_int ldx11,\
                           lapack_complex_double* x21, lapack_int ldx21,\
                           double* theta, lapack_complex_double* u1,\
                           lapack_int ldu1, lapack_complex_double* u2,\
                           lapack_int ldu2, lapack_complex_double* v1t, lapack_int ldv1t  );

    Fptr_NL_LAPACKE_zuncsd2by1 zuncsd2by1;

    zuncsd2by1_obj = new uncsd2by1_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
						   eig_paramslist[idx].m,
						   eig_paramslist[idx].m,
						   eig_paramslist[idx].m, 
						   eig_paramslist[idx].threshold_value);
						   

    idx = Circular_Increment_Index(idx);

    zuncsd2by1_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zuncsd2by1_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zuncsd2by1_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zuncsd2by1_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*zuncsd2by1 library call */
    zuncsd2by1 = (Fptr_NL_LAPACKE_zuncsd2by1)dlsym(zuncsd2by1_obj->hModule, "LAPACKE_zuncsd2by1");
    ASSERT_TRUE(zuncsd2by1 != NULL) << "failed to get the Netlib LAPACKE_zuncsd2by1 symbol";
	
/*Compute zuncsd2by1's  o/p */
    zuncsd2by1_obj->inforef = zuncsd2by1( zuncsd2by1_obj->matrix_layout, zuncsd2by1_obj->jobu1, zuncsd2by1_obj->jobu2, zuncsd2by1_obj->jobv1t,
								zuncsd2by1_obj->m, zuncsd2by1_obj->p,zuncsd2by1_obj->q, zuncsd2by1_obj->x11ref,
								zuncsd2by1_obj->ldx11, zuncsd2by1_obj->x21ref, zuncsd2by1_obj->ldx21, zuncsd2by1_obj->thetaref, zuncsd2by1_obj->u1ref, zuncsd2by1_obj->ldu1,
								zuncsd2by1_obj->u2ref, zuncsd2by1_obj->ldu2, zuncsd2by1_obj->v1tref, zuncsd2by1_obj->ldv1t);

    /* Compute libflame's Lapacke o/p  */
    zuncsd2by1_obj->info = LAPACKE_zuncsd2by1(zuncsd2by1_obj->matrix_layout, zuncsd2by1_obj->jobu1, zuncsd2by1_obj->jobu2, zuncsd2by1_obj->jobv1t,
								zuncsd2by1_obj->m, zuncsd2by1_obj->p,zuncsd2by1_obj->q, zuncsd2by1_obj->x11,
								zuncsd2by1_obj->ldx11, zuncsd2by1_obj->x21, zuncsd2by1_obj->ldx21, zuncsd2by1_obj->theta, zuncsd2by1_obj->u1, 
								zuncsd2by1_obj->ldu1, zuncsd2by1_obj->u2, zuncsd2by1_obj->ldu2, zuncsd2by1_obj->v1t, zuncsd2by1_obj->ldv1t);

    if( zuncsd2by1_obj->info < 0 ) {
        printf( "\n warning: The i:%x22 th argument with libflame \
        LAPACKE_zuncsd2by1 is wrong\n", zuncsd2by1_obj->info );
    }
    if( zuncsd2by1_obj->inforef < 0 ) {
        printf( "The i:%x22 th argument with Netlib LAPACKE_zuncsd2by1 is wrong\n", 
        zuncsd2by1_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zuncsd2by1_obj->diff_x11 =  computeDiff_z( zuncsd2by1_obj->bufsize_x11, zuncsd2by1_obj->x11, zuncsd2by1_obj->x11ref );
	zuncsd2by1_obj->diff_x21 =  computeDiff_z( zuncsd2by1_obj->bufsize_x21, zuncsd2by1_obj->x21, zuncsd2by1_obj->x21ref );
	zuncsd2by1_obj->diff_u1 =  computeDiff_z(zuncsd2by1_obj->bufsize_u1, zuncsd2by1_obj->u1, zuncsd2by1_obj->u1ref );
	zuncsd2by1_obj->diff_u2 =  computeDiff_z( zuncsd2by1_obj->bufsize_u2, zuncsd2by1_obj->u2, zuncsd2by1_obj->u2ref );
	zuncsd2by1_obj->diff_v1t =  computeDiff_z( zuncsd2by1_obj->bufsize_v1t, zuncsd2by1_obj->v1t, zuncsd2by1_obj->v1tref );
	zuncsd2by1_obj->diff_theta =  computeDiff_d( zuncsd2by1_obj->bufsize_theta, zuncsd2by1_obj->theta, zuncsd2by1_obj->thetaref );

}

TEST_F(zuncsd2by1_test, zuncsd2by11) {
    EXPECT_NEAR(0.0, zuncsd2by1_obj->diff_x11, zuncsd2by1_obj->uncsd2by1_threshold);
	EXPECT_NEAR(0.0, zuncsd2by1_obj->diff_x21, zuncsd2by1_obj->uncsd2by1_threshold);
	EXPECT_NEAR(0.0, zuncsd2by1_obj->diff_u1, zuncsd2by1_obj->uncsd2by1_threshold);
	EXPECT_NEAR(0.0, zuncsd2by1_obj->diff_u2, zuncsd2by1_obj->uncsd2by1_threshold);
	EXPECT_NEAR(0.0, zuncsd2by1_obj->diff_v1t, zuncsd2by1_obj->uncsd2by1_threshold);
	EXPECT_NEAR(0.0, zuncsd2by1_obj->diff_theta, zuncsd2by1_obj->uncsd2by1_threshold);
}

TEST_F(zuncsd2by1_test, zuncsd2by12) {
    EXPECT_NEAR(0.0, zuncsd2by1_obj->diff_x11, zuncsd2by1_obj->uncsd2by1_threshold);
	EXPECT_NEAR(0.0, zuncsd2by1_obj->diff_x21, zuncsd2by1_obj->uncsd2by1_threshold);
	EXPECT_NEAR(0.0, zuncsd2by1_obj->diff_u1, zuncsd2by1_obj->uncsd2by1_threshold);
	EXPECT_NEAR(0.0, zuncsd2by1_obj->diff_u2, zuncsd2by1_obj->uncsd2by1_threshold);
	EXPECT_NEAR(0.0, zuncsd2by1_obj->diff_v1t, zuncsd2by1_obj->uncsd2by1_threshold);
	EXPECT_NEAR(0.0, zuncsd2by1_obj->diff_theta, zuncsd2by1_obj->uncsd2by1_threshold);
}

TEST_F(zuncsd2by1_test, zuncsd2by13) {
    EXPECT_NEAR(0.0, zuncsd2by1_obj->diff_x11, zuncsd2by1_obj->uncsd2by1_threshold);
	EXPECT_NEAR(0.0, zuncsd2by1_obj->diff_x21, zuncsd2by1_obj->uncsd2by1_threshold);
	EXPECT_NEAR(0.0, zuncsd2by1_obj->diff_u1, zuncsd2by1_obj->uncsd2by1_threshold);
	EXPECT_NEAR(0.0, zuncsd2by1_obj->diff_u2, zuncsd2by1_obj->uncsd2by1_threshold);
	EXPECT_NEAR(0.0, zuncsd2by1_obj->diff_v1t, zuncsd2by1_obj->uncsd2by1_threshold);
	EXPECT_NEAR(0.0, zuncsd2by1_obj->diff_theta, zuncsd2by1_obj->uncsd2by1_threshold);

}

TEST_F(zuncsd2by1_test, zuncsd2by14) {
    EXPECT_NEAR(0.0, zuncsd2by1_obj->diff_x11, zuncsd2by1_obj->uncsd2by1_threshold);
	EXPECT_NEAR(0.0, zuncsd2by1_obj->diff_x21, zuncsd2by1_obj->uncsd2by1_threshold);
	EXPECT_NEAR(0.0, zuncsd2by1_obj->diff_u1, zuncsd2by1_obj->uncsd2by1_threshold);
	EXPECT_NEAR(0.0, zuncsd2by1_obj->diff_u2, zuncsd2by1_obj->uncsd2by1_threshold);
	EXPECT_NEAR(0.0, zuncsd2by1_obj->diff_v1t, zuncsd2by1_obj->uncsd2by1_threshold);
	EXPECT_NEAR(0.0, zuncsd2by1_obj->diff_theta, zuncsd2by1_obj->uncsd2by1_threshold);
}