#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define unbdb_free() \
if (x11!=NULL)    free(x11); \
if (x11ref!=NULL) free(x11ref);\
if (x12 != NULL) free(x12); \
if (x12ref != NULL) free(x12ref); \
if (x21 != NULL) free(x21); \
if (x21ref!=NULL) free(x21ref);\
if (taup1!=NULL)  free(taup1);\
if (taup1ref!=NULL) free(taup1ref); \
if (taup2!=NULL)  free(taup2);\
if (taup2ref!=NULL) free(taup2ref); \
if (tauq1!=NULL)  free(tauq1);\
if (tauq1ref!=NULL) free(tauq1ref); \
if (tauq2!=NULL)  free(tauq2);\
if (tauq2ref!=NULL) free(tauq2ref); \
if (theta != NULL) free(theta); \
if (thetaref != NULL) free(thetaref); \
if (phi!=NULL) free(phi);\
if (phiref!=NULL) free(phiref);\
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin scomplex_common_parameters  class definition */
class unbdb_scomplex_parameters{

   public:
	int bufsize_x11, bufsize_x12, bufsize_x21, bufsize_x22, bufsize_phi, bufsize_theta;
	int bufsize_taup1, bufsize_taup2, bufsize_tauq1, bufsize_tauq2;
	void *hModule, *dModule;
	float diff_x11, diff_x12, diff_x21, diff_x22, diff_theta, diff_phi;
	float diff_taup1, diff_taup2, diff_tauq1, diff_tauq2;
   /*input parameters */
	int matrix_layout;
	lapack_int p;
	lapack_int m;
	lapack_int q;
	lapack_complex_float* x11, *x12, *x21, *x22;
	char signs;
	char trans;
	lapack_int unbdb_threshold;
	lapack_int ldx11, ldx12, ldx21, ldx22;
	/*Output Parameter*/
	float* theta, *phi;
	float* thetaref, *phiref;
	lapack_complex_float *taup1, *taup2, *tauq1, *tauq2;
	lapack_complex_float *taup1ref, *taup2ref, *tauq1ref, *tauq2ref;
	lapack_complex_float *x11ref, *x12ref, *x21ref, *x22ref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      unbdb_scomplex_parameters (int matrix_layout, char trans, char signs, lapack_int m, lapack_int p, lapack_int q, lapack_int unbdb_threshold);
      ~unbdb_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
unbdb_scomplex_parameters:: unbdb_scomplex_parameters (int matrix_layout_i, char trans_i, char signs_i, lapack_int m_i, lapack_int p_i, lapack_int q_i, lapack_int unbdb_threshold_i)
{
	matrix_layout = matrix_layout_i;
	m = m_i;
	p = p_i;
	q = q_i;
	m = m+p;
	unbdb_threshold = unbdb_threshold_i;
	signs = signs_i;
	trans = trans_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n unbdb scomplex: matrix_layout = %d, trans:%c, signs:%c, m: %d, p:%d, q:%d \n", matrix_layout, trans, signs, m, p, q);
	#endif
	
	/*check for signs*/
	if (signs == 'L')
		signs = 'O';
	
	p < min(m-p,m-q)? q = p: q = min(m-p,m-q);
	
	if (trans == 'C') // Trans is either 'N' or 'T'
		trans = 'N';
	
	if (trans == 'T')
	{
		if (matrix_layout == LAPACK_COL_MAJOR)
		{	
			ldx11 = p;
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{	
			ldx11 = q;
		}else
			EXPECT_TRUE(false) << "matrix_layout invalid";
		ldx12 = m-q;
		ldx21 = q;
		ldx22 = m - q;
	}else if (trans == 'N')
	{
		if (matrix_layout == LAPACK_COL_MAJOR)
		{	
			ldx12 = p;
			ldx21 = m-p;
			ldx22 = m-p;
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{	
			ldx12 = m-q;
			ldx21 = q;
			ldx22 = m - q;
		}else
			EXPECT_TRUE(false) << "matrix_layout invalid";		
		ldx11 = p;
	} else {		
		ldx11 = q;
		ldx12 = m-q;
		ldx21 = q;
		ldx22 = m - q;
	}
	/* bufsizes based on matrix layout*/
	if (matrix_layout == LAPACK_COL_MAJOR)
	{	
		bufsize_x11 = ldx11*q;
		bufsize_x12 = ldx12*(m-q);
		bufsize_x21 = ldx21*q;
		bufsize_x22 = ldx22*(m-q);
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{	
		bufsize_x11 = ldx11*p;
		bufsize_x12 = ldx12*p;
		bufsize_x21 = ldx21*(m-p);
		bufsize_x22 = ldx22*(m - p);
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	bufsize_taup1 = p;
	bufsize_taup2 = m-p;
	bufsize_tauq1 = q;
	bufsize_tauq2 = m-q;
	bufsize_phi = q -1;
	bufsize_theta = q;	
		
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&x11, &x11ref, bufsize_x11);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&x12, &x12ref, bufsize_x12);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&x21, &x21ref, bufsize_x21);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&x22, &x22ref, bufsize_x22);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&taup1, &taup1ref, bufsize_taup1);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&taup2, &taup2ref, bufsize_taup2);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&tauq1, &tauq1ref, bufsize_tauq1);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&tauq2, &tauq2ref, bufsize_tauq2);
	lapacke_gtest_alloc_float_buffer_pair(&theta, &thetaref, bufsize_tauq1);
	lapacke_gtest_alloc_float_buffer_pair(&phi, &phiref, (bufsize_tauq1-1));

	if ((x11==NULL) || (x11ref==NULL) || \
		(x12==NULL) || (x12ref==NULL) || \
		(x21==NULL) || (x21ref==NULL) || \
		(x22==NULL) || (x22ref==NULL) || \
		(taup1==NULL) || (taup1ref==NULL) ||\
		(taup2==NULL) || (taup2ref==NULL) ||\
		(tauq1==NULL) || (tauq1ref==NULL) ||\
		(tauq2==NULL) || (tauq2ref==NULL) ||\
		(theta==NULL) || (thetaref==NULL) ||\
		(phi==NULL) || (phiref==NULL)){
		EXPECT_FALSE( true) << "unbdb_float_parameters object: malloc error.";
		unbdb_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( x11, x11ref, bufsize_x11);
	lapacke_gtest_init_scomplex_buffer_pair_rand( x12, x12ref, bufsize_x12);
	lapacke_gtest_init_scomplex_buffer_pair_rand( x21, x21ref, bufsize_x21);
	lapacke_gtest_init_scomplex_buffer_pair_rand( x22, x22ref, bufsize_x22);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
unbdb_scomplex_parameters :: ~unbdb_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" unbdb_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   unbdb_free();

}
/*  Test fixture class definition */
class cunbdb_test  : public  ::testing::Test {
public:
   unbdb_scomplex_parameters  *cunbdb_obj;
   void SetUp();
   void TearDown () { delete cunbdb_obj;}
};

void cunbdb_test::SetUp(){

	/* LAPACKE cunbdb prototype */
    typedef int (*Fptr_NL_LAPACKE_cunbdb) (int matrix_layout, char trans, char signs,\
                           lapack_int m, lapack_int p, lapack_int q,\
                           lapack_complex_float* x11, lapack_int ldx11,\
                           lapack_complex_float* x12, lapack_int ldx12,\
                           lapack_complex_float* x21, lapack_int ldx21,\
                           lapack_complex_float* x22, lapack_int ldx22,\
                           float* theta, float* phi,\
                           lapack_complex_float* taup1,\
                           lapack_complex_float* taup2,\
                           lapack_complex_float* tauq1,\
                           lapack_complex_float* tauq2);

    Fptr_NL_LAPACKE_cunbdb cunbdb;

    cunbdb_obj = new unbdb_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].trans,
                           eig_paramslist[idx].uplo,
						   eig_paramslist[idx].m,
						   eig_paramslist[idx].m,
						   eig_paramslist[idx].m,
						   eig_paramslist[idx].threshold_value);
						   

    idx = Circular_Increment_Index(idx);

    cunbdb_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cunbdb_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cunbdb_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cunbdb_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*cunbdb library call */
    cunbdb = (Fptr_NL_LAPACKE_cunbdb)dlsym(cunbdb_obj->hModule, "LAPACKE_cunbdb");
    ASSERT_TRUE(cunbdb != NULL) << "failed to get the Netlib LAPACKE_cunbdb symbol";
	
/*Compute cunbdb's  o/p */
    cunbdb_obj->inforef = cunbdb( cunbdb_obj->matrix_layout, cunbdb_obj->trans, cunbdb_obj->signs, cunbdb_obj->m, cunbdb_obj->p,
								cunbdb_obj->q, cunbdb_obj->x11ref, cunbdb_obj->ldx11, cunbdb_obj->x12ref, cunbdb_obj->ldx12, 
								cunbdb_obj->x21ref, cunbdb_obj->ldx21, cunbdb_obj->x22ref, cunbdb_obj->ldx22, cunbdb_obj->thetaref, cunbdb_obj->phiref,
								cunbdb_obj->taup1ref, cunbdb_obj->taup2ref, cunbdb_obj->tauq1ref, cunbdb_obj->tauq2ref);

    /* Compute libflame's Lapacke o/p  */
    cunbdb_obj->info = LAPACKE_cunbdb( cunbdb_obj->matrix_layout, cunbdb_obj->trans, cunbdb_obj->signs, cunbdb_obj->m, cunbdb_obj->p,
								cunbdb_obj->q, cunbdb_obj->x11, cunbdb_obj->ldx11, cunbdb_obj->x12, cunbdb_obj->ldx12, 
								cunbdb_obj->x21, cunbdb_obj->ldx21, cunbdb_obj->x22, cunbdb_obj->ldx22, cunbdb_obj->theta, cunbdb_obj->phi,
								cunbdb_obj->taup1, cunbdb_obj->taup2, cunbdb_obj->tauq1, cunbdb_obj->tauq2);

    if( cunbdb_obj->info < 0 ) {
        printf( "\n warning: The i:%x22 th argument with libflame \
        LAPACKE_cunbdb is wrong\n", cunbdb_obj->info );
    }
    if( cunbdb_obj->inforef < 0 ) {
        printf( "The i:%x22 th argument with Netlib LAPACKE_cunbdb is wrong\n", 
        cunbdb_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cunbdb_obj->diff_x11 =  computeDiff_c( cunbdb_obj->bufsize_x11, cunbdb_obj->x11, cunbdb_obj->x11ref );
	cunbdb_obj->diff_x12 =  computeDiff_c( cunbdb_obj->bufsize_x12, cunbdb_obj->x12, cunbdb_obj->x12ref );
	cunbdb_obj->diff_x21 =  computeDiff_c( cunbdb_obj->bufsize_x21, cunbdb_obj->x21, cunbdb_obj->x21ref );
	cunbdb_obj->diff_x22 =  computeDiff_c( cunbdb_obj->bufsize_x22, cunbdb_obj->x21, cunbdb_obj->x22ref );
	cunbdb_obj->diff_taup1 =  computeDiff_c( cunbdb_obj->bufsize_taup1, cunbdb_obj->taup1, cunbdb_obj->taup1ref );
	cunbdb_obj->diff_taup2 =  computeDiff_c( cunbdb_obj->bufsize_taup2, cunbdb_obj->taup2, cunbdb_obj->taup2ref );
	cunbdb_obj->diff_tauq1 =  computeDiff_c( cunbdb_obj->bufsize_tauq1, cunbdb_obj->tauq1, cunbdb_obj->tauq1ref );
	cunbdb_obj->diff_tauq2 =  computeDiff_c( cunbdb_obj->bufsize_tauq2, cunbdb_obj->tauq2, cunbdb_obj->tauq2ref );
	cunbdb_obj->diff_theta =  computeDiff_s( cunbdb_obj->bufsize_theta, cunbdb_obj->theta, cunbdb_obj->thetaref );
	cunbdb_obj->diff_phi =  computeDiff_s( cunbdb_obj->bufsize_phi, cunbdb_obj->phi, cunbdb_obj->phiref );

}

TEST_F(cunbdb_test, cunbdb1) {
    EXPECT_NEAR(0.0, cunbdb_obj->diff_x11, cunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, cunbdb_obj->diff_x12, cunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, cunbdb_obj->diff_x21, cunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, cunbdb_obj->diff_x22, cunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, cunbdb_obj->diff_taup1, cunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, cunbdb_obj->diff_taup2, cunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, cunbdb_obj->diff_tauq1, cunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, cunbdb_obj->diff_tauq2, cunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, cunbdb_obj->diff_theta, cunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, cunbdb_obj->diff_phi, cunbdb_obj->unbdb_threshold);
}

TEST_F(cunbdb_test, cunbdb2) {
    EXPECT_NEAR(0.0, cunbdb_obj->diff_x11, cunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, cunbdb_obj->diff_x12, cunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, cunbdb_obj->diff_x21, cunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, cunbdb_obj->diff_x22, cunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, cunbdb_obj->diff_taup1, cunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, cunbdb_obj->diff_taup2, cunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, cunbdb_obj->diff_tauq1, cunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, cunbdb_obj->diff_tauq2, cunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, cunbdb_obj->diff_theta, cunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, cunbdb_obj->diff_phi, cunbdb_obj->unbdb_threshold);
}

TEST_F(cunbdb_test, cunbdb3) {
    EXPECT_NEAR(0.0, cunbdb_obj->diff_x11, cunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, cunbdb_obj->diff_x12, cunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, cunbdb_obj->diff_x21, cunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, cunbdb_obj->diff_x22, cunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, cunbdb_obj->diff_taup1, cunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, cunbdb_obj->diff_taup2, cunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, cunbdb_obj->diff_tauq1, cunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, cunbdb_obj->diff_tauq2, cunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, cunbdb_obj->diff_theta, cunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, cunbdb_obj->diff_phi, cunbdb_obj->unbdb_threshold);
}

TEST_F(cunbdb_test, cunbdb4) {
    EXPECT_NEAR(0.0, cunbdb_obj->diff_x11, cunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, cunbdb_obj->diff_x12, cunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, cunbdb_obj->diff_x21, cunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, cunbdb_obj->diff_x22, cunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, cunbdb_obj->diff_taup1, cunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, cunbdb_obj->diff_taup2, cunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, cunbdb_obj->diff_tauq1, cunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, cunbdb_obj->diff_tauq2, cunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, cunbdb_obj->diff_theta, cunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, cunbdb_obj->diff_phi, cunbdb_obj->unbdb_threshold);
}

/* Begin dcomplex_common_parameters  class definition */
class unbdb_dcomplex_parameters{

   public:
	int bufsize_x11, bufsize_x12, bufsize_x21, bufsize_x22, bufsize_phi, bufsize_theta;
	int bufsize_taup1, bufsize_taup2, bufsize_tauq1, bufsize_tauq2;
	void *hModule, *dModule;
	double diff_x11, diff_x12, diff_x21, diff_x22, diff_theta, diff_phi;
	double diff_taup1, diff_taup2, diff_tauq1, diff_tauq2;
   /*input parameters */
	int matrix_layout;
	lapack_int p;
	lapack_int m;
	lapack_int q;
	lapack_int unbdb_threshold;
	lapack_complex_double* x11, *x12, *x21, *x22;
	char signs;
	char trans;
	lapack_int ldx11, ldx12, ldx21, ldx22;
	/*Output Parameter*/
	double* theta, *phi;
	double* thetaref, *phiref;
	lapack_complex_double *taup1, *taup2, *tauq1, *tauq2;
	lapack_complex_double *taup1ref, *taup2ref, *tauq1ref, *tauq2ref;
	lapack_complex_double *x11ref, *x12ref, *x21ref, *x22ref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      unbdb_dcomplex_parameters (int matrix_layout, char trans, char signs, lapack_int m, lapack_int p, lapack_int q, lapack_int unbdb_threshold);
      ~unbdb_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
unbdb_dcomplex_parameters:: unbdb_dcomplex_parameters (int matrix_layout_i, char trans_i, char signs_i, lapack_int m_i, lapack_int p_i, lapack_int q_i, lapack_int unbdb_threshold_i)
{
	matrix_layout = matrix_layout_i;
	m = m_i;
	p = p_i;
	q = q_i;
	m = m+p;
	unbdb_threshold = unbdb_threshold_i;
	signs = signs_i;
	trans = trans_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n unbdb dcomplex: matrix_layout = %d, trans:%c, signs:%c, m: %d, p:%d, q:%d \n", matrix_layout, trans, signs, m, p, q);
	#endif
	
	/*check for signs*/
	if (signs == 'L')
		signs = 'O';
	
	p < min(m-p,m-q)? q = p: q = min(m-p,m-q);

	if (trans == 'C') // Trans is either 'N' or 'T'
		trans = 'N';
	
	
	if (trans == 'T')
	{
		if (matrix_layout == LAPACK_COL_MAJOR)
		{	
			ldx11 = p;
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{	
			ldx11 = q;
		}else
			EXPECT_TRUE(false) << "matrix_layout invalid";
		ldx12 = m-q;
		ldx21 = q;
		ldx22 = m - q;
	}else if (trans == 'N')
	{
		if (matrix_layout == LAPACK_COL_MAJOR)
		{	
			ldx12 = p;
			ldx21 = m-p;
			ldx22 = m-p;
		}else if (matrix_layout == LAPACK_ROW_MAJOR)
		{	
			ldx12 = m-q;
			ldx21 = q;
			ldx22 = m - q;
		}else
			EXPECT_TRUE(false) << "matrix_layout invalid";		
		ldx11 = p;
	} else {		
		ldx11 = q;
		ldx12 = m-q;
		ldx21 = q;
		ldx22 = m - q;
	}
	/* bufsizes based on matrix layout*/
	if (matrix_layout == LAPACK_COL_MAJOR)
	{	
		bufsize_x11 = ldx11*q;
		bufsize_x12 = ldx12*(m-q);
		bufsize_x21 = ldx21*q;
		bufsize_x22 = ldx22*(m-q);
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{	
		bufsize_x11 = ldx11*p;
		bufsize_x12 = ldx12*p;
		bufsize_x21 = ldx21*(m-p);
		bufsize_x22 = ldx22*(m - p);
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	bufsize_taup1 = p;
	bufsize_taup2 = m-p;
	bufsize_tauq1 = q;
	bufsize_tauq2 = m-q;
	bufsize_phi = q -1;
	bufsize_theta = q;

		
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&x11, &x11ref, bufsize_x11);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&x12, &x12ref, bufsize_x12);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&x21, &x21ref, bufsize_x21);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&x22, &x22ref, bufsize_x22);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&taup1, &taup1ref, bufsize_taup1);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&taup2, &taup2ref, bufsize_taup2);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&tauq1, &tauq1ref, bufsize_tauq1);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&tauq2, &tauq2ref, bufsize_tauq2);
	lapacke_gtest_alloc_double_buffer_pair(&theta, &thetaref, bufsize_tauq1);
	lapacke_gtest_alloc_double_buffer_pair(&phi, &phiref, (bufsize_tauq1-1));

	if ((x11==NULL) || (x11ref==NULL) || \
		(x12==NULL) || (x12ref==NULL) || \
		(x21==NULL) || (x21ref==NULL) || \
		(x22==NULL) || (x22ref==NULL) || \
		(taup1==NULL) || (taup1ref==NULL) ||\
		(taup2==NULL) || (taup2ref==NULL) ||\
		(tauq1==NULL) || (tauq1ref==NULL) ||\
		(tauq2==NULL) || (tauq2ref==NULL) ||\
		(theta==NULL) || (thetaref==NULL) ||\
		(phi==NULL) || (phiref==NULL)){
		EXPECT_FALSE( true) << "unbdb_double_parameters object: malloc error.";
		unbdb_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( x11, x11ref, bufsize_x11);
	lapacke_gtest_init_dcomplex_buffer_pair_rand( x12, x12ref, bufsize_x12);
	lapacke_gtest_init_dcomplex_buffer_pair_rand( x21, x21ref, bufsize_x21);
	lapacke_gtest_init_dcomplex_buffer_pair_rand( x22, x22ref, bufsize_x22);

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
unbdb_dcomplex_parameters :: ~unbdb_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" unbdb_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   unbdb_free();

}
/*  Test fixture class definition */
class zunbdb_test  : public  ::testing::Test {
public:
   unbdb_dcomplex_parameters  *zunbdb_obj;
   void SetUp();
   void TearDown () { delete zunbdb_obj;}
};

void zunbdb_test::SetUp(){

	/* LAPACKE zunbdb prototype */
    typedef int (*Fptr_NL_LAPACKE_zunbdb) (int matrix_layout, char trans, char signs,\
                           lapack_int m, lapack_int p, lapack_int q,\
                           lapack_complex_double* x11, lapack_int ldx11,\
                           lapack_complex_double* x12, lapack_int ldx12,\
                           lapack_complex_double* x21, lapack_int ldx21,\
                           lapack_complex_double* x22, lapack_int ldx22,\
                           double* theta, double* phi,\
                           lapack_complex_double* taup1,\
                           lapack_complex_double* taup2,\
                           lapack_complex_double* tauq1,\
                           lapack_complex_double* tauq2);

    Fptr_NL_LAPACKE_zunbdb zunbdb;

    zunbdb_obj = new unbdb_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].trans,
                           eig_paramslist[idx].uplo,
						   eig_paramslist[idx].m,
						   eig_paramslist[idx].m,
						   eig_paramslist[idx].m,
						   eig_paramslist[idx].threshold_value);
						   

    idx = Circular_Increment_Index(idx);

    zunbdb_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zunbdb_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zunbdb_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zunbdb_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*zunbdb library call */
    zunbdb = (Fptr_NL_LAPACKE_zunbdb)dlsym(zunbdb_obj->hModule, "LAPACKE_zunbdb");
    ASSERT_TRUE(zunbdb != NULL) << "failed to get the Netlib LAPACKE_zunbdb symbol";
	
/*Compute zunbdb's  o/p */
    zunbdb_obj->inforef = zunbdb( zunbdb_obj->matrix_layout, zunbdb_obj->trans, zunbdb_obj->signs, zunbdb_obj->m, zunbdb_obj->p,
								zunbdb_obj->q, zunbdb_obj->x11ref, zunbdb_obj->ldx11, zunbdb_obj->x12ref, zunbdb_obj->ldx12, 
								zunbdb_obj->x21ref, zunbdb_obj->ldx21, zunbdb_obj->x22ref, zunbdb_obj->ldx22, zunbdb_obj->thetaref, zunbdb_obj->phiref,
								zunbdb_obj->taup1ref, zunbdb_obj->taup2ref, zunbdb_obj->tauq1ref, zunbdb_obj->tauq2ref);

    /* Compute libflame's Lapacke o/p  */
    zunbdb_obj->info = LAPACKE_zunbdb( zunbdb_obj->matrix_layout, zunbdb_obj->trans, zunbdb_obj->signs, zunbdb_obj->m, zunbdb_obj->p,
								zunbdb_obj->q, zunbdb_obj->x11, zunbdb_obj->ldx11, zunbdb_obj->x12, zunbdb_obj->ldx12, 
								zunbdb_obj->x21, zunbdb_obj->ldx21, zunbdb_obj->x22, zunbdb_obj->ldx22, zunbdb_obj->theta, zunbdb_obj->phi,
								zunbdb_obj->taup1, zunbdb_obj->taup2, zunbdb_obj->tauq1, zunbdb_obj->tauq2);

    if( zunbdb_obj->info < 0 ) {
        printf( "\n warning: The i:%x22 th argument with libflame \
        LAPACKE_zunbdb is wrong\n", zunbdb_obj->info );
    }
    if( zunbdb_obj->inforef < 0 ) {
        printf( "The i:%x22 th argument with Netlib LAPACKE_zunbdb is wrong\n", 
        zunbdb_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zunbdb_obj->diff_x11 =  computeDiff_z( zunbdb_obj->bufsize_x11, zunbdb_obj->x11, zunbdb_obj->x11ref );
	zunbdb_obj->diff_x12 =  computeDiff_z( zunbdb_obj->bufsize_x12, zunbdb_obj->x12, zunbdb_obj->x12ref );
	zunbdb_obj->diff_x21 =  computeDiff_z( zunbdb_obj->bufsize_x21, zunbdb_obj->x21, zunbdb_obj->x21ref );
	zunbdb_obj->diff_x22 =  computeDiff_z( zunbdb_obj->bufsize_x22, zunbdb_obj->x21, zunbdb_obj->x22ref );
	zunbdb_obj->diff_taup1 =  computeDiff_z( zunbdb_obj->bufsize_taup1, zunbdb_obj->taup1, zunbdb_obj->taup1ref );
	zunbdb_obj->diff_taup2 =  computeDiff_z( zunbdb_obj->bufsize_taup2, zunbdb_obj->taup2, zunbdb_obj->taup2ref );
	zunbdb_obj->diff_tauq1 =  computeDiff_z( zunbdb_obj->bufsize_tauq1, zunbdb_obj->tauq1, zunbdb_obj->tauq1ref );
	zunbdb_obj->diff_tauq2 =  computeDiff_z( zunbdb_obj->bufsize_tauq2, zunbdb_obj->tauq2, zunbdb_obj->tauq2ref );
	zunbdb_obj->diff_theta =  computeDiff_d( zunbdb_obj->bufsize_theta, zunbdb_obj->theta, zunbdb_obj->thetaref );
	zunbdb_obj->diff_phi =  computeDiff_d( zunbdb_obj->bufsize_phi, zunbdb_obj->phi, zunbdb_obj->phiref );

}

TEST_F(zunbdb_test, zunbdb1) {
    EXPECT_NEAR(0.0, zunbdb_obj->diff_x11, zunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, zunbdb_obj->diff_x12, zunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, zunbdb_obj->diff_x21, zunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, zunbdb_obj->diff_x22, zunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, zunbdb_obj->diff_taup1, zunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, zunbdb_obj->diff_taup2, zunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, zunbdb_obj->diff_tauq1, zunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, zunbdb_obj->diff_tauq2, zunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, zunbdb_obj->diff_theta, zunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, zunbdb_obj->diff_phi, zunbdb_obj->unbdb_threshold);
}

TEST_F(zunbdb_test, zunbdb2) {
    EXPECT_NEAR(0.0, zunbdb_obj->diff_x11, zunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, zunbdb_obj->diff_x12, zunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, zunbdb_obj->diff_x21, zunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, zunbdb_obj->diff_x22, zunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, zunbdb_obj->diff_taup1, zunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, zunbdb_obj->diff_taup2, zunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, zunbdb_obj->diff_tauq1, zunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, zunbdb_obj->diff_tauq2, zunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, zunbdb_obj->diff_theta, zunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, zunbdb_obj->diff_phi, zunbdb_obj->unbdb_threshold);
}

TEST_F(zunbdb_test, zunbdb3) {
    EXPECT_NEAR(0.0, zunbdb_obj->diff_x11, zunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, zunbdb_obj->diff_x12, zunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, zunbdb_obj->diff_x21, zunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, zunbdb_obj->diff_x22, zunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, zunbdb_obj->diff_taup1, zunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, zunbdb_obj->diff_taup2, zunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, zunbdb_obj->diff_tauq1, zunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, zunbdb_obj->diff_tauq2, zunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, zunbdb_obj->diff_theta, zunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, zunbdb_obj->diff_phi, zunbdb_obj->unbdb_threshold);
}

TEST_F(zunbdb_test, zunbdb4) {
    EXPECT_NEAR(0.0, zunbdb_obj->diff_x11, zunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, zunbdb_obj->diff_x12, zunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, zunbdb_obj->diff_x21, zunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, zunbdb_obj->diff_x22, zunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, zunbdb_obj->diff_taup1, zunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, zunbdb_obj->diff_taup2, zunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, zunbdb_obj->diff_tauq1, zunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, zunbdb_obj->diff_tauq2, zunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, zunbdb_obj->diff_theta, zunbdb_obj->unbdb_threshold);
	EXPECT_NEAR(0.0, zunbdb_obj->diff_phi, zunbdb_obj->unbdb_threshold);
}
