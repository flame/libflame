#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"

#define gesvdq_free() \
if (a!=NULL)    free(a); \
if (aref!=NULL) free(aref);\
if (u != NULL) free(u); \
if (uref != NULL) free(uref); \
if (s != NULL) free(s); \
if (sref != NULL) free(sref); \
if (numrank != NULL) free(numrank);\
if (numrankref != NULL) free(numrankref);\
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin float_common_parameters  class definition */
class gesvdq_float_parameters{

   public:
	int bufsize_a, bufsize_u, bufsize_s, bufsize_v;
	void *hModule, *dModule;
	float diff_a, diff_u, diff_s, diff_v;
   /*input parameters */
	int matrix_layout;
	lapack_int m;
	lapack_int n;
	float* a, *v;
	char joba, jobp, jobr, jobu, jobv;
	lapack_int lda, ldu, ldv;
	/*Output Parameter*/
	float* s;
	float* sref;
	lapack_int *numrank, *numrankref;
	float *u;
	float *aref, *uref, *vref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      gesvdq_float_parameters (int matrix_layout, char joba, char jobp, char jobr, char jobu, char jobv,  lapack_int m, lapack_int n);
      ~gesvdq_float_parameters ();

};

/* Constructor definition  float_common_parameters */
gesvdq_float_parameters:: gesvdq_float_parameters (int matrix_layout_i, char joba_i, char jobp_i, char jobr_i, char jobu_i, char jobv_i,  lapack_int m_i, lapack_int n_i)
{
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	joba = joba_i; //A, H, M, E
	jobp = jobp_i; // P, N
	jobr = jobr_i; // T, N
	jobu = jobu_i; //A, F, S, R, N
	jobv = jobv_i; // A, V, R, N
	
//#if LAPACKE_TEST_VERBOSE
		printf(" \n gesvdq float: matrix_layout = %d, joba:%c, jobp:%c, jobr:%c, jobu:%c, jobv:%c, m:%d, n:%d\n",  matrix_layout,  joba, jobp, jobr, jobu, jobv, m, n);	
	//#endif
	
	/*if ( jobr == 'C')
		jobr = 'N';
		
	if (jobp == 'U')
		jobp == 'P';
	
	if (jobu == 'V')
		jobu = 'S';*/
	
	/* Buffer sizes */
	/*buffer A*/
	lda = m;
	bufsize_a = lda *n;
	/*Buffer s */
	bufsize_s = n;
	/*buffer U*/
	if ((jobu == 'A') ||(jobu  == 'S') ||(jobu == 'U')  ||( jobu == 'R'))
		ldu = m;
	else if (jobu == 'F')
		ldu = n;
	/*Buffer V*/
	ldv = n;
	bufsize_v = ldv *n;	
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&a, &aref, bufsize_a);
	lapacke_gtest_alloc_float_buffer_pair(&u, &uref, bufsize_u);
	lapacke_gtest_alloc_float_buffer_pair(&v, &vref, bufsize_v);
	lapacke_gtest_alloc_float_buffer_pair(&s, &sref, bufsize_s);
	lapacke_gtest_alloc_int_buffer_pair(&numrank, &numrankref, n);

	if ((a==NULL) || (aref==NULL) || \
		(u==NULL) || (uref==NULL) || \
		(v==NULL) || (vref==NULL) ||\
		(s==NULL) || (sref==NULL) ||\
		(numrank==NULL) || (numrankref==NULL)){
		EXPECT_FALSE( true) << "gesvdq_float_parameters object: malloc error.";
		gesvdq_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand( a, aref, bufsize_a);
	lapacke_gtest_init_int_buffer_pair_with_constant(numrank, numrankref, n, 0);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
gesvdq_float_parameters :: ~gesvdq_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" gesvdq_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   gesvdq_free();

}
/*  Test fixture class definition */
class sgesvdq_test  : public  ::testing::Test {
public:
   gesvdq_float_parameters  *sgesvdq_obj;
   void SetUp();
   void TearDown () { delete sgesvdq_obj;}
};

void sgesvdq_test::SetUp(){

	/* LAPACKE sgesvdq prototype */
    typedef int (*Fptr_NL_LAPACKE_sgesvdq) ( int matrix_layout, char joba, char jobp, char jobr, char jobu, char jobv,\
                           lapack_int m, lapack_int n, float* a,\
                           lapack_int lda, float* s, float* u,\
                           lapack_int ldu, float* v,\
                           lapack_int ldv, lapack_int* numrank);

    Fptr_NL_LAPACKE_sgesvdq sgesvdq;

    sgesvdq_obj = new gesvdq_float_parameters ( svd_paramslist[idx].matrix_layout,
                           svd_paramslist[idx].joba_gesvdq,
                           svd_paramslist[idx].jobp_gejsv,
						   svd_paramslist[idx].jobt_gejsv,
						   svd_paramslist[idx].jobu_gesvdq,
						   svd_paramslist[idx].jobv_gesvdq,
						   svd_paramslist[idx].m_gejsv,
						   svd_paramslist[idx].n_gejsv);
						   

    idx = Circular_Increment_Index(idx);

    sgesvdq_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgesvdq_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgesvdq_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgesvdq_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*sgesvdq library call */
    sgesvdq = (Fptr_NL_LAPACKE_sgesvdq)dlsym(sgesvdq_obj->hModule, "LAPACKE_sgesvdq");
    ASSERT_TRUE(sgesvdq != NULL) << "failed to get the Netlib LAPACKE_sgesvdq symbol";
	
/*Compute sgesvdq's  o/p */
    sgesvdq_obj->inforef = sgesvdq( sgesvdq_obj->matrix_layout, sgesvdq_obj->joba, sgesvdq_obj->jobp, sgesvdq_obj->jobr, sgesvdq_obj->jobu,\
								sgesvdq_obj->jobv, sgesvdq_obj->m, sgesvdq_obj->n, sgesvdq_obj->aref, sgesvdq_obj->lda, sgesvdq_obj->sref,\
								sgesvdq_obj->uref, sgesvdq_obj->ldu, sgesvdq_obj->vref, sgesvdq_obj->ldv, sgesvdq_obj->numrankref); 

    /* Compute libflame's Lapacke o/p  */
    sgesvdq_obj->info = LAPACKE_sgesvdq(sgesvdq_obj->matrix_layout, sgesvdq_obj->joba, sgesvdq_obj->jobp, sgesvdq_obj->jobr, sgesvdq_obj->jobu,\
								sgesvdq_obj->jobv, sgesvdq_obj->m, sgesvdq_obj->n, sgesvdq_obj->a, sgesvdq_obj->lda, sgesvdq_obj->s,\
								sgesvdq_obj->u, sgesvdq_obj->ldu, sgesvdq_obj->v, sgesvdq_obj->ldv, sgesvdq_obj->numrank);

    if( sgesvdq_obj->info < 0 ) {
        printf( "\n warning: The i:%x22 th argument with libflame \
        LAPACKE_sgesvdq is wrong\n", sgesvdq_obj->info );
    }
    if( sgesvdq_obj->inforef < 0 ) {
        printf( "The i:%x22 th argument with Netlib LAPACKE_sgesvdq is wrong\n", 
        sgesvdq_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    sgesvdq_obj->diff_a =  computeDiff_s( sgesvdq_obj->bufsize_a, sgesvdq_obj->a, sgesvdq_obj->aref );
	sgesvdq_obj->diff_u =  computeDiff_s( sgesvdq_obj->bufsize_u, sgesvdq_obj->u, sgesvdq_obj->uref );
	sgesvdq_obj->diff_v =  computeDiff_s( sgesvdq_obj->bufsize_v, sgesvdq_obj->v, sgesvdq_obj->vref );
	sgesvdq_obj->diff_s =  computeDiff_s( sgesvdq_obj->bufsize_s, sgesvdq_obj->s, sgesvdq_obj->sref );

}

TEST_F(sgesvdq_test, sgesvdq1) {
    EXPECT_NEAR(0.0, sgesvdq_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, sgesvdq_obj->diff_u, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, sgesvdq_obj->diff_v,  LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, sgesvdq_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
	
}

TEST_F(sgesvdq_test, sgesvdq2) {
    EXPECT_NEAR(0.0, sgesvdq_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, sgesvdq_obj->diff_u, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, sgesvdq_obj->diff_v,  LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, sgesvdq_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgesvdq_test, sgesvdq3) {
    EXPECT_NEAR(0.0, sgesvdq_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, sgesvdq_obj->diff_u, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, sgesvdq_obj->diff_v,  LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, sgesvdq_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgesvdq_test, sgesvdq4) {
    EXPECT_NEAR(0.0, sgesvdq_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, sgesvdq_obj->diff_u, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, sgesvdq_obj->diff_v,  LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, sgesvdq_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
}


/* Begin scomplex_common_parameters  class definition */
class gesvdq_scomplex_parameters{

   public:
	int bufsize_a, bufsize_u, bufsize_s, bufsize_v;
	void *hModule, *dModule;
	float diff_a, diff_u, diff_s, diff_v;
   /*input parameters */
	int matrix_layout;
	lapack_int m;
	lapack_int n;
	lapack_complex_float* a, *v;
	char joba, jobp, jobr, jobu, jobv;
	lapack_int lda, ldu, ldv;
	/*Output Parameter*/
	float* s;
	float* sref;
	lapack_int *numrank, *numrankref;
	lapack_complex_float *u;
	lapack_complex_float *aref, *uref, *vref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      gesvdq_scomplex_parameters (int matrix_layout, char joba, char jobp, char jobr, char jobu, char jobv,  lapack_int m, lapack_int n);
      ~gesvdq_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
gesvdq_scomplex_parameters:: gesvdq_scomplex_parameters (int matrix_layout_i, char joba_i, char jobp_i, char jobr_i, char jobu_i, char jobv_i,  lapack_int m_i, lapack_int n_i)
{
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	joba = joba_i; //A, H, M, E
	jobp = jobp_i; // P, N
	jobr = jobr_i; // T, N
	jobu = jobu_i; //A, F, S, R, N
	jobv = jobv_i; // A, V, R, N
	
	#if LAPACKE_TEST_VERBOSE
		printf(" \n gesvdq scomplex: matrix_layout = %d, joba:%c, jobp:%c, jobr:%c, jobu:%c, jobv:%c, m:%d, n:%d\n", joba, jobp, jobr, jobu, jobv, m, n);	
	#endif
	
	if ( jobr == 'C')
		jobr = 'N';
		
	if (jobp == 'U')
		jobp == 'P';
	
	if (jobu == 'V')
		jobu = 'S';
	
	/* Buffer sizes */
	/*buffer A*/
	lda = m;
	bufsize_a = lda *n;
	/*Buffer s */
	bufsize_s = n;
	/*buffer U*/
	if ((jobu == 'A') ||(jobu  == 'S') ||(jobu == 'U')  ||( jobu == 'R'))
		ldu = m;
	else if (jobu == 'F')
		ldu = n;
	/*Buffer V*/
	ldv = n;
	bufsize_v = ldv *n;	
		
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&a, &aref, bufsize_a);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&u, &uref, bufsize_u);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&v, &vref, bufsize_v);
	lapacke_gtest_alloc_float_buffer_pair(&s, &sref, bufsize_s);
	lapacke_gtest_alloc_int_buffer_pair(&numrank, &numrankref, n);

	if ((a==NULL) || (aref==NULL) || \
		(u==NULL) || (uref==NULL) || \
		(v==NULL) || (vref==NULL) ||\
		(s==NULL) || (sref==NULL) ||\
		(numrank==NULL) || (numrankref==NULL)){
		EXPECT_FALSE( true) << "gesvdq_float_parameters object: malloc error.";
		gesvdq_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, bufsize_a);
	lapacke_gtest_init_int_buffer_pair_with_constant(numrank, numrankref, n, n);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
gesvdq_scomplex_parameters :: ~gesvdq_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" gesvdq_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   gesvdq_free();

}
/*  Test fixture class definition */
class cgesvdq_test  : public  ::testing::Test {
public:
   gesvdq_scomplex_parameters  *cgesvdq_obj;
   void SetUp();
   void TearDown () { delete cgesvdq_obj;}
};

void cgesvdq_test::SetUp(){

	/* LAPACKE cgesvdq prototype */
    typedef int (*Fptr_NL_LAPACKE_cgesvdq) ( int matrix_layout, char joba, char jobp, char jobr, char jobu, char jobv,\
                           lapack_int m, lapack_int n, lapack_complex_float* a,\
                           lapack_int lda, float* s, lapack_complex_float* u,\
                           lapack_int ldu, lapack_complex_float* v,\
                           lapack_int ldv, lapack_int* numrank);

    Fptr_NL_LAPACKE_cgesvdq cgesvdq;

    cgesvdq_obj = new gesvdq_scomplex_parameters ( svd_paramslist[idx].matrix_layout,
                           svd_paramslist[idx].joba_gesvdq,
                           svd_paramslist[idx].jobp_gejsv,
						   svd_paramslist[idx].jobt_gejsv,
						   svd_paramslist[idx].jobu_gesvdq,
						   svd_paramslist[idx].jobv_gesvdq,
						   svd_paramslist[idx].m_gejsv,
						   svd_paramslist[idx].n_gejsv);
						   

    idx = Circular_Increment_Index(idx);

    cgesvdq_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgesvdq_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgesvdq_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgesvdq_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*cgesvdq library call */
    cgesvdq = (Fptr_NL_LAPACKE_cgesvdq)dlsym(cgesvdq_obj->hModule, "LAPACKE_cgesvdq");
    ASSERT_TRUE(cgesvdq != NULL) << "failed to get the Netlib LAPACKE_cgesvdq symbol";
	
/*Compute cgesvdq's  o/p */
    cgesvdq_obj->inforef = cgesvdq( cgesvdq_obj->matrix_layout, cgesvdq_obj->joba, cgesvdq_obj->jobp, cgesvdq_obj->jobr, cgesvdq_obj->jobu,\
								cgesvdq_obj->jobv, cgesvdq_obj->m, cgesvdq_obj->n, cgesvdq_obj->aref, cgesvdq_obj->lda, cgesvdq_obj->sref,\
								cgesvdq_obj->uref, cgesvdq_obj->ldu, cgesvdq_obj->vref, cgesvdq_obj->ldv, cgesvdq_obj->numrankref); 

    /* Compute libflame's Lapacke o/p  */
    cgesvdq_obj->info = LAPACKE_cgesvdq(cgesvdq_obj->matrix_layout, cgesvdq_obj->joba, cgesvdq_obj->jobp, cgesvdq_obj->jobr, cgesvdq_obj->jobu,\
								cgesvdq_obj->jobv, cgesvdq_obj->m, cgesvdq_obj->n, cgesvdq_obj->a, cgesvdq_obj->lda, cgesvdq_obj->s,\
								cgesvdq_obj->u, cgesvdq_obj->ldu, cgesvdq_obj->v, cgesvdq_obj->ldv, cgesvdq_obj->numrank);

    if( cgesvdq_obj->info < 0 ) {
        printf( "\n warning: The i:%x22 th argument with libflame \
        LAPACKE_cgesvdq is wrong\n", cgesvdq_obj->info );
    }
    if( cgesvdq_obj->inforef < 0 ) {
        printf( "The i:%x22 th argument with Netlib LAPACKE_cgesvdq is wrong\n", 
        cgesvdq_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cgesvdq_obj->diff_a =  computeDiff_c( cgesvdq_obj->bufsize_a, cgesvdq_obj->a, cgesvdq_obj->aref );
	cgesvdq_obj->diff_u =  computeDiff_c( cgesvdq_obj->bufsize_u, cgesvdq_obj->u, cgesvdq_obj->uref );
	cgesvdq_obj->diff_v =  computeDiff_c( cgesvdq_obj->bufsize_v, cgesvdq_obj->v, cgesvdq_obj->vref );
	cgesvdq_obj->diff_s =  computeDiff_s( cgesvdq_obj->bufsize_s, cgesvdq_obj->s, cgesvdq_obj->sref );

}
/*
TEST_F(cgesvdq_test, cgesvdq1) {
    EXPECT_NEAR(0.0, cgesvdq_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, cgesvdq_obj->diff_u, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, cgesvdq_obj->diff_v,  LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, cgesvdq_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
	
}

TEST_F(cgesvdq_test, cgesvdq2) {
    EXPECT_NEAR(0.0, cgesvdq_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, cgesvdq_obj->diff_u, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, cgesvdq_obj->diff_v,  LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, cgesvdq_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgesvdq_test, cgesvdq3) {
    EXPECT_NEAR(0.0, cgesvdq_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, cgesvdq_obj->diff_u, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, cgesvdq_obj->diff_v,  LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, cgesvdq_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgesvdq_test, cgesvdq4) {
    EXPECT_NEAR(0.0, cgesvdq_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, cgesvdq_obj->diff_u, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, cgesvdq_obj->diff_v,  LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, cgesvdq_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
}*/

/* Begin dcomplex_common_parameters  class definition */

class gesvdq_dcomplex_parameters{

   public:
	int bufsize_a, bufsize_u, bufsize_s, bufsize_v;
	void *hModule, *dModule;
	double diff_a, diff_u, diff_s, diff_v;
   /*input parameters */
	int matrix_layout;
	lapack_int m;
	lapack_int n;
	lapack_complex_double* a, *v;
	char joba, jobp, jobr, jobu, jobv;
	lapack_int lda, ldu, ldv;
	/*Output Parameter*/
	double* s;
	double* sref;
	lapack_int *numrank, *numrankref;
	lapack_complex_double *u;
	lapack_complex_double *aref, *uref, *vref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      gesvdq_dcomplex_parameters (int matrix_layout, char joba, char jobp, char jobr, char jobu, char jobv,  lapack_int m, lapack_int n);
      ~gesvdq_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
gesvdq_dcomplex_parameters:: gesvdq_dcomplex_parameters (int matrix_layout_i, char joba_i, char jobp_i, char jobr_i, char jobu_i, char jobv_i,  lapack_int m_i, lapack_int n_i)
{
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	joba = joba_i; //A, H, M, E
	jobp = jobp_i; // P, N
	jobr = jobr_i; // T, N
	jobu = jobu_i; //A, F, S, R, N
	jobv = jobv_i; // A, V, R, N
	
	#if LAPACKE_TEST_VERBOSE
		printf(" \n gesvdq dcomplex: matrix_layout = %d, joba:%c, jobp:%c, jobr:%c, jobu:%c, jobv:%c, m:%d, n:%d\n", joba, jobp, jobr, jobu, jobv, m, n);	
	#endif
	
	if ( jobr == 'C')
		jobr = 'N';
		
	if (jobp == 'U')
		jobp == 'P';
	
	if (jobu == 'V')
		jobu = 'S';
	
	/* Buffer sizes */
	/*buffer A*/
	lda = m;
	bufsize_a = lda *n;
	/*Buffer s */
	bufsize_s = n;
	/*buffer U*/
	if ((jobu == 'A') ||(jobu  == 'S') ||(jobu == 'U')  ||( jobu == 'R'))
		ldu = m;
	else if (jobu == 'F')
		ldu = n;
	/*Buffer V*/
	ldv = n;
	bufsize_v = ldv *n;	
		
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&a, &aref, bufsize_a);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&u, &uref, bufsize_u);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&v, &vref, bufsize_v);
	lapacke_gtest_alloc_double_buffer_pair(&s, &sref, bufsize_s);
	lapacke_gtest_alloc_int_buffer_pair(&numrank, &numrankref, n);

	if ((a==NULL) || (aref==NULL) || \
		(u==NULL) || (uref==NULL) || \
		(v==NULL) || (vref==NULL) ||\
		(s==NULL) || (sref==NULL) ||\
		(numrank==NULL) || (numrankref==NULL)){
		EXPECT_FALSE( true) << "gesvdq_double_parameters object: malloc error.";
		gesvdq_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, bufsize_a);
	lapacke_gtest_init_int_buffer_pair_with_constant(numrank, numrankref, n, n);

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
gesvdq_dcomplex_parameters :: ~gesvdq_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" gesvdq_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   gesvdq_free();

}
/*  Test fixture class definition */
class zgesvdq_test  : public  ::testing::Test {
public:
   gesvdq_dcomplex_parameters  *zgesvdq_obj;
   void SetUp();
   void TearDown () { delete zgesvdq_obj;}
};

void zgesvdq_test::SetUp(){

	/* LAPACKE zgesvdq prototype */
    typedef int (*Fptr_NL_LAPACKE_zgesvdq) ( int matrix_layout, char joba, char jobp, char jobr, char jobu, char jobv,\
                           lapack_int m, lapack_int n, lapack_complex_double* a,\
                           lapack_int lda, double* s, lapack_complex_double* u,\
                           lapack_int ldu, lapack_complex_double* v,\
                           lapack_int ldv, lapack_int* numrank);

    Fptr_NL_LAPACKE_zgesvdq zgesvdq;

    zgesvdq_obj = new gesvdq_dcomplex_parameters ( svd_paramslist[idx].matrix_layout,
                           svd_paramslist[idx].joba_gesvdq,
                           svd_paramslist[idx].jobp_gejsv,
						   svd_paramslist[idx].jobt_gejsv,
						   svd_paramslist[idx].jobu_gesvdq,
						   svd_paramslist[idx].jobv_gesvdq,
						   svd_paramslist[idx].m_gejsv,
						   svd_paramslist[idx].n_gejsv);

    idx = Circular_Increment_Index(idx);

    zgesvdq_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgesvdq_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgesvdq_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgesvdq_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*zgesvdq library call */
    zgesvdq = (Fptr_NL_LAPACKE_zgesvdq)dlsym(zgesvdq_obj->hModule, "LAPACKE_zgesvdq");
    ASSERT_TRUE(zgesvdq != NULL) << "failed to get the Netlib LAPACKE_zgesvdq symbol";
	
/*Compute zgesvdq's  o/p */
    zgesvdq_obj->inforef = zgesvdq( zgesvdq_obj->matrix_layout, zgesvdq_obj->joba, zgesvdq_obj->jobp, zgesvdq_obj->jobr, zgesvdq_obj->jobu,\
								zgesvdq_obj->jobv, zgesvdq_obj->m, zgesvdq_obj->n, zgesvdq_obj->aref, zgesvdq_obj->lda, zgesvdq_obj->sref,\
								zgesvdq_obj->uref, zgesvdq_obj->ldu, zgesvdq_obj->vref, zgesvdq_obj->ldv, zgesvdq_obj->numrankref); 

    /* Compute libflame's Lapacke o/p  */
    zgesvdq_obj->info = LAPACKE_zgesvdq(zgesvdq_obj->matrix_layout, zgesvdq_obj->joba, zgesvdq_obj->jobp, zgesvdq_obj->jobr, zgesvdq_obj->jobu,\
								zgesvdq_obj->jobv, zgesvdq_obj->m, zgesvdq_obj->n, zgesvdq_obj->a, zgesvdq_obj->lda, zgesvdq_obj->s,\
								zgesvdq_obj->u, zgesvdq_obj->ldu, zgesvdq_obj->v, zgesvdq_obj->ldv, zgesvdq_obj->numrank);

    if( zgesvdq_obj->info < 0 ) {
        printf( "\n warning: The i:%x22 th argument with libflame \
        LAPACKE_zgesvdq is wrong\n", zgesvdq_obj->info );
    }
    if( zgesvdq_obj->inforef < 0 ) {
        printf( "The i:%x22 th argument with Netlib LAPACKE_zgesvdq is wrong\n", 
        zgesvdq_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zgesvdq_obj->diff_a =  computeDiff_z( zgesvdq_obj->bufsize_a, zgesvdq_obj->a, zgesvdq_obj->aref );
	zgesvdq_obj->diff_u =  computeDiff_z( zgesvdq_obj->bufsize_u, zgesvdq_obj->u, zgesvdq_obj->uref );
	zgesvdq_obj->diff_v =  computeDiff_z( zgesvdq_obj->bufsize_v, zgesvdq_obj->v, zgesvdq_obj->vref );
	zgesvdq_obj->diff_s =  computeDiff_d( zgesvdq_obj->bufsize_s, zgesvdq_obj->s, zgesvdq_obj->sref );

}
/*
TEST_F(zgesvdq_test, zgesvdq1) {
    EXPECT_NEAR(0.0, zgesvdq_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zgesvdq_obj->diff_u, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zgesvdq_obj->diff_v,  LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zgesvdq_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
	
}

TEST_F(zgesvdq_test, zgesvdq2) {
    EXPECT_NEAR(0.0, zgesvdq_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zgesvdq_obj->diff_u, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zgesvdq_obj->diff_v,  LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zgesvdq_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgesvdq_test, zgesvdq3) {
    EXPECT_NEAR(0.0, zgesvdq_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zgesvdq_obj->diff_u, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zgesvdq_obj->diff_v,  LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zgesvdq_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgesvdq_test, zgesvdq4) {
    EXPECT_NEAR(0.0, zgesvdq_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zgesvdq_obj->diff_u, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zgesvdq_obj->diff_v,  LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zgesvdq_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
}*/