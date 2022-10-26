#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"


#define trexc_free() \
if (q!=NULL)    free(q); \
if (qref!=NULL) free(qref);\
if (t!=NULL)  free(t); \
if (tref!=NULL)  free(tref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule);
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
/* Begin float_common_parameters  class definition */
class trexc_float_parameters{

   public:
	int bufsize_q;
	int bufsize_t;
	void *hModule, *dModule;
	float diff_t;
	float diff_q;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	char compq;
	float* q;
	float* t;
	lapack_int ifst,ilst;
	lapack_int ldq;
	lapack_int ldt;
	/*Output Parameter*/
	float *qref,*tref;	
	/*Return Values*/
	lapack_int info, inforef;

   public:
      trexc_float_parameters (int matrix_layout, char compq, lapack_int n, lapack_int ldt, lapack_int ldq, lapack_int ifst, lapack_int ilst );
      ~trexc_float_parameters ();

};

/* Constructor definition  float_common_parameters */
trexc_float_parameters:: trexc_float_parameters (int matrix_layout_i, char compq_i, lapack_int n_i, lapack_int ldt_i, lapack_int ldq_i, lapack_int ifst_i, lapack_int ilst_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	compq = compq_i;
	ldt = ldt_i;
	ldq = ldq_i;
	ifst = ifst_i;
	ilst = ilst_i;
	
	//#if LAPACKE_TEST_VERBOSE
	printf(" \n trexc float:  n: %d, ldt: %d ldq: %d, compq = %c \n",  n, ldt, ldq, compq);
//	#endif
	if (compq == 'U')
		compq = 'N';
	else if (compq == 'V')
		ldq = n;
	
	bufsize_q = (ldq*n);
	bufsize_t = (ldt*n);
	
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&q, &qref, bufsize_q);	
	lapacke_gtest_alloc_float_buffer_pair(&t, &tref, bufsize_t);
	
	if ((q==NULL) || (qref==NULL) || 
		(t == NULL) || (tref == NULL)){
		EXPECT_FALSE( true) << "trexc_float_parameters object: malloc error.";
		trexc_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand( q, qref, bufsize_q);
	lapacke_gtest_init_float_buffer_pair_rand( t, tref, bufsize_t);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
trexc_float_parameters :: ~trexc_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" trexc_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   trexc_free();

}
/*  Test fixture class definition */
class strexc_test  : public  ::testing::Test {
public:
   trexc_float_parameters  *strexc_obj;
   void SetUp();
   void TearDown () { delete strexc_obj; }
};

void strexc_test::SetUp(){

    /* LAPACKE strexc prototype */
    typedef int (*Fptr_NL_LAPACKE_strexc) (int matrix_layout, char compq, lapack_int n, float* t, lapack_int ldt, float* q, lapack_int ldq, lapack_int* ifst, lapack_int* ilst );

    Fptr_NL_LAPACKE_strexc strexc;

    strexc_obj = new trexc_float_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].vect_rd,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].n,
						   eig_paramslist[idx].m,
						   eig_paramslist[idx].k,
						   eig_paramslist[idx].kb);

    idx = Circular_Increment_Index(idx);
	

    strexc_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    strexc_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(strexc_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(strexc_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    strexc = (Fptr_NL_LAPACKE_strexc)dlsym(strexc_obj->hModule, "LAPACKE_strexc");
    ASSERT_TRUE(strexc != NULL) << "failed to get the Netlib LAPACKE_strexc symbol";
    

    strexc_obj->inforef = strexc( strexc_obj->matrix_layout,strexc_obj->compq,  strexc_obj->n,
								strexc_obj->tref, strexc_obj->ldt, strexc_obj->qref, strexc_obj->ldq, 
								&(strexc_obj->ifst), &(strexc_obj->ilst));

    /* Compute libflame's Lapacke o/p  */
    strexc_obj->info = LAPACKE_strexc( strexc_obj->matrix_layout, strexc_obj->compq,  strexc_obj->n,
								strexc_obj->t, strexc_obj->ldt, strexc_obj->q, strexc_obj->ldq,
								&(strexc_obj->ifst), &(strexc_obj->ilst));

    if( strexc_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_strexc is wrong\n", strexc_obj->info );
    }
    if( strexc_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_strexc is wrong\n", 
        strexc_obj->inforef );
    }

    /* Compute diff_qerence between libflame and Netlib o/ps  */
    strexc_obj->diff_t =  computeDiff_s( strexc_obj->bufsize_t, 
                strexc_obj->t, strexc_obj->tref );

}

TEST_F(strexc_test, strexc1) {
    EXPECT_NEAR(0.0, strexc_obj->diff_t, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(strexc_test, strexc2) {
    EXPECT_NEAR(0.0, strexc_obj->diff_t, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(strexc_test, strexc3) {
    EXPECT_NEAR(0.0, strexc_obj->diff_t, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(strexc_test, strexc4) {
    EXPECT_NEAR(0.0, strexc_obj->diff_t, LAPACKE_GTEST_THRESHOLD);
}


/* Begin double_common_parameters  class definition */
class trexc_double_parameters{

   public:
	int bufsize_q;
	int bufsize_t;
	void *hModule, *dModule;
	double diff_t;
	double diff_q;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	char compq;
	double* q;
	double* t;
	lapack_int  ifst, ilst;
	lapack_int ldq;
	lapack_int ldt;
	/*Output Parameter*/
	double *qref,*tref;	
	/*Return Values*/
	lapack_int info, inforef;

   public:
      trexc_double_parameters (int matrix_layout, char compq, lapack_int n, lapack_int ldt, lapack_int ldq, lapack_int ifst, lapack_int ilst);
      ~trexc_double_parameters ();

};

/* Constructor definition  double_common_parameters */
trexc_double_parameters:: trexc_double_parameters (int matrix_layout_i, char compq_i, lapack_int n_i, lapack_int ldt_i, lapack_int ldq_i, lapack_int ifst_i, lapack_int ilst_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	compq = compq_i;
	ldt = ldt_i;
	ldq = ldq_i;
	ifst = ifst_i;
	ilst = ilst_i;
	
	#if LAPACKE_TEST_VERBOSE
	printf(" \n trexc double:  n: %d, ldt: %d ldq: %d, compq = %c \n",  n, ldt, ldq, compq);
	#endif
	if (compq == 'U')
		compq = 'N';
	else if (compq == 'V')
		ldq = n;
	
	bufsize_q = (ldq*n);
	bufsize_t = (ldt*n);
	
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&q, &qref, bufsize_q);
	lapacke_gtest_alloc_double_buffer_pair(&t, &tref, bufsize_t);
	
	if ((q==NULL) || (qref==NULL) || 
		(t == NULL) || (tref == NULL)){
		EXPECT_FALSE( true) << "trexc_double_parameters object: malloc error.";
		trexc_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand( q, qref, bufsize_q);
	lapacke_gtest_init_double_buffer_pair_rand( t, tref, bufsize_t);

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
trexc_double_parameters :: ~trexc_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" trexc_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   trexc_free();

}
/*  Test fixture class definition */
class dtrexc_test  : public  ::testing::Test {
public:
   trexc_double_parameters  *dtrexc_obj;
   void SetUp();
   void TearDown () { delete dtrexc_obj; }
};

void dtrexc_test::SetUp(){

    /* LAPACKE dtrexc prototype */
    typedef int (*Fptr_NL_LAPACKE_dtrexc) (int matrix_layout, char compq, lapack_int n, double* t, lapack_int ldt, double* q, lapack_int ldq, lapack_int* ifst, lapack_int* ilst );

    Fptr_NL_LAPACKE_dtrexc dtrexc;

    dtrexc_obj = new trexc_double_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].vect_rd,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].n,
						   eig_paramslist[idx].m, 
						   eig_paramslist[idx].k,
						   eig_paramslist[idx].kb);

    idx = Circular_Increment_Index(idx);
	

    dtrexc_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dtrexc_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dtrexc_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dtrexc_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    dtrexc = (Fptr_NL_LAPACKE_dtrexc)dlsym(dtrexc_obj->hModule, "LAPACKE_dtrexc");
    ASSERT_TRUE(dtrexc != NULL) << "failed to get the Netlib LAPACKE_dtrexc symbol";
    

    dtrexc_obj->inforef = dtrexc( dtrexc_obj->matrix_layout,dtrexc_obj->compq,  dtrexc_obj->n,
								dtrexc_obj->tref, dtrexc_obj->ldt, dtrexc_obj->qref, dtrexc_obj->ldq, 
								&(dtrexc_obj->ifst), &(dtrexc_obj->ilst));

    /* Compute libflame's Lapacke o/p  */
    dtrexc_obj->info = LAPACKE_dtrexc( dtrexc_obj->matrix_layout,dtrexc_obj->compq,  dtrexc_obj->n,
								dtrexc_obj->t, dtrexc_obj->ldt, dtrexc_obj->q, dtrexc_obj->ldq,
								&(dtrexc_obj->ifst), &(dtrexc_obj->ilst));

    if( dtrexc_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dtrexc is wrong\n", dtrexc_obj->info );
    }
    if( dtrexc_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dtrexc is wrong\n", 
        dtrexc_obj->inforef );
    }

    /* Compute diff_qerence between libflame and Netlib o/ps  */
    dtrexc_obj->diff_t =  computeDiff_d( dtrexc_obj->bufsize_t, 
                dtrexc_obj->t, dtrexc_obj->tref );

}

TEST_F(dtrexc_test, dtrexc1) {
    EXPECT_NEAR(0.0, dtrexc_obj->diff_t, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dtrexc_test, dtrexc2) {
    EXPECT_NEAR(0.0, dtrexc_obj->diff_t, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dtrexc_test, dtrexc3) {
    EXPECT_NEAR(0.0, dtrexc_obj->diff_t, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dtrexc_test, dtrexc4) {
    EXPECT_NEAR(0.0, dtrexc_obj->diff_t, LAPACKE_GTEST_THRESHOLD);
}

/* Begin scomplex_common_parameters  class definition */

class trexc_scomplex_parameters{

   public:
	int bufsize_q;
	int bufsize_t;
	void *hModule, *dModule;
	float diff_t;
	float diff_q;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	char compq;
	lapack_complex_float* q;
	lapack_complex_float* t;
	lapack_int ifst,ilst;
	lapack_int ldq;
	lapack_int ldt;
	/*Output Parameter*/
	lapack_complex_float *qref,*tref;	
	/*Return Values*/
	lapack_int info, inforef;

   public:
      trexc_scomplex_parameters (int matrix_layout, char compq, lapack_int n, lapack_int ldt, lapack_int ldq, lapack_int ifst, lapack_int ilst );
      ~trexc_scomplex_parameters ();

};

/* Constructor definition  scomplex_common_parameters */
trexc_scomplex_parameters:: trexc_scomplex_parameters (int matrix_layout_i, char compq_i, lapack_int n_i, lapack_int ldt_i, lapack_int ldq_i, lapack_int ifst_i, lapack_int ilst_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	compq = compq_i;
	ldt = ldt_i;
	ldq = ldq_i;
	ifst = ifst_i;
	ilst = ilst_i;
	
	#if LAPACKE_TEST_VERBOSE
	printf(" \n trexc scomplex:  n: %d, ldt: %d ldq: %d, compq = %c \n",  n, ldt, ldq, compq);
	#endif
	if (compq == 'U')
		compq = 'N';
	else if (compq == 'V')
		ldq = n;
	
	bufsize_q = (ldq*n);
	bufsize_t = (ldt*n);
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&q, &qref, bufsize_q);	
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&t, &tref, bufsize_t);
	
	if ((q==NULL) || (qref==NULL) || 
		(t == NULL) || (tref == NULL)){
		EXPECT_FALSE( true) << "trexc_scomplex_parameters object: malloc error.";
		trexc_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( q, qref, bufsize_q);
	lapacke_gtest_init_scomplex_buffer_pair_rand( t, tref, bufsize_t);

} /* end of Constructor  */

/* Destructor definition  'scomplex_common_parameters' class  */
trexc_scomplex_parameters :: ~trexc_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" trexc_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   trexc_free();

}
/*  Test fixture class definition */
class ctrexc_test  : public  ::testing::Test {
public:
   trexc_scomplex_parameters  *ctrexc_obj;
   void SetUp();
   void TearDown () { delete ctrexc_obj; }
};

void ctrexc_test::SetUp(){

    /* LAPACKE ctrexc prototype */
    typedef int (*Fptr_NL_LAPACKE_ctrexc) (int matrix_layout, char compq, lapack_int n, lapack_complex_float* t, lapack_int ldt, lapack_complex_float* q, lapack_int ldq, lapack_int ifst, lapack_int ilst );

    Fptr_NL_LAPACKE_ctrexc ctrexc;

    ctrexc_obj = new trexc_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].vect_rd,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].n,
						   eig_paramslist[idx].m,
						   eig_paramslist[idx].k,
						   eig_paramslist[idx].kb);

    idx = Circular_Increment_Index(idx);
	

    ctrexc_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ctrexc_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ctrexc_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ctrexc_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    ctrexc = (Fptr_NL_LAPACKE_ctrexc)dlsym(ctrexc_obj->hModule, "LAPACKE_ctrexc");
    ASSERT_TRUE(ctrexc != NULL) << "failed to get the Netlib LAPACKE_ctrexc symbol";
    

    ctrexc_obj->inforef = ctrexc( ctrexc_obj->matrix_layout,ctrexc_obj->compq,  ctrexc_obj->n,
								ctrexc_obj->tref, ctrexc_obj->ldt, ctrexc_obj->qref, ctrexc_obj->ldq, 
								ctrexc_obj->ifst, ctrexc_obj->ilst);

    /* Compute libflame's Lapacke o/p  */
    ctrexc_obj->info = LAPACKE_ctrexc( ctrexc_obj->matrix_layout, ctrexc_obj->compq,  ctrexc_obj->n,
								ctrexc_obj->t, ctrexc_obj->ldt, ctrexc_obj->q, ctrexc_obj->ldq,
								ctrexc_obj->ifst, ctrexc_obj->ilst);

    if( ctrexc_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ctrexc is wrong\n", ctrexc_obj->info );
    }
    if( ctrexc_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ctrexc is wrong\n", 
        ctrexc_obj->inforef );
    }

    /* Compute diff_qerence between libflame and Netlib o/ps  */
    ctrexc_obj->diff_t =  computeDiff_c( ctrexc_obj->bufsize_t, 
                ctrexc_obj->t, ctrexc_obj->tref );

}

TEST_F(ctrexc_test, ctrexc1) {
    EXPECT_NEAR(0.0, ctrexc_obj->diff_t, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ctrexc_test, ctrexc2) {
    EXPECT_NEAR(0.0, ctrexc_obj->diff_t, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ctrexc_test, ctrexc3) {
    EXPECT_NEAR(0.0, ctrexc_obj->diff_t, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ctrexc_test, ctrexc4) {
    EXPECT_NEAR(0.0, ctrexc_obj->diff_t, LAPACKE_GTEST_THRESHOLD);
}


/* Begin dcomplex_common_parameters  class definition */
class trexc_dcomplex_parameters{

   public:
	int bufsize_q;
	int bufsize_t;
	void *hModule, *dModule;
	double diff_t;
	double diff_q;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	char compq;
	lapack_complex_double* q;
	lapack_complex_double* t;
	lapack_int ifst,ilst;
	lapack_int ldq;
	lapack_int ldt;
	/*Output Parameter*/
	lapack_complex_double *qref,*tref;	
	/*Return Values*/
	lapack_int info, inforef;

   public:
      trexc_dcomplex_parameters (int matrix_layout, char compq, lapack_int n, lapack_int ldt, lapack_int ldq, lapack_int ifst, lapack_int ilst );
      ~trexc_dcomplex_parameters ();

};

/* Constructor definition  dcomplex_common_parameters */
trexc_dcomplex_parameters:: trexc_dcomplex_parameters (int matrix_layout_i, char compq_i, lapack_int n_i, lapack_int ldt_i, lapack_int ldq_i, lapack_int ifst_i, lapack_int ilst_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	compq = compq_i;
	ldt = ldt_i;
	ldq = ldq_i;
	ifst = ifst_i;
	ilst = ilst_i;
	
	#if LAPACKE_TEST_VERBOSE
	printf(" \n trexc dcomplex:  n: %d, ldt: %d ldq: %d, compq = %c \n",  n, ldt, ldq, compq);
	#endif
	if (compq == 'U')
		compq = 'N';
	else if (compq == 'V')
		ldq = n;
	
	bufsize_q = (ldq*n);
	bufsize_t = (ldt*n);
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&q, &qref, bufsize_q);	
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&t, &tref, bufsize_t);
	
	if ((q==NULL) || (qref==NULL) || 
		(t == NULL) || (tref == NULL)){
		EXPECT_FALSE( true) << "trexc_dcomplex_parameters object: malloc error.";
		trexc_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( q, qref, bufsize_q);
	lapacke_gtest_init_dcomplex_buffer_pair_rand( t, tref, bufsize_t);

} /* end of Constructor  */

/* Destructor definition  'dcomplex_common_parameters' class  */
trexc_dcomplex_parameters :: ~trexc_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" trexc_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   trexc_free();

}
/*  Test fixture class definition */
class ztrexc_test  : public  ::testing::Test {
public:
   trexc_dcomplex_parameters  *ztrexc_obj;
   void SetUp();
   void TearDown () { delete ztrexc_obj; }
};

void ztrexc_test::SetUp(){

    /* LAPACKE ztrexc prototype */
    typedef int (*Fptr_NL_LAPACKE_ztrexc) (int matrix_layout, char compq, lapack_int n, lapack_complex_double* t, lapack_int ldt, lapack_complex_double* q, lapack_int ldq, lapack_int ifst, lapack_int ilst );

    Fptr_NL_LAPACKE_ztrexc ztrexc;

    ztrexc_obj = new trexc_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].vect_rd,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].n,
						   eig_paramslist[idx].m,
						   eig_paramslist[idx].k,
						   eig_paramslist[idx].kb);

    idx = Circular_Increment_Index(idx);
	

    ztrexc_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ztrexc_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ztrexc_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ztrexc_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    ztrexc = (Fptr_NL_LAPACKE_ztrexc)dlsym(ztrexc_obj->hModule, "LAPACKE_ztrexc");
    ASSERT_TRUE(ztrexc != NULL) << "failed to get the Netlib LAPACKE_ztrexc symbol";
    

    ztrexc_obj->inforef = ztrexc( ztrexc_obj->matrix_layout,ztrexc_obj->compq,  ztrexc_obj->n,
								ztrexc_obj->tref, ztrexc_obj->ldt, ztrexc_obj->qref, ztrexc_obj->ldq, 
								ztrexc_obj->ifst, ztrexc_obj->ilst);

    /* Compute libflame's Lapacke o/p  */
    ztrexc_obj->info = LAPACKE_ztrexc( ztrexc_obj->matrix_layout, ztrexc_obj->compq,  ztrexc_obj->n,
								ztrexc_obj->t, ztrexc_obj->ldt, ztrexc_obj->q, ztrexc_obj->ldq,
								ztrexc_obj->ifst, ztrexc_obj->ilst);

    if( ztrexc_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ztrexc is wrong\n", ztrexc_obj->info );
    }
    if( ztrexc_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ztrexc is wrong\n", 
        ztrexc_obj->inforef );
    }

    /* Compute diff_qerence between libflame and Netlib o/ps  */
    ztrexc_obj->diff_t =  computeDiff_z( ztrexc_obj->bufsize_t, 
                ztrexc_obj->t, ztrexc_obj->tref );

}

TEST_F(ztrexc_test, ztrexc1) {
    EXPECT_NEAR(0.0, ztrexc_obj->diff_t, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ztrexc_test, ztrexc2) {
    EXPECT_NEAR(0.0, ztrexc_obj->diff_t, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ztrexc_test, ztrexc3) {
    EXPECT_NEAR(0.0, ztrexc_obj->diff_t, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ztrexc_test, ztrexc4) {
    EXPECT_NEAR(0.0, ztrexc_obj->diff_t, LAPACKE_GTEST_THRESHOLD);
}

