#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"


#define geqpf_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (tau!=NULL)  free(tau);\
if (tauref!=NULL) free(tauref); \
if (jpvt!=NULL)  free(jpvt);\
if (jpvtref!=NULL)  free(jpvtref);\
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
/* Begin float_common_parameters  class definition */
class geqpf_float_parameters{

   public:
	int bufsize;
	int jpvt_size;
	void *hModule, *dModule;
	float diff;
	float diff_jpvt;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	float* A;
	lapack_int *jpvt;
	lapack_int lda;
	/*Output Parameter*/
	float* tau;
	float *Aref, *tauref;
	lapack_int *jpvtref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      geqpf_float_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda);
      ~geqpf_float_parameters ();

};

/* Constructor definition  float_common_parameters */
geqpf_float_parameters:: geqpf_float_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n geqpf float:  m: %d, n: %d lda: %d \n",  m, n, lda);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{	lda = m;
		bufsize = lda*n*sizeof(float);		
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{	lda = n;
		bufsize = lda*m*sizeof(float);		
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	jpvt_size = n*sizeof(float);
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_int_buffer_pair(&jpvt, &jpvtref, jpvt_size);
	lapacke_gtest_alloc_float_buffer_pair(&tau, &tauref, (min(m,n)*sizeof(float)));	
	if ((A==NULL) || (Aref==NULL) || \
		(tau==NULL) || (tauref==NULL) || \
		(jpvt == NULL)|| (jpvtref == NULL)){
		EXPECT_FALSE( true) << "geqpf_float_parameters object: malloc error.";
		geqpf_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand( A, Aref, bufsize);
	lapacke_gtest_init_int_buffer_pair_rand( jpvt, jpvtref, jpvt_size);

	/*initialize output matrix by 0 */
	for(i=0;i<(min(m,n));i++) {
		tau[i] = 0;
		tauref[i] = tau[i];
	}

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
geqpf_float_parameters :: ~geqpf_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" geqpf_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   geqpf_free();

}
/*  Test fixture class definition */
class sgeqpf_test  : public  ::testing::Test {
public:
   geqpf_float_parameters  *sgeqpf_obj;
   void SetUp();
   void TearDown () { delete sgeqpf_obj; }
};

void sgeqpf_test::SetUp(){

    /* LAPACKE sgeqpf prototype */
    typedef int (*Fptr_NL_LAPACKE_sgeqpf) (int matrix_layout, lapack_int m,lapack_int n,\
											float *A, lapack_int lda, lapack_int *jpvt, float* tau);

    Fptr_NL_LAPACKE_sgeqpf sgeqpf;

    sgeqpf_obj = new geqpf_float_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);
	

    sgeqpf_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgeqpf_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgeqpf_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgeqpf_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    sgeqpf = (Fptr_NL_LAPACKE_sgeqpf)dlsym(sgeqpf_obj->hModule, "LAPACKE_sgeqpf");
    ASSERT_TRUE(sgeqpf != NULL) << "failed to get the Netlib LAPACKE_sgeqpf symbol";
    

    sgeqpf_obj->inforef = sgeqpf( sgeqpf_obj->matrix_layout, sgeqpf_obj->m,
								sgeqpf_obj->n,sgeqpf_obj->Aref,
								sgeqpf_obj->lda, sgeqpf_obj->jpvtref, sgeqpf_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    sgeqpf_obj->info = LAPACKE_sgeqpf( sgeqpf_obj->matrix_layout, sgeqpf_obj->m,
										sgeqpf_obj->n,sgeqpf_obj->A, 
										sgeqpf_obj->lda, sgeqpf_obj->jpvt, sgeqpf_obj->tau);

    if( sgeqpf_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_sgeqpf is wrong\n", sgeqpf_obj->info );
    }
    if( sgeqpf_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sgeqpf is wrong\n", 
        sgeqpf_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    sgeqpf_obj->diff =  computeDiff_s( sgeqpf_obj->bufsize, 
                sgeqpf_obj->A, sgeqpf_obj->Aref );
	sgeqpf_obj->diff_jpvt = computeDiff_i( sgeqpf_obj->jpvt_size,\
							sgeqpf_obj->jpvt, sgeqpf_obj->jpvtref); 

}

TEST_F(sgeqpf_test, sgeqpf1) {
    EXPECT_NEAR(0.0, sgeqpf_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, sgeqpf_obj->diff_jpvt, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgeqpf_test, sgeqpf2) {
    EXPECT_NEAR(0.0, sgeqpf_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, sgeqpf_obj->diff_jpvt, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgeqpf_test, sgeqpf3) {
    EXPECT_NEAR(0.0, sgeqpf_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, sgeqpf_obj->diff_jpvt, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgeqpf_test, sgeqpf4) {
    EXPECT_NEAR(0.0, sgeqpf_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, sgeqpf_obj->diff_jpvt, LAPACKE_GTEST_THRESHOLD);
}


/* Begin double_common_parameters  class definition */
class geqpf_double_parameters{
   public:
	int bufsize;
	int jpvt_size;
	void *hModule, *dModule;
	double diff;
	double diff_jpvt;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	double* A;
	lapack_int *jpvt;
	lapack_int lda;
	/*Output Parameter*/
	double* tau;
	double *Aref, *tauref;
	lapack_int *jpvtref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      geqpf_double_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda);
      ~geqpf_double_parameters ();

};

/* Constructor definition  double_common_parameters */
geqpf_double_parameters:: geqpf_double_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n geqpf double:  m: %d, n: %d lda: %d \n",  m, n, lda);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{	lda = m;
		bufsize = lda*n*sizeof(double);		
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{	lda = n;
		bufsize = lda*m*sizeof(double);		
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	jpvt_size = n*sizeof(double);
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_int_buffer_pair(&jpvt, &jpvtref, jpvt_size);
	lapacke_gtest_alloc_double_buffer_pair(&tau, &tauref, (min(m,n)*sizeof(double)));	
	if ((A==NULL) || (Aref==NULL) || \
		(tau==NULL) || (tauref==NULL) || \
		(jpvt == NULL)|| (jpvtref == NULL)){
		EXPECT_FALSE( true) << "geqpf_double_parameters object: malloc error.";
		geqpf_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand( A, Aref, bufsize);
	lapacke_gtest_init_int_buffer_pair_rand( jpvt, jpvtref, jpvt_size);
	/*initialize output matrix by 0 */
	for(i=0;i<(min(m,n));i++) {
		tau[i] = 0;
		tauref[i] = tau[i];
	}

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
geqpf_double_parameters :: ~geqpf_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" geqpf_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   geqpf_free();

}
/*  Test fixture class definition */
class dgeqpf_test  : public  ::testing::Test {
public:
   geqpf_double_parameters  *dgeqpf_obj;
   void SetUp();
   void TearDown () { delete dgeqpf_obj; }
};

void dgeqpf_test::SetUp(){

    /* LAPACKE dgeqpf prototype */
    typedef int (*Fptr_NL_LAPACKE_dgeqpf) (int matrix_layout, lapack_int m,lapack_int n,\
											double *A, lapack_int lda, lapack_int *jpvt, double* tau);

    Fptr_NL_LAPACKE_dgeqpf dgeqpf;

    dgeqpf_obj = new geqpf_double_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);
	

    dgeqpf_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgeqpf_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgeqpf_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgeqpf_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    dgeqpf = (Fptr_NL_LAPACKE_dgeqpf)dlsym(dgeqpf_obj->hModule, "LAPACKE_dgeqpf");
    ASSERT_TRUE(dgeqpf != NULL) << "failed to get the Netlib LAPACKE_dgeqpf symbol";
    

    dgeqpf_obj->inforef = dgeqpf( dgeqpf_obj->matrix_layout, dgeqpf_obj->m,
								dgeqpf_obj->n,dgeqpf_obj->Aref,
								dgeqpf_obj->lda, dgeqpf_obj->jpvtref, dgeqpf_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    dgeqpf_obj->info = LAPACKE_dgeqpf( dgeqpf_obj->matrix_layout, dgeqpf_obj->m,
										dgeqpf_obj->n,dgeqpf_obj->A, 
										dgeqpf_obj->lda, dgeqpf_obj->jpvt, dgeqpf_obj->tau);

    if( dgeqpf_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dgeqpf is wrong\n", dgeqpf_obj->info );
    }
    if( dgeqpf_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dgeqpf is wrong\n", 
        dgeqpf_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dgeqpf_obj->diff =  computeDiff_d( dgeqpf_obj->bufsize, 
                dgeqpf_obj->A, dgeqpf_obj->Aref );
	dgeqpf_obj->diff_jpvt = computeDiff_i( dgeqpf_obj->jpvt_size,\
							dgeqpf_obj->jpvt, dgeqpf_obj->jpvtref); 

}

TEST_F(dgeqpf_test, dgeqpf1) {
    EXPECT_NEAR(0.0, dgeqpf_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, dgeqpf_obj->diff_jpvt, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgeqpf_test, dgeqpf2) {
    EXPECT_NEAR(0.0, dgeqpf_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, dgeqpf_obj->diff_jpvt, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgeqpf_test, dgeqpf3) {
    EXPECT_NEAR(0.0, dgeqpf_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, dgeqpf_obj->diff_jpvt, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgeqpf_test, dgeqpf4) {
    EXPECT_NEAR(0.0, dgeqpf_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, dgeqpf_obj->diff_jpvt, LAPACKE_GTEST_THRESHOLD);
}
 