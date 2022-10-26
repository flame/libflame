#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"


#define geqp3_free() \
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
class geqp3_float_parameters{

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
      geqp3_float_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda);
      ~geqp3_float_parameters ();

};

/* Constructor definition  float_common_parameters */
geqp3_float_parameters:: geqp3_float_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n geqp3 float:  m: %d, n: %d lda: %d \n",  m, n, lda);
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
		EXPECT_FALSE( true) << "geqp3_float_parameters object: malloc error.";
		geqp3_free();
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
geqp3_float_parameters :: ~geqp3_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" geqp3_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   geqp3_free();

}
/*  Test fixture class definition */
class sgeqp3_test  : public  ::testing::Test {
public:
   geqp3_float_parameters  *sgeqp3_obj;
   void SetUp();
   void TearDown () { delete sgeqp3_obj; }
};

void sgeqp3_test::SetUp(){

    /* LAPACKE sgeqp3 prototype */
    typedef int (*Fptr_NL_LAPACKE_sgeqp3) (int matrix_layout, lapack_int m,lapack_int n,\
											float *A, lapack_int lda, lapack_int *jpvt, float* tau);

    Fptr_NL_LAPACKE_sgeqp3 sgeqp3;

    sgeqp3_obj = new geqp3_float_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);
	

    sgeqp3_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgeqp3_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgeqp3_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgeqp3_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    sgeqp3 = (Fptr_NL_LAPACKE_sgeqp3)dlsym(sgeqp3_obj->hModule, "LAPACKE_sgeqp3");
    ASSERT_TRUE(sgeqp3 != NULL) << "failed to get the Netlib LAPACKE_sgeqp3 symbol";
    

    sgeqp3_obj->inforef = sgeqp3( sgeqp3_obj->matrix_layout, sgeqp3_obj->m,
								sgeqp3_obj->n,sgeqp3_obj->Aref,
								sgeqp3_obj->lda, sgeqp3_obj->jpvtref, sgeqp3_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    sgeqp3_obj->info = LAPACKE_sgeqp3( sgeqp3_obj->matrix_layout, sgeqp3_obj->m,
										sgeqp3_obj->n,sgeqp3_obj->A, 
										sgeqp3_obj->lda, sgeqp3_obj->jpvt, sgeqp3_obj->tau);

    if( sgeqp3_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_sgeqp3 is wrong\n", sgeqp3_obj->info );
    }
    if( sgeqp3_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sgeqp3 is wrong\n", 
        sgeqp3_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    sgeqp3_obj->diff =  computeDiff_s( sgeqp3_obj->bufsize, 
                sgeqp3_obj->A, sgeqp3_obj->Aref );
	sgeqp3_obj->diff_jpvt = computeDiff_i( sgeqp3_obj->jpvt_size,\
							sgeqp3_obj->jpvt, sgeqp3_obj->jpvtref); 

}

TEST_F(sgeqp3_test, sgeqp31) {
    EXPECT_NEAR(0.0, sgeqp3_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, sgeqp3_obj->diff_jpvt, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgeqp3_test, sgeqp32) {
    EXPECT_NEAR(0.0, sgeqp3_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, sgeqp3_obj->diff_jpvt, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgeqp3_test, sgeqp33) {
    EXPECT_NEAR(0.0, sgeqp3_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, sgeqp3_obj->diff_jpvt, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgeqp3_test, sgeqp34) {
    EXPECT_NEAR(0.0, sgeqp3_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, sgeqp3_obj->diff_jpvt, LAPACKE_GTEST_THRESHOLD);
}


/* Begin double_common_parameters  class definition */
class geqp3_double_parameters{
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
      geqp3_double_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda);
      ~geqp3_double_parameters ();

};

/* Constructor definition  double_common_parameters */
geqp3_double_parameters:: geqp3_double_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n geqp3 double:  m: %d, n: %d lda: %d \n",  m, n, lda);
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
		EXPECT_FALSE( true) << "geqp3_double_parameters object: malloc error.";
		geqp3_free();
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
geqp3_double_parameters :: ~geqp3_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" geqp3_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   geqp3_free();

}
/*  Test fixture class definition */
class dgeqp3_test  : public  ::testing::Test {
public:
   geqp3_double_parameters  *dgeqp3_obj;
   void SetUp();
   void TearDown () { delete dgeqp3_obj; }
};

void dgeqp3_test::SetUp(){

    /* LAPACKE dgeqp3 prototype */
    typedef int (*Fptr_NL_LAPACKE_dgeqp3) (int matrix_layout, lapack_int m,lapack_int n,\
											double *A, lapack_int lda, lapack_int *jpvt, double* tau);

    Fptr_NL_LAPACKE_dgeqp3 dgeqp3;

    dgeqp3_obj = new geqp3_double_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);
	

    dgeqp3_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgeqp3_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgeqp3_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgeqp3_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    dgeqp3 = (Fptr_NL_LAPACKE_dgeqp3)dlsym(dgeqp3_obj->hModule, "LAPACKE_dgeqp3");
    ASSERT_TRUE(dgeqp3 != NULL) << "failed to get the Netlib LAPACKE_dgeqp3 symbol";
    

    dgeqp3_obj->inforef = dgeqp3( dgeqp3_obj->matrix_layout, dgeqp3_obj->m,
								dgeqp3_obj->n,dgeqp3_obj->Aref,
								dgeqp3_obj->lda, dgeqp3_obj->jpvtref, dgeqp3_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    dgeqp3_obj->info = LAPACKE_dgeqp3( dgeqp3_obj->matrix_layout, dgeqp3_obj->m,
										dgeqp3_obj->n,dgeqp3_obj->A, 
										dgeqp3_obj->lda, dgeqp3_obj->jpvt, dgeqp3_obj->tau);

    if( dgeqp3_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dgeqp3 is wrong\n", dgeqp3_obj->info );
    }
    if( dgeqp3_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dgeqp3 is wrong\n", 
        dgeqp3_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dgeqp3_obj->diff =  computeDiff_d( dgeqp3_obj->bufsize, 
                dgeqp3_obj->A, dgeqp3_obj->Aref );
	dgeqp3_obj->diff_jpvt = computeDiff_i( dgeqp3_obj->jpvt_size,\
							dgeqp3_obj->jpvt, dgeqp3_obj->jpvtref); 

}

TEST_F(dgeqp3_test, dgeqp31) {
    EXPECT_NEAR(0.0, dgeqp3_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, dgeqp3_obj->diff_jpvt, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgeqp3_test, dgeqp32) {
    EXPECT_NEAR(0.0, dgeqp3_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, dgeqp3_obj->diff_jpvt, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgeqp3_test, dgeqp33) {
    EXPECT_NEAR(0.0, dgeqp3_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, dgeqp3_obj->diff_jpvt, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgeqp3_test, dgeqp34) {
    EXPECT_NEAR(0.0, dgeqp3_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, dgeqp3_obj->diff_jpvt, LAPACKE_GTEST_THRESHOLD);
}
 