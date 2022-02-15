#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"


#define lag2c_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (sa!=NULL)    free(sa); \
if (saref!=NULL)    free(saref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin dcomplex_common_parameters  class definition */
class lag2c_dcomplex_parameters{

   public:
	int bufsize, bufsize_sa;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n, m;
	lapack_complex_double* A;
	lapack_int lda, ldsa;
	/*Output Parameter*/
	lapack_complex_double *Aref;
	lapack_complex_float *sa , *saref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      lag2c_dcomplex_parameters (int matrix_layout, lapack_int m , lapack_int n);
      ~lag2c_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
lag2c_dcomplex_parameters::lag2c_dcomplex_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i)
{
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n lag2c dcomplex:  m: %d, n: %d\n",  m, n);
	#endif
	
	//if (matrix_layout == LAPACK_COL_MAJOR) {
		lda = m;
		ldsa = n;
		bufsize = lda*n;
		bufsize_sa = ldsa*m;
     
    /*}else if (matrix_layout == LAPACK_ROW_MAJOR) {
		lda = n;
		ldsa = m;
		bufsize = lda*m;
		bufsize_sa = ldsa*n;     
    }else
                EXPECT_TRUE(false) << "matrix_layout invalid";
	*/
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&sa, &saref, bufsize_sa);
	if ((A==NULL) || (Aref==NULL) ||\
		(sa == NULL) || (saref == NULL)){
		EXPECT_FALSE( true) << "lag2c_double_parameters object: malloc error.";
		lag2c_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( A, Aref, bufsize);

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
lag2c_dcomplex_parameters :: ~lag2c_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" lag2c_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lag2c_free();

}
/*  Test fixture class definition */
class zlag2c_test  : public  ::testing::Test {
public:
   lag2c_dcomplex_parameters  *zlag2c_obj;
   void SetUp();
   void TearDown () { delete zlag2c_obj; }
};

void zlag2c_test::SetUp(){

    /* LAPACKE zlag2c prototype */
    typedef int (*Fptr_NL_LAPACKE_zlag2c) (int matrix_layout, lapack_int m, lapack_int n,\
                           const lapack_complex_double* a, lapack_int lda,\
                           lapack_complex_float* sa, lapack_int ldsa);

    Fptr_NL_LAPACKE_zlag2c zlag2c;

    zlag2c_obj = new lag2c_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
						eig_paramslist[idx].m,
						eig_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    zlag2c_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zlag2c_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zlag2c_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zlag2c_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zlag2c = (Fptr_NL_LAPACKE_zlag2c)dlsym(zlag2c_obj->hModule, "LAPACKE_zlag2c");
    ASSERT_TRUE(zlag2c != NULL) << "failed to get the Netlib LAPACKE_zlag2c symbol";
    

    zlag2c_obj->inforef = zlag2c( zlag2c_obj->matrix_layout, zlag2c_obj->m, zlag2c_obj->n,
								(const lapack_complex_double *)zlag2c_obj->Aref,
								zlag2c_obj->lda, zlag2c_obj->saref, zlag2c_obj->ldsa);

    /* Compute libflame's Lapacke o/p  */
    zlag2c_obj->info = LAPACKE_zlag2c( zlag2c_obj->matrix_layout, zlag2c_obj->m, zlag2c_obj->n,
										(const lapack_complex_double*)zlag2c_obj->A, 
										zlag2c_obj->lda, zlag2c_obj->sa, zlag2c_obj->ldsa);

    if( zlag2c_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zlag2c is wrong\n", zlag2c_obj->info );
    }
    if( zlag2c_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zlag2c is wrong\n", 
        zlag2c_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zlag2c_obj->diff =  computeDiff_c( zlag2c_obj->bufsize_sa,
                zlag2c_obj->sa, zlag2c_obj->saref );

}

TEST_F(zlag2c_test, zlag2c1) {
    EXPECT_NEAR(0.0, zlag2c_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zlag2c_test, zlag2c2) {
    EXPECT_NEAR(0.0, zlag2c_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zlag2c_test, zlag2c3) {
    EXPECT_NEAR(0.0, zlag2c_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zlag2c_test, zlag2c4) {
    EXPECT_NEAR(0.0, zlag2c_obj->diff, LAPACKE_EIG_THRESHOLD);
}
