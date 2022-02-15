#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"


#define tgsyl_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (b!=NULL)    free(b); \
if (bref!=NULL) free(bref);\
if (c!=NULL)    free(c); \
if (cref!=NULL) free(cref);\
if (d!=NULL)    free(d); \
if (dref!=NULL) free(dref);\
if (e!=NULL)    free(e); \
if (eref!=NULL) free(eref);\
if (f!=NULL)    free(f); \
if (fref!=NULL) free(fref);\
if( hModule != NULL) dlclose(hModule); \
if( dModule != NULL) dlclose(dModule);

	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
/* Begin float_common_parameters  class definition */
class tgsyl_float_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	float diff_c, diff_f;
	float threshold;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	char trans;
	int ijob;
	float* A;
	float* b;
	float* c, *d, *e, *f;	
	lapack_int lda;
	lapack_int ldb;
	lapack_int ldc, ldd, lde, ldf;
	lapack_int isgn;	
	/*Output Parameter*/
	float scale, scaleref;
	float dif, difref;
	float *Aref;
	float *bref, *cref;
	float *dref, *eref, *fref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      tgsyl_float_parameters (int matrix_layout, char trans, int ijob, lapack_int m,\
					     lapack_int n );
      ~tgsyl_float_parameters ();

};

/* Constructor definition  float_common_parameters */
tgsyl_float_parameters:: tgsyl_float_parameters (int matrix_layout_i, char trans_i, int ijob_i, \
				lapack_int m_i, lapack_int n_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	trans = trans_i;
	ijob = ijob_i;
	n = n_i;
	m = m_i;
	lda = ldd = m;
	ldb = lde = n;
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{
		ldc = m;
		ldf = m;
		bufsize = ldc*n*sizeof(float);
	}
	else if (matrix_layout == LAPACK_ROW_MAJOR)
	{
		ldc = n;
		ldf = n;
		bufsize = ldc*m*sizeof(float);
	}
	else
		EXPECT_TRUE(false) << "matrix_layout invalid";

	#if LAPACKE_TEST_VERBOSE
	printf(" \n tgsyl float:  m: %d, n: %d lda: %d \n",  m, n, lda);
	#endif
	
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&A, &Aref, (lda*m)*sizeof(float));
	lapacke_gtest_alloc_float_buffer_pair(&d, &dref, (ldd*m)*sizeof(float));
	lapacke_gtest_alloc_float_buffer_pair(&b, &bref, (lda*n)*sizeof(float));
	lapacke_gtest_alloc_float_buffer_pair(&e, &eref, (lde*n)*sizeof(float));
	lapacke_gtest_alloc_float_buffer_pair(&c, &cref, bufsize);
	lapacke_gtest_alloc_float_buffer_pair(&f, &fref, bufsize);
	//lapacke_gtest_alloc_float_buffer_pair(&scale, &scaleref, n*sizeof(float));
	
	if ((A==NULL) || (Aref==NULL) || \
		(b==NULL) || (bref==NULL) || \
		(c == NULL) || (cref == NULL) ||
		(e == NULL) || (eref == NULL) ||
		(f == NULL) || (fref == NULL) ||
		(d == NULL) || (dref == NULL)){
		EXPECT_FALSE( true) << "tgsyl_float_parameters object: malloc error.";
		tgsyl_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand_custom_matrix( A, Aref, lda, m, 'U' );
	lapacke_gtest_init_float_buffer_pair_rand_custom_matrix( b, bref, ldb, n, 'U' );
	lapacke_gtest_init_float_buffer_pair_rand( c, cref, bufsize);
	lapacke_gtest_init_float_buffer_pair_rand_custom_matrix( d, dref, ldd, m, 'U' );
	lapacke_gtest_init_float_buffer_pair_rand_custom_matrix( e, eref, lde, n, 'U' );
	lapacke_gtest_init_float_buffer_pair_rand( f, fref, bufsize);

	scale = 0;
	scaleref = 0;

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
tgsyl_float_parameters :: ~tgsyl_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" tgsyl_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   tgsyl_free();

}
/*  Test fixture class definition */
class stgsyl_test  : public  ::testing::Test {
public:
   tgsyl_float_parameters  *stgsyl_obj;
   void SetUp();
   void TearDown () { delete stgsyl_obj; }
};

void stgsyl_test::SetUp(){

    /* LAPACKE stgsyl prototype */
    typedef int (*Fptr_NL_LAPACKE_stgsyl) ( int matrix_layout, char trans,
			lapack_int ijob, lapack_int m, lapack_int n, const float* a,
			lapack_int lda, const float* b, lapack_int ldb, float* c,
			lapack_int ldc, const float* d, lapack_int ldd, const float* e,
			lapack_int lde, float* f, lapack_int ldf, float* scale, float* dif );

    Fptr_NL_LAPACKE_stgsyl stgsyl;

    stgsyl_obj = new tgsyl_float_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].trans,
                           eig_non_sym_paramslist[idx].tgsen_ijob,
						   eig_paramslist[idx].m,
						   eig_paramslist[idx].n
						   );
    stgsyl_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;

    idx = Circular_Increment_Index(idx);
	

    stgsyl_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    stgsyl_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(stgsyl_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(stgsyl_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    stgsyl = (Fptr_NL_LAPACKE_stgsyl)dlsym(stgsyl_obj->hModule, "LAPACKE_stgsyl");
    ASSERT_TRUE(stgsyl != NULL) << "failed to get the Netlib LAPACKE_stgsyl symbol";
    
    stgsyl_obj->inforef = stgsyl( stgsyl_obj->matrix_layout,stgsyl_obj->trans, stgsyl_obj->ijob,
								stgsyl_obj->m, stgsyl_obj->n, stgsyl_obj->Aref, stgsyl_obj->lda,
								stgsyl_obj->bref, stgsyl_obj->ldb, stgsyl_obj->cref, stgsyl_obj->ldc,
								stgsyl_obj->d, stgsyl_obj->ldd, stgsyl_obj->e, stgsyl_obj->lde,
								stgsyl_obj->f, stgsyl_obj->ldf, &stgsyl_obj->scaleref,
								&stgsyl_obj->difref);

    /* Compute libflame's Lapacke o/p  */
    stgsyl_obj->info = LAPACKE_stgsyl( stgsyl_obj->matrix_layout,stgsyl_obj->trans, stgsyl_obj->ijob,\
								stgsyl_obj->m, stgsyl_obj->n, stgsyl_obj->A, stgsyl_obj->lda,
								stgsyl_obj->b, stgsyl_obj->ldb, stgsyl_obj->c, stgsyl_obj->ldc,
								stgsyl_obj->d, stgsyl_obj->ldd, stgsyl_obj->e, stgsyl_obj->lde,
								stgsyl_obj->f, stgsyl_obj->ldf, &stgsyl_obj->scale,
								&stgsyl_obj->dif);

    if( stgsyl_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_stgsyl is wrong\n", stgsyl_obj->info );
    }
    if( stgsyl_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_stgsyl is wrong\n", 
        stgsyl_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    stgsyl_obj->diff_c =  computeDiff_s( stgsyl_obj->bufsize, 
                stgsyl_obj->c, stgsyl_obj->cref );

    stgsyl_obj->diff_f =  computeDiff_s( stgsyl_obj->bufsize, 
                stgsyl_obj->f, stgsyl_obj->fref );

}

TEST_F(stgsyl_test, stgsyl1) {
    EXPECT_NEAR(0.0, stgsyl_obj->diff_c, stgsyl_obj->threshold);
    EXPECT_NEAR(0.0, stgsyl_obj->diff_f, stgsyl_obj->threshold);
    EXPECT_NEAR( stgsyl_obj->scaleref, stgsyl_obj->scale, stgsyl_obj->threshold);
    EXPECT_NEAR( stgsyl_obj->difref, stgsyl_obj->dif, stgsyl_obj->threshold);
}

TEST_F(stgsyl_test, stgsyl2) {
    EXPECT_NEAR(0.0, stgsyl_obj->diff_c, stgsyl_obj->threshold);
    EXPECT_NEAR(0.0, stgsyl_obj->diff_f, stgsyl_obj->threshold);
    EXPECT_NEAR( stgsyl_obj->scaleref, stgsyl_obj->scale, stgsyl_obj->threshold);
    EXPECT_NEAR( stgsyl_obj->difref, stgsyl_obj->dif, stgsyl_obj->threshold);
}

TEST_F(stgsyl_test, stgsyl3) {
    EXPECT_NEAR(0.0, stgsyl_obj->diff_c, stgsyl_obj->threshold);
    EXPECT_NEAR(0.0, stgsyl_obj->diff_f, stgsyl_obj->threshold);
    EXPECT_NEAR( stgsyl_obj->scaleref, stgsyl_obj->scale, stgsyl_obj->threshold);
    EXPECT_NEAR( stgsyl_obj->difref, stgsyl_obj->dif, stgsyl_obj->threshold);
}

TEST_F(stgsyl_test, stgsyl4) {
    EXPECT_NEAR(0.0, stgsyl_obj->diff_c, stgsyl_obj->threshold);
    EXPECT_NEAR(0.0, stgsyl_obj->diff_f, stgsyl_obj->threshold);
    EXPECT_NEAR( stgsyl_obj->scaleref, stgsyl_obj->scale, stgsyl_obj->threshold);
    EXPECT_NEAR( stgsyl_obj->difref, stgsyl_obj->dif, stgsyl_obj->threshold);
}

