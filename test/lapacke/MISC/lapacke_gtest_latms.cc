#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"


#define latms_free() \
if (a!=NULL)    free(a); \
if (aref!=NULL) free(aref);\
if (d!=NULL)    free(d); \
if (dref!=NULL)    free(dref); \
if (iseed!=NULL)  free(iseed);\
if (iseedref!=NULL) free(iseedref); \
if( hModule != NULL) dlclose(hModule); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
/* Begin scomplex_common_parameters  class definition */
class latms_scomplex_parameters{

   public:
	int bufsize, bufsize_d;
	int bufsize_iseed;
	int matrix_layout;
	void *hModule, *dModule;
	float diff, diff_iseed, diff_d;
   /*input parameters */
	lapack_int n, m;
	lapack_int lda;
	lapack_int mode, kl, ku;
	lapack_complex_float * a;
	float *d, cond, dmax;
	char dist, sym, pack;
	/*Output Parameter*/
	lapack_int* iseed, *iseedref;
	lapack_complex_float *aref;
	float *dref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      latms_scomplex_parameters (int matrix_layout, lapack_int m, lapack_int n, char dist, char sym, lapack_int mode,\
									float cond, float dmax, lapack_int kl, lapack_int ku, char pack);
      ~latms_scomplex_parameters ();

};

/* Constructor definition  scomplex_common_parameters */
latms_scomplex_parameters::latms_scomplex_parameters ( int matrix_layout_i, lapack_int m_i, lapack_int n_i, char dist_i,\
								char sym_i, lapack_int mode_i, float cond_i, float dmax_i, lapack_int kl_i, lapack_int ku_i, char pack_i)
{
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	dist = dist_i;
	sym = sym_i;
	mode = mode_i;
	cond = cond_i;
	dmax = dmax_i;
	kl = kl_i;
	ku = ku_i;
	pack = pack_i;
	#if LAPACKE_TEST_VERBOSE
		printf(" \n latms scomplex: matrix_layout: %d , m: %d, n:%d, dist:%c, sym:%c, mode:%d, cond:%d, dmax:%d, kl:%d, ku:%d, pack:%c \n",
			matrix_layout, m, n, dist, sym, mode, cond, dmax, kl, ku, pack);
	#endif
	
	if (dist == 'V')
		dist = 'S';
		
	if (sym  == 'B')
		sym = 'N';
	
	if (pack == 'p')
		pack = 'U';
	else if (pack == 'S')
		pack = 'L';
	
	/*if symmetric or Hermitian matrix kl = ku*/
	if ((sym == 'S') || (sym == 'P'))
		kl = ku;
	/*Based on pack LDA value differs*/
	if ((pack == 'N') ||(pack == 'U') ||(pack == 'L') ||(pack == 'C') ||(pack == 'R'))
	{
		if (matrix_layout == LAPACK_COL_MAJOR)
			lda = m;
		else if (matrix_layout == LAPACK_ROW_MAJOR)
			lda = n;
		
	}else if ((pack == 'B') || (pack == 'Q'))
	{
		lda = fla_min(kl, (m-1));
	}
	
	/*buffer sizes */
	bufsize = lda*n;
	bufsize_iseed = 4;
	bufsize_d = fla_min(m,n);
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&a, &aref, bufsize);
	lapacke_gtest_alloc_int_buffer_pair(&iseed, &iseedref, bufsize_iseed);
	lapacke_gtest_alloc_float_buffer_pair(&d, &dref, bufsize_d);
	if ((a==NULL) || (aref==NULL) ||
		(iseed==NULL) || (iseedref==NULL) ||
		(d==NULL) || (dref==NULL)){
		EXPECT_FALSE( true) << "latms_scomplex_parameters object: malloc error.";
		latms_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_int_buffer_pair_rand(iseed, iseedref, bufsize_iseed);
	lapacke_gtest_init_float_buffer_pair_rand(d, dref, n);

} /* end of Constructor  */

/* Destructor definition  'scomplex_common_parameters' class  */
latms_scomplex_parameters :: ~latms_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" latms_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   latms_free();

}
/*  Test fixture class definition */
class clatms_test  : public  ::testing::Test {
public:
   latms_scomplex_parameters  *clatms_obj;
   void SetUp();
   void TearDown () { delete clatms_obj; }
};

void clatms_test::SetUp(){

    /* LAPACKE clatms prototype */
    typedef int (*Fptr_NL_LAPACKE_clatms) (int matrix_layout, lapack_int m, lapack_int n,\
                           char dist, lapack_int* iseed, char sym, float* d,\
                           lapack_int mode, float cond, float dmax,\
                           lapack_int kl, lapack_int ku, char pack,\
                           lapack_complex_float* a, lapack_int lda );

    Fptr_NL_LAPACKE_clatms clatms;

    clatms_obj = new latms_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
												eig_paramslist[idx].m,
												eig_paramslist[idx].n,
												eig_paramslist[idx].vect_rd,
												eig_paramslist[idx].job,
												eig_paramslist[idx].nrhs,
												eig_paramslist[idx].itype,
												eig_paramslist[idx].n,
												lin_solver_paramslist[idx].kl,
												lin_solver_paramslist[idx].ku,
												eig_paramslist[idx].job);

    idx = Circular_Increment_Index(idx);
	

    clatms_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    clatms_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(clatms_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(clatms_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    clatms = (Fptr_NL_LAPACKE_clatms)dlsym(clatms_obj->hModule, "LAPACKE_clatms");
    ASSERT_TRUE(clatms != NULL) << "failed to get the Netlib LAPACKE_clatms symbol";
    

    clatms_obj->inforef = clatms(clatms_obj->matrix_layout, clatms_obj->m, clatms_obj->n, clatms_obj->dist, clatms_obj->iseedref,\
								clatms_obj->sym, clatms_obj->dref, clatms_obj->mode, clatms_obj->cond, clatms_obj->dmax, clatms_obj->kl,\
								clatms_obj->ku, clatms_obj->pack, clatms_obj->aref, clatms_obj->lda);

    /* Compute libflame's Lapacke o/p  */
    clatms_obj->info = LAPACKE_clatms(clatms_obj->matrix_layout, clatms_obj->m, clatms_obj->n, clatms_obj->dist, clatms_obj->iseed,\
								clatms_obj->sym, clatms_obj->d, clatms_obj->mode, clatms_obj->cond, clatms_obj->dmax, clatms_obj->kl,\
								clatms_obj->ku, clatms_obj->pack, clatms_obj->a, clatms_obj->lda);

    if( clatms_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_clatms is wrong\n", clatms_obj->info );
    }
    if( clatms_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_clatms is wrong\n", 
        clatms_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    clatms_obj->diff =  computeDiff_c( clatms_obj->bufsize,  clatms_obj->a, clatms_obj->aref );
	clatms_obj->diff_d =  computeDiff_s( clatms_obj->bufsize_d,  clatms_obj->d, clatms_obj->dref );

}

TEST_F(clatms_test, clatms1) {
    EXPECT_NEAR(0.0, clatms_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, clatms_obj->diff_d, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(clatms_test, clatms2) {
    EXPECT_NEAR(0.0, clatms_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, clatms_obj->diff_d, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(clatms_test, clatms3) {
    EXPECT_NEAR(0.0, clatms_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, clatms_obj->diff_d, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(clatms_test, clatms4) {
    EXPECT_NEAR(0.0, clatms_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, clatms_obj->diff_d, LAPACKE_GTEST_THRESHOLD);
}


/* Begin dcomplex_common_parameters  class definition */
