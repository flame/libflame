#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"


#define hbgst_free() \
if (ab!=NULL)    free(ab); \
if (abref!=NULL) free(abref);\
if (bb!=NULL)  free(bb);\
if (bbref!=NULL) free(bbref); \
if (x != NULL) free(x); \
if (xref!=NULL) free(xref);\
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule);\
if( bModule != NULL) dlclose(bModule); \
if(lModule != NULL) dlclose(lModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin scomplex_common_parameters  class definition */
class hbgst_scomplex_parameters{

   public:
	int bufsize_a;
	int bufsize_b;
	void *hModule, *dModule;
	void *bModule, *lModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int ka, kb;
	lapack_complex_float* ab;
	char uplo, vect;
	lapack_int ldab,ldbb, ldx ;
	/*Output Parameter*/
	lapack_complex_float* bb, *x;
	lapack_complex_float *abref, *bbref, *xref;
	/*Return Values*/
	lapack_int info, inforef;
	lapack_int info_pbstf, inforef_pbstf;

   public:
      hbgst_scomplex_parameters (int matrix_layout, char vect, char uplo, lapack_int n, lapack_int ka, lapack_int kb, lapack_int ldx);
      ~hbgst_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
hbgst_scomplex_parameters:: hbgst_scomplex_parameters (int matrix_layout_i, char vect_i, char uplo_i, lapack_int n_i, lapack_int ka_i, lapack_int kb_i,lapack_int ldx_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	ldx = ldx_i;
	vect = vect_i;
	ka = ka_i;
	kb = kb_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n hbgst scomplex: matrix_layout = %d, vect:%c, uplo:%c, n: %d, ka:%d, kb:%d,ldx:%d \n", matrix_layout, vect, uplo, n, ka, kb,ldx);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
    {
		ldab = ka+1;
		ldbb = kb+1;
		bufsize_a = ldab*n;
		bufsize_b = ldbb*n;
			
    }else if (matrix_layout == LAPACK_ROW_MAJOR)
    {
		ldab = n;
		ldbb = n;
		bufsize_a = (ldab*(ka + 1));
		bufsize_b = ldbb*(kb + 1);
    }else
    EXPECT_TRUE(false) << "matrix_layout invalid";
	/*correction for the vector value*/
	if (vect == 'U')
		vect = 'N';
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&ab, &abref, bufsize_a);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&bb, &bbref, bufsize_b);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&x, &xref, (ldx*n));
	if ((ab==NULL) || (abref==NULL) || 
		(bb==NULL) || (bbref==NULL) || 
		(x == NULL) || (xref == NULL)){
		EXPECT_FALSE( true) << "hbgst_float_parameters object: malloc error.";
		hbgst_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( ab, abref, bufsize_a);
	lapacke_gtest_init_scomplex_buffer_pair_rand( bb, bbref, bufsize_b);
	//lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix(ab, abref, ldab, bufsize_a, uplo);
	//lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix(bb, bbref, ldbb, bufsize_b, uplo);
	/*initialize output matrix by 0 */
	
} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
hbgst_scomplex_parameters :: ~hbgst_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" hbgst_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   hbgst_free();

}
/*  Test fixture class definition */
class chbgst_test  : public  ::testing::Test {
public:
   hbgst_scomplex_parameters  *chbgst_obj;
   void SetUp();
   void TearDown () { delete chbgst_obj;}
};

void chbgst_test::SetUp(){

	/* LAPACKE cpbstf prototype */
    typedef int (*Fptr_NL_LAPACKE_cpbstf) (int matrix_layout, char uplo, lapack_int n, lapack_int kb, lapack_complex_float *bb, lapack_int ldbb);
	 Fptr_NL_LAPACKE_cpbstf cpbstf;
	 /* LAPACKE chbgst prototype */
    typedef int (*Fptr_NL_LAPACKE_chbgst) (int matrix_layout, char vect, char uplo, lapack_int n, lapack_int ka, lapack_int kb, 
											lapack_complex_float* ab, lapack_int ldab, const lapack_complex_float* bb, lapack_int ldbb,
											lapack_complex_float* x, lapack_int ldx);

    Fptr_NL_LAPACKE_chbgst chbgst;

    chbgst_obj = new hbgst_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
							eig_paramslist[idx].vect_rd,
							eig_paramslist[idx].uplo,
							eig_paramslist[idx].n,
							eig_paramslist[idx].k, 
							eig_paramslist[idx].kb,
							eig_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    chbgst_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    chbgst_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(chbgst_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(chbgst_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*chbgst library call */
    chbgst = (Fptr_NL_LAPACKE_chbgst)dlsym(chbgst_obj->hModule, "LAPACKE_chbgst");
    ASSERT_TRUE(chbgst != NULL) << "failed to get the Netlib LAPACKE_chbgst symbol";

	/*cpbstf library call*/
	chbgst_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    chbgst_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(chbgst_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(chbgst_obj->lModule != NULL) << "Netlib lapacke handle NULL";

    cpbstf = (Fptr_NL_LAPACKE_cpbstf)dlsym(chbgst_obj->lModule, "LAPACKE_cpbstf");
    ASSERT_TRUE(cpbstf != NULL) << "failed to get the Netlib LAPACKE_cpbstf symbol";    
    

    chbgst_obj->inforef_pbstf = cpbstf( chbgst_obj->matrix_layout, chbgst_obj->uplo,
								chbgst_obj->n, chbgst_obj->kb, chbgst_obj->bbref, chbgst_obj->ldbb);

    /* Compute libflame's Lapacke o/p  */
    chbgst_obj->info_pbstf = LAPACKE_cpbstf( chbgst_obj->matrix_layout, chbgst_obj->uplo,
										chbgst_obj->n,chbgst_obj->kb, chbgst_obj->bb, chbgst_obj->ldbb);

    if( chbgst_obj->info_pbstf < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cpbstf is wrong\n", chbgst_obj->info_pbstf );
    }
    if( chbgst_obj->inforef_pbstf < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cpbstf is wrong\n", 
        chbgst_obj->inforef_pbstf );
    }  
/*Compute chbgst's  o/p */
    chbgst_obj->inforef = chbgst( chbgst_obj->matrix_layout, chbgst_obj->vect, chbgst_obj->uplo, chbgst_obj->n, chbgst_obj->ka,
								chbgst_obj->kb, chbgst_obj->abref, chbgst_obj->ldab, (const lapack_complex_float*)chbgst_obj->bbref,
								chbgst_obj->ldbb, chbgst_obj->xref, chbgst_obj->ldx);

    /* Compute libflame's Lapacke o/p  */
    chbgst_obj->info = LAPACKE_chbgst( chbgst_obj->matrix_layout, chbgst_obj->vect, chbgst_obj->uplo, chbgst_obj->n, chbgst_obj->ka,
										chbgst_obj->kb, chbgst_obj->ab, chbgst_obj->ldab, (const lapack_complex_float*)chbgst_obj->bb,
										chbgst_obj->ldbb, chbgst_obj->x, chbgst_obj->ldx);

    if( chbgst_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_chbgst is wrong\n", chbgst_obj->info );
    }
    if( chbgst_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chbgst is wrong\n", 
        chbgst_obj->inforef );
    }

    /* Compute difference between libflame and Netlib o/ps  */
    chbgst_obj->diff =  computeDiff_c( chbgst_obj->bufsize_a, 
                chbgst_obj->ab, chbgst_obj->abref );

}

TEST_F(chbgst_test, chbgst1) {
    EXPECT_NEAR(0.0, chbgst_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(chbgst_test, chbgst2) {
    EXPECT_NEAR(0.0, chbgst_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(chbgst_test, chbgst3) {
    EXPECT_NEAR(0.0, chbgst_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(chbgst_test, chbgst4) {
    EXPECT_NEAR(0.0, chbgst_obj->diff, LAPACKE_EIG_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class hbgst_dcomplex_parameters{

   public:
	int bufsize_a;
	int bufsize_b;
	void *hModule, *dModule;
	void *bModule, *lModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int ka, kb;
	lapack_complex_double* ab;
	char uplo, vect;
	lapack_int ldab,ldbb, ldx ;
	/*Output Parameter*/
	lapack_complex_double* bb, *x;
	lapack_complex_double *abref, *bbref, *xref;
	/*Return Values*/
	lapack_int info, inforef;
	lapack_int info_pbstf, inforef_pbstf;

   public:
      hbgst_dcomplex_parameters (int matrix_layout, char vect, char uplo, lapack_int n, lapack_int ka, lapack_int kb, lapack_int ldx);
      ~hbgst_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
hbgst_dcomplex_parameters:: hbgst_dcomplex_parameters (int matrix_layout_i, char vect_i, char uplo_i, lapack_int n_i, lapack_int ka_i, lapack_int kb_i,lapack_int ldx_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	ldx = ldx_i;
	vect = vect_i;
	ka = ka_i;
	kb = kb_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n hbgst dcomplex: matrix_layout = %d, vect:%c, uplo:%c, n: %d, ka:%d, kb:%d,ldx:%d \n", matrix_layout, vect, uplo, n, ka, kb,ldx);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
    {
		ldab = ka+1;
		ldbb = kb+1;
		bufsize_a = ldab*n;
		bufsize_b = ldbb*n;
			
    }else if (matrix_layout == LAPACK_ROW_MAJOR)
    {
		ldab = n;
		ldbb = n;
		bufsize_a = (ldab*(ka + 1));
		bufsize_b = ldbb*(kb + 1);
    }else
    EXPECT_TRUE(false) << "matrix_layout invalid";
	/*correction for the vector value*/
	if (vect == 'U')
		vect = 'N';
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&ab, &abref, bufsize_a);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&bb, &bbref, bufsize_b);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&x, &xref, (ldx*n));
	if ((ab==NULL) || (abref==NULL) || 
		(bb==NULL) || (bbref==NULL) || 
		(x == NULL) || (xref == NULL)){
		EXPECT_FALSE( true) << "hbgst_double_parameters object: malloc error.";
		hbgst_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( ab, abref, bufsize_a);
	lapacke_gtest_init_dcomplex_buffer_pair_rand( bb, bbref, bufsize_b);
	//lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix(ab, abref, ldab, bufsize_a, uplo);
	//lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix(bb, bbref, ldbb, bufsize_b, uplo);
	/*initialize output matrix by 0 */
	
} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
hbgst_dcomplex_parameters :: ~hbgst_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" hbgst_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   hbgst_free();

}
/*  Test fixture class definition */
class zhbgst_test  : public  ::testing::Test {
public:
   hbgst_dcomplex_parameters  *zhbgst_obj;
   void SetUp();
   void TearDown () { delete zhbgst_obj;}
};

void zhbgst_test::SetUp(){

	/* LAPACKE zpbstf prototype */
    typedef int (*Fptr_NL_LAPACKE_zpbstf) (int matrix_layout, char uplo, lapack_int n, lapack_int kb, lapack_complex_double *bb, lapack_int ldbb);
	 Fptr_NL_LAPACKE_zpbstf zpbstf;
	 /* LAPACKE zhbgst prototype */
    typedef int (*Fptr_NL_LAPACKE_zhbgst) (int matrix_layout, char vect, char uplo, lapack_int n, lapack_int ka, lapack_int kb, 
											lapack_complex_double* ab, lapack_int ldab, const lapack_complex_double* bb, lapack_int ldbb,
											lapack_complex_double* x, lapack_int ldx);

    Fptr_NL_LAPACKE_zhbgst zhbgst;

    zhbgst_obj = new hbgst_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
							eig_paramslist[idx].vect_rd,
							eig_paramslist[idx].uplo,
							eig_paramslist[idx].n,
							eig_paramslist[idx].k, 
							eig_paramslist[idx].kb,
							eig_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    zhbgst_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zhbgst_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zhbgst_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zhbgst_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*zhbgst library call */
    zhbgst = (Fptr_NL_LAPACKE_zhbgst)dlsym(zhbgst_obj->hModule, "LAPACKE_zhbgst");
    ASSERT_TRUE(zhbgst != NULL) << "failed to get the Netlib LAPACKE_zhbgst symbol";

	/*zpbstf library call*/
	zhbgst_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zhbgst_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zhbgst_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zhbgst_obj->lModule != NULL) << "Netlib lapacke handle NULL";

    zpbstf = (Fptr_NL_LAPACKE_zpbstf)dlsym(zhbgst_obj->lModule, "LAPACKE_zpbstf");
    ASSERT_TRUE(zpbstf != NULL) << "failed to get the Netlib LAPACKE_zpbstf symbol";    
    

    zhbgst_obj->inforef_pbstf = zpbstf( zhbgst_obj->matrix_layout, zhbgst_obj->uplo,
								zhbgst_obj->n, zhbgst_obj->kb, zhbgst_obj->bbref, zhbgst_obj->ldbb);

    /* Compute libflame's Lapacke o/p  */
    zhbgst_obj->info_pbstf = LAPACKE_zpbstf( zhbgst_obj->matrix_layout, zhbgst_obj->uplo,
										zhbgst_obj->n,zhbgst_obj->kb, zhbgst_obj->bb, zhbgst_obj->ldbb);

    if( zhbgst_obj->info_pbstf < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zpbstf is wrong\n", zhbgst_obj->info_pbstf );
    }
    if( zhbgst_obj->inforef_pbstf < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zpbstf is wrong\n", 
        zhbgst_obj->inforef_pbstf );
    }  
/*Compute zhbgst's  o/p */
    zhbgst_obj->inforef = zhbgst( zhbgst_obj->matrix_layout, zhbgst_obj->vect, zhbgst_obj->uplo, zhbgst_obj->n, zhbgst_obj->ka,
								zhbgst_obj->kb, zhbgst_obj->abref, zhbgst_obj->ldab, (const lapack_complex_double*)zhbgst_obj->bbref,
								zhbgst_obj->ldbb, zhbgst_obj->xref, zhbgst_obj->ldx);

    /* Compute libflame's Lapacke o/p  */
    zhbgst_obj->info = LAPACKE_zhbgst( zhbgst_obj->matrix_layout, zhbgst_obj->vect, zhbgst_obj->uplo, zhbgst_obj->n, zhbgst_obj->ka,
										zhbgst_obj->kb, zhbgst_obj->ab, zhbgst_obj->ldab, (const lapack_complex_double*)zhbgst_obj->bb,
										zhbgst_obj->ldbb, zhbgst_obj->x, zhbgst_obj->ldx);

    if( zhbgst_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zhbgst is wrong\n", zhbgst_obj->info );
    }
    if( zhbgst_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhbgst is wrong\n", 
        zhbgst_obj->inforef );
    }

    /* Compute difference between libflame and Netlib o/ps  */
    zhbgst_obj->diff =  computeDiff_z( zhbgst_obj->bufsize_a, 
                zhbgst_obj->ab, zhbgst_obj->abref );

}

TEST_F(zhbgst_test, zhbgst1) {
    EXPECT_NEAR(0.0, zhbgst_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zhbgst_test, zhbgst2) {
    EXPECT_NEAR(0.0, zhbgst_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zhbgst_test, zhbgst3) {
    EXPECT_NEAR(0.0, zhbgst_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zhbgst_test, zhbgst4) {
    EXPECT_NEAR(0.0, zhbgst_obj->diff, LAPACKE_EIG_THRESHOLD);
}
