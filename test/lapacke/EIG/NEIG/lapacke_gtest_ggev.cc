#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"
#define LAPACKE_TEST_VERBOSE  (1)

#define ggev_free() \
    if (a!=NULL)        free(a); \
    if (aref!=NULL)     free(aref); \
    if (b!=NULL)        free(b); \
    if (bref!=NULL)     free(bref); \
    if (vsl!=NULL)       free(vsl); \
    if (vslref!=NULL)    free(vslref); \
    if (vsr!=NULL)       free(vsr); \
    if (vsrref!=NULL)    free(vsrref); \
    if (alphar!=NULL)     free(alphar); \
    if (alpharref!=NULL)     free(alpharref); \
    if (alphai!=NULL)     free(alphai); \
    if (alphairef!=NULL)     free(alphairef); \
    if (beta!=NULL)     free(beta); \
    if (betaref!=NULL)     free(betaref); \
    if (select!=NULL)     free(select); \
    if (selectref!=NULL)     free(selectref)


// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin ggev_float_common_parameters  class definition */
class ggev_float_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_a, diff_b;
    float diff_alphai, diff_alphar, diff_beta, diff_vsl, diff_vsr;

    float threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR	
    char jobvsl; // Must be 'N', or 'V'.
    char jobvsr; // Must be 'N', or 'V'.
	char sort; // Must be 'N', or 'S'.
	int *select, *selectref; // Specifies the eigenvectors to be computed.  

    lapack_int n; // The order of the matrices A and B (n≥ 0).
    lapack_int lda; //  The leading dimension of a
    lapack_int ldb; //  The leading dimension of b
    lapack_int ldvsl;
    lapack_int ldvsr;

    /* Input / Output parameters */
    float* a, *aref; // contains the n-by-n general matrix A.
    float* b, *bref; // contains the n-by-n upper triangular matrix B.

    /* Output parameters */
    lapack_int sdim, sdimref;
    float *vsl, *vslref; // left eigenvectors.
    float *vsr, *vsrref; // right eigenvectors

    float* alphar, *alpharref;
	float* alphai, *alphairef;
	float* beta, *betaref;
    /*Return Values */
    int info, inforef;

      ggev_float_parameters (int matrix_layout_i, char jobvsl_i,
						char jobvsr_i, char sort_i, lapack_int n_i );
      ~ggev_float_parameters ();
};

/* Constructor definition  float_common_parameters */
ggev_float_parameters:: ggev_float_parameters (int matrix_layout_i,
            char jobvsl_i, char jobvsr_i, char sort_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    jobvsl = jobvsl_i;
    jobvsr = jobvsr_i;
	sort = sort_i;
	
    n  = n_i;
    sdim = 0;
	sdimref = 0;

    lda = n;
    ldb = n;
    ldvsl = n;
    ldvsr = n;

    hModule = NULL;
    dModule = NULL;

    diff_a = 0;
    diff_b = 0;
    diff_vsl = 0;
    diff_vsr = 0;
    diff_alphai = 0;
    diff_alphar = 0;
    diff_beta = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n ggev float: matrix_layout: %d n: %d  jobvsl: %c \
jobvsr: %c \n", matrix_layout, n, jobvsr, jobvsl);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, lda*n );
    lapacke_gtest_alloc_float_buffer_pair( &b, &bref, ldb*n );
    lapacke_gtest_alloc_float_buffer_pair( &alphar, &alpharref, n );
    lapacke_gtest_alloc_float_buffer_pair( &alphai, &alphairef, n );
    lapacke_gtest_alloc_float_buffer_pair( &beta, &betaref, n );
    lapacke_gtest_alloc_int_buffer_pair( &select, &selectref, n );	
    lapacke_gtest_alloc_float_buffer_pair( &vsl, &vslref, n*ldvsl );
    lapacke_gtest_alloc_float_buffer_pair( &vsr, &vsrref, n*ldvsr );


    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (alphar==NULL) || (alpharref==NULL) || \
        (alphai==NULL) || (alphairef==NULL) || \
        (beta==NULL) || (betaref==NULL) || \
        (vsl==NULL) || (vslref==NULL) || \
        (vsr==NULL) || (vsrref==NULL) || \
		(select==NULL) || (selectref==NULL) ){
       EXPECT_FALSE( true) << "ggev_float_parameters object: malloc error. Exiting ";
       ggev_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, lda*n );
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, ldb*n );
    lapacke_gtest_init_float_buffer_pair_with_constant( vsl, vslref, ldvsl*n, 0 );
    lapacke_gtest_init_float_buffer_pair_with_constant( vsr, vsrref, ldvsr*n, 0 );
    lapacke_gtest_init_float_buffer_pair_with_constant( alphar, alpharref, n, 0);
    lapacke_gtest_init_float_buffer_pair_with_constant( alphai, alphairef, n, 0);
    lapacke_gtest_init_float_buffer_pair_with_constant( beta, betaref, n, 0);
    lapacke_gtest_init_int_buffer_pair_with_constant(select, selectref, n, -1);

   } /* end of Constructor  */
    

/* Destructor definition  'ggev_float_common_parameters' */
ggev_float_parameters :: ~ggev_float_parameters ()
{
   /* De-Allocate memory for the input matrices */
   ggev_free();
} 

//  Test fixture class definition
class sggev_test  : public  ::testing::Test {
public:
   ggev_float_parameters  *sggev_obj;
   void SetUp();  
   void TearDown () { delete sggev_obj; }
};

void sggev_test::SetUp()
{

    /* LAPACKE_sggev prototype */	
    typedef int (*Fptr_NL_LAPACKE_sggev) (int matrix_layout, char jobvl,
		char jobvr, lapack_int n, float* a, lapack_int lda, float* b,
		lapack_int ldb, float* alphar, float* alphai, float* beta,
		float* vl, lapack_int ldvl, float* vr, lapack_int ldvr);
				 
    Fptr_NL_LAPACKE_sggev SGGEV;

    sggev_obj = new  ggev_float_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].jobvsl,
										 eig_non_sym_paramslist[idx].jobvsr,
										 eig_non_sym_paramslist[idx].sort,
                                         eig_paramslist[idx].n );
                                         
    sggev_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    sggev_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sggev_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sggev_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sggev_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SGGEV = (Fptr_NL_LAPACKE_sggev)dlsym(sggev_obj->hModule, "LAPACKE_sggev");
    ASSERT_TRUE(SGGEV != NULL) << "failed to ppt the Netlib LAPACKE_sggev symbol";

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    sggev_obj->inforef = SGGEV(   sggev_obj->matrix_layout,
                                    sggev_obj->jobvsl,
                                    sggev_obj->jobvsr,
                                    sggev_obj->n,
                                    sggev_obj->aref,
                                    sggev_obj->lda,
                                    sggev_obj->bref,
                                    sggev_obj->ldb,
                                    sggev_obj->alpharref,
                                    sggev_obj->alphairef,
                                    sggev_obj->betaref,
									sggev_obj->vslref,
									sggev_obj->ldvsl,
									sggev_obj->vsrref,
									sggev_obj->ldvsr
                                    );

    /* Compute libflame's Lapacke o/p  */
    sggev_obj->inforef =  LAPACKE_sggev (   sggev_obj->matrix_layout,
                                    sggev_obj->jobvsl,
                                    sggev_obj->jobvsr,
                                    sggev_obj->n,
                                    sggev_obj->a,
                                    sggev_obj->lda,
                                    sggev_obj->b,
                                    sggev_obj->ldb,
                                    sggev_obj->alphar,
                                    sggev_obj->alphai,
                                    sggev_obj->beta,
									sggev_obj->vsl,
									sggev_obj->ldvsl,
									sggev_obj->vsr,
									sggev_obj->ldvsr
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    
    sggev_obj->diff_a =  computeDiff_s( (sggev_obj->lda)*(sggev_obj->n), 
                sggev_obj->a, sggev_obj->aref );

    sggev_obj->diff_b =  computeDiff_s( (sggev_obj->ldb)*(sggev_obj->n), 
                sggev_obj->b, sggev_obj->bref );

    sggev_obj->diff_vsl =  computeDiff_s( (sggev_obj->ldvsl)*(sggev_obj->n), 
                sggev_obj->vsl, sggev_obj->vslref );

    sggev_obj->diff_vsr =  computeDiff_s( (sggev_obj->ldvsr)*(sggev_obj->n), 
                sggev_obj->vsr, sggev_obj->vsrref );

    sggev_obj->diff_alphar =  computeDiff_s( sggev_obj->n, 
                sggev_obj->alphar, sggev_obj->alpharref );

    sggev_obj->diff_alphai =  computeDiff_s( sggev_obj->n, 
                sggev_obj->alphai, sggev_obj->alphairef );

    sggev_obj->diff_beta =  computeDiff_s( sggev_obj->n, 
                sggev_obj->beta, sggev_obj->betaref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n ggev float: \n diff_a: %f \n diff_b: %f \n \
diff_vsl: %f \n diff_vsr: %f \n diff_alphar: %f \n\
diff_alphai: %f \n diff_beta: %f \n ",
       sggev_obj->diff_a, sggev_obj->diff_b, sggev_obj->diff_vsl,
	   sggev_obj->diff_vsr, sggev_obj->diff_alphar,
       sggev_obj->diff_alphai, sggev_obj->diff_beta );
#endif
}

TEST_F(sggev_test, sggev1) {
    //EXPECT_NEAR(0.0, sggev_obj->diff_a, sggev_obj->threshold);
    //EXPECT_NEAR(0.0, sggev_obj->diff_b, sggev_obj->threshold);
    EXPECT_NEAR(0.0, sggev_obj->diff_vsl, sggev_obj->threshold);
    EXPECT_NEAR(0.0, sggev_obj->diff_vsr, sggev_obj->threshold);
    EXPECT_NEAR(0.0, sggev_obj->diff_alphar, sggev_obj->threshold);
    EXPECT_NEAR(0.0, sggev_obj->diff_alphai, sggev_obj->threshold);
    EXPECT_NEAR(0.0, sggev_obj->diff_beta, sggev_obj->threshold);
}

TEST_F(sggev_test, sggev2) {
    //EXPECT_NEAR(0.0, sggev_obj->diff_a, sggev_obj->threshold);
    //EXPECT_NEAR(0.0, sggev_obj->diff_b, sggev_obj->threshold);
    EXPECT_NEAR(0.0, sggev_obj->diff_vsl, sggev_obj->threshold);
    EXPECT_NEAR(0.0, sggev_obj->diff_vsr, sggev_obj->threshold);
    EXPECT_NEAR(0.0, sggev_obj->diff_alphar, sggev_obj->threshold);
    EXPECT_NEAR(0.0, sggev_obj->diff_alphai, sggev_obj->threshold);
    EXPECT_NEAR(0.0, sggev_obj->diff_beta, sggev_obj->threshold);
}
TEST_F(sggev_test, sggev3) {
    //EXPECT_NEAR(0.0, sggev_obj->diff_a, sggev_obj->threshold);
    //EXPECT_NEAR(0.0, sggev_obj->diff_b, sggev_obj->threshold);
    EXPECT_NEAR(0.0, sggev_obj->diff_vsl, sggev_obj->threshold);
    EXPECT_NEAR(0.0, sggev_obj->diff_vsr, sggev_obj->threshold);
    EXPECT_NEAR(0.0, sggev_obj->diff_alphar, sggev_obj->threshold);
    EXPECT_NEAR(0.0, sggev_obj->diff_alphai, sggev_obj->threshold);
    EXPECT_NEAR(0.0, sggev_obj->diff_beta, sggev_obj->threshold);
}
TEST_F(sggev_test, sggev4) {
    //EXPECT_NEAR(0.0, sggev_obj->diff_a, sggev_obj->threshold);
    //EXPECT_NEAR(0.0, sggev_obj->diff_b, sggev_obj->threshold);
    EXPECT_NEAR(0.0, sggev_obj->diff_vsl, sggev_obj->threshold);
    EXPECT_NEAR(0.0, sggev_obj->diff_vsr, sggev_obj->threshold);
    EXPECT_NEAR(0.0, sggev_obj->diff_alphar, sggev_obj->threshold);
    EXPECT_NEAR(0.0, sggev_obj->diff_alphai, sggev_obj->threshold);
    EXPECT_NEAR(0.0, sggev_obj->diff_beta, sggev_obj->threshold);
}


/* Begin ggev_double_common_parameters  class definition */
class ggev_double_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    double diff_a, diff_b;
    double diff_alphai, diff_alphar, diff_beta, diff_vsl, diff_vsr;

    float threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR	
    char jobvsl; // Must be 'N', or 'V'.
    char jobvsr; // Must be 'N', or 'V'.
	char sort; // Must be 'N', or 'S'.
	int *select, *selectref; // Specifies the eigenvectors to be computed.  

    lapack_int n; // The order of the matrices A and B (n≥ 0).
    lapack_int lda; //  The leading dimension of a
    lapack_int ldb; //  The leading dimension of b
    lapack_int ldvsl;
    lapack_int ldvsr;

    /* Input / Output parameters */
    double* a, *aref; // contains the n-by-n general matrix A.
    double* b, *bref; // contains the n-by-n upper triangular matrix B.

    /* Output parameters */
    lapack_int sdim, sdimref;
    double *vsl, *vslref; // left eigenvectors.
    double *vsr, *vsrref; // right eigenvectors

    double* alphar, *alpharref;
	double* alphai, *alphairef;
	double* beta, *betaref;
    /*Return Values */
    int info, inforef;

      ggev_double_parameters (int matrix_layout_i, char jobvsl_i,
						char jobvsr_i, char sort_i, lapack_int n_i );
      ~ggev_double_parameters ();
};

/* Constructor definition  double_common_parameters */
ggev_double_parameters:: ggev_double_parameters (int matrix_layout_i,
            char jobvsl_i, char jobvsr_i, char sort_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    jobvsl = jobvsl_i;
    jobvsr = jobvsr_i;
	sort = sort_i;
	
    n  = n_i;
    sdim = 0;
	sdimref = 0;

    lda = n;
    ldb = n;
    ldvsl = n;
    ldvsr = n;

    hModule = NULL;
    dModule = NULL;

    diff_a = 0;
    diff_b = 0;
    diff_vsl = 0;
    diff_vsr = 0;
    diff_alphai = 0;
    diff_alphar = 0;
    diff_beta = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n ggev double: matrix_layout: %d n: %d  jobvsl: %c \
jobvsr: %c \n", matrix_layout, n, jobvsr, jobvsl);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, lda*n );
    lapacke_gtest_alloc_double_buffer_pair( &b, &bref, ldb*n );
    lapacke_gtest_alloc_double_buffer_pair( &alphar, &alpharref, n );
    lapacke_gtest_alloc_double_buffer_pair( &alphai, &alphairef, n );
    lapacke_gtest_alloc_double_buffer_pair( &beta, &betaref, n );
    lapacke_gtest_alloc_int_buffer_pair( &select, &selectref, n );	
    lapacke_gtest_alloc_double_buffer_pair( &vsl, &vslref, n*ldvsl );
    lapacke_gtest_alloc_double_buffer_pair( &vsr, &vsrref, n*ldvsr );


    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (alphar==NULL) || (alpharref==NULL) || \
        (alphai==NULL) || (alphairef==NULL) || \
        (beta==NULL) || (betaref==NULL) || \
        (vsl==NULL) || (vslref==NULL) || \
        (vsr==NULL) || (vsrref==NULL) || \
		(select==NULL) || (selectref==NULL) ){
       EXPECT_FALSE( true) << "ggev_double_parameters object: malloc error. Exiting ";
       ggev_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, lda*n );
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, ldb*n );
    lapacke_gtest_init_double_buffer_pair_with_constant( vsl, vslref, ldvsl*n, 0 );
    lapacke_gtest_init_double_buffer_pair_with_constant( vsr, vsrref, ldvsr*n, 0 );
    lapacke_gtest_init_double_buffer_pair_with_constant( alphar, alpharref, n, 0);
    lapacke_gtest_init_double_buffer_pair_with_constant( alphai, alphairef, n, 0);
    lapacke_gtest_init_double_buffer_pair_with_constant( beta, betaref, n, 0);
    lapacke_gtest_init_int_buffer_pair_with_constant(select, selectref, n, -1);

   } /* end of Constructor  */
    

/* Destructor definition  'ggev_double_common_parameters' */
ggev_double_parameters :: ~ggev_double_parameters ()
{
   /* De-Allocate memory for the input matrices */
   ggev_free();
} 

//  Test fixture class definition
class dggev_test  : public  ::testing::Test {
public:
   ggev_double_parameters  *dggev_obj;
   void SetUp();  
   void TearDown () { delete dggev_obj; }
};

void dggev_test::SetUp()
{

    /* LAPACKE_dggev prototype */	
    typedef int (*Fptr_NL_LAPACKE_dggev) (int matrix_layout, char jobvl,
		char jobvr, lapack_int n, double* a, lapack_int lda, double* b,
		lapack_int ldb, double* alphar, double* alphai, double* beta,
		double* vl, lapack_int ldvl, double* vr, lapack_int ldvr);
				 
    Fptr_NL_LAPACKE_dggev DGGEV;

    dggev_obj = new  ggev_double_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].jobvsl,
										 eig_non_sym_paramslist[idx].jobvsr,
										 eig_non_sym_paramslist[idx].sort,
                                         eig_paramslist[idx].n );
                                         
    dggev_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    dggev_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dggev_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dggev_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dggev_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DGGEV = (Fptr_NL_LAPACKE_dggev)dlsym(dggev_obj->hModule, "LAPACKE_dggev");
    ASSERT_TRUE(DGGEV != NULL) << "failed to ppt the Netlib LAPACKE_dggev symbol";

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    dggev_obj->inforef = DGGEV(   dggev_obj->matrix_layout,
                                    dggev_obj->jobvsl,
                                    dggev_obj->jobvsr,
                                    dggev_obj->n,
                                    dggev_obj->aref,
                                    dggev_obj->lda,
                                    dggev_obj->bref,
                                    dggev_obj->ldb,
                                    dggev_obj->alpharref,
                                    dggev_obj->alphairef,
                                    dggev_obj->betaref,
									dggev_obj->vslref,
									dggev_obj->ldvsl,
									dggev_obj->vsrref,
									dggev_obj->ldvsr
                                    );

    /* Compute libflame's Lapacke o/p  */
    dggev_obj->inforef =  LAPACKE_dggev (   dggev_obj->matrix_layout,
                                    dggev_obj->jobvsl,
                                    dggev_obj->jobvsr,
                                    dggev_obj->n,
                                    dggev_obj->a,
                                    dggev_obj->lda,
                                    dggev_obj->b,
                                    dggev_obj->ldb,
                                    dggev_obj->alphar,
                                    dggev_obj->alphai,
                                    dggev_obj->beta,
									dggev_obj->vsl,
									dggev_obj->ldvsl,
									dggev_obj->vsr,
									dggev_obj->ldvsr
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    
    dggev_obj->diff_a =  computeDiff_d( (dggev_obj->lda)*(dggev_obj->n), 
                dggev_obj->a, dggev_obj->aref );

    dggev_obj->diff_b =  computeDiff_d( (dggev_obj->ldb)*(dggev_obj->n), 
                dggev_obj->b, dggev_obj->bref );

    dggev_obj->diff_vsl =  computeDiff_d( (dggev_obj->ldvsl)*(dggev_obj->n), 
                dggev_obj->vsl, dggev_obj->vslref );

    dggev_obj->diff_vsr =  computeDiff_d( (dggev_obj->ldvsr)*(dggev_obj->n), 
                dggev_obj->vsr, dggev_obj->vsrref );

    dggev_obj->diff_alphar =  computeDiff_d( dggev_obj->n, 
                dggev_obj->alphar, dggev_obj->alpharref );

    dggev_obj->diff_alphai =  computeDiff_d( dggev_obj->n, 
                dggev_obj->alphai, dggev_obj->alphairef );

    dggev_obj->diff_beta = computeDiff_d ( dggev_obj->n, 
                dggev_obj->beta, dggev_obj->betaref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n ggev double: \n diff_a: %f \n diff_b: %f \n \
diff_vsl: %f \n diff_vsr: %f \n diff_alphar: %f \n\
diff_alphai: %f \n diff_beta: %f \n ",
       dggev_obj->diff_a, dggev_obj->diff_b, dggev_obj->diff_vsl,
	   dggev_obj->diff_vsr, dggev_obj->diff_alphar,
       dggev_obj->diff_alphai, dggev_obj->diff_beta );
#endif
}

TEST_F(dggev_test, dggev1) {
    //EXPECT_NEAR(0.0, dggev_obj->diff_a, dggev_obj->threshold);
    //EXPECT_NEAR(0.0, dggev_obj->diff_b, dggev_obj->threshold);
    EXPECT_NEAR(0.0, dggev_obj->diff_vsl, dggev_obj->threshold);
    EXPECT_NEAR(0.0, dggev_obj->diff_vsr, dggev_obj->threshold);
    EXPECT_NEAR(0.0, dggev_obj->diff_alphar, dggev_obj->threshold);
    EXPECT_NEAR(0.0, dggev_obj->diff_alphai, dggev_obj->threshold);
    EXPECT_NEAR(0.0, dggev_obj->diff_beta, dggev_obj->threshold);
}

TEST_F(dggev_test, dggev2) {
    //EXPECT_NEAR(0.0, dggev_obj->diff_a, dggev_obj->threshold);
    //EXPECT_NEAR(0.0, dggev_obj->diff_b, dggev_obj->threshold);
    EXPECT_NEAR(0.0, dggev_obj->diff_vsl, dggev_obj->threshold);
    EXPECT_NEAR(0.0, dggev_obj->diff_vsr, dggev_obj->threshold);
    EXPECT_NEAR(0.0, dggev_obj->diff_alphar, dggev_obj->threshold);
    EXPECT_NEAR(0.0, dggev_obj->diff_alphai, dggev_obj->threshold);
    EXPECT_NEAR(0.0, dggev_obj->diff_beta, dggev_obj->threshold);
}
TEST_F(dggev_test, dggev3) {
    //EXPECT_NEAR(0.0, dggev_obj->diff_a, dggev_obj->threshold);
    //EXPECT_NEAR(0.0, dggev_obj->diff_b, dggev_obj->threshold);
    EXPECT_NEAR(0.0, dggev_obj->diff_vsl, dggev_obj->threshold);
    EXPECT_NEAR(0.0, dggev_obj->diff_vsr, dggev_obj->threshold);
    EXPECT_NEAR(0.0, dggev_obj->diff_alphar, dggev_obj->threshold);
    EXPECT_NEAR(0.0, dggev_obj->diff_alphai, dggev_obj->threshold);
    EXPECT_NEAR(0.0, dggev_obj->diff_beta, dggev_obj->threshold);
}
TEST_F(dggev_test, dggev4) {
    //EXPECT_NEAR(0.0, dggev_obj->diff_a, dggev_obj->threshold);
    //EXPECT_NEAR(0.0, dggev_obj->diff_b, dggev_obj->threshold);
    EXPECT_NEAR(0.0, dggev_obj->diff_vsl, dggev_obj->threshold);
    EXPECT_NEAR(0.0, dggev_obj->diff_vsr, dggev_obj->threshold);
    EXPECT_NEAR(0.0, dggev_obj->diff_alphar, dggev_obj->threshold);
    EXPECT_NEAR(0.0, dggev_obj->diff_alphai, dggev_obj->threshold);
    EXPECT_NEAR(0.0, dggev_obj->diff_beta, dggev_obj->threshold);
}

/* Begin ggev_scomplex_common_parameters  class definition */
class ggev_scomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_a, diff_b;
    float diff_alphai, diff_alphar, diff_beta, diff_vsl, diff_vsr;

    float threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR	
    char jobvsl; // Must be 'N', or 'V'.
    char jobvsr; // Must be 'N', or 'V'.
	char sort; // Must be 'N', or 'S'.
	int *select, *selectref; // Specifies the eigenvectors to be computed.  

    lapack_int n; // The order of the matrices A and B (n≥ 0).
    lapack_int lda; //  The leading dimension of a
    lapack_int ldb; //  The leading dimension of b
    lapack_int ldvsl;
    lapack_int ldvsr;

    /* Input / Output parameters */
    lapack_complex_float* a, *aref; // contains the n-by-n general matrix A.
    lapack_complex_float* b, *bref; // contains the n-by-n upper triangular matrix B.

    /* Output parameters */
    lapack_int sdim, sdimref;
    lapack_complex_float *vsl, *vslref; // left eigenvectors.
    lapack_complex_float *vsr, *vsrref; // right eigenvectors

    lapack_complex_float* alphar, *alpharref;
	lapack_complex_float* alphai, *alphairef;
	lapack_complex_float* beta, *betaref;
    /*Return Values */
    int info, inforef;

      ggev_scomplex_parameters (int matrix_layout_i, char jobvsl_i,
						char jobvsr_i, char sort_i, lapack_int n_i );
      ~ggev_scomplex_parameters ();
};

/* Constructor definition  lapack_complex_float_common_parameters */
ggev_scomplex_parameters:: ggev_scomplex_parameters (int matrix_layout_i,
            char jobvsl_i, char jobvsr_i, char sort_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    jobvsl = jobvsl_i;
    jobvsr = jobvsr_i;
	sort = sort_i;
	
    n  = n_i;
    sdim = 0;
	sdimref = 0;

    lda = n;
    ldb = n;
    ldvsl = n;
    ldvsr = n;

    hModule = NULL;
    dModule = NULL;

    diff_a = 0;
    diff_b = 0;
    diff_vsl = 0;
    diff_vsr = 0;
    diff_alphai = 0;
    diff_alphar = 0;
    diff_beta = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n ggev lapack_complex_float: matrix_layout: %d n: %d  jobvsl: %c \
jobvsr: %c \n", matrix_layout, n, jobvsr, jobvsl);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, lda*n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &b, &bref, ldb*n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &alphar, &alpharref, n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &alphai, &alphairef, n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &beta, &betaref, n );
    lapacke_gtest_alloc_int_buffer_pair( &select, &selectref, n );	
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &vsl, &vslref, n*ldvsl );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &vsr, &vsrref, n*ldvsr );


    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (alphar==NULL) || (alpharref==NULL) || \
        (alphai==NULL) || (alphairef==NULL) || \
        (beta==NULL) || (betaref==NULL) || \
        (vsl==NULL) || (vslref==NULL) || \
        (vsr==NULL) || (vsrref==NULL) || \
		(select==NULL) || (selectref==NULL) ){
       EXPECT_FALSE( true) << "ggev_scomplex_parameters object: malloc error. Exiting ";
       ggev_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, lda*n );
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, ldb*n );
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( vsl, vslref, ldvsl*n, 0 );
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( vsr, vsrref, ldvsr*n, 0 );
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( alphar, alpharref, n, 0);
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( alphai, alphairef, n, 0);
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( beta, betaref, n, 0);
    lapacke_gtest_init_int_buffer_pair_with_constant(select, selectref, n, -1);

   } /* end of Constructor  */
    

/* Destructor definition  'ggev_scomplex_common_parameters' */
ggev_scomplex_parameters :: ~ggev_scomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   ggev_free();
} 

//  Test fixture class definition
class cggev_test  : public  ::testing::Test {
public:
   ggev_scomplex_parameters  *cggev_obj;
   void SetUp();  
   void TearDown () { delete cggev_obj; }
};

void cggev_test::SetUp()
{

    /* LAPACKE_cggev prototype */	
    typedef int (*Fptr_NL_LAPACKE_cggev) (int matrix_layout, char jobvl,
		char jobvr, lapack_int n, lapack_complex_float* a, lapack_int lda,
		lapack_complex_float* b, lapack_int ldb, lapack_complex_float* alpha,
		lapack_complex_float* beta, lapack_complex_float* vl, lapack_int ldvl,
		lapack_complex_float* vr, lapack_int ldvr);
				 
    Fptr_NL_LAPACKE_cggev CGGEV;

    cggev_obj = new  ggev_scomplex_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].jobvsl,
										 eig_non_sym_paramslist[idx].jobvsr,
										 eig_non_sym_paramslist[idx].sort,
                                         eig_paramslist[idx].n );
                                         
    cggev_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    cggev_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cggev_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cggev_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cggev_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CGGEV = (Fptr_NL_LAPACKE_cggev)dlsym(cggev_obj->hModule, "LAPACKE_cggev");
    ASSERT_TRUE(CGGEV != NULL) << "failed to ppt the Netlib LAPACKE_cggev symbol";

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    cggev_obj->inforef = CGGEV(   cggev_obj->matrix_layout,
                                    cggev_obj->jobvsl,
                                    cggev_obj->jobvsr,
                                    cggev_obj->n,
                                    cggev_obj->aref,
                                    cggev_obj->lda,
                                    cggev_obj->bref,
                                    cggev_obj->ldb,
                                    cggev_obj->alpharref,
                                    cggev_obj->betaref,
									cggev_obj->vslref,
									cggev_obj->ldvsl,
									cggev_obj->vsrref,
									cggev_obj->ldvsr
                                    );

    /* Compute libflame's Lapacke o/p  */
    cggev_obj->inforef =  LAPACKE_cggev (   cggev_obj->matrix_layout,
                                    cggev_obj->jobvsl,
                                    cggev_obj->jobvsr,
                                    cggev_obj->n,
                                    cggev_obj->a,
                                    cggev_obj->lda,
                                    cggev_obj->b,
                                    cggev_obj->ldb,
                                    cggev_obj->alphar,
                                    cggev_obj->beta,
									cggev_obj->vsl,
									cggev_obj->ldvsl,
									cggev_obj->vsr,
									cggev_obj->ldvsr
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    
    cggev_obj->diff_a =  computeDiff_c( (cggev_obj->lda)*(cggev_obj->n), 
                cggev_obj->a, cggev_obj->aref );

    cggev_obj->diff_b =  computeDiff_c( (cggev_obj->ldb)*(cggev_obj->n), 
                cggev_obj->b, cggev_obj->bref );

    cggev_obj->diff_vsl =  computeDiff_c( (cggev_obj->ldvsl)*(cggev_obj->n), 
                cggev_obj->vsl, cggev_obj->vslref );

    cggev_obj->diff_vsr =  computeDiff_c( (cggev_obj->ldvsr)*(cggev_obj->n), 
                cggev_obj->vsr, cggev_obj->vsrref );

    cggev_obj->diff_alphar =  computeDiff_c( cggev_obj->n, 
                cggev_obj->alphar, cggev_obj->alpharref );

    cggev_obj->diff_beta =  computeDiff_c( cggev_obj->n, 
                cggev_obj->beta, cggev_obj->betaref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n ggev lapack_complex_float: \n diff_a: %f \n diff_b: %f \n \
diff_vsl: %f \n diff_vsr: %f \n diff_alphar: %f \n\
 \n diff_beta: %f \n ",
       cggev_obj->diff_a, cggev_obj->diff_b, cggev_obj->diff_vsl,
	   cggev_obj->diff_vsr, cggev_obj->diff_alphar,
        cggev_obj->diff_beta );
#endif
}

TEST_F(cggev_test, cggev1) {
    //EXPECT_NEAR(0.0, cggev_obj->diff_a, cggev_obj->threshold);
    //EXPECT_NEAR(0.0, cggev_obj->diff_b, cggev_obj->threshold);
    EXPECT_NEAR(0.0, cggev_obj->diff_vsl, cggev_obj->threshold);
    EXPECT_NEAR(0.0, cggev_obj->diff_vsr, cggev_obj->threshold);
    EXPECT_NEAR(0.0, cggev_obj->diff_alphar, cggev_obj->threshold);
    EXPECT_NEAR(0.0, cggev_obj->diff_beta, cggev_obj->threshold);
}

TEST_F(cggev_test, cggev2) {
    //EXPECT_NEAR(0.0, cggev_obj->diff_a, cggev_obj->threshold);
    //EXPECT_NEAR(0.0, cggev_obj->diff_b, cggev_obj->threshold);
    EXPECT_NEAR(0.0, cggev_obj->diff_vsl, cggev_obj->threshold);
    EXPECT_NEAR(0.0, cggev_obj->diff_vsr, cggev_obj->threshold);
    EXPECT_NEAR(0.0, cggev_obj->diff_alphar, cggev_obj->threshold);
    EXPECT_NEAR(0.0, cggev_obj->diff_beta, cggev_obj->threshold);
}
TEST_F(cggev_test, cggev3) {
    //EXPECT_NEAR(0.0, cggev_obj->diff_a, cggev_obj->threshold);
    //EXPECT_NEAR(0.0, cggev_obj->diff_b, cggev_obj->threshold);
    EXPECT_NEAR(0.0, cggev_obj->diff_vsl, cggev_obj->threshold);
    EXPECT_NEAR(0.0, cggev_obj->diff_vsr, cggev_obj->threshold);
    EXPECT_NEAR(0.0, cggev_obj->diff_alphar, cggev_obj->threshold);
    EXPECT_NEAR(0.0, cggev_obj->diff_beta, cggev_obj->threshold);
}
TEST_F(cggev_test, cggev4) {
    //EXPECT_NEAR(0.0, cggev_obj->diff_a, cggev_obj->threshold);
    //EXPECT_NEAR(0.0, cggev_obj->diff_b, cggev_obj->threshold);
    EXPECT_NEAR(0.0, cggev_obj->diff_vsl, cggev_obj->threshold);
    EXPECT_NEAR(0.0, cggev_obj->diff_vsr, cggev_obj->threshold);
    EXPECT_NEAR(0.0, cggev_obj->diff_alphar, cggev_obj->threshold);
    EXPECT_NEAR(0.0, cggev_obj->diff_beta, cggev_obj->threshold);
}

/* Begin ggev_dcomplex_common_parameters  class definition */
class ggev_dcomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    double diff_a, diff_b;
    double diff_alphai, diff_alphar, diff_beta, diff_vsl, diff_vsr;

    float threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR	
    char jobvsl; // Must be 'N', or 'V'.
    char jobvsr; // Must be 'N', or 'V'.
	char sort; // Must be 'N', or 'S'.
	int *select, *selectref; // Specifies the eigenvectors to be computed.  

    lapack_int n; // The order of the matrices A and B (n≥ 0).
    lapack_int lda; //  The leading dimension of a
    lapack_int ldb; //  The leading dimension of b
    lapack_int ldvsl;
    lapack_int ldvsr;

    /* Input / Output parameters */
    lapack_complex_double* a, *aref; // contains the n-by-n general matrix A.
    lapack_complex_double* b, *bref; // contains the n-by-n upper triangular matrix B.

    /* Output parameters */
    lapack_int sdim, sdimref;
    lapack_complex_double *vsl, *vslref; // left eigenvectors.
    lapack_complex_double *vsr, *vsrref; // right eigenvectors

    lapack_complex_double* alphar, *alpharref;
	lapack_complex_double* alphai, *alphairef;
	lapack_complex_double* beta, *betaref;
    /*Return Values */
    int info, inforef;

      ggev_dcomplex_parameters (int matrix_layout_i, char jobvsl_i,
						char jobvsr_i, char sort_i, lapack_int n_i );
      ~ggev_dcomplex_parameters ();
};

/* Constructor definition  lapack_complex_double_common_parameters */
ggev_dcomplex_parameters:: ggev_dcomplex_parameters (int matrix_layout_i,
            char jobvsl_i, char jobvsr_i, char sort_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    jobvsl = jobvsl_i;
    jobvsr = jobvsr_i;
	sort = sort_i;
	
    n  = n_i;
    sdim = 0;
	sdimref = 0;

    lda = n;
    ldb = n;
    ldvsl = n;
    ldvsr = n;

    hModule = NULL;
    dModule = NULL;

    diff_a = 0;
    diff_b = 0;
    diff_vsl = 0;
    diff_vsr = 0;
    diff_alphai = 0;
    diff_alphar = 0;
    diff_beta = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n ggev lapack_complex_double: matrix_layout: %d n: %d  jobvsl: %c \
jobvsr: %c \n", matrix_layout, n, jobvsr, jobvsl);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, lda*n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &b, &bref, ldb*n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &alphar, &alpharref, n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &alphai, &alphairef, n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &beta, &betaref, n );
    lapacke_gtest_alloc_int_buffer_pair( &select, &selectref, n );	
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &vsl, &vslref, n*ldvsl );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &vsr, &vsrref, n*ldvsr );


    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (alphar==NULL) || (alpharref==NULL) || \
        (alphai==NULL) || (alphairef==NULL) || \
        (beta==NULL) || (betaref==NULL) || \
        (vsl==NULL) || (vslref==NULL) || \
        (vsr==NULL) || (vsrref==NULL) || \
		(select==NULL) || (selectref==NULL) ){
       EXPECT_FALSE( true) << "ggev_dcomplex_parameters object: malloc error. Exiting ";
       ggev_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, lda*n );
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, ldb*n );
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( vsl, vslref, ldvsl*n, 0 );
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( vsr, vsrref, ldvsr*n, 0 );
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( alphar, alpharref, n, 0);
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( alphai, alphairef, n, 0);
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( beta, betaref, n, 0);
    lapacke_gtest_init_int_buffer_pair_with_constant(select, selectref, n, -1);

   } /* end of Constructor  */
    

/* Destructor definition  'ggev_dcomplex_common_parameters' */
ggev_dcomplex_parameters :: ~ggev_dcomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   ggev_free();
} 

//  Test fixture class definition
class zggev_test  : public  ::testing::Test {
public:
   ggev_dcomplex_parameters  *zggev_obj;
   void SetUp();  
   void TearDown () { delete zggev_obj; }
};

void zggev_test::SetUp()
{

    /* LAPACKE_zggev prototype */	
    typedef int (*Fptr_NL_LAPACKE_zggev) (int matrix_layout, char jobvl,
		char jobvr, lapack_int n, lapack_complex_double* a, lapack_int lda,
		lapack_complex_double* b, lapack_int ldb, lapack_complex_double* alpha,
		lapack_complex_double* beta, lapack_complex_double* vl, lapack_int ldvl,
		lapack_complex_double* vr, lapack_int ldvr);
				 
    Fptr_NL_LAPACKE_zggev ZGGEV;

    zggev_obj = new  ggev_dcomplex_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].jobvsl,
										 eig_non_sym_paramslist[idx].jobvsr,
										 eig_non_sym_paramslist[idx].sort,
                                         eig_paramslist[idx].n );
                                         
    zggev_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    zggev_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zggev_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zggev_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zggev_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZGGEV = (Fptr_NL_LAPACKE_zggev)dlsym(zggev_obj->hModule, "LAPACKE_zggev");
    ASSERT_TRUE(ZGGEV != NULL) << "failed to ppt the Netlib LAPACKE_zggev symbol";

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    zggev_obj->inforef = ZGGEV(   zggev_obj->matrix_layout,
                                    zggev_obj->jobvsl,
                                    zggev_obj->jobvsr,
                                    zggev_obj->n,
                                    zggev_obj->aref,
                                    zggev_obj->lda,
                                    zggev_obj->bref,
                                    zggev_obj->ldb,
                                    zggev_obj->alpharref,
                                    zggev_obj->betaref,
									zggev_obj->vslref,
									zggev_obj->ldvsl,
									zggev_obj->vsrref,
									zggev_obj->ldvsr
                                    );

    /* Compute libflame's Lapacke o/p  */
    zggev_obj->inforef =  LAPACKE_zggev (   zggev_obj->matrix_layout,
                                    zggev_obj->jobvsl,
                                    zggev_obj->jobvsr,
                                    zggev_obj->n,
                                    zggev_obj->a,
                                    zggev_obj->lda,
                                    zggev_obj->b,
                                    zggev_obj->ldb,
                                    zggev_obj->alphar,
                                    zggev_obj->beta,
									zggev_obj->vsl,
									zggev_obj->ldvsl,
									zggev_obj->vsr,
									zggev_obj->ldvsr
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    
    zggev_obj->diff_a =  computeDiff_z( (zggev_obj->lda)*(zggev_obj->n), 
                zggev_obj->a, zggev_obj->aref );

    zggev_obj->diff_b =  computeDiff_z( (zggev_obj->ldb)*(zggev_obj->n), 
                zggev_obj->b, zggev_obj->bref );

    zggev_obj->diff_vsl =  computeDiff_z( (zggev_obj->ldvsl)*(zggev_obj->n), 
                zggev_obj->vsl, zggev_obj->vslref );

    zggev_obj->diff_vsr =  computeDiff_z( (zggev_obj->ldvsr)*(zggev_obj->n), 
                zggev_obj->vsr, zggev_obj->vsrref );

    zggev_obj->diff_alphar =  computeDiff_z( zggev_obj->n, 
                zggev_obj->alphar, zggev_obj->alpharref );

    zggev_obj->diff_beta =  computeDiff_z( zggev_obj->n, 
                zggev_obj->beta, zggev_obj->betaref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n ggev lapack_complex_double: \n diff_a: %f \n diff_b: %f \n \
diff_vsl: %f \n diff_vsr: %f \n diff_alphar: %f \n\
 \n diff_beta: %f \n ",
       zggev_obj->diff_a, zggev_obj->diff_b, zggev_obj->diff_vsl,
	   zggev_obj->diff_vsr, zggev_obj->diff_alphar,
        zggev_obj->diff_beta );
#endif
}

TEST_F(zggev_test, zggev1) {
    //EXPECT_NEAR(0.0, zggev_obj->diff_a, zggev_obj->threshold);
    //EXPECT_NEAR(0.0, zggev_obj->diff_b, zggev_obj->threshold);
    EXPECT_NEAR(0.0, zggev_obj->diff_vsl, zggev_obj->threshold);
    EXPECT_NEAR(0.0, zggev_obj->diff_vsr, zggev_obj->threshold);
    EXPECT_NEAR(0.0, zggev_obj->diff_alphar, zggev_obj->threshold);
    EXPECT_NEAR(0.0, zggev_obj->diff_beta, zggev_obj->threshold);
}

TEST_F(zggev_test, zggev2) {
    //EXPECT_NEAR(0.0, zggev_obj->diff_a, zggev_obj->threshold);
    //EXPECT_NEAR(0.0, zggev_obj->diff_b, zggev_obj->threshold);
    EXPECT_NEAR(0.0, zggev_obj->diff_vsl, zggev_obj->threshold);
    EXPECT_NEAR(0.0, zggev_obj->diff_vsr, zggev_obj->threshold);
    EXPECT_NEAR(0.0, zggev_obj->diff_alphar, zggev_obj->threshold);
    EXPECT_NEAR(0.0, zggev_obj->diff_beta, zggev_obj->threshold);
}
TEST_F(zggev_test, zggev3) {
    //EXPECT_NEAR(0.0, zggev_obj->diff_a, zggev_obj->threshold);
    //EXPECT_NEAR(0.0, zggev_obj->diff_b, zggev_obj->threshold);
    EXPECT_NEAR(0.0, zggev_obj->diff_vsl, zggev_obj->threshold);
    EXPECT_NEAR(0.0, zggev_obj->diff_vsr, zggev_obj->threshold);
    EXPECT_NEAR(0.0, zggev_obj->diff_alphar, zggev_obj->threshold);
    EXPECT_NEAR(0.0, zggev_obj->diff_beta, zggev_obj->threshold);
}
TEST_F(zggev_test, zggev4) {
    //EXPECT_NEAR(0.0, zggev_obj->diff_a, zggev_obj->threshold);
    //EXPECT_NEAR(0.0, zggev_obj->diff_b, zggev_obj->threshold);
    EXPECT_NEAR(0.0, zggev_obj->diff_vsl, zggev_obj->threshold);
    EXPECT_NEAR(0.0, zggev_obj->diff_vsr, zggev_obj->threshold);
    EXPECT_NEAR(0.0, zggev_obj->diff_alphar, zggev_obj->threshold);
    EXPECT_NEAR(0.0, zggev_obj->diff_beta, zggev_obj->threshold);
}
