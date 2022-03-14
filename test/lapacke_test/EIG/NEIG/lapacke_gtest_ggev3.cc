#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"
#define LAPACKE_TEST_VERBOSE  (1)

#define ggev3_free() \
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

/* Begin ggev3_float_common_parameters  class definition */
class ggev3_float_parameters{

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

      ggev3_float_parameters (int matrix_layout_i, char jobvsl_i,
						char jobvsr_i, char sort_i, lapack_int n_i );
      ~ggev3_float_parameters ();
};

/* Constructor definition  float_common_parameters */
ggev3_float_parameters:: ggev3_float_parameters (int matrix_layout_i,
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
   printf(" \n ggev3 float: matrix_layout: %d n: %d  jobvsl: %c \
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
       EXPECT_FALSE( true) << "ggev3_float_parameters object: malloc error. Exiting ";
       ggev3_free();
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
    

/* Destructor definition  'ggev3_float_common_parameters' */
ggev3_float_parameters :: ~ggev3_float_parameters ()
{
   /* De-Allocate memory for the input matrices */
   ggev3_free();
} 

//  Test fixture class definition
class sggev3_test  : public  ::testing::Test {
public:
   ggev3_float_parameters  *sggev3_obj;
   void SetUp();  
   void TearDown () { delete sggev3_obj; }
};

void sggev3_test::SetUp()
{

    /* LAPACKE_sggev3 prototype */	
    typedef int (*Fptr_NL_LAPACKE_sggev3) (int matrix_layout, char jobvl,
		char jobvr, lapack_int n, float* a, lapack_int lda, float* b,
		lapack_int ldb, float* alphar, float* alphai, float* beta,
		float* vl, lapack_int ldvl, float* vr, lapack_int ldvr);
				 
    Fptr_NL_LAPACKE_sggev3 SGGEV3;

    sggev3_obj = new  ggev3_float_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].jobvsl,
										 eig_non_sym_paramslist[idx].jobvsr,
										 eig_non_sym_paramslist[idx].sort,
                                         eig_paramslist[idx].n );
                                         
    sggev3_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    sggev3_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sggev3_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sggev3_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sggev3_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SGGEV3 = (Fptr_NL_LAPACKE_sggev3)dlsym(sggev3_obj->hModule, "LAPACKE_sggev3");
    ASSERT_TRUE(SGGEV3 != NULL) << "failed to ppt the Netlib LAPACKE_sggev3 symbol";

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    sggev3_obj->inforef = SGGEV3(   sggev3_obj->matrix_layout,
                                    sggev3_obj->jobvsl,
                                    sggev3_obj->jobvsr,
                                    sggev3_obj->n,
                                    sggev3_obj->aref,
                                    sggev3_obj->lda,
                                    sggev3_obj->bref,
                                    sggev3_obj->ldb,
                                    sggev3_obj->alpharref,
                                    sggev3_obj->alphairef,
                                    sggev3_obj->betaref,
									sggev3_obj->vslref,
									sggev3_obj->ldvsl,
									sggev3_obj->vsrref,
									sggev3_obj->ldvsr
                                    );

    /* Compute libflame's Lapacke o/p  */
    sggev3_obj->inforef =  LAPACKE_sggev3 (   sggev3_obj->matrix_layout,
                                    sggev3_obj->jobvsl,
                                    sggev3_obj->jobvsr,
                                    sggev3_obj->n,
                                    sggev3_obj->a,
                                    sggev3_obj->lda,
                                    sggev3_obj->b,
                                    sggev3_obj->ldb,
                                    sggev3_obj->alphar,
                                    sggev3_obj->alphai,
                                    sggev3_obj->beta,
									sggev3_obj->vsl,
									sggev3_obj->ldvsl,
									sggev3_obj->vsr,
									sggev3_obj->ldvsr
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    
    sggev3_obj->diff_a =  computeDiff_s( (sggev3_obj->lda)*(sggev3_obj->n), 
                sggev3_obj->a, sggev3_obj->aref );

    sggev3_obj->diff_b =  computeDiff_s( (sggev3_obj->ldb)*(sggev3_obj->n), 
                sggev3_obj->b, sggev3_obj->bref );

    sggev3_obj->diff_vsl =  computeDiff_s( (sggev3_obj->ldvsl)*(sggev3_obj->n), 
                sggev3_obj->vsl, sggev3_obj->vslref );

    sggev3_obj->diff_vsr =  computeDiff_s( (sggev3_obj->ldvsr)*(sggev3_obj->n), 
                sggev3_obj->vsr, sggev3_obj->vsrref );

    sggev3_obj->diff_alphar =  computeDiff_s( sggev3_obj->n, 
                sggev3_obj->alphar, sggev3_obj->alpharref );

    sggev3_obj->diff_alphai =  computeDiff_s( sggev3_obj->n, 
                sggev3_obj->alphai, sggev3_obj->alphairef );

    sggev3_obj->diff_beta =  computeDiff_s( sggev3_obj->n, 
                sggev3_obj->beta, sggev3_obj->betaref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n ggev3 float: \n diff_a: %f \n diff_b: %f \n \
diff_vsl: %f \n diff_vsr: %f \n diff_alphar: %f \n\
diff_alphai: %f \n diff_beta: %f \n ",
       sggev3_obj->diff_a, sggev3_obj->diff_b, sggev3_obj->diff_vsl,
	   sggev3_obj->diff_vsr, sggev3_obj->diff_alphar,
       sggev3_obj->diff_alphai, sggev3_obj->diff_beta );
#endif
}

TEST_F(sggev3_test, sggev31) {
    //EXPECT_NEAR(0.0, sggev3_obj->diff_a, sggev3_obj->threshold);
    //EXPECT_NEAR(0.0, sggev3_obj->diff_b, sggev3_obj->threshold);
    EXPECT_NEAR(0.0, sggev3_obj->diff_vsl, sggev3_obj->threshold);
    EXPECT_NEAR(0.0, sggev3_obj->diff_vsr, sggev3_obj->threshold);
    EXPECT_NEAR(0.0, sggev3_obj->diff_alphar, sggev3_obj->threshold);
    EXPECT_NEAR(0.0, sggev3_obj->diff_alphai, sggev3_obj->threshold);
    EXPECT_NEAR(0.0, sggev3_obj->diff_beta, sggev3_obj->threshold);
}

TEST_F(sggev3_test, sggev32) {
    //EXPECT_NEAR(0.0, sggev3_obj->diff_a, sggev3_obj->threshold);
    //EXPECT_NEAR(0.0, sggev3_obj->diff_b, sggev3_obj->threshold);
    EXPECT_NEAR(0.0, sggev3_obj->diff_vsl, sggev3_obj->threshold);
    EXPECT_NEAR(0.0, sggev3_obj->diff_vsr, sggev3_obj->threshold);
    EXPECT_NEAR(0.0, sggev3_obj->diff_alphar, sggev3_obj->threshold);
    EXPECT_NEAR(0.0, sggev3_obj->diff_alphai, sggev3_obj->threshold);
    EXPECT_NEAR(0.0, sggev3_obj->diff_beta, sggev3_obj->threshold);
}
TEST_F(sggev3_test, sggev33) {
    //EXPECT_NEAR(0.0, sggev3_obj->diff_a, sggev3_obj->threshold);
    //EXPECT_NEAR(0.0, sggev3_obj->diff_b, sggev3_obj->threshold);
    EXPECT_NEAR(0.0, sggev3_obj->diff_vsl, sggev3_obj->threshold);
    EXPECT_NEAR(0.0, sggev3_obj->diff_vsr, sggev3_obj->threshold);
    EXPECT_NEAR(0.0, sggev3_obj->diff_alphar, sggev3_obj->threshold);
    EXPECT_NEAR(0.0, sggev3_obj->diff_alphai, sggev3_obj->threshold);
    EXPECT_NEAR(0.0, sggev3_obj->diff_beta, sggev3_obj->threshold);
}
TEST_F(sggev3_test, sggev34) {
    //EXPECT_NEAR(0.0, sggev3_obj->diff_a, sggev3_obj->threshold);
    //EXPECT_NEAR(0.0, sggev3_obj->diff_b, sggev3_obj->threshold);
    EXPECT_NEAR(0.0, sggev3_obj->diff_vsl, sggev3_obj->threshold);
    EXPECT_NEAR(0.0, sggev3_obj->diff_vsr, sggev3_obj->threshold);
    EXPECT_NEAR(0.0, sggev3_obj->diff_alphar, sggev3_obj->threshold);
    EXPECT_NEAR(0.0, sggev3_obj->diff_alphai, sggev3_obj->threshold);
    EXPECT_NEAR(0.0, sggev3_obj->diff_beta, sggev3_obj->threshold);
}


/* Begin ggev3_double_common_parameters  class definition */
class ggev3_double_parameters{

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

      ggev3_double_parameters (int matrix_layout_i, char jobvsl_i,
						char jobvsr_i, char sort_i, lapack_int n_i );
      ~ggev3_double_parameters ();
};

/* Constructor definition  double_common_parameters */
ggev3_double_parameters:: ggev3_double_parameters (int matrix_layout_i,
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
   printf(" \n ggev3 double: matrix_layout: %d n: %d  jobvsl: %c \
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
       EXPECT_FALSE( true) << "ggev3_double_parameters object: malloc error. Exiting ";
       ggev3_free();
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
    

/* Destructor definition  'ggev3_double_common_parameters' */
ggev3_double_parameters :: ~ggev3_double_parameters ()
{
   /* De-Allocate memory for the input matrices */
   ggev3_free();
} 

//  Test fixture class definition
class dggev3_test  : public  ::testing::Test {
public:
   ggev3_double_parameters  *dggev3_obj;
   void SetUp();  
   void TearDown () { delete dggev3_obj; }
};

void dggev3_test::SetUp()
{

    /* LAPACKE_dggev3 prototype */	
    typedef int (*Fptr_NL_LAPACKE_dggev3) (int matrix_layout, char jobvl,
		char jobvr, lapack_int n, double* a, lapack_int lda, double* b,
		lapack_int ldb, double* alphar, double* alphai, double* beta,
		double* vl, lapack_int ldvl, double* vr, lapack_int ldvr);
				 
    Fptr_NL_LAPACKE_dggev3 DGGEV3;

    dggev3_obj = new  ggev3_double_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].jobvsl,
										 eig_non_sym_paramslist[idx].jobvsr,
										 eig_non_sym_paramslist[idx].sort,
                                         eig_paramslist[idx].n );
                                         
    dggev3_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    dggev3_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dggev3_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dggev3_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dggev3_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DGGEV3 = (Fptr_NL_LAPACKE_dggev3)dlsym(dggev3_obj->hModule, "LAPACKE_dggev3");
    ASSERT_TRUE(DGGEV3 != NULL) << "failed to ppt the Netlib LAPACKE_dggev3 symbol";

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    dggev3_obj->inforef = DGGEV3(   dggev3_obj->matrix_layout,
                                    dggev3_obj->jobvsl,
                                    dggev3_obj->jobvsr,
                                    dggev3_obj->n,
                                    dggev3_obj->aref,
                                    dggev3_obj->lda,
                                    dggev3_obj->bref,
                                    dggev3_obj->ldb,
                                    dggev3_obj->alpharref,
                                    dggev3_obj->alphairef,
                                    dggev3_obj->betaref,
									dggev3_obj->vslref,
									dggev3_obj->ldvsl,
									dggev3_obj->vsrref,
									dggev3_obj->ldvsr
                                    );

    /* Compute libflame's Lapacke o/p  */
    dggev3_obj->inforef =  LAPACKE_dggev3 (   dggev3_obj->matrix_layout,
                                    dggev3_obj->jobvsl,
                                    dggev3_obj->jobvsr,
                                    dggev3_obj->n,
                                    dggev3_obj->a,
                                    dggev3_obj->lda,
                                    dggev3_obj->b,
                                    dggev3_obj->ldb,
                                    dggev3_obj->alphar,
                                    dggev3_obj->alphai,
                                    dggev3_obj->beta,
									dggev3_obj->vsl,
									dggev3_obj->ldvsl,
									dggev3_obj->vsr,
									dggev3_obj->ldvsr
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    
    dggev3_obj->diff_a =  computeDiff_d( (dggev3_obj->lda)*(dggev3_obj->n), 
                dggev3_obj->a, dggev3_obj->aref );

    dggev3_obj->diff_b =  computeDiff_d( (dggev3_obj->ldb)*(dggev3_obj->n), 
                dggev3_obj->b, dggev3_obj->bref );

    dggev3_obj->diff_vsl =  computeDiff_d( (dggev3_obj->ldvsl)*(dggev3_obj->n), 
                dggev3_obj->vsl, dggev3_obj->vslref );

    dggev3_obj->diff_vsr =  computeDiff_d( (dggev3_obj->ldvsr)*(dggev3_obj->n), 
                dggev3_obj->vsr, dggev3_obj->vsrref );

    dggev3_obj->diff_alphar =  computeDiff_d( dggev3_obj->n, 
                dggev3_obj->alphar, dggev3_obj->alpharref );

    dggev3_obj->diff_alphai =  computeDiff_d( dggev3_obj->n, 
                dggev3_obj->alphai, dggev3_obj->alphairef );

    dggev3_obj->diff_beta = computeDiff_d ( dggev3_obj->n, 
                dggev3_obj->beta, dggev3_obj->betaref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n ggev3 double: \n diff_a: %f \n diff_b: %f \n \
diff_vsl: %f \n diff_vsr: %f \n diff_alphar: %f \n\
diff_alphai: %f \n diff_beta: %f \n ",
       dggev3_obj->diff_a, dggev3_obj->diff_b, dggev3_obj->diff_vsl,
	   dggev3_obj->diff_vsr, dggev3_obj->diff_alphar,
       dggev3_obj->diff_alphai, dggev3_obj->diff_beta );
#endif
}

TEST_F(dggev3_test, dggev31) {
    //EXPECT_NEAR(0.0, dggev3_obj->diff_a, dggev3_obj->threshold);
    //EXPECT_NEAR(0.0, dggev3_obj->diff_b, dggev3_obj->threshold);
    EXPECT_NEAR(0.0, dggev3_obj->diff_vsl, dggev3_obj->threshold);
    EXPECT_NEAR(0.0, dggev3_obj->diff_vsr, dggev3_obj->threshold);
    EXPECT_NEAR(0.0, dggev3_obj->diff_alphar, dggev3_obj->threshold);
    EXPECT_NEAR(0.0, dggev3_obj->diff_alphai, dggev3_obj->threshold);
    EXPECT_NEAR(0.0, dggev3_obj->diff_beta, dggev3_obj->threshold);
}

TEST_F(dggev3_test, dggev32) {
    //EXPECT_NEAR(0.0, dggev3_obj->diff_a, dggev3_obj->threshold);
    //EXPECT_NEAR(0.0, dggev3_obj->diff_b, dggev3_obj->threshold);
    EXPECT_NEAR(0.0, dggev3_obj->diff_vsl, dggev3_obj->threshold);
    EXPECT_NEAR(0.0, dggev3_obj->diff_vsr, dggev3_obj->threshold);
    EXPECT_NEAR(0.0, dggev3_obj->diff_alphar, dggev3_obj->threshold);
    EXPECT_NEAR(0.0, dggev3_obj->diff_alphai, dggev3_obj->threshold);
    EXPECT_NEAR(0.0, dggev3_obj->diff_beta, dggev3_obj->threshold);
}
TEST_F(dggev3_test, dggev33) {
    //EXPECT_NEAR(0.0, dggev3_obj->diff_a, dggev3_obj->threshold);
    //EXPECT_NEAR(0.0, dggev3_obj->diff_b, dggev3_obj->threshold);
    EXPECT_NEAR(0.0, dggev3_obj->diff_vsl, dggev3_obj->threshold);
    EXPECT_NEAR(0.0, dggev3_obj->diff_vsr, dggev3_obj->threshold);
    EXPECT_NEAR(0.0, dggev3_obj->diff_alphar, dggev3_obj->threshold);
    EXPECT_NEAR(0.0, dggev3_obj->diff_alphai, dggev3_obj->threshold);
    EXPECT_NEAR(0.0, dggev3_obj->diff_beta, dggev3_obj->threshold);
}
TEST_F(dggev3_test, dggev34) {
    //EXPECT_NEAR(0.0, dggev3_obj->diff_a, dggev3_obj->threshold);
    //EXPECT_NEAR(0.0, dggev3_obj->diff_b, dggev3_obj->threshold);
    EXPECT_NEAR(0.0, dggev3_obj->diff_vsl, dggev3_obj->threshold);
    EXPECT_NEAR(0.0, dggev3_obj->diff_vsr, dggev3_obj->threshold);
    EXPECT_NEAR(0.0, dggev3_obj->diff_alphar, dggev3_obj->threshold);
    EXPECT_NEAR(0.0, dggev3_obj->diff_alphai, dggev3_obj->threshold);
    EXPECT_NEAR(0.0, dggev3_obj->diff_beta, dggev3_obj->threshold);
}

/* Begin ggev3_scomplex_common_parameters  class definition */
class ggev3_scomplex_parameters{

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

      ggev3_scomplex_parameters (int matrix_layout_i, char jobvsl_i,
						char jobvsr_i, char sort_i, lapack_int n_i );
      ~ggev3_scomplex_parameters ();
};

/* Constructor definition  lapack_complex_float_common_parameters */
ggev3_scomplex_parameters:: ggev3_scomplex_parameters (int matrix_layout_i,
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
   printf(" \n ggev3 lapack_complex_float: matrix_layout: %d n: %d  jobvsl: %c \
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
       EXPECT_FALSE( true) << "ggev3_scomplex_parameters object: malloc error. Exiting ";
       ggev3_free();
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
    

/* Destructor definition  'ggev3_scomplex_common_parameters' */
ggev3_scomplex_parameters :: ~ggev3_scomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   ggev3_free();
} 

//  Test fixture class definition
class cggev3_test  : public  ::testing::Test {
public:
   ggev3_scomplex_parameters  *cggev3_obj;
   void SetUp();  
   void TearDown () { delete cggev3_obj; }
};

void cggev3_test::SetUp()
{

    /* LAPACKE_cggev3 prototype */	
    typedef int (*Fptr_NL_LAPACKE_cggev3) (int matrix_layout, char jobvl,
		char jobvr, lapack_int n, lapack_complex_float* a, lapack_int lda,
		lapack_complex_float* b, lapack_int ldb, lapack_complex_float* alpha,
		lapack_complex_float* beta, lapack_complex_float* vl, lapack_int ldvl,
		lapack_complex_float* vr, lapack_int ldvr);
				 
    Fptr_NL_LAPACKE_cggev3 CGGEV3;

    cggev3_obj = new  ggev3_scomplex_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].jobvsl,
										 eig_non_sym_paramslist[idx].jobvsr,
										 eig_non_sym_paramslist[idx].sort,
                                         eig_paramslist[idx].n );
                                         
    cggev3_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    cggev3_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cggev3_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cggev3_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cggev3_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CGGEV3 = (Fptr_NL_LAPACKE_cggev3)dlsym(cggev3_obj->hModule, "LAPACKE_cggev3");
    ASSERT_TRUE(CGGEV3 != NULL) << "failed to ppt the Netlib LAPACKE_cggev3 symbol";

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    cggev3_obj->inforef = CGGEV3(   cggev3_obj->matrix_layout,
                                    cggev3_obj->jobvsl,
                                    cggev3_obj->jobvsr,
                                    cggev3_obj->n,
                                    cggev3_obj->aref,
                                    cggev3_obj->lda,
                                    cggev3_obj->bref,
                                    cggev3_obj->ldb,
                                    cggev3_obj->alpharref,
                                    cggev3_obj->betaref,
									cggev3_obj->vslref,
									cggev3_obj->ldvsl,
									cggev3_obj->vsrref,
									cggev3_obj->ldvsr
                                    );

    /* Compute libflame's Lapacke o/p  */
    cggev3_obj->inforef =  LAPACKE_cggev3 (   cggev3_obj->matrix_layout,
                                    cggev3_obj->jobvsl,
                                    cggev3_obj->jobvsr,
                                    cggev3_obj->n,
                                    cggev3_obj->a,
                                    cggev3_obj->lda,
                                    cggev3_obj->b,
                                    cggev3_obj->ldb,
                                    cggev3_obj->alphar,
                                    cggev3_obj->beta,
									cggev3_obj->vsl,
									cggev3_obj->ldvsl,
									cggev3_obj->vsr,
									cggev3_obj->ldvsr
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    
    cggev3_obj->diff_a =  computeDiff_c( (cggev3_obj->lda)*(cggev3_obj->n), 
                cggev3_obj->a, cggev3_obj->aref );

    cggev3_obj->diff_b =  computeDiff_c( (cggev3_obj->ldb)*(cggev3_obj->n), 
                cggev3_obj->b, cggev3_obj->bref );

    cggev3_obj->diff_vsl =  computeDiff_c( (cggev3_obj->ldvsl)*(cggev3_obj->n), 
                cggev3_obj->vsl, cggev3_obj->vslref );

    cggev3_obj->diff_vsr =  computeDiff_c( (cggev3_obj->ldvsr)*(cggev3_obj->n), 
                cggev3_obj->vsr, cggev3_obj->vsrref );

    cggev3_obj->diff_alphar =  computeDiff_c( cggev3_obj->n, 
                cggev3_obj->alphar, cggev3_obj->alpharref );

    cggev3_obj->diff_beta =  computeDiff_c( cggev3_obj->n, 
                cggev3_obj->beta, cggev3_obj->betaref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n ggev3 lapack_complex_float: \n diff_a: %f \n diff_b: %f \n \
diff_vsl: %f \n diff_vsr: %f \n diff_alphar: %f \n\
 \n diff_beta: %f \n ",
       cggev3_obj->diff_a, cggev3_obj->diff_b, cggev3_obj->diff_vsl,
	   cggev3_obj->diff_vsr, cggev3_obj->diff_alphar,
        cggev3_obj->diff_beta );
#endif
}

TEST_F(cggev3_test, cggev31) {
    //EXPECT_NEAR(0.0, cggev3_obj->diff_a, cggev3_obj->threshold);
    //EXPECT_NEAR(0.0, cggev3_obj->diff_b, cggev3_obj->threshold);
    EXPECT_NEAR(0.0, cggev3_obj->diff_vsl, cggev3_obj->threshold);
    EXPECT_NEAR(0.0, cggev3_obj->diff_vsr, cggev3_obj->threshold);
    EXPECT_NEAR(0.0, cggev3_obj->diff_alphar, cggev3_obj->threshold);
    EXPECT_NEAR(0.0, cggev3_obj->diff_beta, cggev3_obj->threshold);
}

TEST_F(cggev3_test, cggev32) {
    //EXPECT_NEAR(0.0, cggev3_obj->diff_a, cggev3_obj->threshold);
    //EXPECT_NEAR(0.0, cggev3_obj->diff_b, cggev3_obj->threshold);
    EXPECT_NEAR(0.0, cggev3_obj->diff_vsl, cggev3_obj->threshold);
    EXPECT_NEAR(0.0, cggev3_obj->diff_vsr, cggev3_obj->threshold);
    EXPECT_NEAR(0.0, cggev3_obj->diff_alphar, cggev3_obj->threshold);
    EXPECT_NEAR(0.0, cggev3_obj->diff_beta, cggev3_obj->threshold);
}
TEST_F(cggev3_test, cggev33) {
    //EXPECT_NEAR(0.0, cggev3_obj->diff_a, cggev3_obj->threshold);
    //EXPECT_NEAR(0.0, cggev3_obj->diff_b, cggev3_obj->threshold);
    EXPECT_NEAR(0.0, cggev3_obj->diff_vsl, cggev3_obj->threshold);
    EXPECT_NEAR(0.0, cggev3_obj->diff_vsr, cggev3_obj->threshold);
    EXPECT_NEAR(0.0, cggev3_obj->diff_alphar, cggev3_obj->threshold);
    EXPECT_NEAR(0.0, cggev3_obj->diff_beta, cggev3_obj->threshold);
}
TEST_F(cggev3_test, cggev34) {
    //EXPECT_NEAR(0.0, cggev3_obj->diff_a, cggev3_obj->threshold);
    //EXPECT_NEAR(0.0, cggev3_obj->diff_b, cggev3_obj->threshold);
    EXPECT_NEAR(0.0, cggev3_obj->diff_vsl, cggev3_obj->threshold);
    EXPECT_NEAR(0.0, cggev3_obj->diff_vsr, cggev3_obj->threshold);
    EXPECT_NEAR(0.0, cggev3_obj->diff_alphar, cggev3_obj->threshold);
    EXPECT_NEAR(0.0, cggev3_obj->diff_beta, cggev3_obj->threshold);
}

/* Begin ggev3_dcomplex_common_parameters  class definition */
class ggev3_dcomplex_parameters{

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

      ggev3_dcomplex_parameters (int matrix_layout_i, char jobvsl_i,
						char jobvsr_i, char sort_i, lapack_int n_i );
      ~ggev3_dcomplex_parameters ();
};

/* Constructor definition  lapack_complex_double_common_parameters */
ggev3_dcomplex_parameters:: ggev3_dcomplex_parameters (int matrix_layout_i,
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
   printf(" \n ggev3 lapack_complex_double: matrix_layout: %d n: %d  jobvsl: %c \
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
       EXPECT_FALSE( true) << "ggev3_dcomplex_parameters object: malloc error. Exiting ";
       ggev3_free();
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
    

/* Destructor definition  'ggev3_dcomplex_common_parameters' */
ggev3_dcomplex_parameters :: ~ggev3_dcomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   ggev3_free();
} 

//  Test fixture class definition
class zggev3_test  : public  ::testing::Test {
public:
   ggev3_dcomplex_parameters  *zggev3_obj;
   void SetUp();  
   void TearDown () { delete zggev3_obj; }
};

void zggev3_test::SetUp()
{

    /* LAPACKE_zggev3 prototype */	
    typedef int (*Fptr_NL_LAPACKE_zggev3) (int matrix_layout, char jobvl,
		char jobvr, lapack_int n, lapack_complex_double* a, lapack_int lda,
		lapack_complex_double* b, lapack_int ldb, lapack_complex_double* alpha,
		lapack_complex_double* beta, lapack_complex_double* vl, lapack_int ldvl,
		lapack_complex_double* vr, lapack_int ldvr);
				 
    Fptr_NL_LAPACKE_zggev3 ZGGEV3;

    zggev3_obj = new  ggev3_dcomplex_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].jobvsl,
										 eig_non_sym_paramslist[idx].jobvsr,
										 eig_non_sym_paramslist[idx].sort,
                                         eig_paramslist[idx].n );
                                         
    zggev3_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    zggev3_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zggev3_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zggev3_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zggev3_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZGGEV3 = (Fptr_NL_LAPACKE_zggev3)dlsym(zggev3_obj->hModule, "LAPACKE_zggev3");
    ASSERT_TRUE(ZGGEV3 != NULL) << "failed to ppt the Netlib LAPACKE_zggev3 symbol";

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    zggev3_obj->inforef = ZGGEV3(   zggev3_obj->matrix_layout,
                                    zggev3_obj->jobvsl,
                                    zggev3_obj->jobvsr,
                                    zggev3_obj->n,
                                    zggev3_obj->aref,
                                    zggev3_obj->lda,
                                    zggev3_obj->bref,
                                    zggev3_obj->ldb,
                                    zggev3_obj->alpharref,
                                    zggev3_obj->betaref,
									zggev3_obj->vslref,
									zggev3_obj->ldvsl,
									zggev3_obj->vsrref,
									zggev3_obj->ldvsr
                                    );

    /* Compute libflame's Lapacke o/p  */
    zggev3_obj->inforef =  LAPACKE_zggev3 (   zggev3_obj->matrix_layout,
                                    zggev3_obj->jobvsl,
                                    zggev3_obj->jobvsr,
                                    zggev3_obj->n,
                                    zggev3_obj->a,
                                    zggev3_obj->lda,
                                    zggev3_obj->b,
                                    zggev3_obj->ldb,
                                    zggev3_obj->alphar,
                                    zggev3_obj->beta,
									zggev3_obj->vsl,
									zggev3_obj->ldvsl,
									zggev3_obj->vsr,
									zggev3_obj->ldvsr
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    
    zggev3_obj->diff_a =  computeDiff_z( (zggev3_obj->lda)*(zggev3_obj->n), 
                zggev3_obj->a, zggev3_obj->aref );

    zggev3_obj->diff_b =  computeDiff_z( (zggev3_obj->ldb)*(zggev3_obj->n), 
                zggev3_obj->b, zggev3_obj->bref );

    zggev3_obj->diff_vsl =  computeDiff_z( (zggev3_obj->ldvsl)*(zggev3_obj->n), 
                zggev3_obj->vsl, zggev3_obj->vslref );

    zggev3_obj->diff_vsr =  computeDiff_z( (zggev3_obj->ldvsr)*(zggev3_obj->n), 
                zggev3_obj->vsr, zggev3_obj->vsrref );

    zggev3_obj->diff_alphar =  computeDiff_z( zggev3_obj->n, 
                zggev3_obj->alphar, zggev3_obj->alpharref );

    zggev3_obj->diff_beta =  computeDiff_z( zggev3_obj->n, 
                zggev3_obj->beta, zggev3_obj->betaref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n ggev3 lapack_complex_double: \n diff_a: %f \n diff_b: %f \n \
diff_vsl: %f \n diff_vsr: %f \n diff_alphar: %f \n\
 \n diff_beta: %f \n ",
       zggev3_obj->diff_a, zggev3_obj->diff_b, zggev3_obj->diff_vsl,
	   zggev3_obj->diff_vsr, zggev3_obj->diff_alphar,
        zggev3_obj->diff_beta );
#endif
}

TEST_F(zggev3_test, zggev31) {
    //EXPECT_NEAR(0.0, zggev3_obj->diff_a, zggev3_obj->threshold);
    //EXPECT_NEAR(0.0, zggev3_obj->diff_b, zggev3_obj->threshold);
    EXPECT_NEAR(0.0, zggev3_obj->diff_vsl, zggev3_obj->threshold);
    EXPECT_NEAR(0.0, zggev3_obj->diff_vsr, zggev3_obj->threshold);
    EXPECT_NEAR(0.0, zggev3_obj->diff_alphar, zggev3_obj->threshold);
    EXPECT_NEAR(0.0, zggev3_obj->diff_beta, zggev3_obj->threshold);
}

TEST_F(zggev3_test, zggev32) {
    //EXPECT_NEAR(0.0, zggev3_obj->diff_a, zggev3_obj->threshold);
    //EXPECT_NEAR(0.0, zggev3_obj->diff_b, zggev3_obj->threshold);
    EXPECT_NEAR(0.0, zggev3_obj->diff_vsl, zggev3_obj->threshold);
    EXPECT_NEAR(0.0, zggev3_obj->diff_vsr, zggev3_obj->threshold);
    EXPECT_NEAR(0.0, zggev3_obj->diff_alphar, zggev3_obj->threshold);
    EXPECT_NEAR(0.0, zggev3_obj->diff_beta, zggev3_obj->threshold);
}
TEST_F(zggev3_test, zggev33) {
    //EXPECT_NEAR(0.0, zggev3_obj->diff_a, zggev3_obj->threshold);
    //EXPECT_NEAR(0.0, zggev3_obj->diff_b, zggev3_obj->threshold);
    EXPECT_NEAR(0.0, zggev3_obj->diff_vsl, zggev3_obj->threshold);
    EXPECT_NEAR(0.0, zggev3_obj->diff_vsr, zggev3_obj->threshold);
    EXPECT_NEAR(0.0, zggev3_obj->diff_alphar, zggev3_obj->threshold);
    EXPECT_NEAR(0.0, zggev3_obj->diff_beta, zggev3_obj->threshold);
}
TEST_F(zggev3_test, zggev34) {
    //EXPECT_NEAR(0.0, zggev3_obj->diff_a, zggev3_obj->threshold);
    //EXPECT_NEAR(0.0, zggev3_obj->diff_b, zggev3_obj->threshold);
    EXPECT_NEAR(0.0, zggev3_obj->diff_vsl, zggev3_obj->threshold);
    EXPECT_NEAR(0.0, zggev3_obj->diff_vsr, zggev3_obj->threshold);
    EXPECT_NEAR(0.0, zggev3_obj->diff_alphar, zggev3_obj->threshold);
    EXPECT_NEAR(0.0, zggev3_obj->diff_beta, zggev3_obj->threshold);
}
