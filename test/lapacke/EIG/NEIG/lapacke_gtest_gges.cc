#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define gges_free() \
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

/* Begin gges_float_common_parameters  class definition */
class gges_float_parameters{

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

      gges_float_parameters (int matrix_layout_i, char jobvsl_i,
						char jobvsr_i, char sort_i, lapack_int n_i );
      ~gges_float_parameters ();
};

/* Constructor definition  float_common_parameters */
gges_float_parameters:: gges_float_parameters (int matrix_layout_i,
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
   printf(" \n gges float: matrix_layout: %d n: %d  jobvsl: %c \
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
       EXPECT_FALSE( true) << "gges_float_parameters object: malloc error. Exiting ";
       gges_free();
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
    

/* Destructor definition  'gges_float_common_parameters' */
gges_float_parameters :: ~gges_float_parameters ()
{
   /* De-Allocate memory for the input matrices */
   gges_free();
} 

//  Test fixture class definition
class sgges_test  : public  ::testing::Test {
public:
   gges_float_parameters  *sgges_obj;
   void SetUp();  
   void TearDown () { delete sgges_obj; }
};

void sgges_test::SetUp()
{

    /* LAPACKE_sgges prototype */	
    typedef int (*Fptr_NL_LAPACKE_sgges) (int matrix_layout, char jobvsl,
		char jobvsr, char sort, LAPACK_S_SELECT3 select, lapack_int n, 
		float* a, lapack_int lda, float* b, lapack_int ldb, lapack_int* sdim,
		float* alphar, float* alphai, float* beta, float* vsl,
		lapack_int ldvsl, float* vsr, lapack_int ldvsr);
				 
    Fptr_NL_LAPACKE_sgges SGGES;

    sgges_obj = new  gges_float_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].jobvsl,
										 eig_non_sym_paramslist[idx].jobvsr,
										 eig_non_sym_paramslist[idx].sort,
                                         eig_paramslist[idx].n );
                                         
    sgges_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    sgges_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgges_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgges_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgges_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SGGES = (Fptr_NL_LAPACKE_sgges)dlsym(sgges_obj->hModule, "LAPACKE_sgges");
    ASSERT_TRUE(SGGES != NULL) << "failed to ppt the Netlib LAPACKE_sgges symbol";

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    sgges_obj->inforef = SGGES(   sgges_obj->matrix_layout,
                                    sgges_obj->jobvsl,
                                    sgges_obj->jobvsr,
                                    sgges_obj->sort,
									(LAPACK_S_SELECT3)sgges_obj->selectref,
                                    sgges_obj->n,
                                    sgges_obj->aref,
                                    sgges_obj->lda,
                                    sgges_obj->bref,
                                    sgges_obj->ldb,
                                    &sgges_obj->sdimref,
                                    sgges_obj->alpharref,
                                    sgges_obj->alphairef,
                                    sgges_obj->betaref,
									sgges_obj->vslref,
									sgges_obj->ldvsl,
									sgges_obj->vsrref,
									sgges_obj->ldvsr
                                    );

    /* Compute libflame's Lapacke o/p  */
    sgges_obj->inforef =  LAPACKE_sgges (   sgges_obj->matrix_layout,
                                    sgges_obj->jobvsl,
                                    sgges_obj->jobvsr,
                                    sgges_obj->sort,
									(LAPACK_S_SELECT3)sgges_obj->select,
                                    sgges_obj->n,
                                    sgges_obj->a,
                                    sgges_obj->lda,
                                    sgges_obj->b,
                                    sgges_obj->ldb,
                                    &sgges_obj->sdim,
                                    sgges_obj->alphar,
                                    sgges_obj->alphai,
                                    sgges_obj->beta,
									sgges_obj->vsl,
									sgges_obj->ldvsl,
									sgges_obj->vsr,
									sgges_obj->ldvsr
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    
    sgges_obj->diff_a =  computeDiff_s( (sgges_obj->lda)*(sgges_obj->n), 
                sgges_obj->a, sgges_obj->aref );

    sgges_obj->diff_b =  computeDiff_s( (sgges_obj->ldb)*(sgges_obj->n), 
                sgges_obj->b, sgges_obj->bref );

    sgges_obj->diff_vsl =  computeDiff_s( (sgges_obj->ldvsl)*(sgges_obj->n), 
                sgges_obj->vsl, sgges_obj->vslref );

    sgges_obj->diff_vsr =  computeDiff_s( (sgges_obj->ldvsr)*(sgges_obj->n), 
                sgges_obj->vsr, sgges_obj->vsrref );

    sgges_obj->diff_alphar =  computeDiff_s( sgges_obj->n, 
                sgges_obj->alphar, sgges_obj->alpharref );

    sgges_obj->diff_alphai =  computeDiff_s( sgges_obj->n, 
                sgges_obj->alphai, sgges_obj->alphairef );

    sgges_obj->diff_beta =  computeDiff_s( sgges_obj->n, 
                sgges_obj->beta, sgges_obj->betaref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n gges float: \n diff_a: %f \n diff_b: %f \n \
diff_vsl: %f \n diff_vsr: %f \n diff_alphar: %f \n\
diff_alphai: %f \n diff_beta: %f \n ",
       sgges_obj->diff_a, sgges_obj->diff_b, sgges_obj->diff_vsl,
	   sgges_obj->diff_vsr, sgges_obj->diff_alphar,
       sgges_obj->diff_alphai, sgges_obj->diff_beta );
#endif
}

TEST_F(sgges_test, sgges1) {
    EXPECT_NEAR(0.0, sgges_obj->diff_vsl, sgges_obj->threshold);
    EXPECT_NEAR(0.0, sgges_obj->diff_vsr, sgges_obj->threshold);
    EXPECT_NEAR(0.0, sgges_obj->diff_alphar, sgges_obj->threshold);
    EXPECT_NEAR(0.0, sgges_obj->diff_alphai, sgges_obj->threshold);
    EXPECT_NEAR(0.0, sgges_obj->diff_beta, sgges_obj->threshold);
    EXPECT_EQ(sgges_obj->sdim, sgges_obj->sdimref);
}

TEST_F(sgges_test, sgges2) {
    EXPECT_NEAR(0.0, sgges_obj->diff_vsl, sgges_obj->threshold);
    EXPECT_NEAR(0.0, sgges_obj->diff_vsr, sgges_obj->threshold);
    EXPECT_NEAR(0.0, sgges_obj->diff_alphar, sgges_obj->threshold);
    EXPECT_NEAR(0.0, sgges_obj->diff_alphai, sgges_obj->threshold);
    EXPECT_NEAR(0.0, sgges_obj->diff_beta, sgges_obj->threshold);
    EXPECT_EQ(sgges_obj->sdim, sgges_obj->sdimref);
}
TEST_F(sgges_test, sgges3) {
    EXPECT_NEAR(0.0, sgges_obj->diff_vsl, sgges_obj->threshold);
    EXPECT_NEAR(0.0, sgges_obj->diff_vsr, sgges_obj->threshold);
    EXPECT_NEAR(0.0, sgges_obj->diff_alphar, sgges_obj->threshold);
    EXPECT_NEAR(0.0, sgges_obj->diff_alphai, sgges_obj->threshold);
    EXPECT_NEAR(0.0, sgges_obj->diff_beta, sgges_obj->threshold);
    EXPECT_EQ(sgges_obj->sdim, sgges_obj->sdimref);
}
TEST_F(sgges_test, sgges4) {
    EXPECT_NEAR(0.0, sgges_obj->diff_vsl, sgges_obj->threshold);
    EXPECT_NEAR(0.0, sgges_obj->diff_vsr, sgges_obj->threshold);
    EXPECT_NEAR(0.0, sgges_obj->diff_alphar, sgges_obj->threshold);
    EXPECT_NEAR(0.0, sgges_obj->diff_alphai, sgges_obj->threshold);
    EXPECT_NEAR(0.0, sgges_obj->diff_beta, sgges_obj->threshold);
    EXPECT_EQ(sgges_obj->sdim, sgges_obj->sdimref);
}


/* Begin gges_double_common_parameters  class definition */
class gges_double_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    double diff_a, diff_b;
    double diff_alphai, diff_alphar, diff_beta, diff_vsl, diff_vsr;

    double threshold;
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

      gges_double_parameters (int matrix_layout_i, char jobvsl_i,
						char jobvsr_i, char sort_i, lapack_int n_i );
      ~gges_double_parameters ();
};

/* Constructor definition  double_common_parameters */
gges_double_parameters:: gges_double_parameters (int matrix_layout_i,
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
   printf(" \n gges double: matrix_layout: %d n: %d  jobvsl: %c \
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
       EXPECT_FALSE( true) << "gges_double_parameters object: malloc error. Exiting ";
       gges_free();
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
    

/* Destructor definition  'gges_double_common_parameters' */
gges_double_parameters :: ~gges_double_parameters ()
{
   /* De-Allocate memory for the input matrices */
   gges_free();
} 

//  Test fixture class definition
class dgges_test  : public  ::testing::Test {
public:
   gges_double_parameters  *dgges_obj;
   void SetUp();  
   void TearDown () { delete dgges_obj; }
};

void dgges_test::SetUp()
{

    /* LAPACKE_dgges prototype */	
    typedef int (*Fptr_NL_LAPACKE_dgges) (int matrix_layout, char jobvsl,
		char jobvsr, char sort, LAPACK_D_SELECT3 select, lapack_int n, 
		double* a, lapack_int lda, double* b, lapack_int ldb, lapack_int* sdim,
		double* alphar, double* alphai, double* beta, double* vsl,
		lapack_int ldvsl, double* vsr, lapack_int ldvsr);
				 
    Fptr_NL_LAPACKE_dgges DGGES;

    dgges_obj = new  gges_double_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].jobvsl,
										 eig_non_sym_paramslist[idx].jobvsr,
										 eig_non_sym_paramslist[idx].sort,
                                         eig_paramslist[idx].n );
                                         
    dgges_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    dgges_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgges_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgges_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgges_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DGGES = (Fptr_NL_LAPACKE_dgges)dlsym(dgges_obj->hModule, "LAPACKE_dgges");
    ASSERT_TRUE(DGGES != NULL) << "failed to ppt the Netlib LAPACKE_dgges symbol";

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    dgges_obj->inforef = DGGES(   dgges_obj->matrix_layout,
                                    dgges_obj->jobvsl,
                                    dgges_obj->jobvsr,
                                    dgges_obj->sort,
									(LAPACK_D_SELECT3)dgges_obj->selectref,
                                    dgges_obj->n,
                                    dgges_obj->aref,
                                    dgges_obj->lda,
                                    dgges_obj->bref,
                                    dgges_obj->ldb,
                                    &dgges_obj->sdimref,
                                    dgges_obj->alpharref,
                                    dgges_obj->alphairef,
                                    dgges_obj->betaref,
									dgges_obj->vslref,
									dgges_obj->ldvsl,
									dgges_obj->vsrref,
									dgges_obj->ldvsr
                                    );

    /* Compute libflame's Lapacke o/p  */
    dgges_obj->inforef =  LAPACKE_dgges (   dgges_obj->matrix_layout,
                                    dgges_obj->jobvsl,
                                    dgges_obj->jobvsr,
                                    dgges_obj->sort,
									(LAPACK_D_SELECT3)dgges_obj->select,
                                    dgges_obj->n,
                                    dgges_obj->a,
                                    dgges_obj->lda,
                                    dgges_obj->b,
                                    dgges_obj->ldb,
                                    &dgges_obj->sdim,
                                    dgges_obj->alphar,
                                    dgges_obj->alphai,
                                    dgges_obj->beta,
									dgges_obj->vsl,
									dgges_obj->ldvsl,
									dgges_obj->vsr,
									dgges_obj->ldvsr
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    
    dgges_obj->diff_a =  computeDiff_d( (dgges_obj->lda)*(dgges_obj->n), 
                dgges_obj->a, dgges_obj->aref );

    dgges_obj->diff_b =  computeDiff_d( (dgges_obj->ldb)*(dgges_obj->n), 
                dgges_obj->b, dgges_obj->bref );

    dgges_obj->diff_vsl =  computeDiff_d( (dgges_obj->ldvsl)*(dgges_obj->n), 
                dgges_obj->vsl, dgges_obj->vslref );

    dgges_obj->diff_vsr =  computeDiff_d( (dgges_obj->ldvsr)*(dgges_obj->n), 
                dgges_obj->vsr, dgges_obj->vsrref );

    dgges_obj->diff_alphar =  computeDiff_d( dgges_obj->n, 
                dgges_obj->alphar, dgges_obj->alpharref );

    dgges_obj->diff_alphai =  computeDiff_d( dgges_obj->n, 
                dgges_obj->alphai, dgges_obj->alphairef );

    dgges_obj->diff_beta =  computeDiff_d( dgges_obj->n, 
                dgges_obj->beta, dgges_obj->betaref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n gges double: \n diff_a: %f \n diff_b: %f \n \
diff_vsl: %f \n diff_vsr: %f \n diff_alphar: %f \n\
diff_alphai: %f \n diff_beta: %f \n ",
       dgges_obj->diff_a, dgges_obj->diff_b, dgges_obj->diff_vsl,
	   dgges_obj->diff_vsr, dgges_obj->diff_alphar,
       dgges_obj->diff_alphai, dgges_obj->diff_beta );
#endif
}

TEST_F(dgges_test, dgges1) {
    EXPECT_NEAR(0.0, dgges_obj->diff_vsl, dgges_obj->threshold);
    EXPECT_NEAR(0.0, dgges_obj->diff_vsr, dgges_obj->threshold);
    EXPECT_NEAR(0.0, dgges_obj->diff_alphar, dgges_obj->threshold);
    EXPECT_NEAR(0.0, dgges_obj->diff_alphai, dgges_obj->threshold);
    EXPECT_NEAR(0.0, dgges_obj->diff_beta, dgges_obj->threshold);
    EXPECT_EQ(dgges_obj->sdim, dgges_obj->sdimref);
}

TEST_F(dgges_test, dgges2) {
    EXPECT_NEAR(0.0, dgges_obj->diff_vsl, dgges_obj->threshold);
    EXPECT_NEAR(0.0, dgges_obj->diff_vsr, dgges_obj->threshold);
    EXPECT_NEAR(0.0, dgges_obj->diff_alphar, dgges_obj->threshold);
    EXPECT_NEAR(0.0, dgges_obj->diff_alphai, dgges_obj->threshold);
    EXPECT_NEAR(0.0, dgges_obj->diff_beta, dgges_obj->threshold);
    EXPECT_EQ(dgges_obj->sdim, dgges_obj->sdimref);
}
TEST_F(dgges_test, dgges3) {
    EXPECT_NEAR(0.0, dgges_obj->diff_vsl, dgges_obj->threshold);
    EXPECT_NEAR(0.0, dgges_obj->diff_vsr, dgges_obj->threshold);
    EXPECT_NEAR(0.0, dgges_obj->diff_alphar, dgges_obj->threshold);
    EXPECT_NEAR(0.0, dgges_obj->diff_alphai, dgges_obj->threshold);
    EXPECT_NEAR(0.0, dgges_obj->diff_beta, dgges_obj->threshold);
    EXPECT_EQ(dgges_obj->sdim, dgges_obj->sdimref);
}
TEST_F(dgges_test, dgges4) {
    EXPECT_NEAR(0.0, dgges_obj->diff_vsl, dgges_obj->threshold);
    EXPECT_NEAR(0.0, dgges_obj->diff_vsr, dgges_obj->threshold);
    EXPECT_NEAR(0.0, dgges_obj->diff_alphar, dgges_obj->threshold);
    EXPECT_NEAR(0.0, dgges_obj->diff_alphai, dgges_obj->threshold);
    EXPECT_NEAR(0.0, dgges_obj->diff_beta, dgges_obj->threshold);
    EXPECT_EQ(dgges_obj->sdim, dgges_obj->sdimref);
}


/* Begin gges_scomplex_common_parameters  class definition */
class gges_scomplex_parameters{

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

      gges_scomplex_parameters (int matrix_layout_i, char jobvsl_i,
						char jobvsr_i, char sort_i, lapack_int n_i );
      ~gges_scomplex_parameters ();
};

/* Constructor definition  lapack_complex_float_common_parameters */
gges_scomplex_parameters:: gges_scomplex_parameters (int matrix_layout_i,
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
   printf(" \n gges lapack_complex_float: matrix_layout: %d n: %d  jobvsl: %c \
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
       EXPECT_FALSE( true) << "gges_scomplex_parameters object: malloc error. Exiting ";
       gges_free();
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
    

/* Destructor definition  'gges_scomplex_common_parameters' */
gges_scomplex_parameters :: ~gges_scomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   gges_free();
} 

//  Test fixture class definition
class cgges_test  : public  ::testing::Test {
public:
   gges_scomplex_parameters  *cgges_obj;
   void SetUp();  
   void TearDown () { delete cgges_obj; }
};

void cgges_test::SetUp()
{

    /* LAPACKE_cgges prototype */	
    typedef int (*Fptr_NL_LAPACKE_cgges) (int matrix_layout, char jobvsl, char jobvsr, char sort, LAPACK_C_SELECT2 select, lapack_int n, lapack_complex_float* a, lapack_int lda, lapack_complex_float* b, lapack_int ldb, lapack_int* sdim, lapack_complex_float* alpha, lapack_complex_float* beta, lapack_complex_float* vsl, lapack_int ldvsl, lapack_complex_float* vsr, lapack_int ldvsr);
				 
    Fptr_NL_LAPACKE_cgges CGGES;

    cgges_obj = new  gges_scomplex_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].jobvsl,
										 eig_non_sym_paramslist[idx].jobvsr,
										 eig_non_sym_paramslist[idx].sort,
                                         eig_paramslist[idx].n );
                                         
    cgges_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    cgges_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgges_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgges_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgges_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CGGES = (Fptr_NL_LAPACKE_cgges)dlsym(cgges_obj->hModule, "LAPACKE_cgges");
    ASSERT_TRUE(CGGES != NULL) << "failed to ppt the Netlib LAPACKE_cgges symbol";
    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    cgges_obj->inforef = CGGES(   cgges_obj->matrix_layout,
                                    cgges_obj->jobvsl,
                                    cgges_obj->jobvsr,
                                    cgges_obj->sort,
									(LAPACK_C_SELECT2)cgges_obj->selectref,
                                    cgges_obj->n,
                                    cgges_obj->aref,
                                    cgges_obj->lda,
                                    cgges_obj->bref,
                                    cgges_obj->ldb,
                                    &cgges_obj->sdimref,
                                    cgges_obj->alpharref,
                                    cgges_obj->betaref,
									cgges_obj->vslref,
									cgges_obj->ldvsl,
									cgges_obj->vsrref,
									cgges_obj->ldvsr
                                    );

    /* Compute libflame's Lapacke o/p  */
    cgges_obj->inforef =  LAPACKE_cgges (   cgges_obj->matrix_layout,
                                    cgges_obj->jobvsl,
                                    cgges_obj->jobvsr,
                                    cgges_obj->sort,
									(LAPACK_C_SELECT2)cgges_obj->select,
                                    cgges_obj->n,
                                    cgges_obj->a,
                                    cgges_obj->lda,
                                    cgges_obj->b,
                                    cgges_obj->ldb,
                                    &cgges_obj->sdim,
                                    cgges_obj->alphar,
                                    cgges_obj->beta,
									cgges_obj->vsl,
									cgges_obj->ldvsl,
									cgges_obj->vsr,
									cgges_obj->ldvsr
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    
    cgges_obj->diff_a =  computeDiff_c( (cgges_obj->lda)*(cgges_obj->n), 
                cgges_obj->a, cgges_obj->aref );

    cgges_obj->diff_b =  computeDiff_c( (cgges_obj->ldb)*(cgges_obj->n), 
                cgges_obj->b, cgges_obj->bref );

    cgges_obj->diff_vsl =  computeDiff_c( (cgges_obj->ldvsl)*(cgges_obj->n), 
                cgges_obj->vsl, cgges_obj->vslref );

    cgges_obj->diff_vsr =  computeDiff_c( (cgges_obj->ldvsr)*(cgges_obj->n), 
                cgges_obj->vsr, cgges_obj->vsrref );

    cgges_obj->diff_alphar =  computeDiff_c( cgges_obj->n, 
                cgges_obj->alphar, cgges_obj->alpharref );

    cgges_obj->diff_beta =  computeDiff_c( cgges_obj->n, 
                cgges_obj->beta, cgges_obj->betaref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n gges lapack_complex_float: \n diff_a: %f \n diff_b: %f \n \
diff_vsl: %f \n diff_vsr: %f \n diff_alphar: %f \n\
diff_beta: %f \n ",
       cgges_obj->diff_a, cgges_obj->diff_b, cgges_obj->diff_vsl,
	   cgges_obj->diff_vsr, cgges_obj->diff_alphar,
       cgges_obj->diff_beta );
#endif
}

TEST_F(cgges_test, cgges1) {
    EXPECT_NEAR(0.0, cgges_obj->diff_vsl, cgges_obj->threshold);
    EXPECT_NEAR(0.0, cgges_obj->diff_vsr, cgges_obj->threshold);
    EXPECT_NEAR(0.0, cgges_obj->diff_alphar, cgges_obj->threshold);
    EXPECT_NEAR(0.0, cgges_obj->diff_beta, cgges_obj->threshold);
    EXPECT_EQ(cgges_obj->sdim, cgges_obj->sdimref);
}

TEST_F(cgges_test, cgges2) {
    EXPECT_NEAR(0.0, cgges_obj->diff_vsl, cgges_obj->threshold);
    EXPECT_NEAR(0.0, cgges_obj->diff_vsr, cgges_obj->threshold);
    EXPECT_NEAR(0.0, cgges_obj->diff_alphar, cgges_obj->threshold);
    EXPECT_NEAR(0.0, cgges_obj->diff_beta, cgges_obj->threshold);
    EXPECT_EQ(cgges_obj->sdim, cgges_obj->sdimref);
}
TEST_F(cgges_test, cgges3) {
    EXPECT_NEAR(0.0, cgges_obj->diff_vsl, cgges_obj->threshold);
    EXPECT_NEAR(0.0, cgges_obj->diff_vsr, cgges_obj->threshold);
    EXPECT_NEAR(0.0, cgges_obj->diff_alphar, cgges_obj->threshold);
    EXPECT_NEAR(0.0, cgges_obj->diff_beta, cgges_obj->threshold);
    EXPECT_EQ(cgges_obj->sdim, cgges_obj->sdimref);
}
TEST_F(cgges_test, cgges4) {
    EXPECT_NEAR(0.0, cgges_obj->diff_vsl, cgges_obj->threshold);
    EXPECT_NEAR(0.0, cgges_obj->diff_vsr, cgges_obj->threshold);
    EXPECT_NEAR(0.0, cgges_obj->diff_alphar, cgges_obj->threshold);
    EXPECT_NEAR(0.0, cgges_obj->diff_beta, cgges_obj->threshold);
    EXPECT_EQ(cgges_obj->sdim, cgges_obj->sdimref);
}

/* Begin gges_dcomplex_common_parameters  class definition */
class gges_dcomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    double diff_a, diff_b;
    double diff_alphai, diff_alphar, diff_beta, diff_vsl, diff_vsr;

    double threshold;
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

      gges_dcomplex_parameters (int matrix_layout_i, char jobvsl_i,
						char jobvsr_i, char sort_i, lapack_int n_i );
      ~gges_dcomplex_parameters ();
};

/* Constructor definition  lapack_complex_double_common_parameters */
gges_dcomplex_parameters:: gges_dcomplex_parameters (int matrix_layout_i,
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
   printf(" \n gges lapack_complex_double: matrix_layout: %d n: %d  jobvsl: %c \
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
       EXPECT_FALSE( true) << "gges_dcomplex_parameters object: malloc error. Exiting ";
       gges_free();
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
    

/* Destructor definition  'gges_dcomplex_common_parameters' */
gges_dcomplex_parameters :: ~gges_dcomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   gges_free();
} 

//  Test fixture class definition
class zgges_test  : public  ::testing::Test {
public:
   gges_dcomplex_parameters  *zgges_obj;
   void SetUp();  
   void TearDown () { delete zgges_obj; }
};

void zgges_test::SetUp()
{

    /* LAPACKE_zgges prototype */	
    typedef int (*Fptr_NL_LAPACKE_zgges) (int matrix_layout, char jobvsl, char jobvsr, char sort, LAPACK_Z_SELECT2 select, lapack_int n, lapack_complex_double* a, lapack_int lda, lapack_complex_double* b, lapack_int ldb, lapack_int* sdim, lapack_complex_double* alpha, lapack_complex_double* beta, lapack_complex_double* vsl, lapack_int ldvsl, lapack_complex_double* vsr, lapack_int ldvsr);
				 
    Fptr_NL_LAPACKE_zgges ZGGES;

    zgges_obj = new  gges_dcomplex_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].jobvsl,
										 eig_non_sym_paramslist[idx].jobvsr,
										 eig_non_sym_paramslist[idx].sort,
                                         eig_paramslist[idx].n );
                                         
    zgges_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    zgges_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgges_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgges_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgges_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZGGES = (Fptr_NL_LAPACKE_zgges)dlsym(zgges_obj->hModule, "LAPACKE_zgges");
    ASSERT_TRUE(ZGGES != NULL) << "failed to ppt the Netlib LAPACKE_zgges symbol";
    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    zgges_obj->inforef = ZGGES(   zgges_obj->matrix_layout,
                                    zgges_obj->jobvsl,
                                    zgges_obj->jobvsr,
                                    zgges_obj->sort,
									(LAPACK_Z_SELECT2)zgges_obj->selectref,
                                    zgges_obj->n,
                                    zgges_obj->aref,
                                    zgges_obj->lda,
                                    zgges_obj->bref,
                                    zgges_obj->ldb,
                                    &zgges_obj->sdimref,
                                    zgges_obj->alpharref,
                                    zgges_obj->betaref,
									zgges_obj->vslref,
									zgges_obj->ldvsl,
									zgges_obj->vsrref,
									zgges_obj->ldvsr
                                    );

    /* Compute libflame's Lapacke o/p  */
    zgges_obj->inforef =  LAPACKE_zgges (   zgges_obj->matrix_layout,
                                    zgges_obj->jobvsl,
                                    zgges_obj->jobvsr,
                                    zgges_obj->sort,
									(LAPACK_Z_SELECT2)zgges_obj->select,
                                    zgges_obj->n,
                                    zgges_obj->a,
                                    zgges_obj->lda,
                                    zgges_obj->b,
                                    zgges_obj->ldb,
                                    &zgges_obj->sdim,
                                    zgges_obj->alphar,
                                    zgges_obj->beta,
									zgges_obj->vsl,
									zgges_obj->ldvsl,
									zgges_obj->vsr,
									zgges_obj->ldvsr
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    
    zgges_obj->diff_a =  computeDiff_z( (zgges_obj->lda)*(zgges_obj->n), 
                zgges_obj->a, zgges_obj->aref );

    zgges_obj->diff_b =  computeDiff_z( (zgges_obj->ldb)*(zgges_obj->n), 
                zgges_obj->b, zgges_obj->bref );

    zgges_obj->diff_vsl =  computeDiff_z( (zgges_obj->ldvsl)*(zgges_obj->n), 
                zgges_obj->vsl, zgges_obj->vslref );

    zgges_obj->diff_vsr =  computeDiff_z( (zgges_obj->ldvsr)*(zgges_obj->n), 
                zgges_obj->vsr, zgges_obj->vsrref );

    zgges_obj->diff_alphar =  computeDiff_z( zgges_obj->n, 
                zgges_obj->alphar, zgges_obj->alpharref );

    zgges_obj->diff_beta =  computeDiff_z( zgges_obj->n, 
                zgges_obj->beta, zgges_obj->betaref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n gges lapack_complex_double: \n diff_a: %f \n diff_b: %f \n \
diff_vsl: %f \n diff_vsr: %f \n diff_alphar: %f \n\
diff_beta: %f \n ",
       zgges_obj->diff_a, zgges_obj->diff_b, zgges_obj->diff_vsl,
	   zgges_obj->diff_vsr, zgges_obj->diff_alphar,
       zgges_obj->diff_beta );
#endif
}

TEST_F(zgges_test, zgges1) {
    EXPECT_NEAR(0.0, zgges_obj->diff_vsl, zgges_obj->threshold);
    EXPECT_NEAR(0.0, zgges_obj->diff_vsr, zgges_obj->threshold);
    EXPECT_NEAR(0.0, zgges_obj->diff_alphar, zgges_obj->threshold);
    EXPECT_NEAR(0.0, zgges_obj->diff_beta, zgges_obj->threshold);
    EXPECT_EQ(zgges_obj->sdim, zgges_obj->sdimref);
}

TEST_F(zgges_test, zgges2) {
    EXPECT_NEAR(0.0, zgges_obj->diff_vsl, zgges_obj->threshold);
    EXPECT_NEAR(0.0, zgges_obj->diff_vsr, zgges_obj->threshold);
    EXPECT_NEAR(0.0, zgges_obj->diff_alphar, zgges_obj->threshold);
    EXPECT_NEAR(0.0, zgges_obj->diff_beta, zgges_obj->threshold);
    EXPECT_EQ(zgges_obj->sdim, zgges_obj->sdimref);
}
TEST_F(zgges_test, zgges3) {
    EXPECT_NEAR(0.0, zgges_obj->diff_vsl, zgges_obj->threshold);
    EXPECT_NEAR(0.0, zgges_obj->diff_vsr, zgges_obj->threshold);
    EXPECT_NEAR(0.0, zgges_obj->diff_alphar, zgges_obj->threshold);
    EXPECT_NEAR(0.0, zgges_obj->diff_beta, zgges_obj->threshold);
    EXPECT_EQ(zgges_obj->sdim, zgges_obj->sdimref);
}
TEST_F(zgges_test, zgges4) {
    EXPECT_NEAR(0.0, zgges_obj->diff_vsl, zgges_obj->threshold);
    EXPECT_NEAR(0.0, zgges_obj->diff_vsr, zgges_obj->threshold);
    EXPECT_NEAR(0.0, zgges_obj->diff_alphar, zgges_obj->threshold);
    EXPECT_NEAR(0.0, zgges_obj->diff_beta, zgges_obj->threshold);
    EXPECT_EQ(zgges_obj->sdim, zgges_obj->sdimref);
}
