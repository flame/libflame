#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"

#define gges3_free() \
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

/* Begin gges3_float_common_parameters  class definition */
class gges3_float_parameters{

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

      gges3_float_parameters (int matrix_layout_i, char jobvsl_i,
						char jobvsr_i, char sort_i, lapack_int n_i );
      ~gges3_float_parameters ();
};

/* Constructor definition  float_common_parameters */
gges3_float_parameters:: gges3_float_parameters (int matrix_layout_i,
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
   printf(" \n gges3 float: matrix_layout: %d n: %d  jobvsl: %c \
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
       EXPECT_FALSE( true) << "gges3_float_parameters object: malloc error. Exiting ";
       gges3_free();
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
    

/* Destructor definition  'gges3_float_common_parameters' */
gges3_float_parameters :: ~gges3_float_parameters ()
{
   /* De-Allocate memory for the input matrices */
   gges3_free();
} 

//  Test fixture class definition
class sgges3_test  : public  ::testing::Test {
public:
   gges3_float_parameters  *sgges3_obj;
   void SetUp();  
   void TearDown () { delete sgges3_obj; }
};

void sgges3_test::SetUp()
{

    /* LAPACKE_sgges3 prototype */	
    typedef int (*Fptr_NL_LAPACKE_sgges3) (int matrix_layout, char jobvsl,
		char jobvsr, char sort, LAPACK_S_SELECT3 select, lapack_int n, 
		float* a, lapack_int lda, float* b, lapack_int ldb, lapack_int* sdim,
		float* alphar, float* alphai, float* beta, float* vsl,
		lapack_int ldvsl, float* vsr, lapack_int ldvsr);
				 
    Fptr_NL_LAPACKE_sgges3 SGGES;

    sgges3_obj = new  gges3_float_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].jobvsl,
										 eig_non_sym_paramslist[idx].jobvsr,
										 eig_non_sym_paramslist[idx].sort,
                                         eig_paramslist[idx].n );
                                         
    sgges3_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    sgges3_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgges3_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgges3_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgges3_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SGGES = (Fptr_NL_LAPACKE_sgges3)dlsym(sgges3_obj->hModule, "LAPACKE_sgges3");
    ASSERT_TRUE(SGGES != NULL) << "failed to ppt the Netlib LAPACKE_sgges3 symbol";

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    sgges3_obj->inforef = SGGES(   sgges3_obj->matrix_layout,
                                    sgges3_obj->jobvsl,
                                    sgges3_obj->jobvsr,
                                    sgges3_obj->sort,
									(LAPACK_S_SELECT3)sgges3_obj->selectref,
                                    sgges3_obj->n,
                                    sgges3_obj->aref,
                                    sgges3_obj->lda,
                                    sgges3_obj->bref,
                                    sgges3_obj->ldb,
                                    &sgges3_obj->sdimref,
                                    sgges3_obj->alpharref,
                                    sgges3_obj->alphairef,
                                    sgges3_obj->betaref,
									sgges3_obj->vslref,
									sgges3_obj->ldvsl,
									sgges3_obj->vsrref,
									sgges3_obj->ldvsr
                                    );

    /* Compute libflame's Lapacke o/p  */
    sgges3_obj->inforef =  LAPACKE_sgges3 (   sgges3_obj->matrix_layout,
                                    sgges3_obj->jobvsl,
                                    sgges3_obj->jobvsr,
                                    sgges3_obj->sort,
									(LAPACK_S_SELECT3)sgges3_obj->select,
                                    sgges3_obj->n,
                                    sgges3_obj->a,
                                    sgges3_obj->lda,
                                    sgges3_obj->b,
                                    sgges3_obj->ldb,
                                    &sgges3_obj->sdim,
                                    sgges3_obj->alphar,
                                    sgges3_obj->alphai,
                                    sgges3_obj->beta,
									sgges3_obj->vsl,
									sgges3_obj->ldvsl,
									sgges3_obj->vsr,
									sgges3_obj->ldvsr
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    
    sgges3_obj->diff_a =  computeDiff_s( (sgges3_obj->lda)*(sgges3_obj->n), 
                sgges3_obj->a, sgges3_obj->aref );

    sgges3_obj->diff_b =  computeDiff_s( (sgges3_obj->ldb)*(sgges3_obj->n), 
                sgges3_obj->b, sgges3_obj->bref );

    sgges3_obj->diff_vsl =  computeDiff_s( (sgges3_obj->ldvsl)*(sgges3_obj->n), 
                sgges3_obj->vsl, sgges3_obj->vslref );

    sgges3_obj->diff_vsr =  computeDiff_s( (sgges3_obj->ldvsr)*(sgges3_obj->n), 
                sgges3_obj->vsr, sgges3_obj->vsrref );

    sgges3_obj->diff_alphar =  computeDiff_s( sgges3_obj->n, 
                sgges3_obj->alphar, sgges3_obj->alpharref );

    sgges3_obj->diff_alphai =  computeDiff_s( sgges3_obj->n, 
                sgges3_obj->alphai, sgges3_obj->alphairef );

    sgges3_obj->diff_beta =  computeDiff_s( sgges3_obj->n, 
                sgges3_obj->beta, sgges3_obj->betaref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n gges3 float: \n diff_a: %f \n diff_b: %f \n \
diff_vsl: %f \n diff_vsr: %f \n diff_alphar: %f \n\
diff_alphai: %f \n diff_beta: %f \n ",
       sgges3_obj->diff_a, sgges3_obj->diff_b, sgges3_obj->diff_vsl,
	   sgges3_obj->diff_vsr, sgges3_obj->diff_alphar,
       sgges3_obj->diff_alphai, sgges3_obj->diff_beta );
#endif
}

TEST_F(sgges3_test, sgges31) {
    //EXPECT_NEAR(0.0, sgges3_obj->diff_a, sgges3_obj->threshold);
    //EXPECT_NEAR(0.0, sgges3_obj->diff_b, sgges3_obj->threshold);
    EXPECT_NEAR(0.0, sgges3_obj->diff_vsl, sgges3_obj->threshold);
    EXPECT_NEAR(0.0, sgges3_obj->diff_vsr, sgges3_obj->threshold);
    EXPECT_NEAR(0.0, sgges3_obj->diff_alphar, sgges3_obj->threshold);
    EXPECT_NEAR(0.0, sgges3_obj->diff_alphai, sgges3_obj->threshold);
    EXPECT_NEAR(0.0, sgges3_obj->diff_beta, sgges3_obj->threshold);
    EXPECT_EQ(sgges3_obj->sdim, sgges3_obj->sdimref);
}

TEST_F(sgges3_test, sgges32) {
    //EXPECT_NEAR(0.0, sgges3_obj->diff_a, sgges3_obj->threshold);
    //EXPECT_NEAR(0.0, sgges3_obj->diff_b, sgges3_obj->threshold);
    EXPECT_NEAR(0.0, sgges3_obj->diff_vsl, sgges3_obj->threshold);
    EXPECT_NEAR(0.0, sgges3_obj->diff_vsr, sgges3_obj->threshold);
    EXPECT_NEAR(0.0, sgges3_obj->diff_alphar, sgges3_obj->threshold);
    EXPECT_NEAR(0.0, sgges3_obj->diff_alphai, sgges3_obj->threshold);
    EXPECT_NEAR(0.0, sgges3_obj->diff_beta, sgges3_obj->threshold);
    EXPECT_EQ(sgges3_obj->sdim, sgges3_obj->sdimref);
}
TEST_F(sgges3_test, sgges33) {
    //EXPECT_NEAR(0.0, sgges3_obj->diff_a, sgges3_obj->threshold);
    //EXPECT_NEAR(0.0, sgges3_obj->diff_b, sgges3_obj->threshold);
    EXPECT_NEAR(0.0, sgges3_obj->diff_vsl, sgges3_obj->threshold);
    EXPECT_NEAR(0.0, sgges3_obj->diff_vsr, sgges3_obj->threshold);
    EXPECT_NEAR(0.0, sgges3_obj->diff_alphar, sgges3_obj->threshold);
    EXPECT_NEAR(0.0, sgges3_obj->diff_alphai, sgges3_obj->threshold);
    EXPECT_NEAR(0.0, sgges3_obj->diff_beta, sgges3_obj->threshold);
    EXPECT_EQ(sgges3_obj->sdim, sgges3_obj->sdimref);
}
TEST_F(sgges3_test, sgges34) {
    //EXPECT_NEAR(0.0, sgges3_obj->diff_a, sgges3_obj->threshold);
    //EXPECT_NEAR(0.0, sgges3_obj->diff_b, sgges3_obj->threshold);
    EXPECT_NEAR(0.0, sgges3_obj->diff_vsl, sgges3_obj->threshold);
    EXPECT_NEAR(0.0, sgges3_obj->diff_vsr, sgges3_obj->threshold);
    EXPECT_NEAR(0.0, sgges3_obj->diff_alphar, sgges3_obj->threshold);
    EXPECT_NEAR(0.0, sgges3_obj->diff_alphai, sgges3_obj->threshold);
    EXPECT_NEAR(0.0, sgges3_obj->diff_beta, sgges3_obj->threshold);
    EXPECT_EQ(sgges3_obj->sdim, sgges3_obj->sdimref);
}


/* Begin gges3_double_common_parameters  class definition */
class gges3_double_parameters{

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

      gges3_double_parameters (int matrix_layout_i, char jobvsl_i,
						char jobvsr_i, char sort_i, lapack_int n_i );
      ~gges3_double_parameters ();
};

/* Constructor definition  double_common_parameters */
gges3_double_parameters:: gges3_double_parameters (int matrix_layout_i,
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
   printf(" \n gges3 double: matrix_layout: %d n: %d  jobvsl: %c \
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
       EXPECT_FALSE( true) << "gges3_double_parameters object: malloc error. Exiting ";
       gges3_free();
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
    

/* Destructor definition  'gges3_double_common_parameters' */
gges3_double_parameters :: ~gges3_double_parameters ()
{
   /* De-Allocate memory for the input matrices */
   gges3_free();
} 

//  Test fixture class definition
class dgges3_test  : public  ::testing::Test {
public:
   gges3_double_parameters  *dgges3_obj;
   void SetUp();  
   void TearDown () { delete dgges3_obj; }
};

void dgges3_test::SetUp()
{

    /* LAPACKE_dgges3 prototype */	
    typedef int (*Fptr_NL_LAPACKE_dgges3) (int matrix_layout, char jobvsl,
		char jobvsr, char sort, LAPACK_D_SELECT3 select, lapack_int n, 
		double* a, lapack_int lda, double* b, lapack_int ldb, lapack_int* sdim,
		double* alphar, double* alphai, double* beta, double* vsl,
		lapack_int ldvsl, double* vsr, lapack_int ldvsr);
				 
    Fptr_NL_LAPACKE_dgges3 DGGES;

    dgges3_obj = new  gges3_double_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].jobvsl,
										 eig_non_sym_paramslist[idx].jobvsr,
										 eig_non_sym_paramslist[idx].sort,
                                         eig_paramslist[idx].n );
                                         
    dgges3_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    dgges3_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgges3_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgges3_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgges3_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DGGES = (Fptr_NL_LAPACKE_dgges3)dlsym(dgges3_obj->hModule, "LAPACKE_dgges3");
    ASSERT_TRUE(DGGES != NULL) << "failed to ppt the Netlib LAPACKE_dgges3 symbol";

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    dgges3_obj->inforef = DGGES(   dgges3_obj->matrix_layout,
                                    dgges3_obj->jobvsl,
                                    dgges3_obj->jobvsr,
                                    dgges3_obj->sort,
									(LAPACK_D_SELECT3)dgges3_obj->selectref,
                                    dgges3_obj->n,
                                    dgges3_obj->aref,
                                    dgges3_obj->lda,
                                    dgges3_obj->bref,
                                    dgges3_obj->ldb,
                                    &dgges3_obj->sdimref,
                                    dgges3_obj->alpharref,
                                    dgges3_obj->alphairef,
                                    dgges3_obj->betaref,
									dgges3_obj->vslref,
									dgges3_obj->ldvsl,
									dgges3_obj->vsrref,
									dgges3_obj->ldvsr
                                    );

    /* Compute libflame's Lapacke o/p  */
    dgges3_obj->inforef =  LAPACKE_dgges3 (   dgges3_obj->matrix_layout,
                                    dgges3_obj->jobvsl,
                                    dgges3_obj->jobvsr,
                                    dgges3_obj->sort,
									(LAPACK_D_SELECT3)dgges3_obj->select,
                                    dgges3_obj->n,
                                    dgges3_obj->a,
                                    dgges3_obj->lda,
                                    dgges3_obj->b,
                                    dgges3_obj->ldb,
                                    &dgges3_obj->sdim,
                                    dgges3_obj->alphar,
                                    dgges3_obj->alphai,
                                    dgges3_obj->beta,
									dgges3_obj->vsl,
									dgges3_obj->ldvsl,
									dgges3_obj->vsr,
									dgges3_obj->ldvsr
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    
    dgges3_obj->diff_a =  computeDiff_d( (dgges3_obj->lda)*(dgges3_obj->n), 
                dgges3_obj->a, dgges3_obj->aref );

    dgges3_obj->diff_b =  computeDiff_d( (dgges3_obj->ldb)*(dgges3_obj->n), 
                dgges3_obj->b, dgges3_obj->bref );

    dgges3_obj->diff_vsl =  computeDiff_d( (dgges3_obj->ldvsl)*(dgges3_obj->n), 
                dgges3_obj->vsl, dgges3_obj->vslref );

    dgges3_obj->diff_vsr =  computeDiff_d( (dgges3_obj->ldvsr)*(dgges3_obj->n), 
                dgges3_obj->vsr, dgges3_obj->vsrref );

    dgges3_obj->diff_alphar =  computeDiff_d( dgges3_obj->n, 
                dgges3_obj->alphar, dgges3_obj->alpharref );

    dgges3_obj->diff_alphai =  computeDiff_d( dgges3_obj->n, 
                dgges3_obj->alphai, dgges3_obj->alphairef );

    dgges3_obj->diff_beta =  computeDiff_d( dgges3_obj->n, 
                dgges3_obj->beta, dgges3_obj->betaref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n gges3 double: \n diff_a: %f \n diff_b: %f \n \
diff_vsl: %f \n diff_vsr: %f \n diff_alphar: %f \n\
diff_alphai: %f \n diff_beta: %f \n ",
       dgges3_obj->diff_a, dgges3_obj->diff_b, dgges3_obj->diff_vsl,
	   dgges3_obj->diff_vsr, dgges3_obj->diff_alphar,
       dgges3_obj->diff_alphai, dgges3_obj->diff_beta );
#endif
}

TEST_F(dgges3_test, dgges31) {
    //EXPECT_NEAR(0.0, dgges3_obj->diff_a, dgges3_obj->threshold);
    //EXPECT_NEAR(0.0, dgges3_obj->diff_b, dgges3_obj->threshold);
    EXPECT_NEAR(0.0, dgges3_obj->diff_vsl, dgges3_obj->threshold);
    EXPECT_NEAR(0.0, dgges3_obj->diff_vsr, dgges3_obj->threshold);
    EXPECT_NEAR(0.0, dgges3_obj->diff_alphar, dgges3_obj->threshold);
    EXPECT_NEAR(0.0, dgges3_obj->diff_alphai, dgges3_obj->threshold);
    EXPECT_NEAR(0.0, dgges3_obj->diff_beta, dgges3_obj->threshold);
    EXPECT_EQ(dgges3_obj->sdim, dgges3_obj->sdimref);
}

TEST_F(dgges3_test, dgges32) {
    //EXPECT_NEAR(0.0, dgges3_obj->diff_a, dgges3_obj->threshold);
    //EXPECT_NEAR(0.0, dgges3_obj->diff_b, dgges3_obj->threshold);
    EXPECT_NEAR(0.0, dgges3_obj->diff_vsl, dgges3_obj->threshold);
    EXPECT_NEAR(0.0, dgges3_obj->diff_vsr, dgges3_obj->threshold);
    EXPECT_NEAR(0.0, dgges3_obj->diff_alphar, dgges3_obj->threshold);
    EXPECT_NEAR(0.0, dgges3_obj->diff_alphai, dgges3_obj->threshold);
    EXPECT_NEAR(0.0, dgges3_obj->diff_beta, dgges3_obj->threshold);
    EXPECT_EQ(dgges3_obj->sdim, dgges3_obj->sdimref);
}
TEST_F(dgges3_test, dgges33) {
    //EXPECT_NEAR(0.0, dgges3_obj->diff_a, dgges3_obj->threshold);
    //EXPECT_NEAR(0.0, dgges3_obj->diff_b, dgges3_obj->threshold);
    EXPECT_NEAR(0.0, dgges3_obj->diff_vsl, dgges3_obj->threshold);
    EXPECT_NEAR(0.0, dgges3_obj->diff_vsr, dgges3_obj->threshold);
    EXPECT_NEAR(0.0, dgges3_obj->diff_alphar, dgges3_obj->threshold);
    EXPECT_NEAR(0.0, dgges3_obj->diff_alphai, dgges3_obj->threshold);
    EXPECT_NEAR(0.0, dgges3_obj->diff_beta, dgges3_obj->threshold);
    EXPECT_EQ(dgges3_obj->sdim, dgges3_obj->sdimref);
}
TEST_F(dgges3_test, dgges34) {
    //EXPECT_NEAR(0.0, dgges3_obj->diff_a, dgges3_obj->threshold);
    //EXPECT_NEAR(0.0, dgges3_obj->diff_b, dgges3_obj->threshold);
    EXPECT_NEAR(0.0, dgges3_obj->diff_vsl, dgges3_obj->threshold);
    EXPECT_NEAR(0.0, dgges3_obj->diff_vsr, dgges3_obj->threshold);
    EXPECT_NEAR(0.0, dgges3_obj->diff_alphar, dgges3_obj->threshold);
    EXPECT_NEAR(0.0, dgges3_obj->diff_alphai, dgges3_obj->threshold);
    EXPECT_NEAR(0.0, dgges3_obj->diff_beta, dgges3_obj->threshold);
    EXPECT_EQ(dgges3_obj->sdim, dgges3_obj->sdimref);
}


/* Begin gges3_scomplex_common_parameters  class definition */
class gges3_scomplex_parameters{

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

      gges3_scomplex_parameters (int matrix_layout_i, char jobvsl_i,
						char jobvsr_i, char sort_i, lapack_int n_i );
      ~gges3_scomplex_parameters ();
};

/* Constructor definition  lapack_complex_float_common_parameters */
gges3_scomplex_parameters:: gges3_scomplex_parameters (int matrix_layout_i,
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
   printf(" \n gges3 lapack_complex_float: matrix_layout: %d n: %d  jobvsl: %c \
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
       EXPECT_FALSE( true) << "gges3_scomplex_parameters object: malloc error. Exiting ";
       gges3_free();
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
    

/* Destructor definition  'gges3_scomplex_common_parameters' */
gges3_scomplex_parameters :: ~gges3_scomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   gges3_free();
} 

//  Test fixture class definition
class cgges3_test  : public  ::testing::Test {
public:
   gges3_scomplex_parameters  *cgges3_obj;
   void SetUp();  
   void TearDown () { delete cgges3_obj; }
};

void cgges3_test::SetUp()
{

    /* LAPACKE_cgges3 prototype */	
    typedef int (*Fptr_NL_LAPACKE_cgges3) (int matrix_layout, char jobvsl, char jobvsr, char sort, LAPACK_C_SELECT2 select, lapack_int n, lapack_complex_float* a, lapack_int lda, lapack_complex_float* b, lapack_int ldb, lapack_int* sdim, lapack_complex_float* alpha, lapack_complex_float* beta, lapack_complex_float* vsl, lapack_int ldvsl, lapack_complex_float* vsr, lapack_int ldvsr);
				 
    Fptr_NL_LAPACKE_cgges3 CGGES;

    cgges3_obj = new  gges3_scomplex_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].jobvsl,
										 eig_non_sym_paramslist[idx].jobvsr,
										 eig_non_sym_paramslist[idx].sort,
                                         eig_paramslist[idx].n );
                                         
    cgges3_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    cgges3_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgges3_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgges3_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgges3_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CGGES = (Fptr_NL_LAPACKE_cgges3)dlsym(cgges3_obj->hModule, "LAPACKE_cgges3");
    ASSERT_TRUE(CGGES != NULL) << "failed to ppt the Netlib LAPACKE_cgges3 symbol";
    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    cgges3_obj->inforef = CGGES(   cgges3_obj->matrix_layout,
                                    cgges3_obj->jobvsl,
                                    cgges3_obj->jobvsr,
                                    cgges3_obj->sort,
									(LAPACK_C_SELECT2)cgges3_obj->selectref,
                                    cgges3_obj->n,
                                    cgges3_obj->aref,
                                    cgges3_obj->lda,
                                    cgges3_obj->bref,
                                    cgges3_obj->ldb,
                                    &cgges3_obj->sdimref,
                                    cgges3_obj->alpharref,
                                    cgges3_obj->betaref,
									cgges3_obj->vslref,
									cgges3_obj->ldvsl,
									cgges3_obj->vsrref,
									cgges3_obj->ldvsr
                                    );

    /* Compute libflame's Lapacke o/p  */
    cgges3_obj->inforef =  LAPACKE_cgges3 (   cgges3_obj->matrix_layout,
                                    cgges3_obj->jobvsl,
                                    cgges3_obj->jobvsr,
                                    cgges3_obj->sort,
									(LAPACK_C_SELECT2)cgges3_obj->select,
                                    cgges3_obj->n,
                                    cgges3_obj->a,
                                    cgges3_obj->lda,
                                    cgges3_obj->b,
                                    cgges3_obj->ldb,
                                    &cgges3_obj->sdim,
                                    cgges3_obj->alphar,
                                    cgges3_obj->beta,
									cgges3_obj->vsl,
									cgges3_obj->ldvsl,
									cgges3_obj->vsr,
									cgges3_obj->ldvsr
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    
    cgges3_obj->diff_a =  computeDiff_c( (cgges3_obj->lda)*(cgges3_obj->n), 
                cgges3_obj->a, cgges3_obj->aref );

    cgges3_obj->diff_b =  computeDiff_c( (cgges3_obj->ldb)*(cgges3_obj->n), 
                cgges3_obj->b, cgges3_obj->bref );

    cgges3_obj->diff_vsl =  computeDiff_c( (cgges3_obj->ldvsl)*(cgges3_obj->n), 
                cgges3_obj->vsl, cgges3_obj->vslref );

    cgges3_obj->diff_vsr =  computeDiff_c( (cgges3_obj->ldvsr)*(cgges3_obj->n), 
                cgges3_obj->vsr, cgges3_obj->vsrref );

    cgges3_obj->diff_alphar =  computeDiff_c( cgges3_obj->n, 
                cgges3_obj->alphar, cgges3_obj->alpharref );

    cgges3_obj->diff_beta =  computeDiff_c( cgges3_obj->n, 
                cgges3_obj->beta, cgges3_obj->betaref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n gges3 lapack_complex_float: \n diff_a: %f \n diff_b: %f \n \
diff_vsl: %f \n diff_vsr: %f \n diff_alphar: %f \n\
diff_beta: %f \n ",
       cgges3_obj->diff_a, cgges3_obj->diff_b, cgges3_obj->diff_vsl,
	   cgges3_obj->diff_vsr, cgges3_obj->diff_alphar,
       cgges3_obj->diff_beta );
#endif
}

TEST_F(cgges3_test, cgges31) {
    //EXPECT_NEAR(0.0, cgges3_obj->diff_a, cgges3_obj->threshold);
    //EXPECT_NEAR(0.0, cgges3_obj->diff_b, cgges3_obj->threshold);
    EXPECT_NEAR(0.0, cgges3_obj->diff_vsl, cgges3_obj->threshold);
    EXPECT_NEAR(0.0, cgges3_obj->diff_vsr, cgges3_obj->threshold);
    EXPECT_NEAR(0.0, cgges3_obj->diff_alphar, cgges3_obj->threshold);
    EXPECT_NEAR(0.0, cgges3_obj->diff_beta, cgges3_obj->threshold);
    EXPECT_EQ(cgges3_obj->sdim, cgges3_obj->sdimref);
}

TEST_F(cgges3_test, cgges32) {
    //EXPECT_NEAR(0.0, cgges3_obj->diff_a, cgges3_obj->threshold);
    //EXPECT_NEAR(0.0, cgges3_obj->diff_b, cgges3_obj->threshold);
    EXPECT_NEAR(0.0, cgges3_obj->diff_vsl, cgges3_obj->threshold);
    EXPECT_NEAR(0.0, cgges3_obj->diff_vsr, cgges3_obj->threshold);
    EXPECT_NEAR(0.0, cgges3_obj->diff_alphar, cgges3_obj->threshold);
    EXPECT_NEAR(0.0, cgges3_obj->diff_beta, cgges3_obj->threshold);
    EXPECT_EQ(cgges3_obj->sdim, cgges3_obj->sdimref);
}
TEST_F(cgges3_test, cgges33) {
    //EXPECT_NEAR(0.0, cgges3_obj->diff_a, cgges3_obj->threshold);
    //EXPECT_NEAR(0.0, cgges3_obj->diff_b, cgges3_obj->threshold);
    EXPECT_NEAR(0.0, cgges3_obj->diff_vsl, cgges3_obj->threshold);
    EXPECT_NEAR(0.0, cgges3_obj->diff_vsr, cgges3_obj->threshold);
    EXPECT_NEAR(0.0, cgges3_obj->diff_alphar, cgges3_obj->threshold);
    EXPECT_NEAR(0.0, cgges3_obj->diff_beta, cgges3_obj->threshold);
    EXPECT_EQ(cgges3_obj->sdim, cgges3_obj->sdimref);
}
TEST_F(cgges3_test, cgges34) {
    //EXPECT_NEAR(0.0, cgges3_obj->diff_a, cgges3_obj->threshold);
    //EXPECT_NEAR(0.0, cgges3_obj->diff_b, cgges3_obj->threshold);
    EXPECT_NEAR(0.0, cgges3_obj->diff_vsl, cgges3_obj->threshold);
    EXPECT_NEAR(0.0, cgges3_obj->diff_vsr, cgges3_obj->threshold);
    EXPECT_NEAR(0.0, cgges3_obj->diff_alphar, cgges3_obj->threshold);
    EXPECT_NEAR(0.0, cgges3_obj->diff_beta, cgges3_obj->threshold);
    EXPECT_EQ(cgges3_obj->sdim, cgges3_obj->sdimref);
}

/* Begin gges3_dcomplex_common_parameters  class definition */
class gges3_dcomplex_parameters{

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

      gges3_dcomplex_parameters (int matrix_layout_i, char jobvsl_i,
						char jobvsr_i, char sort_i, lapack_int n_i );
      ~gges3_dcomplex_parameters ();
};

/* Constructor definition  lapack_complex_double_common_parameters */
gges3_dcomplex_parameters:: gges3_dcomplex_parameters (int matrix_layout_i,
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
   printf(" \n gges3 lapack_complex_double: matrix_layout: %d n: %d  jobvsl: %c \
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
       EXPECT_FALSE( true) << "gges3_dcomplex_parameters object: malloc error. Exiting ";
       gges3_free();
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
    

/* Destructor definition  'gges3_dcomplex_common_parameters' */
gges3_dcomplex_parameters :: ~gges3_dcomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   gges3_free();
} 

//  Test fixture class definition
class zgges3_test  : public  ::testing::Test {
public:
   gges3_dcomplex_parameters  *zgges3_obj;
   void SetUp();  
   void TearDown () { delete zgges3_obj; }
};

void zgges3_test::SetUp()
{

    /* LAPACKE_zgges3 prototype */	
    typedef int (*Fptr_NL_LAPACKE_zgges3) (int matrix_layout, char jobvsl, char jobvsr, char sort, LAPACK_Z_SELECT2 select, lapack_int n, lapack_complex_double* a, lapack_int lda, lapack_complex_double* b, lapack_int ldb, lapack_int* sdim, lapack_complex_double* alpha, lapack_complex_double* beta, lapack_complex_double* vsl, lapack_int ldvsl, lapack_complex_double* vsr, lapack_int ldvsr);
				 
    Fptr_NL_LAPACKE_zgges3 ZGGES;

    zgges3_obj = new  gges3_dcomplex_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].jobvsl,
										 eig_non_sym_paramslist[idx].jobvsr,
										 eig_non_sym_paramslist[idx].sort,
                                         eig_paramslist[idx].n );
                                         
    zgges3_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    zgges3_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgges3_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgges3_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgges3_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZGGES = (Fptr_NL_LAPACKE_zgges3)dlsym(zgges3_obj->hModule, "LAPACKE_zgges3");
    ASSERT_TRUE(ZGGES != NULL) << "failed to ppt the Netlib LAPACKE_zgges3 symbol";
    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    zgges3_obj->inforef = ZGGES(   zgges3_obj->matrix_layout,
                                    zgges3_obj->jobvsl,
                                    zgges3_obj->jobvsr,
                                    zgges3_obj->sort,
									(LAPACK_Z_SELECT2)zgges3_obj->selectref,
                                    zgges3_obj->n,
                                    zgges3_obj->aref,
                                    zgges3_obj->lda,
                                    zgges3_obj->bref,
                                    zgges3_obj->ldb,
                                    &zgges3_obj->sdimref,
                                    zgges3_obj->alpharref,
                                    zgges3_obj->betaref,
									zgges3_obj->vslref,
									zgges3_obj->ldvsl,
									zgges3_obj->vsrref,
									zgges3_obj->ldvsr
                                    );

    /* Compute libflame's Lapacke o/p  */
    zgges3_obj->inforef =  LAPACKE_zgges3 (   zgges3_obj->matrix_layout,
                                    zgges3_obj->jobvsl,
                                    zgges3_obj->jobvsr,
                                    zgges3_obj->sort,
									(LAPACK_Z_SELECT2)zgges3_obj->select,
                                    zgges3_obj->n,
                                    zgges3_obj->a,
                                    zgges3_obj->lda,
                                    zgges3_obj->b,
                                    zgges3_obj->ldb,
                                    &zgges3_obj->sdim,
                                    zgges3_obj->alphar,
                                    zgges3_obj->beta,
									zgges3_obj->vsl,
									zgges3_obj->ldvsl,
									zgges3_obj->vsr,
									zgges3_obj->ldvsr
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    
    zgges3_obj->diff_a =  computeDiff_z( (zgges3_obj->lda)*(zgges3_obj->n), 
                zgges3_obj->a, zgges3_obj->aref );

    zgges3_obj->diff_b =  computeDiff_z( (zgges3_obj->ldb)*(zgges3_obj->n), 
                zgges3_obj->b, zgges3_obj->bref );

    zgges3_obj->diff_vsl =  computeDiff_z( (zgges3_obj->ldvsl)*(zgges3_obj->n), 
                zgges3_obj->vsl, zgges3_obj->vslref );

    zgges3_obj->diff_vsr =  computeDiff_z( (zgges3_obj->ldvsr)*(zgges3_obj->n), 
                zgges3_obj->vsr, zgges3_obj->vsrref );

    zgges3_obj->diff_alphar =  computeDiff_z( zgges3_obj->n, 
                zgges3_obj->alphar, zgges3_obj->alpharref );

    zgges3_obj->diff_beta =  computeDiff_z( zgges3_obj->n, 
                zgges3_obj->beta, zgges3_obj->betaref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n gges3 lapack_complex_double: \n diff_a: %f \n diff_b: %f \n \
diff_vsl: %f \n diff_vsr: %f \n diff_alphar: %f \n\
diff_beta: %f \n ",
       zgges3_obj->diff_a, zgges3_obj->diff_b, zgges3_obj->diff_vsl,
	   zgges3_obj->diff_vsr, zgges3_obj->diff_alphar,
       zgges3_obj->diff_beta );
#endif
}

TEST_F(zgges3_test, zgges31) {
    //EXPECT_NEAR(0.0, zgges3_obj->diff_a, zgges3_obj->threshold);
    //EXPECT_NEAR(0.0, zgges3_obj->diff_b, zgges3_obj->threshold);
    EXPECT_NEAR(0.0, zgges3_obj->diff_vsl, zgges3_obj->threshold);
    EXPECT_NEAR(0.0, zgges3_obj->diff_vsr, zgges3_obj->threshold);
    EXPECT_NEAR(0.0, zgges3_obj->diff_alphar, zgges3_obj->threshold);
    EXPECT_NEAR(0.0, zgges3_obj->diff_beta, zgges3_obj->threshold);
    EXPECT_EQ(zgges3_obj->sdim, zgges3_obj->sdimref);
}

TEST_F(zgges3_test, zgges32) {
    //EXPECT_NEAR(0.0, zgges3_obj->diff_a, zgges3_obj->threshold);
    //EXPECT_NEAR(0.0, zgges3_obj->diff_b, zgges3_obj->threshold);
    EXPECT_NEAR(0.0, zgges3_obj->diff_vsl, zgges3_obj->threshold);
    EXPECT_NEAR(0.0, zgges3_obj->diff_vsr, zgges3_obj->threshold);
    EXPECT_NEAR(0.0, zgges3_obj->diff_alphar, zgges3_obj->threshold);
    EXPECT_NEAR(0.0, zgges3_obj->diff_beta, zgges3_obj->threshold);
    EXPECT_EQ(zgges3_obj->sdim, zgges3_obj->sdimref);
}
TEST_F(zgges3_test, zgges33) {
    //EXPECT_NEAR(0.0, zgges3_obj->diff_a, zgges3_obj->threshold);
    //EXPECT_NEAR(0.0, zgges3_obj->diff_b, zgges3_obj->threshold);
    EXPECT_NEAR(0.0, zgges3_obj->diff_vsl, zgges3_obj->threshold);
    EXPECT_NEAR(0.0, zgges3_obj->diff_vsr, zgges3_obj->threshold);
    EXPECT_NEAR(0.0, zgges3_obj->diff_alphar, zgges3_obj->threshold);
    EXPECT_NEAR(0.0, zgges3_obj->diff_beta, zgges3_obj->threshold);
    EXPECT_EQ(zgges3_obj->sdim, zgges3_obj->sdimref);
}
TEST_F(zgges3_test, zgges34) {
    //EXPECT_NEAR(0.0, zgges3_obj->diff_a, zgges3_obj->threshold);
    //EXPECT_NEAR(0.0, zgges3_obj->diff_b, zgges3_obj->threshold);
    EXPECT_NEAR(0.0, zgges3_obj->diff_vsl, zgges3_obj->threshold);
    EXPECT_NEAR(0.0, zgges3_obj->diff_vsr, zgges3_obj->threshold);
    EXPECT_NEAR(0.0, zgges3_obj->diff_alphar, zgges3_obj->threshold);
    EXPECT_NEAR(0.0, zgges3_obj->diff_beta, zgges3_obj->threshold);
    EXPECT_EQ(zgges3_obj->sdim, zgges3_obj->sdimref);
}
