#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define LAPACKE_TEST_VERBOSE  (1)

#define ggesx_free() \
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

/* Begin ggesx_float_common_parameters  class definition */
class ggesx_float_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_a, diff_b, diff_rconde, diff_rcondv;
    float diff_alphai, diff_alphar, diff_beta, diff_vsl, diff_vsr;

    float threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR	
    char jobvsl; // Must be 'N', or 'V'.
    char jobvsr; // Must be 'N', or 'V'.
	char sort; // Must be 'N', or 'S'.
	char sense; //  Must be 'N', 'E', 'V', or 'B'.
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
	float rconde[2], rcondv[2]; // contain the reciprocal condition numbers 
	float rconderef[2], rcondvref[2]; // contain the reciprocal condition numbers 
    /*Return Values */
    int info, inforef;

      ggesx_float_parameters (int matrix_layout_i, char jobvsl_i,
				char jobvsr_i, char sort_i, char sense_i, lapack_int n_i );
      ~ggesx_float_parameters ();
};

/* Constructor definition  float_common_parameters */
ggesx_float_parameters:: ggesx_float_parameters (int matrix_layout_i,
     char jobvsl_i, char jobvsr_i, char sort_i, char sense_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    jobvsl = jobvsl_i;
    jobvsr = jobvsr_i;
	sort = sort_i;
	sense = sense_i,
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
   printf(" \n ggesx float: matrix_layout: %d n: %d  jobvsl: %c \
jobvsr: %c \t sort: %c \t sense: %c\n", matrix_layout, n, jobvsr, 
jobvsl, sort, sense);
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
       EXPECT_FALSE( true) << "ggesx_float_parameters object: malloc error. Exiting ";
       ggesx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, lda*n );
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, ldb*n );
    lapacke_gtest_init_float_buffer_pair_with_constant( rconde, rconderef, 2, 0 );
    lapacke_gtest_init_float_buffer_pair_with_constant( rcondv, rcondvref, 2, 0 );
    lapacke_gtest_init_float_buffer_pair_with_constant( vsl, vslref, ldvsl*n, 0 );
    lapacke_gtest_init_float_buffer_pair_with_constant( vsr, vsrref, ldvsr*n, 0 );
    lapacke_gtest_init_float_buffer_pair_with_constant( alphar, alpharref, n, 0);
    lapacke_gtest_init_float_buffer_pair_with_constant( alphai, alphairef, n, 0);
    lapacke_gtest_init_float_buffer_pair_with_constant( beta, betaref, n, 0);
    lapacke_gtest_init_int_buffer_pair_with_constant(select, selectref, n, -1);

   } /* end of Constructor  */
    

/* Destructor definition  'ggesx_float_common_parameters' */
ggesx_float_parameters :: ~ggesx_float_parameters ()
{
   /* De-Allocate memory for the input matrices */
   ggesx_free();
} 

//  Test fixture class definition
class sggesx_test  : public  ::testing::Test {
public:
   ggesx_float_parameters  *sggesx_obj;
   void SetUp();  
   void TearDown () { delete sggesx_obj; }
};

void sggesx_test::SetUp()
{

    /* LAPACKE_sggesx prototype */	
    typedef int (*Fptr_NL_LAPACKE_sggesx) (int matrix_layout, char jobvsl, char jobvsr, char sort, LAPACK_S_SELECT3 select, char sense, lapack_int n, float* a, lapack_int lda, float* b, lapack_int ldb, lapack_int* sdim, float* alphar, float* alphai, float* beta, float* vsl, lapack_int ldvsl, float* vsr, lapack_int ldvsr, float* rconde, float* rcondv);
				 
    Fptr_NL_LAPACKE_sggesx SGGESX;

    sggesx_obj = new  ggesx_float_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].jobvsl,
										 eig_non_sym_paramslist[idx].jobvsr,
										 eig_non_sym_paramslist[idx].sort,
										 eig_non_sym_paramslist[idx].sense_ggesx,
                                         eig_paramslist[idx].n );
                                         
    sggesx_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    sggesx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sggesx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sggesx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sggesx_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SGGESX = (Fptr_NL_LAPACKE_sggesx)dlsym(sggesx_obj->hModule, "LAPACKE_sggesx");
    ASSERT_TRUE(SGGESX != NULL) << "failed to ppt the Netlib LAPACKE_sggesx symbol";

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    sggesx_obj->inforef = SGGESX(   sggesx_obj->matrix_layout,
                                    sggesx_obj->jobvsl,
                                    sggesx_obj->jobvsr,
                                    sggesx_obj->sort,
									(LAPACK_S_SELECT3)sggesx_obj->selectref,
                                    sggesx_obj->sense,
                                    sggesx_obj->n,
                                    sggesx_obj->aref,
                                    sggesx_obj->lda,
                                    sggesx_obj->bref,
                                    sggesx_obj->ldb,
                                    &sggesx_obj->sdimref,
                                    sggesx_obj->alpharref,
                                    sggesx_obj->alphairef,
                                    sggesx_obj->betaref,
									sggesx_obj->vslref,
									sggesx_obj->ldvsl,
									sggesx_obj->vsrref,
									sggesx_obj->ldvsr,
									sggesx_obj->rconderef,
									sggesx_obj->rcondvref
                                    );

    /* Compute libflame's Lapacke o/p  */
    sggesx_obj->inforef =  LAPACKE_sggesx (   sggesx_obj->matrix_layout,
                                    sggesx_obj->jobvsl,
                                    sggesx_obj->jobvsr,
                                    sggesx_obj->sort,
									(LAPACK_S_SELECT3)sggesx_obj->select,
                                    sggesx_obj->sense,
                                    sggesx_obj->n,
                                    sggesx_obj->a,
                                    sggesx_obj->lda,
                                    sggesx_obj->b,
                                    sggesx_obj->ldb,
                                    &sggesx_obj->sdim,
                                    sggesx_obj->alphar,
                                    sggesx_obj->alphai,
                                    sggesx_obj->beta,
									sggesx_obj->vsl,
									sggesx_obj->ldvsl,
									sggesx_obj->vsr,
									sggesx_obj->ldvsr,
									sggesx_obj->rconde,
									sggesx_obj->rcondv
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    
    sggesx_obj->diff_a =  computeDiff_s( (sggesx_obj->lda)*(sggesx_obj->n), 
                sggesx_obj->a, sggesx_obj->aref );

    sggesx_obj->diff_b =  computeDiff_s( (sggesx_obj->ldb)*(sggesx_obj->n), 
                sggesx_obj->b, sggesx_obj->bref );

    sggesx_obj->diff_vsl =  computeDiff_s( (sggesx_obj->ldvsl)*(sggesx_obj->n), 
                sggesx_obj->vsl, sggesx_obj->vslref );

    sggesx_obj->diff_vsr =  computeDiff_s( (sggesx_obj->ldvsr)*(sggesx_obj->n), 
                sggesx_obj->vsr, sggesx_obj->vsrref );

    sggesx_obj->diff_alphar =  computeDiff_s( sggesx_obj->n, 
                sggesx_obj->alphar, sggesx_obj->alpharref );

    sggesx_obj->diff_alphai =  computeDiff_s( sggesx_obj->n, 
                sggesx_obj->alphai, sggesx_obj->alphairef );

    sggesx_obj->diff_beta =  computeDiff_s( sggesx_obj->n, 
                sggesx_obj->beta, sggesx_obj->betaref );

    sggesx_obj->diff_rconde =  computeDiff_s( 2, 
                sggesx_obj->rconde, sggesx_obj->rconderef );

    sggesx_obj->diff_rcondv =  computeDiff_s( 2, 
                sggesx_obj->rcondv, sggesx_obj->rcondvref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n ggesx float: \n diff_a: %f \n diff_b: %f \n \
diff_vsl: %f \n diff_vsr: %f \n diff_alphar: %f \n\
diff_alphai: %f \n diff_beta: %f \n ",
       sggesx_obj->diff_a, sggesx_obj->diff_b, sggesx_obj->diff_vsl,
	   sggesx_obj->diff_vsr, sggesx_obj->diff_alphar,
       sggesx_obj->diff_alphai, sggesx_obj->diff_beta );
#endif
}

TEST_F(sggesx_test, sggesx1) {
    //EXPECT_NEAR(0.0, sggesx_obj->diff_a, sggesx_obj->threshold);
    //EXPECT_NEAR(0.0, sggesx_obj->diff_b, sggesx_obj->threshold);
    EXPECT_NEAR(0.0, sggesx_obj->diff_vsl, sggesx_obj->threshold);
    EXPECT_NEAR(0.0, sggesx_obj->diff_vsr, sggesx_obj->threshold);
    EXPECT_NEAR(0.0, sggesx_obj->diff_alphar, sggesx_obj->threshold);
    EXPECT_NEAR(0.0, sggesx_obj->diff_alphai, sggesx_obj->threshold);
    EXPECT_NEAR(0.0, sggesx_obj->diff_beta, sggesx_obj->threshold);
    EXPECT_NEAR(0.0, sggesx_obj->diff_rconde, sggesx_obj->threshold);
    EXPECT_NEAR(0.0, sggesx_obj->diff_rcondv, sggesx_obj->threshold);
    EXPECT_EQ(sggesx_obj->sdim, sggesx_obj->sdimref);
}

TEST_F(sggesx_test, sggesx2) {
    //EXPECT_NEAR(0.0, sggesx_obj->diff_a, sggesx_obj->threshold);
    //EXPECT_NEAR(0.0, sggesx_obj->diff_b, sggesx_obj->threshold);
    EXPECT_NEAR(0.0, sggesx_obj->diff_vsl, sggesx_obj->threshold);
    EXPECT_NEAR(0.0, sggesx_obj->diff_vsr, sggesx_obj->threshold);
    EXPECT_NEAR(0.0, sggesx_obj->diff_alphar, sggesx_obj->threshold);
    EXPECT_NEAR(0.0, sggesx_obj->diff_alphai, sggesx_obj->threshold);
    EXPECT_NEAR(0.0, sggesx_obj->diff_beta, sggesx_obj->threshold);
    EXPECT_NEAR(0.0, sggesx_obj->diff_rconde, sggesx_obj->threshold);
    EXPECT_NEAR(0.0, sggesx_obj->diff_rcondv, sggesx_obj->threshold);
    EXPECT_EQ(sggesx_obj->sdim, sggesx_obj->sdimref);
}
TEST_F(sggesx_test, sggesx3) {
    //EXPECT_NEAR(0.0, sggesx_obj->diff_a, sggesx_obj->threshold);
    //EXPECT_NEAR(0.0, sggesx_obj->diff_b, sggesx_obj->threshold);
    EXPECT_NEAR(0.0, sggesx_obj->diff_vsl, sggesx_obj->threshold);
    EXPECT_NEAR(0.0, sggesx_obj->diff_vsr, sggesx_obj->threshold);
    EXPECT_NEAR(0.0, sggesx_obj->diff_alphar, sggesx_obj->threshold);
    EXPECT_NEAR(0.0, sggesx_obj->diff_alphai, sggesx_obj->threshold);
    EXPECT_NEAR(0.0, sggesx_obj->diff_beta, sggesx_obj->threshold);
    EXPECT_NEAR(0.0, sggesx_obj->diff_rconde, sggesx_obj->threshold);
    EXPECT_NEAR(0.0, sggesx_obj->diff_rcondv, sggesx_obj->threshold);
    EXPECT_EQ(sggesx_obj->sdim, sggesx_obj->sdimref);
}
TEST_F(sggesx_test, sggesx4) {
    //EXPECT_NEAR(0.0, sggesx_obj->diff_a, sggesx_obj->threshold);
    //EXPECT_NEAR(0.0, sggesx_obj->diff_b, sggesx_obj->threshold);
    EXPECT_NEAR(0.0, sggesx_obj->diff_vsl, sggesx_obj->threshold);
    EXPECT_NEAR(0.0, sggesx_obj->diff_vsr, sggesx_obj->threshold);
    EXPECT_NEAR(0.0, sggesx_obj->diff_alphar, sggesx_obj->threshold);
    EXPECT_NEAR(0.0, sggesx_obj->diff_alphai, sggesx_obj->threshold);
    EXPECT_NEAR(0.0, sggesx_obj->diff_beta, sggesx_obj->threshold);
    EXPECT_NEAR(0.0, sggesx_obj->diff_rconde, sggesx_obj->threshold);
    EXPECT_NEAR(0.0, sggesx_obj->diff_rcondv, sggesx_obj->threshold);
    EXPECT_EQ(sggesx_obj->sdim, sggesx_obj->sdimref);
}

/* Begin ggesx_double_common_parameters  class definition */
class ggesx_double_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    double diff_a, diff_b, diff_rconde, diff_rcondv;
    double diff_alphai, diff_alphar, diff_beta, diff_vsl, diff_vsr;

    float threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR	
    char jobvsl; // Must be 'N', or 'V'.
    char jobvsr; // Must be 'N', or 'V'.
	char sort; // Must be 'N', or 'S'.
	char sense; //  Must be 'N', 'E', 'V', or 'B'.
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
	double rconde[2], rcondv[2]; // contain the reciprocal condition numbers 
	double rconderef[2], rcondvref[2]; // contain the reciprocal condition numbers 
    /*Return Values */
    int info, inforef;

      ggesx_double_parameters (int matrix_layout_i, char jobvsl_i,
				char jobvsr_i, char sort_i, char sense_i, lapack_int n_i );
      ~ggesx_double_parameters ();
};

/* Constructor definition  double_common_parameters */
ggesx_double_parameters:: ggesx_double_parameters (int matrix_layout_i,
     char jobvsl_i, char jobvsr_i, char sort_i, char sense_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    jobvsl = jobvsl_i;
    jobvsr = jobvsr_i;
	sort = sort_i;
	sense = sense_i,
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
   printf(" \n ggesx double: matrix_layout: %d n: %d  jobvsl: %c \
jobvsr: %c \t sort: %c \t sense: %c\n", matrix_layout, n, jobvsr, 
jobvsl, sort, sense);
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
       EXPECT_FALSE( true) << "ggesx_double_parameters object: malloc error. Exiting ";
       ggesx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, lda*n );
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, ldb*n );
    lapacke_gtest_init_double_buffer_pair_with_constant( rconde, rconderef, 2, 0 );
    lapacke_gtest_init_double_buffer_pair_with_constant( rcondv, rcondvref, 2, 0 );
    lapacke_gtest_init_double_buffer_pair_with_constant( vsl, vslref, ldvsl*n, 0 );
    lapacke_gtest_init_double_buffer_pair_with_constant( vsr, vsrref, ldvsr*n, 0 );
    lapacke_gtest_init_double_buffer_pair_with_constant( alphar, alpharref, n, 0);
    lapacke_gtest_init_double_buffer_pair_with_constant( alphai, alphairef, n, 0);
    lapacke_gtest_init_double_buffer_pair_with_constant( beta, betaref, n, 0);
    lapacke_gtest_init_int_buffer_pair_with_constant(select, selectref, n, -1);

   } /* end of Constructor  */
    

/* Destructor definition  'ggesx_double_common_parameters' */
ggesx_double_parameters :: ~ggesx_double_parameters ()
{
   /* De-Allocate memory for the input matrices */
   ggesx_free();
} 

//  Test fixture class definition
class dggesx_test  : public  ::testing::Test {
public:
   ggesx_double_parameters  *dggesx_obj;
   void SetUp();  
   void TearDown () { delete dggesx_obj; }
};

void dggesx_test::SetUp()
{

    /* LAPACKE_dggesx prototype */	
    typedef int (*Fptr_NL_LAPACKE_dggesx) (int matrix_layout, char jobvsl, char jobvsr, char sort, LAPACK_D_SELECT3 select, char sense, lapack_int n, double* a, lapack_int lda, double* b, lapack_int ldb, lapack_int* sdim, double* alphar, double* alphai, double* beta, double* vsl, lapack_int ldvsl, double* vsr, lapack_int ldvsr, double* rconde, double* rcondv);
				 
    Fptr_NL_LAPACKE_dggesx DGGESX;

    dggesx_obj = new  ggesx_double_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].jobvsl,
										 eig_non_sym_paramslist[idx].jobvsr,
										 eig_non_sym_paramslist[idx].sort,
										 eig_non_sym_paramslist[idx].sense_ggesx,
                                         eig_paramslist[idx].n );
                                         
    dggesx_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    dggesx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dggesx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dggesx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dggesx_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DGGESX = (Fptr_NL_LAPACKE_dggesx)dlsym(dggesx_obj->hModule, "LAPACKE_dggesx");
    ASSERT_TRUE(DGGESX != NULL) << "failed to ppt the Netlib LAPACKE_dggesx symbol";

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    dggesx_obj->inforef = DGGESX(   dggesx_obj->matrix_layout,
                                    dggesx_obj->jobvsl,
                                    dggesx_obj->jobvsr,
                                    dggesx_obj->sort,
									(LAPACK_D_SELECT3)dggesx_obj->selectref,
                                    dggesx_obj->sense,
                                    dggesx_obj->n,
                                    dggesx_obj->aref,
                                    dggesx_obj->lda,
                                    dggesx_obj->bref,
                                    dggesx_obj->ldb,
                                    &dggesx_obj->sdimref,
                                    dggesx_obj->alpharref,
                                    dggesx_obj->alphairef,
                                    dggesx_obj->betaref,
									dggesx_obj->vslref,
									dggesx_obj->ldvsl,
									dggesx_obj->vsrref,
									dggesx_obj->ldvsr,
									dggesx_obj->rconderef,
									dggesx_obj->rcondvref
                                    );

    /* Compute libflame's Lapacke o/p  */
    dggesx_obj->inforef =  LAPACKE_dggesx (   dggesx_obj->matrix_layout,
                                    dggesx_obj->jobvsl,
                                    dggesx_obj->jobvsr,
                                    dggesx_obj->sort,
									(LAPACK_D_SELECT3)dggesx_obj->select,
                                    dggesx_obj->sense,
                                    dggesx_obj->n,
                                    dggesx_obj->a,
                                    dggesx_obj->lda,
                                    dggesx_obj->b,
                                    dggesx_obj->ldb,
                                    &dggesx_obj->sdim,
                                    dggesx_obj->alphar,
                                    dggesx_obj->alphai,
                                    dggesx_obj->beta,
									dggesx_obj->vsl,
									dggesx_obj->ldvsl,
									dggesx_obj->vsr,
									dggesx_obj->ldvsr,
									dggesx_obj->rconde,
									dggesx_obj->rcondv
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    
    dggesx_obj->diff_a =  computeDiff_d( (dggesx_obj->lda)*(dggesx_obj->n), 
                dggesx_obj->a, dggesx_obj->aref );

    dggesx_obj->diff_b =  computeDiff_d( (dggesx_obj->ldb)*(dggesx_obj->n), 
                dggesx_obj->b, dggesx_obj->bref );

    dggesx_obj->diff_vsl =  computeDiff_d( (dggesx_obj->ldvsl)*(dggesx_obj->n), 
                dggesx_obj->vsl, dggesx_obj->vslref );

    dggesx_obj->diff_vsr =  computeDiff_d( (dggesx_obj->ldvsr)*(dggesx_obj->n), 
                dggesx_obj->vsr, dggesx_obj->vsrref );

    dggesx_obj->diff_alphar =  computeDiff_d( dggesx_obj->n, 
                dggesx_obj->alphar, dggesx_obj->alpharref );

    dggesx_obj->diff_alphai =  computeDiff_d( dggesx_obj->n, 
                dggesx_obj->alphai, dggesx_obj->alphairef );

    dggesx_obj->diff_beta =  computeDiff_d( dggesx_obj->n, 
                dggesx_obj->beta, dggesx_obj->betaref );

    dggesx_obj->diff_rconde =  computeDiff_d( 2, 
                dggesx_obj->rconde, dggesx_obj->rconderef );

    dggesx_obj->diff_rcondv =  computeDiff_d( 2, 
                dggesx_obj->rcondv, dggesx_obj->rcondvref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n ggesx double: \n diff_a: %f \n diff_b: %f \n \
diff_vsl: %f \n diff_vsr: %f \n diff_alphar: %f \n\
diff_alphai: %f \n diff_beta: %f \n ",
       dggesx_obj->diff_a, dggesx_obj->diff_b, dggesx_obj->diff_vsl,
	   dggesx_obj->diff_vsr, dggesx_obj->diff_alphar,
       dggesx_obj->diff_alphai, dggesx_obj->diff_beta );
#endif
}

TEST_F(dggesx_test, dggesx1) {
    //EXPECT_NEAR(0.0, dggesx_obj->diff_a, dggesx_obj->threshold);
    //EXPECT_NEAR(0.0, dggesx_obj->diff_b, dggesx_obj->threshold);
    EXPECT_NEAR(0.0, dggesx_obj->diff_vsl, dggesx_obj->threshold);
    EXPECT_NEAR(0.0, dggesx_obj->diff_vsr, dggesx_obj->threshold);
    EXPECT_NEAR(0.0, dggesx_obj->diff_alphar, dggesx_obj->threshold);
    EXPECT_NEAR(0.0, dggesx_obj->diff_alphai, dggesx_obj->threshold);
    EXPECT_NEAR(0.0, dggesx_obj->diff_beta, dggesx_obj->threshold);
    EXPECT_NEAR(0.0, dggesx_obj->diff_rconde, dggesx_obj->threshold);
    EXPECT_NEAR(0.0, dggesx_obj->diff_rcondv, dggesx_obj->threshold);
    EXPECT_EQ(dggesx_obj->sdim, dggesx_obj->sdimref);
}

TEST_F(dggesx_test, dggesx2) {
    //EXPECT_NEAR(0.0, dggesx_obj->diff_a, dggesx_obj->threshold);
    //EXPECT_NEAR(0.0, dggesx_obj->diff_b, dggesx_obj->threshold);
    EXPECT_NEAR(0.0, dggesx_obj->diff_vsl, dggesx_obj->threshold);
    EXPECT_NEAR(0.0, dggesx_obj->diff_vsr, dggesx_obj->threshold);
    EXPECT_NEAR(0.0, dggesx_obj->diff_alphar, dggesx_obj->threshold);
    EXPECT_NEAR(0.0, dggesx_obj->diff_alphai, dggesx_obj->threshold);
    EXPECT_NEAR(0.0, dggesx_obj->diff_beta, dggesx_obj->threshold);
    EXPECT_NEAR(0.0, dggesx_obj->diff_rconde, dggesx_obj->threshold);
    EXPECT_NEAR(0.0, dggesx_obj->diff_rcondv, dggesx_obj->threshold);
    EXPECT_EQ(dggesx_obj->sdim, dggesx_obj->sdimref);
}
TEST_F(dggesx_test, dggesx3) {
    //EXPECT_NEAR(0.0, dggesx_obj->diff_a, dggesx_obj->threshold);
    //EXPECT_NEAR(0.0, dggesx_obj->diff_b, dggesx_obj->threshold);
    EXPECT_NEAR(0.0, dggesx_obj->diff_vsl, dggesx_obj->threshold);
    EXPECT_NEAR(0.0, dggesx_obj->diff_vsr, dggesx_obj->threshold);
    EXPECT_NEAR(0.0, dggesx_obj->diff_alphar, dggesx_obj->threshold);
    EXPECT_NEAR(0.0, dggesx_obj->diff_alphai, dggesx_obj->threshold);
    EXPECT_NEAR(0.0, dggesx_obj->diff_beta, dggesx_obj->threshold);
    EXPECT_NEAR(0.0, dggesx_obj->diff_rconde, dggesx_obj->threshold);
    EXPECT_NEAR(0.0, dggesx_obj->diff_rcondv, dggesx_obj->threshold);
    EXPECT_EQ(dggesx_obj->sdim, dggesx_obj->sdimref);
}
TEST_F(dggesx_test, dggesx4) {
    //EXPECT_NEAR(0.0, dggesx_obj->diff_a, dggesx_obj->threshold);
    //EXPECT_NEAR(0.0, dggesx_obj->diff_b, dggesx_obj->threshold);
    EXPECT_NEAR(0.0, dggesx_obj->diff_vsl, dggesx_obj->threshold);
    EXPECT_NEAR(0.0, dggesx_obj->diff_vsr, dggesx_obj->threshold);
    EXPECT_NEAR(0.0, dggesx_obj->diff_alphar, dggesx_obj->threshold);
    EXPECT_NEAR(0.0, dggesx_obj->diff_alphai, dggesx_obj->threshold);
    EXPECT_NEAR(0.0, dggesx_obj->diff_beta, dggesx_obj->threshold);
    EXPECT_NEAR(0.0, dggesx_obj->diff_rconde, dggesx_obj->threshold);
    EXPECT_NEAR(0.0, dggesx_obj->diff_rcondv, dggesx_obj->threshold);
    EXPECT_EQ(dggesx_obj->sdim, dggesx_obj->sdimref);
}

/* Begin ggesx_scomplex_common_parameters  class definition */
class ggesx_scomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_a, diff_b, diff_rconde, diff_rcondv;
    float diff_alphai, diff_alphar, diff_beta, diff_vsl, diff_vsr;

    float threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR	
    char jobvsl; // Must be 'N', or 'V'.
    char jobvsr; // Must be 'N', or 'V'.
	char sort; // Must be 'N', or 'S'.
	char sense; //  Must be 'N', 'E', 'V', or 'B'.
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
	float rconde[2], rcondv[2]; // contain the reciprocal condition numbers 
	float rconderef[2], rcondvref[2]; // contain the reciprocal condition numbers 
    /*Return Values */
    int info, inforef;

      ggesx_scomplex_parameters (int matrix_layout_i, char jobvsl_i,
				char jobvsr_i, char sort_i, char sense_i, lapack_int n_i );
      ~ggesx_scomplex_parameters ();
};

/* Constructor definition  lapack_complex_float_common_parameters */
ggesx_scomplex_parameters:: ggesx_scomplex_parameters (int matrix_layout_i,
     char jobvsl_i, char jobvsr_i, char sort_i, char sense_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    jobvsl = jobvsl_i;
    jobvsr = jobvsr_i;
	sort = sort_i;
	sense = sense_i,
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
   printf(" \n ggesx lapack_complex_float: matrix_layout: %d n: %d  jobvsl: %c \
jobvsr: %c \t sort: %c \t sense: %c\n", matrix_layout, n, jobvsr, 
jobvsl, sort, sense);
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
       EXPECT_FALSE( true) << "ggesx_scomplex_parameters object: malloc error. Exiting ";
       ggesx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, lda*n );
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, ldb*n );
    lapacke_gtest_init_float_buffer_pair_with_constant( rconde, rconderef, 2, 0 );
    lapacke_gtest_init_float_buffer_pair_with_constant( rcondv, rcondvref, 2, 0 );
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( vsl, vslref, ldvsl*n, 0 );
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( vsr, vsrref, ldvsr*n, 0 );
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( alphar, alpharref, n, 0);
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( alphai, alphairef, n, 0);
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( beta, betaref, n, 0);
    lapacke_gtest_init_int_buffer_pair_with_constant(select, selectref, n, -1);

   } /* end of Constructor  */
    

/* Destructor definition  'ggesx_scomplex_common_parameters' */
ggesx_scomplex_parameters :: ~ggesx_scomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   ggesx_free();
} 

//  Test fixture class definition
class cggesx_test  : public  ::testing::Test {
public:
   ggesx_scomplex_parameters  *cggesx_obj;
   void SetUp();  
   void TearDown () { delete cggesx_obj; }
};

void cggesx_test::SetUp()
{

    /* LAPACKE_cggesx prototype */	
    typedef int (*Fptr_NL_LAPACKE_cggesx) (int matrix_layout, char jobvsl, char jobvsr, char sort, LAPACK_C_SELECT2 select, char sense, lapack_int n, lapack_complex_float* a, lapack_int lda, lapack_complex_float* b, lapack_int ldb, lapack_int* sdim, lapack_complex_float* alpha, lapack_complex_float* beta, lapack_complex_float* vsl, lapack_int ldvsl, lapack_complex_float* vsr, lapack_int ldvsr, float* rconde, float* rcondv);
				 
    Fptr_NL_LAPACKE_cggesx CGGESX;

    cggesx_obj = new  ggesx_scomplex_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].jobvsl,
										 eig_non_sym_paramslist[idx].jobvsr,
										 eig_non_sym_paramslist[idx].sort,
										 eig_non_sym_paramslist[idx].sense_ggesx,
                                         eig_paramslist[idx].n );
                                         
    cggesx_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    cggesx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cggesx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cggesx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cggesx_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CGGESX = (Fptr_NL_LAPACKE_cggesx)dlsym(cggesx_obj->hModule, "LAPACKE_cggesx");
    ASSERT_TRUE(CGGESX != NULL) << "failed to ppt the Netlib LAPACKE_cggesx symbol";
	
    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    cggesx_obj->inforef = CGGESX(   cggesx_obj->matrix_layout,
                                    cggesx_obj->jobvsl,
                                    cggesx_obj->jobvsr,
                                    cggesx_obj->sort,
									(LAPACK_C_SELECT2)cggesx_obj->selectref,
                                    cggesx_obj->sense,
                                    cggesx_obj->n,
                                    cggesx_obj->aref,
                                    cggesx_obj->lda,
                                    cggesx_obj->bref,
                                    cggesx_obj->ldb,
                                    &cggesx_obj->sdimref,
                                    cggesx_obj->alpharref,
                                    cggesx_obj->betaref,
									cggesx_obj->vslref,
									cggesx_obj->ldvsl,
									cggesx_obj->vsrref,
									cggesx_obj->ldvsr,
									cggesx_obj->rconderef,
									cggesx_obj->rcondvref
                                    );

    /* Compute libflame's Lapacke o/p  */
    cggesx_obj->inforef =  LAPACKE_cggesx (   cggesx_obj->matrix_layout,
                                    cggesx_obj->jobvsl,
                                    cggesx_obj->jobvsr,
                                    cggesx_obj->sort,
									(LAPACK_C_SELECT2)cggesx_obj->select,
                                    cggesx_obj->sense,
                                    cggesx_obj->n,
                                    cggesx_obj->a,
                                    cggesx_obj->lda,
                                    cggesx_obj->b,
                                    cggesx_obj->ldb,
                                    &cggesx_obj->sdim,
                                    cggesx_obj->alphar,
                                    cggesx_obj->beta,
									cggesx_obj->vsl,
									cggesx_obj->ldvsl,
									cggesx_obj->vsr,
									cggesx_obj->ldvsr,
									cggesx_obj->rconde,
									cggesx_obj->rcondv
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    
    cggesx_obj->diff_a =  computeDiff_c( (cggesx_obj->lda)*(cggesx_obj->n), 
                cggesx_obj->a, cggesx_obj->aref );

    cggesx_obj->diff_b =  computeDiff_c( (cggesx_obj->ldb)*(cggesx_obj->n), 
                cggesx_obj->b, cggesx_obj->bref );

    cggesx_obj->diff_vsl =  computeDiff_c( (cggesx_obj->ldvsl)*(cggesx_obj->n), 
                cggesx_obj->vsl, cggesx_obj->vslref );

    cggesx_obj->diff_vsr =  computeDiff_c( (cggesx_obj->ldvsr)*(cggesx_obj->n), 
                cggesx_obj->vsr, cggesx_obj->vsrref );

    cggesx_obj->diff_alphar =  computeDiff_c( cggesx_obj->n, 
                cggesx_obj->alphar, cggesx_obj->alpharref );

    cggesx_obj->diff_beta =  computeDiff_c( cggesx_obj->n, 
                cggesx_obj->beta, cggesx_obj->betaref );

    cggesx_obj->diff_rconde =  computeDiff_s( 2, 
                cggesx_obj->rconde, cggesx_obj->rconderef );

    cggesx_obj->diff_rcondv =  computeDiff_s( 2, 
                cggesx_obj->rcondv, cggesx_obj->rcondvref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n ggesx lapack_complex_float: \n diff_a: %f \n diff_b: %f \n \
diff_vsl: %f \n diff_vsr: %f \n diff_alphar: %f \n\
 diff_beta: %f \n ",
       cggesx_obj->diff_a, cggesx_obj->diff_b, cggesx_obj->diff_vsl,
	   cggesx_obj->diff_vsr, cggesx_obj->diff_alphar,
        cggesx_obj->diff_beta );
#endif
}

TEST_F(cggesx_test, cggesx1) {
    //EXPECT_NEAR(0.0, cggesx_obj->diff_a, cggesx_obj->threshold);
    //EXPECT_NEAR(0.0, cggesx_obj->diff_b, cggesx_obj->threshold);
    EXPECT_NEAR(0.0, cggesx_obj->diff_vsl, cggesx_obj->threshold);
    EXPECT_NEAR(0.0, cggesx_obj->diff_vsr, cggesx_obj->threshold);
    EXPECT_NEAR(0.0, cggesx_obj->diff_alphar, cggesx_obj->threshold);
    EXPECT_NEAR(0.0, cggesx_obj->diff_beta, cggesx_obj->threshold);
    EXPECT_NEAR(0.0, cggesx_obj->diff_rconde, cggesx_obj->threshold);
    EXPECT_NEAR(0.0, cggesx_obj->diff_rcondv, cggesx_obj->threshold);
    EXPECT_EQ(cggesx_obj->sdim, cggesx_obj->sdimref);
}

TEST_F(cggesx_test, cggesx2) {
    //EXPECT_NEAR(0.0, cggesx_obj->diff_a, cggesx_obj->threshold);
    //EXPECT_NEAR(0.0, cggesx_obj->diff_b, cggesx_obj->threshold);
    EXPECT_NEAR(0.0, cggesx_obj->diff_vsl, cggesx_obj->threshold);
    EXPECT_NEAR(0.0, cggesx_obj->diff_vsr, cggesx_obj->threshold);
    EXPECT_NEAR(0.0, cggesx_obj->diff_alphar, cggesx_obj->threshold);
    EXPECT_NEAR(0.0, cggesx_obj->diff_beta, cggesx_obj->threshold);
    EXPECT_NEAR(0.0, cggesx_obj->diff_rconde, cggesx_obj->threshold);
    EXPECT_NEAR(0.0, cggesx_obj->diff_rcondv, cggesx_obj->threshold);
    EXPECT_EQ(cggesx_obj->sdim, cggesx_obj->sdimref);
}
TEST_F(cggesx_test, cggesx3) {
    //EXPECT_NEAR(0.0, cggesx_obj->diff_a, cggesx_obj->threshold);
    //EXPECT_NEAR(0.0, cggesx_obj->diff_b, cggesx_obj->threshold);
    EXPECT_NEAR(0.0, cggesx_obj->diff_vsl, cggesx_obj->threshold);
    EXPECT_NEAR(0.0, cggesx_obj->diff_vsr, cggesx_obj->threshold);
    EXPECT_NEAR(0.0, cggesx_obj->diff_alphar, cggesx_obj->threshold);
    EXPECT_NEAR(0.0, cggesx_obj->diff_beta, cggesx_obj->threshold);
    EXPECT_NEAR(0.0, cggesx_obj->diff_rconde, cggesx_obj->threshold);
    EXPECT_NEAR(0.0, cggesx_obj->diff_rcondv, cggesx_obj->threshold);
    EXPECT_EQ(cggesx_obj->sdim, cggesx_obj->sdimref);
}
TEST_F(cggesx_test, cggesx4) {
    //EXPECT_NEAR(0.0, cggesx_obj->diff_a, cggesx_obj->threshold);
    //EXPECT_NEAR(0.0, cggesx_obj->diff_b, cggesx_obj->threshold);
    EXPECT_NEAR(0.0, cggesx_obj->diff_vsl, cggesx_obj->threshold);
    EXPECT_NEAR(0.0, cggesx_obj->diff_vsr, cggesx_obj->threshold);
    EXPECT_NEAR(0.0, cggesx_obj->diff_alphar, cggesx_obj->threshold);
    EXPECT_NEAR(0.0, cggesx_obj->diff_beta, cggesx_obj->threshold);
    EXPECT_NEAR(0.0, cggesx_obj->diff_rconde, cggesx_obj->threshold);
    EXPECT_NEAR(0.0, cggesx_obj->diff_rcondv, cggesx_obj->threshold);
    EXPECT_EQ(cggesx_obj->sdim, cggesx_obj->sdimref);
}

/* Begin ggesx_dcomplex_common_parameters  class definition */
class ggesx_dcomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    double diff_a, diff_b, diff_rconde, diff_rcondv;
    double diff_alphai, diff_alphar, diff_beta, diff_vsl, diff_vsr;

    float threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR	
    char jobvsl; // Must be 'N', or 'V'.
    char jobvsr; // Must be 'N', or 'V'.
	char sort; // Must be 'N', or 'S'.
	char sense; //  Must be 'N', 'E', 'V', or 'B'.
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
	double rconde[2], rcondv[2]; // contain the reciprocal condition numbers 
	double rconderef[2], rcondvref[2]; // contain the reciprocal condition numbers 
    /*Return Values */
    int info, inforef;

      ggesx_dcomplex_parameters (int matrix_layout_i, char jobvsl_i,
				char jobvsr_i, char sort_i, char sense_i, lapack_int n_i );
      ~ggesx_dcomplex_parameters ();
};

/* Constructor definition  lapack_complex_double_common_parameters */
ggesx_dcomplex_parameters:: ggesx_dcomplex_parameters (int matrix_layout_i,
     char jobvsl_i, char jobvsr_i, char sort_i, char sense_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    jobvsl = jobvsl_i;
    jobvsr = jobvsr_i;
	sort = sort_i;
	sense = sense_i,
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
   printf(" \n ggesx lapack_complex_double: matrix_layout: %d n: %d  jobvsl: %c \
jobvsr: %c \t sort: %c \t sense: %c\n", matrix_layout, n, jobvsr, 
jobvsl, sort, sense);
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
       EXPECT_FALSE( true) << "ggesx_dcomplex_parameters object: malloc error. Exiting ";
       ggesx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, lda*n );
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, ldb*n );
    lapacke_gtest_init_double_buffer_pair_with_constant( rconde, rconderef, 2, 0 );
    lapacke_gtest_init_double_buffer_pair_with_constant( rcondv, rcondvref, 2, 0 );
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( vsl, vslref, ldvsl*n, 0 );
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( vsr, vsrref, ldvsr*n, 0 );
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( alphar, alpharref, n, 0);
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( alphai, alphairef, n, 0);
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( beta, betaref, n, 0);
    lapacke_gtest_init_int_buffer_pair_with_constant(select, selectref, n, -1);

   } /* end of Constructor  */
    

/* Destructor definition  'ggesx_dcomplex_common_parameters' */
ggesx_dcomplex_parameters :: ~ggesx_dcomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   ggesx_free();
} 

//  Test fixture class definition
class zggesx_test  : public  ::testing::Test {
public:
   ggesx_dcomplex_parameters  *zggesx_obj;
   void SetUp();  
   void TearDown () { delete zggesx_obj; }
};

void zggesx_test::SetUp()
{

    /* LAPACKE_zggesx prototype */	
    typedef int (*Fptr_NL_LAPACKE_zggesx) (int matrix_layout, char jobvsl, char jobvsr, char sort, LAPACK_Z_SELECT2 select, char sense, lapack_int n, lapack_complex_double* a, lapack_int lda, lapack_complex_double* b, lapack_int ldb, lapack_int* sdim, lapack_complex_double* alpha, lapack_complex_double* beta, lapack_complex_double* vsl, lapack_int ldvsl, lapack_complex_double* vsr, lapack_int ldvsr, double* rconde, double* rcondv);
				 
    Fptr_NL_LAPACKE_zggesx ZGGESX;

    zggesx_obj = new  ggesx_dcomplex_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].jobvsl,
										 eig_non_sym_paramslist[idx].jobvsr,
										 eig_non_sym_paramslist[idx].sort,
										 eig_non_sym_paramslist[idx].sense_ggesx,
                                         eig_paramslist[idx].n );
                                         
    zggesx_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    zggesx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zggesx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zggesx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zggesx_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZGGESX = (Fptr_NL_LAPACKE_zggesx)dlsym(zggesx_obj->hModule, "LAPACKE_zggesx");
    ASSERT_TRUE(ZGGESX != NULL) << "failed to ppt the Netlib LAPACKE_zggesx symbol";
	
    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    zggesx_obj->inforef = ZGGESX(   zggesx_obj->matrix_layout,
                                    zggesx_obj->jobvsl,
                                    zggesx_obj->jobvsr,
                                    zggesx_obj->sort,
									(LAPACK_Z_SELECT2)zggesx_obj->selectref,
                                    zggesx_obj->sense,
                                    zggesx_obj->n,
                                    zggesx_obj->aref,
                                    zggesx_obj->lda,
                                    zggesx_obj->bref,
                                    zggesx_obj->ldb,
                                    &zggesx_obj->sdimref,
                                    zggesx_obj->alpharref,
                                    zggesx_obj->betaref,
									zggesx_obj->vslref,
									zggesx_obj->ldvsl,
									zggesx_obj->vsrref,
									zggesx_obj->ldvsr,
									zggesx_obj->rconderef,
									zggesx_obj->rcondvref
                                    );

    /* Compute libflame's Lapacke o/p  */
    zggesx_obj->inforef =  LAPACKE_zggesx (   zggesx_obj->matrix_layout,
                                    zggesx_obj->jobvsl,
                                    zggesx_obj->jobvsr,
                                    zggesx_obj->sort,
									(LAPACK_Z_SELECT2)zggesx_obj->select,
                                    zggesx_obj->sense,
                                    zggesx_obj->n,
                                    zggesx_obj->a,
                                    zggesx_obj->lda,
                                    zggesx_obj->b,
                                    zggesx_obj->ldb,
                                    &zggesx_obj->sdim,
                                    zggesx_obj->alphar,
                                    zggesx_obj->beta,
									zggesx_obj->vsl,
									zggesx_obj->ldvsl,
									zggesx_obj->vsr,
									zggesx_obj->ldvsr,
									zggesx_obj->rconde,
									zggesx_obj->rcondv
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    
    zggesx_obj->diff_a =  computeDiff_z( (zggesx_obj->lda)*(zggesx_obj->n), 
                zggesx_obj->a, zggesx_obj->aref );

    zggesx_obj->diff_b =  computeDiff_z( (zggesx_obj->ldb)*(zggesx_obj->n), 
                zggesx_obj->b, zggesx_obj->bref );

    zggesx_obj->diff_vsl =  computeDiff_z( (zggesx_obj->ldvsl)*(zggesx_obj->n), 
                zggesx_obj->vsl, zggesx_obj->vslref );

    zggesx_obj->diff_vsr =  computeDiff_z( (zggesx_obj->ldvsr)*(zggesx_obj->n), 
                zggesx_obj->vsr, zggesx_obj->vsrref );

    zggesx_obj->diff_alphar =  computeDiff_z( zggesx_obj->n, 
                zggesx_obj->alphar, zggesx_obj->alpharref );

    zggesx_obj->diff_beta =  computeDiff_z( zggesx_obj->n, 
                zggesx_obj->beta, zggesx_obj->betaref );

    zggesx_obj->diff_rconde =  computeDiff_d( 2, 
                zggesx_obj->rconde, zggesx_obj->rconderef );

    zggesx_obj->diff_rcondv =  computeDiff_d( 2, 
                zggesx_obj->rcondv, zggesx_obj->rcondvref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n ggesx lapack_complex_double: \n diff_a: %f \n diff_b: %f \n \
diff_vsl: %f \n diff_vsr: %f \n diff_alphar: %f \n\
 diff_beta: %f \n ",
       zggesx_obj->diff_a, zggesx_obj->diff_b, zggesx_obj->diff_vsl,
	   zggesx_obj->diff_vsr, zggesx_obj->diff_alphar,
        zggesx_obj->diff_beta );
#endif
}

TEST_F(zggesx_test, zggesx1) {
    //EXPECT_NEAR(0.0, zggesx_obj->diff_a, zggesx_obj->threshold);
    //EXPECT_NEAR(0.0, zggesx_obj->diff_b, zggesx_obj->threshold);
    EXPECT_NEAR(0.0, zggesx_obj->diff_vsl, zggesx_obj->threshold);
    EXPECT_NEAR(0.0, zggesx_obj->diff_vsr, zggesx_obj->threshold);
    EXPECT_NEAR(0.0, zggesx_obj->diff_alphar, zggesx_obj->threshold);
    EXPECT_NEAR(0.0, zggesx_obj->diff_beta, zggesx_obj->threshold);
    EXPECT_NEAR(0.0, zggesx_obj->diff_rconde, zggesx_obj->threshold);
    EXPECT_NEAR(0.0, zggesx_obj->diff_rcondv, zggesx_obj->threshold);
    EXPECT_EQ(zggesx_obj->sdim, zggesx_obj->sdimref);
}

TEST_F(zggesx_test, zggesx2) {
    //EXPECT_NEAR(0.0, zggesx_obj->diff_a, zggesx_obj->threshold);
    //EXPECT_NEAR(0.0, zggesx_obj->diff_b, zggesx_obj->threshold);
    EXPECT_NEAR(0.0, zggesx_obj->diff_vsl, zggesx_obj->threshold);
    EXPECT_NEAR(0.0, zggesx_obj->diff_vsr, zggesx_obj->threshold);
    EXPECT_NEAR(0.0, zggesx_obj->diff_alphar, zggesx_obj->threshold);
    EXPECT_NEAR(0.0, zggesx_obj->diff_beta, zggesx_obj->threshold);
    EXPECT_NEAR(0.0, zggesx_obj->diff_rconde, zggesx_obj->threshold);
    EXPECT_NEAR(0.0, zggesx_obj->diff_rcondv, zggesx_obj->threshold);
    EXPECT_EQ(zggesx_obj->sdim, zggesx_obj->sdimref);
}
TEST_F(zggesx_test, zggesx3) {
    //EXPECT_NEAR(0.0, zggesx_obj->diff_a, zggesx_obj->threshold);
    //EXPECT_NEAR(0.0, zggesx_obj->diff_b, zggesx_obj->threshold);
    EXPECT_NEAR(0.0, zggesx_obj->diff_vsl, zggesx_obj->threshold);
    EXPECT_NEAR(0.0, zggesx_obj->diff_vsr, zggesx_obj->threshold);
    EXPECT_NEAR(0.0, zggesx_obj->diff_alphar, zggesx_obj->threshold);
    EXPECT_NEAR(0.0, zggesx_obj->diff_beta, zggesx_obj->threshold);
    EXPECT_NEAR(0.0, zggesx_obj->diff_rconde, zggesx_obj->threshold);
    EXPECT_NEAR(0.0, zggesx_obj->diff_rcondv, zggesx_obj->threshold);
    EXPECT_EQ(zggesx_obj->sdim, zggesx_obj->sdimref);
}
TEST_F(zggesx_test, zggesx4) {
    //EXPECT_NEAR(0.0, zggesx_obj->diff_a, zggesx_obj->threshold);
    //EXPECT_NEAR(0.0, zggesx_obj->diff_b, zggesx_obj->threshold);
    EXPECT_NEAR(0.0, zggesx_obj->diff_vsl, zggesx_obj->threshold);
    EXPECT_NEAR(0.0, zggesx_obj->diff_vsr, zggesx_obj->threshold);
    EXPECT_NEAR(0.0, zggesx_obj->diff_alphar, zggesx_obj->threshold);
    EXPECT_NEAR(0.0, zggesx_obj->diff_beta, zggesx_obj->threshold);
    EXPECT_NEAR(0.0, zggesx_obj->diff_rconde, zggesx_obj->threshold);
    EXPECT_NEAR(0.0, zggesx_obj->diff_rcondv, zggesx_obj->threshold);
    EXPECT_EQ(zggesx_obj->sdim, zggesx_obj->sdimref);
}
