#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define gees_free() \
    if (a!=NULL)        free(a); \
    if (aref!=NULL)     free(aref); \
    if (vs!=NULL)       free(vs); \
    if (vsref!=NULL)    free(vsref); \
    if (wr!=NULL)     free(wr); \
    if (wrref!=NULL)     free(wrref); \
    if (wi!=NULL)     free(wi); \
    if (wiref!=NULL)     free(wiref); \
    if (select!=NULL)     free(select); \
    if (selectref!=NULL)     free(selectref)


// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin gees_float_common_parameters  class definition */
class gees_float_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_a;
    float diff_wi, diff_wr, diff_vs;

    float threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
	
    char jobvs; // Must be 'N', or 'V'.
	char sort; // Must be 'N', or 'S'.
	int *select, *selectref; // Specifies the eigenvectors to be computed.  

    lapack_int n; // The order of the matrices A and B (n≥ 0).
    lapack_int lda; //  The leading dimension of a
    lapack_int ldvs;

    /* Input / Output parameters */
    float* a, *aref; // contains the n-by-n general matrix A.

    /* Output parameters */
    lapack_int sdim, sdimref;
    float *vs, *vsref; // Schur vectors.

    //  real and imaginary parts, respectively, of the computed eigenvalues
    float* wr, *wrref;
	float* wi, *wiref;
    /*Return Values */
    int info, inforef;

      gees_float_parameters (int matrix_layout_i, char jobvs_i,
						 char sort_i, lapack_int n_i );
      ~gees_float_parameters ();
};

/* Constructor definition  float_common_parameters */
gees_float_parameters:: gees_float_parameters (int matrix_layout_i,
            char jobvs_i,  char sort_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    jobvs = jobvs_i;
	sort = sort_i;
	
    n  = n_i;
    sdim = 0;
	sdimref = 0;

    lda = n;
    ldvs = n;

    hModule = NULL;
    dModule = NULL;

    diff_a = 0;
    diff_vs = 0;
    diff_wi = 0;
    diff_wr = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n gees float: matrix_layout: %d n: %d  jobvs: %c \
jobvsr: %c \n", matrix_layout, n, jobvsr, jobvs);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, lda*n );
    lapacke_gtest_alloc_float_buffer_pair( &wr, &wrref, n );
    lapacke_gtest_alloc_float_buffer_pair( &wi, &wiref, n );
    lapacke_gtest_alloc_int_buffer_pair( &select, &selectref, n );	
    lapacke_gtest_alloc_float_buffer_pair( &vs, &vsref, n*ldvs );


    if( (a==NULL) || (aref==NULL) ||  \
        (wr==NULL) || (wrref==NULL) || \
        (wi==NULL) || (wiref==NULL) || \
        (vs==NULL) || (vsref==NULL) || \
		(select==NULL) || (selectref==NULL) ){
       EXPECT_FALSE( true) << "gees_float_parameters object: malloc error. Exiting ";
       gees_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, lda*n );
    lapacke_gtest_init_float_buffer_pair_with_constant( vs, vsref, ldvs*n, 0 );
    lapacke_gtest_init_float_buffer_pair_with_constant( wr, wrref, n, 0);
    lapacke_gtest_init_float_buffer_pair_with_constant( wi, wiref, n, 0);
    lapacke_gtest_init_int_buffer_pair_with_constant(select, selectref, n, -1);

   } /* end of Constructor  */
    

/* Destructor definition  'gees_float_common_parameters' */
gees_float_parameters :: ~gees_float_parameters ()
{
   /* De-Allocate memory for the input matrices */
   gees_free();
} 

//  Test fixture class definition
class sgees_test  : public  ::testing::Test {
public:
   gees_float_parameters  *sgees_obj;
   void SetUp();  
   void TearDown () { delete sgees_obj; }
};

void sgees_test::SetUp()
{

    /* LAPACKE_sgees prototype */	
    typedef int (*Fptr_NL_LAPACKE_sgees) (int matrix_layout,
		char jobvs, char sort, LAPACK_S_SELECT2 select, lapack_int n,
		float* a, lapack_int lda, lapack_int* sdim, float* wr,
		float* wi, float* vs, lapack_int ldvs);
				 
    Fptr_NL_LAPACKE_sgees SGEES;

    sgees_obj = new  gees_float_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].jobvsl,
										 eig_non_sym_paramslist[idx].sort_gees,
                                         eig_paramslist[idx].n );
                                         
    sgees_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    sgees_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgees_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgees_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgees_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SGEES = (Fptr_NL_LAPACKE_sgees)dlsym(sgees_obj->hModule, "LAPACKE_sgees");
    ASSERT_TRUE(SGEES != NULL) << "failed to ppt the Netlib LAPACKE_sgees symbol";

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    sgees_obj->inforef = SGEES(   sgees_obj->matrix_layout,
                                    sgees_obj->jobvs,
                                    sgees_obj->sort,
									(LAPACK_S_SELECT2)sgees_obj->selectref,
                                    sgees_obj->n,
                                    sgees_obj->aref,
                                    sgees_obj->lda,
                                    &sgees_obj->sdimref,
                                    sgees_obj->wrref,
                                    sgees_obj->wiref,
									sgees_obj->vsref,
									sgees_obj->ldvs
                                    );

    /* Compute libflame's Lapacke o/p  */
    sgees_obj->inforef =  LAPACKE_sgees (   sgees_obj->matrix_layout,
                                    sgees_obj->jobvs,
                                    sgees_obj->sort,
									(LAPACK_S_SELECT2)sgees_obj->select,
                                    sgees_obj->n,
                                    sgees_obj->a,
                                    sgees_obj->lda,
                                    &sgees_obj->sdim,
                                    sgees_obj->wr,
                                    sgees_obj->wi,
									sgees_obj->vs,
									sgees_obj->ldvs
                                    );
    /* Capture the Netlib, libflame o/p buffers' differences */
    
    sgees_obj->diff_a =  computeDiff_s( (sgees_obj->lda)*(sgees_obj->n), 
                sgees_obj->a, sgees_obj->aref );


    sgees_obj->diff_vs =  computeDiff_s( (sgees_obj->ldvs)*(sgees_obj->n), 
                sgees_obj->vs, sgees_obj->vsref );

    sgees_obj->diff_wr =  computeDiff_s( sgees_obj->n, 
                sgees_obj->wr, sgees_obj->wrref );

    sgees_obj->diff_wi =  computeDiff_s( sgees_obj->n, 
                sgees_obj->wi, sgees_obj->wiref );


#if LAPACKE_TEST_VERBOSE
   printf(" \n gees float: \n diff_a: %f \n \
diff_vs: %f \n diff_wr: %f \n diff_wi: %f \n ",
       sgees_obj->diff_a, sgees_obj->diff_vs,
       sgees_obj->diff_wr, sgees_obj->diff_wi );
#endif
}

TEST_F(sgees_test, sgees1) {
    EXPECT_NEAR(0.0, sgees_obj->diff_vs, sgees_obj->threshold);
    EXPECT_NEAR(0.0, sgees_obj->diff_wr, sgees_obj->threshold);
    EXPECT_NEAR(0.0, sgees_obj->diff_wi, sgees_obj->threshold);
    EXPECT_EQ(sgees_obj->sdim, sgees_obj->sdimref);
}

TEST_F(sgees_test, sgees2) {
    EXPECT_NEAR(0.0, sgees_obj->diff_vs, sgees_obj->threshold);
    EXPECT_NEAR(0.0, sgees_obj->diff_wr, sgees_obj->threshold);
    EXPECT_NEAR(0.0, sgees_obj->diff_wi, sgees_obj->threshold);
    EXPECT_EQ(sgees_obj->sdim, sgees_obj->sdimref);
}
TEST_F(sgees_test, sgees3) {
    EXPECT_NEAR(0.0, sgees_obj->diff_vs, sgees_obj->threshold);
    EXPECT_NEAR(0.0, sgees_obj->diff_wr, sgees_obj->threshold);
    EXPECT_NEAR(0.0, sgees_obj->diff_wi, sgees_obj->threshold);
    EXPECT_EQ(sgees_obj->sdim, sgees_obj->sdimref);
}
TEST_F(sgees_test, sgees4) {
    EXPECT_NEAR(0.0, sgees_obj->diff_vs, sgees_obj->threshold);
    EXPECT_NEAR(0.0, sgees_obj->diff_wr, sgees_obj->threshold);
    EXPECT_NEAR(0.0, sgees_obj->diff_wi, sgees_obj->threshold);
    EXPECT_EQ(sgees_obj->sdim, sgees_obj->sdimref);
}

/* Begin gees_double_common_parameters  class definition */
class gees_double_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    double diff_a;
    double diff_wi, diff_wr, diff_vs;

    double threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
	
    char jobvs; // Must be 'N', or 'V'.
	char sort; // Must be 'N', or 'S'.
	int *select, *selectref; // Specifies the eigenvectors to be computed.  

    lapack_int n; // The order of the matrices A and B (n≥ 0).
    lapack_int lda; //  The leading dimension of a
    lapack_int ldvs;

    /* Input / Output parameters */
    double* a, *aref; // contains the n-by-n general matrix A.

    /* Output parameters */
    lapack_int sdim, sdimref;
    double *vs, *vsref; // Schur vectors.

    //  real and imaginary parts, respectively, of the computed eigenvalues
    double* wr, *wrref;
	double* wi, *wiref;
    /*Return Values */
    int info, inforef;

      gees_double_parameters (int matrix_layout_i, char jobvs_i,
						 char sort_i, lapack_int n_i );
      ~gees_double_parameters ();
};

/* Constructor definition  double_common_parameters */
gees_double_parameters:: gees_double_parameters (int matrix_layout_i,
            char jobvs_i,  char sort_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    jobvs = jobvs_i;
	sort = sort_i;
	
    n  = n_i;
    sdim = 0;
	sdimref = 0;

    lda = n;
    ldvs = n;

    hModule = NULL;
    dModule = NULL;

    diff_a = 0;
    diff_vs = 0;
    diff_wi = 0;
    diff_wr = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n gees double: matrix_layout: %d n: %d  jobvs: %c \
jobvsr: %c \n", matrix_layout, n, jobvsr, jobvs);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, lda*n );
    lapacke_gtest_alloc_double_buffer_pair( &wr, &wrref, n );
    lapacke_gtest_alloc_double_buffer_pair( &wi, &wiref, n );
    lapacke_gtest_alloc_int_buffer_pair( &select, &selectref, n );	
    lapacke_gtest_alloc_double_buffer_pair( &vs, &vsref, n*ldvs );


    if( (a==NULL) || (aref==NULL) ||  \
        (wr==NULL) || (wrref==NULL) || \
        (wi==NULL) || (wiref==NULL) || \
        (vs==NULL) || (vsref==NULL) || \
		(select==NULL) || (selectref==NULL) ){
       EXPECT_FALSE( true) << "gees_double_parameters object: malloc error. Exiting ";
       gees_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, lda*n );
    lapacke_gtest_init_double_buffer_pair_with_constant( vs, vsref, ldvs*n, 0 );
    lapacke_gtest_init_double_buffer_pair_with_constant( wr, wrref, n, 0);
    lapacke_gtest_init_double_buffer_pair_with_constant( wi, wiref, n, 0);
    lapacke_gtest_init_int_buffer_pair_with_constant(select, selectref, n, -1);

   } /* end of Constructor  */
    

/* Destructor definition  'gees_double_common_parameters' */
gees_double_parameters :: ~gees_double_parameters ()
{
   /* De-Allocate memory for the input matrices */
   gees_free();
} 

//  Test fixture class definition
class dgees_test  : public  ::testing::Test {
public:
   gees_double_parameters  *dgees_obj;
   void SetUp();  
   void TearDown () { delete dgees_obj; }
};

void dgees_test::SetUp()
{

    /* LAPACKE_dgees prototype */	
    typedef int (*Fptr_NL_LAPACKE_dgees) (int matrix_layout,
		char jobvs, char sort, LAPACK_D_SELECT2 select, lapack_int n,
		double* a, lapack_int lda, lapack_int* sdim, double* wr,
		double* wi, double* vs, lapack_int ldvs);
				 
    Fptr_NL_LAPACKE_dgees DGEES;

    dgees_obj = new  gees_double_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].jobvsl,
										 eig_non_sym_paramslist[idx].sort_gees,
                                         eig_paramslist[idx].n );
                                         
    dgees_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    dgees_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgees_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgees_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgees_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DGEES = (Fptr_NL_LAPACKE_dgees)dlsym(dgees_obj->hModule, "LAPACKE_dgees");
    ASSERT_TRUE(DGEES != NULL) << "failed to ppt the Netlib LAPACKE_dgees symbol";

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    dgees_obj->inforef = DGEES(   dgees_obj->matrix_layout,
                                    dgees_obj->jobvs,
                                    dgees_obj->sort,
									(LAPACK_D_SELECT2)dgees_obj->selectref,
                                    dgees_obj->n,
                                    dgees_obj->aref,
                                    dgees_obj->lda,
                                    &dgees_obj->sdimref,
                                    dgees_obj->wrref,
                                    dgees_obj->wiref,
									dgees_obj->vsref,
									dgees_obj->ldvs
                                    );

    /* Compute libflame's Lapacke o/p  */
    dgees_obj->inforef =  LAPACKE_dgees (   dgees_obj->matrix_layout,
                                    dgees_obj->jobvs,
                                    dgees_obj->sort,
									(LAPACK_D_SELECT2)dgees_obj->select,
                                    dgees_obj->n,
                                    dgees_obj->a,
                                    dgees_obj->lda,
                                    &dgees_obj->sdim,
                                    dgees_obj->wr,
                                    dgees_obj->wi,
									dgees_obj->vs,
									dgees_obj->ldvs
                                    );
    /* Capture the Netlib, libflame o/p buffers' differences */
    
    dgees_obj->diff_a =  computeDiff_d( (dgees_obj->lda)*(dgees_obj->n), 
                dgees_obj->a, dgees_obj->aref );


    dgees_obj->diff_vs =  computeDiff_d( (dgees_obj->ldvs)*(dgees_obj->n), 
                dgees_obj->vs, dgees_obj->vsref );

    dgees_obj->diff_wr =  computeDiff_d( dgees_obj->n, 
                dgees_obj->wr, dgees_obj->wrref );

    dgees_obj->diff_wi =  computeDiff_d( dgees_obj->n, 
                dgees_obj->wi, dgees_obj->wiref );


#if LAPACKE_TEST_VERBOSE
   printf(" \n gees double: \n diff_a: %f \n \
diff_vs: %f \n diff_wr: %f \n diff_wi: %f \n ",
       dgees_obj->diff_a, dgees_obj->diff_vs,
       dgees_obj->diff_wr, dgees_obj->diff_wi );
#endif
}

TEST_F(dgees_test, dgees1) {
    EXPECT_NEAR(0.0, dgees_obj->diff_vs, dgees_obj->threshold);
    EXPECT_NEAR(0.0, dgees_obj->diff_wr, dgees_obj->threshold);
    EXPECT_NEAR(0.0, dgees_obj->diff_wi, dgees_obj->threshold);
    EXPECT_EQ(dgees_obj->sdim, dgees_obj->sdimref);
}

TEST_F(dgees_test, dgees2) {
    EXPECT_NEAR(0.0, dgees_obj->diff_vs, dgees_obj->threshold);
    EXPECT_NEAR(0.0, dgees_obj->diff_wr, dgees_obj->threshold);
    EXPECT_NEAR(0.0, dgees_obj->diff_wi, dgees_obj->threshold);
    EXPECT_EQ(dgees_obj->sdim, dgees_obj->sdimref);
}
TEST_F(dgees_test, dgees3) {
    EXPECT_NEAR(0.0, dgees_obj->diff_vs, dgees_obj->threshold);
    EXPECT_NEAR(0.0, dgees_obj->diff_wr, dgees_obj->threshold);
    EXPECT_NEAR(0.0, dgees_obj->diff_wi, dgees_obj->threshold);
    EXPECT_EQ(dgees_obj->sdim, dgees_obj->sdimref);
}
TEST_F(dgees_test, dgees4) {
    EXPECT_NEAR(0.0, dgees_obj->diff_vs, dgees_obj->threshold);
    EXPECT_NEAR(0.0, dgees_obj->diff_wr, dgees_obj->threshold);
    EXPECT_NEAR(0.0, dgees_obj->diff_wi, dgees_obj->threshold);
    EXPECT_EQ(dgees_obj->sdim, dgees_obj->sdimref);
}

/* Begin gees_scomplex_common_parameters  class definition */
class gees_scomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_a;
    float diff_wi, diff_wr, diff_vs;

    float threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
	
    char jobvs; // Must be 'N', or 'V'.
	char sort; // Must be 'N', or 'S'.
	int *select, *selectref; // Specifies the eigenvectors to be computed.  

    lapack_int n; // The order of the matrices A and B (n≥ 0).
    lapack_int lda; //  The leading dimension of a
    lapack_int ldvs;

    /* Input / Output parameters */
    lapack_complex_float* a, *aref; // contains the n-by-n general matrix A.

    /* Output parameters */
    lapack_int sdim, sdimref;
    lapack_complex_float *vs, *vsref; // Schur vectors.

    //  real and imaginary parts, respectively, of the computed eigenvalues
    lapack_complex_float* wr, *wrref;
	lapack_complex_float* wi, *wiref;
    /*Return Values */
    int info, inforef;

      gees_scomplex_parameters (int matrix_layout_i, char jobvs_i,
						 char sort_i, lapack_int n_i );
      ~gees_scomplex_parameters ();
};

/* Constructor definition  lapack_complex_float_common_parameters */
gees_scomplex_parameters:: gees_scomplex_parameters (int matrix_layout_i,
            char jobvs_i,  char sort_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    jobvs = jobvs_i;
	sort = sort_i;
	
    n  = n_i;
    sdim = 0;
	sdimref = 0;

    lda = n;
    ldvs = n;

    hModule = NULL;
    dModule = NULL;

    diff_a = 0;
    diff_vs = 0;
    diff_wi = 0;
    diff_wr = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n gees lapack_complex_float: matrix_layout: %d n: %d  jobvs: %c \
jobvsr: %c \n", matrix_layout, n, jobvsr, jobvs);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, lda*n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &wr, &wrref, n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &wi, &wiref, n );
    lapacke_gtest_alloc_int_buffer_pair( &select, &selectref, n );	
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &vs, &vsref, n*ldvs );


    if( (a==NULL) || (aref==NULL) ||  \
        (wr==NULL) || (wrref==NULL) || \
        (wi==NULL) || (wiref==NULL) || \
        (vs==NULL) || (vsref==NULL) || \
		(select==NULL) || (selectref==NULL) ){
       EXPECT_FALSE( true) << "gees_scomplex_parameters object: malloc error. Exiting ";
       gees_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, lda*n );
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( vs, vsref, ldvs*n, 0 );
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( wr, wrref, n, 0);
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( wi, wiref, n, 0);
    lapacke_gtest_init_int_buffer_pair_with_constant(select, selectref, n, -1);

   } /* end of Constructor  */
    

/* Destructor definition  'gees_scomplex_common_parameters' */
gees_scomplex_parameters :: ~gees_scomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   gees_free();
} 

//  Test fixture class definition
class cgees_test  : public  ::testing::Test {
public:
   gees_scomplex_parameters  *cgees_obj;
   void SetUp();  
   void TearDown () { delete cgees_obj; }
};

void cgees_test::SetUp()
{

    /* LAPACKE_cgees prototype */	
    typedef int (*Fptr_NL_LAPACKE_cgees) (int matrix_layout, char jobvs,
		char sort, LAPACK_C_SELECT1 select, lapack_int n, 
		lapack_complex_float* a, lapack_int lda, lapack_int* sdim, 
		lapack_complex_float* w, lapack_complex_float* vs, lapack_int ldvs);
				 
    Fptr_NL_LAPACKE_cgees CGEES;

    cgees_obj = new  gees_scomplex_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].jobvsl,
										 eig_non_sym_paramslist[idx].sort_gees,
                                         eig_paramslist[idx].n );
                                         
    cgees_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    cgees_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgees_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgees_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgees_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CGEES = (Fptr_NL_LAPACKE_cgees)dlsym(cgees_obj->hModule, "LAPACKE_cgees");
    ASSERT_TRUE(CGEES != NULL) << "failed to ppt the Netlib LAPACKE_cgees symbol";

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    cgees_obj->inforef = CGEES(   cgees_obj->matrix_layout,
                                    cgees_obj->jobvs,
                                    cgees_obj->sort,
									(LAPACK_C_SELECT1)cgees_obj->selectref,
                                    cgees_obj->n,
                                    cgees_obj->aref,
                                    cgees_obj->lda,
                                    &cgees_obj->sdimref,
                                    cgees_obj->wrref,
									cgees_obj->vsref,
									cgees_obj->ldvs
                                    );

    /* Compute libflame's Lapacke o/p  */
    cgees_obj->inforef =  LAPACKE_cgees (   cgees_obj->matrix_layout,
                                    cgees_obj->jobvs,
                                    cgees_obj->sort,
									(LAPACK_C_SELECT1)cgees_obj->select,
                                    cgees_obj->n,
                                    cgees_obj->a,
                                    cgees_obj->lda,
                                    &cgees_obj->sdim,
                                    cgees_obj->wr,
									cgees_obj->vs,
									cgees_obj->ldvs
                                    );
    /* Capture the Netlib, libflame o/p buffers' differences */
    
    cgees_obj->diff_a =  computeDiff_c( (cgees_obj->lda)*(cgees_obj->n), 
                cgees_obj->a, cgees_obj->aref );


    cgees_obj->diff_vs =  computeDiff_c( (cgees_obj->ldvs)*(cgees_obj->n), 
                cgees_obj->vs, cgees_obj->vsref );

    cgees_obj->diff_wr =  computeDiff_c( cgees_obj->n, 
                cgees_obj->wr, cgees_obj->wrref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n gees lapack_complex_float: \n diff_a: %f \n \
diff_vs: %f \n diff_wr: %f \n diff_wi: %f \n ",
       cgees_obj->diff_a, cgees_obj->diff_vs,
       cgees_obj->diff_wr, cgees_obj->diff_wi );
#endif
}

TEST_F(cgees_test, cgees1) {
    EXPECT_NEAR(0.0, cgees_obj->diff_vs, cgees_obj->threshold);
    EXPECT_NEAR(0.0, cgees_obj->diff_wr, cgees_obj->threshold);
    EXPECT_EQ(cgees_obj->sdim, cgees_obj->sdimref);
}

TEST_F(cgees_test, cgees2) {
    EXPECT_NEAR(0.0, cgees_obj->diff_vs, cgees_obj->threshold);
    EXPECT_NEAR(0.0, cgees_obj->diff_wr, cgees_obj->threshold);
    EXPECT_EQ(cgees_obj->sdim, cgees_obj->sdimref);
}
TEST_F(cgees_test, cgees3) {
    EXPECT_NEAR(0.0, cgees_obj->diff_vs, cgees_obj->threshold);
    EXPECT_NEAR(0.0, cgees_obj->diff_wr, cgees_obj->threshold);
    EXPECT_EQ(cgees_obj->sdim, cgees_obj->sdimref);
}
TEST_F(cgees_test, cgees4) {
    EXPECT_NEAR(0.0, cgees_obj->diff_vs, cgees_obj->threshold);
    EXPECT_NEAR(0.0, cgees_obj->diff_wr, cgees_obj->threshold);
    EXPECT_EQ(cgees_obj->sdim, cgees_obj->sdimref);
}

/* Begin gees_dcomplex_common_parameters  class definition */
class gees_dcomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    double diff_a;
    double diff_wi, diff_wr, diff_vs;

    double threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
	
    char jobvs; // Must be 'N', or 'V'.
	char sort; // Must be 'N', or 'S'.
	int *select, *selectref; // Specifies the eigenvectors to be computed.  

    lapack_int n; // The order of the matrices A and B (n≥ 0).
    lapack_int lda; //  The leading dimension of a
    lapack_int ldvs;

    /* Input / Output parameters */
    lapack_complex_double* a, *aref; // contains the n-by-n general matrix A.

    /* Output parameters */
    lapack_int sdim, sdimref;
    lapack_complex_double *vs, *vsref; // Schur vectors.

    //  real and imaginary parts, respectively, of the computed eigenvalues
    lapack_complex_double* wr, *wrref;
	lapack_complex_double* wi, *wiref;
    /*Return Values */
    int info, inforef;

      gees_dcomplex_parameters (int matrix_layout_i, char jobvs_i,
						 char sort_i, lapack_int n_i );
      ~gees_dcomplex_parameters ();
};

/* Constructor definition  lapack_complex_double_common_parameters */
gees_dcomplex_parameters:: gees_dcomplex_parameters (int matrix_layout_i,
            char jobvs_i,  char sort_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    jobvs = jobvs_i;
	sort = sort_i;
	
    n  = n_i;
    sdim = 0;
	sdimref = 0;

    lda = n;
    ldvs = n;

    hModule = NULL;
    dModule = NULL;

    diff_a = 0;
    diff_vs = 0;
    diff_wi = 0;
    diff_wr = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n gees lapack_complex_double: matrix_layout: %d n: %d  jobvs: %c \
jobvsr: %c \n", matrix_layout, n, jobvsr, jobvs);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, lda*n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &wr, &wrref, n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &wi, &wiref, n );
    lapacke_gtest_alloc_int_buffer_pair( &select, &selectref, n );	
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &vs, &vsref, n*ldvs );


    if( (a==NULL) || (aref==NULL) ||  \
        (wr==NULL) || (wrref==NULL) || \
        (wi==NULL) || (wiref==NULL) || \
        (vs==NULL) || (vsref==NULL) || \
		(select==NULL) || (selectref==NULL) ){
       EXPECT_FALSE( true) << "gees_dcomplex_parameters object: malloc error. Exiting ";
       gees_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, lda*n );
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( vs, vsref, ldvs*n, 0 );
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( wr, wrref, n, 0);
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( wi, wiref, n, 0);
    lapacke_gtest_init_int_buffer_pair_with_constant(select, selectref, n, -1);

   } /* end of Constructor  */
    

/* Destructor definition  'gees_dcomplex_common_parameters' */
gees_dcomplex_parameters :: ~gees_dcomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   gees_free();
} 

//  Test fixture class definition
class zgees_test  : public  ::testing::Test {
public:
   gees_dcomplex_parameters  *zgees_obj;
   void SetUp();  
   void TearDown () { delete zgees_obj; }
};

void zgees_test::SetUp()
{

    /* LAPACKE_zgees prototype */	
    typedef int (*Fptr_NL_LAPACKE_zgees) (int matrix_layout, char jobvs,
			char sort, LAPACK_Z_SELECT1 select, lapack_int n,
			lapack_complex_double* a, lapack_int lda, lapack_int* sdim,
			lapack_complex_double* w, lapack_complex_double* vs, 
			lapack_int ldvs);
				 
    Fptr_NL_LAPACKE_zgees ZGEES;

    zgees_obj = new  gees_dcomplex_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].jobvsl,
										 eig_non_sym_paramslist[idx].sort_gees,
                                         eig_paramslist[idx].n );
                                         
    zgees_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    zgees_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgees_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgees_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgees_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZGEES = (Fptr_NL_LAPACKE_zgees)dlsym(zgees_obj->hModule, "LAPACKE_zgees");
    ASSERT_TRUE(ZGEES != NULL) << "failed to ppt the Netlib LAPACKE_zgees symbol";

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    zgees_obj->inforef = ZGEES(   zgees_obj->matrix_layout,
                                    zgees_obj->jobvs,
                                    zgees_obj->sort,
									(LAPACK_Z_SELECT1)zgees_obj->selectref,
                                    zgees_obj->n,
                                    zgees_obj->aref,
                                    zgees_obj->lda,
                                    &zgees_obj->sdimref,
                                    zgees_obj->wrref,
									zgees_obj->vsref,
									zgees_obj->ldvs
                                    );

    /* Compute libflame's Lapacke o/p  */
    zgees_obj->inforef =  LAPACKE_zgees (   zgees_obj->matrix_layout,
                                    zgees_obj->jobvs,
                                    zgees_obj->sort,
									(LAPACK_Z_SELECT1)zgees_obj->select,
                                    zgees_obj->n,
                                    zgees_obj->a,
                                    zgees_obj->lda,
                                    &zgees_obj->sdim,
                                    zgees_obj->wr,
									zgees_obj->vs,
									zgees_obj->ldvs
                                    );
    /* Capture the Netlib, libflame o/p buffers' differences */
    
    zgees_obj->diff_a =  computeDiff_z( (zgees_obj->lda)*(zgees_obj->n), 
                zgees_obj->a, zgees_obj->aref );


    zgees_obj->diff_vs =  computeDiff_z( (zgees_obj->ldvs)*(zgees_obj->n), 
                zgees_obj->vs, zgees_obj->vsref );

    zgees_obj->diff_wr =  computeDiff_z( zgees_obj->n, 
                zgees_obj->wr, zgees_obj->wrref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n gees lapack_complex_double: \n diff_a: %f \n \
diff_vs: %f \n diff_wr: %f \n diff_wi: %f \n ",
       zgees_obj->diff_a, zgees_obj->diff_vs,
       zgees_obj->diff_wr, zgees_obj->diff_wi );
#endif
}

TEST_F(zgees_test, zgees1) {
    EXPECT_NEAR(0.0, zgees_obj->diff_vs, zgees_obj->threshold);
    EXPECT_NEAR(0.0, zgees_obj->diff_wr, zgees_obj->threshold);
    EXPECT_EQ(zgees_obj->sdim, zgees_obj->sdimref);
}

TEST_F(zgees_test, zgees2) {
    EXPECT_NEAR(0.0, zgees_obj->diff_vs, zgees_obj->threshold);
    EXPECT_NEAR(0.0, zgees_obj->diff_wr, zgees_obj->threshold);
    EXPECT_EQ(zgees_obj->sdim, zgees_obj->sdimref);
}
TEST_F(zgees_test, zgees3) {
    EXPECT_NEAR(0.0, zgees_obj->diff_vs, zgees_obj->threshold);
    EXPECT_NEAR(0.0, zgees_obj->diff_wr, zgees_obj->threshold);
    EXPECT_EQ(zgees_obj->sdim, zgees_obj->sdimref);
}
TEST_F(zgees_test, zgees4) {
    EXPECT_NEAR(0.0, zgees_obj->diff_vs, zgees_obj->threshold);
    EXPECT_NEAR(0.0, zgees_obj->diff_wr, zgees_obj->threshold);
    EXPECT_EQ(zgees_obj->sdim, zgees_obj->sdimref);
}
