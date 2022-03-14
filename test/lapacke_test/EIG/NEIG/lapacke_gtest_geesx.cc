#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define LAPACKE_TEST_VERBOSE  (1)

#define geesx_free() \
    if (a!=NULL)        free(a); \
    if (aref!=NULL)     free(aref); \
    if (vs!=NULL)       free(vs); \
    if (vsref!=NULL)    free(vsref); \
    if (wr!=NULL)     free(wr); \
    if (wrref!=NULL)     free(wrref); \
    if (wi!=NULL)     free(wi); \
    if (wiref!=NULL)     free(wiref); \
    if (select!=NULL)     free(select); \
    if (selectref!=NULL)     free(selectref); \
    if (rconde!=NULL)   free(rconde); \
    if (rconderef!=NULL)  free(rconderef); \
    if (rcondv!=NULL)   free(rcondv); \
    if (rcondvref!=NULL)  free(rcondvref)


// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin geesx_scomplex_common_parameters  class definition */
class geesx_scomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_a, diff_rconde, diff_rcondv;
    float diff_wi, diff_wr, diff_vs;

    float threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR	
    char jobvs; // Must be 'N', or 'V'.
	char sort; // Must be 'N', or 'S'.
	char sense; //  Must be 'N', 'E', 'V', or 'B'.
	int *select, *selectref; // Specifies the eigenvectors to be computed.  

    lapack_int n; // The order of the matrices A.
    lapack_int lda; //  The leading dimension of a
    lapack_int ldvs;

    /* Input / Output parameters */
    lapack_complex_float* a, *aref; // contains the n-by-n general matrix A.

    /* Output parameters */
    lapack_int sdim, sdimref;
    lapack_complex_float *vs, *vsref; // orthogonal/unitary matrix Z of Schur vectors
    lapack_complex_float* wr, *wrref;
	lapack_complex_float* wi, *wiref;
	float *rconde, *rcondv; // contain the reciprocal condition numbers 
	float *rconderef, *rcondvref; // contain the reciprocal condition numbers 
    /*Return Values */
    int info, inforef;

      geesx_scomplex_parameters (int matrix_layout_i, char jobvs_i,
				    char sort_i, char sense_i, lapack_int n_i );
      ~geesx_scomplex_parameters ();
};

/* Constructor definition  lapack_complex_float_common_parameters */
geesx_scomplex_parameters:: geesx_scomplex_parameters (int matrix_layout_i,
           char jobvs_i, char sort_i, char sense_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    jobvs = jobvs_i;
	sort = sort_i;
	sense = sense_i,
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
   printf(" \n geesx lapack_complex_float: matrix_layout: %d n: %d  jobvs: %c \
sort: %c \t sense: %c\n", matrix_layout, n, jobvs, 
sort, sense);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, lda*n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &wr, &wrref, n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &wi, &wiref, n );
    lapacke_gtest_alloc_int_buffer_pair( &select, &selectref, n );	
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &vs, &vsref, n*ldvs );
    lapacke_gtest_alloc_float_buffer_pair( &rconde, &rconderef, n );
    lapacke_gtest_alloc_float_buffer_pair( &rcondv, &rcondvref, n );


    if( (a==NULL) || (aref==NULL) ||  \
        (wr==NULL) || (wrref==NULL) || \
        (wi==NULL) || (wiref==NULL) || \
        (vs==NULL) || (vsref==NULL) || \
        (rconde==NULL) || (rconderef==NULL) || \
        (rcondv==NULL) || (rcondvref==NULL) || \
		(select==NULL) || (selectref==NULL) ){
       EXPECT_FALSE( true) << "geesx_scomplex_parameters object: malloc error. Exiting ";
       geesx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, lda*n );
    lapacke_gtest_init_float_buffer_pair_with_constant( rconde, rconderef, n, 0 );
    lapacke_gtest_init_float_buffer_pair_with_constant( rcondv, rcondvref, n, 0 );
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( vs, vsref, ldvs*n, 0 );
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( wr, wrref, n, 0);
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( wi, wiref, n, 0);
    lapacke_gtest_init_int_buffer_pair_with_constant(select, selectref, n, -1);

   } /* end of Constructor  */
    

/* Destructor definition  'geesx_scomplex_common_parameters' */
geesx_scomplex_parameters :: ~geesx_scomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   geesx_free();
} 

//  Test fixture class definition
class cgeesx_test  : public  ::testing::Test {
public:
   geesx_scomplex_parameters  *cgeesx_obj;
   void SetUp();  
   void TearDown () { delete cgeesx_obj; }
};

void cgeesx_test::SetUp()
{

    /* LAPACKE_cgeesx prototype */	
    typedef int (*Fptr_NL_LAPACKE_cgeesx) (int matrix_layout, char jobvs,
		char sort, LAPACK_C_SELECT1 select, char sense, lapack_int n,
		lapack_complex_float* a, lapack_int lda, lapack_int* sdim,
		lapack_complex_float* w, lapack_complex_float* vs,
		lapack_int ldvs, float* rconde, float* rcondv);
				 
    Fptr_NL_LAPACKE_cgeesx CGEESX;

    cgeesx_obj = new  geesx_scomplex_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].jobvsl,
										 eig_non_sym_paramslist[idx].sort,
										 eig_non_sym_paramslist[idx].sense_ggesx,
                                         eig_paramslist[idx].n );
                                         
    cgeesx_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    cgeesx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgeesx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgeesx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgeesx_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CGEESX = (Fptr_NL_LAPACKE_cgeesx)dlsym(cgeesx_obj->hModule, "LAPACKE_cgeesx");
    ASSERT_TRUE(CGEESX != NULL) << "failed to ppt the Netlib LAPACKE_cgeesx symbol";

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    cgeesx_obj->inforef = CGEESX(   cgeesx_obj->matrix_layout,
                                    cgeesx_obj->jobvs,
                                    cgeesx_obj->sort,
									(LAPACK_C_SELECT1)cgeesx_obj->selectref,
                                    cgeesx_obj->sense,
                                    cgeesx_obj->n,
                                    cgeesx_obj->aref,
                                    cgeesx_obj->lda,
                                    &cgeesx_obj->sdimref,
                                    cgeesx_obj->wrref,
									cgeesx_obj->vsref,
									cgeesx_obj->ldvs,
									cgeesx_obj->rconderef,
									cgeesx_obj->rcondvref
                                    );

    /* Compute libflame's Lapacke o/p  */
    cgeesx_obj->inforef =  LAPACKE_cgeesx (   cgeesx_obj->matrix_layout,
                                    cgeesx_obj->jobvs,
                                    cgeesx_obj->sort,
									(LAPACK_C_SELECT1)cgeesx_obj->select,
                                    cgeesx_obj->sense,
                                    cgeesx_obj->n,
                                    cgeesx_obj->a,
                                    cgeesx_obj->lda,
                                    &cgeesx_obj->sdim,
                                    cgeesx_obj->wr,
									cgeesx_obj->vs,
									cgeesx_obj->ldvs,
									cgeesx_obj->rconde,
									cgeesx_obj->rcondv
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    
    cgeesx_obj->diff_a =  computeDiff_c( (cgeesx_obj->lda)*(cgeesx_obj->n), 
                cgeesx_obj->a, cgeesx_obj->aref );

    cgeesx_obj->diff_vs =  computeDiff_c( (cgeesx_obj->ldvs)*(cgeesx_obj->n), 
                cgeesx_obj->vs, cgeesx_obj->vsref );

    cgeesx_obj->diff_wr =  computeDiff_c( cgeesx_obj->n, 
                cgeesx_obj->wr, cgeesx_obj->wrref );

    cgeesx_obj->diff_rconde =  computeDiff_s( 2, 
                cgeesx_obj->rconde, cgeesx_obj->rconderef );

    cgeesx_obj->diff_rcondv =  computeDiff_s( 2, 
                cgeesx_obj->rcondv, cgeesx_obj->rcondvref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n geesx lapack_complex_float: \n diff_a: %f \n \
diff_vs: %f \n diff_wr: %f \n diff_wi: %f \n ",cgeesx_obj->diff_a,
cgeesx_obj->diff_vs, cgeesx_obj->diff_wr, cgeesx_obj->diff_wi );
#endif
}

TEST_F(cgeesx_test, cgeesx1) {
    //EXPECT_NEAR(0.0, cgeesx_obj->diff_a, cgeesx_obj->threshold);
    EXPECT_NEAR(0.0, cgeesx_obj->diff_vs, cgeesx_obj->threshold);
    EXPECT_NEAR(0.0, cgeesx_obj->diff_wr, cgeesx_obj->threshold);
    EXPECT_NEAR(0.0, cgeesx_obj->diff_rconde, cgeesx_obj->threshold);
    EXPECT_NEAR(0.0, cgeesx_obj->diff_rcondv, cgeesx_obj->threshold);
    EXPECT_EQ(cgeesx_obj->sdim, cgeesx_obj->sdimref);
}

TEST_F(cgeesx_test, cgeesx2) {
    //EXPECT_NEAR(0.0, cgeesx_obj->diff_a, cgeesx_obj->threshold);
    EXPECT_NEAR(0.0, cgeesx_obj->diff_vs, cgeesx_obj->threshold);
    EXPECT_NEAR(0.0, cgeesx_obj->diff_wr, cgeesx_obj->threshold);
    EXPECT_NEAR(0.0, cgeesx_obj->diff_rconde, cgeesx_obj->threshold);
    EXPECT_NEAR(0.0, cgeesx_obj->diff_rcondv, cgeesx_obj->threshold);
    EXPECT_EQ(cgeesx_obj->sdim, cgeesx_obj->sdimref);
}
TEST_F(cgeesx_test, cgeesx3) {
    //EXPECT_NEAR(0.0, cgeesx_obj->diff_a, cgeesx_obj->threshold);
    EXPECT_NEAR(0.0, cgeesx_obj->diff_vs, cgeesx_obj->threshold);
    EXPECT_NEAR(0.0, cgeesx_obj->diff_wr, cgeesx_obj->threshold);
    EXPECT_NEAR(0.0, cgeesx_obj->diff_rconde, cgeesx_obj->threshold);
    EXPECT_NEAR(0.0, cgeesx_obj->diff_rcondv, cgeesx_obj->threshold);
    EXPECT_EQ(cgeesx_obj->sdim, cgeesx_obj->sdimref);
}
TEST_F(cgeesx_test, cgeesx4) {
    //EXPECT_NEAR(0.0, cgeesx_obj->diff_a, cgeesx_obj->threshold);
    EXPECT_NEAR(0.0, cgeesx_obj->diff_vs, cgeesx_obj->threshold);
    EXPECT_NEAR(0.0, cgeesx_obj->diff_wr, cgeesx_obj->threshold);
    EXPECT_NEAR(0.0, cgeesx_obj->diff_rconde, cgeesx_obj->threshold);
    EXPECT_NEAR(0.0, cgeesx_obj->diff_rcondv, cgeesx_obj->threshold);
    EXPECT_EQ(cgeesx_obj->sdim, cgeesx_obj->sdimref);
}

/* Begin geesx_dcomplex_common_parameters  class definition */
class geesx_dcomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    double diff_a, diff_rconde, diff_rcondv;
    double diff_wi, diff_wr, diff_vs;

    float threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR	
    char jobvs; // Must be 'N', or 'V'.
	char sort; // Must be 'N', or 'S'.
	char sense; //  Must be 'N', 'E', 'V', or 'B'.
	int *select, *selectref; // Specifies the eigenvectors to be computed.  

    lapack_int n; // The order of the matrices A.
    lapack_int lda; //  The leading dimension of a
    lapack_int ldvs;

    /* Input / Output parameters */
    lapack_complex_double* a, *aref; // contains the n-by-n general matrix A.

    /* Output parameters */
    lapack_int sdim, sdimref;
    lapack_complex_double *vs, *vsref; // orthogonal/unitary matrix Z of Schur vectors
    lapack_complex_double* wr, *wrref;
	lapack_complex_double* wi, *wiref;
	double *rconde, *rcondv; // contain the reciprocal condition numbers 
	double *rconderef, *rcondvref; // contain the reciprocal condition numbers 
    /*Return Values */
    int info, inforef;

      geesx_dcomplex_parameters (int matrix_layout_i, char jobvs_i,
				    char sort_i, char sense_i, lapack_int n_i );
      ~geesx_dcomplex_parameters ();
};

/* Constructor definition  lapack_complex_double_common_parameters */
geesx_dcomplex_parameters:: geesx_dcomplex_parameters (int matrix_layout_i,
           char jobvs_i, char sort_i, char sense_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    jobvs = jobvs_i;
	sort = sort_i;
	sense = sense_i,
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
   printf(" \n geesx lapack_complex_double: matrix_layout: %d n: %d  jobvs: %c \
sort: %c \t sense: %c\n", matrix_layout, n, jobvs, 
sort, sense);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, lda*n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &wr, &wrref, n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &wi, &wiref, n );
    lapacke_gtest_alloc_int_buffer_pair( &select, &selectref, n );	
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &vs, &vsref, n*ldvs );
    lapacke_gtest_alloc_double_buffer_pair( &rconde, &rconderef, n );
    lapacke_gtest_alloc_double_buffer_pair( &rcondv, &rcondvref, n );


    if( (a==NULL) || (aref==NULL) ||  \
        (wr==NULL) || (wrref==NULL) || \
        (wi==NULL) || (wiref==NULL) || \
        (vs==NULL) || (vsref==NULL) || \
        (rconde==NULL) || (rconderef==NULL) || \
        (rcondv==NULL) || (rcondvref==NULL) || \
		(select==NULL) || (selectref==NULL) ){
       EXPECT_FALSE( true) << "geesx_dcomplex_parameters object: malloc error. Exiting ";
       geesx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, lda*n );
    lapacke_gtest_init_double_buffer_pair_with_constant( rconde, rconderef, n, 0 );
    lapacke_gtest_init_double_buffer_pair_with_constant( rcondv, rcondvref, n, 0 );
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( vs, vsref, ldvs*n, 0 );
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( wr, wrref, n, 0);
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( wi, wiref, n, 0);
    lapacke_gtest_init_int_buffer_pair_with_constant(select, selectref, n, -1);

   } /* end of Constructor  */
    

/* Destructor definition  'geesx_dcomplex_common_parameters' */
geesx_dcomplex_parameters :: ~geesx_dcomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   geesx_free();
} 

//  Test fixture class definition
class zgeesx_test  : public  ::testing::Test {
public:
   geesx_dcomplex_parameters  *zgeesx_obj;
   void SetUp();  
   void TearDown () { delete zgeesx_obj; }
};

void zgeesx_test::SetUp()
{

    /* LAPACKE_zgeesx prototype */	
    typedef int (*Fptr_NL_LAPACKE_zgeesx) (int matrix_layout, char jobvs,
		char sort, LAPACK_Z_SELECT1 select, char sense, lapack_int n,
		lapack_complex_double* a, lapack_int lda, lapack_int* sdim,
		lapack_complex_double* w, lapack_complex_double* vs,
		lapack_int ldvs, double* rconde, double* rcondv);
				 
    Fptr_NL_LAPACKE_zgeesx ZGEESX;

    zgeesx_obj = new  geesx_dcomplex_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].jobvsl,
										 eig_non_sym_paramslist[idx].sort,
										 eig_non_sym_paramslist[idx].sense_ggesx,
                                         eig_paramslist[idx].n );
                                         
    zgeesx_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    zgeesx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgeesx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgeesx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgeesx_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZGEESX = (Fptr_NL_LAPACKE_zgeesx)dlsym(zgeesx_obj->hModule, "LAPACKE_zgeesx");
    ASSERT_TRUE(ZGEESX != NULL) << "failed to ppt the Netlib LAPACKE_zgeesx symbol";

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    zgeesx_obj->inforef = ZGEESX(   zgeesx_obj->matrix_layout,
                                    zgeesx_obj->jobvs,
                                    zgeesx_obj->sort,
									(LAPACK_Z_SELECT1)zgeesx_obj->selectref,
                                    zgeesx_obj->sense,
                                    zgeesx_obj->n,
                                    zgeesx_obj->aref,
                                    zgeesx_obj->lda,
                                    &zgeesx_obj->sdimref,
                                    zgeesx_obj->wrref,
									zgeesx_obj->vsref,
									zgeesx_obj->ldvs,
									zgeesx_obj->rconderef,
									zgeesx_obj->rcondvref
                                    );

    /* Compute libflame's Lapacke o/p  */
    zgeesx_obj->inforef =  LAPACKE_zgeesx (   zgeesx_obj->matrix_layout,
                                    zgeesx_obj->jobvs,
                                    zgeesx_obj->sort,
									(LAPACK_Z_SELECT1)zgeesx_obj->select,
                                    zgeesx_obj->sense,
                                    zgeesx_obj->n,
                                    zgeesx_obj->a,
                                    zgeesx_obj->lda,
                                    &zgeesx_obj->sdim,
                                    zgeesx_obj->wr,
									zgeesx_obj->vs,
									zgeesx_obj->ldvs,
									zgeesx_obj->rconde,
									zgeesx_obj->rcondv
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    
    zgeesx_obj->diff_a =  computeDiff_z( (zgeesx_obj->lda)*(zgeesx_obj->n), 
                zgeesx_obj->a, zgeesx_obj->aref );

    zgeesx_obj->diff_vs =  computeDiff_z( (zgeesx_obj->ldvs)*(zgeesx_obj->n), 
                zgeesx_obj->vs, zgeesx_obj->vsref );

    zgeesx_obj->diff_wr =  computeDiff_z( zgeesx_obj->n, 
                zgeesx_obj->wr, zgeesx_obj->wrref );

    zgeesx_obj->diff_rconde =  computeDiff_d( 2, 
                zgeesx_obj->rconde, zgeesx_obj->rconderef );

    zgeesx_obj->diff_rcondv =  computeDiff_d( 2, 
                zgeesx_obj->rcondv, zgeesx_obj->rcondvref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n geesx lapack_complex_double: \n diff_a: %f \n \
diff_vs: %f \n diff_wr: %f \n diff_wi: %f \n ",zgeesx_obj->diff_a,
zgeesx_obj->diff_vs, zgeesx_obj->diff_wr, zgeesx_obj->diff_wi );
#endif
}

TEST_F(zgeesx_test, zgeesx1) {
    //EXPECT_NEAR(0.0, zgeesx_obj->diff_a, zgeesx_obj->threshold);
    EXPECT_NEAR(0.0, zgeesx_obj->diff_vs, zgeesx_obj->threshold);
    EXPECT_NEAR(0.0, zgeesx_obj->diff_wr, zgeesx_obj->threshold);
    EXPECT_NEAR(0.0, zgeesx_obj->diff_rconde, zgeesx_obj->threshold);
    EXPECT_NEAR(0.0, zgeesx_obj->diff_rcondv, zgeesx_obj->threshold);
    EXPECT_EQ(zgeesx_obj->sdim, zgeesx_obj->sdimref);
}

TEST_F(zgeesx_test, zgeesx2) {
    //EXPECT_NEAR(0.0, zgeesx_obj->diff_a, zgeesx_obj->threshold);
    EXPECT_NEAR(0.0, zgeesx_obj->diff_vs, zgeesx_obj->threshold);
    EXPECT_NEAR(0.0, zgeesx_obj->diff_wr, zgeesx_obj->threshold);
    EXPECT_NEAR(0.0, zgeesx_obj->diff_rconde, zgeesx_obj->threshold);
    EXPECT_NEAR(0.0, zgeesx_obj->diff_rcondv, zgeesx_obj->threshold);
    EXPECT_EQ(zgeesx_obj->sdim, zgeesx_obj->sdimref);
}
TEST_F(zgeesx_test, zgeesx3) {
    //EXPECT_NEAR(0.0, zgeesx_obj->diff_a, zgeesx_obj->threshold);
    EXPECT_NEAR(0.0, zgeesx_obj->diff_vs, zgeesx_obj->threshold);
    EXPECT_NEAR(0.0, zgeesx_obj->diff_wr, zgeesx_obj->threshold);
    EXPECT_NEAR(0.0, zgeesx_obj->diff_rconde, zgeesx_obj->threshold);
    EXPECT_NEAR(0.0, zgeesx_obj->diff_rcondv, zgeesx_obj->threshold);
    EXPECT_EQ(zgeesx_obj->sdim, zgeesx_obj->sdimref);
}
TEST_F(zgeesx_test, zgeesx4) {
    //EXPECT_NEAR(0.0, zgeesx_obj->diff_a, zgeesx_obj->threshold);
    EXPECT_NEAR(0.0, zgeesx_obj->diff_vs, zgeesx_obj->threshold);
    EXPECT_NEAR(0.0, zgeesx_obj->diff_wr, zgeesx_obj->threshold);
    EXPECT_NEAR(0.0, zgeesx_obj->diff_rconde, zgeesx_obj->threshold);
    EXPECT_NEAR(0.0, zgeesx_obj->diff_rcondv, zgeesx_obj->threshold);
    EXPECT_EQ(zgeesx_obj->sdim, zgeesx_obj->sdimref);
}

#if 0
/* Begin geesx_float_common_parameters  class definition */
class geesx_float_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_a, diff_rconde, diff_rcondv;
    float diff_wi, diff_wr, diff_vs;

    float threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR	
    char jobvs; // Must be 'N', or 'V'.
	char sort; // Must be 'N', or 'S'.
	char sense; //  Must be 'N', 'E', 'V', or 'B'.
	int *select, *selectref; // Specifies the eigenvectors to be computed.  

    lapack_int n; // The order of the matrices A.
    lapack_int lda; //  The leading dimension of a
    lapack_int ldvs;

    /* Input / Output parameters */
    float* a, *aref; // contains the n-by-n general matrix A.

    /* Output parameters */
    lapack_int sdim, sdimref;
    float *vs, *vsref; // orthogonal/unitary matrix Z of Schur vectors
    float* wr, *wrref;
	float* wi, *wiref;
	float *rconde, *rcondv; // contain the reciprocal condition numbers 
	float *rconderef, *rcondvref; // contain the reciprocal condition numbers 
    /*Return Values */
    int info, inforef;

      geesx_float_parameters (int matrix_layout_i, char jobvs_i,
				    char sort_i, char sense_i, lapack_int n_i );
      ~geesx_float_parameters ();
};

/* Constructor definition  float_common_parameters */
geesx_float_parameters:: geesx_float_parameters (int matrix_layout_i,
           char jobvs_i, char sort_i, char sense_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    jobvs = jobvs_i;
	sort = sort_i;
	sense = sense_i,
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
   printf(" \n geesx float: matrix_layout: %d n: %d  jobvs: %c \
sort: %c \t sense: %c\n", matrix_layout, n, jobvs, 
sort, sense);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, lda*n );
    lapacke_gtest_alloc_float_buffer_pair( &wr, &wrref, n );
    lapacke_gtest_alloc_float_buffer_pair( &wi, &wiref, n );
    lapacke_gtest_alloc_int_buffer_pair( &select, &selectref, n );	
    lapacke_gtest_alloc_float_buffer_pair( &vs, &vsref, n*ldvs );
    lapacke_gtest_alloc_float_buffer_pair( &rconde, &rconderef, n );
    lapacke_gtest_alloc_float_buffer_pair( &rcondv, &rcondvref, n );


    if( (a==NULL) || (aref==NULL) ||  \
        (wr==NULL) || (wrref==NULL) || \
        (wi==NULL) || (wiref==NULL) || \
        (vs==NULL) || (vsref==NULL) || \
        (rconde==NULL) || (rconderef==NULL) || \
        (rcondv==NULL) || (rcondvref==NULL) || \
		(select==NULL) || (selectref==NULL) ){
       EXPECT_FALSE( true) << "geesx_float_parameters object: malloc error. Exiting ";
       geesx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, lda*n );
    lapacke_gtest_init_float_buffer_pair_with_constant( rconde, rconderef, n, 0 );
    lapacke_gtest_init_float_buffer_pair_with_constant( rcondv, rcondvref, n, 0 );
    lapacke_gtest_init_float_buffer_pair_with_constant( vs, vsref, ldvs*n, 0 );
    lapacke_gtest_init_float_buffer_pair_with_constant( wr, wrref, n, 0);
    lapacke_gtest_init_float_buffer_pair_with_constant( wi, wiref, n, 0);
    lapacke_gtest_init_int_buffer_pair_with_constant(select, selectref, n, -1);

   } /* end of Constructor  */
    

/* Destructor definition  'geesx_float_common_parameters' */
geesx_float_parameters :: ~geesx_float_parameters ()
{
   /* De-Allocate memory for the input matrices */
   geesx_free();
} 

//  Test fixture class definition
class sgeesx_test  : public  ::testing::Test {
public:
   geesx_float_parameters  *sgeesx_obj;
   void SetUp();  
   void TearDown () { delete sgeesx_obj; }
};

void sgeesx_test::SetUp()
{

    /* LAPACKE_sgeesx prototype */	
    typedef int (*Fptr_NL_LAPACKE_sgeesx) (int matrix_layout, char jobvs,
		char sort, LAPACK_S_SELECT2 select, char sense, lapack_int n,
		float* a, lapack_int lda, lapack_int* sdim, float* wr, float* wi,
		float* vs, lapack_int ldvs, float* rconde, float* rcondv);
				 
    Fptr_NL_LAPACKE_sgeesx SGGESX;

    sgeesx_obj = new  geesx_float_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].jobvsl,
										 eig_non_sym_paramslist[idx].sort,
										 eig_non_sym_paramslist[idx].sense_ggesx,
                                         eig_paramslist[idx].n );
                                         
    sgeesx_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    sgeesx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgeesx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgeesx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgeesx_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SGGESX = (Fptr_NL_LAPACKE_sgeesx)dlsym(sgeesx_obj->hModule, "LAPACKE_sgeesx");
    ASSERT_TRUE(SGGESX != NULL) << "failed to ppt the Netlib LAPACKE_sgeesx symbol";

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    sgeesx_obj->inforef = SGGESX(   sgeesx_obj->matrix_layout,
                                    sgeesx_obj->jobvs,
                                    sgeesx_obj->sort,
									(LAPACK_S_SELECT2)sgeesx_obj->selectref,
                                    sgeesx_obj->sense,
                                    sgeesx_obj->n,
                                    sgeesx_obj->aref,
                                    sgeesx_obj->lda,
                                    &sgeesx_obj->sdimref,
                                    sgeesx_obj->wrref,
                                    sgeesx_obj->wiref,
									sgeesx_obj->vsref,
									sgeesx_obj->ldvs,
									sgeesx_obj->rconderef,
									sgeesx_obj->rcondvref
                                    );

    /* Compute libflame's Lapacke o/p  */
    sgeesx_obj->inforef =  LAPACKE_sgeesx (   sgeesx_obj->matrix_layout,
                                    sgeesx_obj->jobvs,
                                    sgeesx_obj->sort,
									(LAPACK_S_SELECT2)sgeesx_obj->select,
                                    sgeesx_obj->sense,
                                    sgeesx_obj->n,
                                    sgeesx_obj->a,
                                    sgeesx_obj->lda,
                                    &sgeesx_obj->sdim,
                                    sgeesx_obj->wr,
                                    sgeesx_obj->wi,
									sgeesx_obj->vs,
									sgeesx_obj->ldvs,
									sgeesx_obj->rconde,
									sgeesx_obj->rcondv
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    
    sgeesx_obj->diff_a =  computeDiff_s( (sgeesx_obj->lda)*(sgeesx_obj->n), 
                sgeesx_obj->a, sgeesx_obj->aref );

    sgeesx_obj->diff_vs =  computeDiff_s( (sgeesx_obj->ldvs)*(sgeesx_obj->n), 
                sgeesx_obj->vs, sgeesx_obj->vsref );

    sgeesx_obj->diff_wr =  computeDiff_s( sgeesx_obj->n, 
                sgeesx_obj->wr, sgeesx_obj->wrref );

    sgeesx_obj->diff_wi =  computeDiff_s( sgeesx_obj->n, 
                sgeesx_obj->wi, sgeesx_obj->wiref );

    sgeesx_obj->diff_rconde =  computeDiff_s( 2, 
                sgeesx_obj->rconde, sgeesx_obj->rconderef );

    sgeesx_obj->diff_rcondv =  computeDiff_s( 2, 
                sgeesx_obj->rcondv, sgeesx_obj->rcondvref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n geesx float: \n diff_a: %f \n \
diff_vs: %f \n diff_wr: %f \n diff_wi: %f \n ",sgeesx_obj->diff_a,
sgeesx_obj->diff_vs, sgeesx_obj->diff_wr, sgeesx_obj->diff_wi );
#endif
}

TEST_F(sgeesx_test, sgeesx1) {
    //EXPECT_NEAR(0.0, sgeesx_obj->diff_a, sgeesx_obj->threshold);
    EXPECT_NEAR(0.0, sgeesx_obj->diff_vs, sgeesx_obj->threshold);
    EXPECT_NEAR(0.0, sgeesx_obj->diff_wr, sgeesx_obj->threshold);
    EXPECT_NEAR(0.0, sgeesx_obj->diff_wi, sgeesx_obj->threshold);
    EXPECT_NEAR(0.0, sgeesx_obj->diff_rconde, sgeesx_obj->threshold);
    EXPECT_NEAR(0.0, sgeesx_obj->diff_rcondv, sgeesx_obj->threshold);
    EXPECT_EQ(sgeesx_obj->sdim, sgeesx_obj->sdimref);
}

TEST_F(sgeesx_test, sgeesx2) {
    //EXPECT_NEAR(0.0, sgeesx_obj->diff_a, sgeesx_obj->threshold);
    EXPECT_NEAR(0.0, sgeesx_obj->diff_vs, sgeesx_obj->threshold);
    EXPECT_NEAR(0.0, sgeesx_obj->diff_wr, sgeesx_obj->threshold);
    EXPECT_NEAR(0.0, sgeesx_obj->diff_wi, sgeesx_obj->threshold);
    EXPECT_NEAR(0.0, sgeesx_obj->diff_rconde, sgeesx_obj->threshold);
    EXPECT_NEAR(0.0, sgeesx_obj->diff_rcondv, sgeesx_obj->threshold);
    EXPECT_EQ(sgeesx_obj->sdim, sgeesx_obj->sdimref);
}
TEST_F(sgeesx_test, sgeesx3) {
    //EXPECT_NEAR(0.0, sgeesx_obj->diff_a, sgeesx_obj->threshold);
    EXPECT_NEAR(0.0, sgeesx_obj->diff_vs, sgeesx_obj->threshold);
    EXPECT_NEAR(0.0, sgeesx_obj->diff_wr, sgeesx_obj->threshold);
    EXPECT_NEAR(0.0, sgeesx_obj->diff_wi, sgeesx_obj->threshold);
    EXPECT_NEAR(0.0, sgeesx_obj->diff_rconde, sgeesx_obj->threshold);
    EXPECT_NEAR(0.0, sgeesx_obj->diff_rcondv, sgeesx_obj->threshold);
    EXPECT_EQ(sgeesx_obj->sdim, sgeesx_obj->sdimref);
}
TEST_F(sgeesx_test, sgeesx4) {
    //EXPECT_NEAR(0.0, sgeesx_obj->diff_a, sgeesx_obj->threshold);
    EXPECT_NEAR(0.0, sgeesx_obj->diff_vs, sgeesx_obj->threshold);
    EXPECT_NEAR(0.0, sgeesx_obj->diff_wr, sgeesx_obj->threshold);
    EXPECT_NEAR(0.0, sgeesx_obj->diff_wi, sgeesx_obj->threshold);
    EXPECT_NEAR(0.0, sgeesx_obj->diff_rconde, sgeesx_obj->threshold);
    EXPECT_NEAR(0.0, sgeesx_obj->diff_rcondv, sgeesx_obj->threshold);
    EXPECT_EQ(sgeesx_obj->sdim, sgeesx_obj->sdimref);
}
#endif