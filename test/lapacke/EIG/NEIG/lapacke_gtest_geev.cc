#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"
#define LAPACKE_TEST_VERBOSE  (1)

#define geev_free() \
    if (a!=NULL)        free(a); \
    if (aref!=NULL)     free(aref); \
    if (vl!=NULL)       free(vl); \
    if (vlref!=NULL)    free(vlref); \
    if (vr!=NULL)       free(vr); \
    if (vrref!=NULL)    free(vrref); \
    if (wr!=NULL)     free(wr); \
    if (wrref!=NULL)     free(wrref); \
    if (wi!=NULL)     free(wi); \
    if (wiref!=NULL)     free(wiref)


// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin geev_float_common_parameters  class definition */
class geev_float_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_a;
    float diff_wi, diff_wr, diff_vl, diff_vr;

    float threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR	
    char jobvl; // Must be 'N', or 'V'.
    char jobvr; // Must be 'N', or 'V'.

    lapack_int n; // The order of the matrices A and B (n≥ 0).
    lapack_int lda; //  The leading dimension of a
    lapack_int ldvl;
    lapack_int ldvr;

    /* Input / Output parameters */
    float* a, *aref; // contains the n-by-n general matrix A.
    float* b, *bref; // contains the n-by-n upper triangular matrix B.

    /* Output parameters */
    float *vl, *vlref; // left eigenvectors.
    float *vr, *vrref; // right eigenvectors

    float* wr, *wrref;
	float* wi, *wiref;
    /*Return Values */
    int info, inforef;

      geev_float_parameters (int matrix_layout_i, char jobvl_i,
						char jobvr_i,  lapack_int n_i );
      ~geev_float_parameters ();
};

/* Constructor definition  float_common_parameters */
geev_float_parameters:: geev_float_parameters (int matrix_layout_i,
            char jobvl_i, char jobvr_i,  lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    jobvl = jobvl_i;
    jobvr = jobvr_i;
	
    n  = n_i;
    lda = n;
    ldvl = n;
    ldvr = n;

    hModule = NULL;
    dModule = NULL;

    diff_a = 0;
    diff_vl = 0;
    diff_vr = 0;
    diff_wi = 0;
    diff_wr = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n geev float: matrix_layout: %d n: %d  jobvl: %c \
jobvr: %c \n", matrix_layout, n, jobvr, jobvl);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, lda*n );
    lapacke_gtest_alloc_float_buffer_pair( &wr, &wrref, n );
    lapacke_gtest_alloc_float_buffer_pair( &wi, &wiref, n );
    lapacke_gtest_alloc_float_buffer_pair( &vl, &vlref, n*ldvl );
    lapacke_gtest_alloc_float_buffer_pair( &vr, &vrref, n*ldvr );


    if( (a==NULL) || (aref==NULL) ||  \
        (wr==NULL) || (wrref==NULL) || \
        (wi==NULL) || (wiref==NULL) || \
        (vl==NULL) || (vlref==NULL) || \
        (vr==NULL) || (vrref==NULL) ){
       EXPECT_FALSE( true) << "geev_float_parameters object: malloc error. Exiting ";
       geev_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, lda*n );
    lapacke_gtest_init_float_buffer_pair_with_constant( vl, vlref, ldvl*n, 0 );
    lapacke_gtest_init_float_buffer_pair_with_constant( vr, vrref, ldvr*n, 0 );
    lapacke_gtest_init_float_buffer_pair_with_constant( wr, wrref, n, 0);
    lapacke_gtest_init_float_buffer_pair_with_constant( wi, wiref, n, 0);

   } /* end of Constructor  */
    

/* Destructor definition  'geev_float_common_parameters' */
geev_float_parameters :: ~geev_float_parameters ()
{
   /* De-Allocate memory for the input matrices */
   geev_free();
} 

//  Test fixture class definition
class sgeev_test  : public  ::testing::Test {
public:
   geev_float_parameters  *sgeev_obj;
   void SetUp();  
   void TearDown () { delete sgeev_obj; }
};

void sgeev_test::SetUp()
{

    /* LAPACKE_sgeev prototype */	
    typedef int (*Fptr_NL_LAPACKE_sgeev) (int matrix_layout, char jobvl,
		char jobvr, lapack_int n, float* a, lapack_int lda, float* wr,
		float* wi, float* vl, lapack_int ldvl, float* vr, lapack_int ldvr);
				 
    Fptr_NL_LAPACKE_sgeev SGEEV;

    sgeev_obj = new  geev_float_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].jobvsl,
										 eig_non_sym_paramslist[idx].jobvsr,
                                         eig_paramslist[idx].n );
                                         
    sgeev_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    sgeev_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgeev_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgeev_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgeev_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SGEEV = (Fptr_NL_LAPACKE_sgeev)dlsym(sgeev_obj->hModule, "LAPACKE_sgeev");
    ASSERT_TRUE(SGEEV != NULL) << "failed to ppt the Netlib LAPACKE_sgeev symbol";

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    sgeev_obj->inforef = SGEEV(   sgeev_obj->matrix_layout,
                                    sgeev_obj->jobvl,
                                    sgeev_obj->jobvr,
                                    sgeev_obj->n,
                                    sgeev_obj->aref,
                                    sgeev_obj->lda,
                                    sgeev_obj->wrref,
                                    sgeev_obj->wiref,
									sgeev_obj->vlref,
									sgeev_obj->ldvl,
									sgeev_obj->vrref,
									sgeev_obj->ldvr
                                    );

    /* Compute libflame's Lapacke o/p  */
    sgeev_obj->inforef =  LAPACKE_sgeev (   sgeev_obj->matrix_layout,
                                    sgeev_obj->jobvl,
                                    sgeev_obj->jobvr,
                                    sgeev_obj->n,
                                    sgeev_obj->a,
                                    sgeev_obj->lda,
                                    sgeev_obj->wr,
                                    sgeev_obj->wi,
									sgeev_obj->vl,
									sgeev_obj->ldvl,
									sgeev_obj->vr,
									sgeev_obj->ldvr
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    
    sgeev_obj->diff_a =  computeDiff_s( (sgeev_obj->lda)*(sgeev_obj->n), 
                sgeev_obj->a, sgeev_obj->aref );

    sgeev_obj->diff_vl =  computeDiff_s( (sgeev_obj->ldvl)*(sgeev_obj->n), 
                sgeev_obj->vl, sgeev_obj->vlref );

    sgeev_obj->diff_vr =  computeDiff_s( (sgeev_obj->ldvr)*(sgeev_obj->n), 
                sgeev_obj->vr, sgeev_obj->vrref );

    sgeev_obj->diff_wr =  computeDiff_s( sgeev_obj->n, 
                sgeev_obj->wr, sgeev_obj->wrref );

    sgeev_obj->diff_wi =  computeDiff_s( sgeev_obj->n, 
                sgeev_obj->wi, sgeev_obj->wiref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n geev float: \n diff_a: %f \n  \
diff_vl: %f \n diff_vr: %f \n diff_wr: %f \n\
diff_wi: %f \n  ",
       sgeev_obj->diff_a, sgeev_obj->diff_vl,
	   sgeev_obj->diff_vr, sgeev_obj->diff_wr,
       sgeev_obj->diff_wi );
#endif
}

TEST_F(sgeev_test, sgeev1) {
    //EXPECT_NEAR(0.0, sgeev_obj->diff_a, sgeev_obj->threshold);
    EXPECT_NEAR(0.0, sgeev_obj->diff_vl, sgeev_obj->threshold);
    EXPECT_NEAR(0.0, sgeev_obj->diff_vr, sgeev_obj->threshold);
    EXPECT_NEAR(0.0, sgeev_obj->diff_wr, sgeev_obj->threshold);
    EXPECT_NEAR(0.0, sgeev_obj->diff_wi, sgeev_obj->threshold);
}

TEST_F(sgeev_test, sgeev2) {
    //EXPECT_NEAR(0.0, sgeev_obj->diff_a, sgeev_obj->threshold);
    EXPECT_NEAR(0.0, sgeev_obj->diff_vl, sgeev_obj->threshold);
    EXPECT_NEAR(0.0, sgeev_obj->diff_vr, sgeev_obj->threshold);
    EXPECT_NEAR(0.0, sgeev_obj->diff_wr, sgeev_obj->threshold);
    EXPECT_NEAR(0.0, sgeev_obj->diff_wi, sgeev_obj->threshold);
}
TEST_F(sgeev_test, sgeev3) {
    //EXPECT_NEAR(0.0, sgeev_obj->diff_a, sgeev_obj->threshold);
    EXPECT_NEAR(0.0, sgeev_obj->diff_vl, sgeev_obj->threshold);
    EXPECT_NEAR(0.0, sgeev_obj->diff_vr, sgeev_obj->threshold);
    EXPECT_NEAR(0.0, sgeev_obj->diff_wr, sgeev_obj->threshold);
    EXPECT_NEAR(0.0, sgeev_obj->diff_wi, sgeev_obj->threshold);
}
TEST_F(sgeev_test, sgeev4) {
    //EXPECT_NEAR(0.0, sgeev_obj->diff_a, sgeev_obj->threshold);
    EXPECT_NEAR(0.0, sgeev_obj->diff_vl, sgeev_obj->threshold);
    EXPECT_NEAR(0.0, sgeev_obj->diff_vr, sgeev_obj->threshold);
    EXPECT_NEAR(0.0, sgeev_obj->diff_wr, sgeev_obj->threshold);
    EXPECT_NEAR(0.0, sgeev_obj->diff_wi, sgeev_obj->threshold);
}

/* Begin geev_double_common_parameters  class definition */
class geev_double_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    double diff_a;
    double diff_wi, diff_wr, diff_vl, diff_vr;

    float threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR	
    char jobvl; // Must be 'N', or 'V'.
    char jobvr; // Must be 'N', or 'V'.

    lapack_int n; // The order of the matrices A and B (n≥ 0).
    lapack_int lda; //  The leading dimension of a
    lapack_int ldvl;
    lapack_int ldvr;

    /* Input / Output parameters */
    double* a, *aref; // contains the n-by-n general matrix A.
    double* b, *bref; // contains the n-by-n upper triangular matrix B.

    /* Output parameters */
    double *vl, *vlref; // left eigenvectors.
    double *vr, *vrref; // right eigenvectors

    double* wr, *wrref;
	double* wi, *wiref;
    /*Return Values */
    int info, inforef;

      geev_double_parameters (int matrix_layout_i, char jobvl_i,
						char jobvr_i,  lapack_int n_i );
      ~geev_double_parameters ();
};

/* Constructor definition  double_common_parameters */
geev_double_parameters:: geev_double_parameters (int matrix_layout_i,
            char jobvl_i, char jobvr_i,  lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    jobvl = jobvl_i;
    jobvr = jobvr_i;
	
    n  = n_i;
    lda = n;
    ldvl = n;
    ldvr = n;

    hModule = NULL;
    dModule = NULL;

    diff_a = 0;
    diff_vl = 0;
    diff_vr = 0;
    diff_wi = 0;
    diff_wr = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n geev double: matrix_layout: %d n: %d  jobvl: %c \
jobvr: %c \n", matrix_layout, n, jobvr, jobvl);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, lda*n );
    lapacke_gtest_alloc_double_buffer_pair( &wr, &wrref, n );
    lapacke_gtest_alloc_double_buffer_pair( &wi, &wiref, n );
    lapacke_gtest_alloc_double_buffer_pair( &vl, &vlref, n*ldvl );
    lapacke_gtest_alloc_double_buffer_pair( &vr, &vrref, n*ldvr );


    if( (a==NULL) || (aref==NULL) ||  \
        (wr==NULL) || (wrref==NULL) || \
        (wi==NULL) || (wiref==NULL) || \
        (vl==NULL) || (vlref==NULL) || \
        (vr==NULL) || (vrref==NULL) ){
       EXPECT_FALSE( true) << "geev_double_parameters object: malloc error. Exiting ";
       geev_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, lda*n );
    lapacke_gtest_init_double_buffer_pair_with_constant( vl, vlref, ldvl*n, 0 );
    lapacke_gtest_init_double_buffer_pair_with_constant( vr, vrref, ldvr*n, 0 );
    lapacke_gtest_init_double_buffer_pair_with_constant( wr, wrref, n, 0);
    lapacke_gtest_init_double_buffer_pair_with_constant( wi, wiref, n, 0);

   } /* end of Constructor  */
    

/* Destructor definition  'geev_double_common_parameters' */
geev_double_parameters :: ~geev_double_parameters ()
{
   /* De-Allocate memory for the input matrices */
   geev_free();
} 

//  Test fixture class definition
class dgeev_test  : public  ::testing::Test {
public:
   geev_double_parameters  *dgeev_obj;
   void SetUp();  
   void TearDown () { delete dgeev_obj; }
};

void dgeev_test::SetUp()
{

    /* LAPACKE_dgeev prototype */	
    typedef int (*Fptr_NL_LAPACKE_dgeev) (int matrix_layout, char jobvl,
		char jobvr, lapack_int n, double* a, lapack_int lda, double* wr,
		double* wi, double* vl, lapack_int ldvl, double* vr, lapack_int ldvr);
				 
    Fptr_NL_LAPACKE_dgeev DGEEV;

    dgeev_obj = new  geev_double_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].jobvsl,
										 eig_non_sym_paramslist[idx].jobvsr,
                                         eig_paramslist[idx].n );
                                         
    dgeev_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    dgeev_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgeev_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgeev_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgeev_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DGEEV = (Fptr_NL_LAPACKE_dgeev)dlsym(dgeev_obj->hModule, "LAPACKE_dgeev");
    ASSERT_TRUE(DGEEV != NULL) << "failed to ppt the Netlib LAPACKE_dgeev symbol";

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    dgeev_obj->inforef = DGEEV(   dgeev_obj->matrix_layout,
                                    dgeev_obj->jobvl,
                                    dgeev_obj->jobvr,
                                    dgeev_obj->n,
                                    dgeev_obj->aref,
                                    dgeev_obj->lda,
                                    dgeev_obj->wrref,
                                    dgeev_obj->wiref,
									dgeev_obj->vlref,
									dgeev_obj->ldvl,
									dgeev_obj->vrref,
									dgeev_obj->ldvr
                                    );

    /* Compute libflame's Lapacke o/p  */
    dgeev_obj->inforef =  LAPACKE_dgeev (   dgeev_obj->matrix_layout,
                                    dgeev_obj->jobvl,
                                    dgeev_obj->jobvr,
                                    dgeev_obj->n,
                                    dgeev_obj->a,
                                    dgeev_obj->lda,
                                    dgeev_obj->wr,
                                    dgeev_obj->wi,
									dgeev_obj->vl,
									dgeev_obj->ldvl,
									dgeev_obj->vr,
									dgeev_obj->ldvr
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    
    dgeev_obj->diff_a =  computeDiff_d( (dgeev_obj->lda)*(dgeev_obj->n), 
                dgeev_obj->a, dgeev_obj->aref );

    dgeev_obj->diff_vl =  computeDiff_d( (dgeev_obj->ldvl)*(dgeev_obj->n), 
                dgeev_obj->vl, dgeev_obj->vlref );

    dgeev_obj->diff_vr =  computeDiff_d( (dgeev_obj->ldvr)*(dgeev_obj->n), 
                dgeev_obj->vr, dgeev_obj->vrref );

    dgeev_obj->diff_wr =  computeDiff_d( dgeev_obj->n, 
                dgeev_obj->wr, dgeev_obj->wrref );

    dgeev_obj->diff_wi =  computeDiff_d( dgeev_obj->n, 
                dgeev_obj->wi, dgeev_obj->wiref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n geev double: \n diff_a: %f \n  \
diff_vl: %f \n diff_vr: %f \n diff_wr: %f \n\
diff_wi: %f \n  ",
       dgeev_obj->diff_a, dgeev_obj->diff_vl,
	   dgeev_obj->diff_vr, dgeev_obj->diff_wr,
       dgeev_obj->diff_wi );
#endif
}

TEST_F(dgeev_test, dgeev1) {
    //EXPECT_NEAR(0.0, dgeev_obj->diff_a, dgeev_obj->threshold);
    EXPECT_NEAR(0.0, dgeev_obj->diff_vl, dgeev_obj->threshold);
    EXPECT_NEAR(0.0, dgeev_obj->diff_vr, dgeev_obj->threshold);
    EXPECT_NEAR(0.0, dgeev_obj->diff_wr, dgeev_obj->threshold);
    EXPECT_NEAR(0.0, dgeev_obj->diff_wi, dgeev_obj->threshold);
}

TEST_F(dgeev_test, dgeev2) {
    //EXPECT_NEAR(0.0, dgeev_obj->diff_a, dgeev_obj->threshold);
    EXPECT_NEAR(0.0, dgeev_obj->diff_vl, dgeev_obj->threshold);
    EXPECT_NEAR(0.0, dgeev_obj->diff_vr, dgeev_obj->threshold);
    EXPECT_NEAR(0.0, dgeev_obj->diff_wr, dgeev_obj->threshold);
    EXPECT_NEAR(0.0, dgeev_obj->diff_wi, dgeev_obj->threshold);
}
TEST_F(dgeev_test, dgeev3) {
    //EXPECT_NEAR(0.0, dgeev_obj->diff_a, dgeev_obj->threshold);
    EXPECT_NEAR(0.0, dgeev_obj->diff_vl, dgeev_obj->threshold);
    EXPECT_NEAR(0.0, dgeev_obj->diff_vr, dgeev_obj->threshold);
    EXPECT_NEAR(0.0, dgeev_obj->diff_wr, dgeev_obj->threshold);
    EXPECT_NEAR(0.0, dgeev_obj->diff_wi, dgeev_obj->threshold);
}
TEST_F(dgeev_test, dgeev4) {
    //EXPECT_NEAR(0.0, dgeev_obj->diff_a, dgeev_obj->threshold);
    EXPECT_NEAR(0.0, dgeev_obj->diff_vl, dgeev_obj->threshold);
    EXPECT_NEAR(0.0, dgeev_obj->diff_vr, dgeev_obj->threshold);
    EXPECT_NEAR(0.0, dgeev_obj->diff_wr, dgeev_obj->threshold);
    EXPECT_NEAR(0.0, dgeev_obj->diff_wi, dgeev_obj->threshold);
}

/* Begin geev_scomplex_common_parameters  class definition */
class geev_scomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_a;
    float diff_wi, diff_wr, diff_vl, diff_vr;

    float threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR	
    char jobvl; // Must be 'N', or 'V'.
    char jobvr; // Must be 'N', or 'V'.

    lapack_int n; // The order of the matrices A and B (n≥ 0).
    lapack_int lda; //  The leading dimension of a
    lapack_int ldvl;
    lapack_int ldvr;

    /* Input / Output parameters */
    lapack_complex_float* a, *aref; // contains the n-by-n general matrix A.
    lapack_complex_float* b, *bref; // contains the n-by-n upper triangular matrix B.

    /* Output parameters */
    lapack_complex_float *vl, *vlref; // left eigenvectors.
    lapack_complex_float *vr, *vrref; // right eigenvectors

    lapack_complex_float* wr, *wrref;
	lapack_complex_float* wi, *wiref;
    /*Return Values */
    int info, inforef;

      geev_scomplex_parameters (int matrix_layout_i, char jobvl_i,
						char jobvr_i,  lapack_int n_i );
      ~geev_scomplex_parameters ();
};

/* Constructor definition  lapack_complex_float_common_parameters */
geev_scomplex_parameters:: geev_scomplex_parameters (int matrix_layout_i,
            char jobvl_i, char jobvr_i,  lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    jobvl = jobvl_i;
    jobvr = jobvr_i;
	
    n  = n_i;
    lda = n;
    ldvl = n;
    ldvr = n;

    hModule = NULL;
    dModule = NULL;

    diff_a = 0;
    diff_vl = 0;
    diff_vr = 0;
    diff_wi = 0;
    diff_wr = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n geev lapack_complex_float: matrix_layout: %d n: %d  jobvl: %c \
jobvr: %c \n", matrix_layout, n, jobvr, jobvl);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, lda*n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &wr, &wrref, n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &wi, &wiref, n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &vl, &vlref, n*ldvl );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &vr, &vrref, n*ldvr );


    if( (a==NULL) || (aref==NULL) ||  \
        (wr==NULL) || (wrref==NULL) || \
        (wi==NULL) || (wiref==NULL) || \
        (vl==NULL) || (vlref==NULL) || \
        (vr==NULL) || (vrref==NULL) ){
       EXPECT_FALSE( true) << "geev_scomplex_parameters object: malloc error. Exiting ";
       geev_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, lda*n );
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( vl, vlref, ldvl*n, 0 );
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( vr, vrref, ldvr*n, 0 );
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( wr, wrref, n, 0);
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( wi, wiref, n, 0);

   } /* end of Constructor  */
    

/* Destructor definition  'geev_scomplex_common_parameters' */
geev_scomplex_parameters :: ~geev_scomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   geev_free();
} 

//  Test fixture class definition
class cgeev_test  : public  ::testing::Test {
public:
   geev_scomplex_parameters  *cgeev_obj;
   void SetUp();  
   void TearDown () { delete cgeev_obj; }
};

void cgeev_test::SetUp()
{

    /* LAPACKE_cgeev prototype */	
    typedef int (*Fptr_NL_LAPACKE_cgeev) (int matrix_layout, char jobvl,
		char jobvr, lapack_int n, lapack_complex_float* a, lapack_int lda,
		lapack_complex_float* w, lapack_complex_float* vl, lapack_int ldvl,
		lapack_complex_float* vr, lapack_int ldvr);
				 
    Fptr_NL_LAPACKE_cgeev CGEEV;

    cgeev_obj = new  geev_scomplex_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].jobvsl,
										 eig_non_sym_paramslist[idx].jobvsr,
                                         eig_paramslist[idx].n );
                                         
    cgeev_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    cgeev_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgeev_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgeev_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgeev_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CGEEV = (Fptr_NL_LAPACKE_cgeev)dlsym(cgeev_obj->hModule, "LAPACKE_cgeev");
    ASSERT_TRUE(CGEEV != NULL) << "failed to ppt the Netlib LAPACKE_cgeev symbol";

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    cgeev_obj->inforef = CGEEV(   cgeev_obj->matrix_layout,
                                    cgeev_obj->jobvl,
                                    cgeev_obj->jobvr,
                                    cgeev_obj->n,
                                    cgeev_obj->aref,
                                    cgeev_obj->lda,
                                    cgeev_obj->wrref,
									cgeev_obj->vlref,
									cgeev_obj->ldvl,
									cgeev_obj->vrref,
									cgeev_obj->ldvr
                                    );

    /* Compute libflame's Lapacke o/p  */
    cgeev_obj->inforef =  LAPACKE_cgeev (   cgeev_obj->matrix_layout,
                                    cgeev_obj->jobvl,
                                    cgeev_obj->jobvr,
                                    cgeev_obj->n,
                                    cgeev_obj->a,
                                    cgeev_obj->lda,
                                    cgeev_obj->wr,
									cgeev_obj->vl,
									cgeev_obj->ldvl,
									cgeev_obj->vr,
									cgeev_obj->ldvr
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    
    cgeev_obj->diff_a =  computeDiff_c( (cgeev_obj->lda)*(cgeev_obj->n), 
                cgeev_obj->a, cgeev_obj->aref );

    cgeev_obj->diff_vl =  computeDiff_c( (cgeev_obj->ldvl)*(cgeev_obj->n), 
                cgeev_obj->vl, cgeev_obj->vlref );

    cgeev_obj->diff_vr =  computeDiff_c( (cgeev_obj->ldvr)*(cgeev_obj->n), 
                cgeev_obj->vr, cgeev_obj->vrref );

    cgeev_obj->diff_wr =  computeDiff_c( cgeev_obj->n, 
                cgeev_obj->wr, cgeev_obj->wrref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n geev lapack_complex_float: \n diff_a: %f \n  \
diff_vl: %f \n diff_vr: %f \n diff_wr: %f \n\
diff_wi: %f \n  ",
       cgeev_obj->diff_a, cgeev_obj->diff_vl,
	   cgeev_obj->diff_vr, cgeev_obj->diff_wr,
       cgeev_obj->diff_wi );
#endif
}

TEST_F(cgeev_test, cgeev1) {
    //EXPECT_NEAR(0.0, cgeev_obj->diff_a, cgeev_obj->threshold);
    EXPECT_NEAR(0.0, cgeev_obj->diff_vl, cgeev_obj->threshold);
    EXPECT_NEAR(0.0, cgeev_obj->diff_vr, cgeev_obj->threshold);
    EXPECT_NEAR(0.0, cgeev_obj->diff_wr, cgeev_obj->threshold);
}

TEST_F(cgeev_test, cgeev2) {
    //EXPECT_NEAR(0.0, cgeev_obj->diff_a, cgeev_obj->threshold);
    EXPECT_NEAR(0.0, cgeev_obj->diff_vl, cgeev_obj->threshold);
    EXPECT_NEAR(0.0, cgeev_obj->diff_vr, cgeev_obj->threshold);
    EXPECT_NEAR(0.0, cgeev_obj->diff_wr, cgeev_obj->threshold);
}
TEST_F(cgeev_test, cgeev3) {
    //EXPECT_NEAR(0.0, cgeev_obj->diff_a, cgeev_obj->threshold);
    EXPECT_NEAR(0.0, cgeev_obj->diff_vl, cgeev_obj->threshold);
    EXPECT_NEAR(0.0, cgeev_obj->diff_vr, cgeev_obj->threshold);
    EXPECT_NEAR(0.0, cgeev_obj->diff_wr, cgeev_obj->threshold);
}
TEST_F(cgeev_test, cgeev4) {
    //EXPECT_NEAR(0.0, cgeev_obj->diff_a, cgeev_obj->threshold);
    EXPECT_NEAR(0.0, cgeev_obj->diff_vl, cgeev_obj->threshold);
    EXPECT_NEAR(0.0, cgeev_obj->diff_vr, cgeev_obj->threshold);
    EXPECT_NEAR(0.0, cgeev_obj->diff_wr, cgeev_obj->threshold);
}

/* Begin geev_dcomplex_common_parameters  class definition */
class geev_dcomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    double diff_a;
    double diff_wi, diff_wr, diff_vl, diff_vr;

    double threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR	
    char jobvl; // Must be 'N', or 'V'.
    char jobvr; // Must be 'N', or 'V'.

    lapack_int n; // The order of the matrices A and B (n≥ 0).
    lapack_int lda; //  The leading dimension of a
    lapack_int ldvl;
    lapack_int ldvr;

    /* Input / Output parameters */
    lapack_complex_double* a, *aref; // contains the n-by-n general matrix A.
    lapack_complex_double* b, *bref; // contains the n-by-n upper triangular matrix B.

    /* Output parameters */
    lapack_complex_double *vl, *vlref; // left eigenvectors.
    lapack_complex_double *vr, *vrref; // right eigenvectors

    lapack_complex_double* wr, *wrref;
	lapack_complex_double* wi, *wiref;
    /*Return Values */
    int info, inforef;

      geev_dcomplex_parameters (int matrix_layout_i, char jobvl_i,
						char jobvr_i,  lapack_int n_i );
      ~geev_dcomplex_parameters ();
};

/* Constructor definition  lapack_complex_double_common_parameters */
geev_dcomplex_parameters:: geev_dcomplex_parameters (int matrix_layout_i,
            char jobvl_i, char jobvr_i,  lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    jobvl = jobvl_i;
    jobvr = jobvr_i;
	
    n  = n_i;
    lda = n;
    ldvl = n;
    ldvr = n;

    hModule = NULL;
    dModule = NULL;

    diff_a = 0;
    diff_vl = 0;
    diff_vr = 0;
    diff_wi = 0;
    diff_wr = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n geev lapack_complex_double: matrix_layout: %d n: %d  jobvl: %c \
jobvr: %c \n", matrix_layout, n, jobvr, jobvl);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, lda*n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &wr, &wrref, n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &wi, &wiref, n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &vl, &vlref, n*ldvl );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &vr, &vrref, n*ldvr );


    if( (a==NULL) || (aref==NULL) ||  \
        (wr==NULL) || (wrref==NULL) || \
        (wi==NULL) || (wiref==NULL) || \
        (vl==NULL) || (vlref==NULL) || \
        (vr==NULL) || (vrref==NULL) ){
       EXPECT_FALSE( true) << "geev_dcomplex_parameters object: malloc error. Exiting ";
       geev_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, lda*n );
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( vl, vlref, ldvl*n, 0 );
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( vr, vrref, ldvr*n, 0 );
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( wr, wrref, n, 0);
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( wi, wiref, n, 0);

   } /* end of Constructor  */
    

/* Destructor definition  'geev_dcomplex_common_parameters' */
geev_dcomplex_parameters :: ~geev_dcomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   geev_free();
} 

//  Test fixture class definition
class zgeev_test  : public  ::testing::Test {
public:
   geev_dcomplex_parameters  *zgeev_obj;
   void SetUp();  
   void TearDown () { delete zgeev_obj; }
};

void zgeev_test::SetUp()
{

    /* LAPACKE_zgeev prototype */	
    typedef int (*Fptr_NL_LAPACKE_zgeev) (int matrix_layout, char jobvl,
		char jobvr, lapack_int n, lapack_complex_double* a, lapack_int lda,
		lapack_complex_double* w, lapack_complex_double* vl, lapack_int ldvl,
		lapack_complex_double* vr, lapack_int ldvr);
				 
    Fptr_NL_LAPACKE_zgeev ZGEEV;

    zgeev_obj = new  geev_dcomplex_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].jobvsl,
										 eig_non_sym_paramslist[idx].jobvsr,
                                         eig_paramslist[idx].n );
                                         
    zgeev_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    zgeev_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgeev_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgeev_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgeev_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZGEEV = (Fptr_NL_LAPACKE_zgeev)dlsym(zgeev_obj->hModule, "LAPACKE_zgeev");
    ASSERT_TRUE(ZGEEV != NULL) << "failed to ppt the Netlib LAPACKE_zgeev symbol";

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    zgeev_obj->inforef = ZGEEV(   zgeev_obj->matrix_layout,
                                    zgeev_obj->jobvl,
                                    zgeev_obj->jobvr,
                                    zgeev_obj->n,
                                    zgeev_obj->aref,
                                    zgeev_obj->lda,
                                    zgeev_obj->wrref,
									zgeev_obj->vlref,
									zgeev_obj->ldvl,
									zgeev_obj->vrref,
									zgeev_obj->ldvr
                                    );

    /* Compute libflame's Lapacke o/p  */
    zgeev_obj->inforef =  LAPACKE_zgeev (   zgeev_obj->matrix_layout,
                                    zgeev_obj->jobvl,
                                    zgeev_obj->jobvr,
                                    zgeev_obj->n,
                                    zgeev_obj->a,
                                    zgeev_obj->lda,
                                    zgeev_obj->wr,
									zgeev_obj->vl,
									zgeev_obj->ldvl,
									zgeev_obj->vr,
									zgeev_obj->ldvr
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    
    zgeev_obj->diff_a =  computeDiff_z( (zgeev_obj->lda)*(zgeev_obj->n), 
                zgeev_obj->a, zgeev_obj->aref );

    zgeev_obj->diff_vl =  computeDiff_z( (zgeev_obj->ldvl)*(zgeev_obj->n), 
                zgeev_obj->vl, zgeev_obj->vlref );

    zgeev_obj->diff_vr =  computeDiff_z( (zgeev_obj->ldvr)*(zgeev_obj->n), 
                zgeev_obj->vr, zgeev_obj->vrref );

    zgeev_obj->diff_wr =  computeDiff_z( zgeev_obj->n, 
                zgeev_obj->wr, zgeev_obj->wrref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n geev lapack_complex_double: \n diff_a: %f \n  \
diff_vl: %f \n diff_vr: %f \n diff_wr: %f \n\
diff_wi: %f \n  ",
       zgeev_obj->diff_a, zgeev_obj->diff_vl,
	   zgeev_obj->diff_vr, zgeev_obj->diff_wr,
       zgeev_obj->diff_wi );
#endif
}

TEST_F(zgeev_test, zgeev1) {
    //EXPECT_NEAR(0.0, zgeev_obj->diff_a, zgeev_obj->threshold);
    EXPECT_NEAR(0.0, zgeev_obj->diff_vl, zgeev_obj->threshold);
    EXPECT_NEAR(0.0, zgeev_obj->diff_vr, zgeev_obj->threshold);
    EXPECT_NEAR(0.0, zgeev_obj->diff_wr, zgeev_obj->threshold);
}

TEST_F(zgeev_test, zgeev2) {
    //EXPECT_NEAR(0.0, zgeev_obj->diff_a, zgeev_obj->threshold);
    EXPECT_NEAR(0.0, zgeev_obj->diff_vl, zgeev_obj->threshold);
    EXPECT_NEAR(0.0, zgeev_obj->diff_vr, zgeev_obj->threshold);
    EXPECT_NEAR(0.0, zgeev_obj->diff_wr, zgeev_obj->threshold);
}
TEST_F(zgeev_test, zgeev3) {
    //EXPECT_NEAR(0.0, zgeev_obj->diff_a, zgeev_obj->threshold);
    EXPECT_NEAR(0.0, zgeev_obj->diff_vl, zgeev_obj->threshold);
    EXPECT_NEAR(0.0, zgeev_obj->diff_vr, zgeev_obj->threshold);
    EXPECT_NEAR(0.0, zgeev_obj->diff_wr, zgeev_obj->threshold);
}
TEST_F(zgeev_test, zgeev4) {
    //EXPECT_NEAR(0.0, zgeev_obj->diff_a, zgeev_obj->threshold);
    EXPECT_NEAR(0.0, zgeev_obj->diff_vl, zgeev_obj->threshold);
    EXPECT_NEAR(0.0, zgeev_obj->diff_vr, zgeev_obj->threshold);
    EXPECT_NEAR(0.0, zgeev_obj->diff_wr, zgeev_obj->threshold);
}
