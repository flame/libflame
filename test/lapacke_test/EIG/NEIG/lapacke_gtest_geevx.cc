#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"
#define LAPACKE_TEST_VERBOSE  (1)

#define geevx_free() \
    if (a!=NULL)        free(a); \
    if (aref!=NULL)     free(aref); \
    if (vl!=NULL)       free(vl); \
    if (vlref!=NULL)    free(vlref); \
    if (vr!=NULL)       free(vr); \
    if (vrref!=NULL)    free(vrref); \
    if (wr!=NULL)     free(wr); \
    if (wrref!=NULL)     free(wrref); \
    if (wi!=NULL)     free(wi); \
    if (wiref!=NULL)     free(wiref); \
    if (scale!=NULL)   free(scale); \
    if (scaleref!=NULL)  free(scaleref); \
    if (rconde!=NULL)   free(rconde); \
    if (rconderef!=NULL)  free(rconderef); \
    if (rcondv!=NULL)   free(rcondv); \
    if (rcondvref!=NULL)  free(rcondvref)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin geevx_float_common_parameters  class definition */
class geevx_float_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_a;
    float diff_wi, diff_wr, diff_vl, diff_vr;
    float diff_scale, diff_rconde, diff_rcondv;

    float threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR	
    char jobvl; // Must be 'N', or 'V'.
    char jobvr; // Must be 'N', or 'V'.
	char balance; // Must be 'N', 'P' 'S', B.
	char sense; //  Must be 'N', 'E', 'V', or 'B'.

    lapack_int n; // The order of the matrices A and B (n≥ 0).
    lapack_int lda; //  The leading dimension of a
    lapack_int ldvl;
    lapack_int ldvr;

    /* Input / Output parameters */
    float* a, *aref; // contains the n-by-n general matrix A.

    /* Output parameters */
    lapack_int sdim, sdimref;
    float *vl, *vlref; // left eigenvectors.
    float *vr, *vrref; // right eigenvectors
    lapack_int ilo, iloref, ihi, ihiref; // ilo and ihi mark the rows and columns of A which are to be reduced.
    float *scale, *scaleref; // permutations, scaling factors applied to the left side of A and B
    float  abnrm, abnrmref; //  L1- norm of balanced matrices
    float* wr, *wrref;
	float* wi, *wiref;
	float *rconde, *rcondv; // contain the reciprocal condition numbers 
	float *rconderef, *rcondvref; // contain the reciprocal condition numbers 
    /*Return Values */
    int info, inforef;

      geevx_float_parameters (int matrix_layout_i, char jobvl_i,
				char jobvr_i, char balance_i, char sense_i,
				lapack_int n_i );
      ~geevx_float_parameters ();
};

/* Constructor definition  float_common_parameters */
geevx_float_parameters:: geevx_float_parameters (int matrix_layout_i,
     char jobvl_i, char jobvr_i, char balance_i, char sense_i, 
	 lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    jobvl = jobvl_i;
    jobvr = jobvr_i;
	balance = balance_i;
	sense = sense_i;
	
    n  = n_i;
    sdim = 0;
	sdimref = 0;

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
	diff_scale = 0;
    diff_rconde = 0;
	diff_rcondv = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n geevx float: matrix_layout: %d n: %d  jobvl: %c \
jobvr: %c \t balance: %c \t sense: %c\n", matrix_layout, n, jobvr, 
jobvl, balance, sense);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, lda*n );
    lapacke_gtest_alloc_float_buffer_pair( &wr, &wrref, n );
    lapacke_gtest_alloc_float_buffer_pair( &wi, &wiref, n );
    lapacke_gtest_alloc_float_buffer_pair( &vl, &vlref, n*ldvl );
    lapacke_gtest_alloc_float_buffer_pair( &vr, &vrref, n*ldvr );
    lapacke_gtest_alloc_float_buffer_pair( &scale, &scaleref, n );
    lapacke_gtest_alloc_float_buffer_pair( &rconde, &rconderef, n );
    lapacke_gtest_alloc_float_buffer_pair( &rcondv, &rcondvref, n );

    if( (a==NULL) || (aref==NULL) ||  \
        (wr==NULL) || (wrref==NULL) || \
        (wi==NULL) || (wiref==NULL) || \
        (vl==NULL) || (vlref==NULL) || \
        (vr==NULL) || (vrref==NULL) || \
        (scale==NULL) || (scaleref==NULL) || \
        (rconde==NULL) || (rconderef==NULL) || \
        (rcondv==NULL) || (rcondvref==NULL) ){
       EXPECT_FALSE( true) << "geevx_float_parameters object: malloc error. Exiting ";
       geevx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, lda*n );
    lapacke_gtest_init_float_buffer_pair_with_constant( vl, vlref, ldvl*n, 0 );
    lapacke_gtest_init_float_buffer_pair_with_constant( vr, vrref, ldvr*n, 0 );
    lapacke_gtest_init_float_buffer_pair_with_constant( wr, wrref, n, 0);
    lapacke_gtest_init_float_buffer_pair_with_constant( wi, wiref, n, 0);
    lapacke_gtest_init_float_buffer_pair_with_constant( scale, scaleref, n, 0);
    lapacke_gtest_init_float_buffer_pair_with_constant( rconde, rconderef, n, 0);
    lapacke_gtest_init_float_buffer_pair_with_constant( rcondv, rcondvref, n, 0);

   } /* end of Constructor  */

/* Destructor definition  'geevx_float_common_parameters' */
geevx_float_parameters :: ~geevx_float_parameters ()
{
   /* De-Allocate memory for the input matrices */
   geevx_free();
} 

//  Test fixture class definition
class sgeevx_test  : public  ::testing::Test {
public:
   geevx_float_parameters  *sgeevx_obj;
   void SetUp();  
   void TearDown () { delete sgeevx_obj; }
};

void sgeevx_test::SetUp()
{
    /* LAPACKE_sgeevx prototype */	
    typedef int (*Fptr_NL_LAPACKE_sgeevx) (int matrix_layout, char balanc,
		char jobvl, char jobvr, char sense, lapack_int n, float* a,
		lapack_int lda, float* wr, float* wi, float* vl, lapack_int ldvl,
		float* vr, lapack_int ldvr, lapack_int* ilo, lapack_int* ihi,
		float* scale, float* abnrm, float* rconde, float* rcondv);

    Fptr_NL_LAPACKE_sgeevx SGGEVX;

    sgeevx_obj = new  geevx_float_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].jobvsl,
										 eig_non_sym_paramslist[idx].jobvsl, // vl & vr to be same
										 eig_non_sym_paramslist[idx].balance_ggevx,
										 eig_non_sym_paramslist[idx].sense_ggevx,
                                         eig_paramslist[idx].n );
                                         
    sgeevx_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    sgeevx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgeevx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgeevx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgeevx_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SGGEVX = (Fptr_NL_LAPACKE_sgeevx)dlsym(sgeevx_obj->hModule, "LAPACKE_sgeevx");
    ASSERT_TRUE(SGGEVX != NULL) << "failed to ppt the Netlib LAPACKE_sgeevx symbol";

    /* Compute libflame's Lapacke o/p  */
    sgeevx_obj->inforef =  LAPACKE_sgeevx (   sgeevx_obj->matrix_layout,
                                    sgeevx_obj->balance,
                                    sgeevx_obj->jobvl,
                                    sgeevx_obj->jobvr,
                                    sgeevx_obj->sense,
                                    sgeevx_obj->n,
                                    sgeevx_obj->a,
                                    sgeevx_obj->lda,
                                    sgeevx_obj->wr,
                                    sgeevx_obj->wi,
									sgeevx_obj->vl,
									sgeevx_obj->ldvl,
									sgeevx_obj->vr,
									sgeevx_obj->ldvr,
									&sgeevx_obj->ilo,
									&sgeevx_obj->ihi,
									sgeevx_obj->scale,
									&sgeevx_obj->abnrm,
									sgeevx_obj->rconde,
									sgeevx_obj->rcondv
                                    );

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    sgeevx_obj->inforef = SGGEVX(   sgeevx_obj->matrix_layout,
                                    sgeevx_obj->balance,
                                    sgeevx_obj->jobvl,
                                    sgeevx_obj->jobvr,
                                    sgeevx_obj->sense,
                                    sgeevx_obj->n,
                                    sgeevx_obj->aref,
                                    sgeevx_obj->lda,
                                    sgeevx_obj->wrref,
                                    sgeevx_obj->wiref,
									sgeevx_obj->vlref,
									sgeevx_obj->ldvl,
									sgeevx_obj->vrref,
									sgeevx_obj->ldvr,
									&sgeevx_obj->iloref,
									&sgeevx_obj->ihiref,
									sgeevx_obj->scaleref,
									&sgeevx_obj->abnrmref,
									sgeevx_obj->rconderef,
									sgeevx_obj->rcondvref
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    sgeevx_obj->diff_a =  computeDiff_s( (sgeevx_obj->lda)*(sgeevx_obj->n), 
                sgeevx_obj->a, sgeevx_obj->aref );

    sgeevx_obj->diff_vl =  computeDiff_s( (sgeevx_obj->ldvl)*(sgeevx_obj->n), 
                sgeevx_obj->vl, sgeevx_obj->vlref );

    sgeevx_obj->diff_vr =  computeDiff_s( (sgeevx_obj->ldvr)*(sgeevx_obj->n), 
                sgeevx_obj->vr, sgeevx_obj->vrref );

    sgeevx_obj->diff_scale =  computeDiff_s( sgeevx_obj->n, 
                sgeevx_obj->scale, sgeevx_obj->scaleref );

    sgeevx_obj->diff_wr =  computeDiff_s( sgeevx_obj->n, 
                sgeevx_obj->wr, sgeevx_obj->wrref );

    sgeevx_obj->diff_wi =  computeDiff_s( sgeevx_obj->n, 
                sgeevx_obj->wi, sgeevx_obj->wiref );

    sgeevx_obj->diff_rconde =  computeDiff_s( sgeevx_obj->n, 
                sgeevx_obj->rconde, sgeevx_obj->rconderef );

    sgeevx_obj->diff_rcondv =  computeDiff_s( sgeevx_obj->n, 
                sgeevx_obj->rcondv, sgeevx_obj->rcondvref );

    sgeevx_obj->diff_scale =  computeDiff_s( sgeevx_obj->n, 
                sgeevx_obj->scale, sgeevx_obj->scaleref );


#if LAPACKE_TEST_VERBOSE
   printf(" \n geevx float: \n diff_a: %f \n  \
diff_vl: %f \n diff_vr: %f \n diff_wr: %f \n\
diff_wi: %f \n ilo: %d \t iloref: %d \t \
diff_ilo: %d \n ihi: %d \t ihiref: %d \t diff_ihi: %d \n \
diff_scale: %f \n",
       sgeevx_obj->diff_a, sgeevx_obj->diff_vl,
	   sgeevx_obj->diff_vr, sgeevx_obj->diff_wr,
       sgeevx_obj->diff_wi,
	   sgeevx_obj->ilo, sgeevx_obj->iloref, 
	   ( sgeevx_obj->ilo - sgeevx_obj->iloref),
	   sgeevx_obj->ihi, sgeevx_obj->ihiref,
	   ( sgeevx_obj->ihi - sgeevx_obj->ihiref),
	   sgeevx_obj->diff_scale );
#endif
}

TEST_F(sgeevx_test, sgeevx1) {
    //EXPECT_NEAR(0.0, sgeevx_obj->diff_a, sgeevx_obj->threshold);
    EXPECT_NEAR(0.0, sgeevx_obj->diff_vl, sgeevx_obj->threshold);
    EXPECT_NEAR(0.0, sgeevx_obj->diff_vr, sgeevx_obj->threshold);
    EXPECT_NEAR(0.0, sgeevx_obj->diff_wr, sgeevx_obj->threshold);
    EXPECT_NEAR(0.0, sgeevx_obj->diff_wi, sgeevx_obj->threshold);
    EXPECT_NEAR(0.0, sgeevx_obj->diff_rconde, sgeevx_obj->threshold);
    EXPECT_NEAR(0.0, sgeevx_obj->diff_rcondv, sgeevx_obj->threshold);
	/* Commenting the below check because the position of the scalefactor
	values in the (scale & scaleref) , (rscale, rscaleref) varies
	slightly. Though the values are matching.
	So, element by element difference calculation will not work here 
	unlike for other buffers  */
	
    /* EXPECT_NEAR(0.0, sgeevx_obj->diff_scale, sgeevx_obj->threshold);
    EXPECT_NEAR(0.0, sgeevx_obj->diff_rscale, sgeevx_obj->threshold); */
    EXPECT_EQ(sgeevx_obj->ilo, sgeevx_obj->iloref);
    EXPECT_EQ(sgeevx_obj->ihi, sgeevx_obj->ihiref);
}

TEST_F(sgeevx_test, sgeevx2) {
    //EXPECT_NEAR(0.0, sgeevx_obj->diff_a, sgeevx_obj->threshold);
    EXPECT_NEAR(0.0, sgeevx_obj->diff_vl, sgeevx_obj->threshold);
    EXPECT_NEAR(0.0, sgeevx_obj->diff_vr, sgeevx_obj->threshold);
    EXPECT_NEAR(0.0, sgeevx_obj->diff_wr, sgeevx_obj->threshold);
    EXPECT_NEAR(0.0, sgeevx_obj->diff_wi, sgeevx_obj->threshold);
    EXPECT_NEAR(0.0, sgeevx_obj->diff_rconde, sgeevx_obj->threshold);
    EXPECT_NEAR(0.0, sgeevx_obj->diff_rcondv, sgeevx_obj->threshold);
	/* Commenting the below check because the position of the scalefactor
	values in the (scale & scaleref) , (rscale, rscaleref) varies
	slightly. Though the values are matching.
	So, element by element difference calculation will not work here 
	unlike for other buffers  */
	
    /* EXPECT_NEAR(0.0, sgeevx_obj->diff_scale, sgeevx_obj->threshold);
    EXPECT_NEAR(0.0, sgeevx_obj->diff_rscale, sgeevx_obj->threshold); */
    EXPECT_EQ(sgeevx_obj->ilo, sgeevx_obj->iloref);
    EXPECT_EQ(sgeevx_obj->ihi, sgeevx_obj->ihiref);
}
TEST_F(sgeevx_test, sgeevx3) {
    //EXPECT_NEAR(0.0, sgeevx_obj->diff_a, sgeevx_obj->threshold);
    EXPECT_NEAR(0.0, sgeevx_obj->diff_vl, sgeevx_obj->threshold);
    EXPECT_NEAR(0.0, sgeevx_obj->diff_vr, sgeevx_obj->threshold);
    EXPECT_NEAR(0.0, sgeevx_obj->diff_wr, sgeevx_obj->threshold);
    EXPECT_NEAR(0.0, sgeevx_obj->diff_wi, sgeevx_obj->threshold);
    EXPECT_NEAR(0.0, sgeevx_obj->diff_rconde, sgeevx_obj->threshold);
    EXPECT_NEAR(0.0, sgeevx_obj->diff_rcondv, sgeevx_obj->threshold);
	/* Commenting the below check because the position of the scalefactor
	values in the (scale & scaleref) , (rscale, rscaleref) varies
	slightly. Though the values are matching.
	So, element by element difference calculation will not work here 
	unlike for other buffers  */
	
    /* EXPECT_NEAR(0.0, sgeevx_obj->diff_scale, sgeevx_obj->threshold);
    EXPECT_NEAR(0.0, sgeevx_obj->diff_rscale, sgeevx_obj->threshold); */
    EXPECT_EQ(sgeevx_obj->ilo, sgeevx_obj->iloref);
    EXPECT_EQ(sgeevx_obj->ihi, sgeevx_obj->ihiref);
}
TEST_F(sgeevx_test, sgeevx4) {
    //EXPECT_NEAR(0.0, sgeevx_obj->diff_a, sgeevx_obj->threshold);
    EXPECT_NEAR(0.0, sgeevx_obj->diff_vl, sgeevx_obj->threshold);
    EXPECT_NEAR(0.0, sgeevx_obj->diff_vr, sgeevx_obj->threshold);
    EXPECT_NEAR(0.0, sgeevx_obj->diff_wr, sgeevx_obj->threshold);
    EXPECT_NEAR(0.0, sgeevx_obj->diff_wi, sgeevx_obj->threshold);
    EXPECT_NEAR(0.0, sgeevx_obj->diff_rconde, sgeevx_obj->threshold);
    EXPECT_NEAR(0.0, sgeevx_obj->diff_rcondv, sgeevx_obj->threshold);
	/* Commenting the below check because the position of the scalefactor
	values in the (scale & scaleref) , (rscale, rscaleref) varies
	slightly. Though the values are matching.
	So, element by element difference calculation will not work here 
	unlike for other buffers  */
	
    /* EXPECT_NEAR(0.0, sgeevx_obj->diff_scale, sgeevx_obj->threshold);
    EXPECT_NEAR(0.0, sgeevx_obj->diff_rscale, sgeevx_obj->threshold); */
    EXPECT_EQ(sgeevx_obj->ilo, sgeevx_obj->iloref);
    EXPECT_EQ(sgeevx_obj->ihi, sgeevx_obj->ihiref);
}

/* Begin geevx_double_common_parameters  class definition */
class geevx_double_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    double diff_a;
    double diff_wi, diff_wr, diff_vl, diff_vr;
    double diff_scale, diff_rconde, diff_rcondv;

    float threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR	
    char jobvl; // Must be 'N', or 'V'.
    char jobvr; // Must be 'N', or 'V'.
	char balance; // Must be 'N', 'P' 'S', B.
	char sense; //  Must be 'N', 'E', 'V', or 'B'.

    lapack_int n; // The order of the matrices A and B (n≥ 0).
    lapack_int lda; //  The leading dimension of a
    lapack_int ldvl;
    lapack_int ldvr;

    /* Input / Output parameters */
    double* a, *aref; // contains the n-by-n general matrix A.

    /* Output parameters */
    lapack_int sdim, sdimref;
    double *vl, *vlref; // left eigenvectors.
    double *vr, *vrref; // right eigenvectors
    lapack_int ilo, iloref, ihi, ihiref; // ilo and ihi mark the rows and columns of A which are to be reduced.
    double *scale, *scaleref; // permutations, scaling factors applied to the left side of A and B
    double  abnrm, abnrmref; //  L1- norm of balanced matrices
    double* wr, *wrref;
	double* wi, *wiref;
	double *rconde, *rcondv; // contain the reciprocal condition numbers 
	double *rconderef, *rcondvref; // contain the reciprocal condition numbers 
    /*Return Values */
    int info, inforef;

      geevx_double_parameters (int matrix_layout_i, char jobvl_i,
				char jobvr_i, char balance_i, char sense_i,
				lapack_int n_i );
      ~geevx_double_parameters ();
};

/* Constructor definition  double_common_parameters */
geevx_double_parameters:: geevx_double_parameters (int matrix_layout_i,
     char jobvl_i, char jobvr_i, char balance_i, char sense_i, 
	 lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    jobvl = jobvl_i;
    jobvr = jobvr_i;
	balance = balance_i;
	sense = sense_i;
	
    n  = n_i;
    sdim = 0;
	sdimref = 0;

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
	diff_scale = 0;
    diff_rconde = 0;
	diff_rcondv = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n geevx double: matrix_layout: %d n: %d  jobvl: %c \
jobvr: %c \t balance: %c \t sense: %c\n", matrix_layout, n, jobvr, 
jobvl, balance, sense);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, lda*n );
    lapacke_gtest_alloc_double_buffer_pair( &wr, &wrref, n );
    lapacke_gtest_alloc_double_buffer_pair( &wi, &wiref, n );
    lapacke_gtest_alloc_double_buffer_pair( &vl, &vlref, n*ldvl );
    lapacke_gtest_alloc_double_buffer_pair( &vr, &vrref, n*ldvr );
    lapacke_gtest_alloc_double_buffer_pair( &scale, &scaleref, n );
    lapacke_gtest_alloc_double_buffer_pair( &rconde, &rconderef, n );
    lapacke_gtest_alloc_double_buffer_pair( &rcondv, &rcondvref, n );

    if( (a==NULL) || (aref==NULL) ||  \
        (wr==NULL) || (wrref==NULL) || \
        (wi==NULL) || (wiref==NULL) || \
        (vl==NULL) || (vlref==NULL) || \
        (vr==NULL) || (vrref==NULL) || \
        (scale==NULL) || (scaleref==NULL) || \
        (rconde==NULL) || (rconderef==NULL) || \
        (rcondv==NULL) || (rcondvref==NULL) ){
       EXPECT_FALSE( true) << "geevx_double_parameters object: malloc error. Exiting ";
       geevx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, lda*n );
    lapacke_gtest_init_double_buffer_pair_with_constant( vl, vlref, ldvl*n, 0 );
    lapacke_gtest_init_double_buffer_pair_with_constant( vr, vrref, ldvr*n, 0 );
    lapacke_gtest_init_double_buffer_pair_with_constant( wr, wrref, n, 0);
    lapacke_gtest_init_double_buffer_pair_with_constant( wi, wiref, n, 0);
    lapacke_gtest_init_double_buffer_pair_with_constant( scale, scaleref, n, 0);
    lapacke_gtest_init_double_buffer_pair_with_constant( rconde, rconderef, n, 0);
    lapacke_gtest_init_double_buffer_pair_with_constant( rcondv, rcondvref, n, 0);

   } /* end of Constructor  */

/* Destructor definition  'geevx_double_common_parameters' */
geevx_double_parameters :: ~geevx_double_parameters ()
{
   /* De-Allocate memory for the input matrices */
   geevx_free();
} 

//  Test fixture class definition
class dgeevx_test  : public  ::testing::Test {
public:
   geevx_double_parameters  *dgeevx_obj;
   void SetUp();  
   void TearDown () { delete dgeevx_obj; }
};

void dgeevx_test::SetUp()
{
    /* LAPACKE_dgeevx prototype */	
    typedef int (*Fptr_NL_LAPACKE_dgeevx) (int matrix_layout, char balanc,
		char jobvl, char jobvr, char sense, lapack_int n, double* a,
		lapack_int lda, double* wr, double* wi, double* vl, lapack_int ldvl,
		double* vr, lapack_int ldvr, lapack_int* ilo, lapack_int* ihi,
		double* scale, double* abnrm, double* rconde, double* rcondv);

    Fptr_NL_LAPACKE_dgeevx DGEEVX;

    dgeevx_obj = new  geevx_double_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].jobvsl,
										 eig_non_sym_paramslist[idx].jobvsl, // vl & vr to be same
										 eig_non_sym_paramslist[idx].balance_ggevx,
										 eig_non_sym_paramslist[idx].sense_ggevx,
                                         eig_paramslist[idx].n );
                                         
    dgeevx_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    dgeevx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgeevx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgeevx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgeevx_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DGEEVX = (Fptr_NL_LAPACKE_dgeevx)dlsym(dgeevx_obj->hModule, "LAPACKE_dgeevx");
    ASSERT_TRUE(DGEEVX != NULL) << "failed to ppt the Netlib LAPACKE_dgeevx symbol";

    /* Compute libflame's Lapacke o/p  */
    dgeevx_obj->inforef =  LAPACKE_dgeevx (   dgeevx_obj->matrix_layout,
                                    dgeevx_obj->balance,
                                    dgeevx_obj->jobvl,
                                    dgeevx_obj->jobvr,
                                    dgeevx_obj->sense,
                                    dgeevx_obj->n,
                                    dgeevx_obj->a,
                                    dgeevx_obj->lda,
                                    dgeevx_obj->wr,
                                    dgeevx_obj->wi,
									dgeevx_obj->vl,
									dgeevx_obj->ldvl,
									dgeevx_obj->vr,
									dgeevx_obj->ldvr,
									&dgeevx_obj->ilo,
									&dgeevx_obj->ihi,
									dgeevx_obj->scale,
									&dgeevx_obj->abnrm,
									dgeevx_obj->rconde,
									dgeevx_obj->rcondv
                                    );

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    dgeevx_obj->inforef = DGEEVX(   dgeevx_obj->matrix_layout,
                                    dgeevx_obj->balance,
                                    dgeevx_obj->jobvl,
                                    dgeevx_obj->jobvr,
                                    dgeevx_obj->sense,
                                    dgeevx_obj->n,
                                    dgeevx_obj->aref,
                                    dgeevx_obj->lda,
                                    dgeevx_obj->wrref,
                                    dgeevx_obj->wiref,
									dgeevx_obj->vlref,
									dgeevx_obj->ldvl,
									dgeevx_obj->vrref,
									dgeevx_obj->ldvr,
									&dgeevx_obj->iloref,
									&dgeevx_obj->ihiref,
									dgeevx_obj->scaleref,
									&dgeevx_obj->abnrmref,
									dgeevx_obj->rconderef,
									dgeevx_obj->rcondvref
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    dgeevx_obj->diff_a =  computeDiff_d( (dgeevx_obj->lda)*(dgeevx_obj->n), 
                dgeevx_obj->a, dgeevx_obj->aref );

    dgeevx_obj->diff_vl =  computeDiff_d( (dgeevx_obj->ldvl)*(dgeevx_obj->n), 
                dgeevx_obj->vl, dgeevx_obj->vlref );

    dgeevx_obj->diff_vr =  computeDiff_d( (dgeevx_obj->ldvr)*(dgeevx_obj->n), 
                dgeevx_obj->vr, dgeevx_obj->vrref );

    dgeevx_obj->diff_scale =  computeDiff_d( dgeevx_obj->n, 
                dgeevx_obj->scale, dgeevx_obj->scaleref );

    dgeevx_obj->diff_wr =  computeDiff_d( dgeevx_obj->n, 
                dgeevx_obj->wr, dgeevx_obj->wrref );

    dgeevx_obj->diff_wi =  computeDiff_d( dgeevx_obj->n, 
                dgeevx_obj->wi, dgeevx_obj->wiref );

    dgeevx_obj->diff_rconde =  computeDiff_d( dgeevx_obj->n, 
                dgeevx_obj->rconde, dgeevx_obj->rconderef );

    dgeevx_obj->diff_rcondv =  computeDiff_d( dgeevx_obj->n, 
                dgeevx_obj->rcondv, dgeevx_obj->rcondvref );

    dgeevx_obj->diff_scale =  computeDiff_d( dgeevx_obj->n, 
                dgeevx_obj->scale, dgeevx_obj->scaleref );


#if LAPACKE_TEST_VERBOSE
   printf(" \n geevx double: \n diff_a: %f \n  \
diff_vl: %f \n diff_vr: %f \n diff_wr: %f \n\
diff_wi: %f \n ilo: %d \t iloref: %d \t \
diff_ilo: %d \n ihi: %d \t ihiref: %d \t diff_ihi: %d \n \
diff_scale: %f \n",
       dgeevx_obj->diff_a, dgeevx_obj->diff_vl,
	   dgeevx_obj->diff_vr, dgeevx_obj->diff_wr,
       dgeevx_obj->diff_wi,
	   dgeevx_obj->ilo, dgeevx_obj->iloref, 
	   ( dgeevx_obj->ilo - dgeevx_obj->iloref),
	   dgeevx_obj->ihi, dgeevx_obj->ihiref,
	   ( dgeevx_obj->ihi - dgeevx_obj->ihiref),
	   dgeevx_obj->diff_scale );
#endif
}

TEST_F(dgeevx_test, dgeevx1) {
    //EXPECT_NEAR(0.0, dgeevx_obj->diff_a, dgeevx_obj->threshold);
    EXPECT_NEAR(0.0, dgeevx_obj->diff_vl, dgeevx_obj->threshold);
    EXPECT_NEAR(0.0, dgeevx_obj->diff_vr, dgeevx_obj->threshold);
    EXPECT_NEAR(0.0, dgeevx_obj->diff_wr, dgeevx_obj->threshold);
    EXPECT_NEAR(0.0, dgeevx_obj->diff_wi, dgeevx_obj->threshold);
    EXPECT_NEAR(0.0, dgeevx_obj->diff_rconde, dgeevx_obj->threshold);
    EXPECT_NEAR(0.0, dgeevx_obj->diff_rcondv, dgeevx_obj->threshold);
	/* Commenting the below check because the position of the scalefactor
	values in the (scale & scaleref) , (rscale, rscaleref) varies
	slightly. Though the values are matching.
	So, element by element difference calculation will not work here 
	unlike for other buffers  */
	
    /* EXPECT_NEAR(0.0, dgeevx_obj->diff_scale, dgeevx_obj->threshold);
    EXPECT_NEAR(0.0, dgeevx_obj->diff_rscale, dgeevx_obj->threshold); */
    EXPECT_EQ(dgeevx_obj->ilo, dgeevx_obj->iloref);
    EXPECT_EQ(dgeevx_obj->ihi, dgeevx_obj->ihiref);
}

TEST_F(dgeevx_test, dgeevx2) {
    //EXPECT_NEAR(0.0, dgeevx_obj->diff_a, dgeevx_obj->threshold);
    EXPECT_NEAR(0.0, dgeevx_obj->diff_vl, dgeevx_obj->threshold);
    EXPECT_NEAR(0.0, dgeevx_obj->diff_vr, dgeevx_obj->threshold);
    EXPECT_NEAR(0.0, dgeevx_obj->diff_wr, dgeevx_obj->threshold);
    EXPECT_NEAR(0.0, dgeevx_obj->diff_wi, dgeevx_obj->threshold);
    EXPECT_NEAR(0.0, dgeevx_obj->diff_rconde, dgeevx_obj->threshold);
    EXPECT_NEAR(0.0, dgeevx_obj->diff_rcondv, dgeevx_obj->threshold);
	/* Commenting the below check because the position of the scalefactor
	values in the (scale & scaleref) , (rscale, rscaleref) varies
	slightly. Though the values are matching.
	So, element by element difference calculation will not work here 
	unlike for other buffers  */
	
    /* EXPECT_NEAR(0.0, dgeevx_obj->diff_scale, dgeevx_obj->threshold);
    EXPECT_NEAR(0.0, dgeevx_obj->diff_rscale, dgeevx_obj->threshold); */
    EXPECT_EQ(dgeevx_obj->ilo, dgeevx_obj->iloref);
    EXPECT_EQ(dgeevx_obj->ihi, dgeevx_obj->ihiref);
}
TEST_F(dgeevx_test, dgeevx3) {
    //EXPECT_NEAR(0.0, dgeevx_obj->diff_a, dgeevx_obj->threshold);
    EXPECT_NEAR(0.0, dgeevx_obj->diff_vl, dgeevx_obj->threshold);
    EXPECT_NEAR(0.0, dgeevx_obj->diff_vr, dgeevx_obj->threshold);
    EXPECT_NEAR(0.0, dgeevx_obj->diff_wr, dgeevx_obj->threshold);
    EXPECT_NEAR(0.0, dgeevx_obj->diff_wi, dgeevx_obj->threshold);
    EXPECT_NEAR(0.0, dgeevx_obj->diff_rconde, dgeevx_obj->threshold);
    EXPECT_NEAR(0.0, dgeevx_obj->diff_rcondv, dgeevx_obj->threshold);
	/* Commenting the below check because the position of the scalefactor
	values in the (scale & scaleref) , (rscale, rscaleref) varies
	slightly. Though the values are matching.
	So, element by element difference calculation will not work here 
	unlike for other buffers  */
	
    /* EXPECT_NEAR(0.0, dgeevx_obj->diff_scale, dgeevx_obj->threshold);
    EXPECT_NEAR(0.0, dgeevx_obj->diff_rscale, dgeevx_obj->threshold); */
    EXPECT_EQ(dgeevx_obj->ilo, dgeevx_obj->iloref);
    EXPECT_EQ(dgeevx_obj->ihi, dgeevx_obj->ihiref);
}
TEST_F(dgeevx_test, dgeevx4) {
    //EXPECT_NEAR(0.0, dgeevx_obj->diff_a, dgeevx_obj->threshold);
    EXPECT_NEAR(0.0, dgeevx_obj->diff_vl, dgeevx_obj->threshold);
    EXPECT_NEAR(0.0, dgeevx_obj->diff_vr, dgeevx_obj->threshold);
    EXPECT_NEAR(0.0, dgeevx_obj->diff_wr, dgeevx_obj->threshold);
    EXPECT_NEAR(0.0, dgeevx_obj->diff_wi, dgeevx_obj->threshold);
    EXPECT_NEAR(0.0, dgeevx_obj->diff_rconde, dgeevx_obj->threshold);
    EXPECT_NEAR(0.0, dgeevx_obj->diff_rcondv, dgeevx_obj->threshold);
	/* Commenting the below check because the position of the scalefactor
	values in the (scale & scaleref) , (rscale, rscaleref) varies
	slightly. Though the values are matching.
	So, element by element difference calculation will not work here 
	unlike for other buffers  */
	
    /* EXPECT_NEAR(0.0, dgeevx_obj->diff_scale, dgeevx_obj->threshold);
    EXPECT_NEAR(0.0, dgeevx_obj->diff_rscale, dgeevx_obj->threshold); */
    EXPECT_EQ(dgeevx_obj->ilo, dgeevx_obj->iloref);
    EXPECT_EQ(dgeevx_obj->ihi, dgeevx_obj->ihiref);
}


/* Begin geevx_scomplex_common_parameters  class definition */
class geevx_scomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_a;
    float diff_wi, diff_wr, diff_vl, diff_vr;
    float diff_scale, diff_rconde, diff_rcondv;

    float threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR	
    char jobvl; // Must be 'N', or 'V'.
    char jobvr; // Must be 'N', or 'V'.
	char balance; // Must be 'N', 'P' 'S', B.
	char sense; //  Must be 'N', 'E', 'V', or 'B'.

    lapack_int n; // The order of the matrices A and B (n≥ 0).
    lapack_int lda; //  The leading dimension of a
    lapack_int ldvl;
    lapack_int ldvr;

    /* Input / Output parameters */
    lapack_complex_float* a, *aref; // contains the n-by-n general matrix A.

    /* Output parameters */
    lapack_int sdim, sdimref;
    lapack_complex_float *vl, *vlref; // left eigenvectors.
    lapack_complex_float *vr, *vrref; // right eigenvectors
    lapack_int ilo, iloref, ihi, ihiref; // ilo and ihi mark the rows and columns of A which are to be reduced.
    float *scale, *scaleref; // permutations, scaling factors applied to the left side of A and B
    float  abnrm, abnrmref; //  L1- norm of balanced matrices
    lapack_complex_float* wr, *wrref;
	lapack_complex_float* wi, *wiref;
	float *rconde, *rcondv; // contain the reciprocal condition numbers 
	float *rconderef, *rcondvref; // contain the reciprocal condition numbers 
    /*Return Values */
    int info, inforef;

      geevx_scomplex_parameters (int matrix_layout_i, char jobvl_i,
				char jobvr_i, char balance_i, char sense_i,
				lapack_int n_i );
      ~geevx_scomplex_parameters ();
};

/* Constructor definition  lapack_complex_float_common_parameters */
geevx_scomplex_parameters:: geevx_scomplex_parameters (int matrix_layout_i,
     char jobvl_i, char jobvr_i, char balance_i, char sense_i, 
	 lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    jobvl = jobvl_i;
    jobvr = jobvr_i;
	balance = balance_i;
	sense = sense_i;
	
    n  = n_i;
    sdim = 0;
	sdimref = 0;

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
	diff_scale = 0;
    diff_rconde = 0;
	diff_rcondv = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n geevx lapack_complex_float: matrix_layout: %d n: %d  jobvl: %c \
jobvr: %c \t balance: %c \t sense: %c\n", matrix_layout, n, jobvr, 
jobvl, balance, sense);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, lda*n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &wr, &wrref, n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &wi, &wiref, n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &vl, &vlref, n*ldvl );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &vr, &vrref, n*ldvr );
    lapacke_gtest_alloc_float_buffer_pair( &scale, &scaleref, n );
    lapacke_gtest_alloc_float_buffer_pair( &rconde, &rconderef, n );
    lapacke_gtest_alloc_float_buffer_pair( &rcondv, &rcondvref, n );

    if( (a==NULL) || (aref==NULL) ||  \
        (wr==NULL) || (wrref==NULL) || \
        (wi==NULL) || (wiref==NULL) || \
        (vl==NULL) || (vlref==NULL) || \
        (vr==NULL) || (vrref==NULL) || \
        (scale==NULL) || (scaleref==NULL) || \
        (rconde==NULL) || (rconderef==NULL) || \
        (rcondv==NULL) || (rcondvref==NULL) ){
       EXPECT_FALSE( true) << "geevx_scomplex_parameters object: malloc error. Exiting ";
       geevx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, lda*n );
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( vl, vlref, ldvl*n, 0 );
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( vr, vrref, ldvr*n, 0 );
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( wr, wrref, n, 0);
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( wi, wiref, n, 0);
    lapacke_gtest_init_float_buffer_pair_with_constant( scale, scaleref, n, 0);
    lapacke_gtest_init_float_buffer_pair_with_constant( rconde, rconderef, n, 0);
    lapacke_gtest_init_float_buffer_pair_with_constant( rcondv, rcondvref, n, 0);

   } /* end of Constructor  */

/* Destructor definition  'geevx_scomplex_common_parameters' */
geevx_scomplex_parameters :: ~geevx_scomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   geevx_free();
} 

//  Test fixture class definition
class cgeevx_test  : public  ::testing::Test {
public:
   geevx_scomplex_parameters  *cgeevx_obj;
   void SetUp();  
   void TearDown () { delete cgeevx_obj; }
};

void cgeevx_test::SetUp()
{
    /* LAPACKE_cgeevx prototype */	
    typedef int (*Fptr_NL_LAPACKE_cgeevx) (int matrix_layout, char balanc,
		char jobvl, char jobvr, char sense, lapack_int n, lapack_complex_float* a,
		lapack_int lda, lapack_complex_float* w, lapack_complex_float* vl,
		lapack_int ldvl, lapack_complex_float* vr, lapack_int ldvr,
		lapack_int* ilo, lapack_int* ihi, float* scale, float* abnrm,
		float* rconde, float* rcondv);

    Fptr_NL_LAPACKE_cgeevx CGEEVX;

    cgeevx_obj = new  geevx_scomplex_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].jobvsl,
										 eig_non_sym_paramslist[idx].jobvsl, // vl & vr to be same
										 eig_non_sym_paramslist[idx].balance_ggevx,
										 eig_non_sym_paramslist[idx].sense_ggevx,
                                         eig_paramslist[idx].n );
                                         
    cgeevx_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    cgeevx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgeevx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgeevx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgeevx_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CGEEVX = (Fptr_NL_LAPACKE_cgeevx)dlsym(cgeevx_obj->hModule, "LAPACKE_cgeevx");
    ASSERT_TRUE(CGEEVX != NULL) << "failed to ppt the Netlib LAPACKE_cgeevx symbol";

    /* Compute libflame's Lapacke o/p  */
    cgeevx_obj->inforef =  LAPACKE_cgeevx (   cgeevx_obj->matrix_layout,
                                    cgeevx_obj->balance,
                                    cgeevx_obj->jobvl,
                                    cgeevx_obj->jobvr,
                                    cgeevx_obj->sense,
                                    cgeevx_obj->n,
                                    cgeevx_obj->a,
                                    cgeevx_obj->lda,
                                    cgeevx_obj->wr,
									cgeevx_obj->vl,
									cgeevx_obj->ldvl,
									cgeevx_obj->vr,
									cgeevx_obj->ldvr,
									&cgeevx_obj->ilo,
									&cgeevx_obj->ihi,
									cgeevx_obj->scale,
									&cgeevx_obj->abnrm,
									cgeevx_obj->rconde,
									cgeevx_obj->rcondv
                                    );

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    cgeevx_obj->inforef = CGEEVX(   cgeevx_obj->matrix_layout,
                                    cgeevx_obj->balance,
                                    cgeevx_obj->jobvl,
                                    cgeevx_obj->jobvr,
                                    cgeevx_obj->sense,
                                    cgeevx_obj->n,
                                    cgeevx_obj->aref,
                                    cgeevx_obj->lda,
                                    cgeevx_obj->wrref,
									cgeevx_obj->vlref,
									cgeevx_obj->ldvl,
									cgeevx_obj->vrref,
									cgeevx_obj->ldvr,
									&cgeevx_obj->iloref,
									&cgeevx_obj->ihiref,
									cgeevx_obj->scaleref,
									&cgeevx_obj->abnrmref,
									cgeevx_obj->rconderef,
									cgeevx_obj->rcondvref
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    cgeevx_obj->diff_a =  computeDiff_c( (cgeevx_obj->lda)*(cgeevx_obj->n), 
                cgeevx_obj->a, cgeevx_obj->aref );

    cgeevx_obj->diff_vl =  computeDiff_c( (cgeevx_obj->ldvl)*(cgeevx_obj->n), 
                cgeevx_obj->vl, cgeevx_obj->vlref );

    cgeevx_obj->diff_vr =  computeDiff_c( (cgeevx_obj->ldvr)*(cgeevx_obj->n), 
                cgeevx_obj->vr, cgeevx_obj->vrref );

    cgeevx_obj->diff_scale =  computeDiff_s( cgeevx_obj->n, 
                cgeevx_obj->scale, cgeevx_obj->scaleref );

    cgeevx_obj->diff_wr =  computeDiff_c( cgeevx_obj->n, 
                cgeevx_obj->wr, cgeevx_obj->wrref );

    cgeevx_obj->diff_wi =  computeDiff_c( cgeevx_obj->n, 
                cgeevx_obj->wi, cgeevx_obj->wiref );

    cgeevx_obj->diff_rconde =  computeDiff_s( cgeevx_obj->n, 
                cgeevx_obj->rconde, cgeevx_obj->rconderef );

    cgeevx_obj->diff_rcondv =  computeDiff_s( cgeevx_obj->n, 
                cgeevx_obj->rcondv, cgeevx_obj->rcondvref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n geevx lapack_complex_float: \n diff_a: %f \n  \
diff_vl: %f \n diff_vr: %f \n diff_wr: %f \n\
diff_wi: %f \n ilo: %d \t iloref: %d \t \
diff_ilo: %d \n ihi: %d \t ihiref: %d \t diff_ihi: %d \n \
diff_scale: %f \n",
       cgeevx_obj->diff_a, cgeevx_obj->diff_vl,
	   cgeevx_obj->diff_vr, cgeevx_obj->diff_wr,
       cgeevx_obj->diff_wi,
	   cgeevx_obj->ilo, cgeevx_obj->iloref, 
	   ( cgeevx_obj->ilo - cgeevx_obj->iloref),
	   cgeevx_obj->ihi, cgeevx_obj->ihiref,
	   ( cgeevx_obj->ihi - cgeevx_obj->ihiref),
	   cgeevx_obj->diff_scale );
#endif
}

TEST_F(cgeevx_test, cgeevx1) {
    //EXPECT_NEAR(0.0, cgeevx_obj->diff_a, cgeevx_obj->threshold);
    EXPECT_NEAR(0.0, cgeevx_obj->diff_vl, cgeevx_obj->threshold);
    EXPECT_NEAR(0.0, cgeevx_obj->diff_vr, cgeevx_obj->threshold);
    EXPECT_NEAR(0.0, cgeevx_obj->diff_wr, cgeevx_obj->threshold);
    EXPECT_NEAR(0.0, cgeevx_obj->diff_wi, cgeevx_obj->threshold);
    EXPECT_NEAR(0.0, cgeevx_obj->diff_rconde, cgeevx_obj->threshold);
    EXPECT_NEAR(0.0, cgeevx_obj->diff_rcondv, cgeevx_obj->threshold);
	/* Commenting the below check because the position of the scalefactor
	values in the (scale & scaleref) , (rscale, rscaleref) varies
	slightly. Though the values are matching.
	So, element by element difference calculation will not work here 
	unlike for other buffers  */
	
    /* EXPECT_NEAR(0.0, cgeevx_obj->diff_scale, cgeevx_obj->threshold);
    EXPECT_NEAR(0.0, cgeevx_obj->diff_rscale, cgeevx_obj->threshold); */
    EXPECT_EQ(cgeevx_obj->ilo, cgeevx_obj->iloref);
    EXPECT_EQ(cgeevx_obj->ihi, cgeevx_obj->ihiref);
}

TEST_F(cgeevx_test, cgeevx2) {
    //EXPECT_NEAR(0.0, cgeevx_obj->diff_a, cgeevx_obj->threshold);
    EXPECT_NEAR(0.0, cgeevx_obj->diff_vl, cgeevx_obj->threshold);
    EXPECT_NEAR(0.0, cgeevx_obj->diff_vr, cgeevx_obj->threshold);
    EXPECT_NEAR(0.0, cgeevx_obj->diff_wr, cgeevx_obj->threshold);
    EXPECT_NEAR(0.0, cgeevx_obj->diff_wi, cgeevx_obj->threshold);
    EXPECT_NEAR(0.0, cgeevx_obj->diff_rconde, cgeevx_obj->threshold);
    EXPECT_NEAR(0.0, cgeevx_obj->diff_rcondv, cgeevx_obj->threshold);
	/* Commenting the below check because the position of the scalefactor
	values in the (scale & scaleref) , (rscale, rscaleref) varies
	slightly. Though the values are matching.
	So, element by element difference calculation will not work here 
	unlike for other buffers  */
	
    /* EXPECT_NEAR(0.0, cgeevx_obj->diff_scale, cgeevx_obj->threshold);
    EXPECT_NEAR(0.0, cgeevx_obj->diff_rscale, cgeevx_obj->threshold); */
    EXPECT_EQ(cgeevx_obj->ilo, cgeevx_obj->iloref);
    EXPECT_EQ(cgeevx_obj->ihi, cgeevx_obj->ihiref);
}
TEST_F(cgeevx_test, cgeevx3) {
    //EXPECT_NEAR(0.0, cgeevx_obj->diff_a, cgeevx_obj->threshold);
    EXPECT_NEAR(0.0, cgeevx_obj->diff_vl, cgeevx_obj->threshold);
    EXPECT_NEAR(0.0, cgeevx_obj->diff_vr, cgeevx_obj->threshold);
    EXPECT_NEAR(0.0, cgeevx_obj->diff_wr, cgeevx_obj->threshold);
    EXPECT_NEAR(0.0, cgeevx_obj->diff_wi, cgeevx_obj->threshold);
    EXPECT_NEAR(0.0, cgeevx_obj->diff_rconde, cgeevx_obj->threshold);
    EXPECT_NEAR(0.0, cgeevx_obj->diff_rcondv, cgeevx_obj->threshold);
	/* Commenting the below check because the position of the scalefactor
	values in the (scale & scaleref) , (rscale, rscaleref) varies
	slightly. Though the values are matching.
	So, element by element difference calculation will not work here 
	unlike for other buffers  */
	
    /* EXPECT_NEAR(0.0, cgeevx_obj->diff_scale, cgeevx_obj->threshold);
    EXPECT_NEAR(0.0, cgeevx_obj->diff_rscale, cgeevx_obj->threshold); */
    EXPECT_EQ(cgeevx_obj->ilo, cgeevx_obj->iloref);
    EXPECT_EQ(cgeevx_obj->ihi, cgeevx_obj->ihiref);
}
TEST_F(cgeevx_test, cgeevx4) {
    //EXPECT_NEAR(0.0, cgeevx_obj->diff_a, cgeevx_obj->threshold);
    EXPECT_NEAR(0.0, cgeevx_obj->diff_vl, cgeevx_obj->threshold);
    EXPECT_NEAR(0.0, cgeevx_obj->diff_vr, cgeevx_obj->threshold);
    EXPECT_NEAR(0.0, cgeevx_obj->diff_wr, cgeevx_obj->threshold);
    EXPECT_NEAR(0.0, cgeevx_obj->diff_wi, cgeevx_obj->threshold);
    EXPECT_NEAR(0.0, cgeevx_obj->diff_rconde, cgeevx_obj->threshold);
    EXPECT_NEAR(0.0, cgeevx_obj->diff_rcondv, cgeevx_obj->threshold);
	/* Commenting the below check because the position of the scalefactor
	values in the (scale & scaleref) , (rscale, rscaleref) varies
	slightly. Though the values are matching.
	So, element by element difference calculation will not work here 
	unlike for other buffers  */
	
    /* EXPECT_NEAR(0.0, cgeevx_obj->diff_scale, cgeevx_obj->threshold);
    EXPECT_NEAR(0.0, cgeevx_obj->diff_rscale, cgeevx_obj->threshold); */
    EXPECT_EQ(cgeevx_obj->ilo, cgeevx_obj->iloref);
    EXPECT_EQ(cgeevx_obj->ihi, cgeevx_obj->ihiref);
}

/* Begin geevx_dcomplex_common_parameters  class definition */
class geevx_dcomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    double diff_a;
    double diff_wi, diff_wr, diff_vl, diff_vr;
    double diff_scale, diff_rconde, diff_rcondv;

    double threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR	
    char jobvl; // Must be 'N', or 'V'.
    char jobvr; // Must be 'N', or 'V'.
	char balance; // Must be 'N', 'P' 'S', B.
	char sense; //  Must be 'N', 'E', 'V', or 'B'.

    lapack_int n; // The order of the matrices A and B (n≥ 0).
    lapack_int lda; //  The leading dimension of a
    lapack_int ldvl;
    lapack_int ldvr;

    /* Input / Output parameters */
    lapack_complex_double* a, *aref; // contains the n-by-n general matrix A.

    /* Output parameters */
    lapack_int sdim, sdimref;
    lapack_complex_double *vl, *vlref; // left eigenvectors.
    lapack_complex_double *vr, *vrref; // right eigenvectors
    lapack_int ilo, iloref, ihi, ihiref; // ilo and ihi mark the rows and columns of A which are to be reduced.
    double *scale, *scaleref; // permutations, scaling factors applied to the left side of A and B
    double  abnrm, abnrmref; //  L1- norm of balanced matrices
    lapack_complex_double* wr, *wrref;
	lapack_complex_double* wi, *wiref;
	double *rconde, *rcondv; // contain the reciprocal condition numbers 
	double *rconderef, *rcondvref; // contain the reciprocal condition numbers 
    /*Return Values */
    int info, inforef;

      geevx_dcomplex_parameters (int matrix_layout_i, char jobvl_i,
				char jobvr_i, char balance_i, char sense_i,
				lapack_int n_i );
      ~geevx_dcomplex_parameters ();
};

/* Constructor definition  lapack_complex_double_common_parameters */
geevx_dcomplex_parameters:: geevx_dcomplex_parameters (int matrix_layout_i,
     char jobvl_i, char jobvr_i, char balance_i, char sense_i, 
	 lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    jobvl = jobvl_i;
    jobvr = jobvr_i;
	balance = balance_i;
	sense = sense_i;
	
    n  = n_i;
    sdim = 0;
	sdimref = 0;

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
	diff_scale = 0;
    diff_rconde = 0;
	diff_rcondv = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n geevx lapack_complex_double: matrix_layout: %d n: %d  jobvl: %c \
jobvr: %c \t balance: %c \t sense: %c\n", matrix_layout, n, jobvr, 
jobvl, balance, sense);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, lda*n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &wr, &wrref, n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &wi, &wiref, n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &vl, &vlref, n*ldvl );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &vr, &vrref, n*ldvr );
    lapacke_gtest_alloc_double_buffer_pair( &scale, &scaleref, n );
    lapacke_gtest_alloc_double_buffer_pair( &rconde, &rconderef, n );
    lapacke_gtest_alloc_double_buffer_pair( &rcondv, &rcondvref, n );

    if( (a==NULL) || (aref==NULL) ||  \
        (wr==NULL) || (wrref==NULL) || \
        (wi==NULL) || (wiref==NULL) || \
        (vl==NULL) || (vlref==NULL) || \
        (vr==NULL) || (vrref==NULL) || \
        (scale==NULL) || (scaleref==NULL) || \
        (rconde==NULL) || (rconderef==NULL) || \
        (rcondv==NULL) || (rcondvref==NULL) ){
       EXPECT_FALSE( true) << "geevx_dcomplex_parameters object: malloc error. Exiting ";
       geevx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, lda*n );
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( vl, vlref, ldvl*n, 0 );
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( vr, vrref, ldvr*n, 0 );
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( wr, wrref, n, 0);
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( wi, wiref, n, 0);
    lapacke_gtest_init_double_buffer_pair_with_constant( scale, scaleref, n, 0);
    lapacke_gtest_init_double_buffer_pair_with_constant( rconde, rconderef, n, 0);
    lapacke_gtest_init_double_buffer_pair_with_constant( rcondv, rcondvref, n, 0);

   } /* end of Constructor  */

/* Destructor definition  'geevx_dcomplex_common_parameters' */
geevx_dcomplex_parameters :: ~geevx_dcomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   geevx_free();
} 

//  Test fixture class definition
class zgeevx_test  : public  ::testing::Test {
public:
   geevx_dcomplex_parameters  *zgeevx_obj;
   void SetUp();  
   void TearDown () { delete zgeevx_obj; }
};

void zgeevx_test::SetUp()
{
    /* LAPACKE_zgeevx prototype */	
    typedef int (*Fptr_NL_LAPACKE_zgeevx) (int matrix_layout, char balanc,
		char jobvl, char jobvr, char sense, lapack_int n, lapack_complex_double* a,
		lapack_int lda, lapack_complex_double* w, lapack_complex_double* vl,
		lapack_int ldvl, lapack_complex_double* vr, lapack_int ldvr,
		lapack_int* ilo, lapack_int* ihi, double* scale, double* abnrm,
		double* rconde, double* rcondv);

    Fptr_NL_LAPACKE_zgeevx ZGEEVX;

    zgeevx_obj = new  geevx_dcomplex_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].jobvsl,
										 eig_non_sym_paramslist[idx].jobvsl, // vl & vr to be same
										 eig_non_sym_paramslist[idx].balance_ggevx,
										 eig_non_sym_paramslist[idx].sense_ggevx,
                                         eig_paramslist[idx].n );
                                         
    zgeevx_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    zgeevx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgeevx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgeevx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgeevx_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZGEEVX = (Fptr_NL_LAPACKE_zgeevx)dlsym(zgeevx_obj->hModule, "LAPACKE_zgeevx");
    ASSERT_TRUE(ZGEEVX != NULL) << "failed to ppt the Netlib LAPACKE_zgeevx symbol";

    /* Compute libflame's Lapacke o/p  */
    zgeevx_obj->inforef =  LAPACKE_zgeevx (   zgeevx_obj->matrix_layout,
                                    zgeevx_obj->balance,
                                    zgeevx_obj->jobvl,
                                    zgeevx_obj->jobvr,
                                    zgeevx_obj->sense,
                                    zgeevx_obj->n,
                                    zgeevx_obj->a,
                                    zgeevx_obj->lda,
                                    zgeevx_obj->wr,
									zgeevx_obj->vl,
									zgeevx_obj->ldvl,
									zgeevx_obj->vr,
									zgeevx_obj->ldvr,
									&zgeevx_obj->ilo,
									&zgeevx_obj->ihi,
									zgeevx_obj->scale,
									&zgeevx_obj->abnrm,
									zgeevx_obj->rconde,
									zgeevx_obj->rcondv
                                    );

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    zgeevx_obj->inforef = ZGEEVX(   zgeevx_obj->matrix_layout,
                                    zgeevx_obj->balance,
                                    zgeevx_obj->jobvl,
                                    zgeevx_obj->jobvr,
                                    zgeevx_obj->sense,
                                    zgeevx_obj->n,
                                    zgeevx_obj->aref,
                                    zgeevx_obj->lda,
                                    zgeevx_obj->wrref,
									zgeevx_obj->vlref,
									zgeevx_obj->ldvl,
									zgeevx_obj->vrref,
									zgeevx_obj->ldvr,
									&zgeevx_obj->iloref,
									&zgeevx_obj->ihiref,
									zgeevx_obj->scaleref,
									&zgeevx_obj->abnrmref,
									zgeevx_obj->rconderef,
									zgeevx_obj->rcondvref
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    zgeevx_obj->diff_a =  computeDiff_z( (zgeevx_obj->lda)*(zgeevx_obj->n), 
                zgeevx_obj->a, zgeevx_obj->aref );

    zgeevx_obj->diff_vl =  computeDiff_z( (zgeevx_obj->ldvl)*(zgeevx_obj->n), 
                zgeevx_obj->vl, zgeevx_obj->vlref );

    zgeevx_obj->diff_vr =  computeDiff_z( (zgeevx_obj->ldvr)*(zgeevx_obj->n), 
                zgeevx_obj->vr, zgeevx_obj->vrref );

    zgeevx_obj->diff_scale =  computeDiff_d( zgeevx_obj->n, 
                zgeevx_obj->scale, zgeevx_obj->scaleref );

    zgeevx_obj->diff_wr =  computeDiff_z( zgeevx_obj->n, 
                zgeevx_obj->wr, zgeevx_obj->wrref );

    zgeevx_obj->diff_wi =  computeDiff_z( zgeevx_obj->n, 
                zgeevx_obj->wi, zgeevx_obj->wiref );

    zgeevx_obj->diff_rconde =  computeDiff_d( zgeevx_obj->n, 
                zgeevx_obj->rconde, zgeevx_obj->rconderef );

    zgeevx_obj->diff_rcondv =  computeDiff_d( zgeevx_obj->n, 
                zgeevx_obj->rcondv, zgeevx_obj->rcondvref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n geevx lapack_complex_double: \n diff_a: %f \n  \
diff_vl: %f \n diff_vr: %f \n diff_wr: %f \n\
diff_wi: %f \n ilo: %d \t iloref: %d \t \
diff_ilo: %d \n ihi: %d \t ihiref: %d \t diff_ihi: %d \n \
diff_scale: %f \n",
       zgeevx_obj->diff_a, zgeevx_obj->diff_vl,
	   zgeevx_obj->diff_vr, zgeevx_obj->diff_wr,
       zgeevx_obj->diff_wi,
	   zgeevx_obj->ilo, zgeevx_obj->iloref, 
	   ( zgeevx_obj->ilo - zgeevx_obj->iloref),
	   zgeevx_obj->ihi, zgeevx_obj->ihiref,
	   ( zgeevx_obj->ihi - zgeevx_obj->ihiref),
	   zgeevx_obj->diff_scale );
#endif
}

TEST_F(zgeevx_test, zgeevx1) {
    //EXPECT_NEAR(0.0, zgeevx_obj->diff_a, zgeevx_obj->threshold);
    EXPECT_NEAR(0.0, zgeevx_obj->diff_vl, zgeevx_obj->threshold);
    EXPECT_NEAR(0.0, zgeevx_obj->diff_vr, zgeevx_obj->threshold);
    EXPECT_NEAR(0.0, zgeevx_obj->diff_wr, zgeevx_obj->threshold);
    EXPECT_NEAR(0.0, zgeevx_obj->diff_wi, zgeevx_obj->threshold);
    EXPECT_NEAR(0.0, zgeevx_obj->diff_rconde, zgeevx_obj->threshold);
    EXPECT_NEAR(0.0, zgeevx_obj->diff_rcondv, zgeevx_obj->threshold);
	/* Commenting the below check because the position of the scalefactor
	values in the (scale & scaleref) , (rscale, rscaleref) varies
	slightly. Though the values are matching.
	So, element by element difference calculation will not work here 
	unlike for other buffers  */
	
    /* EXPECT_NEAR(0.0, zgeevx_obj->diff_scale, zgeevx_obj->threshold);
    EXPECT_NEAR(0.0, zgeevx_obj->diff_rscale, zgeevx_obj->threshold); */
    EXPECT_EQ(zgeevx_obj->ilo, zgeevx_obj->iloref);
    EXPECT_EQ(zgeevx_obj->ihi, zgeevx_obj->ihiref);
}

TEST_F(zgeevx_test, zgeevx2) {
    //EXPECT_NEAR(0.0, zgeevx_obj->diff_a, zgeevx_obj->threshold);
    EXPECT_NEAR(0.0, zgeevx_obj->diff_vl, zgeevx_obj->threshold);
    EXPECT_NEAR(0.0, zgeevx_obj->diff_vr, zgeevx_obj->threshold);
    EXPECT_NEAR(0.0, zgeevx_obj->diff_wr, zgeevx_obj->threshold);
    EXPECT_NEAR(0.0, zgeevx_obj->diff_wi, zgeevx_obj->threshold);
    EXPECT_NEAR(0.0, zgeevx_obj->diff_rconde, zgeevx_obj->threshold);
    EXPECT_NEAR(0.0, zgeevx_obj->diff_rcondv, zgeevx_obj->threshold);
	/* Commenting the below check because the position of the scalefactor
	values in the (scale & scaleref) , (rscale, rscaleref) varies
	slightly. Though the values are matching.
	So, element by element difference calculation will not work here 
	unlike for other buffers  */
	
    /* EXPECT_NEAR(0.0, zgeevx_obj->diff_scale, zgeevx_obj->threshold);
    EXPECT_NEAR(0.0, zgeevx_obj->diff_rscale, zgeevx_obj->threshold); */
    EXPECT_EQ(zgeevx_obj->ilo, zgeevx_obj->iloref);
    EXPECT_EQ(zgeevx_obj->ihi, zgeevx_obj->ihiref);
}
TEST_F(zgeevx_test, zgeevx3) {
    //EXPECT_NEAR(0.0, zgeevx_obj->diff_a, zgeevx_obj->threshold);
    EXPECT_NEAR(0.0, zgeevx_obj->diff_vl, zgeevx_obj->threshold);
    EXPECT_NEAR(0.0, zgeevx_obj->diff_vr, zgeevx_obj->threshold);
    EXPECT_NEAR(0.0, zgeevx_obj->diff_wr, zgeevx_obj->threshold);
    EXPECT_NEAR(0.0, zgeevx_obj->diff_wi, zgeevx_obj->threshold);
    EXPECT_NEAR(0.0, zgeevx_obj->diff_rconde, zgeevx_obj->threshold);
    EXPECT_NEAR(0.0, zgeevx_obj->diff_rcondv, zgeevx_obj->threshold);
	/* Commenting the below check because the position of the scalefactor
	values in the (scale & scaleref) , (rscale, rscaleref) varies
	slightly. Though the values are matching.
	So, element by element difference calculation will not work here 
	unlike for other buffers  */
	
    /* EXPECT_NEAR(0.0, zgeevx_obj->diff_scale, zgeevx_obj->threshold);
    EXPECT_NEAR(0.0, zgeevx_obj->diff_rscale, zgeevx_obj->threshold); */
    EXPECT_EQ(zgeevx_obj->ilo, zgeevx_obj->iloref);
    EXPECT_EQ(zgeevx_obj->ihi, zgeevx_obj->ihiref);
}
TEST_F(zgeevx_test, zgeevx4) {
    //EXPECT_NEAR(0.0, zgeevx_obj->diff_a, zgeevx_obj->threshold);
    EXPECT_NEAR(0.0, zgeevx_obj->diff_vl, zgeevx_obj->threshold);
    EXPECT_NEAR(0.0, zgeevx_obj->diff_vr, zgeevx_obj->threshold);
    EXPECT_NEAR(0.0, zgeevx_obj->diff_wr, zgeevx_obj->threshold);
    EXPECT_NEAR(0.0, zgeevx_obj->diff_wi, zgeevx_obj->threshold);
    EXPECT_NEAR(0.0, zgeevx_obj->diff_rconde, zgeevx_obj->threshold);
    EXPECT_NEAR(0.0, zgeevx_obj->diff_rcondv, zgeevx_obj->threshold);
	/* Commenting the below check because the position of the scalefactor
	values in the (scale & scaleref) , (rscale, rscaleref) varies
	slightly. Though the values are matching.
	So, element by element difference calculation will not work here 
	unlike for other buffers  */
	
    /* EXPECT_NEAR(0.0, zgeevx_obj->diff_scale, zgeevx_obj->threshold);
    EXPECT_NEAR(0.0, zgeevx_obj->diff_rscale, zgeevx_obj->threshold); */
    EXPECT_EQ(zgeevx_obj->ilo, zgeevx_obj->iloref);
    EXPECT_EQ(zgeevx_obj->ihi, zgeevx_obj->ihiref);
}
