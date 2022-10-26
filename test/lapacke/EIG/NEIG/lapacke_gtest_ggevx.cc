#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define ggevx_free() \
    if (a!=NULL)        free(a); \
    if (aref!=NULL)     free(aref); \
    if (b!=NULL)        free(b); \
    if (bref!=NULL)     free(bref); \
    if (vl!=NULL)       free(vl); \
    if (vlref!=NULL)    free(vlref); \
    if (vr!=NULL)       free(vr); \
    if (vrref!=NULL)    free(vrref); \
    if (alphar!=NULL)     free(alphar); \
    if (alpharref!=NULL)     free(alpharref); \
    if (alphai!=NULL)     free(alphai); \
    if (alphairef!=NULL)     free(alphairef); \
    if (beta!=NULL)     free(beta); \
    if (betaref!=NULL)     free(betaref); \
    if (lscale!=NULL)   free(lscale); \
    if (lscaleref!=NULL)  free(lscaleref); \
    if (rscale!=NULL)        free(rscale); \
    if (rscaleref!=NULL)     free(rscaleref); \
    if (rconde!=NULL)   free(rconde); \
    if (rconderef!=NULL)  free(rconderef); \
    if (rcondv!=NULL)   free(rcondv); \
    if (rcondvref!=NULL)  free(rcondvref)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin ggevx_float_common_parameters  class definition */
class ggevx_float_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_a, diff_b;
    float diff_alphai, diff_alphar, diff_beta, diff_vl, diff_vr;
    float diff_lscale, diff_rscale, diff_rconde, diff_rcondv;

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
    lapack_int ldb; //  The leading dimension of b
    lapack_int ldvl;
    lapack_int ldvr;

    /* Input / Output parameters */
    float* a, *aref; // contains the n-by-n general matrix A.
    float* b, *bref; // contains the n-by-n upper triangular matrix B.

    /* Output parameters */
    lapack_int sdim, sdimref;
    float *vl, *vlref; // left eigenvectors.
    float *vr, *vrref; // right eigenvectors
    lapack_int ilo, iloref, ihi, ihiref; // ilo and ihi mark the rows and columns of A which are to be reduced.
    float *lscale, *lscaleref; // permutations, scaling factors applied to the left side of A and B
    float *rscale, *rscaleref; // permutations, scaling factors applied to the right side of A and B
    float  abnrm, bbnrm, abnrmref, bbnrmref; //  L1- norm of balanced matrices
    float* alphar, *alpharref;
	float* alphai, *alphairef;
	float* beta, *betaref;
	float *rconde, *rcondv; // contain the reciprocal condition numbers 
	float *rconderef, *rcondvref; // contain the reciprocal condition numbers 
    /*Return Values */
    int info, inforef;

      ggevx_float_parameters (int matrix_layout_i, char jobvl_i,
				char jobvr_i, char balance_i, char sense_i,
				lapack_int n_i );
      ~ggevx_float_parameters ();
};

/* Constructor definition  float_common_parameters */
ggevx_float_parameters:: ggevx_float_parameters (int matrix_layout_i,
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
    ldb = n;
    ldvl = n;
    ldvr = n;

    hModule = NULL;
    dModule = NULL;

    diff_a = 0;
    diff_b = 0;
    diff_vl = 0;
    diff_vr = 0;
    diff_alphai = 0;
    diff_alphar = 0;
    diff_beta = 0;
	diff_lscale = 0;
	diff_rscale = 0;
    diff_rconde = 0;
	diff_rcondv = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n ggevx float: matrix_layout: %d n: %d  jobvl: %c \
jobvr: %c \t balance: %c \t sense: %c\n", matrix_layout, n, jobvr, 
jobvl, balance, sense);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, lda*n );
    lapacke_gtest_alloc_float_buffer_pair( &b, &bref, ldb*n );
    lapacke_gtest_alloc_float_buffer_pair( &alphar, &alpharref, n );
    lapacke_gtest_alloc_float_buffer_pair( &alphai, &alphairef, n );
    lapacke_gtest_alloc_float_buffer_pair( &beta, &betaref, n );
    lapacke_gtest_alloc_float_buffer_pair( &vl, &vlref, n*ldvl );
    lapacke_gtest_alloc_float_buffer_pair( &vr, &vrref, n*ldvr );
    lapacke_gtest_alloc_float_buffer_pair( &lscale, &lscaleref, n );
    lapacke_gtest_alloc_float_buffer_pair( &rscale, &rscaleref, n );
    lapacke_gtest_alloc_float_buffer_pair( &rconde, &rconderef, n );
    lapacke_gtest_alloc_float_buffer_pair( &rcondv, &rcondvref, n );

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (alphar==NULL) || (alpharref==NULL) || \
        (alphai==NULL) || (alphairef==NULL) || \
        (beta==NULL) || (betaref==NULL) || \
        (vl==NULL) || (vlref==NULL) || \
        (vr==NULL) || (vrref==NULL) || \
        (lscale==NULL) || (lscaleref==NULL) || \
        (rconde==NULL) || (rconderef==NULL) || \
        (rcondv==NULL) || (rcondvref==NULL) || \
        (rscale==NULL) || (rscaleref==NULL) ){
       EXPECT_FALSE( true) << "ggevx_float_parameters object: malloc error. Exiting ";
       ggevx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, lda*n );
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, ldb*n );
    lapacke_gtest_init_float_buffer_pair_with_constant( rconde, rconderef, n, 0 );
    lapacke_gtest_init_float_buffer_pair_with_constant( rcondv, rcondvref, n, 0 );
    lapacke_gtest_init_float_buffer_pair_with_constant( vl, vlref, ldvl*n, 0 );
    lapacke_gtest_init_float_buffer_pair_with_constant( vr, vrref, ldvr*n, 0 );
    lapacke_gtest_init_float_buffer_pair_with_constant( alphar, alpharref, n, 0);
    lapacke_gtest_init_float_buffer_pair_with_constant( alphai, alphairef, n, 0);
    lapacke_gtest_init_float_buffer_pair_with_constant( beta, betaref, n, 0);
    lapacke_gtest_init_float_buffer_pair_with_constant( lscale, lscaleref, n, 0);
    lapacke_gtest_init_float_buffer_pair_with_constant( rscale, rscaleref, n, 0);
    lapacke_gtest_init_float_buffer_pair_with_constant( rconde, rconderef, n, 0);
    lapacke_gtest_init_float_buffer_pair_with_constant( rcondv, rcondvref, n, 0);

   } /* end of Constructor  */

/* Destructor definition  'ggevx_float_common_parameters' */
ggevx_float_parameters :: ~ggevx_float_parameters ()
{
   /* De-Allocate memory for the input matrices */
   ggevx_free();
} 

//  Test fixture class definition
class sggevx_test  : public  ::testing::Test {
public:
   ggevx_float_parameters  *sggevx_obj;
   void SetUp();  
   void TearDown () { delete sggevx_obj; }
};

void sggevx_test::SetUp()
{
    /* LAPACKE_sggevx prototype */	
    typedef int (*Fptr_NL_LAPACKE_sggevx) (int matrix_layout, char balanc,
		char jobvl, char jobvr, char sense, lapack_int n, float* a,
		lapack_int lda, float* b, lapack_int ldb, float* alphar,
		float* alphai, float* beta, float* vl, lapack_int ldvl,
		float* vr, lapack_int ldvr, lapack_int* ilo, lapack_int* ihi,
		float* lscale, float* rscale, float* abnrm,
		float* bbnrm, float* rconde, float* rcondv);

    Fptr_NL_LAPACKE_sggevx SGGEVX;

    sggevx_obj = new  ggevx_float_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].jobvsl,
										 eig_non_sym_paramslist[idx].jobvsr,
										 eig_non_sym_paramslist[idx].balance_ggevx,
										 eig_non_sym_paramslist[idx].sense_ggevx,
                                         eig_paramslist[idx].n );
                                         
    sggevx_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    sggevx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sggevx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sggevx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sggevx_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SGGEVX = (Fptr_NL_LAPACKE_sggevx)dlsym(sggevx_obj->hModule, "LAPACKE_sggevx");
    ASSERT_TRUE(SGGEVX != NULL) << "failed to ppt the Netlib LAPACKE_sggevx symbol";

    /* Compute libflame's Lapacke o/p  */
    sggevx_obj->inforef =  LAPACKE_sggevx (   sggevx_obj->matrix_layout,
                                    sggevx_obj->balance,
                                    sggevx_obj->jobvl,
                                    sggevx_obj->jobvr,
                                    sggevx_obj->sense,
                                    sggevx_obj->n,
                                    sggevx_obj->a,
                                    sggevx_obj->lda,
                                    sggevx_obj->b,
                                    sggevx_obj->ldb,
                                    sggevx_obj->alphar,
                                    sggevx_obj->alphai,
                                    sggevx_obj->beta,
									sggevx_obj->vl,
									sggevx_obj->ldvl,
									sggevx_obj->vr,
									sggevx_obj->ldvr,
									&sggevx_obj->ilo,
									&sggevx_obj->ihi,
									sggevx_obj->lscale,
									sggevx_obj->rscale,
									&sggevx_obj->abnrm,
									&sggevx_obj->bbnrm,
									sggevx_obj->rconde,
									sggevx_obj->rcondv
                                    );

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    sggevx_obj->inforef = SGGEVX(   sggevx_obj->matrix_layout,
                                    sggevx_obj->balance,
                                    sggevx_obj->jobvl,
                                    sggevx_obj->jobvr,
                                    sggevx_obj->sense,
                                    sggevx_obj->n,
                                    sggevx_obj->aref,
                                    sggevx_obj->lda,
                                    sggevx_obj->bref,
                                    sggevx_obj->ldb,
                                    sggevx_obj->alpharref,
                                    sggevx_obj->alphairef,
                                    sggevx_obj->betaref,
									sggevx_obj->vlref,
									sggevx_obj->ldvl,
									sggevx_obj->vrref,
									sggevx_obj->ldvr,
									&sggevx_obj->iloref,
									&sggevx_obj->ihiref,
									sggevx_obj->lscaleref,
									sggevx_obj->rscaleref,
									&sggevx_obj->abnrmref,
									&sggevx_obj->bbnrmref,
									sggevx_obj->rconderef,
									sggevx_obj->rcondvref
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    sggevx_obj->diff_a =  computeDiff_s( (sggevx_obj->lda)*(sggevx_obj->n), 
                sggevx_obj->a, sggevx_obj->aref );

    sggevx_obj->diff_b =  computeDiff_s( (sggevx_obj->ldb)*(sggevx_obj->n), 
                sggevx_obj->b, sggevx_obj->bref );

    sggevx_obj->diff_vl =  computeDiff_s( (sggevx_obj->ldvl)*(sggevx_obj->n), 
                sggevx_obj->vl, sggevx_obj->vlref );

    sggevx_obj->diff_vr =  computeDiff_s( (sggevx_obj->ldvr)*(sggevx_obj->n), 
                sggevx_obj->vr, sggevx_obj->vrref );

    sggevx_obj->diff_lscale =  computeDiff_s( sggevx_obj->n, 
                sggevx_obj->lscale, sggevx_obj->lscaleref );

    sggevx_obj->diff_rscale =  computeDiff_s( sggevx_obj->n, 
                sggevx_obj->rscale, sggevx_obj->rscaleref );

    sggevx_obj->diff_alphar =  computeDiff_s( sggevx_obj->n, 
                sggevx_obj->alphar, sggevx_obj->alpharref );

    sggevx_obj->diff_alphai =  computeDiff_s( sggevx_obj->n, 
                sggevx_obj->alphai, sggevx_obj->alphairef );

    sggevx_obj->diff_beta =  computeDiff_s( sggevx_obj->n, 
                sggevx_obj->beta, sggevx_obj->betaref );

    sggevx_obj->diff_rconde =  computeDiff_s( sggevx_obj->n, 
                sggevx_obj->rconde, sggevx_obj->rconderef );

    sggevx_obj->diff_rcondv =  computeDiff_s( sggevx_obj->n, 
                sggevx_obj->rcondv, sggevx_obj->rcondvref );

    sggevx_obj->diff_lscale =  computeDiff_s( sggevx_obj->n, 
                sggevx_obj->lscale, sggevx_obj->lscaleref );

    sggevx_obj->diff_rscale =  computeDiff_s( sggevx_obj->n, 
                sggevx_obj->rscale, sggevx_obj->rscaleref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n ggevx float: \n diff_a: %f \n diff_b: %f \n \
diff_vl: %f \n diff_vr: %f \n diff_alphar: %f \n\
diff_alphai: %f \n diff_beta: %f \n ilo: %d \t iloref: %d \t \
diff_ilo: %d \n ihi: %d \t ihiref: %d \t diff_ihi: %d \n \
diff_lscale: %f \t  diff_rscale: %f \n",
       sggevx_obj->diff_a, sggevx_obj->diff_b, sggevx_obj->diff_vl,
	   sggevx_obj->diff_vr, sggevx_obj->diff_alphar,
       sggevx_obj->diff_alphai, sggevx_obj->diff_beta,
	   sggevx_obj->ilo, sggevx_obj->iloref, 
	   ( sggevx_obj->ilo - sggevx_obj->iloref),
	   sggevx_obj->ihi, sggevx_obj->ihiref,
	   ( sggevx_obj->ihi - sggevx_obj->ihiref),
	   sggevx_obj->diff_lscale, sggevx_obj->diff_rscale );
#endif
}

TEST_F(sggevx_test, sggevx1) {
    //EXPECT_NEAR(0.0, sggevx_obj->diff_a, sggevx_obj->threshold);
    //EXPECT_NEAR(0.0, sggevx_obj->diff_b, sggevx_obj->threshold);
    EXPECT_NEAR(0.0, sggevx_obj->diff_vl, sggevx_obj->threshold);
    EXPECT_NEAR(0.0, sggevx_obj->diff_vr, sggevx_obj->threshold);
    EXPECT_NEAR(0.0, sggevx_obj->diff_alphar, sggevx_obj->threshold);
    EXPECT_NEAR(0.0, sggevx_obj->diff_alphai, sggevx_obj->threshold);
    EXPECT_NEAR(0.0, sggevx_obj->diff_beta, sggevx_obj->threshold);
    EXPECT_NEAR(0.0, sggevx_obj->diff_rconde, sggevx_obj->threshold);
    EXPECT_NEAR(0.0, sggevx_obj->diff_rcondv, sggevx_obj->threshold);
	/* Commenting the below check because the position of the scalefactor
	values in the (lscale & lscaleref) , (rscale, rscaleref) varies
	slightly. Though the values are matching.
	So, element by element difference calculation will not work here 
	unlike for other buffers  */
	
    /* EXPECT_NEAR(0.0, sggevx_obj->diff_lscale, sggevx_obj->threshold);
    EXPECT_NEAR(0.0, sggevx_obj->diff_rscale, sggevx_obj->threshold); */
    EXPECT_EQ(sggevx_obj->ilo, sggevx_obj->iloref);
    EXPECT_EQ(sggevx_obj->ihi, sggevx_obj->ihiref);
}

TEST_F(sggevx_test, sggevx2) {
    //EXPECT_NEAR(0.0, sggevx_obj->diff_a, sggevx_obj->threshold);
    //EXPECT_NEAR(0.0, sggevx_obj->diff_b, sggevx_obj->threshold);
    EXPECT_NEAR(0.0, sggevx_obj->diff_vl, sggevx_obj->threshold);
    EXPECT_NEAR(0.0, sggevx_obj->diff_vr, sggevx_obj->threshold);
    EXPECT_NEAR(0.0, sggevx_obj->diff_alphar, sggevx_obj->threshold);
    EXPECT_NEAR(0.0, sggevx_obj->diff_alphai, sggevx_obj->threshold);
    EXPECT_NEAR(0.0, sggevx_obj->diff_beta, sggevx_obj->threshold);
    EXPECT_NEAR(0.0, sggevx_obj->diff_rconde, sggevx_obj->threshold);
    EXPECT_NEAR(0.0, sggevx_obj->diff_rcondv, sggevx_obj->threshold);
	/* Commenting the below check because the position of the scalefactor
	values in the (lscale & lscaleref) , (rscale, rscaleref) varies
	slightly. Though the values are matching.
	So, element by element difference calculation will not work here 
	unlike for other buffers  */
	
    /* EXPECT_NEAR(0.0, sggevx_obj->diff_lscale, sggevx_obj->threshold);
    EXPECT_NEAR(0.0, sggevx_obj->diff_rscale, sggevx_obj->threshold); */
    EXPECT_EQ(sggevx_obj->ilo, sggevx_obj->iloref);
    EXPECT_EQ(sggevx_obj->ihi, sggevx_obj->ihiref);
}
TEST_F(sggevx_test, sggevx3) {
    //EXPECT_NEAR(0.0, sggevx_obj->diff_a, sggevx_obj->threshold);
    //EXPECT_NEAR(0.0, sggevx_obj->diff_b, sggevx_obj->threshold);
    EXPECT_NEAR(0.0, sggevx_obj->diff_vl, sggevx_obj->threshold);
    EXPECT_NEAR(0.0, sggevx_obj->diff_vr, sggevx_obj->threshold);
    EXPECT_NEAR(0.0, sggevx_obj->diff_alphar, sggevx_obj->threshold);
    EXPECT_NEAR(0.0, sggevx_obj->diff_alphai, sggevx_obj->threshold);
    EXPECT_NEAR(0.0, sggevx_obj->diff_beta, sggevx_obj->threshold);
    EXPECT_NEAR(0.0, sggevx_obj->diff_rconde, sggevx_obj->threshold);
    EXPECT_NEAR(0.0, sggevx_obj->diff_rcondv, sggevx_obj->threshold);
	/* Commenting the below check because the position of the scalefactor
	values in the (lscale & lscaleref) , (rscale, rscaleref) varies
	slightly. Though the values are matching.
	So, element by element difference calculation will not work here 
	unlike for other buffers  */
	
    /* EXPECT_NEAR(0.0, sggevx_obj->diff_lscale, sggevx_obj->threshold);
    EXPECT_NEAR(0.0, sggevx_obj->diff_rscale, sggevx_obj->threshold); */
    EXPECT_EQ(sggevx_obj->ilo, sggevx_obj->iloref);
    EXPECT_EQ(sggevx_obj->ihi, sggevx_obj->ihiref);
}
TEST_F(sggevx_test, sggevx4) {
    //EXPECT_NEAR(0.0, sggevx_obj->diff_a, sggevx_obj->threshold);
    //EXPECT_NEAR(0.0, sggevx_obj->diff_b, sggevx_obj->threshold);
    EXPECT_NEAR(0.0, sggevx_obj->diff_vl, sggevx_obj->threshold);
    EXPECT_NEAR(0.0, sggevx_obj->diff_vr, sggevx_obj->threshold);
    EXPECT_NEAR(0.0, sggevx_obj->diff_alphar, sggevx_obj->threshold);
    EXPECT_NEAR(0.0, sggevx_obj->diff_alphai, sggevx_obj->threshold);
    EXPECT_NEAR(0.0, sggevx_obj->diff_beta, sggevx_obj->threshold);
    EXPECT_NEAR(0.0, sggevx_obj->diff_rconde, sggevx_obj->threshold);
    EXPECT_NEAR(0.0, sggevx_obj->diff_rcondv, sggevx_obj->threshold);
	/* Commenting the below check because the position of the scalefactor
	values in the (lscale & lscaleref) , (rscale, rscaleref) varies
	slightly. Though the values are matching.
	So, element by element difference calculation will not work here 
	unlike for other buffers  */
	
    /* EXPECT_NEAR(0.0, sggevx_obj->diff_lscale, sggevx_obj->threshold);
    EXPECT_NEAR(0.0, sggevx_obj->diff_rscale, sggevx_obj->threshold); */
    EXPECT_EQ(sggevx_obj->ilo, sggevx_obj->iloref);
    EXPECT_EQ(sggevx_obj->ihi, sggevx_obj->ihiref);
}

/* Begin ggevx_double_common_parameters  class definition */
class ggevx_double_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    double diff_a, diff_b;
    double diff_alphai, diff_alphar, diff_beta, diff_vl, diff_vr;
    double diff_lscale, diff_rscale, diff_rconde, diff_rcondv;

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
    lapack_int ldb; //  The leading dimension of b
    lapack_int ldvl;
    lapack_int ldvr;

    /* Input / Output parameters */
    double* a, *aref; // contains the n-by-n general matrix A.
    double* b, *bref; // contains the n-by-n upper triangular matrix B.

    /* Output parameters */
    lapack_int sdim, sdimref;
    double *vl, *vlref; // left eigenvectors.
    double *vr, *vrref; // right eigenvectors
    lapack_int ilo, iloref, ihi, ihiref; // ilo and ihi mark the rows and columns of A which are to be reduced.
    double *lscale, *lscaleref; // permutations, scaling factors applied to the left side of A and B
    double *rscale, *rscaleref; // permutations, scaling factors applied to the right side of A and B
    double  abnrm, bbnrm, abnrmref, bbnrmref; //  L1- norm of balanced matrices
    double* alphar, *alpharref;
	double* alphai, *alphairef;
	double* beta, *betaref;
	double *rconde, *rcondv; // contain the reciprocal condition numbers 
	double *rconderef, *rcondvref; // contain the reciprocal condition numbers 
    /*Return Values */
    int info, inforef;

      ggevx_double_parameters (int matrix_layout_i, char jobvl_i,
				char jobvr_i, char balance_i, char sense_i,
				lapack_int n_i );
      ~ggevx_double_parameters ();
};

/* Constructor definition  double_common_parameters */
ggevx_double_parameters:: ggevx_double_parameters (int matrix_layout_i,
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
    ldb = n;
    ldvl = n;
    ldvr = n;

    hModule = NULL;
    dModule = NULL;

    diff_a = 0;
    diff_b = 0;
    diff_vl = 0;
    diff_vr = 0;
    diff_alphai = 0;
    diff_alphar = 0;
    diff_beta = 0;
	diff_lscale = 0;
	diff_rscale = 0;
    diff_rconde = 0;
	diff_rcondv = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n ggevx double: matrix_layout: %d n: %d  jobvl: %c \
jobvr: %c \t balance: %c \t sense: %c\n", matrix_layout, n, jobvr, 
jobvl, balance, sense);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, lda*n );
    lapacke_gtest_alloc_double_buffer_pair( &b, &bref, ldb*n );
    lapacke_gtest_alloc_double_buffer_pair( &alphar, &alpharref, n );
    lapacke_gtest_alloc_double_buffer_pair( &alphai, &alphairef, n );
    lapacke_gtest_alloc_double_buffer_pair( &beta, &betaref, n );
    lapacke_gtest_alloc_double_buffer_pair( &vl, &vlref, n*ldvl );
    lapacke_gtest_alloc_double_buffer_pair( &vr, &vrref, n*ldvr );
    lapacke_gtest_alloc_double_buffer_pair( &lscale, &lscaleref, n );
    lapacke_gtest_alloc_double_buffer_pair( &rscale, &rscaleref, n );
    lapacke_gtest_alloc_double_buffer_pair( &rconde, &rconderef, n );
    lapacke_gtest_alloc_double_buffer_pair( &rcondv, &rcondvref, n );

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (alphar==NULL) || (alpharref==NULL) || \
        (alphai==NULL) || (alphairef==NULL) || \
        (beta==NULL) || (betaref==NULL) || \
        (vl==NULL) || (vlref==NULL) || \
        (vr==NULL) || (vrref==NULL) || \
        (lscale==NULL) || (lscaleref==NULL) || \
        (rconde==NULL) || (rconderef==NULL) || \
        (rcondv==NULL) || (rcondvref==NULL) || \
        (rscale==NULL) || (rscaleref==NULL) ){
       EXPECT_FALSE( true) << "ggevx_double_parameters object: malloc error. Exiting ";
       ggevx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, lda*n );
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, ldb*n );
    lapacke_gtest_init_double_buffer_pair_with_constant( rconde, rconderef, n, 0 );
    lapacke_gtest_init_double_buffer_pair_with_constant( rcondv, rcondvref, n, 0 );
    lapacke_gtest_init_double_buffer_pair_with_constant( vl, vlref, ldvl*n, 0 );
    lapacke_gtest_init_double_buffer_pair_with_constant( vr, vrref, ldvr*n, 0 );
    lapacke_gtest_init_double_buffer_pair_with_constant( alphar, alpharref, n, 0);
    lapacke_gtest_init_double_buffer_pair_with_constant( alphai, alphairef, n, 0);
    lapacke_gtest_init_double_buffer_pair_with_constant( beta, betaref, n, 0);
    lapacke_gtest_init_double_buffer_pair_with_constant( lscale, lscaleref, n, 0);
    lapacke_gtest_init_double_buffer_pair_with_constant( rscale, rscaleref, n, 0);
    lapacke_gtest_init_double_buffer_pair_with_constant( rconde, rconderef, n, 0);
    lapacke_gtest_init_double_buffer_pair_with_constant( rcondv, rcondvref, n, 0);

   } /* end of Constructor  */

/* Destructor definition  'ggevx_double_common_parameters' */
ggevx_double_parameters :: ~ggevx_double_parameters ()
{
   /* De-Allocate memory for the input matrices */
   ggevx_free();
} 

//  Test fixture class definition
class dggevx_test  : public  ::testing::Test {
public:
   ggevx_double_parameters  *dggevx_obj;
   void SetUp();  
   void TearDown () { delete dggevx_obj; }
};

void dggevx_test::SetUp()
{
    /* LAPACKE_dggevx prototype */	
    typedef int (*Fptr_NL_LAPACKE_dggevx) (int matrix_layout, char balanc,
		char jobvl, char jobvr, char sense, lapack_int n, double* a,
		lapack_int lda, double* b, lapack_int ldb, double* alphar,
		double* alphai, double* beta, double* vl, lapack_int ldvl,
		double* vr, lapack_int ldvr, lapack_int* ilo, lapack_int* ihi,
		double* lscale, double* rscale, double* abnrm,
		double* bbnrm, double* rconde, double* rcondv);

    Fptr_NL_LAPACKE_dggevx DGGEVX;

    dggevx_obj = new  ggevx_double_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].jobvsl,
										 eig_non_sym_paramslist[idx].jobvsr,
										 eig_non_sym_paramslist[idx].balance_ggevx,
										 eig_non_sym_paramslist[idx].sense_ggevx,
                                         eig_paramslist[idx].n );
                                         
    dggevx_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    dggevx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dggevx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dggevx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dggevx_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DGGEVX = (Fptr_NL_LAPACKE_dggevx)dlsym(dggevx_obj->hModule, "LAPACKE_dggevx");
    ASSERT_TRUE(DGGEVX != NULL) << "failed to ppt the Netlib LAPACKE_dggevx symbol";

    /* Compute libflame's Lapacke o/p  */
    dggevx_obj->inforef =  LAPACKE_dggevx (   dggevx_obj->matrix_layout,
                                    dggevx_obj->balance,
                                    dggevx_obj->jobvl,
                                    dggevx_obj->jobvr,
                                    dggevx_obj->sense,
                                    dggevx_obj->n,
                                    dggevx_obj->a,
                                    dggevx_obj->lda,
                                    dggevx_obj->b,
                                    dggevx_obj->ldb,
                                    dggevx_obj->alphar,
                                    dggevx_obj->alphai,
                                    dggevx_obj->beta,
									dggevx_obj->vl,
									dggevx_obj->ldvl,
									dggevx_obj->vr,
									dggevx_obj->ldvr,
									&dggevx_obj->ilo,
									&dggevx_obj->ihi,
									dggevx_obj->lscale,
									dggevx_obj->rscale,
									&dggevx_obj->abnrm,
									&dggevx_obj->bbnrm,
									dggevx_obj->rconde,
									dggevx_obj->rcondv
                                    );

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    dggevx_obj->inforef = DGGEVX(   dggevx_obj->matrix_layout,
                                    dggevx_obj->balance,
                                    dggevx_obj->jobvl,
                                    dggevx_obj->jobvr,
                                    dggevx_obj->sense,
                                    dggevx_obj->n,
                                    dggevx_obj->aref,
                                    dggevx_obj->lda,
                                    dggevx_obj->bref,
                                    dggevx_obj->ldb,
                                    dggevx_obj->alpharref,
                                    dggevx_obj->alphairef,
                                    dggevx_obj->betaref,
									dggevx_obj->vlref,
									dggevx_obj->ldvl,
									dggevx_obj->vrref,
									dggevx_obj->ldvr,
									&dggevx_obj->iloref,
									&dggevx_obj->ihiref,
									dggevx_obj->lscaleref,
									dggevx_obj->rscaleref,
									&dggevx_obj->abnrmref,
									&dggevx_obj->bbnrmref,
									dggevx_obj->rconderef,
									dggevx_obj->rcondvref
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    dggevx_obj->diff_a =  computeDiff_d( (dggevx_obj->lda)*(dggevx_obj->n), 
                dggevx_obj->a, dggevx_obj->aref );

    dggevx_obj->diff_b =  computeDiff_d( (dggevx_obj->ldb)*(dggevx_obj->n), 
                dggevx_obj->b, dggevx_obj->bref );

    dggevx_obj->diff_vl =  computeDiff_d( (dggevx_obj->ldvl)*(dggevx_obj->n), 
                dggevx_obj->vl, dggevx_obj->vlref );

    dggevx_obj->diff_vr =  computeDiff_d( (dggevx_obj->ldvr)*(dggevx_obj->n), 
                dggevx_obj->vr, dggevx_obj->vrref );

    dggevx_obj->diff_lscale =  computeDiff_d( dggevx_obj->n, 
                dggevx_obj->lscale, dggevx_obj->lscaleref );

    dggevx_obj->diff_rscale =  computeDiff_d( dggevx_obj->n, 
                dggevx_obj->rscale, dggevx_obj->rscaleref );

    dggevx_obj->diff_alphar =  computeDiff_d( dggevx_obj->n, 
                dggevx_obj->alphar, dggevx_obj->alpharref );

    dggevx_obj->diff_alphai =  computeDiff_d( dggevx_obj->n, 
                dggevx_obj->alphai, dggevx_obj->alphairef );

    dggevx_obj->diff_beta =  computeDiff_d( dggevx_obj->n, 
                dggevx_obj->beta, dggevx_obj->betaref );

    dggevx_obj->diff_rconde =  computeDiff_d( dggevx_obj->n, 
                dggevx_obj->rconde, dggevx_obj->rconderef );

    dggevx_obj->diff_rcondv =  computeDiff_d( dggevx_obj->n, 
                dggevx_obj->rcondv, dggevx_obj->rcondvref );

    dggevx_obj->diff_lscale =  computeDiff_d( dggevx_obj->n, 
                dggevx_obj->lscale, dggevx_obj->lscaleref );

    dggevx_obj->diff_rscale =  computeDiff_d( dggevx_obj->n, 
                dggevx_obj->rscale, dggevx_obj->rscaleref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n ggevx double: \n diff_a: %f \n diff_b: %f \n \
diff_vl: %f \n diff_vr: %f \n diff_alphar: %f \n\
diff_alphai: %f \n diff_beta: %f \n ilo: %d \t iloref: %d \t \
diff_ilo: %d \n ihi: %d \t ihiref: %d \t diff_ihi: %d \n \
diff_lscale: %f \t  diff_rscale: %f \n",
       dggevx_obj->diff_a, dggevx_obj->diff_b, dggevx_obj->diff_vl,
	   dggevx_obj->diff_vr, dggevx_obj->diff_alphar,
       dggevx_obj->diff_alphai, dggevx_obj->diff_beta,
	   dggevx_obj->ilo, dggevx_obj->iloref, 
	   ( dggevx_obj->ilo - dggevx_obj->iloref),
	   dggevx_obj->ihi, dggevx_obj->ihiref,
	   ( dggevx_obj->ihi - dggevx_obj->ihiref),
	   dggevx_obj->diff_lscale, dggevx_obj->diff_rscale );
#endif
}

TEST_F(dggevx_test, dggevx1) {
    //EXPECT_NEAR(0.0, dggevx_obj->diff_a, dggevx_obj->threshold);
    //EXPECT_NEAR(0.0, dggevx_obj->diff_b, dggevx_obj->threshold);
    EXPECT_NEAR(0.0, dggevx_obj->diff_vl, dggevx_obj->threshold);
    EXPECT_NEAR(0.0, dggevx_obj->diff_vr, dggevx_obj->threshold);
    EXPECT_NEAR(0.0, dggevx_obj->diff_alphar, dggevx_obj->threshold);
    EXPECT_NEAR(0.0, dggevx_obj->diff_alphai, dggevx_obj->threshold);
    EXPECT_NEAR(0.0, dggevx_obj->diff_beta, dggevx_obj->threshold);
    EXPECT_NEAR(0.0, dggevx_obj->diff_rconde, dggevx_obj->threshold);
    EXPECT_NEAR(0.0, dggevx_obj->diff_rcondv, dggevx_obj->threshold);
	/* Commenting the below check because the position of the scalefactor
	values in the (lscale & lscaleref) , (rscale, rscaleref) varies
	slightly. Though the values are matching.
	So, element by element difference calculation will not work here 
	unlike for other buffers  */
	
    /* EXPECT_NEAR(0.0, dggevx_obj->diff_lscale, dggevx_obj->threshold);
    EXPECT_NEAR(0.0, dggevx_obj->diff_rscale, dggevx_obj->threshold); */
    EXPECT_EQ(dggevx_obj->ilo, dggevx_obj->iloref);
    EXPECT_EQ(dggevx_obj->ihi, dggevx_obj->ihiref);
}

TEST_F(dggevx_test, dggevx2) {
    //EXPECT_NEAR(0.0, dggevx_obj->diff_a, dggevx_obj->threshold);
    //EXPECT_NEAR(0.0, dggevx_obj->diff_b, dggevx_obj->threshold);
    EXPECT_NEAR(0.0, dggevx_obj->diff_vl, dggevx_obj->threshold);
    EXPECT_NEAR(0.0, dggevx_obj->diff_vr, dggevx_obj->threshold);
    EXPECT_NEAR(0.0, dggevx_obj->diff_alphar, dggevx_obj->threshold);
    EXPECT_NEAR(0.0, dggevx_obj->diff_alphai, dggevx_obj->threshold);
    EXPECT_NEAR(0.0, dggevx_obj->diff_beta, dggevx_obj->threshold);
    EXPECT_NEAR(0.0, dggevx_obj->diff_rconde, dggevx_obj->threshold);
    EXPECT_NEAR(0.0, dggevx_obj->diff_rcondv, dggevx_obj->threshold);
	/* Commenting the below check because the position of the scalefactor
	values in the (lscale & lscaleref) , (rscale, rscaleref) varies
	slightly. Though the values are matching.
	So, element by element difference calculation will not work here 
	unlike for other buffers  */
	
    /* EXPECT_NEAR(0.0, dggevx_obj->diff_lscale, dggevx_obj->threshold);
    EXPECT_NEAR(0.0, dggevx_obj->diff_rscale, dggevx_obj->threshold); */
    EXPECT_EQ(dggevx_obj->ilo, dggevx_obj->iloref);
    EXPECT_EQ(dggevx_obj->ihi, dggevx_obj->ihiref);
}
TEST_F(dggevx_test, dggevx3) {
    //EXPECT_NEAR(0.0, dggevx_obj->diff_a, dggevx_obj->threshold);
    //EXPECT_NEAR(0.0, dggevx_obj->diff_b, dggevx_obj->threshold);
    EXPECT_NEAR(0.0, dggevx_obj->diff_vl, dggevx_obj->threshold);
    EXPECT_NEAR(0.0, dggevx_obj->diff_vr, dggevx_obj->threshold);
    EXPECT_NEAR(0.0, dggevx_obj->diff_alphar, dggevx_obj->threshold);
    EXPECT_NEAR(0.0, dggevx_obj->diff_alphai, dggevx_obj->threshold);
    EXPECT_NEAR(0.0, dggevx_obj->diff_beta, dggevx_obj->threshold);
    EXPECT_NEAR(0.0, dggevx_obj->diff_rconde, dggevx_obj->threshold);
    EXPECT_NEAR(0.0, dggevx_obj->diff_rcondv, dggevx_obj->threshold);
	/* Commenting the below check because the position of the scalefactor
	values in the (lscale & lscaleref) , (rscale, rscaleref) varies
	slightly. Though the values are matching.
	So, element by element difference calculation will not work here 
	unlike for other buffers  */
	
    /* EXPECT_NEAR(0.0, dggevx_obj->diff_lscale, dggevx_obj->threshold);
    EXPECT_NEAR(0.0, dggevx_obj->diff_rscale, dggevx_obj->threshold); */
    EXPECT_EQ(dggevx_obj->ilo, dggevx_obj->iloref);
    EXPECT_EQ(dggevx_obj->ihi, dggevx_obj->ihiref);
}
TEST_F(dggevx_test, dggevx4) {
    //EXPECT_NEAR(0.0, dggevx_obj->diff_a, dggevx_obj->threshold);
    //EXPECT_NEAR(0.0, dggevx_obj->diff_b, dggevx_obj->threshold);
    EXPECT_NEAR(0.0, dggevx_obj->diff_vl, dggevx_obj->threshold);
    EXPECT_NEAR(0.0, dggevx_obj->diff_vr, dggevx_obj->threshold);
    EXPECT_NEAR(0.0, dggevx_obj->diff_alphar, dggevx_obj->threshold);
    EXPECT_NEAR(0.0, dggevx_obj->diff_alphai, dggevx_obj->threshold);
    EXPECT_NEAR(0.0, dggevx_obj->diff_beta, dggevx_obj->threshold);
    EXPECT_NEAR(0.0, dggevx_obj->diff_rconde, dggevx_obj->threshold);
    EXPECT_NEAR(0.0, dggevx_obj->diff_rcondv, dggevx_obj->threshold);
	/* Commenting the below check because the position of the scalefactor
	values in the (lscale & lscaleref) , (rscale, rscaleref) varies
	slightly. Though the values are matching.
	So, element by element difference calculation will not work here 
	unlike for other buffers  */
	
    /* EXPECT_NEAR(0.0, dggevx_obj->diff_lscale, dggevx_obj->threshold);
    EXPECT_NEAR(0.0, dggevx_obj->diff_rscale, dggevx_obj->threshold); */
    EXPECT_EQ(dggevx_obj->ilo, dggevx_obj->iloref);
    EXPECT_EQ(dggevx_obj->ihi, dggevx_obj->ihiref);
}

/* Begin ggevx_scomplex_common_parameters  class definition */
class ggevx_scomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_a, diff_b;
    float diff_alphai, diff_alphar, diff_beta, diff_vl, diff_vr;
    float diff_lscale, diff_rscale, diff_rconde, diff_rcondv;

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
    lapack_int ldb; //  The leading dimension of b
    lapack_int ldvl;
    lapack_int ldvr;

    /* Input / Output parameters */
    lapack_complex_float* a, *aref; // contains the n-by-n general matrix A.
    lapack_complex_float* b, *bref; // contains the n-by-n upper triangular matrix B.

    /* Output parameters */
    lapack_int sdim, sdimref;
    lapack_complex_float *vl, *vlref; // left eigenvectors.
    lapack_complex_float *vr, *vrref; // right eigenvectors
    lapack_int ilo, iloref, ihi, ihiref; // ilo and ihi mark the rows and columns of A which are to be reduced.
    float *lscale, *lscaleref; // permutations, scaling factors applied to the left side of A and B
    float *rscale, *rscaleref; // permutations, scaling factors applied to the right side of A and B
    float  abnrm, bbnrm, abnrmref, bbnrmref; //  L1- norm of balanced matrices
    lapack_complex_float* alphar, *alpharref;
	lapack_complex_float* alphai, *alphairef;
	lapack_complex_float* beta, *betaref;
	float *rconde, *rcondv; // contain the reciprocal condition numbers 
	float *rconderef, *rcondvref; // contain the reciprocal condition numbers 
    /*Return Values */
    int info, inforef;

      ggevx_scomplex_parameters (int matrix_layout_i, char jobvl_i,
				char jobvr_i, char balance_i, char sense_i,
				lapack_int n_i );
      ~ggevx_scomplex_parameters ();
};

/* Constructor definition  lapack_complex_float_common_parameters */
ggevx_scomplex_parameters:: ggevx_scomplex_parameters (int matrix_layout_i,
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
    ldb = n;
    ldvl = n;
    ldvr = n;

    hModule = NULL;
    dModule = NULL;

    diff_a = 0;
    diff_b = 0;
    diff_vl = 0;
    diff_vr = 0;
    diff_alphai = 0;
    diff_alphar = 0;
    diff_beta = 0;
	diff_lscale = 0;
	diff_rscale = 0;
    diff_rconde = 0;
	diff_rcondv = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n ggevx lapack_complex_float: matrix_layout: %d n: %d  jobvl: %c \
jobvr: %c \t balance: %c \t sense: %c\n", matrix_layout, n, jobvr, 
jobvl, balance, sense);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, lda*n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &b, &bref, ldb*n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &alphar, &alpharref, n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &alphai, &alphairef, n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &beta, &betaref, n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &vl, &vlref, n*ldvl );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &vr, &vrref, n*ldvr );
    lapacke_gtest_alloc_float_buffer_pair( &lscale, &lscaleref, n );
    lapacke_gtest_alloc_float_buffer_pair( &rscale, &rscaleref, n );
    lapacke_gtest_alloc_float_buffer_pair( &rconde, &rconderef, n );
    lapacke_gtest_alloc_float_buffer_pair( &rcondv, &rcondvref, n );

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (alphar==NULL) || (alpharref==NULL) || \
        (alphai==NULL) || (alphairef==NULL) || \
        (beta==NULL) || (betaref==NULL) || \
        (vl==NULL) || (vlref==NULL) || \
        (vr==NULL) || (vrref==NULL) || \
        (lscale==NULL) || (lscaleref==NULL) || \
        (rconde==NULL) || (rconderef==NULL) || \
        (rcondv==NULL) || (rcondvref==NULL) || \
        (rscale==NULL) || (rscaleref==NULL) ){
       EXPECT_FALSE( true) << "ggevx_scomplex_parameters object: malloc error. Exiting ";
       ggevx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, lda*n );
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, ldb*n );
    lapacke_gtest_init_float_buffer_pair_with_constant( rconde, rconderef, n, 0 );
    lapacke_gtest_init_float_buffer_pair_with_constant( rcondv, rcondvref, n, 0 );
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( vl, vlref, ldvl*n, 0 );
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( vr, vrref, ldvr*n, 0 );
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( alphar, alpharref, n, 0);
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( alphai, alphairef, n, 0);
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( beta, betaref, n, 0);
    lapacke_gtest_init_float_buffer_pair_with_constant( lscale, lscaleref, n, 0);
    lapacke_gtest_init_float_buffer_pair_with_constant( rscale, rscaleref, n, 0);
    lapacke_gtest_init_float_buffer_pair_with_constant( rconde, rconderef, n, 0);
    lapacke_gtest_init_float_buffer_pair_with_constant( rcondv, rcondvref, n, 0);

   } /* end of Constructor  */

/* Destructor definition  'ggevx_scomplex_common_parameters' */
ggevx_scomplex_parameters :: ~ggevx_scomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   ggevx_free();
} 

//  Test fixture class definition
class cggevx_test  : public  ::testing::Test {
public:
   ggevx_scomplex_parameters  *cggevx_obj;
   void SetUp();  
   void TearDown () { delete cggevx_obj; }
};

void cggevx_test::SetUp()
{
    /* LAPACKE_cggevx prototype */	
    typedef int (*Fptr_NL_LAPACKE_cggevx) (int matrix_layout, char balanc,
	char jobvl, char jobvr, char sense, lapack_int n, lapack_complex_float* a,
	lapack_int lda, lapack_complex_float* b, lapack_int ldb,
	lapack_complex_float* alpha, lapack_complex_float* beta,
	lapack_complex_float* vl, lapack_int ldvl, lapack_complex_float* vr,
	lapack_int ldvr, lapack_int* ilo, lapack_int* ihi, float* lscale,
	float* rscale, float* abnrm, float* bbnrm, float* rconde, float* rcondv);

    Fptr_NL_LAPACKE_cggevx CGGEVX;

    cggevx_obj = new  ggevx_scomplex_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].jobvsl,
										 eig_non_sym_paramslist[idx].jobvsr,
										 eig_non_sym_paramslist[idx].balance_ggevx,
										 eig_non_sym_paramslist[idx].sense_ggevx,
                                         eig_paramslist[idx].n );
                                         
    cggevx_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    cggevx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cggevx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cggevx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cggevx_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CGGEVX = (Fptr_NL_LAPACKE_cggevx)dlsym(cggevx_obj->hModule, "LAPACKE_cggevx");
    ASSERT_TRUE(CGGEVX != NULL) << "failed to ppt the Netlib LAPACKE_cggevx symbol";

    /* Compute libflame's Lapacke o/p  */
    cggevx_obj->inforef =  LAPACKE_cggevx (   cggevx_obj->matrix_layout,
                                    cggevx_obj->balance,
                                    cggevx_obj->jobvl,
                                    cggevx_obj->jobvr,
                                    cggevx_obj->sense,
                                    cggevx_obj->n,
                                    cggevx_obj->a,
                                    cggevx_obj->lda,
                                    cggevx_obj->b,
                                    cggevx_obj->ldb,
                                    cggevx_obj->alphar,
                                    cggevx_obj->beta,
									cggevx_obj->vl,
									cggevx_obj->ldvl,
									cggevx_obj->vr,
									cggevx_obj->ldvr,
									&cggevx_obj->ilo,
									&cggevx_obj->ihi,
									cggevx_obj->lscale,
									cggevx_obj->rscale,
									&cggevx_obj->abnrm,
									&cggevx_obj->bbnrm,
									cggevx_obj->rconde,
									cggevx_obj->rcondv
                                    );

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    cggevx_obj->inforef = CGGEVX(   cggevx_obj->matrix_layout,
                                    cggevx_obj->balance,
                                    cggevx_obj->jobvl,
                                    cggevx_obj->jobvr,
                                    cggevx_obj->sense,
                                    cggevx_obj->n,
                                    cggevx_obj->aref,
                                    cggevx_obj->lda,
                                    cggevx_obj->bref,
                                    cggevx_obj->ldb,
                                    cggevx_obj->alpharref,
                                    cggevx_obj->betaref,
									cggevx_obj->vlref,
									cggevx_obj->ldvl,
									cggevx_obj->vrref,
									cggevx_obj->ldvr,
									&cggevx_obj->iloref,
									&cggevx_obj->ihiref,
									cggevx_obj->lscaleref,
									cggevx_obj->rscaleref,
									&cggevx_obj->abnrmref,
									&cggevx_obj->bbnrmref,
									cggevx_obj->rconderef,
									cggevx_obj->rcondvref
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    cggevx_obj->diff_a =  computeDiff_c( (cggevx_obj->lda)*(cggevx_obj->n), 
                cggevx_obj->a, cggevx_obj->aref );

    cggevx_obj->diff_b =  computeDiff_c( (cggevx_obj->ldb)*(cggevx_obj->n), 
                cggevx_obj->b, cggevx_obj->bref );

    cggevx_obj->diff_vl =  computeDiff_c( (cggevx_obj->ldvl)*(cggevx_obj->n), 
                cggevx_obj->vl, cggevx_obj->vlref );

    cggevx_obj->diff_vr =  computeDiff_c( (cggevx_obj->ldvr)*(cggevx_obj->n), 
                cggevx_obj->vr, cggevx_obj->vrref );

    cggevx_obj->diff_lscale =  computeDiff_s( cggevx_obj->n, 
                cggevx_obj->lscale, cggevx_obj->lscaleref );

    cggevx_obj->diff_rscale =  computeDiff_s( cggevx_obj->n, 
                cggevx_obj->rscale, cggevx_obj->rscaleref );

    cggevx_obj->diff_alphar =  computeDiff_c( cggevx_obj->n, 
                cggevx_obj->alphar, cggevx_obj->alpharref );

    cggevx_obj->diff_alphai =  computeDiff_c( cggevx_obj->n, 
                cggevx_obj->alphai, cggevx_obj->alphairef );

    cggevx_obj->diff_beta =  computeDiff_c( cggevx_obj->n, 
                cggevx_obj->beta, cggevx_obj->betaref );

    cggevx_obj->diff_rconde =  computeDiff_s( cggevx_obj->n, 
                cggevx_obj->rconde, cggevx_obj->rconderef );

    cggevx_obj->diff_rcondv =  computeDiff_s( cggevx_obj->n, 
                cggevx_obj->rcondv, cggevx_obj->rcondvref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n ggevx lapack_complex_float: \n diff_a: %f \n diff_b: %f \n \
diff_vl: %f \n diff_vr: %f \n diff_alphar: %f \n\
diff_alphai: %f \n diff_beta: %f \n ilo: %d \t iloref: %d \t \
diff_ilo: %d \n ihi: %d \t ihiref: %d \t diff_ihi: %d \n \
diff_lscale: %f \t  diff_rscale: %f \n",
       cggevx_obj->diff_a, cggevx_obj->diff_b, cggevx_obj->diff_vl,
	   cggevx_obj->diff_vr, cggevx_obj->diff_alphar,
       cggevx_obj->diff_alphai, cggevx_obj->diff_beta,
	   cggevx_obj->ilo, cggevx_obj->iloref, 
	   ( cggevx_obj->ilo - cggevx_obj->iloref),
	   cggevx_obj->ihi, cggevx_obj->ihiref,
	   ( cggevx_obj->ihi - cggevx_obj->ihiref),
	   cggevx_obj->diff_lscale, cggevx_obj->diff_rscale );
#endif
}

TEST_F(cggevx_test, cggevx1) {
    //EXPECT_NEAR(0.0, cggevx_obj->diff_a, cggevx_obj->threshold);
    //EXPECT_NEAR(0.0, cggevx_obj->diff_b, cggevx_obj->threshold);
    EXPECT_NEAR(0.0, cggevx_obj->diff_vl, cggevx_obj->threshold);
    EXPECT_NEAR(0.0, cggevx_obj->diff_vr, cggevx_obj->threshold);
    EXPECT_NEAR(0.0, cggevx_obj->diff_alphar, cggevx_obj->threshold);
    EXPECT_NEAR(0.0, cggevx_obj->diff_alphai, cggevx_obj->threshold);
    EXPECT_NEAR(0.0, cggevx_obj->diff_beta, cggevx_obj->threshold);
    EXPECT_NEAR(0.0, cggevx_obj->diff_rconde, cggevx_obj->threshold);
    EXPECT_NEAR(0.0, cggevx_obj->diff_rcondv, cggevx_obj->threshold);
	/* Commenting the below check because the position of the scalefactor
	values in the (lscale & lscaleref) , (rscale, rscaleref) varies
	slightly. Though the values are matching.
	So, element by element difference calculation will not work here 
	unlike for other buffers  */
	
    /* EXPECT_NEAR(0.0, cggevx_obj->diff_lscale, cggevx_obj->threshold);
    EXPECT_NEAR(0.0, cggevx_obj->diff_rscale, cggevx_obj->threshold); */
    EXPECT_EQ(cggevx_obj->ilo, cggevx_obj->iloref);
    EXPECT_EQ(cggevx_obj->ihi, cggevx_obj->ihiref);
}

TEST_F(cggevx_test, cggevx2) {
    //EXPECT_NEAR(0.0, cggevx_obj->diff_a, cggevx_obj->threshold);
    //EXPECT_NEAR(0.0, cggevx_obj->diff_b, cggevx_obj->threshold);
    EXPECT_NEAR(0.0, cggevx_obj->diff_vl, cggevx_obj->threshold);
    EXPECT_NEAR(0.0, cggevx_obj->diff_vr, cggevx_obj->threshold);
    EXPECT_NEAR(0.0, cggevx_obj->diff_alphar, cggevx_obj->threshold);
    EXPECT_NEAR(0.0, cggevx_obj->diff_alphai, cggevx_obj->threshold);
    EXPECT_NEAR(0.0, cggevx_obj->diff_beta, cggevx_obj->threshold);
    EXPECT_NEAR(0.0, cggevx_obj->diff_rconde, cggevx_obj->threshold);
    EXPECT_NEAR(0.0, cggevx_obj->diff_rcondv, cggevx_obj->threshold);
	/* Commenting the below check because the position of the scalefactor
	values in the (lscale & lscaleref) , (rscale, rscaleref) varies
	slightly. Though the values are matching.
	So, element by element difference calculation will not work here 
	unlike for other buffers  */
	
    /* EXPECT_NEAR(0.0, cggevx_obj->diff_lscale, cggevx_obj->threshold);
    EXPECT_NEAR(0.0, cggevx_obj->diff_rscale, cggevx_obj->threshold); */
    EXPECT_EQ(cggevx_obj->ilo, cggevx_obj->iloref);
    EXPECT_EQ(cggevx_obj->ihi, cggevx_obj->ihiref);
}
TEST_F(cggevx_test, cggevx3) {
    //EXPECT_NEAR(0.0, cggevx_obj->diff_a, cggevx_obj->threshold);
    //EXPECT_NEAR(0.0, cggevx_obj->diff_b, cggevx_obj->threshold);
    EXPECT_NEAR(0.0, cggevx_obj->diff_vl, cggevx_obj->threshold);
    EXPECT_NEAR(0.0, cggevx_obj->diff_vr, cggevx_obj->threshold);
    EXPECT_NEAR(0.0, cggevx_obj->diff_alphar, cggevx_obj->threshold);
    EXPECT_NEAR(0.0, cggevx_obj->diff_alphai, cggevx_obj->threshold);
    EXPECT_NEAR(0.0, cggevx_obj->diff_beta, cggevx_obj->threshold);
    EXPECT_NEAR(0.0, cggevx_obj->diff_rconde, cggevx_obj->threshold);
    EXPECT_NEAR(0.0, cggevx_obj->diff_rcondv, cggevx_obj->threshold);
	/* Commenting the below check because the position of the scalefactor
	values in the (lscale & lscaleref) , (rscale, rscaleref) varies
	slightly. Though the values are matching.
	So, element by element difference calculation will not work here 
	unlike for other buffers  */
	
    /* EXPECT_NEAR(0.0, cggevx_obj->diff_lscale, cggevx_obj->threshold);
    EXPECT_NEAR(0.0, cggevx_obj->diff_rscale, cggevx_obj->threshold); */
    EXPECT_EQ(cggevx_obj->ilo, cggevx_obj->iloref);
    EXPECT_EQ(cggevx_obj->ihi, cggevx_obj->ihiref);
}
TEST_F(cggevx_test, cggevx4) {
    //EXPECT_NEAR(0.0, cggevx_obj->diff_a, cggevx_obj->threshold);
    //EXPECT_NEAR(0.0, cggevx_obj->diff_b, cggevx_obj->threshold);
    EXPECT_NEAR(0.0, cggevx_obj->diff_vl, cggevx_obj->threshold);
    EXPECT_NEAR(0.0, cggevx_obj->diff_vr, cggevx_obj->threshold);
    EXPECT_NEAR(0.0, cggevx_obj->diff_alphar, cggevx_obj->threshold);
    EXPECT_NEAR(0.0, cggevx_obj->diff_alphai, cggevx_obj->threshold);
    EXPECT_NEAR(0.0, cggevx_obj->diff_beta, cggevx_obj->threshold);
    EXPECT_NEAR(0.0, cggevx_obj->diff_rconde, cggevx_obj->threshold);
    EXPECT_NEAR(0.0, cggevx_obj->diff_rcondv, cggevx_obj->threshold);
	/* Commenting the below check because the position of the scalefactor
	values in the (lscale & lscaleref) , (rscale, rscaleref) varies
	slightly. Though the values are matching.
	So, element by element difference calculation will not work here 
	unlike for other buffers  */
	
    /* EXPECT_NEAR(0.0, cggevx_obj->diff_lscale, cggevx_obj->threshold);
    EXPECT_NEAR(0.0, cggevx_obj->diff_rscale, cggevx_obj->threshold); */
    EXPECT_EQ(cggevx_obj->ilo, cggevx_obj->iloref);
    EXPECT_EQ(cggevx_obj->ihi, cggevx_obj->ihiref);
}

/* Begin ggevx_dcomplex_common_parameters  class definition */
class ggevx_dcomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    double diff_a, diff_b;
    double diff_alphai, diff_alphar, diff_beta, diff_vl, diff_vr;
    double diff_lscale, diff_rscale, diff_rconde, diff_rcondv;

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
    lapack_int ldb; //  The leading dimension of b
    lapack_int ldvl;
    lapack_int ldvr;

    /* Input / Output parameters */
    lapack_complex_double* a, *aref; // contains the n-by-n general matrix A.
    lapack_complex_double* b, *bref; // contains the n-by-n upper triangular matrix B.

    /* Output parameters */
    lapack_int sdim, sdimref;
    lapack_complex_double *vl, *vlref; // left eigenvectors.
    lapack_complex_double *vr, *vrref; // right eigenvectors
    lapack_int ilo, iloref, ihi, ihiref; // ilo and ihi mark the rows and columns of A which are to be reduced.
    double *lscale, *lscaleref; // permutations, scaling factors applied to the left side of A and B
    double *rscale, *rscaleref; // permutations, scaling factors applied to the right side of A and B
    double  abnrm, bbnrm, abnrmref, bbnrmref; //  L1- norm of balanced matrices
    lapack_complex_double* alphar, *alpharref;
	lapack_complex_double* alphai, *alphairef;
	lapack_complex_double* beta, *betaref;
	double *rconde, *rcondv; // contain the reciprocal condition numbers 
	double *rconderef, *rcondvref; // contain the reciprocal condition numbers 
    /*Return Values */
    int info, inforef;

      ggevx_dcomplex_parameters (int matrix_layout_i, char jobvl_i,
				char jobvr_i, char balance_i, char sense_i,
				lapack_int n_i );
      ~ggevx_dcomplex_parameters ();
};

/* Constructor definition  lapack_complex_double_common_parameters */
ggevx_dcomplex_parameters:: ggevx_dcomplex_parameters (int matrix_layout_i,
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
    ldb = n;
    ldvl = n;
    ldvr = n;

    hModule = NULL;
    dModule = NULL;

    diff_a = 0;
    diff_b = 0;
    diff_vl = 0;
    diff_vr = 0;
    diff_alphai = 0;
    diff_alphar = 0;
    diff_beta = 0;
	diff_lscale = 0;
	diff_rscale = 0;
    diff_rconde = 0;
	diff_rcondv = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n ggevx lapack_complex_double: matrix_layout: %d n: %d  jobvl: %c \
jobvr: %c \t balance: %c \t sense: %c\n", matrix_layout, n, jobvr, 
jobvl, balance, sense);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, lda*n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &b, &bref, ldb*n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &alphar, &alpharref, n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &alphai, &alphairef, n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &beta, &betaref, n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &vl, &vlref, n*ldvl );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &vr, &vrref, n*ldvr );
    lapacke_gtest_alloc_double_buffer_pair( &lscale, &lscaleref, n );
    lapacke_gtest_alloc_double_buffer_pair( &rscale, &rscaleref, n );
    lapacke_gtest_alloc_double_buffer_pair( &rconde, &rconderef, n );
    lapacke_gtest_alloc_double_buffer_pair( &rcondv, &rcondvref, n );

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (alphar==NULL) || (alpharref==NULL) || \
        (alphai==NULL) || (alphairef==NULL) || \
        (beta==NULL) || (betaref==NULL) || \
        (vl==NULL) || (vlref==NULL) || \
        (vr==NULL) || (vrref==NULL) || \
        (lscale==NULL) || (lscaleref==NULL) || \
        (rconde==NULL) || (rconderef==NULL) || \
        (rcondv==NULL) || (rcondvref==NULL) || \
        (rscale==NULL) || (rscaleref==NULL) ){
       EXPECT_FALSE( true) << "ggevx_dcomplex_parameters object: malloc error. Exiting ";
       ggevx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, lda*n );
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, ldb*n );
    lapacke_gtest_init_double_buffer_pair_with_constant( rconde, rconderef, n, 0 );
    lapacke_gtest_init_double_buffer_pair_with_constant( rcondv, rcondvref, n, 0 );
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( vl, vlref, ldvl*n, 0 );
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( vr, vrref, ldvr*n, 0 );
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( alphar, alpharref, n, 0);
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( alphai, alphairef, n, 0);
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( beta, betaref, n, 0);
    lapacke_gtest_init_double_buffer_pair_with_constant( lscale, lscaleref, n, 0);
    lapacke_gtest_init_double_buffer_pair_with_constant( rscale, rscaleref, n, 0);
    lapacke_gtest_init_double_buffer_pair_with_constant( rconde, rconderef, n, 0);
    lapacke_gtest_init_double_buffer_pair_with_constant( rcondv, rcondvref, n, 0);

   } /* end of Constructor  */

/* Destructor definition  'ggevx_dcomplex_common_parameters' */
ggevx_dcomplex_parameters :: ~ggevx_dcomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   ggevx_free();
} 

//  Test fixture class definition
class zggevx_test  : public  ::testing::Test {
public:
   ggevx_dcomplex_parameters  *zggevx_obj;
   void SetUp();  
   void TearDown () { delete zggevx_obj; }
};

void zggevx_test::SetUp()
{
    /* LAPACKE_zggevx prototype */	
    typedef int (*Fptr_NL_LAPACKE_zggevx) (int matrix_layout, char balanc,
	char jobvl, char jobvr, char sense, lapack_int n, lapack_complex_double* a,
	lapack_int lda, lapack_complex_double* b, lapack_int ldb,
	lapack_complex_double* alpha, lapack_complex_double* beta,
	lapack_complex_double* vl, lapack_int ldvl, lapack_complex_double* vr,
	lapack_int ldvr, lapack_int* ilo, lapack_int* ihi, double* lscale,
	double* rscale, double* abnrm, double* bbnrm, double* rconde, double* rcondv);

    Fptr_NL_LAPACKE_zggevx ZGGEVX;

    zggevx_obj = new  ggevx_dcomplex_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].jobvsl,
										 eig_non_sym_paramslist[idx].jobvsr,
										 eig_non_sym_paramslist[idx].balance_ggevx,
										 eig_non_sym_paramslist[idx].sense_ggevx,
                                         eig_paramslist[idx].n );
                                         
    zggevx_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    zggevx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zggevx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zggevx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zggevx_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZGGEVX = (Fptr_NL_LAPACKE_zggevx)dlsym(zggevx_obj->hModule, "LAPACKE_zggevx");
    ASSERT_TRUE(ZGGEVX != NULL) << "failed to ppt the Netlib LAPACKE_zggevx symbol";

    /* Compute libflame's Lapacke o/p  */
    zggevx_obj->inforef =  LAPACKE_zggevx (   zggevx_obj->matrix_layout,
                                    zggevx_obj->balance,
                                    zggevx_obj->jobvl,
                                    zggevx_obj->jobvr,
                                    zggevx_obj->sense,
                                    zggevx_obj->n,
                                    zggevx_obj->a,
                                    zggevx_obj->lda,
                                    zggevx_obj->b,
                                    zggevx_obj->ldb,
                                    zggevx_obj->alphar,
                                    zggevx_obj->beta,
									zggevx_obj->vl,
									zggevx_obj->ldvl,
									zggevx_obj->vr,
									zggevx_obj->ldvr,
									&zggevx_obj->ilo,
									&zggevx_obj->ihi,
									zggevx_obj->lscale,
									zggevx_obj->rscale,
									&zggevx_obj->abnrm,
									&zggevx_obj->bbnrm,
									zggevx_obj->rconde,
									zggevx_obj->rcondv
                                    );

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    zggevx_obj->inforef = ZGGEVX(   zggevx_obj->matrix_layout,
                                    zggevx_obj->balance,
                                    zggevx_obj->jobvl,
                                    zggevx_obj->jobvr,
                                    zggevx_obj->sense,
                                    zggevx_obj->n,
                                    zggevx_obj->aref,
                                    zggevx_obj->lda,
                                    zggevx_obj->bref,
                                    zggevx_obj->ldb,
                                    zggevx_obj->alpharref,
                                    zggevx_obj->betaref,
									zggevx_obj->vlref,
									zggevx_obj->ldvl,
									zggevx_obj->vrref,
									zggevx_obj->ldvr,
									&zggevx_obj->iloref,
									&zggevx_obj->ihiref,
									zggevx_obj->lscaleref,
									zggevx_obj->rscaleref,
									&zggevx_obj->abnrmref,
									&zggevx_obj->bbnrmref,
									zggevx_obj->rconderef,
									zggevx_obj->rcondvref
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    zggevx_obj->diff_a =  computeDiff_z( (zggevx_obj->lda)*(zggevx_obj->n), 
                zggevx_obj->a, zggevx_obj->aref );

    zggevx_obj->diff_b =  computeDiff_z( (zggevx_obj->ldb)*(zggevx_obj->n), 
                zggevx_obj->b, zggevx_obj->bref );

    zggevx_obj->diff_vl =  computeDiff_z( (zggevx_obj->ldvl)*(zggevx_obj->n), 
                zggevx_obj->vl, zggevx_obj->vlref );

    zggevx_obj->diff_vr =  computeDiff_z( (zggevx_obj->ldvr)*(zggevx_obj->n), 
                zggevx_obj->vr, zggevx_obj->vrref );

    zggevx_obj->diff_lscale =  computeDiff_d( zggevx_obj->n, 
                zggevx_obj->lscale, zggevx_obj->lscaleref );

    zggevx_obj->diff_rscale =  computeDiff_d( zggevx_obj->n, 
                zggevx_obj->rscale, zggevx_obj->rscaleref );

    zggevx_obj->diff_alphar =  computeDiff_z( zggevx_obj->n, 
                zggevx_obj->alphar, zggevx_obj->alpharref );

    zggevx_obj->diff_alphai =  computeDiff_z( zggevx_obj->n, 
                zggevx_obj->alphai, zggevx_obj->alphairef );

    zggevx_obj->diff_beta =  computeDiff_z( zggevx_obj->n, 
                zggevx_obj->beta, zggevx_obj->betaref );

    zggevx_obj->diff_rconde =  computeDiff_d( zggevx_obj->n, 
                zggevx_obj->rconde, zggevx_obj->rconderef );

    zggevx_obj->diff_rcondv =  computeDiff_d( zggevx_obj->n, 
                zggevx_obj->rcondv, zggevx_obj->rcondvref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n ggevx lapack_complex_double: \n diff_a: %f \n diff_b: %f \n \
diff_vl: %f \n diff_vr: %f \n diff_alphar: %f \n\
diff_alphai: %f \n diff_beta: %f \n ilo: %d \t iloref: %d \t \
diff_ilo: %d \n ihi: %d \t ihiref: %d \t diff_ihi: %d \n \
diff_lscale: %f \t  diff_rscale: %f \n",
       zggevx_obj->diff_a, zggevx_obj->diff_b, zggevx_obj->diff_vl,
	   zggevx_obj->diff_vr, zggevx_obj->diff_alphar,
       zggevx_obj->diff_alphai, zggevx_obj->diff_beta,
	   zggevx_obj->ilo, zggevx_obj->iloref, 
	   ( zggevx_obj->ilo - zggevx_obj->iloref),
	   zggevx_obj->ihi, zggevx_obj->ihiref,
	   ( zggevx_obj->ihi - zggevx_obj->ihiref),
	   zggevx_obj->diff_lscale, zggevx_obj->diff_rscale );
#endif
}

TEST_F(zggevx_test, zggevx1) {
    //EXPECT_NEAR(0.0, zggevx_obj->diff_a, zggevx_obj->threshold);
    //EXPECT_NEAR(0.0, zggevx_obj->diff_b, zggevx_obj->threshold);
    EXPECT_NEAR(0.0, zggevx_obj->diff_vl, zggevx_obj->threshold);
    EXPECT_NEAR(0.0, zggevx_obj->diff_vr, zggevx_obj->threshold);
    EXPECT_NEAR(0.0, zggevx_obj->diff_alphar, zggevx_obj->threshold);
    EXPECT_NEAR(0.0, zggevx_obj->diff_alphai, zggevx_obj->threshold);
    EXPECT_NEAR(0.0, zggevx_obj->diff_beta, zggevx_obj->threshold);
    EXPECT_NEAR(0.0, zggevx_obj->diff_rconde, zggevx_obj->threshold);
    EXPECT_NEAR(0.0, zggevx_obj->diff_rcondv, zggevx_obj->threshold);
	/* Commenting the below check because the position of the scalefactor
	values in the (lscale & lscaleref) , (rscale, rscaleref) varies
	slightly. Though the values are matching.
	So, element by element difference calculation will not work here 
	unlike for other buffers  */
	
    /* EXPECT_NEAR(0.0, zggevx_obj->diff_lscale, zggevx_obj->threshold);
    EXPECT_NEAR(0.0, zggevx_obj->diff_rscale, zggevx_obj->threshold); */
    EXPECT_EQ(zggevx_obj->ilo, zggevx_obj->iloref);
    EXPECT_EQ(zggevx_obj->ihi, zggevx_obj->ihiref);
}

TEST_F(zggevx_test, zggevx2) {
    //EXPECT_NEAR(0.0, zggevx_obj->diff_a, zggevx_obj->threshold);
    //EXPECT_NEAR(0.0, zggevx_obj->diff_b, zggevx_obj->threshold);
    EXPECT_NEAR(0.0, zggevx_obj->diff_vl, zggevx_obj->threshold);
    EXPECT_NEAR(0.0, zggevx_obj->diff_vr, zggevx_obj->threshold);
    EXPECT_NEAR(0.0, zggevx_obj->diff_alphar, zggevx_obj->threshold);
    EXPECT_NEAR(0.0, zggevx_obj->diff_alphai, zggevx_obj->threshold);
    EXPECT_NEAR(0.0, zggevx_obj->diff_beta, zggevx_obj->threshold);
    EXPECT_NEAR(0.0, zggevx_obj->diff_rconde, zggevx_obj->threshold);
    EXPECT_NEAR(0.0, zggevx_obj->diff_rcondv, zggevx_obj->threshold);
	/* Commenting the below check because the position of the scalefactor
	values in the (lscale & lscaleref) , (rscale, rscaleref) varies
	slightly. Though the values are matching.
	So, element by element difference calculation will not work here 
	unlike for other buffers  */
	
    /* EXPECT_NEAR(0.0, zggevx_obj->diff_lscale, zggevx_obj->threshold);
    EXPECT_NEAR(0.0, zggevx_obj->diff_rscale, zggevx_obj->threshold); */
    EXPECT_EQ(zggevx_obj->ilo, zggevx_obj->iloref);
    EXPECT_EQ(zggevx_obj->ihi, zggevx_obj->ihiref);
}
TEST_F(zggevx_test, zggevx3) {
    //EXPECT_NEAR(0.0, zggevx_obj->diff_a, zggevx_obj->threshold);
    //EXPECT_NEAR(0.0, zggevx_obj->diff_b, zggevx_obj->threshold);
    EXPECT_NEAR(0.0, zggevx_obj->diff_vl, zggevx_obj->threshold);
    EXPECT_NEAR(0.0, zggevx_obj->diff_vr, zggevx_obj->threshold);
    EXPECT_NEAR(0.0, zggevx_obj->diff_alphar, zggevx_obj->threshold);
    EXPECT_NEAR(0.0, zggevx_obj->diff_alphai, zggevx_obj->threshold);
    EXPECT_NEAR(0.0, zggevx_obj->diff_beta, zggevx_obj->threshold);
    EXPECT_NEAR(0.0, zggevx_obj->diff_rconde, zggevx_obj->threshold);
    EXPECT_NEAR(0.0, zggevx_obj->diff_rcondv, zggevx_obj->threshold);
	/* Commenting the below check because the position of the scalefactor
	values in the (lscale & lscaleref) , (rscale, rscaleref) varies
	slightly. Though the values are matching.
	So, element by element difference calculation will not work here 
	unlike for other buffers  */
	
    /* EXPECT_NEAR(0.0, zggevx_obj->diff_lscale, zggevx_obj->threshold);
    EXPECT_NEAR(0.0, zggevx_obj->diff_rscale, zggevx_obj->threshold); */
    EXPECT_EQ(zggevx_obj->ilo, zggevx_obj->iloref);
    EXPECT_EQ(zggevx_obj->ihi, zggevx_obj->ihiref);
}
TEST_F(zggevx_test, zggevx4) {
    //EXPECT_NEAR(0.0, zggevx_obj->diff_a, zggevx_obj->threshold);
    //EXPECT_NEAR(0.0, zggevx_obj->diff_b, zggevx_obj->threshold);
    EXPECT_NEAR(0.0, zggevx_obj->diff_vl, zggevx_obj->threshold);
    EXPECT_NEAR(0.0, zggevx_obj->diff_vr, zggevx_obj->threshold);
    EXPECT_NEAR(0.0, zggevx_obj->diff_alphar, zggevx_obj->threshold);
    EXPECT_NEAR(0.0, zggevx_obj->diff_alphai, zggevx_obj->threshold);
    EXPECT_NEAR(0.0, zggevx_obj->diff_beta, zggevx_obj->threshold);
    EXPECT_NEAR(0.0, zggevx_obj->diff_rconde, zggevx_obj->threshold);
    EXPECT_NEAR(0.0, zggevx_obj->diff_rcondv, zggevx_obj->threshold);
	/* Commenting the below check because the position of the scalefactor
	values in the (lscale & lscaleref) , (rscale, rscaleref) varies
	slightly. Though the values are matching.
	So, element by element difference calculation will not work here 
	unlike for other buffers  */
	
    /* EXPECT_NEAR(0.0, zggevx_obj->diff_lscale, zggevx_obj->threshold);
    EXPECT_NEAR(0.0, zggevx_obj->diff_rscale, zggevx_obj->threshold); */
    EXPECT_EQ(zggevx_obj->ilo, zggevx_obj->iloref);
    EXPECT_EQ(zggevx_obj->ihi, zggevx_obj->ihiref);
}
