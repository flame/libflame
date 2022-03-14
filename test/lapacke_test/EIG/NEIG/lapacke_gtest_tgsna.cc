#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define tgsna_free() \
    if (a!=NULL)        free(a); \
    if (aref!=NULL)     free(aref); \
    if (a1!=NULL)        free(a1); \
    if (a1ref!=NULL)     free(a1ref); \
    if (b!=NULL)        free(b); \
    if (bref!=NULL)     free(bref); \
    if (b1!=NULL)        free(b1); \
    if (b1ref!=NULL)     free(b1ref); \
    if (select!=NULL)   free(select); \
    if (selectref!=NULL) free(selectref); \
    if (z!=NULL)        free(z); \
    if (zref!=NULL)     free(zref); \
    if (alphar!=NULL)     free(alphar); \
    if (alpharref!=NULL)     free(alpharref); \
    if (alphai!=NULL)     free(alphai); \
    if (alphairef!=NULL)     free(alphairef); \
    if (beta!=NULL)     free(beta); \
    if (betaref!=NULL)     free(betaref); \
    if (vl!=NULL)       free(vl); \
    if (vlref!=NULL)    free(vlref); \
    if (vr!=NULL)       free(vr); \
    if (s!=NULL)        free(s); \
    if (sref!=NULL)     free(sref); \
    if (dif!=NULL)      free(dif); \
    if (difref!=NULL)   free(difref); \
    if (wr!=NULL)       free(wr); \
    if (wrref!=NULL)    free(wrref); \
    if (wi!=NULL)       free(wi); \
    if (wiref!=NULL)    free(wiref); \
    if (vrref!=NULL)    free(vrref)

#define tgsna_cplx_free() \
    if (a!=NULL)        free(a); \
    if (aref!=NULL)     free(aref); \
    if (b!=NULL)        free(b); \
    if (bref!=NULL)     free(bref); \
    if (b1!=NULL)        free(b1); \
    if (b1ref!=NULL)     free(b1ref); \
    if (select!=NULL)   free(select); \
    if (selectref!=NULL) free(selectref); \
    if (z!=NULL)        free(z); \
    if (zref!=NULL)     free(zref); \
    if (alphar!=NULL)     free(alphar); \
    if (alpharref!=NULL)     free(alpharref); \
    if (alphai!=NULL)     free(alphai); \
    if (alphairef!=NULL)     free(alphairef); \
    if (beta!=NULL)     free(beta); \
    if (betaref!=NULL)     free(betaref); \
    if (vl!=NULL)       free(vl); \
    if (vlref!=NULL)    free(vlref); \
    if (vr!=NULL)       free(vr); \
    if (s!=NULL)        free(s); \
    if (sref!=NULL)     free(sref); \
    if (dif!=NULL)      free(dif); \
    if (difref!=NULL)   free(difref); \
    if (w!=NULL)       	free(w); \
    if (wref!=NULL)     free(wref); \
    if (vrref!=NULL)    free(vrref)


// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin tgsna_float_common_parameters  class definition */
class tgsna_float_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_s, diff_dif;;
    void *hModule, *dModule;
    float threshold;

    /* Intermediate buffers to generate Schur canonical form using gges routine */
    float *z, *zref; // o/p unitary or orthogonal matrix Z of the Schur vectors of H.
    float *wr, *wrref; // holds the eigen values real part
    float *wi, *wiref; // holds the eigen values imaginary part
    float *a1, *a1ref; // The hessenberg matrix / Schur decomposition o/p
    float *b1, *b1ref; // The hessenberg matrix / Schur decomposition o/p
    float* alphar, *alpharref;
	float* alphai, *alphairef;
	float* beta, *betaref;
   
    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char job;  // Must be 'E' or 'V' or 'B'.
    char howmny; // Must be 'A' or 'B' or 'S'.
    char *select, *selectref; // Specifies the eigenvectors to be computed.  

    lapack_int n; // order of matrix
    lapack_int lda, ldb;
    lapack_int ldvl;
    lapack_int ldvr;
    lapack_int mm, mmref; // The number of columns in the arrays VL and/or VR

    float *vl, *vlref; // starting vectors for the inverse iteration for the left eigenvectors.
    float *vr, *vrref; // starting vectors for the inverse iteration for the right eigenvectors
    float *a, *aref; // The hessenberg matrix / Schur decomposition o/p
    float *b, *bref; // The hessenberg matrix / Schur decomposition o/p
    
    /* Output parameters */
    lapack_int m, mref; /* For real variants: number of elements of s and/or dif
	actually used to store the 	estimated condition numbers. */
    float *s, *sref; //  reciprocal condition numbers of the selected eigenvalues
    float *dif, *difref; //  estimated reciprocal condition numbers of the selected right eigenvectors 
    
    /*Return Values */
    int info, inforef;

   public:
      tgsna_float_parameters (int matrix_layout_i, char job_i, char howmny_i,
                              lapack_int n_i );

      ~tgsna_float_parameters ();
};

/* Constructor definition  float_common_parameters */
tgsna_float_parameters:: tgsna_float_parameters (int matrix_layout_i,
                            char job_i, char howmny_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    job = job_i;
    howmny = howmny_i;
    n  = n_i;
   
    lda = n;
	ldb = n;
    ldvl = n;
    ldvr = n;
    
    hModule = NULL;
    dModule = NULL;
    m = 0;
    mref = 0;
    diff_s = 0;
    diff_dif = 0;
	mm = n;
	mmref = n;

#if LAPACKE_TEST_VERBOSE
   printf(" \n tgsna float: matrix_layout: %d n: %d  job: %c \
howmny: %c  \n", matrix_layout, n, job, howmny);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, n*lda );
    lapacke_gtest_alloc_float_buffer_pair( &a1, &a1ref, n*lda );
    lapacke_gtest_alloc_float_buffer_pair( &b, &bref, n*ldb );
    lapacke_gtest_alloc_float_buffer_pair( &b1, &b1ref, n*ldb );
    lapacke_gtest_alloc_float_buffer_pair( &vl, &vlref, n*ldvl );
    lapacke_gtest_alloc_float_buffer_pair( &vr, &vrref, n*ldvr );
    lapacke_gtest_alloc_char_buffer_pair( &select, &selectref, n );
    lapacke_gtest_alloc_float_buffer_pair( &z, &zref, n*n );
    lapacke_gtest_alloc_float_buffer_pair( &wr, &wrref, n );
    lapacke_gtest_alloc_float_buffer_pair( &wi, &wiref, n );
    lapacke_gtest_alloc_float_buffer_pair( &s, &sref, n );
    lapacke_gtest_alloc_float_buffer_pair( &dif, &difref, n );
    lapacke_gtest_alloc_float_buffer_pair( &alphar, &alpharref, n );
    lapacke_gtest_alloc_float_buffer_pair( &alphai, &alphairef, n );
    lapacke_gtest_alloc_float_buffer_pair( &beta, &betaref, n );

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (b1==NULL) || (b1ref==NULL) || \
        (z==NULL) || (zref==NULL) || \
        (vl==NULL) || (vlref==NULL) || \
        (vr==NULL) || (vrref==NULL) || \
        (wr==NULL) || (wrref==NULL) || \
        (wi==NULL) || (wiref==NULL) || \
        (s==NULL) || (sref==NULL) || \
        (dif==NULL) || (difref==NULL) || \
        (alphar==NULL) || (alpharref==NULL) || \
        (alphai==NULL) || (alphairef==NULL) || \
        (beta==NULL) || (betaref==NULL) || \
        (select==NULL) || (selectref==NULL)  ){
       EXPECT_FALSE( true) << "tgsna_float_parameters object: malloc error. Exiting ";
       tgsna_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    /* Works for both general matrix, upper traingualr i/p matrices */
    lapacke_gtest_init_float_buffer_pair_rand_custom_matrix( a, aref, n, lda, 'U');
    lapacke_gtest_init_float_buffer_pair_rand_custom_matrix( b, bref, n, ldb, 'U');
	
	// Additional a1, a1ref buffers to generate eigen values through ggev API
	memcpy(a1, a, n*lda*sizeof(float));
	memcpy(a1ref, a, n*lda*sizeof(float));
	memcpy(b1, b, n*ldb*sizeof(float));
	memcpy(b1ref, b, n*ldb*sizeof(float));

    lapacke_gtest_init_float_buffer_pair_rand(vl, vlref, n*ldvl);
    lapacke_gtest_init_float_buffer_pair_rand(vr, vrref, n*ldvr);
    lapacke_gtest_init_char_buffer_pair_with_constant(select, selectref, n, 0xff);
    lapacke_gtest_init_float_buffer_pair_rand(z, zref, n*n);
    lapacke_gtest_init_float_buffer_pair_with_constant(wr, wrref, n, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(wi, wiref, n, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(s, sref, n, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(dif, difref, n, 0.0);
     lapacke_gtest_init_float_buffer_pair_with_constant( alphar, alpharref, n, 0);
    lapacke_gtest_init_float_buffer_pair_with_constant( alphai, alphairef, n, 0);
    lapacke_gtest_init_float_buffer_pair_with_constant( beta, betaref, n, 0);
   
   } /* end of Constructor  */
    

/* Destructor definition  'tgsna_float_common_parameters' */
tgsna_float_parameters :: ~tgsna_float_parameters ()
{
   /* De-Allocate memory for the input matrices */
   tgsna_free();
} 

//  Test fixture class definition
class stgsna_test  : public  ::testing::Test {
public:
   tgsna_float_parameters  *stgsna_obj;
   void SetUp();  
   void TearDown () { delete stgsna_obj; }
};

void stgsna_test::SetUp()
{
    /* LAPACKE Stgsna prototype */
    typedef int (*Fptr_NL_LAPACKE_stgsna) ( int matrix_layout, 
		char job, char howmny, const lapack_logical* select, 
		lapack_int n, const float* a, lapack_int lda, const float* b, 
		lapack_int ldb, const float* vl, lapack_int ldvl, const float* vr,
		lapack_int ldvr, float* s, float* dif, lapack_int mm, lapack_int* m);
            
    Fptr_NL_LAPACKE_stgsna Stgsna;

    /* local variables needed for calling gges, ggev APIs  */	
    lapack_int sdim;
	float vsl, vsr;

    stgsna_obj = new  tgsna_float_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].job,
                                         eig_non_sym_paramslist[idx].howmny_trsna,
                                         eig_paramslist[idx].n );
    stgsna_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;

    idx = Circular_Increment_Index(idx);

    stgsna_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    stgsna_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(stgsna_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(stgsna_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    Stgsna = (Fptr_NL_LAPACKE_stgsna)dlsym(stgsna_obj->hModule, "LAPACKE_stgsna");
    ASSERT_TRUE(Stgsna != NULL) << "failed to ppt the Netlib LAPACKE_stgsna symbol";

    /*  generate Schur form (S,T) of (A,B) through gges API */	
    stgsna_obj->inforef =  LAPACKE_sgges (   stgsna_obj->matrix_layout,
                                    'N',
                                    'N',
                                    'N',
									(LAPACK_S_SELECT3)stgsna_obj->select,
                                    stgsna_obj->n,
                                    stgsna_obj->a,
                                    stgsna_obj->lda,
                                    stgsna_obj->b,
                                    stgsna_obj->ldb,
                                    &sdim,
                                    stgsna_obj->alphar,
                                    stgsna_obj->alphai,
                                    stgsna_obj->beta,
									&vsl,
									1,
									&vsr,
									1
                                    );

	memcpy(stgsna_obj->aref, stgsna_obj->a, stgsna_obj->n*stgsna_obj->lda*sizeof(float));
	memcpy(stgsna_obj->bref, stgsna_obj->b, stgsna_obj->n*stgsna_obj->ldb*sizeof(float));
	
    stgsna_obj->inforef =  LAPACKE_sggev (   stgsna_obj->matrix_layout,
                                    'V',
                                    'V',
                                    stgsna_obj->n,
                                    stgsna_obj->a1,
                                    stgsna_obj->lda,
                                    stgsna_obj->b1,
                                    stgsna_obj->ldb,
                                    stgsna_obj->alphar,
                                    stgsna_obj->alphai,
                                    stgsna_obj->beta,
									stgsna_obj->vl,
									stgsna_obj->ldvl,
									stgsna_obj->vr,
									stgsna_obj->ldvr
                                    );
									
	memcpy(stgsna_obj->vlref, stgsna_obj->vl, (stgsna_obj->n)*(stgsna_obj->ldvl)*(sizeof(float)) );
	memcpy(stgsna_obj->vrref, stgsna_obj->vr, (stgsna_obj->n)*(stgsna_obj->ldvr)*(sizeof(float)) );

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    stgsna_obj->inforef = Stgsna(   stgsna_obj->matrix_layout,
                                    stgsna_obj->job,
                                    stgsna_obj->howmny, 
                  (lapack_logical *)stgsna_obj->selectref,
                                    stgsna_obj->n,
                                    stgsna_obj->aref, 
                                    stgsna_obj->lda,
                                    stgsna_obj->bref, 
                                    stgsna_obj->ldb,
                                    stgsna_obj->vlref,
                                    stgsna_obj->ldvl,
                                    stgsna_obj->vrref,
                                    stgsna_obj->ldvr,
                                    stgsna_obj->sref,
                                    stgsna_obj->difref,
                                    stgsna_obj->mmref,
                                    &stgsna_obj->mref
                                    );

    /* Compute libflame's Lapacke o/p  */
    stgsna_obj->info = LAPACKE_stgsna(   stgsna_obj->matrix_layout,
                                    stgsna_obj->job,
                                    stgsna_obj->howmny, 
                  (lapack_logical *)stgsna_obj->select,
                                    stgsna_obj->n,
                                    stgsna_obj->a, 
                                    stgsna_obj->lda,
                                    stgsna_obj->b, 
                                    stgsna_obj->ldb,
                                    stgsna_obj->vl,
                                    stgsna_obj->ldvl,
                                    stgsna_obj->vr,
                                    stgsna_obj->ldvr,
                                    stgsna_obj->s,
                                    stgsna_obj->dif,
                                    stgsna_obj->mm,
                                    &stgsna_obj->m
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    /* If howmny = 'A' or 'S', then vl need not be set. */
	if ( stgsna_obj->job != 'V') {
		stgsna_obj->diff_s =  computeDiff_s( stgsna_obj->n, 
                stgsna_obj->s, stgsna_obj->sref );
	}
    
	if ( stgsna_obj->job != 'E') {
        stgsna_obj->diff_dif =  computeDiff_s( stgsna_obj->mm, 
                stgsna_obj->dif, stgsna_obj->difref );
	}

}

TEST_F(stgsna_test, stgsna1) {
    EXPECT_NEAR(0.0, stgsna_obj->diff_dif, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, stgsna_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
	EXPECT_EQ(stgsna_obj->m, stgsna_obj->mref);
}

TEST_F(stgsna_test, stgsna2) {
    EXPECT_NEAR(0.0, stgsna_obj->diff_dif, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, stgsna_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
	EXPECT_EQ(stgsna_obj->m, stgsna_obj->mref);
}

TEST_F(stgsna_test, stgsna3) {
    EXPECT_NEAR(0.0, stgsna_obj->diff_dif, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, stgsna_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
	EXPECT_EQ(stgsna_obj->m, stgsna_obj->mref);
}

TEST_F(stgsna_test, stgsna4) {
    EXPECT_NEAR(0.0, stgsna_obj->diff_dif, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, stgsna_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
	EXPECT_EQ(stgsna_obj->m, stgsna_obj->mref);
}
