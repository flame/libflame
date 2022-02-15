#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"

#define trsna_free() \
    if (t!=NULL)        free(t); \
    if (tref!=NULL)     free(tref); \
    if (select!=NULL)   free(select); \
    if (selectref!=NULL) free(selectref); \
    if (z!=NULL)        free(z); \
    if (zref!=NULL)     free(zref); \
    if (vl!=NULL)       free(vl); \
    if (vlref!=NULL)    free(vlref); \
    if (vr!=NULL)       free(vr); \
    if (s!=NULL)        free(s); \
    if (sref!=NULL)     free(sref); \
    if (sep!=NULL)      free(sep); \
    if (sepref!=NULL)   free(sepref); \
    if (wr!=NULL)       free(wr); \
    if (wrref!=NULL)    free(wrref); \
    if (wi!=NULL)       free(wi); \
    if (wiref!=NULL)    free(wiref); \
    if (vrref!=NULL)    free(vrref)

#define trsna_cplx_free() \
    if (t!=NULL)        free(t); \
    if (tref!=NULL)     free(tref); \
    if (select!=NULL)   free(select); \
    if (selectref!=NULL) free(selectref); \
    if (z!=NULL)        free(z); \
    if (zref!=NULL)     free(zref); \
    if (vl!=NULL)       free(vl); \
    if (vlref!=NULL)    free(vlref); \
    if (vr!=NULL)       free(vr); \
    if (s!=NULL)        free(s); \
    if (sref!=NULL)     free(sref); \
    if (sep!=NULL)      free(sep); \
    if (sepref!=NULL)   free(sepref); \
    if (w!=NULL)       	free(w); \
    if (wref!=NULL)     free(wref); \
    if (vrref!=NULL)    free(vrref)


// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin trsna_float_common_parameters  class definition */
class trsna_float_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_s, diff_sep;;
    void *hModule, *dModule;

    /* Auxiliary buffers to generate Schur canonical form using hseqr routine */
    float *z, *zref; // o/p unitary or orthogonal matrix Z of the Schur vectors of H.
    float *wr, *wrref; // holds the eigen values real part
    float *wi, *wiref; // holds the eigen values imaginary part
    
    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char job;  // Must be 'E' or 'V' or 'B'.
    char howmny; // Must be 'A' or 'B' or 'S'.
    char *select, *selectref; // Specifies the eigenvectors to be computed.  

    lapack_int n; // order of matrix
    lapack_int ldt;
    lapack_int ldvl;
    lapack_int ldvr;
    lapack_int mm; // The number of columns in the arrays VL and/or VR

    float *vl, *vlref; // starting vectors for the inverse iteration for the left eigenvectors.
    float *vr, *vrref; // starting vectors for the inverse iteration for the right eigenvectors
    float *t, *tref; // The hessenberg matrix / Schur decomposition o/p
    float *h, *href; // The hessenberg matrix / Schur decomposition o/p
    
    /* Output parameters */
    lapack_int m, mref; /* For real variants: number of elements of s and/or sep
	actually used to store the 	estimated condition numbers. */
    float *s, *sref; //  reciprocal condition numbers of the selected eigenvalues
    float *sep, *sepref; //  estimated reciprocal condition numbers of the selected right eigenvectors 
    
    /*Return Values */
    int info, inforef;

   public:
      trsna_float_parameters (int matrix_layout_i, char job_i, char howmny_i,
                              lapack_int n_i );

      ~trsna_float_parameters ();
};

/* Constructor definition  float_common_parameters */
trsna_float_parameters:: trsna_float_parameters (int matrix_layout_i,
                            char job_i, char howmny_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    job = job_i;
    howmny = howmny_i;
    n  = n_i;
    mm = n;
   
    ldt = n;
    ldvl = n+1;
    ldvr = n+1;
    
    hModule = NULL;
    dModule = NULL;
    m = 0;
    mref = 0;
    diff_s = 0;
    diff_sep = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n trsna float: matrix_layout: %d n: %d  job: %c \
howmny: %c  \n", matrix_layout, n, job, howmny);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &t, &tref, n*n );
    lapacke_gtest_alloc_float_buffer_pair( &vl, &vlref, n*ldvl );
    lapacke_gtest_alloc_float_buffer_pair( &vr, &vrref, n*ldvr );
    lapacke_gtest_alloc_char_buffer_pair( &select, &selectref, n );
    lapacke_gtest_alloc_float_buffer_pair( &z, &zref, n*n );
    lapacke_gtest_alloc_float_buffer_pair( &wr, &wrref, n );
    lapacke_gtest_alloc_float_buffer_pair( &wi, &wiref, n );
    lapacke_gtest_alloc_float_buffer_pair( &s, &sref, n );
    lapacke_gtest_alloc_float_buffer_pair( &sep, &sepref, n );

    if( (t==NULL) || (tref==NULL) ||  \
        (z==NULL) || (zref==NULL) || \
        (vl==NULL) || (vlref==NULL) || \
        (vr==NULL) || (vrref==NULL) || \
        (wr==NULL) || (wrref==NULL) || \
        (wi==NULL) || (wiref==NULL) || \
        (s==NULL) || (sref==NULL) || \
        (sep==NULL) || (sepref==NULL) || \
        (select==NULL) || (selectref==NULL)  ){
       EXPECT_FALSE( true) << "trsna_float_parameters object: malloc error. Exiting ";
       trsna_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    /* Works for both general matrix, upper traingualr i/p matrices */
    //lapacke_gtest_init_float_buffer_pair_rand_custom_matrix( t, tref, n,ldt, 'U');
    lapacke_gtest_init_float_buffer_pair_rand( t, tref, n*ldt);
    lapacke_gtest_init_float_buffer_pair_rand(vl, vlref, n*ldvl);
    lapacke_gtest_init_float_buffer_pair_rand(vr, vrref, n*ldvr);
    lapacke_gtest_init_char_buffer_pair_with_constant(select, selectref, n, 0xff);
    lapacke_gtest_init_float_buffer_pair_rand(z, zref, n*n);
    lapacke_gtest_init_float_buffer_pair_with_constant(wr, wrref, n, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(wi, wiref, n, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(s, sref, n, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(sep, sepref, n, 0.0);
    
   } /* end of Constructor  */
    

/* Destructor definition  'trsna_float_common_parameters' */
trsna_float_parameters :: ~trsna_float_parameters ()
{
   /* De-Allocate memory for the input matrices */
   trsna_free();
} 

//  Test fixture class definition
class strsna_test  : public  ::testing::Test {
public:
   trsna_float_parameters  *strsna_obj;
   void SetUp();  
   void TearDown () { delete strsna_obj; }
};

void strsna_test::SetUp()
{
    /* LAPACKE STRSNA prototype */
    typedef int (*Fptr_NL_LAPACKE_strsna) ( int matrix_layout, char job,
			char howmny, const lapack_logical* select, lapack_int n,
			const float* t, lapack_int ldt, const float* vl, 
			lapack_int ldvl, const float* vr, lapack_int ldvr, float* s,
			float* sep, lapack_int mm, lapack_int* m);
            
    Fptr_NL_LAPACKE_strsna STRSNA;      
    /* LAPACKE STREVC prototype */
    typedef int (*Fptr_NL_LAPACKE_strevc) ( int matrix_layout, char side,
                char howmny, lapack_logical* select, lapack_int n,
                const float* t, lapack_int ldt, float* vl, lapack_int ldvl,
                float* vr, lapack_int ldvr, lapack_int mm, lapack_int* m );
            
    Fptr_NL_LAPACKE_strevc STREVC;      

    /* LAPACKE SHSEQR prototype */
    typedef int (*Fptr_NL_LAPACKE_shseqr) ( int matrix_layout, char job,
                char compz, lapack_int n, lapack_int ilo, lapack_int ihi,
                float* h, lapack_int ldh, float* wr, float* wi, float* z,
                lapack_int ldz );
    Fptr_NL_LAPACKE_shseqr SHSEQR;



    strsna_obj = new  trsna_float_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].job,
                                         eig_non_sym_paramslist[idx].howmny_trsna,
                                         eig_paramslist[idx].n );
    idx = Circular_Increment_Index(idx);

    strsna_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    strsna_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(strsna_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(strsna_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    STRSNA = (Fptr_NL_LAPACKE_strsna)dlsym(strsna_obj->hModule, "LAPACKE_strsna");
    ASSERT_TRUE(STRSNA != NULL) << "failed to ppt the Netlib LAPACKE_strsna symbol";

    STREVC = (Fptr_NL_LAPACKE_strevc)dlsym(strsna_obj->hModule, "LAPACKE_strevc");
    ASSERT_TRUE(STREVC != NULL) << "failed to ppt the Netlib LAPACKE_strevc symbol";

    SHSEQR = (Fptr_NL_LAPACKE_shseqr)dlsym(strsna_obj->hModule, "LAPACKE_shseqr");
    ASSERT_TRUE(SHSEQR != NULL) << "failed to ppt the Netlib LAPACKE_shseqr symbol";
    

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    strsna_obj->inforef = SHSEQR( strsna_obj->matrix_layout,
                                'S',
                                'I' , 
                                strsna_obj->n,
                                1,
                                strsna_obj->n,
                                strsna_obj->tref,
                                strsna_obj->ldt,
                                strsna_obj->wrref,
                                strsna_obj->wiref,
                                strsna_obj->zref,
                                strsna_obj->n );
                                
    strsna_obj->inforef = STREVC(   strsna_obj->matrix_layout,
                                    'B',
                                    'A', 
                  (lapack_logical *)strsna_obj->selectref,
                                    strsna_obj->n,
                                    strsna_obj->tref, 
                                    strsna_obj->ldt,
                                    strsna_obj->vlref,
                                    strsna_obj->ldvl,
                                    strsna_obj->vrref,
                                    strsna_obj->ldvr,
                                    strsna_obj->mm,
                                    &strsna_obj->mref
                                    );

    strsna_obj->inforef = STRSNA(   strsna_obj->matrix_layout,
                                    strsna_obj->job,
                                    strsna_obj->howmny, 
                  (lapack_logical *)strsna_obj->selectref,
                                    strsna_obj->n,
                                    strsna_obj->tref, 
                                    strsna_obj->ldt,
                                    strsna_obj->vlref,
                                    strsna_obj->ldvl,
                                    strsna_obj->vrref,
                                    strsna_obj->ldvr,
                                    strsna_obj->sref,
                                    strsna_obj->sepref,
                                    strsna_obj->mm,
                                    &strsna_obj->mref
                                    );

    /* Compute libflame's Lapacke o/p  */
    strsna_obj->inforef = LAPACKE_shseqr( strsna_obj->matrix_layout,
                                'S',
                                'I' , 
                                strsna_obj->n,
                                1,
                                strsna_obj->n,
                                strsna_obj->t,
                                strsna_obj->ldt,
                                strsna_obj->wrref,
                                strsna_obj->wiref,
                                strsna_obj->zref,
                                strsna_obj->n );

								
	strsna_obj->info = LAPACKE_strevc( strsna_obj->matrix_layout,
                                    'B',
                                    'A', 
                  (lapack_logical *)strsna_obj->select,
                                    strsna_obj->n,
                                    strsna_obj->t, 
                                    strsna_obj->ldt,
                                    strsna_obj->vl,
                                    strsna_obj->ldvl,
                                    strsna_obj->vr,
                                    strsna_obj->ldvr,
                                    strsna_obj->mm,
                                    &strsna_obj->m
                                    );

#if 1


    strsna_obj->info = LAPACKE_strsna(  strsna_obj->matrix_layout,
                                    strsna_obj->job,
                                    strsna_obj->howmny, 
                  (lapack_logical *)strsna_obj->select,
                                    strsna_obj->n,
                                    strsna_obj->t, 
                                    strsna_obj->ldt,
                                    strsna_obj->vl,
                                    strsna_obj->ldvl,
                                    strsna_obj->vr,
                                    strsna_obj->ldvr,
                                    strsna_obj->s,
                                    strsna_obj->sep,
                                    strsna_obj->mm,
                                    &strsna_obj->m
                                    );
#endif

    /* Capture the Netlib, libflame o/p buffers' differences */
    /* If howmny = 'A' or 'S', then vl need not be set. */
    strsna_obj->diff_s =  computeDiff_s( strsna_obj->n, 
                strsna_obj->s, strsna_obj->sref );
    
    strsna_obj->diff_sep =  computeDiff_s( strsna_obj->n, 
                strsna_obj->sep, strsna_obj->sepref );

}

TEST_F(strsna_test, strsna1) {
    EXPECT_NEAR(0.0, strsna_obj->diff_sep, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, strsna_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(strsna_test, strsna2) {
    EXPECT_NEAR(0.0, strsna_obj->diff_sep, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, strsna_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(strsna_test, strsna3) {
    EXPECT_NEAR(0.0, strsna_obj->diff_sep, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, strsna_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(strsna_test, strsna4) {
    EXPECT_NEAR(0.0, strsna_obj->diff_sep, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, strsna_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
}
