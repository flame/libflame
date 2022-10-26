#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define trsen_free() \
    if (t!=NULL)        free(t); \
    if (tref!=NULL)     free(tref); \
    if (q!=NULL)        free(q); \
    if (qref!=NULL)     free(qref); \
    if (select!=NULL)   free(select); \
    if (selectref!=NULL) free(selectref); \
    if (s!=NULL)        free(s); \
    if (sref!=NULL)     free(sref); \
    if (sep!=NULL)      free(sep); \
    if (sepref!=NULL)   free(sepref); \
    if (wr!=NULL)       free(wr); \
    if (wrref!=NULL)    free(wrref); \
    if (wi!=NULL)       free(wi); \
    if (wiref!=NULL)    free(wiref);

#define trsen_cplx_free() \
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

/* Begin trsen_float_common_parameters  class definition */
class trsen_float_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_s, diff_sep;;
    void *hModule, *dModule;

    /* Intermediate buffers to generate Schur canonical form using hseqr routine 
    float *z, *zref; // o/p unitary or orthogonal matrix Z of the Schur vectors of H.
    float *wr, *wrref; // holds the eigen values real part
    float *wi, *wiref; // holds the eigen values imaginary part*/
    
    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char job;  // Must be 'N' or 'E' or 'V' or 'B'.
    char compq; // Must be 'V' or 'N' .
    char *select, *selectref; // Specifies the eigenvectors to be computed.  

    lapack_int n; // order of matrix
    lapack_int ldt; //   leading dimension of t
    lapack_int ldq; //   leading dimension of q
	
	
    /* input / Output parameters */
    float *t, *tref; // Theupper quasi-triangular n-by-n matrix T, in Schur canonical form
    float *q, *qref; // The hessenberg matrix / Schur decomposition o/p
	
    /* Output parameters */	
    float *w, *wref; // The recorded eigenvalues of R
    float *wi, *wiref; // reordered eigenvalues of R.
    float *wr, *wrref; // reordered eigenvalues of R
    
    lapack_int m, mref; /* For real variants: number of elements of s and/or sep
	actually used to store the 	estimated condition numbers. */
    float *s, *sref; //  reciprocal condition numbers of the selected eigenvalues
    float *sep, *sepref; //  estimated reciprocal condition numbers of the selected right eigenvectors 
    
    /*Return Values */
    int info, inforef;

   public:
      trsen_float_parameters (int matrix_layout_i, char job_i, char compq_i,
                              lapack_int n_i );

      ~trsen_float_parameters ();
};

/* Constructor definition  float_common_parameters */
trsen_float_parameters:: trsen_float_parameters (int matrix_layout_i,
                            char job_i, char compq_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    job = job_i;
    compq = compq_i;
    n  = n_i;
	info = 0;
	inforef = 0;
   
    ldt = n;
    ldq = n;
    
    hModule = NULL;
    dModule = NULL;
    m = 0;
    mref = 0;
    diff_s = 0;
    diff_sep = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n trsen float: matrix_layout: %d n: %d  job: %c \
compq: %c  \n", matrix_layout, n, job, compq);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &t, &tref, n*ldt );
    lapacke_gtest_alloc_float_buffer_pair( &q, &qref, n*ldq );
    lapacke_gtest_alloc_float_buffer_pair( &w, &wref, n );
    lapacke_gtest_alloc_float_buffer_pair( &wr, &wrref, n );
    lapacke_gtest_alloc_float_buffer_pair( &wi, &wiref, n );
	
    lapacke_gtest_alloc_char_buffer_pair( &select, &selectref, n );
    lapacke_gtest_alloc_float_buffer_pair( &s, &sref, n );
    lapacke_gtest_alloc_float_buffer_pair( &sep, &sepref, n );

    if( (t==NULL) || (tref==NULL) ||  \
        (q==NULL) || (qref==NULL) || \
        (w==NULL) || (wref==NULL) || \
        (wr==NULL) || (wrref==NULL) || \
        (wi==NULL) || (wiref==NULL) || \
        (s==NULL) || (sref==NULL) || \
        (sep==NULL) || (sepref==NULL) || \
        (select==NULL) || (selectref==NULL)  ){
       EXPECT_FALSE( true) << "trsen_float_parameters object: malloc error. Exiting ";
       trsen_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    /* Works for both general matrix, upper traingualr i/p matrices */
    lapacke_gtest_init_float_buffer_pair_rand_custom_matrix( t, tref, n,ldt, 'U');
    //lapacke_gtest_init_float_buffer_pair_rand( t, tref, n*ldt);
    lapacke_gtest_init_float_buffer_pair_rand(w, wref, n);
    lapacke_gtest_init_char_buffer_pair_with_constant(select, selectref, n, 0xff);
    lapacke_gtest_init_float_buffer_pair_rand(q, qref, n*ldq);
    lapacke_gtest_init_float_buffer_pair_with_constant(wr, wrref, n, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(wi, wiref, n, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(s, sref, n, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(sep, sepref, n, 0.0);
    
   } /* end of Constructor  */
    

/* Destructor definition  'trsen_float_common_parameters' */
trsen_float_parameters :: ~trsen_float_parameters ()
{
   /* De-Allocate memory for the input matrices */
   trsen_free();
} 

//  Test fixture class definition
class strsen_test  : public  ::testing::Test {
public:
   trsen_float_parameters  *strsen_obj;
   void SetUp();  
   void TearDown () { delete strsen_obj; }
};

void strsen_test::SetUp()
{
    /* LAPACKE STRSEN prototype */
    typedef int (*Fptr_NL_LAPACKE_strsen) ( int matrix_layout, char job,
				char compq, const lapack_logical* select, lapack_int n, 
				float* t, lapack_int ldt, float* q, lapack_int ldq, 
				float* wr, float* wi, lapack_int* m, float* s, float* sep );			
            
    Fptr_NL_LAPACKE_strsen STRSEN;      
    /* LAPACKE STREVC prototype */
    typedef int (*Fptr_NL_LAPACKE_strevc) ( int matrix_layout, char side,
                char compq, lapack_logical* select, lapack_int n,
                const float* t, lapack_int ldt, float* vl, lapack_int ldvl,
                float* vr, lapack_int ldvr, lapack_int mm, lapack_int* m );
            
    Fptr_NL_LAPACKE_strevc STREVC;      

    /* LAPACKE SHSEQR prototype */
    typedef int (*Fptr_NL_LAPACKE_shseqr) ( int matrix_layout, char job,
                char compz, lapack_int n, lapack_int ilo, lapack_int ihi,
                float* h, lapack_int ldh, float* wr, float* wi, float* z,
                lapack_int ldz );
    Fptr_NL_LAPACKE_shseqr SHSEQR;



    strsen_obj = new  trsen_float_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].job_trsen,
                                         eig_non_sym_paramslist[idx].compq,
                                         eig_paramslist[idx].n );
    idx = Circular_Increment_Index(idx);

    strsen_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    strsen_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(strsen_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(strsen_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    STRSEN = (Fptr_NL_LAPACKE_strsen)dlsym(strsen_obj->hModule, "LAPACKE_strsen");
    ASSERT_TRUE(STRSEN != NULL) << "failed to ppt the Netlib LAPACKE_strsen symbol";

    STREVC = (Fptr_NL_LAPACKE_strevc)dlsym(strsen_obj->hModule, "LAPACKE_strevc");
    ASSERT_TRUE(STREVC != NULL) << "failed to ppt the Netlib LAPACKE_strevc symbol";

    SHSEQR = (Fptr_NL_LAPACKE_shseqr)dlsym(strsen_obj->hModule, "LAPACKE_shseqr");
    ASSERT_TRUE(SHSEQR != NULL) << "failed to ppt the Netlib LAPACKE_shseqr symbol";
    

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
	
#if 0	
    strsen_obj->inforef = SHSEQR( strsen_obj->matrix_layout,
                                'S',
                                'I' , 
                                strsen_obj->n,
                                1,
                                strsen_obj->n,
                                strsen_obj->tref,
                                strsen_obj->ldt,
                                strsen_obj->wrref,
                                strsen_obj->wiref,
                                strsen_obj->zref,
                                strsen_obj->n );
                                
    strsen_obj->inforef = STREVC(   strsen_obj->matrix_layout,
                                    'B',
                                    'A', 
                  (lapack_logical *)strsen_obj->selectref,
                                    strsen_obj->n,
                                    strsen_obj->tref, 
                                    strsen_obj->ldt,
                                    strsen_obj->vlref,
                                    strsen_obj->ldvl,
                                    strsen_obj->vrref,
                                    strsen_obj->ldvr,
                                    strsen_obj->mm,
                                    &strsen_obj->mref
                                    );
#endif

    strsen_obj->info = LAPACKE_strsen(  strsen_obj->matrix_layout,
                                    strsen_obj->job,
                                    strsen_obj->compq, 
                  (lapack_logical *)strsen_obj->select,
                                    strsen_obj->n,
                                    strsen_obj->t, 
                                    strsen_obj->ldt,
                                    strsen_obj->q,
                                    strsen_obj->ldq,
                                    strsen_obj->wr,
                                    strsen_obj->wi,
                                    &strsen_obj->m,
                                    strsen_obj->s,
                                    strsen_obj->sep
                                    );

    strsen_obj->inforef = STRSEN(   strsen_obj->matrix_layout,
                                    strsen_obj->job,
                                    strsen_obj->compq, 
                  (lapack_logical *)strsen_obj->selectref,
                                    strsen_obj->n,
                                    strsen_obj->tref, 
                                    strsen_obj->ldt,
                                    strsen_obj->qref,
                                    strsen_obj->ldq,
                                    strsen_obj->wrref,
                                    strsen_obj->wiref,
                                    &strsen_obj->mref,
                                    strsen_obj->sref,
                                    strsen_obj->sepref
                                    );

    /* Compute libflame's Lapacke o/p  */
#if 0
    strsen_obj->inforef = LAPACKE_shseqr( strsen_obj->matrix_layout,
                                'S',
                                'I' , 
                                strsen_obj->n,
                                1,
                                strsen_obj->n,
                                strsen_obj->t,
                                strsen_obj->ldt,
                                strsen_obj->wrref,
                                strsen_obj->wiref,
                                strsen_obj->zref,
                                strsen_obj->n );

								
	strsen_obj->info = LAPACKE_strevc( strsen_obj->matrix_layout,
                                    'B',
                                    'A', 
                  (lapack_logical *)strsen_obj->select,
                                    strsen_obj->n,
                                    strsen_obj->t, 
                                    strsen_obj->ldt,
                                    strsen_obj->vl,
                                    strsen_obj->ldvl,
                                    strsen_obj->vr,
                                    strsen_obj->ldvr,
                                    strsen_obj->mm,
                                    &strsen_obj->m
                                    );
#endif

    /* Capture the Netlib, libflame o/p buffers' differences */
    /* If compq = 'A' or 'S', then vl need not be set. */
    strsen_obj->diff_s =  computeDiff_s( strsen_obj->n, 
                strsen_obj->s, strsen_obj->sref );
    
    strsen_obj->diff_sep =  computeDiff_s( strsen_obj->n, 
                strsen_obj->sep, strsen_obj->sepref );

}

TEST_F(strsen_test, strsen1) {
    EXPECT_NEAR(0.0, strsen_obj->diff_sep, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, strsen_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(strsen_test, strsen2) {
    EXPECT_NEAR(0.0, strsen_obj->diff_sep, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, strsen_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(strsen_test, strsen3) {
    EXPECT_NEAR(0.0, strsen_obj->diff_sep, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, strsen_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(strsen_test, strsen4) {
    EXPECT_NEAR(0.0, strsen_obj->diff_sep, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, strsen_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
}
