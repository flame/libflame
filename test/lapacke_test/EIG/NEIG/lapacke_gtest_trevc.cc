#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define trevc_free() \
    if (t!=NULL)        free(t); \
    if (tref!=NULL)     free(tref); \
    if (select!=NULL)   free(select); \
    if (selectref!=NULL) free(selectref); \
    if (z!=NULL)        free(z); \
    if (zref!=NULL)     free(zref); \
    if (vl!=NULL)       free(vl); \
    if (vlref!=NULL)    free(vlref); \
    if (vr!=NULL)       free(vr); \
    if (wr!=NULL)       free(wr); \
    if (wrref!=NULL)    free(wrref); \
    if (wi!=NULL)       free(wi); \
    if (wiref!=NULL)    free(wiref); \
    if (vrref!=NULL)    free(vrref)

#define trevc_cplx_free() \
    if (t!=NULL)        free(t); \
    if (tref!=NULL)     free(tref); \
    if (select!=NULL)   free(select); \
    if (selectref!=NULL) free(selectref); \
    if (z!=NULL)        free(z); \
    if (zref!=NULL)     free(zref); \
    if (vl!=NULL)       free(vl); \
    if (vlref!=NULL)    free(vlref); \
    if (vr!=NULL)       free(vr); \
    if (w!=NULL)       	free(w); \
    if (wref!=NULL)     free(wref); \
    if (vrref!=NULL)    free(vrref)


// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin trevc_float_common_parameters  class definition */
class trevc_float_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_vl, diff_vr;;
    void *hModule, *dModule;

    /* Intermediate buffers to generate Schur canonical form using hseqr routine */
    float *z, *zref; // o/p unitary or orthogonal matrix Z of the Schur vectors of H.
    float *wr, *wrref; // holds the eigen values real part
    float *wi, *wiref; // holds the eigen values imaginary part
    
    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char side;  // Must be 'R' or 'L' or 'B'.
    char howmny; // Must be 'A' or 'B' or 'S'.
    char *select, *selectref; // Specifies the eigenvectors to be computed.  

    lapack_int n; // order of matrix
    lapack_int ldt;
    lapack_int ldvl;
    lapack_int ldvr;
    lapack_int mm; // The number of columns in the arrays VL and/or VR

    /* Input/ Output parameters */
    float *vl, *vlref; // starting vectors for the inverse iteration for the left eigenvectors.
    float *vr, *vrref; // starting vectors for the inverse iteration for the right eigenvectors
    float *t, *tref; // The hessenberg matrix / Schur decomposition o/p
    float *h, *href; // The hessenberg matrix / Schur decomposition o/p
    
    /* Output parameters */
    lapack_int m, mref; // The number of columns in the arrays VL and/or VR required to store the eigenvectors
    
    /*Return Values */
    int info, inforef;

   public:
      trevc_float_parameters (int matrix_layout_i, char side_i, char howmny_i,
                              lapack_int n_i );

      ~trevc_float_parameters ();
};

/* Constructor definition  float_common_parameters */
trevc_float_parameters:: trevc_float_parameters (int matrix_layout_i,
                            char side_i, char howmny_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    side = side_i;
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
    diff_vl = 0;
    diff_vr = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n trevc float: matrix_layout: %d n: %d  side: %c \
howmny: %c  \n", matrix_layout, n, side, howmny);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &t, &tref, n*n );
    lapacke_gtest_alloc_float_buffer_pair( &vl, &vlref, n*ldvl );
    lapacke_gtest_alloc_float_buffer_pair( &vr, &vrref, n*ldvr );
    lapacke_gtest_alloc_char_buffer_pair( &select, &selectref, n );
    lapacke_gtest_alloc_float_buffer_pair( &z, &zref, n*n );
    lapacke_gtest_alloc_float_buffer_pair( &wr, &wrref, n );
    lapacke_gtest_alloc_float_buffer_pair( &wi, &wiref, n );

    if( (t==NULL) || (tref==NULL) ||  \
        (z==NULL) || (zref==NULL) || \
        (vl==NULL) || (vlref==NULL) || \
        (vr==NULL) || (vrref==NULL) || \
        (wr==NULL) || (wrref==NULL) || \
        (wi==NULL) || (wiref==NULL) || \
        (select==NULL) || (selectref==NULL)  ){
       EXPECT_FALSE( true) << "trevc_float_parameters object: malloc error. Exiting ";
       trevc_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand_custom_matrix( t, tref, n, ldt, 'U');
    lapacke_gtest_init_float_buffer_pair_rand(vl, vlref, n*ldvl);
    lapacke_gtest_init_float_buffer_pair_rand(vr, vrref, n*ldvr);
    lapacke_gtest_init_char_buffer_pair_with_constant(select, selectref, n, 0xff);
    lapacke_gtest_init_float_buffer_pair_rand(z, zref, n*n);
    lapacke_gtest_init_float_buffer_pair_with_constant(wr, wrref, n, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(wi, wiref, n, 0.0);
    
   } /* end of Constructor  */
    

/* Destructor definition  'trevc_float_common_parameters' */
trevc_float_parameters :: ~trevc_float_parameters ()
{
   /* De-Allocate memory for the input matrices */
   trevc_free();
} 

//  Test fixture class definition
class strevc_test  : public  ::testing::Test {
public:
   trevc_float_parameters  *strevc_obj;
   void SetUp();  
   void TearDown () { delete strevc_obj; }
};

void strevc_test::SetUp()
{
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

    strevc_obj = new  trevc_float_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_paramslist[idx].side,
                                         eig_non_sym_paramslist[idx].howmny,
                                         eig_paramslist[idx].n );
    idx = Circular_Increment_Index(idx);

    strevc_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    strevc_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(strevc_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(strevc_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    STREVC = (Fptr_NL_LAPACKE_strevc)dlsym(strevc_obj->hModule, "LAPACKE_strevc");
    ASSERT_TRUE(STREVC != NULL) << "failed to ppt the Netlib LAPACKE_strevc symbol";

    SHSEQR = (Fptr_NL_LAPACKE_shseqr)dlsym(strevc_obj->hModule, "LAPACKE_shseqr");
    ASSERT_TRUE(SHSEQR != NULL) << "failed to ppt the Netlib LAPACKE_shseqr symbol";
    

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    /*  make 'strevc_obj->tref' in Schur canonical form with hseqr API
        NOTE: The test also passes even withour this pre conversion to 
        Schur canonical form    */
    strevc_obj->inforef = SHSEQR( strevc_obj->matrix_layout,
                                'S',
                                'I' , 
                                strevc_obj->n,
                                1,
                                strevc_obj->n,
                                strevc_obj->tref,
                                strevc_obj->ldt,
                                strevc_obj->wrref,
                                strevc_obj->wiref,
                                strevc_obj->zref,
                                strevc_obj->n );
                                
    strevc_obj->inforef = STREVC(   strevc_obj->matrix_layout,
                                    strevc_obj->side,
                                    strevc_obj->howmny, 
                  (lapack_logical *)strevc_obj->selectref,
                                    strevc_obj->n,
                                    strevc_obj->tref, 
                                    strevc_obj->ldt,
                                    strevc_obj->vlref,
                                    strevc_obj->ldvl,
                                    strevc_obj->vrref,
                                    strevc_obj->ldvr,
                                    strevc_obj->mm,
                                    &strevc_obj->mref
                                    );

    /* Compute libflame's Lapacke o/p  */
    /*  make 'strevc_obj->t' in Schur canonical form with hseqr API 
        NOTE: The test also passes even withour this pre conversion to 
        Schur canonical form    */
strevc_obj->info = LAPACKE_shseqr( strevc_obj->matrix_layout,
                                'S',
                                'I' , 
                                strevc_obj->n,
                                1,
                                strevc_obj->n,
                                strevc_obj->t,
                                strevc_obj->ldt,
                                strevc_obj->wr,
                                strevc_obj->wi,
                                strevc_obj->z,
                                strevc_obj->n );

    strevc_obj->info = LAPACKE_strevc(  strevc_obj->matrix_layout,
                                    strevc_obj->side,
                                    strevc_obj->howmny, 
                  (lapack_logical *)strevc_obj->select,
                                    strevc_obj->n,
                                    strevc_obj->t, 
                                    strevc_obj->ldt,
                                    strevc_obj->vl,
                                    strevc_obj->ldvl,
                                    strevc_obj->vr,
                                    strevc_obj->ldvr,
                                    strevc_obj->mm,
                                    &strevc_obj->m
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    /* If howmny = 'A' or 'S', then vl need not be set. */
    if( (strevc_obj->howmny == 'B') && ((strevc_obj->side == 'L') || (strevc_obj->side == 'B')) ){
        strevc_obj->diff_vl =  computeDiff_s( strevc_obj->n * strevc_obj->ldvl, 
                strevc_obj->vl, strevc_obj->vlref );
    }
    
    if( (strevc_obj->side == 'R') || (strevc_obj->side == 'B') ){
        strevc_obj->diff_vr =  computeDiff_s( strevc_obj->n * strevc_obj->n, 
                strevc_obj->vr, strevc_obj->vrref );
    }

}

TEST_F(strevc_test, strevc1) {
    EXPECT_NEAR(0.0, strevc_obj->diff_vr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, strevc_obj->diff_vl, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(strevc_test, strevc2) {
    EXPECT_NEAR(0.0, strevc_obj->diff_vr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, strevc_obj->diff_vl, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(strevc_test, strevc3) {
    EXPECT_NEAR(0.0, strevc_obj->diff_vr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, strevc_obj->diff_vl, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(strevc_test, strevc4) {
    EXPECT_NEAR(0.0, strevc_obj->diff_vr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, strevc_obj->diff_vl, LAPACKE_GTEST_THRESHOLD);
}

/* Begin trevc_double_common_parameters  class definition */
class trevc_double_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    double diff_vl, diff_vr;;
    void *hModule, *dModule;

    /* Intermediate buffers to generate Schur canonical form using hseqr routine */
    double *z, *zref; // o/p unitary or orthogonal matrix Z of the Schur vectors of H.
    double *wr, *wrref; // holds the eigen values real part
    double *wi, *wiref; // holds the eigen values imaginary part
    
    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char side;  // Must be 'R' or 'L' or 'B'.
    char howmny; // Must be 'A' or 'B' or 'S'.
    char *select, *selectref; // Specifies the eigenvectors to be computed.  

    lapack_int n; // order of matrix
    lapack_int ldt;
    lapack_int ldvl;
    lapack_int ldvr;
    lapack_int mm; // The number of columns in the arrays VL and/or VR

    /* Input/ Output parameters */
    double *vl, *vlref; // starting vectors for the inverse iteration for the left eigenvectors.
    double *vr, *vrref; // starting vectors for the inverse iteration for the right eigenvectors
    double *t, *tref; // The hessenberg matrix / Schur decomposition o/p
    double *h, *href; // The hessenberg matrix / Schur decomposition o/p
    
    /* Output parameters */
    lapack_int m, mref; // The number of columns in the arrays VL and/or VR required to store the eigenvectors
    
    /*Return Values */
    int info, inforef;

   public:
      trevc_double_parameters (int matrix_layout_i, char side_i, char howmny_i,
                              lapack_int n_i );

      ~trevc_double_parameters ();
};

/* Constructor definition  double_common_parameters */
trevc_double_parameters:: trevc_double_parameters (int matrix_layout_i,
                            char side_i, char howmny_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    side = side_i;
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
    diff_vl = 0;
    diff_vr = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n trevc double: matrix_layout: %d n: %d  side: %c \
howmny: %c  \n", matrix_layout, n, side, howmny);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &t, &tref, n*n );
    lapacke_gtest_alloc_double_buffer_pair( &vl, &vlref, n*ldvl );
    lapacke_gtest_alloc_double_buffer_pair( &vr, &vrref, n*ldvr );
    lapacke_gtest_alloc_char_buffer_pair( &select, &selectref, n );
    lapacke_gtest_alloc_double_buffer_pair( &z, &zref, n*n );
    lapacke_gtest_alloc_double_buffer_pair( &wr, &wrref, n );
    lapacke_gtest_alloc_double_buffer_pair( &wi, &wiref, n );

    if( (t==NULL) || (tref==NULL) ||  \
        (z==NULL) || (zref==NULL) || \
        (vl==NULL) || (vlref==NULL) || \
        (vr==NULL) || (vrref==NULL) || \
        (wr==NULL) || (wrref==NULL) || \
        (wi==NULL) || (wiref==NULL) || \
        (select==NULL) || (selectref==NULL)  ){
       EXPECT_FALSE( true) << "trevc_double_parameters object: malloc error. Exiting ";
       trevc_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand_custom_matrix( t, tref, n,ldt, 'U');
    lapacke_gtest_init_double_buffer_pair_rand(vl, vlref, n*ldvl);
    lapacke_gtest_init_double_buffer_pair_rand(vr, vrref, n*ldvr);
    lapacke_gtest_init_char_buffer_pair_with_constant(select, selectref, n, 0xff);
    lapacke_gtest_init_double_buffer_pair_rand(z, zref, n*n);
    lapacke_gtest_init_double_buffer_pair_with_constant(wr, wrref, n, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(wi, wiref, n, 0.0);
    
   } /* end of Constructor  */
    

/* Destructor definition  'trevc_double_common_parameters' */
trevc_double_parameters :: ~trevc_double_parameters ()
{
   /* De-Allocate memory for the input matrices */
   trevc_free();
} 

//  Test fixture class definition
class dtrevc_test  : public  ::testing::Test {
public:
   trevc_double_parameters  *dtrevc_obj;
   void SetUp();  
   void TearDown () { delete dtrevc_obj; }
};

void dtrevc_test::SetUp()
{
    /* LAPACKE DTREVC prototype */
    typedef int (*Fptr_NL_LAPACKE_dtrevc) ( int matrix_layout, char side,
                char howmny, lapack_logical* select, lapack_int n,
                const double* t, lapack_int ldt, double* vl, lapack_int ldvl,
                double* vr, lapack_int ldvr, lapack_int mm, lapack_int* m );
            
    Fptr_NL_LAPACKE_dtrevc DTREVC;      

    /* LAPACKE DHSEQR prototype */
    typedef int (*Fptr_NL_LAPACKE_dhseqr) ( int matrix_layout, char job,
                char compz, lapack_int n, lapack_int ilo, lapack_int ihi,
                double* h, lapack_int ldh, double* wr, double* wi, double* z,
                lapack_int ldz );
    Fptr_NL_LAPACKE_dhseqr DHSEQR;

    dtrevc_obj = new  trevc_double_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_paramslist[idx].side,
                                         eig_non_sym_paramslist[idx].howmny,
                                         eig_paramslist[idx].n );
    idx = Circular_Increment_Index(idx);

    dtrevc_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dtrevc_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dtrevc_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dtrevc_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DTREVC = (Fptr_NL_LAPACKE_dtrevc)dlsym(dtrevc_obj->hModule, "LAPACKE_dtrevc");
    ASSERT_TRUE(DTREVC != NULL) << "failed to ppt the Netlib LAPACKE_dtrevc symbol";

    DHSEQR = (Fptr_NL_LAPACKE_dhseqr)dlsym(dtrevc_obj->hModule, "LAPACKE_dhseqr");
    ASSERT_TRUE(DHSEQR != NULL) << "failed to ppt the Netlib LAPACKE_dhseqr symbol";
    

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    /*  make 'dtrevc_obj->tref' in Schur canonical form with hseqr API
        NOTE: The test also passes even withour this pre conversion to 
        Schur canonical form    */
    dtrevc_obj->inforef = DHSEQR( dtrevc_obj->matrix_layout,
                                'S',
                                'I' , 
                                dtrevc_obj->n,
                                1,
                                dtrevc_obj->n,
                                dtrevc_obj->tref,
                                dtrevc_obj->ldt,
                                dtrevc_obj->wrref,
                                dtrevc_obj->wiref,
                                dtrevc_obj->zref,
                                dtrevc_obj->n );
                                
    dtrevc_obj->inforef = DTREVC(   dtrevc_obj->matrix_layout,
                                    dtrevc_obj->side,
                                    dtrevc_obj->howmny, 
                  (lapack_logical *)dtrevc_obj->selectref,
                                    dtrevc_obj->n,
                                    dtrevc_obj->tref, 
                                    dtrevc_obj->ldt,
                                    dtrevc_obj->vlref,
                                    dtrevc_obj->ldvl,
                                    dtrevc_obj->vrref,
                                    dtrevc_obj->ldvr,
                                    dtrevc_obj->mm,
                                    &dtrevc_obj->mref
                                    );

    /* Compute libflame's Lapacke o/p  */
    /*  make 'dtrevc_obj->t' in Schur canonical form with hseqr API 
        NOTE: The test also passes even withour this pre conversion to 
        Schur canonical form    */
dtrevc_obj->info = LAPACKE_dhseqr( dtrevc_obj->matrix_layout,
                                'S',
                                'I' , 
                                dtrevc_obj->n,
                                1,
                                dtrevc_obj->n,
                                dtrevc_obj->t,
                                dtrevc_obj->ldt,
                                dtrevc_obj->wr,
                                dtrevc_obj->wi,
                                dtrevc_obj->z,
                                dtrevc_obj->n );

    dtrevc_obj->info = LAPACKE_dtrevc(  dtrevc_obj->matrix_layout,
                                    dtrevc_obj->side,
                                    dtrevc_obj->howmny, 
                  (lapack_logical *)dtrevc_obj->select,
                                    dtrevc_obj->n,
                                    dtrevc_obj->t, 
                                    dtrevc_obj->ldt,
                                    dtrevc_obj->vl,
                                    dtrevc_obj->ldvl,
                                    dtrevc_obj->vr,
                                    dtrevc_obj->ldvr,
                                    dtrevc_obj->mm,
                                    &dtrevc_obj->m
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    /* If howmny = 'A' or 'S', then vl need not be set. */
    if( (dtrevc_obj->howmny == 'B') && ((dtrevc_obj->side == 'L') || (dtrevc_obj->side == 'B')) ){
        dtrevc_obj->diff_vl =  computeDiff_d( dtrevc_obj->n * dtrevc_obj->ldvl, 
                dtrevc_obj->vl, dtrevc_obj->vlref );
    }
    
    if( (dtrevc_obj->side == 'R') || (dtrevc_obj->side == 'B') ){
        dtrevc_obj->diff_vr =  computeDiff_d( dtrevc_obj->n * dtrevc_obj->ldvr, 
                dtrevc_obj->vr, dtrevc_obj->vrref );
    }

}

TEST_F(dtrevc_test, dtrevc1) {
    EXPECT_NEAR(0.0, dtrevc_obj->diff_vr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dtrevc_obj->diff_vl, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dtrevc_test, dtrevc2) {
    EXPECT_NEAR(0.0, dtrevc_obj->diff_vr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dtrevc_obj->diff_vl, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dtrevc_test, dtrevc3) {
    EXPECT_NEAR(0.0, dtrevc_obj->diff_vr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dtrevc_obj->diff_vl, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dtrevc_test, dtrevc4) {
    EXPECT_NEAR(0.0, dtrevc_obj->diff_vr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dtrevc_obj->diff_vl, LAPACKE_GTEST_THRESHOLD);
}

/* Begin trevc_scomplex_common_parameters  class definition */
class trevc_scomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_vl, diff_vr;;
    void *hModule, *dModule;

    /* Intermediate buffers to generate Schur canonical form using hseqr routine */
    lapack_complex_float *z, *zref; // o/p unitary or orthogonal matrix Z of the Schur vectors of H.
    lapack_complex_float *w, *wref; // holds the complex eigen values
    
    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char side;  // Must be 'R' or 'L' or 'B'.
    char howmny; // Must be 'A' or 'B' or 'S'.
    char *select, *selectref; // Specifies the eigenvectors to be computed.  

    lapack_int n; // order of matrix
    lapack_int ldt;
    lapack_int ldvl;
    lapack_int ldvr;
    lapack_int mm; // The number of columns in the arrays VL and/or VR

    /* Input/ Output parameters */
    lapack_complex_float *vl, *vlref; // starting vectors for the inverse iteration for the left eigenvectors.
    lapack_complex_float *vr, *vrref; // starting vectors for the inverse iteration for the right eigenvectors
    lapack_complex_float *t, *tref; // The hessenberg matrix / Schur decomposition o/p
    lapack_complex_float *h, *href; // The hessenberg matrix / Schur decomposition o/p
    
    /* Output parameters */
    lapack_int m, mref; // The number of columns in the arrays VL and/or VR required to store the eigenvectors
    
    /*Return Values */
    int info, inforef;

   public:
      trevc_scomplex_parameters (int matrix_layout_i, char side_i, char howmny_i,
                              lapack_int n_i );

      ~trevc_scomplex_parameters ();
};

/* Constructor definition  lapack_complex_float_common_parameters */
trevc_scomplex_parameters:: trevc_scomplex_parameters (int matrix_layout_i,
                            char side_i, char howmny_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    side = side_i;
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
    diff_vl = 0;
    diff_vr = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n trevc lapack_complex_float: matrix_layout: %d n: %d  side: %c \
howmny: %c  \n", matrix_layout, n, side, howmny);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &t, &tref, n*n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &vl, &vlref, n*ldvl );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &vr, &vrref, n*ldvr );
    lapacke_gtest_alloc_char_buffer_pair( &select, &selectref, n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &z, &zref, n*n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &w, &wref, n );

    if( (t==NULL) || (tref==NULL) ||  \
        (z==NULL) || (zref==NULL) || \
        (vl==NULL) || (vlref==NULL) || \
        (vr==NULL) || (vrref==NULL) || \
        (w==NULL) || (wref==NULL) || \
        (select==NULL) || (selectref==NULL)  ){
       EXPECT_FALSE( true) << "trevc_scomplex_parameters object: malloc error. Exiting ";
       trevc_cplx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix( t, tref, n,ldt, 'U');
    lapacke_gtest_init_scomplex_buffer_pair_rand(vl, vlref, n*ldvl);
    lapacke_gtest_init_scomplex_buffer_pair_rand(vr, vrref, n*ldvr);
    lapacke_gtest_init_char_buffer_pair_with_constant(select, selectref, n, 0xff);
    lapacke_gtest_init_scomplex_buffer_pair_rand(z, zref, n*n);
    lapacke_gtest_init_scomplex_buffer_pair_with_constant(w, wref, n, 0.0);
    
   } /* end of Constructor  */
    

/* Destructor definition  'trevc_scomplex_common_parameters' */
trevc_scomplex_parameters :: ~trevc_scomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   trevc_cplx_free();
} 

//  Test fixture class definition
class ctrevc_test  : public  ::testing::Test {
public:
   trevc_scomplex_parameters  *ctrevc_obj;
   void SetUp();  
   void TearDown () { delete ctrevc_obj; }
};

void ctrevc_test::SetUp()
{
    /* LAPACKE CTREVC prototype */
    typedef int (*Fptr_NL_LAPACKE_ctrevc) ( int matrix_layout, char side,
                char howmny, lapack_logical* select, lapack_int n,
                const lapack_complex_float* t, lapack_int ldt, lapack_complex_float* vl, lapack_int ldvl,
                lapack_complex_float* vr, lapack_int ldvr, lapack_int mm, lapack_int* m );
            
    Fptr_NL_LAPACKE_ctrevc CTREVC;      

    /* LAPACKE CHSEQR prototype */
    typedef int (*Fptr_NL_LAPACKE_chseqr) ( int matrix_layout, char job,
                char compz, lapack_int n, lapack_int ilo, lapack_int ihi,
                lapack_complex_float* h, lapack_int ldh, lapack_complex_float* w, 
				lapack_complex_float* z, lapack_int ldz );

    Fptr_NL_LAPACKE_chseqr CHSEQR;

    ctrevc_obj = new  trevc_scomplex_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_paramslist[idx].side,
                                         eig_non_sym_paramslist[idx].howmny,
                                         eig_paramslist[idx].n );
    idx = Circular_Increment_Index(idx);

    ctrevc_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ctrevc_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ctrevc_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ctrevc_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CTREVC = (Fptr_NL_LAPACKE_ctrevc)dlsym(ctrevc_obj->hModule, "LAPACKE_ctrevc");
    ASSERT_TRUE(CTREVC != NULL) << "failed to ppt the Netlib LAPACKE_ctrevc symbol";

    CHSEQR = (Fptr_NL_LAPACKE_chseqr)dlsym(ctrevc_obj->hModule, "LAPACKE_chseqr");
    ASSERT_TRUE(CHSEQR != NULL) << "failed to ppt the Netlib LAPACKE_chseqr symbol";
    

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    /*  make 'ctrevc_obj->tref' in Schur canonical form with hseqr API
        NOTE: Complex variants do not need this pre conversion to 
        Schur canonical form    */
    /* ctrevc_obj->inforef = CHSEQR( ctrevc_obj->matrix_layout,
                                'S',
                                'I' , 
                                ctrevc_obj->n,
                                1,
                                ctrevc_obj->n,
                                ctrevc_obj->tref,
                                ctrevc_obj->ldt,
                                ctrevc_obj->wref,
                                ctrevc_obj->zref,
                                ctrevc_obj->n );	*/
                                
    ctrevc_obj->inforef = CTREVC(   ctrevc_obj->matrix_layout,
                                    ctrevc_obj->side,
                                    ctrevc_obj->howmny, 
                  (lapack_logical *)ctrevc_obj->selectref,
                                    ctrevc_obj->n,
                                    ctrevc_obj->tref, 
                                    ctrevc_obj->ldt,
                                    ctrevc_obj->vlref,
                                    ctrevc_obj->ldvl,
                                    ctrevc_obj->vrref,
                                    ctrevc_obj->ldvr,
                                    ctrevc_obj->mm,
                                    &ctrevc_obj->mref
                                    );

    /* Compute libflame's Lapacke o/p  */
    /*  make 'ctrevc_obj->t' in Schur canonical form with hseqr API 
        NOTE: Complex variants do not need this pre conversion to 
        Schur canonical form    */
	/* ctrevc_obj->info = LAPACKE_chseqr( ctrevc_obj->matrix_layout,
                                'S',
                                'I' , 
                                ctrevc_obj->n,
                                1,
                                ctrevc_obj->n,
                                ctrevc_obj->t,
                                ctrevc_obj->ldt,
                                ctrevc_obj->w,
                                ctrevc_obj->z,
                                ctrevc_obj->n );	*/

    ctrevc_obj->info = LAPACKE_ctrevc(  ctrevc_obj->matrix_layout,
                                    ctrevc_obj->side,
                                    ctrevc_obj->howmny, 
                  (lapack_logical *)ctrevc_obj->select,
                                    ctrevc_obj->n,
                                    ctrevc_obj->t, 
                                    ctrevc_obj->ldt,
                                    ctrevc_obj->vl,
                                    ctrevc_obj->ldvl,
                                    ctrevc_obj->vr,
                                    ctrevc_obj->ldvr,
                                    ctrevc_obj->mm,
                                    &ctrevc_obj->m
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    /* If howmny = 'A' or 'S', then vl need not be set. */
    if( (ctrevc_obj->howmny == 'B') && ((ctrevc_obj->side == 'L') || (ctrevc_obj->side == 'B')) ){
        ctrevc_obj->diff_vl =  computeDiff_c( ctrevc_obj->n * ctrevc_obj->ldvl, 
                ctrevc_obj->vl, ctrevc_obj->vlref );
    }
    
    if( (ctrevc_obj->side == 'R') || (ctrevc_obj->side == 'B') ){
        ctrevc_obj->diff_vr =  computeDiff_c( ctrevc_obj->n * ctrevc_obj->ldvr, 
                ctrevc_obj->vr, ctrevc_obj->vrref );
    }

}

TEST_F(ctrevc_test, ctrevc1) {
    EXPECT_NEAR(0.0, ctrevc_obj->diff_vr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, ctrevc_obj->diff_vl, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ctrevc_test, ctrevc2) {
    EXPECT_NEAR(0.0, ctrevc_obj->diff_vr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, ctrevc_obj->diff_vl, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ctrevc_test, ctrevc3) {
    EXPECT_NEAR(0.0, ctrevc_obj->diff_vr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, ctrevc_obj->diff_vl, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ctrevc_test, ctrevc4) {
    EXPECT_NEAR(0.0, ctrevc_obj->diff_vr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, ctrevc_obj->diff_vl, LAPACKE_GTEST_THRESHOLD);
}

/* Begin trevc_dcomplex_common_parameters  class definition */
class trevc_dcomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    double diff_vl, diff_vr;;
    void *hModule, *dModule;

    /* Intermediate buffers to generate Schur canonical form using hseqr routine */
    lapack_complex_double *z, *zref; // o/p unitary or orthogonal matrix Z of the Schur vectors of H.
    lapack_complex_double *w, *wref; // holds the complex eigen values
    
    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char side;  // Must be 'R' or 'L' or 'B'.
    char howmny; // Must be 'A' or 'B' or 'S'.
    char *select, *selectref; // Specifies the eigenvectors to be computed.  

    lapack_int n; // order of matrix
    lapack_int ldt;
    lapack_int ldvl;
    lapack_int ldvr;
    lapack_int mm; // The number of columns in the arrays VL and/or VR

    /* Input/ Output parameters */
    lapack_complex_double *vl, *vlref; // starting vectors for the inverse iteration for the left eigenvectors.
    lapack_complex_double *vr, *vrref; // starting vectors for the inverse iteration for the right eigenvectors
    lapack_complex_double *t, *tref; // The hessenberg matrix / Schur decomposition o/p
    lapack_complex_double *h, *href; // The hessenberg matrix / Schur decomposition o/p
    
    /* Output parameters */
    lapack_int m, mref; // The number of columns in the arrays VL and/or VR required to store the eigenvectors
    
    /*Return Values */
    int info, inforef;

   public:
      trevc_dcomplex_parameters (int matrix_layout_i, char side_i, char howmny_i,
                              lapack_int n_i );

      ~trevc_dcomplex_parameters ();
};

/* Constructor definition  lapack_complex_double_common_parameters */
trevc_dcomplex_parameters:: trevc_dcomplex_parameters (int matrix_layout_i,
                            char side_i, char howmny_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    side = side_i;
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
    diff_vl = 0;
    diff_vr = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n trevc lapack_complex_double: matrix_layout: %d n: %d  side: %c \
howmny: %c  \n", matrix_layout, n, side, howmny);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &t, &tref, n*n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &vl, &vlref, n*ldvl );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &vr, &vrref, n*ldvr );
    lapacke_gtest_alloc_char_buffer_pair( &select, &selectref, n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &z, &zref, n*n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &w, &wref, n );

    if( (t==NULL) || (tref==NULL) ||  \
        (z==NULL) || (zref==NULL) || \
        (vl==NULL) || (vlref==NULL) || \
        (vr==NULL) || (vrref==NULL) || \
        (w==NULL) || (wref==NULL) || \
        (select==NULL) || (selectref==NULL)  ){
       EXPECT_FALSE( true) << "trevc_dcomplex_parameters object: malloc error. Exiting ";
       trevc_cplx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix( t, tref, n,ldt, 'U');
    lapacke_gtest_init_dcomplex_buffer_pair_rand(vl, vlref, n*ldvl);
    lapacke_gtest_init_dcomplex_buffer_pair_rand(vr, vrref, n*ldvr);
    lapacke_gtest_init_char_buffer_pair_with_constant(select, selectref, n, 0xff);
    lapacke_gtest_init_dcomplex_buffer_pair_rand(z, zref, n*n);
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant(w, wref, n, 0.0);
    
   } /* end of Constructor  */
    

/* Destructor definition  'trevc_dcomplex_common_parameters' */
trevc_dcomplex_parameters :: ~trevc_dcomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   trevc_cplx_free();
} 

//  Test fixture class definition
class ztrevc_test  : public  ::testing::Test {
public:
   trevc_dcomplex_parameters  *ztrevc_obj;
   void SetUp();  
   void TearDown () { delete ztrevc_obj; }
};

void ztrevc_test::SetUp()
{
    /* LAPACKE ZTREVC prototype */
    typedef int (*Fptr_NL_LAPACKE_ztrevc) ( int matrix_layout, char side,
                char howmny, lapack_logical* select, lapack_int n,
                const lapack_complex_double* t, lapack_int ldt, lapack_complex_double* vl, lapack_int ldvl,
                lapack_complex_double* vr, lapack_int ldvr, lapack_int mm, lapack_int* m );
            
    Fptr_NL_LAPACKE_ztrevc ZTREVC;      

    /* LAPACKE ZHSEQR prototype */
    typedef int (*Fptr_NL_LAPACKE_zhseqr) ( int matrix_layout, char job,
                char compz, lapack_int n, lapack_int ilo, lapack_int ihi,
                lapack_complex_double* h, lapack_int ldh, lapack_complex_double* w, 
				lapack_complex_double* z, lapack_int ldz );

    Fptr_NL_LAPACKE_zhseqr ZHSEQR;

    ztrevc_obj = new  trevc_dcomplex_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_paramslist[idx].side,
                                         eig_non_sym_paramslist[idx].howmny,
                                         eig_paramslist[idx].n );
    idx = Circular_Increment_Index(idx);

    ztrevc_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ztrevc_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ztrevc_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ztrevc_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZTREVC = (Fptr_NL_LAPACKE_ztrevc)dlsym(ztrevc_obj->hModule, "LAPACKE_ztrevc");
    ASSERT_TRUE(ZTREVC != NULL) << "failed to ppt the Netlib LAPACKE_ztrevc symbol";

    ZHSEQR = (Fptr_NL_LAPACKE_zhseqr)dlsym(ztrevc_obj->hModule, "LAPACKE_zhseqr");
    ASSERT_TRUE(ZHSEQR != NULL) << "failed to ppt the Netlib LAPACKE_zhseqr symbol";
    

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    /*  make 'ztrevc_obj->tref' in Schur canonical form with hseqr API
        NOTE: Complex variants do not need this pre conversion to 
        Schur canonical form    */
    /* ztrevc_obj->inforef = ZHSEQR( ztrevc_obj->matrix_layout,
                                'S',
                                'I' , 
                                ztrevc_obj->n,
                                1,
                                ztrevc_obj->n,
                                ztrevc_obj->tref,
                                ztrevc_obj->ldt,
                                ztrevc_obj->wref,
                                ztrevc_obj->zref,
                                ztrevc_obj->n );	*/
                                
    ztrevc_obj->inforef = ZTREVC(   ztrevc_obj->matrix_layout,
                                    ztrevc_obj->side,
                                    ztrevc_obj->howmny, 
                  (lapack_logical *)ztrevc_obj->selectref,
                                    ztrevc_obj->n,
                                    ztrevc_obj->tref, 
                                    ztrevc_obj->ldt,
                                    ztrevc_obj->vlref,
                                    ztrevc_obj->ldvl,
                                    ztrevc_obj->vrref,
                                    ztrevc_obj->ldvr,
                                    ztrevc_obj->mm,
                                    &ztrevc_obj->mref
                                    );

    /* Compute libflame's Lapacke o/p  */
    /*  make 'ztrevc_obj->t' in Schur canonical form with hseqr API 
        NOTE: Complex variants do not need this pre conversion to 
        Schur canonical form    */
	/* ztrevc_obj->info = LAPACKE_zhseqr( ztrevc_obj->matrix_layout,
                                'S',
                                'I' , 
                                ztrevc_obj->n,
                                1,
                                ztrevc_obj->n,
                                ztrevc_obj->t,
                                ztrevc_obj->ldt,
                                ztrevc_obj->w,
                                ztrevc_obj->z,
                                ztrevc_obj->n );	*/

    ztrevc_obj->info = LAPACKE_ztrevc(  ztrevc_obj->matrix_layout,
                                    ztrevc_obj->side,
                                    ztrevc_obj->howmny, 
                  (lapack_logical *)ztrevc_obj->select,
                                    ztrevc_obj->n,
                                    ztrevc_obj->t, 
                                    ztrevc_obj->ldt,
                                    ztrevc_obj->vl,
                                    ztrevc_obj->ldvl,
                                    ztrevc_obj->vr,
                                    ztrevc_obj->ldvr,
                                    ztrevc_obj->mm,
                                    &ztrevc_obj->m
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    /* If howmny = 'A' or 'S', then vl need not be set. */
    if( (ztrevc_obj->howmny == 'B') && ((ztrevc_obj->side == 'L') || (ztrevc_obj->side == 'B')) ){
        ztrevc_obj->diff_vl =  computeDiff_z( ztrevc_obj->n * ztrevc_obj->ldvl, 
                ztrevc_obj->vl, ztrevc_obj->vlref );
    }
    
    if( (ztrevc_obj->side == 'R') || (ztrevc_obj->side == 'B') ){
        ztrevc_obj->diff_vr =  computeDiff_z( ztrevc_obj->n * ztrevc_obj->ldvr, 
                ztrevc_obj->vr, ztrevc_obj->vrref );
    }

}

TEST_F(ztrevc_test, ztrevc1) {
    EXPECT_NEAR(0.0, ztrevc_obj->diff_vr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, ztrevc_obj->diff_vl, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ztrevc_test, ztrevc2) {
    EXPECT_NEAR(0.0, ztrevc_obj->diff_vr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, ztrevc_obj->diff_vl, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ztrevc_test, ztrevc3) {
    EXPECT_NEAR(0.0, ztrevc_obj->diff_vr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, ztrevc_obj->diff_vl, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ztrevc_test, ztrevc4) {
    EXPECT_NEAR(0.0, ztrevc_obj->diff_vr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, ztrevc_obj->diff_vl, LAPACKE_GTEST_THRESHOLD);
}