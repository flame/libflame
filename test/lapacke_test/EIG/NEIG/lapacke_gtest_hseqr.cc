#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"

#define hseqr_free() \
    if (h!=NULL)        free(h); \
    if (href!=NULL)     free(href); \
    if (z!=NULL)        free(z); \
    if (zref!=NULL)     free(zref); \
    if (wr!=NULL)       free(wr); \
    if (wrref!=NULL)    free(wrref); \
    if (wi!=NULL)       free(wi); \
    if (wiref!=NULL)    free(wiref);

#define hseqr_cplx_free() \
    if (h!=NULL)        free(h); \
    if (href!=NULL)     free(href); \
    if (z!=NULL)        free(z); \
    if (zref!=NULL)     free(zref); \
    if (w!=NULL)       free(w); \
    if (wref!=NULL)    free(wref);

#define LAPACKE_TEST_VERBOSE  (1)
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin hseqr_float_common_parameters  class definition */
class hseqr_float_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_h, diff_z;
    float diff_wr, diff_wi;
    void *hModule, *dModule;
    
    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char job;  // Must be 'N' or 'P' or 'S' or 'B'.
    char compz; // Must be  'N' or 'I' or 'V'.
    lapack_int n; // order of matrix
    lapack_int ilo, iloref;
    lapack_int ihi, ihiref;

    lapack_int ldh; // leading dimension of h
    lapack_int ldz; // leading dimension of z

    /* Input/ Output parameters */
    float *h, *href; // The hessenberg matrix / Schur decomposition o/p
    float *z, *zref; // o/p unitary or orthogonal matrix Z of the Schur vectors of H.
    
    /* Output parameters */
    float *wr, *wrref; // holds the eigen values real part
    float *wi, *wiref; // holds the eigen values imaginary part
    /*Return Values */
    int info, inforef;

   public:
      hseqr_float_parameters (int matrix_layout_i, char job_i, char compz_i,
                                            lapack_int n_i );

      ~hseqr_float_parameters ();
};

/* Constructor definition  float_common_parameters */
hseqr_float_parameters:: hseqr_float_parameters (int matrix_layout_i,
            char job_i,  char compz_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    job = job_i;
    compz = compz_i;
    n  = n_i;
   
    ilo = 1;
    ihi = n;
    ldh = n;
    ldz = n;
    
    hModule = NULL;
    dModule = NULL;
    diff_h = 0;
    diff_z = 0;
    diff_wr = 0;
    diff_wi = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n hseqr float:  n: %d  job: %c compz: %c   \
 matrix_layout: %d \n",  n, job, compz, matrix_layout);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &h, &href, n*ldh );
    lapacke_gtest_alloc_float_buffer_pair( &z, &zref, n*ldz );
    lapacke_gtest_alloc_float_buffer_pair( &wr, &wrref, n );
    lapacke_gtest_alloc_float_buffer_pair( &wi, &wiref, n );

    if( (h==NULL) || (href==NULL) ||  \
        (z==NULL) || (zref==NULL) || \
        (wr==NULL) || (wrref==NULL) || \
        (wi==NULL) || (wiref==NULL)  ){
       EXPECT_FALSE( true) << "hseqr_float_parameters object: malloc error. Exiting ";
       hseqr_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand_custom_matrix( h, href, n,ldh, 'U');
    lapacke_gtest_init_float_buffer_pair_rand(z, zref, n*ldz);
    lapacke_gtest_init_float_buffer_pair_with_constant(wr, wrref, n, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(wi, wiref, n, 0.0);
    
   } /* end of Constructor  */
    

/* Destructor definition  'hseqr_float_common_parameters' */
hseqr_float_parameters :: ~hseqr_float_parameters ()
{
   /* De-Allocate memory for the input matrices */
   hseqr_free();
} 

//  Test fixture class definition
class shseqr_test  : public  ::testing::Test {
public:
   hseqr_float_parameters  *shseqr_obj;
   void SetUp();  
   void TearDown () { delete shseqr_obj; }
};

void shseqr_test::SetUp()
{
    /* LAPACKE SHSEQR prototype */
    typedef int (*Fptr_NL_LAPACKE_shseqr) ( int matrix_layout, char job,
                char compz, lapack_int n, lapack_int ilo, lapack_int ihi,
                float* h, lapack_int ldh, float* wr, float* wi, float* z,
                lapack_int ldz );
            
    /* LAPACKE DGEBAL prototype */
    typedef int (*Fptr_NL_LAPACKE_sgebal) ( int matrix_layout, char job, 
                    lapack_int n, float *A, lapack_int lda, lapack_int* ilo, 
                    lapack_int* ihi, float* scale);
                            
    Fptr_NL_LAPACKE_sgebal SGEBAL;
    Fptr_NL_LAPACKE_shseqr SHSEQR;

    shseqr_obj = new  hseqr_float_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_paramslist[idx].job_seqr,
                                         eig_paramslist[idx].compz,
                                         eig_paramslist[idx].n );
    idx = Circular_Increment_Index(idx);

    shseqr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    shseqr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(shseqr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(shseqr_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SHSEQR = (Fptr_NL_LAPACKE_shseqr)dlsym(shseqr_obj->hModule, "LAPACKE_shseqr");
    ASSERT_TRUE(SHSEQR != NULL) << "failed to ppt the Netlib LAPACKE_shseqr symbol";

    SGEBAL = (Fptr_NL_LAPACKE_sgebal)dlsym(shseqr_obj->hModule, "LAPACKE_sgebal");
    ASSERT_TRUE(SGEBAL != NULL) << "failed to get the Netlib LAPACKE_sgebal symbol";
    
    /* invoke the LAPACKE_sgebal API to balance A 
    shseqr_obj->inforef = SGEBAL( shseqr_obj->matrix_layout, shseqr_obj->job,
                            shseqr_obj->n,shseqr_obj->aref, shseqr_obj->lda, 
                            &shseqr_obj->iloref, &shseqr_obj->ihiref, 
                            shseqr_obj->scaleref);
    
    shseqr_obj->info = LAPACKE_sgebal( shseqr_obj->matrix_layout, shseqr_obj->job,
                            shseqr_obj->n, shseqr_obj->a, shseqr_obj->lda, 
                            &shseqr_obj->ilo, &shseqr_obj->ihi, shseqr_obj->scale);
    
    shseqr_obj->diff_ilo = fabs(shseqr_obj->ilo - shseqr_obj->iloref);
    shseqr_obj->diff_ihi = fabs(shseqr_obj->ihi - shseqr_obj->ihiref);*/

    /* Compute the reference o/p by invoking Netlib-Lapack's API */
    shseqr_obj->inforef = SHSEQR(   shseqr_obj->matrix_layout,
                                    shseqr_obj->job,
                                    shseqr_obj->compz, 
                                    shseqr_obj->n,
                                    shseqr_obj->ilo,
                                    shseqr_obj->ihi,
                                    shseqr_obj->href, 
                                    shseqr_obj->ldh,
                                    shseqr_obj->wrref,
                                    shseqr_obj->wiref,
                                    shseqr_obj->zref,
                                    shseqr_obj->ldz);

    /* Compute libflame's Lapacke o/p  */
    shseqr_obj->info = LAPACKE_shseqr(  shseqr_obj->matrix_layout,
                                        shseqr_obj->job,
                                        shseqr_obj->compz,
                                        shseqr_obj->n,
                                        shseqr_obj->ilo,
                                        shseqr_obj->ihi,
                                        shseqr_obj->h,
                                        shseqr_obj->ldh,
                                        shseqr_obj->wr,
                                        shseqr_obj->wi,
                                        shseqr_obj->z,
                                        shseqr_obj->ldz);
    /* Capture the Netlib, libflame o/p buffers' differences */
    shseqr_obj->diff_h =  computeDiff_s( shseqr_obj->n*shseqr_obj->ldh, 
                shseqr_obj->h, shseqr_obj->href );

    shseqr_obj->diff_z =  computeDiff_s( shseqr_obj->n*shseqr_obj->ldz, 
                shseqr_obj->z, shseqr_obj->zref );

    shseqr_obj->diff_wi =  computeDiff_s( shseqr_obj->n, 
                shseqr_obj->wi, shseqr_obj->wiref );

    shseqr_obj->diff_wr =  computeDiff_s( shseqr_obj->n, 
                shseqr_obj->wr, shseqr_obj->wrref );

}

TEST_F(shseqr_test, shseqr1) {
    EXPECT_NEAR(0.0, shseqr_obj->diff_h, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, shseqr_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, shseqr_obj->diff_wr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, shseqr_obj->diff_wi, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(shseqr_test, shseqr2) {
    EXPECT_NEAR(0.0, shseqr_obj->diff_h, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, shseqr_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, shseqr_obj->diff_wr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, shseqr_obj->diff_wi, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(shseqr_test, shseqr3) {
    EXPECT_NEAR(0.0, shseqr_obj->diff_h, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, shseqr_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, shseqr_obj->diff_wr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, shseqr_obj->diff_wi, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(shseqr_test, shseqr4) {
    EXPECT_NEAR(0.0, shseqr_obj->diff_h, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, shseqr_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, shseqr_obj->diff_wr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, shseqr_obj->diff_wi, LAPACKE_GTEST_THRESHOLD);
}

/* Begin hseqr_double_common_parameters  class definition */
class hseqr_double_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    double diff_h, diff_z;
    double diff_wr, diff_wi;
    void *hModule, *dModule;
    
    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char job;  // Must be 'N' or 'P' or 'S' or 'B'.
    char compz; // Must be  'N' or 'I' or 'V'.
    lapack_int n; // order of matrix
    lapack_int ilo, iloref;
    lapack_int ihi, ihiref;

    lapack_int ldh; // leading dimension of h
    lapack_int ldz; // leading dimension of z

    /* Input/ Output parameters */
    double *h, *href; // The hessenberg matrix / Schur decomposition o/p
    double *z, *zref; // o/p unitary or orthogonal matrix Z of the Schur vectors of H.
    
    /* Output parameters */
    double *wr, *wrref; // holds the eigen values real part
    double *wi, *wiref; // holds the eigen values imaginary part
    /*Return Values */
    int info, inforef;

   public:
      hseqr_double_parameters (int matrix_layout_i, char job_i, char compz_i,
                                            lapack_int n_i );

      ~hseqr_double_parameters ();
};

/* Constructor definition  double_common_parameters */
hseqr_double_parameters:: hseqr_double_parameters (int matrix_layout_i,
            char job_i,  char compz_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    job = job_i;
    compz = compz_i;
    n  = n_i;
   
    ilo = 1;
    ihi = n;
    ldh = n;
    ldz = n;
    
    hModule = NULL;
    dModule = NULL;
    diff_h = 0;
    diff_z = 0;
    diff_wr = 0;
    diff_wi = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n hseqr double:  n: %d  job: %c compz: %c   \
 matrix_layout: %d \n",  n, job, compz, matrix_layout);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &h, &href, n*ldh );
    lapacke_gtest_alloc_double_buffer_pair( &z, &zref, n*ldz );
    lapacke_gtest_alloc_double_buffer_pair( &wr, &wrref, n );
    lapacke_gtest_alloc_double_buffer_pair( &wi, &wiref, n );

    if( (h==NULL) || (href==NULL) ||  \
        (z==NULL) || (zref==NULL) || \
        (wr==NULL) || (wrref==NULL) || \
        (wi==NULL) || (wiref==NULL)  ){
       EXPECT_FALSE( true) << "hseqr_double_parameters object: malloc error. Exiting ";
       hseqr_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand_custom_matrix( h, href, n,ldh, 'U');
    lapacke_gtest_init_double_buffer_pair_rand(z, zref, n*ldz);
    lapacke_gtest_init_double_buffer_pair_with_constant(wr, wrref, n, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(wi, wiref, n, 0.0);
    
   } /* end of Constructor  */
    

/* Destructor definition  'hseqr_double_common_parameters' */
hseqr_double_parameters :: ~hseqr_double_parameters ()
{
   /* De-Allocate memory for the input matrices */
   hseqr_free();
} 

//  Test fixture class definition
class dhseqr_test  : public  ::testing::Test {
public:
   hseqr_double_parameters  *dhseqr_obj;
   void SetUp();  
   void TearDown () { delete dhseqr_obj; }
};

void dhseqr_test::SetUp()
{
    /* LAPACKE DHSEQR prototype */
    typedef int (*Fptr_NL_LAPACKE_dhseqr) ( int matrix_layout, char job,
                char compz, lapack_int n, lapack_int ilo, lapack_int ihi,
                double* h, lapack_int ldh, double* wr, double* wi, double* z,
                lapack_int ldz );
            
    /* LAPACKE DGEBAL prototype */
    typedef int (*Fptr_NL_LAPACKE_sgebal) ( int matrix_layout, char job, 
                    lapack_int n, double *A, lapack_int lda, lapack_int* ilo, 
                    lapack_int* ihi, double* scale);
                            
    Fptr_NL_LAPACKE_sgebal SGEBAL;      
    Fptr_NL_LAPACKE_dhseqr DHSEQR;      

    dhseqr_obj = new  hseqr_double_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_paramslist[idx].job_seqr,
                                         eig_paramslist[idx].compz,
                                         eig_paramslist[idx].n );
    idx = Circular_Increment_Index(idx);

    dhseqr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dhseqr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dhseqr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dhseqr_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DHSEQR = (Fptr_NL_LAPACKE_dhseqr)dlsym(dhseqr_obj->hModule, "LAPACKE_dhseqr");
    ASSERT_TRUE(DHSEQR != NULL) << "failed to ppt the Netlib LAPACKE_dhseqr symbol";

    SGEBAL = (Fptr_NL_LAPACKE_sgebal)dlsym(dhseqr_obj->hModule, "LAPACKE_sgebal");
    ASSERT_TRUE(SGEBAL != NULL) << "failed to get the Netlib LAPACKE_sgebal symbol";
    
    /* invoke the LAPACKE_sgebal API to balance A 
    dhseqr_obj->inforef = SGEBAL( dhseqr_obj->matrix_layout, dhseqr_obj->job,
                            dhseqr_obj->n,dhseqr_obj->aref, dhseqr_obj->lda, 
                            &dhseqr_obj->iloref, &dhseqr_obj->ihiref, 
                            dhseqr_obj->scaleref);
    
    dhseqr_obj->info = LAPACKE_sgebal( dhseqr_obj->matrix_layout, dhseqr_obj->job,
                            dhseqr_obj->n, dhseqr_obj->a, dhseqr_obj->lda, 
                            &dhseqr_obj->ilo, &dhseqr_obj->ihi, dhseqr_obj->scale);
    
    dhseqr_obj->diff_ilo = fabs(dhseqr_obj->ilo - dhseqr_obj->iloref);
    dhseqr_obj->diff_ihi = fabs(dhseqr_obj->ihi - dhseqr_obj->ihiref);*/

    /* Compute the reference o/p by invoking Netlib-Lapack's API */
    dhseqr_obj->inforef = DHSEQR(   dhseqr_obj->matrix_layout,
                                    dhseqr_obj->job,
                                    dhseqr_obj->compz, 
                                    dhseqr_obj->n,
                                    dhseqr_obj->ilo,
                                    dhseqr_obj->ihi,
                                    dhseqr_obj->href, 
                                    dhseqr_obj->ldh,
                                    dhseqr_obj->wrref,
                                    dhseqr_obj->wiref,
                                    dhseqr_obj->zref,
                                    dhseqr_obj->ldz);

    /* Compute libflame's Lapacke o/p  */
    dhseqr_obj->info = LAPACKE_dhseqr(  dhseqr_obj->matrix_layout,
                                        dhseqr_obj->job,
                                        dhseqr_obj->compz,
                                        dhseqr_obj->n,
                                        dhseqr_obj->ilo,
                                        dhseqr_obj->ihi,
                                        dhseqr_obj->h,
                                        dhseqr_obj->ldh,
                                        dhseqr_obj->wr,
                                        dhseqr_obj->wi,
                                        dhseqr_obj->z,
                                        dhseqr_obj->ldz);
    /* Capture the Netlib, libflame o/p buffers' differences */
    dhseqr_obj->diff_h =  computeDiff_d( dhseqr_obj->n*dhseqr_obj->ldh, 
                dhseqr_obj->h, dhseqr_obj->href );

    dhseqr_obj->diff_z =  computeDiff_d( dhseqr_obj->n*dhseqr_obj->ldz, 
                dhseqr_obj->z, dhseqr_obj->zref );

    dhseqr_obj->diff_wi =  computeDiff_d( dhseqr_obj->n, 
                dhseqr_obj->wi, dhseqr_obj->wiref );

    dhseqr_obj->diff_wr =  computeDiff_d( dhseqr_obj->n, 
                dhseqr_obj->wr, dhseqr_obj->wrref );

}

TEST_F(dhseqr_test, dhseqr1) {
    EXPECT_NEAR(0.0, dhseqr_obj->diff_h, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dhseqr_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dhseqr_obj->diff_wr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dhseqr_obj->diff_wi, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dhseqr_test, dhseqr2) {
    EXPECT_NEAR(0.0, dhseqr_obj->diff_h, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dhseqr_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dhseqr_obj->diff_wr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dhseqr_obj->diff_wi, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dhseqr_test, dhseqr3) {
    EXPECT_NEAR(0.0, dhseqr_obj->diff_h, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dhseqr_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dhseqr_obj->diff_wr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dhseqr_obj->diff_wi, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dhseqr_test, dhseqr4) {
    EXPECT_NEAR(0.0, dhseqr_obj->diff_h, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dhseqr_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dhseqr_obj->diff_wr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dhseqr_obj->diff_wi, LAPACKE_GTEST_THRESHOLD);
}

/* Begin hseqr_scomplex_common_parameters  class definition */
class hseqr_scomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_h, diff_z;
    float diff_w;
    void *hModule, *dModule;
    
    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char job;  // Must be 'N' or 'P' or 'S' or 'B'.
    char compz; // Must be  'N' or 'I' or 'V'.
    lapack_int n; // order of matrix
    lapack_int ilo, iloref;
    lapack_int ihi, ihiref;

    lapack_int ldh; // leading dimension of h
    lapack_int ldz; // leading dimension of z

    /* Input/ Output parameters */
    lapack_complex_float *h, *href; // The hessenberg matrix / Schur decomposition o/p
    lapack_complex_float *z, *zref; // o/p unitary or orthogonal matrix Z of the Schur vectors of H.
    
    /* Output parameters */
    lapack_complex_float *w, *wref; // holds the eigen values real part
    /*Return Values */
    int info, inforef;

   public:
      hseqr_scomplex_parameters (int matrix_layout_i, char job_i, char compz_i,
                                            lapack_int n_i );

      ~hseqr_scomplex_parameters ();
};

/* Constructor definition  float_common_parameters */
hseqr_scomplex_parameters:: hseqr_scomplex_parameters (int matrix_layout_i,
            char job_i,  char compz_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    job = job_i;
    compz = compz_i;
    n  = n_i;
   
    ilo = 1;
    ihi = n;
    ldh = n;
    ldz = n;
    
    hModule = NULL;
    dModule = NULL;
    diff_h = 0;
    diff_z = 0;
    diff_w = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n hseqr float:  n: %d  job: %c compz: %c   \
 matrix_layout: %d \n",  n, job, compz, matrix_layout);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &h, &href, n*ldh );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &z, &zref, n*ldz );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &w, &wref, n );

    if( (h==NULL) || (href==NULL) ||  \
        (z==NULL) || (zref==NULL) || \
        (w==NULL) || (wref==NULL)  ){
       EXPECT_FALSE( true) << "hseqr_scomplex_parameters object: malloc error. Exiting ";
       hseqr_cplx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix( h, href, n,ldh, 'U');
    lapacke_gtest_init_scomplex_buffer_pair_rand(z, zref, n*ldz);
    lapacke_gtest_init_scomplex_buffer_pair_with_constant(w, wref, n, 0.0);
    
   } /* end of Constructor  */
    

/* Destructor definition  'hseqr_scomplex_common_parameters' */
hseqr_scomplex_parameters :: ~hseqr_scomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   hseqr_cplx_free();
} 

//  Test fixture class definition
class chseqr_test  : public  ::testing::Test {
public:
   hseqr_scomplex_parameters  *chseqr_obj;
   void SetUp();  
   void TearDown () { delete chseqr_obj; }
};

void chseqr_test::SetUp()
{
    /* LAPACKE CHSEQR prototype */
    typedef int (*Fptr_NL_LAPACKE_chseqr) ( int matrix_layout, char job,
				char compz, lapack_int n, lapack_int ilo, lapack_int ihi,
				lapack_complex_float* h, lapack_int ldh,
				lapack_complex_float* w, lapack_complex_float* z,
				lapack_int ldz );
            
    /* LAPACKE DGEBAL prototype */
    typedef int (*Fptr_NL_LAPACKE_sgebal) ( int matrix_layout, char job, 
                    lapack_int n, float *A, lapack_int lda, lapack_int* ilo, 
                    lapack_int* ihi, float* scale);
                            
    Fptr_NL_LAPACKE_sgebal SGEBAL;      
    Fptr_NL_LAPACKE_chseqr CHSEQR;      

    chseqr_obj = new  hseqr_scomplex_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_paramslist[idx].job_seqr,
                                         eig_paramslist[idx].compz,
                                         eig_paramslist[idx].n );
    idx = Circular_Increment_Index(idx);

    chseqr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    chseqr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(chseqr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(chseqr_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CHSEQR = (Fptr_NL_LAPACKE_chseqr)dlsym(chseqr_obj->hModule, "LAPACKE_chseqr");
    ASSERT_TRUE(CHSEQR != NULL) << "failed to ppt the Netlib LAPACKE_chseqr symbol";

    SGEBAL = (Fptr_NL_LAPACKE_sgebal)dlsym(chseqr_obj->hModule, "LAPACKE_sgebal");
    ASSERT_TRUE(SGEBAL != NULL) << "failed to get the Netlib LAPACKE_sgebal symbol";
    
    /* invoke the LAPACKE_sgebal API to balance A 
    chseqr_obj->inforef = SGEBAL( chseqr_obj->matrix_layout, chseqr_obj->job,
                            chseqr_obj->n,chseqr_obj->aref, chseqr_obj->lda, 
                            &chseqr_obj->iloref, &chseqr_obj->ihiref, 
                            chseqr_obj->scaleref);
    
    chseqr_obj->info = LAPACKE_sgebal( chseqr_obj->matrix_layout, chseqr_obj->job,
                            chseqr_obj->n, chseqr_obj->a, chseqr_obj->lda, 
                            &chseqr_obj->ilo, &chseqr_obj->ihi, chseqr_obj->scale);
    
    chseqr_obj->diff_ilo = fabs(chseqr_obj->ilo - chseqr_obj->iloref);
    chseqr_obj->diff_ihi = fabs(chseqr_obj->ihi - chseqr_obj->ihiref);*/

    /* Compute the reference o/p by invoking Netlib-Lapack's API */
    chseqr_obj->inforef = CHSEQR(   chseqr_obj->matrix_layout,
                                    chseqr_obj->job,
                                    chseqr_obj->compz, 
                                    chseqr_obj->n,
                                    chseqr_obj->ilo,
                                    chseqr_obj->ihi,
                                    chseqr_obj->href, 
                                    chseqr_obj->ldh,
                                    chseqr_obj->wref,
                                    chseqr_obj->zref,
                                    chseqr_obj->ldz);

    /* Compute libflame's Lapacke o/p  */
    chseqr_obj->info = LAPACKE_chseqr(  chseqr_obj->matrix_layout,
                                        chseqr_obj->job,
                                        chseqr_obj->compz,
                                        chseqr_obj->n,
                                        chseqr_obj->ilo,
                                        chseqr_obj->ihi,
                                        chseqr_obj->h,
                                        chseqr_obj->ldh,
                                        chseqr_obj->w,
                                        chseqr_obj->z,
                                        chseqr_obj->ldz);
    /* Capture the Netlib, libflame o/p buffers' differences */
    chseqr_obj->diff_h =  computeDiff_c( chseqr_obj->n*chseqr_obj->ldh, 
                chseqr_obj->h, chseqr_obj->href );

    chseqr_obj->diff_z =  computeDiff_c( chseqr_obj->n*chseqr_obj->ldz, 
                chseqr_obj->z, chseqr_obj->zref );

    chseqr_obj->diff_w =  computeDiff_c( chseqr_obj->n, 
                chseqr_obj->w, chseqr_obj->wref );

}

TEST_F(chseqr_test, chseqr1) {
    EXPECT_NEAR(0.0, chseqr_obj->diff_h, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, chseqr_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, chseqr_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chseqr_test, chseqr2) {
    EXPECT_NEAR(0.0, chseqr_obj->diff_h, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, chseqr_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, chseqr_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chseqr_test, chseqr3) {
    EXPECT_NEAR(0.0, chseqr_obj->diff_h, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, chseqr_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, chseqr_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chseqr_test, chseqr4) {
    EXPECT_NEAR(0.0, chseqr_obj->diff_h, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, chseqr_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, chseqr_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}

/* Begin hseqr_dcomplex_common_parameters  class definition */
class hseqr_dcomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    double diff_h, diff_z;
    double diff_w;
    void *hModule, *dModule;
    
    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char job;  // Must be 'N' or 'P' or 'S' or 'B'.
    char compz; // Must be  'N' or 'I' or 'V'.
    lapack_int n; // order of matrix
    lapack_int ilo, iloref;
    lapack_int ihi, ihiref;

    lapack_int ldh; // leading dimension of h
    lapack_int ldz; // leading dimension of z

    /* Input/ Output parameters */
    lapack_complex_double *h, *href; // The hessenberg matrix / Schur decomposition o/p
    lapack_complex_double *z, *zref; // o/p unitary or orthogonal matrix Z of the Schur vectors of H.
    
    /* Output parameters */
    lapack_complex_double *w, *wref; // holds the eigen values real part
    /*Return Values */
    int info, inforef;

   public:
      hseqr_dcomplex_parameters (int matrix_layout_i, char job_i, char compz_i,
                                            lapack_int n_i );

      ~hseqr_dcomplex_parameters ();
};

/* Constructor definition  double_common_parameters */
hseqr_dcomplex_parameters:: hseqr_dcomplex_parameters (int matrix_layout_i,
            char job_i,  char compz_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    job = job_i;
    compz = compz_i;
    n  = n_i;
   
    ilo = 1;
    ihi = n;
    ldh = n;
    ldz = n;
    
    hModule = NULL;
    dModule = NULL;
    diff_h = 0;
    diff_z = 0;
    diff_w = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n hseqr double:  n: %d  job: %c compz: %c   \
 matrix_layout: %d \n",  n, job, compz, matrix_layout);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &h, &href, n*ldh );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &z, &zref, n*ldz );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &w, &wref, n );

    if( (h==NULL) || (href==NULL) ||  \
        (z==NULL) || (zref==NULL) || \
        (w==NULL) || (wref==NULL)  ){
       EXPECT_FALSE( true) << "hseqr_dcomplex_parameters object: malloc error. Exiting ";
       hseqr_cplx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix( h, href, n,ldh, 'U');
    lapacke_gtest_init_dcomplex_buffer_pair_rand(z, zref, n*ldz);
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant(w, wref, n, 0.0);
    
   } /* end of Constructor  */
    

/* Destructor definition  'hseqr_dcomplex_common_parameters' */
hseqr_dcomplex_parameters :: ~hseqr_dcomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   hseqr_cplx_free();
} 

//  Test fixture class definition
class zhseqr_test  : public  ::testing::Test {
public:
   hseqr_dcomplex_parameters  *zhseqr_obj;
   void SetUp();  
   void TearDown () { delete zhseqr_obj; }
};

void zhseqr_test::SetUp()
{
    /* LAPACKE ZHSEQR prototype */
    typedef int (*Fptr_NL_LAPACKE_zhseqr) ( int matrix_layout, char job,
				char compz, lapack_int n, lapack_int ilo, lapack_int ihi,
				lapack_complex_double* h, lapack_int ldh,
				lapack_complex_double* w, lapack_complex_double* z,
				lapack_int ldz );
            
    /* LAPACKE DGEBAL prototype */
    typedef int (*Fptr_NL_LAPACKE_sgebal) ( int matrix_layout, char job, 
                    lapack_int n, double *A, lapack_int lda, lapack_int* ilo, 
                    lapack_int* ihi, double* scale);
                            
    Fptr_NL_LAPACKE_sgebal SGEBAL;      
    Fptr_NL_LAPACKE_zhseqr ZHSEQR;      

    zhseqr_obj = new  hseqr_dcomplex_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_paramslist[idx].job_seqr,
                                         eig_paramslist[idx].compz,
                                         eig_paramslist[idx].n );
    idx = Circular_Increment_Index(idx);

    zhseqr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zhseqr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zhseqr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zhseqr_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZHSEQR = (Fptr_NL_LAPACKE_zhseqr)dlsym(zhseqr_obj->hModule, "LAPACKE_zhseqr");
    ASSERT_TRUE(ZHSEQR != NULL) << "failed to ppt the Netlib LAPACKE_zhseqr symbol";

    SGEBAL = (Fptr_NL_LAPACKE_sgebal)dlsym(zhseqr_obj->hModule, "LAPACKE_sgebal");
    ASSERT_TRUE(SGEBAL != NULL) << "failed to get the Netlib LAPACKE_sgebal symbol";
    
    /* invoke the LAPACKE_sgebal API to balance A 
    zhseqr_obj->inforef = SGEBAL( zhseqr_obj->matrix_layout, zhseqr_obj->job,
                            zhseqr_obj->n,zhseqr_obj->aref, zhseqr_obj->lda, 
                            &zhseqr_obj->iloref, &zhseqr_obj->ihiref, 
                            zhseqr_obj->scaleref);
    
    zhseqr_obj->info = LAPACKE_sgebal( zhseqr_obj->matrix_layout, zhseqr_obj->job,
                            zhseqr_obj->n, zhseqr_obj->a, zhseqr_obj->lda, 
                            &zhseqr_obj->ilo, &zhseqr_obj->ihi, zhseqr_obj->scale);
    
    zhseqr_obj->diff_ilo = fabs(zhseqr_obj->ilo - zhseqr_obj->iloref);
    zhseqr_obj->diff_ihi = fabs(zhseqr_obj->ihi - zhseqr_obj->ihiref);*/

    /* Compute the reference o/p by invoking Netlib-Lapack's API */
    zhseqr_obj->inforef = ZHSEQR(   zhseqr_obj->matrix_layout,
                                    zhseqr_obj->job,
                                    zhseqr_obj->compz, 
                                    zhseqr_obj->n,
                                    zhseqr_obj->ilo,
                                    zhseqr_obj->ihi,
                                    zhseqr_obj->href, 
                                    zhseqr_obj->ldh,
                                    zhseqr_obj->wref,
                                    zhseqr_obj->zref,
                                    zhseqr_obj->ldz);

    /* Compute libflame's Lapacke o/p  */
    zhseqr_obj->info = LAPACKE_zhseqr(  zhseqr_obj->matrix_layout,
                                        zhseqr_obj->job,
                                        zhseqr_obj->compz,
                                        zhseqr_obj->n,
                                        zhseqr_obj->ilo,
                                        zhseqr_obj->ihi,
                                        zhseqr_obj->h,
                                        zhseqr_obj->ldh,
                                        zhseqr_obj->w,
                                        zhseqr_obj->z,
                                        zhseqr_obj->ldz);
    /* Capture the Netlib, libflame o/p buffers' differences */
    zhseqr_obj->diff_h =  computeDiff_z( zhseqr_obj->n*zhseqr_obj->ldh, 
                zhseqr_obj->h, zhseqr_obj->href );

    zhseqr_obj->diff_z =  computeDiff_z( zhseqr_obj->n*zhseqr_obj->ldz, 
                zhseqr_obj->z, zhseqr_obj->zref );

    zhseqr_obj->diff_w =  computeDiff_z( zhseqr_obj->n, 
                zhseqr_obj->w, zhseqr_obj->wref );

}

TEST_F(zhseqr_test, zhseqr1) {
    EXPECT_NEAR(0.0, zhseqr_obj->diff_h, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zhseqr_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zhseqr_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhseqr_test, zhseqr2) {
    EXPECT_NEAR(0.0, zhseqr_obj->diff_h, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zhseqr_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zhseqr_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhseqr_test, zhseqr3) {
    EXPECT_NEAR(0.0, zhseqr_obj->diff_h, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zhseqr_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zhseqr_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhseqr_test, zhseqr4) {
    EXPECT_NEAR(0.0, zhseqr_obj->diff_h, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zhseqr_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zhseqr_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}


