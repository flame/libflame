#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define gebak_free() \
    if(v != NULL)  		free(v);    \
    if(vref != NULL)  	free(vref); \
    if(a != NULL)  		free(a);    \
    if(aref != NULL)  	free(aref); \
    if(scale != NULL)  	free(scale); \
    if(scaleref!=NULL) 	free(scaleref)
	
#define LAPACKE_TEST_VERBOSE  (1)
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin gebak_float_common_parameters  class definition */
class gebak_float_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_v;
    float diff_ilo, diff_ihi, diff_scale;
    void *hModule, *dModule;
    
    float *a, *aref;
    lapack_int lda;
    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char job;  // Must be 'N' or 'P' or 'S' or 'B'.
    char side; // Must be 'L' or 'R'.
    lapack_int n; // number of rows
    lapack_int ilo, iloref;
    lapack_int ihi, ihiref;
    float *scale, *scaleref; //  permutations / the scaling factor
    lapack_int m; // number of columns
    lapack_int ldv; // leading dimension of v

    /* Input/ Output parameters */
    float *v, *vref; // contains the matrix of left or right eigenvectors

    /*Return Values */
    int info, inforef;

   public:
      gebak_float_parameters (int matrix_layout_i, char job_i, char side_i,
            lapack_int n_i, lapack_int m_i, lapack_int lda_i );

      ~gebak_float_parameters ();
};

/* Constructor definition  float_common_parameters */
gebak_float_parameters:: gebak_float_parameters (int matrix_layout_i,
            char job_i,  char side_i, lapack_int n_i, lapack_int m_i, 
            lapack_int lda_i)
{
    matrix_layout = matrix_layout_i;
    job = job_i;
    side = side_i;
    n  = n_i;
    m = m_i;
    lda = lda_i;
    
    ilo = ihi = 0;
    hModule = NULL;
    dModule = NULL;

    if(matrix_layout==LAPACK_COL_MAJOR){
        ldv = n;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        ldv = m;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }
    
#if LAPACKE_TEST_VERBOSE
   printf(" \n gebak float:  n: %d m: %d job: %c side: %c   \
lda: %d  matrix_layout: %d \n",  n, m, job, side,  
                                          lda, matrix_layout);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, n*lda );
    lapacke_gtest_alloc_float_buffer_pair (&scale, &scaleref, n);   
    lapacke_gtest_alloc_float_buffer_pair( &v, &vref, n*m );

    if( (a==NULL) || (aref==NULL) ||  \
        (v==NULL) || (vref==NULL) || \
        (scale==NULL) || (scaleref==NULL) ){
       EXPECT_FALSE( true) << "gebak_float_parameters object: malloc error. Exiting ";
       gebak_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_float_buffer_pair_with_constant(scale, scaleref, n, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(v, vref, n*m, 0.0);
    
   } /* end of Constructor  */
    

/* Destructor definition  'gebak_float_common_parameters' */
gebak_float_parameters :: ~gebak_float_parameters ()
{
   /* De-Allocate memory for the input matrices */
   gebak_free();
} 

//  Test fixture class definition
class sgebak_test  : public  ::testing::Test {
public:
   gebak_float_parameters  *sgebak_obj;
   void SetUp();  
   void TearDown () { delete sgebak_obj; }
};

void sgebak_test::SetUp()
{
    /* LAPACKE SGEBAK prototype */
    typedef int (*Fptr_NL_LAPACKE_sgebak) ( int matrix_layout, char job,
            char side, lapack_int n, lapack_int ilo, lapack_int ihi,
            const float* scale, lapack_int m, float* v, lapack_int ldv  );
            
    /* LAPACKE DGEBAL prototype */
    typedef int (*Fptr_NL_LAPACKE_sgebal) ( int matrix_layout, char job, 
                    lapack_int n, float *A, lapack_int lda, lapack_int* ilo, 
                    lapack_int* ihi, float* scale);
                            
    Fptr_NL_LAPACKE_sgebal SGEBAL;      
    Fptr_NL_LAPACKE_sgebak SGEBAK;      

    sgebak_obj = new  gebak_float_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_paramslist[idx].job,
                                         eig_paramslist[idx].side,
                                         eig_paramslist[idx].n,
                                         eig_paramslist[idx].m,
                                         eig_paramslist[idx].lda );
    idx = Circular_Increment_Index(idx);

    sgebak_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgebak_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgebak_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgebak_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SGEBAK = (Fptr_NL_LAPACKE_sgebak)dlsym(sgebak_obj->hModule, "LAPACKE_sgebak");
    ASSERT_TRUE(SGEBAK != NULL) << "failed to ppt the Netlib LAPACKE_sgebak symbol";

    SGEBAL = (Fptr_NL_LAPACKE_sgebal)dlsym(sgebak_obj->hModule, "LAPACKE_sgebal");
    ASSERT_TRUE(SGEBAL != NULL) << "failed to get the Netlib LAPACKE_sgebal symbol";
	
    /* invoke the LAPACKE_sgebal API to balance A */
    sgebak_obj->inforef = SGEBAL( sgebak_obj->matrix_layout, sgebak_obj->job,
                            sgebak_obj->n,sgebak_obj->aref, sgebak_obj->lda, 
                            &sgebak_obj->iloref, &sgebak_obj->ihiref, 
                            sgebak_obj->scaleref);
	
    sgebak_obj->info = LAPACKE_sgebal( sgebak_obj->matrix_layout, sgebak_obj->job,
                            sgebak_obj->n, sgebak_obj->a, sgebak_obj->lda, 
                            &sgebak_obj->ilo, &sgebak_obj->ihi, sgebak_obj->scale);
	
    sgebak_obj->diff_ilo = fabs(sgebak_obj->ilo - sgebak_obj->iloref);
    sgebak_obj->diff_ihi = fabs(sgebak_obj->ihi - sgebak_obj->ihiref);
	
    /* Compute the reference o/p by invoking Netlib-Lapack's API */
    sgebak_obj->inforef = SGEBAK( sgebak_obj->matrix_layout, sgebak_obj->job,
                            sgebak_obj->side, sgebak_obj->n, sgebak_obj->iloref,
                            sgebak_obj->ihiref, sgebak_obj->scaleref, 
                            sgebak_obj->m, sgebak_obj->vref, sgebak_obj->ldv );

    /* Compute libflame's Lapacke o/p  */
    sgebak_obj->info = LAPACKE_sgebak( sgebak_obj->matrix_layout, sgebak_obj->job,
                            sgebak_obj->side, sgebak_obj->n, sgebak_obj->ilo,
                            sgebak_obj->ihi, sgebak_obj->scale, 
                            sgebak_obj->m, sgebak_obj->v, sgebak_obj->ldv);

    sgebak_obj->diff_scale =  computeDiff_s( sgebak_obj->n, 
                sgebak_obj->scale, sgebak_obj->scaleref );

    sgebak_obj->diff_v =  computeDiff_s( sgebak_obj->n*sgebak_obj->m, 
                sgebak_obj->v, sgebak_obj->vref );
}

TEST_F(sgebak_test, sgebak1) {
    EXPECT_NEAR(0.0, sgebak_obj->diff_v, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgebak_obj->diff_scale, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgebak_test, sgebak2) {
    EXPECT_NEAR(0.0, sgebak_obj->diff_v, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgebak_obj->diff_scale, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgebak_test, sgebak3) {
    EXPECT_NEAR(0.0, sgebak_obj->diff_v, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgebak_obj->diff_scale, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgebak_test, sgebak4) {
    EXPECT_NEAR(0.0, sgebak_obj->diff_v, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgebak_obj->diff_scale, LAPACKE_GTEST_THRESHOLD);
}

/* Begin gebak_double_parameters  class definition */
class gebak_double_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    double diff_v;
    double diff_ilo, diff_ihi, diff_scale;
    void *hModule, *dModule;
    
    double *a, *aref;
    lapack_int lda;
    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char job;  // Must be 'N' or 'P' or 'S' or 'B'.
    char side; // Must be 'L' or 'R'.
    lapack_int n; // number of rows
    lapack_int ilo, iloref;
    lapack_int ihi, ihiref;
    double *scale, *scaleref; //  permutations / the scaling factor
    lapack_int m; // number of columns
    lapack_int ldv; // leading dimension of v

    /* Input/ Output parameters */
    double *v, *vref; // contains the matrix of left or right eigenvectors

    /*Return Values */
    int info, inforef;

   public:
      gebak_double_parameters (int matrix_layout_i, char job_i, char side_i,
            lapack_int n_i, lapack_int m_i, lapack_int lda_i );

      ~gebak_double_parameters ();
};

/* Constructor definition  double_common_parameters */
gebak_double_parameters:: gebak_double_parameters (int matrix_layout_i,
            char job_i,  char side_i, lapack_int n_i, lapack_int m_i, 
            lapack_int lda_i)
{
    matrix_layout = matrix_layout_i;
    job = job_i;
    side = side_i;
    n  = n_i;
    m = m_i;
    lda = lda_i;
    
    ilo = ihi = 0;
    hModule = NULL;
    dModule = NULL;

    if(matrix_layout==LAPACK_COL_MAJOR){
        ldv = n;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        ldv = m;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }
    
#if LAPACKE_TEST_VERBOSE
   printf(" \n gebak double:  n: %d m: %d job: %c side: %c   \
lda: %d  matrix_layout: %d \n",  n, m, job, side,  
                                          lda, matrix_layout);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, n*lda );
    lapacke_gtest_alloc_double_buffer_pair (&scale, &scaleref, n);   
    lapacke_gtest_alloc_double_buffer_pair( &v, &vref, n*m );

    if( (a==NULL) || (aref==NULL) ||  \
        (v==NULL) || (vref==NULL) || \
        (scale==NULL) || (scaleref==NULL) ){
       EXPECT_FALSE( true) << "gebak_double_parameters object: malloc error. Exiting ";
       gebak_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_double_buffer_pair_with_constant(scale, scaleref, n, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(v, vref, n*m, 0.0);
    
   } /* end of Constructor  */
    

/* Destructor definition  'gebak_double_common_parameters' */
gebak_double_parameters :: ~gebak_double_parameters ()
{
   /* De-Allocate memory for the input matrices */
   gebak_free();
} 

//  Test fixture class definition
class dgebak_test  : public  ::testing::Test {
public:
   gebak_double_parameters  *dgebak_obj;
   void SetUp();  
   void TearDown () { delete dgebak_obj; }
};

void dgebak_test::SetUp()
{
    /* LAPACKE DGEBAK prototype */
    typedef int (*Fptr_NL_LAPACKE_dgebak) ( int matrix_layout, char job,
            char side, lapack_int n, lapack_int ilo, lapack_int ihi,
            const double* scale, lapack_int m, double* v, lapack_int ldv  );
            
    /* LAPACKE DGEBAL prototype */
    typedef int (*Fptr_NL_LAPACKE_dgebal) ( int matrix_layout, char job, 
                    lapack_int n, double *A, lapack_int lda, lapack_int* ilo, 
                    lapack_int* ihi, double* scale);
                            
    Fptr_NL_LAPACKE_dgebal DGEBAL;      
    Fptr_NL_LAPACKE_dgebak DGEBAK;      

    dgebak_obj = new  gebak_double_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_paramslist[idx].job,
                                         eig_paramslist[idx].side,
                                         eig_paramslist[idx].n,
                                         eig_paramslist[idx].m,
                                         eig_paramslist[idx].lda );
    idx = Circular_Increment_Index(idx);

    dgebak_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgebak_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgebak_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgebak_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DGEBAK = (Fptr_NL_LAPACKE_dgebak)dlsym(dgebak_obj->hModule, "LAPACKE_dgebak");
    ASSERT_TRUE(DGEBAK != NULL) << "failed to ppt the Netlib LAPACKE_dgebak symbol";

    DGEBAL = (Fptr_NL_LAPACKE_dgebal)dlsym(dgebak_obj->hModule, "LAPACKE_dgebal");
    ASSERT_TRUE(DGEBAL != NULL) << "failed to get the Netlib LAPACKE_dgebal symbol";
	
    /* invoke the LAPACKE_dgebal API to balance A */
    dgebak_obj->inforef = DGEBAL( dgebak_obj->matrix_layout, dgebak_obj->job,
                            dgebak_obj->n,dgebak_obj->aref, dgebak_obj->lda, 
                            &dgebak_obj->iloref, &dgebak_obj->ihiref, 
                            dgebak_obj->scaleref);
	
    dgebak_obj->info = LAPACKE_dgebal( dgebak_obj->matrix_layout, dgebak_obj->job,
                            dgebak_obj->n, dgebak_obj->a, dgebak_obj->lda, 
                            &dgebak_obj->ilo, &dgebak_obj->ihi, dgebak_obj->scale);
	
    dgebak_obj->diff_ilo = fabs(dgebak_obj->ilo - dgebak_obj->iloref);
    dgebak_obj->diff_ihi = fabs(dgebak_obj->ihi - dgebak_obj->ihiref);
	
    /* Compute the reference o/p by invoking Netlib-Lapack's API */
    dgebak_obj->inforef = DGEBAK( dgebak_obj->matrix_layout, dgebak_obj->job,
                            dgebak_obj->side, dgebak_obj->n, dgebak_obj->iloref,
                            dgebak_obj->ihiref, dgebak_obj->scaleref, 
                            dgebak_obj->m, dgebak_obj->vref, dgebak_obj->ldv );

    /* Compute libflame's Lapacke o/p  */
    dgebak_obj->info = LAPACKE_dgebak( dgebak_obj->matrix_layout, dgebak_obj->job,
                            dgebak_obj->side, dgebak_obj->n, dgebak_obj->ilo,
                            dgebak_obj->ihi, dgebak_obj->scale, 
                            dgebak_obj->m, dgebak_obj->v, dgebak_obj->ldv);

    dgebak_obj->diff_scale =  computeDiff_d( dgebak_obj->n, 
                dgebak_obj->scale, dgebak_obj->scaleref );

    dgebak_obj->diff_v =  computeDiff_d( dgebak_obj->n*dgebak_obj->m, 
                dgebak_obj->v, dgebak_obj->vref );
}

TEST_F(dgebak_test, dgebak1) {
    EXPECT_NEAR(0.0, dgebak_obj->diff_v, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgebak_obj->diff_scale, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgebak_test, dgebak2) {
    EXPECT_NEAR(0.0, dgebak_obj->diff_v, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgebak_obj->diff_scale, LAPACKE_GTEST_THRESHOLD);
}
	
TEST_F(dgebak_test, dgebak3) {
    EXPECT_NEAR(0.0, dgebak_obj->diff_v, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgebak_obj->diff_scale, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgebak_test, dgebak4) {
    EXPECT_NEAR(0.0, dgebak_obj->diff_v, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgebak_obj->diff_scale, LAPACKE_GTEST_THRESHOLD);
}


/* Begin gebak_scomplex_parameters  class definition */
class gebak_scomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_v;
    float diff_ilo, diff_ihi, diff_scale;
    void *hModule, *dModule;
    
    lapack_complex_float *a, *aref;
    lapack_int lda;
    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char job;  // Must be 'N' or 'P' or 'S' or 'B'.
    char side; // Must be 'L' or 'R'.
    lapack_int n; // number of rows
    lapack_int ilo, iloref;
    lapack_int ihi, ihiref;
    float *scale, *scaleref; //  permutations / the scaling factor
    lapack_int m; // number of columns
    lapack_int ldv; // leading dimension of v

    /* Input/ Output parameters */
    lapack_complex_float *v, *vref; // contains the matrix of left or right eigenvectors

    /*Return Values */
    int info, inforef;

   public:
      gebak_scomplex_parameters (int matrix_layout_i, char job_i, char side_i,
            lapack_int n_i, lapack_int m_i, lapack_int lda_i );

      ~gebak_scomplex_parameters ();
};

/* Constructor definition  lapack_scomplex_parameters */
gebak_scomplex_parameters:: gebak_scomplex_parameters (int matrix_layout_i,
            char job_i,  char side_i, lapack_int n_i, lapack_int m_i, 
            lapack_int lda_i)
{
    matrix_layout = matrix_layout_i;
    job = job_i;
    side = side_i;
    n  = n_i;
    m = m_i;
    lda = lda_i;
    
    ilo = ihi = 0;
    hModule = NULL;
    dModule = NULL;

    if(matrix_layout==LAPACK_COL_MAJOR){
        ldv = n;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        ldv = m;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }
    
#if LAPACKE_TEST_VERBOSE
   printf(" \n gebak lapack_complex_float:  n: %d m: %d job: %c side: %c   \
lda: %d  matrix_layout: %d \n",  n, m, job, side,  
                                          lda, matrix_layout);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, n*lda );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &v, &vref, n*m );
    lapacke_gtest_alloc_float_buffer_pair (&scale, &scaleref, n);   

    if( (a==NULL) || (aref==NULL) ||  \
        (v==NULL) || (vref==NULL) || \
        (scale==NULL) || (scaleref==NULL) ){
       EXPECT_FALSE( true) << "gebak_scomplex_parameters object: malloc error. Exiting ";
       gebak_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_float_buffer_pair_with_constant(scale, scaleref, n, 0.0);
    lapacke_gtest_init_scomplex_buffer_pair_with_constant(v, vref, n*m, 0.0);
    
   } /* end of Constructor  */
    

/* Destructor definition  'gebak_scomplex_parameters' */
gebak_scomplex_parameters :: ~gebak_scomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   gebak_free();
} 

//  Test fixture class definition
class cgebak_test  : public  ::testing::Test {
public:
   gebak_scomplex_parameters  *cgebak_obj;
   void SetUp();  
   void TearDown () { delete cgebak_obj; }
};

void cgebak_test::SetUp()
{
    /* LAPACKE CGEBAK prototype */
    typedef int (*Fptr_NL_LAPACKE_cgebak) ( int matrix_layout, char job,
				char side, lapack_int n, lapack_int ilo, lapack_int ihi,
				const float* scale, lapack_int m, lapack_complex_float* v,
				lapack_int ldv  );
            
    /* LAPACKE DGEBAL prototype */
    typedef int (*Fptr_NL_LAPACKE_cgebal) ( int matrix_layout, char job,
				lapack_int n, lapack_complex_float* a, lapack_int lda, 
				lapack_int* ilo, lapack_int* ihi, float* scale);
                            
    Fptr_NL_LAPACKE_cgebal CGEBAL;      
    Fptr_NL_LAPACKE_cgebak CGEBAK;      

    cgebak_obj = new  gebak_scomplex_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_paramslist[idx].job,
                                         eig_paramslist[idx].side,
                                         eig_paramslist[idx].n,
                                         eig_paramslist[idx].m,
                                         eig_paramslist[idx].lda );
    idx = Circular_Increment_Index(idx);

    cgebak_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgebak_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgebak_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgebak_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CGEBAK = (Fptr_NL_LAPACKE_cgebak)dlsym(cgebak_obj->hModule, "LAPACKE_cgebak");
    ASSERT_TRUE(CGEBAK != NULL) << "failed to ppt the Netlib LAPACKE_cgebak symbol";

    CGEBAL = (Fptr_NL_LAPACKE_cgebal)dlsym(cgebak_obj->hModule, "LAPACKE_cgebal");
    ASSERT_TRUE(CGEBAL != NULL) << "failed to get the Netlib LAPACKE_cgebal symbol";
	
    /* invoke the LAPACKE_cgebal API to balance A */
    cgebak_obj->inforef = CGEBAL( cgebak_obj->matrix_layout, cgebak_obj->job,
                            cgebak_obj->n,cgebak_obj->aref, cgebak_obj->lda, 
                            &cgebak_obj->iloref, &cgebak_obj->ihiref, 
                            cgebak_obj->scaleref);
	
    cgebak_obj->info = LAPACKE_cgebal( cgebak_obj->matrix_layout, cgebak_obj->job,
                            cgebak_obj->n, cgebak_obj->a, cgebak_obj->lda, 
                            &cgebak_obj->ilo, &cgebak_obj->ihi, cgebak_obj->scale);
	
    cgebak_obj->diff_ilo = fabs(cgebak_obj->ilo - cgebak_obj->iloref);
    cgebak_obj->diff_ihi = fabs(cgebak_obj->ihi - cgebak_obj->ihiref);
	
    /* Compute the reference o/p by invoking Netlib-Lapack's API */
    cgebak_obj->inforef = CGEBAK( cgebak_obj->matrix_layout, cgebak_obj->job,
                            cgebak_obj->side, cgebak_obj->n, cgebak_obj->iloref,
                            cgebak_obj->ihiref, cgebak_obj->scaleref, 
                            cgebak_obj->m, cgebak_obj->vref, cgebak_obj->ldv );

    /* Compute libflame's Lapacke o/p  */
    cgebak_obj->info = LAPACKE_cgebak( cgebak_obj->matrix_layout, cgebak_obj->job,
                            cgebak_obj->side, cgebak_obj->n, cgebak_obj->ilo,
                            cgebak_obj->ihi, cgebak_obj->scale, 
                            cgebak_obj->m, cgebak_obj->v, cgebak_obj->ldv);

    cgebak_obj->diff_scale =  computeDiff_s( cgebak_obj->n, 
                cgebak_obj->scale, cgebak_obj->scaleref );

    cgebak_obj->diff_v =  computeDiff_c( cgebak_obj->n*cgebak_obj->m, 
                cgebak_obj->v, cgebak_obj->vref );
}

TEST_F(cgebak_test, cgebak1) {
    EXPECT_NEAR(0.0, cgebak_obj->diff_v, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgebak_obj->diff_scale, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgebak_test, cgebak2) {
    EXPECT_NEAR(0.0, cgebak_obj->diff_v, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgebak_obj->diff_scale, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgebak_test, cgebak3) {
    EXPECT_NEAR(0.0, cgebak_obj->diff_v, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgebak_obj->diff_scale, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgebak_test, cgebak4) {
    EXPECT_NEAR(0.0, cgebak_obj->diff_v, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgebak_obj->diff_scale, LAPACKE_GTEST_THRESHOLD);
}

/* Begin gebak_dcomplex_parameters  class definition */
class gebak_dcomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    double diff_v;
    double diff_ilo, diff_ihi, diff_scale;
    void *hModule, *dModule;
    
    lapack_complex_double *a, *aref;
    lapack_int lda;
    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char job;  // Must be 'N' or 'P' or 'S' or 'B'.
    char side; // Must be 'L' or 'R'.
    lapack_int n; // number of rows
    lapack_int ilo, iloref;
    lapack_int ihi, ihiref;
    double *scale, *scaleref; //  permutations / the scaling factor
    lapack_int m; // number of columns
    lapack_int ldv; // leading dimension of v

    /* Input/ Output parameters */
    lapack_complex_double *v, *vref; // contains the matrix of left or right eigenvectors

    /*Return Values */
    int info, inforef;

   public:
      gebak_dcomplex_parameters (int matrix_layout_i, char job_i, char side_i,
            lapack_int n_i, lapack_int m_i, lapack_int lda_i );

      ~gebak_dcomplex_parameters ();
};

/* Constructor definition  lapack_dcomplex_parameters */
gebak_dcomplex_parameters:: gebak_dcomplex_parameters (int matrix_layout_i,
            char job_i,  char side_i, lapack_int n_i, lapack_int m_i, 
            lapack_int lda_i)
{
    matrix_layout = matrix_layout_i;
    job = job_i;
    side = side_i;
    n  = n_i;
    m = m_i;
    lda = lda_i;
    
    ilo = ihi = 0;
    hModule = NULL;
    dModule = NULL;

    if(matrix_layout==LAPACK_COL_MAJOR){
        ldv = n;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        ldv = m;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }
    
#if LAPACKE_TEST_VERBOSE
   printf(" \n gebak lapack_complex_double:  n: %d m: %d job: %c side: %c   \
lda: %d  matrix_layout: %d \n",  n, m, job, side,  
                                          lda, matrix_layout);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, n*lda );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &v, &vref, n*m );
    lapacke_gtest_alloc_double_buffer_pair (&scale, &scaleref, n);   

    if( (a==NULL) || (aref==NULL) ||  \
        (v==NULL) || (vref==NULL) || \
        (scale==NULL) || (scaleref==NULL) ){
       EXPECT_FALSE( true) << "gebak_dcomplex_parameters object: malloc error. Exiting ";
       gebak_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_double_buffer_pair_with_constant(scale, scaleref, n, 0.0);
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant(v, vref, n*m, 0.0);
    
   } /* end of Constructor  */
    

/* Destructor definition  'gebak_dcomplex_parameters' */
gebak_dcomplex_parameters :: ~gebak_dcomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   gebak_free();
} 

//  Test fixture class definition
class zgebak_test  : public  ::testing::Test {
public:
   gebak_dcomplex_parameters  *zgebak_obj;
   void SetUp();  
   void TearDown () { delete zgebak_obj; }
};

void zgebak_test::SetUp()
{
    /* LAPACKE ZGEBAK prototype */
    typedef int (*Fptr_NL_LAPACKE_zgebak) ( int matrix_layout, char job,
				char side, lapack_int n, lapack_int ilo, lapack_int ihi,
				const double* scale, lapack_int m, lapack_complex_double* v,
				lapack_int ldv  );
            
    /* LAPACKE DGEBAL prototype */
    typedef int (*Fptr_NL_LAPACKE_zgebal) ( int matrix_layout, char job,
				lapack_int n, lapack_complex_double* a, lapack_int lda, 
				lapack_int* ilo, lapack_int* ihi, double* scale);
                            
    Fptr_NL_LAPACKE_zgebal ZGEBAL;      
    Fptr_NL_LAPACKE_zgebak ZGEBAK;      

    zgebak_obj = new  gebak_dcomplex_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_paramslist[idx].job,
                                         eig_paramslist[idx].side,
                                         eig_paramslist[idx].n,
                                         eig_paramslist[idx].m,
                                         eig_paramslist[idx].lda );
    idx = Circular_Increment_Index(idx);

    zgebak_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgebak_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgebak_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgebak_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZGEBAK = (Fptr_NL_LAPACKE_zgebak)dlsym(zgebak_obj->hModule, "LAPACKE_zgebak");
    ASSERT_TRUE(ZGEBAK != NULL) << "failed to ppt the Netlib LAPACKE_zgebak symbol";

    ZGEBAL = (Fptr_NL_LAPACKE_zgebal)dlsym(zgebak_obj->hModule, "LAPACKE_zgebal");
    ASSERT_TRUE(ZGEBAL != NULL) << "failed to get the Netlib LAPACKE_zgebal symbol";
	
    /* invoke the LAPACKE_zgebal API to balance A */
    zgebak_obj->inforef = ZGEBAL( zgebak_obj->matrix_layout, zgebak_obj->job,
                            zgebak_obj->n,zgebak_obj->aref, zgebak_obj->lda, 
                            &zgebak_obj->iloref, &zgebak_obj->ihiref, 
                            zgebak_obj->scaleref);
	
    zgebak_obj->info = LAPACKE_zgebal( zgebak_obj->matrix_layout, zgebak_obj->job,
                            zgebak_obj->n, zgebak_obj->a, zgebak_obj->lda, 
                            &zgebak_obj->ilo, &zgebak_obj->ihi, zgebak_obj->scale);
	
    zgebak_obj->diff_ilo = fabs(zgebak_obj->ilo - zgebak_obj->iloref);
    zgebak_obj->diff_ihi = fabs(zgebak_obj->ihi - zgebak_obj->ihiref);
	
    /* Compute the reference o/p by invoking Netlib-Lapack's API */
    zgebak_obj->inforef = ZGEBAK( zgebak_obj->matrix_layout, zgebak_obj->job,
                            zgebak_obj->side, zgebak_obj->n, zgebak_obj->iloref,
                            zgebak_obj->ihiref, zgebak_obj->scaleref, 
                            zgebak_obj->m, zgebak_obj->vref, zgebak_obj->ldv );

    /* Compute libflame's Lapacke o/p  */
    zgebak_obj->info = LAPACKE_zgebak( zgebak_obj->matrix_layout, zgebak_obj->job,
                            zgebak_obj->side, zgebak_obj->n, zgebak_obj->ilo,
                            zgebak_obj->ihi, zgebak_obj->scale, 
                            zgebak_obj->m, zgebak_obj->v, zgebak_obj->ldv);

    zgebak_obj->diff_scale =  computeDiff_d( zgebak_obj->n, 
                zgebak_obj->scale, zgebak_obj->scaleref );

    zgebak_obj->diff_v =  computeDiff_z( zgebak_obj->n*zgebak_obj->m, 
                zgebak_obj->v, zgebak_obj->vref );
}

TEST_F(zgebak_test, zgebak1) {
    EXPECT_NEAR(0.0, zgebak_obj->diff_v, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgebak_obj->diff_scale, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgebak_test, zgebak2) {
    EXPECT_NEAR(0.0, zgebak_obj->diff_v, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgebak_obj->diff_scale, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgebak_test, zgebak3) {
    EXPECT_NEAR(0.0, zgebak_obj->diff_v, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgebak_obj->diff_scale, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgebak_test, zgebak4) {
    EXPECT_NEAR(0.0, zgebak_obj->diff_v, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgebak_obj->diff_scale, LAPACKE_GTEST_THRESHOLD);
}