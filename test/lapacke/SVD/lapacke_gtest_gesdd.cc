#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"
#define LAPACKE_TEST_VERBOSE  (1)
#define gesdd_free() \
    if (a!=NULL)        free(a); \
    if (aref!=NULL)     free(aref); \
    if (s!=NULL)        free(s); \
    if (sref!=NULL)     free(sref); \
    if (u!=NULL)        free(u); \
    if (uref!=NULL)     free(uref); \
    if (vt!=NULL)        free(vt); \
    if (vtref!=NULL)     free(vtref)


// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin gesdd_float_common_parameters  class definition */
class gesdd_float_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_a, diff_s, diff_u, diff_vt;
    float threshold;
    void *hModule, *dModule;
	int min_mn;
    
    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char jobz; // Must be 'A', 'S', 'O', or 'N'. 

    lapack_int m; // The number of rows of the matrix A
    lapack_int n; // The number of columns of the matrix A 
    lapack_int lda; //  The leading dimension of a

    lapack_int ldu; // leading dimension of the output array u . ldu≥ fla_max(1, m)
    lapack_int ldvt; // The leading dimension of the output array vt . ldvt≥ fla_max(1, p) 
    /* Input / Output parameters */
    float* a, *aref; // contains m-by-n matrix A.

    /* Output parameters */
    // Below buffers contain the o/p orthogonal/unitary matrices
    float* u, *uref;
    float* vt, *vtref;
    float* s, *sref;
    
    /*Return Values */
    int info, inforef;

   public:
      gesdd_float_parameters (int matrix_layout_i, char jobz,
                                  lapack_int m, lapack_int n);
      ~gesdd_float_parameters ();
};

/* Constructor definition  gesdd float_common_parameters */
gesdd_float_parameters:: gesdd_float_parameters (int matrix_layout_i,
                    char jobz_i, lapack_int m_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    jobz = jobz_i;
    n = n_i;
    m = m_i;
	min_mn = (m<n)?m:n;
    
	switch (jobz)
	{
		case 'A': 
			ldu = m;
			ldvt = n;
			break;
		case 'N': 
			ldu = 1;
			ldvt = 1;
			break;
		case 'O':
		   if (m >= n)
		   {
		    ldu = 1;
			ldvt = n;
		   } else {
			ldu = m;
			ldvt = 1;
		   }
		   break;
		case 'S':
		   ldu = m;
		   ldvt = n;
		   if ((m >= n) && (matrix_layout == LAPACK_COL_MAJOR))
		   {
				  ldu = n;
				  ldvt = m;
		   }
		   break;
		default:
		   printf( "\n invalid jobz option supplied \n");
		   break;

	}
	
    if (matrix_layout == LAPACK_ROW_MAJOR)
    {
        lda = n;
    }   else
    {
        lda = m;
    }

    // Initialize with default values
    hModule = NULL;
    dModule = NULL;
    diff_a = 0;
    diff_s = 0;
    diff_u = 0;
    diff_vt = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n gesdd float: matrix_layout: %d   jobz: %c \t \
m: %d \t n: %d  \t lda: %d  \t ldvt: %d \t ldu: %d \n",
matrix_layout, jobz_i, m_i, n_i, lda, ldvt, ldu);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, m*n );
    lapacke_gtest_alloc_float_buffer_pair( &s, &sref, min_mn );
    lapacke_gtest_alloc_float_buffer_pair( &u, &uref, m*m );
    lapacke_gtest_alloc_float_buffer_pair( &vt, &vtref, n*n );

    if( (a==NULL) || (aref==NULL) ||  \
        (s==NULL) || (sref==NULL) || \
        (u==NULL) || (uref==NULL) || \
        (vt==NULL) || (vtref==NULL) ){
       EXPECT_FALSE( true) << "gesdd_float_parameters object: malloc error. Exiting ";
       gesdd_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, m*n );
    lapacke_gtest_init_float_buffer_pair_with_constant( u, uref, m*m, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant( vt, vtref, n*n, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant( s, sref, min_mn, 0.0);

   } /* end of Constructor  */

/* Destructor definition  'gesdd_float_common_parameters' */
gesdd_float_parameters :: ~gesdd_float_parameters ()
{
   /* De-Allocate memory for the input matrices */
   gesdd_free();
} 

//  Test fixture class definition
class sgesdd_test  : public  ::testing::Test {
public:
   gesdd_float_parameters  *sgesdd_obj;
   void SetUp();  
   void TearDown () { delete sgesdd_obj; }
};

void sgesdd_test::SetUp()
{
    /* LAPACKE sgesdd prototype */
    typedef int (*Fptr_NL_LAPACKE_sgesdd) ( int matrix_layout, char jobz,
		lapack_int m, lapack_int n, float* a, lapack_int lda, float* s,
		float* u, lapack_int ldu, float* vt, lapack_int ldvt);
            
    Fptr_NL_LAPACKE_sgesdd SGESDD;

    sgesdd_obj = new  gesdd_float_parameters( svd_paramslist[idx].matrix_layout,
                                         svd_paramslist[idx].jobu_gesvd,
                                         svd_paramslist[idx].m,
                                         svd_paramslist[idx].n );

    idx = Circular_Increment_Index(idx);
    sgesdd_obj->threshold = svd_paramslist[idx].svd_threshold;
    sgesdd_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgesdd_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgesdd_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgesdd_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SGESDD = (Fptr_NL_LAPACKE_sgesdd)dlsym(sgesdd_obj->hModule, "LAPACKE_sgesdd");
    ASSERT_TRUE(SGESDD != NULL) << "failed to ppt the Netlib LAPACKE_sgesdd symbol";

    /* Compute libflame's Lapacke o/p  */
    sgesdd_obj->info = LAPACKE_sgesdd(  sgesdd_obj->matrix_layout,
                                    sgesdd_obj->jobz,
                                    sgesdd_obj->m,
                                    sgesdd_obj->n,
                                    sgesdd_obj->a,
                                    sgesdd_obj->lda,
                                    sgesdd_obj->s,
                                    sgesdd_obj->u,
                                    sgesdd_obj->ldu,
                                    sgesdd_obj->vt,
                                    sgesdd_obj->ldvt
                                    );

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    sgesdd_obj->inforef = SGESDD( sgesdd_obj->matrix_layout,
                                    sgesdd_obj->jobz,
                                    sgesdd_obj->m,
                                    sgesdd_obj->n,
                                    sgesdd_obj->aref,
                                    sgesdd_obj->lda,
                                    sgesdd_obj->sref,
                                    sgesdd_obj->uref,
                                    sgesdd_obj->ldu,
                                    sgesdd_obj->vtref,
                                    sgesdd_obj->ldvt
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    sgesdd_obj->diff_a =  computeDiff_s( (sgesdd_obj->m)*(sgesdd_obj->n), 
                sgesdd_obj->a, sgesdd_obj->aref );

    sgesdd_obj->diff_u =  computeDiff_s( (sgesdd_obj->m)*(sgesdd_obj->m), 
                sgesdd_obj->u, sgesdd_obj->uref );

    sgesdd_obj->diff_vt =  computeDiff_s( (sgesdd_obj->n)*(sgesdd_obj->n), 
                sgesdd_obj->vt, sgesdd_obj->vtref );

    sgesdd_obj->diff_s =  computeDiff_s( sgesdd_obj->min_mn, 
                sgesdd_obj->s, sgesdd_obj->sref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n gesdd float: \n diff_a: %f \n diff_s: %f \n \
diff_u: %f \n diff_vt: %f \n  info: %d \t inforef: %d\n ",
       sgesdd_obj->diff_a, sgesdd_obj->diff_s, sgesdd_obj->diff_u,
       sgesdd_obj->diff_vt, sgesdd_obj->info, sgesdd_obj->inforef);
#endif
}

TEST_F(sgesdd_test, sgesdd_1) {
    EXPECT_NEAR(0.0, sgesdd_obj->diff_a, sgesdd_obj->threshold);
    EXPECT_NEAR(0.0, sgesdd_obj->diff_s, sgesdd_obj->threshold);
    EXPECT_NEAR(0.0, sgesdd_obj->diff_u, sgesdd_obj->threshold);
    EXPECT_NEAR(0.0, sgesdd_obj->diff_vt, sgesdd_obj->threshold);
}

TEST_F(sgesdd_test, sgesdd_2) {
    EXPECT_NEAR(0.0, sgesdd_obj->diff_a, sgesdd_obj->threshold);
    EXPECT_NEAR(0.0, sgesdd_obj->diff_s, sgesdd_obj->threshold);
    EXPECT_NEAR(0.0, sgesdd_obj->diff_u, sgesdd_obj->threshold);
    EXPECT_NEAR(0.0, sgesdd_obj->diff_vt, sgesdd_obj->threshold);
}
TEST_F(sgesdd_test, sgesdd_3) {
    EXPECT_NEAR(0.0, sgesdd_obj->diff_a, sgesdd_obj->threshold);
    EXPECT_NEAR(0.0, sgesdd_obj->diff_s, sgesdd_obj->threshold);
    EXPECT_NEAR(0.0, sgesdd_obj->diff_u, sgesdd_obj->threshold);
    EXPECT_NEAR(0.0, sgesdd_obj->diff_vt, sgesdd_obj->threshold);
}
TEST_F(sgesdd_test, sgesdd_4) {
    EXPECT_NEAR(0.0, sgesdd_obj->diff_a, sgesdd_obj->threshold);
    EXPECT_NEAR(0.0, sgesdd_obj->diff_s, sgesdd_obj->threshold);
    EXPECT_NEAR(0.0, sgesdd_obj->diff_u, sgesdd_obj->threshold);
    EXPECT_NEAR(0.0, sgesdd_obj->diff_vt, sgesdd_obj->threshold);
}

/* Begin gesdd_double_common_parameters  class definition */
class gesdd_double_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    double diff_a, diff_s, diff_u, diff_vt;
    float threshold;
    void *hModule, *dModule;
	int min_mn;
    
    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char jobz; // Must be 'A', 'S', 'O', or 'N'. 

    lapack_int m; // The number of rows of the matrix A
    lapack_int n; // The number of columns of the matrix A 
    lapack_int lda; //  The leading dimension of a

    lapack_int ldu; // leading dimension of the output array u . ldu≥ fla_max(1, m)
    lapack_int ldvt; // The leading dimension of the output array vt . ldvt≥ fla_max(1, p) 
    /* Input / Output parameters */
    double* a, *aref; // contains m-by-n matrix A.

    /* Output parameters */
    // Below buffers contain the o/p orthogonal/unitary matrices
    double* u, *uref;
    double* vt, *vtref;
    double* s, *sref;
    
    /*Return Values */
    int info, inforef;

   public:
      gesdd_double_parameters (int matrix_layout_i, char jobz,
                                  lapack_int m, lapack_int n);
      ~gesdd_double_parameters ();
};

/* Constructor definition  gesdd double_common_parameters */
gesdd_double_parameters:: gesdd_double_parameters (int matrix_layout_i,
                    char jobz_i, lapack_int m_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    jobz = jobz_i;
    n = n_i;
    m = m_i;
	min_mn = (m<n)?m:n;
    
	switch (jobz)
	{
		case 'A': 
			ldu = m;
			ldvt = n;
			break;
		case 'N': 
			ldu = 1;
			ldvt = 1;
			break;
		case 'O':
		   if (m >= n)
		   {
		    ldu = 1;
			ldvt = n;
		   } else {
			ldu = m;
			ldvt = 1;
		   }
		   break;
		case 'S':
		   ldu = m;
		   ldvt = n;
		   if ((m >= n) && (matrix_layout == LAPACK_COL_MAJOR))
		   {
				  ldu = n;
				  ldvt = m;
		   }
		   break;
		default:
		   printf( "\n invalid jobz option supplied \n");
		   break;

	}
	
    if (matrix_layout == LAPACK_ROW_MAJOR)
    {
        lda = n;
    }   else
    {
        lda = m;
    }

    // Initialize with default values
    hModule = NULL;
    dModule = NULL;
    diff_a = 0;
    diff_s = 0;
    diff_u = 0;
    diff_vt = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n gesdd double: matrix_layout: %d   jobz: %c \t \
m: %d \t n: %d  \t lda: %d  \t ldvt: %d \t ldu: %d \n",
matrix_layout, jobz_i, m_i, n_i, lda, ldvt, ldu);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, m*n );
    lapacke_gtest_alloc_double_buffer_pair( &s, &sref, min_mn );
    lapacke_gtest_alloc_double_buffer_pair( &u, &uref, m*m );
    lapacke_gtest_alloc_double_buffer_pair( &vt, &vtref, n*n );

    if( (a==NULL) || (aref==NULL) ||  \
        (s==NULL) || (sref==NULL) || \
        (u==NULL) || (uref==NULL) || \
        (vt==NULL) || (vtref==NULL) ){
       EXPECT_FALSE( true) << "gesdd_double_parameters object: malloc error. Exiting ";
       gesdd_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, m*n );
    lapacke_gtest_init_double_buffer_pair_with_constant( u, uref, m*m, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant( vt, vtref, n*n, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant( s, sref, min_mn, 0.0);

   } /* end of Constructor  */

/* Destructor definition  'gesdd_double_common_parameters' */
gesdd_double_parameters :: ~gesdd_double_parameters ()
{
   /* De-Allocate memory for the input matrices */
   gesdd_free();
} 

//  Test fixture class definition
class dgesdd_test  : public  ::testing::Test {
public:
   gesdd_double_parameters  *dgesdd_obj;
   void SetUp();  
   void TearDown () { delete dgesdd_obj; }
};

void dgesdd_test::SetUp()
{
    /* LAPACKE dgesdd prototype */
    typedef int (*Fptr_NL_LAPACKE_dgesdd) ( int matrix_layout, char jobz,
		lapack_int m, lapack_int n, double* a, lapack_int lda, double* s,
		double* u, lapack_int ldu, double* vt, lapack_int ldvt);
            
    Fptr_NL_LAPACKE_dgesdd DGESDD;

    dgesdd_obj = new  gesdd_double_parameters( svd_paramslist[idx].matrix_layout,
                                         svd_paramslist[idx].jobu_gesvd,
                                         svd_paramslist[idx].m,
                                         svd_paramslist[idx].n );

    idx = Circular_Increment_Index(idx);
    dgesdd_obj->threshold = svd_paramslist[idx].svd_threshold;
    dgesdd_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgesdd_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgesdd_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgesdd_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DGESDD = (Fptr_NL_LAPACKE_dgesdd)dlsym(dgesdd_obj->hModule, "LAPACKE_dgesdd");
    ASSERT_TRUE(DGESDD != NULL) << "failed to ppt the Netlib LAPACKE_dgesdd symbol";

    /* Compute libflame's Lapacke o/p  */
    dgesdd_obj->info = LAPACKE_dgesdd(  dgesdd_obj->matrix_layout,
                                    dgesdd_obj->jobz,
                                    dgesdd_obj->m,
                                    dgesdd_obj->n,
                                    dgesdd_obj->a,
                                    dgesdd_obj->lda,
                                    dgesdd_obj->s,
                                    dgesdd_obj->u,
                                    dgesdd_obj->ldu,
                                    dgesdd_obj->vt,
                                    dgesdd_obj->ldvt
                                    );

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    dgesdd_obj->inforef = DGESDD( dgesdd_obj->matrix_layout,
                                    dgesdd_obj->jobz,
                                    dgesdd_obj->m,
                                    dgesdd_obj->n,
                                    dgesdd_obj->aref,
                                    dgesdd_obj->lda,
                                    dgesdd_obj->sref,
                                    dgesdd_obj->uref,
                                    dgesdd_obj->ldu,
                                    dgesdd_obj->vtref,
                                    dgesdd_obj->ldvt
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    dgesdd_obj->diff_a =  computeDiff_d( (dgesdd_obj->m)*(dgesdd_obj->n), 
                dgesdd_obj->a, dgesdd_obj->aref );

    dgesdd_obj->diff_u =  computeDiff_d( (dgesdd_obj->m)*(dgesdd_obj->m), 
                dgesdd_obj->u, dgesdd_obj->uref );

    dgesdd_obj->diff_vt =  computeDiff_d( (dgesdd_obj->n)*(dgesdd_obj->n), 
                dgesdd_obj->vt, dgesdd_obj->vtref );

    dgesdd_obj->diff_s =  computeDiff_d( dgesdd_obj->min_mn, 
                dgesdd_obj->s, dgesdd_obj->sref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n gesdd double: \n diff_a: %f \n diff_s: %f \n \
diff_u: %f \n diff_vt: %f \n  info: %d \t inforef: %d\n ",
       dgesdd_obj->diff_a, dgesdd_obj->diff_s, dgesdd_obj->diff_u,
       dgesdd_obj->diff_vt, dgesdd_obj->info, dgesdd_obj->inforef);
#endif
}

TEST_F(dgesdd_test, dgesdd_1) {
    EXPECT_NEAR(0.0, dgesdd_obj->diff_a, dgesdd_obj->threshold);
    EXPECT_NEAR(0.0, dgesdd_obj->diff_s, dgesdd_obj->threshold);
    EXPECT_NEAR(0.0, dgesdd_obj->diff_u, dgesdd_obj->threshold);
    EXPECT_NEAR(0.0, dgesdd_obj->diff_vt, dgesdd_obj->threshold);
}

TEST_F(dgesdd_test, dgesdd_2) {
    EXPECT_NEAR(0.0, dgesdd_obj->diff_a, dgesdd_obj->threshold);
    EXPECT_NEAR(0.0, dgesdd_obj->diff_s, dgesdd_obj->threshold);
    EXPECT_NEAR(0.0, dgesdd_obj->diff_u, dgesdd_obj->threshold);
    EXPECT_NEAR(0.0, dgesdd_obj->diff_vt, dgesdd_obj->threshold);
}
TEST_F(dgesdd_test, dgesdd_3) {
    EXPECT_NEAR(0.0, dgesdd_obj->diff_a, dgesdd_obj->threshold);
    EXPECT_NEAR(0.0, dgesdd_obj->diff_s, dgesdd_obj->threshold);
    EXPECT_NEAR(0.0, dgesdd_obj->diff_u, dgesdd_obj->threshold);
    EXPECT_NEAR(0.0, dgesdd_obj->diff_vt, dgesdd_obj->threshold);
}
TEST_F(dgesdd_test, dgesdd_4) {
    EXPECT_NEAR(0.0, dgesdd_obj->diff_a, dgesdd_obj->threshold);
    EXPECT_NEAR(0.0, dgesdd_obj->diff_s, dgesdd_obj->threshold);
    EXPECT_NEAR(0.0, dgesdd_obj->diff_u, dgesdd_obj->threshold);
    EXPECT_NEAR(0.0, dgesdd_obj->diff_vt, dgesdd_obj->threshold);
}

/* Begin gesdd_scomplex_common_parameters  class definition */
class gesdd_scomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_a, diff_s, diff_u, diff_vt;
    float threshold;
    void *hModule, *dModule;
	int min_mn;
    
    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char jobz; // Must be 'A', 'S', 'O', or 'N'. 

    lapack_int m; // The number of rows of the matrix A
    lapack_int n; // The number of columns of the matrix A 
    lapack_int lda; //  The leading dimension of a

    lapack_int ldu; // leading dimension of the output array u . ldu≥ fla_max(1, m)
    lapack_int ldvt; // The leading dimension of the output array vt . ldvt≥ fla_max(1, p) 
    /* Input / Output parameters */
    lapack_complex_float* a, *aref; // contains m-by-n matrix A.

    /* Output parameters */
    // Below buffers contain the o/p orthogonal/unitary matrices
    lapack_complex_float* u, *uref;
    lapack_complex_float* vt, *vtref;
    float* s, *sref;
    
    /*Return Values */
    int info, inforef;

   public:
      gesdd_scomplex_parameters (int matrix_layout_i, char jobz,
                                  lapack_int m, lapack_int n);
      ~gesdd_scomplex_parameters ();
};

/* Constructor definition  gesdd lapack_complex_float_common_parameters */
gesdd_scomplex_parameters:: gesdd_scomplex_parameters (int matrix_layout_i,
                    char jobz_i, lapack_int m_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    jobz = jobz_i;
    n = n_i;
    m = m_i;
	min_mn = (m<n)?m:n;
    
	switch (jobz)
	{
		case 'A': 
			ldu = m;
			ldvt = n;
			break;
		case 'N': 
			ldu = 1;
			ldvt = 1;
			break;
		case 'O':
		   if (m >= n)
		   {
		    ldu = 1;
			ldvt = n;
		   } else {
			ldu = m;
			ldvt = 1;
		   }
		   break;
		case 'S':
		   ldu = m;
		   ldvt = n;
		   if ((m >= n) && (matrix_layout == LAPACK_COL_MAJOR))
		   {
				  ldu = n;
				  ldvt = m;
		   }
		   break;
		default:
		   printf( "\n invalid jobz option supplied \n");
		   break;

	}
	
    if (matrix_layout == LAPACK_ROW_MAJOR)
    {
        lda = n;
    }   else
    {
        lda = m;
    }

    // Initialize with default values
    hModule = NULL;
    dModule = NULL;
    diff_a = 0;
    diff_s = 0;
    diff_u = 0;
    diff_vt = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n gesdd lapack_complex_float: matrix_layout: %d   jobz: %c \t \
m: %d \t n: %d  \t lda: %d  \t ldvt: %d \t ldu: %d \n",
matrix_layout, jobz_i, m_i, n_i, lda, ldvt, ldu);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, m*n );
    lapacke_gtest_alloc_float_buffer_pair( &s, &sref, min_mn );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &u, &uref, m*m );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &vt, &vtref, n*n );

    if( (a==NULL) || (aref==NULL) ||  \
        (s==NULL) || (sref==NULL) || \
        (u==NULL) || (uref==NULL) || \
        (vt==NULL) || (vtref==NULL) ){
       EXPECT_FALSE( true) << "gesdd_scomplex_parameters object: malloc error. Exiting ";
       gesdd_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, m*n );
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( u, uref, m*m, 0.0);
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( vt, vtref, n*n, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant( s, sref, min_mn, 0.0);

   } /* end of Constructor  */

/* Destructor definition  'gesdd_scomplex_common_parameters' */
gesdd_scomplex_parameters :: ~gesdd_scomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   gesdd_free();
} 

//  Test fixture class definition
class cgesdd_test  : public  ::testing::Test {
public:
   gesdd_scomplex_parameters  *cgesdd_obj;
   void SetUp();  
   void TearDown () { delete cgesdd_obj; }
};

void cgesdd_test::SetUp()
{
    /* LAPACKE cgesdd prototype */
    typedef int (*Fptr_NL_LAPACKE_cgesdd) ( int matrix_layout, char jobz,
	lapack_int m, lapack_int n, lapack_complex_float* a, lapack_int lda,
	float* s, lapack_complex_float* u, lapack_int ldu,
	lapack_complex_float* vt, lapack_int ldvt);
            
    Fptr_NL_LAPACKE_cgesdd CGESDD;

    cgesdd_obj = new  gesdd_scomplex_parameters( svd_paramslist[idx].matrix_layout,
                                         svd_paramslist[idx].jobu_gesvd,
                                         svd_paramslist[idx].m,
                                         svd_paramslist[idx].n );

    idx = Circular_Increment_Index(idx);
    cgesdd_obj->threshold = svd_paramslist[idx].svd_threshold;
    cgesdd_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgesdd_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgesdd_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgesdd_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CGESDD = (Fptr_NL_LAPACKE_cgesdd)dlsym(cgesdd_obj->hModule, "LAPACKE_cgesdd");
    ASSERT_TRUE(CGESDD != NULL) << "failed to ppt the Netlib LAPACKE_cgesdd symbol";

    /* Compute libflame's Lapacke o/p  */
    cgesdd_obj->info = LAPACKE_cgesdd(  cgesdd_obj->matrix_layout,
                                    cgesdd_obj->jobz,
                                    cgesdd_obj->m,
                                    cgesdd_obj->n,
                                    cgesdd_obj->a,
                                    cgesdd_obj->lda,
                                    cgesdd_obj->s,
                                    cgesdd_obj->u,
                                    cgesdd_obj->ldu,
                                    cgesdd_obj->vt,
                                    cgesdd_obj->ldvt
                                    );

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    cgesdd_obj->inforef = CGESDD( cgesdd_obj->matrix_layout,
                                    cgesdd_obj->jobz,
                                    cgesdd_obj->m,
                                    cgesdd_obj->n,
                                    cgesdd_obj->aref,
                                    cgesdd_obj->lda,
                                    cgesdd_obj->sref,
                                    cgesdd_obj->uref,
                                    cgesdd_obj->ldu,
                                    cgesdd_obj->vtref,
                                    cgesdd_obj->ldvt
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    cgesdd_obj->diff_a =  computeDiff_c( (cgesdd_obj->m)*(cgesdd_obj->n), 
                cgesdd_obj->a, cgesdd_obj->aref );

    cgesdd_obj->diff_u =  computeDiff_c( (cgesdd_obj->m)*(cgesdd_obj->m), 
                cgesdd_obj->u, cgesdd_obj->uref );

    cgesdd_obj->diff_vt =  computeDiff_c( (cgesdd_obj->n)*(cgesdd_obj->n), 
                cgesdd_obj->vt, cgesdd_obj->vtref );

    cgesdd_obj->diff_s =  computeDiff_s( cgesdd_obj->min_mn, 
                cgesdd_obj->s, cgesdd_obj->sref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n gesdd lapack_complex_float: \n diff_a: %f \n diff_s: %f \n \
diff_u: %f \n diff_vt: %f \n  info: %d \t inforef: %d\n ",
       cgesdd_obj->diff_a, cgesdd_obj->diff_s, cgesdd_obj->diff_u,
       cgesdd_obj->diff_vt, cgesdd_obj->info, cgesdd_obj->inforef);
#endif
}

TEST_F(cgesdd_test, cgesdd_1) {
    EXPECT_NEAR(0.0, cgesdd_obj->diff_a, cgesdd_obj->threshold);
    EXPECT_NEAR(0.0, cgesdd_obj->diff_s, cgesdd_obj->threshold);
    EXPECT_NEAR(0.0, cgesdd_obj->diff_u, cgesdd_obj->threshold);
    EXPECT_NEAR(0.0, cgesdd_obj->diff_vt, cgesdd_obj->threshold);
}

TEST_F(cgesdd_test, cgesdd_2) {
    EXPECT_NEAR(0.0, cgesdd_obj->diff_a, cgesdd_obj->threshold);
    EXPECT_NEAR(0.0, cgesdd_obj->diff_s, cgesdd_obj->threshold);
    EXPECT_NEAR(0.0, cgesdd_obj->diff_u, cgesdd_obj->threshold);
    EXPECT_NEAR(0.0, cgesdd_obj->diff_vt, cgesdd_obj->threshold);
}
TEST_F(cgesdd_test, cgesdd_3) {
    EXPECT_NEAR(0.0, cgesdd_obj->diff_a, cgesdd_obj->threshold);
    EXPECT_NEAR(0.0, cgesdd_obj->diff_s, cgesdd_obj->threshold);
    EXPECT_NEAR(0.0, cgesdd_obj->diff_u, cgesdd_obj->threshold);
    EXPECT_NEAR(0.0, cgesdd_obj->diff_vt, cgesdd_obj->threshold);
}
TEST_F(cgesdd_test, cgesdd_4) {
    EXPECT_NEAR(0.0, cgesdd_obj->diff_a, cgesdd_obj->threshold);
    EXPECT_NEAR(0.0, cgesdd_obj->diff_s, cgesdd_obj->threshold);
    EXPECT_NEAR(0.0, cgesdd_obj->diff_u, cgesdd_obj->threshold);
    EXPECT_NEAR(0.0, cgesdd_obj->diff_vt, cgesdd_obj->threshold);
}

/* Begin gesdd_dcomplex_common_parameters  class definition */
class gesdd_dcomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    double diff_a, diff_s, diff_u, diff_vt;
    double threshold;
    void *hModule, *dModule;
	int min_mn;
    
    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char jobz; // Must be 'A', 'S', 'O', or 'N'. 

    lapack_int m; // The number of rows of the matrix A
    lapack_int n; // The number of columns of the matrix A 
    lapack_int lda; //  The leading dimension of a

    lapack_int ldu; // leading dimension of the output array u . ldu≥ fla_max(1, m)
    lapack_int ldvt; // The leading dimension of the output array vt . ldvt≥ fla_max(1, p) 
    /* Input / Output parameters */
    lapack_complex_double* a, *aref; // contains m-by-n matrix A.

    /* Output parameters */
    // Below buffers contain the o/p orthogonal/unitary matrices
    lapack_complex_double* u, *uref;
    lapack_complex_double* vt, *vtref;
    double* s, *sref;
    
    /*Return Values */
    int info, inforef;

   public:
      gesdd_dcomplex_parameters (int matrix_layout_i, char jobz,
                                  lapack_int m, lapack_int n);
      ~gesdd_dcomplex_parameters ();
};

/* Constructor definition  gesdd lapack_complex_double_common_parameters */
gesdd_dcomplex_parameters:: gesdd_dcomplex_parameters (int matrix_layout_i,
                    char jobz_i, lapack_int m_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    jobz = jobz_i;
    n = n_i;
    m = m_i;
	min_mn = (m<n)?m:n;
    
	switch (jobz)
	{
		case 'A': 
			ldu = m;
			ldvt = n;
			break;
		case 'N': 
			ldu = 1;
			ldvt = 1;
			break;
		case 'O':
		   if (m >= n)
		   {
		    ldu = 1;
			ldvt = n;
		   } else {
			ldu = m;
			ldvt = 1;
		   }
		   break;
		case 'S':
		   ldu = m;
		   ldvt = n;
		   if ((m >= n) && (matrix_layout == LAPACK_COL_MAJOR))
		   {
				  ldu = n;
				  ldvt = m;
		   }
		   break;
		default:
		   printf( "\n invalid jobz option supplied \n");
		   break;

	}
	
    if (matrix_layout == LAPACK_ROW_MAJOR)
    {
        lda = n;
    }   else
    {
        lda = m;
    }

    // Initialize with default values
    hModule = NULL;
    dModule = NULL;
    diff_a = 0;
    diff_s = 0;
    diff_u = 0;
    diff_vt = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n gesdd lapack_complex_double: matrix_layout: %d   jobz: %c \t \
m: %d \t n: %d  \t lda: %d  \t ldvt: %d \t ldu: %d \n",
matrix_layout, jobz_i, m_i, n_i, lda, ldvt, ldu);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, m*n );
    lapacke_gtest_alloc_double_buffer_pair( &s, &sref, min_mn );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &u, &uref, m*m );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &vt, &vtref, n*n );

    if( (a==NULL) || (aref==NULL) ||  \
        (s==NULL) || (sref==NULL) || \
        (u==NULL) || (uref==NULL) || \
        (vt==NULL) || (vtref==NULL) ){
       EXPECT_FALSE( true) << "gesdd_dcomplex_parameters object: malloc error. Exiting ";
       gesdd_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, m*n );
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( u, uref, m*m, 0.0);
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( vt, vtref, n*n, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant( s, sref, min_mn, 0.0);

   } /* end of Constructor  */

/* Destructor definition  'gesdd_dcomplex_common_parameters' */
gesdd_dcomplex_parameters :: ~gesdd_dcomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   gesdd_free();
} 

//  Test fixture class definition
class zgesdd_test  : public  ::testing::Test {
public:
   gesdd_dcomplex_parameters  *zgesdd_obj;
   void SetUp();  
   void TearDown () { delete zgesdd_obj; }
};

void zgesdd_test::SetUp()
{
    /* LAPACKE zgesdd prototype */
    typedef int (*Fptr_NL_LAPACKE_zgesdd) ( int matrix_layout, char jobz,
	lapack_int m, lapack_int n, lapack_complex_double* a, lapack_int lda,
	double* s, lapack_complex_double* u, lapack_int ldu,
	lapack_complex_double* vt, lapack_int ldvt);
            
    Fptr_NL_LAPACKE_zgesdd ZGESDD;

    zgesdd_obj = new  gesdd_dcomplex_parameters( svd_paramslist[idx].matrix_layout,
                                         svd_paramslist[idx].jobu_gesvd,
                                         svd_paramslist[idx].m,
                                         svd_paramslist[idx].n );

    idx = Circular_Increment_Index(idx);
    zgesdd_obj->threshold = svd_paramslist[idx].svd_threshold;
    zgesdd_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgesdd_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgesdd_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgesdd_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZGESDD = (Fptr_NL_LAPACKE_zgesdd)dlsym(zgesdd_obj->hModule, "LAPACKE_zgesdd");
    ASSERT_TRUE(ZGESDD != NULL) << "failed to ppt the Netlib LAPACKE_zgesdd symbol";

    /* Compute libflame's Lapacke o/p  */
    zgesdd_obj->info = LAPACKE_zgesdd(  zgesdd_obj->matrix_layout,
                                    zgesdd_obj->jobz,
                                    zgesdd_obj->m,
                                    zgesdd_obj->n,
                                    zgesdd_obj->a,
                                    zgesdd_obj->lda,
                                    zgesdd_obj->s,
                                    zgesdd_obj->u,
                                    zgesdd_obj->ldu,
                                    zgesdd_obj->vt,
                                    zgesdd_obj->ldvt
                                    );

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    zgesdd_obj->inforef = ZGESDD( zgesdd_obj->matrix_layout,
                                    zgesdd_obj->jobz,
                                    zgesdd_obj->m,
                                    zgesdd_obj->n,
                                    zgesdd_obj->aref,
                                    zgesdd_obj->lda,
                                    zgesdd_obj->sref,
                                    zgesdd_obj->uref,
                                    zgesdd_obj->ldu,
                                    zgesdd_obj->vtref,
                                    zgesdd_obj->ldvt
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    zgesdd_obj->diff_a =  computeDiff_z( (zgesdd_obj->m)*(zgesdd_obj->n), 
                zgesdd_obj->a, zgesdd_obj->aref );

    zgesdd_obj->diff_u =  computeDiff_z( (zgesdd_obj->m)*(zgesdd_obj->m), 
                zgesdd_obj->u, zgesdd_obj->uref );

    zgesdd_obj->diff_vt =  computeDiff_z( (zgesdd_obj->n)*(zgesdd_obj->n), 
                zgesdd_obj->vt, zgesdd_obj->vtref );

    zgesdd_obj->diff_s =  computeDiff_d( zgesdd_obj->min_mn, 
                zgesdd_obj->s, zgesdd_obj->sref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n gesdd lapack_complex_double: \n diff_a: %f \n diff_s: %f \n \
diff_u: %f \n diff_vt: %f \n  info: %d \t inforef: %d\n ",
       zgesdd_obj->diff_a, zgesdd_obj->diff_s, zgesdd_obj->diff_u,
       zgesdd_obj->diff_vt, zgesdd_obj->info, zgesdd_obj->inforef);
#endif
}

TEST_F(zgesdd_test, zgesdd_1) {
    EXPECT_NEAR(0.0, zgesdd_obj->diff_a, zgesdd_obj->threshold);
    EXPECT_NEAR(0.0, zgesdd_obj->diff_s, zgesdd_obj->threshold);
    EXPECT_NEAR(0.0, zgesdd_obj->diff_u, zgesdd_obj->threshold);
    EXPECT_NEAR(0.0, zgesdd_obj->diff_vt, zgesdd_obj->threshold);
}

TEST_F(zgesdd_test, zgesdd_2) {
    EXPECT_NEAR(0.0, zgesdd_obj->diff_a, zgesdd_obj->threshold);
    EXPECT_NEAR(0.0, zgesdd_obj->diff_s, zgesdd_obj->threshold);
    EXPECT_NEAR(0.0, zgesdd_obj->diff_u, zgesdd_obj->threshold);
    EXPECT_NEAR(0.0, zgesdd_obj->diff_vt, zgesdd_obj->threshold);
}
TEST_F(zgesdd_test, zgesdd_3) {
    EXPECT_NEAR(0.0, zgesdd_obj->diff_a, zgesdd_obj->threshold);
    EXPECT_NEAR(0.0, zgesdd_obj->diff_s, zgesdd_obj->threshold);
    EXPECT_NEAR(0.0, zgesdd_obj->diff_u, zgesdd_obj->threshold);
    EXPECT_NEAR(0.0, zgesdd_obj->diff_vt, zgesdd_obj->threshold);
}
TEST_F(zgesdd_test, zgesdd_4) {
    EXPECT_NEAR(0.0, zgesdd_obj->diff_a, zgesdd_obj->threshold);
    EXPECT_NEAR(0.0, zgesdd_obj->diff_s, zgesdd_obj->threshold);
    EXPECT_NEAR(0.0, zgesdd_obj->diff_u, zgesdd_obj->threshold);
    EXPECT_NEAR(0.0, zgesdd_obj->diff_vt, zgesdd_obj->threshold);
}