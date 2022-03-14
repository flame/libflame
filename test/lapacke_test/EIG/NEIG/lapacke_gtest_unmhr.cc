#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define LAPACKE_TEST_VERBOSE (1)
#define unmhr_free() \
    if (a!=NULL)        free(a); \
    if (aref!=NULL)     free(aref); \
    if (c!=NULL)        free(c); \
    if (cref!=NULL)     free(cref); \
    if (tau!=NULL)     free(tau); \
    if (tauref!=NULL)  free(tauref)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin unmhr_scomplex_common_parameters  class definition */
class unmhr_scomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_c;
    float threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char trans; // Must be 'N' or 'C'.
    char side	; // Must be 'L' or 'R'.
    lapack_int m,n; // The order of the matrices A and B (n≥ 0).
    lapack_int lda, ldc; //  The leading dimension of a
    lapack_int ilo, ihi;
    lapack_int iloref, ihiref;
    lapack_complex_float* a, *aref; // contains the n-by-n general matrix A.
    lapack_complex_float* tau, *tauref; // contains the elementary reflectors.
	

    /* Output parameters */
    lapack_complex_float* c, *cref; // contains the n-by-n upper triangular matrix B.

    /*Return Values */
    int info, inforef;

   public:
      unmhr_scomplex_parameters (int matrix_layout_i, char side_i, char trans_i,
												lapack_int m_i, lapack_int n_i);
      ~unmhr_scomplex_parameters ();
};

/* Constructor definition  lapack_complex_float_common_parameters */
unmhr_scomplex_parameters:: unmhr_scomplex_parameters (int matrix_layout_i,
				char side_i, char trans_i, lapack_int m_i, lapack_int n_i)
{
	hModule = NULL;
    dModule = NULL;
    diff_c = 0;

    // Initialize 'ilo', 'ihi'.
    ilo = 1;
    iloref = 1;
    ihi = 0;
    ihiref = 0;

    matrix_layout = matrix_layout_i;
    side = side_i;
    trans = trans_i;
    n = n_i;
    m = m_i;
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{
		ldc = m;
	}
	else if (matrix_layout == LAPACK_ROW_MAJOR)
	{
		ldc = n;
	}
	else 
	{
		EXPECT_TRUE(false) << "matrix_layout invalid";
	}
	
    if (side == 'L'){
	   lda = m;
       ihi = m-1;
       ihiref = m-1;
	} else if(side == 'R'){
	   lda = n;
       ihi = n-1;
       ihiref = n-1;
	}

#if LAPACKE_TEST_VERBOSE
   printf(" \n unmhr lapack_complex_float: matrix_layout: %d \t n: %d \t m : %d \t  \
   side: %c  trans: %c \n", matrix_layout, n, m, side, trans );
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, lda*n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &c, &cref, ldc*n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &tau, &tauref, n );

    if( (a==NULL) || (aref==NULL) ||  \
        (c==NULL) || (cref==NULL) || \
        (tau==NULL) || (tauref==NULL) ){
       EXPECT_FALSE( true) << "unmhr_scomplex_parameters object: malloc error. Exiting ";
       unmhr_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, lda*n );
    lapacke_gtest_init_scomplex_buffer_pair_rand( c, cref, ldc*n );
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( tau, tauref, n, 0);
    
   } /* end of Constructor  */
    

/* Destructor definition  'unmhr_scomplex_common_parameters' */
unmhr_scomplex_parameters :: ~unmhr_scomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   //unmhr_free();
    if (a!=NULL)        free(a); 
    if (aref!=NULL)     free(aref); 
    if (c!=NULL)        free(c);
    if (cref!=NULL)     free(cref);
    if (tau!=NULL)     free(tau);
    if (tauref!=NULL)  free(tauref);

} 

//  Test fixture class definition
class cunmhr_test  : public  ::testing::Test {
public:
   unmhr_scomplex_parameters  *cunmhr_obj;
   void SetUp();  
   void TearDown () { delete cunmhr_obj; }
};

void cunmhr_test::SetUp()
{
    /* LAPACKE STRSNA prototype */
    typedef int (*Fptr_NL_LAPACKE_cunmhr) ( int matrix_layout, 
		char side, char trans, lapack_int m, lapack_int n, lapack_int ilo,
		lapack_int ihi, const lapack_complex_float* a, lapack_int lda,
		const lapack_complex_float* tau, lapack_complex_float* c, lapack_int ldc);
                 
    Fptr_NL_LAPACKE_cunmhr Cunmhr;

    cunmhr_obj = new  unmhr_scomplex_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_paramslist[idx].side,
										 eig_non_sym_paramslist[idx].unmhr_trans,
                                         eig_paramslist[idx].n,
                                         eig_paramslist[idx].m );
                                         
    cunmhr_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
	
    cunmhr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cunmhr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cunmhr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cunmhr_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    Cunmhr = (Fptr_NL_LAPACKE_cunmhr)dlsym(cunmhr_obj->hModule, "LAPACKE_cunmhr");
    ASSERT_TRUE(Cunmhr != NULL) << "failed to ppt the Netlib LAPACKE_cunmhr symbol";

    // Compute the inputs to unmhr by invoking gehrd API
    cunmhr_obj->inforef = LAPACKE_cgehrd(   cunmhr_obj->matrix_layout,
                                    cunmhr_obj->n,
                                    cunmhr_obj->ilo,
                                    cunmhr_obj->n,
                                    cunmhr_obj->a,
                                    cunmhr_obj->lda,
                                    cunmhr_obj->tau
                                    );

    cunmhr_obj->inforef = LAPACKE_cgehrd(   cunmhr_obj->matrix_layout,
                                    cunmhr_obj->n,
                                    cunmhr_obj->iloref,
                                    cunmhr_obj->n,
                                    cunmhr_obj->aref,
                                    cunmhr_obj->lda,
                                    cunmhr_obj->tauref
                                    );

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
	lapack_int LAPACKE_cunmhr (int matrix_layout, char side, char trans, lapack_int m, lapack_int n, lapack_int ilo, lapack_int ihi, const lapack_complex_float* a, lapack_int lda, const lapack_complex_float* tau, lapack_complex_float* c, lapack_int ldc);								
    cunmhr_obj->inforef = Cunmhr(   cunmhr_obj->matrix_layout,
                                    cunmhr_obj->side,
                                    cunmhr_obj->trans,
                                    cunmhr_obj->m,
                                    cunmhr_obj->n,
                                    cunmhr_obj->iloref,
                                    cunmhr_obj->ihiref,
                                    cunmhr_obj->aref,
                                    cunmhr_obj->lda,
                                    cunmhr_obj->tauref,
                                    cunmhr_obj->cref,
                                    cunmhr_obj->ldc
                                    );

    /* Compute libflame's Lapacke o/p  */
    cunmhr_obj->info = LAPACKE_cunmhr(  cunmhr_obj->matrix_layout,
                                    cunmhr_obj->side,
                                    cunmhr_obj->trans,
                                    cunmhr_obj->m,
                                    cunmhr_obj->n,
                                    cunmhr_obj->ilo,
                                    cunmhr_obj->ihi,
                                    cunmhr_obj->a,
                                    cunmhr_obj->lda,
                                    cunmhr_obj->tau,
                                    cunmhr_obj->c,
                                    cunmhr_obj->ldc
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    
    cunmhr_obj->diff_c =  computeDiff_c( (cunmhr_obj->ldc)*(cunmhr_obj->n), 
                cunmhr_obj->c, cunmhr_obj->cref );

}

TEST_F(cunmhr_test, cunmhr1) {
    EXPECT_NEAR(0.0, cunmhr_obj->diff_c, cunmhr_obj->threshold);
}

TEST_F(cunmhr_test, cunmhr2) {
    EXPECT_NEAR(0.0, cunmhr_obj->diff_c, cunmhr_obj->threshold);
}

TEST_F(cunmhr_test, cunmhr3) {
    EXPECT_NEAR(0.0, cunmhr_obj->diff_c, cunmhr_obj->threshold);
}

TEST_F(cunmhr_test, cunmhr4) {
    EXPECT_NEAR(0.0, cunmhr_obj->diff_c, cunmhr_obj->threshold);
}


/* Begin unmhr_dcomplex_common_parameters  class definition */
class unmhr_dcomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_c;
    double threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char trans; // Must be 'N' or 'C'.
    char side	; // Must be 'L' or 'R'.
    lapack_int m,n; // The order of the matrices A and B (n≥ 0).
    lapack_int lda, ldc; //  The leading dimension of a
    lapack_int ilo, ihi;
    lapack_int iloref, ihiref;
    lapack_complex_double* a, *aref; // contains the n-by-n general matrix A.
    lapack_complex_double* tau, *tauref; // contains the elementary reflectors.
	

    /* Output parameters */
    lapack_complex_double* c, *cref; // contains the n-by-n upper triangular matrix B.

    /*Return Values */
    int info, inforef;

   public:
      unmhr_dcomplex_parameters (int matrix_layout_i, char side_i, char trans_i,
												lapack_int m_i, lapack_int n_i);
      ~unmhr_dcomplex_parameters ();
};

/* Constructor definition  lapack_complex_double_common_parameters */
unmhr_dcomplex_parameters:: unmhr_dcomplex_parameters (int matrix_layout_i,
				char side_i, char trans_i, lapack_int m_i, lapack_int n_i)
{
	hModule = NULL;
    dModule = NULL;
    diff_c = 0;

    // Initialize 'ilo', 'ihi'.
    ilo = 1;
    iloref = 1;
    ihi = 0;
    ihiref = 0;

    matrix_layout = matrix_layout_i;
    side = side_i;
    trans = trans_i;
    n = n_i;
    m = m_i;
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{
		ldc = m;
	}
	else if (matrix_layout == LAPACK_ROW_MAJOR)
	{
		ldc = n;
	}
	else 
	{
		EXPECT_TRUE(false) << "matrix_layout invalid";
	}
	
    if (side == 'L'){
	   lda = m;
       ihi = m-1;
       ihiref = m-1;
	} else if(side == 'R'){
	   lda = n;
       ihi = n-1;
       ihiref = n-1;
	}

#if LAPACKE_TEST_VERBOSE
   printf(" \n unmhr lapack_complex_double: matrix_layout: %d \t n: %d \t m : %d \t  \
   side: %c  trans: %c \n", matrix_layout, n, m, side, trans );
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, lda*n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &c, &cref, ldc*n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &tau, &tauref, n );

    if( (a==NULL) || (aref==NULL) ||  \
        (c==NULL) || (cref==NULL) || \
        (tau==NULL) || (tauref==NULL) ){
       EXPECT_FALSE( true) << "unmhr_dcomplex_parameters object: malloc error. Exiting ";
       unmhr_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, lda*n );
    lapacke_gtest_init_dcomplex_buffer_pair_rand( c, cref, ldc*n );
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( tau, tauref, n, 0);
    
   } /* end of Constructor  */
    

/* Destructor definition  'unmhr_dcomplex_common_parameters' */
unmhr_dcomplex_parameters :: ~unmhr_dcomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   //unmhr_free();
    if (a!=NULL)        free(a); 
    if (aref!=NULL)     free(aref); 
    if (c!=NULL)        free(c);
    if (cref!=NULL)     free(cref);
    if (tau!=NULL)     free(tau);
    if (tauref!=NULL)  free(tauref);

} 

//  Test fixture class definition
class zunmhr_test  : public  ::testing::Test {
public:
   unmhr_dcomplex_parameters  *zunmhr_obj;
   void SetUp();  
   void TearDown () { delete zunmhr_obj; }
};

void zunmhr_test::SetUp()
{
    /* LAPACKE STRSNA prototype */
    typedef int (*Fptr_NL_LAPACKE_zunmhr) ( int matrix_layout, 
		char side, char trans, lapack_int m, lapack_int n, lapack_int ilo,
		lapack_int ihi, const lapack_complex_double* a, lapack_int lda,
		const lapack_complex_double* tau, lapack_complex_double* c, lapack_int ldc);
                 
    Fptr_NL_LAPACKE_zunmhr zunmhr;

    zunmhr_obj = new  unmhr_dcomplex_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_paramslist[idx].side,
										 eig_non_sym_paramslist[idx].unmhr_trans,
                                         eig_paramslist[idx].n,
                                         eig_paramslist[idx].m );
                                         
    zunmhr_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
	
    zunmhr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zunmhr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zunmhr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zunmhr_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    zunmhr = (Fptr_NL_LAPACKE_zunmhr)dlsym(zunmhr_obj->hModule, "LAPACKE_zunmhr");
    ASSERT_TRUE(zunmhr != NULL) << "failed to ppt the Netlib LAPACKE_zunmhr symbol";

    // Compute the inputs to unmhr by invoking gehrd API
    zunmhr_obj->inforef = LAPACKE_zgehrd(   zunmhr_obj->matrix_layout,
                                    zunmhr_obj->n,
                                    zunmhr_obj->ilo,
                                    zunmhr_obj->n,
                                    zunmhr_obj->a,
                                    zunmhr_obj->lda,
                                    zunmhr_obj->tau
                                    );

    zunmhr_obj->inforef = LAPACKE_zgehrd(   zunmhr_obj->matrix_layout,
                                    zunmhr_obj->n,
                                    zunmhr_obj->iloref,
                                    zunmhr_obj->n,
                                    zunmhr_obj->aref,
                                    zunmhr_obj->lda,
                                    zunmhr_obj->tauref
                                    );

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
	lapack_int LAPACKE_zunmhr (int matrix_layout, char side, char trans, lapack_int m, lapack_int n, lapack_int ilo, lapack_int ihi, const lapack_complex_double* a, lapack_int lda, const lapack_complex_double* tau, lapack_complex_double* c, lapack_int ldc);								
    zunmhr_obj->inforef = zunmhr(   zunmhr_obj->matrix_layout,
                                    zunmhr_obj->side,
                                    zunmhr_obj->trans,
                                    zunmhr_obj->m,
                                    zunmhr_obj->n,
                                    zunmhr_obj->iloref,
                                    zunmhr_obj->ihiref,
                                    zunmhr_obj->aref,
                                    zunmhr_obj->lda,
                                    zunmhr_obj->tauref,
                                    zunmhr_obj->cref,
                                    zunmhr_obj->ldc
                                    );

    /* Compute libflame's Lapacke o/p  */
    zunmhr_obj->info = LAPACKE_zunmhr(  zunmhr_obj->matrix_layout,
                                    zunmhr_obj->side,
                                    zunmhr_obj->trans,
                                    zunmhr_obj->m,
                                    zunmhr_obj->n,
                                    zunmhr_obj->ilo,
                                    zunmhr_obj->ihi,
                                    zunmhr_obj->a,
                                    zunmhr_obj->lda,
                                    zunmhr_obj->tau,
                                    zunmhr_obj->c,
                                    zunmhr_obj->ldc
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    
    zunmhr_obj->diff_c =  computeDiff_z( (zunmhr_obj->ldc)*(zunmhr_obj->n), 
                zunmhr_obj->c, zunmhr_obj->cref );

}

TEST_F(zunmhr_test, zunmhr1) {
    EXPECT_NEAR(0.0, zunmhr_obj->diff_c, zunmhr_obj->threshold);
}

TEST_F(zunmhr_test, zunmhr2) {
    EXPECT_NEAR(0.0, zunmhr_obj->diff_c, zunmhr_obj->threshold);
}

TEST_F(zunmhr_test, zunmhr3) {
    EXPECT_NEAR(0.0, zunmhr_obj->diff_c, zunmhr_obj->threshold);
}

TEST_F(zunmhr_test, zunmhr4) {
    EXPECT_NEAR(0.0, zunmhr_obj->diff_c, zunmhr_obj->threshold);
}
