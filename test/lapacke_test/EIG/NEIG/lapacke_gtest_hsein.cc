#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"

#define hsein_free() \
    if (h!=NULL)        free(h); \
    if (href!=NULL)     free(href); \
    if (select!=NULL)   free(select); \
    if (selectref!=NULL) free(selectref); \
    if (wr!=NULL)       free(wr); \
    if (wrref!=NULL)    free(wrref); \
    if (wi!=NULL)       free(wi); \
    if (wiref!=NULL)    free(wiref);\
    if (vl!=NULL)       free(vl); \
    if (vlref!=NULL)    free(vlref); \
    if (vr!=NULL)       free(vr); \
    if (vrref!=NULL)    free(vrref); \
    if (ifaill!=NULL)    free(ifaill); \
    if (ifaillref!=NULL) free(ifaillref); \
    if (ifailr!=NULL)    free(ifailr); \
    if (ifailrref!=NULL) free(ifailrref)

#define hsein_cplx_free() \
    if (h!=NULL)        free(h); \
    if (href!=NULL)     free(href); \
    if (select!=NULL)   free(select); \
    if (selectref!=NULL) free(selectref); \
    if (w!=NULL)       free(w); \
    if (wref!=NULL)    free(wref); \
    if (vl!=NULL)       free(vl); \
    if (vlref!=NULL)    free(vlref); \
    if (vr!=NULL)       free(vr); \
    if (vrref!=NULL)    free(vrref); \
    if (ifaill!=NULL)    free(ifaill); \
    if (ifaillref!=NULL) free(ifaillref); \
    if (ifailr!=NULL)    free(ifailr); \
    if (ifailrref!=NULL) free(ifailrref)


#define LAPACKE_TEST_VERBOSE  (1)
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin hsein_float_common_parameters  class definition */
class hsein_float_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_wr, diff_wi;;
    void *hModule, *dModule;
    
    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char side;  // Must be 'R' or 'L' or 'B'.
    char eigsrc; // Must be 'Q' or  'N'.
	char initv; // Must be 'N' or 'U'.
	lapack_int ldh;
	float *vl, *vlref; // starting vectors for the inverse iteration for the left eigenvectors.
	lapack_int ldvl;
	float *vr, *vrref; // starting vectors for the inverse iteration for the right eigenvectors
	lapack_int ldvr;
	lapack_int mm; // The number of columns in the arrays VL and/or VR
    lapack_int n; // order of matrix

    /* Input/ Output parameters */
	char *select, *selectref; // Specifies the eigenvectors to be computed.  
	float *h, *href; // The hessenberg matrix / Schur decomposition o/p
	float *wr, *wrref; // the real part of eigenvalues of the matrix H.
    
    /* Output parameters */
	/* Primary o/p is wi, wr - holds the o/p eigen values of i/p matrix 'h' */
    lapack_int m, mref; // The number of columns in the arrays VL and/or VR required to store the eigenvectors
	float *wi, *wiref; // the imaginary part of eigenvalues of the matrix H.
	lapack_int *ifaill, *ifaillref;
	lapack_int *ifailr, *ifailrref;
	
    /*Return Values */
    int info, inforef;

   public:
      hsein_float_parameters (int matrix_layout_i, char side_i, char eigsrc_i,
                              char initv_i, lapack_int n_i, lapack_int mm_i );

      ~hsein_float_parameters ();
};

/* Constructor definition  float_common_parameters */
hsein_float_parameters:: hsein_float_parameters (int matrix_layout_i,
							char side_i, char eigsrc_i, char initv_i,
									lapack_int n_i, lapack_int mm_i)
{
    matrix_layout = matrix_layout_i;
    side = side_i;
    eigsrc = eigsrc_i;
	initv = initv_i;
    n  = n_i;
	mm = n;;
   
    ldh = n;	
	if ( matrix_layout== LAPACK_ROW_MAJOR)
	{
		ldvl = mm;
		ldvr = mm;
	}
	else
	{
		ldvl = n;
		ldvr = n;
	}
	
    hModule = NULL;
    dModule = NULL;
	m = n;
	
	mref = 0;
    diff_wr = 0;
    diff_wi = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n hsein float: matrix_layout: %d n: %d  mm: %d side: %c \
eigsrc: %c  initv: %c \n", matrix_layout, n, mm, side, eigsrc, initv);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &h, &href, n*ldh );
    lapacke_gtest_alloc_float_buffer_pair( &wr, &wrref, n );
    lapacke_gtest_alloc_float_buffer_pair( &wi, &wiref, n );
    lapacke_gtest_alloc_float_buffer_pair( &vl, &vlref, n*mm );
    lapacke_gtest_alloc_float_buffer_pair( &vr, &vrref, n*mm );
    lapacke_gtest_alloc_int_buffer_pair( &ifaill, &ifaillref, mm );
    lapacke_gtest_alloc_int_buffer_pair( &ifailr, &ifailrref, mm );
    lapacke_gtest_alloc_char_buffer_pair( &select, &selectref, n );

    if( (h==NULL) || (href==NULL) ||  \
        (vl==NULL) || (vlref==NULL) || \
        (vr==NULL) || (vrref==NULL) || \
        (ifaill==NULL) || (ifaillref==NULL) || \
        (ifailr==NULL) || (ifailrref==NULL) || \
        (select==NULL) || (selectref==NULL) || \
        (wr==NULL) || (wrref==NULL) || \
        (wi==NULL) || (wiref==NULL)  ){
       EXPECT_FALSE( true) << "hsein_float_parameters object: malloc error. Exiting ";
       hsein_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
	/* Works for both general matrix, upper traingualr i/p matrices */
    lapacke_gtest_init_float_buffer_pair_rand_custom_matrix( h, href, n,ldh, 'U');
    //lapacke_gtest_init_float_buffer_pair_rand( h, href, n*ldh);
    lapacke_gtest_init_float_buffer_pair_rand(wr, wrref, n);
    lapacke_gtest_init_float_buffer_pair_with_constant(wi, wiref, n, 0.0);
    lapacke_gtest_init_float_buffer_pair_rand(vl, vlref, n*mm);
    lapacke_gtest_init_float_buffer_pair_rand(vr, vrref, n*mm);
    lapacke_gtest_init_int_buffer_pair_with_constant(ifailr, ifailrref, mm, 0);
    lapacke_gtest_init_int_buffer_pair_with_constant(ifaill, ifaillref, mm, 0);
    lapacke_gtest_init_char_buffer_pair_with_constant(select, selectref, n, 0xff);
    
   } /* end of Constructor  */
    

/* Destructor definition  'hsein_float_common_parameters' */
hsein_float_parameters :: ~hsein_float_parameters ()
{
   /* De-Allocate memory for the input matrices */
   hsein_free();
} 

//  Test fixture class definition
class shsein_test  : public  ::testing::Test {
public:
   hsein_float_parameters  *shsein_obj;
   void SetUp();  
   void TearDown () { delete shsein_obj; }
};

void shsein_test::SetUp()
{
    /* LAPACKE SHSEIN prototype */
    typedef int (*Fptr_NL_LAPACKE_shsein) ( int matrix_layout, char side,
			char eigsrc, char initv, lapack_logical* select, lapack_int n,
			const float* h, lapack_int ldh, float* wr, const float* wi,
			float* vl, lapack_int ldvl, float* vr, lapack_int ldvr,
			lapack_int mm, lapack_int* m, lapack_int* ifaill,
			lapack_int* ifailr );
            
    Fptr_NL_LAPACKE_shsein SHSEIN;      

    shsein_obj = new  hsein_float_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_paramslist[idx].side,
                                         eig_paramslist[idx].eigsrc,
                                         eig_paramslist[idx].initv,
                                         eig_paramslist[idx].n,
                                         eig_paramslist[idx].kb
										 );
    idx = Circular_Increment_Index(idx);

    shsein_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    shsein_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(shsein_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(shsein_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SHSEIN = (Fptr_NL_LAPACKE_shsein)dlsym(shsein_obj->hModule, "LAPACKE_shsein");
    ASSERT_TRUE(SHSEIN != NULL) << "failed to ppt the Netlib LAPACKE_shsein symbol";

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    shsein_obj->inforef = SHSEIN(   shsein_obj->matrix_layout,
                                    shsein_obj->side,
                                    shsein_obj->eigsrc, 
                                    shsein_obj->initv,
                  (lapack_logical *)shsein_obj->selectref,
                                    shsein_obj->n,
                                    shsein_obj->href, 
                                    shsein_obj->ldh,
                                    shsein_obj->wrref,
                                    shsein_obj->wiref,
                                    shsein_obj->vlref,
                                    shsein_obj->ldvl,
                                    shsein_obj->vrref,
                                    shsein_obj->ldvr,
                                    shsein_obj->mm,
                                    &shsein_obj->mref,
                                    shsein_obj->ifaillref,
                                    shsein_obj->ifailrref
									);

    /* Compute libflame's Lapacke o/p  */
    shsein_obj->info = LAPACKE_shsein(  shsein_obj->matrix_layout,
                                    shsein_obj->side,
                                    shsein_obj->eigsrc, 
                                    shsein_obj->initv,
                  (lapack_logical *)shsein_obj->select,
                                    shsein_obj->n,
                                    shsein_obj->h, 
                                    shsein_obj->ldh,
                                    shsein_obj->wr,
                                    shsein_obj->wi,
                                    shsein_obj->vl,
                                    shsein_obj->ldvl,
                                    shsein_obj->vr,
                                    shsein_obj->ldvr,
                                    shsein_obj->mm,
                                    &shsein_obj->m,
                                    shsein_obj->ifaill,
                                    shsein_obj->ifailr
									);

    /* Capture the Netlib, libflame o/p buffers' differences */
    shsein_obj->diff_wr =  computeDiff_s( shsein_obj->n, 
                shsein_obj->wr, shsein_obj->wrref );

    shsein_obj->diff_wi =  computeDiff_s( shsein_obj->n, 
                shsein_obj->wi, shsein_obj->wiref );

}

TEST_F(shsein_test, shsein1) {
    EXPECT_NEAR(0.0, shsein_obj->diff_wr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, shsein_obj->diff_wi, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(shsein_test, shsein2) {
    EXPECT_NEAR(0.0, shsein_obj->diff_wr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, shsein_obj->diff_wi, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(shsein_test, shsein3) {
    EXPECT_NEAR(0.0, shsein_obj->diff_wr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, shsein_obj->diff_wi, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(shsein_test, shsein4) {
    EXPECT_NEAR(0.0, shsein_obj->diff_wr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, shsein_obj->diff_wi, LAPACKE_GTEST_THRESHOLD);
}

/* Begin hsein_double_common_parameters  class definition */
class hsein_double_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    double diff_wr, diff_wi;;
    void *hModule, *dModule;
    
    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char side;  // Must be 'R' or 'L' or 'B'.
    char eigsrc; // Must be 'Q' or  'N'.
	char initv; // Must be 'N' or 'U'.
	lapack_int ldh;
	double *vl, *vlref; // starting vectors for the inverse iteration for the left eigenvectors.
	lapack_int ldvl;
	double *vr, *vrref; // starting vectors for the inverse iteration for the right eigenvectors
	lapack_int ldvr;
	lapack_int mm; // The number of columns in the arrays VL and/or VR
    lapack_int n; // order of matrix

    /* Input/ Output parameters */
	char *select, *selectref; // Specifies the eigenvectors to be computed.  
	double *h, *href; // The hessenberg matrix / Schur decomposition o/p
	double *wr, *wrref; // the real part of eigenvalues of the matrix H.
    
    /* Output parameters */
	/* Primary o/p is wi, wr - holds the o/p eigen values of i/p matrix 'h' */
    lapack_int m, mref; // The number of columns in the arrays VL and/or VR required to store the eigenvectors
	double *wi, *wiref; // the imaginary part of eigenvalues of the matrix H.
	lapack_int *ifaill, *ifaillref;
	lapack_int *ifailr, *ifailrref;
	
    /*Return Values */
    int info, inforef;

   public:
      hsein_double_parameters (int matrix_layout_i, char side_i, char eigsrc_i,
                              char initv_i, lapack_int n_i, lapack_int mm_i );

      ~hsein_double_parameters ();
};

/* Constructor definition  double_common_parameters */
hsein_double_parameters:: hsein_double_parameters (int matrix_layout_i,
							char side_i, char eigsrc_i, char initv_i,
									lapack_int n_i, lapack_int mm_i)
{
    matrix_layout = matrix_layout_i;
    side = side_i;
    eigsrc = eigsrc_i;
	initv = initv_i;
    n  = n_i;
	mm = n;;
   
    ldh = n;	
	if ( matrix_layout== LAPACK_ROW_MAJOR)
	{
		ldvl = mm;
		ldvr = mm;
	}
	else
	{
		ldvl = n;
		ldvr = n;
	}
	
    hModule = NULL;
    dModule = NULL;
	m = n;
	
	mref = 0;
    diff_wr = 0;
    diff_wi = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n hsein double: matrix_layout: %d n: %d  mm: %d side: %c \
eigsrc: %c  initv: %c \n", matrix_layout, n, mm, side, eigsrc, initv);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &h, &href, n*ldh );
    lapacke_gtest_alloc_double_buffer_pair( &wr, &wrref, n );
    lapacke_gtest_alloc_double_buffer_pair( &wi, &wiref, n );
    lapacke_gtest_alloc_double_buffer_pair( &vl, &vlref, n*mm );
    lapacke_gtest_alloc_double_buffer_pair( &vr, &vrref, n*mm );
    lapacke_gtest_alloc_int_buffer_pair( &ifaill, &ifaillref, mm );
    lapacke_gtest_alloc_int_buffer_pair( &ifailr, &ifailrref, mm );
    lapacke_gtest_alloc_char_buffer_pair( &select, &selectref, n );

    if( (h==NULL) || (href==NULL) ||  \
        (vl==NULL) || (vlref==NULL) || \
        (vr==NULL) || (vrref==NULL) || \
        (ifaill==NULL) || (ifaillref==NULL) || \
        (ifailr==NULL) || (ifailrref==NULL) || \
        (select==NULL) || (selectref==NULL) || \
        (wr==NULL) || (wrref==NULL) || \
        (wi==NULL) || (wiref==NULL)  ){
       EXPECT_FALSE( true) << "hsein_double_parameters object: malloc error. Exiting ";
       hsein_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
	/* Works for both general matrix, upper traingualr i/p matrices */
    lapacke_gtest_init_double_buffer_pair_rand_custom_matrix( h, href, n,ldh, 'U');
    //lapacke_gtest_init_double_buffer_pair_rand( h, href, n*ldh);
    lapacke_gtest_init_double_buffer_pair_rand(wr, wrref, n);
    lapacke_gtest_init_double_buffer_pair_with_constant(wi, wiref, n, 0.0);
    lapacke_gtest_init_double_buffer_pair_rand(vl, vlref, n*mm);
    lapacke_gtest_init_double_buffer_pair_rand(vr, vrref, n*mm);
    lapacke_gtest_init_int_buffer_pair_with_constant(ifailr, ifailrref, mm, 0);
    lapacke_gtest_init_int_buffer_pair_with_constant(ifaill, ifaillref, mm, 0);
    lapacke_gtest_init_char_buffer_pair_with_constant(select, selectref, n, 0xff);
    
   } /* end of Constructor  */
    

/* Destructor definition  'hsein_double_common_parameters' */
hsein_double_parameters :: ~hsein_double_parameters ()
{
   /* De-Allocate memory for the input matrices */
   hsein_free();
} 

//  Test fixture class definition
class dhsein_test  : public  ::testing::Test {
public:
   hsein_double_parameters  *dhsein_obj;
   void SetUp();  
   void TearDown () { delete dhsein_obj; }
};

void dhsein_test::SetUp()
{
    /* LAPACKE DHSEIN prototype */
    typedef int (*Fptr_NL_LAPACKE_dhsein) ( int matrix_layout, char side,
			char eigsrc, char initv, lapack_logical* select, lapack_int n,
			const double* h, lapack_int ldh, double* wr, const double* wi,
			double* vl, lapack_int ldvl, double* vr, lapack_int ldvr,
			lapack_int mm, lapack_int* m, lapack_int* ifaill,
			lapack_int* ifailr );
            
    Fptr_NL_LAPACKE_dhsein DHSEIN;      

    dhsein_obj = new  hsein_double_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_paramslist[idx].side,
                                         eig_paramslist[idx].eigsrc,
                                         eig_paramslist[idx].initv,
                                         eig_paramslist[idx].n,
                                         eig_paramslist[idx].kb
										 );
    idx = Circular_Increment_Index(idx);

    dhsein_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dhsein_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dhsein_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dhsein_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DHSEIN = (Fptr_NL_LAPACKE_dhsein)dlsym(dhsein_obj->hModule, "LAPACKE_dhsein");
    ASSERT_TRUE(DHSEIN != NULL) << "failed to ppt the Netlib LAPACKE_dhsein symbol";

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    dhsein_obj->inforef = DHSEIN(   dhsein_obj->matrix_layout,
                                    dhsein_obj->side,
                                    dhsein_obj->eigsrc, 
                                    dhsein_obj->initv,
                  (lapack_logical *)dhsein_obj->selectref,
                                    dhsein_obj->n,
                                    dhsein_obj->href, 
                                    dhsein_obj->ldh,
                                    dhsein_obj->wrref,
                                    dhsein_obj->wiref,
                                    dhsein_obj->vlref,
                                    dhsein_obj->ldvl,
                                    dhsein_obj->vrref,
                                    dhsein_obj->ldvr,
                                    dhsein_obj->mm,
                                    &dhsein_obj->mref,
                                    dhsein_obj->ifaillref,
                                    dhsein_obj->ifailrref
									);

    /* Compute libflame's Lapacke o/p  */
    dhsein_obj->info = LAPACKE_dhsein(  dhsein_obj->matrix_layout,
                                    dhsein_obj->side,
                                    dhsein_obj->eigsrc, 
                                    dhsein_obj->initv,
                  (lapack_logical *)dhsein_obj->select,
                                    dhsein_obj->n,
                                    dhsein_obj->h, 
                                    dhsein_obj->ldh,
                                    dhsein_obj->wr,
                                    dhsein_obj->wi,
                                    dhsein_obj->vl,
                                    dhsein_obj->ldvl,
                                    dhsein_obj->vr,
                                    dhsein_obj->ldvr,
                                    dhsein_obj->mm,
                                    &dhsein_obj->m,
                                    dhsein_obj->ifaill,
                                    dhsein_obj->ifailr
									);

    /* Capture the Netlib, libflame o/p buffers' differences */
    dhsein_obj->diff_wr =  computeDiff_d( dhsein_obj->n, 
                dhsein_obj->wr, dhsein_obj->wrref );

    dhsein_obj->diff_wi =  computeDiff_d( dhsein_obj->n, 
                dhsein_obj->wi, dhsein_obj->wiref );

}

TEST_F(dhsein_test, dhsein1) {
    EXPECT_NEAR(0.0, dhsein_obj->diff_wr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dhsein_obj->diff_wi, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dhsein_test, dhsein2) {
    EXPECT_NEAR(0.0, dhsein_obj->diff_wr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dhsein_obj->diff_wi, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dhsein_test, dhsein3) {
    EXPECT_NEAR(0.0, dhsein_obj->diff_wr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dhsein_obj->diff_wi, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dhsein_test, dhsein4) {
    EXPECT_NEAR(0.0, dhsein_obj->diff_wr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dhsein_obj->diff_wi, LAPACKE_GTEST_THRESHOLD);
}

/* Begin hsein_scomplex_common_parameters  class definition */
class hsein_scomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_w;
    void *hModule, *dModule;
    
    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char side;  // Must be 'R' or 'L' or 'B'.
    char eigsrc; // Must be 'Q' or  'N'.
	char initv; // Must be 'N' or 'U'.
	lapack_int ldh;
	lapack_complex_float *vl, *vlref; // starting vectors for the inverse iteration for the left eigenvectors.
	lapack_int ldvl;
	lapack_complex_float *vr, *vrref; // starting vectors for the inverse iteration for the right eigenvectors
	lapack_int ldvr;
	lapack_int mm; // The number of columns in the arrays VL and/or VR
    lapack_int n; // order of matrix

    /* Input/ Output parameters */
	char *select, *selectref; // Specifies the eigenvectors to be computed.  
	lapack_complex_float *h, *href; // The hessenberg matrix / Schur decomposition o/p
	lapack_complex_float *w, *wref; // the real part of eigenvalues of the matrix H.
    
    /* Output parameters */
	/* Primary o/p is 'w' - holds the o/p eigen values of i/p matrix 'h' */
    lapack_int m, mref; // The number of columns in the arrays VL and/or VR required to store the eigenvectors
	lapack_int *ifaill, *ifaillref;
	lapack_int *ifailr, *ifailrref;
	
    /*Return Values */
    int info, inforef;

   public:
      hsein_scomplex_parameters (int matrix_layout_i, char side_i, char eigsrc_i,
                              char initv_i, lapack_int n_i, lapack_int mm_i );

      ~hsein_scomplex_parameters ();
};

/* Constructor definition  lapack_complex_float_common_parameters */
hsein_scomplex_parameters:: hsein_scomplex_parameters (int matrix_layout_i,
							char side_i, char eigsrc_i, char initv_i,
									lapack_int n_i, lapack_int mm_i)
{
    matrix_layout = matrix_layout_i;
    side = side_i;
    eigsrc = eigsrc_i;
	initv = initv_i;
    n  = n_i;
	mm = n;;
   
    ldh = n;	
	if ( matrix_layout== LAPACK_ROW_MAJOR)
	{
		ldvl = mm;
		ldvr = mm;
	}
	else
	{
		ldvl = n;
		ldvr = n;
	}
	
    hModule = NULL;
    dModule = NULL;
	m = n;
	
	mref = 0;
    diff_w = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n hsein lapack_complex_float: matrix_layout: %d n: %d  mm: %d side: %c \
eigsrc: %c  initv: %c \n", matrix_layout, n, mm, side, eigsrc, initv);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &h, &href, n*ldh );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &w, &wref, n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &vl, &vlref, n*mm );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &vr, &vrref, n*mm );
    lapacke_gtest_alloc_int_buffer_pair( &ifaill, &ifaillref, mm );
    lapacke_gtest_alloc_int_buffer_pair( &ifailr, &ifailrref, mm );
    lapacke_gtest_alloc_char_buffer_pair( &select, &selectref, n );

    if( (h==NULL) || (href==NULL) ||  \
        (vl==NULL) || (vlref==NULL) || \
        (vr==NULL) || (vrref==NULL) || \
        (ifaill==NULL) || (ifaillref==NULL) || \
        (ifailr==NULL) || (ifailrref==NULL) || \
        (select==NULL) || (selectref==NULL) || \
        (w==NULL) || (wref==NULL)  ){
       EXPECT_FALSE( true) << "hsein_scomplex_parameters object: malloc error. Exiting ";
       hsein_cplx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
	/* Works for both general matrix, upper traingualr i/p matrices */
    lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix( h, href, n,ldh, 'U');
    //lapacke_gtest_init_scomplex_buffer_pair_rand( h, href, n*ldh);
    lapacke_gtest_init_scomplex_buffer_pair_rand(w, wref, n);
    lapacke_gtest_init_scomplex_buffer_pair_rand(vl, vlref, n*mm);
    lapacke_gtest_init_scomplex_buffer_pair_rand(vr, vrref, n*mm);
    lapacke_gtest_init_int_buffer_pair_with_constant(ifailr, ifailrref, mm, 0);
    lapacke_gtest_init_int_buffer_pair_with_constant(ifaill, ifaillref, mm, 0);
    lapacke_gtest_init_char_buffer_pair_with_constant(select, selectref, n, 0xff);
    
   } /* end of Constructor  */
    

/* Destructor definition  'hsein_scomplex_common_parameters' */
hsein_scomplex_parameters :: ~hsein_scomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   hsein_cplx_free();
} 

//  Test fixture class definition
class chsein_test  : public  ::testing::Test {
public:
   hsein_scomplex_parameters  *chsein_obj;
   void SetUp();  
   void TearDown () { delete chsein_obj; }
};

void chsein_test::SetUp()
{
    /* LAPACKE CHSEIN prototype */
    typedef int (*Fptr_NL_LAPACKE_chsein) ( int matrix_layout, char side,
			char eigsrc, char initv, lapack_logical* select, lapack_int n,
			const lapack_complex_float* h, lapack_int ldh, 
			lapack_complex_float* w, lapack_complex_float* vl,
			lapack_int ldvl, lapack_complex_float* vr, lapack_int ldvr,
			lapack_int mm, lapack_int* m, lapack_int* ifaill,
			lapack_int* ifailr );
            
    Fptr_NL_LAPACKE_chsein CHSEIN;      

    chsein_obj = new  hsein_scomplex_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_paramslist[idx].side,
                                         eig_paramslist[idx].eigsrc,
                                         eig_paramslist[idx].initv,
                                         eig_paramslist[idx].n,
                                         eig_paramslist[idx].kb
										 );
    idx = Circular_Increment_Index(idx);

    chsein_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    chsein_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(chsein_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(chsein_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CHSEIN = (Fptr_NL_LAPACKE_chsein)dlsym(chsein_obj->hModule, "LAPACKE_chsein");
    ASSERT_TRUE(CHSEIN != NULL) << "failed to ppt the Netlib LAPACKE_chsein symbol";

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    chsein_obj->inforef = CHSEIN(   chsein_obj->matrix_layout,
                                    chsein_obj->side,
                                    chsein_obj->eigsrc, 
                                    chsein_obj->initv,
                  (lapack_logical *)chsein_obj->selectref,
                                    chsein_obj->n,
                                    chsein_obj->href, 
                                    chsein_obj->ldh,
                                    chsein_obj->wref,
                                    chsein_obj->vlref,
                                    chsein_obj->ldvl,
                                    chsein_obj->vrref,
                                    chsein_obj->ldvr,
                                    chsein_obj->mm,
                                    &chsein_obj->mref,
                                    chsein_obj->ifaillref,
                                    chsein_obj->ifailrref
									);

    /* Compute libflame's Lapacke o/p  */
    chsein_obj->info = LAPACKE_chsein(  chsein_obj->matrix_layout,
                                    chsein_obj->side,
                                    chsein_obj->eigsrc, 
                                    chsein_obj->initv,
                  (lapack_logical *)chsein_obj->select,
                                    chsein_obj->n,
                                    chsein_obj->h, 
                                    chsein_obj->ldh,
                                    chsein_obj->w,
                                    chsein_obj->vl,
                                    chsein_obj->ldvl,
                                    chsein_obj->vr,
                                    chsein_obj->ldvr,
                                    chsein_obj->mm,
                                    &chsein_obj->m,
                                    chsein_obj->ifaill,
                                    chsein_obj->ifailr
									);

    /* Capture the Netlib, libflame o/p buffers' differences */
    chsein_obj->diff_w =  computeDiff_c( chsein_obj->n, 
                chsein_obj->w, chsein_obj->wref );

}

TEST_F(chsein_test, chsein1) {
    EXPECT_NEAR(0.0, chsein_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chsein_test, chsein2) {
    EXPECT_NEAR(0.0, chsein_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chsein_test, chsein3) {
    EXPECT_NEAR(0.0, chsein_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chsein_test, chsein4) {
    EXPECT_NEAR(0.0, chsein_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}


/* Begin hsein_dcomplex_common_parameters  class definition */
class hsein_dcomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    double diff_w;
    void *hModule, *dModule;
    
    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char side;  // Must be 'R' or 'L' or 'B'.
    char eigsrc; // Must be 'Q' or  'N'.
	char initv; // Must be 'N' or 'U'.
	lapack_int ldh;
	lapack_complex_double *vl, *vlref; // starting vectors for the inverse iteration for the left eigenvectors.
	lapack_int ldvl;
	lapack_complex_double *vr, *vrref; // starting vectors for the inverse iteration for the right eigenvectors
	lapack_int ldvr;
	lapack_int mm; // The number of columns in the arrays VL and/or VR
    lapack_int n; // order of matrix

    /* Input/ Output parameters */
	char *select, *selectref; // Specifies the eigenvectors to be computed.  
	lapack_complex_double *h, *href; // The hessenberg matrix / Schur decomposition o/p
	lapack_complex_double *w, *wref; // the real part of eigenvalues of the matrix H.
    
    /* Output parameters */
	/* Primary o/p is 'w' - holds the o/p eigen values of i/p matrix 'h' */
    lapack_int m, mref; // The number of columns in the arrays VL and/or VR required to store the eigenvectors
	lapack_int *ifaill, *ifaillref;
	lapack_int *ifailr, *ifailrref;
	
    /*Return Values */
    int info, inforef;

   public:
      hsein_dcomplex_parameters (int matrix_layout_i, char side_i, char eigsrc_i,
                              char initv_i, lapack_int n_i, lapack_int mm_i );

      ~hsein_dcomplex_parameters ();
};

/* Constructor definition  lapack_complex_double_common_parameters */
hsein_dcomplex_parameters:: hsein_dcomplex_parameters (int matrix_layout_i,
							char side_i, char eigsrc_i, char initv_i,
									lapack_int n_i, lapack_int mm_i)
{
    matrix_layout = matrix_layout_i;
    side = side_i;
    eigsrc = eigsrc_i;
	initv = initv_i;
    n  = n_i;
	mm = n;;
   
    ldh = n;	
	if ( matrix_layout== LAPACK_ROW_MAJOR)
	{
		ldvl = mm;
		ldvr = mm;
	}
	else
	{
		ldvl = n;
		ldvr = n;
	}
	
    hModule = NULL;
    dModule = NULL;
	m = n;
	
	mref = 0;
    diff_w = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n hsein lapack_complex_double: matrix_layout: %d n: %d  mm: %d side: %c \
eigsrc: %c  initv: %c \n", matrix_layout, n, mm, side, eigsrc, initv);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &h, &href, n*ldh );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &w, &wref, n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &vl, &vlref, n*mm );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &vr, &vrref, n*mm );
    lapacke_gtest_alloc_int_buffer_pair( &ifaill, &ifaillref, mm );
    lapacke_gtest_alloc_int_buffer_pair( &ifailr, &ifailrref, mm );
    lapacke_gtest_alloc_char_buffer_pair( &select, &selectref, n );

    if( (h==NULL) || (href==NULL) ||  \
        (vl==NULL) || (vlref==NULL) || \
        (vr==NULL) || (vrref==NULL) || \
        (ifaill==NULL) || (ifaillref==NULL) || \
        (ifailr==NULL) || (ifailrref==NULL) || \
        (select==NULL) || (selectref==NULL) || \
        (w==NULL) || (wref==NULL)  ){
       EXPECT_FALSE( true) << "hsein_dcomplex_parameters object: malloc error. Exiting ";
       hsein_cplx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
	/* Works for both general matrix, upper traingualr i/p matrices */
    lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix( h, href, n,ldh, 'U');
    //lapacke_gtest_init_dcomplex_buffer_pair_rand( h, href, n*ldh);
    lapacke_gtest_init_dcomplex_buffer_pair_rand(w, wref, n);
    lapacke_gtest_init_dcomplex_buffer_pair_rand(vl, vlref, n*mm);
    lapacke_gtest_init_dcomplex_buffer_pair_rand(vr, vrref, n*mm);
    lapacke_gtest_init_int_buffer_pair_with_constant(ifailr, ifailrref, mm, 0);
    lapacke_gtest_init_int_buffer_pair_with_constant(ifaill, ifaillref, mm, 0);
    lapacke_gtest_init_char_buffer_pair_with_constant(select, selectref, n, 0xff);
    
   } /* end of Constructor  */
    

/* Destructor definition  'hsein_dcomplex_common_parameters' */
hsein_dcomplex_parameters :: ~hsein_dcomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   hsein_cplx_free();
} 

//  Test fixture class definition
class zhsein_test  : public  ::testing::Test {
public:
   hsein_dcomplex_parameters  *zhsein_obj;
   void SetUp();  
   void TearDown () { delete zhsein_obj; }
};

void zhsein_test::SetUp()
{
    /* LAPACKE ZHSEIN prototype */
    typedef int (*Fptr_NL_LAPACKE_zhsein) ( int matrix_layout, char side,
			char eigsrc, char initv, lapack_logical* select, lapack_int n,
			const lapack_complex_double* h, lapack_int ldh, 
			lapack_complex_double* w, lapack_complex_double* vl,
			lapack_int ldvl, lapack_complex_double* vr, lapack_int ldvr,
			lapack_int mm, lapack_int* m, lapack_int* ifaill,
			lapack_int* ifailr );
            
    Fptr_NL_LAPACKE_zhsein ZHSEIN;      

    zhsein_obj = new  hsein_dcomplex_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_paramslist[idx].side,
                                         eig_paramslist[idx].eigsrc,
                                         eig_paramslist[idx].initv,
                                         eig_paramslist[idx].n,
                                         eig_paramslist[idx].kb
										 );
    idx = Circular_Increment_Index(idx);

    zhsein_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zhsein_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zhsein_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zhsein_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZHSEIN = (Fptr_NL_LAPACKE_zhsein)dlsym(zhsein_obj->hModule, "LAPACKE_zhsein");
    ASSERT_TRUE(ZHSEIN != NULL) << "failed to ppt the Netlib LAPACKE_zhsein symbol";

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    zhsein_obj->inforef = ZHSEIN(   zhsein_obj->matrix_layout,
                                    zhsein_obj->side,
                                    zhsein_obj->eigsrc, 
                                    zhsein_obj->initv,
                  (lapack_logical *)zhsein_obj->selectref,
                                    zhsein_obj->n,
                                    zhsein_obj->href, 
                                    zhsein_obj->ldh,
                                    zhsein_obj->wref,
                                    zhsein_obj->vlref,
                                    zhsein_obj->ldvl,
                                    zhsein_obj->vrref,
                                    zhsein_obj->ldvr,
                                    zhsein_obj->mm,
                                    &zhsein_obj->mref,
                                    zhsein_obj->ifaillref,
                                    zhsein_obj->ifailrref
									);

    /* Compute libflame's Lapacke o/p  */
    zhsein_obj->info = LAPACKE_zhsein(  zhsein_obj->matrix_layout,
                                    zhsein_obj->side,
                                    zhsein_obj->eigsrc, 
                                    zhsein_obj->initv,
                  (lapack_logical *)zhsein_obj->select,
                                    zhsein_obj->n,
                                    zhsein_obj->h, 
                                    zhsein_obj->ldh,
                                    zhsein_obj->w,
                                    zhsein_obj->vl,
                                    zhsein_obj->ldvl,
                                    zhsein_obj->vr,
                                    zhsein_obj->ldvr,
                                    zhsein_obj->mm,
                                    &zhsein_obj->m,
                                    zhsein_obj->ifaill,
                                    zhsein_obj->ifailr
									);

    /* Capture the Netlib, libflame o/p buffers' differences */
    zhsein_obj->diff_w =  computeDiff_z( zhsein_obj->n, 
                zhsein_obj->w, zhsein_obj->wref );

}

TEST_F(zhsein_test, zhsein1) {
    EXPECT_NEAR(0.0, zhsein_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhsein_test, zhsein2) {
    EXPECT_NEAR(0.0, zhsein_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhsein_test, zhsein3) {
    EXPECT_NEAR(0.0, zhsein_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhsein_test, zhsein4) {
    EXPECT_NEAR(0.0, zhsein_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}


