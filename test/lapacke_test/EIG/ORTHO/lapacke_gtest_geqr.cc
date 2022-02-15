#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"


#define geqr_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (t!=NULL)  free(t);\
if (tref!=NULL) free(tref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
/* Begin float_common_parameters  class definition */
class geqr_float_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	float* A;
	lapack_int lda;
	int tsize;
	/*Output Parameter*/
	float* t;
	float *Aref, *tref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      geqr_float_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda, int tsize);
      ~geqr_float_parameters ();

};

/* Constructor definition  float_common_parameters */
geqr_float_parameters:: geqr_float_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i, int tsize_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	tsize = tsize_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n geqr float:  m: %d, n: %d lda: %d, tsize:%d \n",  m, n, lda, tsize);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
		bufsize = lda*n;
	else if (matrix_layout == LAPACK_ROW_MAJOR)
		bufsize = lda*m;
	else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_float_buffer_pair(&t, &tref, tsize);
	
	if ((A==NULL) || (Aref==NULL) ||\
		(t==NULL) || (tref==NULL)){
		EXPECT_FALSE( true) << "geqr_float_parameters object: malloc error.";
		geqr_free();
		exit(0);
	}	
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand( A, Aref, bufsize);
	

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
geqr_float_parameters :: ~geqr_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" geqr_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   geqr_free();

}
/*  Test fixture class definition */
class sgeqrest  : public  ::testing::Test {
public:
   geqr_float_parameters  *sgeqr_obj;
   void SetUp();
   void TearDown () { delete sgeqr_obj;}
};

void sgeqrest::SetUp(){

    /* LAPACKE sgeqr prototype */
    typedef int (*Fptr_NL_LAPACKE_sgeqr) (int matrix_layout, lapack_int m,lapack_int n, 
											float *A, lapack_int lda, float* t, lapack_int tsize);

    Fptr_NL_LAPACKE_sgeqr sgeqr;

    sgeqr_obj = new geqr_float_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda,
						   eig_paramslist[idx].lda); 
	
	idx = Circular_Increment_Index(idx);

    sgeqr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgeqr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgeqr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgeqr_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    sgeqr = (Fptr_NL_LAPACKE_sgeqr)dlsym(sgeqr_obj->hModule, "LAPACKE_sgeqr");
    ASSERT_TRUE(sgeqr != NULL) << "failed to get the Netlib LAPACKE_sgeqr symbol";
    

    sgeqr_obj->inforef = sgeqr( sgeqr_obj->matrix_layout, sgeqr_obj->m,
								sgeqr_obj->n,sgeqr_obj->Aref,
								sgeqr_obj->lda, sgeqr_obj->tref, sgeqr_obj->tsize);

    /* Compute libflame's Lapacke o/p  */
    sgeqr_obj->info = LAPACKE_sgeqr( sgeqr_obj->matrix_layout, sgeqr_obj->m,
										sgeqr_obj->n,sgeqr_obj->A, 
										sgeqr_obj->lda, sgeqr_obj->t, sgeqr_obj->tsize);
	
	if( sgeqr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_sgeqr is wrong\n", sgeqr_obj->info );
    }
    if( sgeqr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sgeqr is wrong\n", 
        sgeqr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    sgeqr_obj->diff =  computeDiff_s( sgeqr_obj->bufsize, 
                sgeqr_obj->A, sgeqr_obj->Aref );

}

TEST_F(sgeqrest, sgeqr1) {
    EXPECT_NEAR(0.0, sgeqr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgeqrest, sgeqr2) {
    EXPECT_NEAR(0.0, sgeqr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgeqrest, sgeqr3) {
    EXPECT_NEAR(0.0, sgeqr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgeqrest, sgeqr4) {
    EXPECT_NEAR(0.0, sgeqr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin double_common_parameters  class definition */
class geqr_double_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	double* A;
	lapack_int lda;
	int tsize;
	/*Output Parameter*/
	double* t;
	double *Aref, *tref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      geqr_double_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda, int tsize);
      ~geqr_double_parameters ();

};

/* Constructor definition  double_common_parameters */
geqr_double_parameters:: geqr_double_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i, int tsize_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	tsize = tsize_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n geqr double:  m: %d, n: %d lda: %d, tsize:%d \n",  m, n, lda, tsize);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
		bufsize = lda*n;
	else if (matrix_layout == LAPACK_ROW_MAJOR)
		bufsize = lda*m;
	else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_double_buffer_pair(&t, &tref, tsize);
	
	if ((A==NULL) || (Aref==NULL) ||\
		(t==NULL) || (tref==NULL)){
		EXPECT_FALSE( true) << "geqr_double_parameters object: malloc error.";
		geqr_free();
		exit(0);
	}	
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand( A, Aref, bufsize);
	

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
geqr_double_parameters :: ~geqr_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" geqr_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   geqr_free();

}
/*  Test fixture class definition */
class dgeqrest  : public  ::testing::Test {
public:
   geqr_double_parameters  *dgeqr_obj;
   void SetUp();
   void TearDown () { delete dgeqr_obj;}
};

void dgeqrest::SetUp(){

    /* LAPACKE dgeqr prototype */
    typedef int (*Fptr_NL_LAPACKE_dgeqr) (int matrix_layout, lapack_int m,lapack_int n, 
											double *A, lapack_int lda, double* t, lapack_int tsize);

    Fptr_NL_LAPACKE_dgeqr dgeqr;

    dgeqr_obj = new geqr_double_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda,
						   eig_paramslist[idx].lda); 
	
	idx = Circular_Increment_Index(idx);

    dgeqr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgeqr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgeqr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgeqr_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    dgeqr = (Fptr_NL_LAPACKE_dgeqr)dlsym(dgeqr_obj->hModule, "LAPACKE_dgeqr");
    ASSERT_TRUE(dgeqr != NULL) << "failed to get the Netlib LAPACKE_dgeqr symbol";
    

    dgeqr_obj->inforef = dgeqr( dgeqr_obj->matrix_layout, dgeqr_obj->m,
								dgeqr_obj->n,dgeqr_obj->Aref,
								dgeqr_obj->lda, dgeqr_obj->tref, dgeqr_obj->tsize);

    /* Compute libflame's Lapacke o/p  */
    dgeqr_obj->info = LAPACKE_dgeqr( dgeqr_obj->matrix_layout, dgeqr_obj->m,
										dgeqr_obj->n,dgeqr_obj->A, 
										dgeqr_obj->lda, dgeqr_obj->t, dgeqr_obj->tsize);
	
	if( dgeqr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dgeqr is wrong\n", dgeqr_obj->info );
    }
    if( dgeqr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dgeqr is wrong\n", 
        dgeqr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dgeqr_obj->diff =  computeDiff_d( dgeqr_obj->bufsize, 
                dgeqr_obj->A, dgeqr_obj->Aref );

}

TEST_F(dgeqrest, dgeqr1) {
    EXPECT_NEAR(0.0, dgeqr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgeqrest, dgeqr2) {
    EXPECT_NEAR(0.0, dgeqr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgeqrest, dgeqr3) {
    EXPECT_NEAR(0.0, dgeqr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgeqrest, dgeqr4) {
    EXPECT_NEAR(0.0, dgeqr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin scomplex_common_parameters  class definition */
class geqr_scomplex_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_complex_float* A;
	lapack_int lda;
	int tsize;
	/*Output Parameter*/
	lapack_complex_float* t;
	lapack_complex_float *Aref, *tref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      geqr_scomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda, int tsize);
      ~geqr_scomplex_parameters ();

};

/* Constructor definition  scomplex_common_parameters */
geqr_scomplex_parameters:: geqr_scomplex_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i, int tsize_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	tsize = tsize_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n geqr scomplex:  m: %d, n: %d lda: %d, tsize:%d \n",  m, n, lda, tsize);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
		bufsize = lda*n;
	else if (matrix_layout == LAPACK_ROW_MAJOR)
		bufsize = lda*m;
	else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&t, &tref, tsize);
	
	if ((A==NULL) || (Aref==NULL) ||\
		(t==NULL) || (tref==NULL)){
		EXPECT_FALSE( true) << "geqr_scomplex_parameters object: malloc error.";
		geqr_free();
		exit(0);
	}	
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( A, Aref, bufsize);
	

} /* end of Constructor  */

/* Destructor definition  'scomplex_common_parameters' class  */
geqr_scomplex_parameters :: ~geqr_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" geqr_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   geqr_free();

}
/*  Test fixture class definition */
class cgeqrest  : public  ::testing::Test {
public:
   geqr_scomplex_parameters  *cgeqr_obj;
   void SetUp();
   void TearDown () { delete cgeqr_obj;}
};

void cgeqrest::SetUp(){

    /* LAPACKE cgeqr prototype */
    typedef int (*Fptr_NL_LAPACKE_cgeqr) (int matrix_layout, lapack_int m,lapack_int n, 
											lapack_complex_float *A, lapack_int lda, lapack_complex_float* t, lapack_int tsize);

    Fptr_NL_LAPACKE_cgeqr cgeqr;

    cgeqr_obj = new geqr_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda,
						   eig_paramslist[idx].lda); 
	
	idx = Circular_Increment_Index(idx);

    cgeqr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgeqr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgeqr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgeqr_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    cgeqr = (Fptr_NL_LAPACKE_cgeqr)dlsym(cgeqr_obj->hModule, "LAPACKE_cgeqr");
    ASSERT_TRUE(cgeqr != NULL) << "failed to get the Netlib LAPACKE_cgeqr symbol";
    

    cgeqr_obj->inforef = cgeqr( cgeqr_obj->matrix_layout, cgeqr_obj->m,
								cgeqr_obj->n,cgeqr_obj->Aref,
								cgeqr_obj->lda, cgeqr_obj->tref, cgeqr_obj->tsize);

    /* Compute libflame's Lapacke o/p  */
    cgeqr_obj->info = LAPACKE_cgeqr( cgeqr_obj->matrix_layout, cgeqr_obj->m,
										cgeqr_obj->n,cgeqr_obj->A, 
										cgeqr_obj->lda, cgeqr_obj->t, cgeqr_obj->tsize);
	
	if( cgeqr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cgeqr is wrong\n", cgeqr_obj->info );
    }
    if( cgeqr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cgeqr is wrong\n", 
        cgeqr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cgeqr_obj->diff =  computeDiff_c( cgeqr_obj->bufsize, 
                cgeqr_obj->A, cgeqr_obj->Aref );

}

TEST_F(cgeqrest, cgeqr1) {
    EXPECT_NEAR(0.0, cgeqr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgeqrest, cgeqr2) {
    EXPECT_NEAR(0.0, cgeqr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgeqrest, cgeqr3) {
    EXPECT_NEAR(0.0, cgeqr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgeqrest, cgeqr4) {
    EXPECT_NEAR(0.0, cgeqr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class geqr_dcomplex_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_complex_double* A;
	lapack_int lda;
	int tsize;
	/*Output Parameter*/
	lapack_complex_double* t;
	lapack_complex_double *Aref, *tref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      geqr_dcomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda, int tsize);
      ~geqr_dcomplex_parameters ();

};

/* Constructor definition  dcomplex_common_parameters */
geqr_dcomplex_parameters:: geqr_dcomplex_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i, int tsize_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	tsize = tsize_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n geqr dcomplex:  m: %d, n: %d lda: %d, tsize:%d \n",  m, n, lda, tsize);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
		bufsize = lda*n;
	else if (matrix_layout == LAPACK_ROW_MAJOR)
		bufsize = lda*m;
	else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&t, &tref, tsize);
	
	if ((A==NULL) || (Aref==NULL)||\
		(t==NULL) || (tref==NULL)){
		EXPECT_FALSE( true) << "geqr_dcomplex_parameters object: malloc error.";
		geqr_free();
		exit(0);
	}	
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( A, Aref, bufsize);
	

} /* end of Constructor  */

/* Destructor definition  'dcomplex_common_parameters' class  */
geqr_dcomplex_parameters :: ~geqr_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" geqr_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   geqr_free();

}
/*  Test fixture class definition */
class zgeqrest  : public  ::testing::Test {
public:
   geqr_dcomplex_parameters  *zgeqr_obj;
   void SetUp();
   void TearDown () { delete zgeqr_obj;}
};

void zgeqrest::SetUp(){

    /* LAPACKE zgeqr prototype */
    typedef int (*Fptr_NL_LAPACKE_zgeqr) (int matrix_layout, lapack_int m,lapack_int n, 
											lapack_complex_double *A, lapack_int lda, lapack_complex_double* t, lapack_int tsize);

    Fptr_NL_LAPACKE_zgeqr zgeqr;

    zgeqr_obj = new geqr_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda,
						   eig_paramslist[idx].lda); 
	
	idx = Circular_Increment_Index(idx);

    zgeqr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgeqr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgeqr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgeqr_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zgeqr = (Fptr_NL_LAPACKE_zgeqr)dlsym(zgeqr_obj->hModule, "LAPACKE_zgeqr");
    ASSERT_TRUE(zgeqr != NULL) << "failed to get the Netlib LAPACKE_zgeqr symbol";
    

    zgeqr_obj->inforef = zgeqr( zgeqr_obj->matrix_layout, zgeqr_obj->m,
								zgeqr_obj->n,zgeqr_obj->Aref,
								zgeqr_obj->lda, zgeqr_obj->tref, zgeqr_obj->tsize);

    /* Compute libflame's Lapacke o/p  */
    zgeqr_obj->info = LAPACKE_zgeqr( zgeqr_obj->matrix_layout, zgeqr_obj->m,
										zgeqr_obj->n,zgeqr_obj->A, 
										zgeqr_obj->lda, zgeqr_obj->t, zgeqr_obj->tsize);
	
	if( zgeqr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zgeqr is wrong\n", zgeqr_obj->info );
    }
    if( zgeqr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zgeqr is wrong\n", 
        zgeqr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zgeqr_obj->diff =  computeDiff_z( zgeqr_obj->bufsize, 
                zgeqr_obj->A, zgeqr_obj->Aref );

}

TEST_F(zgeqrest, zgeqr1) {
    EXPECT_NEAR(0.0, zgeqr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgeqrest, zgeqr2) {
    EXPECT_NEAR(0.0, zgeqr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgeqrest, zgeqr3) {
    EXPECT_NEAR(0.0, zgeqr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgeqrest, zgeqr4) {
    EXPECT_NEAR(0.0, zgeqr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


