#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"


#define gelq_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (t!=NULL)  free(t);\
if (tref!=NULL) free(tref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
/* Begin float_common_parameters  class definition */
class gelq_float_parameters{

   public:
	int bufsize;
	int bufsize_t;
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
      gelq_float_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda, int tsize);
      ~gelq_float_parameters ();

};

/* Constructor definition  float_common_parameters */
gelq_float_parameters:: gelq_float_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i, int tsize_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	tsize = tsize_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n gelq float:  m: %d, n: %d lda: %d, tsize:%d \n",  m, n, lda, tsize);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR) {
		bufsize = lda*n;
		tsize = bufsize;
	}else if (matrix_layout == LAPACK_ROW_MAJOR) {
		bufsize = lda*m;
		tsize = bufsize;
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";

	
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_float_buffer_pair(&t, &tref, tsize);
	
	if ((A==NULL) || (Aref==NULL) ||\
		(t==NULL) || (tref==NULL)){
		EXPECT_FALSE( true) << "gelq_float_parameters object: malloc error.";
		gelq_free();
		exit(0);
	}	
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand( A, Aref, bufsize);
	

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
gelq_float_parameters :: ~gelq_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" gelq_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   gelq_free();

}
/*  Test fixture class definition */
class sgelqest  : public  ::testing::Test {
public:
   gelq_float_parameters  *sgelq_obj;
   void SetUp();
   void TearDown () { delete sgelq_obj;}
};

void sgelqest::SetUp(){

    /* LAPACKE sgelq prototype */
    typedef int (*Fptr_NL_LAPACKE_sgelq) (int matrix_layout, lapack_int m,lapack_int n, 
											float *A, lapack_int lda, float* t, lapack_int tsize);

    Fptr_NL_LAPACKE_sgelq sgelq;

    sgelq_obj = new gelq_float_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda,
						   eig_paramslist[idx].m); 
	
	idx = Circular_Increment_Index(idx);

    sgelq_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgelq_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgelq_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgelq_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    sgelq = (Fptr_NL_LAPACKE_sgelq)dlsym(sgelq_obj->hModule, "LAPACKE_sgelq");
    ASSERT_TRUE(sgelq != NULL) << "failed to get the Netlib LAPACKE_sgelq symbol";
    

    sgelq_obj->inforef = sgelq( sgelq_obj->matrix_layout, sgelq_obj->m,
								sgelq_obj->n,sgelq_obj->Aref,
								sgelq_obj->lda, sgelq_obj->tref, sgelq_obj->tsize);

    /* Compute libflame's Lapacke o/p  */
    sgelq_obj->info = LAPACKE_sgelq( sgelq_obj->matrix_layout, sgelq_obj->m,
										sgelq_obj->n,sgelq_obj->A, 
										sgelq_obj->lda, sgelq_obj->t, sgelq_obj->tsize);
	
	if( sgelq_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_sgelq is wrong\n", sgelq_obj->info );
    }
    if( sgelq_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sgelq is wrong\n", 
        sgelq_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    sgelq_obj->diff =  computeDiff_s( sgelq_obj->bufsize, 
                sgelq_obj->A, sgelq_obj->Aref );

}

TEST_F(sgelqest, sgelq1) {
    EXPECT_NEAR(0.0, sgelq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgelqest, sgelq2) {
    EXPECT_NEAR(0.0, sgelq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgelqest, sgelq3) {
    EXPECT_NEAR(0.0, sgelq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgelqest, sgelq4) {
    EXPECT_NEAR(0.0, sgelq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin double_common_parameters  class definition */
class gelq_double_parameters{

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
      gelq_double_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda, int tsize);
      ~gelq_double_parameters ();

};

/* Constructor definition  double_common_parameters */
gelq_double_parameters:: gelq_double_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i, int tsize_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	tsize = tsize_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n gelq double:  m: %d, n: %d lda: %d, tsize:%d \n",  m, n, lda, tsize);
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
		EXPECT_FALSE( true) << "gelq_double_parameters object: malloc error.";
		gelq_free();
		exit(0);
	}	
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand( A, Aref, bufsize);
	

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
gelq_double_parameters :: ~gelq_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" gelq_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   gelq_free();

}
/*  Test fixture class definition */
class dgelqest  : public  ::testing::Test {
public:
   gelq_double_parameters  *dgelq_obj;
   void SetUp();
   void TearDown () { delete dgelq_obj;}
};

void dgelqest::SetUp(){

    /* LAPACKE dgelq prototype */
    typedef int (*Fptr_NL_LAPACKE_dgelq) (int matrix_layout, lapack_int m,lapack_int n, 
											double *A, lapack_int lda, double* t, lapack_int tsize);

    Fptr_NL_LAPACKE_dgelq dgelq;

    dgelq_obj = new gelq_double_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].lda,
						   eig_paramslist[idx].lda); 
	
	idx = Circular_Increment_Index(idx);

    dgelq_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgelq_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgelq_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgelq_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    dgelq = (Fptr_NL_LAPACKE_dgelq)dlsym(dgelq_obj->hModule, "LAPACKE_dgelq");
    ASSERT_TRUE(dgelq != NULL) << "failed to get the Netlib LAPACKE_dgelq symbol";
    

    dgelq_obj->inforef = dgelq( dgelq_obj->matrix_layout, dgelq_obj->m,
								dgelq_obj->n,dgelq_obj->Aref,
								dgelq_obj->lda, dgelq_obj->tref, dgelq_obj->tsize);

    /* Compute libflame's Lapacke o/p  */
    dgelq_obj->info = LAPACKE_dgelq( dgelq_obj->matrix_layout, dgelq_obj->m,
										dgelq_obj->n,dgelq_obj->A, 
										dgelq_obj->lda, dgelq_obj->t, dgelq_obj->tsize);
	
	if( dgelq_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dgelq is wrong\n", dgelq_obj->info );
    }
    if( dgelq_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dgelq is wrong\n", 
        dgelq_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dgelq_obj->diff =  computeDiff_d( dgelq_obj->bufsize, 
                dgelq_obj->A, dgelq_obj->Aref );

}

TEST_F(dgelqest, dgelq1) {
    EXPECT_NEAR(0.0, dgelq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgelqest, dgelq2) {
    EXPECT_NEAR(0.0, dgelq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgelqest, dgelq3) {
    EXPECT_NEAR(0.0, dgelq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgelqest, dgelq4) {
    EXPECT_NEAR(0.0, dgelq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin scomplex_common_parameters  class definition */
class gelq_scomplex_parameters{

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
      gelq_scomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda, int tsize);
      ~gelq_scomplex_parameters ();

};

/* Constructor definition  scomplex_common_parameters */
gelq_scomplex_parameters:: gelq_scomplex_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i, int tsize_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	tsize = tsize_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n gelq scomplex:  m: %d, n: %d lda: %d, tsize:%d \n",  m, n, lda, tsize);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR) {
		bufsize = lda*n;
		tsize = bufsize;
	}else if (matrix_layout == LAPACK_ROW_MAJOR) {
		bufsize = lda*m;
		tsize = bufsize;
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";


	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&t, &tref, tsize);
	
	if ((A==NULL) || (Aref==NULL) ||\
		(t==NULL) || (tref==NULL)){
		EXPECT_FALSE( true) << "gelq_scomplex_parameters object: malloc error.";
		gelq_free();
		exit(0);
	}	
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( A, Aref, bufsize);
	

} /* end of Constructor  */

/* Destructor definition  'scomplex_common_parameters' class  */
gelq_scomplex_parameters :: ~gelq_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" gelq_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   gelq_free();

}
/*  Test fixture class definition */
class cgelqest  : public  ::testing::Test {
public:
   gelq_scomplex_parameters  *cgelq_obj;
   void SetUp();
   void TearDown () { delete cgelq_obj;}
};

void cgelqest::SetUp(){

    /* LAPACKE cgelq prototype */
    typedef int (*Fptr_NL_LAPACKE_cgelq) (int matrix_layout, lapack_int m,lapack_int n, 
											lapack_complex_float *A, lapack_int lda, lapack_complex_float* t, lapack_int tsize);

    Fptr_NL_LAPACKE_cgelq cgelq;

    cgelq_obj = new gelq_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].lda,
						   eig_paramslist[idx].lda); 
	
	idx = Circular_Increment_Index(idx);

    cgelq_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgelq_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgelq_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgelq_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    cgelq = (Fptr_NL_LAPACKE_cgelq)dlsym(cgelq_obj->hModule, "LAPACKE_cgelq");
    ASSERT_TRUE(cgelq != NULL) << "failed to get the Netlib LAPACKE_cgelq symbol";
    

    cgelq_obj->inforef = cgelq( cgelq_obj->matrix_layout, cgelq_obj->m,
								cgelq_obj->n,cgelq_obj->Aref,
								cgelq_obj->lda, cgelq_obj->tref, cgelq_obj->tsize);

    /* Compute libflame's Lapacke o/p  */
    cgelq_obj->info = LAPACKE_cgelq( cgelq_obj->matrix_layout, cgelq_obj->m,
										cgelq_obj->n,cgelq_obj->A, 
										cgelq_obj->lda, cgelq_obj->t, cgelq_obj->tsize);
	
	if( cgelq_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cgelq is wrong\n", cgelq_obj->info );
    }
    if( cgelq_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cgelq is wrong\n", 
        cgelq_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cgelq_obj->diff =  computeDiff_c( cgelq_obj->bufsize, 
                cgelq_obj->A, cgelq_obj->Aref );

}

TEST_F(cgelqest, cgelq1) {
    EXPECT_NEAR(0.0, cgelq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgelqest, cgelq2) {
    EXPECT_NEAR(0.0, cgelq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgelqest, cgelq3) {
    EXPECT_NEAR(0.0, cgelq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgelqest, cgelq4) {
    EXPECT_NEAR(0.0, cgelq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class gelq_dcomplex_parameters{

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
      gelq_dcomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda, int tsize);
      ~gelq_dcomplex_parameters ();

};

/* Constructor definition  dcomplex_common_parameters */
gelq_dcomplex_parameters:: gelq_dcomplex_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i, int tsize_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	tsize = tsize_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n gelq dcomplex:  m: %d, n: %d lda: %d, tsize:%d \n",  m, n, lda, tsize);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR) {
		bufsize = lda*n;
		tsize =  bufsize;
	}else if (matrix_layout == LAPACK_ROW_MAJOR) {
		bufsize = lda*m;
		tsize = bufsize;
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";

	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&t, &tref, tsize);
	
	if ((A==NULL) || (Aref==NULL)||\
		(t==NULL) || (tref==NULL)){
		EXPECT_FALSE( true) << "gelq_dcomplex_parameters object: malloc error.";
		gelq_free();
		exit(0);
	}	
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( A, Aref, bufsize);
	

} /* end of Constructor  */

/* Destructor definition  'dcomplex_common_parameters' class  */
gelq_dcomplex_parameters :: ~gelq_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" gelq_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   gelq_free();

}
/*  Test fixture class definition */
class zgelqest  : public  ::testing::Test {
public:
   gelq_dcomplex_parameters  *zgelq_obj;
   void SetUp();
   void TearDown () { delete zgelq_obj;}
};

void zgelqest::SetUp(){

    /* LAPACKE zgelq prototype */
    typedef int (*Fptr_NL_LAPACKE_zgelq) (int matrix_layout, lapack_int m,lapack_int n, 
											lapack_complex_double *A, lapack_int lda, lapack_complex_double* t, lapack_int tsize);

    Fptr_NL_LAPACKE_zgelq zgelq;

    zgelq_obj = new gelq_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda,
						   eig_paramslist[idx].m); 
	
	idx = Circular_Increment_Index(idx);

    zgelq_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgelq_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgelq_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgelq_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zgelq = (Fptr_NL_LAPACKE_zgelq)dlsym(zgelq_obj->hModule, "LAPACKE_zgelq");
    ASSERT_TRUE(zgelq != NULL) << "failed to get the Netlib LAPACKE_zgelq symbol";
    

    zgelq_obj->inforef = zgelq( zgelq_obj->matrix_layout, zgelq_obj->m,
								zgelq_obj->n,zgelq_obj->Aref,
								zgelq_obj->lda, zgelq_obj->tref, zgelq_obj->tsize);

    /* Compute libflame's Lapacke o/p  */
    zgelq_obj->info = LAPACKE_zgelq( zgelq_obj->matrix_layout, zgelq_obj->m,
										zgelq_obj->n,zgelq_obj->A, 
										zgelq_obj->lda, zgelq_obj->t, zgelq_obj->tsize);
	
	if( zgelq_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zgelq is wrong\n", zgelq_obj->info );
    }
    if( zgelq_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zgelq is wrong\n", 
        zgelq_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zgelq_obj->diff =  computeDiff_z( zgelq_obj->bufsize, 
                zgelq_obj->A, zgelq_obj->Aref );

}

TEST_F(zgelqest, zgelq1) {
    EXPECT_NEAR(0.0, zgelq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgelqest, zgelq2) {
    EXPECT_NEAR(0.0, zgelq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgelqest, zgelq3) {
    EXPECT_NEAR(0.0, zgelq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgelqest, zgelq4) {
    EXPECT_NEAR(0.0, zgelq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


