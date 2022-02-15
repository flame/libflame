#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"

#define larfx_free() \
if (v!=NULL)    free(v); \
if (vref!=NULL) free(vref);\
if (work!=NULL)  free(work);\
if (workref!=NULL) free(workref); \
if (c != NULL) free(c); \
if (cref != NULL) free(cref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
/* Begin float_common_parameters  class definition */
class larfx_float_parameters{

   public:
	int bufsize_c;
	int bufsize_w;
	int bufsize_v;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	float* v, tau;
	char side;
	lapack_int ldc;
	/*Output Parameter*/
	float* work, *c;	
	float *vref, *workref, *cref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      larfx_float_parameters (int matrix_layout, char side , lapack_int m, lapack_int n, float tau);;
      ~larfx_float_parameters ();

};

/* Constructor definition  float_common_parameters */
larfx_float_parameters:: larfx_float_parameters (int matrix_layout_i, char side_i , lapack_int m_i, lapack_int n_i, float tau_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	side = side_i;
	tau = tau_i;
	ldc = m;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n larfx float: matrix_layout = %d, side :%c, m:%d, n: %d \n", matrix_layout, side, m, n);
	#endif
	if (side == 'L')
	{
		bufsize_w = n;
		bufsize_v = m;
	}else if (side == 'R') {
		bufsize_w = m;
		bufsize_v = n;
	}else 
		EXPECT_TRUE(false) << "side is invalid";

	if (matrix_layout == LAPACK_COL_MAJOR)
	{		
		bufsize_c = ldc*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{		bufsize_c = ldc*m;
	}else
	{
		EXPECT_TRUE(false) << "matrix_layout invalid";
	}

	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&v, &vref, bufsize_v);
	lapacke_gtest_alloc_float_buffer_pair(&work, &workref, bufsize_w);	
	lapacke_gtest_alloc_float_buffer_pair(&c, &cref, bufsize_c);
	
	if ((v==NULL) || (vref==NULL) ||
		(work==NULL) || (workref==NULL) ||
		(c==NULL) || (cref==NULL))
	{
		EXPECT_FALSE( true) << "larfx_float_parameters object: malloc error.";
		larfx_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand( v, vref, bufsize_v);
	lapacke_gtest_init_float_buffer_pair_rand( work, workref, bufsize_w);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
larfx_float_parameters :: ~larfx_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" larfx_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   larfx_free();

}
/*  Test fixture class definition */
class slarfx_test  : public  ::testing::Test {
public:
   larfx_float_parameters  *slarfx_obj;
   void SetUp();
   void TearDown () { delete slarfx_obj;}
};

void slarfx_test::SetUp(){
	 /* LAPACKE slarfx prototype */
    typedef int (*Fptr_NL_LAPACKE_slarfx) (int matrix_layout , char side , lapack_int m , lapack_int n ,\
	const float * v , float tau , float * c , lapack_int ldc , float * work);

    Fptr_NL_LAPACKE_slarfx slarfx;

    slarfx_obj = new larfx_float_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].side,
						   eig_paramslist[idx].m,
						   eig_paramslist[idx].n,
                           eig_paramslist[idx].kb);
						   

    idx = Circular_Increment_Index(idx);

    slarfx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    slarfx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(slarfx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(slarfx_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*slarfx library call */
    slarfx = (Fptr_NL_LAPACKE_slarfx)dlsym(slarfx_obj->hModule, "LAPACKE_slarfx");
    ASSERT_TRUE(slarfx != NULL) << "failed to get the Netlib LAPACKE_slarfx symbol";
  
/*Compute slarfx's  o/p */
    slarfx_obj->inforef = slarfx( slarfx_obj->matrix_layout, slarfx_obj->side, slarfx_obj->m, slarfx_obj->n,\
	(const float*)slarfx_obj->vref,slarfx_obj->tau, slarfx_obj->cref, slarfx_obj->ldc, slarfx_obj->workref);

    /* Compute libflame's Lapacke o/p  */
    slarfx_obj->info = LAPACKE_slarfx( slarfx_obj->matrix_layout, slarfx_obj->side, slarfx_obj->m, slarfx_obj->n,\
	(const float*)slarfx_obj->v, slarfx_obj->tau, slarfx_obj->c, slarfx_obj->ldc, slarfx_obj->work);
    if( slarfx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame\
        LAPACKE_slarfx is wrong\n", slarfx_obj->info );
    }
    if( slarfx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_slarfx is wrong\n", 
        slarfx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    slarfx_obj->diff =  computeDiff_s( slarfx_obj->bufsize_c, 
                slarfx_obj->c, slarfx_obj->cref );

}

TEST_F(slarfx_test, slarfx1) {
    EXPECT_NEAR(0.0, slarfx_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(slarfx_test, slarfx2) {
    EXPECT_NEAR(0.0, slarfx_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(slarfx_test, slarfx3) {
    EXPECT_NEAR(0.0, slarfx_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin double_common_parameters  class definition */
class larfx_double_parameters{

   public:
	int bufsize_c;
	int bufsize_w;
	int bufsize_v;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	double* v, tau;
	char side;
	lapack_int ldc;
	/*Output Parameter*/
	double* work, *c;	
	double *vref, *workref, *cref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      larfx_double_parameters (int matrix_layout, char side , lapack_int m, lapack_int n, double tau);;
      ~larfx_double_parameters ();

};

/* Constructor definition  double_common_parameters */
larfx_double_parameters:: larfx_double_parameters (int matrix_layout_i, char side_i , lapack_int m_i, lapack_int n_i, double tau_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	side = side_i;
	tau = tau_i;
	ldc = m;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n larfx double: matrix_layout = %d, side :%c, m:%d, n: %d \n", matrix_layout, side, m, n);
	#endif
	if (side == 'L')
	{
		bufsize_w = n;
		bufsize_v = m;
	}else if (side == 'R') {
		bufsize_w = m;
		bufsize_v = n;
	}else 
		EXPECT_TRUE(false) << "side is invalid";

	if (matrix_layout == LAPACK_COL_MAJOR)
	{		
		bufsize_c = ldc*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{		bufsize_c = ldc*m;
	}else
	{
		EXPECT_TRUE(false) << "matrix_layout invalid";
	}

	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&v, &vref, bufsize_v);
	lapacke_gtest_alloc_double_buffer_pair(&work, &workref, bufsize_w);	
	lapacke_gtest_alloc_double_buffer_pair(&c, &cref, bufsize_c);
	
	if ((v==NULL) || (vref==NULL) ||
		(work==NULL) || (workref==NULL) ||
		(c==NULL) || (cref==NULL))
	{
		EXPECT_FALSE( true) << "larfx_double_parameters object: malloc error.";
		larfx_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand( v, vref, bufsize_v);
	lapacke_gtest_init_double_buffer_pair_rand( work, workref, bufsize_w);

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
larfx_double_parameters :: ~larfx_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" larfx_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   larfx_free();

}
/*  Test fixture class definition */
class dlarfx_test  : public  ::testing::Test {
public:
   larfx_double_parameters  *dlarfx_obj;
   void SetUp();
   void TearDown () { delete dlarfx_obj;}
};

void dlarfx_test::SetUp(){
	 /* LAPACKE dlarfx prototype */
    typedef int (*Fptr_NL_LAPACKE_dlarfx) (int matrix_layout , char side , lapack_int m , lapack_int n ,\
	const double * v , double tau , double * c , lapack_int ldc , double * work);

    Fptr_NL_LAPACKE_dlarfx dlarfx;

    dlarfx_obj = new larfx_double_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].side,
						   eig_paramslist[idx].m,
						   eig_paramslist[idx].n,
                           eig_paramslist[idx].kb);
						   

    idx = Circular_Increment_Index(idx);

    dlarfx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dlarfx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dlarfx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dlarfx_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*dlarfx library call */
    dlarfx = (Fptr_NL_LAPACKE_dlarfx)dlsym(dlarfx_obj->hModule, "LAPACKE_dlarfx");
    ASSERT_TRUE(dlarfx != NULL) << "failed to get the Netlib LAPACKE_dlarfx symbol";
  
/*Compute dlarfx's  o/p */
    dlarfx_obj->inforef = dlarfx( dlarfx_obj->matrix_layout, dlarfx_obj->side, dlarfx_obj->m, dlarfx_obj->n,\
	(const double*)dlarfx_obj->vref,dlarfx_obj->tau, dlarfx_obj->cref, dlarfx_obj->ldc, dlarfx_obj->workref);

    /* Compute libflame's Lapacke o/p  */
    dlarfx_obj->info = LAPACKE_dlarfx( dlarfx_obj->matrix_layout, dlarfx_obj->side, dlarfx_obj->m, dlarfx_obj->n,\
	(const double*)dlarfx_obj->v, dlarfx_obj->tau, dlarfx_obj->c, dlarfx_obj->ldc, dlarfx_obj->work);
    if( dlarfx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame\
        LAPACKE_dlarfx is wrong\n", dlarfx_obj->info );
    }
    if( dlarfx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dlarfx is wrong\n", 
        dlarfx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dlarfx_obj->diff =  computeDiff_d( dlarfx_obj->bufsize_c, 
                dlarfx_obj->c, dlarfx_obj->cref );

}

TEST_F(dlarfx_test, dlarfx1) {
    EXPECT_NEAR(0.0, dlarfx_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dlarfx_test, dlarfx2) {
    EXPECT_NEAR(0.0, dlarfx_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dlarfx_test, dlarfx3) {
    EXPECT_NEAR(0.0, dlarfx_obj->diff, LAPACKE_GTEST_THRESHOLD);
}
/* Begin scomplex_common_parameters  class definition */
class larfx_scomplex_parameters{

   public:
	int bufsize_c;
	int bufsize_w;
	int bufsize_v;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_complex_float* v, tau;
	char side;
	lapack_int ldc;
	/*Output Parameter*/
	lapack_complex_float* work, *c;	
	lapack_complex_float *vref, *workref, *cref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      larfx_scomplex_parameters (int matrix_layout, char side , lapack_int m, lapack_int n, lapack_complex_float tau);;
      ~larfx_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
larfx_scomplex_parameters:: larfx_scomplex_parameters (int matrix_layout_i, char side_i , lapack_int m_i, lapack_int n_i, lapack_complex_float tau_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	side = side_i;
	tau = tau_i;
	ldc = m;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n larfx scomplex: matrix_layout = %d, side :%c, m:%d, n: %d \n", matrix_layout, side, m, n);
	#endif
	if (side == 'L')
	{
		bufsize_w = n;
		bufsize_v = m;
	}else if (side == 'R') {
		bufsize_w = m;
		bufsize_v = n;
	}else 
		EXPECT_TRUE(false) << "side is invalid";

	if (matrix_layout == LAPACK_COL_MAJOR)
	{		
		bufsize_c = ldc*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{		bufsize_c = ldc*m;
	}else
	{
		EXPECT_TRUE(false) << "matrix_layout invalid";
	}

	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&v, &vref, bufsize_v);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&work, &workref, bufsize_w);	
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&c, &cref, bufsize_c);
	
	if ((v==NULL) || (vref==NULL) ||
		(work==NULL) || (workref==NULL) ||
		(c==NULL) || (cref==NULL))
	{
		EXPECT_FALSE( true) << "larfx_float_parameters object: malloc error.";
		larfx_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( v, vref, bufsize_v);
	lapacke_gtest_init_scomplex_buffer_pair_rand( work, workref, bufsize_w);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
larfx_scomplex_parameters :: ~larfx_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" larfx_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   larfx_free();

}
/*  Test fixture class definition */
class clarfx_test  : public  ::testing::Test {
public:
   larfx_scomplex_parameters  *clarfx_obj;
   void SetUp();
   void TearDown () { delete clarfx_obj;}
};

void clarfx_test::SetUp(){
	 /* LAPACKE clarfx prototype */
    typedef int (*Fptr_NL_LAPACKE_clarfx) (int matrix_layout , char side , lapack_int m , lapack_int n ,\
	const lapack_complex_float * v , lapack_complex_float tau , lapack_complex_float * c , lapack_int ldc , lapack_complex_float * work);

    Fptr_NL_LAPACKE_clarfx clarfx;

    clarfx_obj = new larfx_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].side,
						   eig_paramslist[idx].m,
						   eig_paramslist[idx].n,
                           eig_paramslist[idx].kb);
						   

    idx = Circular_Increment_Index(idx);

    clarfx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    clarfx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(clarfx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(clarfx_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*clarfx library call */
    clarfx = (Fptr_NL_LAPACKE_clarfx)dlsym(clarfx_obj->hModule, "LAPACKE_clarfx");
    ASSERT_TRUE(clarfx != NULL) << "failed to get the Netlib LAPACKE_clarfx symbol";
  
/*Compute clarfx's  o/p */
    clarfx_obj->inforef = clarfx( clarfx_obj->matrix_layout, clarfx_obj->side, clarfx_obj->m, clarfx_obj->n,\
	(const lapack_complex_float*)clarfx_obj->vref,clarfx_obj->tau, clarfx_obj->cref, clarfx_obj->ldc, clarfx_obj->workref);

    /* Compute libflame's Lapacke o/p  */
    clarfx_obj->info = LAPACKE_clarfx( clarfx_obj->matrix_layout, clarfx_obj->side, clarfx_obj->m, clarfx_obj->n,\
	(const lapack_complex_float*)clarfx_obj->v, clarfx_obj->tau, clarfx_obj->c, clarfx_obj->ldc, clarfx_obj->work);
    if( clarfx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame\
        LAPACKE_clarfx is wrong\n", clarfx_obj->info );
    }
    if( clarfx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_clarfx is wrong\n", 
        clarfx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    clarfx_obj->diff =  computeDiff_c( clarfx_obj->bufsize_c, 
                clarfx_obj->c, clarfx_obj->cref );

}

TEST_F(clarfx_test, clarfx1) {
    EXPECT_NEAR(0.0, clarfx_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(clarfx_test, clarfx2) {
    EXPECT_NEAR(0.0, clarfx_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(clarfx_test, clarfx3) {
    EXPECT_NEAR(0.0, clarfx_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(clarfx_test, clarfx4) {
    EXPECT_NEAR(0.0, clarfx_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class larfx_dcomplex_parameters{

   public:
	int bufsize_c;
	int bufsize_w;
	int bufsize_v;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_complex_double* v, tau;
	char side;
	lapack_int ldc;
	/*Output Parameter*/
	lapack_complex_double* work, *c;	
	lapack_complex_double *vref, *workref, *cref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      larfx_dcomplex_parameters (int matrix_layout, char side , lapack_int m, lapack_int n, lapack_complex_double tau);;
      ~larfx_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
larfx_dcomplex_parameters:: larfx_dcomplex_parameters (int matrix_layout_i, char side_i , lapack_int m_i, lapack_int n_i, lapack_complex_double tau_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	side = side_i;
	tau = tau_i;
	ldc = m;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n larfx dcomplex: matrix_layout = %d, side :%c, m:%d, n: %d \n", matrix_layout, side, m, n);
	#endif
	if (side == 'L')
	{
		bufsize_w = n;
		bufsize_v = m;
	}else if (side == 'R') {
		bufsize_w = m;
		bufsize_v = n;
	}else 
		EXPECT_TRUE(false) << "side is invalid";

	if (matrix_layout == LAPACK_COL_MAJOR)
	{		
		bufsize_c = ldc*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{		bufsize_c = ldc*m;
	}else
	{
		EXPECT_TRUE(false) << "matrix_layout invalid";
	}

	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&v, &vref, bufsize_v);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&work, &workref, bufsize_w);	
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&c, &cref, bufsize_c);
	
	if ((v==NULL) || (vref==NULL) ||
		(work==NULL) || (workref==NULL) ||
		(c==NULL) || (cref==NULL))
	{
		EXPECT_FALSE( true) << "larfx_double_parameters object: malloc error.";
		larfx_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( v, vref, bufsize_v);
	lapacke_gtest_init_dcomplex_buffer_pair_rand( work, workref, bufsize_w);

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
larfx_dcomplex_parameters :: ~larfx_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" larfx_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   larfx_free();

}
/*  Test fixture class definition */
class zlarfx_test  : public  ::testing::Test {
public:
   larfx_dcomplex_parameters  *zlarfx_obj;
   void SetUp();
   void TearDown () { delete zlarfx_obj;}
};

void zlarfx_test::SetUp(){
	 /* LAPACKE zlarfx prototype */
    typedef int (*Fptr_NL_LAPACKE_zlarfx) (int matrix_layout , char side , lapack_int m , lapack_int n ,\
	const lapack_complex_double * v , lapack_complex_double tau , lapack_complex_double * c , lapack_int ldc , lapack_complex_double * work);

    Fptr_NL_LAPACKE_zlarfx zlarfx;

    zlarfx_obj = new larfx_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].side,
						   eig_paramslist[idx].m,
						   eig_paramslist[idx].n,
                           eig_paramslist[idx].kb);
						   

    idx = Circular_Increment_Index(idx);

    zlarfx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zlarfx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zlarfx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zlarfx_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*zlarfx library call */
    zlarfx = (Fptr_NL_LAPACKE_zlarfx)dlsym(zlarfx_obj->hModule, "LAPACKE_zlarfx");
    ASSERT_TRUE(zlarfx != NULL) << "failed to get the Netlib LAPACKE_zlarfx symbol";
  
/*Compute zlarfx's  o/p */
    zlarfx_obj->inforef = zlarfx( zlarfx_obj->matrix_layout, zlarfx_obj->side, zlarfx_obj->m, zlarfx_obj->n,\
	(const lapack_complex_double*)zlarfx_obj->vref,zlarfx_obj->tau, zlarfx_obj->cref, zlarfx_obj->ldc, zlarfx_obj->workref);

    /* Compute libflame's Lapacke o/p  */
    zlarfx_obj->info = LAPACKE_zlarfx( zlarfx_obj->matrix_layout, zlarfx_obj->side, zlarfx_obj->m, zlarfx_obj->n,\
	(const lapack_complex_double*)zlarfx_obj->v, zlarfx_obj->tau, zlarfx_obj->c, zlarfx_obj->ldc, zlarfx_obj->work);
    if( zlarfx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame\
        LAPACKE_zlarfx is wrong\n", zlarfx_obj->info );
    }
    if( zlarfx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zlarfx is wrong\n", 
        zlarfx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zlarfx_obj->diff =  computeDiff_z( zlarfx_obj->bufsize_c, 
                zlarfx_obj->c, zlarfx_obj->cref );

}

TEST_F(zlarfx_test, zlarfx1) {
    EXPECT_NEAR(0.0, zlarfx_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlarfx_test, zlarfx2) {
    EXPECT_NEAR(0.0, zlarfx_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlarfx_test, zlarfx3) {
    EXPECT_NEAR(0.0, zlarfx_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlarfx_test, zlarfx4) {
    EXPECT_NEAR(0.0, zlarfx_obj->diff, LAPACKE_GTEST_THRESHOLD);
}