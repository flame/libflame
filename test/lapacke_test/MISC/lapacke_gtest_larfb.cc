#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"

#define larfb_free() \
if (v!=NULL)    free(v); \
if (vref!=NULL) free(vref);\
if (t!=NULL)  free(t);\
if (tref!=NULL) free(tref); \
if (c != NULL) free(c); \
if (cref != NULL) free(cref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
/* Begin float_common_parameters  class definition */
class larfb_float_parameters{

   public:
	int bufsize_c;
	int bufsize_t;
	int bufsize_v;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int k;
	lapack_int ldc, ldv, ldt;
	float* v;
	char side, trans, direct, storev;
	/*Output Parameter*/
	float* t, *c;
	float *vref, *tref, *cref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      larfb_float_parameters (int matrix_layout, char side, char trans,  lapack_int m, lapack_int n, lapack_int k);
      ~larfb_float_parameters ();

};

/* Constructor definition  float_common_parameters */
larfb_float_parameters:: larfb_float_parameters (int matrix_layout_i, char side_i, char trans_i , lapack_int m_i, lapack_int n_i, lapack_int k_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	k = k_i;
	ldt = k;
	side = side_i;
	trans = trans_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n larfb float: matrix_layout = %d,  m:%d, n: %d, trans:%c \n", matrix_layout, side, m, n, trans);
	#endif

	if (matrix_layout == LAPACK_COL_MAJOR)
	{	/* checking sizes with storev*/	
		storev = 'C';
		if ( side == 'L')
		{	
			ldv = m;
			bufsize_v = ldv*k;
		} else if (side == 'R')
		{
			ldv = n;
			bufsize_v = ldv*k;			
		}
		/* checking sizes with storev*/	
		storev = 'R';
		if ( side == 'L')
		{	
			ldv = k;
			bufsize_v = ldv*m;
		} else if (side == 'R')
		{
			ldv = k;
			bufsize_v = ldv*n;			
		}
		direct = 'F';
		ldc = m;
		bufsize_c = ldc*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{		
		storev = 'R';
		if ( side == 'L')
		{	
			ldv = m;
			bufsize_v = ldv*k;
		} else if (side == 'R')
		{
			ldv = n;
			bufsize_v = ldv*k;			
		}
		/* checking sizes with storev*/	
		storev = 'C';
		if ( side == 'L')
		{	
			ldv = k;
			bufsize_v = ldv*m;
		} else if (side == 'R')
		{
			ldv = k;
			bufsize_v = ldv*n;			
		}
		direct = 'B';
		ldc = n;
		bufsize_c = ldc*m;
	}else
	{
		EXPECT_TRUE(false) << "matrix_layout invalid";
	}
	
	bufsize_t = ldt*k;
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&v, &vref, bufsize_v);
	lapacke_gtest_alloc_float_buffer_pair(&t, &tref, bufsize_t);	
	lapacke_gtest_alloc_float_buffer_pair(&c, &cref, bufsize_c);
	
	if ((v==NULL) || (vref==NULL) ||
		(t==NULL) || (tref==NULL) ||
		(c==NULL) || (cref==NULL))
	{
		EXPECT_FALSE( true) << "larfb_float_parameters object: malloc error.";
		larfb_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand( v, vref, bufsize_v);
	lapacke_gtest_init_float_buffer_pair_rand_custom_matrix( t, tref, ldt, k, 'S');

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
larfb_float_parameters :: ~larfb_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" larfb_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   larfb_free();

}
/*  Test fixture class definition */
class slarfb_test  : public  ::testing::Test {
public:
   larfb_float_parameters  *slarfb_obj;
   void SetUp();
   void TearDown () { delete slarfb_obj;}
};

void slarfb_test::SetUp(){
	 /* LAPACKE slarfb prototype */
    typedef int (*Fptr_NL_LAPACKE_slarfb) (int matrix_layout , char side , char trans , char direct , char storev ,\
	lapack_int m , lapack_int n , lapack_int k , const float * v , lapack_int ldv , const float * t , lapack_int ldt , float * c , lapack_int ldc);

    Fptr_NL_LAPACKE_slarfb slarfb;

    slarfb_obj = new larfb_float_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].side,
						   eig_paramslist[idx].trans,
						   eig_paramslist[idx].m,
						   eig_paramslist[idx].n,
                           eig_paramslist[idx].n);
						   

    idx = Circular_Increment_Index(idx);

    slarfb_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    slarfb_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(slarfb_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(slarfb_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*slarfb library call */
    slarfb = (Fptr_NL_LAPACKE_slarfb)dlsym(slarfb_obj->hModule, "LAPACKE_slarfb");
    ASSERT_TRUE(slarfb != NULL) << "failed to get the Netlib LAPACKE_slarfb symbol";
  
/*Compute slarfb's  o/p */
    slarfb_obj->inforef = slarfb( slarfb_obj->matrix_layout, slarfb_obj->side, slarfb_obj->trans, slarfb_obj->direct, slarfb_obj->storev,\
	slarfb_obj->m, slarfb_obj->n,slarfb_obj->k, slarfb_obj->vref,slarfb_obj->ldv, (const float*)slarfb_obj->tref, slarfb_obj->ldt,
	slarfb_obj->cref, slarfb_obj->ldc);

    /* Compute libflame's Lapacke o/p  */
    slarfb_obj->info = LAPACKE_slarfb( slarfb_obj->matrix_layout, slarfb_obj->side, slarfb_obj->trans, slarfb_obj->direct, slarfb_obj->storev,\
	slarfb_obj->m, slarfb_obj->n, slarfb_obj->k, slarfb_obj->v, slarfb_obj->ldv, (const float*)slarfb_obj->t, slarfb_obj->ldt, \
	slarfb_obj->c, slarfb_obj->ldc);
	
    if( slarfb_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame\
        LAPACKE_slarfb is wrong\n", slarfb_obj->info );
    }
    if( slarfb_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_slarfb is wrong\n", 
        slarfb_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    slarfb_obj->diff =  computeDiff_s( slarfb_obj->bufsize_c, 
                slarfb_obj->c, slarfb_obj->cref );

}

TEST_F(slarfb_test, slarfb1) {
    EXPECT_NEAR(0.0, slarfb_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(slarfb_test, slarfb2) {
    EXPECT_NEAR(0.0, slarfb_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(slarfb_test, slarfb3) {
    EXPECT_NEAR(0.0, slarfb_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin double_common_parameters  class definition */
class larfb_double_parameters{

   public:
	int bufsize_c;
	int bufsize_t;
	int bufsize_v;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int k;
	lapack_int ldc, ldv, ldt;
	double* v;
	char side, trans, direct, storev;
	/*Output Parameter*/
	double* t, *c;
	double *vref, *tref, *cref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      larfb_double_parameters (int matrix_layout, char side, char trans,  lapack_int m, lapack_int n, lapack_int k);
      ~larfb_double_parameters ();

};

/* Constructor definition  double_common_parameters */
larfb_double_parameters:: larfb_double_parameters (int matrix_layout_i, char side_i, char trans_i , lapack_int m_i, lapack_int n_i, lapack_int k_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	k = k_i;
	ldt = k;
	side = side_i;
	trans = trans_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n larfb double: matrix_layout = %d,  m:%d, n: %d \n", matrix_layout, side, m, n);
	#endif

	if (matrix_layout == LAPACK_COL_MAJOR)
	{	/* checking sizes with storev*/	
		storev = 'C';
		if ( side == 'L')
		{	
			ldv = m;
			bufsize_v = ldv*k;
		} else if (side == 'R')
		{
			ldv = n;
			bufsize_v = ldv*k;			
		}
		/* checking sizes with storev*/	
		storev = 'R';
		if ( side == 'L')
		{	
			ldv = k;
			bufsize_v = ldv*m;
		} else if (side == 'R')
		{
			ldv = k;
			bufsize_v = ldv*n;			
		}
		direct = 'F';
		ldc = m;
		bufsize_c = ldc*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{		
		storev = 'R';
		if ( side == 'L')
		{	
			ldv = m;
			bufsize_v = ldv*k;
		} else if (side == 'R')
		{
			ldv = n;
			bufsize_v = ldv*k;			
		}
		/* checking sizes with storev*/	
		storev = 'C';
		if ( side == 'L')
		{	
			ldv = k;
			bufsize_v = ldv*m;
		} else if (side == 'R')
		{
			ldv = k;
			bufsize_v = ldv*n;			
		}
		direct = 'B';
		ldc = n;
		bufsize_c = ldc*m;
	}else
	{
		EXPECT_TRUE(false) << "matrix_layout invalid";
	}
	
	bufsize_t = ldt*k;
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&v, &vref, bufsize_v);
	lapacke_gtest_alloc_double_buffer_pair(&t, &tref, bufsize_t);	
	lapacke_gtest_alloc_double_buffer_pair(&c, &cref, bufsize_c);
	
	if ((v==NULL) || (vref==NULL) ||
		(t==NULL) || (tref==NULL) ||
		(c==NULL) || (cref==NULL))
	{
		EXPECT_FALSE( true) << "larfb_double_parameters object: malloc error.";
		larfb_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand( v, vref, bufsize_v);
	lapacke_gtest_init_double_buffer_pair_rand_custom_matrix( t, tref, ldt, k, 'S');

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
larfb_double_parameters :: ~larfb_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" larfb_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   larfb_free();

}
/*  Test fixture class definition */
class dlarfb_test  : public  ::testing::Test {
public:
   larfb_double_parameters  *dlarfb_obj;
   void SetUp();
   void TearDown () { delete dlarfb_obj;}
};

void dlarfb_test::SetUp(){
	 /* LAPACKE dlarfb prototype */
    typedef int (*Fptr_NL_LAPACKE_dlarfb) (int matrix_layout , char side , char trans , char direct , char storev ,\
	lapack_int m , lapack_int n , lapack_int k , const double * v , lapack_int ldv , const double * t , lapack_int ldt , double * c , lapack_int ldc);

    Fptr_NL_LAPACKE_dlarfb dlarfb;

    dlarfb_obj = new larfb_double_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].side,
						   eig_paramslist[idx].trans,
						   eig_paramslist[idx].m,
						   eig_paramslist[idx].n,
                           eig_paramslist[idx].n);
						   

    idx = Circular_Increment_Index(idx);

    dlarfb_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dlarfb_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dlarfb_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dlarfb_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*dlarfb library call */
    dlarfb = (Fptr_NL_LAPACKE_dlarfb)dlsym(dlarfb_obj->hModule, "LAPACKE_dlarfb");
    ASSERT_TRUE(dlarfb != NULL) << "failed to get the Netlib LAPACKE_dlarfb symbol";
  
/*Compute dlarfb's  o/p */
    dlarfb_obj->inforef = dlarfb( dlarfb_obj->matrix_layout, dlarfb_obj->side, dlarfb_obj->trans, dlarfb_obj->direct, dlarfb_obj->storev,\
	dlarfb_obj->m, dlarfb_obj->n,dlarfb_obj->k, dlarfb_obj->vref,dlarfb_obj->ldv, (const double*)dlarfb_obj->tref, dlarfb_obj->ldt,
	dlarfb_obj->cref, dlarfb_obj->ldc);

    /* Compute libflame's Lapacke o/p  */
    dlarfb_obj->info = LAPACKE_dlarfb( dlarfb_obj->matrix_layout, dlarfb_obj->side, dlarfb_obj->trans, dlarfb_obj->direct, dlarfb_obj->storev,\
	dlarfb_obj->m, dlarfb_obj->n, dlarfb_obj->k, dlarfb_obj->v, dlarfb_obj->ldv, (const double*)dlarfb_obj->t, dlarfb_obj->ldt, \
	dlarfb_obj->c, dlarfb_obj->ldc);
	
    if( dlarfb_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame\
        LAPACKE_dlarfb is wrong\n", dlarfb_obj->info );
    }
    if( dlarfb_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dlarfb is wrong\n", 
        dlarfb_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dlarfb_obj->diff =  computeDiff_d( dlarfb_obj->bufsize_c, 
                dlarfb_obj->c, dlarfb_obj->cref );

}

TEST_F(dlarfb_test, dlarfb1) {
    EXPECT_NEAR(0.0, dlarfb_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dlarfb_test, dlarfb2) {
    EXPECT_NEAR(0.0, dlarfb_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dlarfb_test, dlarfb3) {
    EXPECT_NEAR(0.0, dlarfb_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/*Beging scomplex common parameter class definition*/
class larfb_scomplex_parameters{

   public:
	int bufsize_c;
	int bufsize_t;
	int bufsize_v;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int k;
	lapack_int ldc, ldv, ldt;
	lapack_complex_float* v;
	char side, trans, direct, storev;
	/*Output Parameter*/
	lapack_complex_float* t, *c;
	lapack_complex_float *vref, *tref, *cref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      larfb_scomplex_parameters (int matrix_layout, char side, char trans,  lapack_int m, lapack_int n, lapack_int k);
      ~larfb_scomplex_parameters ();

};

/* Constructor definition  scomplex_common_parameters */
larfb_scomplex_parameters:: larfb_scomplex_parameters (int matrix_layout_i, char side_i, char trans_i , lapack_int m_i, lapack_int n_i, lapack_int k_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	k = k_i;
	ldt = k;
	side = side_i;
	trans = trans_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n larfb scomplex: matrix_layout = %d,  m:%d, n: %d \n", matrix_layout, side, m, n);
	#endif

	if (matrix_layout == LAPACK_COL_MAJOR)
	{	/* checking sizes with storev*/	
		storev = 'C';
		if ( side == 'L')
		{	
			ldv = m;
			bufsize_v = ldv*k;
		} else if (side == 'R')
		{
			ldv = n;
			bufsize_v = ldv*k;			
		}
		/* checking sizes with storev*/	
		storev = 'R';
		if ( side == 'L')
		{	
			ldv = k;
			bufsize_v = ldv*m;
		} else if (side == 'R')
		{
			ldv = k;
			bufsize_v = ldv*n;			
		}
		direct = 'F';
		ldc = m;
		bufsize_c = ldc*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{		
		storev = 'R';
		if ( side == 'L')
		{	
			ldv = m;
			bufsize_v = ldv*k;
		} else if (side == 'R')
		{
			ldv = n;
			bufsize_v = ldv*k;			
		}
		/* checking sizes with storev*/	
		storev = 'C';
		if ( side == 'L')
		{	
			ldv = k;
			bufsize_v = ldv*m;
		} else if (side == 'R')
		{
			ldv = k;
			bufsize_v = ldv*n;			
		}
		direct = 'B';
		ldc = n;
		bufsize_c = ldc*m;
	}else
	{
		EXPECT_TRUE(false) << "matrix_layout invalid";
	}
	
	bufsize_t = ldt*k;
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&v, &vref, bufsize_v);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&t, &tref, bufsize_t);	
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&c, &cref, bufsize_c);
	
	if ((v==NULL) || (vref==NULL) ||
		(t==NULL) || (tref==NULL) ||
		(c==NULL) || (cref==NULL))
	{
		EXPECT_FALSE( true) << "larfb_scomplex_parameters object: malloc error.";
		larfb_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( v, vref, bufsize_v);
	lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix( t, tref, ldt, k, 'S');

} /* end of Constructor  */

/* Destructor definition  'scomplex_common_parameters' class  */
larfb_scomplex_parameters :: ~larfb_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" larfb_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   larfb_free();

}
/*  Test fixture class definition */
class clarfb_test  : public  ::testing::Test {
public:
   larfb_scomplex_parameters  *clarfb_obj;
   void SetUp();
   void TearDown () { delete clarfb_obj;}
};

void clarfb_test::SetUp(){
	 /* LAPACKE clarfb prototype */
    typedef int (*Fptr_NL_LAPACKE_clarfb) (int matrix_layout , char side , char trans , char direct , char storev ,\
	lapack_int m , lapack_int n , lapack_int k , const lapack_complex_float * v , lapack_int ldv , const lapack_complex_float * t , lapack_int ldt , lapack_complex_float * c , lapack_int ldc);

    Fptr_NL_LAPACKE_clarfb clarfb;

    clarfb_obj = new larfb_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].side,
						   eig_paramslist[idx].trans,
						   eig_paramslist[idx].m,
						   eig_paramslist[idx].n,
                           eig_paramslist[idx].n);
						   

    idx = Circular_Increment_Index(idx);

    clarfb_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    clarfb_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(clarfb_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(clarfb_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*clarfb library call */
    clarfb = (Fptr_NL_LAPACKE_clarfb)dlsym(clarfb_obj->hModule, "LAPACKE_clarfb");
    ASSERT_TRUE(clarfb != NULL) << "failed to get the Netlib LAPACKE_clarfb symbol";
  
/*Compute clarfb's  o/p */
    clarfb_obj->inforef = clarfb( clarfb_obj->matrix_layout, clarfb_obj->side, clarfb_obj->trans, clarfb_obj->direct, clarfb_obj->storev,\
	clarfb_obj->m, clarfb_obj->n,clarfb_obj->k, clarfb_obj->vref,clarfb_obj->ldv, (const lapack_complex_float*)clarfb_obj->tref, clarfb_obj->ldt,
	clarfb_obj->cref, clarfb_obj->ldc);

    /* Compute libflame's Lapacke o/p  */
    clarfb_obj->info = LAPACKE_clarfb( clarfb_obj->matrix_layout, clarfb_obj->side, clarfb_obj->trans, clarfb_obj->direct, clarfb_obj->storev,\
	clarfb_obj->m, clarfb_obj->n, clarfb_obj->k, clarfb_obj->v, clarfb_obj->ldv, (const lapack_complex_float*)clarfb_obj->t, clarfb_obj->ldt, \
	clarfb_obj->c, clarfb_obj->ldc);
	
    if( clarfb_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame\
        LAPACKE_clarfb is wrong\n", clarfb_obj->info );
    }
    if( clarfb_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_clarfb is wrong\n", 
        clarfb_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    clarfb_obj->diff =  computeDiff_c( clarfb_obj->bufsize_c, 
                clarfb_obj->c, clarfb_obj->cref );

}

TEST_F(clarfb_test, clarfb1) {
    EXPECT_NEAR(0.0, clarfb_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(clarfb_test, clarfb2) {
    EXPECT_NEAR(0.0, clarfb_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(clarfb_test, clarfb3) {
    EXPECT_NEAR(0.0, clarfb_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/*Beging dcomplex common parameter class definition*/
class larfb_dcomplex_parameters{

   public:
	int bufsize_c;
	int bufsize_t;
	int bufsize_v;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int k;
	lapack_int ldc, ldv, ldt;
	lapack_complex_double* v;
	char side, trans, direct, storev;
	/*Output Parameter*/
	lapack_complex_double* t, *c;
	lapack_complex_double *vref, *tref, *cref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      larfb_dcomplex_parameters (int matrix_layout, char side, char trans,  lapack_int m, lapack_int n, lapack_int k);
      ~larfb_dcomplex_parameters ();

};

/* Constructor definition  dcomplex_common_parameters */
larfb_dcomplex_parameters:: larfb_dcomplex_parameters (int matrix_layout_i, char side_i, char trans_i , lapack_int m_i, lapack_int n_i, lapack_int k_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	k = k_i;
	ldt = k;
	side = side_i;
	trans = trans_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n larfb dcomplex: matrix_layout = %d,  m:%d, n: %d \n", matrix_layout, side, m, n);
	#endif

	if (matrix_layout == LAPACK_COL_MAJOR)
	{	/* checking sizes with storev*/	
		storev = 'C';
		if ( side == 'L')
		{	
			ldv = m;
			bufsize_v = ldv*k;
		} else if (side == 'R')
		{
			ldv = n;
			bufsize_v = ldv*k;			
		}
		/* checking sizes with storev*/	
		storev = 'R';
		if ( side == 'L')
		{	
			ldv = k;
			bufsize_v = ldv*m;
		} else if (side == 'R')
		{
			ldv = k;
			bufsize_v = ldv*n;			
		}
		direct = 'F';
		ldc = m;
		bufsize_c = ldc*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{		
		storev = 'R';
		if ( side == 'L')
		{	
			ldv = m;
			bufsize_v = ldv*k;
		} else if (side == 'R')
		{
			ldv = n;
			bufsize_v = ldv*k;			
		}
		/* checking sizes with storev*/	
		storev = 'C';
		if ( side == 'L')
		{	
			ldv = k;
			bufsize_v = ldv*m;
		} else if (side == 'R')
		{
			ldv = k;
			bufsize_v = ldv*n;			
		}
		direct = 'B';
		ldc = n;
		bufsize_c = ldc*m;
	}else
	{
		EXPECT_TRUE(false) << "matrix_layout invalid";
	}
	
	bufsize_t = ldt*k;
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&v, &vref, bufsize_v);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&t, &tref, bufsize_t);	
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&c, &cref, bufsize_c);
	
	if ((v==NULL) || (vref==NULL) ||
		(t==NULL) || (tref==NULL) ||
		(c==NULL) || (cref==NULL))
	{
		EXPECT_FALSE( true) << "larfb_dcomplex_parameters object: malloc error.";
		larfb_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( v, vref, bufsize_v);
	lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix( t, tref, ldt, k, 'S');

} /* end of Constructor  */

/* Destructor definition  'dcomplex_common_parameters' class  */
larfb_dcomplex_parameters :: ~larfb_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" larfb_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   larfb_free();

}
/*  Test fixture class definition */
class zlarfb_test  : public  ::testing::Test {
public:
   larfb_dcomplex_parameters  *zlarfb_obj;
   void SetUp();
   void TearDown () { delete zlarfb_obj;}
};

void zlarfb_test::SetUp(){
	 /* LAPACKE zlarfb prototype */
    typedef int (*Fptr_NL_LAPACKE_zlarfb) (int matrix_layout , char side , char trans , char direct , char storev ,\
	lapack_int m , lapack_int n , lapack_int k , const lapack_complex_double * v , lapack_int ldv , const lapack_complex_double * t , lapack_int ldt , lapack_complex_double * c , lapack_int ldc);

    Fptr_NL_LAPACKE_zlarfb zlarfb;

    zlarfb_obj = new larfb_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].side,
						   eig_paramslist[idx].trans,
						   eig_paramslist[idx].m,
						   eig_paramslist[idx].n,
                           eig_paramslist[idx].n);
						   

    idx = Circular_Increment_Index(idx);

    zlarfb_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zlarfb_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zlarfb_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zlarfb_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*zlarfb library call */
    zlarfb = (Fptr_NL_LAPACKE_zlarfb)dlsym(zlarfb_obj->hModule, "LAPACKE_zlarfb");
    ASSERT_TRUE(zlarfb != NULL) << "failed to get the Netlib LAPACKE_zlarfb symbol";
  
/*Compute zlarfb's  o/p */
    zlarfb_obj->inforef = zlarfb( zlarfb_obj->matrix_layout, zlarfb_obj->side, zlarfb_obj->trans, zlarfb_obj->direct, zlarfb_obj->storev,\
	zlarfb_obj->m, zlarfb_obj->n,zlarfb_obj->k, zlarfb_obj->vref,zlarfb_obj->ldv, (const lapack_complex_double*)zlarfb_obj->tref, zlarfb_obj->ldt,
	zlarfb_obj->cref, zlarfb_obj->ldc);

    /* Compute libflame's Lapacke o/p  */
    zlarfb_obj->info = LAPACKE_zlarfb( zlarfb_obj->matrix_layout, zlarfb_obj->side, zlarfb_obj->trans, zlarfb_obj->direct, zlarfb_obj->storev,\
	zlarfb_obj->m, zlarfb_obj->n, zlarfb_obj->k, zlarfb_obj->v, zlarfb_obj->ldv, (const lapack_complex_double*)zlarfb_obj->t, zlarfb_obj->ldt, \
	zlarfb_obj->c, zlarfb_obj->ldc);
	
    if( zlarfb_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame\
        LAPACKE_zlarfb is wrong\n", zlarfb_obj->info );
    }
    if( zlarfb_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zlarfb is wrong\n", 
        zlarfb_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zlarfb_obj->diff =  computeDiff_z( zlarfb_obj->bufsize_c, 
                zlarfb_obj->c, zlarfb_obj->cref );

}

TEST_F(zlarfb_test, zlarfb1) {
    EXPECT_NEAR(0.0, zlarfb_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlarfb_test, zlarfb2) {
    EXPECT_NEAR(0.0, zlarfb_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlarfb_test, zlarfb3) {
    EXPECT_NEAR(0.0, zlarfb_obj->diff, LAPACKE_GTEST_THRESHOLD);
}