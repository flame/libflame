#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"

#define lascl_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule); \
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
/* Begin float_common_parameters  class definition */
class lascl_float_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	char type;
	lapack_int n,m, lda, count;
	float* A;
	float cto, cfrom;
	lapack_int kl, ku;
	float* Aref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      lascl_float_parameters (int matrix_layout, char type, lapack_int kl, lapack_int ku,  float cfrom ,float cto, lapack_int m, lapack_int n, lapack_int lda, int count);
      ~lascl_float_parameters ();

};	

/* Constructor definition  float_common_parameters */
lascl_float_parameters:: lascl_float_parameters (int matrix_layout_i, char type_i, lapack_int kl_i, lapack_int ku_i,  float cfrom_i ,
float cto_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i, int count_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	kl = kl_i;
	ku = ku_i;
	type = type_i;
	count = count_i;
	cfrom = cfrom_i;
	cto = cto_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n lascl float: matrix_layout = %d,  n: %d \n", matrix_layout,  n);
	#endif
	
	if (count == 0)
	{
		type  = 'G';
	}else if (count == 3)
	{
		type = 'H';
	}	
		bufsize = lda*n;
		

	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&A, &Aref, bufsize);
	if ((A==NULL) || (Aref==NULL)) {	
		EXPECT_FALSE( true) << "lascl_float_parameters object: malloc error.";
		lascl_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand_custom_matrix( A, Aref, m, n, type);
} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
lascl_float_parameters :: ~lascl_float_parameters ()
{ 
	#if LAPACKE_TEST_VERBOSE
	printf(" lascl_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lascl_free();

}
/*  Test fixture class definition */
class slascl_test  : public  ::testing::Test {
public:
   lascl_float_parameters  *slascl_obj;
   void SetUp();
   void TearDown () { delete slascl_obj;}
};

void slascl_test::SetUp(){

 /* LAPACKE slascl prototype */
    typedef int (*Fptr_NL_LAPACKE_slascl) (int matrix_layout, char type, lapack_int kl, lapack_int ku, float cfrom, float cto,\
	lapack_int m, lapack_int n, float * a, lapack_int lda);

    Fptr_NL_LAPACKE_slascl slascl;

    slascl_obj = new lascl_float_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].uplo,
						   lin_solver_paramslist[idx].kl,
						  lin_solver_paramslist[idx].ku,
						   eig_paramslist[idx].k,
						   eig_paramslist[idx].nb,
						   eig_paramslist[idx].n,
						   eig_paramslist[idx].m,
						   eig_paramslist[idx].lda,
						   idx);

    idx = Circular_Increment_Index(idx);

    slascl_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    slascl_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(slascl_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(slascl_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*slascl library call */
    slascl = (Fptr_NL_LAPACKE_slascl)dlsym(slascl_obj->hModule, "LAPACKE_slascl");
    ASSERT_TRUE(slascl != NULL) << "failed to get the Netlib LAPACKE_slascl symbol";

/*Compute slascl's  o/p */
    slascl_obj->inforef = slascl( slascl_obj->matrix_layout, slascl_obj->type, slascl_obj->kl, slascl_obj->ku,\
	slascl_obj->cfrom, slascl_obj->cto, slascl_obj->m, slascl_obj->n, slascl_obj->Aref, slascl_obj->lda );

    /* Compute libflame's Lapacke o/p  */
	
    slascl_obj->info = LAPACKE_slascl( slascl_obj->matrix_layout, slascl_obj->type, slascl_obj->kl, slascl_obj->ku,\
	slascl_obj->cfrom, slascl_obj->cto, slascl_obj->m, slascl_obj->n, slascl_obj->A, slascl_obj->lda );

    if( slascl_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_slascl is wrong\n", slascl_obj->info );
    }
    if( slascl_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_slascl is wrong\n", 
        slascl_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    slascl_obj->diff =  computeDiff_s( slascl_obj->bufsize, 
                slascl_obj->A, slascl_obj->Aref );

}

TEST_F(slascl_test, slascl1) {
    EXPECT_NEAR(0.0, slascl_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(slascl_test, slascl2) {
    EXPECT_NEAR(0.0, slascl_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(slascl_test, slascl3) {
    EXPECT_NEAR(0.0, slascl_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(slascl_test, slascl4) {
    EXPECT_NEAR(0.0, slascl_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin double_common_parameters  class definition */
class lascl_double_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	char type;
	lapack_int n,m, lda, count;
	double* A;
	double cto, cfrom;
	lapack_int kl, ku;
	double* Aref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      lascl_double_parameters (int matrix_layout, char type, lapack_int kl, lapack_int ku,  double cfrom ,double cto, lapack_int m, lapack_int n, lapack_int lda, lapack_int count);
      ~lascl_double_parameters ();

};	

/* Constructor definition  double_common_parameters */
lascl_double_parameters:: lascl_double_parameters (int matrix_layout_i, char type_i, lapack_int kl_i, lapack_int ku_i,  double cfrom_i ,
double cto_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i, lapack_int count_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	kl = kl_i;
	ku = ku_i;
	type = type_i;
	count = count_i;
	cfrom = cfrom_i;
	cto = cto_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n lascl double: matrix_layout = %d,  n: %d \n", matrix_layout,  n);
	#endif
	
	if (count == 0) {
		type  = 'G';
	}else if (count == 3) {
		type = 'H';
	}
		bufsize = lda*n;
		

	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&A, &Aref, bufsize);
	if ((A==NULL) || (Aref==NULL)) {	
		EXPECT_FALSE( true) << "lascl_double_parameters object: malloc error.";
		lascl_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand_custom_matrix( A, Aref, m, n, type);
} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
lascl_double_parameters :: ~lascl_double_parameters ()
{ 
	#if LAPACKE_TEST_VERBOSE
	printf(" lascl_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lascl_free();

}
/*  Test fixture class definition */
class dlascl_test  : public  ::testing::Test {
public:
   lascl_double_parameters  *dlascl_obj;
   void SetUp();
   void TearDown () { delete dlascl_obj;}
};

void dlascl_test::SetUp(){

 /* LAPACKE dlascl prototype */
    typedef int (*Fptr_NL_LAPACKE_dlascl) (int matrix_layout, char type, lapack_int kl, lapack_int ku, double cfrom, double cto,\
	lapack_int m, lapack_int n, double * a, lapack_int lda);

    Fptr_NL_LAPACKE_dlascl dlascl;

    dlascl_obj = new lascl_double_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].uplo,
						  lin_solver_paramslist[idx].kl,
						  lin_solver_paramslist[idx].ku,
						   eig_paramslist[idx].k,
						   eig_paramslist[idx].nb,
						   eig_paramslist[idx].n,
						   eig_paramslist[idx].m,
						   eig_paramslist[idx].lda,
						   idx);

    idx = Circular_Increment_Index(idx);

    dlascl_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dlascl_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dlascl_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dlascl_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*dlascl library call */
    dlascl = (Fptr_NL_LAPACKE_dlascl)dlsym(dlascl_obj->hModule, "LAPACKE_dlascl");
    ASSERT_TRUE(dlascl != NULL) << "failed to get the Netlib LAPACKE_dlascl symbol";

/*Compute dlascl's  o/p */
    dlascl_obj->inforef = dlascl( dlascl_obj->matrix_layout, dlascl_obj->type, dlascl_obj->kl, dlascl_obj->ku,\
	dlascl_obj->cfrom, dlascl_obj->cto, dlascl_obj->m, dlascl_obj->n, dlascl_obj->Aref, dlascl_obj->lda );

    /* Compute libflame's Lapacke o/p  */
	
    dlascl_obj->info = LAPACKE_dlascl( dlascl_obj->matrix_layout, dlascl_obj->type, dlascl_obj->kl, dlascl_obj->ku,\
	dlascl_obj->cfrom, dlascl_obj->cto, dlascl_obj->m, dlascl_obj->n, dlascl_obj->A, dlascl_obj->lda );

    if( dlascl_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dlascl is wrong\n", dlascl_obj->info );
    }
    if( dlascl_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dlascl is wrong\n", 
        dlascl_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dlascl_obj->diff =  computeDiff_d( dlascl_obj->bufsize, 
                dlascl_obj->A, dlascl_obj->Aref );

}

TEST_F(dlascl_test, dlascl1) {
    EXPECT_NEAR(0.0, dlascl_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dlascl_test, dlascl2) {
    EXPECT_NEAR(0.0, dlascl_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dlascl_test, dlascl3) {
    EXPECT_NEAR(0.0, dlascl_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dlascl_test, dlascl4) {
    EXPECT_NEAR(0.0, dlascl_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin scomplex_common_parameters  class definition */
class lascl_scomplex_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	char type;
	lapack_int n,m, lda, count;
	lapack_complex_float* A;
	float cto, cfrom;
	lapack_int kl, ku;
	lapack_complex_float* Aref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      lascl_scomplex_parameters (int matrix_layout, char type, lapack_int kl, lapack_int ku,  float cfrom ,float cto, lapack_int m, lapack_int n, lapack_int lda, lapack_int count);
      ~lascl_scomplex_parameters ();

};	

/* Constructor definition  scomplex_common_parameters */
lascl_scomplex_parameters:: lascl_scomplex_parameters (int matrix_layout_i, char type_i, lapack_int kl_i, lapack_int ku_i,  float cfrom_i ,
float cto_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i, lapack_int count_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	kl = kl_i;
	ku = ku_i;
	type = type_i;
	count = count_i;
	cfrom = cfrom_i;
	cto = cto_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n lascl scomplex: matrix_layout = %d,  n: %d \n", matrix_layout,  n);
	#endif
	
	if (count == 0) {
		type  = 'G';
	} else if (count == 3) {
		type = 'H';
	}
		bufsize = lda*n;
		

	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize);
	if ((A==NULL) || (Aref==NULL)) {	
		EXPECT_FALSE( true) << "lascl_scomplex_parameters object: malloc error.";
		lascl_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix( A, Aref, m, n, type);
} /* end of Constructor  */

/* Destructor definition  'scomplex_common_parameters' class  */
lascl_scomplex_parameters :: ~lascl_scomplex_parameters ()
{ 
	#if LAPACKE_TEST_VERBOSE
	printf(" lascl_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lascl_free();

}
/*  Test fixture class definition */
class clascl_test  : public  ::testing::Test {
public:
   lascl_scomplex_parameters  *clascl_obj;
   void SetUp();
   void TearDown () { delete clascl_obj;}
};

void clascl_test::SetUp(){

 /* LAPACKE clascl prototype */
    typedef int (*Fptr_NL_LAPACKE_clascl) (int matrix_layout, char type, lapack_int kl, lapack_int ku, float cfrom, float cto,\
	lapack_int m, lapack_int n, lapack_complex_float * a, lapack_int lda);

    Fptr_NL_LAPACKE_clascl clascl;

    clascl_obj = new lascl_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].uplo,
						  lin_solver_paramslist[idx].kl,
						  lin_solver_paramslist[idx].ku,
						   eig_paramslist[idx].k,
						   eig_paramslist[idx].nb,
						   eig_paramslist[idx].n,
						   eig_paramslist[idx].m,
						   eig_paramslist[idx].lda,
						   idx);

    idx = Circular_Increment_Index(idx);

    clascl_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    clascl_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(clascl_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(clascl_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*clascl library call */
    clascl = (Fptr_NL_LAPACKE_clascl)dlsym(clascl_obj->hModule, "LAPACKE_clascl");
    ASSERT_TRUE(clascl != NULL) << "failed to get the Netlib LAPACKE_clascl symbol";

/*Compute clascl's  o/p */
    clascl_obj->inforef = clascl( clascl_obj->matrix_layout, clascl_obj->type, clascl_obj->kl, clascl_obj->ku,\
	clascl_obj->cfrom, clascl_obj->cto, clascl_obj->m, clascl_obj->n, clascl_obj->Aref, clascl_obj->lda );

    /* Compute libflame's Lapacke o/p  */
	
    clascl_obj->info = LAPACKE_clascl( clascl_obj->matrix_layout, clascl_obj->type, clascl_obj->kl, clascl_obj->ku,\
	clascl_obj->cfrom, clascl_obj->cto, clascl_obj->m, clascl_obj->n, clascl_obj->A, clascl_obj->lda );

    if( clascl_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_clascl is wrong\n", clascl_obj->info );
    }
    if( clascl_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_clascl is wrong\n", 
        clascl_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    clascl_obj->diff =  computeDiff_c( clascl_obj->bufsize, 
                clascl_obj->A, clascl_obj->Aref );

}

TEST_F(clascl_test, clascl1) {
    EXPECT_NEAR(0.0, clascl_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(clascl_test, clascl2) {
    EXPECT_NEAR(0.0, clascl_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(clascl_test, clascl3) {
    EXPECT_NEAR(0.0, clascl_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(clascl_test, clascl4) {
    EXPECT_NEAR(0.0, clascl_obj->diff, LAPACKE_GTEST_THRESHOLD);
}



/* Begin dcomplex_common_parameters  class definition */
class lascl_dcomplex_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	char type;
	lapack_int n,m, lda, count;
	lapack_complex_double* A;
	double cto, cfrom;
	lapack_int kl, ku;
	lapack_complex_double* Aref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      lascl_dcomplex_parameters (int matrix_layout, char type, lapack_int kl, lapack_int ku,  double cfrom ,double cto, lapack_int m, lapack_int n, lapack_int lda, lapack_int count);
      ~lascl_dcomplex_parameters ();

};	

/* Constructor definition  dcomplex_common_parameters */
lascl_dcomplex_parameters:: lascl_dcomplex_parameters (int matrix_layout_i, char type_i, lapack_int kl_i, lapack_int ku_i,  double cfrom_i ,
double cto_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i,  lapack_int count_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	kl = kl_i;
	ku = ku_i;
	type = type_i;
	count = count_i;
	cfrom = cfrom_i;
	cto = cto_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n lascl dcomplex: matrix_layout = %d,  n: %d \n", matrix_layout,  n);
	#endif
	
	if (count == 0) {
		type  = 'G';
	} else if (count == 3) {
		type = 'H';
	}
		bufsize = lda*n;
		

	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize);
	if ((A==NULL) || (Aref==NULL)) {	
		EXPECT_FALSE( true) << "lascl_dcomplex_parameters object: malloc error.";
		lascl_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix( A, Aref, m, n, type);
} /* end of Constructor  */

/* Destructor definition  'dcomplex_common_parameters' class  */
lascl_dcomplex_parameters :: ~lascl_dcomplex_parameters ()
{ 
	#if LAPACKE_TEST_VERBOSE
	printf(" lascl_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lascl_free();

}
/*  Test fixture class definition */
class zlascl_test  : public  ::testing::Test {
public:
   lascl_dcomplex_parameters  *zlascl_obj;
   void SetUp();
   void TearDown () { delete zlascl_obj;}
};

void zlascl_test::SetUp(){

 /* LAPACKE zlascl prototype */
    typedef int (*Fptr_NL_LAPACKE_zlascl) (int matrix_layout, char type, lapack_int kl, lapack_int ku, double cfrom, double cto,\
	lapack_int m, lapack_int n, lapack_complex_double * a, lapack_int lda);

    Fptr_NL_LAPACKE_zlascl zlascl;

    zlascl_obj = new lascl_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].uplo,
						  lin_solver_paramslist[idx].kl,
						  lin_solver_paramslist[idx].ku,
						   eig_paramslist[idx].k,
						   eig_paramslist[idx].nb,
						   eig_paramslist[idx].n,
						   eig_paramslist[idx].m,
						   eig_paramslist[idx].lda,
						   idx);

    idx = Circular_Increment_Index(idx);

    zlascl_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zlascl_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zlascl_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zlascl_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*zlascl library call */
    zlascl = (Fptr_NL_LAPACKE_zlascl)dlsym(zlascl_obj->hModule, "LAPACKE_zlascl");
    ASSERT_TRUE(zlascl != NULL) << "failed to get the Netlib LAPACKE_zlascl symbol";

/*Compute zlascl's  o/p */
    zlascl_obj->inforef = zlascl( zlascl_obj->matrix_layout, zlascl_obj->type, zlascl_obj->kl, zlascl_obj->ku,\
	zlascl_obj->cfrom, zlascl_obj->cto, zlascl_obj->m, zlascl_obj->n, zlascl_obj->Aref, zlascl_obj->lda );

    /* Compute libflame's Lapacke o/p  */
	
    zlascl_obj->info = LAPACKE_zlascl( zlascl_obj->matrix_layout, zlascl_obj->type, zlascl_obj->kl, zlascl_obj->ku,\
	zlascl_obj->cfrom, zlascl_obj->cto, zlascl_obj->m, zlascl_obj->n, zlascl_obj->A, zlascl_obj->lda );

    if( zlascl_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zlascl is wrong\n", zlascl_obj->info );
    }
    if( zlascl_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zlascl is wrong\n", 
        zlascl_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zlascl_obj->diff =  computeDiff_z( zlascl_obj->bufsize, 
                zlascl_obj->A, zlascl_obj->Aref );

}

TEST_F(zlascl_test, zlascl1) {
    EXPECT_NEAR(0.0, zlascl_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlascl_test, zlascl2) {
    EXPECT_NEAR(0.0, zlascl_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlascl_test, zlascl3) {
    EXPECT_NEAR(0.0, zlascl_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlascl_test, zlascl4) {
    EXPECT_NEAR(0.0, zlascl_obj->diff, LAPACKE_GTEST_THRESHOLD);
}