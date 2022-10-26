#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define stegr_free() \
if (e!=NULL)    free(e); \
if (eref!=NULL) free(eref);\
if (d!=NULL)    free(d); \
if (dref!=NULL)    free(dref); \
if (w!=NULL)  free(w);\
if (wref!=NULL) free(wref); \
if (z!=NULL)  free(z);\
if (zref!=NULL)  free(zref);\
if (isuppz!=NULL)  free(isuppz);\
if (isuppzref!=NULL)  free(isuppzref);\
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin float_common_parameters  class definition */
class stegr_float_parameters{

   public:
	int bufsize_e;
	int bufsize_z;
	int bufsize_w;
	void *hModule, *dModule;
	float diff_e;
	float diff_d;
	float diff_w;
	float diff_z;
	float vl, vu;
	float abstol;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int ldz;
	lapack_int il, iu;
	lapack_int* isuppz;
	float *z;
	float *e, *d;
	char jobz;
	char range;
	/*Output Parameter*/
	float* w;
	lapack_int m;
	lapack_int* isuppzref;
	float *zref;
	float *eref, *dref;
	float* wref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      stegr_float_parameters (int matrix_layout, char jobz, char range,  lapack_int n, lapack_int m, lapack_int ldz);
      ~stegr_float_parameters ();

};

/* Constructor definition  float_common_parameters */
stegr_float_parameters:: stegr_float_parameters (int matrix_layout_i,  char jobz_i, char range_i, lapack_int n_i, lapack_int m_i, lapack_int ldz_i)
{
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	ldz = ldz_i;
	jobz = jobz_i;
	range = range_i;
	iu = il = vl = vu = 0;
	
	if (jobz == 'U')
		jobz = 'N';
	
	if (range == 'N')
	{	range = 'A';
		m = n;
	}else if (range == 'I')
	{
		iu = n/2;
		il = n/5;
		m = iu-il+1;		
	} else if (range == 'V')
	{	vl = n-5;
		vu= n+5;

	}
	
	#if LAPACKE_TEST_VERBOSE
	printf(" \n stegr float: matrix_layout = %d,  jobz:%c, range:%c, uplo:%c, n: %d , m:%d, ldz: %d, iu:%d, Il:%d \n", matrix_layout, jobz, range, uplo, n, m, ldz, iu, il);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{
		ldz = n;
		bufsize_z = ldz*m;
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
    {
		ldz = m;
		bufsize_z = ldz*n;
    }else
		EXPECT_TRUE(false) << "matrix_layout invalid";

	bufsize_e = n;
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&e, &eref, bufsize_e);
	lapacke_gtest_alloc_float_buffer_pair(&d, &dref, bufsize_e);
	lapacke_gtest_alloc_float_buffer_pair(&z, &zref, bufsize_z);
	lapacke_gtest_alloc_float_buffer_pair(&w, &wref, n);
	lapacke_gtest_alloc_int_buffer_pair(&isuppz, &isuppzref, 2*m);
	
	if ((e==NULL) || (eref==NULL) || 
		(w == NULL) || (wref == NULL) ||
		(z == NULL) || (zref == NULL) ||
		(d == NULL) || (dref == NULL) ||
		(isuppz == NULL) || (isuppzref == NULL)){
		EXPECT_FALSE( true) << "stegr_float_parameters object: malloc error.";
		stegr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	//lapacke_gtest_init_float_buffer_pair_rand_custom_matrix(A, Aref, n, n, uplo);
	lapacke_gtest_init_float_buffer_pair_rand(e, eref, bufsize_e);
	lapacke_gtest_init_float_buffer_pair_rand(d, dref, bufsize_e);
} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
stegr_float_parameters :: ~stegr_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" stegr_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   stegr_free();
  

}
/*  Test fixture class definition */
class sstegr_test  : public  ::testing::Test {
public:
   stegr_float_parameters  *sstegr_obj;
   void SetUp();
   void TearDown () { delete sstegr_obj; }
};

void sstegr_test::SetUp(){

    /* LAPACKE sstegr prototype */
    typedef int (*Fptr_NL_LAPACKE_sstegr) (int matrix_layout, char jobz, char range, lapack_int n, float* d, float* e,\
	float vl, float vu, lapack_int il, lapack_int iu, float abstol, lapack_int* m, float* w, float* z, lapack_int ldz, lapack_int* isuppz);

    Fptr_NL_LAPACKE_sstegr sstegr;

    sstegr_obj = new stegr_float_parameters ( eig_paramslist[idx].matrix_layout,
							eig_paramslist[idx].vect_rd,
							eig_paramslist[idx].compz,
							eig_paramslist[idx].m,
							eig_paramslist[idx].n,
							eig_paramslist[idx].m);

    idx = Circular_Increment_Index(idx);

    sstegr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sstegr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sstegr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sstegr_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    sstegr = (Fptr_NL_LAPACKE_sstegr)dlsym(sstegr_obj->hModule, "LAPACKE_sstegr");
    ASSERT_TRUE(sstegr != NULL) << "failed to get the Netlib LAPACKE_sstegr symbol";
    

    sstegr_obj->inforef = sstegr( sstegr_obj->matrix_layout,  sstegr_obj->jobz, sstegr_obj->range, sstegr_obj->n, 
	sstegr_obj->dref, sstegr_obj->eref, sstegr_obj->vl, sstegr_obj->vu, sstegr_obj->il, sstegr_obj->iu,	sstegr_obj->abstol, 
	&sstegr_obj->m, sstegr_obj->w, sstegr_obj->zref, sstegr_obj->ldz, sstegr_obj->isuppzref);

    /* Compute libflame's Lapacke o/p  */
    sstegr_obj->info = LAPACKE_sstegr( sstegr_obj->matrix_layout, sstegr_obj->jobz, sstegr_obj->range, sstegr_obj->n, 
	sstegr_obj->d, sstegr_obj->e, sstegr_obj->vl, sstegr_obj->vu, sstegr_obj->il, sstegr_obj->iu, sstegr_obj->abstol, &sstegr_obj->m, 
	sstegr_obj->w, sstegr_obj->z, sstegr_obj->ldz, sstegr_obj->isuppz);

    if( sstegr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_sstegr is wrong\n", sstegr_obj->info );
    }
    if( sstegr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sstegr is wrong\n", 
        sstegr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    sstegr_obj->diff_e =  computeDiff_s( sstegr_obj->bufsize_e, sstegr_obj->e, sstegr_obj->eref );
	sstegr_obj->diff_d =  computeDiff_s( sstegr_obj->bufsize_e, sstegr_obj->e, sstegr_obj->eref );
	sstegr_obj->diff_z =  computeDiff_s( sstegr_obj->bufsize_z, sstegr_obj->z, sstegr_obj->zref );
	sstegr_obj->diff_w =  computeDiff_s( sstegr_obj->n, sstegr_obj->w, sstegr_obj->wref );
}

TEST_F(sstegr_test, sstegr1) {
    EXPECT_NEAR(0.0, sstegr_obj->diff_e, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, sstegr_obj->diff_d, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, sstegr_obj->diff_z, LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, sstegr_obj->diff_w, LAPACKE_EIG_THRESHOLD);

}

TEST_F(sstegr_test, sstegr2) {
    EXPECT_NEAR(0.0, sstegr_obj->diff_e, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, sstegr_obj->diff_d, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, sstegr_obj->diff_z, LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, sstegr_obj->diff_w, LAPACKE_EIG_THRESHOLD);
}

TEST_F(sstegr_test, sstegr3) {
    EXPECT_NEAR(0.0, sstegr_obj->diff_e, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, sstegr_obj->diff_d, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, sstegr_obj->diff_z, LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, sstegr_obj->diff_w, LAPACKE_EIG_THRESHOLD);
}

TEST_F(sstegr_test, sstegr4) {
   EXPECT_NEAR(0.0, sstegr_obj->diff_e, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, sstegr_obj->diff_d, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, sstegr_obj->diff_z, LAPACKE_EIG_THRESHOLD);
//	EXPECT_NEAR(0.0, sstegr_obj->diff_w, LAPACKE_EIG_THRESHOLD);
	
}

/* Begin double_common_parameters  class definition */
class stegr_double_parameters{

   public:
	int bufsize_e;
	int bufsize_z;
	int bufsize_w;
	void *hModule, *dModule;
	void *lmodule, *bmodule;
	double diff_e;
	double diff_d;
	double diff_w;
	double diff_z;
	double vl, vu;
	double abstol;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int ldz;
	lapack_int il, iu;
	lapack_int* isuppz;
	double *z;
	double *e, *d;
	char jobz;
	char range;
	/*Output Parameter*/
	double* w;
	lapack_int m;
	lapack_int* isuppzref;
	double *zref;
	double *eref, *dref;
	double* wref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      stegr_double_parameters (int matrix_layout, char jobz, char range,  lapack_int n, lapack_int m, lapack_int ldz);
      ~stegr_double_parameters ();

};

/* Constructor definition  double_common_parameters */
stegr_double_parameters:: stegr_double_parameters (int matrix_layout_i,  char jobz_i, char range_i, lapack_int n_i, lapack_int m_i, lapack_int ldz_i)
{
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	ldz = ldz_i;
	jobz = jobz_i;
	range = range_i;
	iu = il = vl = vu = 0;
	
	if (jobz == 'U')
		jobz = 'N';
	
	if (range == 'N')
	{	range = 'A';
		m = n;
	}else if (range == 'I')
	{
		iu = n/2;
		il = n/5;
		m = iu-il+1;		
	} else if (range == 'V')
	{	vl = n-5;
		vu= n+5;

	}
	
	#if LAPACKE_TEST_VERBOSE
	printf(" \n stegr double: matrix_layout = %d,  jobz:%c, range:%c, uplo:%c, n: %d , m:%d, ldz: %d, iu:%d, Il:%d \n", matrix_layout, jobz, range, uplo, n, m, ldz, iu, il);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{
		ldz = n;
		bufsize_z = ldz*m;
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
    {
		ldz = m;
		bufsize_z = ldz*n;
    }else
		EXPECT_TRUE(false) << "matrix_layout invalid";

	bufsize_e = n;
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&e, &eref, bufsize_e);
	lapacke_gtest_alloc_double_buffer_pair(&d, &dref, bufsize_e);
	lapacke_gtest_alloc_double_buffer_pair(&z, &zref, bufsize_z);
	lapacke_gtest_alloc_double_buffer_pair(&w, &wref, n);
	lapacke_gtest_alloc_int_buffer_pair(&isuppz, &isuppzref, 2*m);
	
	if ((e==NULL) || (eref==NULL) || 
		(w == NULL) || (wref == NULL) ||
		(z == NULL) || (zref == NULL) ||
		(d == NULL) || (dref == NULL) ||
		(isuppz == NULL) || (isuppzref == NULL)){
		EXPECT_FALSE( true) << "stegr_double_parameters object: malloc error.";
		stegr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	//lapacke_gtest_init_double_buffer_pair_rand_custom_matrix(A, Aref, n, n, uplo);
	lapacke_gtest_init_double_buffer_pair_rand(e, eref, bufsize_e);
	lapacke_gtest_init_double_buffer_pair_rand(d, dref, bufsize_e);
} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
stegr_double_parameters :: ~stegr_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" stegr_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   stegr_free();
  

}
/*  Test fixture class definition */
class dstegr_test  : public  ::testing::Test {
public:
   stegr_double_parameters  *dstegr_obj;
   void SetUp();
   void TearDown () { delete dstegr_obj; }
};

void dstegr_test::SetUp(){

    /* LAPACKE dstegr prototype */
    typedef int (*Fptr_NL_LAPACKE_dstegr) (int matrix_layout, char jobz, char range, lapack_int n, double* d, double* e,\
	double vl, double vu, lapack_int il, lapack_int iu, double abstol, lapack_int* m, double* w, double* z, lapack_int ldz, lapack_int* isuppz);

    Fptr_NL_LAPACKE_dstegr dstegr;

    dstegr_obj = new stegr_double_parameters ( eig_paramslist[idx].matrix_layout,
							eig_paramslist[idx].vect_rd,
							eig_paramslist[idx].compz,
							eig_paramslist[idx].m,
							eig_paramslist[idx].n,
							eig_paramslist[idx].m);

    idx = Circular_Increment_Index(idx);

    dstegr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dstegr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dstegr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dstegr_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    dstegr = (Fptr_NL_LAPACKE_dstegr)dlsym(dstegr_obj->hModule, "LAPACKE_dstegr");
    ASSERT_TRUE(dstegr != NULL) << "failed to get the Netlib LAPACKE_dstegr symbol";
    

    dstegr_obj->inforef = dstegr( dstegr_obj->matrix_layout,  dstegr_obj->jobz, dstegr_obj->range, dstegr_obj->n, 
	dstegr_obj->dref, dstegr_obj->eref, dstegr_obj->vl, dstegr_obj->vu, dstegr_obj->il, dstegr_obj->iu,	dstegr_obj->abstol, 
	&dstegr_obj->m, dstegr_obj->w, dstegr_obj->zref, dstegr_obj->ldz, dstegr_obj->isuppzref);

    /* Compute libflame's Lapacke o/p  */
    dstegr_obj->info = LAPACKE_dstegr( dstegr_obj->matrix_layout, dstegr_obj->jobz, dstegr_obj->range, dstegr_obj->n, 
	dstegr_obj->d, dstegr_obj->e, dstegr_obj->vl, dstegr_obj->vu, dstegr_obj->il, dstegr_obj->iu, dstegr_obj->abstol, &dstegr_obj->m, 
	dstegr_obj->w, dstegr_obj->z, dstegr_obj->ldz, dstegr_obj->isuppz);

    if( dstegr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dstegr is wrong\n", dstegr_obj->info );
    }
    if( dstegr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dstegr is wrong\n", 
        dstegr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dstegr_obj->diff_e =  computeDiff_d( dstegr_obj->bufsize_e, dstegr_obj->e, dstegr_obj->eref );
	dstegr_obj->diff_d =  computeDiff_d( dstegr_obj->bufsize_e, dstegr_obj->e, dstegr_obj->eref );
	dstegr_obj->diff_z =  computeDiff_d( dstegr_obj->bufsize_z, dstegr_obj->z, dstegr_obj->zref );
	dstegr_obj->diff_w =  computeDiff_d( dstegr_obj->n, dstegr_obj->w, dstegr_obj->wref );
}

TEST_F(dstegr_test, dstegr1) {
    EXPECT_NEAR(0.0, dstegr_obj->diff_e, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, dstegr_obj->diff_d, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, dstegr_obj->diff_z, LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, dstegr_obj->diff_w, LAPACKE_EIG_THRESHOLD);

}

TEST_F(dstegr_test, dstegr2) {
    EXPECT_NEAR(0.0, dstegr_obj->diff_e, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, dstegr_obj->diff_d, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, dstegr_obj->diff_z, LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, dstegr_obj->diff_w, LAPACKE_EIG_THRESHOLD);
}

TEST_F(dstegr_test, dstegr3) {
    EXPECT_NEAR(0.0, dstegr_obj->diff_e, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, dstegr_obj->diff_d, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, dstegr_obj->diff_z, LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, dstegr_obj->diff_w, LAPACKE_EIG_THRESHOLD);
}

TEST_F(dstegr_test, dstegr4) {
   EXPECT_NEAR(0.0, dstegr_obj->diff_e, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, dstegr_obj->diff_d, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, dstegr_obj->diff_z, LAPACKE_EIG_THRESHOLD);
//	EXPECT_NEAR(0.0, dstegr_obj->diff_w, LAPACKE_EIG_THRESHOLD);
	
}

/* Begin scomplex_common_parameters  class definition */
class stegr_scomplex_parameters{

   public:
	int bufsize_e;
	int bufsize_z;
	int bufsize_w;
	void *hModule, *dModule;
	void *lmodule, *bmodule;
	float diff_e;
	float diff_d;
	float diff_w;
	float diff_z;
	float vl, vu;
	float abstol;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int ldz;
	lapack_int il, iu;
	lapack_int* isuppz;
	lapack_complex_float *z;
	float *e, *d;
	char jobz;
	char range;
	/*Output Parameter*/
	float* w;
	lapack_int m;
	lapack_int* isuppzref;
	lapack_complex_float *zref;
	float *eref, *dref;
	float* wref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      stegr_scomplex_parameters (int matrix_layout, char jobz, char range,  lapack_int n, lapack_int m, lapack_int ldz);
      ~stegr_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
stegr_scomplex_parameters:: stegr_scomplex_parameters (int matrix_layout_i,  char jobz_i, char range_i, lapack_int n_i, lapack_int m_i, lapack_int ldz_i)
{
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	ldz = ldz_i;
	jobz = jobz_i;
	range = range_i;
	iu = il = vl = vu = 0;
	
	if (jobz == 'U')
		jobz = 'N';
	
	if (range == 'N')
	{	range = 'A';
		m = n;
	}else if (range == 'I')
	{
		iu = n/2;
		il = n/5;
		m = iu-il+1;		
	} else if (range == 'V')
	{	vl = n-5;
		vu= n+5;

	}
	
	#if LAPACKE_TEST_VERBOSE
	printf(" \n stegr scomplex: matrix_layout = %d,  jobz:%c, range:%c, uplo:%c, n: %d , m:%d, ldz: %d, iu:%d, Il:%d \n", matrix_layout, jobz, range, uplo, n, m, ldz, iu, il);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{
		ldz = n;
		bufsize_z = ldz*m;
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
    {
		ldz = m;
		bufsize_z = ldz*n;
    }else
		EXPECT_TRUE(false) << "matrix_layout invalid";

	bufsize_e = n;
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&e, &eref, bufsize_e);
	lapacke_gtest_alloc_float_buffer_pair(&d, &dref, bufsize_e);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&z, &zref, bufsize_z);
	lapacke_gtest_alloc_float_buffer_pair(&w, &wref, n);
	lapacke_gtest_alloc_int_buffer_pair(&isuppz, &isuppzref, 2*m);
	
	if ((e==NULL) || (eref==NULL) || 
		(w == NULL) || (wref == NULL) ||
		(z == NULL) || (zref == NULL) ||
		(d == NULL) || (dref == NULL) ||
		(isuppz == NULL) || (isuppzref == NULL)){
		EXPECT_FALSE( true) << "stegr_float_parameters object: malloc error.";
		stegr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	//lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix(A, Aref, n, n, uplo);
	lapacke_gtest_init_float_buffer_pair_rand(e, eref, bufsize_e);
	lapacke_gtest_init_float_buffer_pair_rand(d, dref, bufsize_e);
} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
stegr_scomplex_parameters :: ~stegr_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" stegr_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   stegr_free();
  

}
/*  Test fixture class definition */
class cstegr_test  : public  ::testing::Test {
public:
   stegr_scomplex_parameters  *cstegr_obj;
   void SetUp();
   void TearDown () { delete cstegr_obj; }
};

void cstegr_test::SetUp(){

    /* LAPACKE cstegr prototype */
    typedef int (*Fptr_NL_LAPACKE_cstegr) (int matrix_layout, char jobz, char range, lapack_int n, float* d, float* e,\
	float vl, float vu, lapack_int il, lapack_int iu, float abstol, lapack_int* m, float* w, lapack_complex_float* z, lapack_int ldz, lapack_int* isuppz);

    Fptr_NL_LAPACKE_cstegr cstegr;

    cstegr_obj = new stegr_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
							eig_paramslist[idx].vect_rd,
							eig_paramslist[idx].compz,
							eig_paramslist[idx].m,
							eig_paramslist[idx].n,
							eig_paramslist[idx].m);

    idx = Circular_Increment_Index(idx);

    cstegr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cstegr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cstegr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cstegr_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    cstegr = (Fptr_NL_LAPACKE_cstegr)dlsym(cstegr_obj->hModule, "LAPACKE_cstegr");
    ASSERT_TRUE(cstegr != NULL) << "failed to get the Netlib LAPACKE_cstegr symbol";
    

    cstegr_obj->inforef = cstegr( cstegr_obj->matrix_layout,  cstegr_obj->jobz, cstegr_obj->range, cstegr_obj->n, 
	cstegr_obj->dref, cstegr_obj->eref, cstegr_obj->vl, cstegr_obj->vu, cstegr_obj->il, cstegr_obj->iu,	cstegr_obj->abstol, 
	&cstegr_obj->m, cstegr_obj->w, cstegr_obj->zref, cstegr_obj->ldz, cstegr_obj->isuppzref);

    /* Compute libflame's Lapacke o/p  */
    cstegr_obj->info = LAPACKE_cstegr( cstegr_obj->matrix_layout, cstegr_obj->jobz, cstegr_obj->range, cstegr_obj->n, 
	cstegr_obj->d, cstegr_obj->e, cstegr_obj->vl, cstegr_obj->vu, cstegr_obj->il, cstegr_obj->iu, cstegr_obj->abstol, &cstegr_obj->m, 
	cstegr_obj->w, cstegr_obj->z, cstegr_obj->ldz, cstegr_obj->isuppz);

    if( cstegr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cstegr is wrong\n", cstegr_obj->info );
    }
    if( cstegr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cstegr is wrong\n", 
        cstegr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cstegr_obj->diff_e =  computeDiff_s( cstegr_obj->bufsize_e, cstegr_obj->e, cstegr_obj->eref );
	cstegr_obj->diff_d =  computeDiff_s( cstegr_obj->bufsize_e, cstegr_obj->e, cstegr_obj->eref );
	cstegr_obj->diff_z =  computeDiff_c( cstegr_obj->bufsize_z, cstegr_obj->z, cstegr_obj->zref );
	cstegr_obj->diff_w =  computeDiff_s( cstegr_obj->n, cstegr_obj->w, cstegr_obj->wref );
}

TEST_F(cstegr_test, cstegr1) {
    EXPECT_NEAR(0.0, cstegr_obj->diff_e, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, cstegr_obj->diff_d, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, cstegr_obj->diff_z, LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, cstegr_obj->diff_w, LAPACKE_EIG_THRESHOLD);

}

TEST_F(cstegr_test, cstegr2) {
    EXPECT_NEAR(0.0, cstegr_obj->diff_e, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, cstegr_obj->diff_d, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, cstegr_obj->diff_z, LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, cstegr_obj->diff_w, LAPACKE_EIG_THRESHOLD);
}

TEST_F(cstegr_test, cstegr3) {
    EXPECT_NEAR(0.0, cstegr_obj->diff_e, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, cstegr_obj->diff_d, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, cstegr_obj->diff_z, LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, cstegr_obj->diff_w, LAPACKE_EIG_THRESHOLD);
}

TEST_F(cstegr_test, cstegr4) {
   EXPECT_NEAR(0.0, cstegr_obj->diff_e, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, cstegr_obj->diff_d, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, cstegr_obj->diff_z, LAPACKE_EIG_THRESHOLD);
//	EXPECT_NEAR(0.0, cstegr_obj->diff_w, LAPACKE_EIG_THRESHOLD);
	
}
/* Begin dcomplex_common_parameters  class definition */
class stegr_dcomplex_parameters{

   public:
	int bufsize_e;
	int bufsize_z;
	int bufsize_w;
	void *hModule, *dModule;
	void *lmodule, *bmodule;
	double diff_e;
	double diff_d;
	double diff_w;
	double diff_z;
	double vl, vu;
	double abstol;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int ldz;
	lapack_int il, iu;
	lapack_int* isuppz;
	lapack_complex_double *z;
	double *e, *d;
	char jobz;
	char range;
	/*Output Parameter*/
	double* w;
	lapack_int m;
	lapack_int* isuppzref;
	lapack_complex_double *zref;
	double *eref, *dref;
	double* wref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      stegr_dcomplex_parameters (int matrix_layout, char jobz, char range,  lapack_int n, lapack_int m, lapack_int ldz);
      ~stegr_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
stegr_dcomplex_parameters:: stegr_dcomplex_parameters (int matrix_layout_i,  char jobz_i, char range_i, lapack_int n_i, lapack_int m_i, lapack_int ldz_i)
{
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	ldz = ldz_i;
	jobz = jobz_i;
	range = range_i;
	iu = il = vl = vu = 0;
	
	if (jobz == 'U')
		jobz = 'N';
	
	if (range == 'N')
	{	range = 'A';
		m = n;
	}else if (range == 'I')
	{
		iu = n/2;
		il = n/5;
		m = iu-il+1;		
	} else if (range == 'V')
	{	vl = n-5;
		vu= n+5;

	}
	
	#if LAPACKE_TEST_VERBOSE
	printf(" \n stegr dcomplex: matrix_layout = %d,  jobz:%c, range:%c, uplo:%c, n: %d , m:%d, ldz: %d, iu:%d, Il:%d \n", matrix_layout, jobz, range, uplo, n, m, ldz, iu, il);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{
		ldz = n;
		bufsize_z = ldz*m;
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
    {
		ldz = m;
		bufsize_z = ldz*n;
    }else
		EXPECT_TRUE(false) << "matrix_layout invalid";

	bufsize_e = n;
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&e, &eref, bufsize_e);
	lapacke_gtest_alloc_double_buffer_pair(&d, &dref, bufsize_e);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&z, &zref, bufsize_z);
	lapacke_gtest_alloc_double_buffer_pair(&w, &wref, n);
	lapacke_gtest_alloc_int_buffer_pair(&isuppz, &isuppzref, 2*m);
	
	if ((e==NULL) || (eref==NULL) || 
		(w == NULL) || (wref == NULL) ||
		(z == NULL) || (zref == NULL) ||
		(d == NULL) || (dref == NULL) ||
		(isuppz == NULL) || (isuppzref == NULL)){
		EXPECT_FALSE( true) << "stegr_double_parameters object: malloc error.";
		stegr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	//lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix(A, Aref, n, n, uplo);
	lapacke_gtest_init_double_buffer_pair_rand(e, eref, bufsize_e);
	lapacke_gtest_init_double_buffer_pair_rand(d, dref, bufsize_e);
} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
stegr_dcomplex_parameters :: ~stegr_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" stegr_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   stegr_free();
  

}
/*  Test fixture class definition */
class zstegr_test  : public  ::testing::Test {
public:
   stegr_dcomplex_parameters  *zstegr_obj;
   void SetUp();
   void TearDown () { delete zstegr_obj; }
};

void zstegr_test::SetUp(){

    /* LAPACKE zstegr prototype */
    typedef int (*Fptr_NL_LAPACKE_zstegr) (int matrix_layout, char jobz, char range, lapack_int n, double* d, double* e,\
	double vl, double vu, lapack_int il, lapack_int iu, double abstol, lapack_int* m, double* w, lapack_complex_double* z, lapack_int ldz, lapack_int* isuppz);

    Fptr_NL_LAPACKE_zstegr zstegr;

    zstegr_obj = new stegr_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
							eig_paramslist[idx].vect_rd,
							eig_paramslist[idx].compz,
							eig_paramslist[idx].m,
							eig_paramslist[idx].n,
							eig_paramslist[idx].m);

    idx = Circular_Increment_Index(idx);

    zstegr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zstegr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zstegr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zstegr_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zstegr = (Fptr_NL_LAPACKE_zstegr)dlsym(zstegr_obj->hModule, "LAPACKE_zstegr");
    ASSERT_TRUE(zstegr != NULL) << "failed to get the Netlib LAPACKE_zstegr symbol";
    

    zstegr_obj->inforef = zstegr( zstegr_obj->matrix_layout,  zstegr_obj->jobz, zstegr_obj->range, zstegr_obj->n, 
	zstegr_obj->dref, zstegr_obj->eref, zstegr_obj->vl, zstegr_obj->vu, zstegr_obj->il, zstegr_obj->iu,	zstegr_obj->abstol, 
	&zstegr_obj->m, zstegr_obj->w, zstegr_obj->zref, zstegr_obj->ldz, zstegr_obj->isuppzref);

    /* Compute libflame's Lapacke o/p  */
    zstegr_obj->info = LAPACKE_zstegr( zstegr_obj->matrix_layout, zstegr_obj->jobz, zstegr_obj->range, zstegr_obj->n, 
	zstegr_obj->d, zstegr_obj->e, zstegr_obj->vl, zstegr_obj->vu, zstegr_obj->il, zstegr_obj->iu, zstegr_obj->abstol, &zstegr_obj->m, 
	zstegr_obj->w, zstegr_obj->z, zstegr_obj->ldz, zstegr_obj->isuppz);

    if( zstegr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zstegr is wrong\n", zstegr_obj->info );
    }
    if( zstegr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zstegr is wrong\n", 
        zstegr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zstegr_obj->diff_e =  computeDiff_d( zstegr_obj->bufsize_e, zstegr_obj->e, zstegr_obj->eref );
	zstegr_obj->diff_d =  computeDiff_d( zstegr_obj->bufsize_e, zstegr_obj->e, zstegr_obj->eref );
	zstegr_obj->diff_z =  computeDiff_z( zstegr_obj->bufsize_z, zstegr_obj->z, zstegr_obj->zref );
	zstegr_obj->diff_w =  computeDiff_d( zstegr_obj->n, zstegr_obj->w, zstegr_obj->wref );
}

TEST_F(zstegr_test, zstegr1) {
    EXPECT_NEAR(0.0, zstegr_obj->diff_e, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zstegr_obj->diff_d, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zstegr_obj->diff_z, LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, zstegr_obj->diff_w, LAPACKE_EIG_THRESHOLD);

}

TEST_F(zstegr_test, zstegr2) {
    EXPECT_NEAR(0.0, zstegr_obj->diff_e, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zstegr_obj->diff_d, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zstegr_obj->diff_z, LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, zstegr_obj->diff_w, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zstegr_test, zstegr3) {
    EXPECT_NEAR(0.0, zstegr_obj->diff_e, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zstegr_obj->diff_d, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zstegr_obj->diff_z, LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, zstegr_obj->diff_w, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zstegr_test, zstegr4) {
   EXPECT_NEAR(0.0, zstegr_obj->diff_e, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zstegr_obj->diff_d, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zstegr_obj->diff_z, LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, zstegr_obj->diff_w, LAPACKE_EIG_THRESHOLD);
	
}