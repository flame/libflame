#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"

#define stemr_free() \
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
class stemr_float_parameters{

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
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int ldz,nzc;
	lapack_int il, iu;
	lapack_int* isuppz;
	lapack_logical tryrac;
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
      stemr_float_parameters (int matrix_layout, char jobz, char range,  lapack_int n, lapack_int m, lapack_int ldz);
      ~stemr_float_parameters ();

};

/* Constructor definition  float_common_parameters */
stemr_float_parameters:: stemr_float_parameters (int matrix_layout_i,  char jobz_i, char range_i, lapack_int n_i, lapack_int m_i, lapack_int ldz_i)
{
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	ldz = ldz_i;
	jobz = jobz_i;
	range = range_i;
	tryrac = 0;
	iu = il = vl = vu = 0;
	
	if (jobz == 'U')
		jobz = 'N';
	
	if (range == 'N')
	{	range = 'A';
		m = n;
		nzc = n;
	}else if (range == 'I')
	{
		iu = n/2;
		il = n/5;
		m = iu-il+1;
		nzc = iu-il+1;
	} else if (range == 'V')
	{	vl = n-5;
		vu= n+5;
		nzc = vu;
	}
	
	#if LAPACKE_TEST_VERBOSE
	printf(" \n stemr float: matrix_layout = %d,  jobz:%c, range:%c, uplo:%c, n: %d , m:%d, ldz: %d, iu:%d, Il:%d \n", matrix_layout, jobz, range, uplo, n, m, ldz, iu, il);
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
		EXPECT_FALSE( true) << "stemr_float_parameters object: malloc error.";
		stemr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	//lapacke_gtest_init_float_buffer_pair_rand_custom_matrix(A, Aref, n, n, uplo);
	lapacke_gtest_init_float_buffer_pair_rand(e, eref, bufsize_e);
	lapacke_gtest_init_float_buffer_pair_rand(d, dref, bufsize_e);
} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
stemr_float_parameters :: ~stemr_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" stemr_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   stemr_free();
  

}
/*  Test fixture class definition */
class sstemr_test  : public  ::testing::Test {
public:
   stemr_float_parameters  *sstemr_obj;
   void SetUp();
   void TearDown () { delete sstemr_obj; }
};

void sstemr_test::SetUp(){

    /* LAPACKE sstemr prototype */
    typedef int (*Fptr_NL_LAPACKE_sstemr) (int matrix_layout, char jobz, char range, lapack_int n,\
	 float* d, float* e, float vl, float vu, lapack_int il, lapack_int iu, lapack_int* m, float* w,\
	float* z, lapack_int ldz, lapack_int nzc, lapack_int* isuppz, lapack_logical* tryrac);

    Fptr_NL_LAPACKE_sstemr sstemr;

    sstemr_obj = new stemr_float_parameters ( eig_paramslist[idx].matrix_layout,
							eig_paramslist[idx].vect_rd,
							eig_paramslist[idx].compz,
							eig_paramslist[idx].m,
							eig_paramslist[idx].n,
							eig_paramslist[idx].m);

    idx = Circular_Increment_Index(idx);

    sstemr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sstemr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sstemr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sstemr_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    sstemr = (Fptr_NL_LAPACKE_sstemr)dlsym(sstemr_obj->hModule, "LAPACKE_sstemr");
    ASSERT_TRUE(sstemr != NULL) << "failed to get the Netlib LAPACKE_sstemr symbol";
    

    sstemr_obj->inforef = sstemr( sstemr_obj->matrix_layout,  sstemr_obj->jobz, sstemr_obj->range, sstemr_obj->n, 
	sstemr_obj->dref, sstemr_obj->eref, sstemr_obj->vl, sstemr_obj->vu, sstemr_obj->il, sstemr_obj->iu, 
	&sstemr_obj->m, sstemr_obj->w, sstemr_obj->zref, sstemr_obj->ldz, sstemr_obj->nzc, sstemr_obj->isuppzref, &sstemr_obj->tryrac);

    /* Compute libflame's Lapacke o/p  */
    sstemr_obj->info = LAPACKE_sstemr( sstemr_obj->matrix_layout, sstemr_obj->jobz, sstemr_obj->range, sstemr_obj->n, 
	sstemr_obj->d, sstemr_obj->e, sstemr_obj->vl, sstemr_obj->vu, sstemr_obj->il, sstemr_obj->iu, &sstemr_obj->m, 
	sstemr_obj->w, sstemr_obj->z, sstemr_obj->ldz, sstemr_obj->nzc, sstemr_obj->isuppz, &sstemr_obj->tryrac);

    if( sstemr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_sstemr is wrong\n", sstemr_obj->info );
    }
    if( sstemr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sstemr is wrong\n", 
        sstemr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    sstemr_obj->diff_e =  computeDiff_s( sstemr_obj->bufsize_e, sstemr_obj->e, sstemr_obj->eref );
	sstemr_obj->diff_d =  computeDiff_s( sstemr_obj->bufsize_e, sstemr_obj->e, sstemr_obj->eref );
	sstemr_obj->diff_z =  computeDiff_s( sstemr_obj->bufsize_z, sstemr_obj->z, sstemr_obj->zref );
	sstemr_obj->diff_w =  computeDiff_s( sstemr_obj->n, sstemr_obj->w, sstemr_obj->wref );
}

TEST_F(sstemr_test, sstemr1) {
    EXPECT_NEAR(0.0, sstemr_obj->diff_e, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, sstemr_obj->diff_d, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, sstemr_obj->diff_z, LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, sstemr_obj->diff_w, LAPACKE_EIG_THRESHOLD);

}

TEST_F(sstemr_test, sstemr2) {
    EXPECT_NEAR(0.0, sstemr_obj->diff_e, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, sstemr_obj->diff_d, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, sstemr_obj->diff_z, LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, sstemr_obj->diff_w, LAPACKE_EIG_THRESHOLD);
}

TEST_F(sstemr_test, sstemr3) {
    EXPECT_NEAR(0.0, sstemr_obj->diff_e, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, sstemr_obj->diff_d, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, sstemr_obj->diff_z, LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, sstemr_obj->diff_w, LAPACKE_EIG_THRESHOLD);
}

TEST_F(sstemr_test, sstemr4) {
   EXPECT_NEAR(0.0, sstemr_obj->diff_e, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, sstemr_obj->diff_d, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, sstemr_obj->diff_z, LAPACKE_EIG_THRESHOLD);
//	EXPECT_NEAR(0.0, sstemr_obj->diff_w, LAPACKE_EIG_THRESHOLD);
	
}

/* Begin double_common_parameters  class definition */
class stemr_double_parameters{

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
	lapack_int nzc;
	/*Output Parameter*/
	double* w;
	lapack_int m;
	lapack_int* isuppzref;
	double *zref;
	double *eref, *dref;
	double* wref;
	lapack_logical tryrac;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      stemr_double_parameters (int matrix_layout, char jobz, char range,  lapack_int n, lapack_int m, lapack_int ldz);
      ~stemr_double_parameters ();

};

/* Constructor definition  double_common_parameters */
stemr_double_parameters:: stemr_double_parameters (int matrix_layout_i,  char jobz_i, char range_i, lapack_int n_i, lapack_int m_i, lapack_int ldz_i)
{
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	ldz = ldz_i;
	jobz = jobz_i;
	range = range_i;
	tryrac = 0;
	iu = il = vl = vu = 0;
	
	if (jobz == 'U')
		jobz = 'N';
	
	if (range == 'N')
	{	range = 'A';
		m = n;
		nzc = n;		
	}else if (range == 'I')
	{
		iu = n/2;
		il = n/5;
		m = iu-il+1;
		nzc = iu-il+1;
	} else if (range == 'V')
	{	vl = n-5;
		vu= n+5;
		nzc = vu;
	}
	
	#if LAPACKE_TEST_VERBOSE
	printf(" \n stemr double: matrix_layout = %d,  jobz:%c, range:%c, uplo:%c, n: %d , m:%d, ldz: %d, iu:%d, Il:%d \n", matrix_layout, jobz, range, uplo, n, m, ldz, iu, il);
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
		EXPECT_FALSE( true) << "stemr_double_parameters object: malloc error.";
		stemr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	//lapacke_gtest_init_double_buffer_pair_rand_custom_matrix(A, Aref, n, n, uplo);
	lapacke_gtest_init_double_buffer_pair_rand(e, eref, bufsize_e);
	lapacke_gtest_init_double_buffer_pair_rand(d, dref, bufsize_e);
} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
stemr_double_parameters :: ~stemr_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" stemr_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   stemr_free();
  

}
/*  Test fixture class definition */
class dstemr_test  : public  ::testing::Test {
public:
   stemr_double_parameters  *dstemr_obj;
   void SetUp();
   void TearDown () { delete dstemr_obj; }
};

void dstemr_test::SetUp(){

    /* LAPACKE dstemr prototype */
    typedef int (*Fptr_NL_LAPACKE_dstemr) (int matrix_layout, char jobz, char range, lapack_int n,  double* d, double* e,\
	double vl, double vu, lapack_int il, lapack_int iu, lapack_int* m, double* w, double* z, lapack_int ldz, lapack_int nzc, lapack_int* isuppz, lapack_logical* tryrac);

    Fptr_NL_LAPACKE_dstemr dstemr;

    dstemr_obj = new stemr_double_parameters ( eig_paramslist[idx].matrix_layout,
							eig_paramslist[idx].vect_rd,
							eig_paramslist[idx].compz,
							eig_paramslist[idx].m,
							eig_paramslist[idx].n,
							eig_paramslist[idx].m);

    idx = Circular_Increment_Index(idx);

    dstemr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dstemr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dstemr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dstemr_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    dstemr = (Fptr_NL_LAPACKE_dstemr)dlsym(dstemr_obj->hModule, "LAPACKE_dstemr");
    ASSERT_TRUE(dstemr != NULL) << "failed to get the Netlib LAPACKE_dstemr symbol";
    

    dstemr_obj->inforef = dstemr( dstemr_obj->matrix_layout,  dstemr_obj->jobz, dstemr_obj->range, dstemr_obj->n, 
	dstemr_obj->dref, dstemr_obj->eref, dstemr_obj->vl, dstemr_obj->vu, dstemr_obj->il, dstemr_obj->iu,	 
	&dstemr_obj->m, dstemr_obj->w, dstemr_obj->zref, dstemr_obj->ldz, dstemr_obj->nzc, dstemr_obj->isuppzref, &dstemr_obj->tryrac);

    /* Compute libflame's Lapacke o/p  */
    dstemr_obj->info = LAPACKE_dstemr( dstemr_obj->matrix_layout, dstemr_obj->jobz, dstemr_obj->range, dstemr_obj->n, 
	dstemr_obj->d, dstemr_obj->e, dstemr_obj->vl, dstemr_obj->vu, dstemr_obj->il, dstemr_obj->iu, &dstemr_obj->m, 
	dstemr_obj->w, dstemr_obj->z, dstemr_obj->ldz, dstemr_obj->nzc, dstemr_obj->isuppz, &dstemr_obj->tryrac);

    if( dstemr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dstemr is wrong\n", dstemr_obj->info );
    }
    if( dstemr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dstemr is wrong\n", 
        dstemr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dstemr_obj->diff_e =  computeDiff_d( dstemr_obj->bufsize_e, dstemr_obj->e, dstemr_obj->eref );
	dstemr_obj->diff_d =  computeDiff_d( dstemr_obj->bufsize_e, dstemr_obj->e, dstemr_obj->eref );
	dstemr_obj->diff_z =  computeDiff_d( dstemr_obj->bufsize_z, dstemr_obj->z, dstemr_obj->zref );
	dstemr_obj->diff_w =  computeDiff_d( dstemr_obj->n, dstemr_obj->w, dstemr_obj->wref );
}

TEST_F(dstemr_test, dstemr1) {
    EXPECT_NEAR(0.0, dstemr_obj->diff_e, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, dstemr_obj->diff_d, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, dstemr_obj->diff_z, LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, dstemr_obj->diff_w, LAPACKE_EIG_THRESHOLD);

}

TEST_F(dstemr_test, dstemr2) {
    EXPECT_NEAR(0.0, dstemr_obj->diff_e, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, dstemr_obj->diff_d, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, dstemr_obj->diff_z, LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, dstemr_obj->diff_w, LAPACKE_EIG_THRESHOLD);
}

TEST_F(dstemr_test, dstemr3) {
    EXPECT_NEAR(0.0, dstemr_obj->diff_e, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, dstemr_obj->diff_d, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, dstemr_obj->diff_z, LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, dstemr_obj->diff_w, LAPACKE_EIG_THRESHOLD);
}

TEST_F(dstemr_test, dstemr4) {
   EXPECT_NEAR(0.0, dstemr_obj->diff_e, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, dstemr_obj->diff_d, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, dstemr_obj->diff_z, LAPACKE_EIG_THRESHOLD);
//	EXPECT_NEAR(0.0, dstemr_obj->diff_w, LAPACKE_EIG_THRESHOLD);
	
}

/* Begin scomplex_common_parameters  class definition */
class stemr_scomplex_parameters{

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
	lapack_int nzc;
	/*Output Parameter*/
	float* w;
	lapack_int m;
	lapack_logical tryrac;
	lapack_int* isuppzref;
	lapack_complex_float *zref;
	float *eref, *dref;
	float* wref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      stemr_scomplex_parameters (int matrix_layout, char jobz, char range,  lapack_int n, lapack_int m, lapack_int ldz);
      ~stemr_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
stemr_scomplex_parameters:: stemr_scomplex_parameters (int matrix_layout_i,  char jobz_i, char range_i, lapack_int n_i, lapack_int m_i, lapack_int ldz_i)
{
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	ldz = ldz_i;
	jobz = jobz_i;
	range = range_i;
	tryrac = 0;
	iu = il = vl = vu = 0;
	
	if (jobz == 'U')
		jobz = 'N';
	
	if (range == 'N')
	{	range = 'A';
		m = n;
		nzc = n;
	}else if (range == 'I')
	{
		iu = n/2;
		il = n/5;
		m = iu-il+1;
		nzc = iu-il+1;
	} else if (range == 'V')
	{	vl = n-5;
		vu= n+5;
		nzc = vu;

	}
	
	#if LAPACKE_TEST_VERBOSE
	printf(" \n stemr scomplex: matrix_layout = %d,  jobz:%c, range:%c, uplo:%c, n: %d , m:%d, ldz: %d, iu:%d, Il:%d \n", matrix_layout, jobz, range, uplo, n, m, ldz, iu, il);
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
		EXPECT_FALSE( true) << "stemr_float_parameters object: malloc error.";
		stemr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	//lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix(A, Aref, n, n, uplo);
	lapacke_gtest_init_float_buffer_pair_rand(e, eref, bufsize_e);
	lapacke_gtest_init_float_buffer_pair_rand(d, dref, bufsize_e);
} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
stemr_scomplex_parameters :: ~stemr_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" stemr_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   stemr_free();
  

}
/*  Test fixture class definition */
class cstemr_test  : public  ::testing::Test {
public:
   stemr_scomplex_parameters  *cstemr_obj;
   void SetUp();
   void TearDown () { delete cstemr_obj; }
};

void cstemr_test::SetUp(){

    /* LAPACKE cstemr prototype */
    typedef int (*Fptr_NL_LAPACKE_cstemr) (int matrix_layout, char jobz, char range, lapack_int n,  float* d, float* e,\
	float vl, float vu, lapack_int il, lapack_int iu, lapack_int* m, float* w, lapack_complex_float* z, lapack_int ldz,  lapack_int nzc, lapack_int* isuppz, lapack_logical* tryrac);

    Fptr_NL_LAPACKE_cstemr cstemr;

    cstemr_obj = new stemr_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
							eig_paramslist[idx].vect_rd,
							eig_paramslist[idx].compz,
							eig_paramslist[idx].m,
							eig_paramslist[idx].n,
							eig_paramslist[idx].m);

    idx = Circular_Increment_Index(idx);

    cstemr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cstemr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cstemr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cstemr_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    cstemr = (Fptr_NL_LAPACKE_cstemr)dlsym(cstemr_obj->hModule, "LAPACKE_cstemr");
    ASSERT_TRUE(cstemr != NULL) << "failed to get the Netlib LAPACKE_cstemr symbol";
    

    cstemr_obj->inforef = cstemr( cstemr_obj->matrix_layout,  cstemr_obj->jobz, cstemr_obj->range, cstemr_obj->n, 
	cstemr_obj->dref, cstemr_obj->eref, cstemr_obj->vl, cstemr_obj->vu, cstemr_obj->il, cstemr_obj->iu,
	&cstemr_obj->m, cstemr_obj->w, cstemr_obj->zref, cstemr_obj->ldz, cstemr_obj->nzc, cstemr_obj->isuppzref, &cstemr_obj->tryrac);

    /* Compute libflame's Lapacke o/p  */
    cstemr_obj->info = LAPACKE_cstemr( cstemr_obj->matrix_layout, cstemr_obj->jobz, cstemr_obj->range, cstemr_obj->n, 
	cstemr_obj->d, cstemr_obj->e, cstemr_obj->vl, cstemr_obj->vu, cstemr_obj->il, cstemr_obj->iu,  &cstemr_obj->m, 
	cstemr_obj->w, cstemr_obj->z, cstemr_obj->ldz, cstemr_obj->nzc, cstemr_obj->isuppz, &cstemr_obj->tryrac);

    if( cstemr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cstemr is wrong\n", cstemr_obj->info );
    }
    if( cstemr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cstemr is wrong\n", 
        cstemr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cstemr_obj->diff_e =  computeDiff_s( cstemr_obj->bufsize_e, cstemr_obj->e, cstemr_obj->eref );
	cstemr_obj->diff_d =  computeDiff_s( cstemr_obj->bufsize_e, cstemr_obj->e, cstemr_obj->eref );
	cstemr_obj->diff_z =  computeDiff_c( cstemr_obj->bufsize_z, cstemr_obj->z, cstemr_obj->zref );
	cstemr_obj->diff_w =  computeDiff_s( cstemr_obj->n, cstemr_obj->w, cstemr_obj->wref );
}

TEST_F(cstemr_test, cstemr1) {
    EXPECT_NEAR(0.0, cstemr_obj->diff_e, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, cstemr_obj->diff_d, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, cstemr_obj->diff_z, LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, cstemr_obj->diff_w, LAPACKE_EIG_THRESHOLD);

}

TEST_F(cstemr_test, cstemr2) {
    EXPECT_NEAR(0.0, cstemr_obj->diff_e, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, cstemr_obj->diff_d, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, cstemr_obj->diff_z, LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, cstemr_obj->diff_w, LAPACKE_EIG_THRESHOLD);
}

TEST_F(cstemr_test, cstemr3) {
    EXPECT_NEAR(0.0, cstemr_obj->diff_e, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, cstemr_obj->diff_d, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, cstemr_obj->diff_z, LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, cstemr_obj->diff_w, LAPACKE_EIG_THRESHOLD);
}

TEST_F(cstemr_test, cstemr4) {
   EXPECT_NEAR(0.0, cstemr_obj->diff_e, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, cstemr_obj->diff_d, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, cstemr_obj->diff_z, LAPACKE_EIG_THRESHOLD);
//	EXPECT_NEAR(0.0, cstemr_obj->diff_w, LAPACKE_EIG_THRESHOLD);
	
}
/* Begin dcomplex_common_parameters  class definition */
class stemr_dcomplex_parameters{

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
	lapack_int nzc;
	/*Output Parameter*/
	double* w;
	lapack_logical tryrac;
	lapack_int m;
	lapack_int* isuppzref;
	lapack_complex_double *zref;
	double *eref, *dref;
	double* wref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      stemr_dcomplex_parameters (int matrix_layout, char jobz, char range,  lapack_int n, lapack_int m, lapack_int ldz);
      ~stemr_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
stemr_dcomplex_parameters:: stemr_dcomplex_parameters (int matrix_layout_i,  char jobz_i, char range_i, lapack_int n_i, lapack_int m_i, lapack_int ldz_i)
{
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	ldz = ldz_i;
	jobz = jobz_i;
	range = range_i;
	tryrac = 0;
	iu = il = vl = vu = 0;
	
	if (jobz == 'U')
		jobz = 'N';
	
	if (range == 'N')
	{	range = 'A';
		m = n;
		nzc = n;
	}else if (range == 'I')
	{
		iu = n/2;
		il = n/5;
		m = iu-il+1;
		nzc = iu-il+1;
	} else if (range == 'V')
	{	vl = n-5;
		vu= n+5;
		nzc = vu;
	}
	
	#if LAPACKE_TEST_VERBOSE
	printf(" \n stemr dcomplex: matrix_layout = %d,  jobz:%c, range:%c, uplo:%c, n: %d , m:%d, ldz: %d, iu:%d, Il:%d \n", matrix_layout, jobz, range, uplo, n, m, ldz, iu, il);
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
		EXPECT_FALSE( true) << "stemr_double_parameters object: malloc error.";
		stemr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	//lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix(A, Aref, n, n, uplo);
	lapacke_gtest_init_double_buffer_pair_rand(e, eref, bufsize_e);
	lapacke_gtest_init_double_buffer_pair_rand(d, dref, bufsize_e);
} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
stemr_dcomplex_parameters :: ~stemr_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" stemr_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   stemr_free();
  

}
/*  Test fixture class definition */
class zstemr_test  : public  ::testing::Test {
public:
   stemr_dcomplex_parameters  *zstemr_obj;
   void SetUp();
   void TearDown () { delete zstemr_obj; }
};

void zstemr_test::SetUp(){

    /* LAPACKE zstemr prototype */
    typedef int (*Fptr_NL_LAPACKE_zstemr) (int matrix_layout, char jobz, char range, lapack_int n,  double* d, double* e,\
	double vl, double vu, lapack_int il, lapack_int iu, lapack_int* m, double* w, lapack_complex_double* z, lapack_int ldz, lapack_int nzc, lapack_int* isuppz, lapack_logical* tryrac);

    Fptr_NL_LAPACKE_zstemr zstemr;

    zstemr_obj = new stemr_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
							eig_paramslist[idx].vect_rd,
							eig_paramslist[idx].compz,
							eig_paramslist[idx].m,
							eig_paramslist[idx].n,
							eig_paramslist[idx].m);

    idx = Circular_Increment_Index(idx);

    zstemr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zstemr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zstemr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zstemr_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zstemr = (Fptr_NL_LAPACKE_zstemr)dlsym(zstemr_obj->hModule, "LAPACKE_zstemr");
    ASSERT_TRUE(zstemr != NULL) << "failed to get the Netlib LAPACKE_zstemr symbol";
    

    zstemr_obj->inforef = zstemr( zstemr_obj->matrix_layout,  zstemr_obj->jobz, zstemr_obj->range, zstemr_obj->n, 
	zstemr_obj->dref, zstemr_obj->eref, zstemr_obj->vl, zstemr_obj->vu, zstemr_obj->il, zstemr_obj->iu, 
	&zstemr_obj->m, zstemr_obj->w, zstemr_obj->zref, zstemr_obj->ldz, zstemr_obj->nzc, zstemr_obj->isuppzref, &zstemr_obj->tryrac);

    /* Compute libflame's Lapacke o/p  */
    zstemr_obj->info = LAPACKE_zstemr( zstemr_obj->matrix_layout, zstemr_obj->jobz, zstemr_obj->range, zstemr_obj->n, 
	zstemr_obj->d, zstemr_obj->e, zstemr_obj->vl, zstemr_obj->vu, zstemr_obj->il, zstemr_obj->iu,&zstemr_obj->m, 
	zstemr_obj->w, zstemr_obj->z, zstemr_obj->ldz, zstemr_obj->nzc, zstemr_obj->isuppz, &zstemr_obj->tryrac);

    if( zstemr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zstemr is wrong\n", zstemr_obj->info );
    }
    if( zstemr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zstemr is wrong\n", 
        zstemr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zstemr_obj->diff_e =  computeDiff_d( zstemr_obj->bufsize_e, zstemr_obj->e, zstemr_obj->eref );
	zstemr_obj->diff_d =  computeDiff_d( zstemr_obj->bufsize_e, zstemr_obj->e, zstemr_obj->eref );
	zstemr_obj->diff_z =  computeDiff_z( zstemr_obj->bufsize_z, zstemr_obj->z, zstemr_obj->zref );
	zstemr_obj->diff_w =  computeDiff_d( zstemr_obj->n, zstemr_obj->w, zstemr_obj->wref );
}

TEST_F(zstemr_test, zstemr1) {
    EXPECT_NEAR(0.0, zstemr_obj->diff_e, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zstemr_obj->diff_d, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zstemr_obj->diff_z, LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, zstemr_obj->diff_w, LAPACKE_EIG_THRESHOLD);

}

TEST_F(zstemr_test, zstemr2) {
    EXPECT_NEAR(0.0, zstemr_obj->diff_e, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zstemr_obj->diff_d, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zstemr_obj->diff_z, LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, zstemr_obj->diff_w, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zstemr_test, zstemr3) {
    EXPECT_NEAR(0.0, zstemr_obj->diff_e, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zstemr_obj->diff_d, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zstemr_obj->diff_z, LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, zstemr_obj->diff_w, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zstemr_test, zstemr4) {
   EXPECT_NEAR(0.0, zstemr_obj->diff_e, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zstemr_obj->diff_d, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zstemr_obj->diff_z, LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, zstemr_obj->diff_w, LAPACKE_EIG_THRESHOLD);
	
}