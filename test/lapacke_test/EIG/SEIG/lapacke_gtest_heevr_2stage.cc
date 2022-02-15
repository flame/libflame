#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"

#define heevr_2stage_free() \
if (ap!=NULL)    free(ap); \
if (apref!=NULL) free(apref);\
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

/* Begin scomplex_common_parameters  class definition */
class heevr_2stage_scomplex_parameters{

   public:
	int bufsize_a;
	int bufsize_z;
	void *hModule, *dModule;
	void *lmodule, *bmodule;
	float diff_a;
	float diff_w;
	float diff_z;
	float vl, vu;
	float abstol;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int ldz,lda;
	lapack_int il, iu;
	lapack_int* isuppz;
	lapack_complex_float* ap, *z;
	char uplo;
	char jobz;
	char range;
	/*Output Parameter*/
	float* w;
	lapack_int m;
	lapack_int* isuppzref;
	lapack_complex_float *apref, *zref;
	float* wref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      heevr_2stage_scomplex_parameters (int matrix_layout, char jobz, char range,  char uplo, lapack_int n, lapack_int m, lapack_int ldz);
      ~heevr_2stage_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
heevr_2stage_scomplex_parameters:: heevr_2stage_scomplex_parameters (int matrix_layout_i,  char jobz_i, char range_i, char uplo_i, lapack_int n_i, lapack_int m_i, lapack_int ldz_i)
{
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	ldz = ldz_i;
	uplo = uplo_i;
	jobz = jobz_i;
	range = range_i;
	iu = il = vl = vu = 0;
	
//	if (jobz == 'U')
		jobz = 'N'; //'V' is not supported
	
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
	printf(" \n heevr_2stage scomplex: matrix_layout = %d,  jobz:%c, range:%c, uplo:%c, n: %d , m:%d, ldz: %d, iu:%d, Il:%d \n", matrix_layout, jobz, range, uplo, n, m, ldz, iu, il);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{
		ldz = n;
		bufsize_z = ldz*m;
		lda = n;	
		bufsize_a = lda*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
    {
		ldz = m;
		lda = n;
		bufsize_a = lda*n;
		bufsize_z = ldz*n;
    }else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	ldz =n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&ap, &apref, bufsize_a);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&z, &zref, bufsize_z);
	lapacke_gtest_alloc_float_buffer_pair(&w, &wref, n);
	lapacke_gtest_alloc_int_buffer_pair(&isuppz, &isuppzref, 2*m);
	
	if ((ap==NULL) || (apref==NULL) || 
		(w == NULL) || (wref == NULL) ||
		(z == NULL) || (zref == NULL) ||
		(isuppz == NULL) || (isuppzref == NULL)){
		EXPECT_FALSE( true) << "heevr_2stage_float_parameters object: malloc error.";
		heevr_2stage_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	//lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix(A, Aref, n, n, uplo);
	lapacke_gtest_init_scomplex_buffer_pair_rand(ap, apref, bufsize_a);
} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
heevr_2stage_scomplex_parameters :: ~heevr_2stage_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" heevr_2stage_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   heevr_2stage_free();
  

}
/*  Test fixture class definition */
class cheevr_2stage_test  : public  ::testing::Test {
public:
   heevr_2stage_scomplex_parameters  *cheevr_2stage_obj;
   void SetUp();
   void TearDown () { delete cheevr_2stage_obj; }
};

void cheevr_2stage_test::SetUp(){

    /* LAPACKE cheevr_2stage prototype */
    typedef int (*Fptr_NL_LAPACKE_cheevr_2stage) ( int matrix_layout, char jobz, char range,	char uplo, lapack_int n, lapack_complex_float* ap, 
											lapack_int lda,  float vl, float vu, lapack_int il, lapack_int iu, float abstol, lapack_int* m, 
											float* w, lapack_complex_float* z, lapack_int ldz, lapack_int* isuppz);

    Fptr_NL_LAPACKE_cheevr_2stage cheevr_2stage;

    cheevr_2stage_obj = new heevr_2stage_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
							eig_paramslist[idx].vect_rd,
							eig_paramslist[idx].compz,
							eig_paramslist[idx].uplo,
							eig_paramslist[idx].m,
							eig_paramslist[idx].n,
							eig_paramslist[idx].m);

    idx = Circular_Increment_Index(idx);

    cheevr_2stage_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cheevr_2stage_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cheevr_2stage_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cheevr_2stage_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    cheevr_2stage = (Fptr_NL_LAPACKE_cheevr_2stage)dlsym(cheevr_2stage_obj->hModule, "LAPACKE_cheevr_2stage");
    ASSERT_TRUE(cheevr_2stage != NULL) << "failed to get the Netlib LAPACKE_cheevr_2stage symbol";
    

    cheevr_2stage_obj->inforef = cheevr_2stage( cheevr_2stage_obj->matrix_layout,  cheevr_2stage_obj->jobz, cheevr_2stage_obj->range, cheevr_2stage_obj->uplo, cheevr_2stage_obj->n, 
	cheevr_2stage_obj->apref, cheevr_2stage_obj->lda, cheevr_2stage_obj->vl, cheevr_2stage_obj->vu, cheevr_2stage_obj->il, cheevr_2stage_obj->iu,	cheevr_2stage_obj->abstol, 
	&cheevr_2stage_obj->m, cheevr_2stage_obj->w, cheevr_2stage_obj->zref, cheevr_2stage_obj->ldz, cheevr_2stage_obj->isuppzref);

    /* Compute libflame's Lapacke o/p  */
    cheevr_2stage_obj->info = LAPACKE_cheevr_2stage( cheevr_2stage_obj->matrix_layout, cheevr_2stage_obj->jobz, cheevr_2stage_obj->range, cheevr_2stage_obj->uplo, cheevr_2stage_obj->n, 
	cheevr_2stage_obj->ap, cheevr_2stage_obj->lda, cheevr_2stage_obj->vl, cheevr_2stage_obj->vu, cheevr_2stage_obj->il, cheevr_2stage_obj->iu, cheevr_2stage_obj->abstol, &cheevr_2stage_obj->m, 
	cheevr_2stage_obj->w, cheevr_2stage_obj->z, cheevr_2stage_obj->ldz, cheevr_2stage_obj->isuppz);

    if( cheevr_2stage_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cheevr_2stage is wrong\n", cheevr_2stage_obj->info );
    }
    if( cheevr_2stage_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cheevr_2stage is wrong\n", 
        cheevr_2stage_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cheevr_2stage_obj->diff_a =  computeDiff_c( cheevr_2stage_obj->bufsize_a, cheevr_2stage_obj->ap, cheevr_2stage_obj->apref );
	cheevr_2stage_obj->diff_z =  computeDiff_c( cheevr_2stage_obj->bufsize_z, cheevr_2stage_obj->z, cheevr_2stage_obj->zref );
	cheevr_2stage_obj->diff_w =  computeDiff_s( cheevr_2stage_obj->n, cheevr_2stage_obj->w, cheevr_2stage_obj->wref );
}

TEST_F(cheevr_2stage_test, cheevr_2stage1) {
    EXPECT_NEAR(0.0, cheevr_2stage_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cheevr_2stage_test, cheevr_2stage2) {
    EXPECT_NEAR(0.0, cheevr_2stage_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cheevr_2stage_test, cheevr_2stage3) {
    EXPECT_NEAR(0.0, cheevr_2stage_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cheevr_2stage_test, cheevr_2stage4) {
    EXPECT_NEAR(0.0, cheevr_2stage_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	
}
/* Begin dcomplex_common_parameters  class definition */
class heevr_2stage_dcomplex_parameters{

   public:
	int bufsize_a;
	int bufsize_z;
	void *hModule, *dModule;
	void *lmodule, *bmodule;
	double diff_a;
	double diff_w;
	double diff_z;
	double vl, vu;
	double abstol;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int ldz,lda;
	lapack_int il, iu;
	lapack_int* isuppz;
	lapack_complex_double* ap, *z;
	char uplo;
	char jobz;
	char range;
	/*Output Parameter*/
	double* w;
	lapack_int m;
	lapack_int* isuppzref;
	lapack_complex_double *apref, *zref;
	double* wref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      heevr_2stage_dcomplex_parameters (int matrix_layout, char jobz, char range,  char uplo, lapack_int n, lapack_int m, lapack_int ldz);
      ~heevr_2stage_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
heevr_2stage_dcomplex_parameters:: heevr_2stage_dcomplex_parameters (int matrix_layout_i,  char jobz_i, char range_i, char uplo_i, lapack_int n_i, lapack_int m_i, lapack_int ldz_i)
{
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	ldz = ldz_i;
	uplo = uplo_i;
	jobz = jobz_i;
	range = range_i;
	iu = il =vl = vu = 0;
	
	//if (jobz == 'U')
		jobz = 'N'; //'V' is not supported
	
	if (range == 'N')
	{	range = 'A';
		m = n;
	}else if (range == 'I')
	{
		iu = n/2;
		il = n/5;
		m = iu-il+1;		
	} else if (range == 'V')
	{	
		vl = n-5;
		vu = n+5;
	}
	
	#if LAPACKE_TEST_VERBOSE
	printf(" \n heevr_2stage dcomplex: matrix_layout = %d,  jobz:%c, range:%c, uplo:%c, n: %d , m:%d, ldz: %d, iu:%d, Il:%d \n", matrix_layout, jobz, range, uplo, n, m, ldz, iu, il);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{
		ldz = n;
		bufsize_z = ldz*m;
		lda = n;	
		bufsize_a = lda*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
    {
		ldz = m;
		lda = n;
		bufsize_a = lda*n;
		bufsize_z = ldz*n;
    }else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	ldz =n;
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&ap, &apref, bufsize_a);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&z, &zref, bufsize_z);
	lapacke_gtest_alloc_double_buffer_pair(&w, &wref, n);
	lapacke_gtest_alloc_int_buffer_pair(&isuppz, &isuppzref, 2*m);
	
	if ((ap==NULL) || (apref==NULL) || 
		(w == NULL) || (wref == NULL) ||
		(z == NULL) || (zref == NULL) ||
		(isuppz == NULL) || (isuppzref == NULL)){
		EXPECT_FALSE( true) << "heevr_2stage_double_parameters object: malloc error.";
		heevr_2stage_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	//lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix(A, Aref, n, n, uplo);
	lapacke_gtest_init_dcomplex_buffer_pair_rand(ap, apref, bufsize_a);
} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
heevr_2stage_dcomplex_parameters :: ~heevr_2stage_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" heevr_2stage_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   heevr_2stage_free();
  

}
/*  Test fixture class definition */
class zheevr_2stage_test  : public  ::testing::Test {
public:
   heevr_2stage_dcomplex_parameters  *zheevr_2stage_obj;
   void SetUp();
   void TearDown () { delete zheevr_2stage_obj; }
};

void zheevr_2stage_test::SetUp(){

    /* LAPACKE zheevr_2stage prototype */
    typedef int (*Fptr_NL_LAPACKE_zheevr_2stage) ( int matrix_layout, char jobz, char range, char uplo, lapack_int n, lapack_complex_double* ap, 
											lapack_int lda,  double vl, double vu, lapack_int il, lapack_int iu, double abstol, lapack_int* m, 
											double* w, lapack_complex_double* z, lapack_int ldz, lapack_int* isuppz);

    Fptr_NL_LAPACKE_zheevr_2stage zheevr_2stage;

    zheevr_2stage_obj = new heevr_2stage_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
							eig_paramslist[idx].vect_rd,
							eig_paramslist[idx].compz,
							eig_paramslist[idx].uplo,
							eig_paramslist[idx].m,
							eig_paramslist[idx].n,
							eig_paramslist[idx].m);

    idx = Circular_Increment_Index(idx);

    zheevr_2stage_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zheevr_2stage_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zheevr_2stage_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zheevr_2stage_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zheevr_2stage = (Fptr_NL_LAPACKE_zheevr_2stage)dlsym(zheevr_2stage_obj->hModule, "LAPACKE_zheevr_2stage");
    ASSERT_TRUE(zheevr_2stage != NULL) << "failed to get the Netlib LAPACKE_zheevr_2stage symbol";
    

    zheevr_2stage_obj->inforef = zheevr_2stage( zheevr_2stage_obj->matrix_layout,  zheevr_2stage_obj->jobz, zheevr_2stage_obj->range, zheevr_2stage_obj->uplo, zheevr_2stage_obj->n, 
	zheevr_2stage_obj->apref, zheevr_2stage_obj->lda, zheevr_2stage_obj->vl, zheevr_2stage_obj->vu, zheevr_2stage_obj->il, zheevr_2stage_obj->iu,	zheevr_2stage_obj->abstol, 
	&zheevr_2stage_obj->m, zheevr_2stage_obj->w, zheevr_2stage_obj->zref, zheevr_2stage_obj->ldz, zheevr_2stage_obj->isuppzref);

    /* Compute libflame's Lapacke o/p  */
    zheevr_2stage_obj->info = LAPACKE_zheevr_2stage( zheevr_2stage_obj->matrix_layout, zheevr_2stage_obj->jobz, zheevr_2stage_obj->range, zheevr_2stage_obj->uplo, zheevr_2stage_obj->n, 
	zheevr_2stage_obj->ap, zheevr_2stage_obj->lda, zheevr_2stage_obj->vl, zheevr_2stage_obj->vu, zheevr_2stage_obj->il, zheevr_2stage_obj->iu, zheevr_2stage_obj->abstol, &zheevr_2stage_obj->m, 
	zheevr_2stage_obj->w, zheevr_2stage_obj->z, zheevr_2stage_obj->ldz, zheevr_2stage_obj->isuppz);

    if( zheevr_2stage_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zheevr_2stage is wrong\n", zheevr_2stage_obj->info );
    }
    if( zheevr_2stage_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zheevr_2stage is wrong\n", 
        zheevr_2stage_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zheevr_2stage_obj->diff_a =  computeDiff_z( zheevr_2stage_obj->bufsize_a, zheevr_2stage_obj->ap, zheevr_2stage_obj->apref );
	zheevr_2stage_obj->diff_z =  computeDiff_z( zheevr_2stage_obj->bufsize_z, zheevr_2stage_obj->z, zheevr_2stage_obj->zref );
	zheevr_2stage_obj->diff_w =  computeDiff_d( zheevr_2stage_obj->n, zheevr_2stage_obj->w, zheevr_2stage_obj->wref );
}

TEST_F(zheevr_2stage_test, zheevr_2stage1) {
    EXPECT_NEAR(0.0, zheevr_2stage_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	
}

TEST_F(zheevr_2stage_test, zheevr_2stage2) {
    EXPECT_NEAR(0.0, zheevr_2stage_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	
}

TEST_F(zheevr_2stage_test, zheevr_2stage3) {
    EXPECT_NEAR(0.0, zheevr_2stage_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	
}

TEST_F(zheevr_2stage_test, zheevr_2stage4) {
    EXPECT_NEAR(0.0, zheevr_2stage_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
}
