#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define heevx_2stage_free() \
if (ap!=NULL)    free(ap); \
if (apref!=NULL) free(apref);\
if (w!=NULL)  free(w);\
if (wref!=NULL) free(wref); \
if (z!=NULL)  free(z);\
if (zref!=NULL)  free(zref);\
if (ifail!=NULL)  free(ifail);\
if (ifailref!=NULL)  free(ifailref);\
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin scomplex_common_parameters  class definition */
class heevx_2stage_scomplex_parameters{

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
	lapack_int* ifail;
	lapack_complex_float* ap, *z;
	char uplo;
	char jobz;
	char range;
	/*Output Parameter*/
	float* w;
	lapack_int m;
	lapack_int* ifailref;
	lapack_complex_float *apref, *zref;
	float* wref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      heevx_2stage_scomplex_parameters (int matrix_layout, char jobz, char range,  char uplo, lapack_int n, lapack_int m, lapack_int ldz);
      ~heevx_2stage_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
heevx_2stage_scomplex_parameters:: heevx_2stage_scomplex_parameters (int matrix_layout_i,  char jobz_i, char range_i, char uplo_i, lapack_int n_i, lapack_int m_i, lapack_int ldz_i)
{
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	ldz = ldz_i;
	uplo = uplo_i;
	jobz = jobz_i;
	range = range_i;
	iu = il = vl = vu = 0;
	
	//if (jobz == 'U')
		jobz = 'N'; // V is not supported
	
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
	printf(" \n heevx_2stage scomplex: matrix_layout = %d,  jobz:%c, range:%c, uplo:%c, n: %d , m:%d, ldz: %d, iu:%d, Il:%d \n", matrix_layout, jobz, range, uplo, n, m, ldz, iu, il);
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
	lapacke_gtest_alloc_int_buffer_pair(&ifail, &ifailref, n);
	
	if ((ap==NULL) || (apref==NULL) || 
		(w == NULL) || (wref == NULL) ||
		(z == NULL) || (zref == NULL) ||
		(ifail == NULL) || (ifailref == NULL)){
		EXPECT_FALSE( true) << "heevx_2stage_float_parameters object: malloc error.";
		heevx_2stage_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	//lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix(A, Aref, n, n, uplo);
	lapacke_gtest_init_scomplex_buffer_pair_rand(ap, apref, bufsize_a);
} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
heevx_2stage_scomplex_parameters :: ~heevx_2stage_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" heevx_2stage_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   heevx_2stage_free();
  

}
/*  Test fixture class definition */
class cheevx_2stage_test  : public  ::testing::Test {
public:
   heevx_2stage_scomplex_parameters  *cheevx_2stage_obj;
   void SetUp();
   void TearDown () { delete cheevx_2stage_obj; }
};

void cheevx_2stage_test::SetUp(){

    /* LAPACKE cheevx_2stage prototype */
    typedef int (*Fptr_NL_LAPACKE_cheevx_2stage) ( int matrix_layout, char jobz, char range,	char uplo, lapack_int n, lapack_complex_float* ap, 
											lapack_int lda,  float vl, float vu, lapack_int il, lapack_int iu, float abstol, lapack_int* m, 
											float* w, lapack_complex_float* z, lapack_int ldz, lapack_int* ifail);

    Fptr_NL_LAPACKE_cheevx_2stage cheevx_2stage;

    cheevx_2stage_obj = new heevx_2stage_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
							eig_paramslist[idx].vect_rd,
							eig_paramslist[idx].compz,
							eig_paramslist[idx].uplo,
							eig_paramslist[idx].m,
							eig_paramslist[idx].n,
							eig_paramslist[idx].m);

    idx = Circular_Increment_Index(idx);

    cheevx_2stage_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cheevx_2stage_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cheevx_2stage_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cheevx_2stage_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    cheevx_2stage = (Fptr_NL_LAPACKE_cheevx_2stage)dlsym(cheevx_2stage_obj->hModule, "LAPACKE_cheevx_2stage");
    ASSERT_TRUE(cheevx_2stage != NULL) << "failed to get the Netlib LAPACKE_cheevx_2stage symbol";
    

    cheevx_2stage_obj->inforef = cheevx_2stage( cheevx_2stage_obj->matrix_layout,  cheevx_2stage_obj->jobz, cheevx_2stage_obj->range, cheevx_2stage_obj->uplo, cheevx_2stage_obj->n, 
	cheevx_2stage_obj->apref, cheevx_2stage_obj->lda, cheevx_2stage_obj->vl, cheevx_2stage_obj->vu, cheevx_2stage_obj->il, cheevx_2stage_obj->iu,	cheevx_2stage_obj->abstol, 
	&cheevx_2stage_obj->m, cheevx_2stage_obj->w, cheevx_2stage_obj->zref, cheevx_2stage_obj->ldz, cheevx_2stage_obj->ifailref);

    /* Compute libflame's Lapacke o/p  */
    cheevx_2stage_obj->info = LAPACKE_cheevx_2stage( cheevx_2stage_obj->matrix_layout, cheevx_2stage_obj->jobz, cheevx_2stage_obj->range, cheevx_2stage_obj->uplo, cheevx_2stage_obj->n, 
	cheevx_2stage_obj->ap, cheevx_2stage_obj->lda, cheevx_2stage_obj->vl, cheevx_2stage_obj->vu, cheevx_2stage_obj->il, cheevx_2stage_obj->iu, cheevx_2stage_obj->abstol, &cheevx_2stage_obj->m, 
	cheevx_2stage_obj->w, cheevx_2stage_obj->z, cheevx_2stage_obj->ldz, cheevx_2stage_obj->ifail);

    if( cheevx_2stage_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cheevx_2stage is wrong\n", cheevx_2stage_obj->info );
    }
    if( cheevx_2stage_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cheevx_2stage is wrong\n", 
        cheevx_2stage_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cheevx_2stage_obj->diff_a =  computeDiff_c( cheevx_2stage_obj->bufsize_a, cheevx_2stage_obj->ap, cheevx_2stage_obj->apref );
	cheevx_2stage_obj->diff_z =  computeDiff_c( cheevx_2stage_obj->bufsize_z, cheevx_2stage_obj->z, cheevx_2stage_obj->zref );
	cheevx_2stage_obj->diff_w =  computeDiff_s( cheevx_2stage_obj->n, cheevx_2stage_obj->w, cheevx_2stage_obj->wref );
}

TEST_F(cheevx_2stage_test, cheevx_2stage1) {
    EXPECT_NEAR(0.0, cheevx_2stage_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cheevx_2stage_test, cheevx_2stage2) {
    EXPECT_NEAR(0.0, cheevx_2stage_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cheevx_2stage_test, cheevx_2stage3) {
    EXPECT_NEAR(0.0, cheevx_2stage_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cheevx_2stage_test, cheevx_2stage4) {
    EXPECT_NEAR(0.0, cheevx_2stage_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	
}
/* Begin dcomplex_common_parameters  class definition */
class heevx_2stage_dcomplex_parameters{

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
	lapack_int* ifail;
	lapack_complex_double* ap, *z;
	char uplo;
	char jobz;
	char range;
	/*Output Parameter*/
	double* w;
	lapack_int m;
	lapack_int* ifailref;
	lapack_complex_double *apref, *zref;
	double* wref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      heevx_2stage_dcomplex_parameters (int matrix_layout, char jobz, char range,  char uplo, lapack_int n, lapack_int m, lapack_int ldz);
      ~heevx_2stage_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
heevx_2stage_dcomplex_parameters:: heevx_2stage_dcomplex_parameters (int matrix_layout_i,  char jobz_i, char range_i, char uplo_i, lapack_int n_i, lapack_int m_i, lapack_int ldz_i)
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
		jobz = 'N';// V is not supported
	
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
	printf(" \n heevx_2stage dcomplex: matrix_layout = %d,  jobz:%c, range:%c, uplo:%c, n: %d , m:%d, ldz: %d, iu:%d, Il:%d \n", matrix_layout, jobz, range, uplo, n, m, ldz, iu, il);
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
	lapacke_gtest_alloc_int_buffer_pair(&ifail, &ifailref, n);
	
	if ((ap==NULL) || (apref==NULL) || 
		(w == NULL) || (wref == NULL) ||
		(z == NULL) || (zref == NULL) ||
		(ifail == NULL) || (ifailref == NULL)){
		EXPECT_FALSE( true) << "heevx_2stage_double_parameters object: malloc error.";
		heevx_2stage_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	//lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix(A, Aref, n, n, uplo);
	lapacke_gtest_init_dcomplex_buffer_pair_rand(ap, apref, bufsize_a);
} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
heevx_2stage_dcomplex_parameters :: ~heevx_2stage_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" heevx_2stage_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   heevx_2stage_free();
  

}
/*  Test fixture class definition */
class zheevx_2stage_test  : public  ::testing::Test {
public:
   heevx_2stage_dcomplex_parameters  *zheevx_2stage_obj;
   void SetUp();
   void TearDown () { delete zheevx_2stage_obj; }
};

void zheevx_2stage_test::SetUp(){

    /* LAPACKE zheevx_2stage prototype */
    typedef int (*Fptr_NL_LAPACKE_zheevx_2stage) ( int matrix_layout, char jobz, char range,	char uplo, lapack_int n, lapack_complex_double* ap, 
											lapack_int lda,  double vl, double vu, lapack_int il, lapack_int iu, double abstol, lapack_int* m, 
											double* w, lapack_complex_double* z, lapack_int ldz, lapack_int* ifail);

    Fptr_NL_LAPACKE_zheevx_2stage zheevx_2stage;

    zheevx_2stage_obj = new heevx_2stage_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
							eig_paramslist[idx].vect_rd,
							eig_paramslist[idx].compz,
							eig_paramslist[idx].uplo,
							eig_paramslist[idx].m,
							eig_paramslist[idx].n,
							eig_paramslist[idx].m);

    idx = Circular_Increment_Index(idx);

    zheevx_2stage_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zheevx_2stage_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zheevx_2stage_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zheevx_2stage_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zheevx_2stage = (Fptr_NL_LAPACKE_zheevx_2stage)dlsym(zheevx_2stage_obj->hModule, "LAPACKE_zheevx_2stage");
    ASSERT_TRUE(zheevx_2stage != NULL) << "failed to get the Netlib LAPACKE_zheevx_2stage symbol";
    

    zheevx_2stage_obj->inforef = zheevx_2stage( zheevx_2stage_obj->matrix_layout,  zheevx_2stage_obj->jobz, zheevx_2stage_obj->range, zheevx_2stage_obj->uplo, zheevx_2stage_obj->n, 
	zheevx_2stage_obj->apref, zheevx_2stage_obj->lda, zheevx_2stage_obj->vl, zheevx_2stage_obj->vu, zheevx_2stage_obj->il, zheevx_2stage_obj->iu,	zheevx_2stage_obj->abstol, 
	&zheevx_2stage_obj->m, zheevx_2stage_obj->w, zheevx_2stage_obj->zref, zheevx_2stage_obj->ldz, zheevx_2stage_obj->ifailref);

    /* Compute libflame's Lapacke o/p  */
    zheevx_2stage_obj->info = LAPACKE_zheevx_2stage( zheevx_2stage_obj->matrix_layout, zheevx_2stage_obj->jobz, zheevx_2stage_obj->range, zheevx_2stage_obj->uplo, zheevx_2stage_obj->n, 
	zheevx_2stage_obj->ap, zheevx_2stage_obj->lda, zheevx_2stage_obj->vl, zheevx_2stage_obj->vu, zheevx_2stage_obj->il, zheevx_2stage_obj->iu, zheevx_2stage_obj->abstol, &zheevx_2stage_obj->m, 
	zheevx_2stage_obj->w, zheevx_2stage_obj->z, zheevx_2stage_obj->ldz, zheevx_2stage_obj->ifail);

    if( zheevx_2stage_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zheevx_2stage is wrong\n", zheevx_2stage_obj->info );
    }
    if( zheevx_2stage_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zheevx_2stage is wrong\n", 
        zheevx_2stage_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zheevx_2stage_obj->diff_a =  computeDiff_z( zheevx_2stage_obj->bufsize_a, zheevx_2stage_obj->ap, zheevx_2stage_obj->apref );
	zheevx_2stage_obj->diff_z =  computeDiff_z( zheevx_2stage_obj->bufsize_z, zheevx_2stage_obj->z, zheevx_2stage_obj->zref );
	zheevx_2stage_obj->diff_w =  computeDiff_d( zheevx_2stage_obj->n, zheevx_2stage_obj->w, zheevx_2stage_obj->wref );
}

TEST_F(zheevx_2stage_test, zheevx_2stage1) {
    EXPECT_NEAR(0.0, zheevx_2stage_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	
}

TEST_F(zheevx_2stage_test, zheevx_2stage2) {
    EXPECT_NEAR(0.0, zheevx_2stage_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	
}

TEST_F(zheevx_2stage_test, zheevx_2stage3) {
    EXPECT_NEAR(0.0, zheevx_2stage_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	
}

TEST_F(zheevx_2stage_test, zheevx_2stage4) {
    EXPECT_NEAR(0.0, zheevx_2stage_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
}
