#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define heevx_free() \
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
class heevx_scomplex_parameters{

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
      heevx_scomplex_parameters (int matrix_layout, char jobz, char range,  char uplo, lapack_int n, lapack_int m, lapack_int ldz);
      ~heevx_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
heevx_scomplex_parameters:: heevx_scomplex_parameters (int matrix_layout_i,  char jobz_i, char range_i, char uplo_i, lapack_int n_i, lapack_int m_i, lapack_int ldz_i)
{
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	ldz = ldz_i;
	uplo = uplo_i;
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
	printf(" \n heevx scomplex: matrix_layout = %d,  jobz:%c, range:%c, uplo:%c, n: %d , m:%d, ldz: %d, iu:%d, Il:%d \n", matrix_layout, jobz, range, uplo, n, m, ldz, iu, il);
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

	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&ap, &apref, bufsize_a);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&z, &zref, bufsize_z);
	lapacke_gtest_alloc_float_buffer_pair(&w, &wref, n);
	lapacke_gtest_alloc_int_buffer_pair(&ifail, &ifailref, n);
	
	if ((ap==NULL) || (apref==NULL) || 
		(w == NULL) || (wref == NULL) ||
		(z == NULL) || (zref == NULL) ||
		(ifail == NULL) || (ifailref == NULL)){
		EXPECT_FALSE( true) << "heevx_float_parameters object: malloc error.";
		heevx_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	//lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix(A, Aref, n, n, uplo);
	lapacke_gtest_init_scomplex_buffer_pair_rand(ap, apref, bufsize_a);
} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
heevx_scomplex_parameters :: ~heevx_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" heevx_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   heevx_free();
  

}
/*  Test fixture class definition */
class cheevx_test  : public  ::testing::Test {
public:
   heevx_scomplex_parameters  *cheevx_obj;
   void SetUp();
   void TearDown () { delete cheevx_obj; }
};

void cheevx_test::SetUp(){

    /* LAPACKE cheevx prototype */
    typedef int (*Fptr_NL_LAPACKE_cheevx) ( int matrix_layout, char jobz, char range,	char uplo, lapack_int n, lapack_complex_float* ap, 
											lapack_int lda,  float vl, float vu, lapack_int il, lapack_int iu, float abstol, lapack_int* m, 
											float* w, lapack_complex_float* z, lapack_int ldz, lapack_int* ifail);

    Fptr_NL_LAPACKE_cheevx cheevx;

    cheevx_obj = new heevx_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
							eig_paramslist[idx].vect_rd,
							eig_paramslist[idx].compz,
							eig_paramslist[idx].uplo,
							eig_paramslist[idx].m,
							eig_paramslist[idx].n,
							eig_paramslist[idx].m);

    idx = Circular_Increment_Index(idx);

    cheevx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cheevx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cheevx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cheevx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    cheevx = (Fptr_NL_LAPACKE_cheevx)dlsym(cheevx_obj->hModule, "LAPACKE_cheevx");
    ASSERT_TRUE(cheevx != NULL) << "failed to get the Netlib LAPACKE_cheevx symbol";
    

    cheevx_obj->inforef = cheevx( cheevx_obj->matrix_layout,  cheevx_obj->jobz, cheevx_obj->range, cheevx_obj->uplo, cheevx_obj->n, 
	cheevx_obj->apref, cheevx_obj->lda, cheevx_obj->vl, cheevx_obj->vu, cheevx_obj->il, cheevx_obj->iu,	cheevx_obj->abstol, 
	&cheevx_obj->m, cheevx_obj->w, cheevx_obj->zref, cheevx_obj->ldz, cheevx_obj->ifailref);

    /* Compute libflame's Lapacke o/p  */
    cheevx_obj->info = LAPACKE_cheevx( cheevx_obj->matrix_layout, cheevx_obj->jobz, cheevx_obj->range, cheevx_obj->uplo, cheevx_obj->n, 
	cheevx_obj->ap, cheevx_obj->lda, cheevx_obj->vl, cheevx_obj->vu, cheevx_obj->il, cheevx_obj->iu, cheevx_obj->abstol, &cheevx_obj->m, 
	cheevx_obj->w, cheevx_obj->z, cheevx_obj->ldz, cheevx_obj->ifail);

    if( cheevx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cheevx is wrong\n", cheevx_obj->info );
    }
    if( cheevx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cheevx is wrong\n", 
        cheevx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cheevx_obj->diff_a =  computeDiff_c( cheevx_obj->bufsize_a, cheevx_obj->ap, cheevx_obj->apref );
	cheevx_obj->diff_z =  computeDiff_c( cheevx_obj->bufsize_z, cheevx_obj->z, cheevx_obj->zref );
	cheevx_obj->diff_w =  computeDiff_s( cheevx_obj->n, cheevx_obj->w, cheevx_obj->wref );
}

TEST_F(cheevx_test, cheevx1) {
    EXPECT_NEAR(0.0, cheevx_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cheevx_test, cheevx2) {
    EXPECT_NEAR(0.0, cheevx_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cheevx_test, cheevx3) {
    EXPECT_NEAR(0.0, cheevx_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cheevx_test, cheevx4) {
    EXPECT_NEAR(0.0, cheevx_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	
}
/* Begin dcomplex_common_parameters  class definition */
class heevx_dcomplex_parameters{

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
      heevx_dcomplex_parameters (int matrix_layout, char jobz, char range,  char uplo, lapack_int n, lapack_int m, lapack_int ldz);
      ~heevx_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
heevx_dcomplex_parameters:: heevx_dcomplex_parameters (int matrix_layout_i,  char jobz_i, char range_i, char uplo_i, lapack_int n_i, lapack_int m_i, lapack_int ldz_i)
{
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	ldz = ldz_i;
	uplo = uplo_i;
	jobz = jobz_i;
	range = range_i;
	iu = il =vl = vu = 0;
	
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
	{	
		vl = n-5;
		vu = n+5;
	}
	
	#if LAPACKE_TEST_VERBOSE
	printf(" \n heevx dcomplex: matrix_layout = %d,  jobz:%c, range:%c, uplo:%c, n: %d , m:%d, ldz: %d, iu:%d, Il:%d \n", matrix_layout, jobz, range, uplo, n, m, ldz, iu, il);
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
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&ap, &apref, bufsize_a);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&z, &zref, bufsize_z);
	lapacke_gtest_alloc_double_buffer_pair(&w, &wref, n);
	lapacke_gtest_alloc_int_buffer_pair(&ifail, &ifailref, n);
	
	if ((ap==NULL) || (apref==NULL) || 
		(w == NULL) || (wref == NULL) ||
		(z == NULL) || (zref == NULL) ||
		(ifail == NULL) || (ifailref == NULL)){
		EXPECT_FALSE( true) << "heevx_double_parameters object: malloc error.";
		heevx_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	//lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix(A, Aref, n, n, uplo);
	lapacke_gtest_init_dcomplex_buffer_pair_rand(ap, apref, bufsize_a);
} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
heevx_dcomplex_parameters :: ~heevx_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" heevx_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   heevx_free();
  

}
/*  Test fixture class definition */
class zheevx_test  : public  ::testing::Test {
public:
   heevx_dcomplex_parameters  *zheevx_obj;
   void SetUp();
   void TearDown () { delete zheevx_obj; }
};

void zheevx_test::SetUp(){

    /* LAPACKE zheevx prototype */
    typedef int (*Fptr_NL_LAPACKE_zheevx) ( int matrix_layout, char jobz, char range,	char uplo, lapack_int n, lapack_complex_double* ap, 
											lapack_int lda,  double vl, double vu, lapack_int il, lapack_int iu, double abstol, lapack_int* m, 
											double* w, lapack_complex_double* z, lapack_int ldz, lapack_int* ifail);

    Fptr_NL_LAPACKE_zheevx zheevx;

    zheevx_obj = new heevx_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
							eig_paramslist[idx].vect_rd,
							eig_paramslist[idx].compz,
							eig_paramslist[idx].uplo,
							eig_paramslist[idx].m,
							eig_paramslist[idx].n,
							eig_paramslist[idx].m);

    idx = Circular_Increment_Index(idx);

    zheevx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zheevx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zheevx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zheevx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zheevx = (Fptr_NL_LAPACKE_zheevx)dlsym(zheevx_obj->hModule, "LAPACKE_zheevx");
    ASSERT_TRUE(zheevx != NULL) << "failed to get the Netlib LAPACKE_zheevx symbol";
    

    zheevx_obj->inforef = zheevx( zheevx_obj->matrix_layout,  zheevx_obj->jobz, zheevx_obj->range, zheevx_obj->uplo, zheevx_obj->n, 
	zheevx_obj->apref, zheevx_obj->lda, zheevx_obj->vl, zheevx_obj->vu, zheevx_obj->il, zheevx_obj->iu,	zheevx_obj->abstol, 
	&zheevx_obj->m, zheevx_obj->w, zheevx_obj->zref, zheevx_obj->ldz, zheevx_obj->ifailref);

    /* Compute libflame's Lapacke o/p  */
    zheevx_obj->info = LAPACKE_zheevx( zheevx_obj->matrix_layout, zheevx_obj->jobz, zheevx_obj->range, zheevx_obj->uplo, zheevx_obj->n, 
	zheevx_obj->ap, zheevx_obj->lda, zheevx_obj->vl, zheevx_obj->vu, zheevx_obj->il, zheevx_obj->iu, zheevx_obj->abstol, &zheevx_obj->m, 
	zheevx_obj->w, zheevx_obj->z, zheevx_obj->ldz, zheevx_obj->ifail);

    if( zheevx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zheevx is wrong\n", zheevx_obj->info );
    }
    if( zheevx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zheevx is wrong\n", 
        zheevx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zheevx_obj->diff_a =  computeDiff_z( zheevx_obj->bufsize_a, zheevx_obj->ap, zheevx_obj->apref );
	zheevx_obj->diff_z =  computeDiff_z( zheevx_obj->bufsize_z, zheevx_obj->z, zheevx_obj->zref );
	zheevx_obj->diff_w =  computeDiff_d( zheevx_obj->n, zheevx_obj->w, zheevx_obj->wref );
}

TEST_F(zheevx_test, zheevx1) {
    EXPECT_NEAR(0.0, zheevx_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	
}

TEST_F(zheevx_test, zheevx2) {
    EXPECT_NEAR(0.0, zheevx_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	
}

TEST_F(zheevx_test, zheevx3) {
    EXPECT_NEAR(0.0, zheevx_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	
}

TEST_F(zheevx_test, zheevx4) {
    EXPECT_NEAR(0.0, zheevx_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
}
