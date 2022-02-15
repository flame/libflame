#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"

#define hegvx_free() \
if (ap!=NULL)    free(ap); \
if (apref!=NULL) free(apref);\
if (w!=NULL)  free(w);\
if (wref!=NULL) free(wref); \
if (bp!=NULL)  free(bp);\
if (bpref!=NULL)  free(bpref);\
if (z!=NULL)  free(z);\
if (zref!=NULL)  free(zref);\
if (ifail!=NULL)  free(ifail);\
if (ifailref!=NULL)  free(ifailref);\
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin scomplex_common_parameters  class definition */
class hegvx_scomplex_parameters{

   public:
	int bufsize_a;
	int bufsize_b;
	int bufsize_z;
	void *hModule, *dModule;
	float diff_a;
	float diff_w;
	float diff_b;
	float diff_z;
	float vl, vu;
	float abstol;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int ldz,lda, ldb;	
	lapack_int il, iu;
	lapack_int* ifail;
	lapack_complex_float* ap, *bp, *z;
	char uplo;
	char jobz;
	char range;
	lapack_int itype;
	/*Output Parameter*/
	float* w;
	lapack_int m;
	lapack_int* ifailref;
	lapack_complex_float *apref, *bpref, *zref;
	float* wref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      hegvx_scomplex_parameters (int matrix_layout, lapack_int itype, char jobz, char range,  char uplo, lapack_int n, lapack_int m, lapack_int ldz);
      ~hegvx_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
hegvx_scomplex_parameters:: hegvx_scomplex_parameters (int matrix_layout_i, lapack_int itype_i, char jobz_i, char range_i, char uplo_i, lapack_int n_i, lapack_int m_i, lapack_int ldz_i)
{
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	ldz = ldz_i;
	uplo = uplo_i;
	jobz = jobz_i;
	itype = itype_i;
	range = range_i;
	iu = il =vl = vu = 0;
	
	if (jobz == 'U')
		jobz = 'N';
	
	if (itype == 4)
		itype = 1;
	
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
	printf(" \n hegvx scomplex: matrix_layout = %d, itype:%d, jobz:%c, range:%c, uplo:%c, n: %d , m:%d, ldz: %d, iu:%d, Il:%d \n", matrix_layout, itype, jobz, range, uplo, n, m, ldz, iu, il);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{
		ldz = n;
		bufsize_z = ldz*m;
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
    {
		ldz = m;
		if (range == 'V')
				ldz = n;
		bufsize_z = ldz*n;
    }else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	lda = ldb = n;
	
	bufsize_a =  lda*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&ap, &apref, bufsize_a);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&bp, &bpref, bufsize_a);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&z, &zref, bufsize_z);
	lapacke_gtest_alloc_float_buffer_pair(&w, &wref, n);
	lapacke_gtest_alloc_int_buffer_pair(&ifail, &ifailref,	n);
	
	if ((ap==NULL) || (apref==NULL) || 
		(w == NULL) || (wref == NULL) ||
		(bp == NULL) || (bpref == NULL) ||
		(z == NULL) || (zref == NULL) ||
		(ifail == NULL) || (ifailref == NULL)){
		EXPECT_FALSE( true) << "hegvx_float_parameters object: malloc error.";
		hegvx_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix(ap, apref, bufsize_a, n, uplo);
	//lapacke_gtest_init_scomplex_buffer_pair_rand(ap, apref, bufsize_a);
	lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix(bp, bpref, bufsize_a, n, uplo);
} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
hegvx_scomplex_parameters :: ~hegvx_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" hegvx_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   hegvx_free();
  

}
/*  Test fixture class definition */
class chegvx_test  : public  ::testing::Test {
public:
   hegvx_scomplex_parameters  *chegvx_obj;
   void SetUp();
   void TearDown () { delete chegvx_obj; }
};

void chegvx_test::SetUp(){

    /* LAPACKE chegvx prototype */
    typedef int (*Fptr_NL_LAPACKE_chegvx) ( int matrix_layout, lapack_int itype, char jobz, char range,
										char uplo, lapack_int n, lapack_complex_float* ap, lapack_int lda, lapack_complex_float* bp, lapack_int ldb,
										float vl, float vu,	lapack_int il, lapack_int iu, float abstol, lapack_int* m, float* w, lapack_complex_float* z,
										lapack_int ldz, lapack_int* ifail);

    Fptr_NL_LAPACKE_chegvx chegvx;

    chegvx_obj = new hegvx_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
							eig_paramslist[idx].itype,
							eig_paramslist[idx].vect_rd,
							eig_paramslist[idx].compz,
							eig_paramslist[idx].uplo,
							eig_paramslist[idx].m,
							eig_paramslist[idx].n,
							eig_paramslist[idx].m);

    idx = Circular_Increment_Index(idx);

    chegvx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    chegvx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(chegvx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(chegvx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    chegvx = (Fptr_NL_LAPACKE_chegvx)dlsym(chegvx_obj->hModule, "LAPACKE_chegvx");
    ASSERT_TRUE(chegvx != NULL) << "failed to get the Netlib LAPACKE_chegvx symbol";
    

    chegvx_obj->inforef = chegvx( chegvx_obj->matrix_layout, chegvx_obj->itype, chegvx_obj->jobz, chegvx_obj->range, chegvx_obj->uplo, 
								chegvx_obj->n, chegvx_obj->apref, chegvx_obj->lda, chegvx_obj->bpref, chegvx_obj->ldb, chegvx_obj->vl, chegvx_obj->vu, chegvx_obj->il, chegvx_obj->iu,
								chegvx_obj->abstol, &chegvx_obj->m, chegvx_obj->wref, chegvx_obj->zref, chegvx_obj->ldz, chegvx_obj->ifailref);

    /* Compute libflame's Lapacke o/p  */
    chegvx_obj->info = LAPACKE_chegvx( chegvx_obj->matrix_layout, chegvx_obj->itype, chegvx_obj->jobz, chegvx_obj->range, chegvx_obj->uplo, 
										chegvx_obj->n, chegvx_obj->ap, chegvx_obj->lda, chegvx_obj->bp, chegvx_obj->ldb, chegvx_obj->vl, chegvx_obj->vu, chegvx_obj->il, chegvx_obj->iu,
										chegvx_obj->abstol, &chegvx_obj->m, chegvx_obj->w, chegvx_obj->z, chegvx_obj->ldz, chegvx_obj->ifail);

    if( chegvx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_chegvx is wrong\n", chegvx_obj->info );
    }
    if( chegvx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chegvx is wrong\n", 
        chegvx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    chegvx_obj->diff_a =  computeDiff_c( chegvx_obj->bufsize_a, chegvx_obj->ap, chegvx_obj->apref );
	chegvx_obj->diff_b =  computeDiff_c( chegvx_obj->bufsize_a, chegvx_obj->bp, chegvx_obj->bpref );
	chegvx_obj->diff_z =  computeDiff_c( chegvx_obj->bufsize_z, chegvx_obj->z, chegvx_obj->zref );
	chegvx_obj->diff_w =  computeDiff_s( chegvx_obj->n, chegvx_obj->w, chegvx_obj->wref );
}

TEST_F(chegvx_test, chegvx1) {
    EXPECT_NEAR(0.0, chegvx_obj->diff_a, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, chegvx_obj->diff_w, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, chegvx_obj->diff_b, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, chegvx_obj->diff_z, LAPACKE_EIG_THRESHOLD);
}

TEST_F(chegvx_test, chegvx2) {
    EXPECT_NEAR(0.0, chegvx_obj->diff_a, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, chegvx_obj->diff_w, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, chegvx_obj->diff_b, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, chegvx_obj->diff_z, LAPACKE_EIG_THRESHOLD);
}

TEST_F(chegvx_test, chegvx3) {
    EXPECT_NEAR(0.0, chegvx_obj->diff_a, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, chegvx_obj->diff_w, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, chegvx_obj->diff_b, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, chegvx_obj->diff_z, LAPACKE_EIG_THRESHOLD);
}

TEST_F(chegvx_test, chegvx4) {
    EXPECT_NEAR(0.0, chegvx_obj->diff_a, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, chegvx_obj->diff_w, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, chegvx_obj->diff_b, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, chegvx_obj->diff_z, LAPACKE_EIG_THRESHOLD);
}
/* Begin dcomplex_common_parameters  class definition */
class hegvx_dcomplex_parameters{

   public:
	int bufsize_a;
	int bufsize_b;
	int bufsize_z;
	void *hModule, *dModule;
	double diff_a;
	double diff_w;
	double diff_b;
	double diff_z;
	double vl, vu;
	double abstol;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int ldz, lda, ldb;
	lapack_int il, iu;
	lapack_int* ifail;
	lapack_complex_double* ap, *bp, *z;
	char uplo;
	char jobz;
	char range;
	lapack_int itype;
	/*Output Parameter*/
	double* w;
	lapack_int m;
	lapack_int* ifailref;
	lapack_complex_double *apref, *bpref, *zref;
	double* wref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      hegvx_dcomplex_parameters (int matrix_layout, lapack_int itype, char jobz, char range,  char uplo, lapack_int n, lapack_int m, lapack_int ldz);
      ~hegvx_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
hegvx_dcomplex_parameters:: hegvx_dcomplex_parameters (int matrix_layout_i, lapack_int itype_i, char jobz_i, char range_i, char uplo_i, lapack_int n_i, lapack_int m_i, lapack_int ldz_i)
{
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	ldz = ldz_i;
	uplo = uplo_i;
	jobz = jobz_i;
	itype = itype_i;
	range = range_i;
	iu = il =vl = vu = 0;
	
	if (jobz == 'U')
		jobz = 'N';
	
	if (itype == 4)
		itype = 1;
	
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
		vu = n+5;
	}
	
	#if LAPACKE_TEST_VERBOSE
	printf(" \n hegvx dcomplex: matrix_layout = %d, itype:%d, jobz:%c, range:%c, uplo:%c, n: %d , m:%d, ldz: %d \n", matrix_layout, itype, jobz, range, uplo, n, m, ldz);
	
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{
		ldz = n;
		bufsize_z = ldz*m;
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
    {
		ldz = m;
		if (range == 'V')
				ldz = n;
		bufsize_z = ldz*n;
    }else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	lda =ldb =n;
	bufsize_a =  lda*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&ap, &apref, bufsize_a);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&bp, &bpref, bufsize_a);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&z, &zref, bufsize_z);
	lapacke_gtest_alloc_double_buffer_pair(&w, &wref, n);
	lapacke_gtest_alloc_int_buffer_pair(&ifail, &ifailref,	n);
	
	if ((ap==NULL) || (apref==NULL) || 
		(w == NULL) || (wref == NULL) ||
		(bp == NULL) || (bpref == NULL) ||
		(z == NULL) || (zref == NULL) ||
		(ifail == NULL) || (ifailref == NULL)){
		EXPECT_FALSE( true) << "hegvx_double_parameters object: malloc error.";
		hegvx_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix(ap, apref, bufsize_a, n, uplo);
	//lapacke_gtest_init_dcomplex_buffer_pair_rand(ap, apref, bufsize_a);
	lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix(bp, bpref, bufsize_a, n, uplo);
} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
hegvx_dcomplex_parameters :: ~hegvx_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" hegvx_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   hegvx_free();
  

}
/*  Test fixture class definition */
class zhegvx_test  : public  ::testing::Test {
public:
   hegvx_dcomplex_parameters  *zhegvx_obj;
   void SetUp();
   void TearDown () { delete zhegvx_obj; }
};

void zhegvx_test::SetUp(){

    /* LAPACKE zhegvx prototype */
    typedef int (*Fptr_NL_LAPACKE_zhegvx) ( int matrix_layout, lapack_int itype, char jobz, char range,
										char uplo, lapack_int n, lapack_complex_double* ap, lapack_int lda,  lapack_complex_double* bp, lapack_int ldb, 
										double vl, double vu, lapack_int il, lapack_int iu, double abstol, lapack_int* m, double* w, lapack_complex_double* z,
										lapack_int ldz, lapack_int* ifail);

    Fptr_NL_LAPACKE_zhegvx zhegvx;

    zhegvx_obj = new hegvx_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
							eig_paramslist[idx].itype,
							eig_paramslist[idx].vect_rd,
							eig_paramslist[idx].compz,
							eig_paramslist[idx].uplo,
							eig_paramslist[idx].m,
							eig_paramslist[idx].n,
							eig_paramslist[idx].m);

    idx = Circular_Increment_Index(idx);

    zhegvx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zhegvx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zhegvx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zhegvx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zhegvx = (Fptr_NL_LAPACKE_zhegvx)dlsym(zhegvx_obj->hModule, "LAPACKE_zhegvx");
    ASSERT_TRUE(zhegvx != NULL) << "failed to get the Netlib LAPACKE_zhegvx symbol";
    

    zhegvx_obj->inforef = zhegvx( zhegvx_obj->matrix_layout, zhegvx_obj->itype, zhegvx_obj->jobz, zhegvx_obj->range, zhegvx_obj->uplo, 
								zhegvx_obj->n, zhegvx_obj->apref, zhegvx_obj->lda, zhegvx_obj->bpref, zhegvx_obj->ldb, zhegvx_obj->vl, zhegvx_obj->vu, zhegvx_obj->il, zhegvx_obj->iu,
								zhegvx_obj->abstol, &zhegvx_obj->m, zhegvx_obj->wref, zhegvx_obj->zref, zhegvx_obj->ldz, zhegvx_obj->ifailref);

    /* Compute libflame's Lapacke o/p  */
    zhegvx_obj->info = LAPACKE_zhegvx( zhegvx_obj->matrix_layout, zhegvx_obj->itype, zhegvx_obj->jobz, zhegvx_obj->range, zhegvx_obj->uplo, 
										zhegvx_obj->n, zhegvx_obj->ap, zhegvx_obj->lda, zhegvx_obj->bp, zhegvx_obj->ldb, zhegvx_obj->vl, zhegvx_obj->vu, zhegvx_obj->il, zhegvx_obj->iu,
										zhegvx_obj->abstol, &zhegvx_obj->m, zhegvx_obj->w, zhegvx_obj->z, zhegvx_obj->ldz, zhegvx_obj->ifail);

    if( zhegvx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zhegvx is wrong\n", zhegvx_obj->info );
    }
    if( zhegvx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhegvx is wrong\n", 
        zhegvx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zhegvx_obj->diff_a =  computeDiff_z( zhegvx_obj->bufsize_a, zhegvx_obj->ap, zhegvx_obj->apref );
	zhegvx_obj->diff_b =  computeDiff_z( zhegvx_obj->bufsize_a, zhegvx_obj->bp, zhegvx_obj->bpref );
	zhegvx_obj->diff_z =  computeDiff_z( zhegvx_obj->bufsize_z, zhegvx_obj->z, zhegvx_obj->zref );
	zhegvx_obj->diff_w =  computeDiff_d( zhegvx_obj->n, zhegvx_obj->w, zhegvx_obj->wref );
}

TEST_F(zhegvx_test, zhegvx1) {
    EXPECT_NEAR(0.0, zhegvx_obj->diff_a, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zhegvx_obj->diff_w, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zhegvx_obj->diff_b, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zhegvx_obj->diff_z, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zhegvx_test, zhegvx2) {
    EXPECT_NEAR(0.0, zhegvx_obj->diff_a, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zhegvx_obj->diff_w, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zhegvx_obj->diff_b, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zhegvx_obj->diff_z, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zhegvx_test, zhegvx3) {
    EXPECT_NEAR(0.0, zhegvx_obj->diff_a, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zhegvx_obj->diff_w, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zhegvx_obj->diff_b, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zhegvx_obj->diff_z, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zhegvx_test, zhegvx4) {
    EXPECT_NEAR(0.0, zhegvx_obj->diff_a, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zhegvx_obj->diff_w, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zhegvx_obj->diff_b, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zhegvx_obj->diff_z, LAPACKE_EIG_THRESHOLD);
}
