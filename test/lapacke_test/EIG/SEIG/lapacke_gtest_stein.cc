#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"

#define stein_free() \
if (e!=NULL)    free(e); \
if (eref!=NULL) free(eref);\
if (d!=NULL)    free(d); \
if (dref!=NULL)    free(dref); \
if (w!=NULL)  free(w);\
if (wref!=NULL) free(wref); \
if (z!=NULL)  free(z);\
if (zref!=NULL) free(zref); \
if (ifailv!=NULL)  free(ifailv);\
if (ifailvref!=NULL) free(ifailvref); \
if (iblock!=NULL)  free(iblock);\
if (iblockref!=NULL)  free(iblockref);\
if (isplit!=NULL)  free(isplit);\
if (isplitref!=NULL)  free(isplitref);\
if( hModule != NULL) dlclose(hModule); \
if( dModule != NULL) dlclose(dModule);\
if( bModule != NULL) dlclose(bModule); \
if( lModule != NULL) dlclose(lModule);
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin float_common_parameters  class definition */
class stein_float_parameters{
	public:
		int bufsize_e, bufsize_d;
		int bufsize_w, bufsize_z;
		int bufsize_iblock, bufsize_isplit;
		int bufsize_ifailv;
		void *hModule, *dModule;
		void *lModule, *bModule;
		float diff_e;
		float diff_d;
		float diff_w;
		float diff_z;
		float diff_iblock , diff_isplit;
	/*input parameters */
		int matrix_layout;
		lapack_int n, ldz, m;
		lapack_int il, iu;
		lapack_int* isplit;
		lapack_int *iblock;
		float *e, *d;
		char order;
		char range;
		float vl, vu;
		float abstol;
	/*Output Parameter*/
		float* w;
		float* z, *zref;
		lapack_int m_stebz;
		lapack_int* isplitref;
		lapack_int *iblockref;
		float *eref, *dref;
		float* wref;
		lapack_int* ifailv, *ifailvref;
		lapack_int nsplit;
	/*Return Values*/	
		lapack_int info, inforef;
		lapack_int info_stebz, inforef_stebz;

   public:
      stein_float_parameters (int matrix_layout, char range, char order,  lapack_int n, lapack_int m);
      ~stein_float_parameters ();

};

/* Constructor definition  float_common_parameters */
stein_float_parameters:: stein_float_parameters (int matrix_layout_i, char range_i, char order_i,  lapack_int n_i, lapack_int m_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	range = range_i;
	order = order_i;
	iu = il = vl = vu = 0;
	abstol = 0.0;
	
	order = 'B'; //Order should 'B' for stebz 
	
	if (range == 'A')
	{	
	
	}else if (range == 'I')
	{
		iu = n/2;
		il = n/5;
		m_stebz = iu-il+1;		
	} else if (range == 'V')
	{	vl = n-5;
		vu= n+5;

	}

	#if LAPACKE_TEST_VERBOSE
		printf(" \n stein float: matrix_layout = %d, range:%c, order:%c, n: %d m: %d \n", matrix_layout, range, order,  n, m);
	#endif
	
	/*sizes*/
	bufsize_d = n;
	bufsize_e = n-1;
	bufsize_w = n;
	bufsize_iblock = bufsize_isplit = n;
	bufsize_ifailv = m;
	
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
	
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&e, &eref, bufsize_e);
	lapacke_gtest_alloc_float_buffer_pair(&d, &dref, bufsize_d);
	lapacke_gtest_alloc_int_buffer_pair(&iblock, &iblockref, bufsize_iblock);
	lapacke_gtest_alloc_float_buffer_pair(&w, &wref, bufsize_w);
	lapacke_gtest_alloc_int_buffer_pair(&isplit, &isplitref, bufsize_isplit);
	lapacke_gtest_alloc_int_buffer_pair(&ifailv, &ifailvref, bufsize_ifailv);
	lapacke_gtest_alloc_float_buffer_pair(&z, &zref, bufsize_z);
	
	if ((e==NULL) || (eref==NULL) ||\
		(w == NULL) || (wref == NULL) ||\
		(iblock == NULL) || (iblockref == NULL) ||\
		(d == NULL) || (dref == NULL) ||\
		(isplit == NULL) || (isplitref == NULL)||\
		(ifailv == NULL) || (ifailvref == NULL) ||\
		(z == NULL) || (zref == NULL)){
		EXPECT_FALSE( true) << "stebz_float_parameters object: malloc error.";
		stein_free();
		exit(0);
	}
	
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand(e, eref, bufsize_e);
	lapacke_gtest_init_float_buffer_pair_rand(d, dref, bufsize_d);
	lapacke_gtest_init_int_buffer_pair_with_constant(iblock, iblockref, bufsize_iblock, 1);
	lapacke_gtest_init_int_buffer_pair_with_constant(isplit, isplitref, bufsize_isplit, n);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
stein_float_parameters :: ~stein_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" stein_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   stein_free();

}
/*  Test fixture class definition */
class sstein_test  : public  ::testing::Test {
public:
   stein_float_parameters  *sstein_obj;
   void SetUp();
   void TearDown () { delete sstein_obj; }
};

void sstein_test::SetUp(){

    /* LAPACKE cstebz prototype */
	typedef int (*Fptr_NL_LAPACKE_sstebz) (char range, char order, lapack_int n, float vl,\
                           float vu, lapack_int il, lapack_int iu, float abstol,\
                           const float* d, const float* e, lapack_int* m,\
                           lapack_int* nsplit, float* w, lapack_int* iblock,\
                           lapack_int* isplit  );
						   
	Fptr_NL_LAPACKE_sstebz sstebz;
	/* LAPACKE sstein prototype */
    typedef int (*Fptr_NL_LAPACKE_sstein) (int matrix_layout, lapack_int n, const float* d,\
                           const float* e, lapack_int m, const float* w,\
                           const lapack_int* iblock, const lapack_int* isplit,\
                           float* z, lapack_int ldz,\
                           lapack_int* ifailv);

    Fptr_NL_LAPACKE_sstein sstein;

    sstein_obj = new stein_float_parameters ( svd_paramslist[idx].matrix_layout,
                           svd_paramslist[idx].range_gesvdx,
                           svd_paramslist[idx].jobu,
						   svd_paramslist[idx].n,
						  eig_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);
	
	/*stebz function call*/
	sstein_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sstein_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sstein_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sstein_obj->lModule != NULL) << "Netlib lapacke handle NULL";
	
    sstebz = (Fptr_NL_LAPACKE_sstebz)dlsym(sstein_obj->lModule, "LAPACKE_sstebz");
    ASSERT_TRUE(sstebz != NULL) << "failed to get the Netlib LAPACKE_sstein symbol";
	
	sstein_obj->inforef_stebz = sstebz( sstein_obj->range, sstein_obj->order, sstein_obj->n, sstein_obj->vl, sstein_obj->vu, sstein_obj->il, sstein_obj->iu,\
	sstein_obj->abstol, (const float*)sstein_obj->dref , (const float*)sstein_obj->eref, &sstein_obj->m_stebz, &sstein_obj->nsplit, sstein_obj->wref,\
	sstein_obj->iblockref, sstein_obj->isplitref);

    /* Compute libflame's Lapacke o/p  */
    sstein_obj->info_stebz = LAPACKE_sstebz( sstein_obj->range, sstein_obj->order, sstein_obj->n, sstein_obj->vl, sstein_obj->vu, sstein_obj->il, sstein_obj->iu,\
	sstein_obj->abstol, (const float*)sstein_obj->d, (const float*)sstein_obj->e, &sstein_obj->m_stebz, &sstein_obj->nsplit, sstein_obj->w,\
	sstein_obj->iblock, sstein_obj->isplit);

    if( sstein_obj->info_stebz < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_sstein is wrong\n", sstein_obj->info_stebz );
    }
    if( sstein_obj->inforef_stebz < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sstein is wrong\n", 
        sstein_obj->inforef_stebz );
    }
	
	/*LAPACKE_sstein */

    sstein_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sstein_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sstein_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sstein_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    sstein = (Fptr_NL_LAPACKE_sstein)dlsym(sstein_obj->hModule, "LAPACKE_sstein");
    ASSERT_TRUE(sstein != NULL) << "failed to get the Netlib LAPACKE_sstein symbol";   

		

    sstein_obj->inforef = sstein( sstein_obj->matrix_layout, sstein_obj->n, (const float*)sstein_obj->dref, (const float*)sstein_obj->eref,\
	sstein_obj->m, (const float*)sstein_obj->wref, (const lapack_int*)sstein_obj->iblockref, (const lapack_int*)sstein_obj->isplitref,\
	sstein_obj->zref, sstein_obj->ldz, sstein_obj->ifailvref);

    /* Compute libflame's Lapacke o/p  */
    sstein_obj->info = LAPACKE_sstein( sstein_obj->matrix_layout, sstein_obj->n, (const float*)sstein_obj->d, (const float*)sstein_obj->e,\
	sstein_obj->m, (const float*)sstein_obj->w, (const lapack_int*)sstein_obj->iblock, (const lapack_int*)sstein_obj->isplit,\
	sstein_obj->z, sstein_obj->ldz, sstein_obj->ifailv);

    if( sstein_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_sstein is wrong\n", sstein_obj->info );
    }
    if( sstein_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sstein is wrong\n", 
        sstein_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    sstein_obj->diff_z =  computeDiff_s( sstein_obj->bufsize_z, sstein_obj->z, sstein_obj->zref);

}

TEST_F(sstein_test, sstein1) {
    EXPECT_NEAR(0.0, sstein_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sstein_test, sstein2) {
    EXPECT_NEAR(0.0, sstein_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sstein_test, sstein3) {
    EXPECT_NEAR(0.0, sstein_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sstein_test, sstein4) {
    EXPECT_NEAR(0.0, sstein_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
}

/* Begin double_common_parameters  class definition */
class stein_double_parameters{
	public:
		int bufsize_e, bufsize_d;
		int bufsize_w, bufsize_z;
		int bufsize_iblock, bufsize_isplit;
		int bufsize_ifailv;
		void *hModule, *dModule;
		void *lModule, *bModule;
		double diff_e;
		double diff_d;
		double diff_w;
		double diff_z;
		double diff_iblock , diff_isplit;
	/*input parameters */
		int matrix_layout;
		lapack_int n, ldz, m;
		lapack_int il, iu;
		lapack_int* isplit;
		lapack_int *iblock;
		double *e, *d;
		char order;
		char range;
		double vl, vu;
		double abstol;
	/*Output Parameter*/
		double* w;
		double* z, *zref;
		lapack_int m_stebz;
		lapack_int* isplitref;
		lapack_int *iblockref;
		double *eref, *dref;
		double* wref;
		lapack_int* ifailv, *ifailvref;
		lapack_int nsplit;
	/*Return Values*/	
		lapack_int info, inforef;
		lapack_int info_stebz, inforef_stebz;

   public:
      stein_double_parameters (int matrix_layout, char range, char order,  lapack_int n, lapack_int m);
      ~stein_double_parameters ();

};

/* Constructor definition  double_common_parameters */
stein_double_parameters:: stein_double_parameters (int matrix_layout_i, char range_i, char order_i,  lapack_int n_i, lapack_int m_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	range = range_i;
	order = order_i;
	iu = il = vl = vu = 0;
	abstol = 0.0;
	
	order = 'B'; //Order should 'B' for stebz 
	
		
	if (range == 'A')
	{	
	
	}else if (range == 'I')
	{
		iu = n/2;
		il = n/5;
		m_stebz = iu-il+1;		
	} else if (range == 'V')
	{	vl = n-5;
		vu= n+5;

	}

	#if LAPACKE_TEST_VERBOSE
		printf(" \n stein double: matrix_layout = %d, range:%c, order:5c, n: %d m: %d \n", matrix_layout, range, order,  n, m);
	#endif
	
	/*sizes*/
	bufsize_d = n;
	bufsize_e = n-1;
	bufsize_w = n;
	bufsize_iblock = bufsize_isplit = n;
	bufsize_ifailv = m;
	
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
	
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&e, &eref, bufsize_e);
	lapacke_gtest_alloc_double_buffer_pair(&d, &dref, bufsize_d);
	lapacke_gtest_alloc_int_buffer_pair(&iblock, &iblockref, bufsize_iblock);
	lapacke_gtest_alloc_double_buffer_pair(&w, &wref, bufsize_w);
	lapacke_gtest_alloc_int_buffer_pair(&isplit, &isplitref, bufsize_isplit);
	lapacke_gtest_alloc_int_buffer_pair(&ifailv, &ifailvref, bufsize_ifailv);
	lapacke_gtest_alloc_double_buffer_pair(&z, &zref, bufsize_z);
	
	if ((e==NULL) || (eref==NULL) ||\
		(w == NULL) || (wref == NULL) ||\
		(iblock == NULL) || (iblockref == NULL) ||\
		(d == NULL) || (dref == NULL) ||\
		(isplit == NULL) || (isplitref == NULL)||\
		(ifailv == NULL) || (ifailvref == NULL) ||\
		(z == NULL) || (zref == NULL)){
		EXPECT_FALSE( true) << "stebz_double_parameters object: malloc error.";
		stein_free();
		exit(0);
	}
	
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand(e, eref, bufsize_e);
	lapacke_gtest_init_double_buffer_pair_rand(d, dref, bufsize_d);
	lapacke_gtest_init_int_buffer_pair_with_constant(iblock, iblockref, bufsize_iblock, 1);
	lapacke_gtest_init_int_buffer_pair_with_constant(isplit, isplitref, bufsize_isplit, n);

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
stein_double_parameters :: ~stein_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" stein_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   stein_free();

}
/*  Test fixture class definition */
class dstein_test  : public  ::testing::Test {
public:
   stein_double_parameters  *dstein_obj;
   void SetUp();
   void TearDown () { delete dstein_obj; }
};

void dstein_test::SetUp(){

    /* LAPACKE cstebz prototype */
	typedef int (*Fptr_NL_LAPACKE_dstebz) (char range, char order, lapack_int n, double vl,\
                           double vu, lapack_int il, lapack_int iu, double abstol,\
                           const double* d, const double* e, lapack_int* m,\
                           lapack_int* nsplit, double* w, lapack_int* iblock,\
                           lapack_int* isplit  );
						   
	Fptr_NL_LAPACKE_dstebz dstebz;
	/* LAPACKE dstein prototype */
    typedef int (*Fptr_NL_LAPACKE_dstein) (int matrix_layout, lapack_int n, const double* d,\
                           const double* e, lapack_int m, const double* w,\
                           const lapack_int* iblock, const lapack_int* isplit,\
						   double* z, lapack_int ldz,\
                           lapack_int* ifailv);

    Fptr_NL_LAPACKE_dstein dstein;

    dstein_obj = new stein_double_parameters ( svd_paramslist[idx].matrix_layout,
                           svd_paramslist[idx].range_gesvdx,
                           svd_paramslist[idx].jobu,
						   svd_paramslist[idx].n,
						   eig_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);
	
	/*stebz function call*/
	dstein_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dstein_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dstein_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dstein_obj->lModule != NULL) << "Netlib lapacke handle NULL";
	
    dstebz = (Fptr_NL_LAPACKE_dstebz)dlsym(dstein_obj->lModule, "LAPACKE_dstebz");
    ASSERT_TRUE(dstebz != NULL) << "failed to get the Netlib LAPACKE_dstein symbol";
	
	dstein_obj->inforef_stebz = dstebz( dstein_obj->range, dstein_obj->order, dstein_obj->n, dstein_obj->vl, dstein_obj->vu, dstein_obj->il, dstein_obj->iu,\
	dstein_obj->abstol, (const double*)dstein_obj->dref , (const double*)dstein_obj->eref, &dstein_obj->m_stebz, &dstein_obj->nsplit, dstein_obj->wref,\
	dstein_obj->iblockref, dstein_obj->isplitref);

    /* Compute libflame's Lapacke o/p  */
    dstein_obj->info_stebz = LAPACKE_dstebz( dstein_obj->range, dstein_obj->order, dstein_obj->n, dstein_obj->vl, dstein_obj->vu, dstein_obj->il, dstein_obj->iu,\
	dstein_obj->abstol, (const double*)dstein_obj->d, (const double*)dstein_obj->e, &dstein_obj->m_stebz, &dstein_obj->nsplit, dstein_obj->w,\
	dstein_obj->iblock, dstein_obj->isplit);

    if( dstein_obj->info_stebz < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dstein is wrong\n", dstein_obj->info_stebz );
    }
    if( dstein_obj->inforef_stebz < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dstein is wrong\n", 
        dstein_obj->inforef_stebz );
    }
	
	/*LAPACKE_dstein */

    dstein_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dstein_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dstein_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dstein_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    dstein = (Fptr_NL_LAPACKE_dstein)dlsym(dstein_obj->hModule, "LAPACKE_dstein");
    ASSERT_TRUE(dstein != NULL) << "failed to get the Netlib LAPACKE_dstein symbol";

    

    dstein_obj->inforef = dstein( dstein_obj->matrix_layout, dstein_obj->n, (const double*)dstein_obj->dref, (const double*)dstein_obj->eref,\
	dstein_obj->m, (const double*)dstein_obj->wref, (const lapack_int*)dstein_obj->iblockref, (const lapack_int*)dstein_obj->isplitref,\
	dstein_obj->zref, dstein_obj->ldz, dstein_obj->ifailvref);

    /* Compute libflame's Lapacke o/p  */
    dstein_obj->info = LAPACKE_dstein( dstein_obj->matrix_layout, dstein_obj->n, (const double*)dstein_obj->d, (const double*)dstein_obj->e,\
	dstein_obj->m, (const double*)dstein_obj->w, (const lapack_int*)dstein_obj->iblock, (const lapack_int*)dstein_obj->isplit,\
	dstein_obj->z, dstein_obj->ldz, dstein_obj->ifailv);

    if( dstein_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dstein is wrong\n", dstein_obj->info );
    }
    if( dstein_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dstein is wrong\n", 
        dstein_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dstein_obj->diff_z =  computeDiff_d( dstein_obj->bufsize_z, dstein_obj->z, dstein_obj->zref);

}

TEST_F(dstein_test, dstein1) {
    EXPECT_NEAR(0.0, dstein_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dstein_test, dstein2) {
    EXPECT_NEAR(0.0, dstein_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dstein_test, dstein3) {
    EXPECT_NEAR(0.0, dstein_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dstein_test, dstein4) {
    EXPECT_NEAR(0.0, dstein_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
}

/* Begin scomplex_common_parameters  class definition */
class stein_scomplex_parameters{
	public:
		int bufsize_e, bufsize_d;
		int bufsize_w, bufsize_z;
		int bufsize_iblock, bufsize_isplit;
		int bufsize_ifailv;
		void *hModule, *dModule;
		void *lModule, *bModule;
		float diff_e;
		float diff_d;
		float diff_w;
		float diff_z;
		float diff_iblock , diff_isplit;
	/*input parameters */
		int matrix_layout;
		lapack_int n, ldz, m;
		lapack_int il, iu;
		lapack_int* isplit;
		lapack_int *iblock;
		float *e, *d;
		char order;
		char range;
		float vl, vu;
		float abstol;
	/*Output Parameter*/
		float* w;
		lapack_complex_float* z, *zref;
		lapack_int m_stebz;
		lapack_int* isplitref;
		lapack_int *iblockref;
		float *eref, *dref;
		float* wref;
		lapack_int* ifailv, *ifailvref;
		lapack_int nsplit;
	/*Return Values*/	
		lapack_int info, inforef;
		lapack_int info_stebz, inforef_stebz;

   public:
      stein_scomplex_parameters (int matrix_layout, char range, char order,  lapack_int n, lapack_int m);
      ~stein_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
stein_scomplex_parameters:: stein_scomplex_parameters (int matrix_layout_i, char range_i, char order_i,  lapack_int n_i, lapack_int m_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	range = range_i;
	order = order_i;
	iu = il = vl = vu = 0;
	abstol = 0.0;
	
	order = 'B'; //Order should 'B' for stebz 
		
	if (range == 'A')
	{	
	
	}else if (range == 'I')
	{
		iu = n/2;
		il = n/5;
		m_stebz = iu-il+1;		
	} else if (range == 'V')
	{	vl = n-5;
		vu= n+5;

	}

	#if LAPACKE_TEST_VERBOSE
		printf(" \n stein scomplex: matrix_layout = %d, range:%c, order:%c, n: %d m: %d \n", matrix_layout, range, order,  n, m);
	#endif
	
	/*sizes*/
	bufsize_d = n;
	bufsize_e = n-1;
	bufsize_w = n;
	bufsize_iblock = bufsize_isplit = n;
	bufsize_ifailv = m;
	
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
	
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&e, &eref, bufsize_e);
	lapacke_gtest_alloc_float_buffer_pair(&d, &dref, bufsize_d);
	lapacke_gtest_alloc_int_buffer_pair(&iblock, &iblockref, bufsize_iblock);
	lapacke_gtest_alloc_float_buffer_pair(&w, &wref, bufsize_w);
	lapacke_gtest_alloc_int_buffer_pair(&isplit, &isplitref, bufsize_isplit);
	lapacke_gtest_alloc_int_buffer_pair(&ifailv, &ifailvref, bufsize_ifailv);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&z, &zref, bufsize_z);
	
	if ((e==NULL) || (eref==NULL) ||\
		(w == NULL) || (wref == NULL) ||\
		(iblock == NULL) || (iblockref == NULL) ||\
		(d == NULL) || (dref == NULL) ||\
		(isplit == NULL) || (isplitref == NULL)||\
		(ifailv == NULL) || (ifailvref == NULL) ||\
		(z == NULL) || (zref == NULL)){
		EXPECT_FALSE( true) << "stebz_float_parameters object: malloc error.";
		stein_free();
		exit(0);
	}
	
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand(e, eref, bufsize_e);
	lapacke_gtest_init_float_buffer_pair_rand(d, dref, bufsize_d);
	lapacke_gtest_init_int_buffer_pair_with_constant(iblock, iblockref, bufsize_iblock, 1);
	lapacke_gtest_init_int_buffer_pair_with_constant(isplit, isplitref, bufsize_isplit, n);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
stein_scomplex_parameters :: ~stein_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" stein_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   stein_free();

}
/*  Test fixture class definition */
class cstein_test  : public  ::testing::Test {
public:
   stein_scomplex_parameters  *cstein_obj;
   void SetUp();
   void TearDown () { delete cstein_obj; }
};

void cstein_test::SetUp(){

    /* LAPACKE cstebz prototype */
	typedef int (*Fptr_NL_LAPACKE_sstebz) (char range, char order, lapack_int n, float vl,\
                           float vu, lapack_int il, lapack_int iu, float abstol,\
                           const float* d, const float* e, lapack_int* m,\
                           lapack_int* nsplit, float* w, lapack_int* iblock,\
                           lapack_int* isplit  );
						   
	Fptr_NL_LAPACKE_sstebz sstebz;
	/* LAPACKE cstein prototype */
    typedef int (*Fptr_NL_LAPACKE_cstein) (int matrix_layout, lapack_int n, const float* d,\
                           const float* e, lapack_int m, const float* w,\
                           const lapack_int* iblock, const lapack_int* isplit,\
                           lapack_complex_float* z, lapack_int ldz,\
                           lapack_int* ifailv);

    Fptr_NL_LAPACKE_cstein cstein;

    cstein_obj = new stein_scomplex_parameters ( svd_paramslist[idx].matrix_layout,
                           svd_paramslist[idx].range_gesvdx,
                           svd_paramslist[idx].jobu,
						   svd_paramslist[idx].n,
						  eig_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);
	
	/*stebz function call*/
	cstein_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cstein_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cstein_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cstein_obj->lModule != NULL) << "Netlib lapacke handle NULL";
	
    sstebz = (Fptr_NL_LAPACKE_sstebz)dlsym(cstein_obj->lModule, "LAPACKE_sstebz");
    ASSERT_TRUE(sstebz != NULL) << "failed to get the Netlib LAPACKE_cstein symbol";
	
	cstein_obj->inforef_stebz = sstebz( cstein_obj->range, cstein_obj->order, cstein_obj->n, cstein_obj->vl, cstein_obj->vu, cstein_obj->il, cstein_obj->iu,\
	cstein_obj->abstol, (const float*)cstein_obj->dref , (const float*)cstein_obj->eref, &cstein_obj->m_stebz, &cstein_obj->nsplit, cstein_obj->wref,\
	cstein_obj->iblockref, cstein_obj->isplitref);

    /* Compute libflame's Lapacke o/p  */
    cstein_obj->info_stebz = LAPACKE_sstebz( cstein_obj->range, cstein_obj->order, cstein_obj->n, cstein_obj->vl, cstein_obj->vu, cstein_obj->il, cstein_obj->iu,\
	cstein_obj->abstol, (const float*)cstein_obj->d, (const float*)cstein_obj->e, &cstein_obj->m_stebz, &cstein_obj->nsplit, cstein_obj->w,\
	cstein_obj->iblock, cstein_obj->isplit);

    if( cstein_obj->info_stebz < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cstein is wrong\n", cstein_obj->info_stebz );
    }
    if( cstein_obj->inforef_stebz < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cstein is wrong\n", 
        cstein_obj->inforef_stebz );
    }
	
	/*LAPACKE_cstein */

    cstein_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cstein_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cstein_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cstein_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    cstein = (Fptr_NL_LAPACKE_cstein)dlsym(cstein_obj->hModule, "LAPACKE_cstein");
    ASSERT_TRUE(cstein != NULL) << "failed to get the Netlib LAPACKE_cstein symbol";   

		

    cstein_obj->inforef = cstein( cstein_obj->matrix_layout, cstein_obj->n, (const float*)cstein_obj->dref, (const float*)cstein_obj->eref,\
	cstein_obj->m, (const float*)cstein_obj->wref, (const lapack_int*)cstein_obj->iblockref, (const lapack_int*)cstein_obj->isplitref,\
	cstein_obj->zref, cstein_obj->ldz, cstein_obj->ifailvref);

    /* Compute libflame's Lapacke o/p  */
    cstein_obj->info = LAPACKE_cstein( cstein_obj->matrix_layout, cstein_obj->n, (const float*)cstein_obj->d, (const float*)cstein_obj->e,\
	cstein_obj->m, (const float*)cstein_obj->w, (const lapack_int*)cstein_obj->iblock, (const lapack_int*)cstein_obj->isplit,\
	cstein_obj->z, cstein_obj->ldz, cstein_obj->ifailv);

    if( cstein_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cstein is wrong\n", cstein_obj->info );
    }
    if( cstein_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cstein is wrong\n", 
        cstein_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cstein_obj->diff_z =  computeDiff_c( cstein_obj->bufsize_z, cstein_obj->z, cstein_obj->zref);

}

TEST_F(cstein_test, cstein1) {
    EXPECT_NEAR(0.0, cstein_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cstein_test, cstein2) {
    EXPECT_NEAR(0.0, cstein_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cstein_test, cstein3) {
    EXPECT_NEAR(0.0, cstein_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cstein_test, cstein4) {
    EXPECT_NEAR(0.0, cstein_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class stein_dcomplex_parameters{
	public:
		int bufsize_e, bufsize_d;
		int bufsize_w, bufsize_z;
		int bufsize_iblock, bufsize_isplit;
		int bufsize_ifailv;
		void *hModule, *dModule;
		void *lModule, *bModule;
		double diff_e;
		double diff_d;
		double diff_w;
		double diff_z;
		double diff_iblock , diff_isplit;
	/*input parameters */
		int matrix_layout;
		lapack_int n, ldz, m;
		lapack_int il, iu;
		lapack_int* isplit;
		lapack_int *iblock;
		double *e, *d;
		char order;
		char range;
		double vl, vu;
		double abstol;
	/*Output Parameter*/
		double* w;
		lapack_complex_double* z, *zref;
		lapack_int m_stebz;
		lapack_int* isplitref;
		lapack_int *iblockref;
		double *eref, *dref;
		double* wref;
		lapack_int* ifailv, *ifailvref;
		lapack_int nsplit;
	/*Return Values*/	
		lapack_int info, inforef;
		lapack_int info_stebz, inforef_stebz;

   public:
      stein_dcomplex_parameters (int matrix_layout, char range, char order,  lapack_int n, lapack_int m);
      ~stein_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
stein_dcomplex_parameters:: stein_dcomplex_parameters (int matrix_layout_i, char range_i, char order_i,  lapack_int n_i, lapack_int m_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	range = range_i;
	order = order_i;
	iu = il = vl = vu = 0;
	abstol = 0.0;
	
	order = 'B'; //Order should 'B' for stebz 

	if (range == 'A')
	{	
	
	}else if (range == 'I')
	{
		iu = n/2;
		il = n/5;
		m_stebz = iu-il+1;		
	} else if (range == 'V')
	{	vl = n-5;
		vu= n+5;

	}

	#if LAPACKE_TEST_VERBOSE
		printf(" \n stein dcomplex: matrix_layout = %d, range:%c, order:5c, n: %d m: %d \n", matrix_layout, range, order,  n, m);
	#endif
	
	/*sizes*/
	bufsize_d = n;
	bufsize_e = n-1;
	bufsize_w = n;
	bufsize_iblock = bufsize_isplit = n;
	bufsize_ifailv = m;
	
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
	
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&e, &eref, bufsize_e);
	lapacke_gtest_alloc_double_buffer_pair(&d, &dref, bufsize_d);
	lapacke_gtest_alloc_int_buffer_pair(&iblock, &iblockref, bufsize_iblock);
	lapacke_gtest_alloc_double_buffer_pair(&w, &wref, bufsize_w);
	lapacke_gtest_alloc_int_buffer_pair(&isplit, &isplitref, bufsize_isplit);
	lapacke_gtest_alloc_int_buffer_pair(&ifailv, &ifailvref, bufsize_ifailv);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&z, &zref, bufsize_z);
	
	if ((e==NULL) || (eref==NULL) ||\
		(w == NULL) || (wref == NULL) ||\
		(iblock == NULL) || (iblockref == NULL) ||\
		(d == NULL) || (dref == NULL) ||\
		(isplit == NULL) || (isplitref == NULL)||\
		(ifailv == NULL) || (ifailvref == NULL) ||\
		(z == NULL) || (zref == NULL)){
		EXPECT_FALSE( true) << "stebz_double_parameters object: malloc error.";
		stein_free();
		exit(0);
	}
	
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand(e, eref, bufsize_e);
	lapacke_gtest_init_double_buffer_pair_rand(d, dref, bufsize_d);
	lapacke_gtest_init_int_buffer_pair_with_constant(iblock, iblockref, bufsize_iblock, 1);
	lapacke_gtest_init_int_buffer_pair_with_constant(isplit, isplitref, bufsize_isplit, n);

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
stein_dcomplex_parameters :: ~stein_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" stein_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   stein_free();

}
/*  Test fixture class definition */
class zstein_test  : public  ::testing::Test {
public:
   stein_dcomplex_parameters  *zstein_obj;
   void SetUp();
   void TearDown () { delete zstein_obj; }
};

void zstein_test::SetUp(){

    /* LAPACKE cstebz prototype */
	typedef int (*Fptr_NL_LAPACKE_dstebz) (char range, char order, lapack_int n, double vl,\
                           double vu, lapack_int il, lapack_int iu, double abstol,\
                           const double* d, const double* e, lapack_int* m,\
                           lapack_int* nsplit, double* w, lapack_int* iblock,\
                           lapack_int* isplit  );
						   
	Fptr_NL_LAPACKE_dstebz dstebz;
	/* LAPACKE zstein prototype */
    typedef int (*Fptr_NL_LAPACKE_zstein) (int matrix_layout, lapack_int n, const double* d,\
                           const double* e, lapack_int m, const double* w,\
                           const lapack_int* iblock, const lapack_int* isplit,\
                           lapack_complex_double* z, lapack_int ldz,\
                           lapack_int* ifailv);

    Fptr_NL_LAPACKE_zstein zstein;

    zstein_obj = new stein_dcomplex_parameters ( svd_paramslist[idx].matrix_layout,
                           svd_paramslist[idx].range_gesvdx,
                           svd_paramslist[idx].jobu,
						   svd_paramslist[idx].n,
						   eig_paramslist[idx].nrhs);
						   //svd_paramslist[idx].m);

    idx = Circular_Increment_Index(idx);
	
	/*stebz function call*/
	zstein_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zstein_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zstein_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zstein_obj->lModule != NULL) << "Netlib lapacke handle NULL";
	
    dstebz = (Fptr_NL_LAPACKE_dstebz)dlsym(zstein_obj->lModule, "LAPACKE_dstebz");
    ASSERT_TRUE(dstebz != NULL) << "failed to get the Netlib LAPACKE_zstein symbol";
	
	zstein_obj->inforef_stebz = dstebz( zstein_obj->range, zstein_obj->order, zstein_obj->n, zstein_obj->vl, zstein_obj->vu, zstein_obj->il, zstein_obj->iu,\
	zstein_obj->abstol, (const double*)zstein_obj->dref , (const double*)zstein_obj->eref, &zstein_obj->m_stebz, &zstein_obj->nsplit, zstein_obj->wref,\
	zstein_obj->iblockref, zstein_obj->isplitref);

    /* Compute libflame's Lapacke o/p  */
    zstein_obj->info_stebz = LAPACKE_dstebz( zstein_obj->range, zstein_obj->order, zstein_obj->n, zstein_obj->vl, zstein_obj->vu, zstein_obj->il, zstein_obj->iu,\
	zstein_obj->abstol, (const double*)zstein_obj->d, (const double*)zstein_obj->e, &zstein_obj->m_stebz, &zstein_obj->nsplit, zstein_obj->w,\
	zstein_obj->iblock, zstein_obj->isplit);

    if( zstein_obj->info_stebz < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zstein is wrong\n", zstein_obj->info_stebz );
    }
    if( zstein_obj->inforef_stebz < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zstein is wrong\n", 
        zstein_obj->inforef_stebz );
    }
	
	/*LAPACKE_zstein */

    zstein_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zstein_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zstein_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zstein_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zstein = (Fptr_NL_LAPACKE_zstein)dlsym(zstein_obj->hModule, "LAPACKE_zstein");
    ASSERT_TRUE(zstein != NULL) << "failed to get the Netlib LAPACKE_zstein symbol";

    

    zstein_obj->inforef = zstein( zstein_obj->matrix_layout, zstein_obj->n, (const double*)zstein_obj->dref, (const double*)zstein_obj->eref,\
	zstein_obj->m, (const double*)zstein_obj->wref, (const lapack_int*)zstein_obj->iblockref, (const lapack_int*)zstein_obj->isplitref,\
	zstein_obj->zref, zstein_obj->ldz, zstein_obj->ifailvref);

    /* Compute libflame's Lapacke o/p  */
    zstein_obj->info = LAPACKE_zstein( zstein_obj->matrix_layout, zstein_obj->n, (const double*)zstein_obj->d, (const double*)zstein_obj->e,\
	zstein_obj->m, (const double*)zstein_obj->w, (const lapack_int*)zstein_obj->iblock, (const lapack_int*)zstein_obj->isplit,\
	zstein_obj->z, zstein_obj->ldz, zstein_obj->ifailv);

    if( zstein_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zstein is wrong\n", zstein_obj->info );
    }
    if( zstein_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zstein is wrong\n", 
        zstein_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zstein_obj->diff_z =  computeDiff_z( zstein_obj->bufsize_z, zstein_obj->z, zstein_obj->zref);

}

TEST_F(zstein_test, zstein1) {
    EXPECT_NEAR(0.0, zstein_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zstein_test, zstein2) {
    EXPECT_NEAR(0.0, zstein_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zstein_test, zstein3) {
    EXPECT_NEAR(0.0, zstein_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zstein_test, zstein4) {
    EXPECT_NEAR(0.0, zstein_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
}