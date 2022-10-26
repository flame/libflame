#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define stebz_free() \
if (e!=NULL)    free(e); \
if (eref!=NULL) free(eref);\
if (d!=NULL)    free(d); \
if (dref!=NULL)    free(dref); \
if (w!=NULL)  free(w);\
if (wref!=NULL) free(wref); \
if (iblock!=NULL)  free(iblock);\
if (iblockref!=NULL)  free(iblockref);\
if (isplit!=NULL)  free(isplit);\
if (isplitref!=NULL)  free(isplitref);\
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin float_common_parameters  class definition */
class stebz_float_parameters{

   public:
	int bufsize_e, bufsize_d;
	int bufsize_w;
	int bufsize_iblock, bufsize_isplit;
	void *hModule, *dModule;
	float diff_e;
	float diff_d;
	float diff_w;
	float diff_iblock , diff_isplit;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
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
	lapack_int m;
	lapack_int* isplitref;
	lapack_int *iblockref;
	float *eref, *dref;
	float* wref;
	lapack_int nsplit;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      stebz_float_parameters (int matrix_layout, char range, char order,  lapack_int n);
      ~stebz_float_parameters ();

};

/* Constructor definition  float_common_parameters */
stebz_float_parameters:: stebz_float_parameters (int matrix_layout_i, char range_i, char order_i,  lapack_int n_i)
{
	matrix_layout = matrix_layout_i;
	n = n_i;
	order = order_i;
	range = range_i;
	iu = il = vl = vu = 0;
	abstol = 0.0;
	
	
	/*ORder is either 'B' or 'E'*/
	if (order == 'U')
		order = 'B';
	else 
		order = 'E';
		
	if (range == 'A')
	{	
	
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
	printf(" \n stebz float: matrix_layout = %d,  order :%c, range:%c,  n: %d , iu:%d, Il:%d \n", matrix_layout, order, range,  n, iu, il);
	#endif
	/*sizes*/
	bufsize_d = n;
	bufsize_e = n-1;
	bufsize_w = n;
	bufsize_iblock = bufsize_isplit = n;
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&e, &eref, bufsize_e);
	lapacke_gtest_alloc_float_buffer_pair(&d, &dref, bufsize_d);
	lapacke_gtest_alloc_int_buffer_pair(&iblock, &iblockref, bufsize_iblock);
	lapacke_gtest_alloc_float_buffer_pair(&w, &wref, bufsize_w);
	lapacke_gtest_alloc_int_buffer_pair(&isplit, &isplitref, bufsize_isplit);
	
	if ((e==NULL) || (eref==NULL) || 
		(w == NULL) || (wref == NULL) ||
		(iblock == NULL) || (iblockref == NULL) ||
		(d == NULL) || (dref == NULL) ||
		(isplit == NULL) || (isplitref == NULL)){
		EXPECT_FALSE( true) << "stebz_float_parameters object: malloc error.";
		stebz_free();
		exit(0);
	}
	/* Initialization of input matrices */
	

	lapacke_gtest_init_float_buffer_pair_rand(e, eref, bufsize_e);
	lapacke_gtest_init_float_buffer_pair_rand(d, dref, bufsize_e);
} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
stebz_float_parameters :: ~stebz_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" stebz_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   stebz_free();
  

}
/*  Test fixture class definition */
class sstebz_test  : public  ::testing::Test {
public:
   stebz_float_parameters  *sstebz_obj;
   void SetUp();
   void TearDown () { delete sstebz_obj; }
};

void sstebz_test::SetUp(){

    /* LAPACKE sstebz prototype */
    typedef int (*Fptr_NL_LAPACKE_sstebz) (char range, char order, lapack_int n, float vl,\
                           float vu, lapack_int il, lapack_int iu, float abstol,\
                           const float* d, const float* e, lapack_int* m,\
                           lapack_int* nsplit, float* w, lapack_int* iblock,\
                           lapack_int* isplit );

    Fptr_NL_LAPACKE_sstebz sstebz;

    sstebz_obj = new stebz_float_parameters ( svd_paramslist[idx].matrix_layout,
							svd_paramslist[idx].range_gesvdx,
							svd_paramslist[idx].jobu,
							svd_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    sstebz_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sstebz_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sstebz_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sstebz_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    sstebz = (Fptr_NL_LAPACKE_sstebz)dlsym(sstebz_obj->hModule, "LAPACKE_sstebz");
    ASSERT_TRUE(sstebz != NULL) << "failed to get the Netlib LAPACKE_sstebz symbol";
    

    sstebz_obj->inforef = sstebz( sstebz_obj->range, sstebz_obj->order, sstebz_obj->n, sstebz_obj->vl, sstebz_obj->vu, sstebz_obj->il, sstebz_obj->iu,\
	sstebz_obj->abstol, (const float*)sstebz_obj->dref , (const float*)sstebz_obj->eref, &sstebz_obj->m, &sstebz_obj->nsplit, sstebz_obj->wref,\
	sstebz_obj->iblockref, sstebz_obj->isplitref);

    /* Compute libflame's Lapacke o/p  */
    sstebz_obj->info = LAPACKE_sstebz( sstebz_obj->range, sstebz_obj->order, sstebz_obj->n, sstebz_obj->vl, sstebz_obj->vu, sstebz_obj->il, sstebz_obj->iu,\
	sstebz_obj->abstol, (const float*)sstebz_obj->d , (const float*)sstebz_obj->e, &sstebz_obj->m, &sstebz_obj->nsplit, sstebz_obj->w, sstebz_obj->iblock, sstebz_obj->isplit);

    if( sstebz_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_sstebz is wrong\n", sstebz_obj->info );
    }
    if( sstebz_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sstebz is wrong\n", 
        sstebz_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
  	sstebz_obj->diff_isplit =  computeDiff_i( sstebz_obj->bufsize_isplit, sstebz_obj->isplit, sstebz_obj->isplitref );
	sstebz_obj->diff_iblock =  computeDiff_i( sstebz_obj->bufsize_iblock, sstebz_obj->iblock, sstebz_obj->iblockref );
	sstebz_obj->diff_w =  computeDiff_s( sstebz_obj->bufsize_w, sstebz_obj->w, sstebz_obj->wref );
}

TEST_F(sstebz_test, sstebz1) {
	//EXPECT_NEAR(0.0, sstebz_obj->diff_isplit, LAPACKE_GTEST_THRESHOLD);
	//EXPECT_NEAR(0.0, sstebz_obj->diff_iblock, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, sstebz_obj->diff_w, LAPACKE_GTEST_THRESHOLD);

}

TEST_F(sstebz_test, sstebz2) {
    //EXPECT_NEAR(0.0, sstebz_obj->diff_isplit, LAPACKE_GTEST_THRESHOLD);
	//EXPECT_NEAR(0.0, sstebz_obj->diff_iblock, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, sstebz_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sstebz_test, sstebz3) {
    //EXPECT_NEAR(0.0, sstebz_obj->diff_isplit, LAPACKE_GTEST_THRESHOLD);
//	EXPECT_NEAR(0.0, sstebz_obj->diff_iblock, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, sstebz_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sstebz_test, sstebz4) {
	//EXPECT_NEAR(0.0, sstebz_obj->diff_isplit, LAPACKE_GTEST_THRESHOLD);
//	EXPECT_NEAR(0.0, sstebz_obj->diff_iblock, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, sstebz_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
	
}

/* Begin double_common_parameters  class definition */
class stebz_double_parameters{

   public:
	int bufsize_e, bufsize_d;
	int bufsize_w;
	int bufsize_iblock, bufsize_isplit;
	void *hModule, *dModule;
	double diff_e;
	double diff_d;
	double diff_w;
	double diff_iblock , diff_isplit;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
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
	lapack_int m;
	lapack_int* isplitref;
	lapack_int *iblockref;
	double *eref, *dref;
	double* wref;
	lapack_int nsplit;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      stebz_double_parameters (int matrix_layout, char range, char order,  lapack_int n);
      ~stebz_double_parameters ();

};

/* Constructor definition  double_common_parameters */
stebz_double_parameters:: stebz_double_parameters (int matrix_layout_i, char range_i, char order_i,  lapack_int n_i)
{
	matrix_layout = matrix_layout_i;
	n = n_i;
	order = order_i;
	range = range_i;
	iu = il = vl = vu = 0;
	abstol = 0.0;
	
	
	/*ORder is either 'B' or 'E'*/
	if (order == 'U')
		order = 'B';
	else 
		order = 'E';
		
	if (range == 'A')
	{	
	
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
	printf(" \n stebz double: matrix_layout = %d,  order :%c, range:%c,  n: %d , iu:%d, Il:%d \n", matrix_layout, order, range,  n, iu, il);
	#endif
	/*sizes*/
	bufsize_d = n;
	bufsize_e = n-1;
	bufsize_w = n;
	bufsize_iblock = bufsize_isplit = n;
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&e, &eref, bufsize_e);
	lapacke_gtest_alloc_double_buffer_pair(&d, &dref, bufsize_d);
	lapacke_gtest_alloc_int_buffer_pair(&iblock, &iblockref, bufsize_iblock);
	lapacke_gtest_alloc_double_buffer_pair(&w, &wref, bufsize_w);
	lapacke_gtest_alloc_int_buffer_pair(&isplit, &isplitref, bufsize_isplit);
	
	if ((e==NULL) || (eref==NULL) || 
		(w == NULL) || (wref == NULL) ||
		(iblock == NULL) || (iblockref == NULL) ||
		(d == NULL) || (dref == NULL) ||
		(isplit == NULL) || (isplitref == NULL)){
		EXPECT_FALSE( true) << "stebz_double_parameters object: malloc error.";
		stebz_free();
		exit(0);
	}
	/* Initialization of input matrices */
	

	lapacke_gtest_init_double_buffer_pair_rand(e, eref, bufsize_e);
	lapacke_gtest_init_double_buffer_pair_rand(d, dref, bufsize_e);
} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
stebz_double_parameters :: ~stebz_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" stebz_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   stebz_free();
  

}
/*  Test fixture class definition */
class dstebz_test  : public  ::testing::Test {
public:
   stebz_double_parameters  *dstebz_obj;
   void SetUp();
   void TearDown () { delete dstebz_obj; }
};

void dstebz_test::SetUp(){

    /* LAPACKE dstebz prototype */
    typedef int (*Fptr_NL_LAPACKE_dstebz) (char range, char order, lapack_int n, double vl,\
                           double vu, lapack_int il, lapack_int iu, double abstol,\
                           const double* d, const double* e, lapack_int* m,\
                           lapack_int* nsplit, double* w, lapack_int* iblock,\
                           lapack_int* isplit );

    Fptr_NL_LAPACKE_dstebz dstebz;

    dstebz_obj = new stebz_double_parameters ( svd_paramslist[idx].matrix_layout,
							svd_paramslist[idx].range_gesvdx,
							svd_paramslist[idx].jobu,
							svd_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    dstebz_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dstebz_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dstebz_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dstebz_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    dstebz = (Fptr_NL_LAPACKE_dstebz)dlsym(dstebz_obj->hModule, "LAPACKE_dstebz");
    ASSERT_TRUE(dstebz != NULL) << "failed to get the Netlib LAPACKE_dstebz symbol";
    

    dstebz_obj->inforef = dstebz( dstebz_obj->range, dstebz_obj->order, dstebz_obj->n, dstebz_obj->vl, dstebz_obj->vu, dstebz_obj->il,\
dstebz_obj->iu,	dstebz_obj->abstol, (const double*)dstebz_obj->dref , (const double*)dstebz_obj->eref, &dstebz_obj->m, &dstebz_obj->nsplit,\
	dstebz_obj->wref, dstebz_obj->iblockref, dstebz_obj->isplitref);

    /* Compute libflame's Lapacke o/p  */
    dstebz_obj->info = LAPACKE_dstebz(dstebz_obj->range, dstebz_obj->order, dstebz_obj->n, dstebz_obj->vl, dstebz_obj->vu, dstebz_obj->il,\
	dstebz_obj->iu,	dstebz_obj->abstol, (const double*)dstebz_obj->d , (const double*)dstebz_obj->e,&dstebz_obj->m, &dstebz_obj->nsplit,\
	dstebz_obj->w, dstebz_obj->iblock, dstebz_obj->isplit);

    if( dstebz_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dstebz is wrong\n", dstebz_obj->info );
    }
    if( dstebz_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dstebz is wrong\n", 
        dstebz_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
  	dstebz_obj->diff_isplit =  computeDiff_i( dstebz_obj->bufsize_isplit, dstebz_obj->isplit, dstebz_obj->isplitref );
	dstebz_obj->diff_iblock =  computeDiff_i( dstebz_obj->bufsize_iblock, dstebz_obj->iblock, dstebz_obj->iblockref );
	dstebz_obj->diff_w =  computeDiff_d( dstebz_obj->bufsize_w, dstebz_obj->w, dstebz_obj->wref );
}

TEST_F(dstebz_test, dstebz1) {
	//EXPECT_NEAR(0.0, dstebz_obj->diff_isplit, LAPACKE_GTEST_THRESHOLD);
	//EXPECT_NEAR(0.0, dstebz_obj->diff_iblock, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, dstebz_obj->diff_w, LAPACKE_GTEST_THRESHOLD);

}

TEST_F(dstebz_test, dstebz2) {
//    EXPECT_NEAR(0.0, dstebz_obj->diff_isplit, LAPACKE_GTEST_THRESHOLD);
//	EXPECT_NEAR(0.0, dstebz_obj->diff_iblock, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, dstebz_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dstebz_test, dstebz3) {
  //  EXPECT_NEAR(0.0, dstebz_obj->diff_isplit, LAPACKE_GTEST_THRESHOLD);
	//EXPECT_NEAR(0.0, dstebz_obj->diff_iblock, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, dstebz_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dstebz_test, dstebz4) {
	//EXPECT_NEAR(0.0, dstebz_obj->diff_isplit, LAPACKE_GTEST_THRESHOLD);
//	EXPECT_NEAR(0.0, dstebz_obj->diff_iblock, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, dstebz_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
	
}
