#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"

#define lacn2_free() \
if (x!=NULL)    free(x); \
if (xref!=NULL) free(xref);\
if (v!=NULL)  free(v);\
if (vref!=NULL) free(vref); \
if (est!=NULL)  free(est);\
if (estref!=NULL)  free(estref);\
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule);

	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin float_common_parameters  class definition */
class lacn2_float_parameters{

   public:
   int bufsize_isave;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	lapack_int n;
	float* x;
	lapack_int kase_act, kase_init;
	lapack_int isave[3] = {0};
	lapack_int* isign;
	/*Output Parameter*/
	float* v;	
	float *vref, *xref;
	float *est, *estref;
	lapack_int* isignref;
	lapack_int isaveref[3] = {0};
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      lacn2_float_parameters (lapack_int n, lapack_int kase_act);
      ~lacn2_float_parameters ();

};

/* Constructor definition  float_common_parameters */
lacn2_float_parameters::lacn2_float_parameters (lapack_int n_i, lapack_int kase_act_i)
{	
	n = n_i;
	kase_act = kase_act_i;
	kase_init = 0; // On initial call to routing kase must be 0
	bufsize_isave = 3;
	
	if(kase_act == 3)
		kase_act = 1;

	#if LAPACKE_TEST_VERBOSE
		printf(" \n lacn2 float:  n:%d, kase_act:%d\n",  n, kase_act);
	#endif
	
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&x, &xref, n);
	lapacke_gtest_alloc_float_buffer_pair(&v, &vref, n);
	lapacke_gtest_alloc_float_buffer_pair(&est, &estref, n);
	lapacke_gtest_alloc_int_buffer_pair(&isign, &isignref, n);
	
	if ((x == NULL) || (xref == NULL) ||\
		(est == NULL) || (estref == NULL) ||\
		(isign == NULL) || (isignref == NULL) ||\
		(v == NULL) || (vref == NULL)){
		EXPECT_FALSE( true) << "lacn2_float_parameters object: malloc error.";
		lacn2_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand( x, xref, n);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
lacn2_float_parameters :: ~lacn2_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" lacn2_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lacn2_free();
	if (isign != NULL)    free(isign); 
	if (isignref != NULL) free(isignref);

}
/*  Test fixture class definition */
class slacn2_test  : public  ::testing::Test {
public:
   lacn2_float_parameters  *slacn2_obj;
   void SetUp();
   void TearDown () { delete slacn2_obj; }
};

void slacn2_test::SetUp(){

    /* LAPACKE slacn2 prototype */
    typedef int (*Fptr_NL_LAPACKE_slacn2) ( lapack_int n, float* v, float* x, lapack_int* isgn,
											float* est, lapack_int* kase, lapack_int* isave );

    Fptr_NL_LAPACKE_slacn2 slacn2;

    slacn2_obj = new lacn2_float_parameters (eig_paramslist[idx].n,
						eig_paramslist[idx].itype);

    idx = Circular_Increment_Index(idx);

    slacn2_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    slacn2_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(slacn2_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(slacn2_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    slacn2 = (Fptr_NL_LAPACKE_slacn2)dlsym(slacn2_obj->hModule, "LAPACKE_slacn2");
    ASSERT_TRUE(slacn2 != NULL) << "failed to get the Netlib LAPACKE_slacn2 symbol";
    
	/* Workspace quesry with kase = 0 */
    slacn2_obj->inforef = slacn2(slacn2_obj->n, slacn2_obj->vref, slacn2_obj->xref,  slacn2_obj->isignref,
								slacn2_obj->estref,	&(slacn2_obj->kase_init), slacn2_obj->isaveref);

    /* Compute libflame's Lapacke o/p  */
    slacn2_obj->info = LAPACKE_slacn2(slacn2_obj->n, slacn2_obj->v, slacn2_obj->x,  slacn2_obj->isign,
									slacn2_obj->est, &(slacn2_obj->kase_init), slacn2_obj->isave);

    if( slacn2_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_slacn2 init is wrong\n", slacn2_obj->info );
    }
    if( slacn2_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_slacn2  init is wrong\n", 
        slacn2_obj->inforef );
    }
	/*Actual immediate return call with kase = 1 or 2*/
	slacn2_obj->inforef = slacn2(slacn2_obj->n, slacn2_obj->vref, slacn2_obj->xref, slacn2_obj->isignref,
								slacn2_obj->estref,	&(slacn2_obj->kase_act), slacn2_obj->isaveref);

    /* Compute libflame's Lapacke o/p  */
    slacn2_obj->info = LAPACKE_slacn2(slacn2_obj->n, slacn2_obj->v, slacn2_obj->x, slacn2_obj->isignref,
										slacn2_obj->est, &(slacn2_obj->kase_act), slacn2_obj->isave);

    if( slacn2_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_slacn2 final is wrong\n", slacn2_obj->info );
    }
    if( slacn2_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_slacn2 final is wrong\n", 
        slacn2_obj->inforef );
    }
	
    /* Compute Difference between libflame and Netlib o/ps  */
    slacn2_obj->diff =  computeDiff_s( slacn2_obj->n, 
                slacn2_obj->v, slacn2_obj->vref );

}

TEST_F(slacn2_test, slacn21) {
    EXPECT_NEAR(0.0, slacn2_obj->diff, LAPACKE_GTEST_THRESHOLD);
	
}

TEST_F(slacn2_test, slacn22) {
    EXPECT_NEAR(0.0, slacn2_obj->diff, LAPACKE_GTEST_THRESHOLD);
	
}

TEST_F(slacn2_test, slacn23) {
    EXPECT_NEAR(0.0, slacn2_obj->diff, LAPACKE_GTEST_THRESHOLD);
	
}

TEST_F(slacn2_test, slacn24) {
    EXPECT_NEAR(0.0, slacn2_obj->diff, LAPACKE_GTEST_THRESHOLD);
	
}

/* Begin double_common_parameters  class definition */
class lacn2_double_parameters{

   public:
   int bufsize_isave;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	lapack_int n;
	double* x;
	lapack_int kase_act, kase_init;
	lapack_int isave[3] = {0};
	lapack_int* isign;
	/*Output Parameter*/
	double* v;	
	double *vref, *xref;
	double *est, *estref;
	lapack_int* isignref;
	lapack_int isaveref[3] = {0};
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      lacn2_double_parameters (lapack_int n, lapack_int kase_act);
      ~lacn2_double_parameters ();

};

/* Constructor definition  double_common_parameters */
lacn2_double_parameters::lacn2_double_parameters (lapack_int n_i, lapack_int kase_act_i)
{	
	n = n_i;
	kase_act = kase_act_i;
	kase_init = 0; // On initial call to routing kase must be 0
	bufsize_isave = 3;
	
	if(kase_act == 3)
		kase_act = 1;

	#if LAPACKE_TEST_VERBOSE
		printf(" \n lacn2 double:  n:%d, kase_act:%d\n",  n, kase_act);
	#endif
	
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&x, &xref, n);
	lapacke_gtest_alloc_double_buffer_pair(&v, &vref, n);
	lapacke_gtest_alloc_double_buffer_pair(&est, &estref, n);
	lapacke_gtest_alloc_int_buffer_pair(&isign, &isignref, n);
	
	if ((x == NULL) || (xref == NULL) ||\
		(est == NULL) || (estref == NULL) ||\
		(isign == NULL) || (isignref == NULL) ||\
		(v == NULL) || (vref == NULL)){
		EXPECT_FALSE( true) << "lacn2_double_parameters object: malloc error.";
		lacn2_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand( x, xref, n);

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
lacn2_double_parameters :: ~lacn2_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" lacn2_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lacn2_free();
	if (isign != NULL)    free(isign); 
	if (isignref != NULL) free(isignref);

}
/*  Test fixture class definition */
class dlacn2_test  : public  ::testing::Test {
public:
   lacn2_double_parameters  *dlacn2_obj;
   void SetUp();
   void TearDown () { delete dlacn2_obj; }
};

void dlacn2_test::SetUp(){

    /* LAPACKE dlacn2 prototype */
    typedef int (*Fptr_NL_LAPACKE_dlacn2) ( lapack_int n, double* v, double* x, lapack_int* isgn,
											double* est, lapack_int* kase, lapack_int* isave );

    Fptr_NL_LAPACKE_dlacn2 dlacn2;

    dlacn2_obj = new lacn2_double_parameters (eig_paramslist[idx].n,
						eig_paramslist[idx].itype);

    idx = Circular_Increment_Index(idx);

    dlacn2_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dlacn2_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dlacn2_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dlacn2_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    dlacn2 = (Fptr_NL_LAPACKE_dlacn2)dlsym(dlacn2_obj->hModule, "LAPACKE_dlacn2");
    ASSERT_TRUE(dlacn2 != NULL) << "failed to get the Netlib LAPACKE_dlacn2 symbol";
    
	/* Workspace quesry with kase = 0 */
    dlacn2_obj->inforef = dlacn2(dlacn2_obj->n, dlacn2_obj->vref, dlacn2_obj->xref,  dlacn2_obj->isignref,
								dlacn2_obj->estref,	&(dlacn2_obj->kase_init), dlacn2_obj->isaveref);

    /* Compute libflame's Lapacke o/p  */
    dlacn2_obj->info = LAPACKE_dlacn2(dlacn2_obj->n, dlacn2_obj->v, dlacn2_obj->x,  dlacn2_obj->isign,
									dlacn2_obj->est, &(dlacn2_obj->kase_init), dlacn2_obj->isave);

    if( dlacn2_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dlacn2 init is wrong\n", dlacn2_obj->info );
    }
    if( dlacn2_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dlacn2  init is wrong\n", 
        dlacn2_obj->inforef );
    }
	/*Actual immediate return call with kase = 1 or 2*/
	dlacn2_obj->inforef = dlacn2(dlacn2_obj->n, dlacn2_obj->vref, dlacn2_obj->xref, dlacn2_obj->isignref,
								dlacn2_obj->estref,	&(dlacn2_obj->kase_act), dlacn2_obj->isaveref);

    /* Compute libflame's Lapacke o/p  */
    dlacn2_obj->info = LAPACKE_dlacn2(dlacn2_obj->n, dlacn2_obj->v, dlacn2_obj->x, dlacn2_obj->isignref,
										dlacn2_obj->est, &(dlacn2_obj->kase_act), dlacn2_obj->isave);

    if( dlacn2_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dlacn2 final is wrong\n", dlacn2_obj->info );
    }
    if( dlacn2_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dlacn2 final is wrong\n", 
        dlacn2_obj->inforef );
    }
	
    /* Compute Difference between libflame and Netlib o/ps  */
    dlacn2_obj->diff =  computeDiff_d( dlacn2_obj->n, 
                dlacn2_obj->v, dlacn2_obj->vref );

}

TEST_F(dlacn2_test, dlacn21) {
    EXPECT_NEAR(0.0, dlacn2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dlacn2_test, dlacn22) {
    EXPECT_NEAR(0.0, dlacn2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dlacn2_test, dlacn23) {
    EXPECT_NEAR(0.0, dlacn2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dlacn2_test, dlacn24) {
    EXPECT_NEAR(0.0, dlacn2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin scomplex_common_parameters  class definition */
class lacn2_scomplex_parameters{

   public:
   int bufsize_isave;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	lapack_int n;
	lapack_complex_float* x;
	lapack_int kase_act, kase_init;
	lapack_int isave[3] = {0};
	/*Output Parameter*/
	lapack_complex_float* v;	
	lapack_complex_float *vref, *xref;
	float *est, *estref;
	lapack_int isaveref[3] = {0};
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      lacn2_scomplex_parameters (lapack_int n, lapack_int kase_act);
      ~lacn2_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
lacn2_scomplex_parameters::lacn2_scomplex_parameters (lapack_int n_i, lapack_int kase_act_i)
{	
	n = n_i;
	kase_act = kase_act_i;
	kase_init = 0; // On initial call to routing kase must be 0
	bufsize_isave = 3;
	
	if(kase_act == 3)
		kase_act = 1;

	#if LAPACKE_TEST_VERBOSE
		printf(" \n lacn2 scomplex:  n:%d, kase_act:%d\n",  n, kase_act);
	#endif
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&x, &xref, n);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&v, &vref, n);
	lapacke_gtest_alloc_float_buffer_pair(&est, &estref, n);
	
	
	if ((x == NULL) || (xref == NULL) ||\
		(est == NULL) || (estref == NULL) ||\
		(v == NULL) || (vref == NULL)){
		EXPECT_FALSE( true) << "lacn2_float_parameters object: malloc error.";
		lacn2_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( x, xref, n);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
lacn2_scomplex_parameters :: ~lacn2_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" lacn2_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lacn2_free();

}
/*  Test fixture class definition */
class clacn2_test  : public  ::testing::Test {
public:
   lacn2_scomplex_parameters  *clacn2_obj;
   void SetUp();
   void TearDown () { delete clacn2_obj; }
};

void clacn2_test::SetUp(){

    /* LAPACKE clacn2 prototype */
    typedef int (*Fptr_NL_LAPACKE_clacn2) (lapack_int n, lapack_complex_float* v,
											lapack_complex_float* x,
											float* est, lapack_int* kase, lapack_int* isave);

    Fptr_NL_LAPACKE_clacn2 clacn2;

    clacn2_obj = new lacn2_scomplex_parameters (eig_paramslist[idx].n,
						eig_paramslist[idx].itype);

    idx = Circular_Increment_Index(idx);

    clacn2_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    clacn2_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(clacn2_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(clacn2_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    clacn2 = (Fptr_NL_LAPACKE_clacn2)dlsym(clacn2_obj->hModule, "LAPACKE_clacn2");
    ASSERT_TRUE(clacn2 != NULL) << "failed to get the Netlib LAPACKE_clacn2 symbol";
    
	/* Workspace quesry with kase = 0 */
    clacn2_obj->inforef = clacn2(clacn2_obj->n, clacn2_obj->vref, clacn2_obj->xref, clacn2_obj->estref,
								&(clacn2_obj->kase_init), clacn2_obj->isaveref);

    /* Compute libflame's Lapacke o/p  */
    clacn2_obj->info = LAPACKE_clacn2(clacn2_obj->n, clacn2_obj->v, clacn2_obj->x, clacn2_obj->est,
										&(clacn2_obj->kase_init), clacn2_obj->isave);

    if( clacn2_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_clacn2 init is wrong\n", clacn2_obj->info );
    }
    if( clacn2_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_clacn2  init is wrong\n", 
        clacn2_obj->inforef );
    }
	/*Actual immediate return call with kase = 1 or 2*/
	clacn2_obj->inforef = clacn2(clacn2_obj->n, clacn2_obj->vref, clacn2_obj->xref, clacn2_obj->estref,
								&(clacn2_obj->kase_act), clacn2_obj->isaveref);

    /* Compute libflame's Lapacke o/p  */
    clacn2_obj->info = LAPACKE_clacn2(clacn2_obj->n, clacn2_obj->v, clacn2_obj->x, clacn2_obj->est, 
										&(clacn2_obj->kase_act), clacn2_obj->isave);

    if( clacn2_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_clacn2 final is wrong\n", clacn2_obj->info );
    }
    if( clacn2_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_clacn2 final is wrong\n", 
        clacn2_obj->inforef );
    }
	
    /* Compute Difference between libflame and Netlib o/ps  */
    clacn2_obj->diff =  computeDiff_c( clacn2_obj->n, 
                clacn2_obj->v, clacn2_obj->vref );

}

TEST_F(clacn2_test, clacn21) {
    EXPECT_NEAR(0.0, clacn2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(clacn2_test, clacn22) {
    EXPECT_NEAR(0.0, clacn2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(clacn2_test, clacn23) {
    EXPECT_NEAR(0.0, clacn2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(clacn2_test, clacn24) {
    EXPECT_NEAR(0.0, clacn2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class lacn2_dcomplex_parameters{

   public:
   int bufsize_isave;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	lapack_int n;
	lapack_complex_double* x;
	lapack_int kase_act, kase_init;
	lapack_int isave[3] = {0};
	/*Output Parameter*/
	lapack_complex_double* v;	
	lapack_complex_double *vref, *xref;
	double *est, *estref;
	lapack_int isaveref[3] = {0};
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      lacn2_dcomplex_parameters (lapack_int n, lapack_int kase_act);
      ~lacn2_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
lacn2_dcomplex_parameters::lacn2_dcomplex_parameters (lapack_int n_i, lapack_int kase_act_i)
{	
	n = n_i;
	kase_act = kase_act_i;
	kase_init = 0; // On initial call to routing kase must be 0
	bufsize_isave = 3;
	
	if(kase_act == 3)
		kase_act = 1;

	#if LAPACKE_TEST_VERBOSE
		printf(" \n lacn2 dcomplex:  n:%d, kase_act:%d\n",  n, kase_act);
	#endif
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&x, &xref, n);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&v, &vref, n);
	lapacke_gtest_alloc_double_buffer_pair(&est, &estref, n);
	
	
	if ((x == NULL) || (xref == NULL) ||\
		(est == NULL) || (estref == NULL) ||\
		(v == NULL) || (vref == NULL)){
		EXPECT_FALSE( true) << "lacn2_double_parameters object: malloc error.";
		lacn2_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( x, xref, n);

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
lacn2_dcomplex_parameters :: ~lacn2_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" lacn2_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lacn2_free();

}
/*  Test fixture class definition */
class zlacn2_test  : public  ::testing::Test {
public:
   lacn2_dcomplex_parameters  *zlacn2_obj;
   void SetUp();
   void TearDown () { delete zlacn2_obj; }
};

void zlacn2_test::SetUp(){

    /* LAPACKE zlacn2 prototype */
    typedef int (*Fptr_NL_LAPACKE_zlacn2) (lapack_int n, lapack_complex_double* v,
											lapack_complex_double* x,
											double* est, lapack_int* kase, lapack_int* isave);

    Fptr_NL_LAPACKE_zlacn2 zlacn2;

    zlacn2_obj = new lacn2_dcomplex_parameters (eig_paramslist[idx].n,
						eig_paramslist[idx].itype);

    idx = Circular_Increment_Index(idx);

    zlacn2_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zlacn2_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zlacn2_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zlacn2_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zlacn2 = (Fptr_NL_LAPACKE_zlacn2)dlsym(zlacn2_obj->hModule, "LAPACKE_zlacn2");
    ASSERT_TRUE(zlacn2 != NULL) << "failed to get the Netlib LAPACKE_zlacn2 symbol";
    
	/* Workspace quesry with kase = 0 */
    zlacn2_obj->inforef = zlacn2(zlacn2_obj->n, zlacn2_obj->vref, zlacn2_obj->xref, zlacn2_obj->estref,
								&(zlacn2_obj->kase_init), zlacn2_obj->isaveref);

    /* Compute libflame's Lapacke o/p  */
    zlacn2_obj->info = LAPACKE_zlacn2(zlacn2_obj->n, zlacn2_obj->v, zlacn2_obj->x, zlacn2_obj->est,
										&(zlacn2_obj->kase_init), zlacn2_obj->isave);

    if( zlacn2_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zlacn2 init is wrong\n", zlacn2_obj->info );
    }
    if( zlacn2_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zlacn2  init is wrong\n", 
        zlacn2_obj->inforef );
    }
	/*Actual immediate return call with kase = 1 or 2*/
	zlacn2_obj->inforef = zlacn2(zlacn2_obj->n, zlacn2_obj->vref, zlacn2_obj->xref, zlacn2_obj->estref,
								&(zlacn2_obj->kase_act), zlacn2_obj->isaveref);

    /* Compute libflame's Lapacke o/p  */
    zlacn2_obj->info = LAPACKE_zlacn2(zlacn2_obj->n, zlacn2_obj->v, zlacn2_obj->x, zlacn2_obj->est, 
										&(zlacn2_obj->kase_act), zlacn2_obj->isave);

    if( zlacn2_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zlacn2 final is wrong\n", zlacn2_obj->info );
    }
    if( zlacn2_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zlacn2 final is wrong\n", 
        zlacn2_obj->inforef );
    }
	
    /* Compute Difference between libflame and Netlib o/ps  */
    zlacn2_obj->diff =  computeDiff_z( zlacn2_obj->n, 
                zlacn2_obj->v, zlacn2_obj->vref );

}

TEST_F(zlacn2_test, zlacn21) {
    EXPECT_NEAR(0.0, zlacn2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlacn2_test, zlacn22) {
    EXPECT_NEAR(0.0, zlacn2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlacn2_test, zlacn23) {
    EXPECT_NEAR(0.0, zlacn2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlacn2_test, zlacn24) {
    EXPECT_NEAR(0.0, zlacn2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}