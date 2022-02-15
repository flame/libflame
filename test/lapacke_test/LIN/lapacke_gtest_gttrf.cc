#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"

#define gttrf_free() \
       if (ipiv != NULL) free (ipiv); \
       if (ipivref != NULL) free (ipivref); \
       if (d != NULL)    free (d   ); \
       if (dref != NULL) free (dref); \
       if (du != NULL)    free (du   ); \
       if (duref != NULL) free (duref); \
       if (du2 != NULL)    free (du2   ); \
       if (du2ref != NULL) free (du2ref); \
       if (dl != NULL)    free (dl   ); \
       if (dlref != NULL) free (dlref)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin gttrf_double_parameters  class definition */
class gttrf_double_parameters{

   public:
      /* Input parameters */
      lapack_int n; // No of rows,Columns
    
      /* Input/ Output parameters */
      double *dl,*dlref;// subdiagonal elements of A
      double *d,*dref; // diagonal elements of A
      double *du,*duref; // superdiagonal elements of A

      /* Output parameters */
      double *du2,*du2ref; //
      lapack_int *ipiv,*ipivref; // The pivot indices

      /* Return Values */
      lapack_int info,inforef;

   public: 
      gttrf_double_parameters ( lapack_int n_i);
      ~gttrf_double_parameters (); 
};  /* end of gttrf_double_parameters  class definition */

/* Constructor gttrf_double_parameters definition */
gttrf_double_parameters:: gttrf_double_parameters (lapack_int n_i) {

    n = n_i;
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &d, &dref, n);
    lapacke_gtest_alloc_double_buffer_pair( &dl, &dlref, n-1);
    lapacke_gtest_alloc_double_buffer_pair( &du, &duref, n-1);
    lapacke_gtest_alloc_double_buffer_pair( &du2, &du2ref, n-2);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv,&ipivref,n);

    if((d==NULL) || (dref==NULL) ||  \
       (dl==NULL) || (dlref==NULL) || \
       (du==NULL) || (duref==NULL) ||  \
       (du2==NULL) || (du2ref==NULL) || \
       (ipiv==NULL) || (ipivref==NULL)){
       gttrf_free();
       printf(" gttrf_double_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( d,dref,n);
    lapacke_gtest_init_double_buffer_pair_rand( dl,dlref,n-1);
    lapacke_gtest_init_double_buffer_pair_rand( du,duref,n-1);
    lapacke_gtest_init_double_buffer_pair_rand( du2,du2ref,n-2);

   } /* end of Constructor  */

gttrf_double_parameters:: ~gttrf_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gttrf_double_parameters object: destructor invoked. \n");
#endif
   gttrf_free();
}

//  Test fixture class definition
class dgttrf_test  : public  ::testing::Test {
public:
   gttrf_double_parameters  *dgttrf_obj;
   void SetUp();   
   void TearDown () { delete dgttrf_obj; }
};


void dgttrf_test::SetUp(){

    dgttrf_obj = new  gttrf_double_parameters(lin_solver_paramslist[idx].n);
    idx = Circular_Increment_Index(idx);
    /* LAPACKE dgttrf prototype */
    typedef int (*Fptr_NL_LAPACKE_dgttrf) ( lapack_int n,double *dl,
         double *d,double *du,double *du2,lapack_int *ipiv );

    Fptr_NL_LAPACKE_dgttrf DGTTRF;
    void *hModule,*dModule;


    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit( -1 );
    }

    DGTTRF = (Fptr_NL_LAPACKE_dgttrf)dlsym(hModule,"LAPACKE_dgttrf");
    ASSERT_TRUE(DGTTRF != NULL) << "failed to get the Netlib LAPACKE_dgttrf symbol";

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    dgttrf_obj->inforef = DGTTRF( dgttrf_obj->n,dgttrf_obj->dlref,
                           dgttrf_obj->dref,dgttrf_obj->duref,
                           dgttrf_obj->du2ref,dgttrf_obj->ipivref);

   /* Compute libflame's Lapacke o/p  */
    dgttrf_obj->info     = LAPACKE_dgttrf( dgttrf_obj->n,dgttrf_obj->dl,
                           dgttrf_obj->d,dgttrf_obj->du,dgttrf_obj->du2,
                           dgttrf_obj->ipiv);

    if( dgttrf_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_dgttrf \
                  is wrong\n",dgttrf_obj->info );
    }
    if( dgttrf_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dgttrf is wrong\n",
        dgttrf_obj->inforef );
    }
    dlclose(hModule);
    dlclose(dModule);
}

TEST_F(dgttrf_test, dgttrf1) {
    double diff;
    int ipiv_diff;
    ipiv_diff = computeDiff_i( dgttrf_obj->n,dgttrf_obj->ipiv,dgttrf_obj->ipivref);
    EXPECT_EQ(ipiv_diff, 0) << " pivot computation in dgttrf1 test case failed ";
    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_d( dgttrf_obj->n,dgttrf_obj->d,dgttrf_obj->dref );
    diff +=  computeDiff_d( dgttrf_obj->n-1,dgttrf_obj->du,dgttrf_obj->duref );
    diff +=  computeDiff_d( dgttrf_obj->n-1,dgttrf_obj->dl,dgttrf_obj->dlref );
    diff +=  computeDiff_d( dgttrf_obj->n-2,dgttrf_obj->du2,dgttrf_obj->du2ref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgttrf_test, dgttrf2) {
    double diff;
    int ipiv_diff;
    ipiv_diff = computeDiff_i( dgttrf_obj->n,dgttrf_obj->ipiv,dgttrf_obj->ipivref);
    EXPECT_EQ(ipiv_diff, 0) << " pivot computation in dgttrf1 test case failed ";
    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_d( dgttrf_obj->n,dgttrf_obj->d,dgttrf_obj->dref );
    diff +=  computeDiff_d( dgttrf_obj->n-1,dgttrf_obj->du,dgttrf_obj->duref );
    diff +=  computeDiff_d( dgttrf_obj->n-1,dgttrf_obj->dl,dgttrf_obj->dlref );
    diff +=  computeDiff_d( dgttrf_obj->n-2,dgttrf_obj->du2,dgttrf_obj->du2ref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

/* Begin gttrf_float_parameters  class definition */
class gttrf_float_parameters{

   public:
      /* Input parameters */
      lapack_int n; // No of rows,Columns
    
      /* Input/ Output parameters */
      float *dl,*dlref;// subdiagonal elements of A
      float *d,*dref; // diagonal elements of A
      float *du,*duref; // superdiagonal elements of A

      /* Output parameters */
      float *du2,*du2ref; //
      lapack_int *ipiv,*ipivref; // The pivot indices

      /* Return Values */
      lapack_int info,inforef;

   public: 
      gttrf_float_parameters ( lapack_int n_i);
      ~gttrf_float_parameters (); 
};  /* end of gttrf_float_parameters  class definition */

/* Constructor gttrf_float_parameters definition */
gttrf_float_parameters:: gttrf_float_parameters (lapack_int n_i) {

    n = n_i;
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &d, &dref, n);
    lapacke_gtest_alloc_float_buffer_pair( &dl, &dlref, n-1);
    lapacke_gtest_alloc_float_buffer_pair( &du, &duref, n-1);
    lapacke_gtest_alloc_float_buffer_pair( &du2, &du2ref, n-2);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv,&ipivref,n);

    if((d==NULL) || (dref==NULL) ||  \
       (dl==NULL) || (dlref==NULL) || \
       (du==NULL) || (duref==NULL) ||  \
       (du2==NULL) || (du2ref==NULL) || \
       (ipiv==NULL) || (ipivref==NULL)){
       gttrf_free();
       printf(" gttrf_float_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( d,dref,n);
    lapacke_gtest_init_float_buffer_pair_rand( dl,dlref,n-1);
    lapacke_gtest_init_float_buffer_pair_rand( du,duref,n-1);
    lapacke_gtest_init_float_buffer_pair_rand( du2,du2ref,n-2);

   } /* end of Constructor  */

gttrf_float_parameters:: ~gttrf_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gttrf_float_parameters object: destructor invoked. \n");
#endif
   gttrf_free();
}

//  Test fixture class definition
class sgttrf_test  : public  ::testing::Test {
public:
   gttrf_float_parameters  *sgttrf_obj;
   void SetUp();   
   void TearDown () { delete sgttrf_obj; }
};


void sgttrf_test::SetUp(){

    sgttrf_obj = new  gttrf_float_parameters(lin_solver_paramslist[idx].n);
    idx = Circular_Increment_Index(idx);
    /* LAPACKE sgttrf prototype */
    typedef int (*Fptr_NL_LAPACKE_sgttrf) ( lapack_int n,float *dl,
         float *d,float *du,float *du2,lapack_int *ipiv );

    Fptr_NL_LAPACKE_sgttrf SGTTRF;
    void *hModule,*dModule;


    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit( -1 );
    }

    SGTTRF = (Fptr_NL_LAPACKE_sgttrf)dlsym(hModule,"LAPACKE_sgttrf");
    ASSERT_TRUE(SGTTRF != NULL) << "failed to get the Netlib LAPACKE_sgttrf symbol";

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    sgttrf_obj->inforef = SGTTRF( sgttrf_obj->n,sgttrf_obj->dlref,
                           sgttrf_obj->dref,sgttrf_obj->duref,
                           sgttrf_obj->du2ref,sgttrf_obj->ipivref);

   /* Compute libflame's Lapacke o/p  */
    sgttrf_obj->info     = LAPACKE_sgttrf( sgttrf_obj->n,sgttrf_obj->dl,
                           sgttrf_obj->d,sgttrf_obj->du,sgttrf_obj->du2,
                           sgttrf_obj->ipiv);

    if( sgttrf_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_sgttrf \
                  is wrong\n",sgttrf_obj->info );
    }
    if( sgttrf_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sgttrf is wrong\n",
        sgttrf_obj->inforef );
    }
    dlclose(hModule);
    dlclose(dModule);
}

TEST_F(sgttrf_test, sgttrf1) {
    float diff;
    int ipiv_diff;
    ipiv_diff = computeDiff_i( sgttrf_obj->n,sgttrf_obj->ipiv,sgttrf_obj->ipivref);
    EXPECT_EQ(ipiv_diff, 0) << " pivot computation  failed ";
    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_s( sgttrf_obj->n,sgttrf_obj->d,sgttrf_obj->dref );
    diff +=  computeDiff_s( sgttrf_obj->n-1,sgttrf_obj->du,sgttrf_obj->duref );
    diff +=  computeDiff_s( sgttrf_obj->n-1,sgttrf_obj->dl,sgttrf_obj->dlref );
    diff +=  computeDiff_s( sgttrf_obj->n-2,sgttrf_obj->du2,sgttrf_obj->du2ref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgttrf_test, sgttrf2) {
    float diff;
    int ipiv_diff;
    ipiv_diff = computeDiff_i( sgttrf_obj->n,sgttrf_obj->ipiv,sgttrf_obj->ipivref);
    EXPECT_EQ(ipiv_diff, 0) << " pivot computation  failed ";
    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_s( sgttrf_obj->n,sgttrf_obj->d,sgttrf_obj->dref );
    diff +=  computeDiff_s( sgttrf_obj->n-1,sgttrf_obj->du,sgttrf_obj->duref );
    diff +=  computeDiff_s( sgttrf_obj->n-1,sgttrf_obj->dl,sgttrf_obj->dlref );
    diff +=  computeDiff_s( sgttrf_obj->n-2,sgttrf_obj->du2,sgttrf_obj->du2ref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

/* Begin gttrf_scomplex_parameters  class definition */
class gttrf_scomplex_parameters{

   public:
      /* Input parameters */
      lapack_int n; // No of rows,Columns
    
      /* Input/ Output parameters */
      lapack_complex_float *dl,*dlref;// subdiagonal elements of A
      lapack_complex_float *d,*dref; // diagonal elements of A
      lapack_complex_float *du,*duref; // superdiagonal elements of A

      /* Output parameters */
      lapack_complex_float *du2,*du2ref; //
      lapack_int *ipiv,*ipivref; // The pivot indices

      /* Return Values */
      lapack_int info,inforef;

   public: 
      gttrf_scomplex_parameters ( lapack_int n_i);
      ~gttrf_scomplex_parameters (); 
};  /* end of gttrf_scomplex_parameters  class definition */

/* Constructor gttrf_scomplex_parameters definition */
gttrf_scomplex_parameters:: gttrf_scomplex_parameters (lapack_int n_i) {

    n = n_i;
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &d, &dref, n);
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &dl, &dlref, n-1);
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &du, &duref, n-1);
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &du2, &du2ref, n-2);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv,&ipivref,n);

    if((d==NULL) || (dref==NULL) ||  \
       (dl==NULL) || (dlref==NULL) || \
       (du==NULL) || (duref==NULL) ||  \
       (du2==NULL) || (du2ref==NULL) || \
       (ipiv==NULL) || (ipivref==NULL)){
       gttrf_free();
       printf(" gttrf_scomplex_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( d,dref,n);
    lapacke_gtest_init_scomplex_buffer_pair_rand( dl,dlref,n-1);
    lapacke_gtest_init_scomplex_buffer_pair_rand( du,duref,n-1);
    lapacke_gtest_init_scomplex_buffer_pair_rand( du2,du2ref,n-2);

   } /* end of Constructor  */

gttrf_scomplex_parameters:: ~gttrf_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gttrf_scomplex_parameters object: destructor invoked. \n");
#endif
   gttrf_free();
}

//  Test fixture class definition
class cgttrf_test  : public  ::testing::Test {
public:
   gttrf_scomplex_parameters  *cgttrf_obj;
   void SetUp();   
   void TearDown () { delete cgttrf_obj; }
};


void cgttrf_test::SetUp(){

    cgttrf_obj = new  gttrf_scomplex_parameters(lin_solver_paramslist[idx].n);
    idx = Circular_Increment_Index(idx);
    /* LAPACKE cgttrf prototype */
    typedef int (*Fptr_NL_LAPACKE_cgttrf) ( lapack_int n,lapack_complex_float *dl,
                             lapack_complex_float *d,lapack_complex_float *du,
                                lapack_complex_float *du2,lapack_int *ipiv );

    Fptr_NL_LAPACKE_cgttrf CGTTRF;
    void *hModule,*dModule;


    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit( -1 );
    }

    CGTTRF = (Fptr_NL_LAPACKE_cgttrf)dlsym(hModule,"LAPACKE_cgttrf");
    ASSERT_TRUE(CGTTRF != NULL) << "failed to get the Netlib LAPACKE_cgttrf symbol";

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    cgttrf_obj->inforef = CGTTRF( cgttrf_obj->n,cgttrf_obj->dlref,
                           cgttrf_obj->dref,cgttrf_obj->duref,
                           cgttrf_obj->du2ref,cgttrf_obj->ipivref);

   /* Compute libflame's Lapacke o/p  */
    cgttrf_obj->info     = LAPACKE_cgttrf( cgttrf_obj->n,cgttrf_obj->dl,
                           cgttrf_obj->d,cgttrf_obj->du,cgttrf_obj->du2,
                           cgttrf_obj->ipiv);

    if( cgttrf_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_cgttrf \
                  is wrong\n",cgttrf_obj->info );
    }
    if( cgttrf_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cgttrf is wrong\n",
        cgttrf_obj->inforef );
    }
    dlclose(hModule);
    dlclose(dModule);
}

TEST_F(cgttrf_test, cgttrf1) {
    float diff;    
    int ipiv_diff;
    ipiv_diff = computeDiff_i( cgttrf_obj->n,cgttrf_obj->ipiv,cgttrf_obj->ipivref);
    EXPECT_EQ(ipiv_diff, 0) << " pivot computation  failed ";
    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_c( cgttrf_obj->n,cgttrf_obj->d,cgttrf_obj->dref );
    diff +=  computeDiff_c( cgttrf_obj->n-1,cgttrf_obj->du,cgttrf_obj->duref );
    diff +=  computeDiff_c( cgttrf_obj->n-1,cgttrf_obj->dl,cgttrf_obj->dlref );
    diff +=  computeDiff_c( cgttrf_obj->n-2,cgttrf_obj->du2,cgttrf_obj->du2ref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgttrf_test, cgttrf2) {
    float diff;    
    int ipiv_diff;
    ipiv_diff = computeDiff_i( cgttrf_obj->n,cgttrf_obj->ipiv,cgttrf_obj->ipivref);
    EXPECT_EQ(ipiv_diff, 0) << " pivot computation  failed ";
    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_c( cgttrf_obj->n,cgttrf_obj->d,cgttrf_obj->dref );
    diff +=  computeDiff_c( cgttrf_obj->n-1,cgttrf_obj->du,cgttrf_obj->duref );
    diff +=  computeDiff_c( cgttrf_obj->n-1,cgttrf_obj->dl,cgttrf_obj->dlref );
    diff +=  computeDiff_c( cgttrf_obj->n-2,cgttrf_obj->du2,cgttrf_obj->du2ref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

/* Begin gttrf_dcomplex_parameters  class definition */
class gttrf_dcomplex_parameters{

   public:
      /* Input parameters */
      lapack_int n; // No of rows,Columns
    
      /* Input/ Output parameters */
      lapack_complex_double *dl,*dlref;// subdiagonal elements of A
      lapack_complex_double *d,*dref; // diagonal elements of A
      lapack_complex_double *du,*duref; // superdiagonal elements of A

      /* Output parameters */
      lapack_complex_double *du2,*du2ref; //
      lapack_int *ipiv,*ipivref; // The pivot indices

      /* Return Values */
      lapack_int info,inforef;

   public: 
      gttrf_dcomplex_parameters ( lapack_int n_i);
      ~gttrf_dcomplex_parameters (); 
};  /* end of gttrf_dcomplex_parameters  class definition */

/* Constructor gttrf_dcomplex_parameters definition */
gttrf_dcomplex_parameters:: gttrf_dcomplex_parameters (lapack_int n_i) {

    n = n_i;
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &d, &dref, n);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &dl, &dlref, n-1);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &du, &duref, n-1);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &du2, &du2ref, n-2);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv,&ipivref,n);

    if((d==NULL) || (dref==NULL) ||  \
       (dl==NULL) || (dlref==NULL) || \
       (du==NULL) || (duref==NULL) ||  \
       (du2==NULL) || (du2ref==NULL) || \
       (ipiv==NULL) || (ipivref==NULL)){
       gttrf_free();
       printf(" gttrf_dcomplex_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( d,dref,n);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( dl,dlref,n-1);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( du,duref,n-1);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( du2,du2ref,n-2);

   } /* end of Constructor  */

gttrf_dcomplex_parameters:: ~gttrf_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gttrf_dcomplex_parameters object: destructor invoked. \n");
#endif
   gttrf_free();
}

//  Test fixture class definition
class zgttrf_test  : public  ::testing::Test {
public:
   gttrf_dcomplex_parameters  *zgttrf_obj;
   void SetUp();   
   void TearDown () { delete zgttrf_obj; }
};


void zgttrf_test::SetUp(){

    zgttrf_obj = new  gttrf_dcomplex_parameters(lin_solver_paramslist[idx].n);
    idx = Circular_Increment_Index(idx);
    /* LAPACKE zgttrf prototype */
    typedef int (*Fptr_NL_LAPACKE_zgttrf) ( lapack_int n,lapack_complex_double *dl,
         lapack_complex_double *d,lapack_complex_double *du,lapack_complex_double *du2,lapack_int *ipiv );

    Fptr_NL_LAPACKE_zgttrf ZGTTRF;
    void *hModule,*dModule;


    dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == hModule) || (NULL == dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit( -1 );
    }

    ZGTTRF = (Fptr_NL_LAPACKE_zgttrf)dlsym(hModule,"LAPACKE_zgttrf");
    ASSERT_TRUE(ZGTTRF != NULL) << "failed to get the Netlib LAPACKE_zgttrf symbol";

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    zgttrf_obj->inforef = ZGTTRF( zgttrf_obj->n,zgttrf_obj->dlref,
                           zgttrf_obj->dref,zgttrf_obj->duref,
                           zgttrf_obj->du2ref,zgttrf_obj->ipivref);

   /* Compute libflame's Lapacke o/p  */
    zgttrf_obj->info     = LAPACKE_zgttrf( zgttrf_obj->n,zgttrf_obj->dl,
                           zgttrf_obj->d,zgttrf_obj->du,zgttrf_obj->du2,
                           zgttrf_obj->ipiv);

    if( zgttrf_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zgttrf \
                  is wrong\n",zgttrf_obj->info );
    }
    if( zgttrf_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zgttrf is wrong\n",
        zgttrf_obj->inforef );
    }
    dlclose(hModule);
    dlclose(dModule);
}

TEST_F(zgttrf_test, zgttrf1) {
double diff;    int ipiv_diff;
    ipiv_diff = computeDiff_i( zgttrf_obj->n,zgttrf_obj->ipiv,zgttrf_obj->ipivref);
    EXPECT_EQ(ipiv_diff, 0) << " pivot computation  failed ";
    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_z( zgttrf_obj->n,zgttrf_obj->d,zgttrf_obj->dref );
    diff +=  computeDiff_z( zgttrf_obj->n-1,zgttrf_obj->du,zgttrf_obj->duref );
    diff +=  computeDiff_z( zgttrf_obj->n-1,zgttrf_obj->dl,zgttrf_obj->dlref );
    diff +=  computeDiff_z( zgttrf_obj->n-2,zgttrf_obj->du2,zgttrf_obj->du2ref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgttrf_test, zgttrf2) {
double diff;    int ipiv_diff;
    ipiv_diff = computeDiff_i( zgttrf_obj->n,zgttrf_obj->ipiv,zgttrf_obj->ipivref);
    EXPECT_EQ(ipiv_diff, 0) << " pivot computation  failed ";
    /* Compute Difference in C and CPP buffer */
    diff =  computeDiff_z( zgttrf_obj->n,zgttrf_obj->d,zgttrf_obj->dref );
    diff +=  computeDiff_z( zgttrf_obj->n-1,zgttrf_obj->du,zgttrf_obj->duref );
    diff +=  computeDiff_z( zgttrf_obj->n-1,zgttrf_obj->dl,zgttrf_obj->dlref );
    diff +=  computeDiff_z( zgttrf_obj->n-2,zgttrf_obj->du2,zgttrf_obj->du2ref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}