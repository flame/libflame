#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"
#include <math.h>
#include "lapacke.h"

#define ppcon_free() \
       free (ap   ); \
       free (apref);

class ppcon_double_parameters{
   public:
   /* input params to the API **/
   int  matrix_layout;
   char uplo;
   int  n;
   double  *ap,*apref;
   double anorm;
   double rcond,rcondref;
   int info,inforef;
   public: 
      ppcon_double_parameters (int matrix_layout,char uplo,int n,
                               double anorm );
      ~ppcon_double_parameters ();

}; /* end of ppcon_double_parameters  class definition */

/* Destructor definition **/
ppcon_double_parameters:: ~ppcon_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" ppcon_double_parameters object: destructor invoked. \n");
#endif
   ppcon_free();
}

ppcon_double_parameters::ppcon_double_parameters (int matrix_layout_i,
                             char uplo_i,int n_i,
                             double anormi ){
   int j;
   matrix_layout = matrix_layout_i;
   uplo = uplo_i;
   n = n_i;
   anorm = anormi;

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_double_buffer_pair( &ap, &apref, (n*(n+1)/2));
    if( (ap==NULL) || (apref==NULL) ){
       ppcon_free();
       printf(" ppcon_double_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( ap,apref,(n*(n+1)/2));
    info = 0;
    inforef = 0;
}

TEST(ppcon,dppcon1) {

   /* LAPACKE DPPCON prototype */
   typedef int (*Fptr_NL_LAPACKE_dppcon) ( int matrix_layout,char uplo,
                          lapack_int n,const double* ap, 
                          double anorm,double* rcond  );
                          
   Fptr_NL_LAPACKE_dppcon DPPCON;
   double diff;
   void *hModule,*dModule;

   ppcon_double_parameters dppcon_obj(LAPACK_ROW_MAJOR,
                             'L',512,121.0 );
   dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
   hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);

   if ((NULL == hModule) || (NULL == dModule) )
   {
       printf("Load Library failed. Exiting ....\n");
       exit (0);
   }
   DPPCON = (Fptr_NL_LAPACKE_dppcon)dlsym(hModule,"LAPACKE_dppcon");
   if (NULL == DPPCON)
   {
      printf("Could not get the symbol -LAPACKE_dppcon- . Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
   }
   /* Compute libflame's Lapacke o/p  */
   dppcon_obj.info    = LAPACKE_dppcon( dppcon_obj.matrix_layout,dppcon_obj.uplo,
                                        dppcon_obj.n,(const double*)dppcon_obj.ap,
                                        dppcon_obj.anorm,&dppcon_obj.rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */   
   dppcon_obj.inforef = DPPCON( dppcon_obj.matrix_layout,dppcon_obj.uplo,
                                dppcon_obj.n,(const double*)dppcon_obj.apref,
                                dppcon_obj.anorm,&dppcon_obj.rcondref);
    if( dppcon_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with LAPACKE_dppcon is wrong\n",
                 dppcon_obj.info );
    }
    if( dppcon_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dppcon is wrong\n",
        dppcon_obj.inforef );
    }

    /* Compute Difference in C and CPP buffer */
#if LAPACKE_TEST_VERBOSE
    printf("dppcon_obj.rcond: %lf   dppcon_obj.rcondref: %lf \n",
	                       dppcon_obj.rcond,dppcon_obj.rcondref);
    printf("diff: %lf  \n",diff);
#endif
    EXPECT_NEAR(0.0,abs(dppcon_obj.rcond - dppcon_obj.rcondref),DOUBLE_DIFF_THRESHOLD);
}



class ppcon_float_parameters{
   public:
   /* input params to the API **/
   int  matrix_layout;
   char uplo;
   int  n;
   float  *ap,*apref;
   float anorm;
   float rcond,rcondref;   int info,inforef;
   public: 
      ppcon_float_parameters (int matrix_layout,char uplo,int n,
                             float anorm );
      ~ppcon_float_parameters ();

}; /* end of ppcon_float_parameters  class definition */

/* Destructor definition **/
ppcon_float_parameters:: ~ppcon_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" ppcon_float_parameters object: destructor invoked. \n");
#endif
   ppcon_free();
}

ppcon_float_parameters::ppcon_float_parameters (int matrix_layout_i,
                             char uplo_i,int n_i,float anormi ){
   int j;
   matrix_layout = matrix_layout_i;
   uplo = uplo_i;
   n = n_i;
   anorm = anormi;

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_float_buffer_pair( &ap, &apref, (n*(n+1)/2));
    if( (ap==NULL) || (apref==NULL)){
       ppcon_free();
       printf(" ppcon_float_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( ap,apref,(n*(n+1)/2));
    info = 0;
    inforef = 0;
}

TEST(ppcon,sppcon1) {

   /* LAPACKE SPPCON prototype */
   typedef int (*Fptr_NL_LAPACKE_sppcon) ( int matrix_layout,char uplo,lapack_int n,
                                          const float *a, 
                                          float anorm,float *rcond );
   Fptr_NL_LAPACKE_sppcon SPPCON;
   float diff;
   void *hModule,*dModule;

   ppcon_float_parameters sppcon_obj(LAPACK_COL_MAJOR,
                             'U',1024,0.86 );
   dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
   hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
   if ((NULL == hModule) || (NULL == dModule) )
   {
       printf("Load Library failed. Exiting ....\n");
       exit (0);
   }
   SPPCON = (Fptr_NL_LAPACKE_sppcon)dlsym(hModule,"LAPACKE_sppcon");
   if (NULL == SPPCON)
   {
      printf("Could not get the symbol -LAPACKE_sppcon- . Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
   }
   /* Compute libflame's Lapacke o/p  */
   sppcon_obj.info    = LAPACKE_sppcon( sppcon_obj.matrix_layout,sppcon_obj.uplo,
                                        sppcon_obj.n,
                                        sppcon_obj.ap,
                                        sppcon_obj.anorm,
                                        &sppcon_obj.rcond);

    /* Compute Difference in C and CPP buffer */
#if LAPACKE_TEST_VERBOSE
    printf("sppcon_obj.rcond: %lf   sppcon_obj.rcondref: %lf \n",
	                        sppcon_obj.rcond,sppcon_obj.rcondref);
    printf("diff: %f  \n",diff);
#endif
   EXPECT_NEAR(0.0,abs(sppcon_obj.rcond - sppcon_obj.rcondref),DOUBLE_DIFF_THRESHOLD);
}


class ppcon_scomplex_parameters{
   public:
   /* input params to the API **/
   int  matrix_layout;
   char uplo;
   int  n;
   lapack_complex_float  *ap,*apref;
   float anorm;
   float rcond,rcondref;   int info,inforef;
   public: 
      ppcon_scomplex_parameters (int matrix_layout,char uplo,int n,
                                 float anorm );
      ~ppcon_scomplex_parameters ();

}; /* end of ppcon_float_parameters  class definition */

/* Destructor definition **/
ppcon_scomplex_parameters:: ~ppcon_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" ppcon_float_parameters object: destructor invoked. \n");
#endif
   ppcon_free();
}

ppcon_scomplex_parameters::ppcon_scomplex_parameters (int matrix_layout_i,
                             char uplo_i,int n_i,float anormi ){
   int j;
   matrix_layout = matrix_layout_i;
   uplo = uplo_i;
   n = n_i;
   anorm = anormi;

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &ap, &apref, (n*(n+1)/2));
    if( (ap==NULL) || (apref==NULL)){
       ppcon_free();
       printf(" ppcon_float_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( ap,apref,(n*(n+1)/2));
    info = 0;
    inforef = 0;
}

TEST(ppcon,cppcon1) {

   /* LAPACKE CPPCON prototype */
   typedef int (*Fptr_NL_LAPACKE_cppcon) ( int matrix_layout,char uplo,lapack_int n,
                      const lapack_complex_float *a,
                      float anorm,float *rcond );
   Fptr_NL_LAPACKE_cppcon CPPCON;
   float diff;
   void *hModule,*dModule;

   ppcon_scomplex_parameters cppcon_obj(LAPACK_COL_MAJOR,
                             'U',256,11.0 );
   dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
   hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
   if ((NULL == hModule) || (NULL == dModule) )
   {
       printf("Load Library failed. Exiting ....\n");
       exit (0);
   }
   CPPCON = (Fptr_NL_LAPACKE_cppcon)dlsym(hModule,"LAPACKE_cppcon");
   if (NULL == CPPCON)
   {
      printf("Could not get the symbol -LAPACKE_cppcon- . Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
   }
   /* Compute libflame's Lapacke o/p  */  
   cppcon_obj.info    = LAPACKE_cppcon( cppcon_obj.matrix_layout,cppcon_obj.uplo,
                                        cppcon_obj.n,
                                        cppcon_obj.ap,
                                        cppcon_obj.anorm,
                                        &cppcon_obj.rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */   
   cppcon_obj.inforef = CPPCON( cppcon_obj.matrix_layout,cppcon_obj.uplo,
                                cppcon_obj.n,
                                cppcon_obj.apref,
                                cppcon_obj.anorm,
                                &cppcon_obj.rcondref);

    if( cppcon_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with LAPACKE_cppcon is wrong\n",
                 cppcon_obj.info );
    }
    if( cppcon_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cppcon is wrong\n",
        cppcon_obj.inforef );
    }

#if LAPACKE_TEST_VERBOSE
    printf("cppcon_obj.rcond: %lf   cppcon_obj.rcondref: %lf \n",
	                       cppcon_obj.rcond,cppcon_obj.rcondref);
    printf("diff: %f  \n",diff);
#endif
    EXPECT_NEAR(0.0,abs(cppcon_obj.rcond - cppcon_obj.rcondref),DOUBLE_DIFF_THRESHOLD);
}

class ppcon_dcomplex_parameters{
   public:
   /* input params to the API **/
   int  matrix_layout;
   char uplo;
   int  n;
   lapack_complex_double  *ap,*apref;
   double anorm;
   double rcond,rcondref;   int info,inforef;
   public: 
      ppcon_dcomplex_parameters (int matrix_layout,char uplo,int n,
                                 double anorm );
      ~ppcon_dcomplex_parameters ();

}; /* end of ppcon_float_parameters  class definition */

/* Destructor definition **/
ppcon_dcomplex_parameters:: ~ppcon_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" ppcon_float_parameters object: destructor invoked. \n");
#endif
   ppcon_free();
}

ppcon_dcomplex_parameters::ppcon_dcomplex_parameters (int matrix_layout_i,
                             char uplo_i,int n_i,double anormi ){
   int j;
   matrix_layout = matrix_layout_i;
   uplo = uplo_i;
   n = n_i;
   anorm = anormi;

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &ap, &apref, (n*(n+1)/2));
    if( (ap==NULL) || (apref==NULL) ){
       ppcon_free();
       printf(" ppcon_float_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( ap,apref,(n*(n+1)/2));
    info = 0;
    inforef = 0;
}

TEST(ppcon,zppcon1) {

   /* LAPACKE ZPPCON prototype */
   typedef int (*Fptr_NL_LAPACKE_zppcon) ( int matrix_layout,char uplo,lapack_int n,
                      const lapack_complex_double *a,
                      double anorm,double *rcond );
   Fptr_NL_LAPACKE_zppcon ZPPCON;
   double diff;
   void *hModule,*dModule;

   ppcon_dcomplex_parameters zppcon_obj(LAPACK_ROW_MAJOR,
                             'L',512,167.0 );
   dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
   hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
   if ((NULL == hModule) || (NULL == dModule) )
   {
       printf("Load Library failed. Exiting ....\n");
       exit (0);
   }
   ZPPCON = (Fptr_NL_LAPACKE_zppcon)dlsym(hModule,"LAPACKE_zppcon");
   if (NULL == ZPPCON)
   {
      printf("Could not get the symbol -LAPACKE_zppcon- . Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
   }
   /* Compute libflame's Lapacke o/p  */  
   zppcon_obj.info    = LAPACKE_zppcon( zppcon_obj.matrix_layout,zppcon_obj.uplo,
                                        zppcon_obj.n,
                                        zppcon_obj.ap,
                                        zppcon_obj.anorm,
                                        &zppcon_obj.rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */   
   zppcon_obj.inforef = ZPPCON( zppcon_obj.matrix_layout,zppcon_obj.uplo,
                                zppcon_obj.n,
                                zppcon_obj.apref,
                                zppcon_obj.anorm,
                                &zppcon_obj.rcondref);

    if( zppcon_obj.info < 0 ) {
        printf( "\n warning: The i:%d th argument with LAPACKE_cppcon is wrong\n",
                 zppcon_obj.info );
    }
    if( zppcon_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cppcon is wrong\n",
        zppcon_obj.inforef );
    }

    /* Compute Difference in C and CPP buffer */
#if LAPACKE_TEST_VERBOSE
    printf("zppcon_obj.rcond: %lf   zppcon_obj.rcondref: %lf \n",zppcon_obj.rcond,zppcon_obj.rcondref);
    printf("diff: %lf  \n",diff);
#endif
    EXPECT_NEAR(0.0,abs(zppcon_obj.rcond - zppcon_obj.rcondref),DOUBLE_DIFF_THRESHOLD);

}
