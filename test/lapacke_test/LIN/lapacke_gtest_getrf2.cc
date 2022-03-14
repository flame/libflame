#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"

TEST(getrf2,dgetrf21) {

    /* LAPACKE DGETRF2 prototype */
    typedef int (*Fptr_NL_LAPACKE_dgetrf2) ( int matrix_layout,lapack_int m,lapack_int n,
                                    double* a,lapack_int lda,lapack_int* ipiv );
						    
    Fptr_NL_LAPACKE_dgetrf2 DGETRF2;		
    double diff;
    double_common_parameters dgetrf2_obj(512,1024);
    
    dgetrf2_obj.dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    dgetrf2_obj.hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == dgetrf2_obj.hModule) || (NULL == dgetrf2_obj.dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }
    
    DGETRF2 = (Fptr_NL_LAPACKE_dgetrf2)dlsym(dgetrf2_obj.hModule,"LAPACKE_dgetrf2");
    if (NULL == DGETRF2)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(dgetrf2_obj.hModule);
      dlclose(dgetrf2_obj.dModule);
      exit (-1);
    }
        /* Compute libflame's Lapacke o/p  */
    dgetrf2_obj.info     = LAPACKE_dgetrf2( LAPACK_COL_MAJOR,dgetrf2_obj.m,dgetrf2_obj.n,dgetrf2_obj.A,
	                                             dgetrf2_obj.lda,dgetrf2_obj.ipiv);
	  
        /* Compute the reference o/p by invoking Netlib-Lapack's API */	
        dgetrf2_obj.inforef = DGETRF2( LAPACK_COL_MAJOR,dgetrf2_obj.m,dgetrf2_obj.n,dgetrf2_obj.Aref,
	                                          dgetrf2_obj.lda,dgetrf2_obj.ipivref);

	if( dgetrf2_obj.info < 0 ) {
		printf( "\n warning: The i:%d th argument with libflame LAPACKE_dgetrf2 is wrong\n",
		            dgetrf2_obj.info );
	}
	if( dgetrf2_obj.inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_dgetrf2 is wrong\n",
		dgetrf2_obj.inforef );
	}

	/* Compute Difference in C and CPP buffer */
	diff =  computeDiff_d( dgetrf2_obj.m*dgetrf2_obj.n,dgetrf2_obj.A,dgetrf2_obj.Aref );

	EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}


TEST(getrf2,sgetrf21) {

    /* LAPACKE sgetrf2 prototype */
    typedef int (*Fptr_NL_LAPACKE_sgetrf2) ( int matrix_layout,lapack_int m,lapack_int n,
                                    float* a,lapack_int lda,lapack_int* ipiv );
						    
    Fptr_NL_LAPACKE_sgetrf2 sgetrf2;		
    float diff;    
	
    float_common_parameters sgetrf2_obj(900,900);
    
    sgetrf2_obj.dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);    
    sgetrf2_obj.hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);    
    if ((NULL == sgetrf2_obj.hModule) || (NULL == sgetrf2_obj.dModule) )
    {
    	    printf("Load Library failed. Exiting ....\n");
    	    exit( -1 );
    }
    
    sgetrf2 = (Fptr_NL_LAPACKE_sgetrf2)dlsym(sgetrf2_obj.hModule,"LAPACKE_sgetrf2");
    if (NULL == sgetrf2)
    {
    	  printf("Could not get the symbol. Exiting...\n");
    	  dlclose(sgetrf2_obj.hModule);
    	  dlclose(sgetrf2_obj.dModule);
    	  exit( -1 );
    }

    /* Compute libflame's Lapacke o/p  */
    sgetrf2_obj.info     = LAPACKE_sgetrf2( LAPACK_COL_MAJOR,sgetrf2_obj.m,sgetrf2_obj.n,sgetrf2_obj.A,
                                          sgetrf2_obj.lda,sgetrf2_obj.ipiv);

    /* Compute the reference o/p by invoking Netlib-Lapack's API */	
    sgetrf2_obj.inforef = sgetrf2( LAPACK_COL_MAJOR,sgetrf2_obj.m,sgetrf2_obj.n,sgetrf2_obj.Aref,
                                        sgetrf2_obj.lda,sgetrf2_obj.ipivref);

	if( sgetrf2_obj.info < 0 ) {
		printf( "\n warning: The i:%d th argument with libflame LAPACKE_sgetrf2 is wrong\n",
		            sgetrf2_obj.info );
	}
	if( sgetrf2_obj.inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_sgetrf2 is wrong\n",
		sgetrf2_obj.inforef );
	}

	/* Compute Difference in C and CPP buffer */
	diff =  computeDiff_s( sgetrf2_obj.m*sgetrf2_obj.n,sgetrf2_obj.A,sgetrf2_obj.Aref );
	EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

TEST(getrf2,cgetrf21) {

    /* LAPACKE sgetrf2 prototype */
    typedef int (*Fptr_NL_LAPACKE_cgetrf2) ( int matrix_layout,lapack_int m,lapack_int n,
                                    lapack_complex_float* a,lapack_int lda,lapack_int* ipiv );

    Fptr_NL_LAPACKE_cgetrf2 cgetrf2;		
    float diff;    
	
    scomplex_common_parameters cgetrf2_obj(900,900);
    
    cgetrf2_obj.dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);    
    cgetrf2_obj.hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);    
    if ((NULL == cgetrf2_obj.hModule) || (NULL == cgetrf2_obj.dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit( -1 );
    }
    
    cgetrf2 = (Fptr_NL_LAPACKE_cgetrf2)dlsym(cgetrf2_obj.hModule,"LAPACKE_cgetrf2");
    if (NULL == cgetrf2)
    {
    	printf("Could not get the symbol. Exiting...\n");
    	dlclose(cgetrf2_obj.hModule);
    	dlclose(cgetrf2_obj.dModule);
    	exit( -1 );
    }

    /* Compute libflame's Lapacke o/p  */
    cgetrf2_obj.info     = LAPACKE_cgetrf2( LAPACK_COL_MAJOR,cgetrf2_obj.m,cgetrf2_obj.n,cgetrf2_obj.A,
                                          cgetrf2_obj.lda,cgetrf2_obj.ipiv);

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    cgetrf2_obj.inforef = cgetrf2( LAPACK_COL_MAJOR,cgetrf2_obj.m,cgetrf2_obj.n,cgetrf2_obj.Aref,
                                        cgetrf2_obj.lda,cgetrf2_obj.ipivref);

	if( cgetrf2_obj.info < 0 ) {
		printf( "\n warning: The i:%d th argument with libflame LAPACKE_cgetrf2 is wrong\n",
		            cgetrf2_obj.info );
	}
	if( cgetrf2_obj.inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_cgetrf2 is wrong\n",
		cgetrf2_obj.inforef );
	}

	/* Compute Difference in C and CPP buffer */
	diff =  computeDiff_c( cgetrf2_obj.m*cgetrf2_obj.n,cgetrf2_obj.A,cgetrf2_obj.Aref );
	EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

TEST(getrf2,zgetrf21) {

    /* LAPACKE zgetrf2 prototype */
    typedef int (*Fptr_NL_LAPACKE_zgetrf2) ( int matrix_layout,lapack_int m,lapack_int n,
                             lapack_complex_double* a,lapack_int lda,lapack_int* ipiv );
    Fptr_NL_LAPACKE_zgetrf2 ZGETRF2;		
    float diff;    
	
    dcomplex_common_parameters zgetrf2_obj(900,900);
    
    zgetrf2_obj.dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);    
    zgetrf2_obj.hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);    
    if ((NULL == zgetrf2_obj.hModule) || (NULL == zgetrf2_obj.dModule) )
    {
    	printf("Load Library failed. Exiting ....\n");
    	exit( -1 );
    }
    
    ZGETRF2 = (Fptr_NL_LAPACKE_zgetrf2)dlsym(zgetrf2_obj.hModule,"LAPACKE_zgetrf2");
    if (NULL == ZGETRF2)
    {
    	printf("Could not get the symbol. Exiting...\n");
    	dlclose(zgetrf2_obj.hModule);
    	dlclose(zgetrf2_obj.dModule);
    	exit( -1 );
    }

    /* Compute libflame's Lapacke o/p  */
    zgetrf2_obj.info     = LAPACKE_zgetrf2( LAPACK_COL_MAJOR,zgetrf2_obj.m,zgetrf2_obj.n,
                                          zgetrf2_obj.A,zgetrf2_obj.lda,zgetrf2_obj.ipiv);

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    zgetrf2_obj.inforef = ZGETRF2( LAPACK_COL_MAJOR,zgetrf2_obj.m,zgetrf2_obj.n,
                                 zgetrf2_obj.Aref,zgetrf2_obj.lda,zgetrf2_obj.ipivref);

	if( zgetrf2_obj.info < 0 ) {
		printf( "\n warning: The i:%d th argument with libflame LAPACKE_zgetrf2 is wrong\n",
		            zgetrf2_obj.info );
	}
	if( zgetrf2_obj.inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_zgetrf2 is wrong\n",
		zgetrf2_obj.inforef );
	}

	/* Compute Difference in C and CPP buffer */
	diff =  computeDiff_z( zgetrf2_obj.m*zgetrf2_obj.n,zgetrf2_obj.A,zgetrf2_obj.Aref );
	EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}
