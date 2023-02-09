#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"
#include "lapacke_gtest_gbbrd.h"


gbbrd_double_parameters::gbbrd_double_parameters (int matrix_layout_i,
          char vect_i, int m_i, int n_i, int ncc_i, int kl_i, int ku_i,
          int ldab_i, int ldq_i, int ldpt_i, int ldc_i )
{
	int j, temp;
	
	matrix_layout = matrix_layout_i;
	vect = vect_i;
	m = m_i;
	n = n_i;
	ncc = ncc_i;
	kl = kl_i;
	ku = ku_i;
	ldab = ldab_i;
	ldq = ldq_i;
	ldpt = ldpt_i;
	ldc = ldc_i;
	
/* Memory allocation of the buffers */
#if 1
    lapacke_gtest_alloc_double_buffer_pair( &ab, &abref, (ldab*n));
	lapacke_gtest_alloc_double_buffer_pair( &c,  &cref,  (ldc*ncc));
	
	temp = min(m,n);	
	lapacke_gtest_alloc_double_buffer_pair( &d, &dref, temp);
	lapacke_gtest_alloc_double_buffer_pair( &e, &eref, (temp-1));
	lapacke_gtest_alloc_double_buffer_pair( &q, &qref, (ldq*m));
	lapacke_gtest_alloc_double_buffer_pair( &pt, &ptref, (ldpt*n));
	lapacke_gtest_alloc_double_buffer_pair( &work, &workref, (2*max(m,n)));
#else
	ab    = (double *)malloc((ldab*n)*sizeof(double));
	abref = (double *)malloc((ldab*n)*sizeof(double));

    c    = (double *)malloc((ldc*ncc)*sizeof(double));
	cref = (double *)malloc((ldc*ncc)*sizeof(double));

    temp = min(m,n);
    d    = (double *)malloc(temp*sizeof(double));
	dref = (double *)malloc(temp*sizeof(double));
   
    e    = (double *)malloc((temp-1)*sizeof(double));
	eref = (double *)malloc((temp-1)*sizeof(double));

    q    = (double *)malloc((ldq*m)*sizeof(double));
	qref = (double *)malloc((ldq*m)*sizeof(double));
	
	pt    = (double *)malloc((ldpt*n)*sizeof(double));
	ptref = (double *)malloc((ldpt*n)*sizeof(double));

	j = 2*max(m,n);
    work    = (double *)malloc(j*sizeof(double));
	workref = (double *)malloc(j*sizeof(double));
#endif
	if( (d==NULL) || (dref==NULL) || (e==NULL) || (eref==NULL)|| \
	    (ab==NULL) || (abref==NULL) || (q==NULL) || (qref==NULL)|| \
		(c==NULL) || (cref==NULL) || (work==NULL) || \
		(workref==NULL) || (pt==NULL) || (ptref==NULL) ){

       printf(" gbbrd_double_parameters object: malloc error. Exiting...\n");
       gbbrd_free();
	   exit (-1);
	}

    /* Initialization of input Buffers */
#if 1
    lapacke_gtest_init_double_buffer_pair_rand( ab, abref, (ldab*n));
    lapacke_gtest_init_double_buffer_pair_rand( c, cref, (ldc*ncc));
	lapacke_gtest_init_double_buffer_pair_with_constant( d, dref, temp, 0.0);
	lapacke_gtest_init_double_buffer_pair_with_constant( e, eref, temp-1, 0.0);
	lapacke_gtest_init_double_buffer_pair_with_constant( q, qref, (ldq*m), 0.0);
	lapacke_gtest_init_double_buffer_pair_with_constant( pt, ptref, (ldpt*n), 0.0);
	lapacke_gtest_init_double_buffer_pair_with_constant( work, workref, (2*max(m,n)), 0.0);
#else
		
    for( j = 0; j < (ldab*n); j++ ) {
       ab[j] = ((double) rand()) / ((double) RAND_MAX) - 0.5;
       abref[j] = ab[j];
    }

    for( j = 0; j < (ldc*ncc); j++ ) {
       c[j] = ((double) rand()) / ((double) RAND_MAX) - 0.5;
       cref[j] = c[j];
    }
	/** Below o/p buffers initialized to '0.0' **/  
    for( j = 0; j < temp; j++ ) {
       d[j] = 0.0;
       dref[j] = 0.0;
    }

    for( j = 0; j < temp-1; j++ ) {
       e[j] = 0.0;
       eref[j] = 0.0;
    }

    for( j = 0; j < (ldq*m); j++ ) {
       q[j] = 0.0;
       qref[j] = 0.0;
    }

    for( j = 0; j < (ldpt*n); j++ ) {
       pt[j] = 0.0;
       ptref[j] = 0.0;
    }

    for( j = 0; j < (2*max(m,n)); j++ ) {
       work[j] = 0.0;
       workref[j] = 0.0;
    }
	
#endif

	info = 0;
	inforef = 0;
	
}

/* Destructor definition **/
gbbrd_double_parameters:: ~gbbrd_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gbbrd_double_parameters object: destructor invoked. \n");
#endif
   gbbrd_free();
}


/** constructor for 'gbbrd_float_parameters' class  */
gbbrd_float_parameters::gbbrd_float_parameters (int matrix_layout_i,
          char vect_i, int m_i, int n_i, int ncc_i, int kl_i, int ku_i,
          int ldab_i, int ldq_i, int ldpt_i, int ldc_i )
{
	int j, temp;
	
	matrix_layout = matrix_layout_i;
	vect = vect_i;
	m = m_i;
	n = n_i;
	ncc = ncc_i;
	kl = kl_i;
	ku = ku_i;
	ldab = ldab_i;
	ldq = ldq_i;
	ldpt = ldpt_i;
	ldc = ldc_i;
	
/* Memory allocation of the buffers */
	ab    = (float *)malloc((ldab*n)*sizeof(float));
	abref = (float *)malloc((ldab*n)*sizeof(float));

    c    = (float *)malloc((ldc*ncc)*sizeof(float));
	cref = (float *)malloc((ldc*ncc)*sizeof(float));

    temp = min(m,n);
    d    = (float *)malloc(temp*sizeof(float));
	dref = (float *)malloc(temp*sizeof(float));
   
    e    = (float *)malloc((temp-1)*sizeof(float));
	eref = (float *)malloc((temp-1)*sizeof(float));

    q    = (float *)malloc((ldq*m)*sizeof(float));
	qref = (float *)malloc((ldq*m)*sizeof(float));
	
	pt    = (float *)malloc((ldpt*n)*sizeof(float));
	ptref = (float *)malloc((ldpt*n)*sizeof(float));

	j = 2*max(m,n);
    work    = (float *)malloc(j*sizeof(float));
	workref = (float *)malloc(j*sizeof(float));
	
	if( (d==NULL) || (dref==NULL) || (e==NULL) || (eref==NULL)|| \
	    (ab==NULL) || (abref==NULL) || (q==NULL) || (qref==NULL)|| \
		(c==NULL) || (cref==NULL) || (work==NULL) || \
		(workref==NULL) || (pt==NULL) || (ptref==NULL) ){

       printf(" gbbrd_float_parameters object: malloc error. Exiting...\n");
       gbbrd_free();
	   exit (-1);
	}

    /* Initialization of input Buffers */

    for( j = 0; j < (ldab*n); j++ ) {
       ab[j] = ((float) rand()) / ((float) RAND_MAX) - 0.5;
       abref[j] = ab[j];
    }

    for( j = 0; j < (ldc*ncc); j++ ) {
       c[j] = ((float) rand()) / ((float) RAND_MAX) - 0.5;
       cref[j] = c[j];
    }
	/** Below o/p buffers initialized to '0.0' **/  
    for( j = 0; j < temp; j++ ) {
       d[j] = 0.0;
       dref[j] = 0.0;
    }

    for( j = 0; j < temp-1; j++ ) {
       e[j] = 0.0;
       eref[j] = 0.0;
    }

    for( j = 0; j < (ldq*m); j++ ) {
       q[j] = 0.0;
       qref[j] = 0.0;
    }

    for( j = 0; j < (ldpt*n); j++ ) {
       pt[j] = 0.0;
       ptref[j] = 0.0;
    }

    for( j = 0; j < (2*max(m,n)); j++ ) {
       work[j] = 0.0;
       workref[j] = 0.0;
    }
	info = 0;
	inforef = 0;
	
}

/* Destructor definition **/
gbbrd_float_parameters:: ~gbbrd_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gbbrd_float_parameters object: destructor invoked. \n");
#endif
   gbbrd_free();
}



gbbrd_scomplex_parameters::gbbrd_scomplex_parameters (int matrix_layout_i,
          char vect_i, int m_i, int n_i, int ncc_i, int kl_i, int ku_i,
          int ldab_i, int ldq_i, int ldpt_i, int ldc_i )
{
	int j, temp;
	float dtemp_r, dtemp_i;
	
	matrix_layout = matrix_layout_i;
	vect = vect_i;
	m = m_i;
	n = n_i;
	ncc = ncc_i;
	kl = kl_i;
	ku = ku_i;
	ldab = ldab_i;
	ldq = ldq_i;
	ldpt = ldpt_i;
	ldc = ldc_i;
	
/* Memory allocation of the buffers */
	ab    = (lapack_complex_float *)malloc((ldab*n)*sizeof(lapack_complex_float));
	abref = (lapack_complex_float *)malloc((ldab*n)*sizeof(lapack_complex_float));

    c    = (lapack_complex_float *)malloc((ldc*ncc)*sizeof(lapack_complex_float));
	cref = (lapack_complex_float *)malloc((ldc*ncc)*sizeof(lapack_complex_float));

    temp = min(m,n);
    d    = (float *)malloc(temp*sizeof(float));
	dref = (float *)malloc(temp*sizeof(float));
   
    e    = (float *)malloc((temp-1)*sizeof(float));
	eref = (float *)malloc((temp-1)*sizeof(float));

    q    = (lapack_complex_float *)malloc((ldq*m)*sizeof(lapack_complex_float));
	qref = (lapack_complex_float *)malloc((ldq*m)*sizeof(lapack_complex_float));
	
	pt    = (lapack_complex_float *)malloc((ldpt*n)*sizeof(lapack_complex_float));
	ptref = (lapack_complex_float *)malloc((ldpt*n)*sizeof(lapack_complex_float));

	j = 2*max(m,n);
    work    = (lapack_complex_float *)malloc(j*sizeof(lapack_complex_float));
	workref = (lapack_complex_float *)malloc(j*sizeof(lapack_complex_float));
	
	if( (d==NULL) || (dref==NULL) || (e==NULL) || (eref==NULL)|| \
	    (ab==NULL) || (abref==NULL) || (q==NULL) || (qref==NULL)|| \
		(c==NULL) || (cref==NULL) || (work==NULL) || \
		(workref==NULL) || (pt==NULL) || (ptref==NULL) ){

       printf(" gbbrd_scomplex_parameters object: malloc error. Exiting...\n");
       gbbrd_free();
	   exit (-1);
	}

    /* Initialization of input Buffers */

    for( j = 0; j < (ldab*n); j++ ) {
	   dtemp_r =(((float) rand()) / ((float) RAND_MAX) - 0.5);
	   dtemp_i =(((float) rand()) / ((float) RAND_MAX) - 0.5);
	   ab[j] = lapack_make_complex_float( dtemp_r, dtemp_i);
       abref[j] = lapack_make_complex_float( dtemp_r, dtemp_i);
    }

    for( j = 0; j < (ldc*ncc); j++ ) {
	   dtemp_r =(((float) rand()) / ((float) RAND_MAX) - 0.5);
	   dtemp_i =(((float) rand()) / ((float) RAND_MAX) - 0.5);
	   c[j] = lapack_make_complex_float( dtemp_r, dtemp_i);
       cref[j] = lapack_make_complex_float( dtemp_r, dtemp_i);
    }

	/** Below o/p buffers initialized to '0.0' **/  
    for( j = 0; j < temp; j++ ) {
       d[j] = 0.0;
       dref[j] = 0.0;
    }

    for( j = 0; j < temp-1; j++ ) {
       e[j] = 0.0;
       eref[j] = 0.0;
    }

    for( j = 0; j < (ldq*m); j++ ) {
       q[j] = lapack_make_complex_float( 0.0, 0.0);
       qref[j] = lapack_make_complex_float( 0.0, 0.0);
    }

    for( j = 0; j < (ldpt*n); j++ ) {
       pt[j] = lapack_make_complex_float( 0.0, 0.0);
       ptref[j] = lapack_make_complex_float( 0.0, 0.0);
    }

    for( j = 0; j < (2*max(m,n)); j++ ) {
       work[j] = lapack_make_complex_float( 0.0, 0.0);
       workref[j] = lapack_make_complex_float( 0.0, 0.0);
    }
	info = 0;
	inforef = 0;
	
}

/* Destructor definition **/
gbbrd_scomplex_parameters:: ~gbbrd_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gbbrd_scomplex_parameters object: destructor invoked. \n");
#endif
   gbbrd_free();
}


gbbrd_dcomplex_parameters::gbbrd_dcomplex_parameters (int matrix_layout_i,
          char vect_i, int m_i, int n_i, int ncc_i, int kl_i, int ku_i,
          int ldab_i, int ldq_i, int ldpt_i, int ldc_i )
{
	int j, temp;
	double dtemp_r, dtemp_i;
	
	matrix_layout = matrix_layout_i;
	vect = vect_i;
	m = m_i;
	n = n_i;
	ncc = ncc_i;
	kl = kl_i;
	ku = ku_i;
	ldab = ldab_i;
	ldq = ldq_i;
	ldpt = ldpt_i;
	ldc = ldc_i;
	
/* Memory allocation of the buffers */
	ab    = (lapack_complex_double *)malloc((ldab*n)*sizeof(lapack_complex_double));
	abref = (lapack_complex_double *)malloc((ldab*n)*sizeof(lapack_complex_double));

    c    = (lapack_complex_double *)malloc((ldc*ncc)*sizeof(lapack_complex_double));
	cref = (lapack_complex_double *)malloc((ldc*ncc)*sizeof(lapack_complex_double));

    temp = min(m,n);
    d    = (double *)malloc(temp*sizeof(double));
	dref = (double *)malloc(temp*sizeof(double));
   
    e    = (double *)malloc((temp-1)*sizeof(double));
	eref = (double *)malloc((temp-1)*sizeof(double));

    q    = (lapack_complex_double *)malloc((ldq*m)*sizeof(lapack_complex_double));
	qref = (lapack_complex_double *)malloc((ldq*m)*sizeof(lapack_complex_double));
	
	pt    = (lapack_complex_double *)malloc((ldpt*n)*sizeof(lapack_complex_double));
	ptref = (lapack_complex_double *)malloc((ldpt*n)*sizeof(lapack_complex_double));

	j = 2*max(m,n);
    work    = (lapack_complex_double *)malloc(j*sizeof(lapack_complex_double));
	workref = (lapack_complex_double *)malloc(j*sizeof(lapack_complex_double));
	
	if( (d==NULL) || (dref==NULL) || (e==NULL) || (eref==NULL)|| \
	    (ab==NULL) || (abref==NULL) || (q==NULL) || (qref==NULL)|| \
		(c==NULL) || (cref==NULL) || (work==NULL) || \
		(workref==NULL) || (pt==NULL) || (ptref==NULL) ){

       printf(" gbbrd_dcomplex_parameters object: malloc error. Exiting...\n");
       gbbrd_free();
	   exit (-1);
	}

    /* Initialization of input Buffers */

    for( j = 0; j < (ldab*n); j++ ) {
	   dtemp_r =(((double) rand()) / ((double) RAND_MAX) - 0.5);
	   dtemp_i =(((double) rand()) / ((double) RAND_MAX) - 0.5);
	   ab[j] = lapack_make_complex_double( dtemp_r, dtemp_i);
       abref[j] = lapack_make_complex_double( dtemp_r, dtemp_i);
    }

    for( j = 0; j < (ldc*ncc); j++ ) {
	   dtemp_r =(((double) rand()) / ((double) RAND_MAX) - 0.5);
	   dtemp_i =(((double) rand()) / ((double) RAND_MAX) - 0.5);
	   c[j] = lapack_make_complex_double( dtemp_r, dtemp_i);
       cref[j] = lapack_make_complex_double( dtemp_r, dtemp_i);
    }

	/** Below o/p buffers initialized to '0.0' **/  
    for( j = 0; j < temp; j++ ) {
       d[j] = 0.0;//lapack_make_complex_double( 0.0, 0.0);
       dref[j] = 0.0;//lapack_make_complex_double( 0.0, 0.0);
    }

    for( j = 0; j < temp-1; j++ ) {
       e[j] = 0.0;//lapack_make_complex_double( 0.0, 0.0);
       eref[j] = 0.0;//lapack_make_complex_double( 0.0, 0.0);
    }

    for( j = 0; j < (ldq*m); j++ ) {
       q[j] = lapack_make_complex_double( 0.0, 0.0);
       qref[j] = lapack_make_complex_double( 0.0, 0.0);
    }

    for( j = 0; j < (ldpt*n); j++ ) {
       pt[j] = lapack_make_complex_double( 0.0, 0.0);
       ptref[j] = lapack_make_complex_double( 0.0, 0.0);
    }

    for( j = 0; j < (2*max(m,n)); j++ ) {
       work[j] = lapack_make_complex_double( 0.0, 0.0);
       workref[j] = lapack_make_complex_double( 0.0, 0.0);
    }
	info = 0;
	inforef = 0;
	
}

/* Destructor definition **/
gbbrd_dcomplex_parameters:: ~gbbrd_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gbbrd_dcomplex_parameters object: destructor invoked. \n");
#endif
   gbbrd_free();
}

TEST(gbbrd, dgbbrd1) {
   /* LAPACKE DGESV prototype */
typedef int (*Fptr_NL_LAPACKE_dgbbrd)( int matrix_layout, char vect, int m,
                            int n, int ncc, int kl,
                            int ku, double* ab, int ldab, double* d,
                            double* e, double* q, int ldq, double* pt,
                            int ldpt, double* c, int ldc );
						   
   void *hModule, *dModule;
   Fptr_NL_LAPACKE_dgbbrd DGBBRD;		
   double diff;
   int ipiv_diff;
   
//   gbbrd_double_parameters  dgbbrd_obj (LAPACK_COL_MAJOR,'B', 50, 40, 40, 10, 10,
//                                          25, 50, 40, 50 );
   gbbrd_double_parameters  dgbbrd_obj (LAPACK_COL_MAJOR,'B', 50, 50, 50, 50, 50,
                                          150, 50, 50, 50 );
									  
   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);
  
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
  
   DGBBRD = (Fptr_NL_LAPACKE_dgbbrd)dlsym(hModule, "LAPACKE_dgbbrd");

   if (NULL == DGBBRD)
   {
   	  printf("Could not get the DGBBRD symbol from Netlib Lapack lib. Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }
   
   /* Compute the reference o/p by invoking Netlib-Lapack's API */
   dgbbrd_obj.inforef = DGBBRD( dgbbrd_obj.matrix_layout, dgbbrd_obj.vect, dgbbrd_obj.m,
                            dgbbrd_obj.n, dgbbrd_obj.ncc, dgbbrd_obj.kl,
                            dgbbrd_obj.ku, dgbbrd_obj.abref, dgbbrd_obj.ldab, dgbbrd_obj.dref,
                            dgbbrd_obj.eref, dgbbrd_obj.qref, dgbbrd_obj.ldq, dgbbrd_obj.ptref,
                            dgbbrd_obj.ldpt, dgbbrd_obj.cref, dgbbrd_obj.ldc );

   /* Compute libflame's Lapacke o/p  */
   dgbbrd_obj.info = LAPACKE_dgbbrd( dgbbrd_obj.matrix_layout, dgbbrd_obj.vect, dgbbrd_obj.m,
                            dgbbrd_obj.n, dgbbrd_obj.ncc, dgbbrd_obj.kl,
                            dgbbrd_obj.ku, dgbbrd_obj.ab, dgbbrd_obj.ldab, dgbbrd_obj.d,
                            dgbbrd_obj.e, dgbbrd_obj.q, dgbbrd_obj.ldq, dgbbrd_obj.pt,
                            dgbbrd_obj.ldpt, dgbbrd_obj.c, dgbbrd_obj.ldc );

	if( dgbbrd_obj.info < 0 ) {
		printf( "The i:%d th argument with LibFlame is wrong\n", dgbbrd_obj.info );
		//exit( -1 );
	}
	if( dgbbrd_obj.inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LapackE is wrong\n", dgbbrd_obj.inforef );
		//exit( -1 );
	}

		
  	diff =  computeDiff_d( (dgbbrd_obj.ldab*dgbbrd_obj.n), dgbbrd_obj.ab, dgbbrd_obj.abref );
	EXPECT_NEAR(0.0, diff, DOUBLE_DIFF_THRESHOLD);
}

TEST(gbbrd, sgbbrd1) {
   /* LAPACKE sgbbrd prototype */
typedef int (*Fptr_NL_LAPACKE_sgbbrd)( int matrix_layout, char vect, int m,
                            int n, int ncc, int kl,
                            int ku, float* ab, int ldab, float* d,
                            float* e, float* q, int ldq, float* pt,
                            int ldpt, float* c, int ldc );
						   
   void *hModule, *dModule;
   Fptr_NL_LAPACKE_sgbbrd SGBBRD;		
   float diff;
   int ipiv_diff;
   
   gbbrd_float_parameters  sgbbrd_obj (LAPACK_COL_MAJOR,'B', 50, 50, 50, 50, 40,
                                          150, 50, 50, 50 );
									  
   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);
  
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
  
   SGBBRD = (Fptr_NL_LAPACKE_sgbbrd)dlsym(hModule, "LAPACKE_sgbbrd");

   if (NULL == SGBBRD)
   {
   	  printf("Could not get the SGBBRD symbol from Netlib Lapack lib. Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }
   
   /* Compute the reference o/p by invoking Netlib-Lapack's API */
   sgbbrd_obj.inforef = SGBBRD( sgbbrd_obj.matrix_layout, sgbbrd_obj.vect, sgbbrd_obj.m,
                            sgbbrd_obj.n, sgbbrd_obj.ncc, sgbbrd_obj.kl,
                            sgbbrd_obj.ku, sgbbrd_obj.abref, sgbbrd_obj.ldab, sgbbrd_obj.dref,
                            sgbbrd_obj.eref, sgbbrd_obj.qref, sgbbrd_obj.ldq, sgbbrd_obj.ptref,
                            sgbbrd_obj.ldpt, sgbbrd_obj.cref, sgbbrd_obj.ldc );

   /* Compute libflame's Lapacke o/p  */
   sgbbrd_obj.info = LAPACKE_sgbbrd( sgbbrd_obj.matrix_layout, sgbbrd_obj.vect, sgbbrd_obj.m,
                            sgbbrd_obj.n, sgbbrd_obj.ncc, sgbbrd_obj.kl,
                            sgbbrd_obj.ku, sgbbrd_obj.ab, sgbbrd_obj.ldab, sgbbrd_obj.d,
                            sgbbrd_obj.e, sgbbrd_obj.q, sgbbrd_obj.ldq, sgbbrd_obj.pt,
                            sgbbrd_obj.ldpt, sgbbrd_obj.c, sgbbrd_obj.ldc );

	if( sgbbrd_obj.info < 0 ) {
		printf( "The i:%d th argument with LibFlame is wrong\n", sgbbrd_obj.info );
		//exit( -1 );
	}
	if( sgbbrd_obj.inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LapackE is wrong\n", sgbbrd_obj.inforef );
		//exit( -1 );
	}

		
  	diff =  computeDiff_s( (sgbbrd_obj.ldab*sgbbrd_obj.n), sgbbrd_obj.ab, sgbbrd_obj.abref );
	EXPECT_NEAR(0.0, diff, FLOAT_DIFF_THRESHOLD);
}


TEST(gbbrd, cgbbrd1) {
   /* LAPACKE CGBBRD prototype */
typedef int (*Fptr_NL_LAPACKE_cgbbrd)(int matrix_layout, char vect, lapack_int m,
                            lapack_int n, lapack_int ncc, lapack_int kl,
                            lapack_int ku, lapack_complex_float* ab,
                            lapack_int ldab, float* d, float* e,
                            lapack_complex_float* q, lapack_int ldq,
                            lapack_complex_float* pt, lapack_int ldpt,
                            lapack_complex_float* c, lapack_int ldc );

   void *hModule, *dModule;
   Fptr_NL_LAPACKE_cgbbrd CGBBRD;		
   float diff;
   int ipiv_diff;
   
   gbbrd_scomplex_parameters  cgbbrd_obj (LAPACK_COL_MAJOR,'B', 50, 50, 50, 50, 40,
                                          150, 50, 50, 50 );
									  
   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);
  
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
  
   CGBBRD = (Fptr_NL_LAPACKE_cgbbrd)dlsym(hModule, "LAPACKE_cgbbrd");

   if (NULL == CGBBRD)
   {
   	  printf("Could not get the CGBBRD symbol from Netlib Lapack lib. Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }
   
   /* Compute the reference o/p by invoking Netlib-Lapack's API */
   cgbbrd_obj.inforef = CGBBRD( cgbbrd_obj.matrix_layout, cgbbrd_obj.vect, cgbbrd_obj.m,
                            cgbbrd_obj.n, cgbbrd_obj.ncc, cgbbrd_obj.kl,
                            cgbbrd_obj.ku, cgbbrd_obj.abref, cgbbrd_obj.ldab, cgbbrd_obj.dref,
                            cgbbrd_obj.eref, cgbbrd_obj.qref, cgbbrd_obj.ldq, cgbbrd_obj.ptref,
                            cgbbrd_obj.ldpt, cgbbrd_obj.cref, cgbbrd_obj.ldc );

   /* Compute libflame's Lapacke o/p  */
   cgbbrd_obj.info = LAPACKE_cgbbrd( cgbbrd_obj.matrix_layout, cgbbrd_obj.vect, cgbbrd_obj.m,
                            cgbbrd_obj.n, cgbbrd_obj.ncc, cgbbrd_obj.kl,
                            cgbbrd_obj.ku, cgbbrd_obj.ab, cgbbrd_obj.ldab, cgbbrd_obj.d,
                            cgbbrd_obj.e, cgbbrd_obj.q, cgbbrd_obj.ldq, cgbbrd_obj.pt,
                            cgbbrd_obj.ldpt, cgbbrd_obj.c, cgbbrd_obj.ldc );

	if( cgbbrd_obj.info < 0 ) {
		printf( "The i:%d th argument with LibFlame is wrong\n", cgbbrd_obj.info );
		//exit( -1 );
	}
	if( cgbbrd_obj.inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LapackE is wrong\n", cgbbrd_obj.inforef );
		//exit( -1 );
	}

		
  	diff =  computeDiff_c( (cgbbrd_obj.ldab*cgbbrd_obj.n), cgbbrd_obj.ab, cgbbrd_obj.abref );
	EXPECT_NEAR(0.0, diff, FLOAT_DIFF_THRESHOLD);
}

TEST(gbbrd, zgbbrd1) {
   /* LAPACKE CGBBRD prototype */
typedef int (*Fptr_NL_LAPACKE_zgbbrd)(int matrix_layout, char vect, lapack_int m,
                            lapack_int n, lapack_int ncc, lapack_int kl,
                            lapack_int ku, lapack_complex_double* ab,
                            lapack_int ldab, double* d, double* e,
                            lapack_complex_double* q, lapack_int ldq,
                            lapack_complex_double* pt, lapack_int ldpt,
                            lapack_complex_double* c, lapack_int ldc );

   void *hModule, *dModule;
   Fptr_NL_LAPACKE_zgbbrd ZGBBRD;		
   double diff;
   int ipiv_diff;
   
   gbbrd_dcomplex_parameters  zgbbrd_obj (LAPACK_COL_MAJOR,'B', 50, 50, 50, 50, 40,
                                          150, 50, 50, 50 );
									  
   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);
  
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
  
   ZGBBRD = (Fptr_NL_LAPACKE_zgbbrd)dlsym(hModule, "LAPACKE_zgbbrd");

   if (NULL == ZGBBRD)
   {
   	  printf("Could not get the CGBBRD symbol from Netlib Lapack lib. Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }
   
   /* Compute the reference o/p by invoking Netlib-Lapack's API */
   zgbbrd_obj.inforef = ZGBBRD( zgbbrd_obj.matrix_layout, zgbbrd_obj.vect, zgbbrd_obj.m,
                            zgbbrd_obj.n, zgbbrd_obj.ncc, zgbbrd_obj.kl,
                            zgbbrd_obj.ku, zgbbrd_obj.abref, zgbbrd_obj.ldab, zgbbrd_obj.dref,
                            zgbbrd_obj.eref, zgbbrd_obj.qref, zgbbrd_obj.ldq, zgbbrd_obj.ptref,
                            zgbbrd_obj.ldpt, zgbbrd_obj.cref, zgbbrd_obj.ldc );

   /* Compute libflame's Lapacke o/p  */
   zgbbrd_obj.info = LAPACKE_zgbbrd( zgbbrd_obj.matrix_layout, zgbbrd_obj.vect, zgbbrd_obj.m,
                            zgbbrd_obj.n, zgbbrd_obj.ncc, zgbbrd_obj.kl,
                            zgbbrd_obj.ku, zgbbrd_obj.ab, zgbbrd_obj.ldab, zgbbrd_obj.d,
                            zgbbrd_obj.e, zgbbrd_obj.q, zgbbrd_obj.ldq, zgbbrd_obj.pt,
                            zgbbrd_obj.ldpt, zgbbrd_obj.c, zgbbrd_obj.ldc );

	if( zgbbrd_obj.info < 0 ) {
		printf( "The i:%d th argument with LibFlame is wrong\n", zgbbrd_obj.info );
		//exit( -1 );
	}
	if( zgbbrd_obj.inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LapackE is wrong\n", zgbbrd_obj.inforef );
		//exit( -1 );
	}

		
  	diff =  computeDiff_z( (zgbbrd_obj.ldab*zgbbrd_obj.n), zgbbrd_obj.ab, zgbbrd_obj.abref );
	EXPECT_NEAR(0.0, diff, DOUBLE_DIFF_THRESHOLD);
}