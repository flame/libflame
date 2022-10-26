#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define gebal_free() \
		free(A);	\
		free(Aref); \
		free(scale); \
		free(scaleref); 

/* Begin float_common_parameters  class definition */
class gebal_float_parameters{

   public:
	/* INPUT PARAMETERS */
	int matrix_layout; // Matrix Layout
	lapack_int n; // Number of Column
	float* A;
	lapack_int lda; // lda = n
	char job; // job = 'N'
	/*Output parametes */
	float* scale;
	lapack_int ilo, ihi;
	float *Aref,*scaleref;

	/*Return Values */
	int info, inforef;

   public:
      gebal_float_parameters (int matrix_layout, char job, lapack_int n, lapack_int lda);
      ~gebal_float_parameters ();

};

/* Constructor definition  float_common_parameters */
gebal_float_parameters:: gebal_float_parameters (int matrix_layout_i, char job_i, lapack_int n_i,lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n  = n_i;
	lda = lda_i;
	job = job_i;
	ilo = ihi = 0;

	lapacke_gtest_alloc_float_buffer_pair(&A, &Aref, (lda*n));
	lapacke_gtest_alloc_float_buffer_pair(&scale, &scaleref, n);	
	if ((A==NULL) || (Aref==NULL) || 
		(scale == NULL) || (scaleref == NULL)){
		printf("error of memory allocation. Exiting ...\n");
		gebal_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand(A, Aref, (lda*n));
	for(i = 0; i < n; i++)
	{
		scale[i] = 0;
		scaleref[i] = scale[i];
	}

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
gebal_float_parameters :: ~gebal_float_parameters ()
{
   /* De-Allocate memory for the input matrices */
   gebal_free();
} 

/* Begin double_common_parameters  class definition */
class gebal_double_parameters{

   public:
	/* INPUT PARAMETERS */
	int matrix_layout; // Matrix Layout
	lapack_int n; // Number of Column
	double* A;
	lapack_int lda; // lda = n
	char job; // job = 'N'
	/*Output parametes */
	double* scale;
	lapack_int ilo,  ihi;
	double *Aref,*scaleref;

	/*Return Values */
	int info, inforef;

   public:
      gebal_double_parameters (int matrix_layout, char job, lapack_int n, lapack_int lda);
      ~gebal_double_parameters ();

};

/* Constructor definition  double_common_parameters */
gebal_double_parameters:: gebal_double_parameters (int matrix_layout_i, char job_i, lapack_int n_i,lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n  = n_i;
	lda = lda_i;
	job = job_i;
	ilo = 0;
	ihi = 0;

	lapacke_gtest_alloc_double_buffer_pair(&A, &Aref, (lda*n));
	lapacke_gtest_alloc_double_buffer_pair(&scale, &scaleref, n);	
	if ((A==NULL) || (Aref==NULL) || 
		(scale == NULL) || (scaleref == NULL)){
		printf("error of memory allocation. Exiting ...\n");
		gebal_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand(A, Aref, (lda*n));
	for(i = 0; i < n; i++)
	{
		scale[i] = 0;
		scaleref[i] = scale[i];
	}


} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
gebal_double_parameters :: ~gebal_double_parameters ()
{
   /* De-Allocate memory for the input matrices */
   gebal_free();  
}


/* Begin complex_common_parameters  class definition */
class gebal_complex_parameters{

   public:
	/* INPUT PARAMETERS */
	int matrix_layout; // Matrix Layout
	lapack_int n; // Number of Column
	lapack_complex_float* A;
	lapack_int lda; // lda = n
	char job; // job = 'N'
	/*Output parametes */
	float *scale, *scaleref;
	lapack_int  ilo, ihi;
	lapack_complex_float *Aref;

	/*Return Values */
	int info, inforef;

   public:
      gebal_complex_parameters (int matrix_layout, char job, lapack_int n, lapack_int lda);
      ~gebal_complex_parameters ();

};

/* Constructor definition  complex_common_parameters */
gebal_complex_parameters:: gebal_complex_parameters (int matrix_layout_i, char job_i, lapack_int n_i,lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n  = n_i;
	lda = lda_i;
	job = job_i;
	ilo = 0;
	ihi = 0;

	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, (lda*n));
	lapacke_gtest_alloc_float_buffer_pair(&scale, &scaleref, n);	
	if ((A==NULL) || (Aref==NULL) || 
		(scale == NULL) || (scaleref == NULL)){
		printf("error of memory allocation. Exiting ...\n");
		gebal_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand(A, Aref, (lda*n));
	for(i = 0; i < n; i++)
	{
		scale[i] = 0;
		scaleref[i] = scale[i];
	}

} /* end of Constructor  */

/* Destructor definition  'complex_common_parameters' class  */
gebal_complex_parameters :: ~gebal_complex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   gebal_free();
}


/* Begin complex_double_common_parameters  class definition */
class gebal_complex_double_parameters{

   public:
	/* INPUT PARAMETERS */
	int matrix_layout; // Matrix Layout
	lapack_int n; // Number of Column
	lapack_complex_double * A;
	lapack_int lda; // lda = n
	char job; // job = 'N'
	/*Output parametes */
	double *scale, *scaleref;
	lapack_int ilo, ihi;
	lapack_complex_double *Aref;

	/*Return Values */
	int info, inforef;

   public:
      gebal_complex_double_parameters (int matrix_layout, char job, lapack_int n, lapack_int lda);
      ~gebal_complex_double_parameters ();

};

/* Constructor definition  complex_double_common_parameters */
gebal_complex_double_parameters:: gebal_complex_double_parameters (int matrix_layout_i, char job_i, lapack_int n_i,lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n  = n_i;
	lda = lda_i;
	job = job_i;
	ilo = 0;
	ihi = 0;

	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, (lda*n));
	lapacke_gtest_alloc_double_buffer_pair(&scale, &scaleref, n);	
	if ((A==NULL) || (Aref==NULL) || 
		(scale == NULL) || (scaleref == NULL)){
		printf("error of memory allocation. Exiting ...\n");
		gebal_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand(A, Aref, (lda*n));
	for(i = 0; i < n; i++)
	{
		scale[i] = 0;
		scaleref[i] = scale[i];
	}

} /* end of Constructor  */

/* Destructor definition  'complex_double_common_parameters' class  */
gebal_complex_double_parameters :: ~gebal_complex_double_parameters ()
{
   /* De-Allocate memory for the input matrices */
   gebal_free();
}

/* TESTS --------------------------------------------------------------------------------------*/
TEST(gebal, sgebal) {

   /* LAPACKE DGETRF prototype */
   typedef int (*Fptr_NL_LAPACKE_sgebal) ( int matrix_layout, char job, lapack_int n, float *A, lapack_int lda, 
											lapack_int* ilo, lapack_int* ihi, float* scale);
						   
   Fptr_NL_LAPACKE_sgebal sgebal;		
   float  diff;
   void *hModule, *dModule;

   hModule= NULL; dModule = NULL;   
  
   gebal_float_parameters sgebal_obj(LAPACK_COL_MAJOR, 'B', 256,256); 
   dModule = dlopen("Netlib_lapack_ref_lib/libblas.so.3.9.0", RTLD_NOW | RTLD_GLOBAL);
   
   hModule = dlopen("Netlib_lapack_ref_lib/liblapacke.so.3.9.0", RTLD_NOW);
   
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   
   sgebal = (Fptr_NL_LAPACKE_sgebal)dlsym(hModule, "LAPACKE_sgebal");
   if (NULL == sgebal)
   {
   	  printf("Could not get the symbol. Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (0);
   }
      /* Compute libflame's Lapacke o/p  */
      sgebal_obj.info    = LAPACKE_sgebal( LAPACK_COL_MAJOR, sgebal_obj.job, sgebal_obj.n,sgebal_obj.A, sgebal_obj.lda, 
											&sgebal_obj.ilo, &sgebal_obj.ihi, sgebal_obj.scale);
	  
      /* Compute the reference o/p by invoking Netlib-Lapack's API */	
      sgebal_obj.inforef = sgebal( LAPACK_COL_MAJOR, sgebal_obj.job, sgebal_obj.n,sgebal_obj.Aref, sgebal_obj.lda, 
		     			 &sgebal_obj.ilo, &sgebal_obj.ihi, sgebal_obj.scaleref);

   	/* Check for the exact singularity */
	if( sgebal_obj.info > 0 ) {
		printf( "The diagonal element of the triangular factor of A,\n" );
		printf( "U(%i,%i) is zero, so that A is singular;\n", sgebal_obj.info, sgebal_obj.info );
		printf( "the solution could not be computed. Exiting...  \n" );

		exit( 1 );
	}
	if( sgebal_obj.inforef > 0 ) {
		printf( "The diagonal element of the triangular factor of Aref,\n" );
		printf( "U(%i,%i) is zero, so that Aref is singular;\n", sgebal_obj.inforef, sgebal_obj.inforef );
		printf( "the solution could not be computed. Exiting...  \n" );

		exit( 1 );
	}

	/* Compute Difference in C and CPP buffer */
	diff =  computeDiff_s( sgebal_obj.n*sgebal_obj.lda,  sgebal_obj.A, sgebal_obj.Aref );

	/*
	FLA_Nrm2_external( dgetrf_obj.A, dgetrf_obj.norm );
        FLA_Nrm2_external( dgetrf_obj.Aref, dgetrf_obj. normref );
	diff = FLA_abs( dgetrf_obj. norm,  dgetrf_obj. normref);
	*/
    
    	//EXPECT_EQ (0, diff);
	//EXPECT_FLOAT_EQ (0.0, diff);
	EXPECT_NEAR(0.0, diff, 1.0);
}
/* Test for DOuble */
TEST(gebal, dgebal) {

   /* LAPACKE DGETRF prototype */
   typedef int (*Fptr_NL_LAPACKE_dgebal) ( int matrix_layout, char job, lapack_int n, double *A, lapack_int lda, 
											lapack_int* ilo, lapack_int* ihi, double* scale);
						   
   Fptr_NL_LAPACKE_dgebal dgebal;		
   double diff;
   void *hModule, *dModule;

   hModule= NULL; dModule = NULL;   
  
   gebal_double_parameters dgebal_obj(LAPACK_COL_MAJOR, 'B', 256,256); 
   dModule = dlopen("Netlib_lapack_ref_lib/libblas.so.3.9.0", RTLD_NOW | RTLD_GLOBAL);
   
   hModule = dlopen("Netlib_lapack_ref_lib/liblapacke.so.3.9.0", RTLD_NOW);
   
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   
   dgebal = (Fptr_NL_LAPACKE_dgebal)dlsym(hModule, "LAPACKE_dgebal");
   if (NULL == dgebal)
   {
   	  printf("Could not get the symbol. Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (0);
   }
      /* Compute libflame's Lapacke o/p  */
      dgebal_obj.info    = LAPACKE_dgebal( LAPACK_COL_MAJOR, dgebal_obj.job, dgebal_obj.n, dgebal_obj.A, dgebal_obj.lda, 
											&dgebal_obj.ilo, &dgebal_obj.ihi, dgebal_obj.scale);
	  
      /* Compute the reference o/p by invoking Netlib-Lapack's API */	
      dgebal_obj.inforef = dgebal( LAPACK_COL_MAJOR, dgebal_obj.job, dgebal_obj.n,dgebal_obj.Aref, dgebal_obj.lda, 
		     			 &dgebal_obj.ilo, &dgebal_obj.ihi, dgebal_obj.scaleref);

   	/* Check for the exact singularity */
	if( dgebal_obj.info > 0 ) {
		printf( "The diagonal element of the triangular factor of A,\n" );
		printf( "U(%i,%i) is zero, so that A is singular;\n", dgebal_obj.info, dgebal_obj.info );
		printf( "the solution could not be computed. Exiting...  \n" );

		exit( 1 );
	}
	if( dgebal_obj.inforef > 0 ) {
		printf( "The diagonal element of the triangular factor of Aref,\n" );
		printf( "U(%i,%i) is zero, so that Aref is singular;\n", dgebal_obj.inforef, dgebal_obj.inforef );
		printf( "the solution could not be computed. Exiting...  \n" );

		exit( 1 );
	}

	/* Compute Difference in C and CPP buffer */
	diff =  computeDiff_d( dgebal_obj.n*dgebal_obj.lda,  dgebal_obj.A, dgebal_obj.Aref );

	/*
	FLA_Nrm2_external( dgetrf_obj.A, dgetrf_obj.norm );
        FLA_Nrm2_external( dgetrf_obj.Aref, dgetrf_obj. normref );
	diff = FLA_abs( dgetrf_obj. norm,  dgetrf_obj. normref);
	*/
    
    	//EXPECT_EQ (0, diff);
	//EXPECT_FLOAT_EQ (0.0, diff);
	EXPECT_NEAR(0.0, diff, 1.0);
}

TEST(gebal, cgebal) {

   /* LAPACKE DGETRF prototype */
   typedef int (*Fptr_NL_LAPACKE_cgebal) ( int matrix_layout, char job, lapack_int n, lapack_complex_float* A,
											lapack_int lda, lapack_int* ilo, lapack_int* ihi, float* scale);
						   
   Fptr_NL_LAPACKE_cgebal cgebal;		
   float  diff;
   void *hModule, *dModule;

   hModule= NULL; dModule = NULL;   
  
   gebal_complex_parameters cgebal_obj(LAPACK_COL_MAJOR, 'B', 256,256); 
   dModule = dlopen("Netlib_lapack_ref_lib/libblas.so.3.9.0", RTLD_NOW | RTLD_GLOBAL);
   
   hModule = dlopen("Netlib_lapack_ref_lib/liblapacke.so.3.9.0", RTLD_NOW);
   
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   
   cgebal = (Fptr_NL_LAPACKE_cgebal)dlsym(hModule, "LAPACKE_cgebal");
   if (NULL == cgebal)
   {
   	  printf("Could not get the symbol. Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (0);
   }
      /* Compute libflame's Lapacke o/p  */
      cgebal_obj.info    = LAPACKE_cgebal( LAPACK_COL_MAJOR, cgebal_obj.job, cgebal_obj.n,cgebal_obj.A, cgebal_obj.lda, 
											&cgebal_obj.ilo, &cgebal_obj.ihi, cgebal_obj.scale);
	  
      /* Compute the reference o/p by invoking Netlib-Lapack's API */	
      cgebal_obj.inforef = cgebal( LAPACK_COL_MAJOR, cgebal_obj.job, cgebal_obj.n,cgebal_obj.Aref, cgebal_obj.lda, 
		     			 &cgebal_obj.ilo, &cgebal_obj.ihi, cgebal_obj.scaleref);

   	/* Check for the exact singularity */
	if( cgebal_obj.info > 0 ) {
		printf( "The diagonal element of the triangular factor of A,\n" );
		printf( "U(%i,%i) is zero, so that A is singular;\n", cgebal_obj.info, cgebal_obj.info );
		printf( "the solution could not be computed. Exiting...  \n" );

		exit( 1 );
	}
	if( cgebal_obj.inforef > 0 ) {
		printf( "The diagonal element of the triangular factor of Aref,\n" );
		printf( "U(%i,%i) is zero, so that Aref is singular;\n", cgebal_obj.inforef, cgebal_obj.inforef );
		printf( "the solution could not be computed. Exiting...  \n" );

		exit( 1 );
	}

	/* Compute Difference in C and CPP buffer */
	diff =  computeDiff_c( cgebal_obj.n*cgebal_obj.lda,  cgebal_obj.A, cgebal_obj.Aref );

	/*
	FLA_Nrm2_external( dgetrf_obj.A, dgetrf_obj.norm );
        FLA_Nrm2_external( dgetrf_obj.Aref, dgetrf_obj. normref );
	diff = FLA_abs( dgetrf_obj. norm,  dgetrf_obj. normref);
	*/
    
    	//EXPECT_EQ (0, diff);
	//EXPECT_FLOAT_EQ (0.0, diff);
	EXPECT_NEAR(0.0, diff, 1.0);
}
/***********TEST for DOuble complex **************/
TEST(gebal, zgebal) {

   /* LAPACKE DGETRF prototype */
   typedef int (*Fptr_NL_LAPACKE_zgebal) ( int matrix_layout, char job, lapack_int n, lapack_complex_double *A, lapack_int lda, 
											lapack_int* ilo, lapack_int* ihi, double* scale);
						   
   Fptr_NL_LAPACKE_zgebal zgebal;		
   float  diff;
   void *hModule, *dModule;

   hModule= NULL; dModule = NULL;   
  
   gebal_complex_double_parameters zgebal_obj(LAPACK_COL_MAJOR, 'N', 256,256); 
   dModule = dlopen("Netlib_lapack_ref_lib/libblas.so.3.9.0", RTLD_NOW | RTLD_GLOBAL);
   
   hModule = dlopen("Netlib_lapack_ref_lib/liblapacke.so.3.9.0", RTLD_NOW);
   
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   
   zgebal = (Fptr_NL_LAPACKE_zgebal)dlsym(hModule, "LAPACKE_zgebal");
   if (NULL == zgebal)
   {
   	  printf("Could not get the symbol. Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (0);
   }
      /* Compute libflame's Lapacke o/p  */
      zgebal_obj.info    = LAPACKE_zgebal( LAPACK_COL_MAJOR, zgebal_obj.job, zgebal_obj.n,zgebal_obj.A, zgebal_obj.lda, 
											&zgebal_obj.ilo, &zgebal_obj.ihi, zgebal_obj.scale);
	  
      /* Compute the reference o/p by invoking Netlib-Lapack's API */	
      zgebal_obj.inforef = zgebal( LAPACK_COL_MAJOR, zgebal_obj.job, zgebal_obj.n,zgebal_obj.Aref, zgebal_obj.lda, 
		     			 &zgebal_obj.ilo, &zgebal_obj.ihi, zgebal_obj.scaleref);

   	/* Check for the exact singularity */
	if( zgebal_obj.info > 0 ) {
		printf( "The diagonal element of the triangular factor of A,\n" );
		printf( "U(%i,%i) is zero, so that A is singular;\n", zgebal_obj.info, zgebal_obj.info );
		printf( "the solution could not be computed. Exiting...  \n" );

		exit( 1 );
	}
	if( zgebal_obj.inforef > 0 ) {
		printf( "The diagonal element of the triangular factor of Aref,\n" );
		printf( "U(%i,%i) is zero, so that Aref is singular;\n", zgebal_obj.inforef, zgebal_obj.inforef );
		printf( "the solution could not be computed. Exiting...  \n" );

		exit( 1 );
	}

	/* Compute Difference in C and CPP buffer */
	diff =  computeDiff_z( zgebal_obj.n*zgebal_obj.lda,  zgebal_obj.A, zgebal_obj.Aref );

	/*
	FLA_Nrm2_external( dgetrf_obj.A, dgetrf_obj.norm );
        FLA_Nrm2_external( dgetrf_obj.Aref, dgetrf_obj. normref );
	diff = FLA_abs( dgetrf_obj. norm,  dgetrf_obj. normref);
	*/
    
    	//EXPECT_EQ (0, diff);
	//EXPECT_FLOAT_EQ (0.0, diff);
	EXPECT_NEAR(0.0, diff, 1.0);
}





