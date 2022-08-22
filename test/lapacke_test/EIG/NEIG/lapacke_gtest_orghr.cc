#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"


/* Begin float_common_parameters  class definition */
class orghr_float_parameters{

   public:
	/* INPUT PARAMETERS */
	int matrix_layout; // Matrix Layout
	lapack_int n; // Number of Column
	lapack_int ilo; //
	lapack_int ihi; //
	float* A;       
	lapack_int lda; // lda = n

	/*Output parametes */
	float* tau;
	float *Aref,*tauref;

	/*Return Values */
	int info, inforef;

   public:
      orghr_float_parameters (int matrix_layout, lapack_int n, lapack_int ilo,
		      		lapack_int ihi, lapack_int lda);
      ~orghr_float_parameters ();

};

/* Constructor definition  float_common_parameters */
orghr_float_parameters:: orghr_float_parameters (int matrix_layout_i, lapack_int n_i, lapack_int ilo_i,
	       				lapack_int ihi_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n  = n_i;
	lda = lda_i;
	ilo = ilo_i;
	ihi = ihi_i;

	A = (float *)malloc(lda*n*sizeof(float)) ;
	tau = (float *)malloc((n-1)*sizeof(float)) ;
	if ((A==NULL) || (tau==NULL)){
		printf("error of memory allocation. Exiting ...\n");
		free(A); free(tau);
		exit(0);
	}

	/* Allocation of memory for capturing reference o/ps */
	Aref = (float *)malloc(lda*n*sizeof(float)) ;
	tauref = (float *)malloc((n-1)*sizeof(float)) ;
	if ((Aref==NULL) || (tauref==NULL)){
		printf("error of memory allocation. Exiting ...\n");
		free(Aref); free(tauref);
		exit(0);
	}
	/* Initialization of input matrices */
	for( i = 0; i <(lda*n); i++ ) {
		A[i] = ((float) rand()) / ((float) RAND_MAX) - 0.5;
		Aref[i] = A[i];
	}

	for(i=0;i<(n-1);i++) {
		tau[i] = ((float) rand()) / ((float) RAND_MAX) - 0.5;
		tauref[i] = tau[i];
	}

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
orghr_float_parameters :: ~orghr_float_parameters ()
{
   /* De-Allocate memory for the input matrices */
   if (A!=NULL){
      free(A);
   }
   if (Aref!=NULL){
      free(Aref);
   }
   if (tau!=NULL){
      free(tau);
   }
   if (tauref!=NULL){
      free(tauref);
   }

}

/* Begin double_common_parameters  class definition */
class orghr_double_parameters {

public:
	/*input parametes*/
	int matrix_layout;
	lapack_int n;
	lapack_int ilo;
	lapack_int ihi;
	double* A;
	lapack_int lda;
	/*output parametes*/
	double* tau;
	double* Aref, * tauref;
	/*return values*/
	int info, inforef;

public:
	orghr_double_parameters(int matrix_layout, lapack_int n, lapack_int ilo,
							lapack_int ihi, lapack_int lda);
	~orghr_double_parameters();

};

/* Constructor definition  orghr_double_parameters */
orghr_double_parameters::orghr_double_parameters( int matrix_layout_i, lapack_int n_i, lapack_int ilo_i,
							lapack_int ihi_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i; // set test matrix size
	lda = lda_i;
	ilo = ilo_i;
	ihi = ihi_i;

	A = (double*)malloc(lda * n * sizeof(double));
	tau = (double*)malloc((n - 1) * sizeof(double));
	if ((A == NULL) || (tau == NULL)) {
		printf("error of memory allocation. Exiting ...\n");
		free(A); free(tau);
		exit(0);
	}

	/* Allocation of memory for capturing reference o/ps */
	Aref = (double*)malloc(lda * n * sizeof(double));
	tauref = (double*)malloc((n - 1) * sizeof(double));
	if ((Aref == NULL) || (tauref == NULL)) {
		printf("error of memory allocation. Exiting ...\n");
		free(Aref); free(tauref);
		exit(0);
	}
	/* Initialization of input matrices */
	for (i = 0; i < (lda * n); i++) {
		A[i] = ((double)rand()) / ((double)RAND_MAX) - 0.5;
		Aref[i] = A[i];
	}

	for (i = 0; i < (n - 1); i++) {
		tau[i] = ((double)rand()) / ((double)RAND_MAX) - 0.5;
		tauref[i] = tau[i];
	}

} /* end of Constructor  */

/* Destructor definition  'gethrd_double_parameters' class  */
orghr_double_parameters :: ~orghr_double_parameters()
{
	/* De-Allocate memory for the input matrices */
	if (A != NULL) {
		free(A);
	}
	if (Aref != NULL) {
		free(Aref);
	}
	if (tau != NULL) {
		free(tau);
	}
	if (tauref != NULL) {
		free(tauref);
	}

}
#if 0
/* Begin orghr_complex  class definition */
class orghr_complex_parameters {
public:
	/*input parametes*/
	int matrix_layout;
	lapack_int n;
	lapack_int ilo;
	lapack_int ihi;
	lapack_complex_float * A;
	lapack_int lda;
	/*output parametes*/
	lapack_complex_float * tau;
	lapack_complex_float * Aref, * tauref;
	/*return values*/
	int info, inforef;

public:
	orghr_complex_parameters(int matrix_layout, lapack_int n, lapack_int ilo,
								lapack_int ihi, lapack_int lda);
	~orghr_complex_parameters();

};

/* Constructor definition  orghr_complex_parameters */
orghr_complex_parameters::orghr_complex_parameters(int matrix_layout_i, lapack_int n_i, lapack_int ilo_i,
							lapack_int ihi_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i; // set test matrix size
	lda = lda_i;
	ilo = ilo_i;
	ihi = ihi_i;

	A = (lapack_complex_float*)malloc(lda * n * sizeof(lapack_complex_float));
	tau = (lapack_complex_float*)malloc((n - 1) * sizeof(lapack_complex_float));
	if ((A == NULL) || (tau == NULL)) {
		printf("error of memory allocation. Exiting ...\n");
		free(A); free(tau);
		exit(0);
	}

	/* Allocation of memory for capturing reference o/ps */
	Aref = (lapack_complex_float*)malloc(lda * n * sizeof(lapack_complex_float));
	tauref = (lapack_complex_float*)malloc((n - 1) * sizeof(lapack_complex_float));
	if ((Aref == NULL) || (tauref == NULL)) {
		printf("error of memory allocation. Exiting ...\n");
		free(Aref); free(tauref);
		exit(0);
	}
	/* Initialization of input matrices */
	for (i = 0; i < (lda * n); i++) {
		A[i] = ((lapack_complex_float)rand()) / ((lapack_complex_float)RAND_MAX) - 0.5;
		Aref[i] = A[i];
	}

	for (i = 0; i < (n - 1); i++) {
		tau[i] = ((lapack_complex_float)rand()) / ((lapack_complex_float)RAND_MAX) - 0.5;
		tauref[i] = tau[i];
	}

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
orghr_complex_parameters :: ~orghr_complex_parameters()
{
	/* De-Allocate memory for the input matrices */
	if (A != NULL) {
		free(A);
	}
	if (Aref != NULL) {
		free(Aref);
	}
	if (tau != NULL) {
		free(tau);
	}
	if (tauref != NULL) {
		free(tauref);
	}

}

/* Begin orghr_complex_double_parameters  class definition */
class orghr_complex_double_parameters {

public:
/*input parametes*/
	int matrix_layout;
	lapack_int n;
	lapack_int ilo;
	lapack_int ihi;
	lapack_complex_double* A;
	lapack_int lda;
	/*output parametes*/
	lapack_complex_double* tau;
	lapack_complex_double* Aref, * tauref;
	/*return values*/
	int info, inforef;

public:
	orghr_complex_double_parameters(int matrix_layout, lapack_int n, 
									lapack_int ilo, lapack_int ihi, lapack_int lda);
	~orghr_complex_double_parameters();

};

/* Constructor definition  float_common_parameters */
orghr_complex_double_parameters::orghr_complex_double_parameters(int matrix_layout_i, lapack_int n_i, lapack_int ilo_i,
							lapack_int ihi_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i; // set test matrix size
	lda = lda_i;
	ilo = ilo_i;
	ihi = ihi_i;

	A = (lapack_complex_double*)malloc(lda * n * sizeof(lapack_complex_double));
	tau = (lapack_complex_double*)malloc((n - 1) * sizeof(lapack_complex_double));
	if ((A == NULL) || (tau == NULL)) {
		printf("error of memory allocation. Exiting ...\n");
		free(A); free(tau);
		exit(0);
	}

	/* Allocation of memory for capturing reference o/ps */
	Aref = (lapack_complex_double*)malloc(lda * n * sizeof(lapack_complex_double));
	tauref = (lapack_complex_double*)malloc((n - 1) * sizeof(lapack_complex_double));
	if ((Aref == NULL) || (tauref == NULL)) {
		printf("error of memory allocation. Exiting ...\n");
		free(Aref); free(tauref);
		exit(0);
	}
	/* Initialization of input matrices */
	for (i = 0; i < (lda * n); i++) {
		A[i] = ((lapack_complex_double)rand()) / ((lapack_complex_double)RAND_MAX) - 0.5;
		Aref[i] = A[i];
	}

	for (i = 0; i < (n - 1); i++) {
		tau[i] = ((lapack_complex_double)rand()) / ((lapack_complex_double)RAND_MAX) - 0.5;
		tauref[i] = tau[i];
	}

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
orghr_complex_double_parameters :: ~orghr_complex_double_parameters()
{
	/* De-Allocate memory for the input matrices */
	if (A != NULL) {
		free(A);
	}
	if (Aref != NULL) {
		free(Aref);
	}
	if (tau != NULL) {
		free(tau);
	}
	if (tauref != NULL) {
		free(tauref);
	}

}
#endif

/* TESTS --------------------------------------------------------------------------------------*/
TEST(orghr, sorghr1) {

   /* LAPACKE DGETRF prototype */
   typedef int (*Fptr_NL_LAPACKE_sorghr) ( int matrix_layout, lapack_int n, lapack_int ilo,
                           lapack_int ihi, float *A, lapack_int lda, float* tau );
						   
   Fptr_NL_LAPACKE_sorghr sorghr;		
   float  diff;
   void *hModule, *dModule;

   hModule= NULL; dModule = NULL;
   //ilo = 10;
   //ihi = 10;
  
   orghr_float_parameters sorghr_obj(LAPACK_COL_MAJOR, 128, 10, 10, 128);
   
   dModule = dlopen("Netlib_lapack_ref_lib/libblas.so.3.9.0", RTLD_NOW | RTLD_GLOBAL);
   
   hModule = dlopen("Netlib_lapack_ref_lib/liblapacke.so.3.9.0", RTLD_NOW);
   
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   
   sorghr = (Fptr_NL_LAPACKE_sorghr)dlsym(hModule, "LAPACKE_sorghr");
   if (NULL == sorghr)
   {
   	  printf("Could not get the symbol. Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (0);
   }
      /* Compute libflame's Lapacke o/p  */
      sorghr_obj.info    = LAPACKE_sorghr( LAPACK_COL_MAJOR, sorghr_obj.n, sorghr_obj.ilo,
	                                  sorghr_obj.ihi, sorghr_obj.A, sorghr_obj.lda, sorghr_obj.tau);
	  
      /* Compute the reference o/p by invoking Netlib-Lapack's API */	
      sorghr_obj.inforef = sorghr( LAPACK_COL_MAJOR, sorghr_obj.n, sorghr_obj.ilo, sorghr_obj.ihi,
		     			 sorghr_obj.Aref, sorghr_obj.lda, sorghr_obj.tauref);

   	/* Check for the exact singularity */
	if( sorghr_obj.info > 0 ) {
		printf( "The diagonal element of the triangular factor of A,\n" );
		printf( "U(%i,%i) is zero, so that A is singular;\n", sorghr_obj.info, sorghr_obj.info );
		printf( "the solution could not be computed. Exiting...  \n" );

		exit( 1 );
	}
	if( sorghr_obj.inforef > 0 ) {
		printf( "The diagonal element of the triangular factor of Aref,\n" );
		printf( "U(%i,%i) is zero, so that Aref is singular;\n", sorghr_obj.inforef, sorghr_obj.inforef );
		printf( "the solution could not be computed. Exiting...  \n" );

		exit( 1 );
	}

	/* Compute Difference in C and CPP buffer */
	diff =  computeDiff_s( sorghr_obj.n*sorghr_obj.lda,  sorghr_obj.A, sorghr_obj.Aref );

	/*
	FLA_Nrm2_external( dgetrf_obj.A, dgetrf_obj.norm );
        FLA_Nrm2_external( dgetrf_obj.Aref, dgetrf_obj. normref );
	diff = FLA_abs( dgetrf_obj. norm,  dgetrf_obj. normref);
	*/
    
    	//EXPECT_EQ (0, diff);
	//EXPECT_FLOAT_EQ (0.0, diff);
	EXPECT_NEAR(0.0, diff, 1.0);
}


TEST(orghr, dorghr1) {

	/* LAPACKE Dorghr prototype */
	typedef int (*Fptr_NL_LAPACKE_dorghr) (int matrix_layout, lapack_int n, lapack_int ilo,
		lapack_int ihi, double* A, lapack_int lda, double* tau);

	Fptr_NL_LAPACKE_dorghr dorghr;
	double  diff;
	void* hModule, * dModule;

	hModule = NULL; dModule = NULL;
//	ilo = 10;
	//ihi = 10;

	orghr_double_parameters dorghr_obj(LAPACK_COL_MAJOR, 128, 10, 10, 128);

	dModule = dlopen("Netlib_lapack_ref_lib/libblas.so.3.9.0", RTLD_NOW | RTLD_GLOBAL);

	hModule = dlopen("Netlib_lapack_ref_lib/liblapacke.so.3.9.0", RTLD_NOW);

	if ((NULL == hModule) || (NULL == dModule))
	{
		printf("Load Library failed. Exiting ....\n");
		exit(0);
	}

	dorghr = (Fptr_NL_LAPACKE_dorghr)dlsym(hModule, "LAPACKE_dorghr");
	if (NULL == dorghr)
	{
		printf("Could not get the symbol. Exiting...\n");
		dlclose(hModule);
		dlclose(dModule);
		exit(0);
	}
	/* Compute libflame's Lapacke o/p  */
	dorghr_obj.info = LAPACKE_dorghr(LAPACK_COL_MAJOR, dorghr_obj.n, dorghr_obj.ilo,
		dorghr_obj.ihi, dorghr_obj.A, dorghr_obj.lda, dorghr_obj.tau);

	/* Compute the reference o/p by invoking Netlib-Lapack's API */
	dorghr_obj.inforef = dorghr(LAPACK_COL_MAJOR, dorghr_obj.n, dorghr_obj.ilo, dorghr_obj.ihi,
		dorghr_obj.Aref, dorghr_obj.lda, dorghr_obj.tauref);

	/* Check for the exact singularity */
	if (dorghr_obj.info > 0) {
		printf("The diagonal element of the triangular factor of A,\n");
		printf("U(%i,%i) is zero, so that A is singular;\n", dorghr_obj.info, dorghr_obj.info);
		printf("the solution could not be computed. Exiting...  \n");

		exit(1);
	}
	if (dorghr_obj.inforef > 0) {
		printf("The diagonal element of the triangular factor of Aref,\n");
		printf("U(%i,%i) is zero, so that Aref is singular;\n", dorghr_obj.inforef, dorghr_obj.inforef);
		printf("the solution could not be computed. Exiting...  \n");

		exit(1);
	}

	/* Compute Difference in C and CPP buffer */
	diff = computeDiff_d(dorghr_obj.n * dorghr_obj.lda, dorghr_obj.A, dorghr_obj.Aref);

	/*
	FLA_Nrm2_external( dgetrf_obj.A, dgetrf_obj.norm );
		FLA_Nrm2_external( dgetrf_obj.Aref, dgetrf_obj. normref );
	diff = FLA_abs( dgetrf_obj. norm,  dgetrf_obj. normref);
	*/

	//EXPECT_EQ (0, diff);
//EXPECT_FLOAT_EQ (0.0, diff);
	EXPECT_NEAR(0.0, diff, 1.0);
}
/* Test for Complex */
#if 0
TEST(orghr, corghr1) {

   /* LAPACKE DGETRF prototype */
   typedef int (*Fptr_NL_LAPACKE_corghr) ( int matrix_layout, lapack_int n, lapack_int ilo,
                           lapack_int ihi, lapack_complex_float* A, lapack_int lda, lapack_complex_float* tau );
						   
   Fptr_NL_LAPACKE_corghr corghr;		
   float  diff;
   void *hModule, *dModule;

   hModule= NULL; dModule = NULL;
   //ilo = 10;
   //ihi = 10;
  
   orghr_complex_parameters corghr_obj(LAPACK_COL_MAJOR, 128, 10, 10, 128);
   
   dModule = dlopen("Netlib_lapack_ref_lib/libblas.so.3.9.0", RTLD_NOW | RTLD_GLOBAL);
   
   hModule = dlopen("Netlib_lapack_ref_lib/liblapacke.so.3.9.0", RTLD_NOW);
   
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   
   corghr = (Fptr_NL_LAPACKE_corghr)dlsym(hModule, "LAPACKE_corghr");
   if (NULL == corghr)
   {
   	  printf("Could not get the symbol. Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (0);
   }
      /* Compute libflame's Lapacke o/p  */
      corghr_obj.info    = LAPACKE_corghr( LAPACK_COL_MAJOR, corghr_obj.n, corghr_obj.ilo,
	                                  corghr_obj.ihi, corghr_obj.A, corghr_obj.lda, corghr_obj.tau);
	  
      /* Compute the reference o/p by invoking Netlib-Lapack's API */	
      corghr_obj.inforef = corghr( LAPACK_COL_MAJOR, corghr_obj.n, corghr_obj.ilo, corghr_obj.ihi,
		     			 corghr_obj.Aref, corghr_obj.lda, corghr_obj.tauref);

   	/* Check for the exact singularity */
	if( corghr_obj.info > 0 ) {
		printf( "The diagonal element of the triangular factor of A,\n" );
		printf( "U(%i,%i) is zero, so that A is singular;\n", corghr_obj.info, corghr_obj.info );
		printf( "the solution could not be computed. Exiting...  \n" );

		exit( 1 );
	}
	if( corghr_obj.inforef > 0 ) {
		printf( "The diagonal element of the triangular factor of Aref,\n" );
		printf( "U(%i,%i) is zero, so that Aref is singular;\n", corghr_obj.inforef, corghr_obj.inforef );
		printf( "the solution could not be computed. Exiting...  \n" );

		exit( 1 );
	}

	/* Compute Difference in C and CPP buffer */
	diff =  computeDiff_c( corghr_obj.n*corghr_obj.lda,  corghr_obj.A, corghr_obj.Aref );

	/*
	FLA_Nrm2_external( dgetrf_obj.A, dgetrf_obj.norm );
        FLA_Nrm2_external( dgetrf_obj.Aref, dgetrf_obj. normref );
	diff = FLA_abs( dgetrf_obj. norm,  dgetrf_obj. normref);
	*/
    
    	//EXPECT_EQ (0, diff);
	//EXPECT_FLOAT_EQ (0.0, diff);
	EXPECT_NEAR(0.0, diff, 1.0);
}
/* Test for double Complex */

TEST(orghr, zorghr1) {

	/* LAPACKE zorghr prototype */
	typedef int (*Fptr_NL_LAPACKE_zorghr) (int matrix_layout, lapack_int n, lapack_int ilo,
		lapack_int ihi, lapack_complex_double* A, lapack_int lda, lapack_complex_double* tau);

	Fptr_NL_LAPACKE_zorghr zorghr;
	double  diff;
	void* hModule, * dModule;

	hModule = NULL; dModule = NULL;
//	ilo = 10;
//	ihi = 10;

	orghr_complex_double_parameters zorghr_obj(LAPACK_COL_MAJOR, 128, 10, 10, 128);

	dModule = dlopen("Netlib_lapack_ref_lib/libblas.so.3.9.0", RTLD_NOW | RTLD_GLOBAL);

	hModule = dlopen("Netlib_lapack_ref_lib/liblapacke.so.3.9.0", RTLD_NOW);

	if ((NULL == hModule) || (NULL == dModule))
	{
		printf("Load Library failed. Exiting ....\n");
		exit(0);
	}

	zorghr = (Fptr_NL_LAPACKE_zorghr)dlsym(hModule, "LAPACKE_zorghr");
	if (NULL == zorghr)
	{
		printf("Could not get the symbol. Exiting...\n");
		dlclose(hModule);
		dlclose(dModule);
		exit(0);
	}
	/* Compute libflame's Lapacke o/p  */
	zorghr_obj.info = LAPACKE_zorghr(LAPACK_COL_MAJOR, zorghr_obj.n, zorghr_obj.ilo,
		zorghr_obj.ihi, zorghr_obj.A, zorghr_obj.lda, zorghr_obj.tau);

	/* Compute the reference o/p by invoking Netlib-Lapack's API */
	zorghr_obj.inforef = zorghr(LAPACK_COL_MAJOR, zorghr_obj.n, zorghr_obj.ilo, zorghr_obj.ihi,
		zorghr_obj.Aref, zorghr_obj.lda, zorghr_obj.tauref);

	/* Check for the exact singularity */
	if (zorghr_obj.info > 0) {
		printf("The diagonal element of the triangular factor of A,\n");
		printf("U(%i,%i) is zero, so that A is singular;\n", zorghr_obj.info, zorghr_obj.info);
		printf("the solution could not be computed. Exiting...  \n");

		exit(1);
	}
	if (zorghr_obj.inforef > 0) {
		printf("The diagonal element of the triangular factor of Aref,\n");
		printf("U(%i,%i) is zero, so that Aref is singular;\n", zorghr_obj.inforef, zorghr_obj.inforef);
		printf("the solution could not be computed. Exiting...  \n");

		exit(1);
	}

	/* Compute Difference in C and CPP buffer */
	diff = computeDiff_z(zorghr_obj.n * zorghr_obj.lda, zorghr_obj.A, zorghr_obj.Aref);

	/*
	FLA_Nrm2_external( dgetrf_obj.A, dgetrf_obj.norm );
		FLA_Nrm2_external( dgetrf_obj.Aref, dgetrf_obj. normref );
	diff = FLA_abs( dgetrf_obj. norm,  dgetrf_obj. normref);
	*/
	
	//EXPECT_EQ (0, diff);
//EXPECT_FLOAT_EQ (0.0, diff);
	EXPECT_NEAR(0.0, diff, 1.0);
}
#endif



