#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"


/* Begin float_common_parameters  class definition */
class gehrd_float_parameters{

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
      gehrd_float_parameters (int matrix_layout, lapack_int n, lapack_int ilo,
		      		lapack_int ihi, lapack_int lda);
      ~gehrd_float_parameters ();

};

/* Constructor definition  float_common_parameters */
gehrd_float_parameters:: gehrd_float_parameters (int matrix_layout_i, lapack_int n_i, lapack_int ilo_i,
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
gehrd_float_parameters :: ~gehrd_float_parameters ()
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
class gehrd_double_parameters {

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
	gehrd_double_parameters(int matrix_layout, lapack_int n, lapack_int ilo,
							lapack_int ihi, lapack_int lda);
	~gehrd_double_parameters();

};

/* Constructor definition  gehrd_double_parameters */
gehrd_double_parameters::gehrd_double_parameters( int matrix_layout_i, lapack_int n_i, lapack_int ilo_i,
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
gehrd_double_parameters :: ~gehrd_double_parameters()
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

/* Begin gehrd_complex  class definition */
class gehrd_complex_parameters {
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
	gehrd_complex_parameters(int matrix_layout, lapack_int n, lapack_int ilo,
								lapack_int ihi, lapack_int lda);
	~gehrd_complex_parameters();

};

/* Constructor definition  gehrd_complex_parameters */
gehrd_complex_parameters::gehrd_complex_parameters(int matrix_layout_i, lapack_int n_i, lapack_int ilo_i,
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
gehrd_complex_parameters :: ~gehrd_complex_parameters()
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

/* Begin gehrd_complex_double_parameters  class definition */
class gehrd_complex_double_parameters {

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
	gehrd_complex_double_parameters(int matrix_layout, lapack_int n, 
									lapack_int ilo, lapack_int ihi, lapack_int lda);
	~gehrd_complex_double_parameters();

};

/* Constructor definition  float_common_parameters */
gehrd_complex_double_parameters::gehrd_complex_double_parameters(int matrix_layout_i, lapack_int n_i, lapack_int ilo_i,
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
gehrd_complex_double_parameters :: ~gehrd_complex_double_parameters()
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


/* TESTS --------------------------------------------------------------------------------------*/
TEST(gehrd, sgehrd1) {

   /* LAPACKE DGETRF prototype */
   typedef int (*Fptr_NL_LAPACKE_sgehrd) ( int matrix_layout, lapack_int n, lapack_int ilo,
                           lapack_int ihi, float *A, lapack_int lda, float* tau );
						   
   Fptr_NL_LAPACKE_sgehrd sgehrd;		
   float  diff;
   void *hModule, *dModule;

   hModule= NULL; dModule = NULL;
   //ilo = 10;
   //ihi = 10;
  
   gehrd_float_parameters sgehrd_obj(LAPACK_COL_MAJOR, 128, 10, 10, 128);
   
   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
   
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);
   
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   
   sgehrd = (Fptr_NL_LAPACKE_sgehrd)dlsym(hModule, "LAPACKE_sgehrd");
   if (NULL == sgehrd)
   {
   	  printf("Could not get the symbol. Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (0);
   }
      /* Compute libflame's Lapacke o/p  */
      sgehrd_obj.info    = LAPACKE_sgehrd( LAPACK_COL_MAJOR, sgehrd_obj.n, sgehrd_obj.ilo,
	                                  sgehrd_obj.ihi, sgehrd_obj.A, sgehrd_obj.lda, sgehrd_obj.tau);
	  
      /* Compute the reference o/p by invoking Netlib-Lapack's API */	
      sgehrd_obj.inforef = sgehrd( LAPACK_COL_MAJOR, sgehrd_obj.n, sgehrd_obj.ilo, sgehrd_obj.ihi,
		     			 sgehrd_obj.Aref, sgehrd_obj.lda, sgehrd_obj.tauref);

   	/* Check for the exact singularity */
	if( sgehrd_obj.info > 0 ) {
		printf( "The diagonal element of the triangular factor of A,\n" );
		printf( "U(%i,%i) is zero, so that A is singular;\n", sgehrd_obj.info, sgehrd_obj.info );
		printf( "the solution could not be computed. Exiting...  \n" );

		exit( 1 );
	}
	if( sgehrd_obj.inforef > 0 ) {
		printf( "The diagonal element of the triangular factor of Aref,\n" );
		printf( "U(%i,%i) is zero, so that Aref is singular;\n", sgehrd_obj.inforef, sgehrd_obj.inforef );
		printf( "the solution could not be computed. Exiting...  \n" );

		exit( 1 );
	}

	/* Compute Difference in C and CPP buffer */
	diff =  computeDiff_s( sgehrd_obj.n*sgehrd_obj.lda,  sgehrd_obj.A, sgehrd_obj.Aref );

	/*
	FLA_Nrm2_external( dgetrf_obj.A, dgetrf_obj.norm );
        FLA_Nrm2_external( dgetrf_obj.Aref, dgetrf_obj. normref );
	diff = FLA_abs( dgetrf_obj. norm,  dgetrf_obj. normref);
	*/
    
    	//EXPECT_EQ (0, diff);
	//EXPECT_FLOAT_EQ (0.0, diff);
	EXPECT_NEAR(0.0, diff, 1.0);
}


TEST(gehrd, dgehrd1) {

	/* LAPACKE Dgehrd prototype */
	typedef int (*Fptr_NL_LAPACKE_dgehrd) (int matrix_layout, lapack_int n, lapack_int ilo,
		lapack_int ihi, double* A, lapack_int lda, double* tau);

	Fptr_NL_LAPACKE_dgehrd dgehrd;
	double  diff;
	void* hModule, * dModule;

	hModule = NULL; dModule = NULL;
//	ilo = 10;
	//ihi = 10;

	gehrd_double_parameters dgehrd_obj(LAPACK_COL_MAJOR, 128, 10, 10, 128);

	dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);

	hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

	if ((NULL == hModule) || (NULL == dModule))
	{
		printf("Load Library failed. Exiting ....\n");
		exit(0);
	}

	dgehrd = (Fptr_NL_LAPACKE_dgehrd)dlsym(hModule, "LAPACKE_dgehrd");
	if (NULL == dgehrd)
	{
		printf("Could not get the symbol. Exiting...\n");
		dlclose(hModule);
		dlclose(dModule);
		exit(0);
	}
	/* Compute libflame's Lapacke o/p  */
	dgehrd_obj.info = LAPACKE_dgehrd(LAPACK_COL_MAJOR, dgehrd_obj.n, dgehrd_obj.ilo,
		dgehrd_obj.ihi, dgehrd_obj.A, dgehrd_obj.lda, dgehrd_obj.tau);

	/* Compute the reference o/p by invoking Netlib-Lapack's API */
	dgehrd_obj.inforef = dgehrd(LAPACK_COL_MAJOR, dgehrd_obj.n, dgehrd_obj.ilo, dgehrd_obj.ihi,
		dgehrd_obj.Aref, dgehrd_obj.lda, dgehrd_obj.tauref);

	/* Check for the exact singularity */
	if (dgehrd_obj.info > 0) {
		printf("The diagonal element of the triangular factor of A,\n");
		printf("U(%i,%i) is zero, so that A is singular;\n", dgehrd_obj.info, dgehrd_obj.info);
		printf("the solution could not be computed. Exiting...  \n");

		exit(1);
	}
	if (dgehrd_obj.inforef > 0) {
		printf("The diagonal element of the triangular factor of Aref,\n");
		printf("U(%i,%i) is zero, so that Aref is singular;\n", dgehrd_obj.inforef, dgehrd_obj.inforef);
		printf("the solution could not be computed. Exiting...  \n");

		exit(1);
	}

	/* Compute Difference in C and CPP buffer */
	diff = computeDiff_d(dgehrd_obj.n * dgehrd_obj.lda, dgehrd_obj.A, dgehrd_obj.Aref);

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

TEST(gehrd, cgehrd1) {

   /* LAPACKE DGETRF prototype */
   typedef int (*Fptr_NL_LAPACKE_cgehrd) ( int matrix_layout, lapack_int n, lapack_int ilo,
                           lapack_int ihi, lapack_complex_float* A, lapack_int lda, lapack_complex_float* tau );
						   
   Fptr_NL_LAPACKE_cgehrd cgehrd;		
   float  diff;
   void *hModule, *dModule;

   hModule= NULL; dModule = NULL;
   //ilo = 10;
   //ihi = 10;
  
   gehrd_complex_parameters cgehrd_obj(LAPACK_COL_MAJOR, 128, 10, 10, 128);
   
   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
   
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);
   
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   
   cgehrd = (Fptr_NL_LAPACKE_cgehrd)dlsym(hModule, "LAPACKE_cgehrd");
   if (NULL == cgehrd)
   {
   	  printf("Could not get the symbol. Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (0);
   }
      /* Compute libflame's Lapacke o/p  */
      cgehrd_obj.info    = LAPACKE_cgehrd( LAPACK_COL_MAJOR, cgehrd_obj.n, cgehrd_obj.ilo,
	                                  cgehrd_obj.ihi, cgehrd_obj.A, cgehrd_obj.lda, cgehrd_obj.tau);
	  
      /* Compute the reference o/p by invoking Netlib-Lapack's API */	
      cgehrd_obj.inforef = cgehrd( LAPACK_COL_MAJOR, cgehrd_obj.n, cgehrd_obj.ilo, cgehrd_obj.ihi,
		     			 cgehrd_obj.Aref, cgehrd_obj.lda, cgehrd_obj.tauref);

   	/* Check for the exact singularity */
	if( cgehrd_obj.info > 0 ) {
		printf( "The diagonal element of the triangular factor of A,\n" );
		printf( "U(%i,%i) is zero, so that A is singular;\n", cgehrd_obj.info, cgehrd_obj.info );
		printf( "the solution could not be computed. Exiting...  \n" );

		exit( 1 );
	}
	if( cgehrd_obj.inforef > 0 ) {
		printf( "The diagonal element of the triangular factor of Aref,\n" );
		printf( "U(%i,%i) is zero, so that Aref is singular;\n", cgehrd_obj.inforef, cgehrd_obj.inforef );
		printf( "the solution could not be computed. Exiting...  \n" );

		exit( 1 );
	}

	/* Compute Difference in C and CPP buffer */
	diff =  computeDiff_c( cgehrd_obj.n*cgehrd_obj.lda,  cgehrd_obj.A, cgehrd_obj.Aref );

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

TEST(gehrd, zgehrd1) {

	/* LAPACKE zgehrd prototype */
	typedef int (*Fptr_NL_LAPACKE_zgehrd) (int matrix_layout, lapack_int n, lapack_int ilo,
		lapack_int ihi, lapack_complex_double* A, lapack_int lda, lapack_complex_double* tau);

	Fptr_NL_LAPACKE_zgehrd zgehrd;
	double  diff;
	void* hModule, * dModule;

	hModule = NULL; dModule = NULL;
//	ilo = 10;
//	ihi = 10;

	gehrd_complex_double_parameters zgehrd_obj(LAPACK_COL_MAJOR, 128, 10, 10, 128);

	dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);

	hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

	if ((NULL == hModule) || (NULL == dModule))
	{
		printf("Load Library failed. Exiting ....\n");
		exit(0);
	}

	zgehrd = (Fptr_NL_LAPACKE_zgehrd)dlsym(hModule, "LAPACKE_zgehrd");
	if (NULL == zgehrd)
	{
		printf("Could not get the symbol. Exiting...\n");
		dlclose(hModule);
		dlclose(dModule);
		exit(0);
	}
	/* Compute libflame's Lapacke o/p  */
	zgehrd_obj.info = LAPACKE_zgehrd(LAPACK_COL_MAJOR, zgehrd_obj.n, zgehrd_obj.ilo,
		zgehrd_obj.ihi, zgehrd_obj.A, zgehrd_obj.lda, zgehrd_obj.tau);

	/* Compute the reference o/p by invoking Netlib-Lapack's API */
	zgehrd_obj.inforef = zgehrd(LAPACK_COL_MAJOR, zgehrd_obj.n, zgehrd_obj.ilo, zgehrd_obj.ihi,
		zgehrd_obj.Aref, zgehrd_obj.lda, zgehrd_obj.tauref);

	/* Check for the exact singularity */
	if (zgehrd_obj.info > 0) {
		printf("The diagonal element of the triangular factor of A,\n");
		printf("U(%i,%i) is zero, so that A is singular;\n", zgehrd_obj.info, zgehrd_obj.info);
		printf("the solution could not be computed. Exiting...  \n");

		exit(1);
	}
	if (zgehrd_obj.inforef > 0) {
		printf("The diagonal element of the triangular factor of Aref,\n");
		printf("U(%i,%i) is zero, so that Aref is singular;\n", zgehrd_obj.inforef, zgehrd_obj.inforef);
		printf("the solution could not be computed. Exiting...  \n");

		exit(1);
	}

	/* Compute Difference in C and CPP buffer */
	diff = computeDiff_z(zgehrd_obj.n * zgehrd_obj.lda, zgehrd_obj.A, zgehrd_obj.Aref);

	/*
	FLA_Nrm2_external( dgetrf_obj.A, dgetrf_obj.norm );
		FLA_Nrm2_external( dgetrf_obj.Aref, dgetrf_obj. normref );
	diff = FLA_abs( dgetrf_obj. norm,  dgetrf_obj. normref);
	*/
	
	//EXPECT_EQ (0, diff);
//EXPECT_FLOAT_EQ (0.0, diff);
	EXPECT_NEAR(0.0, diff, 1.0);
}





