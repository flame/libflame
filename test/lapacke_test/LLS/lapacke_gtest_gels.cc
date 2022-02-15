#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"


/* Begin float_common_parameters  class definition */
class gels_float_parameters{

   public:
	/* INPUT PARAMETERS */
	int matrix_layout; // Matrix Layout
	char trans; // 'N' for Linear
	lapack_int m; // Number of rows
	lapack_int n; // Number of ROws
	lapack_int nrhs; // The number of right-hand sides; the number of columns in B
	float* A; // array A 
	float* B; // Array B
	lapack_int lda; // lda = m
	lapack_int ldb; // ldb = n

	/*Output parametes */
	float *Aref,*Bref;

	/*Return Values */
	int info, inforef;

   public:
      gels_float_parameters (int matrix_layout, char trans, lapack_int m, lapack_int n,
		      		lapack_int nrhs, lapack_int lda, lapack_int ldb);
      ~gels_float_parameters ();

};

/* Constructor definition  float_common_parameters */
gels_float_parameters:: gels_float_parameters (int matrix_layout_i, char trans_i, lapack_int m_i, lapack_int n_i,
	       				lapack_int nrhs_i, lapack_int lda_i, lapack_int ldb_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	m  = m_i;
	n = n_i;
	lda = lda_i;
	ldb = ldb_i;
	nrhs = nrhs_i;
	trans = trans_i;	

	A = (float *)malloc(lda*n*sizeof(float)) ;
	B = (float *)malloc(ldb*nrhs*sizeof(float)) ;
	if ((A==NULL) || (B==NULL)){
		printf("error of memory allocation. Exiting ...\n");
		free(A); free(B);
		exit(0);
	}

	/* Allocation of memory for capturing reference o/ps */
	Aref = (float *)malloc(lda*n*sizeof(float)) ;
	Bref = (float *)malloc(ldb*nrhs*sizeof(float)) ;
	if ((Aref==NULL) || (Bref==NULL)){
		printf("error of memory allocation. Exiting ...\n");
		free(Aref); free(Bref);
		exit(0);
	}
	/* Initialization of input matrices */
	for( i = 0; i <(lda*n); i++ ) {
		A[i] = ((float) rand()) / ((float) RAND_MAX) - 0.5;
		Aref[i] = A[i];
	}

	for(i=0;i<ldb*nrhs;i++) {
		B[i] = ((float) rand()) / ((float) RAND_MAX) - 0.5;
		Bref[i] = B[i];
	}

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
gels_float_parameters :: ~gels_float_parameters ()
{
   /* De-Allocate memory for the input matrices */
   if (A!=NULL){
      free(A);
   }
   if (Aref!=NULL){
      free(Aref);
   }
   if (B!=NULL){
      free(B);
   }
   if (Bref!=NULL){
      free(Bref);
   }

}

/* Begin double_common_parameters  class definition */
class gels_double_parameters{

   public:
	/* INPUT PARAMETERS */
	int matrix_layout; // Matrix Layout
	char trans; // 'N' for Linear
	lapack_int m; // Number of rows
	lapack_int n; // Number of ROws
	lapack_int nrhs; // The number of right-hand sides; the number of columns in B
	double* A; // array A 
	double* B; // Array B
	lapack_int lda; // lda = m
	lapack_int ldb; // ldb = n

	/*Output parametes */
	double *Aref,*Bref;

	/*Return Values */
	int info, inforef;

   public:
      gels_double_parameters (int matrix_layout, char trans, lapack_int m, lapack_int n,
		      		lapack_int nrhs, lapack_int lda, lapack_int ldb);
      ~gels_double_parameters ();

};

/* Constructor definition  double_common_parameters */
gels_double_parameters:: gels_double_parameters (int matrix_layout_i, char trans_i, lapack_int m_i, lapack_int n_i,
	       				lapack_int nrhs_i, lapack_int lda_i, lapack_int ldb_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	m  = m_i;
	n = n_i;
	lda = lda_i;
	ldb = ldb_i;
	nrhs = nrhs_i;
	trans = trans_i;	

	A = (double *)malloc(lda*n*sizeof(double)) ;
	B = (double *)malloc(ldb*nrhs*sizeof(double)) ;
	if ((A==NULL) || (B==NULL)){
		printf("error of memory allocation. Exiting ...\n");
		free(A); free(B);
		exit(0);
	}

	/* Allocation of memory for capturing reference o/ps */
	Aref = (double *)malloc(lda*n*sizeof(double)) ;
	Bref = (double *)malloc(ldb*nrhs*sizeof(double)) ;
	if ((Aref==NULL) || (Bref==NULL)){
		printf("error of memory allocation. Exiting ...\n");
		free(Aref); free(Bref);
		exit(0);
	}
	/* Initialization of input matrices */
	for( i = 0; i <(lda*n); i++ ) {
		A[i] = ((double) rand()) / ((double) RAND_MAX) - 0.5;
		Aref[i] = A[i];
	}

	for(i=0;i<ldb*nrhs;i++) {
		B[i] = ((double) rand()) / ((double) RAND_MAX) - 0.5;
		Bref[i] = B[i];
	}

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
gels_double_parameters :: ~gels_double_parameters ()
{
   /* De-Allocate memory for the input matrices */
   if (A!=NULL){
      free(A);
   }
   if (Aref!=NULL){
      free(Aref);
   }
   if (B!=NULL){
      free(B);
   }
   if (Bref!=NULL){
      free(Bref);
   }

}

/* Begin gels_complex  class definition */
class gels_complex_parameters{

   public:
	/* INPUT PARAMETERS */
	int matrix_layout; // Matrix Layout
	char trans; // 'N' for Linear
	lapack_int m; // Number of rows
	lapack_int n; // Number of ROws
	lapack_int nrhs; // The number of right-hand sides; the number of columns in B
	lapack_complex_float* A; // array A 
	lapack_complex_float* B; // Array B
	lapack_int lda; // lda = m
	lapack_int ldb; // ldb = n

	/*Output parametes */
	lapack_complex_float* Aref,*Bref;

	/*Return Values */
	int info, inforef;

   public:
      gels_complex_parameters (int matrix_layout, char trans, lapack_int m, lapack_int n,
		      		lapack_int nrhs, lapack_int lda, lapack_int ldb);
      ~gels_complex_parameters ();

};

/* Constructor definition  complex_common_parameters */
gels_complex_parameters:: gels_complex_parameters (int matrix_layout_i, char trans_i, lapack_int m_i, lapack_int n_i,
	       				lapack_int nrhs_i, lapack_int lda_i, lapack_int ldb_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	m  = m_i;
	n = n_i;
	lda = lda_i;
	ldb = ldb_i;
	nrhs = nrhs_i;
	trans = trans_i;	

	A = (lapack_complex_float *)malloc(lda*n*sizeof(lapack_complex_float)) ;
	B = (lapack_complex_float *)malloc(ldb*nrhs*sizeof(lapack_complex_float)) ;
	if ((A==NULL) || (B==NULL)){
		printf("error of memory allocation. Exiting ...\n");
		free(A); free(B);
		exit(0);
	}

	/* Allocation of memory for capturing reference o/ps */
	Aref = (lapack_complex_float *)malloc(lda*n*sizeof(lapack_complex_float)) ;
	Bref = (lapack_complex_float *)malloc(ldb*nrhs*sizeof(lapack_complex_float)) ;
	if ((Aref==NULL) || (Bref==NULL)){
		printf("error of memory allocation. Exiting ...\n");
		free(Aref); free(Bref);
		exit(0);
	}
	/* Initialization of input matrices */
	for( i = 0; i <(lda*n); i++ ) {
		A[i] = ((lapack_complex_float) rand()) / ((lapack_complex_float) RAND_MAX) - 0.5;
		Aref[i] = A[i];
	}

	for(i=0;i<ldb*nrhs;i++) {
		B[i] = ((lapack_complex_float) rand()) / ((lapack_complex_float) RAND_MAX) - 0.5;
		Bref[i] = B[i];
	}

} /* end of Constructor  */

/* Destructor definition  'lapack_complex_common_parameters' class  */
gels_complex_parameters :: ~gels_complex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   if (A!=NULL){
      free(A);
   }
   if (Aref!=NULL){
      free(Aref);
   }
   if (B!=NULL){
      free(B);
   }
   if (Bref!=NULL){
      free(Bref);
   }

}

/* Begin gels_complex_double_parameters  class definition */
class gels_complex_double_parameters{

   public:
	/* INPUT PARAMETERS */
	int matrix_layout; // Matrix Layout
	char trans; // 'N' for Linear
	lapack_int m; // Number of rows
	lapack_int n; // Number of ROws
	lapack_int nrhs; // The number of right-hand sides; the number of columns in B
	lapack_complex_double* A; // array A 
	lapack_complex_double* B; // Array B
	lapack_int lda; // lda = m
	lapack_int ldb; // ldb = n

	/*Output parametes */
	lapack_complex_double *Aref,*Bref;

	/*Return Values */
	int info, inforef;

   public:
      gels_complex_double_parameters (int matrix_layout, char trans, lapack_int m, lapack_int n,
		      		lapack_int nrhs, lapack_int lda, lapack_int ldb);
      ~gels_complex_double_parameters ();

};

/* Constructor definition  complex_double_common_parameters */
gels_complex_double_parameters:: gels_complex_double_parameters (int matrix_layout_i, char trans_i, lapack_int m_i, lapack_int n_i,
	       				lapack_int nrhs_i, lapack_int lda_i, lapack_int ldb_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	m  = m_i;
	n = n_i;
	lda = lda_i;
	ldb = ldb_i;
	nrhs = nrhs_i;
	trans = trans_i;	

	A = (lapack_complex_double *)malloc(lda*n*sizeof(lapack_complex_double)) ;
	B = (lapack_complex_double *)malloc(ldb*nrhs*sizeof(lapack_complex_double)) ;
	if ((A==NULL) || (B==NULL)){
		printf("error of memory allocation. Exiting ...\n");
		free(A); free(B);
		exit(0);
	}

	/* Allocation of memory for capturing reference o/ps */
	Aref = (lapack_complex_double *)malloc(lda*n*sizeof(lapack_complex_double)) ;
	Bref = (lapack_complex_double *)malloc(ldb*nrhs*sizeof(lapack_complex_double)) ;
	if ((Aref==NULL) || (Bref==NULL)){
		printf("error of memory allocation. Exiting ...\n");
		free(Aref); free(Bref);
		exit(0);
	}
	/* Initialization of input matrices */
	for( i = 0; i <(lda*n); i++ ) {
		A[i] = ((lapack_complex_double) rand()) / ((lapack_complex_double) RAND_MAX) - 0.5;
		Aref[i] = A[i];
	}

	for(i=0;i<ldb*nrhs;i++) {
		B[i] = ((lapack_complex_double) rand()) / ((lapack_complex_double) RAND_MAX) - 0.5;
		Bref[i] = B[i];
	}

} /* end of Constructor  */

/* Destructor definition  'complex_double_common_parameters' class  */
gels_complex_double_parameters :: ~gels_complex_double_parameters ()
{
   /* De-Allocate memory for the input matrices */
   if (A!=NULL){
      free(A);
   }
   if (Aref!=NULL){
      free(Aref);
   }
   if (B!=NULL){
      free(B);
   }
   if (Bref!=NULL){
      free(Bref);
   }

}


/* TESTS --------------------------------------------------------------------------------------*/
TEST(gels, sgels) {

   /* LAPACKE SGELS prototype */
   typedef int (*Fptr_NL_LAPACKE_sgels) ( int matrix_layout, char trans, lapack_int m, lapack_int n, 
                           lapack_int nrhs, float *A, lapack_int lda, float* B, lapack_int ldb );
						   
   Fptr_NL_LAPACKE_sgels sgels;		
   float  diff;
   void *hModule, *dModule;

   hModule= NULL; dModule = NULL;   
  
   gels_float_parameters sgels_obj(LAPACK_COL_MAJOR, 'N', 256, 128, 128, 256, 256);
   
   dModule = dlopen("Netlib_lapack_ref_lib/libblas.so.3.9.0", RTLD_NOW | RTLD_GLOBAL);
   
   hModule = dlopen("Netlib_lapack_ref_lib/liblapacke.so.3.9.0", RTLD_NOW);
   
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   
   sgels = (Fptr_NL_LAPACKE_sgels)dlsym(hModule, "LAPACKE_sgels");
   if (NULL == sgels)
   {
   	  printf("Could not get the symbol. Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (0);
   }
      /* Compute libflame's Lapacke o/p  */
      sgels_obj.info    = LAPACKE_sgels( LAPACK_COL_MAJOR, sgels_obj.trans, sgels_obj.m, sgels_obj.n, 
										sgels_obj.nrhs, sgels_obj.A, sgels_obj.lda, sgels_obj.B, sgels_obj.ldb);
	  
      /* Compute the reference o/p by invoking Netlib-Lapack's API */	
      sgels_obj.inforef = sgels( LAPACK_COL_MAJOR, sgels_obj.trans, sgels_obj.m, sgels_obj.n,sgels_obj.nrhs,
								sgels_obj.Aref, sgels_obj.lda, sgels_obj.Bref, sgels_obj.ldb);

   	/* Check for the exact singularity */
	if( sgels_obj.info > 0 ) {
		printf( "The diagonal element of the triangular factor of A,\n" );
		printf( "U(%i,%i) is zero, so that A is singular;\n", sgels_obj.info, sgels_obj.info );
		printf( "the solution could not be computed. Exiting...  \n" );

		exit( 1 );
	}
	if( sgels_obj.inforef > 0 ) {
		printf( "The diagonal element of the triangular factor of Aref,\n" );
		printf( "U(%i,%i) is zero, so that Aref is singular;\n", sgels_obj.inforef, sgels_obj.inforef );
		printf( "the solution could not be computed. Exiting...  \n" );

		exit( 1 );
	}

	/* Compute Difference in C and CPP buffer */
	diff =  computeDiff_s( sgels_obj.n*sgels_obj.lda,  sgels_obj.A, sgels_obj.Aref );

	/*
	FLA_Nrm2_external( dgetrf_obj.A, dgetrf_obj.norm );
        FLA_Nrm2_external( dgetrf_obj.Aref, dgetrf_obj. normref );
	diff = FLA_abs( dgetrf_obj. norm,  dgetrf_obj. normref);
	*/
    
    	//EXPECT_EQ (0, diff);
	//EXPECT_FLOAT_EQ (0.0, diff);
	EXPECT_NEAR(0.0, diff, 1.0);
}


TEST(gels, dgels) {

	/* LAPACKE Dgels prototype */
	typedef int (*Fptr_NL_LAPACKE_dgels) (int matrix_layout, char trans, lapack_int m, lapack_int n,
		lapack_int nrhs, double* A, lapack_int lda, double* B, lapack_int ldb);

	Fptr_NL_LAPACKE_dgels dgels;
	double  diff;
	void* hModule, * dModule;

	hModule = NULL; dModule = NULL;
//	ilo = 10;
	//ihi = 10;

	gels_double_parameters dgels_obj(LAPACK_COL_MAJOR, 'N', 256, 128, 128, 256, 256);

	dModule = dlopen("Netlib_lapack_ref_lib/libblas.so.3.9.0", RTLD_NOW | RTLD_GLOBAL);

	hModule = dlopen("Netlib_lapack_ref_lib/liblapacke.so.3.9.0", RTLD_NOW);

	if ((NULL == hModule) || (NULL == dModule))
	{
		printf("Load Library failed. Exiting ....\n");
		exit(0);
	}

	dgels = (Fptr_NL_LAPACKE_dgels)dlsym(hModule, "LAPACKE_dgels");
	if (NULL == dgels)
	{
		printf("Could not get the symbol. Exiting...\n");
		dlclose(hModule);
		dlclose(dModule);
		exit(0);
	}
	/* Compute libflame's Lapacke o/p  */
	dgels_obj.info = LAPACKE_dgels(LAPACK_COL_MAJOR, dgels_obj.trans, dgels_obj.m, dgels_obj.n,
		dgels_obj.nrhs, dgels_obj.A, dgels_obj.lda, dgels_obj.B, dgels_obj.ldb);

	/* Compute the reference o/p by invoking Netlib-Lapack's API */
	dgels_obj.inforef = dgels(LAPACK_COL_MAJOR, dgels_obj.trans, dgels_obj.m, dgels_obj.n, dgels_obj.nrhs,
		dgels_obj.Aref, dgels_obj.lda, dgels_obj.Bref, dgels_obj.ldb);

	/* Check for the exact singularity */
	if (dgels_obj.info > 0) {
		printf("The diagonal element of the triangular factor of A,\n");
		printf("U(%i,%i) is zero, so that A is singular;\n", dgels_obj.info, dgels_obj.info);
		printf("the solution could not be computed. Exiting...  \n");

		exit(1);
	}
	if (dgels_obj.inforef > 0) {
		printf("The diagonal element of the triangular factor of Aref,\n");
		printf("U(%i,%i) is zero, so that Aref is singular;\n", dgels_obj.inforef, dgels_obj.inforef);
		printf("the solution could not be computed. Exiting...  \n");

		exit(1);
	}

	/* Compute Difference in C and CPP buffer */
	diff = computeDiff_d(dgels_obj.n * dgels_obj.lda, dgels_obj.A, dgels_obj.Aref);

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

TEST(gels, cgels) {

   /* LAPACKE DGETRF prototype */
   typedef int (*Fptr_NL_LAPACKE_cgels) ( int matrix_layout, char trans, lapack_int m, lapack_int n,
                           lapack_int nrhs, lapack_complex_float* A, lapack_int lda, lapack_complex_float* B, lapack_int ldb );
						   
   Fptr_NL_LAPACKE_cgels cgels;		
   float  diff;
   void *hModule, *dModule;

   hModule= NULL; dModule = NULL;
   //ilo = 10;
   //ihi = 10;
  
   gels_complex_parameters cgels_obj(LAPACK_COL_MAJOR, 'N', 256, 128, 128, 256, 256);
   
   dModule = dlopen("Netlib_lapack_ref_lib/libblas.so.3.9.0", RTLD_NOW | RTLD_GLOBAL);
   
   hModule = dlopen("Netlib_lapack_ref_lib/liblapacke.so.3.9.0", RTLD_NOW);
   
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   
   cgels = (Fptr_NL_LAPACKE_cgels)dlsym(hModule, "LAPACKE_cgels");
   if (NULL == cgels)
   {
   	  printf("Could not get the symbol. Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (0);
   }
      /* Compute libflame's Lapacke o/p  */
      cgels_obj.info    = LAPACKE_cgels( LAPACK_COL_MAJOR, cgels_obj.trans, cgels_obj.m, cgels_obj.n,
	                                  cgels_obj.nrhs, cgels_obj.A, cgels_obj.lda, cgels_obj.B, cgels_obj.ldb);
	  
      /* Compute the reference o/p by invoking Netlib-Lapack's API */	
      cgels_obj.inforef = cgels( LAPACK_COL_MAJOR, cgels_obj.trans, cgels_obj.m, cgels_obj.n, cgels_obj.nrhs,
		     			 cgels_obj.Aref, cgels_obj.lda, cgels_obj.Bref, cgels_obj.ldb);

   	/* Check for the exact singularity */
	if( cgels_obj.info > 0 ) {
		printf( "The diagonal element of the triangular factor of A,\n" );
		printf( "U(%i,%i) is zero, so that A is singular;\n", cgels_obj.info, cgels_obj.info );
		printf( "the solution could not be computed. Exiting...  \n" );

		exit( 1 );
	}
	if( cgels_obj.inforef > 0 ) {
		printf( "The diagonal element of the triangular factor of Aref,\n" );
		printf( "U(%i,%i) is zero, so that Aref is singular;\n", cgels_obj.inforef, cgels_obj.inforef );
		printf( "the solution could not be computed. Exiting...  \n" );

		exit( 1 );
	}

	/* Compute Difference in C and CPP buffer */
	diff =  computeDiff_c( cgels_obj.n*cgels_obj.lda,  cgels_obj.A, cgels_obj.Aref );

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

TEST(gels, zgels) {

	/* LAPACKE zgels prototype */
	typedef int (*Fptr_NL_LAPACKE_zgels) (int matrix_layout, char  trans, lapack_int m, lapack_int n,
		lapack_int nrhs, lapack_complex_double* A, lapack_int lda, lapack_complex_double* B, lapack_int ldb);

	Fptr_NL_LAPACKE_zgels zgels;
	double  diff;
	void* hModule, * dModule;

	hModule = NULL; dModule = NULL;
//	ilo = 10;
//	ihi = 10;

	gels_complex_double_parameters zgels_obj(LAPACK_COL_MAJOR, 'N', 256, 128, 128, 256, 256);

	dModule = dlopen("Netlib_lapack_ref_lib/libblas.so.3.9.0", RTLD_NOW | RTLD_GLOBAL);

	hModule = dlopen("Netlib_lapack_ref_lib/liblapacke.so.3.9.0", RTLD_NOW);

	if ((NULL == hModule) || (NULL == dModule))
	{
		printf("Load Library failed. Exiting ....\n");
		exit(0);
	}

	zgels = (Fptr_NL_LAPACKE_zgels)dlsym(hModule, "LAPACKE_zgels");
	if (NULL == zgels)
	{
		printf("Could not get the symbol. Exiting...\n");
		dlclose(hModule);
		dlclose(dModule);
		exit(0);
	}
	/* Compute libflame's Lapacke o/p  */
	zgels_obj.info = LAPACKE_zgels(LAPACK_COL_MAJOR, zgels_obj.trans, zgels_obj.m, zgels_obj.n,
		zgels_obj.nrhs, zgels_obj.A, zgels_obj.lda, zgels_obj.B, zgels_obj.ldb);

	/* Compute the reference o/p by invoking Netlib-Lapack's API */
	zgels_obj.inforef = zgels(LAPACK_COL_MAJOR, zgels_obj.trans, zgels_obj.m, zgels_obj.n, zgels_obj.nrhs,
		zgels_obj.Aref, zgels_obj.lda, zgels_obj.Bref, zgels_obj.ldb);

	/* Check for the exact singularity */
	if (zgels_obj.info > 0) {
		printf("The diagonal element of the triangular factor of A,\n");
		printf("U(%i,%i) is zero, so that A is singular;\n", zgels_obj.info, zgels_obj.info);
		printf("the solution could not be computed. Exiting...  \n");

		exit(1);
	}
	if (zgels_obj.inforef > 0) {
		printf("The diagonal element of the triangular factor of Aref,\n");
		printf("U(%i,%i) is zero, so that Aref is singular;\n", zgels_obj.inforef, zgels_obj.inforef);
		printf("the solution could not be computed. Exiting...  \n");

		exit(1);
	}

	/* Compute Difference in C and CPP buffer */
	diff = computeDiff_z(zgels_obj.n * zgels_obj.lda, zgels_obj.A, zgels_obj.Aref);

	/*
	FLA_Nrm2_external( dgetrf_obj.A, dgetrf_obj.norm );
		FLA_Nrm2_external( dgetrf_obj.Aref, dgetrf_obj. normref );
	diff = FLA_abs( dgetrf_obj. norm,  dgetrf_obj. normref);
	*/
	
	//EXPECT_EQ (0, diff);
//EXPECT_FLOAT_EQ (0.0, diff);
	EXPECT_NEAR(0.0, diff, 1.0);
}





