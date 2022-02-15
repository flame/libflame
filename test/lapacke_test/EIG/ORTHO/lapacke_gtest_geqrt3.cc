#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"

/* Begin float_common_parameters  class definition */
class geqrt3_float_parameters{

   public:
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int ldt;
	float* A;	
	lapack_int lda;
	/*Output Parameter*/
	float* t;
	float *Aref, *tref;
	/*Return Values*/
	int info, inforef;

   public:
      geqrt3_float_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda, lapack_int ldt);
      ~geqrt3_float_parameters ();

};


/* Begin double_common_parameters  class definition */
class geqrt3_double_parameters {

public:
	/*Input Parameter*/
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int ldt;
	double* A;
	lapack_int lda;
	/*output Parameter*/
	double* t;
	double* Aref, * tref;
	/*Return Values*/
	int info, inforef;

public:
	geqrt3_double_parameters(int matrix_layout, lapack_int m, lapack_int n, lapack_int lda, lapack_int ldt);
	~geqrt3_double_parameters();

};


/* Begin geqrt3_complex  class definition */
class geqrt3_complex_parameters {
public:
	/*Input Parameter*/
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int ldt;
	lapack_complex_float * A;
	lapack_int lda;
	/*Output Parameter*/
	lapack_complex_float * t;
	lapack_complex_float * Aref, * tref;
	/*Return Values*/
	int info, inforef;

public:
	geqrt3_complex_parameters(int matrix_layout, lapack_int m, lapack_int n, lapack_int lda, lapack_int ldt);
	~geqrt3_complex_parameters();

};
/* Begin geqrt3_complex_double_parameters  class definition */
class geqrt3_complex_double_parameters {

public:
	/*Input Parameter*/
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int ldt;
	lapack_complex_double* A;
	lapack_int lda;
	/*Output Parameter*/
	lapack_complex_double* t;
	lapack_complex_double* Aref, * tref;
	/*Return Values*/
	int info, inforef;

public:
	geqrt3_complex_double_parameters(int matrix_layout, lapack_int m, lapack_int n, lapack_int lda, lapack_int ldt);
	~geqrt3_complex_double_parameters();

};


/* Constructor definition  float_common_parameters */
geqrt3_float_parameters:: geqrt3_float_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i, lapack_int ldt_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	ldt = ldt_i;

	A = (float *)malloc(lda*m*sizeof(float)) ;
	t = (float *)malloc(ldt*(min(m,n))*sizeof(float)) ;
	if ((A==NULL) || (t==NULL)){
		printf("error of memory allocation. Exiting ...\n");
		free(A); free(t);
		exit(0);
	}

	/* Allocation of memory for capturing reference o/ps */
	Aref = (float *)malloc(lda*n*sizeof(float)) ;
	tref = (float *)malloc(ldt*(min(m,n))*sizeof(float)) ;
	if ((Aref==NULL) || (t==NULL)){
		printf("error of memory allocation. Exiting ...\n");
		free(Aref); free(tref);
		exit(0);
	}
	/* Initialization of input matrices */
	for( i = 0; i <(lda*n); i++ ) {
		A[i] = ((float) rand()) / ((float) RAND_MAX) - 0.5;
		Aref[i] = A[i];
	}

	for(i=0;i<(ldt*(min(m,n)));i++) {
		t[i] = ((float) rand()) / ((float) RAND_MAX) - 0.5;
		tref[i] = t[i];
	}

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
geqrt3_float_parameters :: ~geqrt3_float_parameters ()
{
   /* De-Allocate memory for the input matrices */
   if (A!=NULL){
      free(A);
   }
   if (Aref!=NULL){
      free(Aref);
   }
   if (t!=NULL){
      free(t);
   }
   if (tref!=NULL){
      free(tref);
   }

}


/* Constructor definition  geqrt3_double_parameters */
geqrt3_double_parameters::geqrt3_double_parameters(int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i, lapack_int ldt_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;	
	ldt = ldt_i;

	A = (double*)malloc(lda * n * sizeof(double));
	t = (double*)malloc(ldt*(min(m,n))* sizeof(double));
	if ((A == NULL) || (t == NULL)) {
		printf("error of memory allocation. Exiting ...\n");
		free(A); free(t);
		exit(0);
	}

	/* Allocation of memory for capturing reference o/ps */
	Aref = (double*)malloc(lda * n * sizeof(double));
	tref = (double*)malloc(ldt*(min(m,n))* sizeof(double));
	if ((Aref == NULL) || (tref == NULL)) {
		printf("error of memory allocation. Exiting ...\n");
		free(Aref); free(tref);
		exit(0);
	}
	/* Initialization of input matrices */
	for (i = 0; i < (lda * n); i++) {
		A[i] = ((double)rand()) / ((double)RAND_MAX) - 0.5;
		Aref[i] = A[i];
	}

	for (i = 0; i < ldt*(min(m,n)); i++) {
		t[i] = ((double)rand()) / ((double)RAND_MAX) - 0.5;
		tref[i] = t[i];
	}

} /* end of Constructor  */

/* Destructor definition  'gethrd_double_parameters' class  */
geqrt3_double_parameters :: ~geqrt3_double_parameters()
{
	/* De-Allocate memory for the input matrices */
	if (A != NULL) {
		free(A);
	}
	if (Aref != NULL) {
		free(Aref);
	}
	if (t != NULL) {
		free(t);
	}
	if (tref != NULL) {
		free(tref);
	}

}


/* Constructor definition  geqrt3_complex_parameters */
geqrt3_complex_parameters::geqrt3_complex_parameters(int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i, lapack_int ldt_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	ldt = ldt_i;

	A = (lapack_complex_float*)malloc(lda * n * sizeof(lapack_complex_float));
	t = (lapack_complex_float*)malloc(ldt*(min(m,n)) * sizeof(lapack_complex_float));
	if ((A == NULL) || (t == NULL)) {
		printf("error of memory allocation. Exiting ...\n");
		free(A); free(t);
		exit(0);
	}

	/* Allocation of memory for capturing reference o/ps */
	Aref = (lapack_complex_float*)malloc(lda * n * sizeof(lapack_complex_float));
	tref = (lapack_complex_float*)malloc(ldt*(min(m,n)) * sizeof(lapack_complex_float));
	if ((Aref == NULL) || (tref == NULL)) {
		printf("error of memory allocation. Exiting ...\n");
		free(Aref); free(tref);
		exit(0);
	}
	/* Initialization of input matrices */
	for (i = 0; i < (lda * n); i++) {
		A[i] = ((lapack_complex_float)rand()) / ((lapack_complex_float)RAND_MAX) - 0.5;
		Aref[i] = A[i];
	}

	for (i = 0; i < ldt*(min(m,n)); i++) {
		t[i] = ((lapack_complex_float)rand()) / ((lapack_complex_float)RAND_MAX) - 0.5;
		tref[i] = t[i];
	}

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
geqrt3_complex_parameters :: ~geqrt3_complex_parameters()
{
	/* De-Allocate memory for the input matrices */
	if (A != NULL) {
		free(A);
	}
	if (Aref != NULL) {
		free(Aref);
	}
	if (t != NULL) {
		free(t);
	}
	if (tref != NULL) {
		free(tref);
	}

}

/* Constructor definition  float_common_parameters */
geqrt3_complex_double_parameters::geqrt3_complex_double_parameters(int matrix_layout_i, lapack_int m_i, lapack_int n_i,lapack_int lda_i, lapack_int ldt_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	ldt = ldt_i;

	A = (lapack_complex_double*)malloc(lda * n * sizeof(lapack_complex_double));
	t = (lapack_complex_double*)malloc(ldt*(min(m,n)) * sizeof(lapack_complex_double));
	if ((A == NULL) || (t == NULL)) {
		printf("error of memory allocation. Exiting ...\n");
		free(A); free(t);
		exit(0);
	}

	/* Allocation of memory for capturing reference o/ps */
	Aref = (lapack_complex_double*)malloc(lda * n * sizeof(lapack_complex_double));
	tref = (lapack_complex_double*)malloc(ldt*(min(m,n)) * sizeof(lapack_complex_double));
	if ((Aref == NULL) || (tref == NULL)) {
		printf("error of memory allocation. Exiting ...\n");
		free(Aref); free(tref);
		exit(0);
	}
	/* Initialization of input matrices */
	for (i = 0; i < (lda * n); i++) {
		A[i] = ((lapack_complex_double)rand()) / ((lapack_complex_double)RAND_MAX) - 0.5;
		Aref[i] = A[i];
	}

	for (i = 0; i < ldt*(min(m,n)); i++) {
		t[i] = ((lapack_complex_double)rand()) / ((lapack_complex_double)RAND_MAX) - 0.5;
		tref[i] = t[i];
	}

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
geqrt3_complex_double_parameters :: ~geqrt3_complex_double_parameters()
{
	/* De-Allocate memory for the input matrices */
	if (A != NULL) {
		free(A);
	}
	if (Aref != NULL) {
		free(Aref);
	}
	if (t != NULL) {
		free(t);
	}
	if (tref != NULL) {
		free(tref);
	}

}


/* TESTS --------------------------------------------------------------------------------------*/
TEST(geqrt3, sgeqrt3) {

   /* LAPACKE DGETRF prototype */
   typedef int (*Fptr_NL_LAPACKE_sgeqrt3) ( int matrix_layout, lapack_int m, lapack_int n, 
                           float *A, lapack_int lda, float* t, lapack_int ldt );
						   
   Fptr_NL_LAPACKE_sgeqrt3 sgeqrt3;	
   float  diff;
   void *hModule, *dModule;

   hModule= NULL; dModule = NULL;   
  
   geqrt3_float_parameters sgeqrt3_obj(LAPACK_COL_MAJOR, 512, 256, 512, 256 );
   
   dModule = dlopen("Netlib_lapack_ref_lib/libblas.so.3.9.0", RTLD_NOW | RTLD_GLOBAL);
   
   hModule = dlopen("Netlib_lapack_ref_lib/liblapacke.so.3.9.0", RTLD_NOW);
   
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   
   sgeqrt3 = (Fptr_NL_LAPACKE_sgeqrt3)dlsym(hModule, "LAPACKE_sgeqrt3");
   if (NULL == sgeqrt3)
   {
   	  printf("Could not get the symbol. Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (0);
   }
      /* Compute libflame's Lapacke o/p  */
      sgeqrt3_obj.info    = LAPACKE_sgeqrt3( LAPACK_COL_MAJOR, sgeqrt3_obj.m,  sgeqrt3_obj.n,
	                                  sgeqrt3_obj.A, sgeqrt3_obj.lda, sgeqrt3_obj.t, sgeqrt3_obj.ldt);
	  
      /* Compute the reference o/p by invoking Netlib-Lapack's API */	
      sgeqrt3_obj.inforef = sgeqrt3( LAPACK_COL_MAJOR, sgeqrt3_obj.m, sgeqrt3_obj.n,
		     			 sgeqrt3_obj.Aref, sgeqrt3_obj.lda, sgeqrt3_obj.tref, sgeqrt3_obj.ldt);

   	/* Check for the exact singularity */
	if( sgeqrt3_obj.info > 0 ) {
		printf( "The diagonal element of the triangular factor of A,\n" );
		printf( "U(%i,%i) is zero, so that A is singular;\n", sgeqrt3_obj.info, sgeqrt3_obj.info );
		printf( "the solution could not be computed. Exiting...  \n" );

		exit( 1 );
	}
	if( sgeqrt3_obj.inforef > 0 ) {
		printf( "The diagonal element of the triangular factor of Aref,\n" );
		printf( "U(%i,%i) is zero, so that Aref is singular;\n", sgeqrt3_obj.inforef, sgeqrt3_obj.inforef );
		printf( "the solution could not be computed. Exiting...  \n" );

		exit( 1 );
	}

	/* Compute Difference in C and CPP buffer */
	diff =  computeDiff_s( sgeqrt3_obj.m*sgeqrt3_obj.n,  sgeqrt3_obj.A, sgeqrt3_obj.Aref );

	/*
	FLA_Nrm2_external( dgetrf_obj.A, dgetrf_obj.norm );
        FLA_Nrm2_external( dgetrf_obj.Aref, dgetrf_obj. normref );
	diff = FLA_abs( dgetrf_obj. norm,  dgetrf_obj. normref);
	*/
    
    	//EXPECT_EQ (0, diff);
	//EXPECT_FLOAT_EQ (0.0, diff);
	EXPECT_NEAR(0.0, diff, 1.0);
}


TEST(geqrt3, dgeqrt3) {

	/* LAPACKE Dgeqrt3 prototype */
	typedef int (*Fptr_NL_LAPACKE_dgeqrt3) (int matrix_layout, lapack_int m, lapack_int n,
											double* A, lapack_int lda, double* t, lapack_int ldt);

	Fptr_NL_LAPACKE_dgeqrt3 dgeqrt3;
	double  diff;
	void* hModule, * dModule;

	hModule = NULL; dModule = NULL;	

	geqrt3_double_parameters dgeqrt3_obj(LAPACK_COL_MAJOR,  512, 256, 512, 256);

	dModule = dlopen("Netlib_lapack_ref_lib/libblas.so.3.9.0", RTLD_NOW | RTLD_GLOBAL);

	hModule = dlopen("Netlib_lapack_ref_lib/liblapacke.so.3.9.0", RTLD_NOW);

	if ((NULL == hModule) || (NULL == dModule))
	{
		printf("Load Library failed. Exiting ....\n");
		exit(0);
	}

	dgeqrt3 = (Fptr_NL_LAPACKE_dgeqrt3)dlsym(hModule, "LAPACKE_dgeqrt3");
	if (NULL == dgeqrt3)
	{
		printf("Could not get the symbol. Exiting...\n");
		dlclose(hModule);
		dlclose(dModule);
		exit(0);
	}
	/* Compute libflame's Lapacke o/p  */
	dgeqrt3_obj.info = LAPACKE_dgeqrt3(LAPACK_COL_MAJOR, dgeqrt3_obj.m, dgeqrt3_obj.n,
									dgeqrt3_obj.A, dgeqrt3_obj.lda, dgeqrt3_obj.t, dgeqrt3_obj.ldt);

	/* Compute the reference o/p by invoking Netlib-Lapack's API */
	dgeqrt3_obj.inforef = dgeqrt3(LAPACK_COL_MAJOR, dgeqrt3_obj.m, dgeqrt3_obj.n,
								dgeqrt3_obj.Aref, dgeqrt3_obj.lda, dgeqrt3_obj.tref, dgeqrt3_obj.ldt);

	/* Check for the exact singularity */
	if (dgeqrt3_obj.info > 0) {
		printf("The diagonal element of the triangular factor of A,\n");
		printf("U(%i,%i) is zero, so that A is singular;\n", dgeqrt3_obj.info, dgeqrt3_obj.info);
		printf("the solution could not be computed. Exiting...  \n");

		exit(1);
	}
	if (dgeqrt3_obj.inforef > 0) {
		printf("The diagonal element of the triangular factor of Aref,\n");
		printf("U(%i,%i) is zero, so that Aref is singular;\n", dgeqrt3_obj.inforef, dgeqrt3_obj.inforef);
		printf("the solution could not be computed. Exiting...  \n");

		exit(1);
	}

	/* Compute Difference in C and CPP buffer */
	diff = computeDiff_d(dgeqrt3_obj.m * dgeqrt3_obj.n, dgeqrt3_obj.A, dgeqrt3_obj.Aref);

	/*
	FLA_Nrm2_external( dgetrf_obj.A, dgetrf_obj.norm );
		FLA_Nrm2_external( dgetrf_obj.Aref, dgetrf_obj. normref );
	diff = FLA_abs( dgetrf_obj. norm,  dgetrf_obj. normref);
	*/

	//EXPECT_EQ (0, diff);
//EXPECT_FLOAT_EQ (0.0, diff);
	EXPECT_NEAR(0.0, diff, 1.0);
}
/*TEST FOR COMPLEX */


TEST(geqrt3, cgeqrt3) {

   /* LAPACKE DGETRF prototype */
   typedef int (*Fptr_NL_LAPACKE_cgeqrt3) ( int matrix_layout, lapack_int m, lapack_int n,
                           lapack_complex_float *A, lapack_int lda, lapack_complex_float* t, lapack_int ldt );
						   
   Fptr_NL_LAPACKE_cgeqrt3 cgeqrt3;	
   float  diff;
   void *hModule, *dModule;

   hModule= NULL; dModule = NULL;   
  
   geqrt3_complex_parameters cgeqrt3_obj(LAPACK_COL_MAJOR,  512, 256, 512, 256);
   
   dModule = dlopen("Netlib_lapack_ref_lib/libblas.so.3.9.0", RTLD_NOW | RTLD_GLOBAL);
   
   hModule = dlopen("Netlib_lapack_ref_lib/liblapacke.so.3.9.0", RTLD_NOW);
   
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   
   cgeqrt3 = (Fptr_NL_LAPACKE_cgeqrt3)dlsym(hModule, "LAPACKE_cgeqrt3");
   if (NULL == cgeqrt3)
   {
   	  printf("Could not get the symbol. Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (0);
   }
      /* Compute libflame's Lapacke o/p  */
      cgeqrt3_obj.info    = LAPACKE_cgeqrt3( LAPACK_COL_MAJOR, cgeqrt3_obj.m, cgeqrt3_obj.n,	
											cgeqrt3_obj.A, cgeqrt3_obj.lda, cgeqrt3_obj.t, cgeqrt3_obj.ldt);
	  
      /* Compute the reference o/p by invoking Netlib-Lapack's API */	
      cgeqrt3_obj.inforef = cgeqrt3( LAPACK_COL_MAJOR, cgeqrt3_obj.m, cgeqrt3_obj.n,
		     			 cgeqrt3_obj.Aref, cgeqrt3_obj.lda, cgeqrt3_obj.tref, cgeqrt3_obj.ldt);

   	/* Check for the exact singularity */
	if( cgeqrt3_obj.info > 0 ) {
		printf( "The diagonal element of the triangular factor of A,\n" );
		printf( "U(%i,%i) is zero, so that A is singular;\n", cgeqrt3_obj.info, cgeqrt3_obj.info );
		printf( "the solution could not be computed. Exiting...  \n" );

		exit( 1 );
	}
	if( cgeqrt3_obj.inforef > 0 ) {
		printf( "The diagonal element of the triangular factor of Aref,\n" );
		printf( "U(%i,%i) is zero, so that Aref is singular;\n", cgeqrt3_obj.inforef, cgeqrt3_obj.inforef );
		printf( "the solution could not be computed. Exiting...  \n" );

		exit( 1 );
	}

	/* Compute Difference in C and CPP buffer */
	diff =  computeDiff_c( cgeqrt3_obj.m*cgeqrt3_obj.n,  cgeqrt3_obj.A, cgeqrt3_obj.Aref );

	/*
	FLA_Nrm2_external( dgetrf_obj.A, dgetrf_obj.norm );
        FLA_Nrm2_external( dgetrf_obj.Aref, dgetrf_obj. normref );
	diff = FLA_abs( dgetrf_obj. norm,  dgetrf_obj. normref);
	*/
    
    	//EXPECT_EQ (0, diff);
	//EXPECT_FLOAT_EQ (0.0, diff);
	EXPECT_NEAR(0.0, diff, 1.0);
}


TEST(geqrt3, zgeqrt3) {

	/* LAPACKE zgeqrt3 prototype */
	typedef int (*Fptr_NL_LAPACKE_zgeqrt3) (int matrix_layout, lapack_int m, lapack_int n, 
										lapack_complex_double* A, lapack_int lda, lapack_complex_double* t, lapack_int ldt);

	Fptr_NL_LAPACKE_zgeqrt3 zgeqrt3;
	double  diff;
	void* hModule, * dModule;

	hModule = NULL; dModule = NULL;	

	geqrt3_complex_double_parameters zgeqrt3_obj(LAPACK_COL_MAJOR,  512, 256, 512, 256);

	dModule = dlopen("Netlib_lapack_ref_lib/libblas.so.3.9.0", RTLD_NOW | RTLD_GLOBAL);

	hModule = dlopen("Netlib_lapack_ref_lib/liblapacke.so.3.9.0", RTLD_NOW);

	if ((NULL == hModule) || (NULL == dModule))
	{
		printf("Load Library failed. Exiting ....\n");
		exit(0);
	}

	zgeqrt3 = (Fptr_NL_LAPACKE_zgeqrt3)dlsym(hModule, "LAPACKE_zgeqrt3");
	if (NULL == zgeqrt3)
	{
		printf("Could not get the symbol. Exiting...\n");
		dlclose(hModule);
		dlclose(dModule);
		exit(0);
	}
	/* Compute libflame's Lapacke o/p  */
	zgeqrt3_obj.info = LAPACKE_zgeqrt3(LAPACK_COL_MAJOR, zgeqrt3_obj.m, zgeqrt3_obj.n,
									zgeqrt3_obj.A, zgeqrt3_obj.lda, zgeqrt3_obj.t, zgeqrt3_obj.ldt);

	/* Compute the reference o/p by invoking Netlib-Lapack's API */
	zgeqrt3_obj.inforef = zgeqrt3(LAPACK_COL_MAJOR, zgeqrt3_obj.m, zgeqrt3_obj.n,
									zgeqrt3_obj.Aref, zgeqrt3_obj.lda, zgeqrt3_obj.tref, zgeqrt3_obj.ldt);

	/* Check for the exact singularity */
	if (zgeqrt3_obj.info > 0) {
		printf("The diagonal element of the triangular factor of A,\n");
		printf("U(%i,%i) is zero, so that A is singular;\n", zgeqrt3_obj.info, zgeqrt3_obj.info);
		printf("the solution could not be computed. Exiting...  \n");

		exit(1);
	}
	if (zgeqrt3_obj.inforef > 0) {
		printf("The diagonal element of the triangular factor of Aref,\n");
		printf("U(%i,%i) is zero, so that Aref is singular;\n", zgeqrt3_obj.inforef, zgeqrt3_obj.inforef);
		printf("the solution could not be computed. Exiting...  \n");

		exit(1);
	}

	/* Compute Difference in C and CPP buffer */
	diff = computeDiff_z(zgeqrt3_obj.m * zgeqrt3_obj.n, zgeqrt3_obj.A, zgeqrt3_obj.Aref);

	/*
	FLA_Nrm2_external( dgetrf_obj.A, dgetrf_obj.norm );
		FLA_Nrm2_external( dgetrf_obj.Aref, dgetrf_obj. normref );
	diff = FLA_abs( dgetrf_obj. norm,  dgetrf_obj. normref);
	*/

	//EXPECT_EQ (0, diff);
//EXPECT_FLOAT_EQ (0.0, diff);
	EXPECT_NEAR(0.0, diff, 1.0);
}


