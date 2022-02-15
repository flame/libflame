#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"

/* Begin float_common_parameters  class definition */
class geqrt2_float_parameters{

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
      geqrt2_float_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda, lapack_int ldt);
      ~geqrt2_float_parameters ();

};


/* Begin double_common_parameters  class definition */
class geqrt2_double_parameters {

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
	geqrt2_double_parameters(int matrix_layout, lapack_int m, lapack_int n, lapack_int lda, lapack_int ldt);
	~geqrt2_double_parameters();

};


/* Begin geqrt2_complex  class definition */
class geqrt2_complex_parameters {
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
	geqrt2_complex_parameters(int matrix_layout, lapack_int m, lapack_int n, lapack_int lda, lapack_int ldt);
	~geqrt2_complex_parameters();

};
/* Begin geqrt2_complex_double_parameters  class definition */
class geqrt2_complex_double_parameters {

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
	geqrt2_complex_double_parameters(int matrix_layout, lapack_int m, lapack_int n, lapack_int lda, lapack_int ldt);
	~geqrt2_complex_double_parameters();

};


/* Constructor definition  float_common_parameters */
geqrt2_float_parameters:: geqrt2_float_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i, lapack_int ldt_i)
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
geqrt2_float_parameters :: ~geqrt2_float_parameters ()
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


/* Constructor definition  geqrt2_double_parameters */
geqrt2_double_parameters::geqrt2_double_parameters(int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i, lapack_int ldt_i)
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
geqrt2_double_parameters :: ~geqrt2_double_parameters()
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


/* Constructor definition  geqrt2_complex_parameters */
geqrt2_complex_parameters::geqrt2_complex_parameters(int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i, lapack_int ldt_i)
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
geqrt2_complex_parameters :: ~geqrt2_complex_parameters()
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
geqrt2_complex_double_parameters::geqrt2_complex_double_parameters(int matrix_layout_i, lapack_int m_i, lapack_int n_i,lapack_int lda_i, lapack_int ldt_i)
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
geqrt2_complex_double_parameters :: ~geqrt2_complex_double_parameters()
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
TEST(geqrt2, sgeqrt2) {

   /* LAPACKE DGETRF prototype */
   typedef int (*Fptr_NL_LAPACKE_sgeqrt2) ( int matrix_layout, lapack_int m, lapack_int n, 
                           float *A, lapack_int lda, float* t, lapack_int ldt );
						   
   Fptr_NL_LAPACKE_sgeqrt2 sgeqrt2;	
   float  diff;
   void *hModule, *dModule;

   hModule= NULL; dModule = NULL;   
  
   geqrt2_float_parameters sgeqrt2_obj(LAPACK_COL_MAJOR, 512, 256, 512, 256 );
   
   dModule = dlopen("Netlib_lapack_ref_lib/libblas.so.3.9.0", RTLD_NOW | RTLD_GLOBAL);
   
   hModule = dlopen("Netlib_lapack_ref_lib/liblapacke.so.3.9.0", RTLD_NOW);
   
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   
   sgeqrt2 = (Fptr_NL_LAPACKE_sgeqrt2)dlsym(hModule, "LAPACKE_sgeqrt2");
   if (NULL == sgeqrt2)
   {
   	  printf("Could not get the symbol. Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (0);
   }
      /* Compute libflame's Lapacke o/p  */
      sgeqrt2_obj.info    = LAPACKE_sgeqrt2( LAPACK_COL_MAJOR, sgeqrt2_obj.m,  sgeqrt2_obj.n,
	                                  sgeqrt2_obj.A, sgeqrt2_obj.lda, sgeqrt2_obj.t, sgeqrt2_obj.ldt);
	  
      /* Compute the reference o/p by invoking Netlib-Lapack's API */	
      sgeqrt2_obj.inforef = sgeqrt2( LAPACK_COL_MAJOR, sgeqrt2_obj.m, sgeqrt2_obj.n,
		     			 sgeqrt2_obj.Aref, sgeqrt2_obj.lda, sgeqrt2_obj.tref, sgeqrt2_obj.ldt);

   	/* Check for the exact singularity */
	if( sgeqrt2_obj.info > 0 ) {
		printf( "The diagonal element of the triangular factor of A,\n" );
		printf( "U(%i,%i) is zero, so that A is singular;\n", sgeqrt2_obj.info, sgeqrt2_obj.info );
		printf( "the solution could not be computed. Exiting...  \n" );

		exit( 1 );
	}
	if( sgeqrt2_obj.inforef > 0 ) {
		printf( "The diagonal element of the triangular factor of Aref,\n" );
		printf( "U(%i,%i) is zero, so that Aref is singular;\n", sgeqrt2_obj.inforef, sgeqrt2_obj.inforef );
		printf( "the solution could not be computed. Exiting...  \n" );

		exit( 1 );
	}

	/* Compute Difference in C and CPP buffer */
	diff =  computeDiff_s( sgeqrt2_obj.m*sgeqrt2_obj.n,  sgeqrt2_obj.A, sgeqrt2_obj.Aref );

	/*
	FLA_Nrm2_external( dgetrf_obj.A, dgetrf_obj.norm );
        FLA_Nrm2_external( dgetrf_obj.Aref, dgetrf_obj. normref );
	diff = FLA_abs( dgetrf_obj. norm,  dgetrf_obj. normref);
	*/
    
    	//EXPECT_EQ (0, diff);
	//EXPECT_FLOAT_EQ (0.0, diff);
	EXPECT_NEAR(0.0, diff, 1.0);
}


TEST(geqrt2, dgeqrt2) {

	/* LAPACKE Dgeqrt2 prototype */
	typedef int (*Fptr_NL_LAPACKE_dgeqrt2) (int matrix_layout, lapack_int m, lapack_int n,
											double* A, lapack_int lda, double* t, lapack_int ldt);

	Fptr_NL_LAPACKE_dgeqrt2 dgeqrt2;
	double  diff;
	void* hModule, * dModule;

	hModule = NULL; dModule = NULL;	

	geqrt2_double_parameters dgeqrt2_obj(LAPACK_COL_MAJOR,  512, 256, 512, 256);

	dModule = dlopen("Netlib_lapack_ref_lib/libblas.so.3.9.0", RTLD_NOW | RTLD_GLOBAL);

	hModule = dlopen("Netlib_lapack_ref_lib/liblapacke.so.3.9.0", RTLD_NOW);

	if ((NULL == hModule) || (NULL == dModule))
	{
		printf("Load Library failed. Exiting ....\n");
		exit(0);
	}

	dgeqrt2 = (Fptr_NL_LAPACKE_dgeqrt2)dlsym(hModule, "LAPACKE_dgeqrt2");
	if (NULL == dgeqrt2)
	{
		printf("Could not get the symbol. Exiting...\n");
		dlclose(hModule);
		dlclose(dModule);
		exit(0);
	}
	/* Compute libflame's Lapacke o/p  */
	dgeqrt2_obj.info = LAPACKE_dgeqrt2(LAPACK_COL_MAJOR, dgeqrt2_obj.m, dgeqrt2_obj.n,
									dgeqrt2_obj.A, dgeqrt2_obj.lda, dgeqrt2_obj.t, dgeqrt2_obj.ldt);

	/* Compute the reference o/p by invoking Netlib-Lapack's API */
	dgeqrt2_obj.inforef = dgeqrt2(LAPACK_COL_MAJOR, dgeqrt2_obj.m, dgeqrt2_obj.n,
								dgeqrt2_obj.Aref, dgeqrt2_obj.lda, dgeqrt2_obj.tref, dgeqrt2_obj.ldt);

	/* Check for the exact singularity */
	if (dgeqrt2_obj.info > 0) {
		printf("The diagonal element of the triangular factor of A,\n");
		printf("U(%i,%i) is zero, so that A is singular;\n", dgeqrt2_obj.info, dgeqrt2_obj.info);
		printf("the solution could not be computed. Exiting...  \n");

		exit(1);
	}
	if (dgeqrt2_obj.inforef > 0) {
		printf("The diagonal element of the triangular factor of Aref,\n");
		printf("U(%i,%i) is zero, so that Aref is singular;\n", dgeqrt2_obj.inforef, dgeqrt2_obj.inforef);
		printf("the solution could not be computed. Exiting...  \n");

		exit(1);
	}

	/* Compute Difference in C and CPP buffer */
	diff = computeDiff_d(dgeqrt2_obj.m * dgeqrt2_obj.n, dgeqrt2_obj.A, dgeqrt2_obj.Aref);

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


TEST(geqrt2, cgeqrt2) {

   /* LAPACKE DGETRF prototype */
   typedef int (*Fptr_NL_LAPACKE_cgeqrt2) ( int matrix_layout, lapack_int m, lapack_int n,
                           lapack_complex_float *A, lapack_int lda, lapack_complex_float* t, lapack_int ldt );
						   
   Fptr_NL_LAPACKE_cgeqrt2 cgeqrt2;	
   float  diff;
   void *hModule, *dModule;

   hModule= NULL; dModule = NULL;   
  
   geqrt2_complex_parameters cgeqrt2_obj(LAPACK_COL_MAJOR,  512, 256, 512, 256);
   
   dModule = dlopen("Netlib_lapack_ref_lib/libblas.so.3.9.0", RTLD_NOW | RTLD_GLOBAL);
   
   hModule = dlopen("Netlib_lapack_ref_lib/liblapacke.so.3.9.0", RTLD_NOW);
   
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   
   cgeqrt2 = (Fptr_NL_LAPACKE_cgeqrt2)dlsym(hModule, "LAPACKE_cgeqrt2");
   if (NULL == cgeqrt2)
   {
   	  printf("Could not get the symbol. Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (0);
   }
      /* Compute libflame's Lapacke o/p  */
      cgeqrt2_obj.info    = LAPACKE_cgeqrt2( LAPACK_COL_MAJOR, cgeqrt2_obj.m, cgeqrt2_obj.n,	
											cgeqrt2_obj.A, cgeqrt2_obj.lda, cgeqrt2_obj.t, cgeqrt2_obj.ldt);
	  
      /* Compute the reference o/p by invoking Netlib-Lapack's API */	
      cgeqrt2_obj.inforef = cgeqrt2( LAPACK_COL_MAJOR, cgeqrt2_obj.m, cgeqrt2_obj.n,
		     			 cgeqrt2_obj.Aref, cgeqrt2_obj.lda, cgeqrt2_obj.tref, cgeqrt2_obj.ldt);

   	/* Check for the exact singularity */
	if( cgeqrt2_obj.info > 0 ) {
		printf( "The diagonal element of the triangular factor of A,\n" );
		printf( "U(%i,%i) is zero, so that A is singular;\n", cgeqrt2_obj.info, cgeqrt2_obj.info );
		printf( "the solution could not be computed. Exiting...  \n" );

		exit( 1 );
	}
	if( cgeqrt2_obj.inforef > 0 ) {
		printf( "The diagonal element of the triangular factor of Aref,\n" );
		printf( "U(%i,%i) is zero, so that Aref is singular;\n", cgeqrt2_obj.inforef, cgeqrt2_obj.inforef );
		printf( "the solution could not be computed. Exiting...  \n" );

		exit( 1 );
	}

	/* Compute Difference in C and CPP buffer */
	diff =  computeDiff_c( cgeqrt2_obj.m*cgeqrt2_obj.n,  cgeqrt2_obj.A, cgeqrt2_obj.Aref );

	/*
	FLA_Nrm2_external( dgetrf_obj.A, dgetrf_obj.norm );
        FLA_Nrm2_external( dgetrf_obj.Aref, dgetrf_obj. normref );
	diff = FLA_abs( dgetrf_obj. norm,  dgetrf_obj. normref);
	*/
    
    	//EXPECT_EQ (0, diff);
	//EXPECT_FLOAT_EQ (0.0, diff);
	EXPECT_NEAR(0.0, diff, 1.0);
}


TEST(geqrt2, zgeqrt2) {

	/* LAPACKE zgeqrt2 prototype */
	typedef int (*Fptr_NL_LAPACKE_zgeqrt2) (int matrix_layout, lapack_int m, lapack_int n, 
										lapack_complex_double* A, lapack_int lda, lapack_complex_double* t, lapack_int ldt);

	Fptr_NL_LAPACKE_zgeqrt2 zgeqrt2;
	double  diff;
	void* hModule, * dModule;

	hModule = NULL; dModule = NULL;	

	geqrt2_complex_double_parameters zgeqrt2_obj(LAPACK_COL_MAJOR,  512, 256, 512, 256);

	dModule = dlopen("Netlib_lapack_ref_lib/libblas.so.3.9.0", RTLD_NOW | RTLD_GLOBAL);

	hModule = dlopen("Netlib_lapack_ref_lib/liblapacke.so.3.9.0", RTLD_NOW);

	if ((NULL == hModule) || (NULL == dModule))
	{
		printf("Load Library failed. Exiting ....\n");
		exit(0);
	}

	zgeqrt2 = (Fptr_NL_LAPACKE_zgeqrt2)dlsym(hModule, "LAPACKE_zgeqrt2");
	if (NULL == zgeqrt2)
	{
		printf("Could not get the symbol. Exiting...\n");
		dlclose(hModule);
		dlclose(dModule);
		exit(0);
	}
	/* Compute libflame's Lapacke o/p  */
	zgeqrt2_obj.info = LAPACKE_zgeqrt2(LAPACK_COL_MAJOR, zgeqrt2_obj.m, zgeqrt2_obj.n,
									zgeqrt2_obj.A, zgeqrt2_obj.lda, zgeqrt2_obj.t, zgeqrt2_obj.ldt);

	/* Compute the reference o/p by invoking Netlib-Lapack's API */
	zgeqrt2_obj.inforef = zgeqrt2(LAPACK_COL_MAJOR, zgeqrt2_obj.m, zgeqrt2_obj.n,
									zgeqrt2_obj.Aref, zgeqrt2_obj.lda, zgeqrt2_obj.tref, zgeqrt2_obj.ldt);

	/* Check for the exact singularity */
	if (zgeqrt2_obj.info > 0) {
		printf("The diagonal element of the triangular factor of A,\n");
		printf("U(%i,%i) is zero, so that A is singular;\n", zgeqrt2_obj.info, zgeqrt2_obj.info);
		printf("the solution could not be computed. Exiting...  \n");

		exit(1);
	}
	if (zgeqrt2_obj.inforef > 0) {
		printf("The diagonal element of the triangular factor of Aref,\n");
		printf("U(%i,%i) is zero, so that Aref is singular;\n", zgeqrt2_obj.inforef, zgeqrt2_obj.inforef);
		printf("the solution could not be computed. Exiting...  \n");

		exit(1);
	}

	/* Compute Difference in C and CPP buffer */
	diff = computeDiff_z(zgeqrt2_obj.m * zgeqrt2_obj.n, zgeqrt2_obj.A, zgeqrt2_obj.Aref);

	/*
	FLA_Nrm2_external( dgetrf_obj.A, dgetrf_obj.norm );
		FLA_Nrm2_external( dgetrf_obj.Aref, dgetrf_obj. normref );
	diff = FLA_abs( dgetrf_obj. norm,  dgetrf_obj. normref);
	*/

	//EXPECT_EQ (0, diff);
//EXPECT_FLOAT_EQ (0.0, diff);
	EXPECT_NEAR(0.0, diff, 1.0);
}


