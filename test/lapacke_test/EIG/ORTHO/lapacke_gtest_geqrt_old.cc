#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"

#define geqrt_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (tau!=NULL)  free(tau);\
if (tauref!=NULL) free(tauref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin float_common_parameters  class definition */
class geqrt_float_parameters{

   public:
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int nb;
	lapack_int ldt;
	float* A;	
	lapack_int lda;
	/*Output Parameter*/
	float* t;
	float *Aref, *tref;
	/*Return Values*/
	int info, inforef;

   public:
      geqrt_float_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int nb, lapack_int lda, lapack_int ldt);
      ~geqrt_float_parameters ();

};


/* Begin double_common_parameters  class definition */
class geqrt_double_parameters {

public:
	/*Input Parameter*/
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int nb;
	lapack_int ldt;
	double* A;
	lapack_int lda;
	/*output Parameter*/
	double* t;
	double* Aref, * tref;
	/*Return Values*/
	int info, inforef;

public:
	geqrt_double_parameters(int matrix_layout, lapack_int m, lapack_int n, lapack_int nb, lapack_int lda, lapack_int ldt);
	~geqrt_double_parameters();

};


/* Begin geqrt_complex  class definition */
class geqrt_complex_parameters {
public:
	/*Input Parameter*/
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int nb;
	lapack_int ldt;
	lapack_complex_float * A;
	lapack_int lda;
	/*Output Parameter*/
	lapack_complex_float * t;
	lapack_complex_float * Aref, * tref;
	/*Return Values*/
	int info, inforef;

public:
	geqrt_complex_parameters(int matrix_layout, lapack_int m, lapack_int n, lapack_int nb, lapack_int lda, lapack_int ldt);
	~geqrt_complex_parameters();

};
/* Begin geqrt_complex_double_parameters  class definition */
class geqrt_complex_double_parameters {

public:
	/*Input Parameter*/
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int nb;
	lapack_int ldt;
	lapack_complex_double* A;
	lapack_int lda;
	/*Output Parameter*/
	lapack_complex_double* t;
	lapack_complex_double* Aref, * tref;
	/*Return Values*/
	int info, inforef;

public:
	geqrt_complex_double_parameters(int matrix_layout, lapack_int m, lapack_int n, lapack_int nb, lapack_int lda, lapack_int ldt);
	~geqrt_complex_double_parameters();

};


/* Constructor definition  float_common_parameters */
geqrt_float_parameters:: geqrt_float_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int nb_i, lapack_int lda_i, lapack_int ldt_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	nb = nb_i;
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
geqrt_float_parameters :: ~geqrt_float_parameters ()
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


/* Constructor definition  geqrt_double_parameters */
geqrt_double_parameters::geqrt_double_parameters(int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int nb_i, lapack_int lda_i, lapack_int ldt_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	nb = nb_i;
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
geqrt_double_parameters :: ~geqrt_double_parameters()
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


/* Constructor definition  geqrt_complex_parameters */
geqrt_complex_parameters::geqrt_complex_parameters(int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int nb_i, lapack_int lda_i, lapack_int ldt_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	ldt = ldt_i;
	nb = nb_i;

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
geqrt_complex_parameters :: ~geqrt_complex_parameters()
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
geqrt_complex_double_parameters::geqrt_complex_double_parameters(int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int nb_i,lapack_int lda_i, lapack_int ldt_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	ldt = ldt_i;
	nb = nb_i;

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
geqrt_complex_double_parameters :: ~geqrt_complex_double_parameters()
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
TEST(geqrt, sgeqrt) {

   /* LAPACKE DGETRF prototype */
   typedef int (*Fptr_NL_LAPACKE_sgeqrt) ( int matrix_layout, lapack_int m, lapack_int n, lapack_int nb, 
                           float *A, lapack_int lda, float* t, lapack_int ldt );
						   
   Fptr_NL_LAPACKE_sgeqrt sgeqrt;	
   float  diff;
   void *hModule, *dModule;

   hModule= NULL; dModule = NULL;   
  
   geqrt_float_parameters sgeqrt_obj(LAPACK_COL_MAJOR, 512, 256, 256, 512, 256 );
   
   dModule = dlopen("Netlib_lapack_ref_lib/libblas.so.3.9.0", RTLD_NOW | RTLD_GLOBAL);
   
   hModule = dlopen("Netlib_lapack_ref_lib/liblapacke.so.3.9.0", RTLD_NOW);
   
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   
   sgeqrt = (Fptr_NL_LAPACKE_sgeqrt)dlsym(hModule, "LAPACKE_sgeqrt");
   if (NULL == sgeqrt)
   {
   	  printf("Could not get the symbol. Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (0);
   }
      /* Compute libflame's Lapacke o/p  */
      sgeqrt_obj.info    = LAPACKE_sgeqrt( LAPACK_COL_MAJOR, sgeqrt_obj.m,  sgeqrt_obj.n, sgeqrt_obj.nb,
	                                  sgeqrt_obj.A, sgeqrt_obj.lda, sgeqrt_obj.t, sgeqrt_obj.ldt);
	  
      /* Compute the reference o/p by invoking Netlib-Lapack's API */	
      sgeqrt_obj.inforef = sgeqrt( LAPACK_COL_MAJOR, sgeqrt_obj.m, sgeqrt_obj.n, sgeqrt_obj.nb,
		     			 sgeqrt_obj.Aref, sgeqrt_obj.lda, sgeqrt_obj.tref, sgeqrt_obj.ldt);

   	/* Check for the exact singularity */
	if( sgeqrt_obj.info > 0 ) {
		printf( "The diagonal element of the triangular factor of A,\n" );
		printf( "U(%i,%i) is zero, so that A is singular;\n", sgeqrt_obj.info, sgeqrt_obj.info );
		printf( "the solution could not be computed. Exiting...  \n" );

		exit( 1 );
	}
	if( sgeqrt_obj.inforef > 0 ) {
		printf( "The diagonal element of the triangular factor of Aref,\n" );
		printf( "U(%i,%i) is zero, so that Aref is singular;\n", sgeqrt_obj.inforef, sgeqrt_obj.inforef );
		printf( "the solution could not be computed. Exiting...  \n" );

		exit( 1 );
	}

	/* Compute Difference in C and CPP buffer */
	diff =  computeDiff_s( sgeqrt_obj.m*sgeqrt_obj.n,  sgeqrt_obj.A, sgeqrt_obj.Aref );

	/*
	FLA_Nrm2_external( dgetrf_obj.A, dgetrf_obj.norm );
        FLA_Nrm2_external( dgetrf_obj.Aref, dgetrf_obj. normref );
	diff = FLA_abs( dgetrf_obj. norm,  dgetrf_obj. normref);
	*/
    
    	//EXPECT_EQ (0, diff);
	//EXPECT_FLOAT_EQ (0.0, diff);
	EXPECT_NEAR(0.0, diff, 1.0);
}


TEST(geqrt, dgeqrt) {

	/* LAPACKE Dgeqrt prototype */
	typedef int (*Fptr_NL_LAPACKE_dgeqrt) (int matrix_layout, lapack_int m, lapack_int n, lapack_int nb,	
											double* A, lapack_int lda, double* t, lapack_int ldt);

	Fptr_NL_LAPACKE_dgeqrt dgeqrt;
	double  diff;
	void* hModule, * dModule;

	hModule = NULL; dModule = NULL;	

	geqrt_double_parameters dgeqrt_obj(LAPACK_COL_MAJOR,  512, 256, 256, 512, 256);

	dModule = dlopen("Netlib_lapack_ref_lib/libblas.so.3.9.0", RTLD_NOW | RTLD_GLOBAL);

	hModule = dlopen("Netlib_lapack_ref_lib/liblapacke.so.3.9.0", RTLD_NOW);

	if ((NULL == hModule) || (NULL == dModule))
	{
		printf("Load Library failed. Exiting ....\n");
		exit(0);
	}

	dgeqrt = (Fptr_NL_LAPACKE_dgeqrt)dlsym(hModule, "LAPACKE_dgeqrt");
	if (NULL == dgeqrt)
	{
		printf("Could not get the symbol. Exiting...\n");
		dlclose(hModule);
		dlclose(dModule);
		exit(0);
	}
	/* Compute libflame's Lapacke o/p  */
	dgeqrt_obj.info = LAPACKE_dgeqrt(LAPACK_COL_MAJOR, dgeqrt_obj.m, dgeqrt_obj.n, dgeqrt_obj.nb,
									dgeqrt_obj.A, dgeqrt_obj.lda, dgeqrt_obj.t, dgeqrt_obj.ldt);

	/* Compute the reference o/p by invoking Netlib-Lapack's API */
	dgeqrt_obj.inforef = dgeqrt(LAPACK_COL_MAJOR, dgeqrt_obj.m, dgeqrt_obj.n, dgeqrt_obj.nb,
								dgeqrt_obj.Aref, dgeqrt_obj.lda, dgeqrt_obj.tref, dgeqrt_obj.ldt);

	/* Check for the exact singularity */
	if (dgeqrt_obj.info > 0) {
		printf("The diagonal element of the triangular factor of A,\n");
		printf("U(%i,%i) is zero, so that A is singular;\n", dgeqrt_obj.info, dgeqrt_obj.info);
		printf("the solution could not be computed. Exiting...  \n");

		exit(1);
	}
	if (dgeqrt_obj.inforef > 0) {
		printf("The diagonal element of the triangular factor of Aref,\n");
		printf("U(%i,%i) is zero, so that Aref is singular;\n", dgeqrt_obj.inforef, dgeqrt_obj.inforef);
		printf("the solution could not be computed. Exiting...  \n");

		exit(1);
	}

	/* Compute Difference in C and CPP buffer */
	diff = computeDiff_d(dgeqrt_obj.m * dgeqrt_obj.n, dgeqrt_obj.A, dgeqrt_obj.Aref);

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


TEST(geqrt, cgeqrt) {

   /* LAPACKE DGETRF prototype */
   typedef int (*Fptr_NL_LAPACKE_cgeqrt) ( int matrix_layout, lapack_int m, lapack_int n, lapack_int nb,
                           lapack_complex_float *A, lapack_int lda, lapack_complex_float* t, lapack_int ldt );
						   
   Fptr_NL_LAPACKE_cgeqrt cgeqrt;	
   float  diff;
   void *hModule, *dModule;

   hModule= NULL; dModule = NULL;   
  
   geqrt_complex_parameters cgeqrt_obj(LAPACK_COL_MAJOR,  512, 256, 256, 512, 256);
   
   dModule = dlopen("Netlib_lapack_ref_lib/libblas.so.3.9.0", RTLD_NOW | RTLD_GLOBAL);
   
   hModule = dlopen("Netlib_lapack_ref_lib/liblapacke.so.3.9.0", RTLD_NOW);
   
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   
   cgeqrt = (Fptr_NL_LAPACKE_cgeqrt)dlsym(hModule, "LAPACKE_cgeqrt");
   if (NULL == cgeqrt)
   {
   	  printf("Could not get the symbol. Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (0);
   }
      /* Compute libflame's Lapacke o/p  */
      cgeqrt_obj.info    = LAPACKE_cgeqrt( LAPACK_COL_MAJOR, cgeqrt_obj.m, cgeqrt_obj.n, cgeqrt_obj.nb,	
											cgeqrt_obj.A, cgeqrt_obj.lda, cgeqrt_obj.t, cgeqrt_obj.ldt);
	  
      /* Compute the reference o/p by invoking Netlib-Lapack's API */	
      cgeqrt_obj.inforef = cgeqrt( LAPACK_COL_MAJOR, cgeqrt_obj.m, cgeqrt_obj.n, cgeqrt_obj.nb,
		     			 cgeqrt_obj.Aref, cgeqrt_obj.lda, cgeqrt_obj.tref, cgeqrt_obj.ldt);

   	/* Check for the exact singularity */
	if( cgeqrt_obj.info > 0 ) {
		printf( "The diagonal element of the triangular factor of A,\n" );
		printf( "U(%i,%i) is zero, so that A is singular;\n", cgeqrt_obj.info, cgeqrt_obj.info );
		printf( "the solution could not be computed. Exiting...  \n" );

		exit( 1 );
	}
	if( cgeqrt_obj.inforef > 0 ) {
		printf( "The diagonal element of the triangular factor of Aref,\n" );
		printf( "U(%i,%i) is zero, so that Aref is singular;\n", cgeqrt_obj.inforef, cgeqrt_obj.inforef );
		printf( "the solution could not be computed. Exiting...  \n" );

		exit( 1 );
	}

	/* Compute Difference in C and CPP buffer */
	diff =  computeDiff_c( cgeqrt_obj.m*cgeqrt_obj.n,  cgeqrt_obj.A, cgeqrt_obj.Aref );

	/*
	FLA_Nrm2_external( dgetrf_obj.A, dgetrf_obj.norm );
        FLA_Nrm2_external( dgetrf_obj.Aref, dgetrf_obj. normref );
	diff = FLA_abs( dgetrf_obj. norm,  dgetrf_obj. normref);
	*/
    
    	//EXPECT_EQ (0, diff);
	//EXPECT_FLOAT_EQ (0.0, diff);
	EXPECT_NEAR(0.0, diff, 1.0);
}


TEST(geqrt, zgeqrt) {

	/* LAPACKE zgeqrt prototype */
	typedef int (*Fptr_NL_LAPACKE_zgeqrt) (int matrix_layout, lapack_int m, lapack_int n,lapack_int nb,
										lapack_complex_double* A, lapack_int lda, lapack_complex_double* t, lapack_int ldt);

	Fptr_NL_LAPACKE_zgeqrt zgeqrt;
	double  diff;
	void* hModule, * dModule;

	hModule = NULL; dModule = NULL;	

	geqrt_complex_double_parameters zgeqrt_obj(LAPACK_COL_MAJOR,  512, 256, 256, 512, 256);

	dModule = dlopen("Netlib_lapack_ref_lib/libblas.so.3.9.0", RTLD_NOW | RTLD_GLOBAL);

	hModule = dlopen("Netlib_lapack_ref_lib/liblapacke.so.3.9.0", RTLD_NOW);

	if ((NULL == hModule) || (NULL == dModule))
	{
		printf("Load Library failed. Exiting ....\n");
		exit(0);
	}

	zgeqrt = (Fptr_NL_LAPACKE_zgeqrt)dlsym(hModule, "LAPACKE_zgeqrt");
	if (NULL == zgeqrt)
	{
		printf("Could not get the symbol. Exiting...\n");
		dlclose(hModule);
		dlclose(dModule);
		exit(0);
	}
	/* Compute libflame's Lapacke o/p  */
	zgeqrt_obj.info = LAPACKE_zgeqrt(LAPACK_COL_MAJOR, zgeqrt_obj.m, zgeqrt_obj.n, zgeqrt_obj.nb,
									zgeqrt_obj.A, zgeqrt_obj.lda, zgeqrt_obj.t, zgeqrt_obj.ldt);

	/* Compute the reference o/p by invoking Netlib-Lapack's API */
	zgeqrt_obj.inforef = zgeqrt(LAPACK_COL_MAJOR, zgeqrt_obj.m, zgeqrt_obj.n, zgeqrt_obj.nb, 
									zgeqrt_obj.Aref, zgeqrt_obj.lda, zgeqrt_obj.tref, zgeqrt_obj.ldt);

	/* Check for the exact singularity */
	if (zgeqrt_obj.info > 0) {
		printf("The diagonal element of the triangular factor of A,\n");
		printf("U(%i,%i) is zero, so that A is singular;\n", zgeqrt_obj.info, zgeqrt_obj.info);
		printf("the solution could not be computed. Exiting...  \n");

		exit(1);
	}
	if (zgeqrt_obj.inforef > 0) {
		printf("The diagonal element of the triangular factor of Aref,\n");
		printf("U(%i,%i) is zero, so that Aref is singular;\n", zgeqrt_obj.inforef, zgeqrt_obj.inforef);
		printf("the solution could not be computed. Exiting...  \n");

		exit(1);
	}

	/* Compute Difference in C and CPP buffer */
	diff = computeDiff_z(zgeqrt_obj.m * zgeqrt_obj.n, zgeqrt_obj.A, zgeqrt_obj.Aref);

	/*
	FLA_Nrm2_external( dgetrf_obj.A, dgetrf_obj.norm );
		FLA_Nrm2_external( dgetrf_obj.Aref, dgetrf_obj. normref );
	diff = FLA_abs( dgetrf_obj. norm,  dgetrf_obj. normref);
	*/

	//EXPECT_EQ (0, diff);
//EXPECT_FLOAT_EQ (0.0, diff);
	EXPECT_NEAR(0.0, diff, 1.0);
}


