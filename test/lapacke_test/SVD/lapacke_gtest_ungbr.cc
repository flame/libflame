#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"

#define ungbr_free() \
        free(A);  \
        free(Aref); \
        free(tau); \
        free(tauref);

#define M_ROWS 512
#define N_COLS 256
#define K_REFL 128

/* Begin ungbr_complex  class definition */
class ungbr_complex_parameters {
public:
	/*Input Parameter*/
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int k;
	char vect;
	lapack_complex_float * A;
	lapack_complex_float * tau;
	lapack_int lda;
	/*Output Parameter*/	
	lapack_complex_float * Aref, * tauref;
	/*Return Values*/
	int info, inforef;

public:
	ungbr_complex_parameters(int matrix_layout, char vect, lapack_int m, lapack_int n, lapack_int k, lapack_int lda);
	~ungbr_complex_parameters();

};
/* Begin ungbr_complex_double_parameters  class definition */
class ungbr_complex_double_parameters {

public:
	/*Input Parameter*/
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int k;
	char vect;
	lapack_complex_double* A;
	lapack_complex_double* tau;
	lapack_int lda;
	/*Output Parameter*/	
	lapack_complex_double* Aref, * tauref;
	/*Return Values*/
	int info, inforef;

public:
	ungbr_complex_double_parameters(int matrix_layout, char vect, lapack_int m, lapack_int n, lapack_int k, lapack_int lda);
	~ungbr_complex_double_parameters();

};

/* Constructor definition  ungbr_complex_parameters */
ungbr_complex_parameters::ungbr_complex_parameters(int matrix_layout_i, char vect_i, lapack_int m_i, lapack_int n_i, lapack_int k_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	k  = k_i;
	vect = vect_i;

	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, (lda*n));
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&tau, &tauref, k);
	if ((A == NULL) || (tau == NULL) || (Aref == NULL) || (tauref == NULL)) {
		printf("error of memory allocation. Exiting ...\n");
		ungbr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand(A, Aref, (lda*n));
	lapacke_gtest_init_scomplex_buffer_pair_rand(tau, tauref, k);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
ungbr_complex_parameters :: ~ungbr_complex_parameters()
{
	/* De-Allocate memory for the input matrices */
	ungbr_free();

}

/* Constructor definition  float_common_parameters */
ungbr_complex_double_parameters::ungbr_complex_double_parameters(int matrix_layout_i, char vect_i, lapack_int m_i, lapack_int n_i, lapack_int k_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	k  = k_i;
	vect = vect_i;

	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, (lda*n));
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&tau, &tauref, k);
	if ((A == NULL) || (tau == NULL) || (Aref == NULL) || (tauref == NULL)) {
		printf("error of memory allocation. Exiting ...\n");
		ungbr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand(A, Aref, (lda*n));
	lapacke_gtest_init_dcomplex_buffer_pair_rand(tau, tauref, k);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
ungbr_complex_double_parameters :: ~ungbr_complex_double_parameters()
{
	/* De-Allocate memory for the input matrices */
	ungbr_free();
}


/* TESTS --------------------------------------------------------------------------------------*/

/*TEST FOR COMPLEX */


TEST(ungbr, cungbr) {

   /* LAPACKE DGETRF prototype */
   typedef int (*Fptr_NL_LAPACKE_cungbr) ( int matrix_layout, char vect, lapack_int m, lapack_int n, lapack_int k,
                           lapack_complex_float *A, lapack_int lda, const lapack_complex_float* tau );
						   
   Fptr_NL_LAPACKE_cungbr cungbr;	
   float  diff;
   void *hModule, *dModule;

   hModule= NULL; dModule = NULL;   
  
   ungbr_complex_parameters cungbr_obj(LAPACK_COL_MAJOR, 'Q', M_ROWS, N_COLS, K_REFL, M_ROWS);
   
   dModule = dlopen("Netlib_lapack_ref_lib/libblas.so.3.9.0", RTLD_NOW | RTLD_GLOBAL);
   
   hModule = dlopen("Netlib_lapack_ref_lib/liblapacke.so.3.9.0", RTLD_NOW);
   
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   
   cungbr = (Fptr_NL_LAPACKE_cungbr)dlsym(hModule, "LAPACKE_cungbr");
   if (NULL == cungbr)
   {
   	  printf("Could not get the symbol. Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (0);
   }
      /* Compute libflame's Lapacke o/p  */
      cungbr_obj.info    = LAPACKE_cungbr( LAPACK_COL_MAJOR, cungbr_obj.vect, cungbr_obj.m, cungbr_obj.n, cungbr_obj.k,	
											cungbr_obj.A, cungbr_obj.lda, cungbr_obj.tau);
	  
      /* Compute the reference o/p by invoking Netlib-Lapack's API */	
      cungbr_obj.inforef = cungbr( LAPACK_COL_MAJOR, cungbr_obj.vect, cungbr_obj.m, cungbr_obj.n, cungbr_obj.k,
		     			 cungbr_obj.Aref, cungbr_obj.lda, cungbr_obj.tauref);

   	/* Check for the exact singularity */
	if( cungbr_obj.info > 0 ) {
		printf( "The diagonal element of the triangular factor of A,\n" );
		printf( "U(%i,%i) is zero, so that A is singular;\n", cungbr_obj.info, cungbr_obj.info );
		printf( "the solution could not be computed. Exiting...  \n" );

		exit( 1 );
	}
	if( cungbr_obj.inforef > 0 ) {
		printf( "The diagonal element of the triangular factor of Aref,\n" );
		printf( "U(%i,%i) is zero, so that Aref is singular;\n", cungbr_obj.inforef, cungbr_obj.inforef );
		printf( "the solution could not be computed. Exiting...  \n" );

		exit( 1 );
	}

	/* Compute Difference in C and CPP buffer */
	diff =  computeDiff_c( cungbr_obj.m*cungbr_obj.n,  cungbr_obj.A, cungbr_obj.Aref );

	/*
	FLA_Nrm2_external( dgetrf_obj.A, dgetrf_obj.norm );
        FLA_Nrm2_external( dgetrf_obj.Aref, dgetrf_obj. normref );
	diff = FLA_abs( dgetrf_obj. norm,  dgetrf_obj. normref);
	*/
    
    	//EXPECT_EQ (0, diff);
	//EXPECT_FLOAT_EQ (0.0, diff);
	EXPECT_NEAR(0.0, diff, 1.0);
}


TEST(ungbr, zungbr) {

	/* LAPACKE zungbr prototype */
	typedef int (*Fptr_NL_LAPACKE_zungbr) (int matrix_layout, char vect, lapack_int m, lapack_int n, lapack_int k,
										lapack_complex_double* A, lapack_int lda, const lapack_complex_double* tau);

	Fptr_NL_LAPACKE_zungbr zungbr;
	double  diff;
	void* hModule, * dModule;

	hModule = NULL; dModule = NULL;	

	ungbr_complex_double_parameters zungbr_obj(LAPACK_COL_MAJOR, 'Q', M_ROWS, N_COLS, K_REFL, M_ROWS);

	dModule = dlopen("Netlib_lapack_ref_lib/libblas.so.3.9.0", RTLD_NOW | RTLD_GLOBAL);

	hModule = dlopen("Netlib_lapack_ref_lib/liblapacke.so.3.9.0", RTLD_NOW);

	if ((NULL == hModule) || (NULL == dModule))
	{
		printf("Load Library failed. Exiting ....\n");
		exit(0);
	}

	zungbr = (Fptr_NL_LAPACKE_zungbr)dlsym(hModule, "LAPACKE_zungbr");
	if (NULL == zungbr)
	{
		printf("Could not get the symbol. Exiting...\n");
		dlclose(hModule);
		dlclose(dModule);
		exit(0);
	}
	/* Compute libflame's Lapacke o/p  */
	zungbr_obj.info = LAPACKE_zungbr(LAPACK_COL_MAJOR, zungbr_obj.vect, zungbr_obj.m, zungbr_obj.n, zungbr_obj.k,
									zungbr_obj.A, zungbr_obj.lda, zungbr_obj.tau);

	/* Compute the reference o/p by invoking Netlib-Lapack's API */
	zungbr_obj.inforef = zungbr(LAPACK_COL_MAJOR, zungbr_obj.vect, zungbr_obj.m, zungbr_obj.n, zungbr_obj.k,
									zungbr_obj.Aref, zungbr_obj.lda, zungbr_obj.tauref);

	/* Check for the exact singularity */
	if (zungbr_obj.info > 0) {
		printf("The diagonal element of the triangular factor of A,\n");
		printf("U(%i,%i) is zero, so that A is singular;\n", zungbr_obj.info, zungbr_obj.info);
		printf("the solution could not be computed. Exiting...  \n");

		exit(1);
	}
	if (zungbr_obj.inforef > 0) {
		printf("The diagonal element of the triangular factor of Aref,\n");
		printf("U(%i,%i) is zero, so that Aref is singular;\n", zungbr_obj.inforef, zungbr_obj.inforef);
		printf("the solution could not be computed. Exiting...  \n");

		exit(1);
	}

	/* Compute Difference in C and CPP buffer */
	diff = computeDiff_z(zungbr_obj.m * zungbr_obj.n, zungbr_obj.A, zungbr_obj.Aref);

	/*
	FLA_Nrm2_external( dgetrf_obj.A, dgetrf_obj.norm );
		FLA_Nrm2_external( dgetrf_obj.Aref, dgetrf_obj. normref );
	diff = FLA_abs( dgetrf_obj. norm,  dgetrf_obj. normref);
	*/

	//EXPECT_EQ (0, diff);
//EXPECT_FLOAT_EQ (0.0, diff);
	EXPECT_NEAR(0.0, diff, 1.0);
}


