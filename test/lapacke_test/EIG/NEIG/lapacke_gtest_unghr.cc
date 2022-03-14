#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define unghr_free() \
        free(A);  \
        free(Aref); \
        free(tau); \
        free(tauref);

#define N_ORDER 512
#define ILO 128
#define IHI 256

/* Begin unghr_complex  class definition */
class unghr_complex_parameters {
public:
	/*Input Parameter*/
	int matrix_layout;
	lapack_int n;
	lapack_int ilo;
	lapack_int ihi;
	lapack_complex_float * A;
	lapack_complex_float * tau;
	lapack_int lda;
	/*Output Parameter*/	
	lapack_complex_float * Aref, * tauref;
	/*Return Values*/
	int info, inforef;

public:
	unghr_complex_parameters(int matrix_layout, lapack_int n, lapack_int ilo, lapack_int ihi, lapack_int lda);
	~unghr_complex_parameters();

};
/* Begin unghr_complex_double_parameters  class definition */
class unghr_complex_double_parameters {

public:
	/*Input Parameter*/
	int matrix_layout;
	lapack_int n;
	lapack_int ilo;
	lapack_int ihi;
	lapack_complex_double* A;
	lapack_complex_double* tau;
	lapack_int lda;
	/*Output Parameter*/	
	lapack_complex_double* Aref, * tauref;
	/*Return Values*/
	int info, inforef;

public:
	unghr_complex_double_parameters(int matrix_layout, lapack_int n, lapack_int ilo, lapack_int ihi, lapack_int lda);
	~unghr_complex_double_parameters();

};

/* Constructor definition  unghr_complex_parameters */
unghr_complex_parameters::unghr_complex_parameters(int matrix_layout_i, lapack_int n_i, lapack_int ilo_i, lapack_int ihi_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	ilo = ilo_i;
	lda = lda_i;
	ihi  = ihi_i;

	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, (lda*n));
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&tau, &tauref, (n-1));
	if ((A == NULL) || (tau == NULL) || (Aref == NULL) || (tauref == NULL)) {
		printf("error of memory allocation. Exiting ...\n");
		unghr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand(A, Aref, (lda*n));
	lapacke_gtest_init_scomplex_buffer_pair_rand(tau, tauref, (n-1));

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
unghr_complex_parameters :: ~unghr_complex_parameters()
{
	/* De-Allocate memory for the input matrices */
	unghr_free();

}

/* Constructor definition  float_common_parameters */
unghr_complex_double_parameters::unghr_complex_double_parameters(int matrix_layout_i, lapack_int n_i, lapack_int ilo_i, lapack_int ihi_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	ilo = ilo_i;
	lda = lda_i;
	ihi  = ihi_i;

	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, (lda*n));
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&tau, &tauref, (n-1));
	if ((A == NULL) || (tau == NULL) || (Aref == NULL) || (tauref == NULL)) {
		printf("error of memory allocation. Exiting ...\n");
		unghr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand(A, Aref, (lda*n));
	lapacke_gtest_init_dcomplex_buffer_pair_rand(tau, tauref, (n-1));

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
unghr_complex_double_parameters :: ~unghr_complex_double_parameters()
{
	/* De-Allocate memory for the input matrices */
	unghr_free();
}


/* TESTS --------------------------------------------------------------------------------------*/

/*TEST FOR COMPLEX */


TEST(unghr, cunghr) {

   /* LAPACKE DGETRF prototype */
   typedef int (*Fptr_NL_LAPACKE_cunghr) ( int matrix_layout, lapack_int n, lapack_int ilo, lapack_int ihi,
                           lapack_complex_float *A, lapack_int lda, const lapack_complex_float* tau );
						   
   Fptr_NL_LAPACKE_cunghr cunghr;	
   float  diff;
   void *hModule, *dModule;

   hModule= NULL; dModule = NULL;   
  
   unghr_complex_parameters cunghr_obj(LAPACK_COL_MAJOR, N_ORDER, ILO, IHI, N_ORDER);
   
   dModule = dlopen("Netlib_lapack_ref_lib/libblas.so.3.9.0", RTLD_NOW | RTLD_GLOBAL);
   
   hModule = dlopen("Netlib_lapack_ref_lib/liblapacke.so.3.9.0", RTLD_NOW);
   
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   
   cunghr = (Fptr_NL_LAPACKE_cunghr)dlsym(hModule, "LAPACKE_cunghr");
   if (NULL == cunghr)
   {
   	  printf("Could not get the symbol. Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (0);
   }
      /* Compute libflame's Lapacke o/p  */
      cunghr_obj.info    = LAPACKE_cunghr( LAPACK_COL_MAJOR, cunghr_obj.n, cunghr_obj.ilo, cunghr_obj.ihi,	
											cunghr_obj.A, cunghr_obj.lda, cunghr_obj.tau);
	  
      /* Compute the reference o/p by invoking Netlib-Lapack's API */	
      cunghr_obj.inforef = cunghr( LAPACK_COL_MAJOR, cunghr_obj.n, cunghr_obj.ilo, cunghr_obj.ihi,
		     			 cunghr_obj.Aref, cunghr_obj.lda, cunghr_obj.tauref);

   	/* Check for the exact singularity */
	if( cunghr_obj.info > 0 ) {
		printf( "The diagonal element of the triangular factor of A,\n" );
		printf( "U(%i,%i) is zero, so that A is singular;\n", cunghr_obj.info, cunghr_obj.info );
		printf( "the solution could not be computed. Exiting...  \n" );

		exit( 1 );
	}
	if( cunghr_obj.inforef > 0 ) {
		printf( "The diagonal element of the triangular factor of Aref,\n" );
		printf( "U(%i,%i) is zero, so that Aref is singular;\n", cunghr_obj.inforef, cunghr_obj.inforef );
		printf( "the solution could not be computed. Exiting...  \n" );

		exit( 1 );
	}

	/* Compute Difference in C and CPP buffer */
	diff =  computeDiff_c( cunghr_obj.lda*cunghr_obj.n,  cunghr_obj.A, cunghr_obj.Aref );

	/*
	FLA_Nrm2_external( dgetrf_obj.A, dgetrf_obj.norm );
        FLA_Nrm2_external( dgetrf_obj.Aref, dgetrf_obj. normref );
	diff = FLA_abs( dgetrf_obj. norm,  dgetrf_obj. normref);
	*/
    
    	//EXPECT_EQ (0, diff);
	//EXPECT_FLOAT_EQ (0.0, diff);
	EXPECT_NEAR(0.0, diff, 1.0);
}


TEST(unghr, zunghr) {

	/* LAPACKE zunghr prototype */
	typedef int (*Fptr_NL_LAPACKE_zunghr) (int matrix_layout, lapack_int n, lapack_int ilo, lapack_int ihi,
										lapack_complex_double* A, lapack_int lda, const lapack_complex_double* tau);

	Fptr_NL_LAPACKE_zunghr zunghr;
	double  diff;
	void* hModule, * dModule;

	hModule = NULL; dModule = NULL;	

	unghr_complex_double_parameters zunghr_obj(LAPACK_COL_MAJOR, N_ORDER, ILO, IHI, N_ORDER);

	dModule = dlopen("Netlib_lapack_ref_lib/libblas.so.3.9.0", RTLD_NOW | RTLD_GLOBAL);

	hModule = dlopen("Netlib_lapack_ref_lib/liblapacke.so.3.9.0", RTLD_NOW);

	if ((NULL == hModule) || (NULL == dModule))
	{
		printf("Load Library failed. Exiting ....\n");
		exit(0);
	}

	zunghr = (Fptr_NL_LAPACKE_zunghr)dlsym(hModule, "LAPACKE_zunghr");
	if (NULL == zunghr)
	{
		printf("Could not get the symbol. Exiting...\n");
		dlclose(hModule);
		dlclose(dModule);
		exit(0);
	}
	/* Compute libflame's Lapacke o/p  */
	zunghr_obj.info = LAPACKE_zunghr(LAPACK_COL_MAJOR, zunghr_obj.n, zunghr_obj.ilo, zunghr_obj.ihi,
									zunghr_obj.A, zunghr_obj.lda, zunghr_obj.tau);

	/* Compute the reference o/p by invoking Netlib-Lapack's API */
	zunghr_obj.inforef = zunghr(LAPACK_COL_MAJOR, zunghr_obj.n, zunghr_obj.ilo, zunghr_obj.ihi,
									zunghr_obj.Aref, zunghr_obj.lda, zunghr_obj.tauref);

	/* Check for the exact singularity */
	if (zunghr_obj.info > 0) {
		printf("The diagonal element of the triangular factor of A,\n");
		printf("U(%i,%i) is zero, so that A is singular;\n", zunghr_obj.info, zunghr_obj.info);
		printf("the solution could not be computed. Exiting...  \n");

		exit(1);
	}
	if (zunghr_obj.inforef > 0) {
		printf("The diagonal element of the triangular factor of Aref,\n");
		printf("U(%i,%i) is zero, so that Aref is singular;\n", zunghr_obj.inforef, zunghr_obj.inforef);
		printf("the solution could not be computed. Exiting...  \n");

		exit(1);
	}

	/* Compute Difference in C and CPP buffer */
	diff = computeDiff_z(zunghr_obj.lda * zunghr_obj.n, zunghr_obj.A, zunghr_obj.Aref);

	/*
	FLA_Nrm2_external( dgetrf_obj.A, dgetrf_obj.norm );
		FLA_Nrm2_external( dgetrf_obj.Aref, dgetrf_obj. normref );
	diff = FLA_abs( dgetrf_obj. norm,  dgetrf_obj. normref);
	*/

	//EXPECT_EQ (0, diff);
//EXPECT_FLOAT_EQ (0.0, diff);
	EXPECT_NEAR(0.0, diff, 1.0);
}


