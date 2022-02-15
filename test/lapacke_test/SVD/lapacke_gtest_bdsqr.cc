#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"

class bdsqr_double_parameters{
	
   public:
   /* input params to the API **/
   int 	matrix_layout;
   char uplo;
   int 	n;
   int 	ncvt;
   int 	nru;
   int 	ncc;
   double *d, *dref;
   double *e, *eref;
   double *vt, *vtref;
   int 	ldvt;
   double *u, *uref;
   int 	ldu;
   double *c, *cref;
   int 	ldc;
   double *work, *workref;
   int info, inforef;

   public: 
      bdsqr_double_parameters (int matrix_layout, char uplo, int n,
                            int ncvt, int nru, int ncc, int ldvt,
                            int ldu, int ldc );
      ~bdsqr_double_parameters ();

}; /* end of bdsqr_double_parameters  class definition */

bdsqr_double_parameters::bdsqr_double_parameters (int matrix_layout_i,
   char uplo_i, int n_i, int ncvt_i, int nru_i, int ncc_i, int ldvt_i,
   int 	ldu_i, int ldc_i )
{
	int j;
	int work_size = 4*(n_i-1);
	matrix_layout = matrix_layout_i;
	uplo = uplo_i;
	n = n_i;
	ncvt = ncvt_i;
	nru = nru_i;
	ncc = ncc_i;
	ldvt = ldvt_i;
	ldu = ldu_i;
	ldc = ldc_i;
	
/* Memory allocation of the buffers */
    d    = (double *)malloc(n*sizeof(double));
	dref = (double *)malloc(n*sizeof(double));
   
    e    = (double *)malloc((n-1)*sizeof(double));
	eref = (double *)malloc((n-1)*sizeof(double));
  
	vt    = (double *)malloc((ldvt*ncvt)*sizeof(double));
	vtref = (double *)malloc((ldvt*ncvt)*sizeof(double));

    u    = (double *)malloc((ldu*n)*sizeof(double));
	uref = (double *)malloc((ldu*n)*sizeof(double));
	
    c    = (double *)malloc((ldc*ncc)*sizeof(double));
	cref = (double *)malloc((ldc*ncc)*sizeof(double));

    work    = (double *)malloc(work_size*sizeof(double));
	workref = (double *)malloc(work_size*sizeof(double));
	
	if( (d==NULL) || (dref==NULL) || (e==NULL) || (eref==NULL)|| \
	    (vt==NULL) || (vtref==NULL) || (u==NULL) || (uref==NULL)|| \
		(c==NULL) || (cref==NULL) || (work==NULL) || \
		(workref==NULL) ){

       printf(" bdsqr_double_parameters object: malloc error. Exiting...\n");
       free (d   );
       free (dref);
       free (e   );
       free (eref);
       free (vt   );
       free (vtref);
       free (u   );
       free (uref);
       free (c   );
       free (cref);
	   free (work  );
       free (workref);

	   exit (-1);
	}

    /* Initialization of input Buffers */
	  
    for( j = 0; j < n; j++ ) {
       d[j] = ((double) rand()) / ((double) RAND_MAX) - 0.5;
       dref[j] = d[j];
    }

    for( j = 0; j < n-1; j++ ) {
       e[j] = ((double) rand()) / ((double) RAND_MAX) - 0.5;
       eref[j] = e[j];
    }

    for( j = 0; j < (ldvt*ncvt); j++ ) {
       vt[j] = ((double) rand()) / ((double) RAND_MAX) - 0.5;
       vtref[j] = vt[j];
    }

    for( j = 0; j < (ldu*n); j++ ) {
       u[j] = ((double) rand()) / ((double) RAND_MAX) - 0.5;
       uref[j] = u[j];
    }

    for( j = 0; j < (ldc*ncc); j++ ) {
       c[j] = ((double) rand()) / ((double) RAND_MAX) - 0.5;
       cref[j] = c[j];
    }
	
	info = 0;
	inforef = 0;
	
}

/* Destructor definition **/
bdsqr_double_parameters:: ~bdsqr_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
       printf(" bdsqr_double_parameters object: destructor invoked. \n");
#endif
       free (d   );
       free (dref);
       free (e   );
       free (eref);
       free (vt   );
       free (vtref);
       free (u   );
       free (uref);
       free (c   );
       free (cref);
	   free (work  );
       free (workref);

}

TEST(bdsqr, dbdsqr1) {
   /* LAPACKE dbdsqr prototype */
typedef int (*Fptr_NL_LAPACKE_dbdsqr)( int matrix_layout, char uplo, int n,
                            int ncvt, int nru, int ncc,
                            double* d, double* e, double* vt, int ldvt,
                            double* u, int ldu, double* c,
                            int ldc );
						   
   void *hModule, *dModule;
   Fptr_NL_LAPACKE_dbdsqr DBDSQR;		
   double diff;
   printf(" NL 1 \n");

   bdsqr_double_parameters dbdsqr_obj(LAPACK_COL_MAJOR, 'u', 50, 20, 20, 20, 50,
                                      20, 50 );
   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);
  
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
  
   DBDSQR = (Fptr_NL_LAPACKE_dbdsqr)dlsym(hModule, "LAPACKE_dbdsqr");

   if (NULL == DBDSQR)
   {
   	  printf("Could not get the dbdsqr symbol. Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }
   dbdsqr_obj.inforef = DBDSQR( dbdsqr_obj.matrix_layout, dbdsqr_obj.uplo, dbdsqr_obj.n,
                            dbdsqr_obj.ncvt, dbdsqr_obj.nru, dbdsqr_obj.ncc,
                            dbdsqr_obj.dref, dbdsqr_obj.eref, dbdsqr_obj.vtref, 
							dbdsqr_obj.ldvt, dbdsqr_obj.uref, dbdsqr_obj.ldu, 
							dbdsqr_obj.cref, dbdsqr_obj.ldc );

   dbdsqr_obj.info = LAPACKE_dbdsqr( /*dbdsqr_obj.matrix_layout,*/LAPACK_COL_MAJOR, dbdsqr_obj.uplo, dbdsqr_obj.n,
                            dbdsqr_obj.ncvt, dbdsqr_obj.nru, dbdsqr_obj.ncc,
                            dbdsqr_obj.d, dbdsqr_obj.e, dbdsqr_obj.vt, 
							dbdsqr_obj.ldvt, dbdsqr_obj.u, dbdsqr_obj.ldu, 
							dbdsqr_obj.c, dbdsqr_obj.ldc );

  /* info check validation is not needed here  */

  	diff =  computeDiff_d( (4*(dbdsqr_obj.n-1)), dbdsqr_obj.work, dbdsqr_obj.workref );
	EXPECT_NEAR(0.0, diff, DOUBLE_DIFF_THRESHOLD);
}



/**   Float variant of the API test case   **/

class bdsqr_float_parameters{
	
   public:
   /* input params to the API **/
   int 	matrix_layout;
   char uplo;
   int 	n;
   int 	ncvt;
   int 	nru;
   int 	ncc;
   float *d, *dref;
   float *e, *eref;
   float *vt, *vtref;
   int 	ldvt;
   float *u, *uref;
   int 	ldu;
   float *c, *cref;
   int 	ldc;
   float *work, *workref;
   int info, inforef;

   public: 
      bdsqr_float_parameters (int matrix_layout, char uplo, int n,
                            int ncvt, int nru, int ncc, int ldvt,
                            int ldu, int ldc );
      ~bdsqr_float_parameters ();

}; /* end of bdsqr_float_parameters  class definition */

bdsqr_float_parameters::bdsqr_float_parameters (int matrix_layout_i,
   char uplo_i, int n_i, int ncvt_i, int nru_i, int ncc_i, int ldvt_i,
   int 	ldu_i, int ldc_i )
{
	int j;
	float diff;
	int work_size = 4*(n_i-1);
	matrix_layout = matrix_layout_i;
	uplo = uplo_i;
	n = n_i;
	ncvt = ncvt_i;
	nru = nru_i;
	ncc = ncc_i;
	ldvt = ldvt_i;
	ldu = ldu_i;
	ldc = ldc_i;
	
/* Memory allocation of the buffers */
    d    = (float *)malloc(n*sizeof(float));
	dref = (float *)malloc(n*sizeof(float));
   
    e    = (float *)malloc((n-1)*sizeof(float));
	eref = (float *)malloc((n-1)*sizeof(float));
  
	vt    = (float *)malloc((ldvt*ncvt)*sizeof(float));
	vtref = (float *)malloc((ldvt*ncvt)*sizeof(float));

    u    = (float *)malloc((ldu*n)*sizeof(float));
	uref = (float *)malloc((ldu*n)*sizeof(float));
	
    c    = (float *)malloc((ldc*ncc)*sizeof(float));
	cref = (float *)malloc((ldc*ncc)*sizeof(float));

    work    = (float *)malloc(work_size*sizeof(float));
	workref = (float *)malloc(work_size*sizeof(float));
	
	if( (d==NULL) || (dref==NULL) || (e==NULL) || (eref==NULL)|| \
	    (vt==NULL) || (vtref==NULL) || (u==NULL) || (uref==NULL)|| \
		(c==NULL) || (cref==NULL) || (work==NULL) || \
		(workref==NULL) ){

       printf(" bdsqr_float_parameters object: malloc error. Exiting...\n");
       free (d   );
       free (dref);
       free (e   );
       free (eref);
       free (vt   );
       free (vtref);
       free (u   );
       free (uref);
       free (c   );
       free (cref);
	   free (work  );
       free (workref);

	   exit (-1);
	}

    /* Initialization of input Buffers */
	  
    for( j = 0; j < n; j++ ) {
       d[j] = ((float) rand()) / ((float) RAND_MAX) - 0.5;
       dref[j] = d[j];
    }

    for( j = 0; j < n-1; j++ ) {
       e[j] = ((float) rand()) / ((float) RAND_MAX) - 0.5;
       eref[j] = e[j];
    }

    for( j = 0; j < (ldvt*ncvt); j++ ) {
       vt[j] = ((float) rand()) / ((float) RAND_MAX) - 0.5;
       vtref[j] = vt[j];
    }

    for( j = 0; j < (ldu*n); j++ ) {
       u[j] = ((float) rand()) / ((float) RAND_MAX) - 0.5;
       uref[j] = u[j];
    }

    for( j = 0; j < (ldc*ncc); j++ ) {
       c[j] = ((float) rand()) / ((float) RAND_MAX) - 0.5;
       cref[j] = c[j];
    }
	
	info = 0;
	inforef = 0;
	
}

/* Destructor definition **/
bdsqr_float_parameters:: ~bdsqr_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
       printf(" bdsqr_float_parameters object: destructor invoked. \n");
#endif
       free (d   );
       free (dref);
       free (e   );
       free (eref);
       free (vt   );
       free (vtref);
       free (u   );
       free (uref);
       free (c   );
       free (cref);
	   free (work  );
       free (workref);

}

TEST(bdsqr, sbdsqr1) {

   /* LAPACKE sbdsqr prototype */
typedef int (*Fptr_NL_LAPACKE_sbdsqr)( int matrix_layout, char uplo, int n,
                            int ncvt, int nru, int ncc,
                            float* d, float* e, float* vt, int ldvt,
                            float* u, int ldu, float* c,
                            int ldc );
						   
   void *hModule, *dModule;
   Fptr_NL_LAPACKE_sbdsqr SBDSQR;		
   float diff;

//   bdsqr_float_parameters sbdsqr_obj(LAPACK_COL_MAJOR, 'u', 50, 20, 20, 20, 50,
//                                      20, 50 );
   bdsqr_float_parameters sbdsqr_obj(LAPACK_COL_MAJOR, 'u', 30, 20, 20, 20, 30,
                                      20, 30 );

   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
   
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

 
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
  
   SBDSQR = (Fptr_NL_LAPACKE_sbdsqr)dlsym(hModule, "LAPACKE_sbdsqr");

   if (NULL == SBDSQR)
   {
   	  printf("Could not get the sbdsqr symbol. Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }
   sbdsqr_obj.inforef = SBDSQR( sbdsqr_obj.matrix_layout, sbdsqr_obj.uplo, sbdsqr_obj.n,
                            sbdsqr_obj.ncvt, sbdsqr_obj.nru, sbdsqr_obj.ncc,
                            sbdsqr_obj.dref, sbdsqr_obj.eref, sbdsqr_obj.vtref, 
							sbdsqr_obj.ldvt, sbdsqr_obj.uref, sbdsqr_obj.ldu, 
							sbdsqr_obj.cref, sbdsqr_obj.ldc );

   sbdsqr_obj.info = LAPACKE_sbdsqr( sbdsqr_obj.matrix_layout, sbdsqr_obj.uplo, sbdsqr_obj.n,
                            sbdsqr_obj.ncvt, sbdsqr_obj.nru, sbdsqr_obj.ncc,
                            sbdsqr_obj.d, sbdsqr_obj.e, sbdsqr_obj.vt, 
							sbdsqr_obj.ldvt, sbdsqr_obj.u, sbdsqr_obj.ldu, 
							sbdsqr_obj.c, sbdsqr_obj.ldc );

  /* info check validation is not needed here  */

  	diff =  computeDiff_s( (4*(sbdsqr_obj.n-1)), sbdsqr_obj.work, sbdsqr_obj.workref );
	EXPECT_NEAR(0.0, diff, FLOAT_DIFF_THRESHOLD);
}
