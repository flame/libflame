#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"

#define gtsvx_free() \
  if (d != NULL)    free (d  ); \
  if (dref != NULL) free (dref); \
  if (b != NULL)    free (b   ); \
  if (bref != NULL) free (bref); \
  if (x != NULL)    free (x  ); \
  if (xref != NULL) free (xref); \
  if (df != NULL)    free (df  ); \
  if (dfref != NULL) free (dfref); \
  if (dl != NULL)    free (dl  ); \
  if (dlref != NULL) free (dlref); \
  if (du != NULL)    free (du  ); \
  if (duref != NULL) free (duref); \
  if (dlf != NULL)    free (dlf  ); \
  if (dlfref != NULL) free (dlfref); \
  if (duf != NULL)    free (duf  ); \
  if (dufref != NULL) free (dufref); \
  if (du2 != NULL)    free (du2  ); \
  if (du2ref != NULL) free (du2ref); \
  if (ferr != NULL)    free (ferr  ); \
  if (ferrref != NULL) free (ferrref); \
  if (berr != NULL)    free (berr  ); \
  if (berrref != NULL) free (berrref); \
  if (ipiv != NULL)    free (ipiv  ); \
  if (ipivref != NULL) free (ipivref); \
  if( hModule != NULL) dlclose(hModule); \
  if(dModule != NULL) dlclose(dModule)


// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin gtsvx_double_parameters  class definition */
class gtsvx_double_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      // capture difference between ref o/p & libflame lapacke o/p.
      double diff_berr, diff_ferr, diff_xerr;
      double diff_dlf, diff_df, diff_duf, diff_du2;
      int diff_piv;
      void *hModule, *dModule;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char fact; //  Must be 'F', 'N', or 'E'.
      char trans; //  Must be 'N' , 'T' or 'C'.

      lapack_int n; // The order of A (Square matrix).
      lapack_int nrhs; // The number of right-hand sides
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_int ldx; //  leading dimension of 'x'
      double *b, *bref; //right-hand sides for the systems of equations.
      double *d, *dref; //size (n), contains the diagonal elements of A.
      double* dl, *dlref; // size (n -1), contains the subdiagonal elements of A.
      double* du, *duref; // size (n -1), contains the superdiagonal elements of A..

      /* Input/ Output parameters */
      double *df, *dfref; // 'n' diagonal elements of the upper triangular matrix
      double *dlf, *dlfref; // the (n -1) multipliers that define the matrix L
      double *duf, *dufref; // the (n -1) elements of the first superdiagonal of U
      double* du2, *du2ref; // the (n-2) elements of the second superdiagonal of U.
      
      lapack_int *ipiv, *ipivref; // pivot buffer
      
      /* Output parameters */
      char   equed, equedref; //  Must be 'N', 'R', 'C', or 'B'.
      double rpivot, rpivotref; // reciprocal pivot growth factor.
      double rcond, rcondref; // estimate of the reciprocal condition number of the matrix A.
      double* x, *xref;
      double* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      double* berr, *berrref; // component-wise relative backward error for each solution vector xj.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      gtsvx_double_parameters ( int matrix_layout_i, char fact_i, char trans_i,
                               char equed_i, lapack_int n_i, lapack_int nrhs_i);
              
      ~gtsvx_double_parameters (); 
};  /* end of gtsvx_double_parameters  class definition */


/* Constructor gtsvx_double_parameters definition */
gtsvx_double_parameters:: gtsvx_double_parameters ( int matrix_layout_i, 
                char fact_i, char trans_i, char equed_i, lapack_int n_i,
                  lapack_int nrhs_i) {
                                    
    matrix_layout = matrix_layout_i;
    fact= fact_i;
    trans = trans_i;
    // equed acts as input when " fact = 'F' " else output
    equed = equed_i;
    equedref = equed_i;
    n = n_i;
    nrhs = nrhs_i;
    
    b_bufsize = n*nrhs;
    x_bufsize = b_bufsize;

    if(matrix_layout==LAPACK_COL_MAJOR){
        ldb = n;
        ldx = n;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        ldb = nrhs;
        ldx = nrhs;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }
#if LAPACKE_TEST_VERBOSE
   printf(" \n gtsvx Double: matrix_layout: %d  n: %d, fact: %c trans: %c   \
ldb: %d nrhs: %d  ldaf: %d  ldx: %d \n", matrix_layout, n, fact, trans,  
                                          ldb, nrhs);
#endif

    hModule = NULL;
    dModule = NULL;
    diff_berr = 0;
    diff_ferr = 0;
    diff_xerr = 0;
    diff_dlf = 0;
    diff_df = 0;
    diff_duf = 0;
    diff_du2 = 0;
    diff_piv = 0;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_double_buffer_pair( &d, &dref, n    );
    lapacke_gtest_alloc_double_buffer_pair( &dl, &dlref, n-1  );
    lapacke_gtest_alloc_double_buffer_pair( &du, &duref, n-1  );
    lapacke_gtest_alloc_double_buffer_pair( &dlf, &dlfref, n-1  );
    lapacke_gtest_alloc_double_buffer_pair( &df, &dfref, n    );
    lapacke_gtest_alloc_double_buffer_pair( &duf, &dufref, n-1  );
    lapacke_gtest_alloc_double_buffer_pair( &du2, &du2ref, n-2  );
 
    lapacke_gtest_alloc_double_buffer_pair( &x, &xref, x_bufsize);
    lapacke_gtest_alloc_double_buffer_pair( &ferr, &ferrref, nrhs);
    lapacke_gtest_alloc_double_buffer_pair( &berr, &berrref, nrhs);
    
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (d==NULL) || (dref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (x==NULL) || (xref==NULL) || \
        (df==NULL) || (dfref==NULL) || \
        (dl==NULL) || (dlref==NULL) || \
        (du==NULL) || (duref==NULL) || \
        (dlf==NULL) || (dlfref==NULL) || \
        (duf==NULL) || (dufref==NULL) || \
        (du2==NULL) || (du2ref==NULL) || \
        (ferr==NULL) || (ferrref==NULL) || \
        (berr==NULL) || (berrref==NULL) || \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "pbtrs_double_parameters object: malloc error.";
       gtsvx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( d, dref, n);
    lapacke_gtest_init_double_buffer_pair_rand( dl, dlref, n-1);
    lapacke_gtest_init_double_buffer_pair_rand( du, duref, n-1);
    lapacke_gtest_init_double_buffer_pair_rand( dlf, dlfref, n-1);
    lapacke_gtest_init_double_buffer_pair_rand( df, dfref, n);
    lapacke_gtest_init_double_buffer_pair_rand( duf, dufref, n-1);
    lapacke_gtest_init_double_buffer_pair_rand( du2, du2ref, n-2);


    
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, b_bufsize);
    
    lapacke_gtest_init_double_buffer_pair_with_constant(x, xref, x_bufsize, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);
    
   } /* end of Constructor  */

gtsvx_double_parameters:: ~gtsvx_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gtsvx_double_parameters object: destructor invoked. \n");
#endif
   gtsvx_free();
}

//  Test fixture class definition
class dgtsvx_test  : public  ::testing::Test {
public:
   gtsvx_double_parameters  *dgtsvx_obj;
   void SetUp();  
   void TearDown () { delete dgtsvx_obj; }
};


void dgtsvx_test::SetUp(){

    /* LAPACKE DGTSVX prototype */
    typedef int (*Fptr_NL_LAPACKE_dgtsvx) (int matrix_layout, char fact,
            char trans, lapack_int n, lapack_int nrhs, const double* dl,
            const double* d, const double* du, double* dlf, double* df,
            double* duf, double* du2, lapack_int* ipiv, const double* b,
            lapack_int ldb, double* x, lapack_int ldx, double* rcond,
            double* ferr, double* berr);

    Fptr_NL_LAPACKE_dgtsvx DGTSVX;


     /* LAPACKE DGETRF prototype */
    typedef int (*Fptr_NL_LAPACKE_dgttrf) ( lapack_int n , double *dl,
            double *d , double *du , double *du2 , lapack_int *ipiv );

    Fptr_NL_LAPACKE_dgttrf DGTTRF;


    dgtsvx_obj = new gtsvx_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_driver_paramslist[idx].fact,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].equed,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    dgtsvx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgtsvx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgtsvx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgtsvx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    DGTSVX = (Fptr_NL_LAPACKE_dgtsvx)dlsym(dgtsvx_obj->hModule, "LAPACKE_dgtsvx");
    ASSERT_TRUE(DGTSVX != NULL) << "failed to get the Netlib LAPACKE_dgtsvx symbol";
    
    /* if  'fact' == 'F', then the i/p matrix A should be in 
       factored form.
       so, invoke the getrf API to compute the factorized A.  */
    if(dgtsvx_obj->fact == 'F') {
        DGTTRF = (Fptr_NL_LAPACKE_dgttrf)dlsym(dgtsvx_obj->hModule,"LAPACKE_dgttrf");
        ASSERT_TRUE(DGTTRF != NULL) << "failed to get the Netlib LAPACKE_dgttrf symbol";
            
        dgtsvx_obj->inforef = DGTTRF( dgtsvx_obj->n,
                                      dgtsvx_obj->dlref, dgtsvx_obj->dref,
                                      dgtsvx_obj->duref,
                                      dgtsvx_obj->du2ref,
                                      dgtsvx_obj->ipivref);
                               
        dgtsvx_obj->info = LAPACKE_dgttrf( dgtsvx_obj->n,
                                      dgtsvx_obj->dl, dgtsvx_obj->d,
                                      dgtsvx_obj->du,
                                      dgtsvx_obj->du2,
                                      dgtsvx_obj->ipiv);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    dgtsvx_obj->inforef = DGTSVX( dgtsvx_obj->matrix_layout, dgtsvx_obj->fact,
                                  dgtsvx_obj->trans, dgtsvx_obj->n,
                                  dgtsvx_obj->nrhs,
                                  dgtsvx_obj->dlref, dgtsvx_obj->dref, 
                                  dgtsvx_obj->duref, dgtsvx_obj->dlfref,
                                  dgtsvx_obj->dfref, dgtsvx_obj->dufref,
                                  dgtsvx_obj->du2ref, dgtsvx_obj->ipivref,
                                  dgtsvx_obj->bref, dgtsvx_obj->ldb,
                                  dgtsvx_obj->xref, dgtsvx_obj->ldx,
                                  &dgtsvx_obj->rcondref, 
                                  dgtsvx_obj->ferrref,
                                  dgtsvx_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    dgtsvx_obj->info = LAPACKE_dgtsvx( dgtsvx_obj->matrix_layout, dgtsvx_obj->fact,
                                  dgtsvx_obj->trans, dgtsvx_obj->n,
                                  dgtsvx_obj->nrhs,
                                  dgtsvx_obj->dl, dgtsvx_obj->d, 
                                  dgtsvx_obj->du, dgtsvx_obj->dlf,
                                  dgtsvx_obj->df, dgtsvx_obj->duf,
                                  dgtsvx_obj->du2, dgtsvx_obj->ipiv,
                                  dgtsvx_obj->b, dgtsvx_obj->ldb,
                                  dgtsvx_obj->x, dgtsvx_obj->ldx,
                                  &dgtsvx_obj->rcond, 
                                  dgtsvx_obj->ferr,
                                  dgtsvx_obj->berr);

    if( dgtsvx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dgtsvx is wrong\n", dgtsvx_obj->info );
    }
    if( dgtsvx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dgtsvx is wrong\n", 
        dgtsvx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    if( (dgtsvx_obj->inforef ==0) || (dgtsvx_obj->inforef == dgtsvx_obj->n+1) ){
       dgtsvx_obj->diff_xerr =  computeDiff_d( dgtsvx_obj->x_bufsize, 
                 dgtsvx_obj->x, dgtsvx_obj->xref );
    }

    dgtsvx_obj->diff_dlf =  computeDiff_d( dgtsvx_obj->n - 1, 
                dgtsvx_obj->dlf, dgtsvx_obj->dlfref );

    dgtsvx_obj->diff_df =  computeDiff_d( dgtsvx_obj->n, 
                dgtsvx_obj->df, dgtsvx_obj->dfref );

    dgtsvx_obj->diff_dlf =  computeDiff_d( dgtsvx_obj->n - 1, 
                dgtsvx_obj->duf, dgtsvx_obj->dufref );

    dgtsvx_obj->diff_du2 =  computeDiff_d( dgtsvx_obj->n - 2, 
                dgtsvx_obj->du2, dgtsvx_obj->du2ref );

    dgtsvx_obj->diff_piv =  computeDiff_i( dgtsvx_obj->n, 
                dgtsvx_obj->ipiv, dgtsvx_obj->ipivref );

    dgtsvx_obj->diff_berr =  computeDiff_d( dgtsvx_obj->nrhs, 
                dgtsvx_obj->berr, dgtsvx_obj->berrref );
                
    dgtsvx_obj->diff_ferr =  computeDiff_d( dgtsvx_obj->nrhs, 
                dgtsvx_obj->ferr, dgtsvx_obj->ferrref );
}

TEST_F(dgtsvx_test, dgtsvx1) {
    EXPECT_NEAR(0.0, dgtsvx_obj->diff_dlf, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgtsvx_obj->diff_df, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgtsvx_obj->diff_duf, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgtsvx_obj->diff_du2, LAPACKE_GTEST_THRESHOLD);
    EXPECT_EQ(0, dgtsvx_obj->diff_piv);
    EXPECT_NEAR(0.0, dgtsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgtsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgtsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgtsvx_test, dgtsvx2) {
    EXPECT_NEAR(0.0, dgtsvx_obj->diff_dlf, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgtsvx_obj->diff_df, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgtsvx_obj->diff_duf, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgtsvx_obj->diff_du2, LAPACKE_GTEST_THRESHOLD);
    EXPECT_EQ(0, dgtsvx_obj->diff_piv);
    EXPECT_NEAR(0.0, dgtsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgtsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgtsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgtsvx_test, dgtsvx3) {
    EXPECT_NEAR(0.0, dgtsvx_obj->diff_dlf, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgtsvx_obj->diff_df, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgtsvx_obj->diff_duf, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgtsvx_obj->diff_du2, LAPACKE_GTEST_THRESHOLD);
    EXPECT_EQ(0, dgtsvx_obj->diff_piv);
    EXPECT_NEAR(0.0, dgtsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgtsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgtsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgtsvx_test, dgtsvx4) {
    EXPECT_NEAR(0.0, dgtsvx_obj->diff_dlf, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgtsvx_obj->diff_df, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgtsvx_obj->diff_duf, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgtsvx_obj->diff_du2, LAPACKE_GTEST_THRESHOLD);
    EXPECT_EQ(0, dgtsvx_obj->diff_piv);
    EXPECT_NEAR(0.0, dgtsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgtsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgtsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

/* Begin gtsvx_float_parameters  class definition */
class gtsvx_float_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      // capture difference between ref o/p & libflame lapacke o/p.
      float diff_berr, diff_ferr, diff_xerr;
      float diff_dlf, diff_df, diff_duf, diff_du2;
      int diff_piv;
      void *hModule, *dModule;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char fact; //  Must be 'F', 'N', or 'E'.
      char trans; //  Must be 'N' , 'T' or 'C'.

      lapack_int n; // The order of A (Square matrix).
      lapack_int nrhs; // The number of right-hand sides
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_int ldx; //  leading dimension of 'x'
      float *b, *bref; //right-hand sides for the systems of equations.
      float *d, *dref; //size (n), contains the diagonal elements of A.
      float* dl, *dlref; // size (n -1), contains the subdiagonal elements of A.
      float* du, *duref; // size (n -1), contains the superdiagonal elements of A..

      /* Input/ Output parameters */
      float *df, *dfref; // 'n' diagonal elements of the upper triangular matrix
      float *dlf, *dlfref; // the (n -1) multipliers that define the matrix L
      float *duf, *dufref; // the (n -1) elements of the first superdiagonal of U
      float* du2, *du2ref; // the (n-2) elements of the second superdiagonal of U.
      
      lapack_int *ipiv, *ipivref; // pivot buffer
      
      /* Output parameters */
      char   equed, equedref; //  Must be 'N', 'R', 'C', or 'B'.
      float rpivot, rpivotref; // reciprocal pivot growth factor.
      float rcond, rcondref; // estimate of the reciprocal condition number of the matrix A.
      float* x, *xref;
      float* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      float* berr, *berrref; // component-wise relative backward error for each solution vector xj.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      gtsvx_float_parameters ( int matrix_layout_i, char fact_i, char trans_i,
                               char equed_i, lapack_int n_i, lapack_int nrhs_i);
              
      ~gtsvx_float_parameters (); 
};  /* end of gtsvx_float_parameters  class definition */


/* Constructor gtsvx_float_parameters definition */
gtsvx_float_parameters:: gtsvx_float_parameters ( int matrix_layout_i, 
                char fact_i, char trans_i, char equed_i, lapack_int n_i,
                  lapack_int nrhs_i) {
                                    
    matrix_layout = matrix_layout_i;
    fact= fact_i;
    trans = trans_i;
    // equed acts as input when " fact = 'F' " else output
    equed = equed_i;
    equedref = equed_i;
    n = n_i;
    nrhs = nrhs_i;
    
    b_bufsize = n*nrhs;
    x_bufsize = b_bufsize;

    if(matrix_layout==LAPACK_COL_MAJOR){
        ldb = n;
        ldx = n;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        ldb = nrhs;
        ldx = nrhs;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }
#if LAPACKE_TEST_VERBOSE
   printf(" \n gtsvx Double: matrix_layout: %d  n: %d, fact: %c trans: %c   \
ldb: %d nrhs: %d  ldaf: %d  ldx: %d \n", matrix_layout, n, fact, trans,  
                                          ldb, nrhs);
#endif

    hModule = NULL;
    dModule = NULL;
    diff_berr = 0;
    diff_ferr = 0;
    diff_xerr = 0;
    diff_dlf = 0;
    diff_df = 0;
    diff_duf = 0;
    diff_du2 = 0;
    diff_piv = 0;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_float_buffer_pair( &d, &dref, n    );
    lapacke_gtest_alloc_float_buffer_pair( &dl, &dlref, n-1  );
    lapacke_gtest_alloc_float_buffer_pair( &du, &duref, n-1  );
    lapacke_gtest_alloc_float_buffer_pair( &dlf, &dlfref, n-1  );
    lapacke_gtest_alloc_float_buffer_pair( &df, &dfref, n    );
    lapacke_gtest_alloc_float_buffer_pair( &duf, &dufref, n-1  );
    lapacke_gtest_alloc_float_buffer_pair( &du2, &du2ref, n-2  );
 
    lapacke_gtest_alloc_float_buffer_pair( &x, &xref, x_bufsize);
    lapacke_gtest_alloc_float_buffer_pair( &ferr, &ferrref, nrhs);
    lapacke_gtest_alloc_float_buffer_pair( &berr, &berrref, nrhs);
    
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (d==NULL) || (dref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (x==NULL) || (xref==NULL) || \
        (df==NULL) || (dfref==NULL) || \
        (dl==NULL) || (dlref==NULL) || \
        (du==NULL) || (duref==NULL) || \
        (dlf==NULL) || (dlfref==NULL) || \
        (duf==NULL) || (dufref==NULL) || \
        (du2==NULL) || (du2ref==NULL) || \
        (ferr==NULL) || (ferrref==NULL) || \
        (berr==NULL) || (berrref==NULL) || \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "pbtrs_float_parameters object: malloc error.";
       gtsvx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( d, dref, n);
    lapacke_gtest_init_float_buffer_pair_rand( dl, dlref, n-1);
    lapacke_gtest_init_float_buffer_pair_rand( du, duref, n-1);
    lapacke_gtest_init_float_buffer_pair_rand( dlf, dlfref, n-1);
    lapacke_gtest_init_float_buffer_pair_rand( df, dfref, n);
    lapacke_gtest_init_float_buffer_pair_rand( duf, dufref, n-1);
    lapacke_gtest_init_float_buffer_pair_rand( du2, du2ref, n-2);


    
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, b_bufsize);
    
    lapacke_gtest_init_float_buffer_pair_with_constant(x, xref, x_bufsize, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);
    
   } /* end of Constructor  */

gtsvx_float_parameters:: ~gtsvx_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gtsvx_float_parameters object: destructor invoked. \n");
#endif
   gtsvx_free();
}

//  Test fixture class definition
class sgtsvx_test  : public  ::testing::Test {
public:
   gtsvx_float_parameters  *sgtsvx_obj;
   void SetUp();  
   void TearDown () { delete sgtsvx_obj; }
};


void sgtsvx_test::SetUp(){

    /* LAPACKE SGTSVX prototype */
    typedef int (*Fptr_NL_LAPACKE_sgtsvx) (int matrix_layout, char fact,
            char trans, lapack_int n, lapack_int nrhs, const float* dl,
            const float* d, const float* du, float* dlf, float* df,
            float* duf, float* du2, lapack_int* ipiv, const float* b,
            lapack_int ldb, float* x, lapack_int ldx, float* rcond,
            float* ferr, float* berr);

    Fptr_NL_LAPACKE_sgtsvx SGTSVX;


     /* LAPACKE DGETRF prototype */
    typedef int (*Fptr_NL_LAPACKE_sgttrf) ( lapack_int n , float *dl,
            float *d , float *du , float *du2 , lapack_int *ipiv );

    Fptr_NL_LAPACKE_sgttrf SGTTRF;


    sgtsvx_obj = new gtsvx_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_driver_paramslist[idx].fact,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].equed,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    sgtsvx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgtsvx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgtsvx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgtsvx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    SGTSVX = (Fptr_NL_LAPACKE_sgtsvx)dlsym(sgtsvx_obj->hModule, "LAPACKE_sgtsvx");
    ASSERT_TRUE(SGTSVX != NULL) << "failed to get the Netlib LAPACKE_sgtsvx symbol";
    
    /* if  'fact' == 'F', then the i/p matrix A should be in 
       factored form.
       so, invoke the getrf API to compute the factorized A.  */
    if(sgtsvx_obj->fact == 'F') {
        SGTTRF = (Fptr_NL_LAPACKE_sgttrf)dlsym(sgtsvx_obj->hModule,"LAPACKE_sgttrf");
        ASSERT_TRUE(SGTTRF != NULL) << "failed to get the Netlib LAPACKE_sgttrf symbol";
            
        sgtsvx_obj->inforef = SGTTRF( sgtsvx_obj->n,
                                      sgtsvx_obj->dlref, sgtsvx_obj->dref,
                                      sgtsvx_obj->duref,
                                      sgtsvx_obj->du2ref,
                                      sgtsvx_obj->ipivref);
                               
        sgtsvx_obj->info = LAPACKE_sgttrf( sgtsvx_obj->n,
                                      sgtsvx_obj->dl, sgtsvx_obj->d,
                                      sgtsvx_obj->du,
                                      sgtsvx_obj->du2,
                                      sgtsvx_obj->ipiv);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    sgtsvx_obj->inforef = SGTSVX( sgtsvx_obj->matrix_layout, sgtsvx_obj->fact,
                                  sgtsvx_obj->trans, sgtsvx_obj->n,
                                  sgtsvx_obj->nrhs,
                                  sgtsvx_obj->dlref, sgtsvx_obj->dref, 
                                  sgtsvx_obj->duref, sgtsvx_obj->dlfref,
                                  sgtsvx_obj->dfref, sgtsvx_obj->dufref,
                                  sgtsvx_obj->du2ref, sgtsvx_obj->ipivref,
                                  sgtsvx_obj->bref, sgtsvx_obj->ldb,
                                  sgtsvx_obj->xref, sgtsvx_obj->ldx,
                                  &sgtsvx_obj->rcondref, 
                                  sgtsvx_obj->ferrref,
                                  sgtsvx_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    sgtsvx_obj->info = LAPACKE_sgtsvx( sgtsvx_obj->matrix_layout, sgtsvx_obj->fact,
                                  sgtsvx_obj->trans, sgtsvx_obj->n,
                                  sgtsvx_obj->nrhs,
                                  sgtsvx_obj->dl, sgtsvx_obj->d, 
                                  sgtsvx_obj->du, sgtsvx_obj->dlf,
                                  sgtsvx_obj->df, sgtsvx_obj->duf,
                                  sgtsvx_obj->du2, sgtsvx_obj->ipiv,
                                  sgtsvx_obj->b, sgtsvx_obj->ldb,
                                  sgtsvx_obj->x, sgtsvx_obj->ldx,
                                  &sgtsvx_obj->rcond, 
                                  sgtsvx_obj->ferr,
                                  sgtsvx_obj->berr);

    if( sgtsvx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_sgtsvx is wrong\n", sgtsvx_obj->info );
    }
    if( sgtsvx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sgtsvx is wrong\n", 
        sgtsvx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    if( (sgtsvx_obj->inforef ==0) || (sgtsvx_obj->inforef == sgtsvx_obj->n+1) ){
       sgtsvx_obj->diff_xerr =  computeDiff_s( sgtsvx_obj->x_bufsize, 
                 sgtsvx_obj->x, sgtsvx_obj->xref );
    }

    sgtsvx_obj->diff_dlf =  computeDiff_s( sgtsvx_obj->n - 1, 
                sgtsvx_obj->dlf, sgtsvx_obj->dlfref );

    sgtsvx_obj->diff_df =  computeDiff_s( sgtsvx_obj->n, 
                sgtsvx_obj->df, sgtsvx_obj->dfref );

    sgtsvx_obj->diff_dlf =  computeDiff_s( sgtsvx_obj->n - 1, 
                sgtsvx_obj->duf, sgtsvx_obj->dufref );

    sgtsvx_obj->diff_du2 =  computeDiff_s( sgtsvx_obj->n - 2, 
                sgtsvx_obj->du2, sgtsvx_obj->du2ref );

    sgtsvx_obj->diff_piv =  computeDiff_i( sgtsvx_obj->n, 
                sgtsvx_obj->ipiv, sgtsvx_obj->ipivref );

    sgtsvx_obj->diff_berr =  computeDiff_s( sgtsvx_obj->nrhs, 
                sgtsvx_obj->berr, sgtsvx_obj->berrref );
                
    sgtsvx_obj->diff_ferr =  computeDiff_s( sgtsvx_obj->nrhs, 
                sgtsvx_obj->ferr, sgtsvx_obj->ferrref );
}

TEST_F(sgtsvx_test, sgtsvx1) {
    EXPECT_NEAR(0.0, sgtsvx_obj->diff_dlf, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgtsvx_obj->diff_df, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgtsvx_obj->diff_duf, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgtsvx_obj->diff_du2, LAPACKE_GTEST_THRESHOLD);
    EXPECT_EQ(0, sgtsvx_obj->diff_piv);
    EXPECT_NEAR(0.0, sgtsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgtsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgtsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgtsvx_test, sgtsvx2) {
    EXPECT_NEAR(0.0, sgtsvx_obj->diff_dlf, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgtsvx_obj->diff_df, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgtsvx_obj->diff_duf, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgtsvx_obj->diff_du2, LAPACKE_GTEST_THRESHOLD);
    EXPECT_EQ(0, sgtsvx_obj->diff_piv);
    EXPECT_NEAR(0.0, sgtsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgtsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgtsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgtsvx_test, sgtsvx3) {
    EXPECT_NEAR(0.0, sgtsvx_obj->diff_dlf, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgtsvx_obj->diff_df, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgtsvx_obj->diff_duf, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgtsvx_obj->diff_du2, LAPACKE_GTEST_THRESHOLD);
    EXPECT_EQ(0, sgtsvx_obj->diff_piv);
    EXPECT_NEAR(0.0, sgtsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgtsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgtsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgtsvx_test, sgtsvx4) {
    EXPECT_NEAR(0.0, sgtsvx_obj->diff_dlf, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgtsvx_obj->diff_df, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgtsvx_obj->diff_duf, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgtsvx_obj->diff_du2, LAPACKE_GTEST_THRESHOLD);
    EXPECT_EQ(0, sgtsvx_obj->diff_piv);
    EXPECT_NEAR(0.0, sgtsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgtsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgtsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

/* Begin gtsvx_scomplex_parameters  class definition */
class gtsvx_scomplex_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      // capture difference between ref o/p & libflame lapacke o/p.
      float diff_berr, diff_ferr, diff_xerr;
      float diff_dlf, diff_df, diff_duf, diff_du2;
      int diff_piv;
      void *hModule, *dModule;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char fact; //  Must be 'F', 'N', or 'E'.
      char trans; //  Must be 'N' , 'T' or 'C'.

      lapack_int n; // The order of A (Square matrix).
      lapack_int nrhs; // The number of right-hand sides
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_int ldx; //  leading dimension of 'x'
      lapack_complex_float *b, *bref; //right-hand sides for the systems of equations.
      lapack_complex_float *d, *dref; //size (n), contains the diagonal elements of A.
      lapack_complex_float* dl, *dlref; // size (n -1), contains the subdiagonal elements of A.
      lapack_complex_float* du, *duref; // size (n -1), contains the superdiagonal elements of A..

      /* Input/ Output parameters */
      lapack_complex_float *df, *dfref; // 'n' diagonal elements of the upper triangular matrix
      lapack_complex_float *dlf, *dlfref; // the (n -1) multipliers that define the matrix L
      lapack_complex_float *duf, *dufref; // the (n -1) elements of the first superdiagonal of U
      lapack_complex_float* du2, *du2ref; // the (n-2) elements of the second superdiagonal of U.
      
      lapack_int *ipiv, *ipivref; // pivot buffer
      
      /* Output parameters */
      char   equed, equedref; //  Must be 'N', 'R', 'C', or 'B'.
      lapack_complex_float rpivot, rpivotref; // reciprocal pivot growth factor.
      float rcond, rcondref; // estimate of the reciprocal condition number of the matrix A.
      lapack_complex_float* x, *xref;
      float* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      float* berr, *berrref; // component-wise relative backward error for each solution vector xj.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      gtsvx_scomplex_parameters ( int matrix_layout_i, char fact_i, char trans_i,
                               char equed_i, lapack_int n_i, lapack_int nrhs_i);
              
      ~gtsvx_scomplex_parameters (); 
};  /* end of gtsvx_scomplex_parameters  class definition */


/* Constructor gtsvx_scomplex_parameters definition */
gtsvx_scomplex_parameters:: gtsvx_scomplex_parameters ( int matrix_layout_i, 
                char fact_i, char trans_i, char equed_i, lapack_int n_i,
                  lapack_int nrhs_i) {
                                    
    matrix_layout = matrix_layout_i;
    fact= fact_i;
    trans = trans_i;
    // equed acts as input when " fact = 'F' " else output
    equed = equed_i;
    equedref = equed_i;
    n = n_i;
    nrhs = nrhs_i;
    
    b_bufsize = n*nrhs;
    x_bufsize = b_bufsize;

    if(matrix_layout==LAPACK_COL_MAJOR){
        ldb = n;
        ldx = n;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        ldb = nrhs;
        ldx = nrhs;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }
#if LAPACKE_TEST_VERBOSE
   printf(" \n gtsvx Double: matrix_layout: %d  n: %d, fact: %c trans: %c   \
ldb: %d nrhs: %d  ldaf: %d  ldx: %d \n", matrix_layout, n, fact, trans,  
                                          ldb, nrhs);
#endif

    hModule = NULL;
    dModule = NULL;
    diff_berr = 0;
    diff_ferr = 0;
    diff_xerr = 0;
    diff_dlf = 0;
    diff_df = 0;
    diff_duf = 0;
    diff_du2 = 0;
    diff_piv = 0;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &d, &dref, n    );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &dl, &dlref, n-1  );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &du, &duref, n-1  );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &dlf, &dlfref, n-1  );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &df, &dfref, n    );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &duf, &dufref, n-1  );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &du2, &du2ref, n-2  ); 
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &x, &xref, x_bufsize);
    
    lapacke_gtest_alloc_float_buffer_pair( &ferr, &ferrref, nrhs);
    lapacke_gtest_alloc_float_buffer_pair( &berr, &berrref, nrhs);
    
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (d==NULL) || (dref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (x==NULL) || (xref==NULL) || \
        (df==NULL) || (dfref==NULL) || \
        (dl==NULL) || (dlref==NULL) || \
        (du==NULL) || (duref==NULL) || \
        (dlf==NULL) || (dlfref==NULL) || \
        (duf==NULL) || (dufref==NULL) || \
        (du2==NULL) || (du2ref==NULL) || \
        (ferr==NULL) || (ferrref==NULL) || \
        (berr==NULL) || (berrref==NULL) || \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "pbtrs_scomplex_parameters object: malloc error.";
       gtsvx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( d, dref, n);
    lapacke_gtest_init_scomplex_buffer_pair_rand( dl, dlref, n-1);
    lapacke_gtest_init_scomplex_buffer_pair_rand( du, duref, n-1);
    lapacke_gtest_init_scomplex_buffer_pair_rand( dlf, dlfref, n-1);
    lapacke_gtest_init_scomplex_buffer_pair_rand( df, dfref, n);
    lapacke_gtest_init_scomplex_buffer_pair_rand( duf, dufref, n-1);
    lapacke_gtest_init_scomplex_buffer_pair_rand( du2, du2ref, n-2);


    
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
    
    lapacke_gtest_init_scomplex_buffer_pair_with_constant(x, xref, x_bufsize, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);
    
   } /* end of Constructor  */

gtsvx_scomplex_parameters:: ~gtsvx_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gtsvx_scomplex_parameters object: destructor invoked. \n");
#endif
   gtsvx_free();
}

//  Test fixture class definition
class cgtsvx_test  : public  ::testing::Test {
public:
   gtsvx_scomplex_parameters  *cgtsvx_obj;
   void SetUp();  
   void TearDown () { delete cgtsvx_obj; }
};


void cgtsvx_test::SetUp(){

    /* LAPACKE CGTSVX prototype */
    typedef int (*Fptr_NL_LAPACKE_cgtsvx) (int matrix_layout, char fact,
            char trans, lapack_int n, lapack_int nrhs, const lapack_complex_float* dl,
            const lapack_complex_float* d, const lapack_complex_float* du, lapack_complex_float* dlf, lapack_complex_float* df,
            lapack_complex_float* duf, lapack_complex_float* du2, lapack_int* ipiv, const lapack_complex_float* b,
            lapack_int ldb, lapack_complex_float* x, lapack_int ldx, float* rcond,
            float* ferr, float* berr);

    Fptr_NL_LAPACKE_cgtsvx CGTSVX;


     /* LAPACKE DGETRF prototype */
    typedef int (*Fptr_NL_LAPACKE_cgttrf) ( lapack_int n , lapack_complex_float *dl,
            lapack_complex_float *d , lapack_complex_float *du , lapack_complex_float *du2 , lapack_int *ipiv );

    Fptr_NL_LAPACKE_cgttrf CGTTRF;


    cgtsvx_obj = new gtsvx_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_driver_paramslist[idx].fact,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].equed,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    cgtsvx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgtsvx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgtsvx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgtsvx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    CGTSVX = (Fptr_NL_LAPACKE_cgtsvx)dlsym(cgtsvx_obj->hModule, "LAPACKE_cgtsvx");
    ASSERT_TRUE(CGTSVX != NULL) << "failed to get the Netlib LAPACKE_cgtsvx symbol";
    
    /* if  'fact' == 'F', then the i/p matrix A should be in 
       factored form.
       so, invoke the getrf API to compute the factorized A.  */
    if(cgtsvx_obj->fact == 'F') {
        CGTTRF = (Fptr_NL_LAPACKE_cgttrf)dlsym(cgtsvx_obj->hModule,"LAPACKE_cgttrf");
        ASSERT_TRUE(CGTTRF != NULL) << "failed to get the Netlib LAPACKE_cgttrf symbol";
            
        cgtsvx_obj->inforef = CGTTRF( cgtsvx_obj->n,
                                      cgtsvx_obj->dlref, cgtsvx_obj->dref,
                                      cgtsvx_obj->duref,
                                      cgtsvx_obj->du2ref,
                                      cgtsvx_obj->ipivref);
                               
        cgtsvx_obj->info = LAPACKE_cgttrf( cgtsvx_obj->n,
                                      cgtsvx_obj->dl, cgtsvx_obj->d,
                                      cgtsvx_obj->du,
                                      cgtsvx_obj->du2,
                                      cgtsvx_obj->ipiv);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    cgtsvx_obj->inforef = CGTSVX( cgtsvx_obj->matrix_layout, cgtsvx_obj->fact,
                                  cgtsvx_obj->trans, cgtsvx_obj->n,
                                  cgtsvx_obj->nrhs,
                                  cgtsvx_obj->dlref, cgtsvx_obj->dref, 
                                  cgtsvx_obj->duref, cgtsvx_obj->dlfref,
                                  cgtsvx_obj->dfref, cgtsvx_obj->dufref,
                                  cgtsvx_obj->du2ref, cgtsvx_obj->ipivref,
                                  cgtsvx_obj->bref, cgtsvx_obj->ldb,
                                  cgtsvx_obj->xref, cgtsvx_obj->ldx,
                                  &cgtsvx_obj->rcondref, 
                                  cgtsvx_obj->ferrref,
                                  cgtsvx_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    cgtsvx_obj->info = LAPACKE_cgtsvx( cgtsvx_obj->matrix_layout, cgtsvx_obj->fact,
                                  cgtsvx_obj->trans, cgtsvx_obj->n,
                                  cgtsvx_obj->nrhs,
                                  cgtsvx_obj->dl, cgtsvx_obj->d, 
                                  cgtsvx_obj->du, cgtsvx_obj->dlf,
                                  cgtsvx_obj->df, cgtsvx_obj->duf,
                                  cgtsvx_obj->du2, cgtsvx_obj->ipiv,
                                  cgtsvx_obj->b, cgtsvx_obj->ldb,
                                  cgtsvx_obj->x, cgtsvx_obj->ldx,
                                  &cgtsvx_obj->rcond, 
                                  cgtsvx_obj->ferr,
                                  cgtsvx_obj->berr);

    if( cgtsvx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cgtsvx is wrong\n", cgtsvx_obj->info );
    }
    if( cgtsvx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cgtsvx is wrong\n", 
        cgtsvx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    if( (cgtsvx_obj->inforef ==0) || (cgtsvx_obj->inforef == cgtsvx_obj->n+1) ){
       cgtsvx_obj->diff_xerr =  computeDiff_c( cgtsvx_obj->x_bufsize, 
                 cgtsvx_obj->x, cgtsvx_obj->xref );
    }

    cgtsvx_obj->diff_dlf =  computeDiff_c( cgtsvx_obj->n - 1, 
                cgtsvx_obj->dlf, cgtsvx_obj->dlfref );

    cgtsvx_obj->diff_df =  computeDiff_c( cgtsvx_obj->n, 
                cgtsvx_obj->df, cgtsvx_obj->dfref );

    cgtsvx_obj->diff_dlf =  computeDiff_c( cgtsvx_obj->n - 1, 
                cgtsvx_obj->duf, cgtsvx_obj->dufref );

    cgtsvx_obj->diff_du2 =  computeDiff_c( cgtsvx_obj->n - 2, 
                cgtsvx_obj->du2, cgtsvx_obj->du2ref );

    cgtsvx_obj->diff_piv =  computeDiff_i( cgtsvx_obj->n, 
                cgtsvx_obj->ipiv, cgtsvx_obj->ipivref );

    cgtsvx_obj->diff_berr =  computeDiff_s( cgtsvx_obj->nrhs, 
                cgtsvx_obj->berr, cgtsvx_obj->berrref );
                
    cgtsvx_obj->diff_ferr =  computeDiff_s( cgtsvx_obj->nrhs, 
                cgtsvx_obj->ferr, cgtsvx_obj->ferrref );
}

TEST_F(cgtsvx_test, cgtsvx1) {
    EXPECT_NEAR(0.0, cgtsvx_obj->diff_dlf, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgtsvx_obj->diff_df, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgtsvx_obj->diff_duf, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgtsvx_obj->diff_du2, LAPACKE_GTEST_THRESHOLD);
    EXPECT_EQ(0, cgtsvx_obj->diff_piv);
    EXPECT_NEAR(0.0, cgtsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgtsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgtsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgtsvx_test, cgtsvx2) {
    EXPECT_NEAR(0.0, cgtsvx_obj->diff_dlf, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgtsvx_obj->diff_df, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgtsvx_obj->diff_duf, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgtsvx_obj->diff_du2, LAPACKE_GTEST_THRESHOLD);
    EXPECT_EQ(0, cgtsvx_obj->diff_piv);
    EXPECT_NEAR(0.0, cgtsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgtsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgtsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgtsvx_test, cgtsvx3) {
    EXPECT_NEAR(0.0, cgtsvx_obj->diff_dlf, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgtsvx_obj->diff_df, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgtsvx_obj->diff_duf, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgtsvx_obj->diff_du2, LAPACKE_GTEST_THRESHOLD);
    EXPECT_EQ(0, cgtsvx_obj->diff_piv);
    EXPECT_NEAR(0.0, cgtsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgtsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgtsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgtsvx_test, cgtsvx4) {
    EXPECT_NEAR(0.0, cgtsvx_obj->diff_dlf, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgtsvx_obj->diff_df, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgtsvx_obj->diff_duf, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgtsvx_obj->diff_du2, LAPACKE_GTEST_THRESHOLD);
    EXPECT_EQ(0, cgtsvx_obj->diff_piv);
    EXPECT_NEAR(0.0, cgtsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgtsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgtsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

/* Begin gtsvx_dcomplex_parameters  class definition */
class gtsvx_dcomplex_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      // capture difference between ref o/p & libflame lapacke o/p.
      double diff_berr, diff_ferr, diff_xerr;
      double diff_dlf, diff_df, diff_duf, diff_du2;
      int diff_piv;
      void *hModule, *dModule;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char fact; //  Must be 'F', 'N', or 'E'.
      char trans; //  Must be 'N' , 'T' or 'C'.

      lapack_int n; // The order of A (Square matrix).
      lapack_int nrhs; // The number of right-hand sides
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_int ldx; //  leading dimension of 'x'
      lapack_complex_double *b, *bref; //right-hand sides for the systems of equations.
      lapack_complex_double *d, *dref; //size (n), contains the diagonal elements of A.
      lapack_complex_double* dl, *dlref; // size (n -1), contains the subdiagonal elements of A.
      lapack_complex_double* du, *duref; // size (n -1), contains the superdiagonal elements of A..

      /* Input/ Output parameters */
      lapack_complex_double *df, *dfref; // 'n' diagonal elements of the upper triangular matrix
      lapack_complex_double *dlf, *dlfref; // the (n -1) multipliers that define the matrix L
      lapack_complex_double *duf, *dufref; // the (n -1) elements of the first superdiagonal of U
      lapack_complex_double* du2, *du2ref; // the (n-2) elements of the second superdiagonal of U.
      
      lapack_int *ipiv, *ipivref; // pivot buffer
      
      /* Output parameters */
      char   equed, equedref; //  Must be 'N', 'R', 'C', or 'B'.
      lapack_complex_double rpivot, rpivotref; // reciprocal pivot growth factor.
      double rcond, rcondref; // estimate of the reciprocal condition number of the matrix A.
      lapack_complex_double* x, *xref;
      double* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      double* berr, *berrref; // component-wise relative backward error for each solution vector xj.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      gtsvx_dcomplex_parameters ( int matrix_layout_i, char fact_i, char trans_i,
                               char equed_i, lapack_int n_i, lapack_int nrhs_i);
              
      ~gtsvx_dcomplex_parameters (); 
};  /* end of gtsvx_dcomplex_parameters  class definition */


/* Constructor gtsvx_dcomplex_parameters definition */
gtsvx_dcomplex_parameters:: gtsvx_dcomplex_parameters ( int matrix_layout_i, 
                char fact_i, char trans_i, char equed_i, lapack_int n_i,
                  lapack_int nrhs_i) {
                                    
    matrix_layout = matrix_layout_i;
    fact= fact_i;
    trans = trans_i;
    // equed acts as input when " fact = 'F' " else output
    equed = equed_i;
    equedref = equed_i;
    n = n_i;
    nrhs = nrhs_i;
    
    b_bufsize = n*nrhs;
    x_bufsize = b_bufsize;

    if(matrix_layout==LAPACK_COL_MAJOR){
        ldb = n;
        ldx = n;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        ldb = nrhs;
        ldx = nrhs;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }
#if LAPACKE_TEST_VERBOSE
   printf(" \n gtsvx Double: matrix_layout: %d  n: %d, fact: %c trans: %c   \
ldb: %d nrhs: %d  ldaf: %d  ldx: %d \n", matrix_layout, n, fact, trans,  
                                          ldb, nrhs);
#endif

    hModule = NULL;
    dModule = NULL;
    diff_berr = 0;
    diff_ferr = 0;
    diff_xerr = 0;
    diff_dlf = 0;
    diff_df = 0;
    diff_duf = 0;
    diff_du2 = 0;
    diff_piv = 0;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &d, &dref, n    );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &dl, &dlref, n-1  );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &du, &duref, n-1  );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &dlf, &dlfref, n-1  );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &df, &dfref, n    );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &duf, &dufref, n-1  );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &du2, &du2ref, n-2  ); 
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &x, &xref, x_bufsize);
    
    lapacke_gtest_alloc_double_buffer_pair( &ferr, &ferrref, nrhs);
    lapacke_gtest_alloc_double_buffer_pair( &berr, &berrref, nrhs);
    
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (d==NULL) || (dref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (x==NULL) || (xref==NULL) || \
        (df==NULL) || (dfref==NULL) || \
        (dl==NULL) || (dlref==NULL) || \
        (du==NULL) || (duref==NULL) || \
        (dlf==NULL) || (dlfref==NULL) || \
        (duf==NULL) || (dufref==NULL) || \
        (du2==NULL) || (du2ref==NULL) || \
        (ferr==NULL) || (ferrref==NULL) || \
        (berr==NULL) || (berrref==NULL) || \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "pbtrs_dcomplex_parameters object: malloc error.";
       gtsvx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( d, dref, n);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( dl, dlref, n-1);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( du, duref, n-1);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( dlf, dlfref, n-1);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( df, dfref, n);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( duf, dufref, n-1);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( du2, du2ref, n-2);


    
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
    
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant(x, xref, x_bufsize, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);
    
   } /* end of Constructor  */

gtsvx_dcomplex_parameters:: ~gtsvx_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gtsvx_dcomplex_parameters object: destructor invoked. \n");
#endif
   gtsvx_free();
}

//  Test fixture class definition
class zgtsvx_test  : public  ::testing::Test {
public:
   gtsvx_dcomplex_parameters  *zgtsvx_obj;
   void SetUp();  
   void TearDown () { delete zgtsvx_obj; }
};


void zgtsvx_test::SetUp(){

    /* LAPACKE ZGTSVX prototype */
    typedef int (*Fptr_NL_LAPACKE_zgtsvx) (int matrix_layout, char fact,
            char trans, lapack_int n, lapack_int nrhs, const lapack_complex_double* dl,
            const lapack_complex_double* d, const lapack_complex_double* du, lapack_complex_double* dlf, lapack_complex_double* df,
            lapack_complex_double* duf, lapack_complex_double* du2, lapack_int* ipiv, const lapack_complex_double* b,
            lapack_int ldb, lapack_complex_double* x, lapack_int ldx, double* rcond,
            double* ferr, double* berr);

    Fptr_NL_LAPACKE_zgtsvx ZGTSVX;


     /* LAPACKE DGETRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zgttrf) ( lapack_int n , lapack_complex_double *dl,
            lapack_complex_double *d , lapack_complex_double *du , lapack_complex_double *du2 , lapack_int *ipiv );

    Fptr_NL_LAPACKE_zgttrf ZGTTRF;


    zgtsvx_obj = new gtsvx_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_driver_paramslist[idx].fact,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].equed,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    zgtsvx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgtsvx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgtsvx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgtsvx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    ZGTSVX = (Fptr_NL_LAPACKE_zgtsvx)dlsym(zgtsvx_obj->hModule, "LAPACKE_zgtsvx");
    ASSERT_TRUE(ZGTSVX != NULL) << "failed to get the Netlib LAPACKE_zgtsvx symbol";
    
    /* if  'fact' == 'F', then the i/p matrix A should be in 
       factored form.
       so, invoke the getrf API to compute the factorized A.  */
    if(zgtsvx_obj->fact == 'F') {
        ZGTTRF = (Fptr_NL_LAPACKE_zgttrf)dlsym(zgtsvx_obj->hModule,"LAPACKE_zgttrf");
        ASSERT_TRUE(ZGTTRF != NULL) << "failed to get the Netlib LAPACKE_zgttrf symbol";
            
        zgtsvx_obj->inforef = ZGTTRF( zgtsvx_obj->n,
                                      zgtsvx_obj->dlref, zgtsvx_obj->dref,
                                      zgtsvx_obj->duref,
                                      zgtsvx_obj->du2ref,
                                      zgtsvx_obj->ipivref);
                               
        zgtsvx_obj->info = LAPACKE_zgttrf( zgtsvx_obj->n,
                                      zgtsvx_obj->dl, zgtsvx_obj->d,
                                      zgtsvx_obj->du,
                                      zgtsvx_obj->du2,
                                      zgtsvx_obj->ipiv);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    zgtsvx_obj->inforef = ZGTSVX( zgtsvx_obj->matrix_layout, zgtsvx_obj->fact,
                                  zgtsvx_obj->trans, zgtsvx_obj->n,
                                  zgtsvx_obj->nrhs,
                                  zgtsvx_obj->dlref, zgtsvx_obj->dref, 
                                  zgtsvx_obj->duref, zgtsvx_obj->dlfref,
                                  zgtsvx_obj->dfref, zgtsvx_obj->dufref,
                                  zgtsvx_obj->du2ref, zgtsvx_obj->ipivref,
                                  zgtsvx_obj->bref, zgtsvx_obj->ldb,
                                  zgtsvx_obj->xref, zgtsvx_obj->ldx,
                                  &zgtsvx_obj->rcondref, 
                                  zgtsvx_obj->ferrref,
                                  zgtsvx_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    zgtsvx_obj->info = LAPACKE_zgtsvx( zgtsvx_obj->matrix_layout, zgtsvx_obj->fact,
                                  zgtsvx_obj->trans, zgtsvx_obj->n,
                                  zgtsvx_obj->nrhs,
                                  zgtsvx_obj->dl, zgtsvx_obj->d, 
                                  zgtsvx_obj->du, zgtsvx_obj->dlf,
                                  zgtsvx_obj->df, zgtsvx_obj->duf,
                                  zgtsvx_obj->du2, zgtsvx_obj->ipiv,
                                  zgtsvx_obj->b, zgtsvx_obj->ldb,
                                  zgtsvx_obj->x, zgtsvx_obj->ldx,
                                  &zgtsvx_obj->rcond, 
                                  zgtsvx_obj->ferr,
                                  zgtsvx_obj->berr);

    if( zgtsvx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zgtsvx is wrong\n", zgtsvx_obj->info );
    }
    if( zgtsvx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zgtsvx is wrong\n", 
        zgtsvx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    if( (zgtsvx_obj->inforef ==0) || (zgtsvx_obj->inforef == zgtsvx_obj->n+1) ){
       zgtsvx_obj->diff_xerr =  computeDiff_z( zgtsvx_obj->x_bufsize, 
                 zgtsvx_obj->x, zgtsvx_obj->xref );
    }

    zgtsvx_obj->diff_dlf =  computeDiff_z( zgtsvx_obj->n - 1, 
                zgtsvx_obj->dlf, zgtsvx_obj->dlfref );

    zgtsvx_obj->diff_df =  computeDiff_z( zgtsvx_obj->n, 
                zgtsvx_obj->df, zgtsvx_obj->dfref );

    zgtsvx_obj->diff_dlf =  computeDiff_z( zgtsvx_obj->n - 1, 
                zgtsvx_obj->duf, zgtsvx_obj->dufref );

    zgtsvx_obj->diff_du2 =  computeDiff_z( zgtsvx_obj->n - 2, 
                zgtsvx_obj->du2, zgtsvx_obj->du2ref );

    zgtsvx_obj->diff_piv =  computeDiff_i( zgtsvx_obj->n, 
                zgtsvx_obj->ipiv, zgtsvx_obj->ipivref );

    zgtsvx_obj->diff_berr =  computeDiff_d( zgtsvx_obj->nrhs, 
                zgtsvx_obj->berr, zgtsvx_obj->berrref );
                
    zgtsvx_obj->diff_ferr =  computeDiff_d( zgtsvx_obj->nrhs, 
                zgtsvx_obj->ferr, zgtsvx_obj->ferrref );
}

TEST_F(zgtsvx_test, zgtsvx1) {
    EXPECT_NEAR(0.0, zgtsvx_obj->diff_dlf, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgtsvx_obj->diff_df, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgtsvx_obj->diff_duf, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgtsvx_obj->diff_du2, LAPACKE_GTEST_THRESHOLD);
    EXPECT_EQ(0, zgtsvx_obj->diff_piv);
    EXPECT_NEAR(0.0, zgtsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgtsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgtsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgtsvx_test, zgtsvx2) {
    EXPECT_NEAR(0.0, zgtsvx_obj->diff_dlf, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgtsvx_obj->diff_df, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgtsvx_obj->diff_duf, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgtsvx_obj->diff_du2, LAPACKE_GTEST_THRESHOLD);
    EXPECT_EQ(0, zgtsvx_obj->diff_piv);
    EXPECT_NEAR(0.0, zgtsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgtsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgtsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgtsvx_test, zgtsvx3) {
    EXPECT_NEAR(0.0, zgtsvx_obj->diff_dlf, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgtsvx_obj->diff_df, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgtsvx_obj->diff_duf, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgtsvx_obj->diff_du2, LAPACKE_GTEST_THRESHOLD);
    EXPECT_EQ(0, zgtsvx_obj->diff_piv);
    EXPECT_NEAR(0.0, zgtsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgtsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgtsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgtsvx_test, zgtsvx4) {
    EXPECT_NEAR(0.0, zgtsvx_obj->diff_dlf, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgtsvx_obj->diff_df, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgtsvx_obj->diff_duf, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgtsvx_obj->diff_du2, LAPACKE_GTEST_THRESHOLD);
    EXPECT_EQ(0, zgtsvx_obj->diff_piv);
    EXPECT_NEAR(0.0, zgtsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgtsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgtsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}
