#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"


#define gtrfs_free() \
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

/* Begin gtrfs_double_parameters  class definition */
class gtrfs_double_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      // capture difference between ref o/p & libflame lapacke o/p.
      double diff_berr, diff_ferr, diff_xerr;
      void *hModule, *dModule;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
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
      double* x, *xref;
      double* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      double* berr, *berrref; // component-wise relative backward error for each solution vector xj.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      gtrfs_double_parameters ( int matrix_layout_i, char trans_i,
                               char equed_i, lapack_int n_i, lapack_int nrhs_i);
              
      ~gtrfs_double_parameters (); 
};  /* end of gtrfs_double_parameters  class definition */


/* Constructor gtrfs_double_parameters definition */
gtrfs_double_parameters:: gtrfs_double_parameters ( int matrix_layout_i, 
                char trans_i, char equed_i, lapack_int n_i,
                  lapack_int nrhs_i) {
                                    
    matrix_layout = matrix_layout_i;
    trans = trans_i;
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
   printf(" \n gtrfs Double: matrix_layout: %d  n: %d,  trans: %c   \
ldb: %d nrhs: %d  ldaf: %d  ldx: %d \n", matrix_layout, n, trans,  
                                          ldb, nrhs);
#endif

    hModule = NULL;
    dModule = NULL;
    diff_berr = 0;
    diff_ferr = 0;
    diff_xerr = 0;

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
       gtrfs_free();
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

gtrfs_double_parameters:: ~gtrfs_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gtrfs_double_parameters object: destructor invoked. \n");
#endif
   gtrfs_free();
}

//  Test fixture class definition
class dgtrfs_test  : public  ::testing::Test {
public:
   gtrfs_double_parameters  *dgtrfs_obj;
   void SetUp();  
   void TearDown () { delete dgtrfs_obj; }
};


void dgtrfs_test::SetUp(){

    /* LAPACKE DGTSVX prototype */
    typedef int (*Fptr_NL_LAPACKE_dgtrfs) (int matrix_layout, 
            char trans, lapack_int n, lapack_int nrhs, const double* dl,
            const double* d, const double* du, double* dlf, double* df,
            double* duf, double* du2, lapack_int* ipiv, const double* b,
            lapack_int ldb, double* x, lapack_int ldx, 
            double* ferr, double* berr);

    Fptr_NL_LAPACKE_dgtrfs DGTSVX;


     /* LAPACKE DGETRF prototype */
    typedef int (*Fptr_NL_LAPACKE_dgttrf) ( lapack_int n , double *dl,
            double *d , double *du , double *du2 , lapack_int *ipiv );

    Fptr_NL_LAPACKE_dgttrf DGTTRF;


    dgtrfs_obj = new gtrfs_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].equed,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    dgtrfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgtrfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgtrfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgtrfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    DGTSVX = (Fptr_NL_LAPACKE_dgtrfs)dlsym(dgtrfs_obj->hModule, "LAPACKE_dgtrfs");
    ASSERT_TRUE(DGTSVX != NULL) << "failed to get the Netlib LAPACKE_dgtrfs symbol";
    
    /* invoke the getrf API to compute the factorized A.  */
        DGTTRF = (Fptr_NL_LAPACKE_dgttrf)dlsym(dgtrfs_obj->hModule,"LAPACKE_dgttrf");
        ASSERT_TRUE(DGTTRF != NULL) << "failed to get the Netlib LAPACKE_dgttrf symbol";
            
        dgtrfs_obj->inforef = DGTTRF( dgtrfs_obj->n,
                                      dgtrfs_obj->dlref, dgtrfs_obj->dref,
                                      dgtrfs_obj->duref,
                                      dgtrfs_obj->du2ref,
                                      dgtrfs_obj->ipivref);
                               
        dgtrfs_obj->info = LAPACKE_dgttrf( dgtrfs_obj->n,
                                      dgtrfs_obj->dl, dgtrfs_obj->d,
                                      dgtrfs_obj->du,
                                      dgtrfs_obj->du2,
                                      dgtrfs_obj->ipiv);

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    dgtrfs_obj->inforef = DGTSVX( dgtrfs_obj->matrix_layout,
                                  dgtrfs_obj->trans, dgtrfs_obj->n,
                                  dgtrfs_obj->nrhs,
                                  dgtrfs_obj->dlref, dgtrfs_obj->dref, 
                                  dgtrfs_obj->duref, dgtrfs_obj->dlfref,
                                  dgtrfs_obj->dfref, dgtrfs_obj->dufref,
                                  dgtrfs_obj->du2ref, dgtrfs_obj->ipivref,
                                  dgtrfs_obj->bref, dgtrfs_obj->ldb,
                                  dgtrfs_obj->xref, dgtrfs_obj->ldx,
                                  dgtrfs_obj->ferrref,
                                  dgtrfs_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    dgtrfs_obj->info = LAPACKE_dgtrfs( dgtrfs_obj->matrix_layout,
                                  dgtrfs_obj->trans, dgtrfs_obj->n,
                                  dgtrfs_obj->nrhs,
                                  dgtrfs_obj->dl, dgtrfs_obj->d, 
                                  dgtrfs_obj->du, dgtrfs_obj->dlf,
                                  dgtrfs_obj->df, dgtrfs_obj->duf,
                                  dgtrfs_obj->du2, dgtrfs_obj->ipiv,
                                  dgtrfs_obj->b, dgtrfs_obj->ldb,
                                  dgtrfs_obj->x, dgtrfs_obj->ldx,
                                  dgtrfs_obj->ferr,
                                  dgtrfs_obj->berr);

    if( dgtrfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dgtrfs is wrong\n", dgtrfs_obj->info );
    }
    if( dgtrfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dgtrfs is wrong\n", 
        dgtrfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
       dgtrfs_obj->diff_xerr =  computeDiff_d( dgtrfs_obj->x_bufsize, 
                 dgtrfs_obj->x, dgtrfs_obj->xref );

    dgtrfs_obj->diff_berr =  computeDiff_d( dgtrfs_obj->nrhs, 
                dgtrfs_obj->berr, dgtrfs_obj->berrref );
                
    dgtrfs_obj->diff_ferr =  computeDiff_d( dgtrfs_obj->nrhs, 
                dgtrfs_obj->ferr, dgtrfs_obj->ferrref );
}

TEST_F(dgtrfs_test, dgtrfs1) {
    EXPECT_NEAR(0.0, dgtrfs_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgtrfs_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgtrfs_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgtrfs_test, dgtrfs2) {
    EXPECT_NEAR(0.0, dgtrfs_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgtrfs_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgtrfs_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgtrfs_test, dgtrfs3) {
    EXPECT_NEAR(0.0, dgtrfs_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgtrfs_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgtrfs_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgtrfs_test, dgtrfs4) {
    EXPECT_NEAR(0.0, dgtrfs_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgtrfs_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dgtrfs_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

/* Begin gtrfs_float_parameters  class definition */
class gtrfs_float_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      // capture difference between ref o/p & libflame lapacke o/p.
      float diff_berr, diff_ferr, diff_xerr;
      void *hModule, *dModule;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
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
      float* x, *xref;
      float* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      float* berr, *berrref; // component-wise relative backward error for each solution vector xj.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      gtrfs_float_parameters ( int matrix_layout_i, char trans_i,
                               char equed_i, lapack_int n_i, lapack_int nrhs_i);
              
      ~gtrfs_float_parameters (); 
};  /* end of gtrfs_float_parameters  class definition */


/* Constructor gtrfs_float_parameters definition */
gtrfs_float_parameters:: gtrfs_float_parameters ( int matrix_layout_i, 
                char trans_i, char equed_i, lapack_int n_i,
                  lapack_int nrhs_i) {
                                    
    matrix_layout = matrix_layout_i;
    trans = trans_i;
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
   printf(" \n gtrfs Double: matrix_layout: %d  n: %d,  trans: %c   \
ldb: %d nrhs: %d  ldaf: %d  ldx: %d \n", matrix_layout, n, trans,  
                                          ldb, nrhs);
#endif

    hModule = NULL;
    dModule = NULL;
    diff_berr = 0;
    diff_ferr = 0;
    diff_xerr = 0;

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
       gtrfs_free();
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

gtrfs_float_parameters:: ~gtrfs_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gtrfs_float_parameters object: destructor invoked. \n");
#endif
   gtrfs_free();
}

//  Test fixture class definition
class sgtrfs_test  : public  ::testing::Test {
public:
   gtrfs_float_parameters  *sgtrfs_obj;
   void SetUp();  
   void TearDown () { delete sgtrfs_obj; }
};


void sgtrfs_test::SetUp(){

    /* LAPACKE SGTSVX prototype */
    typedef int (*Fptr_NL_LAPACKE_sgtrfs) (int matrix_layout, 
            char trans, lapack_int n, lapack_int nrhs, const float* dl,
            const float* d, const float* du, float* dlf, float* df,
            float* duf, float* du2, lapack_int* ipiv, const float* b,
            lapack_int ldb, float* x, lapack_int ldx, 
            float* ferr, float* berr);

    Fptr_NL_LAPACKE_sgtrfs SGTSVX;


     /* LAPACKE DGETRF prototype */
    typedef int (*Fptr_NL_LAPACKE_sgttrf) ( lapack_int n , float *dl,
            float *d , float *du , float *du2 , lapack_int *ipiv );

    Fptr_NL_LAPACKE_sgttrf SGTTRF;


    sgtrfs_obj = new gtrfs_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].equed,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    sgtrfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgtrfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgtrfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgtrfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    SGTSVX = (Fptr_NL_LAPACKE_sgtrfs)dlsym(sgtrfs_obj->hModule, "LAPACKE_sgtrfs");
    ASSERT_TRUE(SGTSVX != NULL) << "failed to get the Netlib LAPACKE_sgtrfs symbol";
    
    /* if  'fact' == 'F', then the i/p matrix A should be in 
       factored form.
       so, invoke the getrf API to compute the factorized A.  */
        SGTTRF = (Fptr_NL_LAPACKE_sgttrf)dlsym(sgtrfs_obj->hModule,"LAPACKE_sgttrf");
        ASSERT_TRUE(SGTTRF != NULL) << "failed to get the Netlib LAPACKE_sgttrf symbol";
            
        sgtrfs_obj->inforef = SGTTRF( sgtrfs_obj->n,
                                      sgtrfs_obj->dlref, sgtrfs_obj->dref,
                                      sgtrfs_obj->duref,
                                      sgtrfs_obj->du2ref,
                                      sgtrfs_obj->ipivref);
                               
        sgtrfs_obj->info = LAPACKE_sgttrf( sgtrfs_obj->n,
                                      sgtrfs_obj->dl, sgtrfs_obj->d,
                                      sgtrfs_obj->du,
                                      sgtrfs_obj->du2,
                                      sgtrfs_obj->ipiv);

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    sgtrfs_obj->inforef = SGTSVX( sgtrfs_obj->matrix_layout,
                                  sgtrfs_obj->trans, sgtrfs_obj->n,
                                  sgtrfs_obj->nrhs,
                                  sgtrfs_obj->dlref, sgtrfs_obj->dref, 
                                  sgtrfs_obj->duref, sgtrfs_obj->dlfref,
                                  sgtrfs_obj->dfref, sgtrfs_obj->dufref,
                                  sgtrfs_obj->du2ref, sgtrfs_obj->ipivref,
                                  sgtrfs_obj->bref, sgtrfs_obj->ldb,
                                  sgtrfs_obj->xref, sgtrfs_obj->ldx,
                                  sgtrfs_obj->ferrref,
                                  sgtrfs_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    sgtrfs_obj->info = LAPACKE_sgtrfs( sgtrfs_obj->matrix_layout,
                                  sgtrfs_obj->trans, sgtrfs_obj->n,
                                  sgtrfs_obj->nrhs,
                                  sgtrfs_obj->dl, sgtrfs_obj->d, 
                                  sgtrfs_obj->du, sgtrfs_obj->dlf,
                                  sgtrfs_obj->df, sgtrfs_obj->duf,
                                  sgtrfs_obj->du2, sgtrfs_obj->ipiv,
                                  sgtrfs_obj->b, sgtrfs_obj->ldb,
                                  sgtrfs_obj->x, sgtrfs_obj->ldx,
                                  sgtrfs_obj->ferr,
                                  sgtrfs_obj->berr);

    if( sgtrfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_sgtrfs is wrong\n", sgtrfs_obj->info );
    }
    if( sgtrfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sgtrfs is wrong\n", 
        sgtrfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
       sgtrfs_obj->diff_xerr =  computeDiff_s( sgtrfs_obj->x_bufsize, 
                 sgtrfs_obj->x, sgtrfs_obj->xref );

    sgtrfs_obj->diff_berr =  computeDiff_s( sgtrfs_obj->nrhs, 
                sgtrfs_obj->berr, sgtrfs_obj->berrref );
                
    sgtrfs_obj->diff_ferr =  computeDiff_s( sgtrfs_obj->nrhs, 
                sgtrfs_obj->ferr, sgtrfs_obj->ferrref );
}

TEST_F(sgtrfs_test, sgtrfs1) {
    EXPECT_NEAR(0.0, sgtrfs_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgtrfs_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgtrfs_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgtrfs_test, sgtrfs2) {
    EXPECT_NEAR(0.0, sgtrfs_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgtrfs_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgtrfs_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgtrfs_test, sgtrfs3) {
    EXPECT_NEAR(0.0, sgtrfs_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgtrfs_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgtrfs_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgtrfs_test, sgtrfs4) {
    EXPECT_NEAR(0.0, sgtrfs_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgtrfs_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgtrfs_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

/* Begin gtrfs_scomplex_parameters  class definition */
class gtrfs_scomplex_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      // capture difference between ref o/p & libflame lapacke o/p.
      float diff_berr, diff_ferr, diff_xerr;
      void *hModule, *dModule;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
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
      lapack_complex_float* x, *xref;
      float* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      float* berr, *berrref; // component-wise relative backward error for each solution vector xj.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      gtrfs_scomplex_parameters ( int matrix_layout_i, char trans_i,
                               char equed_i, lapack_int n_i, lapack_int nrhs_i);
              
      ~gtrfs_scomplex_parameters (); 
};  /* end of gtrfs_scomplex_parameters  class definition */


/* Constructor gtrfs_scomplex_parameters definition */
gtrfs_scomplex_parameters:: gtrfs_scomplex_parameters ( int matrix_layout_i, 
                char trans_i, char equed_i, lapack_int n_i,
                  lapack_int nrhs_i) {
                                    
    matrix_layout = matrix_layout_i;
    trans = trans_i;
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
   printf(" \n gtrfs Double: matrix_layout: %d  n: %d,  trans: %c   \
ldb: %d nrhs: %d  ldaf: %d  ldx: %d \n", matrix_layout, n, trans,  
                                          ldb, nrhs);
#endif

    hModule = NULL;
    dModule = NULL;
    diff_berr = 0;
    diff_ferr = 0;
    diff_xerr = 0;

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
       gtrfs_free();
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

gtrfs_scomplex_parameters:: ~gtrfs_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gtrfs_scomplex_parameters object: destructor invoked. \n");
#endif
   gtrfs_free();
}

//  Test fixture class definition
class cgtrfs_test  : public  ::testing::Test {
public:
   gtrfs_scomplex_parameters  *cgtrfs_obj;
   void SetUp();  
   void TearDown () { delete cgtrfs_obj; }
};


void cgtrfs_test::SetUp(){

    /* LAPACKE CGTSVX prototype */
    typedef int (*Fptr_NL_LAPACKE_cgtrfs) (int matrix_layout, 
            char trans, lapack_int n, lapack_int nrhs, const lapack_complex_float* dl,
            const lapack_complex_float* d, const lapack_complex_float* du, lapack_complex_float* dlf, lapack_complex_float* df,
            lapack_complex_float* duf, lapack_complex_float* du2, lapack_int* ipiv, const lapack_complex_float* b,
            lapack_int ldb, lapack_complex_float* x, lapack_int ldx, 
            float* ferr, float* berr);

    Fptr_NL_LAPACKE_cgtrfs CGTSVX;


     /* LAPACKE DGETRF prototype */
    typedef int (*Fptr_NL_LAPACKE_cgttrf) ( lapack_int n , lapack_complex_float *dl,
            lapack_complex_float *d , lapack_complex_float *du , lapack_complex_float *du2 , lapack_int *ipiv );

    Fptr_NL_LAPACKE_cgttrf CGTTRF;


    cgtrfs_obj = new gtrfs_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].equed,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    cgtrfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgtrfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgtrfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgtrfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    CGTSVX = (Fptr_NL_LAPACKE_cgtrfs)dlsym(cgtrfs_obj->hModule, "LAPACKE_cgtrfs");
    ASSERT_TRUE(CGTSVX != NULL) << "failed to get the Netlib LAPACKE_cgtrfs symbol";
    
    /* if  'fact' == 'F', then the i/p matrix A should be in 
       factored form.
       so, invoke the getrf API to compute the factorized A.  */
        CGTTRF = (Fptr_NL_LAPACKE_cgttrf)dlsym(cgtrfs_obj->hModule,"LAPACKE_cgttrf");
        ASSERT_TRUE(CGTTRF != NULL) << "failed to get the Netlib LAPACKE_cgttrf symbol";
            
        cgtrfs_obj->inforef = CGTTRF( cgtrfs_obj->n,
                                      cgtrfs_obj->dlref, cgtrfs_obj->dref,
                                      cgtrfs_obj->duref,
                                      cgtrfs_obj->du2ref,
                                      cgtrfs_obj->ipivref);
                               
        cgtrfs_obj->info = LAPACKE_cgttrf( cgtrfs_obj->n,
                                      cgtrfs_obj->dl, cgtrfs_obj->d,
                                      cgtrfs_obj->du,
                                      cgtrfs_obj->du2,
                                      cgtrfs_obj->ipiv);

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    cgtrfs_obj->inforef = CGTSVX( cgtrfs_obj->matrix_layout,
                                  cgtrfs_obj->trans, cgtrfs_obj->n,
                                  cgtrfs_obj->nrhs,
                                  cgtrfs_obj->dlref, cgtrfs_obj->dref, 
                                  cgtrfs_obj->duref, cgtrfs_obj->dlfref,
                                  cgtrfs_obj->dfref, cgtrfs_obj->dufref,
                                  cgtrfs_obj->du2ref, cgtrfs_obj->ipivref,
                                  cgtrfs_obj->bref, cgtrfs_obj->ldb,
                                  cgtrfs_obj->xref, cgtrfs_obj->ldx,
                                  cgtrfs_obj->ferrref,
                                  cgtrfs_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    cgtrfs_obj->info = LAPACKE_cgtrfs( cgtrfs_obj->matrix_layout,
                                  cgtrfs_obj->trans, cgtrfs_obj->n,
                                  cgtrfs_obj->nrhs,
                                  cgtrfs_obj->dl, cgtrfs_obj->d, 
                                  cgtrfs_obj->du, cgtrfs_obj->dlf,
                                  cgtrfs_obj->df, cgtrfs_obj->duf,
                                  cgtrfs_obj->du2, cgtrfs_obj->ipiv,
                                  cgtrfs_obj->b, cgtrfs_obj->ldb,
                                  cgtrfs_obj->x, cgtrfs_obj->ldx,
                                  cgtrfs_obj->ferr,
                                  cgtrfs_obj->berr);

    if( cgtrfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cgtrfs is wrong\n", cgtrfs_obj->info );
    }
    if( cgtrfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cgtrfs is wrong\n", 
        cgtrfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
       cgtrfs_obj->diff_xerr =  computeDiff_c( cgtrfs_obj->x_bufsize, 
                 cgtrfs_obj->x, cgtrfs_obj->xref );

    cgtrfs_obj->diff_berr =  computeDiff_s( cgtrfs_obj->nrhs, 
                cgtrfs_obj->berr, cgtrfs_obj->berrref );
                
    cgtrfs_obj->diff_ferr =  computeDiff_s( cgtrfs_obj->nrhs, 
                cgtrfs_obj->ferr, cgtrfs_obj->ferrref );
}

TEST_F(cgtrfs_test, cgtrfs1) {
    EXPECT_NEAR(0.0, cgtrfs_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgtrfs_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgtrfs_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgtrfs_test, cgtrfs2) {
    EXPECT_NEAR(0.0, cgtrfs_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgtrfs_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgtrfs_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgtrfs_test, cgtrfs3) {
    EXPECT_NEAR(0.0, cgtrfs_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgtrfs_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgtrfs_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgtrfs_test, cgtrfs4) {
    EXPECT_NEAR(0.0, cgtrfs_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgtrfs_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cgtrfs_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

/* Begin gtrfs_dcomplex_parameters  class definition */
class gtrfs_dcomplex_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      // capture difference between ref o/p & libflame lapacke o/p.
      double diff_berr, diff_ferr, diff_xerr;
      void *hModule, *dModule;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
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
      lapack_complex_double* x, *xref;
      double* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      double* berr, *berrref; // component-wise relative backward error for each solution vector xj.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      gtrfs_dcomplex_parameters ( int matrix_layout_i, char trans_i,
                               char equed_i, lapack_int n_i, lapack_int nrhs_i);
              
      ~gtrfs_dcomplex_parameters (); 
};  /* end of gtrfs_dcomplex_parameters  class definition */


/* Constructor gtrfs_dcomplex_parameters definition */
gtrfs_dcomplex_parameters:: gtrfs_dcomplex_parameters ( int matrix_layout_i, 
                char trans_i, char equed_i, lapack_int n_i,
                  lapack_int nrhs_i) {
                                    
    matrix_layout = matrix_layout_i;
    trans = trans_i;
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
   printf(" \n gtrfs Double: matrix_layout: %d  n: %d,  trans: %c   \
ldb: %d nrhs: %d  ldaf: %d  ldx: %d \n", matrix_layout, n, trans,  
                                          ldb, nrhs);
#endif

    hModule = NULL;
    dModule = NULL;
    diff_berr = 0;
    diff_ferr = 0;
    diff_xerr = 0;

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
       gtrfs_free();
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

gtrfs_dcomplex_parameters:: ~gtrfs_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gtrfs_dcomplex_parameters object: destructor invoked. \n");
#endif
   gtrfs_free();
}

//  Test fixture class definition
class zgtrfs_test  : public  ::testing::Test {
public:
   gtrfs_dcomplex_parameters  *zgtrfs_obj;
   void SetUp();  
   void TearDown () { delete zgtrfs_obj; }
};


void zgtrfs_test::SetUp(){

    /* LAPACKE ZGTSVX prototype */
    typedef int (*Fptr_NL_LAPACKE_zgtrfs) (int matrix_layout, 
            char trans, lapack_int n, lapack_int nrhs, const lapack_complex_double* dl,
            const lapack_complex_double* d, const lapack_complex_double* du, lapack_complex_double* dlf, lapack_complex_double* df,
            lapack_complex_double* duf, lapack_complex_double* du2, lapack_int* ipiv, const lapack_complex_double* b,
            lapack_int ldb, lapack_complex_double* x, lapack_int ldx, 
            double* ferr, double* berr);

    Fptr_NL_LAPACKE_zgtrfs ZGTSVX;


     /* LAPACKE DGETRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zgttrf) ( lapack_int n , lapack_complex_double *dl,
            lapack_complex_double *d , lapack_complex_double *du , lapack_complex_double *du2 , lapack_int *ipiv );

    Fptr_NL_LAPACKE_zgttrf ZGTTRF;


    zgtrfs_obj = new gtrfs_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].equed,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    zgtrfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgtrfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgtrfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgtrfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    ZGTSVX = (Fptr_NL_LAPACKE_zgtrfs)dlsym(zgtrfs_obj->hModule, "LAPACKE_zgtrfs");
    ASSERT_TRUE(ZGTSVX != NULL) << "failed to get the Netlib LAPACKE_zgtrfs symbol";
    
    /* if  'fact' == 'F', then the i/p matrix A should be in 
       factored form.
       so, invoke the getrf API to compute the factorized A.  */
        ZGTTRF = (Fptr_NL_LAPACKE_zgttrf)dlsym(zgtrfs_obj->hModule,"LAPACKE_zgttrf");
        ASSERT_TRUE(ZGTTRF != NULL) << "failed to get the Netlib LAPACKE_zgttrf symbol";
            
        zgtrfs_obj->inforef = ZGTTRF( zgtrfs_obj->n,
                                      zgtrfs_obj->dlref, zgtrfs_obj->dref,
                                      zgtrfs_obj->duref,
                                      zgtrfs_obj->du2ref,
                                      zgtrfs_obj->ipivref);
                               
        zgtrfs_obj->info = LAPACKE_zgttrf( zgtrfs_obj->n,
                                      zgtrfs_obj->dl, zgtrfs_obj->d,
                                      zgtrfs_obj->du,
                                      zgtrfs_obj->du2,
                                      zgtrfs_obj->ipiv);

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    zgtrfs_obj->inforef = ZGTSVX( zgtrfs_obj->matrix_layout, 
                                  zgtrfs_obj->trans, zgtrfs_obj->n,
                                  zgtrfs_obj->nrhs,
                                  zgtrfs_obj->dlref, zgtrfs_obj->dref, 
                                  zgtrfs_obj->duref, zgtrfs_obj->dlfref,
                                  zgtrfs_obj->dfref, zgtrfs_obj->dufref,
                                  zgtrfs_obj->du2ref, zgtrfs_obj->ipivref,
                                  zgtrfs_obj->bref, zgtrfs_obj->ldb,
                                  zgtrfs_obj->xref, zgtrfs_obj->ldx,
                                  zgtrfs_obj->ferrref,
                                  zgtrfs_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    zgtrfs_obj->info = LAPACKE_zgtrfs( zgtrfs_obj->matrix_layout, 
                                  zgtrfs_obj->trans, zgtrfs_obj->n,
                                  zgtrfs_obj->nrhs,
                                  zgtrfs_obj->dl, zgtrfs_obj->d, 
                                  zgtrfs_obj->du, zgtrfs_obj->dlf,
                                  zgtrfs_obj->df, zgtrfs_obj->duf,
                                  zgtrfs_obj->du2, zgtrfs_obj->ipiv,
                                  zgtrfs_obj->b, zgtrfs_obj->ldb,
                                  zgtrfs_obj->x, zgtrfs_obj->ldx,
                                  zgtrfs_obj->ferr,
                                  zgtrfs_obj->berr);

    if( zgtrfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zgtrfs is wrong\n", zgtrfs_obj->info );
    }
    if( zgtrfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zgtrfs is wrong\n", 
        zgtrfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
       zgtrfs_obj->diff_xerr =  computeDiff_z( zgtrfs_obj->x_bufsize, 
                 zgtrfs_obj->x, zgtrfs_obj->xref );

    zgtrfs_obj->diff_berr =  computeDiff_d( zgtrfs_obj->nrhs, 
                zgtrfs_obj->berr, zgtrfs_obj->berrref );
                
    zgtrfs_obj->diff_ferr =  computeDiff_d( zgtrfs_obj->nrhs, 
                zgtrfs_obj->ferr, zgtrfs_obj->ferrref );
}

TEST_F(zgtrfs_test, zgtrfs1) {
    EXPECT_NEAR(0.0, zgtrfs_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgtrfs_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgtrfs_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgtrfs_test, zgtrfs2) {
    EXPECT_NEAR(0.0, zgtrfs_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgtrfs_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgtrfs_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgtrfs_test, zgtrfs3) {
    EXPECT_NEAR(0.0, zgtrfs_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgtrfs_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgtrfs_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgtrfs_test, zgtrfs4) {
    EXPECT_NEAR(0.0, zgtrfs_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgtrfs_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zgtrfs_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}
