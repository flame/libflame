#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define gesvxx_free() \
  if (a != NULL)    free (a  ); \
  if (aref != NULL) free (aref); \
  if (b != NULL)    free (b   ); \
  if (bref != NULL) free (bref); \
  if (x != NULL)    free (x  ); \
  if (xref != NULL) free (xref); \
  if (af != NULL)    free (af  ); \
  if (afref != NULL) free (afref); \
  if (r != NULL)    free (r  ); \
  if (rref != NULL) free (rref); \
  if (c != NULL)    free (c  ); \
  if (cref != NULL) free (cref); \
  if (err_bnds_norm != NULL) free (err_bnds_norm); \
  if (err_bnds_normref != NULL) free (err_bnds_normref); \
  if (err_bnds_comp != NULL) free (err_bnds_comp); \
  if (err_bnds_compref != NULL) free (err_bnds_compref); \
  if (rpvgrw != NULL)    free (rpvgrw  ); \
  if (rpvgrwref != NULL) free (rpvgrwref); \
  if (berr != NULL)    free (berr  ); \
  if (berrref != NULL) free (berrref); \
  if (ipiv != NULL)    free (ipiv  ); \
  if (ipivref != NULL) free (ipivref); \
  if( hModule != NULL) dlclose(hModule); \
  if(dModule != NULL) dlclose(dModule)


// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin gesvxx_float_parameters  class definition */
class gesvxx_float_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      float diff; // capture difference between ref o/p & libflame lapacke o/p.
      float diff_berr, diff_rpvgrw, diff_xerr;
      void *hModule, *dModule;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char fact; //  Must be 'F', 'N', or 'E'.
      char trans; //  Must be 'N' , 'T' or 'C'.

      lapack_int n; // The order of A (Square matrix).
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_int ldaf;  //  leading dimension of 'af'
      lapack_int ldx; //  leading dimension of 'x'
      lapack_int n_err_bnds; // Number of error bounds
      lapack_int nparams; // number of user given parameters. if 0, then default params are used.
      float* params; // user given prameters. 
      
      /* Input/ Output parameters */
      float *b, *bref; //right-hand sides for the systems of equations.
      float *a, *aref; //The array ab contains the matrix A
      float *af, *afref; //contains the factored form of the matrix A
      float* r, *rref; // the row scale factors for A.
      float* c, *cref; // column scale factors for A.
      lapack_int *ipiv, *ipivref; // pivot buffer
      
      /* Output parameters */
      char   equed, equedref; //  Must be 'N', 'R', 'C', or 'B'.
      float rpivot, rpivotref; // reciprocal pivot growth factor.
      float rcond, rcondref; // estimate of the reciprocal condition number of the matrix A.
      float* x, *xref;
      float* rpvgrw, *rpvgrwref; // estimated forward error bound for each solution vector x.
      float* berr, *berrref; // component-wise relative backward error for each solution vector xj.
      float *err_bnds_norm, *err_bnds_normref; //  error bounds & condition numbers corresponding to the normwise relative error
      float *err_bnds_comp, *err_bnds_compref; //  various error bounds and condition numbers corresponding to the componentwise relative error
      
      /* Return Values */
      lapack_int info, inforef;

   public: 
      gesvxx_float_parameters ( int matrix_layout_i, char fact_i, char trans_i,
                               char equed_i, lapack_int n_i, lapack_int nrhs_i,
                               int nparams_i, int nerrbnds_i);
              
      ~gesvxx_float_parameters (); 
};  /* end of gesvxx_float_parameters  class definition */


/* Constructor gesvxx_float_parameters definition */
gesvxx_float_parameters:: gesvxx_float_parameters ( int matrix_layout_i, 
                char fact_i, char trans_i, char equed_i, lapack_int n_i,
                  lapack_int nrhs_i, int nparams_i, int nerrbnds_i) {
                                    
    matrix_layout = matrix_layout_i;
    fact= fact_i;
    trans = trans_i;
    // equed acts as input when " fact = 'F' " else output
    equed = equed_i;
    equedref = equed_i;
    nparams = nparams_i;
    n_err_bnds = nerrbnds_i;
    n = n_i;
    nrhs = nrhs_i;
    lda = n;
    ldaf = n;
    
#if LAPACKE_TEST_VERBOSE
   printf(" \n gesvxx Double:  n: %d, fact: %c trans: %c  lda: %d  \
ldb: %d nrhs: %d  ldaf: %d  ldx: %d \n",  n, fact, trans, lda, 
                                          ldb, nrhs);
#endif
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

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_float_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_float_buffer_pair( &x, &xref, x_bufsize);
    lapacke_gtest_alloc_float_buffer_pair( &af, &afref, (n*n));
    lapacke_gtest_alloc_float_buffer_pair( &r, &rref, n);
    lapacke_gtest_alloc_float_buffer_pair( &c, &cref, n);
    lapacke_gtest_alloc_float_buffer_pair( &rpvgrw, &rpvgrwref, nrhs);
    lapacke_gtest_alloc_float_buffer_pair( &berr, &berrref, nrhs);
    lapacke_gtest_alloc_float_buffer_pair( &err_bnds_norm, &err_bnds_normref, nrhs*n_err_bnds);
    lapacke_gtest_alloc_float_buffer_pair( &err_bnds_comp, &err_bnds_compref, nrhs*n_err_bnds);
    
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (x==NULL) || (xref==NULL) || \
        (af==NULL) || (afref==NULL) || \
        (r==NULL) || (rref==NULL) || \
        (c==NULL) || (cref==NULL) || \
        (rpvgrw==NULL) || (rpvgrwref==NULL) || \
        (berr==NULL) || (berrref==NULL) || \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "pbtrs_float_parameters object: malloc error.";
       gesvxx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, n*n);
    memcpy(af, a, (n*n*sizeof(float)));
    memcpy(afref, a, (n*n*sizeof(float)));
    
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, b_bufsize);
    
    lapacke_gtest_init_float_buffer_pair_with_constant(x, xref, x_bufsize, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(rpvgrw, rpvgrwref, nrhs, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(r, rref, n, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(c, cref, n, 0.0);

    
   } /* end of Constructor  */

gesvxx_float_parameters:: ~gesvxx_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gesvxx_float_parameters object: destructor invoked. \n");
#endif
   gesvxx_free();
}

//  Test fixture class definition
class sgesvxx_test  : public  ::testing::Test {
public:
   gesvxx_float_parameters  *sgesvxx_obj;
   void SetUp();  
   void TearDown () { delete sgesvxx_obj; }
};


void sgesvxx_test::SetUp(){

    /* LAPACKE SGESVXX prototype */
    typedef int (*Fptr_NL_LAPACKE_sgesvxx) (int matrix_layout,
                                            char fact,
                                            char trans,
                                            lapack_int n,
                                            lapack_int nrhs,
                                            float* a,
                                            lapack_int lda,
                                            float* af,
                                            lapack_int ldaf,
                                            lapack_int* ipiv,
                                            char* equed,
                                            float* r,
                                            float* c,
                                            float* b,
                                            lapack_int ldb,
                                            float* x,
                                            lapack_int ldx,
                                            float* rcond,
                                            float* rpvgrw,
                                            float* berr,
                                            lapack_int n_err_bnds,
                                            float* err_bnds_norm,
                                            float* err_bnds_comp,
                                            lapack_int nparams,
                                            const float* params);
    Fptr_NL_LAPACKE_sgesvxx SGESVXX;

     /* LAPACKE SGETRF prototype */
    typedef int (*Fptr_NL_LAPACKE_sgetrf) ( int matrix_layout,lapack_int m,lapack_int n,
                                    float* a,lapack_int lda,lapack_int* ipiv );

    Fptr_NL_LAPACKE_sgetrf SGETRF;


    sgesvxx_obj = new gesvxx_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].fact,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].equed,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs,
                           lin_driver_paramslist[idx].nparams,
                           lin_driver_paramslist[idx].nerrbnds);

    idx = Circular_Increment_Index(idx);

    sgesvxx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgesvxx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgesvxx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgesvxx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    SGESVXX = (Fptr_NL_LAPACKE_sgesvxx)dlsym(sgesvxx_obj->hModule, "LAPACKE_sgesvxx");
    ASSERT_TRUE(SGESVXX != NULL) << "failed to get the Netlib LAPACKE_sgesvxx symbol";
    
    /* if  'fact' == 'F', then the i/p matrix A should be in 
       factored form.
       so, invoke the getrf API to compute the factorized A.  */
    if(sgesvxx_obj->fact == 'F') {
        SGETRF = (Fptr_NL_LAPACKE_sgetrf)dlsym(sgesvxx_obj->hModule,"LAPACKE_sgetrf");
        ASSERT_TRUE(SGETRF != NULL) << "failed to get the Netlib LAPACKE_sgetrf symbol";
            

        sgesvxx_obj->inforef = SGETRF( sgesvxx_obj->matrix_layout,
                                      sgesvxx_obj->n, sgesvxx_obj->n,
                                      sgesvxx_obj->afref,
                                      sgesvxx_obj->lda,
                                      sgesvxx_obj->ipivref);
                               
        sgesvxx_obj->info = LAPACKE_sgetrf( sgesvxx_obj->matrix_layout,
                                           sgesvxx_obj->n, sgesvxx_obj->n,
                                           sgesvxx_obj->af,
                                           sgesvxx_obj->lda,
                                           sgesvxx_obj->ipiv);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */
    sgesvxx_obj->inforef = SGESVXX(     sgesvxx_obj->matrix_layout,
                                    sgesvxx_obj->fact,
                                    sgesvxx_obj->trans,
                                    sgesvxx_obj->n,
                                    sgesvxx_obj->nrhs,
                                    sgesvxx_obj->aref,
                                    sgesvxx_obj->lda, 
                                    sgesvxx_obj->afref,
                                    sgesvxx_obj->ldaf,
                                    sgesvxx_obj->ipivref,
                                    &sgesvxx_obj->equedref,
                                    sgesvxx_obj->rref,
                                    sgesvxx_obj->cref,
                                    sgesvxx_obj->bref,
                                    sgesvxx_obj->ldb,
                                    sgesvxx_obj->xref,
                                    sgesvxx_obj->ldx,
                                    &sgesvxx_obj->rcondref, 
                                    sgesvxx_obj->rpvgrwref,
                                    sgesvxx_obj->berrref,
                                    sgesvxx_obj->n_err_bnds,
                                    sgesvxx_obj->err_bnds_normref,
                                    sgesvxx_obj->err_bnds_compref,
                                    sgesvxx_obj->nparams,
                                    sgesvxx_obj->params );

    /* Compute libflame's Lapacke o/p  */
    sgesvxx_obj->info = LAPACKE_sgesvxx ( sgesvxx_obj->matrix_layout,
//    sgesvxx_obj->info =  SGESVXX( sgesvxx_obj->matrix_layout,
                                    sgesvxx_obj->fact,
                                    sgesvxx_obj->trans,
                                    sgesvxx_obj->n,
                                    sgesvxx_obj->nrhs,
                                    sgesvxx_obj->a,
                                    sgesvxx_obj->lda, 
                                    sgesvxx_obj->af,
                                    sgesvxx_obj->ldaf,
                                    sgesvxx_obj->ipiv,
                                    &sgesvxx_obj->equed,
                                    sgesvxx_obj->r,
                                    sgesvxx_obj->c,
                                    sgesvxx_obj->b,
                                    sgesvxx_obj->ldb,
                                    sgesvxx_obj->x,
                                    sgesvxx_obj->ldx,
                                    &sgesvxx_obj->rcond, 
                                    sgesvxx_obj->rpvgrw,
                                    sgesvxx_obj->berr,
                                    sgesvxx_obj->n_err_bnds,
                                    sgesvxx_obj->err_bnds_norm,
                                    sgesvxx_obj->err_bnds_comp,
                                    sgesvxx_obj->nparams,
                                    sgesvxx_obj->params);
                                    
    if( sgesvxx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_sgesvxx is wrong\n", sgesvxx_obj->info );
    }
    if( sgesvxx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sgesvxx is wrong\n", 
        sgesvxx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    sgesvxx_obj->diff =  computeDiff_s( sgesvxx_obj->b_bufsize, 
                sgesvxx_obj->b, sgesvxx_obj->bref );

    sgesvxx_obj->diff_xerr =  computeDiff_s( sgesvxx_obj->x_bufsize, 
                sgesvxx_obj->x, sgesvxx_obj->xref );

    sgesvxx_obj->diff_berr =  computeDiff_s( sgesvxx_obj->nrhs, 
                sgesvxx_obj->berr, sgesvxx_obj->berrref );
                
    sgesvxx_obj->diff_rpvgrw =  computeDiff_s( sgesvxx_obj->nrhs, 
                sgesvxx_obj->rpvgrw, sgesvxx_obj->rpvgrwref );
}

TEST_F(sgesvxx_test, sgesvxx1) {
    EXPECT_NEAR(0.0, sgesvxx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgesvxx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgesvxx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgesvxx_obj->diff_rpvgrw, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (sgesvxx_obj->rpivot - sgesvxx_obj->rpivotref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgesvxx_test, sgesvxx2) {
    EXPECT_NEAR(0.0, sgesvxx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgesvxx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgesvxx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgesvxx_obj->diff_rpvgrw, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (sgesvxx_obj->rpivot - sgesvxx_obj->rpivotref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgesvxx_test, sgesvxx3) {
    EXPECT_NEAR(0.0, sgesvxx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgesvxx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgesvxx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgesvxx_obj->diff_rpvgrw, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (sgesvxx_obj->rpivot - sgesvxx_obj->rpivotref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgesvxx_test, sgesvxx4) {
    EXPECT_NEAR(0.0, sgesvxx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgesvxx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgesvxx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sgesvxx_obj->diff_rpvgrw, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (sgesvxx_obj->rpivot - sgesvxx_obj->rpivotref),
                LAPACKE_GTEST_THRESHOLD);
}

