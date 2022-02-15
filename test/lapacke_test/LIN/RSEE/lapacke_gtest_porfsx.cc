#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"

#define porfsx_free() \
  if (a != NULL)    free (a  ); \
  if (aref != NULL) free (aref); \
  if (b != NULL)    free (b   ); \
  if (bref != NULL) free (bref); \
  if (s != NULL)    free (s   ); \
  if (sref != NULL) free (sref); \
  if (x != NULL)    free (x  ); \
  if (xref != NULL) free (xref); \
  if (af != NULL)    free (af  ); \
  if (afref != NULL) free (afref); \
  if (err_bnds_norm != NULL)    free (err_bnds_norm  ); \
  if (err_bnds_normref != NULL) free (err_bnds_normref); \
  if (err_bnds_comp != NULL)    free (err_bnds_comp  ); \
  if (err_bnds_compref != NULL) free (err_bnds_compref); \
  if (ferr != NULL)    free (ferr  ); \
  if (ferrref != NULL) free (ferrref); \
  if (berr != NULL)    free (berr  ); \
  if (berrref != NULL) free (berrref); \
  if( hModule != NULL) dlclose(hModule); \
  if(dModule != NULL) dlclose(dModule)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin porfsx_float_parameters  class definition */
class porfsx_float_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      float diff_berr, diff_xerr;
      float diff_err_bnds_norm, diff_err_bnds_comp;

      void *hModule, *dModule;
      float threshold;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      char equed; // Must be 'N' or 'Y'.
      lapack_int n_err_bnds; // Number of error bounds to return for each RHS
      lapack_int nparams;

      lapack_int n; // The order of A (Square matrix).
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_int ldaf;  //  leading dimension of 'af'
      lapack_int ldx; //  leading dimension of 'x'
      float *s, *sref; // scale factors for A.

      /* Input/ Output parameters */
      float *b, *bref; //right-hand sides for the systems of equations.
      float *a, *aref; //The array ab contains the matrix A
      float *af, *afref; //contains the factored form of the matrix A

      /* Output parameters */
      float* x, *xref;
      float* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      float* berr, *berrref; // component-wise relative backward error for each solution vector xj.
      float rcond, rcondref;
      float* err_bnds_norm, *err_bnds_normref;
      float* err_bnds_comp, *err_bnds_compref;
      float params[3], *paramsref[3]; 
      /* Return Values */
      lapack_int info, inforef;

   public: 
      porfsx_float_parameters ( int matrix_layout_i, char uplo_i,
                               lapack_int n_i, lapack_int nrhs_i,
                               char equed_i, lapack_int n_err_bnds_i,
                               lapack_int nparams_i);
              
      ~porfsx_float_parameters (); 
};  /* end of porfsx_float_parameters  class definition */


/* Constructor porfsx_float_parameters definition */
porfsx_float_parameters:: porfsx_float_parameters ( int matrix_layout_i, 
                         char uplo_i, lapack_int n_i, lapack_int nrhs_i,
          char equed_i, lapack_int n_err_bnds_i, lapack_int nparams_i) {
                                    
    matrix_layout = matrix_layout_i;
    uplo = uplo_i;
    n = n_i;
    nrhs = nrhs_i;
    equed = equed_i;
    n_err_bnds = n_err_bnds_i;
    nparams = nparams_i;

    lda = n;
    ldaf = n;
    
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

    diff_err_bnds_norm = 0;
    diff_err_bnds_comp = 0;
    diff_berr = 0;
    diff_xerr = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n porfsx float:  matrix_layout_i: %d n: %d, uplo: %c  lda: %d  \
ldb: %d nrhs: %d  threshold: %f \n", matrix_layout_i, n, uplo, lda,
ldb, nrhs, lin_solver_paramslist[idx].solver_threhold );
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_float_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_float_buffer_pair( &x, &xref, x_bufsize);
    lapacke_gtest_alloc_float_buffer_pair( &af, &afref, (n*n));
    lapacke_gtest_alloc_float_buffer_pair( &ferr, &ferrref, nrhs);
    lapacke_gtest_alloc_float_buffer_pair( &berr, &berrref, nrhs);
    lapacke_gtest_alloc_float_buffer_pair( &s, &sref, n);
    lapacke_gtest_alloc_float_buffer_pair( &err_bnds_norm, &err_bnds_normref, nrhs*n_err_bnds);
    lapacke_gtest_alloc_float_buffer_pair( &err_bnds_comp, &err_bnds_compref, nrhs*n_err_bnds);
    
    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (s==NULL) || (sref==NULL) || \
        (x==NULL) || (xref==NULL) || \
        (af==NULL) || (afref==NULL) || \
        (ferr==NULL) || (ferrref==NULL) || \
        (err_bnds_norm==NULL) || (err_bnds_normref==NULL) || \
        (err_bnds_comp==NULL) || (err_bnds_compref==NULL) || \
        (berr==NULL) || (berrref==NULL) ){
       EXPECT_FALSE( true) << "pbtrs_float_parameters object: malloc error.";
       porfsx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand_custom_matrix ( a, aref, n, n, uplo);
    memcpy(af, a, (n*n*sizeof(float)));
    memcpy(afref, a, (n*n*sizeof(float)));
    
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, b_bufsize);
    memcpy(x, b, (b_bufsize*sizeof(float)));
    memcpy(xref, b, (b_bufsize*sizeof(float)));
    
    lapacke_gtest_init_float_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(err_bnds_norm, err_bnds_normref, nrhs*n_err_bnds, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(err_bnds_comp, err_bnds_compref, nrhs*n_err_bnds, 0.0);
    
    // TODO: This is a kind of hard coding the array. Need to modify for generalized version.
    lapacke_gtest_init_float_buffer_pair_with_constant(s, sref, n, 2.0);

    
   } /* end of Constructor  */

porfsx_float_parameters:: ~porfsx_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" porfsx_float_parameters object: destructor invoked. \n");
#endif
   porfsx_free();
}

//  Test fixture class definition
class sporfsx_test  : public  ::testing::Test {
public:
   porfsx_float_parameters  *sporfsx_obj;
   void SetUp();  
   void TearDown () { delete sporfsx_obj; }
};


void sporfsx_test::SetUp(){

    /* LAPACKE SPORFSX prototype */
    typedef int (*Fptr_NL_LAPACKE_sporfsx) (int matrix_layout, char uplo,
        char equed, lapack_int n, lapack_int nrhs, const float* a,
        lapack_int lda, const float* af, lapack_int ldaf, const float* s,
        const float* b, lapack_int ldb, float* x, lapack_int ldx,
        float* rcond, float* berr, lapack_int n_err_bnds,
        float* err_bnds_norm, float* err_bnds_comp, lapack_int nparams,
        float* params);

    Fptr_NL_LAPACKE_sporfsx SPORFSX;

     /* LAPACKE SPOTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_spotrf) ( int matrix_layout , char uplo ,
                                lapack_int n , float * a , lapack_int lda );

    Fptr_NL_LAPACKE_spotrf SPOTRF;

    /* LAPACKE SPOTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_spotrs) (int matrix_layout, char uplo,
                        lapack_int n, lapack_int nrhs, const float * a,
                          lapack_int lda, float * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_spotrs SPOTRS;

    sporfsx_obj = new porfsx_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs,
                           lin_solver_paramslist[idx].equed_porfsx,
                           lin_solver_paramslist[idx].n_err_bnds_porfsx,
                           lin_solver_paramslist[idx].nparams_porfsx
                           );

    sporfsx_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
    idx = Circular_Increment_Index(idx);

    sporfsx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sporfsx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sporfsx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sporfsx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    SPORFSX = (Fptr_NL_LAPACKE_sporfsx)dlsym(sporfsx_obj->hModule, "LAPACKE_sporfsx");
    ASSERT_TRUE(SPORFSX != NULL) << "failed to pot the Netlib LAPACKE_sporfsx symbol";

    SPOTRS = (Fptr_NL_LAPACKE_spotrs)dlsym(sporfsx_obj->hModule, "LAPACKE_spotrs");
    ASSERT_TRUE(SPOTRS != NULL) << "failed to get the Netlib LAPACKE_spotrs symbol";
    
    SPOTRF = (Fptr_NL_LAPACKE_spotrf)dlsym(sporfsx_obj->hModule,"LAPACKE_spotrf");
    ASSERT_TRUE(SPOTRF != NULL) << "failed to pot the Netlib LAPACKE_spotrf symbol";
    /*  Generate the i/ps for porfsx by calling the sequence of potrf, potrs APIs  */
            
    sporfsx_obj->inforef = SPOTRF( sporfsx_obj->matrix_layout,
                                  sporfsx_obj->uplo, sporfsx_obj->n,
                                  sporfsx_obj->afref,
                                  sporfsx_obj->lda);
                                  
    sporfsx_obj->inforef = SPOTRS( sporfsx_obj->matrix_layout,
                                  sporfsx_obj->uplo, 
                                  sporfsx_obj->n,
                                  sporfsx_obj->nrhs,
                                  (const float *)sporfsx_obj->afref, 
                                  sporfsx_obj->lda,
                                  sporfsx_obj->xref,
                                  sporfsx_obj->ldb);
                          
    // invoking netlib potrf for generation of libflame i/p as well.                       
    sporfsx_obj->info = SPOTRF( sporfsx_obj->matrix_layout,
                               sporfsx_obj->uplo, sporfsx_obj->n,
                               sporfsx_obj->af,
                               sporfsx_obj->lda);

    sporfsx_obj->info = SPOTRS( sporfsx_obj->matrix_layout,
                               sporfsx_obj->uplo, sporfsx_obj->n,
                               sporfsx_obj->nrhs,
                               (const float *)sporfsx_obj->af, 
                                sporfsx_obj->lda,
                                sporfsx_obj->x,
                                sporfsx_obj->ldb );

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    sporfsx_obj->inforef = SPORFSX( sporfsx_obj->matrix_layout,
                                  sporfsx_obj->uplo,
                                  sporfsx_obj->equed,
                                  sporfsx_obj->n,
                                  sporfsx_obj->nrhs,
                                  sporfsx_obj->aref, sporfsx_obj->lda, 
                                  sporfsx_obj->afref, sporfsx_obj->ldaf,
                                  sporfsx_obj->sref,
                                  sporfsx_obj->bref, sporfsx_obj->ldb,
                                  sporfsx_obj->xref, sporfsx_obj->ldx,
                                  &sporfsx_obj->rcondref,
                                  sporfsx_obj->berrref,
                                  sporfsx_obj->n_err_bnds,
                                  sporfsx_obj->err_bnds_normref,
                                  sporfsx_obj->err_bnds_compref,
                                  sporfsx_obj->nparams,
                                  sporfsx_obj->params
                                  );

    /* Compute libflame's Lapacke o/p  */
    //sporfsx_obj->info = SPORFSX( sporfsx_obj->matrix_layout,
    sporfsx_obj->info = LAPACKE_sporfsx( sporfsx_obj->matrix_layout,
                                  sporfsx_obj->uplo,
                                  sporfsx_obj->equed,
                                  sporfsx_obj->n,
                                  sporfsx_obj->nrhs,
                                  sporfsx_obj->a, sporfsx_obj->lda, 
                                  sporfsx_obj->af, sporfsx_obj->ldaf,
                                  sporfsx_obj->s,
                                  sporfsx_obj->b, sporfsx_obj->ldb,
                                  sporfsx_obj->x, sporfsx_obj->ldx,
                                  &sporfsx_obj->rcond,
                                  sporfsx_obj->berr,
                                  sporfsx_obj->n_err_bnds,
                                  sporfsx_obj->err_bnds_norm,
                                  sporfsx_obj->err_bnds_comp,
                                  sporfsx_obj->nparams,
                                  sporfsx_obj->params
                                  );

    if( sporfsx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_sporfsx is wrong\n", sporfsx_obj->info );
    }
    if( sporfsx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sporfsx is wrong\n", 
        sporfsx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    sporfsx_obj->diff_xerr =  computeDiff_s( sporfsx_obj->x_bufsize, 
                sporfsx_obj->x, sporfsx_obj->xref );

    sporfsx_obj->diff_berr =  computeDiff_s( sporfsx_obj->nrhs, 
                sporfsx_obj->berr, sporfsx_obj->berrref );

    sporfsx_obj->diff_err_bnds_norm =  computeDiff_s( sporfsx_obj->nrhs*sporfsx_obj->n_err_bnds, 
                sporfsx_obj->err_bnds_norm, sporfsx_obj->err_bnds_normref );

    sporfsx_obj->diff_err_bnds_comp =  computeDiff_s( sporfsx_obj->nrhs*sporfsx_obj->n_err_bnds, 
                sporfsx_obj->err_bnds_comp, sporfsx_obj->err_bnds_compref );

}

TEST_F(sporfsx_test, sporfsx1) {
    EXPECT_NEAR(0.0, sporfsx_obj->diff_xerr, sporfsx_obj->threshold);
    EXPECT_NEAR(0.0, sporfsx_obj->diff_berr, sporfsx_obj->threshold);
    EXPECT_NEAR(0.0, sporfsx_obj->diff_err_bnds_norm, sporfsx_obj->threshold);
    EXPECT_NEAR(0.0, sporfsx_obj->diff_err_bnds_comp, sporfsx_obj->threshold);
}

TEST_F(sporfsx_test, sporfsx2) {
    EXPECT_NEAR(0.0, sporfsx_obj->diff_xerr, sporfsx_obj->threshold);
    EXPECT_NEAR(0.0, sporfsx_obj->diff_berr, sporfsx_obj->threshold);
    EXPECT_NEAR(0.0, sporfsx_obj->diff_err_bnds_norm, sporfsx_obj->threshold);
    EXPECT_NEAR(0.0, sporfsx_obj->diff_err_bnds_comp, sporfsx_obj->threshold);
}

TEST_F(sporfsx_test, sporfsx3) {
    EXPECT_NEAR(0.0, sporfsx_obj->diff_xerr, sporfsx_obj->threshold);
    EXPECT_NEAR(0.0, sporfsx_obj->diff_berr, sporfsx_obj->threshold);
    EXPECT_NEAR(0.0, sporfsx_obj->diff_err_bnds_norm, sporfsx_obj->threshold);
    EXPECT_NEAR(0.0, sporfsx_obj->diff_err_bnds_comp, sporfsx_obj->threshold);
}

TEST_F(sporfsx_test, sporfsx4) {
    EXPECT_NEAR(0.0, sporfsx_obj->diff_xerr, sporfsx_obj->threshold);
    EXPECT_NEAR(0.0, sporfsx_obj->diff_berr, sporfsx_obj->threshold);
    EXPECT_NEAR(0.0, sporfsx_obj->diff_err_bnds_norm, sporfsx_obj->threshold);
    EXPECT_NEAR(0.0, sporfsx_obj->diff_err_bnds_comp, sporfsx_obj->threshold);
}