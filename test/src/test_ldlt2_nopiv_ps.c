/*
    Copyright (c) 2020 Advanced Micro Devices, Inc.  All rights reserved.
    Sep 24, 2020
*/

#include "FLAME.h"
#include "test_libflame.h"
#include "test_ldlt2_nopiv_ps.h"
#include "test_common.h"

#define THRESH 30
#define THRESH2 300

// static variables
static char* op_str                   = "Partial / Incomplete LDLT(2) factorization without Pivoting";

extern float  slansy_( char *, char *, integer *, float  *, integer *, float  * );
extern double dlansy_( char *, char *, integer *, double *, integer *, double * );
extern float  clansy_( char *, char *, integer *, scomplex *, integer *, float * );
extern double zlansy_( char *, char *, integer *, dcomplex *, integer *, double * );

extern int sspffrt2_( float  *ap, integer *n, integer * ncolm, float  *work, float  *work2 );
extern int dspffrt2_( double *ap, integer *n, integer * ncolm, double *work, double *work2 );
extern int cspffrt2_( scomplex *ap, integer *n, integer * ncolm, scomplex *work, scomplex *work2 );
extern int zspffrt2_( dcomplex *ap, integer *n, integer * ncolm, dcomplex *work, dcomplex *work2 );

extern float slamch_(char *cmach);
extern double dlamch_(char *cmach);

/*************************************************************************/
/*                     SINGLE PRECISION ROUTINES                         */
/*************************************************************************/

float diff_norm_s( float *ad, float *fod, integer n, integer ni )
{
   float *L, *D, *Lt;
   float *ra1, *ra2, *rad;
   float czero = 0.0;
   float cone = 1.0;
   float snrm, anrm;
   float rnrm, eps;

   integer di, i, j;

   char *tf = "No Transpose";
   char *nm = "1";
   char *ul = "lower";

   L   = (float *) malloc( n * n * sizeof( float ) );
   D   = (float *) malloc( n * n * sizeof( float ) );
   Lt  = (float *) malloc( n * n * sizeof( float ) );
   ra1 = (float *) malloc( n * n * sizeof( float ) );
   ra2 = (float *) malloc( n * n * sizeof( float ) );
   rad = (float *) malloc( n * n * sizeof( float ) );

   set_identity_s(L,  n, n, n);
   set_identity_s(D,  n, n, n);
   set_identity_s(Lt, n, n, n);

   /* Form the original matrix from parts */

   /* Unpack L, D and L' */
   di = 0;
   for( i = 0; i < ni; i++ )
   {
      /* Form L */
      for( j = i; j < n; j++ )
      {
         L[i * n + j] = ad[j - i + di];
      }
      L[i * n + i] = 1.0 / ad[di];

      /* Form L' */
      for( j = i; j < n; j++ )
      {
         Lt[j * n + i] = ad[j - i + di];
      }
      Lt[i * n + i] = 1.0 / ad[di];

      D[i * n + i] = ad[di];

      di += n - i;
   }

   for( ; i < n; i++ )
   {
      for( j = i; j < n; j++ )
      {
         D[j + i * n] = ad[di + j - i];
         D[i + j * n] = ad[di + j - i];
      }

      di += n - i;
   }

   /* Product of L and D */
   sgemm_( tf, tf, &n, &n, &n, &cone, L, &n, D, &n, &czero, ra1, &n );

   /* Product of L*D and L' */
   sgemm_( tf, tf, &n, &n, &n, &cone, ra1, &n, Lt, &n, &czero, ra2, &n );

   /* Calculate difference matrix */
   for( i = 0; i < n; i++ )
   {
      for( j = 0; j < n; j++ )
      {
         rad[i * n + j] = ra2[i * n + j] - fod[i * n + j];
      }
   }
   eps = slamch_("Epsilon");
   anrm = slansy_( nm, ul, &n, fod, &n, L );
   rnrm = slansy_( nm, ul, &n, rad, &n, L );

   /* Calculate the final norm */
   snrm = rnrm / n / anrm / eps;

   free( L );
   free( D );
   free( Lt );
   free( ra1 );
   free( ra2 );
   free( rad );

   return snrm;
}

void test_ldlt2_nopiv_ps_s( test_params_t *params )
{
   integer    n, ni;
   float  *od, *ad;
   float  *fod;
   float  snrm;

   int    i, n_repeats;
   double time, perf;
   double time_min;
   char   *strpass = "PASS";
   char   *strmarg = "MARGINAL";
   char   *strfail = "FAILURE";
   char   *pcode;

   for( n = params->p_first; n <= params->p_max; n += params->p_inc )
   {
      integer pn = n * (n + 1) / 2;
      /* allocate packed matrices */
      od = (float *) malloc( pn * sizeof( float ) );
      ad = (float *) malloc( pn * sizeof( float ) );

      /* full sized matrix to store original */
      fod = (float *) malloc( n * n * sizeof( float ) );

      /* Initialize with random values */
      rand_sym_matrix_s( fod, n, n, n );

      /* Pack the symmetric matrix into input */
      pack_matrix_lt_s( fod, od, n, n );
      pack_matrix_lt_s( fod, ad, n, n );

      pcode = strpass;
      snrm = 0.0;
      perf = 0.0;
      if( (integer) params->p_nfact == -2 )
      {
         n_repeats = (params->n_repeats < 5) ? params->n_repeats : 5;
         for( ni = 0; ni <= n; ni++ )
         {
            time_min   = 1e9;
            for( i = 0; i < n_repeats; i++ )
            {
               copy_matrix_s(od, ad, pn, 1, 1, 1);

               time = FLA_Clock();
               sspffrt2_(ad, &n, &ni, NULL, NULL);
               time = FLA_Clock() - time;
               time_min = min( time_min, time );

               /* Calculate error norm */
               snrm = diff_norm_s( ad, fod, n, ni );

               if( snrm >= THRESH2 )
               {
                  pcode = strfail;
                  break;
               }
               else if( snrm >= THRESH )
               {
                  pcode = strmarg;
                  break;
               }
            }
            /* Calculate performance in GFLOPS */
            perf = ni/6.0f*(2.0*ni*ni-6.0*ni*n+3.0*ni+6.0*n*n-6.0*n+7) / time_min / FLOPS_PER_UNIT_PERF;

            libfla_test_print_result_info("SPFFRT2",
                                          (char *) 's',
                                          "c",
                                          n,
                                          perf,
                                          time_min,
                                          snrm,
                                          pcode,
                                          ni);
            if( pcode[0] == 'M' )
               break;
         }
      }
      else
      {
         time_min   = 1e9;
         snrm = 0.0f;
         if( (integer) params->p_nfact == -1 )
         {
            ni = rand() % n;
         }
         else
         {
            ni = params->p_nfact;
         }

         for( i = 0; i < params->n_repeats; i++ )
         {
            copy_matrix_s(od, ad, pn, 1, 1, 1);

            time = FLA_Clock();
            sspffrt2_(ad, &n, &ni, NULL, NULL);
            time = FLA_Clock() - time;
            time_min = min( time_min, time );

            /* Calculate error norm */
            snrm = diff_norm_s( ad, fod, n, ni );

            if( snrm >= THRESH2 )
            {
               pcode = strfail;
               break;
            }
            else if( snrm >= THRESH )
            {
               pcode = strmarg;
               break;
            }
         }

         /* Calculate performance in GFLOPS */
         perf = ni/6.0f*(2.0*ni*ni-6.0*ni*n+3.0*ni+6.0*n*n-6.0*n+7) / time_min / FLOPS_PER_UNIT_PERF;

         libfla_test_print_result_info("SPFFRT2",
                                       (char *) 's',
                                       "c",
                                       n,
                                       perf,
                                       time_min,
                                       snrm,
                                       pcode,
                                       ni );
      }

      free(od);
      free(ad);
      free(fod);
   }
   return;
}

/*************************************************************************/
/*                     DOUBLE PRECISION ROUTINES                         */
/*************************************************************************/

double diff_norm_d( double *ad, double *fod, integer n, integer ni )
{
   double *L, *D, *Lt;
   double *ra1, *ra2, *rad;
   double czero = 0.0;
   double cone = 1.0;
   double dnrm, anrm;
   double rnrm, eps;

   integer di, i, j;

   char *tf = "No Transpose";
   char *nm = "1";
   char *ul = "lower";

   L   = (double *) malloc( n * n * sizeof( double ) );
   D   = (double *) malloc( n * n * sizeof( double ) );
   Lt  = (double *) malloc( n * n * sizeof( double ) );
   ra1 = (double *) malloc( n * n * sizeof( double ) );
   ra2 = (double *) malloc( n * n * sizeof( double ) );
   rad = (double *) malloc( n * n * sizeof( double ) );

   set_identity_d(L,  n, n, n);
   set_identity_d(D,  n, n, n);
   set_identity_d(Lt, n, n, n);

   /* Form the original matrix from parts */

   /* Unpack L, D and L' */
   di = 0;
   for( i = 0; i < ni; i++ )
   {
      /* Form L */
      for( j = i; j < n; j++ )
      {
         L[i * n + j] = ad[j - i + di];
      }
      L[i * n + i] = 1.0 / ad[di];

      /* Form L' */
      for( j = i; j < n; j++ )
      {
         Lt[j * n + i] = ad[j - i + di];
      }
      Lt[i * n + i] = 1.0 / ad[di];

      D[i * n + i] = ad[di];

      di += n - i;
   }

   for( ; i < n; i++ )
   {
      for( j = i; j < n; j++ )
      {
         D[j + i * n] = ad[di + j - i];
         D[i + j * n] = ad[di + j - i];
      }

      di += n - i;
   }

   /* Product of L and D */
   dgemm_( tf, tf, &n, &n, &n, &cone, L, &n, D, &n, &czero, ra1, &n );

   /* Product of L*D and L' */
   dgemm_( tf, tf, &n, &n, &n, &cone, ra1, &n, Lt, &n, &czero, ra2, &n );

   /* Calculate difference matrix */
   for( i = 0; i < n; i++ )
   {
      for( j = 0; j < n; j++ )
      {
         rad[i * n + j] = ra2[i * n + j] - fod[i * n + j];
      }
   }
   eps = dlamch_("Epsilon");
   anrm = dlansy_( nm, ul, &n, fod, &n, L );
   rnrm = dlansy_( nm, ul, &n, rad, &n, L );

   /* Calculate the final norm */
   dnrm = rnrm / n / anrm / eps;

   free( L );
   free( D );
   free( Lt );
   free( ra1 );
   free( ra2 );
   free( rad );

   return dnrm;
}

void test_ldlt2_nopiv_ps_d( test_params_t *params )
{
   integer n, ni;
   double *od, *ad;
   double *work;
   double *fod;
   double dnrm;

   int    i, n_repeats;
   double time, perf;
   double time_min;
   char   *strpass = "PASS";
   char   *strmarg = "MARGINAL";
   char   *strfail = "FAILURE";
   char   *pcode;

   for( n = params->p_first; n <= params->p_max; n += params->p_inc )
   {
      integer pn = n * (n + 1) / 2;

      work = (double *) malloc( 2 * n * sizeof( double ) );

      /* allocate packed matrices */
      od = (double *) malloc( pn * sizeof( double ) );
      ad = (double *) malloc( pn * sizeof( double ) );

      /* full sized matrix to store original */
      fod = (double *) malloc( n * n * sizeof( double ) );

      /* Initialize with random values */
      rand_sym_matrix_d( fod, n, n, n );

      /* Pack the symmetric matrix into input */
      pack_matrix_lt_d( fod, od, n, n );
      pack_matrix_lt_d( fod, ad, n, n );

      pcode = strpass;
      dnrm = 0.0;
      perf = 0.0;
      if( (integer) params->p_nfact == -2 )
      {        
         n_repeats = (params->n_repeats < 5) ? params->n_repeats : 5;
         for( ni = 0; ni <= n; ni++ )
         {
            time_min   = 1e9;
            for( i = 0; i < n_repeats; i++ )
            {
               copy_matrix_d(od, ad, pn, 1, 1, 1);

               time = FLA_Clock();
               dspffrt2_(ad, &n, &ni, work, NULL);
               time = FLA_Clock() - time;
               time_min = min( time_min, time );

               /* Calculate error norm */
               dnrm = diff_norm_d( ad, fod, n, ni );

               if( dnrm >= THRESH2 )
               {
                  pcode = strfail;
                  break;
               }
               else if( dnrm >= THRESH )
               {
                  pcode = strmarg;
                  break;
               }
            }
            /* Calculate performance in GFLOPS */
            perf = ni/6.0*(2.0*ni*ni-6.0*ni*n+3.0*ni+6.0*n*n-6.0*n+7) / time_min / FLOPS_PER_UNIT_PERF;

            libfla_test_print_result_info("SPFFRT2",
                                          (char *) 'd',
                                          "c",
                                          n,
                                          perf,
                                          time_min,
                                          dnrm,
                                          pcode,
                                          ni);
            if( pcode[0] == 'M' )
               break;
         }
      }
      else
      {
         time_min   = 1e9;
         dnrm = 0.0;
         if( (integer) params->p_nfact == -1 )
         {
            ni = rand() % n;
         }
         else
         {
            ni = params->p_nfact;
         }

         for( i = 0; i < params->n_repeats; i++ )
         {
            copy_matrix_d(od, ad, pn, 1, 1, 1);

            time = FLA_Clock();
            dspffrt2_(ad, &n, &ni, work, NULL);
            time = FLA_Clock() - time;
            time_min = min( time_min, time );

            /* Calculate error norm */
            dnrm = diff_norm_d( ad, fod, n, ni );

            if( dnrm >= THRESH2 )
            {
               pcode = strfail;
               break;
            }
            else if( dnrm >= THRESH )
            {
               pcode = strmarg;
               break;
            }
         }

         /* Calculate performance in GFLOPS */
         perf = ni/6.0f*(2.0*ni*ni-6.0*ni*n+3.0*ni+6.0*n*n-6.0*n+7) / time_min / FLOPS_PER_UNIT_PERF;

         libfla_test_print_result_info("SPFFRT2",
                                       (char *) 'd',
                                       "c",
                                       n,
                                       perf,
                                       time_min,
                                       dnrm,
                                       pcode,
                                       ni);
      }

      free(od);
      free(ad);

      free(fod);
   }
   return;
}

/*************************************************************************/
/*                    COMPLEX PRECISION ROUTINES                         */
/*************************************************************************/

float diff_norm_c( scomplex *ad, scomplex *fod, integer n, integer ni )
{
   scomplex *L, *D, *Lt;
   scomplex *ra1, *ra2, *rad;
   scomplex czero;
   scomplex cone;
   float snrm, anrm, *work;
   double rnrm, eps;

   integer di, i, j;

   char *tf = "No Transpose";
   char *nm = "1";
   char *ul = "lower";

   czero.real = 0.0;
   czero.imag = 0.0;
   cone.real = 1.0;
   cone.imag = 0.0;

   L   = (scomplex *) malloc( n * n * sizeof( scomplex ) );
   D   = (scomplex *) malloc( n * n * sizeof( scomplex ) );
   Lt  = (scomplex *) malloc( n * n * sizeof( scomplex ) );
   ra1 = (scomplex *) malloc( n * n * sizeof( scomplex ) );
   ra2 = (scomplex *) malloc( n * n * sizeof( scomplex ) );
   rad = (scomplex *) malloc( n * n * sizeof( scomplex ) );
   work = (float *) malloc( n * sizeof( float ) );

   set_identity_c(L,  n, n, n);
   set_identity_c(D,  n, n, n);
   set_identity_c(Lt, n, n, n);

   /* Form the original matrix from parts */

   /* Unpack L, D and L' */
   di = 0;
   for( i = 0; i < ni; i++ )
   {
      /* Form L */
      for( j = i; j < n; j++ )
      {
         L[i * n + j].real = ad[j - i + di].real;
         L[i * n + j].imag = ad[j - i + di].imag;
      }
      c_div_t(&L[i * n + i], &cone, &ad[di]);

      /* Form L' */
      for( j = i; j < n; j++ )
      {
         Lt[j * n + i].real = ad[j - i + di].real;
         Lt[j * n + i].imag = ad[j - i + di].imag;
      }
      c_div_t(&Lt[i * n + i], &cone, &ad[di]);

      D[i * n + i].real = ad[di].real;
      D[i * n + i].imag = ad[di].imag;

      di += n - i;
   }

   for( ; i < n; i++ )
   {
      for( j = i; j < n; j++ )
      {
         D[j + i * n].real = ad[di + j - i].real;
         D[j + i * n].imag = ad[di + j - i].imag;

         D[i + j * n].real = ad[di + j - i].real;
         D[i + j * n].imag = ad[di + j - i].imag;
      }

      di += n - i;
   }

   /* Product of L and D */
   cgemm_( tf, tf, &n, &n, &n, &cone, L, &n, D, &n, &czero, ra1, &n );

   /* Product of L*D and L' */
   cgemm_( tf, tf, &n, &n, &n, &cone, ra1, &n, Lt, &n, &czero, ra2, &n );

   /* Calculate difference matrix */
   for( i = 0; i < n; i++ )
   {
      for( j = 0; j < n; j++ )
      {
         rad[i * n + j].real = ra2[i * n + j].real - fod[i * n + j].real;
         rad[i * n + j].imag = ra2[i * n + j].imag - fod[i * n + j].imag;
      }
   }
   eps = slamch_("Epsilon");
   anrm = clansy_( nm, ul, &n, fod, &n, work );
   rnrm = clansy_( nm, ul, &n, rad, &n, work );

   /* Calculate the final norm */
   snrm = rnrm / n / anrm / eps;

   free( L );
   free( D );
   free( Lt );
   free( ra1 );
   free( ra2 );
   free( rad );
   free( work );

   return snrm;
}

void test_ldlt2_nopiv_ps_c( test_params_t *params )
{
   integer n, ni;
   scomplex *od, *ad;
   scomplex *fod;
   float  snrm;

   int    i, n_repeats;
   double time, perf;
   double time_min;
   char   *strpass = "PASS";
   char   *strmarg = "MARGINAL";
   char   *strfail = "FAILURE";
   char   *pcode;

   for( n = params->p_first; n <= params->p_max; n += params->p_inc )
   {
      integer pn = n * (n + 1) / 2;
      /* allocate packed matrices */
      od = (scomplex *) malloc( pn * sizeof( scomplex ) );
      ad = (scomplex *) malloc( pn * sizeof( scomplex ) );

      /* full sized matrix to store original */
      fod = (scomplex *) malloc( n * n * sizeof( scomplex ) );

      /* Initialize with random values */
      rand_sym_matrix_c( fod, n, n, n );

      /* Pack the symmetric matrix into input */
      pack_matrix_lt_c( fod, od, n, n );
      pack_matrix_lt_c( fod, ad, n, n );

      pcode = strpass;
      snrm = 0.0;
      perf = 0.0;
      if( (integer) params->p_nfact == -2 )
      {        
         n_repeats = (params->n_repeats < 5) ? params->n_repeats : 5;
         for( ni = 0; ni <= n; ni++ )
         {
            time_min   = 1e9;
            for( i = 0; i < n_repeats; i++ )
            {
               copy_matrix_c(od, ad, pn, 1, 1, 1);

               time = FLA_Clock();
               cspffrt2_(ad, &n, &ni, NULL, NULL);
               time = FLA_Clock() - time;
               time_min = min( time_min, time );

               /* Calculate error norm */
               snrm = diff_norm_c( ad, fod, n, ni );

               if( snrm >= THRESH2 )
               {
                  pcode = strfail;
                  break;
               }
               else if( snrm >= THRESH )
               {
                  pcode = strmarg;
                  break;
               }
            }

            /* Calculate performance in GFLOPS */
            perf = ni/3.0f*(4.0*ni*ni-12.0*ni*n+9.0*ni+12.0*n*n-18.0*n+8) / time_min / FLOPS_PER_UNIT_PERF;

            libfla_test_print_result_info("SPFFRT2",
                                          (char *) 'c',
                                          "c",
                                          n,
                                          perf,
                                          time_min,
                                          snrm,
                                          pcode,
                                          ni); 
            if( pcode[0] == 'M' )
               break;
         }
      }
      else
      {
         time_min   = 1e9;
         snrm = 0.0f;
         if( (integer) params->p_nfact == -1 )
         {
            ni = rand() % n;
         }
         else
         {
            ni = params->p_nfact;
         }

         for( i = 0; i < params->n_repeats; i++ )
         {
            copy_matrix_c(od, ad, pn, 1, 1, 1);

            time = FLA_Clock();
            cspffrt2_(ad, &n, &ni, NULL, NULL);
            time = FLA_Clock() - time;
            time_min = min( time_min, time );

            /* Calculate error norm */
            snrm = diff_norm_c( ad, fod, n, ni );

            if( snrm >= THRESH2 )
            {
               pcode = strfail;
               break;
            }
            else if( snrm >= THRESH )
            {
               pcode = strmarg;
               break;
            }
         }

         /* Calculate performance in GFLOPS */
         perf = ni/3.0f*(4.0*ni*ni-12.0*ni*n+9.0*ni+12.0*n*n-18.0*n+8) / time_min / FLOPS_PER_UNIT_PERF;

         libfla_test_print_result_info("SPFFRT2",
                                       (char *) 'c',
                                       "c",
                                       n,
                                       perf,
                                       time_min,
                                       snrm,
                                       pcode,
                                       ni);
      }

      free(od);
      free(ad);

      free(fod);
   }
   return;
}

/*************************************************************************/
/*               DOUBLE COMPLEX PRECISION ROUTINES                       */
/*************************************************************************/

double diff_norm_z( dcomplex *ad, dcomplex *fod, integer n, integer ni )
{
   dcomplex *L, *D, *Lt;
   dcomplex *ra1, *ra2, *rad;
   dcomplex czero;
   dcomplex cone;
   double dnrm, anrm;
   double rnrm, eps, *work;

   integer di, i, j;

   char *tf = "No Transpose";
   char *nm = "1";
   char *ul = "lower";

   czero.real = 0.0;
   czero.imag = 0.0;
   cone.real = 1.0;
   cone.imag = 0.0;

   L   = (dcomplex *) malloc( n * n * sizeof( dcomplex ) );
   D   = (dcomplex *) malloc( n * n * sizeof( dcomplex ) );
   Lt  = (dcomplex *) malloc( n * n * sizeof( dcomplex ) );
   ra1 = (dcomplex *) malloc( n * n * sizeof( dcomplex ) );
   ra2 = (dcomplex *) malloc( n * n * sizeof( dcomplex ) );
   rad = (dcomplex *) malloc( n * n * sizeof( dcomplex ) );
   work = (double *) malloc( n * sizeof( double ) );

   set_identity_z(L,  n, n, n);
   set_identity_z(D,  n, n, n);
   set_identity_z(Lt, n, n, n);

   /* Form the original matrix from parts */

   /* Unpack L, D and L' */
   di = 0;
   for( i = 0; i < ni; i++ )
   {
      /* Form L */
      for( j = i; j < n; j++ )
      {
         L[i * n + j].real = ad[j - i + di].real;
         L[i * n + j].imag = ad[j - i + di].imag;
      }
      z_div_t(&L[i * n + i], &cone, &ad[di]);

      /* Form L' */
      for( j = i; j < n; j++ )
      {
         Lt[j * n + i].real = ad[j - i + di].real;
         Lt[j * n + i].imag = ad[j - i + di].imag;
      }
      z_div_t(&Lt[i * n + i], &cone, &ad[di]);

      D[i * n + i].real = ad[di].real;
      D[i * n + i].imag = ad[di].imag;

      di += n - i;
   }

   for( ; i < n; i++ )
   {
      for( j = i; j < n; j++ )
      {
         D[j + i * n].real = ad[di + j - i].real;
         D[j + i * n].imag = ad[di + j - i].imag;

         D[i + j * n].real = ad[di + j - i].real;
         D[i + j * n].imag = ad[di + j - i].imag;
      }

      di += n - i;
   }

   /* Product of L and D */
   zgemm_( tf, tf, &n, &n, &n, &cone, L, &n, D, &n, &czero, ra1, &n );

   /* Product of L*D and L' */
   zgemm_( tf, tf, &n, &n, &n, &cone, ra1, &n, Lt, &n, &czero, ra2, &n );

   /* Calculate difference matrix */
   for( i = 0; i < n; i++ )
   {
      for( j = 0; j < n; j++ )
      {
         rad[i * n + j].real = ra2[i * n + j].real - fod[i * n + j].real;
         rad[i * n + j].imag = ra2[i * n + j].imag - fod[i * n + j].imag;
      }
   }
   eps = dlamch_("Epsilon");
   anrm = zlansy_( nm, ul, &n, fod, &n, work );
   rnrm = zlansy_( nm, ul, &n, rad, &n, work );

   /* Calculate the final norm */
   dnrm = rnrm / n / anrm / eps;

   free( L );
   free( D );
   free( Lt );
   free( ra1 );
   free( ra2 );
   free( rad );
   free( work );

   return dnrm;
}

void test_ldlt2_nopiv_ps_z( test_params_t *params )
{
   integer n, ni;
   dcomplex *od, *ad;
   dcomplex *fod, *work, *work1;
   double dnrm;

   int    i, n_repeats;
   double time, perf;
   double time_min;
   char   *strpass = "PASS";
   char   *strmarg = "MARGINAL";
   char   *strfail = "FAILURE";
   char   *pcode;

   for( n = params->p_first; n <= params->p_max; n += params->p_inc )
   {
      integer pn = n * (n + 1) / 2;
      /* allocate packed matrices */
      od = (dcomplex *) malloc( pn * sizeof( dcomplex ) );
      ad = (dcomplex *) malloc( pn * sizeof( dcomplex ) );
      work = (dcomplex *) malloc( n * sizeof( dcomplex ) );
      work1 = (dcomplex *) malloc( n * sizeof( dcomplex ) );

      /* full sized matrix to store original */
      fod = (dcomplex *) malloc( n * n * sizeof( dcomplex ) );

      /* Initialize with random values */
      rand_sym_matrix_z( fod, n, n, n );

      /* Pack the symmetric matrix into input */
      pack_matrix_lt_z( fod, od, n, n );
      pack_matrix_lt_z( fod, ad, n, n );

      pcode = strpass;
      dnrm = 0.0;
      perf = 0.0;
      if( (integer) params->p_nfact == -2 )
      {        
         n_repeats = (params->n_repeats < 5) ? params->n_repeats : 5;
         for( ni = 0; ni <= n; ni++ )
         {
            time_min   = 1e9;
            for( i = 0; i < n_repeats; i++ )
            {
               copy_matrix_z(od, ad, pn, 1, 1, 1);

               time = FLA_Clock();
               zspffrt2_(ad, &n, &ni, work, work1);
               time = FLA_Clock() - time;
               time_min = min( time_min, time );

               /* Calculate error norm */
               dnrm = diff_norm_z( ad, fod, n, ni );

               if( dnrm >= THRESH2 )
               {
                  pcode = strfail;
                  break;
               }
               else if( dnrm >= THRESH )
               {
                  pcode = strmarg;
                  break;
               }
            }
            /* Calculate performance in GFLOPS */
            perf = ni/3.0f*(4.0*ni*ni-12.0*ni*n+9.0*ni+12.0*n*n-18.0*n+8) / time_min / FLOPS_PER_UNIT_PERF;

            libfla_test_print_result_info("SPFFRT2",
                                          (char *) 'z',
                                          "c",
                                          n,
                                          perf,
                                          time_min,
                                          dnrm,
                                          pcode,
                                          ni);
            if( pcode[0] == 'M' )
               break;
         }
      }
      else
      {
         time_min   = 1e9;
         dnrm = 0.0f;
         if( (integer) params->p_nfact == -1 )
         {
            ni = rand() % n;
         }
         else
         {
            ni = params->p_nfact;
         }

         for( i = 0; i < params->n_repeats; i++ )
         {
            copy_matrix_z(od, ad, pn, 1, 1, 1);

            time = FLA_Clock();
            zspffrt2_(ad, &n, &ni, work, work1);
            time = FLA_Clock() - time;
            time_min = min( time_min, time );

            /* Calculate error norm */
            dnrm = diff_norm_z( ad, fod, n, ni );

            if( dnrm >= THRESH2 )
            {
               pcode = strfail;
               break;
            }
            else if( dnrm >= THRESH )
            {
               pcode = strmarg;
               break;
            }
         }

         /* Calculate performance in GFLOPS */
         perf = ni/3.0f*(4.0*ni*ni-12.0*ni*n+9.0*ni+12.0*n*n-18.0*n+8) / time_min / FLOPS_PER_UNIT_PERF;

         libfla_test_print_result_info("SPFFRT2",
                                       (char *) 'z',
                                       "c",
                                       n,
                                       perf,
                                       time_min,
                                       dnrm,
                                       pcode,
                                       ni);
      }

      free(od);
      free(ad);

      free(fod);
      free(work);
      free(work1);
   }
   return;
}

void libfla_test_ldlt2_nopiv_ps( test_params_t *params, test_op_t op )
{
   int i;

   libfla_test_output_info( "--- %s ---\n", op_str );
   libfla_test_output_info( "\n" );
   libfla_test_output_info( "%3sAPI%28s DATA_TYPE%4s SIZE%1s FLOPS%2s TIME(s)%6s ERROR%5s STATUS\n", "", "", "", "", "", "", "" );
   libfla_test_output_info( "%3s====%28s==========%4s====%1s=======%2s========%5s==========%2s========\n", "", "", "", "", "", "", "" );

   if ( op.fla_blk_ext == ENABLE )
   {
      for( i = 0; i < params->n_datatypes; ++i )
      {
         switch( params->datatype[i] )
         {
            case FLA_FLOAT:
               test_ldlt2_nopiv_ps_s( params );
               break;
            case FLA_DOUBLE:
               test_ldlt2_nopiv_ps_d( params );
               break;
            case FLA_COMPLEX:
               test_ldlt2_nopiv_ps_c( params );
               break;
            case FLA_DOUBLE_COMPLEX:
               test_ldlt2_nopiv_ps_z( params );
               break;
            default:
               printf("FLA_ERROR: Unknown Datatype Encountered\n");
               break;
         }
      }
   }
   libfla_test_output_info( "\n" );
}


