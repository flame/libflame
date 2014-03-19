/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#ifndef BLIS1_MACRO_DEFS_H
#define BLIS1_MACRO_DEFS_H

// --- Constants ---------------------------------------------------------------

#define BLIS1_NO_INTRINSICS  0
#define BLIS1_SSE_INTRINSICS 3

// --- boolean ---

#undef FALSE
#define FALSE 0

#undef TRUE
#define TRUE 1

/*
// --- trans ---

#define BLIS1_NO_TRANSPOSE      'n'
#define BLIS1_TRANSPOSE         't'
#define BLIS1_CONJ_NO_TRANSPOSE 'c'
#define BLIS1_CONJ_TRANSPOSE    'h'

// --- conj ---

#define BLIS1_NO_CONJUGATE      'n'
#define BLIS1_CONJUGATE         'c'

// --- uplo ---

#define BLIS1_LOWER_TRIANGULAR  'l'
#define BLIS1_UPPER_TRIANGULAR  'u'

// --- side ---

#define BLIS1_LEFT              'l'
#define BLIS1_RIGHT             'r'

// --- diag ---

#define BLIS1_NONUNIT_DIAG      'n'
#define BLIS1_UNIT_DIAG         'u'
#define BLIS1_ZERO_DIAG         'z'
*/

// --- Functional macros -------------------------------------------------------

// --- Type-agnostic ---

// min, max, abs

#define bl1_min( a, b )  ( (a) < (b) ? (a) : (b) )
#define bl1_max( a, b )  ( (a) > (b) ? (a) : (b) )
#define bl1_abs( a )     ( (a) <= 0 ? -(a) : (a) )

// fmin, fmax, fabs

#define bl1_fmin( a, b ) bl1_min( a, b )
#define bl1_fmax( a, b ) bl1_max( a, b )
#define bl1_fabs( a )    ( (a) <= 0.0 ? -(a) : (a) )

// fminabs, fmaxabs
#define bl1_fminabs( a, b ) \
\
    bl1_fmin( bl1_fabs( a ), \
              bl1_fabs( b ) )

#define bl1_fmaxabs( a, b ) \
\
    bl1_fmax( bl1_fabs( a ), \
              bl1_fabs( b ) )

// --- Type-dependent ---

// --- neg1 ---

// void bl1_sneg1( float* x );
#define bl1_sneg1( x ) \
*(x)     *= -1.0F;

// void bl1_dneg1( double* x );
#define bl1_dneg1( x ) \
*(x)     *= -1.0;

// void bl1_cneg1( scomplex* x );
#define bl1_cneg1( x ) \
(x)->real *= -1.0F; \
(x)->imag *= -1.0F;

// void bl1_zneg1( dcomplex* x );
#define bl1_zneg1( x ) \
(x)->real *= -1.0; \
(x)->imag *= -1.0;

// --- neg2 ---

// void bl1_sneg2( float* x, float* y );
#define bl1_sneg2( x, y ) \
*(y)      = -1.0F * *(x);

// void bl1_dneg2( double* x, double* y );
#define bl1_dneg2( x, y ) \
*(y)      = -1.0  * *(x);

// void bl1_cneg2( scomplex* x, scomplex* y );
#define bl1_cneg2( x, y ) \
(y)->real = -1.0F * (x)->real; \
(y)->imag = -1.0F * (x)->imag;

// void bl1_zneg2( dcomplex* x, dcomplex* y );
#define bl1_zneg2( x, y ) \
(y)->real = -1.0  * (x)->real; \
(y)->imag = -1.0  * (x)->imag;

// --- sqrte ---

// void bl1_ssqrte( float* alpha, int* error );
#define bl1_ssqrte( alpha, error ) \
if ( *(alpha)      <= 0.0F || isnan( *(alpha) ) ) {  *(error) = FLA_FAILURE; } \
else { *(alpha)      =  ( float ) sqrt( *(alpha) );  *(error) = FLA_SUCCESS; }

// void bl1_dsqrte( double* alpha, int* error );
#define bl1_dsqrte( alpha, error ) \
if ( *(alpha)      <= 0.0 || isnan( *(alpha) ) ) {   *(error) = FLA_FAILURE; } \
else { *(alpha)      = ( double ) sqrt( *(alpha) );  *(error) = FLA_SUCCESS; }

// void bl1_csqrte( scomplex* alpha, int* error );
#define bl1_csqrte( alpha, error ) \
if ( (alpha)->real <= 0.0F || isnan( (alpha)->real) ) \
{                     *(error) = FLA_FAILURE; } \
else { \
(alpha)->real =  ( float ) sqrt( (alpha)->real ); \
(alpha)->imag = 0.0F; *(error) = FLA_SUCCESS; }

// void bl1_zsqrte( dcomplex* alpha, int* error );
#define bl1_zsqrte( alpha, error ) \
if ( (alpha)->real <= 0.0 || isnan( (alpha)->real) )  \
{                     *(error) = FLA_FAILURE; } \
else { \
(alpha)->real = ( double ) sqrt( (alpha)->real ); \
(alpha)->imag = 0.0;  *(error) = FLA_SUCCESS; }

// --- absval2 ---

// void bl1_sabsval2( float* alpha, float* absval );
#define bl1_sabsval2( alpha, absval ) \
*(absval) = ( float ) fabs( ( double ) *(alpha) );

// void bl1_dabsval2( double* alpha, double* absval );
#define bl1_dabsval2( alpha, absval ) \
*(absval) = fabs( *(alpha) );

// void bl1_cabsval2( scomplex* x, scomplex* a );
#define bl1_cabsval2( x, a ) \
{ \
	float  s   = bl1_fmaxabs( (x)->real, (x)->imag ); \
	float  mag = sqrtf( s ) * \
	             sqrtf( ( (x)->real / s ) * (x)->real + \
	                    ( (x)->imag / s ) * (x)->imag ); \
	(a)->real   = mag; \
	(a)->imag   = 0.0F; \
}

// void bl1_csabsval2( scomplex* x, float* a );
#define bl1_csabsval2( x, a ) \
{ \
	float  s   = bl1_fmaxabs( (x)->real, (x)->imag ); \
	float  mag = sqrtf( s ) * \
	             sqrtf( ( (x)->real / s ) * (x)->real + \
	                    ( (x)->imag / s ) * (x)->imag ); \
	*(a)       = mag; \
}

// void bl1_zabsval2( dcomplex* x, dcomplex* a );
#define bl1_zabsval2( x, a ) \
{ \
	double s   = bl1_fmaxabs( (x)->real, (x)->imag ); \
	double mag = sqrt( s ) * \
	             sqrt( ( (x)->real / s ) * (x)->real + \
	                   ( (x)->imag / s ) * (x)->imag ); \
	(a)->real   = mag; \
	(a)->imag   = 0.0; \
}

// void bl1_zdabsval2( dcomplex* x, double* a );
#define bl1_zdabsval2( x, a ) \
{ \
	double s   = bl1_fmaxabs( (x)->real, (x)->imag ); \
	double mag = sqrt( s ) * \
	             sqrt( ( (x)->real / s ) * (x)->real + \
	                   ( (x)->imag / s ) * (x)->imag ); \
	*(a)       = mag; \
}


// --- absqr ---

// void bl1_sabsqr( float* alpha );
#define bl1_sabsqr( alpha ) \
*(alpha) = *(alpha) * *(alpha);

// void bl1_dabsqr( double* alpha );
#define bl1_dabsqr( alpha ) \
*(alpha) = *(alpha) * *(alpha);

// void bl1_cabsqr( scomplex* alpha );
#define bl1_cabsqr( alpha ) \
(alpha)->real = (alpha)->real * (alpha)->real + (alpha)->imag * (alpha)->imag; \
(alpha)->imag = 0.0F;

// void bl1_zabsqr( dcomplex* alpha );
#define bl1_zabsqr( alpha ) \
(alpha)->real = (alpha)->real * (alpha)->real + (alpha)->imag * (alpha)->imag; \
(alpha)->imag = 0.0;

// --- invscals ---

// void bl1_sinvscals( float* a, float* y );
#define bl1_sinvscals( a, y ) \
*(y) = *(y) / *(a);

// void bl1_dinvscals( double* a, double* y );
#define bl1_dinvscals( a, y ) \
*(y) = *(y) / *(a);

// void bl1_csinvscals( float* a, scomplex* y );
#define bl1_csinvscals( a, y ) \
{ \
(y)->real = (y)->real / *(a); \
(y)->imag = (y)->imag / *(a); \
}

// void bl1_cinvscals( scomplex* a, scomplex* y );
#define bl1_cinvscals( a, y ) \
{ \
	float  s     = bl1_fmaxabs( (a)->real, (a)->imag ); \
	float  ar_s  = (a)->real / s; \
	float  ai_s  = (a)->imag / s; \
	float  yrt   = (y)->real; \
	float  temp  = ( ar_s * (a)->real + ai_s * (a)->imag ); \
	(y)->real    = ( (yrt)     * ar_s + (y)->imag * ai_s ) / temp; \
	(y)->imag    = ( (y)->imag * ar_s - (yrt)     * ai_s ) / temp; \
}

// void bl1_zdinvscals( double* a, dcomplex* y );
#define bl1_zdinvscals( a, y ) \
{ \
(y)->real = (y)->real / *(a); \
(y)->imag = (y)->imag / *(a); \
}

// void bl1_zinvscals( dcomplex* a, dcomplex* y );
#define bl1_zinvscals( a, y ) \
{ \
	double s     = bl1_fmaxabs( (a)->real, (a)->imag ); \
	double ar_s  = (a)->real / s; \
	double ai_s  = (a)->imag / s; \
	double yrt   = (y)->real; \
	double temp  = ( ar_s * (a)->real + ai_s * (a)->imag ); \
	(y)->real    = ( (yrt)     * ar_s + (y)->imag * ai_s ) / temp; \
	(y)->imag    = ( (y)->imag * ar_s - (yrt)     * ai_s ) / temp; \
}

// --- div3 ---

// void bl1_sdiv3( float* x, float* y, float* a );
#define bl1_sdiv3( x, y, a ) \
*(a) = *(x) / *(y);

// void bl1_ddiv3( double* x, double* y, double* a );
#define bl1_ddiv3( x, y, a ) \
*(a) = *(x) / *(y);

// void bl1_cdiv3( scomplex* x, scomplex* y, scomplex* a );
// a = x / y;
#define bl1_cdiv3( x, y, a ) \
{ \
	*a = *x; \
	bl1_cinvscals( y, a ); \
}

// void bl1_zdiv3( dcomplex* x, dcomplex* y, dcomplex* a );
#define bl1_zdiv3( x, y, a ) \
{ \
	*a = *x; \
	bl1_zinvscals( y, a ); \
}

// --- add3 ---

// void bl1_sadd3( float* x, float* y, float* a );
#define bl1_sadd3( x, y, a ) \
*(a) = *(x) + *(y);

// void bl1_dadd3( double* x, double* y, double* a );
#define bl1_dadd3( x, y, a ) \
*(a) = *(x) + *(y);

// void bl1_cadd3( scomplex* x, scomplex* y, scomplex* a );
#define bl1_cadd3( x, y, a ) \
{ \
(a)->real = (x)->real + (y)->real; \
(a)->imag = (x)->imag + (y)->imag; \
}

// void bl1_zadd3( dcomplex* x, dcomplex* y, dcomplex* a );
#define bl1_zadd3( x, y, a ) \
{ \
(a)->real = (x)->real + (y)->real; \
(a)->imag = (x)->imag + (y)->imag; \
}

// --- copys ---

// void bl1_scopys( conj1_t conj, float* x, float* y );
#define bl1_scopys( conj, x, y ) \
*(y) = *(x);

// void bl1_dcopys( conj1_t conj, double* x, double* y );
#define bl1_dcopys( conj, x, y ) \
*(y) = *(x);

// void bl1_ccopys( conj1_t conj, scomplex* x, scomplex* y );
#define bl1_ccopys( conj, x, y ) \
*(y) = *(x); \
if ( bl1_is_conj( conj ) ) (y)->imag *= -1.0F;

// void bl1_zcopys( conj1_t conj, dcomplex* x, dcomplex* y );
#define bl1_zcopys( conj, x, y ) \
*(y) = *(x); \
if ( bl1_is_conj( conj ) ) (y)->imag *= -1.0;

// --- scals ---

// void bl1_sscals( float* a, float* y );
#define bl1_sscals( a, y ) \
*(y) = *(a) * *(y);

// void bl1_dscals( double* a, double* y );
#define bl1_dscals( a, y ) \
*(y) = *(a) * *(y);

// void bl1_csscals( float* a, scomplex* y );
#define bl1_csscals( a, y ) \
{ \
(y)->real = *(a) * (y)->real; \
(y)->imag = *(a) * (y)->imag; \
}

// void bl1_cscals( scomplex* a, scomplex* y );
#define bl1_cscals( a, y ) \
{ \
float tempr = (a)->real * (y)->real - (a)->imag * (y)->imag; \
float tempi = (a)->imag * (y)->real + (a)->real * (y)->imag; \
(y)->real = tempr; \
(y)->imag = tempi; \
}

// void bl1_zdscals( double* a, dcomplex* y );
#define bl1_zdscals( a, y ) \
{ \
(y)->real = *(a) * (y)->real; \
(y)->imag = *(a) * (y)->imag; \
}

// void bl1_zscals( dcomplex* a, dcomplex* y );
#define bl1_zscals( a, y ) \
{ \
double tempr = (a)->real * (y)->real - (a)->imag * (y)->imag; \
double tempi = (a)->imag * (y)->real + (a)->real * (y)->imag; \
(y)->real = tempr; \
(y)->imag = tempi; \
}

// --- mult3 ---

// void bl1_smult3( float* x, float* y, float* a );
#define bl1_smult3( x, y, a ) \
*(a) = *(x) * *(y);

// void bl1_dmult3( double* x, double* y, double* a );
#define bl1_dmult3( x, y, a ) \
*(a) = *(x) * *(y);

// void bl1_cmult3( scomplex* x, scomplex* y, scomplex* a );
#define bl1_cmult3( x, y, a ) \
{ \
float tempr = (x)->real * (y)->real - (x)->imag * (y)->imag; \
float tempi = (x)->imag * (y)->real + (x)->real * (y)->imag; \
(a)->real = tempr; \
(a)->imag = tempi; \
}

// void bl1_zmult3( dcomplex* x, dcomplex* y, dcomplex* a );
#define bl1_zmult3( x, y, a ) \
{ \
double tempr = (x)->real * (y)->real - (x)->imag * (y)->imag; \
double tempi = (x)->imag * (y)->real + (x)->real * (y)->imag; \
(a)->real = tempr; \
(a)->imag = tempi; \
}

// --- mult4 ---

// void bl1_smult4( float* alpha, float* x, float* y1, float* y2 );
#define bl1_smult4( alpha, x, y1, y2 ) \
*(y2) = *(y1) + *(alpha) * *(x);

// void bl1_dmult4( double* alpha, double* x, double* y1, double* y2 );
#define bl1_dmult4( alpha, x, y1, y2 ) \
*(y2) = *(y1) + *(alpha) * *(x);

// void bl1_cmult4( scomplex* alpha, scomplex* x, scomplex* y1, scomplex* y2 );
#define bl1_cmult4( alpha, x, y1, y2 ) \
{ \
(y2)->real = (y1)->real + (alpha)->real * (x)->real - (alpha)->imag * (x)->imag; \
(y2)->imag = (y1)->imag + (alpha)->imag * (x)->real + (alpha)->real * (x)->imag; \
}

// void bl1_zmult4( dcomplex* alpha, dcomplex* x, dcomplex* y1, dcomplex* y2 );
#define bl1_zmult4( alpha, x, y1, y2 ) \
{ \
(y2)->real = (y1)->real + (alpha)->real * (x)->real - (alpha)->imag * (x)->imag; \
(y2)->imag = (y1)->imag + (alpha)->imag * (x)->real + (alpha)->real * (x)->imag; \
}

// --- conjs ---

// void bl1_sconjs( float* a );
#define bl1_sconjs( a ) \
;

// void bl1_dconjs( double* a );
#define bl1_dconjs( a ) \
;

// void bl1_cconjs( scomplex* a );
#define bl1_cconjs( a ) \
(a)->imag *= -1.0F;

// void bl1_zconjs( dcomplex* a );
#define bl1_zconjs( a ) \
(a)->imag *= -1.0;

// --- copyconj ---

// void bl1_scopyconj( float* x, float* y );
#define bl1_scopyconj( x, y ) \
*(y) = *(x);

// void bl1_dcopyconj( double* x, double* y );
#define bl1_dcopyconj( x, y ) \
*(y) = *(x);

// void bl1_ccopyconj( scomplex* x, scomplex* y );
#define bl1_ccopyconj( x, y ) \
(y)->real =         (x)->real; \
(y)->imag = -1.0F * (x)->imag;

// void bl1_zcopyconj( dcomplex* x, dcomplex* y );
#define bl1_zcopyconj( x, y ) \
(y)->real =         (x)->real; \
(y)->imag = -1.0  * (x)->imag;

// --- eq1 ---

// void bl1_seq1( float* alpha );
#define bl1_seq1( alpha ) \
  ( *alpha == 1.0F )

// void bl1_deq1( double* alpha );
#define bl1_deq1( alpha ) \
  ( *alpha == 1.0 )

// void bl1_ceq1( scomplex* alpha );
#define bl1_ceq1( alpha ) \
  ( (alpha)->real == 1.0F && (alpha)->imag == 0.0F )

// void bl1_zeq1( dcomplex* alpha );
#define bl1_zeq1( alpha ) \
  ( (alpha)->real == 1.0 && (alpha)->imag == 0.0 )

// --- Swapping/toggle macros --------------------------------------------------

// --- swap_pointers ---

#define bl1_sswap_pointers( a, b ) \
{ \
float* temp = (a); \
(a) = (b); \
(b) = temp; \
}

#define bl1_dswap_pointers( a, b ) \
{ \
double* temp = (a); \
(a) = (b); \
(b) = temp; \
}

#define bl1_cswap_pointers( a, b ) \
{ \
void* temp = (a); \
(a) = (b); \
(b) = temp; \
}

#define bl1_zswap_pointers( a, b ) \
{ \
void* temp = (a); \
(a) = (b); \
(b) = temp; \
}

// --- swap_ints ---

#define bl1_swap_ints( a, b ) \
{ \
int temp = (a); \
(a) = (b); \
(b) = temp; \
}

// --- swap_trans ---

#define bl1_swap_trans( a, b ) \
{ \
trans1_t temp = (a); \
(a) = (b); \
(b) = temp; \
}

// --- swap_conj ---

#define bl1_swap_conj( a, b ) \
{ \
conj1_t temp = (a); \
(a) = (b); \
(b) = temp; \
}

// --- toggle_side ---

#define bl1_toggle_side( side ) \
{ \
if ( bl1_is_left( side ) ) side = BLIS1_RIGHT; \
else                       side = BLIS1_LEFT; \
}

// --- toggle_uplo ---

#define bl1_toggle_uplo( uplo ) \
{ \
if ( bl1_is_lower( uplo ) ) uplo = BLIS1_UPPER_TRIANGULAR; \
else                        uplo = BLIS1_LOWER_TRIANGULAR; \
}

// --- toggle_trans ---
#define bl1_toggle_trans( trans ) \
{ \
if      ( bl1_is_notrans( trans ) )     trans = BLIS1_TRANSPOSE; \
else if ( bl1_is_trans( trans ) )       trans = BLIS1_NO_TRANSPOSE; \
else if ( bl1_is_conjnotrans( trans ) ) trans = BLIS1_CONJ_TRANSPOSE; \
else                                    trans = BLIS1_CONJ_NO_TRANSPOSE; \
}

// --- toggle_conjtrans ---
#define bl1_toggle_conjtrans( trans ) \
{ \
if      ( bl1_is_notrans( trans ) )     trans = BLIS1_CONJ_TRANSPOSE; \
else                                    trans = BLIS1_NO_TRANSPOSE; \
}

// --- toggle_conj ---

#define bl1_toggle_conj( conj ) \
{ \
if ( bl1_is_conj( conj ) ) conj = BLIS1_NO_CONJUGATE; \
else                       conj = BLIS1_CONJUGATE; \
}

#endif // #ifndef BLIS1_MACRO_DEFS_H
