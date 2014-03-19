/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#ifndef BLIS1_TYPE_DEFS_H
#define BLIS1_TYPE_DEFS_H

// --- Basic type definitions -------------------------------------------------

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

#define BLIS1_TRANS_BEGIN 100
#define BLIS1_UPLO_BEGIN  200
#define BLIS1_SIDE_BEGIN  300
#define BLIS1_DIAG_BEGIN  400
#define BLIS1_CONJ_BEGIN  500

typedef enum
{
	BLIS1_NO_TRANSPOSE = BLIS1_TRANS_BEGIN,
	BLIS1_TRANSPOSE,
	BLIS1_CONJ_NO_TRANSPOSE,
	BLIS1_CONJ_TRANSPOSE
} trans1_t;

typedef enum
{
	BLIS1_LOWER_TRIANGULAR = BLIS1_UPLO_BEGIN,
	BLIS1_UPPER_TRIANGULAR
} uplo1_t;

typedef enum
{
	BLIS1_LEFT = BLIS1_SIDE_BEGIN,
	BLIS1_RIGHT
} side1_t;

typedef enum
{
	BLIS1_NONUNIT_DIAG = BLIS1_DIAG_BEGIN,
	BLIS1_UNIT_DIAG,
	BLIS1_ZERO_DIAG
} diag1_t;

typedef enum
{
	BLIS1_NO_CONJUGATE = BLIS1_CONJ_BEGIN,
	BLIS1_CONJUGATE
} conj1_t;





// --- Intrinsic/assembly definitions ----------------------------------------

/*
#ifndef BLIS1_FROM_LIBFLAME
    #error "NOT using blis from libflame"
#else
    #error "using blis from libflame"
#endif
*/

/*
#if BLIS1_VECTOR_INTRINSIC_TYPE == BLIS1_SSE_INTRINSICS
    #error "using sse in blis"
#elif  BLIS1_VECTOR_INTRINSIC_TYPE == BLIS1_NO_INTRINSICS
    #error "NOT using sse in blis"
#else
    #error "undefined case!"
#endif
*/

// Only define vector intrinsics types if they are not already provided by
// libflame.
#ifndef BLIS1_FROM_LIBFLAME

#if BLIS1_VECTOR_INTRINSIC_TYPE == BLIS1_SSE_INTRINSICS

#include "pmmintrin.h"
typedef union
{
    __m128d v; 
    double  d[2];
} v2df_t;
#endif

#endif


// --- Complex type definitions -----------------------------------------------

// Only define complex types if they are not already provided by libflame.
//#ifndef BLIS1_ENABLE_USE_OF_LIBFLAME_TYPES
#ifndef BLIS1_FROM_LIBFLAME

typedef struct scomplex
{
  float real, imag;
} scomplex;

typedef struct dcomplex
{
  double real, imag;
} dcomplex;

#endif


#endif // BLIS1_TYPE_DEFS_H
