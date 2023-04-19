/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/



#ifndef FLASH_LAPACK2FLASH_UTIL_DEFS_H
#define FLASH_LAPACK2FLASH_UTIL_DEFS_H

// --- LAPACK interface-------------- ------------------------------
//#define FLA_LAPACK2FLAME_SUPPORT_COMPLEX

// --- LAPACK fortran datatype macros ------------------------------
/*
typedef real          lapack_stype;
typedef doublereal    lapack_dtype;
typedef complex       lapack_ctype;
typedef doublecomplex lapack_ztype;

typedef real          lapack_sreal;
typedef doublereal    lapack_dreal;
typedef real          lapack_creal;
typedef doublereal    lapack_zreal;
*/

typedef float    lapack_stype;
typedef double   lapack_dtype;
typedef scomplex lapack_ctype;
typedef dcomplex lapack_ztype;

typedef float    lapack_sreal;
typedef double   lapack_dreal;
typedef float    lapack_creal;
typedef double   lapack_zreal;

#define LAPACK_strans "T"
#define LAPACK_dtrans "T"
#define LAPACK_ctrans "C"
#define LAPACK_ztrans "C"

#define PREFIX2LAPACK_TRANS(prefix) LAPACK_ ## prefix ## trans

#define PREFIX2LAPACK_TYPEDEF(prefix) lapack_ ## prefix ## type
#define PREFIX2LAPACK_REALDEF(prefix) lapack_ ## prefix ## real

#define FLAME_stype FLA_FLOAT
#define FLAME_dtype FLA_DOUBLE
#define FLAME_ctype FLA_COMPLEX
#define FLAME_ztype FLA_DOUBLE_COMPLEX

#define FLAME_srealtype FLA_FLOAT
#define FLAME_drealtype FLA_DOUBLE
#define FLAME_crealtype FLA_FLOAT
#define FLAME_zrealtype FLA_DOUBLE

#define PREFIX2FLAME_DATATYPE(prefix) FLAME_ ## prefix ## type
#define PREFIX2FLAME_REALTYPE(prefix) FLAME_ ## prefix ## realtype

#define s_IS_COMPLEX 0
#define d_IS_COMPLEX 0
#define c_IS_COMPLEX 1
#define z_IS_COMPLEX 1

#define IS_COMPLEX_PREFIX(prefix)  prefix ## _IS_COMPLEX

#define FLAME_is_szero( val ) (*val == 0.f)
#define FLAME_is_dzero( val ) (*val == 0.0)
#define FLAME_is_czero( val ) (val->real == 0.f && val->imag == 0.f)
#define FLAME_is_zzero( val ) (val->real == 0.0 && val->imag == 0.0)

#define PREFIX2FLAME_IS_ZERO(prefix, val) FLAME_is_ ## prefix ## zero(val)

extern int FLAME_invert_stau( FLA_Obj t );
extern int FLAME_invert_dtau( FLA_Obj t );
extern int FLAME_invert_ctau( FLA_Obj t );
extern int FLAME_invert_ztau( FLA_Obj t );

#define PREFIX2FLAME_INVERT_TAU(prefix, val) FLAME_invert_ ## prefix ## tau(val)

extern int FLAME_QR_piv_preorder( FLA_Obj A, int *jpiv_lapack, int *jpiv_fla );

FLA_Bool FLASH_Check_offload_to_gpu( dim_t blocksize, int m, int n, unsigned int crossover );
void FLASH_Toggle_gpu_offload( FLA_Bool toggle );

#endif
