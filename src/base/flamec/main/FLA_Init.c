/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"


static TLS_CLASS_SPEC FLA_Bool FLA_initialized = FALSE;

FLA_Obj TLS_CLASS_SPEC FLA_THREE = {};
FLA_Obj TLS_CLASS_SPEC FLA_TWO = {};
FLA_Obj TLS_CLASS_SPEC FLA_ONE = {};
FLA_Obj TLS_CLASS_SPEC FLA_ONE_HALF = {};
FLA_Obj TLS_CLASS_SPEC FLA_ZERO = {};
FLA_Obj TLS_CLASS_SPEC FLA_MINUS_ONE_HALF = {};
FLA_Obj TLS_CLASS_SPEC FLA_MINUS_ONE = {};
FLA_Obj TLS_CLASS_SPEC FLA_MINUS_TWO = {};
FLA_Obj TLS_CLASS_SPEC FLA_MINUS_THREE = {};

FLA_Obj TLS_CLASS_SPEC FLA_EPSILON = {};
FLA_Obj TLS_CLASS_SPEC FLA_SAFE_MIN = {};
FLA_Obj TLS_CLASS_SPEC FLA_SAFE_MIN_SQUARE = {};
FLA_Obj TLS_CLASS_SPEC FLA_SAFE_INV_MIN = {};
FLA_Obj TLS_CLASS_SPEC FLA_SAFE_INV_MIN_SQUARE = {};
FLA_Obj TLS_CLASS_SPEC FLA_UNDERFLOW_THRES = {};
FLA_Obj TLS_CLASS_SPEC FLA_OVERFLOW_THRES = {};
FLA_Obj TLS_CLASS_SPEC FLA_UNDERFLOW_SQUARE_THRES = {};
FLA_Obj TLS_CLASS_SPEC FLA_OVERFLOW_SQUARE_THRES = {};

const TLS_CLASS_SPEC float    fzero = 0.0f;
const TLS_CLASS_SPEC double   dzero = 0.0;
const TLS_CLASS_SPEC scomplex czero = { 0.0f, 0.0f };
const TLS_CLASS_SPEC dcomplex zzero = { 0.0 , 0.0  };

struct _LF_VERSION
{
	char version[1000];
	int once;
};
typedef struct _LF_VERSION LF_VERSION;

#define VERSION_MAKE_STR(x) _VERSION_MAKE_STR(x)
#define _VERSION_MAKE_STR(x) #x

/* *************************************************************************

   FLA_Init()

 *************************************************************************** */

void FLA_Init()
{
  if ( FLA_initialized == TRUE ) return;
  
  FLA_initialized = TRUE;

  FLA_Error_messages_init();

  FLA_Memory_leak_counter_init();

  FLA_Init_constants();

  FLA_Cntl_init();

#if FLA_VECTOR_INTRINSIC_TYPE == FLA_SSE_INTRINSICS
  _MM_SET_FLUSH_ZERO_MODE( _MM_FLUSH_ZERO_ON );
#endif

#ifdef FLA_ENABLE_SUPERMATRIX
  FLASH_Queue_init();
#endif
}

/* *************************************************************************

  FLA_Finalize()

 *************************************************************************** */

void FLA_Finalize()
{
  if ( FLA_initialized == FALSE ) return;

  FLA_initialized = FALSE;

  FLA_Finalize_constants();

  FLA_Cntl_finalize();

#ifdef FLA_ENABLE_SUPERMATRIX
  FLASH_Queue_finalize();
#endif

  FLA_Memory_leak_counter_finalize();
}

/* *************************************************************************

  FLA_Init_safe()

 *************************************************************************** */

void FLA_Init_safe( FLA_Error* init_result )
{
  if ( FLA_Initialized() )
  {
    *init_result = FLA_FAILURE;
  }
  else
  {
    FLA_Init();
    *init_result = FLA_SUCCESS;
  }
}

/* *************************************************************************

  FLA_Finalize_safe()

 *************************************************************************** */

void FLA_Finalize_safe( FLA_Error init_result )
{
  if ( init_result == FLA_SUCCESS )
    FLA_Finalize();
}

/* *************************************************************************

   FLA_Initialized()

 *************************************************************************** */

FLA_Bool FLA_Initialized( void )
{
  return FLA_initialized;
}

/* *************************************************************************

   FLA_Init_constants()

 *************************************************************************** */

void FLA_Init_constants()
{
  FLA_Obj_create_constant(  3.0, &FLA_THREE );
  FLA_Obj_create_constant(  2.0, &FLA_TWO );
  FLA_Obj_create_constant(  1.0, &FLA_ONE );
  FLA_Obj_create_constant(  0.5, &FLA_ONE_HALF );
  FLA_Obj_create_constant(  0.0, &FLA_ZERO );
  FLA_Obj_create_constant( -0.5, &FLA_MINUS_ONE_HALF );
  FLA_Obj_create_constant( -1.0, &FLA_MINUS_ONE );
  FLA_Obj_create_constant( -2.0, &FLA_MINUS_TWO );
  FLA_Obj_create_constant( -3.0, &FLA_MINUS_THREE );


  { 
    float  
      eps_f, 
      sfmin_f = FLT_MIN, sfmin_f2,
      small_f = ( 1.0F / FLT_MAX ), 
      under_f = FLT_MIN, 
      over_f  = FLT_MAX;

    double 
      eps_d, 
      sfmin_d = DBL_MIN, sfmin_d2,
      small_d = ( 1.0  / DBL_MAX ), 
      under_d = DBL_MIN,
      over_d  = DBL_MAX;

    if ( FLT_ROUNDS == 1 )
    {
      eps_f = FLT_EPSILON*0.5F;
      eps_d = DBL_EPSILON*0.5;
    }
    else 
    {
      eps_f = FLT_EPSILON;
      eps_d = DBL_EPSILON;
    }

    if ( small_f >= sfmin_f ) sfmin_f = small_f * ( 1.0F + eps_f );
    if ( small_d >= sfmin_d ) sfmin_d = small_d * ( 1.0  + eps_d );

    sfmin_f  = sfmin_f/eps_f;
    sfmin_d  = sfmin_d/eps_d;

    sfmin_f2 = sqrt( sfmin_f );
    sfmin_d2 = sqrt( sfmin_d );

    FLA_Obj_create_constant_ext( eps_f,           eps_d,           &FLA_EPSILON );

    FLA_Obj_create_constant_ext( sfmin_f,         sfmin_d,         &FLA_SAFE_MIN );
    FLA_Obj_create_constant_ext( 1.0F/sfmin_f,    1.0/sfmin_d,     &FLA_SAFE_INV_MIN );

    FLA_Obj_create_constant_ext( sfmin_f2,        sfmin_d2,        &FLA_SAFE_MIN_SQUARE );
    FLA_Obj_create_constant_ext( 1.0F/sfmin_f2,   1.0/sfmin_d2,    &FLA_SAFE_INV_MIN_SQUARE );
    
    FLA_Obj_create_constant_ext( under_f,         under_d,         &FLA_UNDERFLOW_THRES );
    FLA_Obj_create_constant_ext( over_f,          over_d,          &FLA_OVERFLOW_THRES );

    FLA_Obj_create_constant_ext( sqrt( under_f ), sqrt( under_d ), &FLA_UNDERFLOW_SQUARE_THRES );
    FLA_Obj_create_constant_ext( sqrt( over_f ),  sqrt( over_d ),  &FLA_OVERFLOW_SQUARE_THRES );
  } 
}

/* *************************************************************************

   FLA_Finalize_constants()

 *************************************************************************** */

void FLA_Finalize_constants()
{
  FLA_Obj_free( &FLA_THREE );
  FLA_Obj_free( &FLA_TWO );
  FLA_Obj_free( &FLA_ONE );
  FLA_Obj_free( &FLA_ONE_HALF );
  FLA_Obj_free( &FLA_ZERO );
  FLA_Obj_free( &FLA_MINUS_ONE_HALF );
  FLA_Obj_free( &FLA_MINUS_ONE );
  FLA_Obj_free( &FLA_MINUS_TWO );
  FLA_Obj_free( &FLA_MINUS_THREE );

  FLA_Obj_free( &FLA_EPSILON );
  FLA_Obj_free( &FLA_SAFE_MIN );
  FLA_Obj_free( &FLA_SAFE_MIN_SQUARE );
  FLA_Obj_free( &FLA_SAFE_INV_MIN );
  FLA_Obj_free( &FLA_SAFE_INV_MIN_SQUARE );
  FLA_Obj_free( &FLA_UNDERFLOW_THRES );
  FLA_Obj_free( &FLA_OVERFLOW_THRES );  
  FLA_Obj_free( &FLA_UNDERFLOW_SQUARE_THRES );
  FLA_Obj_free( &FLA_OVERFLOW_SQUARE_THRES );
}

char*     FLA_Get_AOCL_Version( void )
{
     static TLS_CLASS_SPEC LF_VERSION lflibversion;
     if (lflibversion.once)
     {
	    return lflibversion.version;
     }

     char lfmainversion[] = "AOCL-LAPACK ";
     char* lfversion = lflibversion.version;
     char lapackversion[] = ", supports LAPACK 3.11.0";
     int length, i;

     length = 0;
     for (i = 0; lfmainversion[length] != '\0'; ++i, ++length) 
     {
	 lfversion[length] = lfmainversion[i];
     }

#ifdef FLA_LIBFLAME_VERSION
#ifdef FLA_ENABLE_WINDOWS_BUILD
     char configlfversion[] = VERSION_MAKE_STR(FLA_LIBFLAME_VERSION); //Quotions in CMake are a problem. Hence strinifying
#else
     char configlfversion[] = FLA_LIBFLAME_VERSION;
#endif
     for (i = 0; configlfversion[i] != '\0'; ++i, ++length) 
     {
	 lfversion[length] = configlfversion[i];
     }
#endif

     for (i = 0; lapackversion[i] != '\0'; ++i, ++length) 
     {
	 lfversion[length] = lapackversion[i];
     }

     lfversion[length] = '\0';

     lflibversion.once = 1;

     return lflibversion.version;
}
