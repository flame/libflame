/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/
#include "FLAME.h"

int main( int argc, char** argv ) {
  FLA_Error    init_result; 

  FLA_Init_safe( &init_result );          

  FLA_Obj_fshow(stdout, "FLA_EPSILON",                FLA_EPSILON,                "%6.4e", "--");
  FLA_Obj_fshow(stdout, "FLA_SAFE_MIN",               FLA_SAFE_MIN,               "%6.4e", "--");
  FLA_Obj_fshow(stdout, "FLA_SAFE_INV_MIN",           FLA_SAFE_INV_MIN,           "%6.4e", "--");
  FLA_Obj_fshow(stdout, "FLA_UNDERFLOW_THRES",        FLA_UNDERFLOW_THRES,        "%6.4e", "--");
  FLA_Obj_fshow(stdout, "FLA_OVERFLOW_THRES",         FLA_OVERFLOW_THRES,         "%6.4e", "--");
  FLA_Obj_fshow(stdout, "FLA_UNDERFLOW_SQUARE_THRES", FLA_UNDERFLOW_SQUARE_THRES, "%6.4e", "--");
  FLA_Obj_fshow(stdout, "FLA_OVERFLOW_SQUARE_THRES",  FLA_OVERFLOW_SQUARE_THRES,  "%6.4e", "--");

  FLA_Finalize_safe( init_result );     
}
