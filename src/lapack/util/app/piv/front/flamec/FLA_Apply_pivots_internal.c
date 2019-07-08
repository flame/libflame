/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_THREAD_SAFE_INTERFACES
FLA_Error FLA_Apply_pivots_internal_ts( FLA_cntl_init_s *FLA_cntl_init_i, FLA_Side side, FLA_Trans trans, FLA_Obj p, FLA_Obj A, fla_appiv_t* cntl )
{
   FLA_Error r_val = FLA_SUCCESS;

   if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
        FLA_Cntl_variant( cntl ) == FLA_SUBPROBLEM )
   {
      if ( FLASH_Queue_get_enabled_ts(FLA_cntl_init_i) )
      {
         // Enqueue
         ENQUEUE_FLASH_Apply_pivots_macro( side, trans, *FLASH_OBJ_PTR_AT_TS( FLA_cntl_init_i, p ), A, cntl );
      }
      else
      {
         // Execute leaf
         r_val = FLA_Apply_pivots_macro_task_ts( FLA_cntl_init_i, side, trans, *FLASH_OBJ_PTR_AT_TS( FLA_cntl_init_i, p ), A, cntl );
      }
   }
   else
   {
      // Parameter combinations
      if ( trans == FLA_NO_TRANSPOSE )
      {
         if      ( side == FLA_LEFT )
         {
            r_val = FLA_Apply_pivots_ln_ts( FLA_cntl_init_i, p, A, cntl );
         }
         else if ( side == FLA_RIGHT )
         {
            r_val = FLA_Apply_pivots_rn_ts( FLA_cntl_init_i, p, A, cntl );
         }
      }
      else if ( trans == FLA_TRANSPOSE )
      {
         if      ( side == FLA_LEFT )
         {
            r_val = FLA_Apply_pivots_lt_ts( FLA_cntl_init_i, p, A, cntl );
         }
         else if ( side == FLA_RIGHT )
         {
            r_val = FLA_Apply_pivots_rt_ts( FLA_cntl_init_i, p, A, cntl );
         }
      }
   }   

   return r_val;
}
#endif

FLA_Error FLA_Apply_pivots_internal( FLA_Side side, FLA_Trans trans, FLA_Obj p, FLA_Obj A, fla_appiv_t* cntl )
{
   FLA_Error r_val = FLA_SUCCESS;

   if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
        FLA_Cntl_variant( cntl ) == FLA_SUBPROBLEM )
   {
      if ( FLASH_Queue_get_enabled( ) )
      {
         // Enqueue
         ENQUEUE_FLASH_Apply_pivots_macro( side, trans, *FLASH_OBJ_PTR_AT( p ), A, cntl );
      }
      else
      {
         // Execute leaf
         r_val = FLA_Apply_pivots_macro_task( side, trans, *FLASH_OBJ_PTR_AT( p ), A, cntl );
      }
   }
   else
   {
      // Parameter combinations
      if ( trans == FLA_NO_TRANSPOSE )
      {
         if      ( side == FLA_LEFT )
         {
            r_val = FLA_Apply_pivots_ln( p, A, cntl );
         }
         else if ( side == FLA_RIGHT )
         {
            r_val = FLA_Apply_pivots_rn( p, A, cntl );
         }
      }
      else if ( trans == FLA_TRANSPOSE )
      {
         if      ( side == FLA_LEFT )
         {
            r_val = FLA_Apply_pivots_lt( p, A, cntl );
         }
         else if ( side == FLA_RIGHT )
         {
            r_val = FLA_Apply_pivots_rt( p, A, cntl );
         }
      }
   }   

   return r_val;
}
