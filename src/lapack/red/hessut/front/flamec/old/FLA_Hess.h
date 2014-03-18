
#include "FLAME.h"

FLA_Error FLA_Hess( FLA_Obj A, FLA_Obj t, int ilo, int ihi );
FLA_Error FLA_Hess_internal( FLA_Obj A, FLA_Obj t, int ilo, int ihi, fla_hess_t* cntl );
FLA_Error FLA_Hess_blk_var1( FLA_Obj A, FLA_Obj t, int ilo, int ihi, fla_hess_t* cntl );
FLA_Error FLA_Hess_blk_var2( FLA_Obj A, FLA_Obj t, int ilo, int ihi, fla_hess_t* cntl );
FLA_Error FLA_Hess_step_unb_var1( FLA_Obj A, FLA_Obj t, FLA_Obj B, FLA_Obj C, int nb_alg );
FLA_Error FLA_Hess_step_unb_var2( FLA_Obj A, FLA_Obj t, FLA_Obj B, FLA_Obj C, int nb_alg );

void MyFLA_extract_column_block( FLA_Obj A, int ini, int n, FLA_Obj * B );
void MyFLA_extract_row_block( FLA_Obj A, int ini, int n, FLA_Obj * B );
void MyFLA_int_swap( int * a, int * b );
void MyFLA_v_int_swap( int n, int * a, int * b );

