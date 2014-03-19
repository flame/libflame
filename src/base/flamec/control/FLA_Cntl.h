/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

// Shared definitions

typedef struct FLA_Blocksize_s
{
	dim_t s;
	dim_t d;
	dim_t c;
	dim_t z;
} fla_blocksize_t;

#define FLA_SUBPROBLEM                  0
#define FLA_UNBLOCKED_EXTERN           10
#define FLA_BLOCKED_EXTERN             13

#define FLA_UNB_VAR_OFFSET             40
#define FLA_OPT_VAR_OFFSET             80
#define FLA_BLK_VAR_OFFSET            120
#define FLA_BLF_VAR_OFFSET            160

#define FLA_UNBLOCKED_VARIANT1        (FLA_UNB_VAR_OFFSET+1)
#define FLA_UNBLOCKED_VARIANT2        (FLA_UNB_VAR_OFFSET+2)
#define FLA_UNBLOCKED_VARIANT3        (FLA_UNB_VAR_OFFSET+3)
#define FLA_UNBLOCKED_VARIANT4        (FLA_UNB_VAR_OFFSET+4)
#define FLA_UNBLOCKED_VARIANT5        (FLA_UNB_VAR_OFFSET+5)
#define FLA_UNBLOCKED_VARIANT6        (FLA_UNB_VAR_OFFSET+6)
#define FLA_UNBLOCKED_VARIANT7        (FLA_UNB_VAR_OFFSET+7)
#define FLA_UNBLOCKED_VARIANT8        (FLA_UNB_VAR_OFFSET+8)
#define FLA_UNBLOCKED_VARIANT9        (FLA_UNB_VAR_OFFSET+9)
#define FLA_UNBLOCKED_VARIANT10       (FLA_UNB_VAR_OFFSET+10)

#define FLA_UNB_OPT_VARIANT1          (FLA_OPT_VAR_OFFSET+1)
#define FLA_UNB_OPT_VARIANT2          (FLA_OPT_VAR_OFFSET+2)
#define FLA_UNB_OPT_VARIANT3          (FLA_OPT_VAR_OFFSET+3)
#define FLA_UNB_OPT_VARIANT4          (FLA_OPT_VAR_OFFSET+4)
#define FLA_UNB_OPT_VARIANT5          (FLA_OPT_VAR_OFFSET+5)
#define FLA_UNB_OPT_VARIANT6          (FLA_OPT_VAR_OFFSET+6)
#define FLA_UNB_OPT_VARIANT7          (FLA_OPT_VAR_OFFSET+7)
#define FLA_UNB_OPT_VARIANT8          (FLA_OPT_VAR_OFFSET+8)
#define FLA_UNB_OPT_VARIANT9          (FLA_OPT_VAR_OFFSET+9)
#define FLA_UNB_OPT_VARIANT10         (FLA_OPT_VAR_OFFSET+10)

#define FLA_BLOCKED_VARIANT1          (FLA_BLK_VAR_OFFSET+1)
#define FLA_BLOCKED_VARIANT2          (FLA_BLK_VAR_OFFSET+2)
#define FLA_BLOCKED_VARIANT3          (FLA_BLK_VAR_OFFSET+3)
#define FLA_BLOCKED_VARIANT4          (FLA_BLK_VAR_OFFSET+4)
#define FLA_BLOCKED_VARIANT5          (FLA_BLK_VAR_OFFSET+5)
#define FLA_BLOCKED_VARIANT6          (FLA_BLK_VAR_OFFSET+6)
#define FLA_BLOCKED_VARIANT7          (FLA_BLK_VAR_OFFSET+7)
#define FLA_BLOCKED_VARIANT8          (FLA_BLK_VAR_OFFSET+8)
#define FLA_BLOCKED_VARIANT9          (FLA_BLK_VAR_OFFSET+9)
#define FLA_BLOCKED_VARIANT10         (FLA_BLK_VAR_OFFSET+10)
#define FLA_BLOCKED_VARIANT11         (FLA_BLK_VAR_OFFSET+11)
#define FLA_BLOCKED_VARIANT12         (FLA_BLK_VAR_OFFSET+12)
#define FLA_BLOCKED_VARIANT13         (FLA_BLK_VAR_OFFSET+13)
#define FLA_BLOCKED_VARIANT14         (FLA_BLK_VAR_OFFSET+14)
#define FLA_BLOCKED_VARIANT15         (FLA_BLK_VAR_OFFSET+15)
#define FLA_BLOCKED_VARIANT16         (FLA_BLK_VAR_OFFSET+16)
#define FLA_BLOCKED_VARIANT17         (FLA_BLK_VAR_OFFSET+17)
#define FLA_BLOCKED_VARIANT18         (FLA_BLK_VAR_OFFSET+18)
#define FLA_BLOCKED_VARIANT19         (FLA_BLK_VAR_OFFSET+19)
#define FLA_BLOCKED_VARIANT20         (FLA_BLK_VAR_OFFSET+20)

#define FLA_BLK_FUS_VARIANT1          (FLA_BLF_VAR_OFFSET+1)
#define FLA_BLK_FUS_VARIANT2          (FLA_BLF_VAR_OFFSET+2)
#define FLA_BLK_FUS_VARIANT3          (FLA_BLF_VAR_OFFSET+3)
#define FLA_BLK_FUS_VARIANT4          (FLA_BLF_VAR_OFFSET+4)
#define FLA_BLK_FUS_VARIANT5          (FLA_BLF_VAR_OFFSET+5)
#define FLA_BLK_FUS_VARIANT6          (FLA_BLF_VAR_OFFSET+6)
#define FLA_BLK_FUS_VARIANT7          (FLA_BLF_VAR_OFFSET+7)
#define FLA_BLK_FUS_VARIANT8          (FLA_BLF_VAR_OFFSET+8)
#define FLA_BLK_FUS_VARIANT9          (FLA_BLF_VAR_OFFSET+9)
#define FLA_BLK_FUS_VARIANT10         (FLA_BLF_VAR_OFFSET+10)

#define FLA_Cntl_matrix_type( cntl )  cntl->matrix_type
#define FLA_Cntl_blocksize( cntl )    cntl->blocksize
#define FLA_Cntl_variant( cntl )      cntl->variant

void FLA_Cntl_obj_free( void* cntl );


// Include the control tree definitions for each class of operation.
#include "FLA_Cntl_blas1.h"
#include "FLA_Cntl_blas2.h"
#include "FLA_Cntl_blas3.h"
#include "FLA_Cntl_lapack.h"

