/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#ifndef FLASH_MACRO_DEFS_H
#define FLASH_MACRO_DEFS_H

#define FLASH_OBJ_PTR_AT( A )  ( (FLA_Obj *) FLA_Obj_buffer_at_view( A ) )

#define FLA_FLAT_TO_HIER 4000
#define FLA_HIER_TO_FLAT 4001

#endif
