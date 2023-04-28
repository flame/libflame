/*

    Copyright (C) 2023, Advanced Micro Devices, Inc.

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#ifndef FLASH_GET_SET_CONTROLS_H
#define FLASH_GET_SET_CONTROLS_H

#ifdef FLA_ENABLE_LAPACK2FLASH // Start lapack2flash

FLA_Error      FLASH_set_preferred_blocksize( dim_t blocksize );
dim_t          FLASH_get_preferred_blocksize( void );
FLA_Error      FLASH_set_n_preferred_threads( unsigned int threads );
unsigned int   FLASH_get_n_preferred_threads( void );
FLA_Error      FLASH_set_depth( dim_t depth );
dim_t          FLASH_get_depth( void );
FLA_Error      FLASH_set_tile_offload( unsigned int tiles );
unsigned int   FLASH_get_tile_offload( void );

#endif // End lapack2flash

#endif // End header