/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Apply_pivots_macro_external( FLA_Side side, FLA_Trans trans, FLA_Obj p, FLA_Obj A )
{
   int          i, j;
   int          ipiv;
   int*         buf_p    = ( int* ) FLA_Obj_buffer_at_view( p );
   FLA_Obj*     blocks   = FLASH_OBJ_PTR_AT( A );
   int          m_blocks = FLA_Obj_length( A );
   int          m_A      = FLA_Obj_length( *blocks );
   int          n_A      = FLA_Obj_width( *blocks );
   FLA_Datatype datatype = FLA_Obj_datatype( A );

#ifdef FLA_ENABLE_WINDOWS_BUILD
   int* m  = ( int* ) _alloca( m_blocks * sizeof( int ) );
   int* cs = ( int* ) _alloca( m_blocks * sizeof( int ) );
#else
   int* m  = ( int* ) malloc( m_blocks * sizeof( int ) );
   int* cs = ( int* ) malloc( m_blocks * sizeof( int ) );
   //int m[m_blocks];
   //int cs[m_blocks];
#endif

   if ( side != FLA_LEFT || trans != FLA_NO_TRANSPOSE )
      FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );

   switch ( datatype )
   {
      case FLA_FLOAT:
      {
#ifdef FLA_ENABLE_WINDOWS_BUILD
         float** buffer = ( float** ) _alloca( m_blocks * sizeof( float* ) );
#else
         float** buffer = ( float** ) malloc( m_blocks * sizeof( float* ) );
         //float*  buffer[m_blocks];
#endif
         for ( i = 0; i < m_blocks; i++ )
         {
            buffer[i] = ( float* ) FLA_Obj_buffer_at_view( blocks[i] );
            
            m[i] = FLA_Obj_length( blocks[i] );
            cs[i] = FLA_Obj_col_stride( blocks[i] );
         }
         
         for ( j = 0; j < m_A; j++ )
         {
            ipiv = buf_p[j] + j;
            
            if ( ipiv != j )
            {
               i = 0;
               
               while ( ipiv >= m[i] )
               {
                  ipiv = ipiv - m[i];
                  i++;
               }

               bl1_sswapv( n_A,
                           buffer[0] + j,    cs[0],
                           buffer[i] + ipiv, cs[i] );
            }
         }
#ifdef FLA_ENABLE_WINDOWS_BUILD
#else
         free( buffer ); 
#endif
         break;
      }
      case FLA_DOUBLE:
      {
#ifdef FLA_ENABLE_WINDOWS_BUILD
         double** buffer = ( double** ) _alloca( m_blocks * sizeof( double* ) );
#else
         double** buffer = ( double** ) malloc( m_blocks * sizeof( double* ) );
         //double*  buffer[m_blocks];
#endif
         for ( i = 0; i < m_blocks; i++ )
         {
            buffer[i] = ( double* ) FLA_Obj_buffer_at_view( blocks[i] );
            
            m[i] = FLA_Obj_length( blocks[i] );
            cs[i] = FLA_Obj_col_stride( blocks[i] );
         }
         
         for ( j = 0; j < m_A; j++ )
         {
            ipiv = buf_p[j] + j;
            
            if ( ipiv != j )
            {
               i = 0;
               
               while ( ipiv >= m[i] )
               {
                  ipiv = ipiv - m[i];
                  i++;
               }

               bl1_dswapv( n_A,
                           buffer[0] + j,    cs[0],
                           buffer[i] + ipiv, cs[i] );
            }
         }
#ifdef FLA_ENABLE_WINDOWS_BUILD
#else
         free( buffer ); 
#endif
         break;
      }
      case FLA_COMPLEX:
      {
#ifdef FLA_ENABLE_WINDOWS_BUILD
         scomplex** buffer = ( scomplex** ) _alloca( m_blocks * sizeof( scomplex* ) );
#else
         scomplex** buffer = ( scomplex** ) malloc( m_blocks * sizeof( scomplex* ) );
         //scomplex*  buffer[m_blocks];
#endif
         for ( i = 0; i < m_blocks; i++ )
         {
            buffer[i] = ( scomplex* ) FLA_Obj_buffer_at_view( blocks[i] );
            
            m[i] = FLA_Obj_length( blocks[i] );
            cs[i] = FLA_Obj_col_stride( blocks[i] );
         }
         
         for ( j = 0; j < m_A; j++ )
         {
            ipiv = buf_p[j] + j;
            
            if ( ipiv != j )
            {
               i = 0;
               
               while ( ipiv >= m[i] )
               {
                  ipiv = ipiv - m[i];
                  i++;
               }

               bl1_cswapv( n_A,
                           buffer[0] + j,    cs[0],
                           buffer[i] + ipiv, cs[i] );
            }
         }
#ifdef FLA_ENABLE_WINDOWS_BUILD
#else
         free( buffer ); 
#endif
         break;
      }
      case FLA_DOUBLE_COMPLEX:
      {
#ifdef FLA_ENABLE_WINDOWS_BUILD
         dcomplex** buffer = ( dcomplex** ) _alloca( m_blocks * sizeof( dcomplex* ) );
#else
         dcomplex** buffer = ( dcomplex** ) malloc( m_blocks * sizeof( dcomplex* ) );
         //dcomplex*  buffer[m_blocks];
#endif
         for ( i = 0; i < m_blocks; i++ )
         {
            buffer[i] = ( dcomplex* ) FLA_Obj_buffer_at_view( blocks[i] );
            
            m[i] = FLA_Obj_length( blocks[i] );
            cs[i] = FLA_Obj_col_stride( blocks[i] );
         }
         
         for ( j = 0; j < m_A; j++ )
         {
            ipiv = buf_p[j] + j;
            
            if ( ipiv != j )
            {
               i = 0;
               
               while ( ipiv >= m[i] )
               {
                  ipiv = ipiv - m[i];
                  i++;
               }

               bl1_zswapv( n_A,
                           buffer[0] + j,    cs[0],
                           buffer[i] + ipiv, cs[i] );
            }
         }
#ifdef FLA_ENABLE_WINDOWS_BUILD
#else
         free( buffer ); 
#endif
         break;
      }
   }

#ifdef FLA_ENABLE_WINDOWS_BUILD
#else
   free( m );
   free( cs );
#endif

   return FLA_SUCCESS;
}

