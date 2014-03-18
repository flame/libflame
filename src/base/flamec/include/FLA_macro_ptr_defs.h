
#include "FLA_type_defs.h"

// --- Pointer-accessing FLAME macro definitions ------------------------------------

#define FLA_CONSTANT_I_OFFSET  0
#define FLA_CONSTANT_S_OFFSET  ( sizeof(double) )
#define FLA_CONSTANT_D_OFFSET  ( sizeof(double) + sizeof(double) )
#define FLA_CONSTANT_C_OFFSET  ( sizeof(double) + sizeof(double) + sizeof(double) )
#define FLA_CONSTANT_Z_OFFSET  ( sizeof(double) + sizeof(double) + sizeof(double) + sizeof( scomplex ) )
#define FLA_CONSTANT_SIZE      ( sizeof(double) + sizeof(double) + sizeof(double) + sizeof( scomplex ) + sizeof( dcomplex ) )

#define FLA_INT_PTR( x ) \
  ( ((x).base)->datatype == FLA_CONSTANT ? \
    ( ( int * )      ( ( ( char * )     ((x).base)->buffer ) + FLA_CONSTANT_I_OFFSET             ) ) : \
                     ( ( ( int * )      ((x).base)->buffer ) + ( size_t ) (x).offn * ((x).base)->cs + \
                                                               ( size_t ) (x).offm * ((x).base)->rs ) )

#define FLA_FLOAT_PTR( x ) \
  ( ((x).base)->datatype == FLA_CONSTANT ? \
    ( ( float * )    ( ( ( char * )     ((x).base)->buffer ) + FLA_CONSTANT_S_OFFSET             ) ) : \
                     ( ( ( float * )    ((x).base)->buffer ) + ( size_t ) (x).offn * ((x).base)->cs + \
                                                               ( size_t ) (x).offm * ((x).base)->rs ) )

#define FLA_DOUBLE_PTR( x ) \
  ( ((x).base)->datatype == FLA_CONSTANT ? \
    ( ( double * )   ( ( ( char * )     ((x).base)->buffer ) + FLA_CONSTANT_D_OFFSET             ) ) : \
                     ( ( ( double * )   ((x).base)->buffer ) + ( size_t ) (x).offn * ((x).base)->cs + \
                                                               ( size_t ) (x).offm * ((x).base)->rs ) )

#define FLA_COMPLEX_PTR( x ) \
  ( ((x).base)->datatype == FLA_CONSTANT ? \
    ( ( scomplex * ) ( ( ( char * )     ((x).base)->buffer ) + FLA_CONSTANT_C_OFFSET             ) ) : \
                     ( ( ( scomplex * ) ((x).base)->buffer ) + ( size_t ) (x).offn * ((x).base)->cs + \
                                                               ( size_t ) (x).offm * ((x).base)->rs ) )

#define FLA_DOUBLE_COMPLEX_PTR( x ) \
  ( ((x).base)->datatype == FLA_CONSTANT ? \
    ( ( dcomplex * ) ( ( ( char * )     ((x).base)->buffer ) + FLA_CONSTANT_Z_OFFSET             ) ) : \
                     ( ( ( dcomplex * ) ((x).base)->buffer ) + ( size_t ) (x).offn * ((x).base)->cs + \
                                                               ( size_t ) (x).offm * ((x).base)->rs ) )

