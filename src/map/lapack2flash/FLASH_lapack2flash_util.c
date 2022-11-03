/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_LAPACK2FLASH

// Init
#define F77_fla_init                        F77_FUNC( fla_init, FLA_INIT )
#define F77_fla_finalize                    F77_FUNC( fla_finalize , FLA_FINALIZE )
#define F77_fla_initialized                 F77_FUNC( fla_initialized, FLA_INITIALIZED )
#define F77_fla_memory_leak_counter_status  F77_FUNC( fla_memory_leak_counter_status, FLA_MEMORY_LEAK_COUNTER_STATUS )
#define F77_fla_memory_leak_counter_set     F77_FUNC( fla_memory_leak_counter_set, FLA_MEMORY_LEAK_COUNTER_SET )
#define F77_fla_obj_show                    F77_FUNC( fla_obj_show, FLA_OBJ_SHOW )

static dim_t        flash_blocksize = 1024;
static unsigned int flash_n_threads = 1;
static dim_t        flash_depth = 1;

void F77_fla_init()
{
    FLA_Init();
}
void F77_fla_finalize()
{
    FLA_Finalize();
}
void F77_fla_initialized( int* ok )
{
    *ok = (FLA_Initialized() ? 1 : 0);
}
void F77_fla_memory_leak_counter_status( int* stat )
{
    *stat = (FLA_Memory_leak_counter_status() ? 1 : 0);
}
void F77_fla_memory_leak_counter_set( int* stat )
{
    FLA_Memory_leak_counter_set( (*stat ? TRUE : FALSE) );
}
void F77_fla_obj_show( char* prefix, int* m, int* n, void* buffer, int* ldim )
{
     FLA_Error    init_result;
     FLA_Datatype datatype;
     FLA_Obj      A;

     switch( *prefix ) {
     case 'i':
     case 'I': datatype = FLA_INT;            break;
     case 's':
     case 'S': datatype = FLA_FLOAT;          break;
     case 'd':
     case 'D': datatype = FLA_DOUBLE;         break;
     case 'c':
     case 'C': datatype = FLA_COMPLEX;        break;
     case 'z':
     case 'Z': datatype = FLA_DOUBLE_COMPLEX; break;
     default:
       fprintf(stderr, "Invalid prefix %c, where i,s,d,c,z are allowed.\n", *prefix);
       FLA_Abort();
     }

     FLA_Init_safe( &init_result );

     FLA_Obj_create_without_buffer( datatype, *m, *n, &A );
     FLA_Obj_attach_buffer( buffer, 1, *ldim, &A );
     FLA_Obj_fshow( stdout, 
                    "= F77_FLA_OBJ_SHOW =", A, "% 6.4e", 
                    "=-=-=-=-=-=-=-=-=-=-\n");
     FLA_Obj_free_without_buffer( &A );

     FLA_Finalize_safe( init_result );
}

// Transform tau.
int FLAME_invert_stau( FLA_Obj t )
{
    dim_t  m    = FLA_Obj_vector_dim( t );
    dim_t  inc  = FLA_Obj_vector_inc( t );
    float* buff = FLA_Obj_buffer_at_view( t );
    float  one  = 1.0F;
    float  zero = 0.0F;
    float* chi;
    int    i;

    for ( i = 0; i < m; ++i )
    {
        chi = buff + i*inc;
        if ( *chi != zero )
            *chi = ( one / *chi );
    }
    return 0;
}
int FLAME_invert_dtau( FLA_Obj t )
{
    dim_t   m    = FLA_Obj_vector_dim( t );
    dim_t   inc  = FLA_Obj_vector_inc( t );
    double* buff = FLA_Obj_buffer_at_view( t );
    double  one  = 1.0;
    double  zero = 0.0;
    double* chi;
    int     i;

    for ( i = 0; i < m; ++i )
    {
        chi = buff + i*inc;
        if ( *chi != zero )
            *chi = ( one / *chi );
    }
    return 0;
}
int FLAME_invert_ctau( FLA_Obj t )
{
    dim_t     m    = FLA_Obj_vector_dim( t );
    dim_t     inc  = FLA_Obj_vector_inc( t );
    scomplex* buff = FLA_Obj_buffer_at_view( t );
    float     one  = 1.0F;
    float     conjsign = one; // if conjugate -one;
    float     zero = 0.0F;
    float     temp, s, xr_s, xi_s;
    scomplex* chi;
    int       i;

    for ( i = 0; i < m; ++i )
    {
        chi  = buff + i*inc;
        s    = bl1_fmaxabs( chi->real, chi->imag );
        if ( s != zero )
        {
            xr_s = chi->real / s;
            xi_s = chi->imag / s;
            temp = xr_s * chi->real + xi_s * chi->imag;
            chi->real =            xr_s / temp;
            chi->imag = conjsign * xi_s / temp;
        }
    }
    return 0;
}
int FLAME_invert_ztau( FLA_Obj t )
{
    dim_t     m    = FLA_Obj_vector_dim( t );
    dim_t     inc  = FLA_Obj_vector_inc( t );
    dcomplex* buff = FLA_Obj_buffer_at_view( t );
    double    one  = 1.0;
    double    conjsign = one; // if conjugate -one;
    double    zero = 0.0;
    double    temp, s, xr_s, xi_s;
    dcomplex* chi;
    int       i;

    for ( i = 0; i < m; ++i )
    {
        chi  = buff + i*inc;
        s    = bl1_fmaxabs( chi->real, chi->imag );
        if ( s != zero )
        {
            xr_s = chi->real / s;
            xi_s = chi->imag / s;
            temp = xr_s * chi->real + xi_s * chi->imag;
            chi->real =            xr_s / temp;
            chi->imag = conjsign * xi_s / temp;
        }
    }
    return 0;
}

/*
int FLAME_QR_piv_preorder( FLA_Obj A, int *jpiv_lapack, int *jpiv_fla )
{
    FLA_Datatype datatype = FLA_Obj_datatype( A );
    dim_t m_A  = FLA_Obj_length( A );
    dim_t n_A  = FLA_Obj_width( A );
    dim_t cs_A = FLA_Obj_col_stride( A );
    dim_t rs_A = FLA_Obj_row_stride( A );
    int j, head;

    for ( j=0; j<n_A; ++j )
        jpiv_lapack[j] = (j+1);
    return 0;

    switch ( datatype )
    {
    case FLA_FLOAT:
    {
        float* buff_A = FLA_FLOAT_PTR( A );

        head = 0;
        for ( j=0; j<n_A; ++j )
        {
            // For nonzero lapack pivots...
            if ( jpiv_lapack[j] != 0 )
            {
                // If the pivot is not same as the head...
                if ( j != head )
                {
                    // Reorder A by swapping columns.
                    bl1_sswapmt( BLIS1_NO_TRANSPOSE, m_A, 1,
                                 &buff_A[j*   cs_A], rs_A, cs_A,
                                 &buff_A[head*cs_A], rs_A, cs_A );
                    // Swap pivot indices
                    jpiv_lapack[j]    = jpiv_lapack[head];
                    jpiv_lapack[head] = (j+1);

                    // Set FLAME pivots minus one;
                    // this will ingore column pivoting in performing FLA_QR_UT_piv.
                    jpiv_fla[j]       = -1;
                }
                else
                {
                    jpiv_lapack[j] = (j+1);
                }
                // Move forward the head.
                ++head;
            }
            else
            {
                jpiv_lapack[j] = (j+1);
            }
        }
        break;
    }
    case FLA_DOUBLE:
    {
        //double* buff_A = FLA_DOUBLE_PTR( A );
        break;
    }
    case FLA_COMPLEX:
    {
        //scomplex* buff_A = FLA_COMPLEX_PTR( A );
        break;
    }
    case FLA_DOUBLE_COMPLEX:
    {
        //dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );
        break;
    }
    }
    return 0;
}
*/

/*---------------------Get/Set Helper---------------------*/

FLA_Error FLASH_set_preferred_blocksize( dim_t blocksize )
{
    flash_blocksize = blocksize;

    return FLA_SUCCESS;
}

dim_t FLASH_get_preferred_blocksize( void )
{
    return flash_blocksize;
}

FLA_Error FLASH_set_n_preferred_threads( unsigned int threads )
{
    flash_n_threads = threads;

    #ifdef FLA_ENABLE_HIP
        int device_count = 1;
        FLASH_Queue_available_devices_hip( &device_count );

        // Pass number of threads unless it's more than available devices
        FLASH_Queue_set_num_threads( min( threads, device_count) );
    #else
        FLASH_Queue_set_num_threads( threads );
    #endif

    return FLA_SUCCESS;
}

unsigned int FLASH_get_n_preferred_threads( void )
{
    // Returning the number of preferred threads, not actual number.
    return flash_n_threads;
}

FLA_Error FLASH_set_depth( dim_t depth )
{
    flash_depth = depth;

    return FLA_SUCCESS;
}

dim_t FLASH_get_depth( void )
{
    return flash_depth;
}

#endif
