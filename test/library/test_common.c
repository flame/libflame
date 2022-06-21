/*
    Copyright (c) 2022 Advanced Micro Devices, Inc.  All rights reserved.
*/

#include "test_common.h"

// Global variables
integer i_zero = 0, i_one = 1, i_n_one = -1;
float s_zero = 0, s_one = 1, s_n_one = -1;
double d_zero = 0, d_one = 1, d_n_one = -1;
scomplex c_zero = {0,0}, c_one = {1,0}, c_n_one = {-1,0};
dcomplex z_zero = {0,0}, z_one = {1,0}, z_n_one = {-1,0};

/* create vector of given datatype*/
void create_vector(integer datatype, void **A, integer M)
{
    *A = NULL;

    switch(datatype)
    {
        case INTEGER:
        {
            *A = (integer *)malloc(M * sizeof(integer));
            break;
        }

        case FLOAT:
        {
            *A = (float *)malloc(M * sizeof(float));
            break;
        }

        case DOUBLE:
        {
            *A = (double *)malloc(M * sizeof(double));
            break;
        }

        case COMPLEX:
        {
            *A = (scomplex *)malloc(M * sizeof(scomplex));
            break;
        }

        case DOUBLE_COMPLEX:
        {
            *A = (dcomplex *)malloc(M * sizeof(dcomplex));
            break;
        }
    }

    if(*A == NULL)
    {
        fprintf( stderr, "malloc() returned NULL pointer\n");
        abort();
    }

    return;
}


void create_realtype_vector(integer datatype, void **A, integer M)
{
    *A = NULL;

    if(datatype == FLOAT || datatype == COMPLEX)
        *A = (float *)malloc(M * sizeof(float));
    else
        *A = (double *)malloc(M * sizeof(double));

    if(*A == NULL)
    {
        fprintf( stderr, "malloc() returned NULL pointer\n");
        abort();
    }

    return;
}

/* free vector */
void free_vector(void *A)
{
    if(!A)
        return;

    free(A);
}

/* Initialize vector with random values */
void rand_vector(integer datatype, void *A, integer M, integer LDA)
{
    integer i;

    switch( datatype )
    {
        case FLOAT:
        {
            for( i = 0; i < M; i++ )
            {
                ((float *)A)[i * LDA] = SRAND();
            }
            break;
        }
        case DOUBLE:
        {
            for( i = 0; i < M; i++ )
            {
                ((double *)A)[i * LDA] = DRAND();
            }
            break;
        }
        case COMPLEX:
        {
            for( i = 0; i < M; i++ )
            {
                ((scomplex *)A)[i * LDA].real = SRAND();
                ((scomplex *)A)[i * LDA].imag = SRAND();
            }
            break;
        }
        case DOUBLE_COMPLEX:
        {
            for( i = 0; i < M; i++ )
            {
                ((dcomplex *)A)[i * LDA].real = DRAND();
                ((dcomplex *)A)[i * LDA].imag = DRAND();
            }
            break;
        }
    }

    return;
}


/* Copy a vector */
void copy_vector(integer datatype, integer M, void *A, integer LDA, void *B, integer LDB)
{
    switch( datatype )
    {
        case INTEGER:
        {
            integer i;

            for( i = 0; i < M; i++ )
            {
                ((integer *)B)[ i * LDB ] = ((integer *)A)[ i * LDA ];
            }
            break;
        }
        case FLOAT:
        {
            scopy_(&M, A, &LDA, B, &LDB);
            break;
        }
        case DOUBLE:
        {
            dcopy_(&M, A, &LDA, B, &LDB);
            break;
        }
        case COMPLEX:
        {
            ccopy_(&M, A, &LDA, B, &LDB);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            zcopy_(&M, A, &LDA, B, &LDB);
            break;
        }
    }

    return;
}

void copy_realtype_vector(integer datatype, integer M, void *A, integer LDA, void *B, integer LDB)
{
    if(datatype == FLOAT || datatype == COMPLEX)
        scopy_(&M, A, &LDA, B, &LDB);
    else
        dcopy_(&M, A, &LDA, B, &LDB);

    return;
}


/* create matrix of given datatype*/
void create_matrix(integer datatype, void **A, integer M, integer N)
{
    *A = NULL;

    switch(datatype)
    {
        case INTEGER:
        {
            *A = (integer *)malloc(M * N * sizeof(integer));
            break;
        }

        case FLOAT:
        {
            *A = (float *)malloc(M * N * sizeof(float));
            break;
        }

        case DOUBLE:
        {
            *A = (double *)malloc(M * N * sizeof(double));
            break;
        }

        case COMPLEX:
        {
            *A = (scomplex *)malloc(M * N * sizeof(scomplex));
            break;
        }

        case DOUBLE_COMPLEX:
        {
            *A = (dcomplex *)malloc(M * N * sizeof(dcomplex));
            break;
        }
    }

    if(*A == NULL)
    {
        fprintf( stderr, "malloc() returned NULL pointer\n");
        abort();
    }

    return;
}


void create_realtype_matrix(integer datatype, void **A, integer M, integer N)
{
    *A = NULL;

    if(datatype == FLOAT || datatype == COMPLEX)
        *A = (float *)malloc(M * N * sizeof(float));
    else
        *A = (double *)malloc(M * N * sizeof(double));

    if(*A == NULL)
    {
        fprintf( stderr, "malloc() returned NULL pointer\n");
        abort();
    }

    return;
}

void* get_m_ptr(integer datatype, void *A, integer M, integer N, integer LDA)
{
    void *mat = NULL;

    switch(datatype)
    {
        case FLOAT:
        {
            mat = ((float *)A) + M + N * LDA;
            break;
        }
        case DOUBLE:
        {
            mat = ((double *)A) + M + N * LDA;
            break;
        }
        case COMPLEX:
        {
            mat = ((scomplex *)A) + M + N * LDA;
            break;
        }
        case DOUBLE_COMPLEX:
        {
            mat = ((dcomplex *)A) + M + N * LDA;
            break;
        }
    }

    return mat;
}



/* free matrix */
void free_matrix(void *A)
{
    if(!A)
        return;

    free(A);
}


/* Initialize matrix with random values */
void rand_matrix(integer datatype, void *A, integer M, integer N, integer LDA)
{
    integer i, j;

    switch( datatype )
    {
        case FLOAT:
        {
            for( i = 0; i < N; i++ )
            {
                for( j = 0; j < M; j++ )
                {
                    ((float *)A)[i * LDA + j] = SRAND();
                }
            }
            break;
        }
        case DOUBLE:
        {
            for( i = 0; i < N; i++ )
            {
                for( j = 0; j < M; j++ )
                {
                    ((double *)A)[i * LDA + j] = DRAND();
                }
            }
            break;
        }
        case COMPLEX:
        {
            for( i = 0; i < N; i++ )
            {
                for( j = 0; j < M; j++ )
                {
                    ((scomplex *)A)[i * LDA + j].real = SRAND();
                    ((scomplex *)A)[i * LDA + j].imag = SRAND();
                }
            }
            break;
        }
        case DOUBLE_COMPLEX:
        {
            for( i = 0; i < N; i++ )
            {
                for( j = 0; j < M; j++ )
                {
                    ((dcomplex *)A)[i * LDA + j].real = DRAND();
                    ((dcomplex *)A)[i * LDA + j].imag = DRAND();
                }
            }
            break;
        }
    }

    return;
}

/* Initialize symmetric matrix with random values */
void rand_sym_matrix(integer datatype, void *A, integer M, integer N, integer LDA)
{
    integer i, j;

    switch( datatype )
    {
        case FLOAT:
        {
            for( i = 0; i < N; i++ )
            {
                for( j = i; j < M; j++ )
                {
                    ((float *)A)[i * LDA + j] = SRAND();
          ((float *)A)[j * LDA + i] = ((float *)A)[i * LDA + j];
                }
            }
            break;
        }
        case DOUBLE:
        {
            for( i = 0; i < N; i++ )
            {
                for( j = i; j < M; j++ )
                {
                    ((double *)A)[i * LDA + j] = DRAND();
          ((double *)A)[j * LDA + i] = ((double *)A)[i * LDA + j];
                }
            }
            break;
        }
        case COMPLEX:
        {
            for( i = 0; i < N; i++ )
            {
                for( j = i; j < M; j++ )
                {
                    ((scomplex *)A)[i * LDA + j].real = SRAND();
                    ((scomplex *)A)[i * LDA + j].imag = SRAND();
          ((scomplex *)A)[j * LDA + i].real = ((scomplex *)A)[i * LDA + j].real;
                    ((scomplex *)A)[j * LDA + i].imag = ((scomplex *)A)[i * LDA + j].imag;
                }
            }
            break;
        }
        case DOUBLE_COMPLEX:
        {
            for( i = 0; i < N; i++ )
            {
                for( j = i; j < M; j++ )
                {
                    ((dcomplex *)A)[i * LDA + j].real = DRAND();
                    ((dcomplex *)A)[i * LDA + j].imag = DRAND();
          ((dcomplex *)A)[j * LDA + i].real = ((dcomplex *)A)[i * LDA + j].real;
                    ((dcomplex *)A)[j * LDA + i].imag = ((dcomplex *)A)[i * LDA + j].imag;
                }
            }
            break;
        }
    }

    return;
}


/* Copy a matrix */
void copy_matrix(integer datatype, char *uplo, integer M, integer N, void *A, integer LDA, void *B, integer LDB)
{
    switch( datatype )
    {
        case INTEGER:
        {
            integer i, j;

            for( i = 0; i < N; i++ )
            {
                for( j = 0; j < M; j++ )
                {
                    ((integer *)B)[ i * LDB + j ] = ((integer *)A)[ i * LDA + j ];
                }
            }
            break;
        }
        case FLOAT:
        {
            slacpy_(uplo, &M, &N, A, &LDA, B, &LDB);
            break;
        }
        case DOUBLE:
        {
            dlacpy_(uplo, &M, &N, A, &LDA, B, &LDB);
            break;
        }
        case COMPLEX:
        {
            clacpy_(uplo, &M, &N, A, &LDA, B, &LDB);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            zlacpy_(uplo, &M, &N, A, &LDA, B, &LDB);
            break;
        }
    }

    return;
}


void copy_realtype_matrix(integer datatype, char *uplo, integer M, integer N, void *A, integer LDA, void *B, integer LDB)
{
    if(datatype == FLOAT || datatype == COMPLEX)
        slacpy_(uplo, &M, &N, A, &LDA, B, &LDB);
    else
        dlacpy_(uplo, &M, &N, A, &LDA, B, &LDB);

    return;
}



/* Initialize a matrix with zeros */
void reset_matrix(integer datatype, integer M, integer N, void *A, integer LDA)
{
    integer i, j;

    switch( datatype )
    {
        case INTEGER:
        {
            for( i = 0; i < N; i++ )
            {
                for( j = 0; j < M; j++ )
                {
                    ((integer *)A)[ i * LDA + j ] = 0.f;
                }
            }
            break;
        }

        case FLOAT:
        {
            slaset_("A", &M, &N, &s_zero, &s_zero, A, &LDA);
            break;
        }

        case DOUBLE:
        {
            dlaset_("A", &M, &N, &d_zero, &d_zero, A, &LDA);
            break;
        }

        case COMPLEX:
        {
            claset_("A", &M, &N, &c_zero, &c_zero, A, &LDA);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            zlaset_("A", &M, &N, &z_zero, &z_zero, A, &LDA);
            break;
        }
    }

    return;
}


/* Set a matrix to identity */
void set_identity_matrix(integer datatype, integer M, integer N, void *A, integer LDA)
{

    switch( datatype )
    {
        case FLOAT:
        {
            slaset_("A", &M, &N, &s_zero, &s_one, A, &LDA);
            break;
        }

        case DOUBLE:
        {
            dlaset_("A", &M, &N, &d_zero, &d_one, A, &LDA);
            break;
        }

        case COMPLEX:
        {
            claset_("A", &M, &N, &c_zero, &c_one, A, &LDA);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            zlaset_("A", &M, &N, &z_zero, &z_one, A, &LDA);
            break;
        }
    }

    return;
}

void z_div_t(dcomplex *cp, dcomplex *ap, dcomplex *bp)
{
    dcomplex a = *ap;
    dcomplex b = *bp;
    double temp;

    temp = b.real * b.real + b.imag * b.imag;
    if(!temp)
    {
        fprintf( stderr, "z_div_t : temp is zero. Abort\n");
        abort();
    }

    cp->real = ( a.real * b.real + a.imag * b.imag ) / temp;
    cp->imag = ( a.imag * b.real - a.real * b.imag ) / temp;
}


/* Division of complex types */
void c_div_t(scomplex *cp, scomplex *ap, scomplex *bp)
{
    scomplex a = *ap;
    scomplex b = *bp;
    float temp;

    temp = b.real * b.real + b.imag * b.imag;
    if(!temp)
    {
        fprintf( stderr, "z_div_t : temp is zero. Abort\n");
        abort();
    }

    cp->real = ( a.real * b.real + a.imag * b.imag ) / temp;
    cp->imag = ( a.imag * b.real - a.real * b.imag ) / temp;
}

/* work value calculation */
integer get_work_value( integer datatype, void *work )
{
    integer value;

    if(!work)
        return 0;

    switch(datatype)
    {
        case INTEGER:
        {
            value = (*(integer*)work);
            break;
        }
        case FLOAT:
        {
            value = (integer) (*(float*)work);
            break;
        }
        case DOUBLE:
        {
            value = (integer) (*(double*)work);
            break;
        }
        case COMPLEX:
        {
            value = (integer) (((scomplex *)work)->real);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            value = (integer) (((dcomplex *)work)->real);
            break;
        }
        default:
        {
            value = 0;
            break;
        }
    }
    return value;
}

void diagmv( integer datatype, integer m, integer n, void* x, integer incx, void* a, integer a_rs, integer a_cs )
{
    integer inca, lda;
    integer n_iter;
    integer n_elem;
    integer j;

    if(m == 0 || n == 0)
        return;

    // Initialize with optimal values for column-major storage.
    inca   = a_rs;
    lda    = a_cs;
    n_iter = n;
    n_elem = m;

    switch(datatype)
    {
        case FLOAT:
        {
            float *a_begin;
            for ( j = 0; j < n_iter; j++ )
            {
                a_begin = (float *)a + j*lda;
                scalv( datatype, n_elem, x, incx, a_begin, inca );
            }
            break;
        }

        case DOUBLE:
        {
            double *a_begin;
            for ( j = 0; j < n_iter; j++ )
            {
                a_begin = (double *)a + j*lda;
                scalv( datatype, n_elem, x, incx, a_begin, inca );
            }
            break;
        }

        case COMPLEX:
        {
            scomplex *a_begin;
            for ( j = 0; j < n_iter; j++ )
            {
                a_begin = (scomplex *)a + j*lda;
                scalv( datatype, n_elem, x, incx, a_begin, inca );
            }
            break;
        }

        case DOUBLE_COMPLEX:
        {
            dcomplex *a_begin;
            for ( j = 0; j < n_iter; j++ )
            {
                a_begin = (dcomplex *)a + j*lda;
                scalv( datatype, n_elem, x, incx, a_begin, inca );
            }
            break;
        }
    }
}

void scalv( integer datatype, integer n, void* x, integer incx, void* y, integer incy )
{
    integer i;

    switch(datatype)
    {
        case FLOAT:
        {
            float *chi, *psi;
            for ( i = 0; i < n; ++i )
            {
                chi = (float *)x + i*incx;
                psi = (float *)y + i*incy;

                (*psi) = (*chi) * (*psi);
            }
            break;
        }

        case DOUBLE:
        {
            double *chi, *psi;
            for ( i = 0; i < n; ++i )
            {
                chi = (double *)x + i*incx;
                psi = (double *)y + i*incy;

                (*psi) = (*chi) * (*psi);
            }
            break;
        }

        case COMPLEX:
        {
            float *chi;
            scomplex *psi;

            for ( i = 0; i < n; ++i )
            {
                chi = (float *)x + i*incx;
                psi = (scomplex *)y + i*incy;

                psi->real = (*chi) * (psi)->real;
                psi->imag = (*chi) * (psi)->imag;
            }
            break;
        }

        case DOUBLE_COMPLEX:
        {
            double *chi;
            dcomplex *psi;
            for ( i = 0; i < n; ++i )
            {
                chi = (double *)x + i*incx;
                psi = (dcomplex *)y + i*incy;

                psi->real = (*chi) * (psi)->real;
                psi->imag = (*chi) * (psi)->imag;
            }
            break;
        }
    }
}

void set_transpose(integer datatype, char *uplo, char *trans_A, char *trans_B)
{
    if(*uplo == 'L')
    {
        *trans_A = 'N';
        *trans_B = 'C';
    }
    else
    {
	*trans_A = 'C';
        *trans_B = 'N';
    }
}

void rand_spd_matrix(integer datatype, char *uplo, void **A, integer m,integer lda)
{
    void *sample = NULL;
    void *buff_A = NULL, *buff_B = NULL;
    void *I = NULL;
    char trans_A, trans_B;

    create_matrix(datatype, &sample, m, m);
    create_matrix(datatype, &buff_A, m, m);
    create_matrix(datatype, &buff_B, m, m);

    reset_matrix(datatype, m, m, buff_A, m);
    reset_matrix(datatype, m, m, buff_B, m);

    create_matrix(datatype, &I, m, m);
    set_identity_matrix(datatype, m, m, I, m);

    /* Generate random symmetric matrix */
    rand_sym_matrix(datatype, sample, m, m, lda);

    /* Based on uplo set the transpose flag */
    set_transpose(datatype, uplo, &trans_A, &trans_B);

    copy_matrix(datatype, uplo, m, m, sample, m, buff_A, lda);
    copy_matrix(datatype, uplo, m, m, sample, m, buff_B, lda);

    switch(datatype)
    {
        case FLOAT:
        {
            float beta = m;
            sgemm_(&trans_A, &trans_B, &m, &m, &m, &s_one, buff_A, &m, buff_B, &m, &beta, I, &m);
            break;
        }
        case DOUBLE:
        {
            double beta = m;
            dgemm_(&trans_A, &trans_B, &m, &m, &m, &d_one, buff_A, &m, buff_B, &m, &beta, I, &m);
            break;
        }
        case COMPLEX:
        {
            scomplex beta = {m,0};
            cgemm_(&trans_A, &trans_B, &m, &m, &m, &c_one, buff_A, &m, buff_B, &m, &beta, I, &m);
            break;
        }
        case DOUBLE_COMPLEX:
        {
             dcomplex beta = {m,0};
             zgemm_(&trans_A, &trans_B, &m, &m, &m, &z_one, buff_A, &m, buff_B, &m, &beta, I, &m);
             break;
        }
    }
    copy_matrix(datatype, "full", m, m, I, m, *A, lda);

    /* free buffers */
    free_matrix(sample);
    free_matrix(buff_A);
    free_matrix(buff_B);
    free_matrix(I);

    return;
}

void diagonalize_vector(integer datatype, void* s, void* sigma, integer m, integer n, integer LDA)
{
    integer incr, i, j, min_m_n;

    incr = m + 1;
    min_m_n = min(m, n);

    reset_matrix(datatype, m, n, sigma, m);

    switch( datatype )
    {
        case FLOAT:
        {
            scopy_(&min_m_n, s, &i_one, sigma, &incr);		
            break;
        }
        case DOUBLE:
        {
            dcopy_(&min_m_n, s, &i_one, sigma, &incr);		
            break;
        }
        case COMPLEX:
        {
            for( i = 0; i < n; i++ )
            {
                for( j = i; j < m; j++ )
                {
                    if(i == j)
                        ((scomplex *)sigma)[i * LDA + j].real = ((float *)s)[i];
                }
            }
            break;
        }
        case DOUBLE_COMPLEX:
        {
            for( i = 0; i < n; i++ )
            {
                for( j = i; j < m; j++ )
                {
                    if(i == j)
                        ((dcomplex *)sigma)[i * LDA + j].real = ((double *)s)[i];
                }
            }
            break;
        }
    }
    return;
}

/* Generate random Hermitian matrix */
void rand_hermitian_matrix(integer datatype, integer n, void** A, integer lda)
{
    void *B = NULL, *Herm = NULL;

    create_matrix(datatype, &B, n, n);
    create_matrix(datatype, &Herm, n, n);

    reset_matrix(datatype, n, n, B, n);
    reset_matrix(datatype, n, n, Herm, n);

    rand_matrix(datatype, B, n, n, n);

    switch(datatype)
    {
        case COMPLEX:
        {
            cgemm_("N", "C", &n, &n, &n, &c_one, B, &n, B, &n, &c_zero, Herm, &n);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            zgemm_("N", "C", &n, &n, &n, &z_one, B, &n, B, &n, &z_zero, Herm, &n);
            break;
        }
    }
    copy_matrix(datatype, "full", n, n, Herm, lda, *A, lda);
    free_matrix(B);
    free_matrix(Herm);
    return;
}
