/*
    Copyright (C) 2022-2023, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_common.h"
#include "test_prototype.h"

/* Generates Orthogonal matrix from ORGTR() after SYTRD() call. */
void invoke_sytrd(integer datatype, char *uplo, char compz, integer n, void *A, integer lda, void *D, void *E, integer *info)
{
    void *tau = NULL, *work = NULL;
    integer lwork = -1;

    create_vector(datatype, &tau, n-1);
    switch (datatype)
    {
        case FLOAT:
        {
            create_vector(datatype, &work, 1);
            
            fla_lapack_ssytrd(uplo, &n, NULL, &lda, NULL, NULL, tau, work, &lwork, info);

            /* Get work size */
            lwork = get_work_value( datatype, work );
            free_vector(work);
            create_vector(datatype, &work, lwork);

            /* ssytrd_ to form symmetric tridiagonal matrix */
            fla_lapack_ssytrd(uplo, &n, A, &lda, D, E, tau, work, &lwork, info);

            free_vector(work);

            lwork = -1;
            create_vector(datatype, &work, 1);
            fla_lapack_sorgtr(uplo, &n, NULL, &lda, tau, work, &lwork, info);

            /* Get work size */
            lwork = get_work_value( datatype, work );
            free_vector(work);
            create_vector(datatype, &work, lwork);

            /* sorgtr_ to generate orthogonal matrix */
            fla_lapack_sorgtr(uplo, &n, A, &lda, tau, work, &lwork, info);

            free_vector(work);
            break;
        }
        case DOUBLE:
        {
            create_vector(datatype, &work, 1);

            fla_lapack_dsytrd(uplo, &n, NULL, &lda, NULL, NULL, tau, work, &lwork, info);

            /* Get work size */
            lwork = get_work_value( datatype, work );
            free_vector(work);
            create_vector(datatype, &work, lwork);

            /* dsytrd_ to form symmetric tridiagonal matrix */
            fla_lapack_dsytrd(uplo, &n, A, &lda, D, E, tau, work, &lwork, info);

            free_vector(work);

            lwork = -1;
            create_vector(datatype, &work, 1);
            fla_lapack_dorgtr(uplo, &n, NULL, &lda, tau, work, &lwork, info);

            /* Get work size */
            lwork = get_work_value( datatype, work );
            free_vector(work);
            create_vector(datatype, &work, lwork);

            /* dorgtr_ to generate orthogonal matrix */
            fla_lapack_dorgtr(uplo, &n, A, &lda, tau, work, &lwork, info);

            free_vector(work);
            break;
        }
        case COMPLEX:
        {
            create_vector(datatype, &work, 1);

            fla_lapack_chetrd(uplo, &n, NULL, &lda, NULL, NULL, tau, work, &lwork, info);

            /* Get work size */
            lwork = get_work_value( datatype, work );
            free_vector(work);
            create_vector(datatype, &work, lwork);

            /* chetrd_ to form symmetric tridiagonal matrix */
            fla_lapack_chetrd(uplo, &n, A, &lda, D, E, tau, work, &lwork, info);

            free_vector(work);

            lwork = -1;
            create_vector(datatype, &work, 1);
            fla_lapack_cungtr(uplo, &n, NULL, &lda, tau, work, &lwork, info);

            /* Get work size */
            lwork = get_work_value( datatype, work );
            free_vector(work);
            create_vector(datatype, &work, lwork);

            /* cungtr_ to generate orthogonal matrix */
            fla_lapack_cungtr(uplo, &n, A, &lda, tau, work, &lwork, info);

            free_vector(work);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            create_vector(datatype, &work, 1);

            fla_lapack_zhetrd(uplo, &n, NULL, &lda, NULL, NULL, tau, work, &lwork, info);

            /* Get work size */
            lwork = get_work_value( datatype, work );
            free_vector(work);
            create_vector(datatype, &work, lwork);

            /* zhetrd_ to form symmetric tridiagonal matrix */
            fla_lapack_zhetrd(uplo, &n, A, &lda, D, E, tau, work, &lwork, info);

            free_vector(work);

            lwork = -1;
            create_vector(datatype, &work, 1);
            fla_lapack_zungtr(uplo, &n, NULL, &lda, tau, work, &lwork, info);

            /* Get work size */
            lwork = get_work_value( datatype, work );
            free_vector(work);
            create_vector(datatype, &work, lwork);

            /* zungtr_ to generate orthogonal matrix */
            fla_lapack_zungtr(uplo, &n, A, &lda, tau, work, &lwork, info);

            free_vector(work);
            break;
        }
    }
    /* Free buffers */
    free_vector(tau);
}
