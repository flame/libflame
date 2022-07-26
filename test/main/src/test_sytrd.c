/*
    Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
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
            
            ssytrd_(uplo, &n, NULL, &lda, NULL, NULL, tau, work, &lwork, info);

            /* Get work size */
            lwork = get_work_value( datatype, work );
            free_vector(work);
            create_vector(datatype, &work, lwork);

            /* ssytrd_ to form symmetric tridiagonal matrix */
            ssytrd_(uplo, &n, A, &lda, D, E, tau, work, &lwork, info);

            free_vector(work);

            lwork = -1;
            create_vector(datatype, &work, 1);
            sorgtr_(uplo, &n, NULL, &lda, tau, work, &lwork, info);

            /* Get work size */
            lwork = get_work_value( datatype, work );
            free_vector(work);
            create_vector(datatype, &work, lwork);

            /* sorgtr_ to generate orthogonal matrix */
            sorgtr_(uplo, &n, A, &lda, tau, work, &lwork, info);

            free_vector(work);
            break;
        }
        case DOUBLE:
        {
            create_vector(datatype, &work, 1);

            dsytrd_(uplo, &n, NULL, &lda, NULL, NULL, tau, work, &lwork, info);

            /* Get work size */
            lwork = get_work_value( datatype, work );
            free_vector(work);
            create_vector(datatype, &work, lwork);

            /* dsytrd_ to form symmetric tridiagonal matrix */
            dsytrd_(uplo, &n, A, &lda, D, E, tau, work, &lwork, info);

            free_vector(work);

            lwork = -1;
            create_vector(datatype, &work, 1);
            dorgtr_(uplo, &n, NULL, &lda, tau, work, &lwork, info);

            /* Get work size */
            lwork = get_work_value( datatype, work );
            free_vector(work);
            create_vector(datatype, &work, lwork);

            /* dorgtr_ to generate orthogonal matrix */
            dorgtr_(uplo, &n, A, &lda, tau, work, &lwork, info);

            free_vector(work);
            break;
        }
        case COMPLEX:
        {
            create_vector(datatype, &work, 1);

            chetrd_(uplo, &n, NULL, &lda, NULL, NULL, tau, work, &lwork, info);

            /* Get work size */
            lwork = get_work_value( datatype, work );
            free_vector(work);
            create_vector(datatype, &work, lwork);

            /* chetrd_ to form symmetric tridiagonal matrix */
            chetrd_(uplo, &n, A, &lda, D, E, tau, work, &lwork, info);

            free_vector(work);

            lwork = -1;
            create_vector(datatype, &work, 1);
            cungtr_(uplo, &n, NULL, &lda, tau, work, &lwork, info);

            /* Get work size */
            lwork = get_work_value( datatype, work );
            free_vector(work);
            create_vector(datatype, &work, lwork);

            /* cungtr_ to generate orthogonal matrix */
            cungtr_(uplo, &n, A, &lda, tau, work, &lwork, info);

            free_vector(work);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            create_vector(datatype, &work, 1);

            zhetrd_(uplo, &n, NULL, &lda, NULL, NULL, tau, work, &lwork, info);

            /* Get work size */
            lwork = get_work_value( datatype, work );
            free_vector(work);
            create_vector(datatype, &work, lwork);

            /* zhetrd_ to form symmetric tridiagonal matrix */
            zhetrd_(uplo, &n, A, &lda, D, E, tau, work, &lwork, info);

            free_vector(work);

            lwork = -1;
            create_vector(datatype, &work, 1);
            zungtr_(uplo, &n, NULL, &lda, tau, work, &lwork, info);

            /* Get work size */
            lwork = get_work_value( datatype, work );
            free_vector(work);
            create_vector(datatype, &work, lwork);

            /* zungtr_ to generate orthogonal matrix */
            zungtr_(uplo, &n, A, &lda, tau, work, &lwork, info);

            free_vector(work);
            break;
        }
    }
    /* Free buffers */
    free_vector(tau);
}
