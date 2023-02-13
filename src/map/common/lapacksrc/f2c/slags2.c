/* ../netlib/slags2.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b SLAGS2 computes 2-by-2 orthogonal matrices U, V, and Q, and applies them to matrices A and B su ch that the rows of the transformed A and B are parallel. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SLAGS2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slags2. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slags2. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slags2. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SLAGS2( UPPER, A1, A2, A3, B1, B2, B3, CSU, SNU, CSV, */
/* SNV, CSQ, SNQ ) */
/* .. Scalar Arguments .. */
/* LOGICAL UPPER */
/* REAL A1, A2, A3, B1, B2, B3, CSQ, CSU, CSV, SNQ, */
/* $ SNU, SNV */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLAGS2 computes 2-by-2 orthogonal matrices U, V and Q, such */
/* > that if ( UPPER ) then */
/* > */
/* > U**T *A*Q = U**T *( A1 A2 )*Q = ( x 0 ) */
/* > ( 0 A3 ) ( x x ) */
/* > and */
/* > V**T*B*Q = V**T *( B1 B2 )*Q = ( x 0 ) */
/* > ( 0 B3 ) ( x x ) */
/* > */
/* > or if ( .NOT.UPPER ) then */
/* > */
/* > U**T *A*Q = U**T *( A1 0 )*Q = ( x x ) */
/* > ( A2 A3 ) ( 0 x ) */
/* > and */
/* > V**T*B*Q = V**T*( B1 0 )*Q = ( x x ) */
/* > ( B2 B3 ) ( 0 x ) */
/* > */
/* > The rows of the transformed A and B are parallel, where */
/* > */
/* > U = ( CSU SNU ), V = ( CSV SNV ), Q = ( CSQ SNQ ) */
/* > ( -SNU CSU ) ( -SNV CSV ) ( -SNQ CSQ ) */
/* > */
/* > Z**T denotes the transpose of Z. */
/* > */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPPER */
/* > \verbatim */
/* > UPPER is LOGICAL */
/* > = .TRUE.: the input matrices A and B are upper triangular. */
/* > = .FALSE.: the input matrices A and B are lower triangular. */
/* > \endverbatim */
/* > */
/* > \param[in] A1 */
/* > \verbatim */
/* > A1 is REAL */
/* > \endverbatim */
/* > */
/* > \param[in] A2 */
/* > \verbatim */
/* > A2 is REAL */
/* > \endverbatim */
/* > */
/* > \param[in] A3 */
/* > \verbatim */
/* > A3 is REAL */
/* > On entry, A1, A2 and A3 are elements of the input 2-by-2 */
/* > upper (lower) triangular matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] B1 */
/* > \verbatim */
/* > B1 is REAL */
/* > \endverbatim */
/* > */
/* > \param[in] B2 */
/* > \verbatim */
/* > B2 is REAL */
/* > \endverbatim */
/* > */
/* > \param[in] B3 */
/* > \verbatim */
/* > B3 is REAL */
/* > On entry, B1, B2 and B3 are elements of the input 2-by-2 */
/* > upper (lower) triangular matrix B. */
/* > \endverbatim */
/* > */
/* > \param[out] CSU */
/* > \verbatim */
/* > CSU is REAL */
/* > \endverbatim */
/* > */
/* > \param[out] SNU */
/* > \verbatim */
/* > SNU is REAL */
/* > The desired orthogonal matrix U. */
/* > \endverbatim */
/* > */
/* > \param[out] CSV */
/* > \verbatim */
/* > CSV is REAL */
/* > \endverbatim */
/* > */
/* > \param[out] SNV */
/* > \verbatim */
/* > SNV is REAL */
/* > The desired orthogonal matrix V. */
/* > \endverbatim */
/* > */
/* > \param[out] CSQ */
/* > \verbatim */
/* > CSQ is REAL */
/* > \endverbatim */
/* > */
/* > \param[out] SNQ */
/* > \verbatim */
/* > SNQ is REAL */
/* > The desired orthogonal matrix Q. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup realOTHERauxiliary */
/* ===================================================================== */
/* Subroutine */
int slags2_(logical *upper, real *a1, real *a2, real *a3, real *b1, real *b2, real *b3, real *csu, real *snu, real *csv, real * snv, real *csq, real *snq)
{
    /* System generated locals */
    real r__1;
    /* Local variables */
    real a, b, c__, d__, r__, s1, s2, ua11, ua12, ua21, ua22, vb11, vb12, vb21, vb22, csl, csr, snl, snr, aua11, aua12, aua21, aua22, avb11, avb12, avb21, avb22, ua11r, ua22r, vb11r, vb22r;
    extern /* Subroutine */
    int slasv2_(real *, real *, real *, real *, real * , real *, real *, real *, real *), slartg_(real *, real *, real *, real *, real *);
    /* -- LAPACK auxiliary routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Parameters .. */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    if (*upper)
    {
        /* Input matrices A and B are upper triangular matrices */
        /* Form matrix C = A*adj(B) = ( a b ) */
        /* ( 0 d ) */
        a = *a1 * *b3;
        d__ = *a3 * *b1;
        b = *a2 * *b1 - *a1 * *b2;
        /* The SVD of real 2-by-2 triangular C */
        /* ( CSL -SNL )*( A B )*( CSR SNR ) = ( R 0 ) */
        /* ( SNL CSL ) ( 0 D ) ( -SNR CSR ) ( 0 T ) */
        slasv2_(&a, &b, &d__, &s1, &s2, &snr, &csr, &snl, &csl);
        if (f2c_abs(csl) >= f2c_abs(snl) || f2c_abs(csr) >= f2c_abs(snr))
        {
            /* Compute the (1,1) and (1,2) elements of U**T *A and V**T *B, */
            /* and (1,2) element of |U|**T *|A| and |V|**T *|B|. */
            ua11r = csl * *a1;
            ua12 = csl * *a2 + snl * *a3;
            vb11r = csr * *b1;
            vb12 = csr * *b2 + snr * *b3;
            aua12 = f2c_abs(csl) * f2c_abs(*a2) + f2c_abs(snl) * f2c_abs(*a3);
            avb12 = f2c_abs(csr) * f2c_abs(*b2) + f2c_abs(snr) * f2c_abs(*b3);
            /* zero (1,2) elements of U**T *A and V**T *B */
            if (f2c_abs(ua11r) + f2c_abs(ua12) != 0.f)
            {
                if (aua12 / (f2c_abs(ua11r) + f2c_abs(ua12)) <= avb12 / (f2c_abs(vb11r) + f2c_abs(vb12)))
                {
                    r__1 = -ua11r;
                    slartg_(&r__1, &ua12, csq, snq, &r__);
                }
                else
                {
                    r__1 = -vb11r;
                    slartg_(&r__1, &vb12, csq, snq, &r__);
                }
            }
            else
            {
                r__1 = -vb11r;
                slartg_(&r__1, &vb12, csq, snq, &r__);
            }
            *csu = csl;
            *snu = -snl;
            *csv = csr;
            *snv = -snr;
        }
        else
        {
            /* Compute the (2,1) and (2,2) elements of U**T *A and V**T *B, */
            /* and (2,2) element of |U|**T *|A| and |V|**T *|B|. */
            ua21 = -snl * *a1;
            ua22 = -snl * *a2 + csl * *a3;
            vb21 = -snr * *b1;
            vb22 = -snr * *b2 + csr * *b3;
            aua22 = f2c_abs(snl) * f2c_abs(*a2) + f2c_abs(csl) * f2c_abs(*a3);
            avb22 = f2c_abs(snr) * f2c_abs(*b2) + f2c_abs(csr) * f2c_abs(*b3);
            /* zero (2,2) elements of U**T*A and V**T*B, and then swap. */
            if (f2c_abs(ua21) + f2c_abs(ua22) != 0.f)
            {
                if (aua22 / (f2c_abs(ua21) + f2c_abs(ua22)) <= avb22 / (f2c_abs(vb21) + f2c_abs(vb22)))
                {
                    r__1 = -ua21;
                    slartg_(&r__1, &ua22, csq, snq, &r__);
                }
                else
                {
                    r__1 = -vb21;
                    slartg_(&r__1, &vb22, csq, snq, &r__);
                }
            }
            else
            {
                r__1 = -vb21;
                slartg_(&r__1, &vb22, csq, snq, &r__);
            }
            *csu = snl;
            *snu = csl;
            *csv = snr;
            *snv = csr;
        }
    }
    else
    {
        /* Input matrices A and B are lower triangular matrices */
        /* Form matrix C = A*adj(B) = ( a 0 ) */
        /* ( c d ) */
        a = *a1 * *b3;
        d__ = *a3 * *b1;
        c__ = *a2 * *b3 - *a3 * *b2;
        /* The SVD of real 2-by-2 triangular C */
        /* ( CSL -SNL )*( A 0 )*( CSR SNR ) = ( R 0 ) */
        /* ( SNL CSL ) ( C D ) ( -SNR CSR ) ( 0 T ) */
        slasv2_(&a, &c__, &d__, &s1, &s2, &snr, &csr, &snl, &csl);
        if (f2c_abs(csr) >= f2c_abs(snr) || f2c_abs(csl) >= f2c_abs(snl))
        {
            /* Compute the (2,1) and (2,2) elements of U**T *A and V**T *B, */
            /* and (2,1) element of |U|**T *|A| and |V|**T *|B|. */
            ua21 = -snr * *a1 + csr * *a2;
            ua22r = csr * *a3;
            vb21 = -snl * *b1 + csl * *b2;
            vb22r = csl * *b3;
            aua21 = f2c_abs(snr) * f2c_abs(*a1) + f2c_abs(csr) * f2c_abs(*a2);
            avb21 = f2c_abs(snl) * f2c_abs(*b1) + f2c_abs(csl) * f2c_abs(*b2);
            /* zero (2,1) elements of U**T *A and V**T *B. */
            if (f2c_abs(ua21) + f2c_abs(ua22r) != 0.f)
            {
                if (aua21 / (f2c_abs(ua21) + f2c_abs(ua22r)) <= avb21 / (f2c_abs(vb21) + f2c_abs(vb22r)))
                {
                    slartg_(&ua22r, &ua21, csq, snq, &r__);
                }
                else
                {
                    slartg_(&vb22r, &vb21, csq, snq, &r__);
                }
            }
            else
            {
                slartg_(&vb22r, &vb21, csq, snq, &r__);
            }
            *csu = csr;
            *snu = -snr;
            *csv = csl;
            *snv = -snl;
        }
        else
        {
            /* Compute the (1,1) and (1,2) elements of U**T *A and V**T *B, */
            /* and (1,1) element of |U|**T *|A| and |V|**T *|B|. */
            ua11 = csr * *a1 + snr * *a2;
            ua12 = snr * *a3;
            vb11 = csl * *b1 + snl * *b2;
            vb12 = snl * *b3;
            aua11 = f2c_abs(csr) * f2c_abs(*a1) + f2c_abs(snr) * f2c_abs(*a2);
            avb11 = f2c_abs(csl) * f2c_abs(*b1) + f2c_abs(snl) * f2c_abs(*b2);
            /* zero (1,1) elements of U**T*A and V**T*B, and then swap. */
            if (f2c_abs(ua11) + f2c_abs(ua12) != 0.f)
            {
                if (aua11 / (f2c_abs(ua11) + f2c_abs(ua12)) <= avb11 / (f2c_abs(vb11) + f2c_abs(vb12)))
                {
                    slartg_(&ua12, &ua11, csq, snq, &r__);
                }
                else
                {
                    slartg_(&vb12, &vb11, csq, snq, &r__);
                }
            }
            else
            {
                slartg_(&vb12, &vb11, csq, snq, &r__);
            }
            *csu = snr;
            *snu = csr;
            *csv = snl;
            *snv = csl;
        }
    }
    return 0;
    /* End of SLAGS2 */
}
/* slags2_ */
