/* tstiee.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "FLA_f2c.h"


/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__10 = 10;
static integer c__2 = 2;
static integer c__3 = 3;
static integer c__4 = 4;
static integer c__11 = 11;
/*
   static integer c__0 = 0;
   static real c_b227 = 0.f;
   static real c_b228 = 1.f;
*/

/* Main program */
int MAIN__(void)
{
    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen),
            e_wsle(void);

    /* Local variables */
    integer ieeeok;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *,
                           integer *, integer *);

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 6, 0, 0, 0 };
    static cilist io___2 = { 0, 6, 0, 0, 0 };
    static cilist io___3 = { 0, 6, 0, 0, 0 };
    static cilist io___5 = { 0, 6, 0, 0, 0 };
    static cilist io___6 = { 0, 6, 0, 0, 0 };
    static cilist io___7 = { 0, 6, 0, 0, 0 };
    static cilist io___8 = { 0, 6, 0, 0, 0 };
    static cilist io___9 = { 0, 6, 0, 0, 0 };
    static cilist io___10 = { 0, 6, 0, 0, 0 };
    static cilist io___11 = { 0, 6, 0, 0, 0 };
    static cilist io___12 = { 0, 6, 0, 0, 0 };
    static cilist io___13 = { 0, 6, 0, 0, 0 };
    static cilist io___14 = { 0, 6, 0, 0, 0 };
    static cilist io___15 = { 0, 6, 0, 0, 0 };
    static cilist io___16 = { 0, 6, 0, 0, 0 };
    static cilist io___17 = { 0, 6, 0, 0, 0 };
    static cilist io___18 = { 0, 6, 0, 0, 0 };
    static cilist io___19 = { 0, 6, 0, 0, 0 };



    /*  -- LAPACK test routine (version 3.2) -- */
    /*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
    /*     November 2006 */

    /*     .. External Functions .. */
    /*     .. */
    /*     .. Local Scalars .. */
    /*     .. */
    /*     .. Executable Statements .. */

    s_wsle(&io___1);
    do_lio(&c__9, &c__1, "We are about to check whether infinity arithmetic",
           (ftnlen)49);
    e_wsle();
    s_wsle(&io___2);
    do_lio(&c__9, &c__1, "can be trusted.  If this test hangs, set", (ftnlen)
           40);
    e_wsle();
    s_wsle(&io___3);
    do_lio(&c__9, &c__1, "ILAENV = 0 for ISPEC = 10 in LAPACK/SRC/ilaenv.f", (
               ftnlen)48);
    e_wsle();

    ieeeok = ilaenv_(&c__10, "ILAENV", "N", &c__1, &c__2, &c__3, &c__4);
    s_wsle(&io___5);
    e_wsle();

    if (ieeeok == 0)
    {
        s_wsle(&io___6);
        do_lio(&c__9, &c__1, "Infinity arithmetic did not perform per the ie"
               "ee spec", (ftnlen)53);
        e_wsle();
    }
    else
    {
        s_wsle(&io___7);
        do_lio(&c__9, &c__1, "Infinity arithmetic performed as per the ieee "
               "spec.", (ftnlen)51);
        e_wsle();
        s_wsle(&io___8);
        do_lio(&c__9, &c__1, "However, this is not an exhaustive test and do"
               "es not", (ftnlen)52);
        e_wsle();
        s_wsle(&io___9);
        do_lio(&c__9, &c__1, "guarantee that infinity arithmetic meets the", (
                   ftnlen)44);
        do_lio(&c__9, &c__1, " ieee spec.", (ftnlen)11);
        e_wsle();
    }

    s_wsle(&io___10);
    e_wsle();
    s_wsle(&io___11);
    do_lio(&c__9, &c__1, "We are about to check whether NaN arithmetic", (
               ftnlen)44);
    e_wsle();
    s_wsle(&io___12);
    do_lio(&c__9, &c__1, "can be trusted.  If this test hangs, set", (ftnlen)
           40);
    e_wsle();
    s_wsle(&io___13);
    do_lio(&c__9, &c__1, "ILAENV = 0 for ISPEC = 11 in LAPACK/SRC/ilaenv.f", (
               ftnlen)48);
    e_wsle();
    ieeeok = ilaenv_(&c__11, "ILAENV", "N", &c__1, &c__2, &c__3, &c__4);

    s_wsle(&io___14);
    e_wsle();
    if (ieeeok == 0)
    {
        s_wsle(&io___15);
        do_lio(&c__9, &c__1, "NaN arithmetic did not perform per the ieee sp"
               "ec", (ftnlen)48);
        e_wsle();
    }
    else
    {
        s_wsle(&io___16);
        do_lio(&c__9, &c__1, "NaN arithmetic performed as per the ieee", (
                   ftnlen)40);
        do_lio(&c__9, &c__1, " spec.", (ftnlen)6);
        e_wsle();
        s_wsle(&io___17);
        do_lio(&c__9, &c__1, "However, this is not an exhaustive test and do"
               "es not", (ftnlen)52);
        e_wsle();
        s_wsle(&io___18);
        do_lio(&c__9, &c__1, "guarantee that NaN arithmetic meets the", (
                   ftnlen)39);
        do_lio(&c__9, &c__1, " ieee spec.", (ftnlen)11);
        e_wsle();
    }
    s_wsle(&io___19);
    e_wsle();

    return 0;
} /* MAIN__ */

/* Main program alias */ int main_ ()
{
    MAIN__ ();
    return 0;
}
