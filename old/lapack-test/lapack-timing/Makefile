#######################################################################
#  This makefile runs the timing programs for the linear equation routines
#  and the eigenvalue routines in LAPACK.  The timing output files
#  are grouped as follows:
#
#       SLINTIM,SEIGTIM  -- Single precision real time routines
#       CLINTIM,CEIGTIM  -- Single precision complex time routines
#       DLINTIM,DEIGTIM  -- Double precision real time routines
#       ZLINTIM,ZEIGTIM  -- Double precision complex time routines
#
#  Timing programs can be executed for all or some of the four different
#  precisions.  Enter 'make' followed by one or more of the data types
#  desired.
#
#  By default, 'make' alone executes the timings with SMALL data sets.
#
#  Some examples:
#       make single
#       make single complex
#       make single double complex complex16
#  Alternatively, the command
#       make
#  without any arguments runs all eight test programs with the SMALL
#  timing data files.
#
#  If you wish to use the LARGE data sets for your timings, you
#  need to type
#       make large
#  for all tests, which would be equivalent to typing
#       make single_large double_large complex_large complex16_large
#  which could also be executed separately, for example, as
#       make single_large
#
#  The executable files are called:
#       xlintims, xlintimd, xlintimc, and xlintimz for LIN
#       xeigtims, xeigtimd, xeigtimc, and xeigtimz for EIG
#  and exist in the current directory level.
#
#  To remove the output files after the tests have been run, enter
#       make clean
#
#  To re-run specific tests after a make, enter (for example):
#       'rm stime.out; make'  or:
#       'make stime.out' or:
#       'touch stime.in; make' (to re-run the timings.)
#
#       'rm *time.out; make'  (to re-run all the timings.)
#
#######################################################################

include ../make.inc

all:     single complex double complex16
large:	 single_large complex_large double_large complex16_large
small:   single complex double complex16

SEIGTIM= sgeptim.out \
         sneptim.out \
         sseptim.out \
         ssvdtim.out

SEIGTM2= SGEPTIM.out \
         SNEPTIM.out \
         SSEPTIM.out \
         SSVDTIM.out

CEIGTIM= cgeptim.out \
         cneptim.out \
         cseptim.out \
         csvdtim.out

CEIGTM2= CGEPTIM.out \
         CNEPTIM.out \
         CSEPTIM.out \
         CSVDTIM.out

DEIGTIM= dgeptim.out \
         dneptim.out \
         dseptim.out \
         dsvdtim.out

DEIGTM2= DGEPTIM.out \
         DNEPTIM.out \
         DSEPTIM.out \
         DSVDTIM.out

ZEIGTIM= zgeptim.out \
         zneptim.out \
         zseptim.out \
         zsvdtim.out

ZEIGTM2= ZGEPTIM.out \
         ZNEPTIM.out \
         ZSEPTIM.out \
         ZSVDTIM.out

SLINTIM= stime.out \
	 sband.out \
	 stime2.out

SLINTM2= STIME.out \
	 SBAND.out \
	 STIME2.out

CLINTIM= ctime.out \
	 cband.out \
	 ctime2.out

CLINTM2= CTIME.out \
	 CBAND.out \
	 CTIME2.out

DLINTIM= dtime.out \
	 dband.out \
	 dtime2.out

DLINTM2= DTIME.out \
	 DBAND.out \
	 DTIME2.out

ZLINTIM= ztime.out \
	 zband.out \
	 ztime2.out

ZLINTM2= ZTIME.out \
	 ZBAND.out \
	 ZTIME2.out

single:         $(SLINTIM) $(SEIGTIM)
complex:        $(CLINTIM) $(CEIGTIM)
double:         $(DLINTIM) $(DEIGTIM)
complex16:      $(ZLINTIM) $(ZEIGTIM)

single_large:     $(SLINTM2) $(SEIGTM2)
complex_large:    $(CLINTM2) $(CEIGTM2)
double_large:     $(DLINTM2) $(DEIGTM2)
complex16_large:  $(ZLINTM2) $(ZEIGTM2)

#
# ======== SINGLE LIN TIMINGS ===========================

stime.out: stime.in xlintims
	@echo Timing square REAL LAPACK linear equation routines
	./xlintims < stime.in > $@ 2>&1

STIME.out: STIME.in xlintims
	@echo Timing square REAL LAPACK linear equation routines
	./xlintims < STIME.in > $@ 2>&1

sband.out: sband.in xlintims
	@echo Timing banded REAL LAPACK linear equation routines
	./xlintims < sband.in > $@ 2>&1

SBAND.out: SBAND.in xlintims
	@echo Timing banded REAL LAPACK linear equation routines
	./xlintims < SBAND.in > $@ 2>&1

stime2.out: stime2.in xlintims
	@echo Timing rectangular REAL LAPACK linear equation routines
	./xlintims < stime2.in > $@ 2>&1

STIME2.out: STIME2.in xlintims
	@echo Timing rectangular REAL LAPACK linear equation routines
	./xlintims < STIME2.in > $@ 2>&1
#
# ======== COMPLEX LIN TIMINGS ==========================

ctime.out: ctime.in xlintimc
	@echo Timing square COMPLEX LAPACK linear equation routines
	./xlintimc < ctime.in > $@ 2>&1

CTIME.out: CTIME.in xlintimc
	@echo Timing square COMPLEX LAPACK linear equation routines
	./xlintimc < CTIME.in > $@ 2>&1

cband.out: cband.in xlintimc
	@echo Timing banded COMPLEX LAPACK linear equation routines
	./xlintimc < cband.in > $@ 2>&1

CBAND.out: CBAND.in xlintimc
	@echo Timing banded COMPLEX LAPACK linear equation routines
	./xlintimc < CBAND.in > $@ 2>&1

ctime2.out: ctime2.in xlintimc
	@echo Timing rectangular COMPLEX LAPACK linear equation routines
	./xlintimc < ctime2.in > $@ 2>&1

CTIME2.out: CTIME2.in xlintimc
	@echo Timing rectangular COMPLEX LAPACK linear equation routines
	./xlintimc < CTIME2.in > $@ 2>&1
#
# ======== DOUBLE LIN TIMINGS ===========================

dtime.out: dtime.in xlintimd
	@echo Timing square DOUBLE PRECISION LAPACK linear equation routines
	./xlintimd < dtime.in > $@ 2>&1

DTIME.out: DTIME.in xlintimd
	@echo Timing square DOUBLE PRECISION LAPACK linear equation routines
	./xlintimd < DTIME.in > $@ 2>&1

dband.out: dband.in xlintimd
	@echo Timing banded DOUBLE PRECISION LAPACK linear equation routines
	./xlintimd < dband.in > $@ 2>&1

DBAND.out: dband.in xlintimd
	@echo Timing banded DOUBLE PRECISION LAPACK linear equation routines
	./xlintimd < DBAND.in > $@ 2>&1

dtime2.out: dtime2.in xlintimd
	@echo Timing rectangular DOUBLE PRECISION LAPACK linear equation routines
	./xlintimd < dtime2.in > $@ 2>&1

DTIME2.out: DTIME2.in xlintimd
	@echo Timing rectangular DOUBLE PRECISION LAPACK linear equation routines
	./xlintimd < DTIME2.in > $@ 2>&1
#
# ======== COMPLEX16 LIN TIMINGS ========================

ztime.out: ztime.in xlintimz
	@echo Timing square COMPLEX16 LAPACK linear equation routines
	./xlintimz < ztime.in > $@ 2>&1

ZTIME.out: ztime.in xlintimz
	@echo Timing square COMPLEX16 LAPACK linear equation routines
	./xlintimz < ZTIME.in > $@ 2>&1

zband.out: zband.in xlintimz
	@echo Timing banded COMPLEX16 LAPACK linear equation routines
	./xlintimz < zband.in > $@ 2>&1

ZBAND.out: ZBAND.in xlintimz
	@echo Timing banded COMPLEX16 LAPACK linear equation routines
	./xlintimz < ZBAND.in > $@ 2>&1

ztime2.out: ztime2.in xlintimz
	@echo Timing rectangular COMPLEX16 LAPACK linear equation routines
	./xlintimz < ztime2.in > $@ 2>&1

ZTIME2.out: ZTIME2.in xlintimz
	@echo Timing rectangular COMPLEX16 LAPACK linear equation routines
	./xlintimz < ZTIME2.in > $@ 2>&1
#
#
# ======== SINGLE EIG TIMINGS ===========================
#
sgeptim.out: sgeptim.in xeigtims
	@echo GEP: Timing REAL Generalized Nonsymmetric Eigenvalue Problem routines
	./xeigtims < sgeptim.in > $@ 2>&1

SGEPTIM.out: SGEPTIM.in xeigtims
	@echo GEP: Timing REAL Generalized Nonsymmetric Eigenvalue Problem routines
	./xeigtims < SGEPTIM.in > $@ 2>&1

sneptim.out: sneptim.in xeigtims
	@echo NEP: Timing REAL Nonsymmetric Eigenvalue Problem routines
	./xeigtims < sneptim.in > $@ 2>&1

SNEPTIM.out: SNEPTIM.in xeigtims
	@echo NEP: Timing REAL Nonsymmetric Eigenvalue Problem routines
	./xeigtims < SNEPTIM.in > $@ 2>&1

sseptim.out: sseptim.in xeigtims
	@echo SEP: Timing REAL Symmetric Eigenvalue Problem routines
	./xeigtims < sseptim.in > $@ 2>&1

SSEPTIM.out: SSEPTIM.in xeigtims
	@echo SEP: Timing REAL Symmetric Eigenvalue Problem routines
	./xeigtims < SSEPTIM.in > $@ 2>&1

ssvdtim.out: ssvdtim.in xeigtims
	@echo SVD: Timing REAL Singular Value Decomposition routines
	./xeigtims < ssvdtim.in > $@ 2>&1

SSVDTIM.out: SSVDTIM.in xeigtims
	@echo SVD: Timing REAL Singular Value Decomposition routines
	./xeigtims < SSVDTIM.in > $@ 2>&1
#
# ======== COMPLEX EIG TIMINGS ===========================
#
cgeptim.out: cgeptim.in xeigtimc
	@echo GEP: Timing COMPLEX Generalized Nonsymmetric Eigenvalue Problem routines
	./xeigtimc < cgeptim.in > $@ 2>&1

CGEPTIM.out: CGEPTIM.in xeigtimc
	@echo GEP: Timing COMPLEX Generalized Nonsymmetric Eigenvalue Problem routines
	./xeigtimc < cgeptim.in > $@ 2>&1

cneptim.out: cneptim.in xeigtimc
	@echo NEP: Timing COMPLEX Nonsymmetric Eigenvalue Problem routines
	./xeigtimc < cneptim.in > $@ 2>&1

CNEPTIM.out: CNEPTIM.in xeigtimc
	@echo NEP: Timing COMPLEX Nonsymmetric Eigenvalue Problem routines
	./xeigtimc < CNEPTIM.in > $@ 2>&1

cseptim.out: cseptim.in xeigtimc
	@echo SEP: Timing COMPLEX Symmetric Eigenvalue Problem routines
	./xeigtimc < cseptim.in > $@ 2>&1

CSEPTIM.out: CSEPTIM.in xeigtimc
	@echo SEP: Timing COMPLEX Symmetric Eigenvalue Problem routines
	./xeigtimc < CSEPTIM.in > $@ 2>&1

csvdtim.out: csvdtim.in xeigtimc
	@echo SVD: Timing COMPLEX Singular Value Decomposition routines
	./xeigtimc < csvdtim.in > $@ 2>&1

CSVDTIM.out: CSVDTIM.in xeigtimc
	@echo SVD: Timing COMPLEX Singular Value Decomposition routines
	./xeigtimc < CSVDTIM.in > $@ 2>&1
#
# ======== DOUBLE EIG TIMINGS ===========================
#
dgeptim.out: dgeptim.in xeigtimd
	@echo GEP: Timing DOUBLE PRECISION Generalized Nonsymmetric Eigenvalue Problem routines
	./xeigtimd < dgeptim.in > $@ 2>&1

DGEPTIM.out: DGEPTIM.in xeigtimd
	@echo GEP: Timing DOUBLE PRECISION Generalized Nonsymmetric Eigenvalue Problem routines
	./xeigtimd < dgeptim.in > $@ 2>&1

dneptim.out: dneptim.in xeigtimd
	@echo NEP: Timing DOUBLE PRECISION Nonsymmetric Eigenvalue Problem routines
	./xeigtimd < dneptim.in > $@ 2>&1

DNEPTIM.out: DNEPTIM.in xeigtimd
	@echo NEP: Timing DOUBLE PRECISION Nonsymmetric Eigenvalue Problem routines
	./xeigtimd < DNEPTIM.in > $@ 2>&1

dseptim.out: dseptim.in xeigtimd
	@echo SEP: Timing DOUBLE PRECISION Symmetric Eigenvalue Problem routines
	./xeigtimd < dseptim.in > $@ 2>&1

DSEPTIM.out: DSEPTIM.in xeigtimd
	@echo SEP: Timing DOUBLE PRECISION Symmetric Eigenvalue Problem routines
	./xeigtimd < DSEPTIM.in > $@ 2>&1

dsvdtim.out: dsvdtim.in xeigtimd
	@echo SVD: Timing DOUBLE PRECISION Singular Value Decomposition routines
	./xeigtimd < dsvdtim.in > $@ 2>&1

DSVDTIM.out: DSVDTIM.in xeigtimd
	@echo SVD: Timing DOUBLE PRECISION Singular Value Decomposition routines
	./xeigtimd < DSVDTIM.in > $@ 2>&1
#
# ======== COMPLEX16 EIG TIMINGS ===========================
#
zgeptim.out: zgeptim.in xeigtimz
	@echo GEP: Timing COMPLEX16 Generalized Nonsymmetric Eigenvalue Problem routines
	./xeigtimz < zgeptim.in > $@ 2>&1

ZGEPTIM.out: ZGEPTIM.in xeigtimz
	@echo GEP: Timing COMPLEX16 Generalized Nonsymmetric Eigenvalue Problem routines
	./xeigtimz < zgeptim.in > $@ 2>&1

zneptim.out: zneptim.in xeigtimz
	@echo NEP: Timing COMPLEX16 Nonsymmetric Eigenvalue Problem routines
	./xeigtimz < zneptim.in > $@ 2>&1

ZNEPTIM.out: ZNEPTIM.in xeigtimz
	@echo NEP: Timing COMPLEX16 Nonsymmetric Eigenvalue Problem routines
	./xeigtimz < ZNEPTIM.in > $@ 2>&1

zseptim.out: zseptim.in xeigtimz
	@echo SEP: Timing COMPLEX16 Symmetric Eigenvalue Problem routines
	./xeigtimz < zseptim.in > $@ 2>&1

ZSEPTIM.out: ZSEPTIM.in xeigtimz
	@echo SEP: Timing COMPLEX16 Symmetric Eigenvalue Problem routines
	./xeigtimz < ZSEPTIM.in > $@ 2>&1

zsvdtim.out: zsvdtim.in xeigtimz
	@echo SVD: Timing COMPLEX16 Singular Value Decomposition routines
	./xeigtimz < zsvdtim.in > $@ 2>&1

ZSVDTIM.out: ZSVDTIM.in xeigtimz
	@echo SVD: Timing COMPLEX16 Singular Value Decomposition routines
	./xeigtimz < ZSVDTIM.in > $@ 2>&1
# ==============================================================================

xlintims:
	cd LIN/LINSRC ; $(MAKE) single
	cd LIN ; $(MAKE) single

xlintimc:
	cd LIN/LINSRC ; $(MAKE) complex
	cd LIN ; $(MAKE) complex

xlintimd:
	cd LIN/LINSRC ; $(MAKE) double
	cd LIN ; $(MAKE) double

xlintimz:
	cd LIN/LINSRC ; $(MAKE) complex16
	cd LIN ; $(MAKE) complex16


xeigtims:
	cd EIG/EIGSRC ; $(MAKE) single
	cd EIG ; $(MAKE) single

xeigtimc:
	cd EIG/EIGSRC ; $(MAKE) complex
	cd EIG ; $(MAKE) complex

xeigtimd:
	cd EIG/EIGSRC ; $(MAKE) double
	cd EIG ; $(MAKE) double

xeigtimz:
	cd EIG/EIGSRC ; $(MAKE) complex16
	cd EIG ; $(MAKE) complex16

clean:
	rm -f *.out core

cleanup:
	rm -f x* *.out core
