/*
    Copyright (C) 2022-2023, Advanced Micro Devices, Inc. All rights reserved.
*/

#ifndef TEST_PROTOTYPE_H
#define TEST_PROTOTYPE_H

/* Rename API as per API_CALLING_CONVENTION */
#if (UPPER_)
    
#define fla_lapack_sladiv SLADIV_
#define fla_lapack_dladiv DLADIV_

#define fla_lapack_slacpy SLACPY_
#define fla_lapack_dlacpy DLACPY_
#define fla_lapack_clacpy CLACPY_
#define fla_lapack_zlacpy ZLACPY_

#define fla_lapack_slaset SLASET_
#define fla_lapack_dlaset DLASET_
#define fla_lapack_claset CLASET_
#define fla_lapack_zlaset ZLASET_

#define fla_lapack_slange SLANGE_
#define fla_lapack_dlange DLANGE_
#define fla_lapack_clange CLANGE_
#define fla_lapack_zlange ZLANGE_

#define fla_lapack_slaswp SLASWP_
#define fla_lapack_dlaswp DLASWP_
#define fla_lapack_claswp CLASWP_
#define fla_lapack_zlaswp ZLASWP_

#define fla_lapack_slamch SLAMCH_
#define fla_lapack_dlamch DLAMCH_

#define fla_lapack_sorgrq SORGRQ_
#define fla_lapack_dorgrq DORGRQ_
#define fla_lapack_cungrq CUNGRQ_
#define fla_lapack_zungrq ZUNGRQ_

#define fla_lapack_sgeev SGEEV_
#define fla_lapack_dgeev DGEEV_
#define fla_lapack_cgeev CGEEV_
#define fla_lapack_zgeev ZGEEV_

#define fla_lapack_sgeevx SGEEVX_
#define fla_lapack_dgeevx DGEEVX_
#define fla_lapack_cgeevx CGEEVX_
#define fla_lapack_zgeevx ZGEEVX_

#define fla_lapack_sgesdd SGESDD_
#define fla_lapack_dgesdd DGESDD_
#define fla_lapack_cgesdd CGESDD_
#define fla_lapack_zgesdd ZGESDD_

#define fla_lapack_sgesvd SGESVD_
#define fla_lapack_dgesvd DGESVD_
#define fla_lapack_cgesvd CGESVD_
#define fla_lapack_zgesvd ZGESVD_

#define fla_lapack_sgeqrf SGEQRF_
#define fla_lapack_dgeqrf DGEQRF_
#define fla_lapack_cgeqrf CGEQRF_
#define fla_lapack_zgeqrf ZGEQRF_

#define fla_lapack_sgerqf SGERQF_
#define fla_lapack_dgerqf DGERQF_
#define fla_lapack_cgerqf CGERQF_
#define fla_lapack_zgerqf ZGERQF_

#define fla_lapack_sgerq2 SGERQ2_
#define fla_lapack_dgerq2 DGERQ2_
#define fla_lapack_cgerq2 CGERQ2_
#define fla_lapack_zgerq2 ZGERQ2_

#define fla_lapack_sgeqp3 SGEQP3_
#define fla_lapack_dgeqp3 DGEQP3_
#define fla_lapack_cgeqp3 CGEQP3_
#define fla_lapack_zgeqp3 ZGEQP3_

#define fla_lapack_spotrf SPOTRF_
#define fla_lapack_dpotrf DPOTRF_
#define fla_lapack_cpotrf CPOTRF_
#define fla_lapack_zpotrf ZPOTRF_

#define fla_lapack_spotrs SPOTRS_
#define fla_lapack_dpotrs DPOTRS_
#define fla_lapack_cpotrs CPOTRS_
#define fla_lapack_zpotrs ZPOTRS_

#define fla_lapack_sgetrf SGETRF_
#define fla_lapack_dgetrf DGETRF_
#define fla_lapack_cgetrf CGETRF_
#define fla_lapack_zgetrf ZGETRF_

#define fla_lapack_sgetri SGETRI_
#define fla_lapack_dgetri DGETRI_
#define fla_lapack_cgetri CGETRI_
#define fla_lapack_zgetri ZGETRI_

#define fla_lapack_sgetrs SGETRS_
#define fla_lapack_dgetrs DGETRS_
#define fla_lapack_cgetrs CGETRS_
#define fla_lapack_zgetrs ZGETRS_

#define fla_lapack_ssyevd SSYEVD_
#define fla_lapack_dsyevd DSYEVD_
#define fla_lapack_cheevd CHEEVD_
#define fla_lapack_zheevd ZHEEVD_

#define fla_lapack_sggevx SGGEVX_
#define fla_lapack_dggevx DGGEVX_
#define fla_lapack_cggevx CGGEVX_
#define fla_lapack_zggevx ZGGEVX_

#define fla_lapack_sgesv SGESV_
#define fla_lapack_dgesv DGESV_
#define fla_lapack_cgesv CGESV_
#define fla_lapack_zgesv ZGESV_

#define fla_lapack_slapmt SLAPMT_
#define fla_lapack_dlapmt DLAPMT_
#define fla_lapack_clapmt CLAPMT_
#define fla_lapack_zlapmt ZLAPMT_

#define fla_lapack_sggev SGGEV_
#define fla_lapack_dggev DGGEV_
#define fla_lapack_cggev CGGEV_
#define fla_lapack_zggev ZGGEV_

#define fla_lapack_sorgtr SORGTR_
#define fla_lapack_dorgtr DORGTR_
#define fla_lapack_cungtr CUNGTR_
#define fla_lapack_zungtr ZUNGTR_

#define fla_lapack_sorgqr SORGQR_
#define fla_lapack_dorgqr DORGQR_
#define fla_lapack_cungqr CUNGQR_
#define fla_lapack_zungqr ZUNGQR_

#define fla_lapack_ssteqr SSTEQR_
#define fla_lapack_dsteqr DSTEQR_
#define fla_lapack_csteqr CSTEQR_
#define fla_lapack_zsteqr ZSTEQR_

#define fla_lapack_sstevd SSTEVD_
#define fla_lapack_dstevd DSTEVD_

#define fla_lapack_ssytrd SSYTRD_
#define fla_lapack_dsytrd DSYTRD_
#define fla_lapack_chetrd CHETRD_
#define fla_lapack_zhetrd ZHETRD_

#define fla_lapack_sstedc SSTEDC_
#define fla_lapack_dstedc DSTEDC_
#define fla_lapack_cstedc CSTEDC_
#define fla_lapack_zstedc ZSTEDC_

#define fla_lapack_shseqr SHSEQR_
#define fla_lapack_dhseqr DHSEQR_
#define fla_lapack_chseqr CHSEQR_
#define fla_lapack_zhseqr ZHSEQR_

#define fla_lapack_sgehrd SGEHRD_
#define fla_lapack_dgehrd DGEHRD_
#define fla_lapack_cgehrd CGEHRD_
#define fla_lapack_zgehrd ZGEHRD_

#define fla_lapack_sorghr SORGHR_
#define fla_lapack_dorghr DORGHR_
#define fla_lapack_cunghr CUNGHR_
#define fla_lapack_zunghr ZUNGHR_

#define fla_lapack_sstedc SSYEV_
#define fla_lapack_dstedc DSYEV_
#define fla_lapack_cstedc CHEEV_
#define fla_lapack_zstedc ZHEEV_

#define fla_lapack_sspffrt2 SSPFFRT2_
#define fla_lapack_dspffrt2 DSPFFRT2_
#define fla_lapack_cspffrt2 CSPFFRT2_
#define fla_lapack_zspffrt2 ZSPFFRT2_

#define fla_lapack_sspffrtx SSPFFRTX_
#define fla_lapack_dspffrtx DSPFFRTX_
#define fla_lapack_cspffrtx CSPFFRTX_
#define fla_lapack_zspffrtx ZSPFFRTX_

#elif (UPPER)

#define fla_lapack_sladiv SLADIV
#define fla_lapack_dladiv DLADIV

#define fla_lapack_slacpy SLACPY
#define fla_lapack_dlacpy DLACPY
#define fla_lapack_clacpy CLACPY
#define fla_lapack_zlacpy ZLACPY

#define fla_lapack_slaset SLASET
#define fla_lapack_dlaset DLASET
#define fla_lapack_claset CLASET
#define fla_lapack_zlaset ZLASET

#define fla_lapack_slange SLANGE
#define fla_lapack_dlange DLANGE
#define fla_lapack_clange CLANGE
#define fla_lapack_zlange ZLANGE

#define fla_lapack_slaswp SLASWP
#define fla_lapack_dlaswp DLASWP
#define fla_lapack_claswp CLASWP
#define fla_lapack_zlaswp ZLASWP

#define fla_lapack_slamch SLAMCH
#define fla_lapack_dlamch DLAMCH

#define fla_lapack_sorgrq SORGRQ
#define fla_lapack_dorgrq DORGRQ
#define fla_lapack_cungrq CUNGRQ
#define fla_lapack_zungrq ZUNGRQ

#define fla_lapack_sgeev SGEEV
#define fla_lapack_dgeev DGEEV
#define fla_lapack_cgeev CGEEV
#define fla_lapack_zgeev ZGEEV

#define fla_lapack_sgeevx SGEEVX
#define fla_lapack_dgeevx DGEEVX
#define fla_lapack_cgeevx CGEEVX
#define fla_lapack_zgeevx ZGEEVX

#define fla_lapack_sgesdd SGESDD
#define fla_lapack_dgesdd DGESDD
#define fla_lapack_cgesdd CGESDD
#define fla_lapack_zgesdd ZGESDD

#define fla_lapack_sgesvd SGESVD
#define fla_lapack_dgesvd DGESVD
#define fla_lapack_cgesvd CGESVD
#define fla_lapack_zgesvd ZGESVD

#define fla_lapack_sgeqrf SGEQRF
#define fla_lapack_dgeqrf DGEQRF
#define fla_lapack_cgeqrf CGEQRF
#define fla_lapack_zgeqrf ZGEQRF

#define fla_lapack_sgerqf SGERQF
#define fla_lapack_dgerqf DGERQF
#define fla_lapack_cgerqf CGERQF
#define fla_lapack_zgerqf ZGERQF

#define fla_lapack_sgerq2 SGERQ2
#define fla_lapack_dgerq2 DGERQ2
#define fla_lapack_cgerq2 CGERQ2
#define fla_lapack_zgerq2 ZGERQ2

#define fla_lapack_sgeqp3 SGEQP3
#define fla_lapack_dgeqp3 DGEQP3
#define fla_lapack_cgeqp3 CGEQP3
#define fla_lapack_zgeqp3 ZGEQP3

#define fla_lapack_spotrf SPOTRF
#define fla_lapack_dpotrf DPOTRF
#define fla_lapack_cpotrf CPOTRF
#define fla_lapack_zpotrf ZPOTRF

#define fla_lapack_spotrs SPOTRS
#define fla_lapack_dpotrs DPOTRS
#define fla_lapack_cpotrs CPOTRS
#define fla_lapack_zpotrs ZPOTRS

#define fla_lapack_sgetrf SGETRF
#define fla_lapack_dgetrf DGETRF
#define fla_lapack_cgetrf CGETRF
#define fla_lapack_zgetrf ZGETRF

#define fla_lapack_sgetri SGETRI
#define fla_lapack_dgetri DGETRI
#define fla_lapack_cgetri CGETRI
#define fla_lapack_zgetri ZGETRI

#define fla_lapack_sgetrs SGETRS
#define fla_lapack_dgetrs DGETRS
#define fla_lapack_cgetrs CGETRS
#define fla_lapack_zgetrs ZGETRS

#define fla_lapack_ssyevd SSYEVD
#define fla_lapack_dsyevd DSYEVD
#define fla_lapack_cheevd CHEEVD
#define fla_lapack_zheevd ZHEEVD

#define fla_lapack_sggevx SGGEVX
#define fla_lapack_dggevx DGGEVX
#define fla_lapack_cggevx CGGEVX
#define fla_lapack_zggevx ZGGEVX

#define fla_lapack_sgesv SGESV
#define fla_lapack_dgesv DGESV
#define fla_lapack_cgesv CGESV
#define fla_lapack_zgesv ZGESV

#define fla_lapack_slapmt SLAPMT
#define fla_lapack_dlapmt DLAPMT
#define fla_lapack_clapmt CLAPMT
#define fla_lapack_zlapmt ZLAPMT

#define fla_lapack_sggev SGGEV
#define fla_lapack_dggev DGGEV
#define fla_lapack_cggev CGGEV
#define fla_lapack_zggev ZGGEV

#define fla_lapack_sorgtr SORGTR
#define fla_lapack_dorgtr DORGTR
#define fla_lapack_cungtr CUNGTR
#define fla_lapack_zungtr ZUNGTR

#define fla_lapack_sorgqr SORGQR
#define fla_lapack_dorgqr DORGQR
#define fla_lapack_cungqr CUNGQR
#define fla_lapack_zungqr ZUNGQR

#define fla_lapack_ssteqr SSTEQR
#define fla_lapack_dsteqr DSTEQR
#define fla_lapack_csteqr CSTEQR
#define fla_lapack_zsteqr ZSTEQR

#define fla_lapack_sstevd SSTEVD
#define fla_lapack_dstevd DSTEVD

#define fla_lapack_ssytrd SSYTRD
#define fla_lapack_dsytrd DSYTRD
#define fla_lapack_chetrd CHETRD
#define fla_lapack_zhetrd ZHETRD

#define fla_lapack_sstedc SSTEDC
#define fla_lapack_dstedc DSTEDC
#define fla_lapack_cstedc CSTEDC
#define fla_lapack_zstedc ZSTEDC

#define fla_lapack_shseqr SHSEQR
#define fla_lapack_dhseqr DHSEQR
#define fla_lapack_chseqr CHSEQR
#define fla_lapack_zhseqr ZHSEQR

#define fla_lapack_sgehrd SGEHRD
#define fla_lapack_dgehrd DGEHRD
#define fla_lapack_cgehrd CGEHRD
#define fla_lapack_zgehrd ZGEHRD

#define fla_lapack_sorghr SORGHR
#define fla_lapack_dorghr DORGHR
#define fla_lapack_cunghr CUNGHR
#define fla_lapack_zunghr ZUNGHR

#define fla_lapack_sstedc SSYEV
#define fla_lapack_dstedc DSYEV
#define fla_lapack_cstedc CHEEV
#define fla_lapack_zstedc ZHEEV

#define fla_lapack_sspffrt2 SSPFFRT2
#define fla_lapack_dspffrt2 DSPFFRT2
#define fla_lapack_cspffrt2 CSPFFRT2
#define fla_lapack_zspffrt2 ZSPFFRT2

#define fla_lapack_sspffrtx SSPFFRTX
#define fla_lapack_dspffrtx DSPFFRTX
#define fla_lapack_cspffrtx CSPFFRTX
#define fla_lapack_zspffrtx ZSPFFRTX

#elif (LOWER)

#define fla_lapack_sladiv sladiv
#define fla_lapack_dladiv dladiv

#define fla_lapack_slacpy slacpy
#define fla_lapack_dlacpy dlacpy
#define fla_lapack_clacpy clacpy
#define fla_lapack_zlacpy zlacpy

#define fla_lapack_slaset slaset
#define fla_lapack_dlaset dlaset
#define fla_lapack_claset claset
#define fla_lapack_zlaset zlaset

#define fla_lapack_slange slange
#define fla_lapack_dlange dlange
#define fla_lapack_clange clange
#define fla_lapack_zlange zlange

#define fla_lapack_slaswp slaswp
#define fla_lapack_dlaswp dlaswp
#define fla_lapack_claswp claswp
#define fla_lapack_zlaswp zlaswp

#define fla_lapack_slamch slamch
#define fla_lapack_dlamch dlamch

#define fla_lapack_sorgrq sorgrq
#define fla_lapack_dorgrq dorgrq
#define fla_lapack_cungrq cungrq
#define fla_lapack_zungrq zungrq

#define fla_lapack_sgeev sgeev
#define fla_lapack_dgeev dgeev
#define fla_lapack_cgeev cgeev
#define fla_lapack_zgeev zgeev

#define fla_lapack_sgeevx sgeevx
#define fla_lapack_dgeevx dgeevx
#define fla_lapack_cgeevx cgeevx
#define fla_lapack_zgeevx zgeevx

#define fla_lapack_sgesdd sgesdd
#define fla_lapack_dgesdd dgesdd
#define fla_lapack_cgesdd cgesdd
#define fla_lapack_zgesdd zgesdd

#define fla_lapack_sgesvd sgesvd
#define fla_lapack_dgesvd dgesvd
#define fla_lapack_cgesvd cgesvd
#define fla_lapack_zgesvd zgesvd

#define fla_lapack_sgeqrf sgeqrf
#define fla_lapack_dgeqrf dgeqrf
#define fla_lapack_cgeqrf cgeqrf
#define fla_lapack_zgeqrf zgeqrf

#define fla_lapack_sgerqf sgerqf
#define fla_lapack_dgerqf dgerqf
#define fla_lapack_cgerqf cgerqf
#define fla_lapack_zgerqf zgerqf

#define fla_lapack_sgerq2 sgerq2
#define fla_lapack_dgerq2 dgerq2
#define fla_lapack_cgerq2 cgerq2
#define fla_lapack_zgerq2 zgerq2

#define fla_lapack_sgeqp3 sgeqp3
#define fla_lapack_dgeqp3 dgeqp3
#define fla_lapack_cgeqp3 cgeqp3
#define fla_lapack_zgeqp3 zgeqp3

#define fla_lapack_spotrf spotrf
#define fla_lapack_dpotrf dpotrf
#define fla_lapack_cpotrf cpotrf
#define fla_lapack_zpotrf zpotrf

#define fla_lapack_spotrs spotrs
#define fla_lapack_dpotrs dpotrs
#define fla_lapack_cpotrs cpotrs
#define fla_lapack_zpotrs zpotrs

#define fla_lapack_sgetrf sgetrf
#define fla_lapack_dgetrf dgetrf
#define fla_lapack_cgetrf cgetrf
#define fla_lapack_zgetrf zgetrf

#define fla_lapack_sgetri sgetri
#define fla_lapack_dgetri dgetri
#define fla_lapack_cgetri cgetri
#define fla_lapack_zgetri zgetri

#define fla_lapack_sgetrs sgetrs
#define fla_lapack_dgetrs dgetrs
#define fla_lapack_cgetrs cgetrs
#define fla_lapack_zgetrs zgetrs

#define fla_lapack_ssyevd ssyevd
#define fla_lapack_dsyevd dsyevd
#define fla_lapack_cheevd cheevd
#define fla_lapack_zheevd zheevd

#define fla_lapack_sggevx sggevx
#define fla_lapack_dggevx dggevx
#define fla_lapack_cggevx cggevx
#define fla_lapack_zggevx zggevx

#define fla_lapack_sgesv sgesv
#define fla_lapack_dgesv dgesv
#define fla_lapack_cgesv cgesv
#define fla_lapack_zgesv zgesv

#define fla_lapack_slapmt slapmt
#define fla_lapack_dlapmt dlapmt
#define fla_lapack_clapmt clapmt
#define fla_lapack_zlapmt zlapmt

#define fla_lapack_sggev sggev
#define fla_lapack_dggev dggev
#define fla_lapack_cggev cggev
#define fla_lapack_zggev zggev

#define fla_lapack_sorgtr sorgtr
#define fla_lapack_dorgtr dorgtr
#define fla_lapack_cungtr cungtr
#define fla_lapack_zungtr zungtr

#define fla_lapack_sorgqr sorgqr
#define fla_lapack_dorgqr dorgqr
#define fla_lapack_cungqr cungqr
#define fla_lapack_zungqr zungqr

#define fla_lapack_ssteqr ssteqr
#define fla_lapack_dsteqr dsteqr
#define fla_lapack_csteqr csteqr
#define fla_lapack_zsteqr zsteqr

#define fla_lapack_sstevd sstevd
#define fla_lapack_dstevd dstevd

#define fla_lapack_ssytrd ssytrd
#define fla_lapack_dsytrd dsytrd
#define fla_lapack_chetrd chetrd
#define fla_lapack_zhetrd zhetrd

#define fla_lapack_sstedc sstedc
#define fla_lapack_dstedc dstedc
#define fla_lapack_cstedc cstedc
#define fla_lapack_zstedc zstedc

#define fla_lapack_shseqr shseqr
#define fla_lapack_dhseqr dhseqr
#define fla_lapack_chseqr chseqr
#define fla_lapack_zhseqr zhseqr

#define fla_lapack_sgehrd sgehrd
#define fla_lapack_dgehrd dgehrd
#define fla_lapack_cgehrd cgehrd
#define fla_lapack_zgehrd zgehrd

#define fla_lapack_sorghr sorghr
#define fla_lapack_dorghr dorghr
#define fla_lapack_cunghr cunghr
#define fla_lapack_zunghr zunghr

#define fla_lapack_sstedc ssyev
#define fla_lapack_dstedc dsyev
#define fla_lapack_cstedc cheev
#define fla_lapack_zstedc zheev

#define fla_lapack_sspffrt2 sspffrt2
#define fla_lapack_dspffrt2 dspffrt2
#define fla_lapack_cspffrt2 cspffrt2
#define fla_lapack_zspffrt2 zspffrt2

#define fla_lapack_sspffrtx sspffrtx
#define fla_lapack_dspffrtx dspffrtx
#define fla_lapack_cspffrtx cspffrtx
#define fla_lapack_zspffrtx zspffrtx

#else

#define fla_lapack_sladiv sladiv_
#define fla_lapack_dladiv dladiv_

#define fla_lapack_slacpy slacpy_
#define fla_lapack_dlacpy dlacpy_
#define fla_lapack_clacpy clacpy_
#define fla_lapack_zlacpy zlacpy_

#define fla_lapack_slaset slaset_
#define fla_lapack_dlaset dlaset_
#define fla_lapack_claset claset_
#define fla_lapack_zlaset zlaset_

#define fla_lapack_slange slange_
#define fla_lapack_dlange dlange_
#define fla_lapack_clange clange_
#define fla_lapack_zlange zlange_

#define fla_lapack_slaswp slaswp_
#define fla_lapack_dlaswp dlaswp_
#define fla_lapack_claswp claswp_
#define fla_lapack_zlaswp zlaswp_

#define fla_lapack_slamch slamch_
#define fla_lapack_dlamch dlamch_

#define fla_lapack_sorgrq sorgrq_
#define fla_lapack_dorgrq dorgrq_
#define fla_lapack_cungrq cungrq_
#define fla_lapack_zungrq zungrq_

#define fla_lapack_sgeev sgeev_
#define fla_lapack_dgeev dgeev_
#define fla_lapack_cgeev cgeev_
#define fla_lapack_zgeev zgeev_

#define fla_lapack_sgeevx sgeevx_
#define fla_lapack_dgeevx dgeevx_
#define fla_lapack_cgeevx cgeevx_
#define fla_lapack_zgeevx zgeevx_

#define fla_lapack_sgesdd sgesdd_
#define fla_lapack_dgesdd dgesdd_
#define fla_lapack_cgesdd cgesdd_
#define fla_lapack_zgesdd zgesdd_

#define fla_lapack_sgesvd sgesvd_
#define fla_lapack_dgesvd dgesvd_
#define fla_lapack_cgesvd cgesvd_
#define fla_lapack_zgesvd zgesvd_

#define fla_lapack_sgeqrf sgeqrf_
#define fla_lapack_dgeqrf dgeqrf_
#define fla_lapack_cgeqrf cgeqrf_
#define fla_lapack_zgeqrf zgeqrf_

#define fla_lapack_sgerqf sgerqf_
#define fla_lapack_dgerqf dgerqf_
#define fla_lapack_cgerqf cgerqf_
#define fla_lapack_zgerqf zgerqf_

#define fla_lapack_sgerq2 sgerq2_
#define fla_lapack_dgerq2 dgerq2_
#define fla_lapack_cgerq2 cgerq2_
#define fla_lapack_zgerq2 zgerq2_

#define fla_lapack_sgeqp3 sgeqp3_
#define fla_lapack_dgeqp3 dgeqp3_
#define fla_lapack_cgeqp3 cgeqp3_
#define fla_lapack_zgeqp3 zgeqp3_

#define fla_lapack_spotrf spotrf_
#define fla_lapack_dpotrf dpotrf_
#define fla_lapack_cpotrf cpotrf_
#define fla_lapack_zpotrf zpotrf_

#define fla_lapack_spotrs spotrs_
#define fla_lapack_dpotrs dpotrs_
#define fla_lapack_cpotrs cpotrs_
#define fla_lapack_zpotrs zpotrs_

#define fla_lapack_sgetrf sgetrf_
#define fla_lapack_dgetrf dgetrf_
#define fla_lapack_cgetrf cgetrf_
#define fla_lapack_zgetrf zgetrf_

#define fla_lapack_sgetri sgetri_
#define fla_lapack_dgetri dgetri_
#define fla_lapack_cgetri cgetri_
#define fla_lapack_zgetri zgetri_

#define fla_lapack_sgetrs sgetrs_
#define fla_lapack_dgetrs dgetrs_
#define fla_lapack_cgetrs cgetrs_
#define fla_lapack_zgetrs zgetrs_

#define fla_lapack_ssyevd ssyevd_
#define fla_lapack_dsyevd dsyevd_
#define fla_lapack_cheevd cheevd_
#define fla_lapack_zheevd zheevd_

#define fla_lapack_sggevx sggevx_
#define fla_lapack_dggevx dggevx_
#define fla_lapack_cggevx cggevx_
#define fla_lapack_zggevx zggevx_

#define fla_lapack_sgesv sgesv_
#define fla_lapack_dgesv dgesv_
#define fla_lapack_cgesv cgesv_
#define fla_lapack_zgesv zgesv_

#define fla_lapack_slapmt slapmt_
#define fla_lapack_dlapmt dlapmt_
#define fla_lapack_clapmt clapmt_
#define fla_lapack_zlapmt zlapmt_

#define fla_lapack_sggev sggev_
#define fla_lapack_dggev dggev_
#define fla_lapack_cggev cggev_
#define fla_lapack_zggev zggev_

#define fla_lapack_sorgtr sorgtr_
#define fla_lapack_dorgtr dorgtr_
#define fla_lapack_cungtr cungtr_
#define fla_lapack_zungtr zungtr_

#define fla_lapack_sorgqr sorgqr_
#define fla_lapack_dorgqr dorgqr_
#define fla_lapack_cungqr cungqr_
#define fla_lapack_zungqr zungqr_

#define fla_lapack_ssteqr ssteqr_
#define fla_lapack_dsteqr dsteqr_
#define fla_lapack_csteqr csteqr_
#define fla_lapack_zsteqr zsteqr_

#define fla_lapack_sstevd sstevd_
#define fla_lapack_dstevd dstevd_

#define fla_lapack_ssytrd ssytrd_
#define fla_lapack_dsytrd dsytrd_
#define fla_lapack_chetrd chetrd_
#define fla_lapack_zhetrd zhetrd_

#define fla_lapack_sstedc sstedc_
#define fla_lapack_dstedc dstedc_
#define fla_lapack_cstedc cstedc_
#define fla_lapack_zstedc zstedc_

#define fla_lapack_shseqr shseqr_
#define fla_lapack_dhseqr dhseqr_
#define fla_lapack_chseqr chseqr_
#define fla_lapack_zhseqr zhseqr_

#define fla_lapack_sgehrd sgehrd_
#define fla_lapack_dgehrd dgehrd_
#define fla_lapack_cgehrd cgehrd_
#define fla_lapack_zgehrd zgehrd_

#define fla_lapack_sorghr sorghr_
#define fla_lapack_dorghr dorghr_
#define fla_lapack_cunghr cunghr_
#define fla_lapack_zunghr zunghr_

#define fla_lapack_ssyev ssyev_
#define fla_lapack_dsyev dsyev_
#define fla_lapack_cheev cheev_
#define fla_lapack_zheev zheev_

#define fla_lapack_sspffrt2 sspffrt2_
#define fla_lapack_dspffrt2 dspffrt2_
#define fla_lapack_cspffrt2 cspffrt2_
#define fla_lapack_zspffrt2 zspffrt2_

#define fla_lapack_sspffrtx sspffrtx_
#define fla_lapack_dspffrtx dspffrtx_
#define fla_lapack_cspffrtx cspffrtx_
#define fla_lapack_zspffrtx zspffrtx_

#endif /*if UPPER_*/

/* These functions are API invoking functions used in other API test codes */
extern void invoke_getrf(integer datatype, integer* m, integer* n, void* a, integer* lda, integer* ipiv, integer* info);
extern void invoke_potrf(char* uplo, integer datatype, integer* m, void* a, integer* lda, integer* info);
extern void invoke_geqrf(integer datatype, integer* m, integer* n, void* a, integer* lda, void* tau, void* work, integer* lwork, integer* info);
/* Generates Orthogonal matrix from ORGTR() after SYTRD() call. */
extern void invoke_sytrd(integer datatype, char *uplo, char compz, integer n, void *A, integer lda, void *D, void *E, integer *info);

/* --------BLAS APIs -------- */

extern int saxpy_(integer *n, void *a, void *x, integer *incx, void *y, integer *incy);
extern int daxpy_(integer *n, void *a, void *x, integer *incx, void *y, integer *incy);
extern int caxpy_(integer *n, void *a, void *x, integer *incx, void *y, integer *incy);
extern int zaxpy_(integer *n, void *a, void *x, integer *incx, void *y, integer *incy);

extern float snrm2_(integer *n, void *x, integer *incx);
extern double dnrm2_(integer *n, void *x, integer *incx);
extern float scnrm2_(integer *n, void *x, integer *incx);
extern double dznrm2_(integer *n, void *x, integer *incx);

extern int scopy_(integer *n, void *x, integer *incx, void *y, integer *incy);
extern int dcopy_(integer *n, void *x, integer *incx, void *y, integer *incy);
extern int ccopy_(integer *n, void *x, integer *incx, void *y, integer *incy);
extern int zcopy_(integer *n, void *x, integer *incx, void *y, integer *incy);

extern int sgemv_(char *trans, integer *m, integer *n, void *alpha, void *a, integer *lda, void *x, integer *incx, void *beta, void *y, integer *incy);
extern int dgemv_(char *trans, integer *m, integer *n, void *alpha, void *a, integer *lda, void *x, integer *incx, void *beta, void *y, integer *incy);
extern int cgemv_(char *trans, integer *m, integer *n, void *alpha, void *a, integer *lda, void *x, integer *incx, void *beta, void *y, integer *incy);
extern int zgemv_(char *trans, integer *m, integer *n, void *alpha, void *a, integer *lda, void *x, integer *incx, void *beta, void *y, integer *incy);

extern int sgemm_(char *transa, char *transb, integer *m, integer * n, integer *k, void *alpha, void *a, integer *lda, void *b, integer * ldb, void *beta, void *c__, integer *ldc);
extern int dgemm_(char *transa, char *transb, integer *m, integer * n, integer *k, void *alpha, void *a, integer *lda, void *b, integer * ldb, void *beta, void *c__, integer *ldc);
extern int cgemm_(char *transa, char *transb, integer *m, integer * n, integer *k, void *alpha, void *a, integer *lda, void *b, integer * ldb, void *beta, void *c__, integer *ldc);
extern int zgemm_(char *transa, char *transb, integer *m, integer * n, integer *k, void *alpha, void *a, integer *lda, void *b, integer * ldb, void *beta, void *c__, integer *ldc);

extern void dzgemm_(char* transa, char* transb, integer* m, integer* n, integer* k, void* alpha, void* a, integer* lda, void* b, integer* ldb, void* beta, void* c, integer* ldc);

/* --------LAPACK APIs ---------*/

/* performs complex division */
extern int fla_lapack_sladiv(float *a, float *b, float *c__, float *d__, float *p, float *q);
extern int fla_lapack_dladiv(double *a, double *b, double *c__, double *d__, double *p, double *q);

extern float fla_lapack_slamch(char *);
extern double fla_lapack_dlamch(char *);

extern int fla_lapack_claswp(integer* n, scomplex* a, integer* lda, integer* k1, integer* k2, integer* ipiv, integer* incx);
extern int fla_lapack_dlaswp(integer* n, double* a, integer* lda, integer* k1, integer* k2, integer* ipiv, integer* incx);
extern int fla_lapack_slaswp(integer* n, float* a, integer* lda, integer* k1, integer* k2, integer* ipiv, integer* incx);
extern int fla_lapack_zlaswp(integer* n, dcomplex* a, integer* lda, integer* k1, integer* k2, integer* ipiv, integer* incx);

extern float fla_lapack_slange(char* norm, integer* m, integer* n, void* a, integer* lda, void* work);
extern double fla_lapack_dlange(char* norm, integer* m, integer* n, void* a, integer* lda, void* work);
extern float fla_lapack_clange(char* norm, integer* m, integer* n, void* a, integer* lda, void* work);
extern double fla_lapack_zlange(char* norm, integer* m, integer* n, void* a, integer* lda, void* work);

extern int fla_lapack_slaset(char *uplo, integer *m, integer *n, void *alpha, void *beta, void *a, integer *lda);
extern int fla_lapack_dlaset(char *uplo, integer *m, integer *n, void *alpha, void *beta, void *a, integer *lda);
extern int fla_lapack_claset(char *uplo, integer *m, integer *n, void *alpha, void *beta, void *a, integer *lda);
extern int fla_lapack_zlaset(char *uplo, integer *m, integer *n, void *alpha, void *beta, void *a, integer *lda);

extern int fla_lapack_slacpy(char* uplo, integer* m, integer* n, void* a, integer* lda, void* b, integer* ldb);
extern int fla_lapack_dlacpy(char* uplo, integer* m, integer* n, void* a, integer* lda, void* b, integer* ldb);
extern int fla_lapack_clacpy(char* uplo, integer* m, integer* n, void* a, integer* lda, void* b, integer* ldb);
extern int fla_lapack_zlacpy(char* uplo, integer* m, integer* n, void* a, integer* lda, void* b, integer* ldb);

extern int fla_lapack_sorgrq(integer* m, integer* n, integer* k, void* a, integer* lda, void* tau, void* work, integer* lwork, integer* info);
extern int fla_lapack_dorgrq(integer* m, integer* n, integer* k, void* a, integer* lda, void* tau, void* work, integer* lwork, integer* info);
extern int fla_lapack_cungrq(integer* m, integer* n, integer* k, void* a, integer* lda, void* tau, void* work, integer* lwork, integer* info);
extern int fla_lapack_zungrq(integer* m, integer* n, integer* k, void* a, integer* lda, void* tau, void* work, integer* lwork, integer* info);

/* Eigen value and Eigen vectors*/
extern int fla_lapack_sgeevx(char* balanc, char* jobvl, char* jobvr, char* sense, integer* n, void* a, integer* lda, void* wr, void* wi, void* vl, integer* ldvl, void* vr, integer* ldvr, integer* ilo, integer* ihi, void* scale, void* abnrm, void* rconde, void* rcondv, void* work, integer* lwork, integer* iwork, integer* info);
extern int fla_lapack_dgeevx(char* balanc, char* jobvl, char* jobvr, char* sense, integer* n, void* a, integer* lda, void* wr, void* wi, void* vl, integer* ldvl, void* vr, integer* ldvr, integer* ilo, integer* ihi, void* scale, void* abnrm, void* rconde, void* rcondv, void* work, integer* lwork, integer* iwork, integer* info);
extern int fla_lapack_cgeevx(char* balanc, char* jobvl, char* jobvr, char* sense, integer* n, void* a, integer* lda, void* w, void* vl, integer* ldvl, void* vr, integer* ldvr, integer* ilo, integer* ihi, void* scale, void* abnrm, void* rconde, void* rcondv, void* work, integer* lwork, void* rwork, integer* info);
extern int fla_lapack_zgeevx(char* balanc, char* jobvl, char* jobvr, char* sense, integer* n, void* a, integer* lda, void* w, void* vl, integer* ldvl, void* vr, integer* ldvr, integer* ilo, integer* ihi, void* scale, void* abnrm, void* rconde, void* rcondv, void* work, integer* lwork, void* rwork, integer* info);

int fla_lapack_sgeev(char* jobvl, char* jobvr, integer* n, void* a, integer* lda, void* wr, void* wi, void* vl, integer* ldvl, void* vr, integer* ldvr, void* work, integer* lwork, integer* info);
int fla_lapack_dgeev(char* jobvl, char* jobvr, integer* n, void* a, integer* lda, void* wr, void* wi, void* vl, integer* ldvl, void* vr, void* ldvr, void* work, integer* lwork, integer* info);
int fla_lapack_cgeev(char* jobvl, char* jobvr, integer* n, void* a, integer* lda, void* w, void* vl, integer* ldvl, void* vr, integer* ldvr, void* work, integer* lwork, void* rwork, integer* info);
int fla_lapack_zgeev(char* jobvl, char* jobvr, integer* n, void* a, integer* lda, void* w, void* vl, integer* ldvl, void* vr, integer* ldvr, void* work, integer* lwork, void* rwork, integer* info);

/* Singular Value Decomposition */
extern int fla_lapack_sgesdd(char* jobz, integer* m, integer* n, void* a, integer* lda, void* s, void* u, integer* ldu, void* vt, integer* ldvt, void* work, integer* lwork, integer* iwork, integer* info);
extern int fla_lapack_dgesdd(char* jobz, integer* m, integer* n, void* a, integer* lda, void* s, void* u, integer* ldu, void* vt, integer* ldvt, void* work, integer* lwork, integer* iwork, integer* info);
extern int fla_lapack_cgesdd(char* jobz, integer* m, integer* n, void* a, integer* lda, void* s, void* u, integer* ldu, void* vt, integer* ldvt, void* work, integer* lwork, void* rwork, integer* iwork, integer* info);
extern int fla_lapack_zgesdd(char* jobz, integer* m, integer* n, void* a, integer* lda, void* s, void* u, integer* ldu, void* vt, integer* ldvt, void* work, integer* lwork, void* rwork, integer* iwork, integer* info);

extern int fla_lapack_sgesvd(char* jobu, char* jobvt, integer* m, integer* n, void* a, integer* lda, void* s, void* u, integer* ldu, void* vt, integer* ldvt, void* work, integer* lwork, integer* info);
extern int fla_lapack_dgesvd(char* jobu, char *jobvt, integer* m, integer* n, void* a, integer* lda, void* s, void* u, integer* ldu, void* vt, integer* ldvt, void* work, integer* lwork, integer* info);
extern int fla_lapack_cgesvd(char* jobu, char* jobvt, integer* m, integer* n, void* a, integer* lda, void* s, void* u, integer* ldu, void* vt, integer* ldvt, void* work, integer* lwork, integer* rwork, integer* info);
extern int fla_lapack_zgesvd(char* jobu, char* jobvt, integer* m, integer* n, void* a, integer* lda, void* s, void* u, integer* ldu, void* vt, integer* ldvt, void* work, integer* lwork, integer* rwork, integer* info);

/* QR factorization */
extern int fla_lapack_sgeqrf(integer *m, integer *n, void *a, integer *lda, void *tau, void *work, integer *lwork, integer *info);
extern int fla_lapack_dgeqrf(integer *m, integer *n, void *a, integer *lda, void *tau, void *work, integer *lwork, integer *info);
extern int fla_lapack_cgeqrf(integer *m, integer *n, void *a, integer *lda, void *tau, void *work, integer *lwork, integer *info);
extern int fla_lapack_zgeqrf(integer *m, integer *n, void *a, integer *lda, void *tau, void *work, integer *lwork, integer *info);

/* RQ factorization APIs */
extern int fla_lapack_sgerqf(integer *m, integer *n, void *a, integer *lda, void *tau, void *work, integer *lwork, integer *info);
extern int fla_lapack_dgerqf(integer *m, integer *n, void *a, integer *lda, void *tau, void *work, integer *lwork, integer *info);
extern int fla_lapack_cgerqf(integer *m, integer *n, void *a, integer *lda, void *tau, void *work, integer *lwork, integer *info);
extern int fla_lapack_zgerqf(integer *m, integer *n, void *a, integer *lda, void *tau, void *work, integer *lwork, integer *info);

/* RQ factorization with unblocked algorithm APIs*/
extern int fla_lapack_sgerq2(integer *m, integer *n, void *a, integer *lda, void *tau, void *work, integer *info);
extern int fla_lapack_dgerq2(integer *m, integer *n, void *a, integer *lda, void *tau, void *work, integer *info);
extern int fla_lapack_cgerq2(integer *m, integer *n, void *a, integer *lda, void *tau, void *work, integer *info);
extern int fla_lapack_zgerq2(integer *m, integer *n, void *a, integer *lda, void *tau, void *work, integer *info);

/* QR factorization with column pivoting */
extern int fla_lapack_sgeqp3(integer *m, integer *n, void *a, integer *lda, integer *jpvt, void *tau, void *work, integer *lwork, integer *info);
extern int fla_lapack_dgeqp3(integer *m, integer *n, void *a, integer *lda, integer *jpvt, void *tau, void *work, integer *lwork, integer *info);
extern int fla_lapack_cgeqp3(integer *m, integer *n, void *a, integer *lda, integer *jpvt, void *tau, void *work, integer *lwork, void *rwork, integer *info);
extern int fla_lapack_zgeqp3(integer *m, integer *n, void *a, integer *lda, integer *jpvt, void *tau, void *work, integer *lwork, void *rwork, integer *info);

/* Cholesky factorization APIS*/
extern int fla_lapack_spotrf(char* uplo, integer* n, void* a, integer* lda, integer* info);
extern int fla_lapack_dpotrf(char* uplo, integer* n, void* a, integer* lda, integer* info);
extern int fla_lapack_cpotrf(char* uplo, integer* n, void* a, integer* lda, integer* info);
extern int fla_lapack_zpotrf(char* uplo, integer* n, void* a, integer* lda, integer* info);

extern int fla_lapack_spotrs(char* uplo, integer* n, integer* nrhs, void* a, integer* lda, void* b, integer* ldb, integer* info);
extern int fla_lapack_dpotrs(char* uplo, integer* n, integer* nrhs, void* a, integer* lda, void* b, integer* ldb, integer* info);
extern int fla_lapack_cpotrs(char* uplo, integer* n, integer* nrhs, void* a, integer* lda, void* b, integer* ldb, integer* info);
extern int fla_lapack_zpotrs(char* uplo, integer* n, integer* nrhs, void* a, integer* lda, void* b, integer* ldb, integer* info);

/* LU factorization of a general m by n matrix */
extern int fla_lapack_dgetrf(integer* m, integer* n, void* a, integer* lda, integer* ipiv, integer* info);
extern int fla_lapack_sgetrf(integer* m, integer* n, void* a, integer* lda, integer* ipiv, integer* info);
extern int fla_lapack_cgetrf(integer* m, integer* n, void* a, integer* lda, integer* ipiv, integer* info);
extern int fla_lapack_zgetrf(integer* m, integer* n, void* a, integer* lda, integer* ipiv, integer* info);

/* LU factorization of a general m by n matrix */
extern int fla_lapack_dgetri(integer* n, void* a, integer* lda, integer* ipiv, void *work, integer *lwork, integer* info);
extern int fla_lapack_sgetri(integer* n, void* a, integer* lda, integer* ipiv, void* work, integer* lwork, integer* info);
extern int fla_lapack_cgetri(integer* n, void* a, integer* lda, integer* ipiv, void* work, integer* lwork, integer* info);
extern int fla_lapack_zgetri(integer* n, void* a, integer* lda, integer* ipiv, void* work, integer* lwork, integer* info);

/* LU factorization of a general m by n matrix */
extern int fla_lapack_dgetrs(char* trans, integer* n, integer* nrhs, void* a, integer* lda, integer* ipiv, void* b, integer* ldb, integer* info);
extern int fla_lapack_sgetrs(char* trans, integer* n, integer* nrhs, void* a, integer* lda, integer* ipiv, void* b, integer* ldb, integer* info);
extern int fla_lapack_cgetrs(char* trans, integer* n, integer* nrhs, void* a, integer* lda, integer* ipiv, void* b, integer* ldb, integer* info);
extern int fla_lapack_zgetrs(char* trans, integer* n, integer* nrhs, void* a, integer* lda, integer* ipiv, void* b, integer* ldb, integer* info);

/* Computation of Eigen Values and Eigen Vectors */
extern int fla_lapack_ssyevd(char* jobz, char* uplo, integer* n, void* a, integer* lda, void* w, void* work, integer* lwork, integer* iwork, integer* liwork, integer* info);
extern int fla_lapack_dsyevd(char* jobz, char* uplo, integer* n, void* a, integer* lda, void* w, void* work, integer* lwork, integer* iwork, integer* liwork, integer* info);
extern int fla_lapack_cheevd(char* jobz, char* uplo, integer* n, void* a, integer* lda, void* w, void* work, integer* lwork, void* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info);
extern int fla_lapack_zheevd(char* jobz, char* uplo, integer* n, void* a, integer* lda, void* w, void* work, integer* lwork, void* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info);

/* Computation of Eigen Values and Eigen Vectors */
extern int fla_lapack_sggevx(char* balanc, char* jobvl, char* jobvr, char* sense, integer* n, void* a, integer* lda, void* b, integer* ldb, void* alphar, void* alphai, void* beta, void* vl, integer* ldvl, void* vr, integer* ldvr, integer* ilo, integer* ihi, void* lscale, void* rscale, void* abnrm, void* bbnrm, void* rconde, void* rcondv, void* work, integer* lwork, integer* iwork, integer* bwork, integer* info);
extern int fla_lapack_dggevx(char* balanc, char* jobvl, char* jobvr, char* sense, integer* n, void* a, integer* lda, void* b, integer* ldb, void* alphar, void* alphai, void* beta, void* vl, integer* ldvl, void* vr, integer* ldvr, integer* ilo, integer* ihi, void* lscale, void* rscale, void* abnrm, void* bbnrm, void* rconde, void* rcondv, void* work, integer* lwork, integer* iwork, integer* bwork, integer* info);
extern int fla_lapack_cggevx(char* balanc, char* jobvl, char* jobvr, char* sense, integer* n, void* a, integer* lda, void* b, integer* ldb, void* alpha, void* beta, void* vl, integer* ldvl, void* vr, integer* ldvr, integer* ilo, integer* ihi, void* lscale, void* rscale, void* abnrm, void* bbnrm, void* rconde, void* rcondv, void* work, integer* lwork, void* rwork, integer* iwork, integer* bwork, integer* info);
extern int fla_lapack_zggevx(char* balanc, char* jobvl, char* jobvr, char* sense, integer* n, void* a, integer* lda, void* b, integer* ldb, void* alpha, void* beta, void* vl, integer* ldvl, void* vr, integer* ldvr, integer* ilo, integer* ihi, void* lscale, void* rscale, void* abnrm, void* bbnrm, void* rconde, void* rcondv, void* work, integer* lwork, void* rwork, integer* iwork, integer* bwork, integer* info);

extern int fla_lapack_sgesv(integer* n, integer* nrhs, void* a, integer* lda, integer* ipiv, void* b, integer* ldb, integer* info);
extern int fla_lapack_dgesv(integer* n, integer* nrhs, void* a, integer* lda, integer* ipiv, void* b, integer* ldb, integer* info);
extern int fla_lapack_cgesv(integer* n, integer* nrhs, void* a, integer* lda, integer* ipiv, void* b, integer* ldb, integer* info);
extern int fla_lapack_zgesv(integer* n, integer* nrhs, void* a, integer* lda, integer* ipiv, void* b, integer* ldb, integer* info);

extern int fla_lapack_slapmt(integer* forwrd, integer* m, integer* n, void* x, integer* ldx, integer* k);
extern int fla_lapack_dlapmt(integer* forwrd, integer* m, integer* n, void* x, integer* ldx, integer* k);
extern int fla_lapack_clapmt(integer* forwrd, integer* m, integer* n, void* x, integer* ldx, integer* k);
extern int fla_lapack_zlapmt(integer* forwrd, integer* m, integer* n, void* x, integer* ldx, integer* k);

/* Computation of Eigen Values and Eigen Vectors */
extern int fla_lapack_dggev(char* jobvl, char* jobvr, integer* n, void* a, integer* lda, void* b, integer* ldb, void* alphar, void* alphai, void* beta, void* vl, integer* ldvl, void* vr, integer* ldvr, void* work, integer* lwork, integer* info);
extern int fla_lapack_sggev(char* jobvl, char* jobvr, integer* n, void* a, integer* lda, void* b, integer* ldb, void* alphar, void* alphai, void* beta, void* vl, integer* ldvl, void* vr, integer* ldvr, void* work, integer* lwork, integer* info);
extern int fla_lapack_cggev(char* jobvl, char* jobvr, integer* n, void* a, integer* lda, void* b, integer* ldb, void* alpha, void* beta, void* vl, integer* ldvl, void* vr, integer* ldvr, void* work, integer* lwork, void* rwork, integer* info);
extern int fla_lapack_zggev(char* jobvl, char* jobvr, integer* n, void* a, integer* lda, void* b, integer* ldb, void* alpha, void* beta, void* vl, integer* ldvl, void* vr, integer* ldvr, void* work, integer* lwork, void* rwork, integer* info);

/* Generates a real orthogonal matrix Q which is defined as the product of n-1 elementary reflectors of order N.*/
extern int fla_lapack_sorgtr(char* uplo, integer* m, void* a, integer* lda, void* tau, void* work, integer* lwork, integer* info);
extern int fla_lapack_dorgtr(char* uplo, integer* m, void* a, integer* lda, void* tau, void* work, integer* lwork, integer* info);
extern int fla_lapack_cungtr(char* uplo, integer* m, void* a, integer* lda, void* tau, void* work, integer* lwork, integer* info);
extern int fla_lapack_zungtr(char* uplo, integer* m, void* a, integer* lda, void* tau, void* work, integer* lwork, integer* info);

extern int fla_lapack_sorgqr(integer* m, integer* n, integer* k, void* a, integer* lda, void* tau, void* work, integer* lwork, integer* info);
extern int fla_lapack_dorgqr(integer* m, integer* n, integer* k, void* a, integer* lda, void* tau, void* work, integer* lwork, integer* info);
extern int fla_lapack_cungqr(integer* m, integer* n, integer* k, void* a, integer* lda, void* tau, void* work, integer* lwork, integer* info);
extern int fla_lapack_zungqr(integer* m, integer* n, integer* k, void* a, integer* lda, void* tau, void* work, integer* lwork, integer* info);

extern int fla_lapack_ssteqr(char* compz, integer* n, void* d, void* e, void* z, integer* ldz, void* work, integer* info);
extern int fla_lapack_dsteqr(char* compz, integer* n, void* d, void* e, void* z, integer* ldz, void* work, integer* info);
extern int fla_lapack_csteqr(char* compz, integer* n, void* d, void* e, void* z, integer* ldz, void* work, integer* info);
extern int fla_lapack_zsteqr(char* compz, integer* n, void* d, void* e, void* z, integer* ldz, void* work, integer* info);

extern int fla_lapack_sstevd(char* jobz, integer* n, void* d, void* e, void* z, integer* ldz, void* work, integer* lwork, integer* iwork, integer* liwork, integer* info);
extern int fla_lapack_dstevd(char* jobz, integer* n, void* d, void* e, void* z, integer* ldz, void* work, integer* lwork, integer* iwork, integer* liwork, integer* info);

/* Reduces real symmetric matrix A to real symmetric tridiagonal form T
   by an orthogonal similarity transformation*/
extern int fla_lapack_ssytrd(char* uplo, integer* n, void* a, integer* lda, void* d, void* e, void* tau, void* work, integer* lwork, integer* info);
extern int fla_lapack_dsytrd(char* uplo, integer* n, void* a, integer* lda, void* d, void* e, void* tau, void* work, integer* lwork, integer* info);
extern int fla_lapack_chetrd(char* uplo, integer* n, void* a, integer* lda, void* d, void* e, void* tau, void* work, integer* lwork, integer* info);
extern int fla_lapack_zhetrd(char* uplo, integer* n, void* a, integer* lda, void* d, void* e, void* tau, void* work, integer* lwork, integer* info);

/* Computation of Eigen Values and Eigen Vectors */
extern int fla_lapack_sstedc(char* compz, integer* n, void* d, void* e, void* z, integer* ldz, void* work, integer* lwork, integer* iwork, integer* liwork, integer* info);
extern int fla_lapack_dstedc(char* compz, integer* n, void* d, void* e, void* z, integer* ldz, void* work, integer* lwork, integer* iwork, integer* liwork, integer* info);
extern int fla_lapack_cstedc(char* compz, integer* n, void* d, void* e, void* z, integer* ldz, void* work, integer* lwork, void* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info);
extern int fla_lapack_zstedc(char* compz, integer* n, void* d, void* e, void* z, integer* ldz, void* work, integer* lwork, void* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info);

/* Computation of the Eigen Values of a Hessenberg matrix */
extern int fla_lapack_shseqr(char* job, char* compz, integer* n, integer* ilo, integer* ihi, void* h, integer* ldh, void* wr, void* wi, void* z, integer* ldz, void* work, integer* lwork, integer* info);
extern int fla_lapack_dhseqr(char* job, char* compz, integer* n, integer* ilo, integer* ihi, void* h, integer* ldh, void* wr, void* wi, void* z, integer* ldz, void* work, integer* lwork, integer* info);
extern int fla_lapack_chseqr(char* job, char* compz, integer* n, integer* ilo, integer* ihi, void* h, integer* ldh, void* w, void* z, integer* ldz, void* work, integer* lwork, integer* info);
extern int fla_lapack_zhseqr(char* job, char* compz, integer* n, integer* ilo, integer* ihi, void* h, integer* ldh, void* w, void* z, integer* ldz, void* work, integer* lwork, integer* info);

extern int fla_lapack_sgehrd(integer* n, integer* ilo, integer* ihi, void* a, integer* lda, void* tau, void* work, integer* lwork, integer* info);
extern int fla_lapack_dgehrd(integer* n, integer* ilo, integer* ihi, void* a, integer* lda, void* tau, void* work, integer* lwork, integer* info);
extern int fla_lapack_cgehrd(integer* n, integer* ilo, integer* ihi, void* a, integer* lda, void* tau, void* work, integer* lwork, integer* info);
extern int fla_lapack_zgehrd(integer* n, integer* ilo, integer* ihi, void* a, integer* lda, void* tau, void* work, integer* lwork, integer* info);

extern int fla_lapack_sorghr(integer *n, integer *ilo, integer *ihi, void *a, integer *lda, void *tau, void *work, integer *lwork, integer *info);
extern int fla_lapack_dorghr(integer *n, integer *ilo, integer *ihi, void *a, integer *lda, void *tau, void *work, integer *lwork, integer *info);
extern int fla_lapack_cunghr(integer *n, integer *ilo, integer *ihi, void* a, integer *lda, void *tau, void *work, integer *lwork, integer *info);
extern int fla_lapack_zunghr(integer *n, integer *ilo, integer *ihi, void* a, integer *lda, void *tau, void *work, integer *lwork, integer *info);

/* Computation of Eigen Values and Eigen Vectors*/
extern int fla_lapack_ssyev(char *jobz, char *uplo, integer *n, void *a, integer *lda, void *w, void *work, integer *lwork, integer *info);
extern int fla_lapack_dsyev(char *jobz, char *uplo, integer *n, void *a, integer *lda, void *w, void *work, integer *lwork, integer *info);
extern int fla_lapack_cheev(char *jobz, char *uplo, integer *n, void *a, integer *lda, void *w, void *work, integer *lwork, void *rwork, integer *info);
extern int fla_lapack_zheev(char *jobz, char *uplo, integer *n, void *a, integer *lda, void *w, void *work, integer *lwork, void *rwork, integer *info);

/* Computes the partial LDLT factorization of a symmetric matrix using packed storage */
extern int fla_lapack_sspffrt2(void  *ap, integer *n, integer *ncolm, void  *work, void  *work2);
extern int fla_lapack_dspffrt2(void  *ap, integer *n, integer *ncolm, void  *work, void  *work2);
extern int fla_lapack_cspffrt2(void  *ap, integer *n, integer *ncolm, void  *work, void  *work2);
extern int fla_lapack_zspffrt2(void  *ap, integer *n, integer *ncolm, void  *work, void  *work2);

extern int fla_lapack_sspffrtx(void  *ap, integer *n, integer *ncolm, void  *work, void  *work2);
extern int fla_lapack_dspffrtx(void  *ap, integer *n, integer *ncolm, void  *work, void  *work2);
extern int fla_lapack_cspffrtx(void  *ap, integer *n, integer *ncolm, void  *work, void  *work2);
extern int fla_lapack_zspffrtx(void  *ap, integer *n, integer *ncolm, void  *work, void  *work2);

#endif  // TEST_PROTOTYPE_H
