##Copyright (C) 2020, Advanced Micro Devices, Inc.##
if(INCLUDE_LAPACKE)
target_sources("${PROJECT_NAME}"
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_c_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_cgb_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_cgb_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_cge_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_cge_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_cgg_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_cgg_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_cgt_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_chb_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_chb_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_che_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_che_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_chp_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_chp_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_chs_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_chs_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_cpb_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_cpb_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_cpf_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_cpf_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_cpo_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_cpo_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_cpp_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_cpp_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_cpt_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_csp_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_csp_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_cst_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_csy_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_csy_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_ctb_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_ctb_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_ctf_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_ctf_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_ctp_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_ctp_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_ctr_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_ctr_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_d_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_dgb_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_dgb_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_dge_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_dge_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_dgg_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_dgg_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_dgt_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_dhs_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_dhs_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_dpb_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_dpb_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_dpf_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_dpf_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_dpo_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_dpo_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_dpp_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_dpp_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_dpt_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_dsb_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_dsb_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_dsp_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_dsp_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_dst_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_dsy_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_dsy_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_dtb_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_dtb_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_dtf_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_dtf_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_dtp_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_dtp_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_dtr_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_dtr_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_lsame.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_make_complex_double.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_make_complex_float.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_s_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_sgb_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_sgb_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_sge_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_sge_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_sgg_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_sgg_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_sgt_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_shs_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_shs_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_spb_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_spb_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_spf_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_spf_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_spo_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_spo_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_spp_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_spp_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_spt_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_ssb_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_ssb_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_ssp_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_ssp_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_sst_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_ssy_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_ssy_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_stb_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_stb_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_stf_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_stf_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_stp_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_stp_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_str_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_str_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_xerbla.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_z_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_zgb_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_zgb_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_zge_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_zge_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_zgg_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_zgg_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_zgt_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_zhb_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_zhb_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_zhe_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_zhe_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_zhp_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_zhp_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_zhs_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_zhs_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_zpb_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_zpb_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_zpf_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_zpf_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_zpo_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_zpo_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_zpp_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_zpp_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_zpt_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_zsp_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_zsp_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_zst_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_zsy_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_zsy_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_ztb_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_ztb_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_ztf_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_ztf_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_ztp_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_ztp_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_ztr_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_ztr_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_ctz_nancheck.c 
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_dtz_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_stz_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_ztz_nancheck.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_ctz_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_dtz_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_stz_trans.c
${CMAKE_CURRENT_SOURCE_DIR}/lapacke_ztz_trans.c
)
else(NOT INCLUDE_LAPACKE)

set(UTILS
lapacke_c_nancheck.c            lapacke_ctr_trans.c             lapacke_make_complex_float.c    lapacke_zgb_nancheck.c
lapacke_cgb_nancheck.c          lapacke_d_nancheck.c            lapacke_s_nancheck.c            lapacke_zgb_trans.c
lapacke_cgb_trans.c             lapacke_dgb_nancheck.c          lapacke_sgb_nancheck.c          lapacke_zge_nancheck.c
lapacke_cge_nancheck.c          lapacke_dgb_trans.c             lapacke_sgb_trans.c             lapacke_zge_trans.c
lapacke_cge_trans.c             lapacke_dge_nancheck.c          lapacke_sge_nancheck.c          lapacke_zgg_nancheck.c
lapacke_cgg_nancheck.c          lapacke_dge_trans.c             lapacke_sge_trans.c             lapacke_zgg_trans.c
lapacke_cgg_trans.c             lapacke_dgg_nancheck.c          lapacke_sgg_nancheck.c          lapacke_zgt_nancheck.c
lapacke_cgt_nancheck.c          lapacke_dgg_trans.c             lapacke_sgg_trans.c             lapacke_zhb_nancheck.c
lapacke_chb_nancheck.c          lapacke_dgt_nancheck.c          lapacke_sgt_nancheck.c          lapacke_zhb_trans.c
lapacke_chb_trans.c             lapacke_dhs_nancheck.c          lapacke_shs_nancheck.c          lapacke_zhe_nancheck.c
lapacke_che_nancheck.c          lapacke_dhs_trans.c             lapacke_shs_trans.c             lapacke_zhe_trans.c
lapacke_che_trans.c             lapacke_dpb_nancheck.c          lapacke_spb_nancheck.c          lapacke_zhp_nancheck.c
lapacke_chp_nancheck.c          lapacke_dpb_trans.c             lapacke_spb_trans.c             lapacke_zhp_trans.c
lapacke_chp_trans.c             lapacke_dpf_nancheck.c          lapacke_spf_nancheck.c          lapacke_zhs_nancheck.c
lapacke_chs_nancheck.c          lapacke_dpf_trans.c             lapacke_spf_trans.c             lapacke_zhs_trans.c
lapacke_chs_trans.c             lapacke_dpo_nancheck.c          lapacke_spo_nancheck.c          lapacke_zpb_nancheck.c
lapacke_cpb_nancheck.c          lapacke_dpo_trans.c             lapacke_spo_trans.c             lapacke_zpb_trans.c
lapacke_cpb_trans.c             lapacke_dpp_nancheck.c          lapacke_spp_nancheck.c          lapacke_zpf_nancheck.c
lapacke_cpf_nancheck.c          lapacke_dpp_trans.c             lapacke_spp_trans.c             lapacke_zpf_trans.c
lapacke_cpf_trans.c             lapacke_dpt_nancheck.c          lapacke_spt_nancheck.c          lapacke_zpo_nancheck.c
lapacke_cpo_nancheck.c          lapacke_dsb_nancheck.c          lapacke_ssb_nancheck.c          lapacke_zpo_trans.c
lapacke_cpo_trans.c             lapacke_dsb_trans.c             lapacke_ssb_trans.c             lapacke_zpp_nancheck.c
lapacke_cpp_nancheck.c          lapacke_dsp_nancheck.c          lapacke_ssp_nancheck.c          lapacke_zpp_trans.c
lapacke_cpp_trans.c             lapacke_dsp_trans.c             lapacke_ssp_trans.c             lapacke_zpt_nancheck.c
lapacke_cpt_nancheck.c          lapacke_dst_nancheck.c          lapacke_sst_nancheck.c          lapacke_zsp_nancheck.c
lapacke_csp_nancheck.c          lapacke_dsy_nancheck.c          lapacke_ssy_nancheck.c          lapacke_zsp_trans.c
lapacke_csp_trans.c             lapacke_dsy_trans.c             lapacke_ssy_trans.c             lapacke_zst_nancheck.c
lapacke_cst_nancheck.c          lapacke_dtb_nancheck.c          lapacke_stb_nancheck.c          lapacke_zsy_nancheck.c
lapacke_csy_nancheck.c          lapacke_dtb_trans.c             lapacke_stb_trans.c             lapacke_zsy_trans.c
lapacke_csy_trans.c             lapacke_dtf_nancheck.c          lapacke_stf_nancheck.c          lapacke_ztb_nancheck.c
lapacke_ctb_nancheck.c          lapacke_dtf_trans.c             lapacke_stf_trans.c             lapacke_ztb_trans.c
lapacke_ctb_trans.c             lapacke_dtp_nancheck.c          lapacke_stp_nancheck.c          lapacke_ztf_nancheck.c
lapacke_ctf_nancheck.c          lapacke_dtp_trans.c             lapacke_stp_trans.c             lapacke_ztf_trans.c
lapacke_ctf_trans.c             lapacke_dtr_nancheck.c          lapacke_str_nancheck.c          lapacke_ztp_nancheck.c
lapacke_ctp_nancheck.c          lapacke_dtr_trans.c             lapacke_str_trans.c             lapacke_ztp_trans.c
lapacke_ctp_trans.c             lapacke_lsame.c                 lapacke_xerbla.c                lapacke_ztr_nancheck.c
lapacke_ctr_nancheck.c          lapacke_make_complex_double.c   lapacke_z_nancheck.c            lapacke_ztr_trans.c
lapacke_ctz_nancheck.c          lapacke_dtz_nancheck.c          lapacke_stz_nancheck.c          lapacke_ztz_nancheck.c
lapacke_ctz_trans.c             lapacke_dtz_trans.c             lapacke_stz_trans.c             lapacke_ztz_trans.c
)
endif()
