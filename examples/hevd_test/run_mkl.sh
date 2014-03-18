#!/bin/bash

execfile="test_Hevd_lv_components.mkl.x"
op="hevd"
impl="mkl102"
range="2h_3k_2h"
params="b512_k32"
dt="z"

dist="cluster"     ; nohup ./${execfile} < input_${dist} > output_${op}_${impl}_${range}_${params}_${dist}_${dt}.m ;
dist="geometric"   ; nohup ./${execfile} < input_${dist} > output_${op}_${impl}_${range}_${params}_${dist}_${dt}.m ;
dist="inverse"     ; nohup ./${execfile} < input_${dist} > output_${op}_${impl}_${range}_${params}_${dist}_${dt}.m ;
dist="logarithmic" ; nohup ./${execfile} < input_${dist} > output_${op}_${impl}_${range}_${params}_${dist}_${dt}.m ;
dist="random"      ; nohup ./${execfile} < input_${dist} > output_${op}_${impl}_${range}_${params}_${dist}_${dt}.m ;
dist="linear"      ; nohup ./${execfile} < input_${dist} > output_${op}_${impl}_${range}_${params}_${dist}_${dt}.m ;


