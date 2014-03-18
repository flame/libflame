#!/bin/bash

execfile="test_Svd_uv_components.netlib.x"
op="svd"
impl="net331"
range="2h_3k_2h"
params="b512_k32"
shape="m1p_n1p"
dt="z"

dist="cluster"     ; nohup ./${execfile} < input_${dist} > output_${op}_${impl}_${range}_${params}_${shape}_${dist}_${dt}.m ;
dist="geometric"   ; nohup ./${execfile} < input_${dist} > output_${op}_${impl}_${range}_${params}_${shape}_${dist}_${dt}.m ;
dist="inverse"     ; nohup ./${execfile} < input_${dist} > output_${op}_${impl}_${range}_${params}_${shape}_${dist}_${dt}.m ;
dist="logarithmic" ; nohup ./${execfile} < input_${dist} > output_${op}_${impl}_${range}_${params}_${shape}_${dist}_${dt}.m ;
dist="random"      ; nohup ./${execfile} < input_${dist} > output_${op}_${impl}_${range}_${params}_${shape}_${dist}_${dt}.m ;
dist="linear"      ; nohup ./${execfile} < input_${dist} > output_${op}_${impl}_${range}_${params}_${shape}_${dist}_${dt}.m ;


