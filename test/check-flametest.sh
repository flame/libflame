#!/bin/sh
#
#  LIBFLAME
#  An object-based framework for developing high-performance BLAS-like
#  libraries.
#
#  Copyright (C) 2020, Advanced Micro Devices, Inc. All rights reserved.
#
#  Redistribution and use in source and binary forms, with or without
#  modification, are permitted provided that the following conditions are
#  met:
#   - Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
#   - Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in the
#     documentation and/or other materials provided with the distribution.
#   - Neither the name(s) of the copyright holder(s) nor the names of its
#     contributors may be used to endorse or promote products derived
#     from this software without specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
#  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
#  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
#  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
#  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#

script_name=${0##*/}

ansi_red="\033[0;31m"       #prints the message on console in red
ansi_green="\033[0;32m"     #prints the message on console in green
ansi_normal="\033[0m"

passmsg="All LIBFLAME tests passed!"
failmsg0="Some tests have failed."
failmsg1="Please see test/output.test for details!"

grep -q FAILURE $1

if [ $? -eq 0 ]; then
	printf "${ansi_red}""${script_name}: ${failmsg0}""${ansi_normal}\n"
	printf "${ansi_red}""${script_name}: ${failmsg1}""${ansi_normal}\n"
	exit 0
else
	printf "${ansi_green}""${script_name}: ${passmsg}""${ansi_normal}\n"
	exit 0
fi

