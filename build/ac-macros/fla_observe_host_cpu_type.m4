AC_DEFUN([FLA_OBSERVE_HOST_CPU_TYPE],
[
	AC_REQUIRE([AC_CANONICAL_HOST])
	
	case $host in
		dnl Intel Pentium-based class of processors, as well as those processors
		dnl such as AMD Athlons, Durons, etc that fall under the i?86 family.
		i386*-*-* | i586*-*-* | i686*-*-*)
			fla_host_cpu=$host_cpu
			fla_c_compiler_list="gcc icc cc"
			fla_mpicc_compiler_list="mpicc"
			fla_ar=ar
			fla_ranlib=ranlib
		;;
		dnl Intel EM64T or AMD Opteron/Athlon64 processors.
		x86_64*-*-*)
			fla_host_cpu=$host_cpu
			fla_c_compiler_list="gcc icc pathcc cc"
			fla_mpicc_compiler_list="mpicc"
			fla_ar=ar
			fla_ranlib=ranlib
		;;
		dnl Intel Itanium processors.
		ia64*-*-*)
			fla_host_cpu=$host_cpu
			fla_c_compiler_list="icc gcc cc"
			fla_mpicc_compiler_list="mpicc"
			fla_ar=ar
			fla_ranlib=ranlib
		;;
		dnl NEC SX systems.
		sx*-nec-superux*)
			fla_host_cpu=$host_cpu
			fla_c_compiler_list="sxcc"
			fla_mpicc_compiler_list="mpicc"
			fla_ar=sxar
			fla_ranlib=ranlib
		;;
		dnl IBM POWER/AIX systems.
		powerpc*-ibm-aix*)
			fla_host_cpu=$host_cpu
			fla_c_compiler_list="xlc"
			fla_mpicc_compiler_list="mpcc mpicc"
			fla_ar=ar
			fla_ranlib=ranlib
		;;
		dnl PowerPC/Cell systems.
		powerpc64-*-linux-gnu)
			if test "$fla_enable_cell_spu_parallelism" = "yes" ; then
				fla_host_cpu=$host_cpu
				fla_c_compiler_list="ppu-gcc"
				fla_mpicc_compiler_list="mpicc"
				fla_ar=ppu-ar
				fla_ranlib=ppu-ranlib
			else
				fla_host_cpu=$host_cpu
				fla_c_compiler_list="gcc xlc"
				fla_mpicc_compiler_list="mpicc"
				fla_ar=ar
				fla_ranlib=ranlib
			fi
		;;
		dnl For all other proessors, use a basic search path.
		*)
			fla_host_cpu=$host_cpu
			fla_c_compiler_list="gcc cc" 
			fla_mpicc_compiler_list="mpicc"
			fla_ar=ar
			fla_ranlib=ranlib
		;;
	esac
	
	dnl Substitute the cpu type into the autoconf output files
	AC_SUBST(fla_host_cpu)

])
