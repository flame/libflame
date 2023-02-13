#!/bin/bash

print_usage()
{
    local script_name

    # Echo usage info
    echo " "
    echo " "${script_name}
    echo " "
    echo " Usage:"
    echo "   ${script_name} (no arguments)"
    echo " "

    # Exit with non-zero exit status
    exit 1
}

main()
{
	script_name=${0##*/}

	location="../fortran"

	# Check the number of argements
	if [ $# != "0" ]; then
		print_usage
	fi

	tmp_filename="tmp_file.f"

	xblas_leaves="cla_gbrfsx_extended.f cla_gerfsx_extended.f cla_herfsx_extended.f cla_porfsx_extended.f cla_syrfsx_extended.f dla_gbrfsx_extended.f dla_gerfsx_extended.f dla_porfsx_extended.f dla_syrfsx_extended.f sla_gbrfsx_extended.f sla_gerfsx_extended.f sla_porfsx_extended.f sla_syrfsx_extended.f zla_gbrfsx_extended.f zla_gerfsx_extended.f zla_herfsx_extended.f zla_porfsx_extended.f zla_syrfsx_extended.f"

	xblas_dependents="cgbrfsx.f cgbsvxx.f cgerfsx.f cgesvxx.f cherfsx.f chesvxx.f cporfsx.f cposvxx.f csyrfsx.f csysvxx.f dgbrfsx.f dgbsvxx.f dgerfsx.f dgesvxx.f dporfsx.f dposvxx.f dsyrfsx.f dsysvxx.f sgbrfsx.f sgbsvxx.f sgerfsx.f sgesvxx.f sporfsx.f sposvxx.f ssyrfsx.f ssysvxx.f zgbrfsx.f zgbsvxx.f zgerfsx.f zgesvxx.f zherfsx.f zhesvxx.f zporfsx.f zposvxx.f zsyrfsx.f zsysvxx.f"

	xblas_all="${xblas_leaves} ${xblas_dependents}"

	echo "Adding macros to disable xblas:"

	cd ${location}

	for filename in ${xblas_all}; do

		echo "  ${filename}"

		echo "#ifdef FLA_ENABLE_XBLAS" > ${tmp_filename}
		cat ${filename} >> ${tmp_filename}
		echo "#endif" >> ${tmp_filename}

		mv ${tmp_filename} ${filename}
	done

	return 0
}

main "$@"
