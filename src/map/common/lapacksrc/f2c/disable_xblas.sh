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

	# Check the number of argements
	if [ $# != "0" ]; then
		print_usage
	fi

	tmp_filename="tmp_file.c"

	xblas_leaves="cla_gbrfsx_extended.c cla_gerfsx_extended.c cla_herfsx_extended.c cla_porfsx_extended.c cla_syrfsx_extended.c dla_gbrfsx_extended.c dla_gerfsx_extended.c dla_porfsx_extended.c dla_syrfsx_extended.c sla_gbrfsx_extended.c sla_gerfsx_extended.c sla_porfsx_extended.c sla_syrfsx_extended.c zla_gbrfsx_extended.c zla_gerfsx_extended.c zla_herfsx_extended.c zla_porfsx_extended.c zla_syrfsx_extended.c"

	xblas_dependents="cgbrfsx.c cgbsvxx.c cgerfsx.c cgesvxx.c cherfsx.c chesvxx.c cporfsx.c cposvxx.c csyrfsx.c csysvxx.c dgbrfsx.c dgbsvxx.c dgerfsx.c dgesvxx.c dporfsx.c dposvxx.c dsyrfsx.c dsysvxx.c sgbrfsx.c sgbsvxx.c sgerfsx.c sgesvxx.c sporfsx.c sposvxx.c ssyrfsx.c ssysvxx.c zgbrfsx.c zgbsvxx.c zgerfsx.c zgesvxx.c zherfsx.c zhesvxx.c zporfsx.c zposvxx.c zsyrfsx.c zsysvxx.c"

	xblas_all="${xblas_leaves} ${xblas_dependents}"

	echo "Adding cpp macros to disable xblas:"

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
