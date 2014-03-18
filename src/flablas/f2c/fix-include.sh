#!/bin/bash

for cfile in *.c; do

	cp ${cfile} ${cfile}.temp
	cat ${cfile}.temp | sed -e "s/include \"f2c.h\"/include \"FLA_f2c.h\"/g" \
	                  > ${cfile}
	rm ${cfile}.temp
done

