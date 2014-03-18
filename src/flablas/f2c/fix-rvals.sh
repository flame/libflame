#!/bin/bash

blasfile="cdotc.c"
cp ${blasfile} ${blasfile}.temp
cat ${blasfile}.temp | sed -e "s/\/\* Complex \*\/ VOID/complex/g" \
                     | sed -e "s/complex \* ret_val,/\/\*complex \* ret_val,\*\//g" \
                     | sed -e "s/^{/{\n    complex ret_val;\n/g" \
                     | sed -e "s/ret_val->r/ret_val\.r/g" \
                     | sed -e "s/ret_val->i/ret_val\.i/g" \
                     | sed -e "s/return/return ret_val/g" \
                     > ${blasfile}
rm ${blasfile}.temp

blasfile="cdotu.c"
cp ${blasfile} ${blasfile}.temp
cat ${blasfile}.temp | sed -e "s/\/\* Complex \*\/ VOID/complex/g" \
                     | sed -e "s/complex \* ret_val,/\/\*complex \* ret_val,\*\//g" \
                     | sed -e "s/^{/{\n    complex ret_val;\n/g" \
                     | sed -e "s/ret_val->r/ret_val\.r/g" \
                     | sed -e "s/ret_val->i/ret_val\.i/g" \
                     | sed -e "s/return/return ret_val/g" \
                     > ${blasfile}
rm ${blasfile}.temp

blasfile="zdotc.c"
cp ${blasfile} ${blasfile}.temp
cat ${blasfile}.temp | sed -e "s/\/\* Double Complex \*\/ VOID/doublecomplex/g" \
                     | sed -e "s/doublecomplex \* ret_val,/\/\*doublecomplex \* ret_val,\*\//g" \
                     | sed -e "s/^{/{\n    doublecomplex ret_val;\n/g" \
                     | sed -e "s/ret_val->r/ret_val\.r/g" \
                     | sed -e "s/ret_val->i/ret_val\.i/g" \
                     | sed -e "s/return/return ret_val/g" \
                     > ${blasfile}
rm ${blasfile}.temp

blasfile="zdotu.c"
cp ${blasfile} ${blasfile}.temp
cat ${blasfile}.temp | sed -e "s/\/\* Double Complex \*\/ VOID/doublecomplex/g" \
                     | sed -e "s/doublecomplex \* ret_val,/\/\*doublecomplex \* ret_val,\*\//g" \
                     | sed -e "s/^{/{\n    doublecomplex ret_val;\n/g" \
                     | sed -e "s/ret_val->r/ret_val\.r/g" \
                     | sed -e "s/ret_val->i/ret_val\.i/g" \
                     | sed -e "s/return/return ret_val/g" \
                     > ${blasfile}
rm ${blasfile}.temp

