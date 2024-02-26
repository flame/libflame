#!/bin/bash

lcov --rc lcov_branch_coverage=1 --capture --directory . --output-file coverage.info;
lcov --remove --rc lcov_branch_coverage=1 coverage.info '/usr/*' '*/test/*' \
                    '*src/base/flamec/hierarchy*' '*src/base/flamec/old*' '*src/base/flamec/alt*' '*src/base/flamec/supermatrix*' \
                    '*src/base/flamec/blis*' 'base/flamec/wrappers/blas*' '*src/blas/*' '*src/flablas/*' \
                    '*src/aocl_dtl/*' '*src/lapacke/LAPACKE/example*'   \
                    -o filtered_coverage.info

genhtml --rc genhtml_branch_coverage=1 --title "LIBFLAME CODE COVERAGE REPORT" filtered_coverage.info --prefix $PWD --function-coverage --branch-coverage --legend --output-directory out;
cd out; pushd &lt;index.html;  python3 -m http.server 9999; popd;