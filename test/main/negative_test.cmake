set(NEGATIVE_TEST_CASES "ggevx sdcz A V V B 10 10 10 10 10 -1 1 --einfo=-1"
            "ggevx sdcz B A V B 10 10 10 10 10 -1 1 --einfo=-2"
            "ggevx sdcz B V A B 10 10 10 10 10 -1 1 --einfo=-3"
            "ggevx sdcz B V V A 10 10 10 10 10 -1 1 --einfo=-4"
            "ggevx sdcz B V V B -10 10 10 10 10 -1 1 --einfo=-5"
            "ggevx sdcz B V V B 10 -10 10 10 10 -1 1 --einfo=-7"
            "ggevx sdcz B V V B 10 10 -10 10 10 -1 1 --einfo=-9"
            "ggevx sd B V V B 10 10 10 -10 10 -1 1 --einfo=-14"
            "ggevx cz B V V B 10 10 10 -10 10 -1 1 --einfo=-13"
            "ggevx sd B V V B 10 10 10 10 -10 -1 1 --einfo=-16"
            "ggevx cz B V V B 10 10 10 10 -10 -1 1 --einfo=-15"
            "gesv sdcz -10 10 10 10 1 --einfo=-1"
            "gesv sdcz 10 -10 10 10 1 --einfo=-2"
            "gesv sdcz 10 10 -10 10 1 --einfo=-4"
            "gesv sdcz 10 10 10 -10 1 --einfo=-7"
            "geqrf sdcz -10 10 10 -1 1 --einfo=-1"
            "geqrf sdcz 10 -10 10 -1 1 --einfo=-2"
            "geqrf sdcz 10 10 -10 -1 1 --einfo=-4"
            "geqrf sdcz 10 10 10 1 1 --einfo=-7"
            "gerqf sdcz -10 10 10 -1 1 --einfo=-1"
            "gerqf sdcz 10 -10 10 -1 1 --einfo=-2"
            "gerqf sdcz 10 10 -10 -1 1 --einfo=-4"
            "gerqf sdcz 10 10 10 1 1 --einfo=-7"
            "gerq2 sdcz -10 10 10 1 --einfo=-1"
            "gerq2 sdcz 10 -10 10 1 --einfo=-2"
            "gerq2 sdcz 10 10 -10 1 --einfo=-4"
            "gelqf sdcz -10 10 10 -1 1 --einfo=-1"
            "gelqf sdcz 10 -10 10 -1 1 --einfo=-2"
            "gelqf sdcz 10 10 -10 -1 1 --einfo=-4"
            "gelqf sdcz 10 10 10 1 1 --einfo=-7"
            "potrf sdcz U -10 10 1 --einfo=-2"
            "potrf sdcz X 10 10 1 --einfo=-1"
            "potrf sdcz U 10 -10 1 --einfo=-4"
            "getrf sdcz -10 10 10 1 --einfo=-1"
            "getrf sdcz 10 -10 10 1 --einfo=-2"
            "getrf sdcz 10 10 -10 1 --einfo=-4"
            "getri sdcz -10 10 -1 1 --einfo=-1"
            "getri sdcz 10 -10 -1 1 --einfo=-3"
            "getri sdcz 10 10 1 1 --einfo=-6"
            "getrs sdcz A 10 10 10 10 1 --einfo=-1"
            "getrs sdcz C -10 10 10 10 1 --einfo=-2"
            "getrs sdcz C 10 -10 10 10 1 --einfo=-3"
            "getrs sdcz C 10 10 -10 10 1 --einfo=-5"
            "getrs sdcz C 10 10 10 -10 1 --einfo=-8"
            "potrs sdcz A 10 10 10 10 1 --einfo=-1"
            "potrs sdcz U -10 10 10 10 1 --einfo=-2"
            "potrs sdcz U 10 -10 10 10 1 --einfo=-3"
            "potrs sdcz U 10 10 -10 10 1 --einfo=-5"
            "potrs sdcz U 10 10 10 -10 1 --einfo=-7"
            "orgqr sdcz 10 -10 10 -1 1 --einfo=-2"
            "orgqr sdcz 10 10 -10 -1 1 --einfo=-5"
            "orgqr sdcz 10 10 10 1 1 --einfo=-8"
            "gesdd sdcz A -10 10 10 10 10 -1 1 --einfo=-2"
            "gesdd sdcz A 10 -10 10 10 10 -1 1 --einfo=-3"
            "gesdd sdcz P 10 10 10 10 10 -1 1 --einfo=-1"
            "gesdd sdcz A 10 10 -10 10 10 -1 1 --einfo=-5"
            "gesdd sdcz A 10 10 10 -10 10 -1 1 --einfo=-8"
            "gesdd sdcz A 10 10 10 10 -10 -1 1 --einfo=-10"
            "gesdd sdcz A 10 10 10 10 10 1 1 --einfo=-12"
            "syevd sdcz V U -10 10 -1 -1 -1 1 --einfo=-3"
            "syevd sdcz A U 10 10 -1 -1 -1 1 --einfo=-1"
            "syevd sdcz V A 10 10 -1 -1 -1 1 --einfo=-2"
            "syevd sdcz V U 10 -10 -1 -1 -1 1 --einfo=-5"
            "syevd sd V U 10 10 1 1 -1 1 --einfo=-8"
            "syevd sd V U 10 10 1000 1 -1 1 --einfo=-10"
            "syevd cz V U 10 10 1 1 1 1 --einfo=-8"
            "syevd cz V U 10 10 1000 1 1 1 --einfo=-10"
            "syevd cz V U 10 10 1000 1 1000 1 --einfo=-12"
            "gesvd sdcz X A 10 10 10 10 10 -1 1 --einfo=-1"
            "gesvd sdcz A X 10 10 10 10 10 -1 1 --einfo=-2"
            "gesvd sdcz A A -10 10 10 10 10 -1 1 --einfo=-3"
            "gesvd sdcz A A 10 -10 10 10 10 -1 1 --einfo=-4"
            "gesvd sdcz A A 10 10 -10 10 10 -1 1 --einfo=-6"
            "gesvd sdcz A A 10 10 10 -10 10 -1 1 --einfo=-9"
            "gesvd sdcz A A 10 10 10 10 -10 -1 1 --einfo=-11"
            "gesvd sdcz A A 10 10 10 10 10 1 1 --einfo=-13"
            "geevx sdcz A V V B 10 10 10 10 -1 1 --einfo=-1"
            "geevx sdcz B A V V 10 10 10 10 -1 1 --einfo=-2"
            "geevx sdcz B V A V 10 10 10 10 -1 1 --einfo=-3"
            "geevx sdcz B V V A 10 10 10 10 -1 1 --einfo=-4"
            "geevx sdcz B V V B -10 10 10 10 -1 1 --einfo=-5"
            "geevx sdcz B V V B 10 -10 10 10 -1 1 --einfo=-7"
            "geevx sd B V V B 10 10 -10 10 -1 1 --einfo=-11"
            "geevx cz B V V B 10 10 -10 10 -1 1 --einfo=-10"
            "geevx sd B V V B 10 10 10 -10 -1 1 --einfo=-13"
            "geevx cz B V V B 10 10 10 -10 -1 1 --einfo=-12"
            "geevx sd B V V B 10 10 10 10 1 1 --einfo=-21"
            "geevx cz B V V B 10 10 10 10 1 1 --einfo=-20"
            "geev sdcz A V 10 10 10 10 -1 1 --einfo=-1"
            "geev sdcz V A 10 10 10 10 -1 1 --einfo=-2"
            "geev sdcz V V -10 10 10 10 -1 1 --einfo=-3"
            "geev sdcz V V 10 -10 10 10 -1 1 --einfo=-5"
            "geev sd V V 10 10 -10 10 -1 1 --einfo=-9"
            "geev cz V V 10 10 -10 10 -1 1 --einfo=-8"
            "geev sd V V 10 10 10 -10 -1 1 --einfo=-11"
            "geev cz V V 10 10 10 -10 -1 1 --einfo=-10"
            "geev sd V V 10 10 10 10 1 1 --einfo=-13"
            "geev cz V V 10 10 10 10 1 1 --einfo=-12"
            "geqp3 sdcz -10 10 10 -1 1 --einfo=-1"
            "geqp3 sdcz 10 -10 10 -1 1 --einfo=-2"
            "geqp3 sdcz 10 10 -10 -1 1 --einfo=-4"
            "geqp3 sdcz 10 10 10 1 1 --einfo=-8"
            "ggev sdcz V V -10 10 10 10 10 -1 1 --einfo=-3"
            "ggev sdcz V V 10 -10 10 10 10 -1 1 --einfo=-5"
            "ggev sd V V 10 10 10 -10 10 -1 1 --einfo=-12"
            "ggev cz V V 10 10 10 -10 10 -1 1 --einfo=-11"
            "ggev sd V V 10 10 10 10 -10 -1 1 --einfo=-14"
            "ggev cz V V 10 10 10 10 -10 -1 1 --einfo=-13"
            "ggev sd V V 10 10 10 10 10 1 1 --einfo=-16"
            "ggev cz V V 10 10 10 10 10 1 1 --einfo=-15"
            "steqr sdcz A 10 10 1 --einfo=-1"
            "steqr sdcz V -10 10 1 --einfo=-2"
            "steqr sdcz V 10 -10 1 --einfo=-6"
            "stevd sd A 10 10 -1 -1 1 --einfo=-1"
            "stevd sd V -10 10 -1 -1 1 --einfo=-2"
            "stevd sd V 10 -10 -1 -1 1 --einfo=-6"
            "stevd sd V 10 10 1 1 1 --einfo=-8"
            "stevd sd V 10 10 1000 1 1 --einfo=-10"
            "stedc sdcz A 10 10 -1 -1 -1 1 --einfo=-1"
            "stedc sdcz V -10 10 -1 -1 -1 1 --einfo=-2"
            "stedc sdcz V 10 -10 -1 -1 -1 1 --einfo=-6"
            "stedc sd V 10 10 1 1 -1 1 --einfo=-8"
            "stedc cz V 10 10 1000 1 1 1 --einfo=-10"
            "hseqr sdcz V V 10 2 5 10 10 -1 1 --einfo=-1"
            "hseqr sdcz S S 10 2 5 10 10 -1 1 --einfo=-2"
            "hseqr sdcz S V -10 2 5 10 10 -1 1 --einfo=-3"
            "hseqr sdcz S V 10 11 5 10 10 -1 1 --einfo=-4"
            "hseqr sdcz S V 10 2 12 10 10 -1 1 --einfo=-5"
            "hseqr sdcz S V 10 2 5 -10 10 -1 1 --einfo=-7"
            "hseqr sd S V 10 2 5 10 -10 -1 1 --einfo=-11"
            "hseqr cz S V 10 2 5 10 -10 -1 1 --einfo=-10"
            "hseqr sd S V 10 2 5 10 10 1 1 --einfo=-13"
            "hseqr cz S V 10 2 5 10 10 1 1 --einfo=-12"
            "syev sdcz A U 10 10 -1 1 --einfo=-1"
            "syev sdcz V A 10 10 -1 1 --einfo=-2"
            "syev sdcz V U -10 10 -1 1 --einfo=-3"
            "syev sdcz V U 10 -10 -1 1 --einfo=-5"
            "syev sdcz V U 10 10 1 1 --einfo=-8"
            "gehrd sdcz -10 2 5 10 -1 1 --einfo=-1"
            "gehrd sdcz 10 11 5 10 -1 1 --einfo=-2"
            "gehrd sdcz 10 2 12 10 -1 1 --einfo=-3"
            "gehrd sdcz 10 2 5 -10 -1 1 --einfo=-5"
            "gehrd sdcz 10 2 5 10 1 1 --einfo=-8"
            "gghrd sdcz A V 10 2 5 10 10 10 10 1 --einfo=-1"
            "gghrd sdcz V A 10 2 5 10 10 10 10 1 --einfo=-2"
            "gghrd sdcz V V -10 2 5 10 10 10 10 1 --einfo=-3"
            "gghrd sdcz V V 10 0 5 10 10 10 10 1 --einfo=-4"
            "gghrd sdcz V V 10 2 11 10 10 10 10 1 --einfo=-5"
            "gghrd sdcz V V 10 2 5 -10 10 10 10 1 --einfo=-7"
            "gghrd sdcz V V 10 2 5 10 -10 10 10 1 --einfo=-9"
            "gghrd sdcz V V 10 2 5 10 10 10 -10 1 --einfo=-13"
            "hgeqz sdcz A V V 10 2 4 10 10 10 10 -1 1 --einfo=-1"
            "hgeqz sdcz S A V 10 2 4 10 10 10 10 -1 1 --einfo=-2"
            "hgeqz sdcz S V V -10 2 4 10 10 10 10 -1 1 --einfo=-4"
            "hgeqz sdcz S V V 10 0 4 10 10 10 10 -1 1 --einfo=-5"
            "hgeqz sdcz S V V 10 2 12 10 10 10 10 -1 1 --einfo=-6"
            "hgeqz sdcz S V V 10 2 4 -10 10 10 10 -1 1 --einfo=-8"
            "hgeqz sdcz S V V 10 2 4 10 -10 10 10 -1 1 --einfo=-10"
            "hgeqz sd S V V 10 2 4 10 10 -10 10 -1 1 --einfo=-15"
            "hgeqz cz S V V 10 2 4 10 10 -10 10 -1 1 --einfo=-14"
            "hgeqz sd S V V 10 2 4 10 10 10 -10 -1 1 --einfo=-17"
            "hgeqz cz S V V 10 2 4 10 10 10 -10 -1 1 --einfo=-16"
            "hgeqz sd S V V 10 2 4 10 10 10 10 1 1 --einfo=-19"
            "hgeqz cz S V V 10 2 4 10 10 10 10 1 1 --einfo=-18"
            "org2r sdcz 10 -10 10 1 --einfo=-2"
            "org2r sdcz 10 10 -10 1 --einfo=-5")

set(TEST_NUM 1)
foreach(neg_test_cases IN LISTS NEGATIVE_TEST_CASES)
    # this line splits entire string into separate arguments as ctest requres the arguments to be passed separately rather than a single string
    string(REPLACE " " ";" COMMANDLINE_PARAMS ${neg_test_cases})
    set(TEST_NAME NEGATIVE_TEST_CASE_${TEST_NUM} )
    add_test(${TEST_NAME} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${PROJECT_NAME} ${COMMANDLINE_PARAMS})
    set_tests_properties(${TEST_NAME} PROPERTIES FAIL_REGULAR_EXPRESSION "FAIL;No test was run, give valid arguments")
MATH(EXPR TEST_NUM "${TEST_NUM}+1")
endforeach()
