
@echo off

echo. 
echo  Makefile
echo. 
echo  Field G. Van Zee
echo.  
echo  nmake Makefile for building libflame for Microsoft Windows. nmake targets
echo  may be invoked after running the configure.cmd script. Valid targets are:
echo. 
echo    all          - Invoke the lib and dll targets.
echo    lib          - Build libflame as a static library.
echo    dll          - Build libflame as a dynamically-linked library.
echo    help         - Output help and usage information.
echo    clean        - Invoke clean-log and clean-build targets.
echo    clean-log    - Remove any log files present.
echo    clean-config - Remove all products of configure.cmd. Namely, remove the
echo                   config, include, and src directories.
echo    clean-build  - Remove all products of the compilation portion of the build
echo                   process. Namely, remove the obj, lib, and dll directories.
echo    distclean    - Invoke clean-log, clean-config, and clean-build targets.
echo.
echo  The Makefile also recognizes configuration options corresponding to the
echo  following Makefile variables:
echo.
echo    VERBOSE               - When defined, nmake outputs the actual commands
echo                            executed instead of more concise one-line progress
echo                            indicators. (Undefined by default.)
echo.
echo  Typically, these options are specified by commenting or uncommenting the
echo  corresponding lines in the Makefile. However, if the Makefile currently does
echo  not define one of the options, and you wish to enable the corresponding
echo  feature without editing the Makefile, you may define the variable at the
echo  command line when nmake is invoked. For example, you may enable verboseness
echo  while invoking the lib target as follows:
echo.
echo    nmake lib VERBOSE=1
echo.
