
@echo off

:ENVIRONMENT
	set GEN_CHECK_REV_FILE=.\build\gen-check-rev-file.py
	set GATHER_SRC=.\build\gather-src-for-windows.py
	set GEN_CONFIG_FILE=.\build\gen-config-file.py
	set CONFIG_DEFS_TEMPL=.\build\config.mk.in
	set SRC_TREE_DIR=..\src
	set TOP_BUILD_DIR=.

:PARAMS
	if "%1"=="" (goto USAGE)
	if "%2"=="" (goto USAGE)
	if "%3"=="" (goto USAGE)

	set ARCH=%1
	set BUILD=%2
	set CCOMPILER=%3
	
:TASK_UNIT
	echo %0: Checking/updating revision file.
	%GEN_CHECK_REV_FILE% -v
	echo %0: Gathering source files into local flat directories.
	%GATHER_SRC% %SRC_TREE_DIR% %TOP_BUILD_DIR%
	echo %0: Creating configure definitions file.
	%GEN_CONFIG_FILE% %TOP_BUILD_DIR% %ARCH% %BUILD% %CCOMPILER% %CONFIG_DEFS_TEMPL%
	echo %0: Configuration and setup complete. You may now run nmake. 

	goto END

:USAGE
	echo. 
	echo  configure.cmd
	echo. 
	echo  A wrapper script for various configuration and setup scripts that need
	echo. to be run before nmake when building libflame for Microsoft Windows.
	echo. 
	echo  USAGE:
	echo     %0 [arch] [build] [cc]
	echo.
	echo        arch     -- The architecture string to build.
	echo                    Supported values: {x86,x64}
	echo        build    -- The kind of build.
	echo                    Supported values: {debug,release}
	echo        cc       -- The C compiler to use.
	echo                    Supported values: {icl,cl}
	echo. 
	echo  examples:
	echo     %0 x86 debug icl
	echo     %0 x64 release cl
	echo.

:END
