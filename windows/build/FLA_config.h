
// --- General options ---------------------------------------------------------

// Determines whether to enable various segments of code identified as
// providing non-critical functionality.
#define FLA_ENABLE_NON_CRITICAL_CODE

// Determines whether the LAPACK compatibility layer is included in libflame.
// NOTE: If lapack2flame is enabled, external-lapack-for-subproblems MUST
// be disabled!
//#define FLA_ENABLE_LAPACK2FLAME

// Determines whether to enable external LAPACK for small subproblems.
// NOTE: If external-lapack-for-subproblems is enabled, (a) lapack2flame MUST
// be disabled, AND (b) external-lapack-interfaces MUST be enabled.
//#define FLA_ENABLE_EXTERNAL_LAPACK_FOR_SUBPROBLEMS

// Determines whether to enable interfaces to external LAPACK routines.
// NOTE: If external-lapack-interfaces is enabled, an LAPACK library will be
// required at link-time.
//#define FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES

// Determines whether to use control trees to select a reasonable FLAME
// variant and blocksize when level-3 BLAS front-ends are invoked.
#define FLA_ENABLE_BLAS3_FRONT_END_CNTL_TREES

// Determines whether to modify various segments of code needed for
// integrating libflame into Windows.
#define FLA_ENABLE_WINDOWS_BUILD

// Determines whether to define a portable FLA_Clock() in terms of
// gettimeofday() from time.h.
#define FLA_ENABLE_PORTABLE_TIMER


// --- Runtime error checking and debugging ------------------------------------

// Determines whether to enable internal runtime consistency checks of
// function parameters and return values.
#define FLA_ENABLE_INTERNAL_ERROR_CHECKING

// Encodes the default level of internal error checking chosen at
// configure-time.
#define FLA_INTERNAL_ERROR_CHECKING_LEVEL 2   // full error checking
//#define FLA_INTERNAL_ERROR_CHECKING_LEVEL 1   // minimal error checking
//#define FLA_INTERNAL_ERROR_CHECKING_LEVEL 0   // no error checking

// Determines whether to enable the FLA_malloc()/FLA_free() memory counter
// by default.
#define FLA_ENABLE_MEMORY_LEAK_COUNTER


// --- Multithreading and SuperMatrix ------------------------------------------

// Determines whether thread-specific blocks of code should be compiled.
#define FLA_ENABLE_MULTITHREADING

// Encodes the type of multithreading chosen at configure-time.
#define FLA_MULTITHREADING_MODEL 1   // OpenMP
// #define FLA_MULTITHREADING_MODEL 2   // POSIX threads

// Determines whether SuperMatrix-specific blocks of code should be compiled.
#define FLA_ENABLE_SUPERMATRIX


// --- BLAS and blocksizes -----------------------------------------------------

// Determines whether to enable CBLAS interfaces instead of Fortran-77
// interfaces to the BLAS.
// #define FLA_ENABLE_CBLAS_INTERFACES

// Determines whether to enable interfaces to internal/low-level libgoto
// symbols.
// #define FLA_ENABLE_GOTO_INTERFACES

// Sets the default blocksize in the k dimension (used only if
// libgoto interfaces are disabled).
// #define FLA_DEFAULT_K_BLOCKSIZE 128

// Sets the default blocksize in the m dimension (used only if
// libgoto interfaces are disabled).
// #define FLA_DEFAULT_M_BLOCKSIZE 128

// Sets the default blocksize in the n dimension (used only if
// libgoto interfaces are disabled).
// #define FLA_DEFAULT_N_BLOCKSIZE 128


// --- Memory alignment --------------------------------------------------------

// Determines whether memory is aligned to user-requested boundaries.
// #define FLA_ENABLE_MEMORY_ALIGNMENT

// Sets the byte boundary used to align the starting address of all memory
// allocated dynamically through libflame. Only used if
// FLA_ENABLE_MEMORY_ALIGNMENT is defined.
// #define FLA_MEMORY_ALIGNMENT_BOUNDARY 8

// Determines whether to enable code that will increase FLA_Obj leading
// dimensions to ensure that matrix columns adhere to the alignment specified
// by FLA_MEMORY_ALIGNMENT_BOUNDARY.
// #define FLA_ENABLE_LDIM_ALIGNMENT


// --- Fortran-77 compatibility ------------------------------------------------

// Determines whether the Fortran name-mangling suffix was determined at
// configure-time. This option is not used in Windows.
// #define FLA_ENABLE_AUTODETECT_F77_UNDERSCORING

// Determines whether the Fortran 77 compiler appends an underscore to symbol
// names. Not used in Windows.
// #define FLA_F77_UNDERSCORE

// Determines whether the Fortran 77 compiler appends an extra underscore to
// symbol names that already contain at least one underscore. Not used in
// Windows.
// #define FLA_F77_EXTRA_UNDERSCORE

// Determines whether invocations to the BLAS within libflame are converted to
// uppercase symbols.
#define FLA_ENABLE_UPPERCASE_BLAS

// Determines whether invocations to LAPACK within libflame are converted to
// uppercase symbols.
#define FLA_ENABLE_UPPERCASE_LAPACK


// --- Experimental/unsupported/broken options ---------------------------------

// Determines whether GPU-specific blocks of code should be compiled.
// #define FLA_ENABLE_GPU

