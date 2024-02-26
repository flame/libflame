AOCL Progress

The AOCL libraries may be used to perform lengthy computations 
(for example, matrix multiplications, solver involving large matrices). 
These operations/computations may go on for hours. 
AOCL progress feature provides mechanism for the application to check how far the computations have progressed.
The AOCL libraries periodically updates the application with progress made so far via a callback function.

Usages:

The Application needs to define a aocl_fla_progress function or callback function in specific format 
and register this callback function with the Libflame library.

int aocl_fla_progress(
char* api,
integer lenapi,
integer *progress,
integer *current_thread,
integer *total_threads
)

or

The callback function prototype must be as 
defined below, however, the function name can be changed as per user preference.


int test_progress(
char* api,
integer lenapi,
integer *progress,
integer *current_thread,
integer *total_threads
)

The table below explains various parameters
parameters           |        Purpose
---------------------------------------------------------------------
api                  |   Name of the API which is currently running
lapi                 |   Length of API/Operation string
progress	           |   Linear progress made in current thread so far
current_thread	     |   Current thread id
total_threads	       |   Total number of threads in current team

Note: In case of single-threaded AOCL-LAPACK, the values of "current_thread" and "total_threads" are set to 0 and 1 respectively.
As a result, the callback function cannot be used to monitor the thread ID and thread count of the application.

Callback Registration: 

The callback function must be registered with library for it to report the progress. Each library has its own callback registration function. Registration is done by calling:
aocl_fla_set_progress(test_progress);   // for libflame

Example:

int aocl_fla_progress(char* api,integer lenapi,integer *progress,integer *current_thread,integer *total_threads)
{
  char buf[BUFLEN];
  if( lenapi >= BUFLEN ) lenapi = BUFLEN-1;
  strncpy( buf, api, lenapi );
  buf[lenapi] = '\0';
  printf( "In AOCL FLA  Progress thread  %lld", at API  %s, progress  %lld total threads= %lld\n", *current_thread, buf, *progress,*total_threads );
  return 0;

}

or

int test_progress(char* api,integer lenapi,integer *progress,integer *current_thread,integer *total_threads)
{
  char buf[BUFLEN];
  if( lenapi >= BUFLEN ) lenapi = BUFLEN-1;
  strncpy( buf, api, lenapi );
  buf[lenapi] = '\0';
  printf( "In AOCL Progress thread  %lld", at API  %s, progress  %lld total threads= %lld\n", *current_thread, buf, *progress,*total_threads );
  return 0;

}

Register the callback with:
aocl_fla_set_progress(test_progress);

/* NOTE- Application must avoid updating 'lenapi' value and ensure it is within the string length of the API name. */

This will results in output in following format (output truncated):

/libflame/test/aoclflaprogress$ ./test_libFLAME_aocl.x
in main m = 8006 ,n=8006,lda=8006
In AOCL Progress thread  0, at API  DGETRF, progress 8 total threads= 1
In AOCL Progress thread  0, at API  DGETRF, progress 16 total threads= 1
In AOCL Progress thread  0, at API  DGETRF, progress 24 total threads= 1
In AOCL Progress thread  0, at API  DGETRF, progress 32 total threads= 1
In AOCL Progress thread  0, at API  DGETRF, progress 40 total threads= 1
In AOCL Progress thread  0, at API  DGETRF, progress 48 total threads= 1
In AOCL Progress thread  0, at API  DGETRF, progress 56 total threads= 1
In AOCL Progress thread  0, at API  DGETRF, progress 64 total threads= 1
In AOCL Progress thread  0, at API  DGETRF, progress 72 total threads= 1
In AOCL Progress thread  0, at API  DGETRF, progress 80 total threads= 1
In AOCL Progress thread  0, at API  DGETRF, progress 88 total threads= 1
In AOCL Progress thread  0, at API  DGETRF, progress 96 total threads= 1
In AOCL Progress thread  0, at API  DGETRF, progress 104 total threads= 1
In AOCL Progress thread  0, at API  DGETRF, progress 112 total threads= 1
In AOCL Progress thread  0, at API  DGETRF, progress 120 total threads= 1
In AOCL Progress thread  0, at API  DGETRF, progress 128 total threads= 1
In AOCL Progress thread  0, at API  DGETRF, progress 136 total threads= 1
In AOCL Progress thread  0, at API  DGETRF, progress 144 total threads= 1
In AOCL Progress thread  0, at API  DGETRF, progress 152 total threads= 1
In AOCL Progress thread  0, at API  DGETRF, progress 160 total threads= 1
....

Limitations And Caveats:

    On Windows,aocl_fla_progress is not supported. Use Callback registration method.
