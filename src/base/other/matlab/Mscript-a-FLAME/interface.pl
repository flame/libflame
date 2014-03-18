#!/lusr/bin/perl

# Usage:
# 1) interface.pl
#      FLA_foo nrhs nint
# 2) interface.pl < funcList.txt
#
# This script is used to create interface mex files that enable FLAME computational
# routines written in C to be called from Matlab.  For each FLAME function, the
# following things must be specified.
# FLA_foo       Name of the FLAME function being interfaced.
# nrhs          Number of arguments.
# nint          Number of integer (attribute) arguments.
#
# The script produces a file named FLA_foo.c, that contains the inteface code for
# each function FLA_foo that is specified. The above parameters needed for each
# function can be specified interactively or listed in a file, funcList.txt.  Each
# function must be listed on its own line:
# FLA_foo1 nrhs1 nint1
# FLA_foo2 nrhs2 nint2
# ...
#

while (<>) {
  ($func, $nrhs, $nint) = /^\s*(\S+)\s+(\d+)\s+(\d+)\s*$/;

  next if $func eq "";

  if ($nint > $nrhs) {
    die "$func: nint cannot be greater than nrhs";
  }

  # For each function create a C interface file.
  open(FH, ">$func.c") || die "Cannot open file for output.";

  printf FH "/* Mex file created by interface.pl for the FLAME C function %s() */\n", $func;
  printf FH "\n";
  printf FH "#include \"mex.h\"\n";
  printf FH "#include \"FLA_Matlab2C.h\"\n";
  printf FH "\n";
  printf FH "extern double FLA_Clock();\n";
  printf FH "\n";
  printf FH "#define NRHS %d\n", $nrhs;
  printf FH "#define NINT %d\n", $nint;
  printf FH "\n";
  printf FH "#define NOBJ (NRHS-NINT)\n";
  printf FH "\n";
  printf FH "\n";
  printf FH "void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {\n";
  printf FH "  int attr[NINT];\n";
  printf FH "  FLA_Obj obj[NOBJ];\n";
  printf FH "  double *dtime;\n";
  printf FH "\n";
  printf FH "  FLA_Init();\n";
  printf FH "\n";
  printf FH "  /* Check if the number of arguments supplied is correct */\n";
  printf FH "  FLA_M2C_CheckNumArgs(NRHS, nrhs);\n";
  printf FH "\n";
  printf FH "  /* Convert Matlab arguments into the appropriate FLAME C arguments */\n";
  printf FH "  FLA_M2C_ConvertArgs(NRHS, prhs, NINT, attr, obj);\n";
  printf FH "\n";
  printf FH "  /* If an extra argument is supplied, collect timing informaion in it. */\n";
  printf FH "  if (nrhs == NRHS+1)\n";
  printf FH "    dtime = FLA_M2C_ConvertDoublePtr(prhs[NRHS]);\n";
  printf FH "\n";

  printf FH "  /* Now call the C FLAME function, timing it if the extra argument is given. */\n";
  printf FH "\n";
  printf FH "  if (nrhs == NRHS+1)\n";
  printf FH "    *dtime = FLA_Clock();\n";
  printf FH "\n";

  printf FH "  %s(", $func;
  for ($i=0; $i<$nint; $i++) {
    printf FH "attr[%d]", $i;
    if ($i < ($nrhs-1)) {
      printf FH ", ";
    }
  }
  for ($i=$i; $i<$nrhs; $i++) {
    printf FH "obj[%d]", $i-$nint;
    if ($i < ($nrhs-1)) {
      printf FH ", ";
    }
  }
  printf FH ");\n";

  printf FH "\n";
  printf FH "  if (nrhs == NRHS+1)\n";
  printf FH "    *dtime = FLA_Clock() - *dtime;\n";
  printf FH "\n";
  printf FH "  FLA_Finalize();\n";
  printf FH "\n";
  printf FH "}\n";
  printf FH "\n";

  close(FH);
}
