#include "test_lapack2flame.h"

#define FLOAT_NAME    "single"
#define DOUBLE_NAME   "double"
#define SCOMPLEX_NAME "scomplex"
#define DCOMPLEX_NAME "dcomplex"
#define INTEGER_NAME  "integer"

char* datatype2str( FLA_Datatype datatype )
{
    switch (datatype)
    {
    case FLA_FLOAT:
        return FLOAT_NAME;
    case FLA_DOUBLE:
        return DOUBLE_NAME;
    case FLA_COMPLEX:
        return SCOMPLEX_NAME;
    case FLA_DOUBLE_COMPLEX:
        return DCOMPLEX_NAME;
    case FLA_INT:
        return INTEGER_NAME;
    default:
        fprintf(stderr, ">>>> Datatype: %d\n", datatype);
    }
    return "Underined datatype";
}

FLA_Datatype str2datatype( char* name )
{
    char *p = name;

    for ( ; *p; ++p) *p = tolower(*p);

    if      (!strcmp(name, FLOAT_NAME))    return FLA_FLOAT;
    else if (!strcmp(name, DOUBLE_NAME))   return FLA_DOUBLE;
    else if (!strcmp(name, SCOMPLEX_NAME)) return FLA_COMPLEX;
    else if (!strcmp(name, DCOMPLEX_NAME)) return FLA_DOUBLE_COMPLEX;
    else if (!strcmp(name, INTEGER_NAME))  return FLA_INT;
    return 0;
}

/* char* transpose2str( FLA_Trans trans ) */
/* { */
/*     switch (datatype) */
/*     { */
/*     case FLA_FLOAT: */
/*         return FLOAT_NAME; */
/*     case FLA_DOUBLE: */
/*         return DOUBLE_NAME; */
/*     case FLA_COMPLEX: */
/*         return SCOMPLEX_NAME; */
/*     case FLA_DOUBLE_COMPLEX: */
/*         return DCOMPLEX_NAME; */
/*     case FLA_INT: */
/*         return INTEGER_NAME; */
/*     default: */
/*         fprintf(stderr, ">>>> Datatype: %d\n", datatype); */
/*     } */
/*     return "Underined datatype"; */
/* } */
