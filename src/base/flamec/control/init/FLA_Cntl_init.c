
#include "FLAME.h"

void FLA_Cntl_init()
{
  FLA_Cntl_init_flamec();
  FLA_Cntl_init_flash();
}

void FLA_Cntl_finalize()
{
  FLA_Cntl_finalize_flamec();
  FLA_Cntl_finalize_flash();
}

