/* This gives the effect of

	subroutine exit(rc)
	integer*4 rc
	stop
	end

 * with the added side effect of supplying rc as the program's exit code.
 */

#include "FLA_f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

  void exit_(integer *rc)
  {
    exit(*rc);
  }

#ifdef __cplusplus
}
#endif
