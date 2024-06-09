/* Include files */

#include "untitled_cgxe.h"
#include "m_YOartArQux2lGbKyl3KxRB.h"
#include "m_0UplWFobhjrEmZfPJtxyI.h"

unsigned int cgxe_untitled_method_dispatcher(SimStruct* S, int_T method, void
  * data)
{
  if (ssGetChecksum0(S) == 2441808678 &&
      ssGetChecksum1(S) == 1256922200 &&
      ssGetChecksum2(S) == 2201851376 &&
      ssGetChecksum3(S) == 527969350) {
    method_dispatcher_YOartArQux2lGbKyl3KxRB(S, method, data);
    return 1;
  }

  if (ssGetChecksum0(S) == 3811608255 &&
      ssGetChecksum1(S) == 1569159052 &&
      ssGetChecksum2(S) == 22600739 &&
      ssGetChecksum3(S) == 2271534791) {
    method_dispatcher_0UplWFobhjrEmZfPJtxyI(S, method, data);
    return 1;
  }

  return 0;
}
