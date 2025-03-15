#include "mex.h"
#include "flops.h"

void addflops(unsigned fl)
{
  mxArray *flopcount = mexGetVariable("global","flopcount");
  if(flopcount && !mxIsEmpty(flopcount)) {
    *mxGetPr(flopcount) += fl;
    mexPutVariable("global","flopcount",flopcount);
  }
#endif
}


