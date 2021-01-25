#include "mexutil.h"

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  int m,n,i,sz;
  double *p;
  if(nrhs != 2) mexErrMsgTxt("xones(m,n)");
  if(mxGetNumberOfElements(prhs[0]) != 1)
    mexErrMsgTxt("m is not scalar");
  if(mxGetNumberOfElements(prhs[1]) != 1)
    mexErrMsgTxt("n is not scalar");
  m = (int)*mxGetPr(prhs[0]);
  n = (int)*mxGetPr(prhs[1]);
  sz = m*n;
  plhs[0] = mxCreateDoubleMatrixE(m,n,mxREAL);
  p = (double*)mxGetPr(plhs[0]);
#endif
  for(i=0;i<sz;i++) *p++ = 1;
}


