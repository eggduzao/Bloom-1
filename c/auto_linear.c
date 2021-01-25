#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], 
		 int nrhs, const mxArray *prhs[])
{
  double *p;
  if(nrhs != 2) {
    mexErrMsgTxt("usage: sameobject(a,b)");
  }
  plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
  p = mxGetPr(plhs[0]);
  *p = (mxGetData(prhs[0]) == mxGetData(prhs[1]));
}


