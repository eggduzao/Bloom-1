#include "mex.h"
#include "util.h"

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  mwSize len, i, bins;
  double *indata, *outdata;

  if((nrhs < 1) || (nrhs > 2))
    mexErrMsgTxt("Usage: h = int_hist(x, n)");

  indata = mxGetPr(prhs[0]);
  len = mxGetNumberOfElements(prhs[0]);

  if(mxIsSparse(prhs[0]))
    mexErrMsgTxt("Cannot handle sparse matrices.  Sorry.");

  if(nrhs == 2) {
    if(mxGetNumberOfElements(prhs[1]) != 1) mexErrMsgTxt("n is not scalar.");
    bins = *mxGetPr(prhs[1]);
  } else {
    bins = indata[0];
    for(i=0;i<len;i++) {
      if(indata[i] > bins) bins = indata[i];
    }
  }

  plhs[0] = mxCreateDoubleMatrix(1, bins, mxREAL);
  outdata = mxGetPr(plhs[0]);

  for(i=0;i<len;i++) {
    int v = (int)(*indata++) - 1;
    if((v < 0) || (v >= bins))
      mexErrMsgTxt("value out of bounds");
    outdata[v]++;
  }
}


