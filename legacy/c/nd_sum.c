#include "mex.h"
#include "flops.h"

void ndsum(double *dest, double *src, mwSize ndim, mwSize *size,
	   mwSize total, mwSize *mask) 
{

  mwSize *masked_size,*cum_masked_size,*advance,*rewind;
  mwSize *subs = mxCalloc(ndim,sizeof(mwSize));
  mwSize i,j;
  masked_size = mxCalloc(ndim,sizeof(mwSize));
  for(i=0;i<ndim;i++) {
    masked_size[i] = mask[i] ? 1 : size[i];
  }

  cum_masked_size = mxCalloc(ndim,sizeof(mwSize));
  advance = mxCalloc(ndim,sizeof(mwSize));
  rewind = mxCalloc(ndim,sizeof(mwSize));
  cum_masked_size[0] = 1;
  for(i=0;i<ndim-1;i++) {
    cum_masked_size[i+1] = cum_masked_size[i]*masked_size[i];
  }
  for(i=0;i<ndim;i++) {
    advance[i] = cum_masked_size[i]*!mask[i];
    rewind[i] = cum_masked_size[i]*(masked_size[i]-1);
  }

  for(j=0;j<total;j++) {
    *dest += *src++;
    for(i=0;i<ndim;i++) {
      if(subs[i] == size[i]-1) {
	subs[i] = 0;
	dest -= rewind[i];
      }
      else {
	subs[i]++;
	dest += advance[i];
	break;
      }
    }
  }
  mxFree(rewind);
  mxFree(advance);
  mxFree(cum_masked_size);
  mxFree(masked_size);
  mxFree(subs);
}

void mexFunction(int nlhs, mxArray *plhs[], 
		 int nrhs, const mxArray *prhs[])
{
  mwSize i,j;
  mwSize ndim = mxGetNumberOfDimensions(prhs[0]);
  const mwSize *size = mxGetDimensions(prhs[0]);
  double *src = mxGetPr(prhs[0]), *dest;
  double *dim = mxGetPr(prhs[1]);
  mwSize nremove = mxGetNumberOfElements(prhs[1]);
  mwSize dest_ndim = ndim - nremove;
  unsigned *dest_size;
  unsigned *mask;
  mwSize total = mxGetNumberOfElements(prhs[0]), dest_total;

  if(dest_ndim == 0) {
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    dest = mxGetPr(plhs[0]);
    for(i=0;i<total;i++) {
      *dest += src[i];
    }
    addflops(total-1);
    return;
  }
  mask = mxCalloc(ndim,sizeof(mwSize));
  for(i=0;i<nremove;i++) mask[(mwSize)dim[i]-1] = 1;
  dest_size = mxCalloc(dest_ndim,sizeof(mwSize));
  j = 0;
  for(i=0;i<ndim;i++) {
    if(!mask[i]) dest_size[j++] = size[i];
  }
  plhs[0] = mxCreateNumericArray(dest_ndim, dest_size, mxDOUBLE_CLASS, mxREAL);
  dest_total = mxGetNumberOfElements(plhs[0]);
  dest = mxGetPr(plhs[0]);
  addflops(total-dest_total);
  if(dest_total == total) {
    memcpy(dest,src,total*sizeof(double));
    return;
  }
  ndsum(dest,src,ndim,size,total,mask);
}


