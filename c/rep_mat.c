#include "mexutil.h"
#include <string.h>

#define ALWAYS_2D 1

void memrep(char *dest, mwSize chunk, mwSize rep)
{
  if(chunk >= 1024) {
    mwSize i;
    char *p = dest;
    for(i=1;i<rep;i++) {
      p += chunk;
      memcpy(p, dest, chunk);
    }
  } else {
    if(rep == 1) return;
    memcpy(dest + chunk, dest, chunk); 
    if(rep & 1) {
      dest += chunk;
      memcpy(dest + chunk, dest, chunk);
    }
    memrep(dest, chunk<<1, rep>>1);
  }
}

void repmat(char *dest, const char *src, int ndim, mwSize *destdimsize, 
	    mwSize *dimsize, const mwSize *dims, mwSize *rep) 
{
  int d = ndim-1;
  mwSize i;
  mwSize chunk;
  if(d == 0) {
    chunk = dimsize[0];
    memcpy(dest,src,chunk);
  }
  else {
    for(i=0;i<dims[d];i++) {
      repmat(dest + i*destdimsize[d-1], src + i*dimsize[d-1], 
	     ndim-1, destdimsize, dimsize, dims, rep);
    }
    chunk = destdimsize[d-1]*dims[d];
  }
  memrep(dest,chunk,rep[d]);
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  const mxArray *srcmat;
  int ndim, eltsize;
  mwSize *dimsize;
  const mwSize *dims;
  int ndimdest;
  mwSize *destdims, *destdimsize;
  char *src, *dest;
  mwSize *rep;
  int i,nrep;
  int extra_rep = 1;
  int empty;
	double *outp, *inp;
	mwSize m,n,numel;

  if(nrhs < 2) mexErrMsgTxt("Usage: repmat(A, [M N ...])");
  srcmat = prhs[0];
	if(0) {
		m = mxGetM(srcmat);
		n = mxGetN(srcmat);		
		plhs[0] = mxCreateDoubleMatrixE(m, n, mxREAL);
		outp = mxGetPr(plhs[0]);
		inp = mxGetPr(srcmat);
		numel = mxGetNumberOfElements(srcmat);
		memcpy(outp, inp, numel*sizeof(double));
		return;
	}

  if(!mxIsNumeric(srcmat) || mxIsSparse(srcmat) || mxIsCell(srcmat) || mxIsStruct(srcmat)) {
    mexCallMATLAB(nlhs,plhs,nrhs,(mxArray**)prhs,"xrepmat");return;
  }
  ndim = mxGetNumberOfDimensions(srcmat);
  dims = mxGetDimensions(srcmat);
  eltsize = mxGetElementSize(srcmat);

  dimsize = (mwSize*)mxCalloc(ndim, sizeof(mwSize));
  dimsize[0] = eltsize*dims[0];
  for(i=1;i<ndim;i++) dimsize[i] = dimsize[i-1]*dims[i];

  ndimdest = ndim;
  if(nrhs == 2) {
    nrep = mxGetN(prhs[1]);
    if(nrep > ndimdest) ndimdest = nrep;
    rep = (mwSize*)mxCalloc(ndimdest, sizeof(mwSize));
    for(i=0;i<nrep;i++) {
      double repv = mxGetPr(prhs[1])[i];
      rep[i] = (mwSize)repv;
    }
#if ALWAYS_2D
    if(nrep == 1) {
      nrep = 2;
      rep[1] = rep[0];
    }
#endif
  }
  else {
    int ri=0;
    nrep = 0;
    for(i=0;i<nrhs-1;i++) {
      nrep += mxGetNumberOfElements(prhs[i+1]);
    }
    if(nrep > ndimdest) ndimdest = nrep;
    rep = (mwSize*)mxCalloc(ndimdest, sizeof(mwSize));
    for(i=0;i<nrhs-1;i++) {
      double *p = mxGetPr(prhs[i+1]);
      int j, sz = mxGetNumberOfElements(prhs[i+1]);
      for(j=0;j<sz;j++) rep[ri++] = (mwSize)p[j];
    }
  }
  for(i=nrep;i<ndimdest;i++) rep[i] = 1;

  destdims = (mwSize*)mxCalloc(ndimdest, sizeof(mwSize));
  for(i=0;i<ndim;i++) destdims[i] = dims[i]*rep[i];
  for(;i<ndimdest;i++) { 
    destdims[i] = rep[i];
    extra_rep *= rep[i];
  }
  destdimsize = (mwSize*)mxCalloc(ndim, sizeof(mwSize));
  destdimsize[0] = eltsize*destdims[0];
  for(i=1;i<ndim;i++) destdimsize[i] = destdimsize[i-1]*destdims[i];

  plhs[0] = mxCreateNumericArrayE(ndimdest, destdims, mxGetClassID(srcmat), 
				  mxIsComplex(srcmat)?mxCOMPLEX:mxREAL);

  empty = 0;
  for (i=0; i < nrep; i++) {
    if (rep[i]==0) 
      empty = 1;
  }
  if (empty) 
    return;

  src = (char*)mxGetData(srcmat);
  dest = (char*)mxGetData(plhs[0]);
  repmat(dest,src,ndim,destdimsize,dimsize,dims,rep);
  if(ndimdest > ndim) memrep(dest,destdimsize[ndim-1],extra_rep);
  if(mxIsComplex(srcmat)) {
    src = (char*)mxGetPi(srcmat);
    dest = (char*)mxGetPi(plhs[0]);
    repmat(dest,src,ndim,destdimsize,dimsize,dims,rep);
    if(ndimdest > ndim) memrep(dest,destdimsize[ndim-1],extra_rep);
  }
}


