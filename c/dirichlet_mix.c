#define _USE_MATH_DEFINES 1
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include "util.h"

#ifdef _MSC_VER
#define finite _finite
#define isnan _isnan
#endif

#pragma data_seg(".seed")
static long ix = 101;
static long iy = 1001;
static long iz = 10001;
static double RandN_previous = 0;
static int RandN_usePrevious = 0;
#pragma data_seg()
#pragma comment(linker,"/section:.seed,rws")

#if 1

double Rand(void)
{
  static float u;
  
  ix = 171*(ix % 177)-2*(ix/177);
  iy = 172*(iy % 176)-2*(iy/176);
  iz = 170*(iz % 178)-2*(iz/178);
  
  if (ix<0) ix = ix + 30269;
  if (iy<0) iy = iy + 30307;
  if (iz<0) iz = iz + 30323;
  
  u = ((float) ix)/30269 +
                ((float) iy)/30307 + ((float) iz)/30323;
  u -= (float)(int)u;
  return(u);
}
#else

double Rand(void)
{
  mxArray *plhs[1];
  if(mexCallMATLAB(1,plhs,0,NULL,"rand")) {
    mexErrMsgTxt("mexCallMATLAB(rand) failed");
  }
  return mxGetPr(plhs[0])[0];
}
#endif

void ResetSeed(void)
{
  SetSeed(101,1001,10001);
}

void SetSeed(long new_ix, long new_iy, long new_iz)
{
  ix = new_ix;
  iy = new_iy;
  iz = new_iz;
  RandN_usePrevious = 0;
}

void GetSeed(long *ix_out, long *iy_out, long *iz_out)
{
  *ix_out = ix;
  *iy_out = iy;
  *iz_out = iz;
  RandN_usePrevious = 0;
}

double RandN(void)
{
  double x,y,radius;
  if(RandN_usePrevious) {
    RandN_usePrevious = 0;
    return RandN_previous;
  }
  do {
    x = 2*Rand()-1;
    y = 2*Rand()-1;
    radius = (x*x)+(y*y);
  } while((radius >= 1.0) || (radius == 0.0));

  radius = sqrt(-2*log(radius)/radius);
  x *= radius;
  y *= radius;
  RandN_previous = y;
  RandN_usePrevious = 1;
  return x;
}

double GammaRand(double a)
{

  double boost, d, c, v;
  if(a < 1) {
    boost = exp(log(Rand())/a);
    a++;
  } 
  else boost = 1;
  d = a-1.0/3; c = 1.0/sqrt(9*d);
  while(1) {
    double x,u;
    do {
      x = RandN();
      v = 1+c*x;
    } while(v <= 0);
    v = v*v*v;
    x = x*x;
    u = Rand();
    if((u < 1-.0331*x*x) || 
       (log(u) < 0.5*x + d*(1-v+log(v)))) break;
  }
  return( boost*d*v );
}

double BetaRand(double a, double b)
{
  double g = GammaRand(a);
  return g/(g + GammaRand(b));
}

int BinoRand(double p, int n)
{
  int r = 0;
  if(isnan(p)) return 0;
  if(p < DBL_EPSILON) return 0;
  if(p >= 1-DBL_EPSILON) return n;
  if((p > 0.5) && (n < 15)) {
    int i;
    for(i=0;i<n;i++) {
      if(Rand() < p) r++;
    }
    return r;
  }
  if(n*p < 10) {
    double q = -log(1-p), e = -log(Rand()), s;
    r = n;
    for(s = e/r; s <= q; s += e/r) {
      r--;
      if(r == 0) break;
      e = -log(Rand());
    }
    r = n-r;
    return r;
  }
  if (1) {
    int i = (int)(p*(n+1));
    double b = BetaRand(i, n+1-i);
    if(b <= p) r = i + BinoRand((p-b)/(1-b), n-i);
    else r = i - 1 - BinoRand((b-p)/b, i-1);
    return r;
  }
}


