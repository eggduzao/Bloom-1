#include "mex.h"

#if 1

#ifdef UNDERSCORE_LAPACK_CALL
int dtrsm_(char *side, char *uplo, char *transa, char *diag, 
					 int *m, int *n, double *alpha, double *a, int *lda, 
					 double *b, int *ldb)
#else
	int dtrsm(char *side, char *uplo, char *transa, char *diag, 
						int *m, int *n, double *alpha, double *a, int *lda, 
						double *b, int *ldb)
#endif
{
  int i,j,k;
#define A(I,J) a[(I) + (J)*(*lda)]
#define B(I,J) b[(I) + (J)*(*ldb)]

  if(*uplo == 'U') {
    for(j=0;j<*n;j++) {
      for(k=*m-1;k>=0;k--) {
				if(B(k,j) != 0.) {
					B(k,j) /= A(k,k);
					for(i=0;i<k;i++) {
						B(i,j) -= B(k,j) * A(i,k);
					}
				}
      }
    }
  } else {
    for(j=0;j<*n;j++) {
      for(k=0;k<*m;k++) {
				if(B(k,j) != 0.) {
					B(k,j) /= A(k,k);
					for(i=k+1;i<*m;i++) {
						B(i,j) -= B(k,j) * A(i,k);
					}
				}
      }
    }
  }
	return 0;
}

#else

typedef int logical;

logical lsame_(char *ca, char *cb)
{
  return(*ca == *cb);
}

int xerbla_(char *srname, int *info)
{
  mexErrMsgTxt(srname);
}

#ifdef UNDERSCORE_LAPACK_CALL
int dtrsm_(char *side, char *uplo, char *transa, char *diag, 
	  int *m, int *n, double *alpha, double *a, int *lda, 
	  double *b, int *ldb)
#else
int dtrsm(char *side, char *uplo, char *transa, char *diag, 
	  int *m, int *n, double *alpha, double *a, int *lda, 
	  double *b, int *ldb)
#endif
{

    int a_dim1, a_offset, b_dim1, b_offset;

    static int info;
    static double temp;
    static int i, j, k;
    static logical lside;
    static int nrowa;
    static logical upper;
    static logical nounit;

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]

    lside = lsame_(side, "L");
    if (lside) {
	nrowa = *m;
    } else {
	nrowa = *n;
    }
    nounit = lsame_(diag, "N");
    upper = lsame_(uplo, "U");

    info = 0;
    if (! lside && ! lsame_(side, "R")) {
	info = 1;
    } else if (! upper && ! lsame_(uplo, "L")) {
	info = 2;
    } else if (! lsame_(transa, "N") && ! lsame_(transa, "T") 
	    && ! lsame_(transa, "C")) {
	info = 3;
    } else if (! lsame_(diag, "U") && ! lsame_(diag, "N")) {
	info = 4;
    } else if (*m < 0) {
	info = 5;
    } else if (*n < 0) {
	info = 6;
    } else if (*lda < nrowa) {
	info = 9;
    } else if (*ldb < *m) {
	info = 11;
    }
    if (info != 0) {
	xerbla_("DTRSM ", &info);
	return 0;
    }

    if (*n == 0) {
	return 0;
    }

    if (*alpha == 0.) {
	for (j = 1; j <= *n; ++j) {
	    for (i = 1; i <= *m; ++i) {
		B(i,j) = 0.;
	    }
	}
	return 0;
    }

    if (lside) {
	if (lsame_(transa, "N")) {

	    if (upper) {
		for (j = 1; j <= *n; ++j) {
		    if (*alpha != 1.) {
			for (i = 1; i <= *m; ++i) {
			    B(i,j) = *alpha * B(i,j);
			}
		    }
		    for (k = *m; k >= 1; --k) {
		      if (B(k,j) != 0.) {
			printf("%d %d %g %g\n",k,j,B(k,j),A(k,k));
			B(k,j) /= A(k,k);
			for (i = 1; i <= k-1; ++i) {
			  B(i,j) -= B(k,j) * A(i,k);
			}
		      }
		    }
		}
	    } else {
		for (j = 1; j <= *n; ++j) {
		    if (*alpha != 1.) {
			for (i = 1; i <= *m; ++i) {
			    B(i,j) = *alpha * B(i,j);
			}
		    }
		    for (k = 1; k <= *m; ++k) {
			if (B(k,j) != 0.) {
			    if (nounit) {
				B(k,j) /= A(k,k);
			    }
			    for (i = k + 1; i <= *m; ++i) {
				B(i,j) -= B(k,j) * A(i,k);
			    }
			}
		    }
		}
	    }
	} else {

	    if (upper) {
		for (j = 1; j <= *n; ++j) {
		    for (i = 1; i <= *m; ++i) {
			temp = *alpha * B(i,j);
			for (k = 1; k <= i-1; ++k) {
			    temp -= A(k,i) * B(k,j);
			}
			if (nounit) {
			    temp /= A(i,i);
			}
			B(i,j) = temp;
		    }
		}
	    } else {
		for (j = 1; j <= *n; ++j) {
		    for (i = *m; i >= 1; --i) {
			temp = *alpha * B(i,j);
			for (k = i + 1; k <= *m; ++k) {
			    temp -= A(k,i) * B(k,j);
			}
			if (nounit) {
			    temp /= A(i,i);
			}
			B(i,j) = temp;
		    }
		}
	    }
	}
    } else {
	if (lsame_(transa, "N")) {

	    if (upper) {
		for (j = 1; j <= *n; ++j) {
		    if (*alpha != 1.) {
			for (i = 1; i <= *m; ++i) {
			    B(i,j) = *alpha * B(i,j);
			}
		    }
		    for (k = 1; k <= j-1; ++k) {
			if (A(k,j) != 0.) {
			    for (i = 1; i <= *m; ++i) {
				B(i,j) -= A(k,j) * B(i,k);
			    }
			}
		    }
		    if (nounit) {
			temp = 1. / A(j,j);
			for (i = 1; i <= *m; ++i) {
			    B(i,j) = temp * B(i,j);
			}
		    }
		}
	    } else {
		for (j = *n; j >= 1; --j) {
		    if (*alpha != 1.) {
			for (i = 1; i <= *m; ++i) {
			    B(i,j) = *alpha * B(i,j);
			}
		    }
		    for (k = j + 1; k <= *n; ++k) {
			if (A(k,j) != 0.) {
			    for (i = 1; i <= *m; ++i) {
				B(i,j) -= A(k,j) * B(i,k);
			    }
			}
		    }
		    if (nounit) {
			temp = 1. / A(j,j);
			for (i = 1; i <= *m; ++i) {
			    B(i,j) = temp * B(i,j);
			}
		    }
		}
	    }
	} else {

	    if (upper) {
		for (k = *n; k >= 1; --k) {
		    if (nounit) {
			temp = 1. / A(k,k);
			for (i = 1; i <= *m; ++i) {
			    B(i,k) = temp * B(i,k);
			}
		    }
		    for (j = 1; j <= k-1; ++j) {
			if (A(j,k) != 0.) {
			    temp = A(j,k);
			    for (i = 1; i <= *m; ++i) {
				B(i,j) -= temp * B(i,k);
			    }
			}
		    }
		    if (*alpha != 1.) {
			for (i = 1; i <= *m; ++i) {
			    B(i,k) = *alpha * B(i,k);
			}
		    }
		}
	    } else {
		for (k = 1; k <= *n; ++k) {
		    if (nounit) {
			temp = 1. / A(k,k);
			for (i = 1; i <= *m; ++i) {
			    B(i,k) = temp * B(i,k);
			}
		    }
		    for (j = k + 1; j <= *n; ++j) {
			if (A(j,k) != 0.) {
			    temp = A(j,k);
			    for (i = 1; i <= *m; ++i) {
				B(i,j) -= temp * B(i,k);
			    }
			}
		    }
		    if (*alpha != 1.) {
			for (i = 1; i <= *m; ++i) {
			    B(i,k) = *alpha * B(i,k);
			}
		    }
		}
	    }
	}
    }

    return 0;

}
#endif


