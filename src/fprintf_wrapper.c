#include "stdio.h"

int fprintf_wrapper(const char* fn,
  int *i1,int *i3, int *setis1, int *cutv1,
  float *distv1, int *ivec21, int *invvec1, int *cutv12,
  int *vv11, int *vv21, int *vv31, int *vv12, int *vv22,
  int *vv32, int *vv13, int *vv23, int *vv33,
  int *ii1, int *jj1, int *kk1,
  int *newlornot) {
  FILE *fp;
  int retv;
//4(i10,1x),g12.5,1x,1000(i10,1x)
  fp = fopen(fn,"a");
  if (*newlornot > 0) {
    retv = fprintf(fp,"%d %d %d %d %g %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",
    *i1,*i3,*setis1,*cutv1,*distv1,*ivec21,*invvec1,*cutv12,*vv11,*vv21,*vv31,
    *vv12,*vv22,*vv32,*vv13,*vv23,*vv33,*ii1,*jj1,*kk1);
  } else {
    retv = fprintf(fp,"%d %d %d %d %g %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d",
    *i1,*i3,*setis1,*cutv1,*distv1,*ivec21,*invvec1,*cutv12,*vv11,*vv21,*vv31,
    *vv12,*vv22,*vv32,*vv13,*vv23,*vv33,*ii1,*jj1,*kk1);
  }
  fclose(fp);
  return retv;
//
}
