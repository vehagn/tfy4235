#include <stdio.h>
#include <stdlib.h>

#define max(a,b)        (a < b ? b : a)
#define abs(a)          (a < 0 ? -a : a)
#define isPowerOfTwo(n) (n && !(n & (n - 1)))
#define C0  0.4829629131445341
#define C1  0.8365163037378079
#define C2  0.2241438680420134
#define C3 -0.1294095225512604

int IBM = 1234511;
double rinv = 0.5/( ((long int)1 << 31) -  1);

double ranf(){
  IBM *= 16807;
  return (IBM*rinv + 0.5);
}


void daub4(double * a, const int n, int isign){
  if (n >= 4) {
    if (isign == 1) {
      int i,j,nh,nh1;
      double * wksp = (double * )malloc(sizeof(double)*n);

      nh = n / 2;
      nh1 = nh + 1;
      for (j = 0, i = 0; j < n-3; j += 2, i++) {
        wksp[i]    = C0*a[j] + C1*a[j+1] + C2*a[j+2] + C3*a[j+3];
        wksp[i+nh] = C3*a[j] - C2*a[j+1] + C1*a[j+2] - C0*a[j+3];
      }
      wksp[i]    = C0*a[n-2] + C1*a[n-1] + C2*a[0] + C3*a[1];
      wksp[i+nh] = C3*a[n-2] - C2*a[n-1] + C1*a[0] - C0*a[1];

      for (i = 0; i < n; i++) {
        a[i] = wksp[i];
      }
      free(wksp);
    } else if (isign == -1) {
      int i,j,nh,nh1;
      double * wksp = (double * )malloc(sizeof(double)*n);

      nh = n / 2;
      nh1 = nh + 1;
      wksp[0] = C2*a[nh-1]+C1*a[n-1]+C0*a[0]+C3*a[nh1-1];
      wksp[1] = C3*a[nh-1]-C0*a[n-1]+C1*a[0]-C2*a[nh1-1];
      for (i = 0, j = 2; i < nh-1; i++, j += 2) {
        wksp[j]   = C2*a[i] + C1*a[i+nh] + C0*a[i+1] + C3*a[i+nh1];
        wksp[j+1] = C3*a[i] - C0*a[i+nh] + C1*a[i+1] - C2*a[i+nh1];
      }

      for (i = 0; i < n; i++) {
        a[i] = wksp[i];
      }
      free(wksp);
    } else {
      fprintf(stderr,"Invalid transformation type\n Normal: 1\n Inverse: -1\n Current: %d\n",isign);
    }
  }
}

void wt1(double * a, const int n, int isign) {
  if (n >= 4 && isPowerOfTwo(n)) {
    if (isign == 1) {
      double nn = n;
      while (nn >= 4) {
        daub4(a,nn,isign);
        nn /= 2;
      }
    } else if (isign == -1) {
      double nn = 4;
      while (nn <= n) {
        daub4(a,nn,isign);
        nn *= 2;
      }
    } else {
      fprintf(stderr,"Invalid transformation type\n Normal: 1\n Inverse: -1\n Current: %d\n",isign);
    }
  } else {
    fprintf(stderr,"Too few data points, or not power of two\n");
  }
}

int main() {
  int i;
  for (i = 0; i < 1000; i++) {
    ranf();
  }

  int n = (1 << 15);    // n = 2^15
  double filter = 0.001;

  double *newData = (double *)malloc(sizeof(double)*n);
  double *transformedData = (double *)malloc(sizeof(double)*n);
  double *oldData = (double *)malloc(sizeof(double)*n);

  for (i = 1; i < n; i++) {
    double r = 2.0*(int)(ranf() + 0.5) - 1.0;   // -1.0 or 1.0
    newData[i] = newData[i-1] + r;
    oldData[i] = newData[i];
  }


  wt1(newData,n,1);
  double amax = -1.0e+8;
  for (i = 0; i < n; i++) {
    amax = max(amax,abs(newData[i]));
    transformedData[i] = newData[i];
  }
  double adel = amax*filter;
  int numlow = 0;

  for (i = 0; i < n; i++) {
    if (abs(newData[i]) < adel){
      newData[i] = 0.0;
      numlow++;
    }
  }
  printf("%d %d %f\%\n",numlow,n,(double)numlow/n);
  wt1(newData,n,-1);

  FILE * outfile = fopen("out.d","w");
  fprintf(outfile,"# i a[i] b[i] \n# i: index\n# a: compressed random walk\n# b: original random walk\n# \n");
  for (i = 0; i < n; i++) {
    fprintf(outfile,"%d %f %f %f\n",i,newData[i],oldData[i],transformedData[i]);
  }
  fclose(outfile);

  free(oldData);
  free(newData);
  free(transformedData);
  return 0;
}
