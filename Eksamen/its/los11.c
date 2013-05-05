#include <stdio.h>
#include <time.h>


int IBM = 123451;

double rinv = 0.5/( ((long int)1 << 31) - 1.0);

double ranf() {
  IBM *= 16807;
  return (IBM*rinv + 0.5);
}

int calcEnergy(int * spins, int ** couplings, int N){
  int E = 0, i,j;

  for (j = 0; j < N; j++) {
    for (i = 0; i < N; i++) {
      int k = j*N + i;
      int up = N*((j+1)%N) + i;
      int right = j*N + (i+1)%N;
      E -= spins[k]*couplings[k][0]*spins[right] + spins[k]*couplings[k][1]*spins[up];
    }
  }

  return E;
}

int main() {
  int i,j;

  const int N = 5;
  double tstart = 10;
  double ct = 0.999;

  IBM = 2*time(0) + 1;

  for (i = 0; i < 1000; i++) {
    ranf();
  }

  int * spins = (int*)malloc(sizeof(int)*N*N);
  int * up = (int*)malloc(sizeof(int)*N*N);
  int * down = (int*)malloc(sizeof(int)*N*N);
  int * right = (int*)malloc(sizeof(int)*N*N);
  int * left = (int*)malloc(sizeof(int)*N*N);
  int ** couplings = (int**)malloc(sizeof(int*)*N*N);

  for (j = 0; j < N; j++) {
    for (i = 0; i < N; i++) {
      int k = j*N + i;
      spins[k] = (ranf() < 0.5 ? -1 : 1);
      couplings[k] = (int*)malloc(sizeof(int)*2);
      couplings[k][0] = (ranf() < 0.5 ? -1 : 1);  // right coupling
      couplings[k][1] = (ranf() < 0.5 ? -1 : 1);  // up coupling

      //Initializing neighbor arrays for simplified energy calculations.
      up[k] = k + N;
      down[k] = k - N;
      left[k] = k - 1;
      right[k] = k + 1;

      if (j == 0) {
        down[k] += N*N;
      } else if (j == (N - 1)) {
        up[k] -= N*N;
      }
      if (i == 0) {
        left[k] += N;
      } else if (i == (N - 1)) {
        right[k] -= N;
      }
    }
  }

  int E = calcEnergy(spins,couplings,N);

  printf("E=%d\n",E);

  double temp = tstart;//tstart;

  for (j = 0; j < 10000; j++){  // number of generations/sweeps
    for (i = 0; i < N*N; i++){  // full sweep
      int x = (int)(ranf()*N);
      int y = (int)(ranf()*N);

      int k = y*N + x;

      int cLeft = couplings[left[k]][0];
      int cRight = couplings[k][0];
      int cDown = couplings[down[k]][1];
      int cUp = couplings[k][1];

      int sLeft = spins[left[k]];
      int sRight = spins[right[k]];
      int sDown = spins[down[k]];
      int sUp = spins[up[k]];

      int dE = 2*spins[k]*(cRight*sRight + cUp*sUp + cLeft*sLeft + cDown*sDown);

      double r = ranf();
      if (r < exp(-dE/temp)){
        E += dE;
        spins[k] *= -1;
      }
    }
    temp *= ct;
    printf("E=%d, T = %f\n",E,temp);
  }

  return 0;
}
