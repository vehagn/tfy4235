#include<stdio.h>
#include<stdlib.h>
//#include<math.h>
#include<time.h>

double pres = 2.0e-4;
int N = 1024;
int i;
int iterations;

int main (int argc, char** argv){

    if (argc > 1){
        N = atoi(argv[1]);
    }

    double* c = (double*)malloc(sizeof(double)*N);
    double* volt = (double*)malloc(sizeof(double)*(N+1));
    double* r = (double*)malloc(sizeof(double)*(N-1));
    double* p = (double*)malloc(sizeof(double)*(N));
    double* Ap = (double*)malloc(sizeof(double)*(N));
    double alpha;

    srand((unsigned)time(NULL));
    for (i = 0; i < N; i++)
        c[i] = ((double)rand()/(double)RAND_MAX);

    for (i = 0; i<N+1;i++)
        volt[i] = i/(double)N;

    double rTr=0, pTAp;
    for (i=1;i<N;i++){
        r[i-1] = -c[i-1]*(volt[i-1]-volt[i]) - c[i]*(volt[i+1] - volt[i]);
        p[i] = r[i-1];
        rTr += r[i-1]*r[i-1];
    }

    printf("rTr: %f\n",rTr);
    p[0] = 0; r[N] = 0;

    double beta_n, beta_d, beta;
    while (rTr > pres*pres){
        pTAp = 0;
        beta_n = 0; beta_d = 0;
        for (i=1;i<N;i++){
           Ap[i-1] = c[i-1]*(p[i-1]-p[i]) + c[i]*(p[i+1]-p[i]);
           pTAp += p[i]*Ap[i-1];
        }

        alpha = rTr/pTAp;

        for (i=1;i<N;i++){
           volt[i] += alpha*p[i];
           r[i-1] -= alpha*Ap[i-1];
           beta_n += r[i-1]*r[i-1];
        }

        beta_d = rTr;
        beta = beta_n/beta_d;

        for (i=1;i<N;i++){
           p[i] = r[i-1] + beta*p[i];
        }
        rTr = beta_n;
        iterations++;
    }
    printf("Iterations: %i \n",iterations);
    return 0;
}
