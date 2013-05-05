#include<stdio.h>
#include<stdlib.h>
//#include<math.h>
#include<time.h>

int IBM = 1927;  //SEED
const int N = 1024;
double rinv = 0.5/((long int)1 << 31);
double pres = 2.0e-4;
double rps, rpn, am, bm;

double * cond, * resi, * volt, * vnew, * p, * r, * ap, *gs, *sor, *cg;


double randIBM(){          //returns random numbrs on the unit interval
    //IBM*=16807;
    //return IBM*rinv + 0.5;
    return ((double)rand()/(double)RAND_MAX);
}

int main(){
    srand((unsigned)time(NULL));
    int i;
		int n = 0;

    cond = (double*)malloc(sizeof(double)*N);
    volt = (double*)malloc(sizeof(double)*(N+1));
    p = (double*)malloc(sizeof(double)*(N+1));
    r = (double*)malloc(sizeof(double)*(N-1));
    ap = (double*)malloc(sizeof(double)*(N-1));

		for(i=0;i<1000;i++) //throwaway random numbers
        randIBM();

    for (i = 0; i < N; i++)
        cond[i] = randIBM();

    for(i=0;i<N+1;i++)
        volt[i] = (double)i/N;

		//initialization of work vectors
    rps = 0;
    for (i=1;i<N;i++){
        p[i] = -cond[i-1]*(volt[i-1] - volt[i]) - cond[i]*(volt[i+1]-volt[i]);
        r[i-1] = p[i];
        rps += r[i-1]*r[i-1];
    }
    p[0] = 0;
    p[N] = 0;

    while ((rps > pres*pres)&&(5)){
        //printf("rps: %f \n",rps);
        am = 0;
				rpn = 0;

        for (i=1;i<N;i++){
            ap[i-1] = cond[i-1]*(p[i-1] - p[i]) + cond[i]*(p[i+1] - p[i]);
            am += p[i]*ap[i-1];
        }

        am = rps/am;

        for (i=1;i<N;i++){
            volt[i] += am*p[i];
        }

        for (i=1;i<N;i++){
            r[i-1] -= am*ap[i-1];
        }

        for (i=0;i<N-1;i++){
            rpn += r[i]*r[i];
        }

        bm = rpn/rps;
        for (i=1;i<N;i++){
            p[i] = r[i-1] + bm*p[i];
        }

				rps = rpn;
        n++;
    }

    printf("Conjugate gradient: %d iterations\n",n);

    free(cond);
    free(volt);
    free(p);
    free(r);
    free(ap);

    return 0;
}
