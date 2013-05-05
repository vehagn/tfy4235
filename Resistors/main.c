#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

double pres = 2.0e-4;
double eps;
double tmp;
int iterations;
int N = 1024;
int i;

int main (int argc, char** argv){

    if (argc > 1){
        N = atoi(argv[1]);
    }
    double* cond = (double*)malloc(sizeof(double)*N);
    double* v_old = (double*)malloc(sizeof(double)*(N+1));
    double* v_new = (double*)malloc(sizeof(double)*(N+1));

    double* jacobi = (double*)malloc(sizeof(double)*(N+1));
    double* gauss = (double*)malloc(sizeof(double)*(N+1));
    double* sor = (double*)malloc(sizeof(double)*(N+1));

    srand((unsigned)time(NULL));
    for (i = 0; i < N; i++)
        cond[i] = ((double)rand()/(double)RAND_MAX);

    printf("Jacobi:\n");
    for (i = 0; i<N+1;i++)
        v_old[i] = i/(double)N;

    iterations = 0; eps = 1;
    while (eps > pres){
        for (i=1;i<N;i++){
            v_new[i] = (cond[i-1]*v_old[i-1] + cond[i]*v_old[i+1])/(cond[i-1] + cond[i]);
        }
        for(i=1;i<N;i++){
            v_old[i] = v_new[i];
        }
        eps = 0;
        for (i=1;i<N;i++){
            tmp = (cond[i-1]*v_old[i-1] + cond[i]*v_old[i+1] - (cond[i-1] + cond[i])*v_old[i]);
            eps += tmp*tmp;
        }
        eps = sqrt(eps);
        iterations++;
    }
    printf("%i iterations.\n",iterations);
    for (i=0;i<=N;i++)
        jacobi[i] = v_new[i];



    printf("Gauss-Seidel:\n");
    for (i = 0; i<N+1;i++)
        v_new[i] = i/(double)N;

    iterations = 0; eps = 1;
    while (eps > pres){
        for (i=1;i<N;i++){
            v_new[i] = (cond[i-1]*v_new[i-1] + cond[i]*v_new[i+1])/(cond[i-1] + cond[i]);
        }
        eps = 0;
        for (i=1;i<N;i++){
            tmp = (cond[i-1]*v_new[i-1] + cond[i]*v_new[i+1] - (cond[i-1] + cond[i])*v_new[i]);
            eps += tmp*tmp;
        }
        eps = sqrt(eps);
        iterations++;
    }
    printf("%i iterations.\n",iterations);
    for (i=0;i<=N;i++)
        gauss[i] = v_new[i];


    printf("SOR:\n");
    for (i = 0; i<N+1;i++)
        v_old[i] = i/(double)N;

    iterations = 0; eps = 1;
    double omega = 1.5;

    while (eps > pres){
        for (i=0;i<N-1;i++){
            v_new[i] = v_old[i+1];
        }
        for (i=1;i<N;i++){
            v_old[i] = (cond[i-1]*v_old[i-1] + cond[i]*v_old[i+1])/(cond[i-1] + cond[i]);
        }
        for (i=1;i<N;i++){
            v_old[i] = omega*v_old[i] + (1.0-omega)*v_new[i-1];
        }
        eps = 0;
        for (i=1;i<N;i++){
            tmp = (cond[i-1]*v_old[i-1] + cond[i]*v_old[i+1] - (cond[i-1] + cond[i])*v_old[i]);
            eps += tmp*tmp;
        }
        eps = sqrt(eps);
        iterations++;
    }
    printf("%i iterations.\n",iterations);
    for (i=0;i<=N;i++)
        sor[i] = v_new[i];

    printf("\nv[%d]: jac=%.9f gs=%.9f sor=%.9f\n",N/2,jacobi[N/2],gauss[N/2],sor[N/2]);

    free(cond);
    free(v_old);
    free(v_new);

    free(jacobi);
    free(gauss);
    free(sor);

    return 0;
}
