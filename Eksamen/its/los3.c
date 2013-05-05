#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int IBM = 4711;
int N = 1024;
double rinv = 0.5/(1 << 31);
double omega = 1.5;
double pres = 2.0e-5;

double * cond, * resi, * volt, * vnew, *jac, *gs, *sor;

double randIBMf(){  	  	//Returns random numbers between 0 and 1
    IBM*=16807;
    return IBM*rinv + 0.5;
}

int main(){
    int i;


		/* INITIALIZATION */
    for(i=0;i<1000000;i++)
        randIBMf();

    cond = (double*)malloc(sizeof(double)*N);
    volt = (double*)malloc(sizeof(double)*(N+1));
    vnew = (double*)malloc(sizeof(double)*(N+1));

    jac = (double*)malloc(sizeof(double)*(N+1));
    gs = (double*)malloc(sizeof(double)*(N+1));
    sor = (double*)malloc(sizeof(double)*(N+1));

    for (i = 0; i < N; i++)
        cond[i] = randIBMf();


		/* JACOBI */
    for(i=0;i<N+1;i++)
        volt[i] = (double)i/N;

    printf("Jacobi: ");
    int n = 0;
    double eps = 2*pres;
    while (eps > pres){
        for(i=1;i<N;i++){ 
            vnew[i] = (cond[i-1]*volt[i-1] + cond[i]*volt[i+1])/(cond[i-1] + cond[i]);
        }
        eps = 0;
        double tmp;
        for(i=1;i<N;i++){
            volt[i] = vnew[i];
        }
        for(i=1;i<N;i++){
            tmp = (cond[i-1]*volt[i-1] + cond[i]*volt[i+1] - (cond[i-1] + cond[i])*volt[i]);
            eps += tmp*tmp;
        }
        eps = sqrt(eps);
        n++;
    }
    printf("%d iterations.\n",n);
 
   	for (i=0;i<N+1;i++){
      jac[i] = volt[i];
    }
 

		/* GAUSS-SEIDEL */
    for(i=0;i<N+1;i++)
        volt[i] = (double)i/N;

    printf("Gauss-Seidel: ");
    n = 0;
    eps = 2*pres;
    while (eps >= pres){
        for(i=1;i<N;i++){ 
            volt[i] = (cond[i-1]*volt[i-1] + cond[i]*volt[i+1])/(cond[i-1] + cond[i]);
        }
        eps = 0;
        double tmp;
        for(i=1;i<N;i++){
            tmp = (cond[i-1]*volt[i-1] + cond[i]*volt[i+1] - (cond[i-1] + cond[i])*volt[i]);
            eps += tmp*tmp;
        }
        eps = sqrt(eps);
        n++;
    }
    printf("%d iterations.\n",n);

    for (i=0;i<N+1;i++){
      gs[i] = volt[i];
    }


		/* Successive Over-Relaxation */
    for(i=0;i<N+1;i++)
        volt[i] = (double)i/N;

    printf("SOR, omega = %.3f: ",omega);

    n = 0;
    eps = 2*pres;
    while (eps >= pres){
        for(i=0;i<N-1;i++){
            vnew[i] = volt[i+1];
        }
        for(i=1;i<N;i++){ // GS step
            volt[i] = (cond[i-1]*volt[i-1] + cond[i]*volt[i+1])/(cond[i-1] + cond[i]);
        }
        for(i=1;i<N;i++){ 
            volt[i] = omega*volt[i] + (1.0-omega)*vnew[i-1];
        }
        eps = 0;
        double tmp;
        for(i=1;i<N;i++){
            tmp = (cond[i-1]*volt[i-1] + cond[i]*volt[i+1] - (cond[i-1] + cond[i])*volt[i]);
            eps += tmp*tmp;
        }
        eps = sqrt(eps);
        n++;
    }
    printf("%d iterations.\n",n);
    for (i=0;i<N+1;i++){
      sor[i] = volt[i];
    }

    printf("\n");
  
		/* Final output, to compare the convergence */
    printf("v[%d]: jac=%.9f gs=%.9f sor=%.9f\n",N/2,jac[N/2],gs[N/2],sor[N/2]);

    free(cond);
    free(volt);
    free(vnew);

		free(jac);
		free(gs);
		free(sor);

    return 0;
}
