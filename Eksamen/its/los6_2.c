#include <stdio.h>
#include <math.h>

int IBM = 12345;    // "random" seed
const double rinv = 0.5 / ( (long int)1 << 31 );  // scaling, 2^32
const double PI = 4.0*atan(1.0);
const double PI2 = 8.0*atan(1.0);

double randIBMf() {
    IBM *= 16807;
    return ( IBM*rinv + 0.5 );  // double on [0,1)
}

void randIBMn(double * y1, double * y2) {
    double x1 = randIBMf();
    double x2 = randIBMf();

    double r = sqrt(-2.0*log(x1));
    double t = PI2*x2;
    (*y1) = r*cos(t);
    (*y2) = r*sin(t);
}

double max(double x, double y) {
    return (x > y ? x : y);
}
double max3(double a, double b, double c) {
    return max(a,max(b,c));
}

int main() {
    int i,j,N;
    int M = 10000;

    FILE * outfile = fopen("out_second.d","w");
    
    fprintf(outfile,"# output file for problem 2 of problem 6\n");
    fprintf(outfile,"# average maximum value of N random gaussian numbers\n");
    fprintf(outfile,"# averaged over M = %d batches\n",M);
    fprintf(outfile,"# format: N | sqrt(log(N)) (guess) | xavg\n\n");

    for (N = 10; N <= 100000; N*= 10) {
        printf("N = %d",N);

        double xavg = 0;
        double xlar = 0;

        for (j = 0; j < M; j++) {
            xlar = 0;
            double y1,y2,y;
            for (i = 0; i < N/2; i++) {
                randIBMn(&y1,&y2);
                xlar = max3(xlar,y1,y2);
            }
    
            xavg += xlar;
        }

        xavg = xavg/M;

        fprintf(outfile, "%d %f %f\n",N,sqrt(log(N)),xavg);

        printf(" ... done\n");
    }

    fclose(outfile);
}
