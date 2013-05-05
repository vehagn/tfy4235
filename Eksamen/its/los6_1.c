#include <stdio.h>

int IBM = 12345;    // "random" seed
double rinv = 0.5 / ( (long int)1 << 31 );  // scaling, 2^32

double randIBMf() {
    IBM *= 16807;
    return ( IBM*rinv + 0.5 );  // double on [0,1)
}

int main() {
    double a, b = randIBMf();       // b = r_0
    int n = 0, N = 10000;

    double xmax = 1.0, ymax = 1.0;  // change these for different ranges

    FILE * outfile = fopen("out_first.d","w");

    fprintf(outfile,"# output file for problem 1 of problem set 6\n");
    fprintf(outfile,"# %d points restricted to the [0:%.3f], [0:%.3f] quadrant\n",N, xmax, ymax);
    fprintf(outfile,"# format: x-coordinate | y-coordinate\n\n");

    while (n < N) {
        a = b;
        b = randIBMf();

        if (a < xmax && b < ymax) {
            fprintf(outfile,"%f %f\n",a,b);
            n++;
        }
    }

    fclose(outfile);
}
