#include <stdio.h>
#include <math.h>
#include <time.h>

double PI5 = 20.0*atan(1.0);
int IBM = 1311;
double rinv = 0.5/2147483647.0;
const double R1 = 1.0/4294967296.0;  // 1.0/2^32
const double R2 = 1.0/4294967295.0;  // 1.0/(2.0^32-1)

int randIBM(){
    IBM *= 16807;
    return IBM;
}

double ranf(){
    IBM *= 16807;
    return rinv*IBM + 0.5;
}

double intToDouble(int x, double a, double b){
    double c = 0.5*(a+b);
    double d = (b-a)*R1;    //can be precomputed given a and b is known
    return c + d*x;
}
double intToDouble2(int x, double a, double b){
    double c = 0.5*(a+b) + 0.5*(b-a)*R2;
    double d = (b-a)*R2;    //can be precomputed given a and b is known
    return c + d*x;
}

void printBinary(int x){
    int i;
    for (i = 31; i >= 0; i--){
        printf("%d",!!(x&(1<<i)));
        if (!(i & 7))
            printf(" ");
    }
    printf("\n");
}

/*
Explanation swapbits:
 - the bits of z is 1 for the bits we want to swap
 - the bits of x^y is 1 where the bits of x and y differ
 - the bits of w=(x^y)&z is therefore 1 for the bits of x
   and y that we want to change
 - if a <= b: z = z(a,b) = 2^(b+1) - 2^(a), else: z = z(b,a)

*/
void swapbits(int * x1,int * y1,int a,int b){ //
    int x = *x1;
    int y = *y1;
    int z = (a <= b) ? (1 << b+1) - (1 << a) : (1 << a+1) - (1 << b);
    int w = (x^y)&z;
    *x1 = x ^ w;
    *y1 = y ^ w;
}

double f(double x){
    return cos(x)*exp(-1*x*x);
}

int main(){
    /* initialize */
    int i,j,k;
    IBM = 2*time(0)+1;

    for (i = 0; i < 1000; i++)
        randIBM();

    int N = 100, Ngen = 100, imax, ipop[N], ipopn[N];
    double pvlg = 0.75, qvlg = 0.0, rmut = 0.02, xmax, fmax, xpop[N], fpop[N];

    FILE * outfile = fopen("out.d","w");
    FILE * outfile2 = fopen("out2.d","w");

    for (i = 0; i < N; i++){
        ipop[i] = randIBM();
    }

    for (i = 0; i < Ngen; i++){
        /* find maximum and print to file */
        fmax = -10000.0;
        xmax = 0.0;
        for (j = 0; j < N; j++){
            xpop[j] = intToDouble(ipop[j],-PI5,PI5);
            fpop[j] = f(xpop[j]);
            fprintf(outfile2,"%d %g %g\n",i,xpop[j],fpop[j]);

            if (fpop[j] > fmax){
                imax = j;
                fmax = fpop[j];
                xmax = xpop[j];
            }
        }

        printf("%d: most fit individual: i=%d -> %g, fitness=%g\n",i,imax,xmax,fmax);
        fprintf(outfile,"%d %f %f\n",i,xmax,fmax);

        /* selection */
        for (j = 0; j < N; j++){
            int jpop = N*ranf();       // random pair pair, first
            int kpop = N*ranf();       // random pair, second

            int imin = jpop;
            int imax = kpop;

            double fjpop = fpop[jpop]; // fitness first
            double fkpop = fpop[kpop]; // fitness second

            if (fjpop > fkpop){
                imin = kpop;
                imax = jpop;
            }

            // most fit individual chosen with probability pvlg
            // least fit chosen with probability 1 - pvlg
            int ivlg = (int)(pvlg + ranf()); // [pvlg, pvlg+1.0]
            ipopn[j] = (ivlg ? ipop[imax] : ipop[imin]);
        }
        for (j = 0; j < N; j++){
            ipop[j] = ipopn[j];
        }

        /* two-point crossover */
        for (j = 0; j < N; j++){
            int ivlg = (int)(qvlg + ranf());
            if (ivlg) {                 // probability qvlg to be true
                int kpop = N*ranf();    // choose random individual
                int jstr = 32*ranf();   // choose first random bit
                int kstr = 32*ranf();   // choose second random bit
                swapbits(&ipop[j],&ipop[kpop],jstr,kstr);
            }
        }

        /* mutation */
        for (j = 0; j < N; j++){
            int mut = 0;
            for (k = 0; k < 32; k++){
                int imut = (int)(rmut + ranf());
                mut = mut ^ (imut << k);
            }
            ipop[j] = ipop[j] ^ mut; // flip mutated bits
        }
    }

    fclose(outfile);
    fclose(outfile2);

    return 0;
}
