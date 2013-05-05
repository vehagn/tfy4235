#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int IBM = 12345;
double rinv = 0.5/((long int)1 << 31);
const double PI2 = 8.0*atan(1.0);

double minx = -5.0;
double maxx = 5.0;
int bins = 200;
long int N = 1000000;
int corrLen = 1000;
double stepSize = 0.5;

double ranf(){
    IBM *= 16807;
    return rinv*IBM + 0.5;
}

void rang(double * n1, double * n2){
    double x1 = ranf();
    double x2 = ranf();
    double r = sqrt(-2.0*log(x1));
    *n1 = r*cos(PI2*x2);
    *n2 = r*sin(PI2*x2);
}

// SEE METROPOLIS FOR COMMENTS.
void boxMuller(){
    double n1,n2;
    int i,j,bin;
    double dbin = (maxx - minx)/bins;   //bin width

    int * hist;
    double * corr;
    double * xold;

    hist = (int*)malloc(sizeof(int)*bins);
    corr = (double*)malloc(sizeof(double)*corrLen);
    xold = (double*)malloc(sizeof(double)*corrLen);
    // initialize 
    for (i = 0; i < bins; i++){
        hist[i] = 0;
    }
    for (i = 0; i < corrLen; i++){
        xold[i] = 0;
        corr[i] = 0;
    }

    for (i = 0; i  < N/2; i++){
        rang(&n1,&n2);
        
        if (n1 < minx)
            bin = 0;
        else if (n1 >= maxx)
            bin = bins-1;
        else
            bin = (n1 - minx)/dbin;
        
        hist[bin]++;

        for (j = 0; j < corrLen; j++){
            corr[j] += n1*xold[j];
        }

        for (j = 0; j < corrLen-1; j++){
            xold[j] = xold[j+1];
        }

        xold[corrLen-1] = n1;

        if (n2 < minx)
            bin = 0;
        else if (n2 >= maxx)
            bin = bins-1;
        else
            bin = (n2 - minx)/dbin;
        
        hist[bin]++;

        for (j = 0; j < corrLen; j++){
            corr[j] += n2*xold[j];
        }

        for (j = 0; j < corrLen-1; j++){
            xold[j] = xold[j+1];
        }

        xold[corrLen-1] = n2;

        if (i % 100000 == 0)
            printf("i = %d\n",i);

    }

    FILE * outfile = fopen("boxmuller.d","w");
    
    for (i = 0; i < bins; i++){
        fprintf(outfile,"%d %f %d %f\n",i,minx + i*dbin, hist[i],hist[i]/(double)N);
    }
    fclose(outfile);

    outfile = fopen("boxmuller_corr.d","w");

    for (i = 0; i < corrLen; i++){
        fprintf(outfile,"%d %f %f\n",i,corr[i], corr[i]/((N-i)));
    }

    fclose(outfile);
    free(hist);
}

void metropolis(){
    double n1,n2;
    int i,j,bin;
    double dbin = (maxx - minx)/bins;

    int * hist;
    double * corr, * xold;

    hist = (int*)malloc(sizeof(int)*bins);
    corr = (double*)malloc(sizeof(double)*corrLen);
    xold = (double*)malloc(sizeof(double)*corrLen);

    // initialize histogram bins
    for (i = 0; i < bins; i++){
        hist[i] = 0;
    }
    // initialize correlation array and array for temporarily saving values (empty)
    for (i = 0; i < corrLen; i++){
        xold[i] = 0;
        corr[i] = 0;
    }

    // initial position and probability
    double x1 = 0, x2;
    double p1 = exp(-0.5*x1*x1), p2;
    double r;
    for (i = 0; i  < N; i++){
        // MC step
        r = ranf() - 0.5;       // between -0.5 and 0.5
        x2 = x1 + r*stepSize;   // trial position

        p2 = exp(-0.5*x2*x2);   // and corresponding probability

        if (p2/p1 > 1){         // if the trial position is more probable, move
            x1 = x2;
            p1 = p2;
        } else {                // else, pick a random number
            r = ranf();
            if ((p2/p1)>r){     // if the ratio is larger than the random number, move
                x1 = x2;
                p1 = p2;
            }                   // else stay
        }


        // calculate histogram bin correspondning to position
        if (x1 < minx)              // lower limit
            bin = 0;
        else if (x1 >= maxx)        // upper limit
            bin = bins-1;
        else
            bin = (x1 - minx)/dbin;
        
        hist[bin]++;

        // Calculate correlation contribution from x1
        for (j = 0; j < corrLen; j++){
            corr[j] += x1*xold[j];
        }

        // move every old position one step to the right
        for (j = corrLen-1; j > 0; j--){
            xold[j] = xold[j-1];
        }
        
        // insert new position in the first position
        xold[0] = x1;
        
        if (i % 100000 == 0)
            printf("i = %d\n",i);
    }

    // histogram file
    FILE * outfile = fopen("metropolis.d","w");
    for (i = 0; i < bins; i++){
        fprintf(outfile,"%d %f %d %f\n",i,minx + i*dbin, hist[i],hist[i]/(double)N);
    }
    fclose(outfile);

    // correlation file
    outfile = fopen("metropolis_corr.d","w");
    for (i = 0; i < corrLen; i++){
        fprintf(outfile,"%d %f %f\n",i,corr[i], corr[i]/(N-i));
    }
    fclose(outfile);

    free(hist);
    free(corr);
    free(xold);
}

int main() {
    boxMuller();
    metropolis();    
}
