#include<stdio.h>
#include<stdlib.h>
#include<sys/time.h>
#include<math.h>

// RNG Variables and constants
const double RINV = 0.5/(1 << 31);      // 0.5/(2^31)
int IBM;

// Variables and constants
int N;
int dtmax;

int mini;
double minr;
double t_tmp;
double time_init, time_sim, time_out;

// Arrays
double * r;     // Array for storing fitness values
int * mins;     // Array for storing which species was least fit


double walltime ( void ) {  // Returns current time, with microsecond precision
    static struct timeval t;
    gettimeofday ( &t, NULL );
    return ( t.tv_sec + 1e-6 * t.tv_usec );
}

int randIBM() {      // Random integer between -2^31 and 2^31 - 1
    IBM*=16807;     // google: LCG - Linear Congruential Generator
    return IBM;
}

double randIBMf() {  // Random number between 0 and 1, [0,1)
    return randIBM()*RINV + 0.5;
}

int main (int argc, char ** argv) {
    int i,j;

    // INITIALIZATIONS
	t_tmp = walltime();
    for (i = 0; i < argc; i++) {    
        printf("%s ",argv[i]);
    }
    printf("\n");
    
    if (argc == 1) {
        printf("No arguments from command line. Using default values\n");
        N = 64;
        dtmax = 10000;
    } else if (argc == 3) {
        N = atoi(argv[1]);      // First argument, number of species
        dtmax = atoi(argv[2]);  // Second argument, number of iterations
    } else {
        printf("Invalid number of arguments from command line, aborting.\n");
        exit(-1);
    }
    printf(" N = %d\n",N);
    printf(" dtmax = %d\n\n",dtmax);

    FILE * f1 = fopen("out.d","w");

    // Initializing RNG
    IBM = 2*time(0) + 1;
    for (i = 0; i < 1000; i++)
        randIBM();

    r = (double*)malloc(sizeof(double)*N);
    mins = (int*)malloc(sizeof(int)*dtmax);
    
    // Initialize fitness array
    for(i = 0; i < N; i++){
        r[i] = randIBMf();
    }
    time_init = walltime() - t_tmp;

    // MAIN LOOPS
    t_tmp = walltime();
    for(j = 0;j < dtmax; j++){
        minr = 1;
        for(i=0;i<N;i++){
            if (r[i] < minr){
                minr = r[i];
                mini = i;
            }
        }

        mins[j] = mini;
        int right = mini + 1;
        int left = mini - 1;

        if(mini == N-1)
            right = 0;
        if(mini == 0)
            left = N-1;

        r[mini] = randIBMf();
        r[left] = randIBMf();
        r[right] = randIBMf();
    }
    time_sim = walltime() - t_tmp;

    // PRINT TO FILE
    t_tmp = walltime();
    for (i = 0; i < dtmax; i++){
        fprintf(f1,"%d    %d\n",mins[i], i);
    }
    time_out = walltime() - t_tmp;

    // TIMINGS
    printf("time init: %f\n",time_init);
    printf("time sim: %f\n",time_sim);
    printf("time out: %f\n",time_out);

    // FREE MEMORY
    free(r);
    free(mins);
    fclose(f1);
    return 0;
}

