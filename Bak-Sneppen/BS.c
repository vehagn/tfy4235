// Bak-Sneppen.cpp : Defines the entry point for the console application.
//

#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <ctime>
#include <cmath>

const int species = (1<<9)+2;
const int generations = (1<<17);

// Process has done i out of n rounds,
// and we want a bar of width w and resolution r.
static inline void loadBar(int x, int n, int r, int w)
{
    // Only update r times.
    if ( x % (n/r) != 0 ) return;
 
    // Calculuate the ratio of complete-to-incomplete.
    float ratio = x/(float)n;
    int   c     = ratio * w;
 
    // Show the percentage complete.
    printf("%3d%% [", (int)(ratio*100) );
 
    // Show the load bar.
    for (int x=0; x<c; x++)
       printf("=");
 
    for (int x=c; x<w; x++)
       printf(" ");
 
    // ANSI Control codes to go back to the
    // previous line and clear it.
    printf("]\n33[F33[J\r");
}


int main(int argc, char* argv[])
{
	double *fitness = (double*) malloc((species*generations+1)*sizeof(double));
	double *lifetime = (double*) malloc((species*generations+1)*sizeof(double));
	int *change = (int*) malloc(generations*sizeof(int));
	int minIndex;
	double min;

	printf("Memory: %lu MB\n", (2*species*generations*sizeof(double))>>20);

	FILE *fitfile;
	fitfile = fopen("fitness.dat","w");
	FILE *lifefile;	
	lifefile = fopen("lifetime.dat","w");
	FILE *changefile;
	changefile = fopen("change.dat","w");

	srand((unsigned)time(NULL));

	for(int i = 0; i < species; i++ ){
		fitness[i] = ((double)rand()/(double)RAND_MAX);
		lifetime[i] = 0;
	}
	
	int row;
	printf("Calculating: ");
	for(int i = 0; i < generations-1; i++){
		loadBar(i,generations,100,20);
		row = species*i;
		min = 1;
		fitness[row] = fitness[row + species-1];
		fitness[row + species-1] = fitness[row + 1];
		for(int j = 0+1; j < species-1; j++){
			lifetime[row + j] += 1/(double)generations;
			if(min > fitness[row + j]){
				min = fitness[row + j];
				minIndex = j;
			}
		}
		change[i] = minIndex;

		fitness[row + minIndex-1] = ((double)rand()/(double)RAND_MAX);
		fitness[row + minIndex] = ((double)rand()/(double)RAND_MAX);
		fitness[row + minIndex+1] = ((double)rand()/(double)RAND_MAX);

		lifetime[row + minIndex-1] = 0;
		lifetime[row + minIndex] = 0;
		lifetime[row + minIndex+1] = 0;
		
		for(int k = 0; k < species; k++ ){
			fitness[(i+1)*species + k] = fitness[row + k];
			lifetime[(i+1)*species + k] = lifetime[row + k];
		}
	}

	printf("\nCalculations done, saving results...");

	for(int i = 0; i < generations; i++){
		loadBar(i,generations,100,20);
		row = species*i;
		for(int j = 0+1; j < species-1; j++){
			//fprintf(fitfile,"%1.16f\t",fitness[row + j]);
			fprintf(lifefile,"%1.16f\t",lifetime[row + j]);
			}
		//fprintf(fitfile,"\n");
		fprintf(lifefile,"\n");
		fprintf(changefile,"%i\t%i\t\n",change[i]-1,i);
		fprintf(changefile,"%i\t%i\t\n",change[i],i);
		fprintf(changefile,"%i\t%i\t\n",change[i]+1,i);
	}
	
	printf("...Results saved!\n");
	fclose(fitfile); fclose(lifefile); fclose(changefile);
	return 0;
}

