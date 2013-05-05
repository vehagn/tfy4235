#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int IBM = 4171;
double rinv = 0.5/( ((long int)1 << 31));

int n = 100, nh = 100, nlw = 100;
int itermax = 200, isamp = 0, nsamp = 500;
double egen, egan;

double eigenmax = 20.0, eigenmin = -20.0;
double avelambg = 0, avelamsm = 0;
double biglam, smalam;

double * vec, * vecp;
double * histMin, * histMax, * histLW;
double ** M, ** N, ** P;

double randIBMf();
void newSample();

int main(int argc, char *argv[]) {
    int i;
    for(i=0;i<1000;i++){
        randIBMf();
    }

    egen = (eigenmax - eigenmin)/(nh - 1);

    histMin = (double *)malloc(sizeof(double)*nh);
    histMax = (double *)malloc(sizeof(double)*nh);
    histLW = (double *)malloc(sizeof(double)*nlw);

    for(i=0;i<nh;i++){
        histMin[i] = 0;
        histMax[i] = 0;
    }
    for(i=0;i<nlw;i++){
        histLW[i] = 0;
    }

    M = (double **)malloc(sizeof(double*)*n);
    N = (double **)malloc(sizeof(double*)*n);
    P = (double **)malloc(sizeof(double*)*n);
    for(i=0;i<n;i++){
        M[i] = (double*)malloc(sizeof(double)*n);
        N[i] = (double*)malloc(sizeof(double)*n);
        P[i] = (double*)malloc(sizeof(double)*n);
    }

    vec = (double*)malloc(sizeof(double)*n);
    vecp = (double*)malloc(sizeof(double)*n);

    while(isamp < nsamp){
				newSample();
        printf("isamp = %d\n",isamp);
    }
		
    avelambg /= nsamp;
    avelamsm /= nsamp;
    egan = (avelambg-avelamsm)/(nlw-1);

		FILE * out1 = fopen("out1.d","w");
		FILE * out2 = fopen("out2.d","w");

    // Save bin, eigenvalue of bin, histogram of minimums, and of maximums
		for (i = 0; i < nh; i++){
				fprintf(out1,"%d %f %f %f\n",i, eigenmin + i*egen,histMin[i], histMax[i]);
		}
    // Save bin, eigenvalue of bin, and of density of states (cummulative distribution)
		for (i = 0; i < nlw; i++){
				fprintf(out2,"%d %f %f\n",i, avelamsm + i*egan , histLW[i]/nsamp);
		}
		fclose(out1);
		fclose(out2);
}

void newSample(){
    int i,j,iter;
    double sum;

    // Initializing initial random matrix
    for(i=0;i<n;i++){
        for(j=i;j<n;j++){
            M[i][j] = 2*randIBMf() - 1;
            M[j][i] = M[i][j];
        }
        M[i][i] += 1;     // adding a constant of 1 to diagonal
    }
    for(i=0;i<n;i++){
        vec[i] = 1;       // initializing initial vector, must be non-zero
    }

    for(iter=0;iter<itermax;iter++){
        for(i=0;i<n;i++){ 
            vecp[i] = 0;  // vector containing result of matrix-vector product
            for(j=0;j<n;j++){
                vecp[i] += M[i][j]*vec[j];
            }
        }
        sum = 0;
        for(i=0;i<n;i++){
            sum += vecp[i]*vecp[i];   // dot product
        }
        sum = 1.0/sqrt(sum);
        for(i=0;i<n;i++){
            vec[i]=vecp[i]*sum;       // vector for next iteration
        }

        biglam = 0;                   // variable storing the largest eigenvalue
        for(i=0;i<n;i++){
            biglam += vec[i]*vecp[i]; // dot product
        }

        biglam -= 1;                  // subtracting the constant of 1
    }

    int ih = (biglam - eigenmin)/egen;// calculating histogram bin
    histMax[ih]++;                    // corresponding to biglam

    for(i=0;i<n;i++){
        M[i][i] -= 1;                 // subtracting constant from diagonal
    }
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            N[i][j] = -M[i][j];       // New matrix for finding the smallest eigenvalue
        }
        N[i][i] += biglam;            // Adding constant to diagonal
    }
    for(i=0;i<n;i++){
        vec[i] = 1;                   // Initial vector
    }

    for(iter=0;iter<itermax;iter++){
        for(i=0;i<n;i++){
            vecp[i] = 0;
            for(j=0;j<n;j++){
                vecp[i] += N[i][j]*vec[j];
            }
        }
        double sum = 0;
        for(i=0;i<n;i++){
            sum += vecp[i]*vecp[i];
        }
        sum = 1.0/sqrt(sum);

        for(i=0;i<n;i++){
            vec[i]=vecp[i]*sum;
        }

        smalam = 0;
        for(i=0;i<n;i++){
            smalam += vec[i]*vecp[i];
        }
        smalam = biglam - smalam;       // Transforming from largest to smallest
    }
    ih = (smalam - eigenmin)/egen;
    histMin[ih] += 1;

    avelambg += biglam;                 // Average value contribution
    avelamsm += smalam;

    egan = (biglam - smalam)/nlw;

    int ilw,k;
    double alam;

    for(ilw = 0; ilw < nlw;ilw++){

        alam = smalam + egan*ilw;         // current eigenvalue scaled and transformed
        for(i=0;i<n;i++){
            for(j=0;j<n;j++){
                N[i][j] = M[i][j];
            }
        }

        int nshift = 1;
        double ann;
        for(k=n-1;k>1;k--){
            ann = 1.0/(N[k][k] - alam);   // denominator
            if (ann < 0)
                nshift++;

            for(i=0;i<k;i++){
                for(j=0;j<k;j++){
                    P[i][j] = N[i][j] - N[i][k]*N[k][j]*ann;  // new temp matrix.
                }
            }
            for(i=0;i<k;i++){
                for(j=0;j<k;j++){
                    N[i][j] = P[i][j];    // save new matrix
                }
            }
        }
        histLW[ilw] += nshift;            // increment DOS
    }

    isamp++;
}


double randIBMf(){
    IBM *= 16807;
    return (IBM*rinv + 0.5);
}

