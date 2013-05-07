#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <memory.h>
#include <ctime>

const double R = 1.;
const double PI = 4.*atan(1);
const double XI = R/10;

typedef struct{
	double xpos;
	double ypos;
	double r;
} Disc;

double rand0to1();
int randi(int max);

double discDistanceSq(Disc a, Disc b);
int collisionCheck(Disc a, Disc b, double dx, double dy);

void placeDiscs(Disc a[], double Lx, double Ly, int N);
void bSort(int sort[], Disc a[], int start, int end);
void jiggleDisc(Disc a[], int xsort[], int ysort[], int N, double Lx, double Ly);
void pushDisc(Disc a[], int xsort[], int ysort[], int n, int N);
int checkLeaks(Disc a[], int xsort[], int ysort[], int N, double Lx, double Ly);

void printDiscs(Disc a[], int N, int M, int num);

int main (int argc, char** argv){
	srand((unsigned)time(NULL));
	int N = 4;			//Number of discs
	double Lx = 10, Ly = 10;	//Dimensions
	
	if (argc > 1){
		N = atoi(argv[1]);	
		if (argc > 3){
			Lx = atof(argv[2]);
			Ly = atof(argv[3]);
		}
	}
	int M;

	printf("Density: %f\n",N*PI*R*R/(Lx*Ly));
	Disc discs[2*N];
	placeDiscs(discs, Lx, Ly, N);
	printf("Discs placed!\n");

	int xsort[N], ysort[N];
	for (int i = 0; i < N; i++){
		xsort[i] = i;
		ysort[i] = i;
	}

	bSort(xsort, discs, 0, N);
	bSort(ysort, discs, 0, N);	
	M = N + checkLeaks(discs, xsort, ysort, N, Lx, Ly);	

	/*printf("xpos:\t");
	for (int i = 0; i < N; i++){
		printf("%4.2f\t",discs[xsort[i]].xpos);
	}
	printf("\nypos:\t");
	for (int i = 0; i < N; i++){
		printf("%4.2f\t",discs[ysort[i]].ypos);
	}
	printf("\ncoords:\t");
	for (int i = 0; i < N; i++)
		printf("%i:(%4.2f,%4.2f) ",xsort[i],discs[xsort[i]].xpos,discs[xsort[i]].ypos);
	printf("\n");
	
	for (int i = 0; i < N; i++){
		for (int j = i; j < N; j++){
			printf("Distance Sq %i %i: %f\n",i,j,discDistanceSq(discs[i],discs[j]));		
		}
	}*/
	
	int collisions = 0;
	for (int i = 0; i < N; i++){
		for (int j = i+1; j < N; j++){
			collisions += collisionCheck(discs[i],discs[j],0,0);
		}
	}
	printf("Collisions %i\n",collisions);
	
	printDiscs(discs, N, M, 0);
	for(int i = 0; i < 50000; i++){
		jiggleDisc(discs, xsort, ysort, N, Lx, Ly);
	}	

	printDiscs(discs, N, M, 1);

	for (int i = 0; i < N; i++){
		for (int j = i+1; j < N; j++){
			collisions += collisionCheck(discs[i],discs[j],0,0);
		}
	}
	printf("Collisions %i\n",collisions);
	return 0;
}

double rand0to1(){
	return ((double)rand()/(double)RAND_MAX);
}
int randi(int max){
	return rand()%max;
}

double discDistanceSq(Disc a, Disc b){
	return (a.xpos - b.xpos)*(a.xpos - b.xpos) + (a.ypos - b.ypos)*(a.ypos - b.ypos);
}

int collisionCheck(Disc a, Disc b, double dx, double dy){
	if (((a.xpos+dx-b.xpos) < (a.r+b.r)) || ((a.ypos+dy-b.ypos) < (a.r+b.r))){
		if (((a.xpos+dx-b.xpos)*(a.xpos+dx-b.xpos)+(a.ypos+dy-b.ypos)*(a.ypos+dy-b.ypos)) < (a.r+b.r)*(a.r+b.r)){
			//printf("%f\n",(((a.xpos+dx-b.xpos)*(a.xpos+dx-b.xpos)+(a.ypos+dy-b.ypos)*(a.ypos+dy-b.ypos))));
			return 1;
		}else {return 0;}
	}else {return 0;}
}

void placeDiscs(Disc a[], double Lx, double Ly, int N){
	a[0].r = R;
	a[0].xpos = Lx*rand0to1();
	a[0].ypos = Ly*rand0to1();
	for (int i = 1; i < N; i++){
		a[i].r =R;	
		a[i].xpos = Lx*rand0to1();
		a[i].ypos = Ly*rand0to1();
		for (int j = 0; j < i; j++){
		/*	if (((a[i].xpos - a[j].xpos)*(a[i].xpos - a[j].xpos) + 
				(a[i].ypos - a[j].ypos)*(a[i].ypos - a[j].ypos)) < (a[i].r+a[j].r)*(a[i].r+a[j].r)){*/
			if (collisionCheck(a[i],a[j],0,0)){
				i--;
				//printf("COLLISION!!!\t %i\t %4.2f\t %4.2f\n",i,a[i].xpos,a[i].ypos);
				break;
			}else if (a[i].xpos + 2*R > Lx){
				if (collisionCheck(a[i],a[j],-Lx,0)){
					i--;
					//printf("R collision\n");
					break;
				}
			}else if (a[i].xpos - 2*R < 0){
				if (collisionCheck(a[i],a[j],Lx,0)){
					i--;
					//printf("L collision\n");
					break;
				}
			}else if (a[i].ypos + 2*R > Ly){
				if (collisionCheck(a[i],a[j],0,-Ly)){
					i--;
					//printf("T collision\n");
					break;	
				}else if (a[i].xpos + 2*R > Lx){
					if (collisionCheck(a[i],a[j],-Lx,Ly)){
						i--;
						//printf("T R collision\n");
						break;
					}
				}else if (a[i].xpos - 2*R < 0){
					if (collisionCheck(a[i],a[j],Lx,Ly)){					
						i--;
						//printf("T L collision\n");
						break;
					}
				}
			}else if (a[i].ypos - 2*R < 0){
				if (collisionCheck(a[i],a[j],0,Ly)){
					i--;
					//printf("B collision\n");
					break;
				}else if (a[i].xpos + 2*R > Lx){
					if (collisionCheck(a[i],a[j],-Lx,Ly)){
						i--;
						//printf("B R collision\n");
						break;
					}
				}else if (a[i].xpos - 2*R < 0){
					if (collisionCheck(a[i],a[j],Lx,Ly)){					
						i--;
						//printf("B T collision\n");
						break;
					}
				}
			} 
		}
	}	 

}


void bSort(int xsort[], Disc a[], int start, int N){
	int temp;
	int j;
	for (int i = start; i < N-1; i++){
		j = i;
		while (a[xsort[j]].xpos > a[xsort[j+1]].xpos){
			temp = xsort[j];
			xsort[j] = xsort[j+1];
			xsort[j+1] = temp;
			j = j-1;
			if (j < start) {break;}
		}
	}
}

int checkLeaks(Disc a[], int xsort[], int ysort[], int N, double Lx, double Ly){
	double xdispl = 0, ydispl = 0;
	int lleaks = 0, rleaks = N-1;
	int bleaks = 0, tleaks = N-1;
	int leaks = 0;
	int M = N-1;
	//-----------------------------------
	while(a[xsort[lleaks]].xpos < R){ //check left side
		lleaks++;
	}
	leaks += lleaks;
	//printf("lleaks %i\n",lleaks);
	while (lleaks > 0){
		--lleaks;
		a[++M].r = a[xsort[lleaks]].r;
		a[M].xpos = a[xsort[lleaks]].xpos + Lx;
		a[M].ypos = a[xsort[lleaks]].ypos;
		//printf("Left leak index: %i\n",xsort[lleaks]);	
	}
	//------------------------------------
	while(a[xsort[rleaks]].xpos > Lx-R){
		rleaks--;
	}
	rleaks = N-1 - rleaks;
	leaks += rleaks;
	//printf("rleaks %i\n",rleaks);
	while(rleaks > 0){
		a[++M].r = a[xsort[N-rleaks]].r;
		a[M].xpos = a[xsort[N-rleaks]].xpos - Lx;
		a[M].ypos = a[xsort[N-rleaks]].ypos;
		//printf("Right leak index: %i\n",xsort[N-rleaks]);
		--rleaks;	
	}
	//-----------------------------------
	while(a[ysort[bleaks]].ypos < R){
		bleaks++;
	}
	//printf("bleaks %i\n",bleaks);
	leaks += bleaks;
	while (bleaks > 0){
		--bleaks;
		a[++M].r = a[ysort[bleaks]].r;
		a[M].xpos = a[ysort[bleaks]].xpos;
		a[M].ypos = a[ysort[bleaks]].ypos + Ly;
		//printf("Bottom leak index: %i\n",ysort[bleaks]);
	}
	//-----------------------------------
	while(a[ysort[tleaks]].ypos > Ly-R){
		tleaks--;
	}
	tleaks = N-1 - tleaks;
	leaks += tleaks;
	//printf("tleaks %i\n",tleaks);
	while(tleaks > 0){
		a[++M].r = a[ysort[N-tleaks]].r;
		a[M].xpos = a[ysort[N-tleaks]].xpos;
		a[M].ypos = a[ysort[N-tleaks]].ypos - Ly;
		//printf("Top leaks index: %i \n",ysort[N-tleaks]);
		--tleaks;
	}
	printf("Leaks: %i\n",leaks);
	

	return leaks;
}

void jiggleDisc(Disc a[], int xsort[], int ysort[], int N, double Lx, double Ly){
	double r = XI*rand0to1();
	double theta = 2*PI*rand0to1();
	
	double dx = r*cos(theta);
	double dy = r*sin(theta);

	int x = randi(N);
	int y = 0;
	//Find corresponding y-sort position of x-sort pos.
	while(a[xsort[x]].ypos > a[ysort[y]].ypos)
		y++;
	
	int negx = 0, posx = 0;

	double dist = 0;
	while (dist < 2*R + r){
		negx++;
		dist = a[xsort[x]].xpos-a[xsort[(x-negx+2*N)%N]].xpos;
		if (dist < 0){dist += Lx;}
		//printf("negx: %f\n",dist);
	}
	dist = 0;
	while (dist < 2*R + r){
		posx++;
		dist = a[xsort[(x+posx)%N]].xpos-a[xsort[x]].xpos;
		if (dist < 0){dist += Lx;}
		//printf("posx: %f\n",dist);
	}
	
	int negy = 0, posy = 0;

	dist = 0;
	while (dist < 2*R + r){
		negy++;
		dist = a[ysort[y]].ypos-a[ysort[(y-negy+2*N)%N]].ypos;
		if (dist < 0){dist += Ly;}
		//printf("negy: %f\n",dist);
	}
	dist = 0;
	while (dist < 2*R + r){
		posy++;
		dist = a[ysort[(y+posy)%N]].ypos-a[ysort[y]].ypos;
		if (dist < 0){dist += Ly;}
		//printf("posy: %f\n",dist);
	}

	bool collision = false;
	//printf("test for collisions\n");
	for (int i = x-negx; i < x; i++){
		if (collisionCheck(a[xsort[(i+2*N)%N]],a[xsort[x]],dx,dy))
			collision = true;
	}
	for (int i = x+1; i <= x+posx; i++){
		if (collisionCheck(a[xsort[i%N]],a[xsort[x]],dx,dy))
			collision = true;
	}
	
	for (int i = y-negy; i < y; i++){
		if (collisionCheck(a[ysort[(i+2*N)%N]],a[ysort[y]],dx,dy))
			collision = true;
	}
	for (int i = y+1; i <= y+posy; i++){
		if (collisionCheck(a[ysort[i%N]],a[ysort[y]],dx,dy))
			collision = true;
	}
	if (!collision){
		a[xsort[x]].xpos += dx;
		a[xsort[x]].ypos += dy;
		if ((x+posx<N) && ((x-negx) > 0)){bSort(xsort,a,x-negx,x+posx);}
		else {bSort(xsort,a,0,N);}
		if ((y+posy<N) && ((y-negy) > 0)){bSort(ysort,a,y-negy,x+posy);}
		else {bSort(xsort,a,0,N);}
		//printf("JIGGLE!\n");
	}else{}
}



void pushDisc(Disc a[], int xsort[], int ysort[], int n, int N){

}

void printDiscs(Disc a[], int N, int M, int num){
	FILE *discFile;
	FILE *ghostFile;
	char fname[16];
	sprintf(fname,"discs%i.dat",num);
	discFile = fopen(fname,"w");
	for (int i = 0; i < N; i++){
		fprintf(discFile,"%f\t%f\t%f\n",a[i].xpos,a[i].ypos,a[i].r);
	}
	
	ghostFile = fopen("ghosts.dat","w");
	for (int i = N; i < M; i++){
		fprintf(ghostFile,"%f\t%f\t%f\n",a[i].xpos,a[i].ypos,a[i].r);
	}
}
