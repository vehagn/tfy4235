#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <memory.h>
#include <ctime>

const double R = .5;
const double PI = 4.*atan(1);
const double XI = 2*R;

typedef struct{
	double xpos;
	double ypos;
	double r;
} Disc;

double rand0to1();
int randi(int max);

double discDistanceSq(Disc a, Disc b);
double distance1D(Disc a, Disc b);
int collisionCheck(Disc a, Disc b, double dx, double dy);

void placeDiscs(Disc a[], double Lx, double Ly, int N);
void bSort(int xsort[], int ysort[], Disc a[], int end);
void jiggleDisc(Disc a[], int xsort[], int ysort[], int N, int M, double Lx, double Ly);
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

	bSort(xsort, ysort, discs, N);	
	M = N + checkLeaks(discs, xsort, ysort, N, Lx, Ly);	
	int collisions = 0;
	for (int i = 0; i < N; i++){
		for (int j = i+1; j < M; j++){
			if (collisionCheck(discs[i],discs[j],0,0)){
				collisions++;
				printf("%i,%i\t(%4.2f,%4.2f)\t(%4.2f,%4.2f)\n",i,j,discs[i].xpos,discs[i].ypos,discs[j].xpos,discs[j].ypos);
			}
		}
	}
	printf("Collisions %i\n",collisions);
	
	printf("M: %i\n",M);
	printDiscs(discs, N, M, 0);

	int jiggles = 5000;
	for(int i = 0; i < jiggles; i++){
		//printf("%i \n",i);
		if(!(i%(jiggles/100))){printf("%2.0f \n",100.*(i/(double)jiggles));}
		jiggleDisc(discs, xsort, ysort, N, M, Lx, Ly);
	}	
	
	bSort(xsort, ysort, discs, N);	
	M = N + checkLeaks(discs, xsort, ysort, N, Lx, Ly);	
	collisions = 0;
	for (int i = 0; i < N; i++){
		for (int j = i+1; j < M; j++){
			if (collisionCheck(discs[i],discs[j],0,0)){
				collisions++;
				printf("(%4.2f,%4.2f)\t(%4.2f,%4.2f)\n",discs[i].xpos,discs[i].ypos,discs[j].xpos,discs[j].ypos);
			}
		}
	}
	printf("Collisions %i\n",collisions);
	
	printf("\nM: %i\n",M);
	printDiscs(discs, N, M, 1);
	
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
	if (((a.xpos+dx-b.xpos) < 2*R) || ((a.ypos+dy-b.ypos) < 2*R)){
		if (((a.xpos+dx-b.xpos)*(a.xpos+dx-b.xpos)+(a.ypos+dy-b.ypos)*(a.ypos+dy-b.ypos)) < (a.r+b.r)*(a.r+b.r)){return 1;}
		else {return 0;}
	}else {return 0;}
}

void placeDiscs(Disc a[], double Lx, double Ly, int N){
	a[0].r = R;
	a[0].xpos = Lx*rand0to1();
	a[0].ypos = Ly*rand0to1();
	int M = N;
	for (int i = 1; i < N; i++){
		a[i].r =R;	
		a[i].xpos = Lx*rand0to1();
		a[i].ypos = Ly*rand0to1();
		for (int j = 0; j < i; j++){
			//General collision
			if (collisionCheck(a[i],a[j],0,0)){i--;break;}
			//Side collision
			if (a[i].xpos + 2*R > Lx){
				if (collisionCheck(a[i],a[j],-Lx,0)){i--;break;}
			}else if (a[i].xpos - 2*R < 0){
				if (collisionCheck(a[i],a[j],Lx,0)){i--;break;}
			}else if (a[i].ypos + 2*R > Ly){
				if (collisionCheck(a[i],a[j],0,-Ly)){i--;break;}
			}else if (a[i].ypos - 2*R < 0){
				if (collisionCheck(a[i],a[j],0,Ly)){i--;break;}
			}
			//Corner collision
			if (((a[i].xpos +2*R) > Lx)&&((a[i].ypos + 2*R) > Ly)){
				if (collisionCheck(a[i],a[j],-Lx,0)){i--;break;}
				if (collisionCheck(a[i],a[j],0,-Ly)){i--;break;}
				if (collisionCheck(a[i],a[j],-Lx,-Ly)){i--;break;}
			}else if(((a[i].xpos +2*R) > Lx)&&((a[i].ypos - 2*R) < 0)){
				if (collisionCheck(a[i],a[j],-Lx,0)){i--;break;}
				if (collisionCheck(a[i],a[j],0,Ly)){i--;break;}
				if (collisionCheck(a[i],a[j],-Lx,Ly)){i--;break;}
			}else if (((a[i].xpos -2*R) < 0)&&((a[i].ypos - 2*R) < 0)){
				if (collisionCheck(a[i],a[j],Lx,0)){i--;break;}
				if (collisionCheck(a[i],a[j],0,Ly)){i--;break;}
				if (collisionCheck(a[i],a[j],Lx,Ly)){i--;break;}
			}else if (((a[i].xpos -2*R) < 0)&&((a[i].ypos + 2*R) > Ly)){
				if (collisionCheck(a[i],a[j],Lx,0)){i--;break;}
				if (collisionCheck(a[i],a[j],0,-Ly)){i--;break;}
				if (collisionCheck(a[i],a[j],Lx,-Ly)){i--;break;}
			}
		}
	}	 

}


void bSort(int xsort[], int ysort[], Disc a[], int N){
	int temp;
	int j;
	for (int i = 0; i < N-1; i++){
		j = i;
		while (a[xsort[j]].xpos > a[xsort[j+1]].xpos){
			temp = xsort[j];
			xsort[j] = xsort[j+1];
			xsort[j+1] = temp;
			j = j-1;
			if (j < 0) {break;}
		}
		j = i;
		while (a[ysort[j]].ypos > a[ysort[j+1]].ypos){
			temp = ysort[j];
			ysort[j] = ysort[j+1];
			ysort[j+1] = temp;
			j = j-1;
			if (j < 0) {break;}
		}
	}
}

int checkLeaks(Disc a[], int xsort[], int ysort[], int N, double Lx, double Ly){
	double xdispl = 0, ydispl = 0;
	int lleaks = 0, rleaks = N-1;
	int bleaks = 0, tleaks = N-1;
	int cleaks = 0; // corner leaks
	int leaks = 0;
	int M = N-1;
	//-----------------------------------
	while(a[xsort[lleaks]].xpos < R){ //check left side
		lleaks++;
	}
	leaks += lleaks;
	while (lleaks > 0){
		--lleaks;
		a[++M].r = a[xsort[lleaks]].r;
		a[M].xpos = a[xsort[lleaks]].xpos + Lx;
		a[M].ypos = a[xsort[lleaks]].ypos;
		if (a[M].ypos < R){
			a[++M].r = a[xsort[lleaks]].r;
			a[M].xpos = a[xsort[lleaks]].xpos + Lx;
			a[M].ypos = a[xsort[lleaks]].ypos + Ly;
			cleaks++;
		}else if (a[M].ypos > Ly-R){
			a[++M].r = a[xsort[lleaks]].r;
			a[M].xpos = a[xsort[lleaks]].xpos + Lx;
			a[M].ypos = a[xsort[lleaks]].ypos - Ly;
			cleaks++;
		}	
	}
	//------------------------------------
	while(a[xsort[rleaks]].xpos > Lx-R){
		rleaks--;
	}
	rleaks = N-1 - rleaks;
	leaks += rleaks;
	while(rleaks > 0){
		a[++M].r = a[xsort[N-rleaks]].r;
		a[M].xpos = a[xsort[N-rleaks]].xpos - Lx;
		a[M].ypos = a[xsort[N-rleaks]].ypos;
		if (a[M].ypos < R){
			a[++M].r = a[xsort[N-rleaks]].r;
			a[M].xpos = a[xsort[N-rleaks]].xpos - Lx;
			a[M].ypos = a[xsort[N-rleaks]].ypos + Ly;
			cleaks++;
		}else if (a[M].ypos > Ly-R){
			a[++M].r = a[xsort[N-rleaks]].r;
			a[M].xpos = a[xsort[N-rleaks]].xpos - Lx;
			a[M].ypos = a[xsort[N-rleaks]].ypos - Ly;
			cleaks++;
		}	
		
		--rleaks;	
	}
	//-----------------------------------
	while(a[ysort[bleaks]].ypos < R){
		bleaks++;
	}
	leaks += bleaks;
	while (bleaks > 0){
		--bleaks;
		a[++M].r = a[ysort[bleaks]].r;
		a[M].xpos = a[ysort[bleaks]].xpos;
		a[M].ypos = a[ysort[bleaks]].ypos + Ly;
	}
	//-----------------------------------
	while(a[ysort[tleaks]].ypos > Ly-R){
		tleaks--;
	}
	tleaks = N-1 - tleaks;
	leaks += tleaks;
	while(tleaks > 0){
		a[++M].r = a[ysort[N-tleaks]].r;
		a[M].xpos = a[ysort[N-tleaks]].xpos;
		a[M].ypos = a[ysort[N-tleaks]].ypos - Ly;
		--tleaks;
	}
	leaks += cleaks;
	return leaks;
}

void jiggleDisc(Disc a[], int xsort[], int ysort[], int N, int M, double Lx, double Ly){
	double r = XI*rand0to1();
	double theta = 2*PI*rand0to1();
	
	double dx = r*cos(theta);
	double dy = r*sin(theta);

	int x = randi(N);
	int y = 0;
	//printf("Find y pos\n");
	//Find corresponding y-sort position of x-sort pos.
	while(a[xsort[x]].ypos > a[ysort[y]].ypos)
		y++;
		
	double oldx = a[xsort[x]].xpos;
	double oldy = a[xsort[x]].ypos;

	a[xsort[x]].xpos += dx;
	a[xsort[x]].ypos += dy;

	if (a[xsort[x]].xpos > Lx){a[xsort[x]].xpos -= Lx;}
	if (a[xsort[x]].xpos < 0){a[xsort[x]].xpos += Lx;}
	if (a[xsort[x]].ypos > Ly){a[xsort[x]].ypos -= Ly;}
	if (a[xsort[x]].ypos < 0){a[xsort[x]].ypos += Ly;}	

	for (int i = 0; i < N; i++){
		for (int j = 0; j < i; j++){
			//General collision
			if (collisionCheck(a[i],a[j],0,0)){a[xsort[x]].xpos=oldx;a[xsort[x]].ypos=oldy;return;}
			//Side collision
			if (a[i].xpos + 2*R > Lx){
				if (collisionCheck(a[i],a[j],-Lx,0)){a[xsort[x]].xpos=oldx;a[xsort[x]].ypos=oldy;return;}
			}else if (a[i].xpos - 2*R < 0){
				if (collisionCheck(a[i],a[j],Lx,0)){a[xsort[x]].xpos=oldx;a[xsort[x]].ypos=oldy;return;}
			}else if (a[i].ypos + 2*R > Ly){
				if (collisionCheck(a[i],a[j],0,-Ly)){a[xsort[x]].xpos=oldx;a[xsort[x]].ypos=oldy;return;}
			}else if (a[i].ypos - 2*R < 0){
				if (collisionCheck(a[i],a[j],0,Ly)){a[xsort[x]].xpos=oldx;a[xsort[x]].ypos=oldy;return;}
			}
			//Corner collision
			if (((a[i].xpos +2*R) > Lx)&&((a[i].ypos + 2*R) > Ly)){
				if (collisionCheck(a[i],a[j],-Lx,0)){a[xsort[x]].xpos=oldx;a[xsort[x]].ypos=oldy;return;}
				if (collisionCheck(a[i],a[j],0,-Ly)){a[xsort[x]].xpos=oldx;a[xsort[x]].ypos=oldy;return;}
				if (collisionCheck(a[i],a[j],-Lx,-Ly)){a[xsort[x]].xpos=oldx;a[xsort[x]].ypos=oldy;return;}
			}else if(((a[i].xpos +2*R) > Lx)&&((a[i].ypos - 2*R) < 0)){
				if (collisionCheck(a[i],a[j],-Lx,0)){a[xsort[x]].xpos=oldx;a[xsort[x]].ypos=oldy;return;}
				if (collisionCheck(a[i],a[j],0,Ly)){a[xsort[x]].xpos=oldx;a[xsort[x]].ypos=oldy;return;}
				if (collisionCheck(a[i],a[j],-Lx,Ly)){a[xsort[x]].xpos=oldx;a[xsort[x]].ypos=oldy;return;}
			}else if (((a[i].xpos -2*R) < 0)&&((a[i].ypos - 2*R) < 0)){
				if (collisionCheck(a[i],a[j],Lx,0)){a[xsort[x]].xpos=oldx;a[xsort[x]].ypos=oldy;return;}
				if (collisionCheck(a[i],a[j],0,Ly)){a[xsort[x]].xpos=oldx;a[xsort[x]].ypos=oldy;return;}
				if (collisionCheck(a[i],a[j],Lx,Ly)){a[xsort[x]].xpos=oldx;a[xsort[x]].ypos=oldy;return;}
			}else if (((a[i].xpos -2*R) < 0)&&((a[i].ypos + 2*R) > Ly)){
				if (collisionCheck(a[i],a[j],Lx,0)){a[xsort[x]].xpos=oldx;a[xsort[x]].ypos=oldy;return;}
				if (collisionCheck(a[i],a[j],0,-Ly)){a[xsort[x]].xpos=oldx;a[xsort[x]].ypos=oldy;return;}
				if (collisionCheck(a[i],a[j],Lx,-Ly)){a[xsort[x]].xpos=oldx;a[xsort[x]].ypos=oldy;return;}
			} 
		}
	}
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
	fclose(discFile);
	
	sprintf(fname,"ghosts%i.dat",num);
	ghostFile = fopen(fname,"w");
	for (int i = N; i < M; i++){
		fprintf(ghostFile,"%f\t%f\t%f\n",a[i].xpos,a[i].ypos,a[i].r);
		//printf("%f\t%f\t%f\n",a[i].xpos,a[i].ypos,a[i].r);
	}
	fclose(ghostFile);
}
