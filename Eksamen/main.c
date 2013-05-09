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
double distance1D(Disc a, Disc b);
int collisionCheck(Disc a, Disc b, double dx, double dy, double r);

void placeDiscs(Disc a[], double Lx, double Ly, int N);
void bSort(int xsort[], int ysort[], Disc a[], int end);
int jiggleDisc(Disc a[], int xsort[], int ysort[], int N, int M, double Lx, double Ly, int **map);
void pushDisc(Disc a[], int xsort[], int ysort[], int N, double l, double Lx, double Ly, int **map);
int checkLeaks(Disc a[], int xsort[], int ysort[], int N, double Lx, double Ly);

int checkConcentric(Disc a[], Disc b[], N, double r_max, int m);

void printDiscs(Disc a[], int N, int M, int num);
void printMap(int **map, int iLx, int iLy);

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
	int iLx = (int)Lx, iLy = (int)Ly;
	
	int **map;
	map = (int **)malloc(iLx*sizeof(int*));
	map[0] = (int *)malloc(iLx*iLy*sizeof(double));
	for (int i=1; i < iLx; i++){
		map[i] = map[i-1] + iLy;
	}

	for (int i = 0; i < iLx; i++){
		for (int j = 0; j < iLy; j++)
			map[i][j] = 0;
	}

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
			if (collisionCheck(discs[i],discs[j],0,0,0)){
				collisions++;
				printf("%i,%i\t(%4.2f,%4.2f)\t(%4.2f,%4.2f)\n",i,j,discs[i].xpos,discs[i].ypos,discs[j].xpos,discs[j].ypos);
			}
		}
	}
	printf("Collisions %i\n",collisions);
	
	printf("M: %i\n",M);
	printDiscs(discs, N, M, 0);

	int jiggles = 0 << 18 ;
	for(int i = 0; i < jiggles; i++){
		//printf("%i \n",i);
		if(!(i%(jiggles/100))){printf("%03.0f \n",100.*(i/(double)jiggles));}
		jiggleDisc(discs, xsort, ysort, N, M, Lx, Ly, map);
	}	
	
	int pushes = 1 << 20; double l = 10*R;
	for(int i = 0; i < pushes; i++){
		if(!(i%(pushes/100))){printf("%03.0f\n",100.*(i/(double)pushes));}
		pushDisc(discs, xsort, ysort, N, l, Lx, Ly, map); 
	}

	bSort(xsort, ysort, discs, N);	
	M = N + checkLeaks(discs, xsort, ysort, N, Lx, Ly);	
	collisions = 0;
	for (int i = 0; i < N; i++){
		for (int j = i+1; j < M; j++){
			if (collisionCheck(discs[i],discs[j],0,0,0)){
				collisions++;
	//			printf("(%4.2f,%4.2f)\t(%4.2f,%4.2f)\n",discs[i].xpos,discs[i].ypos,discs[j].xpos,discs[j].ypos);
			}
		}
	}
	printf("Collisions %i\n",collisions);
	
	printf("M: %i\n",M);
	printDiscs(discs, N, M, 1);
	printMap(map,iLx,iLy);	

	/*for(int i = 0; i < N; i++){
		printf("%4.2f\t",discs[xsort[i]].xpos);
	}*/
	printf("\n");
	return 0;
}

double rand0to1(){
	return ((double)rand()/(double)RAND_MAX);
}
int randi(int max){
	return rand()%max;
}

int checkConcentric(Disc a[], Disc b[], N, double r_max, int m){
	int inside[m];
	for (int c = 0; c < m; c++){
		inside[c] = 0;
		r = (c/m)*r_max;
		for (int i = 0; i < N; i++){
			for (int j = 0; j < i; j++){
				//General collision
				if (collisionCheck(a[i],a[j],0,0,r)){inside++;break}
				//Side collision
				if (a[i].xpos + 2*R > Lx){
					if (concentricCheck(a[i],a[j],-Lx,0,r)){inside++;break}
				}else if (a[i].xpos - 2*R < 0){
					if (concentricCheck(a[i],a[j],Lx,0,r)){inside++;break}
				}else if (a[i].ypos + 2*R > Ly){
					if (concentricCheck(a[i],a[j],0,-Ly,r)){inside++;break;}
				}else if (a[i].ypos - 2*R < 0){
					if (concentricCheck(a[i],a[j],0,Ly,r)){inside++;break;}
				}
				//Corner collision
				if (((a[i].xpos +2*R) > Lx)&&((a[i].ypos + 2*R) > Ly)){
					if (concentricCheck(a[i],a[j],-Lx,0,r)){inside++;break;}
					if (concentricCheck(a[i],a[j],0,-Ly,r)){inside++;break;}
					if (concentricCheck(a[i],a[j],-Lx,-Ly,r)){inside++;break;}
				}else if(((a[i].xpos +2*R) > Lx)&&((a[i].ypos - 2*R) < 0)){
					if (concentricCheck(a[i],a[j],-Lx,0,r)){inside++;break;}
					if (concentricCheck(a[i],a[j],0,Ly,r)){inside++;break;}
					if (concentricCheck(a[i],a[j],-Lx,Ly,r)){inside++;break;}
				}else if (((a[i].xpos -2*R) < 0)&&((a[i].ypos - 2*R) < 0)){
					if (concentricCheck(a[i],a[j],Lx,0,r)){inside++;break;}
					if (concentricCheck(a[i],a[j],0,Ly,r)){inside++;break;}
					if (concentricCheck(a[i],a[j],Lx,Ly,r)){inside++;break;}
				}else if (((a[i].xpos -2*R) < 0)&&((a[i].ypos + 2*R) > Ly)){
					if (concentricCheck(a[i],a[j],Lx,0,r)){inside++;break;}
					if (concentricCheck(a[i],a[j],0,-Ly,r)){inside++;break;}
					if (concentricCheck(a[i],a[j],Lx,-Ly,r)){inside++;break;}
				}
			}
		}	 
	}
				
}

int concentricCheck(Disc a, Disc b, double dx, double dy, double r){
	if (((a.xpos-b.xpos)*(a.xpos-b.xpos)+(a.ypos-b.ypos)*(a.ypos-b.ypos)) < (a.r+r)(a.r+r)){return 1;}
	else {return 0;}

int collisionCheck(Disc a, Disc b, double dx, double dy, double r){
	if (((a.xpos+dx-b.xpos) < 2*(R+r)) || ((a.ypos+dy-b.ypos) < 2*(R+r))){
		if (((a.xpos+dx-b.xpos)*(a.xpos+dx-b.xpos)+(a.ypos+dy-b.ypos)*(a.ypos+dy-b.ypos)) < (a.r+b.r+2*r)*(a.r+b.r+2*r)){return 1;}
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
			if (collisionCheck(a[i],a[j],0,0,0)){i--;break;}
			//Side collision
			if (a[i].xpos + 2*R > Lx){
				if (collisionCheck(a[i],a[j],-Lx,0,0)){i--;break;}
			}else if (a[i].xpos - 2*R < 0){
				if (collisionCheck(a[i],a[j],Lx,0,0)){i--;break;}
			}else if (a[i].ypos + 2*R > Ly){
				if (collisionCheck(a[i],a[j],0,-Ly,0)){i--;break;}
			}else if (a[i].ypos - 2*R < 0){
				if (collisionCheck(a[i],a[j],0,Ly,0)){i--;break;}
			}
			//Corner collision
			if (((a[i].xpos +2*R) > Lx)&&((a[i].ypos + 2*R) > Ly)){
				if (collisionCheck(a[i],a[j],-Lx,0,0)){i--;break;}
				if (collisionCheck(a[i],a[j],0,-Ly,0)){i--;break;}
				if (collisionCheck(a[i],a[j],-Lx,-Ly,0)){i--;break;}
			}else if(((a[i].xpos +2*R) > Lx)&&((a[i].ypos - 2*R) < 0)){
				if (collisionCheck(a[i],a[j],-Lx,0,0)){i--;break;}
				if (collisionCheck(a[i],a[j],0,Ly,0)){i--;break;}
				if (collisionCheck(a[i],a[j],-Lx,Ly,0)){i--;break;}
			}else if (((a[i].xpos -2*R) < 0)&&((a[i].ypos - 2*R) < 0)){
				if (collisionCheck(a[i],a[j],Lx,0,0)){i--;break;}
				if (collisionCheck(a[i],a[j],0,Ly,0)){i--;break;}
				if (collisionCheck(a[i],a[j],Lx,Ly,0)){i--;break;}
			}else if (((a[i].xpos -2*R) < 0)&&((a[i].ypos + 2*R) > Ly)){
				if (collisionCheck(a[i],a[j],Lx,0,0)){i--;break;}
				if (collisionCheck(a[i],a[j],0,-Ly,0)){i--;break;}
				if (collisionCheck(a[i],a[j],Lx,-Ly,0)){i--;break;}
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

int jiggleDisc(Disc a[], int xsort[], int ysort[], int N, int M, double Lx, double Ly, int** map){
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
			if (collisionCheck(a[i],a[j],0,0,0)){a[xsort[x]].xpos=oldx;a[xsort[x]].ypos=oldy;return 0;}
			//Side collision
			if (a[i].xpos + 2*R > Lx){
				if (collisionCheck(a[i],a[j],-Lx,0,0)){a[xsort[x]].xpos=oldx;a[xsort[x]].ypos=oldy;return 0;}
			}else if (a[i].xpos - 2*R < 0){
				if (collisionCheck(a[i],a[j],Lx,0,0)){a[xsort[x]].xpos=oldx;a[xsort[x]].ypos=oldy;return 0;}
			}else if (a[i].ypos + 2*R > Ly){
				if (collisionCheck(a[i],a[j],0,-Ly,0)){a[xsort[x]].xpos=oldx;a[xsort[x]].ypos=oldy;return 0;}
			}else if (a[i].ypos - 2*R < 0){
				if (collisionCheck(a[i],a[j],0,Ly,0)){a[xsort[x]].xpos=oldx;a[xsort[x]].ypos=oldy;return 0;}
			}
			//Corner collision
			if (((a[i].xpos +2*R) > Lx)&&((a[i].ypos + 2*R) > Ly)){
				if (collisionCheck(a[i],a[j],-Lx,0,0)){a[xsort[x]].xpos=oldx;a[xsort[x]].ypos=oldy;return 0;}
				if (collisionCheck(a[i],a[j],0,-Ly,0)){a[xsort[x]].xpos=oldx;a[xsort[x]].ypos=oldy;return 0;}
				if (collisionCheck(a[i],a[j],-Lx,-Ly,0)){a[xsort[x]].xpos=oldx;a[xsort[x]].ypos=oldy;return 0;}
			}else if(((a[i].xpos +2*R) > Lx)&&((a[i].ypos - 2*R) < 0)){
				if (collisionCheck(a[i],a[j],-Lx,0,0)){a[xsort[x]].xpos=oldx;a[xsort[x]].ypos=oldy;return 0;}
				if (collisionCheck(a[i],a[j],0,Ly,0)){a[xsort[x]].xpos=oldx;a[xsort[x]].ypos=oldy;return 0;}
				if (collisionCheck(a[i],a[j],-Lx,Ly,0)){a[xsort[x]].xpos=oldx;a[xsort[x]].ypos=oldy;return 0;}
			}else if (((a[i].xpos -2*R) < 0)&&((a[i].ypos - 2*R) < 0)){
				if (collisionCheck(a[i],a[j],Lx,0,0)){a[xsort[x]].xpos=oldx;a[xsort[x]].ypos=oldy;return 0;}
				if (collisionCheck(a[i],a[j],0,Ly,0)){a[xsort[x]].xpos=oldx;a[xsort[x]].ypos=oldy;return 0;}
				if (collisionCheck(a[i],a[j],Lx,Ly,0)){a[xsort[x]].xpos=oldx;a[xsort[x]].ypos=oldy;return 0;}
			}else if (((a[i].xpos -2*R) < 0)&&((a[i].ypos + 2*R) > Ly)){
				if (collisionCheck(a[i],a[j],Lx,0,0)){a[xsort[x]].xpos=oldx;a[xsort[x]].ypos=oldy;return 0;}
				if (collisionCheck(a[i],a[j],0,-Ly,0)){a[xsort[x]].xpos=oldx;a[xsort[x]].ypos=oldy;return 0;}
				if (collisionCheck(a[i],a[j],Lx,-Ly,0)){a[xsort[x]].xpos=oldx;a[xsort[x]].ypos=oldy;return 0;}
			} 
		}
	}
	map[(int)a[xsort[x]].xpos][(int)a[xsort[x]].ypos]++; 
	return 1;
}



void pushDisc(Disc a[], int xsort[], int ysort[], int N, double l, double Lx, double Ly, int** map){
	double moved = 0;
	double move = 0, xdist = 0, ydist = 0, dist, mindist = Lx*Ly;
	int n = randi(N);
	int m = (n+1)%N;
	int k = m;
	double above, below, RR;
	double kxdist, kydist; 
	RR = a[xsort[n]].r + a[xsort[m]].r;
	
	while (moved <= l){
		//mindist = Lx*Ly;	
		while (xdist < 2*R+l){
			xdist = a[xsort[m]].xpos - a[xsort[n]].xpos;
			if (xdist < 0){xdist += Lx;}
			above = a[xsort[m]].ypos - a[xsort[n]].ypos;
			if (above < 0){above += Ly;}
			below = a[xsort[n]].ypos - a[xsort[m]].ypos;
			if (below < 0){below += Ly;}
			ydist = (above<below)?(above):(below);
			dist = xdist - RR*sin(acos(ydist/RR));
			if (dist <= 0){
				mindist = Lx*Ly;
				n = m;
			}else if ((mindist > dist)&&(ydist<RR)){
				mindist = dist;
				kxdist = xdist;
				kydist = ydist;				
				k = m;
			}
			m = (m+1)%N;
		}
		//mindist = Lx*Ly;
		/*if (((kxdist >= l-moved)||(kydist >= l-moved))&&(dist<RR)){
			a[xsort[n]].xpos += (l - moved);
			if (a[xsort[n]].xpos > Lx){a[xsort[n]].xpos -= Lx;}
			printf("free\n");
			bSort(xsort, ysort, a, N);
			xdist = 0;
			break;
		}*/
		move = mindist;
	//	printf("%i %i %i %4.2f %4.2f %4.2f\n", n, m, k, move, kxdist, kydist);
		move = (move>l-moved)?(l-moved):(move);
		if (move == mindist){map[(int)a[xsort[k]].xpos][(int)a[xsort[k]].ypos]++;}
	//	printf("%i %i %i %4.2f\n", n, m, k, move);
		a[xsort[n]].xpos += move;
		if (a[xsort[n]].xpos > Lx){a[xsort[n]].xpos -= Lx;}
		moved += move;
	//	printf("moved %4.2f\n",moved);
		if (moved >= l){break;}
		n = k;
		m = (k+1)%N;
		k = m;
		xdist = 0; dist = 0; mindist = Lx*Ly;
		bSort(xsort, ysort, a, N);
	}
	bSort(xsort, ysort, a, N);
	moved = 0; ydist = 0; mindist = Lx*Ly;
	n = randi(N); m = (n+1)%N; k = m;
	while (moved <= l){
		//mindist = Lx*Ly;	
		while (ydist < 2*R+l){
			ydist = a[ysort[m]].ypos - a[ysort[n]].ypos;
			if (ydist < 0){ydist += Ly;}
			above = a[ysort[m]].xpos - a[ysort[n]].xpos;
			if (above < 0){above += Lx;}
			below = a[ysort[n]].xpos - a[ysort[m]].xpos;
			if (below < 0){below += Lx;}
			xdist = (above<below)?(above):(below);
			dist = ydist - RR*sin(acos(xdist/RR));
			if (dist <= 0){
				mindist = Lx*Ly;
				n = m;
			}else if ((mindist > dist)&&(xdist<RR)){
				mindist = dist;
				k = m;
			}
			m = (m+1)%N;
		}
		move = mindist; //kxdist - sqrt(RR*RR-kydist*kydist);
		move = (move>l-moved)?(l-moved):(move);
		if (move == mindist){map[(int)a[ysort[k]].xpos][(int)a[ysort[k]].ypos]++;}
		a[ysort[n]].ypos += move;
		if (a[ysort[n]].ypos > Ly){a[ysort[n]].ypos -= Ly;}
		moved += move;
		if (moved >= l){break;}
		n = k;
		m = (k+1)%N;
		k = m;
		ydist = 0; dist = 0; mindist = Lx*Ly;
		bSort(xsort, ysort, a, N);
	}
	bSort(xsort, ysort, a, N);
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

void printMap(int ** map, int iLx, int iLy){
	FILE *mapFile;
	char fname[16];
	sprintf(fname,"map.dat");
	mapFile = fopen(fname,"w");
	for (int i = 0; i < iLx; i++){
		for (int j=0; j < iLy; j++){
			fprintf(mapFile,"%i\t",map[i][j]);
		}
		fprintf(mapFile,"\n");
	}
	fclose(mapFile);
}
