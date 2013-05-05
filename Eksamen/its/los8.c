#include <stdio.h>
#include <math.h>

int main(){
    int i;
    double PI2 = 8.0*atan(1.0);
    double omega2 = 2.0, beta = 1.0;
    double x0 = 0.0, dx = 0.001, Y0 = 1.0, u0 = 0.0; // OBS: y0 defined in math.h as a Bessel function
    double xmax = 2*(PI2/sqrt(omega2));              // 2 periods for the harmonic oscillator


    double xn = x0,yn = Y0,un = u0,tmp;
    double dd = -omega2*dx, ee;

    FILE * outfile = fopen("out_harm_explicit.d","w");
    fprintf(outfile,"#Output file for harmonic oscillator using explicit euler\n");
    fprintf(outfile,"# xn | yn | un | dx=%f omega2=%f\n",dx,omega2);
    fprintf(outfile,"%f %f %f\n",xn,yn,un);
    for (i = 0; xn < xmax; i++){
        tmp = yn + un*dx;
        un = un + yn*dd;
        yn = tmp;
        xn += dx;
        fprintf(outfile,"%f %f %f\n",xn,yn,un);
    }
    fclose(outfile);


    xn = x0, yn = Y0, un = u0,tmp;
    dd = -omega2*dx, ee = 1.0/(1.0 - dd*dx); // 1.0/(1.0 + dx^2 omega^2)

    outfile = fopen("out_harm_implicit.d","w");
    fprintf(outfile,"#Output file for harmonic oscillator using implicit euler\n");
    fprintf(outfile,"# xn | yn | un | dx=%f omega2=%f\n",dx,omega2);
    fprintf(outfile,"%f %f %f\n",xn,yn,un);
    for (i = 0; xn < xmax; i++){
        tmp = (yn + un*dx)*ee;
        un = (un + yn*dd)*ee;
        yn = tmp;
        xn += dx;
        fprintf(outfile,"%f %f %f\n",xn,yn,un);
    }
    fclose(outfile);


    xn = x0, yn = Y0, un = u0,tmp;
    dd = -omega2*dx, ee = -beta*dx;

    outfile = fopen("out_nonlin_explicit.d","w");
    fprintf(outfile,"#Output file for nonlinear oscillator using explicit euler\n");
    fprintf(outfile,"# xn | yn | un | dx=%f omega2=%f\n",dx,omega2);
    fprintf(outfile,"%f %f %f\n",xn,yn,un);
    for (i = 0; xn < xmax; i++){
        tmp = yn + un*dx;
        un = un + yn*dd + ee*yn*yn*yn;
        yn = tmp;
        xn += dx;
        fprintf(outfile,"%f %f %f\n",xn,yn,un);
    }
    fclose(outfile);
    return 0;
}
