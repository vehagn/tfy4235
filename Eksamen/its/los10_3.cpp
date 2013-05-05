#include <cmath>
#include <ctime>
#include <bitset>
#include <iostream>
#include <fstream>

const double PI5 = 20.0*atan(1.0);

using std::bitset;
using std::ofstream;
using std::cout;
using std::endl;

class RNG {
private:
    int seed;
    int IBM;
    double rinv;
public:
    RNG(int seed = 12345);
    int rand();
    double ranf();
};

/*
template argument of bitset is of type "size_t",
which is a portable type of unsigned integers
see: http://www.cplusplus.com/reference/bitset/bitset/bitset/

only using "unsigned int" works in windows, but fails in linux
*/
template <size_t T>
double toDouble(bitset<T> &, double, double);
template <size_t T>
void swapbits(bitset<T> &, bitset<T> &, int, int);
double f(double);

int main(){

    // Initializing variables and parameters
    RNG * rng = new RNG(2*time(0) + 1);

    const size_t numbits = 32; // number of bits for each individual
    const int N = 100;         // number of individuals
    const int Ngen = 100;      // number of generations
    int imax = 0;
    double fmax = 0.0;
    double xmax = 0.0;
    const double pvlg = 0.75;   // best fit selection probability
    const double qvlg = 0.5;    // crossover probability
    const double rmut = 0.02;   // per-bit mutation probability
    bitset<numbits> ipop[N];    // array containing individual states
    bitset<numbits> ipopn[N];   // temporary array for selected individuals
    double fpop[N];             // array containing fitness for each individual

    // operning files, out.d contains values for all individuals for all generations
    // out2.d contains only maximum values
    ofstream outfile("out.d");
    ofstream outfile2("out2.d");

    // initializing states
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < numbits; j++) {
            ipop[i][j] = static_cast<int>(rng->ranf() + 0.5);
        }
    }

    // main loop
    for (int gen = 0; gen < Ngen; gen++){
        /* Finding maximum and write to files */
        fmax = -10000.0;
        xmax = 0.0;
        for (int i = 0; i < N; i++) {
            double x = toDouble(ipop[i],-PI5, PI5);
            double ff = f(x);
            fpop[i] = ff;
            outfile << gen << " " << x << " " << ff << endl;

            if (ff > fmax) {
                xmax = x;
                imax = i;
                fmax = ff;
            }
        }
        outfile2 << gen << " " << xmax << " " << fmax << endl;
        cout << "Current max fit: gen = " << gen << " " << fmax << ", i = " << imax << ", x = " << xmax << endl;

        /* Selection */
        for (int i = 0; i < N; i++) {
            int j = static_cast<int>( N*rng->ranf() );
            int k = static_cast<int>( N*rng->ranf() );
            double fj = fpop[j];
            double fk = fpop[k];
            int ivlg = static_cast<int>(pvlg + rng->ranf());

            if (fj > fk) {  // choose maximum of pair with probability pvlg
                ipopn[i] = (ivlg ? ipop[j] : ipop[k]);
            } else {
                ipopn[i] = (ivlg ? ipop[k] : ipop[j]);
            }
            // compact: ipopn[i] = (fj > fk ? (ivlg ? ipop[j] : ipop[k]) : (ivlg ? ipop[k] : ipop[j]) );
        }
        for (int i = 0; i < N; i++) {
            ipop[i] = ipopn[i];
        }

        /* Two-point crossover */
        for (int i = 0; i < N; i++){
            int ivlg = static_cast<int>(qvlg + rng->ranf());
            if (ivlg) {                                            // probability qvlg to be true
                int j = static_cast<int>( N*rng->ranf() );         // choose random individual
                int a = static_cast<int>( numbits*rng->ranf() );   // choose first random bit
                int b = static_cast<int>( numbits*rng->ranf() );   // choose second random bit
                swapbits(ipop[i],ipop[j],a,b);
            }
        }

        /* Mutation */
        for (int i = 0; i < N; i++) {
            bitset<numbits> mut(0);
            for (int j = 0; j < numbits; j++) {
                int imut = static_cast<int>( rmut + rng->ranf() );
                mut[j] = imut;
            }
            ipop[i] ^= mut;
        }
    }

    outfile.close();
    outfile2.close();

    delete rng;

    return 0;
}

RNG::RNG(int seed) {
    this->seed = seed;
    this->IBM = seed;
    this->rinv = 0.5/2147483647.0;
}
int RNG::rand() {
    this->IBM *= 16807;
    return (this->IBM);
}
double RNG::ranf() {
    this->IBM *= 16807;
    return (this->IBM*this->rinv + 0.5);
}

template <size_t T>
double toDouble(bitset<T> & x, double a, double b){
    static double R1 = 1.0/pow(2.0,T);  // 1.0/2^32
    double y = 0.0, powr = 1.0;
    for (unsigned int i = 0; i < x.size(); i++) {
        y += x[i]*powr;
        powr *= 2.0;
    }
    return (y*R1*(b-a) + a);
}

template <size_t T>
void swapbits(bitset<T> & x, bitset<T> & y, int a, int b){ //
    bitset<T> z(0), w(0);
    for (int i = a; i <= b; i++)
        z[i] = 1;
    w = (x^y)&z;
    x = x^w;
    y = y^w;
}

double f(double x){
    return cos(x)*exp(-x*x);
}
