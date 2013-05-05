#include <stdio.h>

void print_integer(int in) {
  int i;
  for (i = 31; i >= 0; i--){
    printf("%d",!!((1<<i)&in));
    if (i % 8 == 0)
      printf(" ");
  }
	printf("\n");
}

void print_long_integer(long int in) {
  int i;
  for (i = 63; i >= 0; i--){
    printf("%d",!!(((long int)1<<i)&in));
    if (i % 8 == 0)
      printf(" ");
  }
	printf("\n");
}

void print_float(float * in) {
  int i;
  int * c = (int*)in;
  for (i = 31; i >= 0; i--){
    printf("%d",!!((1<<i)&(*c)));
    if (i == 31 || i == 23)
      printf(" ");
  }
	printf("\n");

}

void print_double(double * in) {
  int i;
  long int * c = (long int*)in;
  
  for (i = 63; i >= 0; i--){
    printf("%d",!!(((long int)1 << i)&(*c)));
    if (i == 63 || i == 52)
      printf(" ");
  }
	printf("\n");

}
int main(){
  int i;
  
  printf("Printing integer:\n");
  print_integer(123);
  
  printf("Printing long integer:\n");
  print_long_integer(123);

  printf("Printing float:\n");
	float f1 = 1.2;
  print_float(&f1);
 
  printf("Printing double:\n");

	double d1 = 1.2;
  print_double(&d1);
	printf("\n");

  printf("Printing float using union:\n");
  union my_union {
    int ii;
    float ff;
  } f;
  
  f.ff = 1.2;
 
  for (i = 8*sizeof(f.ff) - 1; i >= 0; i--){
    printf("%d",!!((1<<i)&(f.ii)));
    if (i == 31 || i == 23)
      printf(" ");
  }
  printf("\n");
  
  printf("Printing double using union:\n");
  union my_union2 {
    unsigned long long int ii;
    double ff;
  } g;

  g.ff = 1.2;

  for (i = 8*sizeof(g.ff) - 1; i >= 0; i--){
    printf("%d",!!(((long int)1 << i)&(g.ii)));
    if (i == 63 || i == 52)
      printf(" ");
  }
  printf("\n");

  printf("Finding max precision:\n");

  double  pres = 1.0;
  
  while (1.0 + pres > 1.0) {
    pres = pres * 0.5;
    printf("pres = %.18lf  ",pres);
    printf("1.0 + pres = %.18lf\n",pres+1.0);
  }
  pres = pres*2.0;
  printf("pres=%.18g\n",pres);
  printf("1.0 + pres=%.18g\n",pres+1.0);
  
  return 0;
}
