#include "nw.h"


void char2index(const char* str, int* indexs, int len){
  int i;
  for(i = 0;i<len;i++){
    switch(str[i]){
    case 'A':indexs[i]=0;break;
    case 'C':indexs[i]=1;break;
    case 'G':indexs[i]=2;break;
    case 'T':indexs[i]=3;break;
    case 'N':indexs[i]=4;break;
    case '-':indexs[i]=5;break;
    }
  }
}

#define gap 5

//Replace function call with substitution matrix A,C,G,T,N,-
//Remember to convert x and y to integer arrays
void nw(int* D, int Dx, int Dy, int* mat, int mx, int my, char *xstr, int xL, char *ystr, int yL){
  int* x = (int*)malloc(sizeof(int)*(xL));
  int* y = (int*)malloc(sizeof(int)*(yL));
  char2index(xstr,x,xL);
  char2index(ystr,y,yL);
  D[0] = 0;
  int i,j;
  for(j = 1 ; j<Dy; j++){
    D[j] = j * mat[gap*my+y[j-1]];
  }
  for(i = 1 ; i<Dx; i++){
    D[i*Dy] = i * mat[x[i-1]*my+gap];
  }
  for(i = 1; i<Dx;i++){
    for(j = 1; j<Dy;j++){
      D[i*Dy+j] = (int) MAX(MAX( (double)D[(i-1)*Dy+(j-1)] + (double)mat[x[i-1]*my+y[j-1]],
				 (double)D[(i-1)*Dy+j]   + (double)mat[x[i-1]*my+gap]),
			    (double)D[i*Dy+(j-1)] + (double)mat[gap*my+y[j-1]]);
    }
  }
  free(x);
  free(y);
}


/* int main(){ */
/*   int test1 = matchCost('-','-'); */
/*   assert (test1 == 1); */
/*   printf("Test1 succeeded\n"); */
/*   int test2 = matchCost('-','A'); */
/*   assert (test2 == -1); */
/*   printf("Test2 succeeded\n"); */
/*   /\* int* test3 = needlemanWunsch("ACGT","ACGT", matchCost); *\/ */
/*   /\* //printf("Test3 value %d\n",test3[15]); *\/ */
/*   /\* int i,j; *\/ */
/*   /\* for(i = 0 ; i < 5; i++){ *\/ */
/*   /\*   for(j = 0 ; j < 5; j++){ *\/ */
/*   /\*     printf("%d\t",test3[j+5*i]); *\/ */
/*   /\*   } *\/ */
/*   /\*   printf("\n"); *\/ */
/*   /\* } *\/ */
/*   /\* assert (test3[24] == 4); *\/ */
/*   /\* printf("Test3 succeeded\n"); *\/ */
/*   return 0; */
/* } */
