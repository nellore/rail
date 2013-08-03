#include <math.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#define MIN(X,Y) ({(X) < (Y) ? (X) : (Y);})
#define MAX(X,Y) ({(X) > (Y) ? (X) : (Y);})



void nw(int* D, int Dx, int Dy, 
	int* mat, int mx, int my, 
	char *xstr, int xL,
	char *ystr, int yL);

