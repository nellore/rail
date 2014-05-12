#include "nw.h"

#define GAP 5

void char2index(const char* str, int* indexs, int len) {
	int i;
	for(i = 0; i < len; i++) {
		switch(str[i]) {
			case 'A': case 'a': { indexs[i] = 0;   break; }
			case 'C': case 'c': { indexs[i] = 1;   break; }
			case 'G': case 'g': { indexs[i] = 2;   break; }
			case 'T': case 't': { indexs[i] = 3;   break; }
			case '-':           { indexs[i] = GAP; break; }
			default:            { indexs[i] = 4;   break; }
		}
	}
}

#define MAX3(x, y, z) MAX((x), MAX((y), (z)))

/**
 * Global alignment with a .  Note we're taking the MAX of the three contributions.
 */
void nw(int* D, int Dx, int Dy, int* mat, int mx, int my, char *xstr, int xL, char *ystr, int yL){
	int* x = (int*)malloc(sizeof(int)*(xL));
	int* y = (int*)malloc(sizeof(int)*(yL));
	char2index(xstr, x, xL); // convert string x to integer array
	char2index(ystr, y, yL); // convert string y to integer array
	D[0] = 0;
	int i,j;
	for(j = 1; j < Dy; j++) {
		D[j] = j * mat[GAP*my + y[j-1]];
	}
	for(i = 1; i < Dx; i++) {
		D[i*Dy] = i * mat[x[i-1]*my + GAP];
	}
	for(i = 1; i < Dx; i++) {
		for(j = 1; j < Dy; j++) {
			D[i*Dy+j] = (int) MAX3(
				D[(i-1) * Dy + (j-1)] + mat[x[i-1]*my + y[j-1]],
				D[(i-1) * Dy +    j ] + mat[x[i-1]*my + GAP],
				D[i     * Dy + (j-1)] + mat[GAP*my    + y[j-1]]);
		}
	}
	free(x);
	free(y);
}
