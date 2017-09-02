//
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

const double THRESH	= 1E-6;
const int DIM		= 2;
const int DEP		= 3;
const double B		= 0.95;
const double r 		= 0.04;  		
const double w 		= 1.0; 	
const double p0		= 0.1;	
const double p1		= 0.9;
const double A 		= 50.0; // Artifical problem size

int M;
int MAXITER;
int iteration;

int converged(double (*parr1)[M], double (*parr2)[M], double (*pdiff)[M], const int size)
{
	double maxd1 = 0, maxd2 = 0;

	for (int i = 0; i < size; i++)
	{
		pdiff[0][i] = fabs (parr1[0][i] - parr2[0][i]);
		pdiff[1][i] = fabs (parr1[1][i] - parr2[1][i]);
		if ( pdiff[0][i] > maxd1 )
			maxd1 = pdiff[0][i];
		if ( pdiff[1][i] > maxd2 )
			maxd2 = pdiff[1][i];
	}
	return ((maxd1 < THRESH) && (maxd2 < THRESH)) ? 1 : 0; 
}

void printres(double (*parr)[M], double** ppol[], const int size)
{
	printf("J[0]\tk[0]\tkp[0]\tc[0]\tJ[1]\tk[1]\tkp[1]\tc[1]\n");
	for (int i = 0; i < size; i++)
		printf("%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n",
		 parr[0][i], ppol[0][0][i], ppol[0][1][i], ppol[0][2][i], 
		 parr[1][i], ppol[1][0][i], ppol[1][1][i], ppol[1][2][i]);
}

int main(int argc, char** argv)
{
	if (argc == 1) 
	{
		printf("Usage: \n");
		printf("%s GridSize MaxIter\n", argv[0]);
		return 0;
	}
	
	M  = atoi (argv[1]);
	MAXITER  = atoi (argv[2]);

	/////
	double** pol[DIM];
	for (int i = 0; i < DIM; i++)
	{
		pol[i]	 = (double **) malloc(DEP * sizeof(double));
		for (int j = 0; j < DEP; j++)
			pol[i][j] = (double *) malloc(M * sizeof(double));
	}

	double *k  = calloc(M, sizeof(double));
	//double *kp = calloc(M, sizeof(double));
	for (int i = 1; i < M; i++)
		k[i] = k[i-1] + A / (M-1);	// Initialize grid
	
	// Alt. malloc: double (*arr)[COLS] = malloc( ROWS * sizeof *arr)
	double (*J_new)[M] = malloc(sizeof (*J_new) * DIM);
	double (*J_old)[M] = malloc(sizeof (*J_old) * DIM);
	double (*diff)[M]  = malloc(sizeof (*diff) * DIM);

	for (int i = 0; i<DIM; i++)
	{
		for (int j=0; j<M; j++)
		{
			J_new[i][j] = J_old[i][j] = 0;
		}
	}
	for (int i = 0; i < DIM; i++)
		for (int j = 0; j < DEP; j++)
			for (int k = 0; k < M; k++)
				pol[i][j][k] = 0;
	/////
		
	for (iteration = 1; iteration <= MAXITER; iteration++)
	{
		memcpy(J_old, J_new, sizeof(*J_new) * DIM);
		for (int i = 0; i < M; i++)
		{
			double c0, c1;
			c0 		   	 = (1+r) * k[i] + w*0 - k[0];
			c1 		   	 = (1+r) * k[i] + w*1 - k[0];
			J_new[0][i]	 = sqrt(c0) + B * ( p0 * J_old[0][0] + p1*J_old[1][0] );			
			J_new[1][i]	 = sqrt(c1) + B * ( p0 * J_old[0][0] + p1*J_old[1][0] );
			for (int j = 0; j < M; j++)
			{	
				double c0, c1, t0, t1;
			 	c0 = (1+r) * k[i] + w*0 - k[j];
				t0 = sqrt(c0) + B * ( p0 * J_old[0][j] + p1 * J_old[1][j] );
				c1 = (1+r) * k[i] + w*1 - k[j];
				t1 = sqrt(c1) + B * ( p0 * J_old[0][j] + p1 * J_old[1][j] );

				if ( c0 >= 0 && t0 >= J_new[0][i] )
					J_new[0][i] = t0;
				if ( c1 >= 0 && t1 >= J_new[1][i] )
					J_new[1][i] = t1;
			}
		}
	}
	if ( converged(J_old, J_new, diff, M) ) 
		{
			printres(J_new, pol, M);
			return 0;
		}
	else
		printf("Did not converge\n");
}

