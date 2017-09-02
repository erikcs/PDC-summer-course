// ~ 1.9 sec
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

const double THRESH	= 1E-6;
const int DIM		= 2;
const int DEP		= 3;

const char* hmsg = "Usage: $ ./dasgn GridSize B r w p0 p1\n"
				   "Example: ./dasgn 10 0.95 0.04 1.0 0.1 0.9\n";
// Parameters
int M; 	
double B;		
double r;  		
double w; 	
double p0;	
double p1;

int converged(double (*parr1)[], double (*parr2)[], double (*pdiff)[], const int size);
void printres(double (*parr)[], double** ppol[], const int size);

int main(int argc, char** argv)
{
	if (argc < 7) 
	{
		printf("%s", hmsg);
		return 0;
	}
	
	M  = atoi (argv[1]);
	B  = atof (argv[2]);
	r  = atof (argv[3]);
	w  = atof (argv[4]);
	p0 = atof (argv[5]);
	p1 = atof (argv[6]);
	
	double** pol[DIM];
	for (int i = 0; i < DIM; i++)
	{
		pol[i]	 = (double **) malloc(DEP * sizeof(double));
		for (int j = 0; j < DEP; j++)
			pol[i][j] = (double *) malloc(M * sizeof(double));
	}

	double *k  = calloc(M, sizeof(double));
	double *kp = calloc(M, sizeof(double));
	for (int i = 1; i < M; i++)
		k[i] = kp[i] = k[i-1] + 5.0 / (M-1);	// Initialize grid
	
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

	for (;;)
	{
		memcpy(J_old, J_new, sizeof(*J_new) * DIM);
		for (int i = 0; i < M; i++)
		{
			for (int e = 0; e < 2; e++)
			{
				double c   	 = (1+r) * k[i] + w*e - kp[0];
				J_new[e][i]	 = sqrt(c) + B * ( p0 * J_old[0][0] + p1*J_old[1][0] );
				pol[e][0][i] = k[i];
				pol[e][1][i] = kp[0];
				pol[e][2][i] = c;				
			}
			for (int j = 0; j < M; j++)
			{
				for (int e = 0; e < 2; e++)
				{
					double c = (1+r) * k[i] + w*e - kp[j];
					if ( c >= 0)
					{
						double tmp = sqrt(c) + B * ( p0 * J_old[0][j] + p1 * J_old[1][j] );
						if ( tmp > J_new[e][i] )
						{
							J_new[e][i] = tmp;
							pol[e][0][i] = k[i];
							pol[e][1][i] = kp[j];
							pol[e][2][i] = c;
						}
					}
				}
			}
		}
		if ( converged(J_old, J_new, diff, M) ) 
		{
			printres(J_new, pol, M);
			return 0;
		}
	}
}

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

