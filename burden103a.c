#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define N 2
#define h 1e-6

double f1(double x[N]) {return (4*pow(x[0], 2) - 20*x[0] + 0.25*pow(x[1], 2) +8);}

double f2(double x[N]) {return (0.5*x[0]*pow(x[1], 2) + 2*x[0] - 5*x[1] + 8);}

double df(double f(), double x[N], int i)
{
	double xant, dy, dx;
	
	xant = x[i];
	x[i] +=h;
	dy=f(x);
	x[i]-=2*h;
	dx=dy-f(x);
	x[i]=xant;
	return (dx/(2*h));
}

void imprime(double **M, int var){
	
	int i, j;
	
	for(i=0;i<var;i++) {
		for(j=0;j<var+1;j++) { 
			printf("%5.2lf\t",M[i][j]);
		}	  
		puts("");
	}
}

void pivoteamento(double **M, int dim)
{
	double aux1, lamb;
	int i, j, k, aux2, o;
	
	
	for(i=0; i<dim; i++)
	{
		aux2=i;
		
		for(j=i+1; j<dim; j++)
		{
			if(fabs(M[aux2][aux2])<fabs(M[j][aux2]))
				aux2 = j;
		}
		
		
		o=i;
		
		if(aux2>i)
		{
			for(k=o; k<dim+1; k++)
			{
				aux1=M[aux2][k];
				M[aux2][k] = M[o][k];
				M[o][k] = aux1;
			}
			
		}
		
			
		for(j=i; j<dim-1; j++)
		{
			lamb=M[j+1][i]/M[i][i];
			for(k=i; k<dim+1; k++)
				M[j+1][k]-=(lamb) * M[i][k];
		}
	}
}

void subsreversa (double **M, double *raizes, int dim)
{
	int i, j;
	double soma=0;
	
	
	for(i=dim-1; i>=0; i--)
	{
		raizes[i] = M[i][dim];
		
		for(j=dim; j>=i+1; j--)
			raizes[i]-=M[i][j]*raizes[j];
			
		raizes[i]/=M[i][i];
	}
}

int main(int argc, char **argv)
{
	double **S;
	double **J;
	double xa[N] ={0, 0}, F[N];
	double *raizes;
	double eps=1e-5, eps2, tol, norm, norma;
	double (*equacao[N])() = {f1, f2}; //vetor de funções
	int i, j;
	
	eps2 = eps*eps;
	raizes = malloc(N*sizeof(double ));
	S = malloc( N*sizeof(double *));
	for( i = 0 ; i < N ; i++ )
		S[i] = (double *) malloc((N+1)*sizeof(double));
	
	J = malloc( N*sizeof(double *));
	for( i = 0 ; i < N ; i++ )
		J[i] = (double *) malloc((N+1)*sizeof(double));
	
	do
	{
		

		norma=norm=0;
		
		for(i=0; i<N; i++)
		{
			F[i] = (-1)* equacao[i](xa);
			norma += xa[i]*xa[i]; 
		
			for(j=0; j<N; j++)
				J[i][j]=df(equacao[i], xa, j);
		}
		
		for(i=0; i<N; i++)
		{
			for(j=0; j<N; j++)
				S[i][j]=J[i][j];
		}
		
		for(i=0; i<N; i++)
			S[i][N] = F[i];
			
		imprime(S, N);
		
		pivoteamento(S, N);
		subsreversa(S, raizes, N);
		
		
		
		for(i=0; i<N; i++)
		{
			xa[i]+=raizes[i];
			norm += xa[i]*xa[i];
		}
		
		puts("\n");
		
		tol = (norm-norma)/norm;
		
	}while(fabs(tol)>eps2);
	
	printf("\n\nSolução\n\n");
	
	for(i=0; i<N; i++)
		printf("%lf\t", xa[i]);
	
	puts("\n");
	
	return 0;
}
