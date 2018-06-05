#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define N 3
#define err 0.05

double f1(double x[N]) {return (15*x[1]+pow(x[1], 2)-4*x[2]-13);}

double f2(double x[N]) {return (pow(x[0], 2)+10*x[1]-x[2]-11);}

double f3(double x[N]) {return (pow(x[1], 3)-25*x[2]+22);}

double g(double x[N]) {return (pow(f1(x), 2) + pow(f2(x), 2) + pow(f3(x), 2));}

double df(double f(), double x[N], int i)
{
	double xant, dy, dx;
	
	xant = x[i];
	x[i] +=err;
	dy=f(x);
	x[i]-=2*err;
	dx=dy-f(x);
	x[i]=xant;
	return (dx/(2*err));
}

double* gradiente(double **M, double F[N])
{
	int i, j, k;
	double *grad;
	
	
	grad = malloc(N*sizeof(double));
	
	for(i=0; i<N; i++)
	{
		for(j=0; j<N; j++)
			grad[i] += -2*M[j][i]*F[j];
	}
	
	return grad;
}

double h(double *grad, double x[N], double a)
{
	int i;
	double xo[N]={0, 0, 0}, t;
	
	for(i=0; i<N; i++)
		xo[i]=x[i]-a*grad[i];
	
	t=g(xo);
	
	return t;
}

double alpha(double h1, double h3, double a2) {return ((h3*a2-h1)/(2*h3));}

void imprime(double **M, int var){
	
	int i, j;
	
	for(i=0;i<var;i++) {
		for(j=0;j<var;j++) { 
			printf("%5.2lf\t",M[i][j]);
		}	  
		puts("");
	}
}

int main()
{
	FILE *fp;
	double **J;
	double xa[N] ={0.1, 0.1, -0.1}, F[N];
	double *grad;
	double eps=1e-5, eps2, tol, norm, norma, a1, a2, a3, h1, h2, h3, amed, g1, g2, g3, normgrad=0;
	double (*equacao[N])() = {f1, f2, f3}; //vetor de funções
	int i, j, cont=0;
	
	eps2 = eps*eps;
	
	J = malloc( N*sizeof(double *));
	for( i = 0 ; i < N ; i++ )
		J[i] = (double *) malloc(N*sizeof(double));
		
	
	do
	{
		norma=norm=normgrad=0;
		
		for(i=0; i<N; i++)
		{
			F[i] = (-1)* equacao[i](xa);
			norma += xa[i]*xa[i];
			
			for(j=0; j<N; j++)
				J[i][j]=df(equacao[i], xa, j);
		}
		
		
		grad=gradiente(J, F);
		
		
		
		for(i=0; i<N; i++)
			normgrad+=grad[i]*grad[i];
		
		normgrad=sqrt(normgrad);
			
		for(i=0; i<N; i++)
			grad[i]/=normgrad;
		
		
		a1=0;
		a3=1;
		
		
		while(h(grad, xa, a3)>h(grad, xa, a1))
			a3=a3/2;
			
		if((a3-a1)<eps2)
		{
			a3=a1;
			a1=1;
			a2=a1/2;
		}
		
		else
			a2=a3/2;
			
			
		g1=h(grad, xa, a1);
		
		g3=h(grad, xa, a3);
			
		g2=h(grad, xa, a2);
		
		h1=(g2-g1)/(a2-a1);
		h2=(g3-g2)/(a3-a2);
		h3=(h2-h1)/(a3-a1);
		
		
		amed=alpha(h1, h3, a2);
		
		
		
		for(i=0; i<N; i++)
		{
			xa[i]-=amed*grad[i];
			norm+=xa[i]*xa[i];
		}
		
		tol= (norm-norma)/norm;
		g1=g2=g3=0;
		cont++;
		
	}while(fabs(tol)>eps2);
	
	printf("\n\nSolução\n\n");
	
	for(i=0; i<N; i++)
		printf("%lf\t", xa[i]);
	
	puts("\n\n");
	
	return 0;
}
		
