#include <stdlib.h>
#include <stdio.h>

#define MSIZE      400   /* maximum number of 0-1 variables */
int n, c, z;
int p[MSIZE][MSIZE];
int w[MSIZE];
int x[MSIZE];


void readmatrix(char *str,int p[MSIZE][MSIZE] ,int dim)
{
	FILE *fp;
	fp=fopen(str,"r");
	for(int i=0;i<dim;i++)
	{
		for(int j=0;j<dim;j++)
		{
			if (!fscanf(fp, "%d", &p[i][j]))
			break;
		}
	}
	fclose(fp);
}
void printmatrix(int p[MSIZE][MSIZE] ,int dim)
{
	for(int i=0;i<dim;i++)
	{
		for(int j=0;j<dim;j++)
		{
			printf("%d\t",p[i][j]);
		}
		printf("\n");
	}
}

void main()
{
    int den,no;
    printf("n: denominator: numerator:\n");
    scanf("%d %d %d",&n,&den,&no);
	//n=98;
	//c=2*n/7;
    c=no*n/den;
	int val;
	readmatrix("p2.txt",p,n);
	for (int i=0;i<n;i++)
		w[i]=1;

	//printmatrix(p,n);
	val=quadknap(n,c,p,w,x);

	printf("val=%d\n",val);
	for (int i=0;i<n;i++) printf("%d ",x[i]);
}
