#include <stdlib.h>
#include <stdio.h>

#define MSIZE      400   /* maximum number of 0-1 variables */

int n, c, z;
int p[MSIZE][MSIZE];
int w[MSIZE];
int x[MSIZE];



void readmatrix(char *str,int dim)
{
	FILE *file;
	file=fopen(str,"r");
	printf("%d\n",file);
	for(int i=0;i<dim;i++)
	{
		for(int j=0;j<dim;j++)
		{
			//if (!fscanf(file, "%d", &p[i][j])) 
			//break;			
		}
	}
	//fclose(file);
}
void printmatrix(int dim)
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
	int val;
	c=2;
	n=8;
	
/*	for (int i=0;i<n;i++)
		for (int j=0;j<n;j++)
			
			if (i==j)
				p[i][i]=0;
			else
				if (i+j==1)
					p[i][j]=1;
				else
				{
					if (i+j==2)
						p[i][j]==2;
					else
						p[i][j]==3;
				}*/
	
/*	p[1][0]=1;
	p[0][1]=1;
	p[2][0]=2;
	p[0][2]=2;
	p[2][1]=-3;
	p[1][2]=-3;
	*/
	for (int i=0;i<n;i++)
		w[i]=1;	
	
	readmatrix("p2.txt",n);
	
	
	
	printmatrix(n);
	/*

	val=quadknap(n,c,p,w,x);
	printf("val=%d\n",val);
	printf("x=");
	for (int i=0;i<3;i++) printf("%d ",x[i]);
	putchar('\n');

	for (int i=0;i<3;i++)
	{
		for (int j=0;j<n;j++)
			printf("%d ",p[i][j]);
		printf("\n");
	}
	*/
}