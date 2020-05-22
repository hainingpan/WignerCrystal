#include <stdio.h>
#include <stdlib.h>
#include "adder.h"

double adder(double *in1, double *in2, int len)
{
  double sum=0;
  for (int i=0;i<len;i++)
  {
	  sum+=in1[i]+in2[i];
  }
  return sum;
}