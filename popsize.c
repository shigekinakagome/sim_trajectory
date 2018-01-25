/*
 *  popsize.c
 *  
 */

#include <stdio.h>
#include <math.h>

int na = 12000 ;
int nb = 1200 ;
int nc = 12000 ;
int t1 = 500 ;
int t2 = 1100 ;
int t3 = 1600 ;
int npres = 120000 ;

	int
popsize(int t)
{
double r ;
int pop ;

	
	if( t <= t1 ) pop = npres ;
	else if (t <= t2 ) pop = na ;
	else if ( t <= t3 ) pop = nb ;
	else pop = nc ;
	return(pop);
}

