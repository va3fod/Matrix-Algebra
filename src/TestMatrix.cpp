#pragma once
#include "matrix.h"
#include <iostream>

double angle(Matrix& VEC1, Matrix& VEC2)
{
	double argument = 0;

	double scalar = VEC1 ^ VEC2; // dot product
	double abs1 = VEC1.mag();
	double abs2 = VEC2.mag();

	double dum = abs1 * abs2;
	const double EPS = 1.e-6;
	double ang = 0;
	
	// check if the dimensions of the vectors are greater then 3 and print a warning if they are
	if (VEC1.num_row > 3 || VEC2.num_col > 3)
	{
		std::cout << "Warning: The vectors dimensions greater than 3" << std::endl;
	}
	
	// check if the dimensions 2 vectors have the same numbers of rows and cols
	if (VEC1.num_row != VEC2.num_row || VEC1.num_col != VEC2.num_col)
	{
		std::cout << "Error: The vectors dimensions are not the same" << std::endl;
		return -1;
	}
	
	if (abs1 * abs2 > EPS)
	{
		argument = scalar / dum;
	}
	else
	{
		argument = 1.0;
	}

	if (argument > 1.0) argument = 1.0;

	if (argument < -1.0) argument = -1.0;

	ang = acos(argument)*180/M_PI;
	std::cout << "angle is " << ang << std::endl;
	
	return ang;
}

int main()
{
    // sample use of Matrix class
	Matrix v1(3, 1, 0);
	Matrix v2(3, 1, 0);
	double ang = 0;

	v1[0][0] = 0;
	v1[1][0] = 1;
	v1[2][0] = 2;
	
	
	
	v2[0][0] = 1;
	v2[1][0] = 3;
	v2[2][0] = 1;
	
	
		
		
	ang = angle(v1, v2);
	
	

/*
	AA[0][0] = 4;
	AA[0][1] = 1;
	AA[0][2] = 1;
    
	AA[1][0] = 1;
	AA[1][1] = 3;
	AA[1][2] = 2;
    
	AA[2][0] = 1;
	AA[2][1] = 2;
    AA[2][2] = 3;
	
	bb[0][0] = 2;
	bb[1][0] = 1;
	bb[2][0] = 1;
	bb.print("bb");
	
	brez = (AA*bb)  + (bb^bb) ;
	brez.print("brez");
	
	BB = ~AA;

	AA.print("AA");
	std::cout << "BB" << BB << std::endl;
	
	// Ax=b;
	xx = AA.inv() * bb;
	xx.print("xx = A.inv() * b");

    // print the inverse of AA
	AA.inv().print("inverse of AA");
	*/

    std::cin.get();
}

