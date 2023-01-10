#pragma once
#include "matrix.h"
#include <iostream>

int main()
{
    // sample use of Matrix class
    Matrix AA(3, 3, 1);
	Matrix BB = AA;
	Matrix bb(3, 1, 1);
	Matrix brez(3, 1, 1);
	Matrix xx;

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
	
	brez = AA*bb;
	brez.print("brez");
	
	brez = (AA*bb)  + (bb^bb) ;
	brez.print("(AA*bb)  + (bb^bb)");
	//std::cout << "brez = \n" << brez << std::endl;
	
	BB = ~AA;

	AA.print("AA");
	std::cout << "BB = \n" << BB << std::endl;
	
	// Ax=b;
	xx = AA.inv() * bb;
	xx.print("xx = A.inv() * b");

    // print the inverse of AA
	AA.inv().print("inverse of AA");
	

    std::cin.get();
}

