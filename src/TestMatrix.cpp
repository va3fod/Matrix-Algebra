#pragma once
#include "matrix.h"
#include <iostream>

int main()
{
    // sample use of Matrix class
    Matrix AA(3, 3, 1);

	AA[0][0] = 4;
	AA[0][1] = 1;
	AA[0][2] = 1;
    
	AA[1][0] = 1;
	AA[1][1] = 3;
	AA[1][2] = 2;
    
	AA[2][0] = 1;
	AA[2][1] = 2;
    AA[2][2] = 3;
	std::cout << "AA = " << std::endl;
	std::cout << AA << std::endl;
	

	Matrix bb(3, 1, 1);    
	Matrix brez(3, 1, 1);
	
	bb[0][0] = 2;
	bb[1][0] = 1;
	bb[2][0] = 1;

	// 
	brez = 2* bb;
	std::cout << "brez = \n" <<brez<< std::endl;
	brez = bb * 2 ;
	std::cout << "brez = \n" << brez << std::endl;
	
	
	bb.print("bb");
    
	Matrix xx(3, 1, -99); // default to -99  just for testing

	xx = conjugateGradient(AA, bb);

	//xx = AA.inv() * bb;
    
	xx.print("xx");

    // print the inverse of AA
	AA.inv().print("inverse of AA");
	

    std::cin.get();
}

