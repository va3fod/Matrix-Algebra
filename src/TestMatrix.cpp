#pragma once
#include "matrix.h"
#include <iostream>

int main()
{
	Matrix AA;
	Matrix bb(3, 1, 7);
	Matrix brez;
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

	bb = AA.getColVec(1);
	
	bb[0][0] = 2;
	bb[1][0] = 1;
	bb[2][0] = 1;

	
	bb.print("bb");
	
	brez =-bb;
	brez.print("brez");
	
	// Ax=b;
	xx = AA.inv() * bb;
	xx.print("xx = A.inv() * b");

    // print AA
	AA.print("AA");

	Matrix subAA;
	subAA = AA.getMatrix(0, 1, 2,2);
	subAA.print("subAA");

	subAA = AA(0, 1, 2, 2);
	subAA.print("subAA");



    std::cin.get();
}

