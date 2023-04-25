#pragma once
#include "matrix.h"
#include <iostream>

int main()
{
	Matrix AA,ainv;
	Matrix bb(3, 1, 7);
	Matrix aa(1, 3, 0);
	Matrix brez;
	Matrix xx;

	Matrix AB(AA);

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
	aa = AA.getRowVec(1);
	//print aa
	aa.print("aa");
	
	bb[0][0] = 2;
	bb[1][0] = 1;
	bb[2][0] = 1;

	
	bb.print("bb");
	
	brez =-bb;
	brez.print("brez");
	
	// Ax=b;
	xx = AA.inv() * bb;
	xx.print("xx = A.inv() * b");
	AA.inv(ainv);
	ainv = AA.inv();

    // print AA
	AA.print("AA");

	Matrix subAA;
	subAA = AA.getMatrix(0, 1, 2,2);
	subAA.print("subAA");

	subAA = AA(0, 1, 2, 2);
	subAA.print("subAA");

	// declare a column vector of 4 elements and initialize it to 1
	Matrix cc(4, 1, 1);
	// populate "cc" with the random values between 0 and 10
	cc[0][0] = 3;
	cc[1][0] = 1;
	cc[2][0] = 2;
	cc[3][0] = 4;
	// print "cc"
	cc.print("cc");

	// calculate the unitvector of "cc"
	Matrix ee;
	//cc.unitvec(ee);
	ee = cc.unitvec();
	// print "ee"
	ee.print("unit vec for ""cc"" is = ");


	// declare a row vector of 4 elements and initialize it to 1
	Matrix dd(1, 4, 1);
	// populate "cc" with the random values between 0 and 10
	dd[0][0] = 3;
	dd[0][1] = 1;
	dd[0][2] = 2;
	dd[0][3] = 4;
	// print "dd"
	dd.print("dd");

	// calculate the unitvector of "dd"
	Matrix ff;
	//dd.unitvec(ff);
	ff = dd.unitvec();
	// print "ff"
	ff.print("unit vec for ""dd"" is = ");




    std::cin.get();
}

