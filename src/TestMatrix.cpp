#pragma once
#include "matrix.h"
#include <iostream>

int main()
{
	Matrix AA, ainv;
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

	// calculate the unitvector of "ff"
	Matrix ff;
	ff = dd.unitvec();
	// print "ff"
	ff.print("unit vec for ""dd"" is = ");

	double dotprodvalue = 0;
	dotprodvalue = cc ^ dd;
	cout << "dot product of 4d vectors is = " << dotprodvalue << endl;

	// add a test for the dot product of two vectors
	Matrix gg(4, 1, 1);
	gg[0][0] = 3;
	gg[1][0] = 1;
	gg[2][0] = 2;
	gg[3][0] = 4;
	// print "gg"
	gg.print("gg");

	// Additional unit tests for operators not used yet

	// Test for addition operator
	Matrix hh = cc + gg;
	hh.print("cc + gg");

	// Test for subtraction operator
	Matrix ii = cc - gg;
	ii.print("cc - gg");

	// Test for scalar multiplication operator
	Matrix jj = cc * 2;
	jj.print("cc * 2");

	// Test for scalar division operator
	Matrix kk = cc / 2;
	kk.print("cc / 2");

	// Test for matrix multiplication operator
	Matrix ll = AA * AB;
	ll.print("AA * AB");

	// insert here the code to test the rest of the operators

	// Test for element-wise multiplication operator
	Matrix nn = cc.mul(gg);
	nn.print("cc element-wise mul gg");

	// Test for negation operator
	Matrix pp = -cc;
	pp.print("negation of cc");

	// Test for equality operator
	bool isEqual = (cc == gg);
	cout << "cc == gg: " << isEqual << endl;

	// Test for inequality operator
	bool isNotEqual = (cc != gg);
	cout << "cc != gg: " << isNotEqual << endl;


	// add code below to make sure we unit test the rest of the operators and functions and methods
	
    // Test for element-wise division operator
    Matrix oo = cc;
    for (int i = 0; i < cc.num_row; ++i) {
    for (int j = 0; j < cc.num_col; ++j) {
    oo[i][j] = cc[i][j] / gg[i][j];
    }
    }
    oo.print("cc element-wise div gg");
	oo.print("cc element-wise div gg");

	// Test for negation operator (already tested as pp = -cc)

	// Test for assignment operator
	Matrix assignTest;
	assignTest = cc;
	assignTest.print("Assignment operator test (assignTest = cc)");

	// Test for scalar addition operator
	Matrix scalarAdd = cc + 5;
	scalarAdd.print("cc + 5");

	// Test for scalar subtraction operator
	Matrix scalarSub = cc - 2;
	scalarSub.print("cc - 2");

	// Test for scalar multiplication assignment operator
	Matrix scalarMulAssign = cc;
	scalarMulAssign *= 3;
	scalarMulAssign.print("cc *= 3");

	// Test for scalar addition assignment operator
	Matrix scalarAddAssign = cc;
	scalarAddAssign += 2;
	scalarAddAssign.print("cc += 2");

	// Test for scalar subtraction assignment operator
	Matrix scalarSubAssign = cc;
	scalarSubAssign -= 1;
	scalarSubAssign.print("cc -= 1");

	// Test for addition assignment operator
	Matrix addAssign = cc;
	addAssign += gg;
	addAssign.print("cc += gg");

	// Test for subtraction assignment operator
	Matrix subAssign = cc;
	subAssign -= gg;
	subAssign.print("cc -= gg");

	// Test for multiplication assignment operator
	Matrix mulAssign = cc;
	mulAssign *= gg;
	mulAssign.print("cc *= gg");

	// Test for subMatrix
	Matrix subMat = AA.subMatrix(0, 0);
	subMat.print("Submatrix of AA (removing row 0, col 0)");

	// Test for buildDiag
	Matrix diagMat;
	cc.buildDiag(diagMat);
	diagMat.print("Diagonal matrix from cc");

	// Test for adjoint
	Matrix adjMat;
	AA.adjoint(adjMat);
	adjMat.print("Adjoint of AA");

	// Test for identity
	Matrix identityMat(3, 3);
	identityMat.identity();
	identityMat.print("Identity matrix 3x3");

	// Test for polar2cartesian and cartesian2polar (example values)
	Matrix polarMat(1, 3);
	polarMat.polar2cartesian(5, 0.5, 0.3);
	Matrix cartMat = polarMat.cartesian2polar();
	cartMat.print("Cartesian to polar");

	// Test for skew_sym
	Matrix skewMat = cc.skew_sym();
	skewMat.print("Skew-symmetric matrix from cc");

	// Test for getIndex
	int idx = cc.getIndex(2, 0);
	std::cout << "Index of (2,0) in cc: " << idx << std::endl;

	// Test for getElem and setElem
	cc.setElem(0, 0, 9.0);
	double elem = cc.getElem(0, 0);
	std::cout << "Element at (0,0) in cc after setElem: " << elem << std::endl;

	// Test for mirrorData
	cc.mirrorData();
	std::cout << "Data vector after mirrorData: ";
	for (auto v : cc.data) std::cout << v << " ";
	std::cout << std::endl;







	std::cin.get();
}
