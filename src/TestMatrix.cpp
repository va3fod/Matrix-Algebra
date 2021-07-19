#pragma once
#include "matrix.h"
#include <iostream>

// testing the Matrix algebra class... VA3FOD
//*********************************************
int main()
{
    // sample use of Matrix class
  
    /*Matrix V1(3, 1, 3);
    Matrix V2(3, 1,99);

    V1.unitvec3(V2);

    V1.setElem(2, 0, 99);
    V2.setMat(2);
    double m = V2.mag();

    V1.print("V1");
    V2.print("V2");
    */

  /*  Matrix AA(3, 3, 1), BB(2, 3, 2);
    
    AA.setElem(0, 0, -3);
    AA.setElem(0, 1, 2);
    AA.setElem(0, 2, -5);

    AA.setElem(1, 0, -1);
    AA.setElem(1, 1, 0);
    AA.setElem(1, 2, -2);

    AA.setElem(2, 0, 3);
    AA.setElem(2, 1, -4);
    AA.setElem(2, 2, 1);
    AA.print("A");

    double det = AA.det();
  
    AA.adjoint(BB);
    BB.print("BB");

    Matrix vec1;
    vec1.print("vec1");
    AA.getCol(vec1, 1);
    vec1.print("vec1");
    AA.getRow(vec1, 1);
    vec1.print("vec1");

    Matrix vec;
    vec.print("vec");
    vec.polar2cartesian(1, 0.3, 0.5);
    vec.print("vec");

    Matrix back(1,2);
    back = vec.cartesian2polar();
    back.print("back");

    Matrix vecc(3, 1, 1);
    vecc.pd[2][0] = 10;

    Matrix skeww(3, 4, 1);
    skeww = vecc.skew_sym();
    skeww.print("skeww");

    Matrix invMat(2, 3);
    AA.inverse(invMat);
    invMat.print("invMat");*/


    Matrix A(4, 4, 2);

    A.print("A is = ");
    int rankvalue = A.rank();
    cout << "rank is = "<< rankvalue << endl;

    
    Matrix B=A;
    
    A += A;
    A.print("A = ");
    B.print("B = ");
    bool eq = A == B ? 1 : 0;

    A = B * A + 2 - B*3;
    A.print("A = ");

 /*   Matrix C(2, 3, 2), D(3, 4, 2);
    Matrix Prod;

    C.print("C = ");
    D.print("D = ");

    Prod = C * D;
    Prod.print("Prd matrix = ");*/

   /* Matrix v1(3, 1, 2);
    Matrix v2(3, 1, -3);

    Matrix v3;
    v3 = v1 * v2;
    v3.print("v3=");*/
    
    Matrix v1(4, 1, 2),v2(1,4,3);
    double dotProd = v1 ^ v2;
    cout << "dot prod = " << dotProd << endl;

    std::cin.get();
}

//*********************************************************
