#pragma once
#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>

using std::cout;
using std::endl;

// Matrix class --->>> VA3FOD
class Matrix
{
private:
	
	const int  MATDIM = 3; // default
	const double EPS_MATRIX = 1.e-7;
	
	int MatTemp = 0; // for memory management
	void AllocateMemory(int row, int col);
	void deAllocateMemory(void);
	bool CheckVectors(Matrix& V1, Matrix& V2, int* length, bool *rowsDim);// returns true if the 2 vectors are identical (dimensions) 
	void CheckDimensions(int r, int c);

	int num_row = MATDIM; 
	int num_col = MATDIM; 

public:

	double** pd=NULL;  //  data
	
	Matrix();
	Matrix(const double &value);
	Matrix(int r, int c);
	Matrix(int r, int c,const double &value);
	Matrix(double* pMat, int row, int col);
	Matrix(Matrix &other); 
	
	~Matrix();

	// Matrix resize
	void resize(int r, int c);

	//Returns magnitude value of a vector
	const double mag(void);

	//Assigns a value to a matrix element
	void setElem(int r,int c,double value);

	//set all matrix elements to the user defined value
	void setMat(double val);

	//Returns the value at location r, c of matrix
	const double getElem(int r, int c) const;

	//Returns sequential index, counting elements row wise, of Matrix -->> a11, a12, a13, a21, etc
	const int getIndex(int r, int c) const;

	//Returns vec 3x1 of Matrix col c
	void getCol(Matrix& out, int c) const;

	//Returns vec 3x1 of from MAtrix row r
	void getRow(Matrix& out, int r) const;

	//Displays Matrix on console
	void print(int r, int c) const;
	void print(void) const;
	void print(const std::string MatName);

	//Returns a nxn diagonal matrix from nx1 or 1xn vector
	void buildDiag(Matrix& out);

	//Returns the adjoint
	void adjoint(Matrix& out);
	
	// Polar coordinates to cartesian transformation
	// |V1|             | cos(elevation)*cos(azimuth)|
	// |V2| = magnitude*|cos(elevation)*sin(azimuth) |
	// |V3|		        |	  -sin(elevation)        |
	void polar2cartesian(double mag,double az,double elev); // azimuth and elevation angles given in radians  , // to be moved to a sim util library type of space... as a function
	
	// Cartesian coordinates to polar transformation
	// magnitude = POLAR(0,0) = |V|
	// azimuth   = POLAR(1,0) = atan2(V2,V1)
	// elevation = POLAR(2,0) = atan2(-V3,sqrt(V1^2+V2^2)
	Matrix& cartesian2polar(void); // to be moved to a sim util library type of space... as a function

	//Returns the skew-symmetric matrix from a 3-dim vector VEC
	//			| 0 -c  b|		|a|
	//			| c  0 -a| <--	|b|
	//			|-b  a  0|		|c|
	Matrix& skew_sym(void);

	//Returns the sub matrix after row r and col c have been ommitted	
	Matrix& sub_matrix(int r, int c);

	//Returns unit vector from 3x1 vector
	void unitvec3(Matrix& out);

	// Calculate the trace of a square matrix
	double trace(void);

	//Calculated the determinant
	double det(void);

	// Calculate the Rank of a Matrix
	int rank(void);

	//Returns an identity matrix of size r x c
	void identity(void);

	//Returns the inverse of a square matrix AMAT
	void inverse(Matrix& out);


	//********************************************************* Operator declarations *********************************************
	//Inequality relational operator, returns true or false
	bool operator!=(Matrix &b);

	//Scalar multiplication assignment operator (scalar element by element multiplication)
	Matrix& operator*(double b);
		
	//Multiplication operator, returns matrix product 
	Matrix& operator*(Matrix& b);

	//Scalar multiplication assignment operator (scalar element by element multiplication)
	void operator*=(double b);

	//Multiplication assignment operator 
	void operator*=(Matrix &b);

	//Scalar Addition operator (scalar element by element addition)
	Matrix& operator+(double b);

	//Addition operator, returns matrix addition
	Matrix& operator+(Matrix &b);

	//Scalar addition assignment operator (scalar element by element addition)
	void operator+=(double b);	
	
	//Addition assignment operator
	void operator+=(Matrix &b);

	//Scalar substraction operator (scalar element by element substraction)
	Matrix& operator-(double b);

	//Substraction operator, returns matrix substraction
	Matrix& operator-(Matrix &b);

	//Scalar substraction assignment operator (scalar element by element substraction)
	void operator-=(double b);

	//Subtraction assignment operator
	void operator-=(Matrix &b);

	//Assignment operator 
	void operator=(Matrix &b);

	//Equality relational operator, returns true or false
	bool operator==(Matrix &b);

	//Returns the scalar(dot) product of two vectors 
	double operator^(Matrix &b);

	//Alternate transposes Aij=>Aji
	Matrix& operator~();
};

