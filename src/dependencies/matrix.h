#pragma once
#include <iostream>
#include <cmath>
#include <vector>

using std::cout;
using std::endl;
using std::max;

// Matrix class --->>> VA3FOD
class Matrix
{
private:
	
	static const int  MATDIM = 3; // default
	static constexpr double EPS_MATRIX = 1.e-7;
	static constexpr double pi = 3.14159265358979323846;
	
	int MatTemp = 0; // for memory management
	void AllocateMemory(const int row,const int col);
	void deAllocateMemory(void);
	bool CheckVectors(Matrix& V1, Matrix& V2, int* length, bool *rowsDim);// returns true if the 2 vectors are identical (dimensions) 
	bool CheckVectors(const Matrix& V1,const  Matrix& V2, int* length, bool* rowsDim);// returns true if the 2 vectors are identical (dimensions) 
	void CheckDimensions(const int r, const int c);

public:

	double** pd=NULL;  //  data
	
	// declare a vector type to copy the data from "pd" to it
	std::vector<double> data;
	
	// mirror copy the data from pd double pointer to the data vector for easier access for debugging.
	void mirrorData(void);

	int num_row = MATDIM;
	int num_col = MATDIM;
	int MatrixSize = num_row * num_col;
	
	Matrix();
	Matrix(const double &value);
	Matrix(const int r, const int c);
	Matrix(const int r, const int c,const double &value);
	Matrix(const double* pMat, const int row, const int col);
	Matrix(Matrix &other);  // copy constructor
	Matrix(const Matrix& other);// copy constructor
	
	~Matrix();

	// Matrix resize
	void resize(int r, int c);

	//Returns magnitude value of a vector
	const double mag(void);

	//Assigns a value to a matrix element
	void setElem(int r,int c,const double &value);

	//set all matrix elements to the user defined value
	void setMat(const double &val);

	//Returns the value at location r, c of matrix
	const double getElem(int r, int c) const;

	//Returns sequential index, counting elements row wise, of Matrix -->> a11, a12, a13, a21, etc
	const int getIndex(int r, int c) const;

	//Returns vec nx1 of Matrix col c
	void getColVec(Matrix& out, int c) const;
	Matrix& getColVec(int c) const;

	//Returns vec nx1 of from Matrix row r
	void getRowVec(Matrix& out, int r) const;
	Matrix& getRowVec(int r) const;
	
	// Returns the size of the matrix
	const int size(void) ;

	//Displays Matrix on console
	void print(int r, int c) const;
	void print(void) const;
	void print(const std::string MatName);
	
	friend std::ostream& operator<<(std::ostream& out, const Matrix& m)
	{
		out << "(" << m.num_row << "," << m.num_col << ") = \n";
		for (int i = 0; i < m.num_row; i++)
		{
			for (int j = 0; j < m.num_col; j++)
			{
				out << m.pd[i][j] << " ";
			}
			out << std::endl;
		}
		return out;
	}

	//Returns a nxn diagonal matrix from nx1 or 1xn vector
	void buildDiag(Matrix& out);

	//Returns the adjoint
	void adjoint(Matrix& out);
	
	void polar2cartesian(double mag, double az, double elev);
	Matrix& cartesian2polar(void);

	//Returns the skew-symmetric matrix from a 3-dim vector VEC
	//			| 0 -c  b|		|a|
	//			| c  0 -a| <--	|b|
	//			|-b  a  0|		|c|
	Matrix& skew_sym(void);

	//Returns the sub matrix after row r and col c have been ommitted	
	Matrix& subMatrix(int r, int c);

	Matrix& getMatrix(int rowStart, int colStart, int rowEnd, int colEnd);

	//Returns unit vector from nx1 vector
	void unitvec(Matrix& out);
	Matrix& unitvec(void);

	// Calculate the trace of a square matrix
	double trace(void);

	//Calculated the determinant
	double det(void);

	// Calculate the Rank of a Matrix
	int rank(void);

	//Returns an identity matrix of size r x c
	void identity(void);

	//Returns the inverse of a square matrix AMAT
	void inv(Matrix& out);
	Matrix & inv(void);

	//Inequality relational operator, returns true or false
	bool operator!=(Matrix &b);

	//Scalar multiplication assignment operator (scalar element by element multiplication)
	Matrix& operator*(double b);
		
	//Multiplication operator, returns matrix product 
	Matrix& operator*(Matrix& b);

	// define the multiplication operator for a "double * Matrix" expression
	friend Matrix& operator*(const double& a, Matrix& b);

	// define the subtraction operator for a "-Matrix" expression
	Matrix& operator-(void);
	
	// define the "/" operator for a matrix division by a scalar
	Matrix& operator/(const double& a);

	//Scalar multiplication assignment operator (scalar element by element multiplication)
	void operator*=(double b);

	//Multiplication assignment operator 
	void operator*=(Matrix &b);

	//Scalar Addition operator (scalar element by element addition)
	Matrix& operator+(double b);

	//Addition operator, returns matrix addition
	Matrix& operator+(Matrix &b);

	// define the addition operator for a "double + Matrix" expression
	friend Matrix& operator+(const double& a, Matrix& b);

	//Scalar addition assignment operator (scalar element by element addition)
	void operator+=(double b);	
	
	//Addition assignment operator
	void operator+=(Matrix &b);

	//Scalar substraction operator (scalar element by element substraction)
	Matrix& operator-(double b);

	//Substraction operator, returns matrix substraction
	Matrix& operator-(Matrix &b);
	
	// define the addition operator for a "double - Matrix" expression
	friend Matrix& operator-(const double& a, Matrix& b);

	//Scalar substraction assignment operator (scalar element by element substraction)
	void operator-=(double b);

	//Subtraction assignment operator
	void operator-=(Matrix &b);

	//Assignment operator 
	void operator=(Matrix &b);
	void operator=(double b);
	void operator=(const Matrix& b);

	//Equality relational operator, returns true or false
	bool operator==(Matrix &b) const;

	//Returns the scalar(dot) product of two vectors 
	double operator^(Matrix &b);

	//Alternate transposes Aij=>Aji
	Matrix& operator~();

	// define operator to access Matrix elements as object[i][j] 
	double* operator[](int i);

	// define an operator to get a submatrix from a matrix like subAA = AA[0,1][2, 2];
	Matrix& operator()(int rowStart, int colStart, int rowEnd, int colEnd);

};
