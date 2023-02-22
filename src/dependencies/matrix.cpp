#pragma once
#include "matrix.h"

///////////////////////////////////////////////////////////////////////////////
// Matrix constructors
Matrix::Matrix():num_row(MATDIM),num_col(MATDIM)
{	
	// allocate here memory for the matrix for a default size 3x3
	AllocateMemory(num_row, num_col);
}
Matrix::Matrix(const double &value)
{
	// allocate here memory for the matrix for a default size 3x3
	AllocateMemory(num_row, num_col);

	setMat(value);
}
Matrix::Matrix(int r, int c)
{
	CheckDimensions(r, c);

	// allocate here memory for the matrix 
	AllocateMemory(num_row, num_col);
}

Matrix::Matrix(int r, int c, const double &value)
{
	CheckDimensions(r, c);

	// allocate here memory for the matrix
	AllocateMemory(num_row, num_col);

	setMat(value);
}

Matrix::Matrix(double* pMat, int row,int col):num_row(row), num_col(col)
{
	AllocateMemory(num_row, num_col);

	for (int ii = 0; ii < row; ii++)
	{
		for (int jj = 0; jj < col; jj++)
		{
			pd[ii][jj] = pMat[jj + ii * num_col];
		}
	}
}

Matrix::Matrix(Matrix& other) // that is the default copy constructor
{
	// this is needed to avoid making a shallow copy.
	* this = other;
	//cout << "default copy constructor called" << endl;
}

Matrix::Matrix(const Matrix& other) // that is the default copy constructor
{
	// this is needed to avoid making a shallow copy.
	*this = other;
	//cout << "default copy constructor called" << endl;
}

void Matrix::CheckDimensions(int r, int c)
{
	if (r > 0)
		num_row = r;
	else
		num_row = MATDIM;

	if (c > 0)
		num_col = c;
	else
		num_col = MATDIM;
}

//Default destructor
Matrix::~Matrix()
{
	deAllocateMemory();
}

void Matrix::deAllocateMemory(void)
{
	if (pd)
	{
		for (int jj = 0; jj < num_row; jj++)
		{
			delete[] pd[jj];
			pd[jj] = NULL;
		}
		delete[] pd;
		pd = NULL;
	}

	num_row = -1;
	num_col = -1;
}

void Matrix::AllocateMemory(const int row, const int col)
{
	pd = NULL;
	size();
	
	// allocate here memory for the matrix via the double pointer
	pd = new double* [row];
	for (int ii = 0; ii < row; ii++)
	{
		pd[ii] = new double[col];
		for (int jj = 0; jj < col; jj++)
		{
			pd[ii][jj] = 0; // default data to zeros;
		}
	}
}

void Matrix::resize(int r, int c)
{
	deAllocateMemory(); // clean up

	if (r > 0)
		num_row = r;
	else
		num_row = MATDIM;

	if (c > 0)
		num_col = c;
	else
		num_col = MATDIM;

	// allocate here memory for the matrix via the double pointer
	AllocateMemory(num_row, num_col);
}

bool Matrix::CheckVectors(Matrix& V1, Matrix& V2, int* length, bool* rowsDim)
{
	if ((V1.num_row == V2.num_row) && (V1.num_col == V2.num_col) && (V1.num_row == 1 || V1.num_col == 1))
	{
		*length = num_row > num_col ? num_row : num_col;
		*rowsDim = num_row > num_col ? true : false;
		return true;
	}
	else
	{
		*length = 0;
		*rowsDim = true;
		cout <<"vectors dimentions are not the same" << endl;
		return false;
	}
}

bool Matrix::CheckVectors(const Matrix& V1, const  Matrix& V2, int* length, bool* rowsDim)
{
	if ((V1.num_row == V2.num_row) && (V1.num_col == V2.num_col) && (V1.num_row == 1 || V1.num_col == 1))
	{
		*length = num_row > num_col ? num_row : num_col;
		*rowsDim = num_row > num_col ? true : false;
		return true;
	}
	else
	{
		*length = 0;
		*rowsDim = true;
		cout << "vectors dimentions are not the same" << endl;
		return false;
	}
}

//Absolute value of vector
const double Matrix::mag(void) 
{
	double d=0.0;
	if (num_row == 1 && num_col>=1)
	{
		for (int i = 0; i < num_col; i++)
			d += pd[0][i] * pd[0][i];
	}
	
	if (num_col==1 && num_row>1)
	{
		for (int i = 0; i < num_row; i++)
			d += pd[i][0] * pd[i][0];
	}
	
	d=sqrt(d);

	if (num_row > 1 && num_col > 1)
	{
		cout << " This is not a vector. Can not calculate the magnitude" << endl;
	}

	return d;
}

//Adjoint matrix (same as det procedure however the matrix element is not multiplied into each cofactor)
void Matrix::adjoint(Matrix &out)
{	
	out.resize(num_row,num_col);

	for (int i = 0; i < num_row; i++)
	{
		for (int j = 0; j < num_col; j++)
		{
			if (((i + j) % 2) != 0)
			{
				out.pd[i][j] = (-1.0) * sub_matrix(i, j).det();
			}
			else
			{
				out.pd[i][j] = sub_matrix(i, j).det();
			}
		}
	}

	out = ~out;
}

///////////////////////////////////////////////////////////////////////////////
//Returns a nxn diagonal matrix from nx1 or 1xn vector
void Matrix::buildDiag(Matrix &out)
{
	int max_size = num_row > num_col ? num_row : num_col;
	out.resize(max_size, max_size);

	for (int i = 0; i < max_size; i++)
	{
		out.pd[i][i] = pd[i][0];
	}
}
//////////////////////////////////////////////////////////////////////////////
//Assigns a value to a matrix element
void Matrix::setElem(int r, int c, const double &value)
{
	if (r < num_row && c < num_col)
	{
		pd[r][c] = value;
	}
	else
	{
		cout << "Index out of bounds" << endl;
	}
}
///////////////////////////////////////////////////////////////////////////////
void Matrix::setMat(const double &val)
{
	//all elements values are set to "val"
	for (int i = 0; i < num_row; i++)
		for (int j = 0; j < num_col; j++)
		{
			pd[i][j] = val;
		}
}

///////////////////////////////////////////////////////////////////////////////
//Returns column vector of col c
void Matrix::getColVec(Matrix &out, int c) const
{
	out.resize(num_row, 1);

	for (int i = 0; i < num_row; i++)
	{
		out.pd[i][0] = pd[i][c];
	}
}

Matrix& Matrix::getColVec(int c) const
{
	Matrix* pMatTmp = new Matrix(num_row,1);
	//out.resize(num_row, 1);

	for (int i = 0; i < num_row; i++)
	{
		pMatTmp->pd[i][0] = pd[i][c];
	}

	pMatTmp->MatTemp = 1;

	return *pMatTmp;
}

///////////////////////////////////////////////////////////////////////////////
//Returns vector of row r
void Matrix::getRowVec(Matrix &out,int r) const
{
	out.resize(1, num_col);

	for (int i = 0; i < num_col; i++)
	{
		out.pd[0][i] = pd[r][i];
	}
}

///////////////////////////////////////////////////////////////////////////////
//Example: print on console your matrix or vector
void Matrix::print(int r,int c) const
{
	cout<<endl;

	if (r < num_row && c < num_col)
	{
		for (int i = 0; i < r; i++)
		{
			for (int j = 0; j < c; j++)
			{
				cout << pd[i][j] << "\t";
			}
			cout << endl;
		}
	}
	else
	{
		cout << "indx out of bounds" << endl;
	}
	
}	
void Matrix::print(void) const
{
	cout << endl;
	
	for (int i = 0; i < num_row; i++)
	{
		for (int j = 0; j < num_col; j++)
		{
			cout << pd[i][j] << "\t";
		}
		cout << endl;
	}
}
void Matrix::print(const std::string MatName)
{
	cout << endl;
	cout << MatName << "(" << num_row << "," << num_col << ") = ";
	print();

}
///////////////////////////////////////////////////////////////////////////////
//Returns sequential index counting elements row wise
const int Matrix::getIndex(int r, int c) const
{
	int index=0;

	if (r < num_row && c < num_col)
	{
		// r and c are C style indexes but the num_col is the actual number of columns, not the index
		index = r * num_col + c + 1;
		return index;
	}
	else
	{
		cout << "Index out of bounds" << endl;
		return -1;
	}
}

///////////////////////////////////////////////////////////////////////////////
//Returns the value at row r col c of Matrix
const double Matrix::getElem(int r, int c) const
{
	if(r<num_row && c<num_col)
		return pd[r][c];		
	else
	{
		cout << "Index out of bounds" << endl;
		return 0.0;
	}
}

///////////////////////////////////////////////////////////////////////////////
// Converts square matrix in an identity matrix
void Matrix::identity(void)
{
	if (num_row==num_col)
	{
		for (int i = 0; i < num_row; i++)
		{
			pd[i][i] = 1.0;
		}
	}
	else
	{
		cout << "Matrix is not square" << endl;
	}
}

///////////////////////////////////////////////////////////////////////////////
//Returns the inverse of a square matrix A
//Inv  A^-1 = (1/det(A)) * adj(A)
void Matrix::inverse(Matrix &out)
{
	double d=0;
	
	if (num_col!=num_row)
	{
		cout << "Matrix is not square"<<endl;
	}
	else
	{
		d=det();
		if (abs(d) > EPS_MATRIX)
		{
			d = 1 / d;
			adjoint(out);
			out *= d;
		}
		else
		{
			cout << "Determinant is zero, inverse does not exist" << endl;
		}
		
	}
}
//Returns the inverse of a square matrix A
//Inv  A^-1 = (1/det(A)) * adj(A)
Matrix & Matrix::inv(void)
{
	double d = 0;
	Matrix* pMatTmp = new Matrix(num_row, num_col);

	if (num_col != num_row)
	{
		cout << "Matrix is not square" << endl;
	}
	else
	{
		d = det();
		if (abs(d) > EPS_MATRIX)
		{
			d = 1 / d;
			adjoint(*pMatTmp);
			*pMatTmp *= d;
		}
		else
		{
			cout << "Determinant is zero, inverse does not exist" << endl;
		}
	}

	pMatTmp->MatTemp = 1;

	return *pMatTmp;
}

///////////////////////////////////////////////////////////////////////////////
//Returns the Determinant using the recursive procedure
double Matrix::det(void)
{
	double result = 0.0;

	if (num_row != num_col)
	{
		cout << "Matrix is not square." << endl;
		return 0;
	}

	//base case of a single matrix element
	if ((num_col == 1) && (num_row == 1))
		return pd[0][0];

	//second base case of a 2x2 matrix
	else if ((num_col == 2) && (num_row == 2))
		return pd[0][0] * pd[1][1] - pd[0][1] * pd[1][0];
	else
	{
		int i = 0;
		for (int j = 0; j < num_col; j++)
		{
			//use cofactors and submatricies to finish for NxN
			if (((i + j) % 2) != 0)
			{
				//(-1)^(i+j)=>i+j is odd
				result += (-1.0) * sub_matrix(i, j).det() * pd[i][j];
			}
			else
			{
				//(-1)^(i+j)=>i+j is even
				result += sub_matrix(i, j).det() * pd[i][j];
			}
		}
	}
	return result;
}

// Calculate the trace of a square Matrix
double Matrix::trace(void)
{
	double trace = 0;
	if (num_col == num_row)
	{
		for (int i = 0; i < num_row; i++)
			trace = trace + pd[i][i];
	}

	return trace;
}

// Calculate the rank of a matrix
int Matrix::rank(void)
{
	Matrix Aech(num_row,num_col,0);

	for (int i = 0; i < num_row; i++)
		for (int j = 0; j < num_col; j++)
			Aech.pd[i][j] = pd[i][j];
	
	int r_pivot = 0;

		for (int c = 0;c<num_col;c++)
		{
			int pivot_found = true;

			if (c < num_row)
			{
				r_pivot = c;
			}
			
			if (abs(Aech.pd[r_pivot][c]) < EPS_MATRIX)
			{
				pivot_found = false;
				for (int rp = r_pivot + 1; rp < num_row - 1; rp++)
				{
					if (pd[rp][c] > EPS_MATRIX && !pivot_found)
					{
					
						pivot_found = true;
						Matrix row;
						Aech.getRowVec(row, r_pivot);
					//	row.print("row to be swapped");
						for (int temp = 0; temp < num_col; temp++)
						{
							Aech.pd[r_pivot][temp] = Aech.pd[rp][temp];
							Aech.pd[rp][temp] = row.pd[0][temp];
						}
						//Aech.print("rows swapped");
					}
				}
			}

			for (int r = r_pivot+1;r<num_row;r++)
			{				
					//Make it zero formula = a(r,c)+a(r-1,c) * ( -a(r,c)/a(r-1,c) )
					for (int c_cur = 0;c_cur<num_col;c_cur++)
					{
						if (c != c_cur ) //&& abs(Aech.pd[r][c_cur])> EPS_MATRIX
						{
							Aech.pd[r][c_cur] = Aech.pd[r][c_cur] + Aech.pd[r_pivot][c_cur] * (-Aech.pd[r][c] / Aech.pd[r_pivot][c]);
						}
							
					}
					Aech.pd[r][c] = Aech.pd[r][c] + Aech.pd[r_pivot][c] * (-Aech.pd[r][c] / Aech.pd[r_pivot][c]); // one time only for the special col that contains the pivot.
				
			//	Aech.print("Aech lin");// for debug purposes
			} //end row loop
		} //end col loop

		// calculate rank by sweeping the Aech matrix here..
		int rank = 0;
		for (int r = 0; r < num_row; r++)
		{
			bool elemnonzero = false;

			for (int c = 0; c < num_col; c++)
			{
				if (abs(Aech.pd[r][c])> EPS_MATRIX)
				{
					elemnonzero = true;
					rank++;
					break;
				}
			}
			
		}
			
	return rank;
}
///////////////////////////////////////////////////////////////////////////////
//Calculates Cartesian from polar coordinates
//|V1|             | cos(elevation)*cos(azimuth)|
//|V2| = magnitude*|cos(elevation)*sin(azimuth) |
//|V3|		       |	  -sin(elevation)       |
void Matrix::polar2cartesian(double mag, double az, double elev)
{
	resize(3, 1);
	setMat(0);
	setElem(0, 0, abs(mag) * (cos(elev) * cos(az)));
	setElem(1, 0, abs(mag) * (cos(elev) * sin(az)));
	setElem(2, 0, abs(mag) * (sin(elev) * (-1.0)));
}

///////////////////////////////////////////////////////////////////////////////
//Returns polar from cartesian coordinates
// magnitude = POLAR(0,0) = |V|
// azimuth   = POLAR(1,0) = atan2(V2,V1)
// elevation = POLAR(2,0) = atan2(-V3,sqrt(V1^2+V2^2)
Matrix& Matrix::cartesian2polar(void)
{
	double d = 0.0; // magnitude
	double azimuth = 0.0;
	double elevation = 0.0;
	double denom = 0;

	int max_size = num_row > num_col ? num_row : num_col;
	Matrix* POLAR = new Matrix(max_size, 1);

	double v1 = getElem(0, 0);
	double v2 = getElem(1, 0);
	double v3 = getElem(2, 0);

	for (int i = 0; i < num_row; i++)
		d += pd[i][0] * pd[i][0];
	d = sqrt(d);

	azimuth = atan2(v2, v1);

	denom = sqrt(v1 * v1 + v2 * v2);
	if (denom > 0.)
		elevation = atan2(-v3, denom);
	else
	{
		if (v3 > 0) elevation = -pi / 2.0;
		if (v3 < 0) elevation = pi / 2.0;
		if (v3 == 0) elevation = 0.;
	}

	POLAR->setMat(0);

	POLAR->pd[0][0] = d;
	POLAR->pd[1][0] = azimuth;
	POLAR->pd[2][0] = elevation;

	POLAR->MatTemp = 1;

	return *POLAR;
}

//Inequality relational operator
bool Matrix::operator!=(Matrix &b)
{
	if ((num_col != b.num_col) || (num_row != b.num_row))
	{
		cout << "invalid dimensions for inequality operation" << endl;
		return true;
	}
		
	for (int i=0;i<num_row;i++)
		for(int j=0;j<num_col;j++)
		{
			//check to see if abs dif value differ by more than EPS
			if (abs(pd[i][j] - b.pd[i][j]) > EPS_MATRIX)
			{
				return true;
			}
		}

	return false;
}

//Scalar multiplication assignment operator
//Note: scalar must be the second operand
Matrix& Matrix::operator*(double b)
{
	Matrix* pMatTmp = new Matrix(num_row, num_col);

	for (int i = 0; i < num_row; i++)
		for (int j = 0; j < num_col; j++)
		{
			pMatTmp->pd[i][j] = pd[i][j] * b;
		}

	pMatTmp->MatTemp = 1;

	return *pMatTmp;
}

// define the multiplication operator for a "double * Matrix" expression
Matrix& operator*(const double& a, Matrix& b)
{
	Matrix* pMatTmp = new Matrix(b.num_row, b.num_col);

	for (int i = 0; i < b.num_row; i++)
		for (int j = 0; j < b.num_col; j++)
		{
			pMatTmp->pd[i][j] = b.pd[i][j] * a;
		}

	pMatTmp->MatTemp = 1;

	return *pMatTmp;
}

Matrix& Matrix::operator-(void)
{
	Matrix* pMatTmp = new Matrix(num_row, num_col);

	for (int i = 0; i < num_row; i++)
		for (int j = 0; j < num_col; j++)
		{
			pMatTmp->pd[i][j] = -pd[i][j];
		}

	pMatTmp->MatTemp = 1;

	return *pMatTmp;
}

///////////////////////////////////////////////////////////////////////////////
//Multiplication operator
Matrix& Matrix::operator*(Matrix& b)
{
	double temp=0.0;
	int dim = 0;
	bool rowsDim = 0;
	
	Matrix* pMatTmp=NULL;

	//check for proper dimensions
	if (num_col==b.num_row) 
	{
		pMatTmp = new Matrix(num_row, b.num_col);
		pMatTmp->MatTemp = 1;
		// normal matrix or vector multiplication, element by element
		
		for (int j=0;j<num_row;j++)
		{
			for (int k = 0; k < b.num_col; k++)
			{
				temp=0.0;
				for(int i=0;i<b.num_row;i++)
				{
					temp+= pd[j][i]*b.pd[i][k];
				}
				pMatTmp->pd[j][k] = temp;
			}
		}				
	}
	else if (CheckVectors(*this,b,&dim,&rowsDim) && dim==3)
	{
		pMatTmp = new Matrix(num_row,num_col);
		pMatTmp->MatTemp = 1;

		// cross product for 3x1 (1x3) vectors
		pMatTmp->pd[0][0] = pd[rowsDim ? 1 : 0][rowsDim ? 0 : 1] * b.pd[rowsDim ? 2 : 0][rowsDim ? 0 : 2] - pd[rowsDim ? 2 : 0][rowsDim ? 0 : 2] * b.pd[rowsDim ? 1 : 0][rowsDim ? 0 : 1];
		pMatTmp->pd[rowsDim ? 1 : 0][rowsDim ? 0 : 1] = pd[rowsDim ? 2 : 0][rowsDim ? 0 : 2] * b.pd[0][0] - pd[0][0] * b.pd[rowsDim ? 2 : 0][rowsDim ? 0 : 2];
		pMatTmp->pd[rowsDim ? 2 : 0][rowsDim ? 0 : 2] = pd[0][0] * b.pd[rowsDim ? 1 : 0][rowsDim ? 0 : 1] - pd[rowsDim ? 1 : 0][rowsDim ? 0 : 1] * b.pd[0][0];
	}
	else
	{
		cout << "invalid dimensions for operation * " << endl;
	}
	
	return *pMatTmp;
}

// define the "/" operator for a matrix division by a scalar
Matrix& Matrix::operator/(const double& a)
{
	Matrix* pMatTmp = new Matrix(num_row, num_col);

	for (int i = 0; i < num_row; i++)
		for (int j = 0; j < num_col; j++)
		{
			// add here a protection for division by zero (a==0)
			if (abs(a) < EPS_MATRIX)
			{
				cout << "division by zero" << endl;
				return *this;
			}
			else pMatTmp->pd[i][j] = pd[i][j] / a;
			
		}

	pMatTmp->MatTemp = 1;

	return *pMatTmp;
}

///////////////////////////////////////////////////////////////////////////////
//Scalar multiplication assignment operator (scalar element by element multiplication)
void Matrix::operator*=(double b)
{
	for (int i=0;i<num_row;i++)
		for (int j = 0; j < num_col; j++)
		{
			pd[i][j] = pd[i][j] * b;
		}	
}

///////////////////////////////////////////////////////////////////////////////
//Multiplication operator
void Matrix::operator*=(Matrix &b)
{
	//create resultant matrix
	Matrix RESULT(num_row, b.num_col);
	double temp=0.0;
	
	//check dimensions
	if (num_col!=b.num_row)
	{
		cout << "invalid dimensions for operation *= " << endl;
	}
	else
	{
		for (int k=0; k<b.num_col;k++)
		{
			for (int j=0;j<num_row;j++)
			{
				temp=0.0;
				for(int i=0;i<num_col;i++)
				{
					temp+= pd[j][i]*b.pd[i][k];
				}
				
				RESULT.pd[j][k] = temp;
			}
		}			

		deAllocateMemory();
		num_col = RESULT.num_col;
		num_row = RESULT.num_row;
		AllocateMemory(num_col, num_row);

		for (int i=0;i<num_row;i++)
			for (int j = 0; j < num_col; j++)
			{
				pd[i][j] = RESULT.pd[i][j];
			}
	}
}
///////////////////////////////////////////////////////////////////////////////
//Scalar Addition operator (scalar element by element addition)
//Note: scalar must be the second operand
Matrix& Matrix::operator+(double b)
{
	Matrix* pMatTmp = new Matrix(num_row, num_col);

	for (int i=0;i<num_row;i++)
		for (int j = 0; j < num_col; j++)
		{
			pMatTmp->pd[i][j] = pd[i][j] + b;
		}
			
	pMatTmp->MatTemp = 1;

	return *pMatTmp;
}
// define the addition operator for a "double + Matrix" expression
Matrix& operator+(const double& a, Matrix& b)
{
	Matrix* pMatTmp = new Matrix(b.num_row, b.num_col);

	for (int i = 0; i < b.num_row; i++)
		for (int j = 0; j < b.num_col; j++)
		{
			pMatTmp->pd[i][j] = b.pd[i][j] + a;
		}

	pMatTmp->MatTemp = 1;

	return *pMatTmp;
}
///////////////////////////////////////////////////////////////////////////////
//Addition operator, returns matrix addition
Matrix& Matrix::operator+(Matrix &b)
{
	Matrix* pMatTmp=nullptr;

	if (num_row == b.num_row && num_col == b.num_col)
	{
		pMatTmp = new Matrix(num_row, num_col);

		for (int i = 0; i < num_row; i++)
			for (int j = 0; j < num_col; j++)
			{
				pMatTmp->pd[i][j] = pd[i][j] + b.pd[i][j];
			}

		pMatTmp->MatTemp = 1;	
	}
	else 
	{
		cout << "Invalid Matrix dimensions in 'operator +" << endl; 
	}
	return *pMatTmp;
}
///////////////////////////////////////////////////////////////////////////////
//Scalar addition assignment operator (scalar element by element addition)
void Matrix::operator+=(double b)
{
	for (int i=0;i<num_row;i++)
		for (int j = 0; j < num_col; j++)
		{
			pd[i][j] = pd[i][j] + b;
		}	
}
///////////////////////////////////////////////////////////////////////////////
//Addition assignment operator
void Matrix::operator+=(Matrix &b)
{
	if ((num_col!=b.num_col)&&(num_row!=b.num_row))
		{cout<<"Invalid Matrix dimensions in 'operator +="<<endl;}
	else
	{
	for (int i=0;i<num_row;i++)
		for (int j = 0; j < num_col; j++)
		{
			pd[i][j] = pd[i][j] + b.pd[i][j];
		}	
	}
}
///////////////////////////////////////////////////////////////////////////////
//Scalar substraction operator (scalar element by element substraction)
//scalar must be the second operand
Matrix& Matrix::operator-(double b)
{
	Matrix* pMatTmp = new Matrix(num_row, num_col);

	for (int i=0;i<num_row;i++)
		for (int j = 0; j < num_col; j++)
		{
			pMatTmp->pd[i][j] = pd[i][j] - b;
		}
			
	pMatTmp->MatTemp = 1;

	return *pMatTmp;
}
// define the subtraction operator for a "double - Matrix" expression
Matrix& operator-(const double& a, Matrix& b)
{
	Matrix* pMatTmp = new Matrix(b.num_row, b.num_col);

	for (int i = 0; i < b.num_row; i++)
		for (int j = 0; j < b.num_col; j++)
		{
			pMatTmp->pd[i][j] = b.pd[i][j] - a;
		}

	pMatTmp->MatTemp = 1;

	return *pMatTmp;
}
///////////////////////////////////////////////////////////////////////////////
//Substraction operator, returns matrix substraction
Matrix& Matrix::operator-(Matrix &b)
{
	Matrix* pMatTmp = nullptr;

	if (num_row == b.num_row && num_col == b.num_col)
	{
		pMatTmp = new Matrix(num_row, num_col);

		for (int i = 0; i < num_row; i++)
			for (int j = 0; j < num_col; j++)
			{
				pMatTmp->pd[i][j] = pd[i][j] - b.pd[i][j];
			}

		pMatTmp->MatTemp = 1;
	}
	else
	{
		cout << "Invalid Matrix dimensions in 'operator -" << endl;
	}

	return *pMatTmp;
}

///////////////////////////////////////////////////////////////////////////////
//Scalar substraction assignment operator (scalar element by element substraction)
void Matrix::operator-=(double b)
{
	for (int i=0;i<num_row;i++)
		for (int j = 0; j < num_col; j++)
		{
			pd[i][j] = pd[i][j] - b;
		}
}

///////////////////////////////////////////////////////////////////////////////
//Subtraction assignment operator
void Matrix::operator-=(Matrix &b)
{
	if ((num_col!=b.num_col)&&(num_row!=b.num_row))
	{
		cout << "Invalid Matrix dimension in operator -=" << endl;
	}
	else
	{
		for (int i=0;i<num_row;i++)
			for (int j = 0; j < num_col; j++)
			{
				pd[i][j] = pd[i][j] - b.pd[i][j];
			}
	}
}

///////////////////////////////////////////////////////////////////////////////
void Matrix::operator=(Matrix &b)
{
	if (&b)
	{
		resize(b.num_row, b.num_col);

		for (int ii = 0; ii < num_row; ii++)
		{
			for (int jj = 0; jj < num_col; jj++)
			{
				pd[ii][jj] = b.pd[ii][jj];
			}
		}

		if (b.MatTemp)
		{
			if (b.pd)
			{
				delete& b;
			}
		}
	}
		
}
void Matrix::operator=(double b)
{
	resize(1, 1);

	pd[0][0] = b;
		
}

void Matrix::operator=(const Matrix& b)
{
	if (&b)
	{
		resize(b.num_row, b.num_col);

		for (int ii = 0; ii < num_row; ii++)
		{
			for (int jj = 0; jj < num_col; jj++)
			{
				pd[ii][jj] = b.pd[ii][jj];
			}
		}

		if (b.MatTemp)
		{
			if (b.pd)
			{
				delete& b;
			}
		}
	}
}
///////////////////////////////////////////////////////////////////////////////
//Equality relational operator of 2 identical matrixes
bool Matrix::operator==(Matrix &b)
{
	//check dimensions
	if ((num_col != b.num_col) || (num_row != b.num_row))
	{
		cout << "invalid dimensions for the equality operation" << endl;
		return false;
	}
		
	for (int i=0;i<num_row;i++)
		for(int j=0;j<num_col;j++)
		{
			//check to see if values differ by more or less than EPS
			if (abs(pd[i][j]-b.pd[i][j])> EPS_MATRIX)
				return false; 
		}

	return true;
}

//Example: value = v1 ^ v2   = v1x*v2x+v1y*v2y+v1z*v2z * ...
// Dot product of 2 vectors, of the same dimensions (1xn dot 1xn or nx1 dot nx1 or nx1 dot 1xn or 1xn dot nx1)
double Matrix::operator^(Matrix& b)
{
	double temp = 0.0;
	int dim = 0;
	bool rowsDim = 1; // 1 means the rows are > cols
	
	bool v1row = false;
	bool v2row = false;

	if (num_col == 1)
	{
		v1row = true;
	}
	
	if (b.num_col == 1)
	{
		v2row = true;
	}

	if ((num_row == 1 || num_col==1) && (b.num_row == 1 || b.num_col==1))
	{
		for (int i = 0; i < MatrixSize; i++)
		{
			temp += pd[v1row ? i : 0][v1row ? 0 : i] * b.pd[v2row ? i : 0][v2row ? 0 : i];
		}
		
		// print dot product
		cout << "Dot product: " << temp << endl;
	}
	else
	{
		cout << "invalid dimensions for the dot production operation" << endl;
	}
			
	return temp;
}

//Transpose Aij=>Aji
Matrix& Matrix::operator~()
{
	Matrix* pMatTmp = new Matrix(num_row, num_col);

	for (int i=0;i<num_row;i++)
		for (int j = 0; j < num_col; j++)
		{
			pMatTmp->pd[j][i] = pd[i][j];
		}
			
	pMatTmp->MatTemp = 1;

	return *pMatTmp;
}

double* Matrix::operator[](int i)
{
	// add checks to make sure we are not accessing index out of bounds for "pd"
	if (i < num_row)
	{
		return pd[i];
	}
	return nullptr;
}

//Returns the skew-symmetric matrix from a 3-dim vector
//			| 0 -c  b|		|a|
//			| c  0 -a| <--	|b|
//			|-b  a  0|		|c|
Matrix& Matrix::skew_sym(void)
{
	Matrix* pMatTmp = new Matrix(3, 3);

	pMatTmp->setMat(0);

	if ((num_row!=3) && (num_col!=1)) 
	{
		cout << "Invalid Matrix dimension" << endl;
	}
	else
	{
		pMatTmp->pd[0][1] = -pd[2][0];
		pMatTmp->pd[1][2] = -pd[0][0];
		pMatTmp->pd[2][0] = -pd[1][0];
		pMatTmp->pd[1][0] = pd[2][0];
		pMatTmp->pd[0][2] = pd[1][0];
		pMatTmp->pd[2][1] = pd[0][0];
	}

	pMatTmp->MatTemp = 1;

	return *pMatTmp;
}

/////////////////////////////////////////////////////////////////////////////
//Returns the sub matrix after row r and col c have been ommitted
Matrix& Matrix::sub_matrix(int r, int c)
{ 
	Matrix* pMatTmp = new Matrix(num_row - 1, num_col - 1);

	if (r < num_row && c < num_col)
	{
		pMatTmp->setMat(0);
		for (int i = 0; i < num_row; i++)
		{
			//skip the row to be removed
			if (i != r)
				if (i < r)
					
					for (int j = 0; j < num_col; j++)
					{
						if (j != c)
							if (j > c)
							{
								//adjust index for the column that has been removed	
								pMatTmp->pd[i][j - 1] = pd[i][j];
							}

							else
							{
								pMatTmp->pd[i][j] = pd[i][j];
							}
					}
				else
					//adjust indices for the row that was removed
					for (int j = 0; j < num_col; j++)
					{
						if (j != c)
							if (j > c)
							{
								//adjust index for the column that has been removed	
								pMatTmp->pd[i - 1][j - 1] = pd[i][j];
							}
							else
							{
								pMatTmp->pd[i - 1][j] = pd[i][j];
							}
					}
		}
	}
	
	pMatTmp->MatTemp = 1;
	return *pMatTmp;
}
Matrix& Matrix::getMatrix(int rowStart, int colStart, int rowEnd, int colEnd)
{
	// extract and return a submatrix from the original Matrix, with elements between start row and start column and ending to the end row and end column.
	// add protection to this code to make sure we are not accessing index out of bounds for "pd"
	if (rowStart < 0 || rowStart >= num_row || rowEnd < 0 || rowEnd >= num_row || colStart < 0 || colStart >= num_col || colEnd < 0 || colEnd >= num_col)
	{
		cout << "Invalid Matrix dimension" << endl;
		return *this;
	}
	
	Matrix* pMatTmp = new Matrix(rowEnd - rowStart + 1, colEnd - colStart + 1,0);
	for (int i = rowStart; i <= rowEnd; i++)
	{
		for (int j = colStart; j <= colEnd; j++)
		{
			pMatTmp->pd[i - rowStart][j - colStart] = pd[i][j];
		}
	}
	



	pMatTmp->MatTemp = 1;
	return *pMatTmp;
}
/////////////////////////////////////////////////////////////////////////////
//Builds unit vector from 3x1 vector
void Matrix::unitvec3(Matrix &out)
{
	double magnitude = mag();
	
	if (magnitude > EPS_MATRIX)
	{
		double first = pd[0][0];
		double second = pd[1][0];
		double third = pd[2][0];

		out.pd[0][0] = first / mag();
		out.pd[1][0] = second / mag();
		out.pd[2][0] = third / mag();
	}
	else
	{
		cout << "zero vector" << endl;
	}
}

const int Matrix::size(void) 
{
	MatrixSize = num_row * num_col;
	return num_row * num_col;
}

 
