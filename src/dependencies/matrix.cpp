#pragma once
#include "matrix.h"

///////////////////////////////////////////////////////////////////////////////
// constructors
Matrix::Matrix():num_row(MATDIM),num_col(MATDIM)
{	
	// allocate here memory for the matrix for a default size 3x3
	AllocateMemory(num_row, num_col);
	
}
Matrix::Matrix(double value)
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

Matrix::Matrix(int r, int c, double value)
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
	* this = other;
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

void Matrix::AllocateMemory(int row, int col)
{
	pd = NULL;
	
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
	deAllocateMemory();

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
		return false;
	}
}
///////////////////////////////////////////////////////////////////////////////
//Absolute value of vector. If the object instance is a matrix, specify the column index
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

///////////////////////////////////////////////////////////////////////////////
//Adjoint matrix (same as det procedure however the matrix element
//is NOT multiplied into each cofactor)
///////////////////////////////////////////////////////////////////////////////
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

///////////////////////////////////////////////////////////////////////////////
//Calculates Cartesian from polar coordinates
//|V1|             | cos(elevation)*cos(azimuth)|
//|V2| = magnitude*|cos(elevation)*sin(azimuth) |
//|V3|		       |	  -sin(elevation)       |
//
//Example: vector.cart_from_pol(dvbe,psivg,thtvg); 	
///////////////////////////////////////////////////////////////////////////////
void Matrix::polar2cartesian(double mag,double az,double elev)
{
	resize(3, 1);
	setMat(0);
	setElem(0,0, abs(mag) *(cos(elev)*cos(az)));
	setElem(1,0, abs(mag) *(cos(elev)*sin(az)));
	setElem(2,0, abs(mag) *(sin(elev)*(-1.0)));
}

///////////////////////////////////////////////////////////////////////////////
//Returns polar from cartesian coordinates
// magnitude = POLAR(0,0) = |V|
// azimuth   = POLAR(1,0) = atan2(V2,V1)
// elevation = POLAR(2,0) = atan2(-V3,sqrt(V1^2+V2^2)
//Example: POLAR = vector.pol_from_cart();
///////////////////////////////////////////////////////////////////////////////
Matrix& Matrix::cartesian2polar(void)
{
	double d = 0.0; // magnitude
	double azimuth = 0.0;
	double elevation = 0.0;
	double denom=0;

	int max_size = num_row > num_col ? num_row : num_col;
	Matrix* POLAR = new Matrix(max_size, 1);
	const double	PI = 3.1415926536;

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
		if (v3 > 0) elevation = -PI / 2.;
		if (v3 < 0) elevation = PI / 2.;
		if (v3 == 0) elevation = 0.;
	}

	POLAR->setMat(0);

	POLAR->pd[0][0] = d;
	POLAR->pd[1][0] = azimuth;
	POLAR->pd[2][0] = elevation;

	POLAR->MatTemp = 1;

	return *POLAR;
}

//////////////////////////////////////////////////////////////////////////////
//Assigns a value to a matrix element
void Matrix::setElem(int r, int c, double value)
{
	if (r < num_row && c < num_col)
	{
		pd[r][c] = value;
	}
	else
	{
		cout << "index out of bounds" << endl;
	}
}
///////////////////////////////////////////////////////////////////////////////
void Matrix::setMat(double val)
{
	//all values are set to "val"
	for (int i = 0; i < num_row; i++)
		for (int j = 0; j < num_col; j++)
		{
			pd[i][j] = val;
		}
}

///////////////////////////////////////////////////////////////////////////////
//Returns column vector of col c
///////////////////////////////////////////////////////////////////////////////
void Matrix::getCol(Matrix &out, int c) const
{
	out.resize(num_row, 1);

	for (int i = 0; i < num_row; i++)
	{
		out.pd[i][0] = pd[i][c];
	}
}

///////////////////////////////////////////////////////////////////////////////
//Returns vector of row r
//Example: vector = MAT.row_vector(2);
///////////////////////////////////////////////////////////////////////////////
void Matrix::getRow(Matrix &out,int r) const
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
	cout << "\""<<MatName <<"\""<< " is :";
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
		cout << "indx out of bounds" << endl;
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
		cout << "indx out of bounds" << endl;
		return 0.0;
	}
}

///////////////////////////////////////////////////////////////////////////////
// Converts sqaure matrix in an identity matrix
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
		cout << "matrix is not square" << endl;
	}
}

///////////////////////////////////////////////////////////////////////////////
//Returns the inverse of a square matrix AMAT
//Inversion  A^-1=(1/det(A))*Adj(A)
//Example: BMAT = AMAT.inverse();
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
		d=1/d;
		adjoint(out);
		out*=d;
	}
}

///////////////////////////////////////////////////////////////////////////////
//Returns the Determinant using the recursive procedure
double Matrix::det(void)
{
	double result = 0.0;

	if (num_row != num_col)
	{
		cout << "matrix is not square." << endl;
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
						Aech.getRow(row, r_pivot);
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
//Inequality relational operator
bool Matrix::operator!=(Matrix &b)
{
	if ((num_col != b.num_col) || (num_row != b.num_row))
	{
		return true;
	}
		
	for (int i=0;i<num_row;i++)
		for(int j=0;j<num_col;j++)
		{
			//check to see if values differ by more or less than EPS
			if (abs(pd[i][j] - b.pd[i][j]) > EPS_MATRIX)
			{
				return true;
			}
		}

	return false;
}

///////////////////////////////////////////////////////////////////////////////
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
		cout << "invalid dimensions for operation" << endl;
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
//Example: CMAT = AMAT + b; 
///////////////////////////////////////////////////////////////////////////////
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

///////////////////////////////////////////////////////////////////////////////
//Addition operator, returns matrix addition
///////////////////////////////////////////////////////////////////////////////
Matrix& Matrix::operator+(Matrix &b)
{
	Matrix* pMatTmp = new Matrix(num_row, num_col);

	for (int i = 0; i < num_row; i++)
		for (int j = 0; j < num_col; j++)
		{
			pMatTmp->pd[i][j] = pd[i][j] + b.pd[i][j];
		}

	pMatTmp->MatTemp = 1;

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
void Matrix::operator+=(Matrix &b)// Modified here
{
	if ((num_col!=b.num_col)&&(num_row!=b.num_row))
		{cout<<"Invalid Matrix dimension in 'operator +="<<endl;}
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
//Note: scalar must be the second operand
//Example: CMAT = AMAT - b;
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

///////////////////////////////////////////////////////////////////////////////
//Substraction operator, returns matrix substraction
Matrix& Matrix::operator-(Matrix &b)
{
	Matrix* pMatTmp = new Matrix(num_row, num_col);

	if ((num_col!=b.num_col)&&(num_row!=b.num_row))
	{
		cout<<"Invalid Matrix dimension in operator-"<<endl;
	}
	else
	{
		for (int i=0;i<num_row;i++)
			for (int j = 0; j < num_col; j++)
			{
				pMatTmp->pd[i][j] = pd[i][j] - b.pd[i][j];
			}
	}

	pMatTmp->MatTemp = 1;

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
				delete &b;
			}
		}
}

///////////////////////////////////////////////////////////////////////////////
//Equality relational operator
bool Matrix::operator==(Matrix &b)
{
	//check dimensions
	if ((num_col!=b.num_col) || (num_row != b.num_row))
		return false;
	
	for (int i=0;i<num_row;i++)
		for(int j=0;j<num_col;j++)
		{
			//check to see if values differ by more or less than EPS
			if (abs(pd[i][j]-b.pd[i][j])> EPS_MATRIX)
				return false; 
		}

	return true;
}

///////////////////////////////////////////////////////////////////////////////
//Example: value = AMAT ^ BMAT;   = ax*bx+ay*by+az*bz
// Dot product of 2 vectors
double Matrix::operator^(Matrix &b)
{
	double temp=0.0;
	int dim = 0;
	bool rowsDim = 1; // 1 means the rows are > cols

	//check dimensions
	if (CheckVectors(*this,b,&dim, &rowsDim))
	{
		for (int i = 0; i < dim; i++)
		{
				temp += pd[rowsDim?i:0][rowsDim ? 0 : i] * b.pd[rowsDim ? i : 0][rowsDim ? 0 : i];
		}
	}
	else
	{
		cout << "invalid dimenstion for the operation requested";
	}
			
	return temp;
}

///////////////////////////////////////////////////////////////////////////////
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

///////////////////////////////////////////////////////////////////////////////
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

///////////////////////////////////////////////////////////////////////////////
//Returns the sub matrix after row r and col c have been ommitted
//Example: BMAT = AMAT.sub_matrix(1,3); 
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
					// DO NOT adjust indices NO row has been removed
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

///////////////////////////////////////////////////////////////////////////////
//Builds unit vector from 3x1 vector
///////////////////////////////////////////////////////////////////////////////
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

