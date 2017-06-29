/***********************************************************************************/
/*  Copyright 2011 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
/***********************************************************************************/
/* This file is part of MeteoIO.
    MeteoIO is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MeteoIO is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with MeteoIO.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <meteoio/dataClasses/Matrix.h>
#include <meteoio/IOUtils.h>
#include <meteoio/IOExceptions.h>

#include <time.h> //needed for random()
#include <cmath> //needed for fabs()
#include <iostream>
#include <iomanip>

namespace mio {

const double Matrix::epsilon = 1e-9; //for considering a determinant to be zero, etc
const double Matrix::epsilon_mtr = 1e-6; //for comparing two matrix

Matrix::Matrix(const int& rows, const int& cols) : vecData(), ncols(0), nrows(0) {
	if (rows<0 || cols<0) {
		std::ostringstream tmp;
		tmp << "Trying construct a matrix with negative dimensions: ";
		tmp << "(" << rows << "," << cols << ")";
		throw IOException(tmp.str(), AT);
	}
	resize((unsigned)rows,(unsigned)cols);
}

Matrix::Matrix(const size_t& n, const double& init) : vecData(n*n, 0.), ncols(n), nrows(n) {
	for (size_t ii=1; ii<=n; ii++) operator()(ii,ii) = init;
}

void Matrix::identity(const size_t& n, const double& init) {
	resize(n,n,0.);
	for (size_t ii=1; ii<=n; ii++) operator()(ii,ii) = init;
}

void Matrix::resize(const size_t& rows, const size_t& cols) {
	clear();
	vecData.resize(rows*cols);
	nrows = rows;
	ncols = cols;
}

void Matrix::resize(const size_t& rows, const size_t& cols, const double& init) {
	clear();
	vecData.resize(rows*cols, init);
	ncols = cols;
	nrows = rows;
}

void Matrix::size(size_t& rows, size_t& cols) const{
	rows=nrows;
	cols=ncols;
}

void Matrix::clear() {
	vecData.clear();
	nrows=ncols=0;
}

void Matrix::random(const double& range) {
	srand( static_cast<unsigned>(time(0)) );

	for (size_t ii=0; ii<vecData.size(); ii++)
		vecData[ii] = (double)rand()/(double)RAND_MAX*range;
}

double& Matrix::operator ()(const size_t& i, const size_t& j) {
#ifndef NOSAFECHECKS
	if ((i<1) || (i > nrows) || (j<1) || (j > ncols)) {
		std::ostringstream ss;
		ss << "Trying to access matrix[" << i << "," << j << "]";
		throw IndexOutOfBoundsException(ss.str(), AT);
	}
#endif
	return vecData[(j-1) + (i-1)*ncols];
}

double Matrix::operator ()(const size_t& i, const size_t& j) const {
#ifndef NOSAFECHECKS
	if ((i<1) || (i > nrows) || (j<1) || (j > ncols)) {
		std::ostringstream ss;
		ss << "Trying to access matrix[" << i << "," << j << "]";
		throw IndexOutOfBoundsException(ss.str(), AT);
	}
#endif
	return vecData[(j-1) + (i-1)*ncols];
}

const std::string Matrix::toString() const {
	std::ostringstream os;
	const size_t wd=6;
	os << "\n┌ ";
	for (size_t jj=1; jj<=(ncols*(wd+1)); jj++)
		os << " ";
	os << " ┐\n";
	for (size_t ii=1; ii<=nrows; ii++) {
		os << "│ ";
		for (size_t jj=1; jj<=ncols; jj++) {
			os << std::setw(wd) << std::fixed << std::setprecision(2) << operator()(ii,jj) << " ";
		}
		os << " │\n";
	}
	os << "└ ";
	for (size_t jj=1; jj<=(ncols*(wd+1)); jj++)
		os << " ";
	os << " ┘\n";
	return os.str();
}

bool Matrix::operator==(const Matrix& in) const {
	size_t in_nrows, in_ncols;
	in.size(in_nrows, in_ncols);

	if (nrows!=in_nrows || ncols!=in_ncols)
		return false;

	for (size_t ii=0; ii<vecData.size(); ii++)
		if ( !IOUtils::checkEpsilonEquality( vecData[ii] , in.vecData[ii], epsilon_mtr) ) return false;

	return true;
}

bool Matrix::operator!=(const Matrix& in) const {
	return !(*this==in);
}

Matrix& Matrix::operator+=(const Matrix& rhs) {
	//check dimensions compatibility
	if (nrows!=rhs.nrows || ncols!=rhs.ncols) {
		std::ostringstream tmp;
		tmp << "Trying to add two matrix with incompatible dimensions: ";
		tmp << "(" << nrows << "," << ncols << ") * ";
		tmp << "(" << rhs.nrows << "," << rhs.ncols << ")";
		throw IOException(tmp.str(), AT);
	}

	//fill sum matrix
	for (size_t ii=0; ii<vecData.size(); ii++)
		vecData[ii] += rhs.vecData[ii];

	return *this;
}

const Matrix Matrix::operator+(const Matrix& rhs) const {
	Matrix result = *this;
	result += rhs; //already implemented

	return result;
}

Matrix& Matrix::operator+=(const double& rhs) {
	//fill sum matrix
	for (size_t ii=0; ii<vecData.size(); ii++)
		vecData[ii] += rhs;

	return *this;
}

const Matrix Matrix::operator+(const double& rhs) const {
	Matrix result = *this;
	result += rhs; //already implemented

	return result;
}

Matrix& Matrix::operator-=(const Matrix& rhs) {
	//check dimensions compatibility
	if (nrows!=rhs.nrows || ncols!=rhs.ncols) {
		std::ostringstream tmp;
		tmp << "Trying to substract two matrix with incompatible dimensions: ";
		tmp << "(" << nrows << "," << ncols << ") * ";
		tmp << "(" << rhs.nrows << "," << rhs.ncols << ")";
		throw IOException(tmp.str(), AT);
	}

	//fill sum matrix
	for (size_t ii=0; ii<vecData.size(); ii++)
		vecData[ii] -= rhs.vecData[ii];

	return *this;
}

const Matrix Matrix::operator-(const Matrix& rhs) const {
	Matrix result = *this;
	result -= rhs; //already implemented

	return result;
}

Matrix& Matrix::operator-=(const double& rhs) {
	*this += -rhs;

	return *this;
}

const Matrix Matrix::operator-(const double& rhs) const {
	Matrix result = *this;
	result += -rhs; //already implemented

	return result;
}

Matrix& Matrix::operator*=(const Matrix& rhs) {
	//check dimensions compatibility
	if (ncols!=rhs.nrows) {
		std::ostringstream tmp;
		tmp << "Trying to multiply two matrix with incompatible dimensions: ";
		tmp << "(" << nrows << "," << ncols << ") * ";
		tmp << "(" << rhs.nrows << "," << rhs.ncols << ")";
		throw IOException(tmp.str(), AT);
	}

	//create new matrix
	Matrix result(nrows, rhs.ncols);

	//fill product matrix
	for (size_t i=1; i<=result.nrows; i++) {
		for (size_t j=1; j<=result.ncols; j++) {
			double sum=0.;
			for (size_t idx=1; idx<=ncols; idx++) {
				sum += operator()(i,idx) * rhs(idx,j);
			}
			result(i,j) = sum;
		}
	}

	*this = result;
	return *this;
}

const Matrix Matrix::operator*(const Matrix& rhs) const {
	Matrix result = *this;
	result *= rhs; //already implemented

	return result;
}

Matrix& Matrix::operator*=(const double& rhs) {
	for (size_t ii=0; ii<vecData.size(); ii++)
		vecData[ii] *= rhs;

	return *this;
}

const Matrix Matrix::operator*(const double& rhs) const {
	Matrix result = *this;
	result *= rhs; //already implemented

	return result;
}

Matrix& Matrix::operator/=(const double& rhs) {
	*this *= (1./rhs);
	return *this;
}

const Matrix Matrix::operator/(const double& rhs) const {
	Matrix result = *this;
	result *= 1./rhs; //already implemented

	return result;
}

double Matrix::scalar(const Matrix& m) {
	return m.scalar();
}

double Matrix::scalar() const {
	if (ncols!=1 || nrows!=1) {
		std::ostringstream tmp;
		tmp << "Trying to get scalar value of a non (1x1) matrix ";
		tmp << "(" << nrows << "," << ncols << ") !";
		throw IOException(tmp.str(), AT);
	}
	return operator()(1,1);
}

double Matrix::dot(const Matrix& A, const Matrix& B) {
	size_t Acols, Arows, Bcols, Brows;
	A.size(Arows, Acols);
	B.size(Brows, Bcols);

	if (Acols!=1 || Bcols!=1) {
		std::ostringstream tmp;
		tmp << "Trying to get dot product of non vector matrix ";
		tmp << "(" << Arows << "," << Acols << ") · ";
		tmp << "(" << Brows << "," << Bcols << ") · ";
		throw IOException(tmp.str(), AT);
	}
	if (Arows!=Brows) {
		std::ostringstream tmp;
		tmp << "Trying to get dot product of incompatible matrix ";
		tmp << "(" << Arows << "," << Acols << ") · ";
		tmp << "(" << Brows << "," << Bcols << ") · ";
		throw IOException(tmp.str(), AT);
	}

	double sum=0.;
	for (size_t i=1; i<=Arows; i++) {
		sum += A(i,1)*B(i,1);
	}

	return sum;
}

Matrix Matrix::T(const Matrix& m) {
	return m.getT();
}

Matrix Matrix::getT() const {
	Matrix result(ncols, nrows);
	for (size_t i=1; i<=result.nrows; i++) {
		for (size_t j=1; j<=result.ncols; j++) {
			result(i,j) = operator()(j,i);
		}
	}
	return result;
}

void Matrix::T() {
	Matrix tmp(*this);
	*this = tmp.getT();
}

double Matrix::det() const {
	if (nrows!=ncols) {
		std::ostringstream tmp;
		tmp << "Trying to calculate the determinant of a non-square matrix ";
		tmp << "(" << nrows << "," << ncols << ") !";
		throw IOException(tmp.str(), AT);
	}
	Matrix L,U;
	if (LU(L,U)==false) return 0.;

	double product=1.;
	for (size_t i=1; i<=nrows; i++) product *= U(i,i);

	return product;
}

bool Matrix::LU(Matrix& L, Matrix& U) const {
//Dolittle algorithm, cf http://math.fullerton.edu/mathews/numerical/linear/dol/dol.html
//HACK: there is no permutation matrix, so it might not be able to give a decomposition...
	if (nrows!=ncols) {
		std::ostringstream tmp;
		tmp << "Trying to calculate the LU decomposition of a non-square matrix ";
		tmp << "(" << nrows << "," << ncols << ") !";
		throw IOException(tmp.str(), AT);
	}

	const size_t n = nrows;
	U.clear();
	U = *this;
	L.identity(n, 1.); //initialized as identity matrix, then populated
	const Matrix& A = *this;

	for (size_t k=1; k<=n; k++) {
		//compute U elements
		for (size_t j=1; j<k; j++) {
			U(k,j) = 0.;
		}
		for (size_t j=k; j<=n; j++) {
			double sum=0.;
			for (size_t m=1; m<=(k-1); m++) sum += L(k,m)*U(m,j);
			U(k,j) = A(k,j) - sum;
		}

		if ( k<n && IOUtils::checkEpsilonEquality(U(k,k), 0., epsilon) ) return false; //we can not compute L
		//compute L elements
		for (size_t i=k+1; i<=n; i++) {
			double sum=0.;
			for (size_t m=1; m<=(k-1); m++) sum += L(i,m)*U(m,k);
			L(i,k) = (A(i,k) - sum) / U(k,k);
		}
	}
	return true;
}

Matrix Matrix::getInv() const {
//This uses an LU decomposition followed by backward and forward solving for the inverse
//See for example Press, William H.; Flannery, Brian P.; Teukolsky, Saul A.; Vetterling, William T. (1992), "LU Decomposition and Its Applications", Numerical Recipes in FORTRAN: The Art of Scientific Computing (2nd ed.), Cambridge University Press, pp. 34–42
	if (nrows!=ncols) {
		std::ostringstream tmp;
		tmp << "Trying to invert a non-square matrix ";
		tmp << "(" << nrows << "," << ncols << ") !";
		throw IOException(tmp.str(), AT);
	}
	const size_t n = nrows;

	Matrix U;
	Matrix L;
	if (LU(L, U)==false) {
		throw IOException("LU decomposition of given matrix not possible", AT);
	}

	//we solve AX=I with X=A-1. Since A=LU, then LUX = I
	//we start by forward solving LY=I with Y=UX
	Matrix Y(n, n);
	for (size_t i=1; i<=n; i++) {
		if (IOUtils::checkEpsilonEquality(L(i,i), 0., epsilon)) {
			throw IOException("The given matrix can not be inversed", AT);
		}
		Y(i,i) = 1./L(i,i); //j==i
		for (size_t j=1; j<i; j++) { //j<i
			double sum=0.;
			for (size_t k=i-1; k>=1; k--) { //equivalent to 1 -> i-1
				sum += L(i,k) * Y(k,j);
			}
			Y(i,j) = -1./L(i,i) * sum;
		}
		for (size_t j=i+1; j<=n; j++) { //j>i
			Y(i,j) = 0.;
		}
	}

	//now, we backward solve UX=Y
	Matrix X(n,n);
	for (size_t i=n; i>=1; i--) { //lines
		if (IOUtils::checkEpsilonEquality(U(i,i), 0., epsilon)) { //HACK: actually, only U(n,n) needs checking
			throw IOException("The given matrix is singular and can not be inversed", AT);
		}
		for (size_t j=1; j<=n; j++) { //lines
			double sum=0.;
			for (size_t k=i+1; k<=n; k++) {
				sum += U(i,k) * X(k,j);
			}
			X(i,j) = (Y(i,j) - sum) / U(i,i);
		}
	}

	return X;
}

bool Matrix::inv() {
//same as getInv() const but we write the final result on top of the input matrix
	if (nrows!=ncols) {
		std::ostringstream tmp;
		tmp << "Trying to invert a non-square matrix ";
		tmp << "(" << nrows << "," << ncols << ") !";
		throw IOException(tmp.str(), AT);
	}
	const size_t n = nrows;

	Matrix U;
	Matrix L;
	if (LU(L, U)==false) {
		return false;
	}

	//we solve AX=I with X=A-1. Since A=LU, then LUX = I
	//we start by forward solving LY=I with Y=UX
	Matrix Y(n, n);
	for (size_t i=1; i<=n; i++) {
		if (IOUtils::checkEpsilonEquality(L(i,i), 0., epsilon)) {
			return false;
		}
		Y(i,i) = 1./L(i,i); //j==i
		for (size_t j=1; j<i; j++) { //j<i
			double sum=0.;
			for (size_t k=i-1; k>=1; k--) { //equivalent to 1 -> i-1
				sum += L(i,k) * Y(k,j);
			}
			Y(i,j) = -1./L(i,i) * sum;
		}
		for (size_t j=i+1; j<=n; j++) { //j>i
			Y(i,j) = 0.;
		}
	}

	//now, we backward solve UX=Y
	Matrix& X = *this; //we write the solution over the input matrix
	for (size_t i=n; i>=1; i--) { //lines
		if (IOUtils::checkEpsilonEquality(U(i,i), 0., epsilon)) { //actually, only U(n,n) needs checking
			return false;
		}
		for (size_t j=1; j<=n; j++) { //lines
			double sum=0.;
			for (size_t k=i+1; k<=n; k++) {
				sum += U(i,k) * X(k,j);
			}
			X(i,j) = (Y(i,j) - sum) / U(i,i);
		}
	}

	return true;
}

bool Matrix::solve(const Matrix& A, const Matrix& B, Matrix& X) {
//This uses an LU decomposition followed by backward and forward solving for A·X=B
	size_t Anrows,Ancols, Bnrows, Bncols;
	A.size(Anrows, Ancols);
	if (Anrows!=Ancols) {
		std::ostringstream tmp;
		tmp << "Trying to solve A·X=B with A non square matrix ";
		tmp << "(" << Anrows << "," << Ancols << ") !";
		throw IOException(tmp.str(), AT);
	}
	B.size(Bnrows, Bncols);
	if (Anrows!=Bnrows)  {
		std::ostringstream tmp;
		tmp << "Trying to solve A·X=B with A and B of incompatible dimensions ";
		tmp << "(" << Anrows << "," << Ancols << ") and (";
		tmp << "(" << Bnrows << "," << Bncols << ") !";
		throw IOException(tmp.str(), AT);
	}
	const size_t n = Anrows;
	const size_t m = Bncols;

	Matrix U;
	Matrix L;
	if (A.LU(L, U)==false) {
		return false;
	}

	//we solve AX=B. Since A=LU, then LUX = B
	//we start by forward solving LY=B with Y=UX
	Matrix Y(n, m);
	for (size_t i=1; i<=n; i++) {
		if (IOUtils::checkEpsilonEquality(L(i,i), 0., epsilon)) {
			return false;
		}
		for (size_t j=1; j<=m; j++) {
			double sum=0.;
			for (size_t k=1; k<i; k++) {
				sum += L(i,k) * Y(k,j);
			}
			Y(i,j) = (B(i,j) - sum) / L(i,i);
		}
	}

	//now, we backward solve UX=Y
	X.resize(n,m); //we need to ensure that X has the correct dimensions
	for (size_t i=n; i>=1; i--) { //lines
		if (IOUtils::checkEpsilonEquality(U(i,i), 0., epsilon)) { //actually, only U(n,n) needs checking
			//singular matrix
			return false;
		}
		for (size_t j=1; j<=m; j++) {
			double sum = 0.;
			for (size_t k=i+1; k<=n; k++) {
				sum += U(i,k) * X(k,j);
			}
			X(i,j) = (Y(i,j) - sum) / U(i,i);
		}
	}

	return true;
}

Matrix Matrix::solve(const Matrix& A, const Matrix& B) {
//This uses an LU decomposition followed by backward and forward solving for A·X=B
	Matrix X;
	if (!solve(A, B, X))
		throw IOException("Matrix inversion failed!", AT);
	return X;
}

bool Matrix::TDMA_solve(const Matrix& A, const Matrix& B, Matrix& X)
{ //Thomas algorithm for tridiagonal matrix solving of A·X=B
	size_t Anrows,Ancols, Bnrows, Bncols;
	A.size(Anrows, Ancols);
	if (Anrows!=Ancols) {
		std::ostringstream tmp;
		tmp << "Trying to solve A·X=B with A non square matrix ";
		tmp << "(" << Anrows << "," << Ancols << ") !";
		throw IOException(tmp.str(), AT);
	}
	B.size(Bnrows, Bncols);
	if (Anrows!=Bnrows)  {
		std::ostringstream tmp;
		tmp << "Trying to solve A·X=B with A and B of incompatible dimensions ";
		tmp << "(" << Anrows << "," << Ancols << ") and (";
		tmp << "(" << Bnrows << "," << Bncols << ") !";
		throw IOException(tmp.str(), AT);
	}
	if (Bncols!=1) {
		std::ostringstream tmp;
		tmp << "Trying to solve A·X=B but B is not a vector! It is ";
		tmp << "(" << Bnrows << "," << Bncols << ") !";
		throw IOException(tmp.str(), AT);
	}

	const size_t n = Anrows;
	std::vector<double> b(n+1), c(n+1), v(n+1); //so we can keep the same index as for the matrix

	b[1] = A(1,1); v[1] = B(1,1); //otherwise they would never be defined
	for (size_t i=2; i<=n; i++) {
		if (IOUtils::checkEpsilonEquality(b[i-1], 0., epsilon))
			return false;
		const double b_i = A(i,i);
		const double v_i = B(i,1);
		const double a_i = A(i,i-1);
		const double m = a_i / b[i-1];
		c[i-1] = A(i-1,i);
		b[i] = b_i - m * c[i-1];
		v[i] = v_i - m * v[i-1];
	}

	X.resize(n,1); //we need to ensure that X has the correct dimensions
	X(n,1) = v[n] / b[n];
	for (size_t i=n-1; i>=1; i--) {
		X(i,1) = ( v[i] - c[i]*X(i+1,1) ) / b[i];
	}

	return true;
}

Matrix Matrix::TDMA_solve(const Matrix& A, const Matrix& B) {
//This uses the Thomas algorithm for tridiagonal matrix solving of A·X=B
	Matrix X;
	if (TDMA_solve(A, B, X))
		return X;
	else
		throw IOException("Matrix inversion failed!", AT);
}

bool Matrix::isIdentity() const {
	if (nrows!=ncols) {
		std::ostringstream tmp;
		tmp << "A non-square matrix ";
		tmp << "(" << nrows << "," << ncols << ") can not be the identity matrix!";
		throw IOException(tmp.str(), AT);
	}

	bool is_identity=true;
	for (size_t i=1; i<=nrows; i++) {
		for (size_t j=1; j<=ncols; j++) {
			const double val = operator()(i,j);
			if (i!=j) {
				if (IOUtils::checkEpsilonEquality(val,0.,epsilon_mtr)==false) {
					is_identity=false;
					break;
				}
			} else {
				if (IOUtils::checkEpsilonEquality(val,1.,epsilon_mtr)==false) {
					is_identity=false;
					break;
				}
			}
		}
	}

	return is_identity;
}

bool Matrix::isIdentity(const Matrix& A) {
	return A.isIdentity();
}

void Matrix::partialPivoting(std::vector<size_t>& pivot_idx) {
	pivot_idx.clear();

	//bad luck: if a row has several elements that are max of their columns,
	//we don't optimize its position. Ie: we can end up with a small element
	//on the diagonal
	for (size_t j=1; j<=ncols; j++) {
		const size_t old_i = j;
		const size_t new_i = findMaxInCol(j);
		if (new_i!=j) { //ie: pivoting needed
			swapRows(old_i, new_i);
			pivot_idx.push_back(new_i);
		} else
			pivot_idx.push_back(old_i);
	}
}

void Matrix::partialPivoting() {
	std::vector<size_t> pivot_idx;
	partialPivoting(pivot_idx);
}

void Matrix::maximalPivoting() {
	std::vector<size_t> pivot_idx;
	Matrix tmp( *this );

	for (size_t i=1; i<=nrows; i++) {
		const double scale = operator()(i,findMaxInRow(i));
		for (size_t j=1; j<=ncols; j++) {
			operator()(i,j) /= scale;
		}
	}
	tmp.partialPivoting(pivot_idx);

	//pivot on original matrix //HACK: not finished yet!
	throw IOException("Method not implemented yet!!", AT);
}

/*void Matrix::bidiagonalize() {
	//Matrix e(1,ncols);
	std::vector<double> e(ncols+1); //so we remain compatible with matrix index
	double g=0., x=0.;

	for (size_t i=1; i<=ncols; i++) {
		e[i]=g; s=0.; l=i+1;
		for (size_t j=i; j<=m; j++) s += ( operator()(i,j)*operator()(i,j) );
	}
}*/

//return the index of the line containing the highest absolute value at column col
size_t Matrix::findMaxInCol(const size_t &col) {
	size_t row_idx = 0;
	double max_val=0.;

	for (size_t i=1; i<=nrows; i++) {
		const double val = fabs( operator()(i,col) );
		if ( val>max_val) {
			max_val=val;
			row_idx=i;
		}
	}
	return row_idx;
}

//return the index of the line containing the highest absolute value at column col
size_t Matrix::findMaxInRow(const size_t &row) {
	size_t col_idx = 0;
	double max_val=0.;

	for (size_t j=1; j<=ncols; j++) {
		const double val = fabs( operator()(row,j) );
		if ( val>max_val) {
			max_val=val;
			col_idx=j;
		}
	}
	return col_idx;
}


void Matrix::swapRows(const size_t &i1, const size_t &i2) {
	for (size_t j=1; j<=ncols; j++) {
		const double tmp = operator()(i2,j);
		operator()(i2,j) = operator()(i1,j);
		operator()(i1,j) = tmp;
	}
}

} //end namespace
