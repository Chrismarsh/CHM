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
#include <meteoio/meteoStats/libfit1D.h>
#include <meteoio/meteoStats/libinterpol1D.h>
#include <meteoio/IOExceptions.h>
#include <meteoio/MathOptim.h>

#include <cmath>
#include <algorithm>
#include <iomanip>

using namespace std;

namespace mio {

Fit1D::Fit1D(const regression& regType, const std::vector<double>& in_X, const std::vector<double>& in_Y, const bool& updatefit) : model(NULL) 
{
	const bool status = setModel(regType, in_X, in_Y, updatefit);
	if (updatefit && status==false) 
		throw NoDataException("The provided data was insufficient when constructing the regression model '"+model->getName()+"'", AT);
}

Fit1D::Fit1D(const std::string& regType, const std::vector<double>& in_X, const std::vector<double>& in_Y, const bool& updatefit) : model(NULL) 
{
	const bool status = setModel(regType, in_X, in_Y, updatefit);
	if (updatefit && status==false) 
		throw NoDataException("The provided data was insufficient when constructing the regression model '"+model->getName()+"'", AT);
}

Fit1D::Fit1D(const Fit1D& i_fit) : model(NULL) { //HACK: the pointer could not be valid anymore
	*this = i_fit;
}

Fit1D& Fit1D::operator=(const Fit1D& source) { //HACK: the pointer could not be valid anymore
	if (this != &source) {
		model = new SimpleLinear; //this is only for memory allocation
		*model = *(source.model); //copy what is pointed to
	}
	return *this;
}

bool Fit1D::setModel(const std::string& i_regType, const std::vector<double>& in_X, const std::vector<double>& in_Y, const bool& updatefit)
{
	regression regType;
	if (i_regType=="ZERO") regType=ZERO;
	else if (i_regType=="SIMPLE_LINEAR") regType=SIMPLE_LINEAR;
	else if (i_regType=="NOISY_LINEAR") regType=NOISY_LINEAR;
	else if (i_regType=="LINVARIO") regType=LINVARIO;
	else if (i_regType=="EXPVARIO") regType=EXPVARIO;
	else if (i_regType=="SPHERICVARIO") regType=SPHERICVARIO;
	else if (i_regType=="RATQUADVARIO") regType=RATQUADVARIO;
	else if (i_regType=="LINEARLS") regType=LINEARLS;
	else if (i_regType=="QUADRATIC") regType=QUADRATIC;
	else {
		throw IOException("The regression algorithm '"+i_regType+"' is not implemented" , AT);
	}

	return setModel(regType, in_X, in_Y, updatefit);
}

bool Fit1D::setModel(const regression& regType, const std::vector<double>& in_X, const std::vector<double>& in_Y, const bool& updatefit)
{
	if (model!=NULL) delete model;

	if (regType==ZERO) model=new Zero;
	if (regType==SIMPLE_LINEAR) model=new SimpleLinear;
	if (regType==NOISY_LINEAR) model=new NoisyLinear;
	if (regType==LINVARIO) model=new LinVario;
	if (regType==EXPVARIO) model=new ExpVario;
	if (regType==SPHERICVARIO) model=new SphericVario;
	if (regType==RATQUADVARIO) model=new RatQuadVario;
	if (regType==LINEARLS) model=new LinearLS;
	if (regType==QUADRATIC) model=new Quadratic;

	//remove nodata points
	std::vector<double> X, Y;
	for (size_t ii=0; ii<in_X.size(); ii++) {
		if (in_X[ii]!=IOUtils::nodata && in_Y[ii]!=IOUtils::nodata) {
			X.push_back( in_X[ii] );
			Y.push_back( in_Y[ii] );
		}
	}

	model->setData(X, Y);
	if (updatefit)
		return fit();
	else
		return true;
}

//////////////////////////////////////////////////////////////
// regression models
void SimpleLinear::setData(const std::vector<double>& in_X, const std::vector<double>& in_Y) {
	X = in_X;
	Y = in_Y;

	fit_ready = false;
}

double SimpleLinear::f(const double& x) const {
	return Lambda.at(0)*x + Lambda.at(1);
}

bool SimpleLinear::fit()
{
	if (!checkInputs())
		return false;

	Lambda.clear();
	double a,b,r;
	std::string mesg;
	std::ostringstream ss;

	if (fixed_lapse_rate==IOUtils::nodata) {
		Interpol1D::LinRegression(X, Y, a, b, r, mesg);
		ss << mesg << "Computed regression with " << regname << " model - r=" << std::setprecision(2) << r;
	} else {
		a = fixed_lapse_rate;
		if (a!=0.) {
			Interpol1D::LinRegression(X, Y, a, b, r, mesg, true);
			ss << mesg << "Computed regression with " << regname << " model ";
			ss << "(fixed lapse rate=" << a << ") - r=" << std::setprecision(2) << r;
		} else {
			a=0.;
			b=0.;
			ss << mesg << "Computed regression with " << regname << " model";
		}
	}
	Lambda.push_back(a);
	Lambda.push_back(b);
	infoString = ss.str();
	fit_ready = true;
	return true;
}

bool NoisyLinear::fit()
{
	if (!checkInputs())
		return false;

	Lambda.clear();
	double a,b,r;
	std::string mesg;
	std::ostringstream ss;

	if (fixed_lapse_rate==IOUtils::nodata) {
		Interpol1D::NoisyLinRegression(X, Y, a, b, r, mesg);
		ss << mesg  << "Computed regression with " << regname << " model - r=" << std::setprecision(2) << r;
	} else {
		a = fixed_lapse_rate;
		if (a!=0.) {
			Interpol1D::NoisyLinRegression(X, Y, a, b, r, mesg, true);
			ss << mesg  << "Computed regression with " << regname << " model ";
			ss << "(fixed lapse rate=" << a << ") - r=" << std::setprecision(2) << r;
		} else {
			a=0.;
			b=0.;
			ss << mesg << "Computed regression with " << regname << " model";
		}
	}
	Lambda.push_back(a);
	Lambda.push_back(b);
	infoString = ss.str();
	fit_ready = true;
	return true;
}

//regression models using the standard least square algorithm
double SphericVario::f(const double& x) const {
	//c0>=0, cs>=0, as>=0
	const double c0 = Lambda.at(0);
	const double cs = Lambda.at(1);
	const double as = Lambda.at(2);

	if (x==0) return 0;

	const double abs_x = fabs(x);
	if (abs_x>0 && abs_x<=as) {
		const double val = abs_x/as;
		const double y = c0 + cs * ( 1.5*val - 0.5*Optim::pow3(val) );
		return y;
	} else {
		return (c0+cs);
	}
}

void SphericVario::setDefaultGuess() {
	Lambda.push_back( *min_element(Y.begin(), Y.end()) );
	Lambda.push_back( *max_element(Y.begin(), Y.end()) );
	Lambda.push_back( *max_element(X.begin(), X.end()) );
}

double LinVario::f(const double& x) const {
	//c0>=0, b1>=0
	const double c0 = Lambda.at(0);
	const double bl = Lambda.at(1);

	if (x==0) {
		return 0;
	} else {
		const double y = c0 + bl * abs(x);
		return y;
	}
}

void LinVario::setDefaultGuess() {
	double xzero=X[0];
	size_t xzero_idx=0;
	for (size_t i=1; i<X.size(); i++) {
		if (abs(X[i])<xzero) { xzero=X[i]; xzero_idx=i;}
	}
	const double slope = Interpol1D::arithmeticMean( Interpol1D::derivative(X, Y) );
	Lambda.push_back( Y[xzero_idx] );
	Lambda.push_back( slope );
}

double ExpVario::f(const double& x) const {
	//c0>=0, ce>=0, ae>=0
	const double c0 = Lambda.at(0);
	const double ce = Lambda.at(1);
	const double ae = Lambda.at(2);

	if (x==0) {
		return 0;
	} else {
		const double y = c0 + ce * (1. - exp(-abs(x)/ae) );
		return y;
	}
}

void ExpVario::setDefaultGuess() {
	double xzero=X[0];
	size_t xzero_idx=0;
	for (size_t i=1; i<X.size(); i++) {
		if (abs(X[i])<xzero) { xzero=X[i]; xzero_idx=i;}
	}
	Lambda.push_back( Y[xzero_idx] );
	Lambda.push_back( Y.back() - Y[xzero_idx] );
	Lambda.push_back( 1. );
}

double RatQuadVario::f(const double& x) const {
	//c0>=0, cr>=0, ar>=0
	const double c0 = Lambda.at(0);
	const double cr = Lambda.at(1);
	const double ar = Lambda.at(2);

	if (x==0) {
		return 0;
	} else {
		const double y = c0 + cr*x*x / (1. + x*x/ar);
		return y;
	}
}

void RatQuadVario::setDefaultGuess() {
	double xzero=X[0];
	size_t xzero_idx=0;
	for (size_t i=1; i<X.size(); i++) {
		if (abs(X[i])<xzero) { xzero=X[i]; xzero_idx=i;}
	}
	Lambda.push_back( Y[xzero_idx] );
	Lambda.push_back( *( std::max_element( Y.begin(), Y.end() ) ) );
	Lambda.push_back( 1. );
}

double LinearLS::f(const double& x) const {
	const double y = Lambda.at(0)*x + Lambda.at(1); //Lambda is a vector
	return y;
}

void LinearLS::setDefaultGuess() {
	double xzero=X[0];
	size_t xzero_idx=0;
	for (size_t i=1; i<X.size(); i++) {
		if (abs(X[i])<xzero) { xzero=X[i]; xzero_idx=i;}
	}

	const double slope = Interpol1D::arithmeticMean( Interpol1D::derivative(X, Y) );
	Lambda.push_back( slope );
	Lambda.push_back( Y[xzero_idx] );
}

double Quadratic::f(const double& x) const {
	const double y = Lambda.at(0)*x*x + Lambda.at(1)*x + Lambda.at(2); //Lambda is a vector
	return y;
}

void Quadratic::setDefaultGuess() {
	std::vector<double> der( Interpol1D::derivative(X, Y) );
	const double acc = 0.5 * Interpol1D::arithmeticMean( Interpol1D::derivative(X, der) );
	double xzero=der[0];
	size_t xzero_idx=0;
	for (size_t i=1; i<der.size(); i++) {
		if (abs(der[i])<xzero) { xzero=der[i]; xzero_idx=i;}
	}

	Lambda.push_back( acc ); //0
	Lambda.push_back( der[xzero_idx] ); //1
	if (acc>0.)
		Lambda.push_back( *( std::min_element( Y.begin(), Y.end() ) ) );
	else
		Lambda.push_back( *( std::max_element( Y.begin(), Y.end() ) ) );
}

} //namespace
