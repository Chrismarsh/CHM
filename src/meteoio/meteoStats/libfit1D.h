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
#ifndef __LIBFIT1D_H__
#define __LIBFIT1D_H__

#include <meteoio/IOExceptions.h>
#include <meteoio/meteoStats/libfit1DCore.h>

#include <vector>

namespace mio {

class Zero : public FitModel {
	public:
		Zero() {fit_ready = true; nParam = 0; min_nb_pts = 0; regname = "Zero";}
		void setData(const std::vector<double>& /*in_X*/, const std::vector<double>& /*in_Y*/) { }
		bool fit() { return true;}
		double f(const double& /*x*/) const {return 0.;}
};

class SimpleLinear : public FitModel {
	public:
		SimpleLinear() : fixed_lapse_rate(IOUtils::nodata) {fit_ready = false; nParam = 2; min_nb_pts = 2; regname = "SimpleLinear";}
		void setData(const std::vector<double>& in_X, const std::vector<double>& in_Y);
		bool fit();
		double f(const double& x) const;
		void setLapseRate(const double& in_lapse_rate) {fixed_lapse_rate = in_lapse_rate; fit_ready = false; min_nb_pts=1;}
	protected:
		bool checkInputs();
		double fixed_lapse_rate;
};

class NoisyLinear : public SimpleLinear {
	public:
		NoisyLinear() {fit_ready = false; nParam = 2; min_nb_pts = 2; regname = "NoisyLinear";}
		bool fit();
};

class SphericVario : public FitLeastSquare {
	public:
		SphericVario() {fit_ready = false; nParam = 3; min_nb_pts = 4; regname = "SphericVario";}
		void setDefaultGuess();
		double f(const double& x) const;
};

class LinVario : public FitLeastSquare {
	public:
		LinVario() {fit_ready = false; nParam = 2; min_nb_pts = 3; regname = "LinVario";}
		void setDefaultGuess();
		double f(const double& x) const;
};

class ExpVario : public FitLeastSquare {
	public:
		ExpVario() {fit_ready = false; nParam = 3; min_nb_pts = 4; regname = "ExpVario";}
		void setDefaultGuess();
		double f(const double& x) const;
};

class RatQuadVario : public FitLeastSquare {
	public:
		RatQuadVario() {fit_ready = false; nParam = 3; min_nb_pts = 4; regname = "RatQuadVario";}
		void setDefaultGuess();
		double f(const double& x) const;
};

class LinearLS : public FitLeastSquare {
	public:
		LinearLS() {fit_ready = false; nParam = 2; min_nb_pts = 3; regname = "LinearLS";}
		void setDefaultGuess();
		double f(const double& x) const;
};

class Quadratic : public FitLeastSquare {
	public:
		Quadratic() {fit_ready = false; nParam = 3; min_nb_pts = 4; regname = "Quadratic";}
		void setDefaultGuess();
		double f(const double& x) const;
};

/**
 * @class Fit1D
 * @brief A class to perform 1D regressions
 * It works on a time serie and uses either ad-hoc methods or matrix arithmetic to perform an arbitrary fit.
 * Currently, the following models are supported:
 * - Specific fits:
 *    - SimpleLinear
 *    - NoisyLinear
 * - Least Square fits:
 *    - SphericVario
 *    - LinVario
 *    - ExpVario
 *    - RatQuadVario
 *    - LinearLS
 *    - Quadratic
 *
 * The various variogram models can be found in
 * <i>"Statistics for spatial data"</i>, Noel A. C. Cressie, John Wiley & Sons, revised edition, 1993, pp63.
 *
 *
 * @ingroup stats
 * @author Mathias Bavay
 * @date   2011-01-20
 */
class Fit1D {
 	public:
		///Keywords for regression model
		typedef enum REGRESSION {
			ZERO, ///< always return zero (this is a way to disable detrending)
			SIMPLE_LINEAR, ///< basic, cheap linear fit
			NOISY_LINEAR, ///< same as SIMPLE_LINEAR but trying to remove outliers
			LINVARIO, ///< linear variogram
			EXPVARIO, ///< exponential variogram
			SPHERICVARIO, ///< spherical variogram
			RATQUADVARIO, ///< rational quadratic variogram
			LINEARLS, ///< linear, using least squares
			QUADRATIC ///< quadratic
		} regression;

		/**
		* @brief Empty Constructor. The model must be set afterwards.
		* If the model has not been set before calling other methods, a NULL pointer exception will be thrown.
		*/
		Fit1D() : model(NULL) {}

		/**
		* @brief Constructor.
		* @param regType regression model to use
		* @param in_X vector of data points abscissae
		* @param in_Y vector of data points ordinates
		* @param updatefit should the fit be redone? (default=true, otherwise you must manually call fit())
		*/
		Fit1D(const regression& regType, const std::vector<double>& in_X, const std::vector<double>& in_Y, const bool& updatefit=true);

		/**
		* @brief Constructor for user provided model type.
		* @param regType regression model to use
		* @param in_X vector of data points abscissae
		* @param in_Y vector of data points ordinates
		* @param updatefit should the fit be redone? (default=true, otherwise you must manually call fit())
		*/
		Fit1D(const std::string& regType, const std::vector<double>& in_X, const std::vector<double>& in_Y, const bool& updatefit=true);

		/**
		* @brief Copy constructor.
		* @param i_fit Object to copy
		*/
		Fit1D(const Fit1D& i_fit);

		~Fit1D() {delete model;}

		/**
		* @brief Set or reset the regression model.
		* @param i_regType regression model to use
		* @param in_X vector of data points abscissae
		* @param in_Y vector of data points ordinates
		* @param updatefit should the fit be redone? (default=true, otherwise you must manually call fit())
		* @return false if could not compute the parameters. if !updatefit, always return true
		*/
		bool setModel(const regression& i_regType, const std::vector<double>& in_X, const std::vector<double>& in_Y, const bool& updatefit=true);

		/**
		* @brief Set or reset the regression model.
		* @param i_regType regression model to use
		* @param in_X vector of data points abscissae
		* @param in_Y vector of data points ordinates
		* @param updatefit should the fit be redone? (default=true, otherwise you must manually call fit())
		* @return false if could not compute the parameters. if !updatefit, always return true
		*/
		bool setModel(const std::string& i_regType, const std::vector<double>& in_X, const std::vector<double>& in_Y, const bool& updatefit=true);

		/**
		* @brief Provide a set of initial values for the model parameters.
		* The model can be used right after providing the guesses, and it would use those guesses as parameters,
		* thus allowing the user to force his model parameters.
		* @param lambda_in one initial value per model parameter
		*/
		void setGuess(const std::vector<double>& lambda_in) {model->setGuess(lambda_in);}

		/**
		* @brief Set a forced lapse rate for linear regressions
		* This will throw an exception for all other regression models!
		* @param lapse_rate lapse rate to set
		*/
		void setLapseRate(const double& lapse_rate) {model->setLapseRate(lapse_rate);}

		/**
		* @brief Compute the regression parameters
		* @return false if could not compute the parameters
		*/
		bool fit() {return model->fit();}

		/**
		* @brief Calculate a value using the computed least square fit.
		* The fit has to be computed before.
		* @param x abscissa
		* @return f(x) using the computed least square fit
		*/
		double f(const double& x) const {return model->f(x);}

		/**
		* @brief Calculate the parameters of the fit.
		* The fit has to be computed before.
		* @param coefficients vector containing the coefficients
		*/
		void getParams(std::vector<double>& coefficients) const {model->getParams(coefficients);}

		/**
		* @brief Return the name of the fit model.
		* @return model name
		*/
		std::string getName() const {return model->getName();}

		/**
		* @brief Return a string of information about the fit.
		* The fit has to be computed before.
		* @return info string
		*/
		std::string getInfo() const {return model->getInfo();}

		/**
		* @brief Set the information string.
		* This is useful to append some extra information to the information string.
		* This should be called <b>after</b> computing the fit (otherwise it will be overwritten).
		* @param info string
		*/
		void setInfo(const std::string& info) {model->setInfo(info);}

		Fit1D& operator =(const Fit1D& source);

		/**
		* @brief Calculate a value using the computed least square fit.
		* The fit has to be computed before.
		* @param x abscissa
		* @return f(x) using the computed least square fit
		*/
		double operator ()(const double& x) const {return model->f(x);}

		std::string toString() const {return model->toString();}

	private:
		FitModel *model;
};

} //end namespace

#endif
