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
#ifndef __LIBFIT1DCORE_H__
#define __LIBFIT1DCORE_H__

#include <meteoio/IOExceptions.h>
#include <meteoio/IOUtils.h>
#include <meteoio/dataClasses/Matrix.h>

#include <vector>

namespace mio {

//class to use as an interface
class FitModel {
	public:
		FitModel() : Lambda(), X(), Y(), infoString(), regname(), nPts(0), nParam(0), min_nb_pts(0), fit_ready(false) {}
		virtual ~FitModel() {}
		virtual void setData(const std::vector<double>& in_X, const std::vector<double>& in_Y) = 0;
		void setGuess(const std::vector<double>& lambda_in);
		virtual void setLapseRate(const double& /*lapse_rate*/) {throw InvalidArgumentException("Lapse rates can only be forced for linear regressions!", AT);}
		virtual bool fit() = 0;
		virtual double f(const double& x) const = 0;
		void getParams(std::vector<double>& o_coefficients) const;
		std::string getName() const {return regname;}
		std::string getInfo() const;
		void setInfo(const std::string& info) {infoString=info;}
		FitModel& operator =(const FitModel& source);
		std::string toString() const;
	protected:
		virtual bool checkInputs() {return true;}

		std::vector<double> Lambda; //parameters of the fit
		std::vector<double> X; //X of input data set to fit
		std::vector<double> Y; //Y of input data set to fit
		std::string infoString;
		std::string regname; //model name
		size_t nPts; //number of data points
		size_t nParam; //number of parameters
		size_t min_nb_pts; //minimum number of data points
		bool fit_ready; //set to false if fit() must be called before using the fit
};

/**
 * @class FitLeastSquare
 * @brief A class to perform non-linear least square fitting.
 * It works on a time serie and uses matrix arithmetic to perform an arbitrary fit
 * (see http://mathworld.wolfram.com/NonlinearLeastSquaresFitting.html).
 *
 * @ingroup stats
 * @author Mathias Bavay
 * @date   2011-01-20
 */
class FitLeastSquare : public FitModel {
 	public:
		FitLeastSquare() {}
		void setData(const std::vector<double>& in_X, const std::vector<double>& in_Y);
		bool fit();
		virtual double f(const double& x) const = 0;

	protected:
		virtual void setDefaultGuess(); //set defaults guess values. Called by setData

	private:
		bool checkInputs();
		void initLambda();
		void initDLambda(Matrix& dLambda) const;
		double getDelta(const double& var) const;
		double DDer(const double& x, const size_t& index);
		bool computeFit();

		static const double lambda_init; //initial default guess
		static const double delta_init_abs; //initial delta, absolute
		static const double delta_init_rel; //initial delta, relative
		static const double eps_conv; //convergence criteria
		static const unsigned int max_iter; //maximum number of iterations
};

} //end namespace

#endif
