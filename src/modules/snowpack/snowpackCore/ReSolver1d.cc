/*
 *  SNOWPACK stand-alone
 *
 *  Copyright WSL Institute for Snow and Avalanche Research SLF, DAVOS, SWITZERLAND
*/
/*  This file is part of Snowpack.
    Snowpack is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Snowpack is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Snowpack.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <snowpack/snowpackCore/ReSolver1d.h>
#include <snowpack/snowpackCore/Snowpack.h>
#ifdef CLAPACK
	// Matching C data types with FORTRAN data types (taken from f2c.h):
	typedef long int integer;
	typedef double doublereal;

	// Declare the function interfaces with the LAPACK library (taken from clapack.h):
	extern "C" {
		/* Subroutine */ int dgesvd_(char *jobu, char *jobvt, integer *m, integer *n,
		doublereal *a, integer *lda, doublereal *s, doublereal *u, integer *
		ldu, doublereal *vt, integer *ldvt, doublereal *work, integer *lwork,
		integer *info);

		/* Subroutine */ int dgesdd_(char *jobz, integer *m, integer *n, doublereal *
		a, integer *lda, doublereal *s, doublereal *u, integer *ldu,
		doublereal *vt, integer *ldvt, doublereal *work, integer *lwork,
		integer *iwork, integer *info);

		/* Subroutine */ int dgtsv_(integer *n, integer *nrhs, doublereal *dl,
		doublereal *d__, doublereal *du, doublereal *b, integer *ldb, integer
		*info);
	}
#endif
#include <string.h>


using namespace std;
using namespace mio;

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wsign-conversion"
#endif

ReSolver1d::ReSolver1d(const SnowpackConfig& cfg)
           : surfacefluxrate(0.), soilsurfacesourceflux(0.), variant(),
             iwatertransportmodel_snow(BUCKET), iwatertransportmodel_soil(BUCKET),
             watertransportmodel_snow("BUCKET"), watertransportmodel_soil("BUCKET"), BottomBC(FREEDRAINAGE),
             sn_dt(IOUtils::nodata), useSoilLayers(false), water_layer(false)
{
	cfg.getValue("VARIANT", "SnowpackAdvanced", variant);

	// Defines whether soil layers are used
	cfg.getValue("SNP_SOIL", "Snowpack", useSoilLayers);

	//To build a thin top rain-water layer over a thin top ice layer, rocks, roads etc.
	cfg.getValue("WATER_LAYER", "SnowpackAdvanced", water_layer);

	//Calculation time step in seconds as derived from CALCULATION_STEP_LENGTH
	double calculation_step_length = cfg.get("CALCULATION_STEP_LENGTH", "Snowpack");
	sn_dt = M_TO_S(calculation_step_length);

	//Water transport model snow
	cfg.getValue("WATERTRANSPORTMODEL_SNOW", "SnowpackAdvanced", watertransportmodel_snow);
	iwatertransportmodel_snow=UNDEFINED;
	if (watertransportmodel_snow=="BUCKET") {
		iwatertransportmodel_snow=BUCKET;
	} else if (watertransportmodel_snow=="NIED") {
		iwatertransportmodel_snow=NIED;
	} else if (watertransportmodel_snow=="RICHARDSEQUATION") {
		iwatertransportmodel_snow=RICHARDSEQUATION;
	}

	//Water transport model soil
	cfg.getValue("WATERTRANSPORTMODEL_SOIL", "SnowpackAdvanced", watertransportmodel_soil);
	iwatertransportmodel_soil=UNDEFINED;
	if (watertransportmodel_soil=="BUCKET") {
		iwatertransportmodel_soil=BUCKET;
	} else if (watertransportmodel_soil=="NIED") {
		iwatertransportmodel_soil=NIED;
	} else if (watertransportmodel_soil=="RICHARDSEQUATION") {
		iwatertransportmodel_soil=RICHARDSEQUATION;
	}

	//Set lower boundary condition
	std::string tmp_lb_cond_waterflux;
	cfg.getValue("LB_COND_WATERFLUX", "SnowpackAdvanced", tmp_lb_cond_waterflux);
	if (tmp_lb_cond_waterflux=="DIRICHLET") {
		BottomBC=DIRICHLET;
	} else if (tmp_lb_cond_waterflux=="WATERTABLE") {
		BottomBC=WATERTABLE;
	} else if (tmp_lb_cond_waterflux=="FREEDRAINAGE") {
		BottomBC=FREEDRAINAGE;
	} else if (tmp_lb_cond_waterflux=="GRAVITATIONALDRAINAGE") {
		BottomBC=GRAVITATIONALDRAINAGE;
	} else if (tmp_lb_cond_waterflux=="SEEPAGE") {
		BottomBC=SEEPAGEBOUNDARY;
	}
}


/**
 * @brief Calculating pressure head from water content \n
 * The following function calculates the pressure head belonging to a given water content \n
 * @author Nander Wever
 * @param theta Water content (m^3/m^3)
 * @param theta_r Residual water content (m^3/m^3)
 * @param theta_s Saturated water content (m^3/m^3)
 * @param alpha Van Genuchten parameter
 * @param m Van Genuchten parameter
 * @param n Van Genuchten parameter
 * @param Sc Cut off value
 * @param h_e Air entry pressure
 * @param h_d Dry limit of pressure head
 */
double ReSolver1d::fromTHETAtoH(double theta, double theta_r, double theta_s, double alpha, double m, double n, double Sc, double h_e, double h_d)
{
	//Inverse of Van Genuchten (1980), Equation 21:
	double returnvalue;
	if (theta<=theta_r) {
		returnvalue=h_d;
	} else {
		if (theta >= theta_s) {
			returnvalue=h_e;
		} else {
			returnvalue=-1.*(1./alpha)*pow( (pow(Sc*((theta-theta_r)/(theta_s-theta_r)), (-1./m)) - 1.), (1./n));
		}
	}
	return returnvalue;
}


/**
 * @brief Calculating pressure head from water content when ice is present \n
 * The following function calculates the pressure head belonging to a given water content when ice is present \n
 * @author Nander Wever
 * @param theta Water content (m^3/m^3)
 * @param theta_r Residual water content (m^3/m^3)
 * @param theta_s Saturated water content (m^3/m^3)
 * @param alpha Van Genuchten parameter
 * @param m Van Genuchten parameter
 * @param n Van Genuchten parameter
 * @param Sc Cut off value
 * @param h_e Air entry pressure
 * @param h_d Dry limit of pressure head
 * @param theta_i Ice content (m^3/m^3)
 */
double ReSolver1d::fromTHETAtoHforICE(double theta, double theta_r, double theta_s, double alpha, double m, double n, double Sc, double h_e, double h_d, double theta_i)
{
	//To have same return value as fromTHETAtoH, call this function with theta_i==0.
	return fromTHETAtoH(theta+(theta_i*(Constants::density_ice/Constants::density_water)), theta_r, theta_s, alpha, m, n, Sc, h_e, h_d);
}


/**
 * @brief Calculating volumetric water content from pressure head \n
 * The following function calculates the volumetric water content belonging to a given pressure head \n
 * @author Nander Wever
 * @param h Pressure head (m)
 * @param theta_r Residual water content (m^3/m^3)
 * @param theta_s Saturated water content (m^3/m^3)
 * @param alpha Van Genuchten parameter
 * @param m Van Genuchten parameter
 * @param n Van Genuchten parameter
 * @param Sc Cut off value
 * @param h_e Air entry pressure
 * @param h_d Dry limit of pressure head
 */
double ReSolver1d::fromHtoTHETA(double h, double theta_r, double theta_s, double alpha, double m, double n, double Sc, double h_e)
{
	double returnvalue=0.;
	//Van Genuchten (1980), Equation 21:
	if (h>h_e) {		//Saturation
		returnvalue=theta_s;
	} else {
		returnvalue=theta_r+( (theta_s-theta_r)*(1./Sc)*pow(1.+pow((alpha*fabs(h)),n),(-1.*m)) );
	}
	return returnvalue;
}


/**
 * @brief Calculating volumetric water content from pressure head when ice is present \n
 * The following function calculates the volumetric water content belonging to a given pressure head when ice is present \n
 * @author Nander Wever
 * @param h Pressure head (m)
 * @param theta_r Residual water content (m^3/m^3)
 * @param theta_s Saturated water content (m^3/m^3)
 * @param alpha Van Genuchten parameter
 * @param m Van Genuchten parameter
 * @param n Van Genuchten parameter
 * @param Sc Cut off value
 * @param h_e Air entry pressure
 * @param h_d Dry limit of pressure head
 * @param theta_i Ice content (m^3/m^3)
 */
double ReSolver1d::fromHtoTHETAforICE(double h, double theta_r, double theta_s, double alpha, double m, double n, double Sc, double h_e, double theta_i)
{
	//To have same return value as fromHtoTHETA, call this function with theta_i==0.
	return fromHtoTHETA(h, theta_r, theta_s, alpha, m, n, Sc, h_e)-(theta_i*(Constants::density_ice/Constants::density_water));
}


/**
 * @brief Calculate air entry pressure head \n
 * Air entry pressure head in [m] that corresponds to a maximum pore size (using Young-Laplace Equation).\n
 * This is a required value for specifying water retention curves, see Ippisch et al. (2006).\n
 * @author Nander Wever
 * @param MaximumPoreSize Maximum pore size (diameter, not radius!) [m]
 * @param Temperature Temperature for determining surface tension [K]
 */
double ReSolver1d::AirEntryPressureHead(double MaximumPoreSize, double Temperature)
{
	//Surface tension is dependent on the temperature. Most simulations will be in the temperature range of -20 - +20 degC.
	//Source: http://en.wikipedia.org/wiki/Surface_tension
	//Surface tension of water in N/m.
	const double SurfaceTension = (Temperature > 293.)? 0.07197 : 0.07564; //Value for 25 degC vs for 0 degC
	const double delta_P=-1.*(2.*SurfaceTension)/(0.5*MaximumPoreSize);
	const double air_entry_head=delta_P/(Constants::density_water*Constants::g);

	return air_entry_head;
}


/**
 * @brief Solving system of equations using Thomas Algorithm \n
 * The following function solves a tridiagonal system of equations using Thomas Algorithm \n
 * @author http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
 * @param n number of equations
 * @param a sub-diagonal (means it is the diagonal below the main diagonal) -- indexed from 0..n-2
 * @param b the main diagonal
 * @param c sup-diagonal (means it is the diagonal above the main diagonal) -- indexed from 0..n-2
 * @param v right part
 * @param x the solution
 */
int ReSolver1d::TDMASolver (int n, double *a, double *b, double *c, double *v, double *x)
{
	// See: http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
	// This solver is very rapid, but has the problem that when elements of the matrix get very small, relative to each other, precision problems propagate.
	// A better option is to use the DGTSV solver from LAPACK, as it does partial pivoting, although even that is not always enough.
	// See for explanation: http://en.wikipedia.org/wiki/Pivot_element#Partial_and_complete_pivoting
        /**
         * n - number of equations
         * a - sub-diagonal (means it is the diagonal below the main diagonal) -- indexed from 0..n-2
         * b - the main diagonal
         * c - sup-diagonal (means it is the diagonal above the main diagonal) -- indexed from 0..n-2
         * v - right part
         * x - the solution
	 * Return value: 0 = succes, otherwise: error
         */
	if (b[n-1]==0.) return -1;		// This will cause division by 0, so return with error code.

	for (int i = 1; i < n; i++) {
		if (b[i-1]==0.) return -1;	// This will cause division by 0, so return with error code.
		double m = a[i-1]/b[i-1];
		b[i] = b[i] - m * c[i - 1];
		v[i] = v[i] - m*v[i-1];
	}

	x[n-1] = v[n-1]/b[n-1];

	for (int i = n - 2; i >= 0; --i)
		x[i] = (v[i] - c[i] * x[i+1]) / b[i];

	return 0;
}

/**
 * @brief Solving system of equations using matrix inversion \n
 * The following function solves a tridiagonal system of equations using \n
 * Moore-Penrose matrix inversion, using SVD, giving a pseudo-inverse of a. \n
 * @author Nander Wever
 * @param m number of rows in matrix
 * @param n number of columns in matrix
 * @param lda leading dimension of matrix
 * @param a pointer to top-left corner of matrix to inverse
 */
#ifdef CLAPACK
int ReSolver1d::pinv(int m, int n, int lda, double *a)
// NOTE: inverse matrix is returned in "a" as well. Make sure to send a copy to pinv if you want to keep the original matrix.
//
// Returned status by dgesvd_/dgesdd_:
//  INFO    (output) INTEGER
//          = 0:  successful exit.
//          < 0:  if INFO = -i, the i-th argument had an illegal value.
//          > 0:  For DGESVD: DBDSQR did not converge, INFO specifies how many. Superdiagonals of an intermediate bidiagonal form B did not converge to zero.
//		  For DGESDD: DBDSDC did not converge, updating process failed. The algorithm failed to compute an singular value. The update process of divide and conquer failed.
{
	//Switch for dgesvd/dgesdd
	bool useOptimezedSVD=true;	//True: dgesdd is used, which is faster, but requires more memory, compared to dgesvd. Note that when dgesdd failes, the function tries dgesvd. When that fails as well, the program is terminated.
					//I encountered some numerical issues with dgesdd, depending on settings in feenable (likely over- or underflow). So setting this switch to false gives a safe mode.
					//Note: there are some bugreports for dgesdd.

	//1D-Array for singular values:
	int nSV = m < n ? m : n;	//nSV: number of singular values
	double *s = (double *)calloc(nSV * sizeof *s, sizeof(double));

	//2D-Arrays for matrices U and Vt:
	double *u = (double *)calloc((m*m) * sizeof(*u), sizeof(double));
	double *vt = (double *)calloc((n*n) * sizeof(*vt), sizeof(double));

	//Setup workspace
	double workspaceSize;
	double *work = &workspaceSize;	//From the documentation, it's not clear to me whether with lwork=-1, work can be considered as an array, or will always return 1 value. If it might return an array, not enough memory is allocated here.
	int lwork = -1;
	int *iwork = (int *)calloc(8*nSV, sizeof(double));
	int info = 0;
	char jobu='A';
	char jobvt='A';

	//Call dgesvd_/dgesdd_ with lwork = -1 to determine optimal workspace size:
	if (useOptimezedSVD) {
		dgesdd_(&jobu, (integer*) &m, (integer*) &n, a, (integer*) &lda, s, u, (integer*) &m, vt, (integer*) &n, work, (integer*) &lwork, (integer*) iwork, (integer*) &info);
	} else {
		dgesvd_(&jobu, &jobvt, (integer*) &m, (integer*) &n, a, (integer*) &lda, s, u, (integer*) &m, vt, (integer*) &n, work, (integer*) &lwork, (integer*) &info);
	}
	if (info!=0) {
		printf ("ERROR in ReSolver1d.cc: SVD [determining optimal workspace size failed]: info = %d\n", info);
		throw;
	}

	//Now, workspace size is returned in work.
	lwork = int(work[0]);

	//Then allocate work again for the determined workspace size
	work = (double *)calloc(lwork * sizeof(*work), sizeof(double));

	// Call dgesdd_ to do the actual computation:
	if (useOptimezedSVD) {
		dgesdd_(&jobu, (integer*) &m, (integer*) &n, a, (integer*) &lda, s, u, (integer*) &m, vt, (integer*) &n, work, (integer*) &lwork, (integer*) iwork, (integer*) &info);
		if (info>0) {		//if info>0, there is a convergence problem, probably because the matrix is ill-conditioned. We can try dgesvd before giving up.
			printf ("ERROR in ReSolver1d.cc: DGESDD failed [info = %d]. Falling back on DGESVD.\n", info);
			dgesvd_(&jobu, &jobvt, (integer*) &m, (integer*) &n, a, (integer*) &lda, s, u, (integer*) &m, vt, (integer*) &n, work, (integer*) &lwork, (integer*) &info);
			//Now, workspace size is returned in work.
			lwork = int(work[0]);

			//Then allocate work again for the determined workspace size
			free(work);								//First free the previously allocated work.
			work = (double *)calloc(lwork * sizeof(*work), sizeof(double));		//Reallocate work.

			dgesvd_(&jobu, &jobvt, (integer*) &m, (integer*) &n, a, (integer*) &lda, s, u, (integer*) &m, vt, (integer*) &n, work, (integer*) &lwork, (integer*) &info);
		}
	} else {
		dgesvd_(&jobu, &jobvt, (integer*) &m, (integer*) &n, a, (integer*) &lda, s, u, (integer*) &m, vt, (integer*) &n, work, (integer*) &lwork, (integer*) &info);
	}
	if (info!=0) {
		printf ("ERROR in ReSolver1d.cc: DGESVD failed [info = %d]. No solution found with currently used time step.\n", info);
		return -1;
	}

	//dgesvd/dgesdd gave us:
	//A=U.S.Vt
	//We now need:
	//A'=V.(S').Ut
	//Because S is a diagonal matrix, S' is just (1/S)
	//Be careful: comparing the results with for example MATLAB is not directly possible.
	//Singular value decomposition is not unique, e.g. sign can differ for U and Vt, but due to multiplication, this sign reversal cancels out.

	std::vector<double> dinv(m*n, 0.);
	int i, j, k;

	//Calculate pseudo-inverse: Step 1: U.S':
	for (i = m-1; i >= 0; i--) {
		for (j = n-1; j >= 0; j--) {
			//dinv[j*m+i]=0.;							//This line is not necessary anymore, as we reset array to zero where it is declared.
			for (k = n-1; k >= 0; k--) {						//Do matrix multiplications
				if (j==k) {							//j==k: this is because s is actually a diagonal matrix, and therefore, represented as a vector.
					if(s[k]!=0.) {						//NANDER TODO HACK: I'm not happy with this solution, but how to circumvent underflow? I just want 0. in that case.
						dinv[j*m+i]+=(vt[i*n+k]*(double(1.)/s[k]));	//Note: vt needs to be transposed AND is in FORTRAN matrix notation, therefore, the indices appear normal again.
					}
				}
			}
		}
	}
	//Step 2: (U.S')*VT
	for (i = m-1; i >= 0; i--) {
		for (j = n-1; j >= 0; j--) {
			a[j*m+i]=0;
			for (k = n-1; k >= 0; k--) {				//Do matrix multiplications
				a[j*m+i]+=(dinv[k*m+j]*u[k*n+i]);		//Note: u needs to be transposed AND is in FORTRAN matrix notation.
			}
		}
	}

	//Free memory
	free(work);
	free(iwork);
	free(s);
	free(u);
	free(vt);

	return 0;
}
#else
int ReSolver1d::pinv(int /*m*/, int /*n*/, int /*lda*/, double */*a*/) {
	//Nothing so far
	throw IOException("This private method is not implemented when the LAPACK and BLAS libraries are not installed", AT);
}
#endif


/**
 * @brief Set soil parameters for a given soil type \n
 * Set soil parameters for a given soil type \n
 * @author Nander Wever
 * @param type Soil type
 * @param theta_r Residual water content
 * @param theta_soil Volumetric soil content
 * @param alpha Van Genuchten parameter
 * @param m Van Genuchten parameter
 * @param n Van Genuchten parameter
 * @param ksat Saturated hydraulic conductivity
 * @param he Air entry pressure
 */
void ReSolver1d::SetSoil(SoilTypes type, double *theta_r, double *theta_s, double *alpha, double *m, double *n, double *ksat, double *he)
{
	double MaximumPoreSize=0.;	//Maximum pore size (diameter) in [m]

	//Set van Genuchten parameters
	switch (type) {
		case ORGANIC:
			//Organic: Nemes (2001), Development of Soil Hydraulic Pedotransfer Functions on a European scale: Their Usefulness in the Assessment of Soil Quality.
			*theta_r=0.01;
			*theta_s=0.766;
			*alpha=1.3;
			*n=1.2039;
			*ksat=8.000/(365.*24.*60.*60.);
			MaximumPoreSize=0.005;
			break;

		//ROSETTA Class Average Hydraulic Parameters: http://ars.usda.gov/Services/docs.htm?docid=8955
		case CLAY:
			*theta_r=0.098;
			*theta_s=0.459;
			*n=1.253;
			*alpha=1.496;
			*ksat=0.14757/(24.*60.*60.);
			MaximumPoreSize=0.005;
			break;

		case CLAYLOAM:
			*theta_r=0.079;
			*theta_s=0.442;
			*n=1.416;
			*alpha=1.581;
			*ksat=0.0818/(24.*60.*60.);
			MaximumPoreSize=0.005;
			break;

		case LOAM:
			*theta_r=0.061;
			*theta_s=0.399;
			*alpha=1.11;
			*n=1.47;
			*ksat=0.02947/(24.*60.*60.);
			MaximumPoreSize=0.005;
			break;

		case LOAMYSAND:
			*theta_r=0.049;
			*theta_s=0.39;
			*n=1.746;
			*alpha=3.475;
			*ksat=1.052/(24.*60.*60.);
			MaximumPoreSize=0.005;
			break;

		case SAND:
			*theta_r=0.053;
			*theta_s=0.375;
			*n=3.177;
			*alpha=3.524;
			*ksat=6.427/(24.*60.*60.);
			MaximumPoreSize=0.005;
			break;

		case SANDYCLAY:
			*theta_r=0.117;
			*theta_s=0.385;
			*n=1.208;
			*alpha=3.342;
			*ksat=0.1135/(24.*60.*60.);
			MaximumPoreSize=0.005;
			break;

		case SANDYCLAYLOAM:
			*theta_r=0.063;
			*theta_s=0.384;
			*n=1.330;
			*alpha=2.109;
			*ksat=0.1318/(24.*60.*60.);
			MaximumPoreSize=0.005;
			break;

		case SANDYLOAM:
			*theta_r=0.039;
			*theta_s=0.387;
			*n=1.4488;
			*alpha=2.667;
			*ksat=0.3828/(24.*60.*60.);
			MaximumPoreSize=0.005;
			break;

		case SILT:
			*theta_r=0.050;
			*theta_s=0.489;
			*n=1.6788;
			*alpha=0.6577;
			*ksat=0.4375/(24.*60.*60.);
			MaximumPoreSize=0.005;
			break;

		case SILTYCLAY:
			*theta_r=0.111;
			*theta_s=0.481;
			*n=1.321;
			*alpha=1.622;
			*ksat=0.09616/(24.*60.*60.);
			MaximumPoreSize=0.005;
			break;

		case SILTYCLAYLOAM:
			*theta_r=0.090;
			*theta_s=0.482;
			*n=1.5205;
			*alpha=0.8395;
			*ksat=0.1112/(24.*60.*60.);
			MaximumPoreSize=0.005;
			break;

		case SILTLOAM:
			*theta_r=0.065;
			*theta_s=0.439;
			*n=1.6634;
			*alpha=0.5058;
			*ksat=0.1824/(24.*60.*60.);
			MaximumPoreSize=0.005;
			break;

		case WFJGRAVELSAND: //Gravel/sand
			*theta_r=0.01;
			*theta_s=0.35;
			*n=4.5;
			*alpha=3.5;
			*ksat=0.000003171; //Equal to 100 m/year, for clean sand and silty sand, according to: http://web.ead.anl.gov/resrad/datacoll/conuct.htm
			MaximumPoreSize=0.005;
			break;
	}

	*he=AirEntryPressureHead(MaximumPoreSize, 273.);
	*m=(*n-1.)/(*n);

	return;
}


/**
 * @brief Solve Richards Equation \n
 * Solve Richards Equation \n
 * @author Nander Wever
 * @param Xdata
 * @param Sdata
 */
void ReSolver1d::SolveRichardsEquation(SnowStation& Xdata, SurfaceFluxes& Sdata)
{
// Main references used to write this code, as a reference to understand the working of this code:
// - Celia, M, A., Bouloutas, E.T., Zabra, R.L. (1990) A general mass-conservative numerical solution for the unsaturated flow equation Water Resources Research (26:7), 1483-1496.
//   Describes the main part of the solver (Picard iteration of the mixed form of the Richards equation, see eq. 17 in that paper)
// - Bastos de Vasconcellos, C.A., Amorim, J.C.C. (2001) in Proceedings of COBEM 2001, Fluid Mechanics, Vol. 8, 139-148.
//   Is using Celia (1990), but I mention it here because it provides an excellent overview in matrix-notation of how to exactly implement the Picard-iteration and fill in the matrices and vectors, also in case fo Dirichlet or Neumann boundary conditions.
//   It will certainly help to understand why the matrices are filled the way they are (see Equations 18, 19, 20 and 21). Note: in Eq. 19, the coefficient for alpha is missing a minus-sign and in Eq. 20, ignore the alpha_1 coefficient, because that one does not exist.
// - Hari Prasad, K.S., Mohan Kumar, M.S and Sekhar, M. (2001) Modelling flow through unsaturated zones: Sensitivity to unsaturated soil properties. Sadhana vol. 26, no 6. 517–528. DOI: 10.1007/BF02703457
//   Provides an explicit description of how to implement Dirichlet boundary conditions.
// - Huang, K., Mohanty, B.P., van Genuchten, M.Th. (1996) A new convergence criterion for the modified Picard iteration method to solve variably saturated flow equation Journal of Hydrology (178), 69-91.
//   Describes a better convergence criterion, based on theta. It is implemented here, but it is not always behaving stable, likely when the convergence criterion is set too low.
// - McCord, J.T. (1991) Application of second-type boundaries in unsaturated flow modelling Water Resources Research (27:12), 3257-3260.
//   Describes the boundary conditions applied here (Neumann and Dirichlet).
// - Paniconi, C., Putti, M. (1994) A comparison of Picard and Newton iteration in the numerical solution of multidimensional variably saturated flow problems. Water Resources Research (30:12) 3357-3374.
//   Describes the time stepping mechanism implemented here.
// - Schaap, M.G. and van Genuchten, M.Th. (2006) A modified Mualem-van Genuchten formulation for improved description of the hydraulic conductivity near saturation. Vadoze Zone Journal (5) 27-34,
//   Describes the small changes made to the van Genuchten parameterization of soil, to better deal with saturated conditions (theta_m and h_s). This, I deactivated, as a better solution seems to be:
// - Ippisch, O., Vogel, H.-J. and Bastian, P. (2006) Validity limits for the van Genuchten-Mualem model and implications for parameter estimation and numerical simulation. Adv. in Water Res. (29), 1780-1789.
//   Describes inconsistenties in near saturation and suggest using an air entry pressure he.
// - Zhang, X., Bengough, A.G., Crawford, J.W. and Young, I.M. (2002) Efficient methods for solving water flow in variably saturated soils under prescribed flux infiltration. Journal of Hydrology (260) 75-87.
//   Describes an more efficient way of treating the gravity term, and describes a way of estimating the hydraulic conductivity at the nodes.
// - Li, C.W. (1993) A Simplified Newton Iteration Method With Linear Finite Elements for Transient Unsaturated Flow. Water Resources Research (29:4) 965-971.
//   Describes how the hydraulic conductivity at the interface nodes can be approximated by integration.
// - Yamaguchi, S., Katsushima, T., Sato, A. and Kumakura, T. (2010) Water retention curve of snow with different grain sizes. Cold Regions Science and Technology (64:2) 87-93.
//   Describes van Genuchten parameters (alpha and n) for snow, based on measurements. Experiments for around 550 kg/m^3 dense snow.
// - Yamaguchi, S., Watanabe, K., Katsushima, T., Sato, A., Kumakura, T. (2012) Dependence of the water retention curve of snow on snow characteristics. Annals of Glaciology 53(61). doi: 10.3189/2012AoG61A001
//   Update of the Yamaguchi (2010) paper: describes van Genuchten parameters (alpha and n) for snow, based on measurements. Experiments for a range of snow densities.
// - Hirashima, H., Yamaguchi, S., Sato, A., and Lehning, M. (2010) Numerical modeling of liquid water movement through layered snow based on new measurements of the water retention curve. Cold Regions Science and Technology (64:2), 94-103.
//   Describes van Genuchten parameters for snow, based on Yamaguchi (2010).
// - Rathfelder, K and Abriola, L. (1994) Mass conservative numerical solutions of the head-based Richards equation. Water Resources Research (30:9) 2579-2586.
//   Describes an implementation of variable grid spacing for solving Richards Equation in 1D.


// PROBLEM SOLVER GUIDE:
// KNOWN ISSUES:
//	- When using Richars-Equation, the new energy conservative PhaseChange-schemes may cause snow temperatures to be above 273.15K. As long as they are transient, it should not considered
//        to be a problem. Future optimization here may be possible. It's likely related to the fact that when solving Richards-Equation, basically every snow layer has some amount of water in it,
//        albeit very tiny. But this causes some difficulties in determining whether snow is wet or dry, so whether the nodes are at melting temperature.
//      - In case of floating point exceptions: ReSolver1d has some problems when (in CMake) DEBUG_ARITHM is set to ON. You can safely set it to OFF, as the code detects for
//        illegal operations itself and takes appropriate measures, like choosing another solver or reducing the time step.
//	- In case of non-convergence of the solver: Numerical problems were found when the SNOWPACK time step is larger than 15 minutes. For example caused by the settling routine,
//	  which is based on 15 minute time steps. So before digging further in the problem, make sure you run SNOWPACK with 15 minute time steps.
//      - Evaporation from soil in dry limit. This cause numerical troubles, but it is also not realistic to force a certain amount of evaporation from a near-dry soil (the water
//	  is just not there!). Set LIMITEDFLUXEVAPORATION or LIMITEDFLUX as top boundary condition to be safe.
//      - Infiltration in soil in wet limit. This can cause numerical trouble, but it is also not realistic. You cannot put more water in the domain then there is room for.
//	  So for example: never use 10 cm of soil with DIRICHLET lower boundary condition and NEUMANN on top. Set LIMITEDFLUXINFILTRATION or LIMITEDFLUX as lower boundary condition to be safe.
//
// SOME DEEPER LYING PROBLEMS:
// -  In case of floating point exceptions with DEBUG_ARITHM set to OFF:
//      - either pressure head is blowing up to extremely high values, before rewinds can be done. This can happen if you change MAX_ALLOWED_DELTA_H, but also make sure that the time step is still larger than the minimum allowed time step.
//        If no rewinds can be done anymore, because the time step is already too small, the solver is just trying until it is killed.
// -  In case of other exceptions, see the error message on screen. Messages from the solver: if info-value is positive, the number indicates quite often the element that is having illegal value, like nan, or inf. This can be a result from very small time steps (check dt), or else it is a bug.
// -  Problems can be caused when the state of the snowcover changes rapidly between calls to ReSolver1d. For example, it was found that when SNOWPACK was run with 60 minute
//    time steps, the settling was so fast that elements collapsed, producing saturated snow layers. Then ReSolver1d was not able to solve the equation anymore with reasonable time steps.
// -  A lot of problems arise from non convergence in the solver and very small time steps. A common reason for this is filling of the model domain. Clear examples:
//    - 10cm soil with saturated lower boundary condition. The soil will saturate quickly. No melt water can infilitrate the soil anymore, but starts ponding. In reality, such cases would lead to a water layer, or overland flow.
//    - Strong infiltration rates (like high precipitation rates)
//    - High values in the source term
//    In these cases the model just cannot do anything with all the water...
// -  The model is still running and producing sensible results, but is having a very small time step.
//      First of all, I found that only numerous subsequent timesteps of <1s. should not be accepted and there is something "wrong" with the code. Else, just leave the code running. Note: you cannot really say wrong, as it is still finding converging solutions.
//      Try to understand why the time step rewinds are done, maybe by even enabling the Level3 output. And using less to search the output with "less" and search for "rewind".
//	  - Is there a massbalance problem? Then probably the fluxes at the top and/or bottom are not correctly calculated, or there is a real bug.
//	    Massbalance errors also arise when k_ip12(i-1) != k_im12(i)! The nodal values should always be the same for both the upper and lower node!
//	  - Is there a very strong gradient in pressure head, for example at the new snow layer? What is the value for h_d, is it very small? Then maybe limit the range over which the Van Genuchten parameters can vary (limiting grain size for example for snow).
//
// TODO IN FUTURE DEVELOPMENT
// -  Implement a strategy what to do with the rejected infilitrating water in case of LIMITEDFLUX and LIMITEDFLUXINFILTRATION. Either built-up a water layer (theta[WATER]==1) on top (real ponding),
//    or write it out in a kind of overland flow variable.

	//Initializations
	enum RunCases{UNIFORMSOIL, IMISDEFAULT, WFJ, CDP, SNOFILE};

	//
	// BEGIN OF SETTINGS
	//
	const RunCases runcase = SNOFILE;					//Defines what the soil looks like. Recommended: SNOFILE, soil type based on grain size in sno file.
	const BoundaryConditions TopBC = LIMITEDFLUX;				//Bottom boundary condition (recommended choice is LIMITEDFLUX, so too much evaporation from dry soil or snow or too much infilitration in wet soil is prohibited).
		//In case you select one of the LIMITEDFLUX options, specify whether these are only for soil, for snow or for both:
		const bool LIMITEDFLUXEVAPORATION_soil=true;
		const bool LIMITEDFLUXEVAPORATION_snow=true;
		const bool LIMITEDFLUXINFILTRATION_soil=true;
		const bool LIMITEDFLUXINFILTRATION_snow=true;
		const bool LIMITEDFLUXINFILTRATION_snowsoil=true;		//This switch allows to limit the infiltration flux from snow into soil, when the snowpack is solved with the Bucket or NIED water transport scheme.
	const bool AllowSoilFreezing=true;					//true: soil may freeze. false: all ice will be removed (if any ice present) and no ice will form.
	const bool ApplyIceImpedance=false;					//Apply impedance on hydraulic conductivity in case of soil freezing. See: Zhao et al. (1997) and Hansson et al. (2004)  [Dall'Amicao, 2011].
	const VanGenuchten_ModelTypesSnow VGModelTypeSnow=YAMAGUCHI2012;	//(Recommended: YAMAGUCHI2012) Set a VanGenuchten model for snow (relates pressure head to theta and vice versa)
	const bool alpine3d=false;						//Flag for alpine3d simulations. Note: this flag is not necessary to set, but it will enforce some precautions to provide extra numerical stability (at the cost of small mass balance violations).

	//Setting some program flow variables
	const bool SafeMode=false;			//Enable safemode only when necessary, for example in operational runs or Alpine3D simulations. It rescues simulations that do not converge, at the cost of violating mass balance.

	const bool WriteOutNumerics_Level0=false;	//true: after time step, some summarizing numerics information is printed
	const bool WriteOutNumerics_Level1=false;	//true: per iteration, some basic numerics information is printed (not so much output, but still is useful for debugging)
	const bool WriteOutNumerics_Level2=false;	//true: only initial conditions and subsequent convergence information is printed during each time step (it's quite some output, but helps debugging)
	const bool WriteOutNumerics_Level3=false;	//true: per iteration, highly detailed layer and solver is printed (it's really a lot of output, but helps debugging the most difficult problems)

	const bool AllowDrySnowLayers=false;		//true: snow layers can be dry (theta=0.) false: snow layers will always be at least theta=theta_r. The necessary water is created by melting ice.
	const bool AllowDrySoilLayers=false;		//true: soil layers can be dry (theta=0.) false: soil layers will always be at least theta=theta_r. The necessary water is created out of nothing!
							//  Will destroy mass balance! Setting this flag true for soil layers is not recommended as it will not be really physical plausible.

//Set the defaults based on whether CLAPACK is available
#ifdef CLAPACK
	SOLVERS PreferredSolver=DGTSV;			//Choose either DGESVD, DGTSV or TDMA.
	const bool AllowSwitchSolver=true;		//If true: solver DGESVD will be tried when DGTSV fails. There is a trade-off here between matrix inversion and smaller time steps. In case of many layers (>100), DGESVD can become very slow, so in case DGTSV does not find a solution, it may be more efficient to
	//take a smaller time step than to do full matrix inversion.
#else
	SOLVERS PreferredSolver=TDMA;			//Without CLAPACK, only TDMA is available.
#endif
							// DGTSV : (recommended) This function does matrix inversion giving the knowledge matrix A is a tridiagonal matrix (fast). Does partial pivoting to stabelize numerical solutions. But partial pivoting is not always enough.
							// DGESVD: (recommended only when stability issues are encounterd) This function does full matrix inversion (slow, but most reliable and maybe useful for finding 2D solutions one day).
							// TDMA  : (recommended only when libraries BLAS and LAPACK are not available) This function does matrix inversion giving the knowledge matrix A is a tridiagonal matrix (really fast). Does no (partial) pivoting at all, so big risks of numerical troubles.

	SOLVERS ActiveSolver=PreferredSolver;		//Set the ActiveSolver to the PreferredSolver. This is because the code tries to prevent "difficult" matrices to be solved by the DGTSV or TDMA algorithm, so we should be able to switch temporarily to another solver.

	//Set parameterization for hydraulic conductivity
	const K_Parameterizations K_PARAM=CALONNE;		// Implemented choices: SHIMIZU, CALONNE, based on Shimizu (1970) and Calonne (2012).
	//Set how the hydraulic conductivity at the interface nodes should be calculated.
	const K_AverageTypes K_AVERAGETYPE=ARITHMETICMEAN;	// Implemented choices: ARITHMETICMEAN (recommended), HARMONICMEAN, GEOMETRICMEAN, MINIMUMVALUE, UPSTREAM





	//
	// END OF SETTINGS
	// WARNING: Below this line, changes to initializations are likely to break the code!
	//

	//Setting convergence criteria and numerical limits
	const double REQUIRED_ACCURACY_H=1E-3;		//Required accuracy for the Richard solver: this is for the delta h convergence criterion
	const double REQUIRED_ACCURACY_THETA=1E-5;	//Required accuracy for the Richard solver: this is for the delta theta convergence criterion. It is recommended to adjust PhaseChange::RE_theta_r in PhaseChanges.cc in case this value is changed.
							//Huang et al. (1996) proposes 0.0001 here (=1E-4). 1E-4 causes some mass balance problems. Therefore, it is set to 1E-5.
	const double convergencecriterionthreshold=0.99;//Based on this value of theta_dim, either theta-based convergence is chosen, or h-based. Note we need to make this destinction, beacuse theta-based does not work close to saturation or with ponding.
	const double MAX_ALLOWED_DELTA_H=1E32;		//Set an upper threshold for the delta_h[i] that is allowed. The idea is that when delta_h for an iteration is too large, we have a too large time step and a rewind is necessary.
	const int INCR_ITER=5;				//Number of iterations for the Richard solver after which time step is increased.
	const int DECR_ITER=10;				//Number of iterations for the Richard solver after which time step is decreased.
	const int MAX_ITER=15;				//Maximum number of iterations for the Richard solver.
	const double MIN_VAL_TIMESTEP=1E-12;		//Minimum time step allowed in Richards solver. Don't set this too low (let's say 1E-40), becuase the calculations are then done at the limits of the floating point precision.
	const double MAX_VAL_TIMESTEP=900.;		//Maximum time step allowed in Richards solver.
	const double MAX_VAL_TIMESTEP_FOR_SNOW=900.;	//Maximum time step allowed in Richards solver when there are snow layers in the domain.
	const double BS_MAX_ITER=5000;			//Maximum allowed number of iterations in the soil-freezing algorithm.
	const double SF_epsilon=1E-4;			//Required accuracy for the root finding algorithm when solving soil freezing/thawing.

	//Initializing and defining Richards solver time domain
	const double snowpack_dt = sn_dt;		//Time step of SNOWPACK (in seconds)
	double dt=10.;					//Set the initial time step for the Richard solver (in seconds). This time step should be smaller or equal to the SNOWPACK time step.
	bool boolFirstFunctionCall;			//true: first execution of this function, false: not the first execution of this function
	if (Xdata.ReSolver_dt>0.) {			//Retrieve last dt used in last performed time step. Note Xdata.SoliNumerics_dt<0 when initialized
		boolFirstFunctionCall=false;		//Subsequent call to ReSolver1d
		dt=Xdata.ReSolver_dt;			//Set time step to last used time step
 	} else {					//else it is the first time this function is called.
		boolFirstFunctionCall=true;		//Set this flag to true, so we know that no previous pressure head information is available, and we can only work with theta.
	}
	double TimeAdvance=0.;				//Time advance of the Richards solver


	//Initializing and defining Richards solver space domain
	const size_t nN=Xdata.getNumberOfNodes();	//Number of nodes
	const size_t nE=nN-1;				//Number of layers
	const double cos_sl = Xdata.cos_sl;		//Slope cosinus, giving cos_sl=1 for flat field.
	vector<ElementData>& EMS = Xdata.Edata;		//Create reference to SNOWPACK elements.
	vector<NodeData>& NDS = Xdata.Ndata;		//Create reference to SNOWPACK nodes.
	std::vector<double> dz(nE, 0.);			//Layer height (in meters)
	std::vector<double> z(nE, 0.);			//Height above the surface (so -1 is 1m below surface)
	std::vector<double> dz_up(nE, 0.);		//Distance to upper node (in meters)
	std::vector<double> dz_down(nE, 0.);		//Distance to lower node (in meters)
	std::vector<double> dz_(nE, 0.);		//Layer distance for the finite differences, see Rathfelder (2004).
	int uppernode=-1;				//Upper node of Richards solver domain
	int lowernode=-1;				//Lower node of Richards solver domain
	std::vector<int>SnowpackElement(nE,0);		//Dictionary between snowpack domain and Richards solver domain. SnowpackElement[j]=i means layer j in Richards solver is layer i in snowpack domain.
							//Then, using EMS[SnowpackElement[j]], we can refer to the SNOWPACK domain from the Richards solver domain.
	int toplayer;					//highest layer (top of snowpack, or top of soil in case of no soil)
	const int nsoillayers_snowpack=int(Xdata.SoilNode);	//where does the soil start? Note, when toplayer is set to nsoillayers_snowpack, only soil is treated with Richards equation. HACK, TODO: remove type inconstency in comparison
							//Note here that Xdata.SoilNode denotes first element as 0, so Xdata.SoilNode=4 denotes 4 soil layers.
	int nsoillayers_richardssolver=0;
	if(iwatertransportmodel_snow != RICHARDSEQUATION) {	//RE only for soil
		toplayer=nsoillayers_snowpack;			//toplayer=nE: all layers are treated by richards equation, toplayer=nsoillayers_snowpack: only soil.
	} else {						//RE for both snow and soil
		toplayer=int(nE);				//toplayer=nE: all layers are treated by richards equation, toplayer=nsoillayers_snowpack: only soil. HACK, TODO: remove type inconstency in comparison
	}
	if(toplayer==0) return;				//Nothing to do here!

	//Initializations of the convergence criteria
	int trigger_layer_accuracy=-1;			//At which layer the accuracy was not reached.
	double track_accuracy_h=0.;			//This variable tracks the accuracy of convergence for all h-convergence based layers.
	double track_accuracy_theta=0.;			//This variable tracks the accuracy of convergence for all theta-convergence based layers.
	double max_delta_h=0.;				//Tracks max_delta_h, to determine if our time step is too large. Note: this is different from checking the convergence criterion. This is just to check the time step. If a too large time step is used, big values of delta_h may arise, which cause pressure head to hit the singularities for dry and wet soil, and causes problems with the power functions in the Von Genuchten-Mualem model.
	int track_trigger_layer_accuracy=-1;		//This variable tracks the layer were the accuracy is smallest (this means: most difficulty in converging).
	bool boolConvergence=false;			//true: convergence is reached, false: convergence not reached
	double mass1=0, mass2=0, massbalanceerror=0.;	//Mass balance check variables.
	const double maxallowedmassbalanceerror=1E-10;	//This value is carefully chosen. It should be considered together with REQUIRED_ACCURACY_THETA and REQUIRED_ACCURACY_H
	double massbalanceerror_sum=0.;			//Sum of mass balance error over snowpack time step.


	//Initializations for summarizing statistics and some supporting variables, like indices, counters, etc.
	double accuracy=0.;				//Keeps track of reached accuracy.
	int niter=0;					//Counts iterations within one time step of the Richards solver
	int niter_snowpack_dt=0;			//Counts iterations within one time step of the SNOWPACK time domain
	int niter_nrewinds=0;				//Counts number of rewinds (i.e. a solution was not found and it is tried again with a smaller time step)
	int niter_seqrewinds=0;				//Counts number of sequential rewinds. We then decrease the time step more, when we encounter sequential rewinds.
	int seq_safemode=0;				//Counts the number of sequential SafeMode actions
	//Numerical performance statistics
	double stats_min_dt=MAX_VAL_TIMESTEP;		//Minimum RE time step in SNOWPACK time step, initialized in a way that the comparison will always work.
	double stats_max_dt=0.;				//Maximum RE time step in SNOWPACK time step, initialized in a way that the comparison will always work.
	int stats_nrewinds=0;				//Number of rewinds over the SNOWPACK time step.
	int stats_niters=0;				//Number of iterations over the SNOWPACK time step, excluding the ones before a rewind.
	int stats_nsteps=0;				//Number of time steps in the RE solver over the SNOWPACK time step.
	size_t bs_stats_totiter=0;			//Soil freezing/thawing solver: total number of iterations over all layers over the SNOWPACK time step,
	size_t bs_stats_maxiter=0;			//Soil freezing/thawing solver: maximum number of iterations in a certain layers over the SNOWPACK time step.
	//Counters, etc.
	int i, j, k;					//Counters for layers
	const size_t nmemstates=1;			//Number of memory states, used to store changes of delta_h between iterations. Currently not used, but possible use is to check if delta_h is blowing up.
	int memstate=0;					//Keeping track of the current memory index
	double h_d=0.;					//Lower limit for pressure head: definition of "completely dry". This value will determined later on.


	//Initializing and declaring boundary conditions and flux variables
	BoundaryConditions aTopBC;			//Actual applied top boundary condition (can only be either Dirichlet or Neumann, as all the others can be translated in an application of either one of those two.)
	BoundaryConditions aBottomBC;			//Actual applied bottom boundary condition (can only be either Dirichlet or Neumann, as all the others can be translated in an application of either one of those two.)
	double htop=0., TopFluxRate=0.;			//Dirichlet (constant head) and Neumann (constant flux) upper boundary values respectively.
	double h_d_uppernode=0.;			//Used for LIMITEDFLUXEVAPORATION boundary condition.
	double hbottom=0., BottomFluxRate=0.;		//Dirichlet (constant head) and Neumann (constant flux) lower boundary values respectively.
	double actualtopflux=0;				//Stores the actual applied flux through top (positive is inflow).
	double actualtopfluxcheck=0.;			//Stores the water change in the top element + outflow to lower layers to derive the input flux at the surface.
	double refusedtopflux=0;			//Stores the difference in flux that was requested, but could not be applied
	double actualbottomflux=0;			//Stores the actual flux through the bottom (positive is outflow).
	double snowsoilinterfaceflux1=0.;		//Stores the actual flux through the soil-snow interface (positive is flow into soil).
	double snowsoilinterfaceflux2=0.;		//Stores the actual flux through the soil-snow interface (positive is flow into soil).
	double snowsoilinterfaceflux_before=0.;		//Stores the actual flux through the soil-snow interface at the beginning of a time step (positive is flow into soil).
	double snowsoilinterfaceflux_after=0.;		//Stores the actual flux through the soil-snow interface at the end of a time step (positive is flow into soil).
	double totalsourcetermflux=0.;			//Stores the total applied source term flux (it's a kind of boundary flux, but then in the middle of the domain).

	//Declare all numerical arrays and matrices:
	std::vector< std::vector<double> > delta_h(nmemstates, std::vector<double> (nE,0.));	//Change in pressure head per iteration
	std::vector<double> delta_h_dt(nE, 0.);		//Change in pressure head per time step.
	std::vector<double> delta_theta(nE, 0.);	//Change in volumetric water content per iteration
	std::vector<double> delta_theta_dt(nE, 0.);	//Change in volumetric water content per time step.
	std::vector<double> delta_theta_i(nE, 0.);	//Change in volumetric ice content per iteration
	std::vector<double> delta_theta_i_dt(nE, 0.);	//Change in volumetric ice content per time step.
	std::vector<double> delta_Te(nE, 0.);		//Change in element temperature per time step due to soil freezing/thawing.
	std::vector<double> delta_Te_i(nE, 0.);		//Change in element temperature per iteration time step due to soil freezing/thawing.
	std::vector<double> delta_Te_adv(nE, 0.);	//Change in element temperature per time step due to heat advection by the water flow.
	std::vector<double> delta_Te_adv_i(nE, 0.);	//Change in element temperature per iteration time step due to heat advection by the water flow.

	//std::vector<std::vector<double> > a(nE, std::vector<double> (nE, 0));	//Left hand side matrix. Note, we write immediately to ainv! But this is kept in to understand the original code.
	std::vector<double> ainv(nE*nE, 0.);			//Inverse of A, written down as a 1D array instead of a 2D array, with the translation: a[i][j]=ainv[i*nsoillayers_snowpack+j]
	std::vector<double> ad(nE, 0.);				//The diagonal of matrix A, used for DGTSV
	std::vector<double> adu(nE, 0.);			//The upper second diagonal of matrix A, used for DGTSV
	std::vector<double> adl(nE, 0.);			//The lower second diagonal of matrix A, used for DGTSV

	std::vector<double> k_np1_m_ip12(nE, 0.);		//Hydraulic conductivity at the upper interface node
	std::vector<double> k_np1_m_im12(nE, 0.);		//Hydraulic conductivity at the lower interface node
	std::vector<double> h_np1_m(nE, 0.);			//Pressure head at beginning of an iteration.
	std::vector<double> h_n(nE, 0.);			//Pressure head at beginning of time step dt. Used to determine delta_h_dt, to better forecast value for next time step.
	std::vector<double> s(nE, 0.);				//Source/sink in terms of theta [m^3/m^3/s].
	std::vector<double> C(nE, 0.);				//Water capacity function. Specific moisture capacity (dtheta/dh), see Celia et al., (1990).
	std::vector<double> ksat(nE, 0.);			//Soil property. Saturation hydraulic conductivity.
	std::vector<double> K(nE, 0.);				//Hydraulic conductivity function
	std::vector<double> impedance(nE, 0.);			//Impedance factor due to ice formation in matrix (see Dall'Amico, 2011);
	std::vector<double> Se(nE, 0.);				//Effective saturation, sometimes called dimensionless volumetric water content.
	std::vector<double> term_up(nE, 0.);			//Variable to support construction of the R.H.S. (R_mpfd in Celia et al., 1990).
	std::vector<double> term_down(nE, 0.);			//Variable to support construction of the R.H.S. (R_mpfd in Celia et al., 1990).
	std::vector<double> r_mpfd(nE, 0.);			//R_mpfd (see Celia et al, 1990).
	std::vector<double> r_mpfd2(nE, 0.);			//Copy of R_mpfd, used for DGTSV. Note: R_mpfd2 is overwritten by DGTSV, so we need a copy.
	std::vector<double> deficit_vector(nE, 0.);		//Deficit vector of solution
	double deficit_vector_norm=0.;				//2nd norm of deficit vector
	std::vector<double> h_np1_mp1(nE, 0.);			//Pressure head for the solution time step in the next iteration
	std::vector<double> theta_np1_m(nE, 0.);		//Theta for the solution time step in the current iteration.
	std::vector<double> theta_np1_mp1(nE, 0.);		//Theta for the solution time step in the next iteration.
	std::vector<double> theta_n(nE, 0.);			//Theta at the current time step.
	std::vector<double> theta_d(nE, 0.);			//There is a singularity for dry soils, at theta=theta_r. There h -> OO. So we limit this. We define a pressure head that we consider "dry soil" (h_d) and then we calculate what theta belongs to this h_d.
	std::vector<double> alpha(nE, 0.);			//Soil property in Van Genuchten model. [m^-1]
	std::vector<double> n(nE, 0.);				//Soil property in Van Genuchten model.
	std::vector<double> m(nE, 0.);				//Soil property in Van Genuchten model.
	std::vector<double> h_e(nE, 0.);			//Soil property, air entry pressure, see Ippisch (2006) for details.
	std::vector<double> theta_r(nE, 0.);			//Soil property, residual water content.
	std::vector<double> theta_s(nE, 0.);			//Soil property, saturation water content.
	std::vector<double> Sc(nE, 0.);				//Saturation at cut-off point h_e (see Ippisch et al (2006)).
	std::vector<double> theta_i_n(nE, 0.);			//Soil state, ice content at the beginning of the time step. Volumetric water content and NOT liquid water equivalent!
	std::vector<double> theta_i_np1_m(nE, 0.);		//Soil state, ice content at the beginning of the current iteration. Volumetric water content and NOT liquid water equivalent!
	std::vector<double> theta_i_np1_mp1(nE, 0.);		//Soil state, ice content at the next iteration. Volumetric water content and NOT liquid water equivalent!

	std::vector<bool> activelayer(nE, true);		//true: layer is active participating in matrix flow (Richards equation). false: layer is inactive (too dry, or ice layer)
	std::vector<double> dT(nE, 0.);				//Stores the energy needed to create theta_r from the ice matrix.
	std::vector<double> snowpackBACKUPTHETAICE(nE, 0.);	//Backup array for the initial SNOWPACK theta ice
	std::vector<double> snowpackBACKUPTHETAWATER(nE, 0.);	//Backup array for the initial SNOWPACK theta water
	std::vector<double> wateroverflow(nE, 0.);		//Array for all the water that is >theta_s (m^3/m^3)]. This water is just thrown away in the model and is a leak in the mass balance.

        //For soil freezing/thawing
	std::vector<double> T_melt(nE, 0.);			//Contains the freezing point depression due to unsaturated conditions (K)
	const double T_0=Constants::freezing_tk;		//Freezing temperature of water at atmospheric pressure (K)
	const double delF=Constants::lh_fusion;			//Heat associated with freezing


	//Prevent buffering on the stdout when we write debugging output. In case of exceptions (program crashes), we don't loose any output which is still in the buffer and we can better track what went wrong.
	if(WriteOutNumerics_Level1==true || WriteOutNumerics_Level2==true || WriteOutNumerics_Level3==true) {
		setvbuf(stdout, NULL, _IONBF, 0);
	}


#ifdef DEBUG_ARITHM
	if(boolFirstFunctionCall==true) {
		// Warn users for compilation with DEBUG_ARITHM = ON. The solver uses isnan and isinf to check if the time step is too large.
		prn_msg( __FILE__, __LINE__, "wrn", Date(), "SNOWPACK has been compiled with DEBUG_ARITHM set to ON. This will likely result in a \"Floating point exception\" when using Richards equation solver. It is strongly recommended to set DEBUG_ARITHM to OFF!");
	}
#endif

	//Backup SNOWPACK state
	for (i = 0; i<toplayer; i++) {		//Cycle through all SNOWPACK layers
		//Do the basic check if there is not too much ice.
		if(EMS[i].theta[ICE]>1.) {
			std::cout << "WARNING: very strange, theta[ICE]>1 at layer " << i << "/" << nE << " (" << EMS[i].theta[ICE] << ")\n";
			EMS[i].theta[ICE]=1.;
		}
		//Do the basic check of the sum of element contents.
		const double sum=EMS[i].theta[AIR] + EMS[i].theta[WATER] + EMS[i].theta[ICE] + EMS[i].theta[SOIL];
		if((sum>1.+Constants::eps || sum<1.-Constants::eps) && (boolFirstFunctionCall!=true)) {
			std::cout << "WARNING: very strange, sum of element contents != 1. (but " << sum << ") at layer " << i << "/" << nE << ". Values scaled.\n";
			//Note: we do not scale theta[SOIL].
			const double correction_factor=(EMS[i].theta[AIR] + EMS[i].theta[WATER] + EMS[i].theta[ICE])/(1.-EMS[i].theta[SOIL]);
			EMS[i].theta[AIR]/=correction_factor;
			wateroverflow[i]+=(EMS[i].theta[ICE]-(EMS[i].theta[ICE]/correction_factor))*(Constants::density_ice/Constants::density_water);	//We just throw away the ice, without considering melting it.
			EMS[i].theta[ICE]/=correction_factor;
			wateroverflow[i]+=EMS[i].theta[WATER]-(EMS[i].theta[WATER]/correction_factor);
			EMS[i].theta[WATER]/=correction_factor;
		} else {
			if(boolFirstFunctionCall==true) {
				EMS[i].theta[AIR]=1.-EMS[i].theta[SOIL]-EMS[i].theta[ICE]-EMS[i].theta[WATER];
			}
		}

		if (WriteOutNumerics_Level2==true) 
			std::cout << "RECEIVING at layer " << i << ": sum=" << sum << std::fixed << std::setprecision(15) <<  " air=" << EMS[i].theta[AIR] << " ice=" << EMS[i].theta[ICE] << " soil=" << EMS[i].theta[SOIL] << " water=" << EMS[i].theta[WATER] << " Te=" << EMS[i].Te << "\n"<< std::setprecision(6) ;

		//In case we don't want to allow soil to freeze, melt all ice that is there:
		if(AllowSoilFreezing==false && i<nsoillayers_snowpack && EMS[i].theta[ICE]>0.) {	//If we are in soil, and have ice, but don't allow soil freezing, melt all the ice.
			if(i>0) {				//Because we now work with Dirichlet BC at the lower boundary, we cannot increase water content. So then, increase air content.
				EMS[i].theta[WATER]+=EMS[i].theta[ICE]*(Constants::density_ice/Constants::density_water);
			} else {
				EMS[i].theta[AIR]+=EMS[i].theta[ICE];
			}
			const double deltaT=(-1.*EMS[i].theta[ICE]) / ((EMS[i].c[TEMPERATURE] * EMS[i].Rho) / ( Constants::density_ice * Constants::lh_fusion ));
			EMS[i].Te+=deltaT;

			if(i==int(nE)-1 && i>=0) {		//HACK, TODO: remove type inconstency in comparison
				NDS[i+1].T+=deltaT;
				NDS[i].T+=deltaT;
			}

			EMS[i].Qmf += (-1.*EMS[i].theta[ICE] * Constants::density_ice * Constants::lh_fusion) / snowpack_dt;	// Units: [W m-3]
			EMS[i].theta[ICE]=0.;
			//And now update state properties.
			EMS[i].Rho = (EMS[i].theta[ICE] * Constants::density_ice) + (EMS[i].theta[WATER] * Constants::density_water) + (EMS[i].theta[SOIL] * EMS[i].soil[SOIL_RHO]);
			EMS[i].M=EMS[i].L*EMS[i].Rho;
		}

		//Make backup of incoming values for theta[ICE] and theta[WATER]. This is used in case we allow dry snow layers  *AND*  MIN_VAL_THETA_SNOWPACK > 0., to determine original water content and how much water was added to the domain.
		snowpackBACKUPTHETAICE[i]=EMS[i].theta[ICE];
		snowpackBACKUPTHETAWATER[i]=EMS[i].theta[WATER];
	}


	//Domain initialization (this needs to be done every time step, as snowpack layers will settle and thereby change height.
	//i=0 is bottom layer, and toplayer-1 is top layer.
	double heightshift=0.;			//heightshift can be used to shift the vertical domain up and down.
	double totalheight=0.;			//tracking the total height of the column
	j=0;					//Reset Richards-solver domain layer count
	for (i = 0; i<toplayer; i++) {		//Cycle over all SNOWPACK layers
		totalheight+=EMS[i].L;
		if(i>nsoillayers_snowpack-1) { 					//We are in snow, and create n sublayers per snow layer
			const int nsublayers=1;					//Specify number of sublayers
			for(k=0; k<nsublayers; k++) {				//Cycle through sublayers
				dz[j]=(EMS[i].L/nsublayers);			//Divide snowpack layer into sublayers
				if(j==0) {	//Lowest element
					z[j]=.5*dz[j]-heightshift;
				} else {
					z[j]=z[j-1]+(dz[j-1]/2.+(dz[j])/2.)-heightshift;
				}
				SnowpackElement[j]=i;				//Fill dictionary
				j++;						//Increase richards solver domain layer count
			}
		} else {
			const int nsublayers=1;					//Specify number of sublayers
			for(k=0; k<nsublayers; k++) {
				dz[j]=EMS[i].L/nsublayers;
				if(j==0) {	//Lowest element
					z[j]=.5*dz[j]-heightshift;
				} else {
					z[j]=z[j-1]+(dz[j-1]/2.+(dz[j])/2.)-heightshift;
				}
				SnowpackElement[j]=i;				//Fill dictionary
				nsoillayers_richardssolver++;			//Count the number of soillayers in richards solver domain
				j++;						//Increase richards solver domain layer count
			}
		}
	}
	uppernode=j-1;	//j increased 1 too much
	lowernode=0;	//lower node is just the bottom


	//Additional domain initialization: determine grid cell sizes, and node distances.
	//See for additional details on finite differences scheme with varying grid cell size: Rathfelder (1994).
	double tmpheight1=0., tmpheight2=0.;
	for (j=lowernode; j<=uppernode; j++) {
		//Distance to lower node
		if(j!=lowernode) {
			dz_down[j]=z[j]-z[j-1];
		}
		tmpheight1+=dz_down[j];
		//Distance to upper node
		if(j!=uppernode) {
			dz_up[j]=z[j+1]-z[j];
		}
		tmpheight2+=dz_up[j];
		//Mean distance
		//dz_[j]=0.5*(dz_down[j]+dz_up[j]);	//This is the definition of dz_ by Rathfelder (2004). However, it does not work, results in mass balance errors.
		dz_[j]=dz[j];				//This works.
		if(WriteOutNumerics_Level3==true) 
			std::cout << "DOMAIN: node " << j << " -- z=" << z[j] << std::fixed << std::setprecision(15) << " dz=" << dz[j] << " dz_up=" << dz_up[j] << " dz_down=" << dz_down[j] << " dz_=" << dz_[j] << "\n" << std::setprecision(6);
	}
	dz_down[lowernode]=totalheight-tmpheight1;
	dz_up[uppernode]=totalheight-tmpheight2;
	if(WriteOutNumerics_Level2==true) std::cout << "SLOPE: " << std::fixed << std::setprecision(15) << cos_sl << "\n" << std::setprecision(6);


	//Now set van Genuchten parameter for each layer
	h_d=0.;									//Set definition of pressure head of completely dry to zero, we will determine it in the next loop.
	double tmpheight=0.;
	for (i=uppernode; i >= 0; i--) {					//Go from top to bottom in Richard solver domain
		if ( SnowpackElement[i] >= nsoillayers_snowpack) {		//Snow, assuming that the use of sublayers (higher resolution) is only used in snow. TODO: this has to be rewritten more nicely!!
			const double max_allowed_ice=0.95;			//An ice pore space of 5% is a reasonable value: K. M. Golden et al. The Percolation Phase Transition in Sea Ice, Science 282, 2238 (1998), doi: 10.1126/science.282.5397.2238
			if(EMS[SnowpackElement[i]].theta[ICE]>max_allowed_ice) {
				//Pure ice layers are a problem for Richards equation (of course...), so we limit the volumetric ice content to 99 %.
				const double tmp_missing_theta=(EMS[SnowpackElement[i]].theta[ICE]-max_allowed_ice)*(Constants::density_ice/Constants::density_water);	//Not too dry (original)
				dT[SnowpackElement[i]]+=tmp_missing_theta*(Constants::density_water/Constants::density_ice) / ((EMS[SnowpackElement[i]].c[TEMPERATURE] * EMS[SnowpackElement[i]].Rho) / ( Constants::density_ice * Constants::lh_fusion ));
				std::cout << "[W] ReSolver1d.cc: ICE LAYER --> WATER CREATED (" << tmp_missing_theta << "): i=" << i << " --- dT=" << dT[SnowpackElement[i]] << " T=" << EMS[SnowpackElement[i]].Te << " theta[WATER]=" << EMS[SnowpackElement[i]].theta[WATER] << " theta[ICE]=" << EMS[SnowpackElement[i]].theta[ICE] << "\n";
				EMS[SnowpackElement[i]].theta[WATER]+=0.99*tmp_missing_theta;	//Here, we make a small mass balance error, but it should prevent fully saturated layers
				EMS[SnowpackElement[i]].theta[ICE]-=tmp_missing_theta*(Constants::density_water/Constants::density_ice);
				EMS[SnowpackElement[i]].theta[AIR]=1.-EMS[SnowpackElement[i]].theta[ICE]-EMS[SnowpackElement[i]].theta[WATER];
			}
			//Scaling theta_r between 0 and 0.02:
			const double TuningFactor=0.75;				//Tuning factor for scaling
			//Increase theta_r in case of wetting:
			theta_r[i]=MAX(0., MIN(0.02, MAX(EMS[SnowpackElement[i]].theta_r, TuningFactor*EMS[SnowpackElement[i]].theta[WATER])));
			//Decrease theta_r in case of refreezing:
			theta_r[i]=MAX(0., MIN(theta_r[i], EMS[SnowpackElement[i]].theta[WATER]-(REQUIRED_ACCURACY_THETA/10.)));

			theta_s[i]=(1. - EMS[SnowpackElement[i]].theta[ICE])*(Constants::density_ice/Constants::density_water);

			//Make ice layers inactive:
			if(theta_s[i]<Constants::eps2) {
				if(WriteOutNumerics_Level3==true) 
					std::cout << "WARNING: layer " << i << "/" << nE << "   theta_s= " << theta_s[i] << "   EMS[i].theta[ICE]=" << EMS[SnowpackElement[i]].theta[ICE] << "\n";
				std::cout << "WARNING: layer " << i << "/" << nE << " is ice layer! (theta_s[i]=" << theta_s[i] << "; EMS[i].theta[ICE]=" << EMS[SnowpackElement[i]].theta[ICE] << ")\n";
				activelayer[i]=false;
			} else {
				activelayer[i]=true;
			}
			//Make sure theta_r << theta_s
			if(theta_s[i]<theta_r[i]+0.01) {
				if(WriteOutNumerics_Level1==true) std::cout << "WARNING: layer " << i << "/" << nE << "  theta_s= " << theta_s[i] << "   theta_r=" << theta_r[i] << "\n";
				theta_r[i]=theta_s[i]/4.;
			}
			//Set air entry pressure
			h_e[i]=AirEntryPressureHead(0.005, 273.);

			//Note: rg is in mm, and it is the radius (confirmed by Charles, see DataClasses.h)
			const double tmprg=EMS[SnowpackElement[i]].rg;	//Backup original grain size value

			switch ( VGModelTypeSnow ) {	//Set Van Genuchten parameters for snow, depending on the chosen model for snow.

			case YAMAGUCHI2012:
				{
					//Calculate ratio density/grain size (see Yamaguchi (2012)):
					double tmp_rho_d=(EMS[SnowpackElement[i]].theta[ICE]*Constants::density_ice)/( (2.*EMS[SnowpackElement[i]].rg) / 1000.);

					//Limit tmp_rho_d to reasonable values, so alpha and especially n remain in numerically stable bounds:
					tmp_rho_d=MAX(2000., tmp_rho_d);
					alpha[i]=4.4E6*pow(tmp_rho_d, -0.98);	//See Eq. 6 in Yamaguchi (2012).
					n[i]=1.+2.7E-3*pow(tmp_rho_d, 0.61);	//See Eq. 7 in Yamaguchi (2012).
					break;
				}

			case YAMAGUCHI2010:
				{
					//Limit grain size, to stay within the bounds of the Van Genuchten parameterizations for snow.
					const double GRAINRADIUSLOWERTHRESHOLD=0.0;		//Lower threshold
					const double GRAINRADIUSUPPERTHRESHOLD=2.0;		//Upper threshold. 2.02 is value for n>1, which is required.
					//Now limit grain sizes
					if(EMS[SnowpackElement[i]].rg>GRAINRADIUSUPPERTHRESHOLD) EMS[SnowpackElement[i]].rg=GRAINRADIUSUPPERTHRESHOLD;
					if(EMS[SnowpackElement[i]].rg<GRAINRADIUSLOWERTHRESHOLD) EMS[SnowpackElement[i]].rg=GRAINRADIUSLOWERTHRESHOLD;

					//Note: rg is in mm, and it is the radius (confirmed by Charles, see DataClasses.h)
					alpha[i]=7.3*(2.*EMS[SnowpackElement[i]].rg)+1.9;			//See Eq. 12 (note d is defined as diameter in mm!) in Yamaguchi (2010).
					n[i]=-3.3*(2.*EMS[SnowpackElement[i]].rg)+14.4;				//See Eq. 11 (note d is defined as diameter in mm!) in Yamaguchi (2010).
					break;
				}

			case YAMAGUCHI2010_ADAPTED:
				{
					//Limit grain size, the parameterizations still hold, but high values of alpha and small values of n are causing numerical troubles.
					const double GRAINRADIUSLOWERTHRESHOLD=0.0;		//Lower threshold
					const double GRAINRADIUSUPPERTHRESHOLD=4.0;		//Upper threshold
					//Now limit grain sizes
					if(EMS[SnowpackElement[i]].rg>GRAINRADIUSUPPERTHRESHOLD) EMS[SnowpackElement[i]].rg=GRAINRADIUSUPPERTHRESHOLD;
					if(EMS[SnowpackElement[i]].rg<GRAINRADIUSLOWERTHRESHOLD) EMS[SnowpackElement[i]].rg=GRAINRADIUSLOWERTHRESHOLD;

					alpha[i]=7.3*(2.*EMS[SnowpackElement[i]].rg)+1.9;			//See Eq. 12 (note d is defined as diameter in mm!) in Yamaguchi (2010).
					//Instead of the linear fit in Yamaguchi (2010), Hirashima (2011) approximated the data with a power law fit, valid for the whole range of grain sizes:
					n[i]=15.68*exp(-0.46*(2.*EMS[SnowpackElement[i]].rg)) + 1.;		//Hirashima (2011), Eq. 17
					break;
				}

			case DAANEN:
				{
					const double GRAINRADIUSLOWERTHRESHOLD=0.0;		//Equal to Yamaguchi adapted
					const double GRAINRADIUSUPPERTHRESHOLD=4.0;		//Equal to Yamaguchi adapted
					//Now limit grain sizes
					if(EMS[SnowpackElement[i]].rg>GRAINRADIUSUPPERTHRESHOLD) EMS[SnowpackElement[i]].rg=GRAINRADIUSUPPERTHRESHOLD;
					if(EMS[SnowpackElement[i]].rg<GRAINRADIUSLOWERTHRESHOLD) EMS[SnowpackElement[i]].rg=GRAINRADIUSLOWERTHRESHOLD;

					alpha[i]=30.*(2.*EMS[SnowpackElement[i]].rg)+12.;
					n[i]=0.800*(2.*EMS[SnowpackElement[i]].rg)+3.;
					break;
				}

			}

			const double tmp_dynamic_viscosity_water=0.001792;				//In Pa/s, from WaterTransport code by Hirashima: 0.001792

			switch ( K_PARAM ) {	//Set saturated hydraulic conductivity

			case SHIMIZU:
				//This formulation for ksat is proposed by Shimizu (1970), and is valid up to 450 kg/m^3. See Equation 5 in Jordan, 1999 + conversion from hydraulic permeability to hydraulic conductivity.
				if(EMS[SnowpackElement[i]].theta[ICE] * Constants::density_ice>450.) {
					ksat[i]=0.077 * (2.*EMS[SnowpackElement[i]].rg / 1000.)*(2.*EMS[SnowpackElement[i]].rg / 1000.) * exp(-0.0078 * 450.) * (Constants::g * Constants::density_water) / tmp_dynamic_viscosity_water;
				} else {
					ksat[i]=0.077 * (2.*EMS[SnowpackElement[i]].rg / 1000.)*(2.*EMS[SnowpackElement[i]].rg / 1000.) * exp(-0.0078 * EMS[SnowpackElement[i]].theta[ICE] * Constants::density_ice) * (Constants::g * Constants::density_water) / tmp_dynamic_viscosity_water;
				}
				break;

			case CALONNE:
				//See: Calonne et al., 3-D image-based numerical computations of snow permeability: links to specific surface area, density, and microstructural anisotropy, TC, 2012.
				ksat[i]=0.75 * (EMS[SnowpackElement[i]].ogs / 1000.)*(EMS[SnowpackElement[i]].ogs / 1000.) * exp(-0.013 * EMS[SnowpackElement[i]].theta[ICE] * Constants::density_ice) * (Constants::g * Constants::density_water) / tmp_dynamic_viscosity_water;
				break;

			}

			//Restore original grain size value from backup
			EMS[SnowpackElement[i]].rg=tmprg;
		} else {  				//Soil
			tmpheight+=dz[i];		//This is only done in soil, so we have a relative reference only for a soil, not for snow.

			switch ( runcase ) {
			case UNIFORMSOIL:
				//Uniform soil
				SetSoil(WFJGRAVELSAND, &theta_r[i], &theta_s[i], &alpha[i], &m[i], &n[i], &ksat[i], &h_e[i]);
				//SetSoil(SAND, &theta_r[i], &theta_s[i], &alpha[i], &m[i], &n[i], &ksat[i], &h_e[i]);
				//SetSoil(SANDYLOAM, &theta_r[i], &theta_s[i], &alpha[i], &m[i], &n[i], &ksat[i], &h_e[i]);
				break;
			case IMISDEFAULT:
				//Default case (IMIS):
				if(tmpheight<=0.25001) {
					//Silt loam
					//SetSoil(ORGANIC, &theta_r[i], &theta_s[i], &alpha[i], &m[i], &n[i], &ksat[i], &h_e[i]);
					//SetSoil(SILTLOAM, &theta_r[i], &theta_s[i], &alpha[i], &m[i], &n[i], &ksat[i], &h_e[i]);
					SetSoil(SANDYLOAM, &theta_r[i], &theta_s[i], &alpha[i], &m[i], &n[i], &ksat[i], &h_e[i]);
				} else { //Gravel/sand
					if(tmpheight<1.001) {
						SetSoil(SAND, &theta_r[i], &theta_s[i], &alpha[i], &m[i], &n[i], &ksat[i], &h_e[i]);
					} else {
						SetSoil(WFJGRAVELSAND, &theta_r[i], &theta_s[i], &alpha[i], &m[i], &n[i], &ksat[i], &h_e[i]);
					}
				}
				break;
			case WFJ:
				//Case WFJ:
				SetSoil(WFJGRAVELSAND, &theta_r[i], &theta_s[i], &alpha[i], &m[i], &n[i], &ksat[i], &h_e[i]);
				break;
			case CDP:
				//Case Col de Porte
				SetSoil(SANDYLOAM, &theta_r[i], &theta_s[i], &alpha[i], &m[i], &n[i], &ksat[i], &h_e[i]);
				break;
			case SNOFILE:
				if(EMS[SnowpackElement[i]].rg < 0.5) {
					SetSoil(ORGANIC, &theta_r[i], &theta_s[i], &alpha[i], &m[i], &n[i], &ksat[i], &h_e[i]);
				} else if (EMS[SnowpackElement[i]].rg < 1.) {
					SetSoil(CLAY, &theta_r[i], &theta_s[i], &alpha[i], &m[i], &n[i], &ksat[i], &h_e[i]);
				} else if (EMS[SnowpackElement[i]].rg < 2.) {
					SetSoil(CLAYLOAM, &theta_r[i], &theta_s[i], &alpha[i], &m[i], &n[i], &ksat[i], &h_e[i]);
				} else if (EMS[SnowpackElement[i]].rg < 3.) {
					SetSoil(LOAM, &theta_r[i], &theta_s[i], &alpha[i], &m[i], &n[i], &ksat[i], &h_e[i]);
				} else if (EMS[SnowpackElement[i]].rg < 4.) {
					SetSoil(LOAMYSAND, &theta_r[i], &theta_s[i], &alpha[i], &m[i], &n[i], &ksat[i], &h_e[i]);
				} else if (EMS[SnowpackElement[i]].rg < 5.) {
					SetSoil(SAND, &theta_r[i], &theta_s[i], &alpha[i], &m[i], &n[i], &ksat[i], &h_e[i]);
				} else if (EMS[SnowpackElement[i]].rg < 6.) {
					SetSoil(SANDYCLAY, &theta_r[i], &theta_s[i], &alpha[i], &m[i], &n[i], &ksat[i], &h_e[i]);
				} else if (EMS[SnowpackElement[i]].rg < 7.) {
					SetSoil(SANDYCLAYLOAM, &theta_r[i], &theta_s[i], &alpha[i], &m[i], &n[i], &ksat[i], &h_e[i]);
				} else if (EMS[SnowpackElement[i]].rg < 8.) {
					SetSoil(SANDYLOAM, &theta_r[i], &theta_s[i], &alpha[i], &m[i], &n[i], &ksat[i], &h_e[i]);
				} else if (EMS[SnowpackElement[i]].rg < 9.) {
					SetSoil(SILT, &theta_r[i], &theta_s[i], &alpha[i], &m[i], &n[i], &ksat[i], &h_e[i]);
				} else if (EMS[SnowpackElement[i]].rg < 10.) {
					SetSoil(SILTYCLAY, &theta_r[i], &theta_s[i], &alpha[i], &m[i], &n[i], &ksat[i], &h_e[i]);
				} else if (EMS[SnowpackElement[i]].rg < 11.) {
					SetSoil(SILTYCLAYLOAM, &theta_r[i], &theta_s[i], &alpha[i], &m[i], &n[i], &ksat[i], &h_e[i]);
				} else if (EMS[SnowpackElement[i]].rg < 12.) {
					SetSoil(SILTLOAM, &theta_r[i], &theta_s[i], &alpha[i], &m[i], &n[i], &ksat[i], &h_e[i]);
				} else {
					SetSoil(WFJGRAVELSAND, &theta_r[i], &theta_s[i], &alpha[i], &m[i], &n[i], &ksat[i], &h_e[i]);
				}
				break;
			}
			//I encountered the following problem: fully saturated soil and freezing water: there is not enough place to store the ice!!!
			//In the old snowpack code, this problem was solved by keeping the increase in volume when all the water in the element would freeze, free as theta[AIR].
			//However, this will not work in the Richards, as theta[WATER] is varying per time step. So we keep free a volume as if the soil is saturated AND will freeze:
			EMS[SnowpackElement[i]].theta[SOIL]=1.-((Constants::density_water/Constants::density_ice)*theta_s[i]);	//Determine the soil content based on the pore space
		}

		//Calculate m:
		m[i]=(n[i]-1.)/n[i];

		//Calculate saturation at cut-off point h_e (see Ippisch et al (2006)).
		Sc[i]=pow((1.+pow(alpha[i]*fabs(h_e[i]), n[i])), -1.*m[i]);

		//Get theta_i_n
		if(EMS[i].theta[SOIL]>Constants::eps2) {	//Only for soil
			theta_i_n[i]=EMS[SnowpackElement[i]].theta[ICE];

			//Get T_melt that suffices partitioning pressure head into part for ice and part for water
			if(theta_i_n[i]>0.) {			//If there is ice in soil, calculate freezing point depression.
				const double hw0=fromTHETAtoH(EMS[SnowpackElement[i]].theta[WATER]+(EMS[SnowpackElement[i]].theta[ICE]*(Constants::density_ice/Constants::density_water)), theta_r[i], theta_s[i], alpha[i], m[i], n[i], Sc[i], h_e[i], h_d);
				T_melt[i]=T_0+((Constants::g*T_0)/delF)*hw0;
			} else {
				T_melt[i]=T_0;			//Melting point is just the standard melting point.
			}
		} else {					//For snow
			theta_i_n[i]=0.;			//This sounds strange for snow, but the idea is that ice in snow functions as soil in soil (being the matrix)
			T_melt[i]=T_0;				//For snow, we currently don't have anything with freezing point depression, as we have in soil.
		}

		//Determine what pressure head should be considered "dry".
		//Explanation: cold dry new snow layers are initialized with this value. We need to make sure that ALL the other layers have at least a higher pressure head when they contain at least a little bit of water. Else, various numerical troubles arise.
		//In case the value is too high, we get fluxes out of the completely dry snow layer, and a too low value causes many numerical difficulties as in that case, it represents a much stronger gradient in pressure head than necessary (many iterations and small time steps).
		//So we check for each particular layer what pressure head is associated with a theta[WATER] that is a smaller deviation from theta_r then the solver will resolve.
		const double tmp_head=fromTHETAtoHforICE(theta_r[i]+(REQUIRED_ACCURACY_THETA/10.), theta_r[i], theta_s[i], alpha[i], m[i], n[i], Sc[i], h_e[i], h_d, theta_i_n[i]);
		if(h_d>tmp_head) h_d=tmp_head;
		if(i==uppernode) h_d_uppernode=tmp_head;	//We store this value in order to use it for the LIMITEDFLUXEVAPORATION
		if (WriteOutNumerics_Level3==true) 
			std::cout << "H_D at " << i << ": " << std::scientific << tmp_head << std::fixed << " [alpha: " << alpha[i] << "; m: " << m[i] << "; n: " << n[i] << "; Sc: " << Sc[i] << "; h_e: " << h_e[i] << "\n";
	}
	if (WriteOutNumerics_Level2==true) std::cout << "MIN_HEAD: " << std::scientific << h_d << std::fixed << "\n";


	//Coupling of SNOWPACK domain to RE-solver domain. Makes sure the EMS.theta[XXX] are within the limits specified by the Van Genuchten parameterizations.
        for (i = uppernode; i >= lowernode; i--) {	//Cycle over all Richards solver domain layers
  		//Now calculate the theta that should be considered "dry soil".
		theta_d[i]=fromHtoTHETAforICE(h_d, theta_r[i], theta_s[i], alpha[i], m[i], n[i], Sc[i], h_e[i], 0.);

		//Now check if this case is not too extreme
		if(theta_d[i]<theta_r[i]+(REQUIRED_ACCURACY_THETA/1000.)) {
			theta_d[i]=theta_r[i]+(REQUIRED_ACCURACY_THETA/1000.);
		}

		//Now make sure that the water content in SNOWPACK's ElementData matches the soil settings (not too wet, not too dry):
		// 1) Not too wet
		if(EMS[SnowpackElement[i]].theta[SOIL]>Constants::eps2) {		//For soil
			if(EMS[SnowpackElement[i]].theta[WATER]+(EMS[SnowpackElement[i]].theta[ICE]*(Constants::density_ice/Constants::density_water)) > theta_s[i]) {
				wateroverflow[i]+=(EMS[SnowpackElement[i]].theta[WATER]+(EMS[SnowpackElement[i]].theta[ICE]*(Constants::density_ice/Constants::density_water))-theta_s[i]);
				EMS[SnowpackElement[i]].theta[WATER]=theta_s[i]-(EMS[SnowpackElement[i]].theta[ICE]*(Constants::density_ice/Constants::density_water));
			}
		} else {								//For snow
  			if(EMS[SnowpackElement[i]].theta[WATER] > theta_s[i]) {
				wateroverflow[i]+=EMS[SnowpackElement[i]].theta[WATER]-theta_s[i];
				EMS[SnowpackElement[i]].theta[WATER]=theta_s[i];
			}
		}

		// 2) Not too dry
		if(EMS[SnowpackElement[i]].theta[SOIL]>Constants::eps2) {		//For soil
			if(AllowDrySoilLayers==true) {
				if(not( (EMS[SnowpackElement[i]].theta[WATER]+(EMS[SnowpackElement[i]].theta[ICE]*(Constants::density_ice/Constants::density_water))) < theta_r[i]) && (EMS[SnowpackElement[i]].theta[WATER]+(EMS[SnowpackElement[i]].theta[ICE]*(Constants::density_ice/Constants::density_water))) < theta_d[i]) {
					wateroverflow[i]+=( (EMS[SnowpackElement[i]].theta[WATER]+(EMS[SnowpackElement[i]].theta[ICE]*(Constants::density_ice/Constants::density_water))) - theta_d[i]);
					EMS[SnowpackElement[i]].theta[WATER]=theta_d[i]-(EMS[SnowpackElement[i]].theta[ICE]*(Constants::density_ice/Constants::density_water));
					activelayer[i]=true;
				} else {
					if( (EMS[SnowpackElement[i]].theta[WATER]+(EMS[SnowpackElement[i]].theta[ICE]*(Constants::density_ice/Constants::density_water))) < theta_r[i]) {
						activelayer[i]=false;
					} else {
						activelayer[i]=true;
					}
				}
			} else {
				if( (EMS[SnowpackElement[i]].theta[WATER]+(EMS[SnowpackElement[i]].theta[ICE]*(Constants::density_ice/Constants::density_water))) < theta_d[i]) {
					wateroverflow[i]+=( (EMS[SnowpackElement[i]].theta[WATER]+(EMS[SnowpackElement[i]].theta[ICE]*(Constants::density_ice/Constants::density_water))) - theta_d[i]);
					EMS[SnowpackElement[i]].theta[WATER]=theta_d[i]-(EMS[SnowpackElement[i]].theta[ICE]*(Constants::density_ice/Constants::density_water));
				}
				activelayer[i]=true;
			}
		} else {
			if(AllowDrySnowLayers==true) {
				//For snow, we have to melt ice to create theta_r!!
				if(not(EMS[SnowpackElement[i]].theta[WATER]<theta_r[i]) && EMS[SnowpackElement[i]].theta[WATER]<theta_d[i]) {
					const double tmp_missing_theta=(theta_d[i]-EMS[SnowpackElement[i]].theta[WATER]);	//Not too dry (original)

					//Note: we do the book keeping for wateroverflow[i] not here, but afterwards, when we know how much water was used.

					//The energy required for this melt is stored. We apply it only if the layer will get enough water to have this water passed on to the rest of SNOWPACK.
					//Note, the plus sign is because we can pass the same SNOWPACK layer several times, when it is divided into sublayers. Then, for the next sublayer, the water is already there, so only 0 is added.
					dT[SnowpackElement[i]]+=tmp_missing_theta*(Constants::density_water/Constants::density_ice) / ((EMS[SnowpackElement[i]].c[TEMPERATURE] * EMS[SnowpackElement[i]].Rho) / ( Constants::density_ice * Constants::lh_fusion ));

					EMS[SnowpackElement[i]].theta[WATER]+=tmp_missing_theta;
					EMS[SnowpackElement[i]].theta[ICE]-=tmp_missing_theta*(Constants::density_water/Constants::density_ice);
					activelayer[i]=true;
				} else {
					if(EMS[SnowpackElement[i]].theta[WATER]<theta_r[i]) {
						activelayer[i]=false;
					} else {
						activelayer[i]=true;
					}
				}
			} else {
				//For snow, we have to melt ice to create theta_r!!
				if(theta_d[i]>0.) {	//Do we have room for melt water?
					if(EMS[SnowpackElement[i]].theta[WATER]<theta_d[i]) {
						const double tmp_missing_theta=(theta_d[i]-EMS[SnowpackElement[i]].theta[WATER]);	//Not too dry (original)
						dT[SnowpackElement[i]]+=tmp_missing_theta*(Constants::density_water/Constants::density_ice) / ((EMS[SnowpackElement[i]].c[TEMPERATURE] * EMS[SnowpackElement[i]].Rho) / ( Constants::density_ice * Constants::lh_fusion ));
						if (WriteOutNumerics_Level1==true) 
							std::cout << "WATER CREATED (" << tmp_missing_theta << "): i=" << i << " --- dT=" << dT[SnowpackElement[i]] << " T=" << EMS[SnowpackElement[i]].Te << "  theta[WATER]=" << EMS[SnowpackElement[i]].theta[WATER] << " theta[ICE]=" << EMS[SnowpackElement[i]].theta[ICE] << "\n";
						EMS[SnowpackElement[i]].theta[WATER]+=tmp_missing_theta;
						EMS[SnowpackElement[i]].theta[ICE]-=tmp_missing_theta*(Constants::density_water/Constants::density_ice);
					}
					activelayer[i]=true;
				} else {		//Else we have a pure ice layer
					activelayer[i]=false;
				}
			}
		}


		//Now copy the EMS water content into the working arrays to solve Richards-equation (so this is the important part were this function is coupled to the rest of SNOWPACK).
		if(activelayer[i]==true) {
			// Now calculate initial pressure head:
			h_n[i]=fromTHETAtoHforICE(EMS[SnowpackElement[i]].theta[WATER], theta_r[i], theta_s[i], alpha[i], m[i], n[i], Sc[i], h_e[i], h_d, theta_i_n[i]);
			theta_n[i]=fromHtoTHETAforICE(h_n[i], theta_r[i], theta_s[i], alpha[i], m[i], n[i], Sc[i], h_e[i], theta_i_n[i]);	//This is the current theta, which we determine from h_n[i].
		} else {
			theta_n[i]=EMS[SnowpackElement[i]].theta[WATER];
			h_n[i]=h_d;
		}

		//Determine source/sink term
		s[i]=0.;								//Reset source/sink term

		//Add wateroverflow (Remember: units wateroverflow [m^3/m^3]):
		if(alpine3d==false && (wateroverflow[i]>0 || SafeMode==false)) {	//In SafeMode, we don't allow the negative wateroverflow to be used as sink term, as a negative wateroverflow is caused by initialization of very dry snow layers, so the sink term would basically be a sink term in very dry conditions, which is numerically unstable.
			if(i==uppernode) {
				surfacefluxrate+=(wateroverflow[i]*dz[i])/sn_dt;
				wateroverflow[i]=0.;
			} else {
				if((wateroverflow[i]*dz[i])/sn_dt < ksat[i+1]) {	//Check if influx is not too large
					s[i]+=wateroverflow[i]/sn_dt;			//These terms mainly are caused by merging elements, where the mass of the upper element is added to the lower one. This can lead to too much water in a certain element. We add this as a source term.
					wateroverflow[i]=0.;				//Since we have put wateroverflow in source/sink term, it's not an overflow anymore for this layer.
				} else {						//Else limit influx and throw other water away... So this is a water hole in the model. I suggest making a variable MS_LATERALRUNOFF to track this water.
					s[i]+=(ksat[i]/dz[i]);
					wateroverflow[i]-=(ksat[i]/dz[i])*sn_dt;
				}
			}
		}

		//Now add soilsurfacesourceflux (in case RemoveElements removed the lowest snow element):
		if(soilsurfacesourceflux>0. && i==nsoillayers_richardssolver) {		//We assign source flux in the lowest snow element if the source flux is >0. This can only be the case when we use RE for snow, so we don't have to check for this.
			//Remember: soilsurfacesourceflux=[m^3/m^2/s]
			s[i]+=soilsurfacesourceflux/dz[i];				//Soilsurfacesourceflux>0. if we remove the first snow element above the soil AND there are more snow layers (else it is a surfaceflux) AND we use RE for snow.
		}

		//Add source/sink term from other parts of SNOWPACK (in particular Canopy.cc)
		s[i]+=EMS[i].lwc_source/sn_dt;
		EMS[i].lwc_source=0.;		// Now that we used the variable, reset it.

		//To now the flux of water in/out of the model domain due to the source/sink term.
		totalsourcetermflux+=s[i]*dz[i];
	}


	//Initialize upper boundary in case of Dirichlet
	if(TopBC==DIRICHLET) {
		aTopBC=DIRICHLET;
		htop=h_n[uppernode];
	}

	//Initialize lower boundary in case of Dirichlet
	if(BottomBC==DIRICHLET) {
		hbottom=h_n[lowernode];
	}

	//Initialize lower boundary in case of WATERTABLE: saturated
	if(BottomBC==WATERTABLE) {
		aBottomBC=DIRICHLET;
		hbottom=h_e[lowernode];
		h_n[lowernode]=hbottom;

		wateroverflow[lowernode]+=(theta_n[lowernode]);	//First we remove all water from the lowest element
		theta_n[lowernode]=fromHtoTHETAforICE(h_n[lowernode], theta_r[lowernode], theta_s[lowernode], alpha[lowernode], m[lowernode], n[lowernode], Sc[lowernode], h_e[lowernode], theta_i_n[lowernode]);
		wateroverflow[lowernode]-=(theta_n[lowernode]);	//Then we add the saturated boundary water content from the lowest element.
	}


	//Note: there are 2 iterations. First, the iteration starts to match the Richards solver time step to the SNOWPACK time step. Simple example: assume SNOWPACK time step is 15 minutes and
	//Richards solver time step is 1 minute, there should be 15 iterations to match the solution to the SNOWPACK time step.
	//Then, for each time step of the Richard solver, iterations are necessary to find the solution to the equation.
	int nsteps=0;			//Counts the number of time steps in the Richards solver.
	bool StopLoop=false;		//Will switch to true when the integrated time step matches the SNOWPACK time step.
	bool DoRewindFlag=false;	//Will switch to true when the time step has to be redone with a smaller time step.

	//Determine mass at beginning of snowpack time step.
	mass1=0.;
	for (i = uppernode; i >= lowernode; i--) {
		mass1+=(theta_n[i]+(theta_i_n[i]*(Constants::density_ice/Constants::density_water)))*dz[i];
	}

	do
	{
		if(DoRewindFlag==false) {		//Only if we are not doing a rewind, we should increase the number of steps (else it basically is the same time step).
			nsteps++;			//Increase the number of steps
			niter_nrewinds=0;		//Reset rewind counter
		}

		Xdata.ReSolver_dt=dt;			//Store the last used time step.
		if ((TimeAdvance+dt)>=snowpack_dt) {	//If our time step is so large that the integrated time step will exceed the SNOWPACK time step, we limit the dt for the current time step...
			dt=snowpack_dt-TimeAdvance;	//...so it matches exactly the SNOWPACK time step.
			StopLoop=true;			//And we set the switch to stop the Richards solver.
		}
		TimeAdvance+=dt;			//Update the total time in this time step. This variable is used to match SNOWPACK time steps.

		//Prepare for next time step:
		niter=0;				//reset iter counter
		accuracy=0.;				//reset accuracy.


		//Set Solver
		ActiveSolver=PreferredSolver;		//We set the active solver to the preferred solver


		//Initialize values for the first iteration (iteration m)
		for (i = uppernode; i >= lowernode; i--) {
			// Note, it is not possible to do an educated guess. The guess should be mass-conservative, which is very difficult to achieve.
			h_np1_m[i]=h_n[i];
			theta_np1_m[i]=theta_n[i];
			theta_i_np1_m[i]=theta_i_n[i];
		}

		//Write out initial water content
		if( (WriteOutNumerics_Level1==true && nsteps==1 && niter_nrewinds==0 ) || (WriteOutNumerics_Level2==true)) {
			for (i = uppernode; i >= lowernode; i--) {
				const string is_active = (activelayer[i])? "true" : "false" ;
				std::cout << "ITER: " << niter << " i: " << i << " active? " << is_active << std::setprecision(15) << "; h_n: " << h_n[i] << " (h_np1: " << h_np1_m[i] << ") theta: " << theta_n[i] << std::setprecision(6) << "(" << theta_r[i] << "-" << theta_s[i] << ") ice: " << EMS[SnowpackElement[i]].theta[ICE] << "/" << theta_i_n[i] << " (vg_params: " << alpha[i] << " " << m[i] << " " << n[i] << ")\n";
			}
		}

		DoRewindFlag=false;							//Reset DoRewindFlag. We do it now, just before starting the solver, so we can use the status of this flag to initialize the solver properly.
		boolConvergence=false;							//Default is no convergence, until proven otherwise

		while (boolConvergence==false && DoRewindFlag==false) {			//In theory, this can create an endless loop, but for this, I put a throw in the code when no convergence is achieved, because then the situation is hopeless anyway.
			niter++;
			niter_snowpack_dt++;
			memstate++;
			int solver_result=0;

			//Prepare matrices
			//Update state properties
			for (i = uppernode; i >= lowernode; i--) {
				if(activelayer[i]==true) {
					//Calculate theta from h
					theta_np1_m[i]=fromHtoTHETAforICE(h_np1_m[i], theta_r[i], theta_s[i], alpha[i], m[i], n[i], Sc[i], h_e[i], theta_i_np1_m[i]);

					//Determine Se
					if(h_np1_m[i]<h_e[i]) {			//See Ippisch (2006)
						//Calculate dimensionless saturation
						Se[i] = ((theta_np1_m[i] + (theta_i_np1_m[i]*(Constants::density_ice/Constants::density_water)) - theta_r[i])/(theta_s[i] - theta_r[i]));
						if(Se[i]<0.) {
							//The formulation of Se[i] as used here may lead to very small negative values for Se. These are corrected here.
							if(Se[i]<-1E-12) std::cout << "WARNING: Se[" << i << "]=" << std::scientific << Se[i] << std::fixed << ".\n";	//This points towards a more serious problem, so give a warning...
							Se[i]=0.;
						}
					} else {					//In case of saturation:
						Se[i]=1.;
					}

					//Determine hydraulic conductivity
					if(Se[i]<1.) {
					  	//Compute the hydraulic conductivity (see Ippisch, 2006)
						K[i]=ksat[i]*sqrt(Se[i])*pow((1.-(pow(1.-pow(Se[i]*Sc[i],(1./m[i])),m[i])))/(1.-pow(1.-pow(Sc[i],(1./m[i])), m[i])),2.);
					} else {
						K[i] = ksat[i];
					}

					//Applying ice impedance on K
					if(ApplyIceImpedance==true) {
						const double omega=7.;		//See Zhao et al. (1997) and Hansson et al. (2004)  [Dall'Amicao, 2011].
						if( SnowpackElement[i] < nsoillayers_snowpack && theta_i_np1_m[i]>0. && K[i]>0. ) {	//Only for soil and when there is ice in the soil
							//q=theta_i_np1_m[i]/(theta_s[i]-theta_r[i]);					//This is how Dall'Amico presents it, but it is based on Hanssen (2004), who defines it as:
							const double q = (theta_i_np1_m[i]*(Constants::density_ice/Constants::density_water))/((theta_np1_m[i]+(theta_i_np1_m[i]*(Constants::density_ice/Constants::density_water)))-theta_r[i]);		//Hanssen (2004).
							impedance[i]=pow(10., -1.*omega*q);
						} else {
							impedance[i]=1.;
						}
						K[i]*=impedance[i];
					}

			 		//Calculate the specific moisture capacity (which is derivative d.theta/d.h)
					if(Se[i]<1.)	{	//No saturation
						C[i]=alpha[i]*n[i]*m[i]*((theta_s[i]-theta_r[i])/Sc[i])*(pow((alpha[i]*fabs(h_np1_m[i])), (n[i]-1.)))*(pow(1.+pow((alpha[i]*fabs(h_np1_m[i])), n[i]), (-1.*m[i]-1.)));
						if(isnan(C[i])) solver_result=-1;
					} else {		//Saturation
						C[i]=0.;
					}
				} else {	//If not an active layer
					K[i]=0.;
					C[i]=0.;
				}
				if(WriteOutNumerics_Level3==true) std::cout << "HYDPROPS: i=" << i << std::scientific << " Se=" << Se[i] << " C=" << C[i] << " K=" << K[i] << ".\n" << std::fixed;
			}

			for (i = lowernode; i <= uppernode; i++) {
				//Determine K at interface nodes
				// 1) Determine k_np1_m_ip12
				if (i!=uppernode && activelayer[i+1]==true) {
					//For the rest of the domain, we might have heterogeneous soils, so we have to approximate the hydraulic conductivity at the interface nodes.

					switch (K_AVERAGETYPE) {
						case ARITHMETICMEAN:
						{
							k_np1_m_ip12[i]=.5*(K[i]+K[i+1]);
							break;
						}

						case GEOMETRICMEAN:
						{
							k_np1_m_ip12[i]=sqrt(K[i]*K[i+1]);
							break;
						}

						case HARMONICMEAN:
						{
							if(K[i]>0. && K[i+1]>0.) {
								k_np1_m_ip12[i]=2./(1./K[i]+1./K[i+1]);
							} else {
								k_np1_m_ip12[i]=0.;
							}
							break;
						}

						case MINIMUMVALUE:
						{
							if(K[i]>K[i+1]) {
								k_np1_m_ip12[i]=K[i+1];
							} else {
								k_np1_m_ip12[i]=K[i];
							}
							break;
						}

						case UPSTREAM:
						{
							if( ((h_np1_m[i+1]-h_np1_m[i])/dz_down[i+1]) - cos_sl > 0.) {
								k_np1_m_ip12[i]=K[i];
							} else {
								k_np1_m_ip12[i]=K[i+1];
							}
							break;
						}
					}
				} else {
					//For the boundaries, we neglect gradients in K. This corresponds to the specified fluid flux boundary condition (Equation 4 of McCord, WRR, 1991).
					if(i==uppernode) {
						k_np1_m_ip12[i]=K[i];
					} else {
						k_np1_m_ip12[i]=0.;
					}
				}

				// 2) Determine k_np1_m_im12
				if (i!=lowernode && activelayer[i-1]==true) {
					// The following statement needs to be true, else you won't have mass balance in the solver!
					k_np1_m_im12[i]=k_np1_m_ip12[i-1];
				} else {
					//For the boundaries, we neglect gradients in K. This corresponds to the specified fluid flux boundary condition (Equation 4 of McCord, WRR, 1991).
					if(i==lowernode) {
						k_np1_m_im12[i]=K[i];
					} else {
						k_np1_m_im12[i]=0.;
					}
				}
				if(WriteOutNumerics_Level3==true) 
					std::cout << "HYDRCONDUCT: node " << i << std::scientific << " -- K=" << K[i] << " k_np1_m_im12=" << k_np1_m_im12[i] << " k_np1_m_ip12=" << k_np1_m_ip12[i] << std::fixed << "   impedance: " << impedance[i] << "\n";
			}

			//Determine which and how boundary conditions should be applied:
			if (TopBC==DIRICHLET) {
				aTopBC=DIRICHLET;				//Set Dirichlet BC
				TopFluxRate=0.;					//Dirichlet BC, so no flux
				theta_np1_m[uppernode]=theta_n[uppernode];
			} else if (TopBC==NEUMANN) {
				//Note: TopFluxRate is defined as gradient over pressure head. For influx, pressure head is increasing with increasing height, so TopFluxRate is positive.
				//Units: surfacefluxrate=[m^3/m^2/s]
				aTopBC=NEUMANN;					//Set Neumann BC
				TopFluxRate=surfacefluxrate;			//Flux for Neumann BC
			} else if (TopBC==LIMITEDFLUXEVAPORATION || TopBC==LIMITEDFLUXINFILTRATION || TopBC==LIMITEDFLUX) {
				//Now check if the topflux is not too big or small, giving positive pressure heads. For example: during heavy rain, the rain rate can be much more than handled by the soil. The upper layer will blow up the model in this case, as it cannot deal with all the incoming water. So the fluxes should not exceed dry or saturated conditions.
				aTopBC=NEUMANN;					// Limited flux is technically just Neumann, but with limited fluxes.
				if(niter==1) TopFluxRate=surfacefluxrate;	// Initial guess for Neumann BC
				// Now reduce flux when necessary:
  				if((TopBC == LIMITEDFLUXINFILTRATION || TopBC == LIMITEDFLUX) && (TopFluxRate>0.) && (
				     (LIMITEDFLUXINFILTRATION_soil==true && int(nsoillayers_snowpack)==int(nE))
				        || (LIMITEDFLUXINFILTRATION_snowsoil==true && int(nsoillayers_snowpack)<int(nE) && toplayer==nsoillayers_snowpack)
				           || (LIMITEDFLUXINFILTRATION_snow==true && int(nsoillayers_snowpack)<int(nE)))) {
					// Influx condition
					// Determine the limiting flux:
					const double flux_compare =														//The limiting flux is:
					        (dz[uppernode]*(theta_s[uppernode] - (theta_np1_m[uppernode] + theta_i_np1_m[uppernode]))/dt)					// net flux that would lead to saturation of the top layer
					                + ((uppernode>0) ? k_np1_m_im12[uppernode]*(((h_np1_m[uppernode]-h_np1_m[uppernode-1])/dz_down[uppernode]) + cos_sl) : 0.);	// plus what could leave below

					// For alpine3d simulations, we are stricter for the sake of stability: we also don't allow a positive influx when there is ponding inside the model domain:
					if(alpine3d==true) {
						bool isPonding=false;
						for(int jj=lowernode; jj<=uppernode; jj++) {
							if(h_np1_m[jj]>h_e[jj]) isPonding=true;
						}
						if(isPonding==true) TopFluxRate=0.;
					}

					if((0.999*flux_compare) < TopFluxRate) {		//Limit flux if necessary. Note: we multiply flux_compare with 0.999 because flux_compare can be
						TopFluxRate=MAX(0., (0.999*flux_compare));	//regarded as the asymptotic case from which we want to stay away a little.
					}
				}
  				if((TopBC == LIMITEDFLUXEVAPORATION || TopBC == LIMITEDFLUX) && (TopFluxRate<0.) && ((LIMITEDFLUXEVAPORATION_soil==true && (int(nsoillayers_snowpack)==int(nE) || toplayer==nsoillayers_snowpack)) || (LIMITEDFLUXEVAPORATION_snow==true && int(nsoillayers_snowpack)<int(nE)))) {
					// Outflux condition
					const double head_compare=h_d_uppernode;
					const double flux_compare=k_np1_m_ip12[uppernode]*(((head_compare-h_np1_m[uppernode])/dz_up[uppernode]) + cos_sl);
					if(flux_compare > TopFluxRate) {
						TopFluxRate=MIN(0., flux_compare);
					}
				}
			} else if (TopBC==WATERTABLE) {
				std::cout << "ERROR in ReSolver1d.cc: WATERTABLE cannot be applied as top boundary condition (doesn't make sense)!\n";
				throw;
			} else if (TopBC==FREEDRAINAGE) {
				std::cout << "ERROR in ReSolver1d.cc: FREEDRAINAGE cannot be applied as top boundary condition (doesn't make sense)!\n";
				throw;
			} else if (TopBC==SEEPAGEBOUNDARY) {
				std::cout << "ERROR in ReSolver1d.cc: SEEPAGEBOUNDARY cannot be applied as top boundary condition (doesn't make sense)!\n";
				throw;
			} else if (TopBC==GRAVITATIONALDRAINAGE) {
				std::cout << "ERROR in ReSolver1d.cc: GRAVITATIONALDRAINAGE cannot be applied as top boundary condition (doesn't make sense)!\n";
				throw;
			}


			if (BottomBC==DIRICHLET) {
				aBottomBC=DIRICHLET;		//Set Dirichlet BC.
				BottomFluxRate=0.;		//Dirichlet BC, so no prescribed flux.
				theta_np1_m[lowernode]=theta_n[lowernode];
			} else if (BottomBC==WATERTABLE) {
				aBottomBC=DIRICHLET;		//Water table is a Dirichlet BC.
				BottomFluxRate=0.;		//Dirichlet BC, so no prescribed flux.
				theta_np1_m[lowernode]=theta_n[lowernode];
			} else if (BottomBC==NEUMANN) {
				aBottomBC=NEUMANN;		//Set Neumann BC.
				//Note: BottomFluxRate is defined as gradient over pressure head. For outflux (drainage), pressure head is increasing with increasing height, so BottomFluxRate is positive.
				BottomFluxRate=0.0000005;	//Flux for Neumann BC.
			} else if (BottomBC==FREEDRAINAGE) {
				//First calculate flux between lowest and lowest+1 element.
				const double tmpgrad=((h_np1_m[lowernode+1]-h_np1_m[lowernode])/dz_up[lowernode]);	//Note: flux would be (tmpgrad * K).
				if((tmpgrad+cos_sl) < 0.) {
					//In this case, we would create influx at lower boundary, which does not work with FREEDRAINAGE.
					//Then set zero flux:
					aBottomBC=NEUMANN;
					BottomFluxRate=0.;
				} else {
					aBottomBC=NEUMANN;
					//Now, prescribe flux at lower boundary equivalent to tmpgrad
					BottomFluxRate=(tmpgrad+cos_sl)*k_np1_m_im12[lowernode];
				}
			} else if (BottomBC==SEEPAGEBOUNDARY) {
				//Neumann with flux=0 in case of unsaturated
				//Dirichlet with h_bottom=0 in case of saturated
				if(h_n[lowernode+1]<0.) {
					aBottomBC=NEUMANN;
					BottomFluxRate=0.;
				} else {
					aBottomBC=DIRICHLET;
					hbottom=0.;
					BottomFluxRate=0.;
				}
			} else if (BottomBC==GRAVITATIONALDRAINAGE) {
				// See: Xubin Zeng and Mark Decker (2008). Improving the Numerical Solution of Soil Moisture–Based Richards Equation for Land Models with a Deep or Shallow Water Table
				// http://dx.doi.org/10.1175/2008JHM1011.1
				aBottomBC=NEUMANN;
				BottomFluxRate=k_np1_m_im12[lowernode];
			} else if (BottomBC==LIMITEDFLUX) {
				//Probably also not necessary.
				std::cout << "ERROR in ReSolver1d.cc: No implementation for LIMITEDFLUX lower boundary condition. Either choose a saturated DIRICHLET (lower boundary in water table), or choose GRAVITATIONAL or FREEDRAINAGE (lower boundary not in water table).\n";
				throw;
			} else if (BottomBC==LIMITEDFLUXEVAPORATION) {
				std::cout << "ERROR in ReSolver1d.cc: LIMITEDFLUXEVAPORATION cannot be applied as bottom boundary condition (doesn't make sense)!\n";
				throw;
			} else if (BottomBC==LIMITEDFLUXINFILTRATION) {
				std::cout << "ERROR in ReSolver1d.cc: LIMITEDFLUXINFILTRATION cannot be applied as bottom boundary condition (doesn't make sense)!\n";
				throw;
			}


			if (WriteOutNumerics_Level2==true) 
				std::cout << "BOUNDARYTOPFLUX: [ BC: " << TopBC << "] " << std::scientific << TopFluxRate << " " << surfacefluxrate << " " << theta_n[lowernode] << " " << K[lowernode] << " " << ((h_np1_mp1[lowernode])+(((TopFluxRate/k_np1_m_im12[lowernode])-1.)*dz_down[lowernode])) << " " << h_np1_mp1[lowernode] << " " << k_np1_m_im12[lowernode] << " " << (TopFluxRate/k_np1_m_im12[lowernode]) << "\n" << std::fixed;
			if (WriteOutNumerics_Level2==true) 
				std::cout << "NUMERICS: BCTOP: " << TopBC << std::scientific << "  TOPFLUXRATE = " << TopFluxRate << " SURFACEFLUXRATE = " << surfacefluxrate << "\n" << std::fixed;


			if (niter==1) {
				if (int(nsoillayers_snowpack)<int(nE)) {	//We have snow layers
					// See McCord (1996). snowsoilinterfaceflux > 0 means influx!
					snowsoilinterfaceflux_before=((((h_n[nsoillayers_richardssolver]-h_n[nsoillayers_richardssolver-1])/dz_up[nsoillayers_richardssolver-1])+cos_sl)*k_np1_m_ip12[nsoillayers_richardssolver-1]*dt);
				}
			}


			//Solve equation
			std::fill(ainv.begin(), ainv.end(), 0.);	//This is very important: with inverting the matrix, it may become non-tridiagonal! So we have to explicitly set its elements to 0, because some of the for-loops only touch the tridiagonal part of the matrix.
			for (i = uppernode; i >= lowernode; i--) {
				j=i;	//As matrix A is tridiagonal, so it can be filled very efficiently. However, I keep the notation of i and j, so it's better understood how the structure of A is. I only evaluate i==j.
				//This part is for the DGESVD/DGESDD solver, which uses full matrix a (ainv). We always need them, because in case DGTSV fails, we should be able to fall back on DGESVD/DGESDD:
				if(i==j) {
					//Set up the matrix diagonal
					ainv[j*(uppernode+1)+i]=(1./dt)*C[i];

					//The following two lines assume Neumann boundary conditions (for upper and lowernode, one of the terms drop out). If Dirichlet is used, this will be corrected later.
					if(i!=lowernode) ainv[j*(uppernode+1)+i]+=(1./(dz_[i]))*((k_np1_m_im12[i]/(dz_down[i])));
					if(i!=uppernode) ainv[j*(uppernode+1)+i]+=(1./(dz_[i]))*((k_np1_m_ip12[i]/(dz_up[i])));

					//Correct diagonal in case of Dirichlet
					if(aTopBC==DIRICHLET && i==uppernode) {
						ainv[i*(uppernode+1)+i]=1.;
					}
					if(aBottomBC==DIRICHLET && i==lowernode) {
						ainv[i*(uppernode+1)+i]=1.;
					}

					//Set up the matrix upper and lower diagonals
					if(i!=lowernode) ainv[i*(uppernode+1)+(i-1)]=(-1./(dz_[i]))*((k_np1_m_im12[i]/(dz_down[i])));
					if(i!=uppernode) ainv[i*(uppernode+1)+(i+1)]=(-1./(dz_[i]))*((k_np1_m_ip12[i]/(dz_up[i])));

					//Correct upper and lower diagonals in case of Dirichlet
					if(aTopBC==DIRICHLET && i==uppernode) {
						ainv[(i-1)*(uppernode+1)+i]=0.;
						ainv[i*(uppernode+1)+(i-1)]=0.;
					}
					if(aBottomBC==DIRICHLET && i==lowernode) {
						ainv[(i+1)*(uppernode+1)+i]=0.;
						ainv[i*(uppernode+1)+(i+1)]=0.;
					}
				}

				//This part is for the DGTSV or TDMA solver, that uses the fact that A is a tridiagonal matrix, so we only have to specify the diagonals and subdiagonals.
				if(ActiveSolver==DGTSV || ActiveSolver==TDMA ) {
					if(i==j) {
						//Set up the matrix diagonal
						ad[i]=(1./dt)*C[i];

						//The following two lines assume Neumann boundary conditions (for upper and lowernode, one of the terms drop out). If Dirichlet is used, this will be corrected later.
						if(i!=lowernode) ad[i]+=(1./(dz_[i]))*((k_np1_m_im12[i]/(dz_down[i])));
						if(i!=uppernode) ad[i]+=(1./(dz_[i]))*((k_np1_m_ip12[i]/(dz_up[i])));

						//Correct diagonal in case of Dirichlet
						if(aTopBC==DIRICHLET && i==uppernode) {
							ad[i]=1.;
						}
						if(aBottomBC==DIRICHLET && i==lowernode) {
							ad[i]=1.;
						}

						//Set up the matrix upper and lower diagonals
						if(i!=lowernode) adl[i-1]=-(1./(dz_[i]))*((k_np1_m_im12[i]/(dz_down[i])));
						if(i!=uppernode) adu[i]=-(1./(dz_[i]))*((k_np1_m_ip12[i]/(dz_up[i])));

						//Correct diagonals in case of Dirichlet
						if(aTopBC==DIRICHLET && i==uppernode) {
							adu[i-1]=0.;
							adl[i-1]=0.;
						}
						if(aBottomBC==DIRICHLET && i==lowernode) {
							adu[i]=0.;
							adl[i]=0.;
						}
					}
				}

				//We copy here the matrix to the ainv, which is passed to the SVD-routine later on. This ainv is altered externally, that's why we need a copy.
				//ainv[j*(uppernode+1)+i]=a[i][j];

				//Determine R:
				term_up[i]=0.;
				term_down[i]=0.;

				//Fill R.H.S. vector
				//Note: the gravity term is not explicitly in Celia et al (1990). It is just z[i], as pressure head should already be scaled by rho_water * g. Then it is taken outside the nabla, by using the chain rule.
				if(i==uppernode) {
					if(aTopBC==NEUMANN) {		//Neumann, following Equation 4 in McCord, WRR (1991).
						term_up[i]=(TopFluxRate)*dz_up[i] - cos_sl*(dz_up[i]*k_np1_m_ip12[i]);
					} else {	//Dirichlet
						term_up[i]=0.;
					}
				} else {
					if(activelayer[i+1]==true) {
						term_up[i]=k_np1_m_ip12[i]*(h_np1_m[i+1]-h_np1_m[i]);
					} else {
						//Analogue to Neumann at top:
						term_up[i]=(0.)*dz_up[i] - cos_sl*(dz_up[i]*k_np1_m_ip12[i]);
					}
				}
				if(i==lowernode) {
					if(aBottomBC == NEUMANN) {	//Neumann, following Equation 4 in McCord, WRR (1991).
						term_down[i]=(BottomFluxRate)*dz_down[i] - cos_sl*(dz_down[i]*k_np1_m_im12[i]);
					} else {			//Dirichlet
						term_down[i]=0.;
					}
				} else {
					if(activelayer[i-1]==true) {
						term_down[i]=k_np1_m_im12[i]*(h_np1_m[i]-h_np1_m[i-1]);
					} else {
						//Analogue to Neumann at top:
						term_down[i]=(0.)*dz_down[i] - cos_sl*(dz_down[i]*k_np1_m_im12[i]);
					}
				}

				//RHS eq. 17 in Celia et al. (1990):
				r_mpfd[i]=(1./(dz_[i]))*((term_up[i]/dz_up[i])-(term_down[i]/dz_down[i])) + cos_sl*((k_np1_m_ip12[i]-k_np1_m_im12[i])/(dz_[i])) - (1./dt)*((theta_np1_m[i]-theta_n[i]) + (theta_i_np1_m[i]-theta_i_n[i])*(Constants::density_ice/Constants::density_water)) + s[i];

				// r_mpfd is an approximation of how far one is away from the solution. So in case of Dirichlet boundaries, we are *at* the solution:
				if(aTopBC==DIRICHLET) r_mpfd[uppernode]=0.;
				if(aBottomBC==DIRICHLET) r_mpfd[lowernode]=0.;

				r_mpfd2[i]=r_mpfd[i];			// We make a copy for use with DGTSV and TDMA solvers.
				if(WriteOutNumerics_Level3==true) {
					std::cout << "SOLVER: i=" << i << std::scientific << " - r_mpfd=" << r_mpfd[i] << " term_up=" << term_up[i] << " term_down=" << term_down[i] << " a=" << ainv[i*(uppernode+1)+i]/*a[i][i]*/ << "adl=" << adl[i] << " adu=" << adu[i] << " [" << K[i] << " - " << C[i] << "]\n" << std::fixed;
				}
			}

			//Before solving the system of equations, reset convergence tracking variables:
			track_accuracy_h=0.;
			track_accuracy_theta=0.;
			track_trigger_layer_accuracy=-1;
			accuracy=-1.;				//-1 is a flag. accuracy can only be positive, so when it is negative, we know that no layer did NOT converged yet.
			trigger_layer_accuracy=-1;		//-1 is a flag. when it is negative, we know that no layer was NOT converged yet.
			int trigger_layer_blowup=-1;		//-1 is a flag. when it is negative, we know that no layer was NOT converged yet.
			max_delta_h=0.;
			boolConvergence=true;			//We initialize it as true, and set it to false when necessary.
			mass2=0.;

			//Now call the designated solver.
			if(solver_result==0) {
				if (ActiveSolver==TDMA) {
					// Note: TDMA is very rapid, but has the problem that when elements in the matrix differ order of magnitudes, rounding errors can occur that destroy accuracy.
					// For this reason, it is better to use DGTSV solver, which does partial pivoting to prevent this. See: http://en.wikipedia.org/wiki/Pivot_element#Partial_and_complete_pivoting
					const int matrixdimensions=(uppernode-lowernode)+1;
					solver_result=TDMASolver(matrixdimensions, &adl[0], &ad[0], &adu[0], &r_mpfd[0], &r_mpfd2[0]);
				}

				if(ActiveSolver==DGTSV) {
#ifdef CLAPACK
					// Solver for Tridiagonal matrices, with partial pivoting.
					int info=0;
					const int matrixdimensions=(uppernode-lowernode)+1;
					const int vectordimensions=1;
					dgtsv_( (integer*) &matrixdimensions, (integer*) &vectordimensions, &adl[0], &ad[0], &adu[0], &r_mpfd2[0], (integer*) &matrixdimensions, (integer*) &info );

					if(info!=0) {
						//= 0: successful exit
						//< 0: if INFO = -i, the i-th argument had an illegal value
						//> 0: if INFO = i, U(i,i) is exactly zero, and the solution
						//    has not been computed.  The factorization has not been
						//    completed unless i = N.
						if(AllowSwitchSolver==true) {
							if(WriteOutNumerics_Level0==true) std::cout << "ERROR in ReSolver1d.cc: DGTSV failed [info = " << info << "]. Trying DGESVD/DGESDD...\n";
							ActiveSolver=DGESVD;
						} else {
							if(WriteOutNumerics_Level0==true) std::cout << "ERROR in ReSolver1d.cc: DGTSV failed [info = " << info << "]. Trying with smaller time step...\n";
							solver_result=-1;
						}
					}
#else
					throw InvalidArgumentException("you cannot use solver DGTSV when libraries BLAS and LAPACK are not installed. Either install these libraries, or choose solver TDMA", AT);
#endif
				}

				if(ActiveSolver==DGESVD) {
#ifdef CLAPACK
					//Do Moore-Penrose matrix inversion, using singular value decomposition (SVD), so we can write: H = A' * R
					solver_result=pinv((uppernode-lowernode)+1, (uppernode-lowernode)+1, (uppernode-lowernode)+1, &ainv[0]);
#else
					throw InvalidArgumentException("you cannot use solver DGESVD when libraries BLAS and LAPACK are not installed. Either install these libraries, or choose solver TDMA", AT);
#endif
				}


				//Apply new iteration solution

				//This is a little bit complicated. The problem started when we did soil freezing. If then suddenly an isnan is detected somewhere in the model domain, some part of the soil is already through the phasechange function, other parts not (maybe).
				//It is difficult to revert this soil freezing, so therefore, we need first to loop over i to determine the complete solution vector delta_h, and then an other loop over i to apply the new solution.
				//However, if a proper way to revert soil freezing is made, this extra loop can be removed.
				for (i = uppernode; i >= lowernode; i--) {
					//Determine delta h:
					if(ActiveSolver==DGESVD) {
						delta_h[memstate%nmemstates][i]=0.;
						//Note: after inverting, ainv is non tridiagonal, so we have to loop over all elements.
						for (k = uppernode; k >= lowernode; k--) {
						//for (k = MIN(uppernode, i+1); k >= MAX(lowernode, i-1); k--) {
							delta_h[memstate%nmemstates][i]+=ainv[i*(uppernode+1)+k]*r_mpfd[k];
						}
					} else {	//In case of DGTSV, solution is returned in r_mpfd2, overwriting original content.
						delta_h[memstate%nmemstates][i]=r_mpfd2[i];
					}
					if(isnan(delta_h[memstate%nmemstates][i])==true || isinf(delta_h[memstate%nmemstates][i])==true) {
						solver_result=-1;
					}
				}
			}


			//Apply Dirichlet BCs:
			if(aTopBC==DIRICHLET) {
				h_np1_mp1[uppernode]=htop;
				delta_h[memstate%nmemstates][uppernode]=0.;
				delta_theta[uppernode]=0.;
			}
			if(aBottomBC==DIRICHLET) {
				h_np1_mp1[lowernode]=hbottom;
				delta_h[memstate%nmemstates][lowernode]=0.;
				delta_theta[lowernode]=0.;
			}

			for (i = uppernode; i >= lowernode; i--) {
				if(activelayer[i]==true) {
					//Keep track of the maximum delta h, to detect possible model blow-ups.
					if(fabs(delta_h[memstate%nmemstates][i])>max_delta_h) {	// If change is too big and we are allowed to do a rewind, don't check for accuracy
						//delta_h[memstate%nmemstates][i]=0.;
						trigger_layer_blowup=i;
						max_delta_h=fabs(delta_h[memstate%nmemstates][i]);
						h_np1_mp1[i]=h_np1_m[i];
						theta_np1_mp1[i]=theta_np1_m[i];
						delta_theta_i[i]=0.;
						delta_theta[i]=1E10;			//Set delta_theta[i] to any value, to make sure the convergence test will fail.
					}

					//if(not(max_delta_h>MAX_ALLOWED_DELTA_H) || niter>MAX_ITER+1) {	//If it is, there is a big chance the call to fromHtoTHETA will fail because of floating point exceptions (overflows). In this case, we will force a rewind later on, so the solution does not matter anymore.
														//The second clause means we cannot increase time step anymore, so we should just try the solution.
					if(solver_result!=-1) {
						//Apply solution
						h_np1_mp1[i]=h_np1_m[i]+delta_h[memstate%nmemstates][i];

						//Calculate theta
						theta_np1_mp1[i]=fromHtoTHETAforICE(h_np1_mp1[i], theta_r[i], theta_s[i], alpha[i], m[i], n[i], Sc[i], h_e[i], theta_i_np1_m[i]);

						//Calculate temperature change of soil layers to reflect heat advected by the flowing water
						if(i<nsoillayers_richardssolver) {
							//Calculate the fluxes from above and below for this layer
							const double tmp_flux_above = (i<toplayer-1) ? (((((h_np1_m[i+1]+delta_h[memstate%nmemstates][i+1])-(h_np1_m[i]+delta_h[memstate%nmemstates][i]))/dz_up[i])+cos_sl)*k_np1_m_ip12[i]*dt) : 0;			//Units: [m^3/m^2]
							const double tmp_flux_below = (i>0) ? (((((h_np1_m[i]+delta_h[memstate%nmemstates][i])-(h_np1_m[i-1]+delta_h[memstate%nmemstates][i-1]))/dz_up[i-1])+cos_sl)*k_np1_m_ip12[i-1]*dt) : 0;				//Units: [m^3/m^2]

							//Calculate intermediate state variables of this layer
							const double tmp_theta_air = 1. - theta_i_n[i] - (theta_np1_mp1[i] + (theta_i_np1_m[i]-theta_i_n[i])*(Constants::density_ice/Constants::density_water)) - EMS[SnowpackElement[i]].theta[SOIL];					//Units: [m^3 m^-3]
							const double tmp_rho = (Constants::density_ice * theta_i_n[i] + Constants::density_water * (theta_np1_mp1[i] + (theta_i_np1_m[i]-theta_i_n[i])*(Constants::density_ice/Constants::density_water)) + EMS[SnowpackElement[i]].soil[SOIL_RHO] * EMS[SnowpackElement[i]].theta[SOIL]);	//Units: [kg m-3]
							const double tmp_c_p = (Constants::density_air * tmp_theta_air * Constants::specific_heat_air							//Units: [J kg-1 K-1]
										+ Constants::density_ice * theta_i_n[i] * Constants::specific_heat_ice
										+ Constants::density_water * (theta_np1_mp1[i] + (theta_i_np1_m[i]-theta_i_n[i])*(Constants::density_ice/Constants::density_water)) * Constants::specific_heat_water
										+ EMS[SnowpackElement[i]].soil[SOIL_RHO] * EMS[SnowpackElement[i]].theta[SOIL] * EMS[SnowpackElement[i]].soil[SOIL_C]
										) / tmp_rho;
							delta_Te_adv_i[i]=0.;
							if (tmp_flux_above>0.) {	//Positve flux from above (= influx in current layer)
								//Advected heat
								const double tmp_adv_heat = ((EMS[SnowpackElement[i+1]].Te + delta_Te_adv[i+1] + delta_Te[i+1]) - (EMS[SnowpackElement[i]].Te + delta_Te_adv[i] + delta_Te[i])) * Constants::density_water * tmp_flux_above * Constants::specific_heat_water;	//Units [J/m^2]
								delta_Te_adv_i[i] = (tmp_adv_heat) / (tmp_c_p * tmp_rho * EMS[SnowpackElement[i]].L);
							}
							if (tmp_flux_below<0.) {	//Negative flux from below (=influx in current layer)
								//Advected heat
								const double tmp_adv_heat = ((EMS[SnowpackElement[i-1]].Te + delta_Te_adv[i-1] + delta_Te[i-1]) - (EMS[SnowpackElement[i]].Te + delta_Te_adv[i] + delta_Te[i])) * Constants::density_water * (-1.*tmp_flux_below) * Constants::specific_heat_water;	//Units [J/m^2]
								//In rare cases, we may have inflow from above AND below, so we add (+=) the temperature change due to heat advection
								delta_Te_adv_i[i] += (tmp_adv_heat) / (tmp_c_p * tmp_rho * EMS[SnowpackElement[i]].L);
							}

							//Repartition ice/water based on new head
							if(AllowSoilFreezing==true) {
								size_t BS_iter=0;			//Counting the number of iterations
								const double hw0=h_np1_mp1[i];
								T_melt[i]=T_0+((Constants::g*T_0)/delF)*hw0;
								// Bisection-Secant method, see wikipedia: http://en.wikipedia.org/wiki/False_position_method
								//   fromHtoTHETA(hw0+(Constants::lh_fusion/(Constants::g*T_melt[i]))*(EMS[SnowpackElement[i]].Te-T_melt[i]), theta_r[i], theta_s[i], alpha[i], m[i], n[i], Sc[i], h_e[i])
								//      +
								//   (theta_i_np1_mp1[i]*(Constants::density_ice/Constants::density_water))
								//      -
								//   fromHtoTHETA(hw0, theta_r[i], theta_s[i], alpha[i], m[i], n[i], Sc[i], h_e[i]);
								//      = 0.
								// Solving this equation for theta_i_np1_mp1[i] (which influences theta_np1_mp1 and Te)

								// So the new liquid water content basically is the same equation, but we have to adapt EMS[SnowpackElement[i]].Te to the amount of ice we create (delta_i).
								if((theta_i_np1_m[i] > 0. && (EMS[SnowpackElement[i]].Te + delta_Te_adv[i] + delta_Te_adv_i[i] + delta_Te[i]) > T_melt[i]) || (EMS[SnowpackElement[i]].Te + delta_Te_adv[i] + delta_Te_adv_i[i] + delta_Te[i]) < T_melt[i]) {

									if(WriteOutNumerics_Level2==true) {
										const double tmp_T = EMS[SnowpackElement[i]].Te + delta_Te_adv[i] + delta_Te_adv_i[i] + delta_Te[i] + delta_Te_i[i];
										std::cout << "BEFORE [" << i << std::fixed << std::setprecision(15) << "]; theta_w: " << theta_np1_mp1[i] << " theta_i_np1_m: " << theta_i_np1_m[i] << " theta_s: " << theta_s[i] << std::setprecision(3) << "  T: " << tmp_T << std::setprecision(8) << "  rho: " << tmp_rho << "  cp: " << tmp_c_p << " ColdC: " << tmp_rho * tmp_c_p * tmp_T * EMS[SnowpackElement[i]].L << "\n" << std::setprecision(6);
									}

									//Determine maximum possible change in ice content, which should be between 0, and theta_water > theta_d (all possible water freezes). Then maximum ice content is determined based on the temperature difference between element and T_melt.
									//const double max_delta_ice=(MIN((theta_np1_mp1[i]-theta_d[i])*(Constants::density_ice/Constants::density_water), MAX(0., T_melt[i]-(EMS[SnowpackElement[i]].Te/*+delta_Te*/)) * ((EMS[SnowpackElement[i]].c[TEMPERATURE] * EMS[i].Rho) / ( Constants::density_ice * Constants::lh_fusion ))));
									double max_delta_ice;
									if((EMS[SnowpackElement[i]].Te + delta_Te_adv[i] + delta_Te_adv_i[i] + delta_Te[i]) > T_melt[i]) {
										// Melt: either all ice will disappear, or a fraction based on available energy
										max_delta_ice=-1.*theta_i_n[i];
									} else {
										// Freeze: either all available water will freeze, or a fraction based on available energy.
										max_delta_ice=(theta_np1_mp1[i]-0.)*(Constants::density_water/Constants::density_ice);
									}

									bool BS_converged=false;
									double ak=0., bk=0., ck=0.;	//These are values for changes in ice content.
									double delta_Te_ak=0., delta_Te_bk=0., delta_Te_ck=0., delta_w_ak=0., delta_w_bk=0., delta_w_ck=0.;
									double ck1=0, delta_Te_ck1=0., delta_w_ck1=0.;
									if(max_delta_ice>0.) {
										ak=0.;
										bk=max_delta_ice;
									} else {
										ak=max_delta_ice;
										bk=0.;
									}
									// Deal with special cases:
									// 1) So much energy available that all ice will melt (note: this case will not be properly solved by Bisection-Secant method.)
									if((T_melt[i]-(EMS[SnowpackElement[i]].Te + delta_Te_adv[i] + delta_Te_adv_i[i] + delta_Te[i])) * ((tmp_c_p * tmp_rho) / ( Constants::density_ice * Constants::lh_fusion )) < -1.*theta_i_n[i] && BS_converged==false) {
										ck=-1.*theta_i_np1_m[i];
										delta_w_ck=-1.*(ck*(Constants::density_ice/Constants::density_water));
										delta_Te_ck=((theta_i_np1_m[i] - theta_i_n[i]) + ck) / ((tmp_c_p * tmp_rho) / ( Constants::density_ice * Constants::lh_fusion ));	//Change in element temperature associated with change in ice content
										if(WriteOutNumerics_Level3==true) {
											const double tmp_T = EMS[SnowpackElement[i]].Te + delta_Te_adv[i] + delta_Te_adv_i[i] + delta_Te[i] + delta_Te_i[i] + delta_Te_ck;
											std::cout << "BS_ITER [" << BS_iter << std::scientific << "], case 2: a=" << ak << " b=" << bk << " c=" << ck << " (max: " << max_delta_ice << ") " << delta_w_ck << " " << tmp_T << " " << T_melt[i] << ": fa: " << (delta_w_ak + ak*(Constants::density_ice/Constants::density_water)) << " fb: " << (delta_w_bk + bk*(Constants::density_ice/Constants::density_water)) << " fc: " << (delta_w_ck + ck*(Constants::density_ice/Constants::density_water)) << "\n" << std::fixed;
										}
										BS_converged=true;
									}
									// 2) Very small temperature difference or very small possible change in ice content
									if(fabs(ak-bk)<SF_epsilon && BS_converged==false) {
										// In this case it is possible that we should melt some ice in order to prevent theta[WATER] to get negative (drainage case):
										ck=0.;
										if(theta_np1_mp1[i]<0.) {
											delta_w_ck=-1.*theta_np1_mp1[i];					//Make sure water gets 0.
											ck=theta_np1_mp1[i]*(Constants::density_water/Constants::density_ice);	//Necessary change in theta[ICE]
											delta_Te_ck=((theta_i_np1_m[i] - theta_i_n[i]) + ck) / ((tmp_c_p * tmp_rho) / ( Constants::density_ice * Constants::lh_fusion ));	//Change in element temperature associated with change in ice content
										}
										if(WriteOutNumerics_Level3==true) {
											const double tmp_T = EMS[SnowpackElement[i]].Te + delta_Te_adv[i] + delta_Te_adv_i[i] + delta_Te[i] + delta_Te_i[i] + delta_Te_ck;
											std::cout << "BS_ITER [" << BS_iter << std::scientific << "], case 1: a=" << ak << " b=" << bk << " c=" << ck << " (max: " << max_delta_ice << ") " << delta_w_ck << " " << tmp_T << " " << T_melt[i] << ": fa: " << (delta_w_ak + ak*(Constants::density_ice/Constants::density_water)) << " fb: " << (delta_w_bk + bk*(Constants::density_ice/Constants::density_water)) << " fc: " << (delta_w_ck + ck*(Constants::density_ice/Constants::density_water)) << "\n" << std::fixed;
										}
										BS_converged=true;
									}
									while (BS_converged==false && BS_iter < BS_MAX_ITER) {
										BS_iter++;
										delta_Te_ak=((theta_i_np1_m[i] - theta_i_n[i]) + ak) / ((tmp_c_p * tmp_rho) / ( Constants::density_ice * Constants::lh_fusion ));			//Change in element temperature associated with change in ice content
										delta_Te_bk=((theta_i_np1_m[i] - theta_i_n[i]) + bk) / ((tmp_c_p * tmp_rho) / ( Constants::density_ice * Constants::lh_fusion ));			//Change in element temperature associated with change in ice content
										delta_w_ak=(fromHtoTHETA(hw0+(Constants::lh_fusion/(Constants::g*T_melt[i]))*MIN(0., (EMS[SnowpackElement[i]].Te + delta_Te_adv[i] + delta_Te_adv_i[i] + delta_Te[i] + delta_Te_ak)-T_melt[i]), theta_r[i], theta_s[i], alpha[i], m[i], n[i], Sc[i], h_e[i])) - theta_np1_mp1[i];
										delta_w_bk=(fromHtoTHETA(hw0+(Constants::lh_fusion/(Constants::g*T_melt[i]))*MIN(0., (EMS[SnowpackElement[i]].Te + delta_Te_adv[i] + delta_Te_adv_i[i] + delta_Te[i] + delta_Te_bk)-T_melt[i]), theta_r[i], theta_s[i], alpha[i], m[i], n[i], Sc[i], h_e[i])) - theta_np1_mp1[i];
										//Now calculate bisect
										ck1=(ak+bk)/2.;
										delta_Te_ck1=((theta_i_np1_m[i] - theta_i_n[i]) + ck1) / ((tmp_c_p * tmp_rho) / ( Constants::density_ice * Constants::lh_fusion ));			//Change in element temperature associated with change in ice content
										delta_w_ck1=(fromHtoTHETA(hw0+(Constants::lh_fusion/(Constants::g*T_melt[i]))*MIN(0., (EMS[SnowpackElement[i]].Te + delta_Te_adv[i] + delta_Te_adv_i[i] + delta_Te[i] + delta_Te_ck1)-T_melt[i]), theta_r[i], theta_s[i], alpha[i], m[i], n[i], Sc[i], h_e[i])) - theta_np1_mp1[i];
										//Now check secant
										ck=((delta_w_bk + bk*(Constants::density_ice/Constants::density_water))*ak  -  (delta_w_ak + ak*(Constants::density_ice/Constants::density_water))*bk)  /  ((delta_w_bk + bk*(Constants::density_ice/Constants::density_water)) - (delta_w_ak + ak*(Constants::density_ice/Constants::density_water)));
										delta_Te_ck=((theta_i_np1_m[i] - theta_i_n[i]) + ck) / ((tmp_c_p * tmp_rho) / ( Constants::density_ice * Constants::lh_fusion ));			//Change in element temperature associated with change in ice content
										delta_w_ck=(fromHtoTHETA(hw0+(Constants::lh_fusion/(Constants::g*T_melt[i]))*MIN(0., (EMS[SnowpackElement[i]].Te + delta_Te_adv[i] + delta_Te_adv_i[i] + delta_Te[i] + delta_Te_ck)-T_melt[i]), theta_r[i], theta_s[i], alpha[i], m[i], n[i], Sc[i], h_e[i])) - theta_np1_mp1[i];
										//Now check if bisect or secant is a better approximation
										if(fabs(delta_w_ck + ck*(Constants::density_ice/Constants::density_water))>fabs(delta_w_ck1+ck1*(Constants::density_ice/Constants::density_water))) {
											ck=ck1;
											delta_Te_ck=delta_Te_ck1;
											delta_w_ck=delta_w_ck1;
										}
										if(WriteOutNumerics_Level3==true) {
											const double tmp_T = EMS[SnowpackElement[i]].Te + delta_Te_adv[i] + delta_Te_adv_i[i] + delta_Te[i] + delta_Te_i[i] + delta_Te_ck;
											std::cout << "BS_ITER [" << BS_iter << std::scientific << "]: a=" << ak << " b=" << bk << " c=" << ck << " (max: " << max_delta_ice << ") " << delta_w_ck << " " << tmp_T << " " << T_melt[i] << ": fa: " << (delta_w_ak + ak*(Constants::density_ice/Constants::density_water)) << " fb: " << (delta_w_bk + bk*(Constants::density_ice/Constants::density_water)) << " fc: " << (delta_w_ck + ck*(Constants::density_ice/Constants::density_water)) << "\n" << std::fixed;
										}
										//Now check if convergence is achieved
										if(fabs(delta_w_ck + ck*(Constants::density_ice/Constants::density_water)) < SF_epsilon) {
											delta_w_ck=-1.*(ck*(Constants::density_ice/Constants::density_water));	//Make delta in water equal to ice, so we keep mass-balance.
											BS_converged=true;
										} else if(fabs(delta_w_ak + ak*(Constants::density_ice/Constants::density_water)) < SF_epsilon) {
											ck=ak;
											delta_w_ck=-1.*(ck*(Constants::density_ice/Constants::density_water));	//Make delta in water equal to ice, so we keep mass-balance.
											delta_Te_ck=delta_Te_ak;
											BS_converged=true;
										} else if(fabs(delta_w_bk + bk*(Constants::density_ice/Constants::density_water)) < SF_epsilon) {
											ck=bk;
											delta_w_ck=-1.*(ck*(Constants::density_ice/Constants::density_water));	//Make delta in water equal to ice, so we keep mass-balance.
											delta_Te_ck=delta_Te_bk;
											BS_converged=true;
										} else {
											//And determine whether to update the left or right point
											if((delta_w_ck + ck*(Constants::density_ice/Constants::density_water)) * (delta_w_ak + ak*(Constants::density_ice/Constants::density_water)) > 0.) {	//Multiply to check if same sign
												ak=ck;
											} else {
												bk=ck;
											}
										}
									}
									if(BS_converged==false) {
										if(WriteOutNumerics_Level0==true) std::cout << "[W] ReSolver1d.cc: Bisect-Secant method failed to converge in soil freezing with dt = " << dt << ".\n";
										if(WriteOutNumerics_Level1==true) {
											const double tmp_T = EMS[SnowpackElement[i]].Te + delta_Te_adv[i] + delta_Te_adv_i[i] + delta_Te[i] + delta_Te_i[i] + delta_Te_ck;
											std::cout << "  -- BS_ITER [" << BS_iter << std::scientific << "]: a=" << ak << " b=" << bk << " c=" << ck << " (max: " << max_delta_ice << ") " << delta_w_ck << " " << tmp_T << " " << T_melt[i] << ": fa: " << (delta_w_ak + ak*(Constants::density_ice/Constants::density_water)) << " fb: " << (delta_w_bk + bk*(Constants::density_ice/Constants::density_water)) << " fc: " << (delta_w_ck + ck*(Constants::density_ice/Constants::density_water)) << "\n" << std::fixed;
											std::cout << "  -- " << std::setprecision(15) << T_melt[i] << " " << EMS[SnowpackElement[i]].Te << " " << delta_Te_adv[i] << " " << delta_Te_adv_i[i] << " " << delta_Te[i] << "   " << EMS[SnowpackElement[i]].theta[WATER] << " " << EMS[SnowpackElement[i]].theta[ICE] << "\n" << std::setprecision(6);
										}
										max_delta_h=2.*MAX_ALLOWED_DELTA_H;
										solver_result=-1;
									} else {
										//Final solution
										const double tmp_delta_i=ck;
										const double tmp_delta_w=delta_w_ck;
										const double tmp_delta_Te=delta_Te_ck;
										//Apply final solution
										delta_Te_i[i]=tmp_delta_Te;
										theta_i_np1_mp1[i]=theta_i_np1_m[i]+tmp_delta_i;
										theta_np1_mp1[i]+=tmp_delta_w;
									}
								} else {
									theta_i_np1_mp1[i]=0.;
									theta_np1_mp1[i]=fromHtoTHETAforICE(h_np1_mp1[i], theta_r[i], theta_s[i], alpha[i], m[i], n[i], Sc[i], h_e[i], 0.);
								}
								//Update BS-solver statistics
								bs_stats_totiter+=BS_iter;
								if(BS_iter>bs_stats_maxiter) bs_stats_maxiter=BS_iter;
								if(WriteOutNumerics_Level2==true) 
									std::cout << "AFTER [" << i << std::setprecision(15) << "]: theta_w: " << theta_np1_mp1[i] << " theta_i_np1_m: " << theta_i_np1_mp1[i] << " theta_s:" << theta_s[i] << std::setprecision(3) << "  T: " << EMS[SnowpackElement[i]].Te + delta_Te_adv[i] + delta_Te_adv_i[i] + delta_Te[i] + delta_Te_i[i] << " (niter=" << BS_iter << ")\n" << std::setprecision(6);
							} //END OF REPARTITIONING ICE/WATER
						}

						delta_theta[i]=theta_np1_mp1[i]-theta_np1_m[i];
						delta_theta_i[i]=theta_i_np1_mp1[i]-theta_i_np1_m[i];
					} else {
						// Solver failed, trigger rewind
						max_delta_h=2.*MAX_ALLOWED_DELTA_H;
					}
				} else {
					theta_np1_mp1[i]=theta_n[i];
					h_np1_mp1[i]=h_n[i];
					delta_theta[i]=0.;
					delta_theta_i[i]=0.;
					delta_h[memstate%nmemstates][i]=0.;
				}

				//Update mass balance
				mass2+=(theta_np1_mp1[i]+(theta_i_np1_mp1[i]*(Constants::density_ice/Constants::density_water)))*dz[i];

				if(WriteOutNumerics_Level2==true) {
					std::cout << "ITER: " << niter << " i: " << i << std::scientific << std::setprecision(10) << " ---  h1: " << h_np1_m[i] << " d_h: " << delta_h[memstate%nmemstates][i] << " h2: " << h_np1_mp1[i] << std::setprecision(12) << " --- theta1: " << theta_np1_m[i] << " d_theta: " << delta_theta[i] << " theta2: " << theta_np1_mp1[i] << "\n" << std::setprecision(6) << std::fixed;
				}


				// Note: boundaries are not considered for determination of the accuracy in case of Dirichlet (of course!).

				//Absolute accuracy in h: This seems to produce the best (in sense of most stable) behaviour.
				//Note: because we want to be able to very accurately determine the flux over the snow-soil boundary (for MS_SNOW_RUNOFF) and model boundaries (for MS_SOIL_RUNOFF),
				//we ALWAYS have to assess the accuracy in head in this region! If we don't do this, then in case of dry soil layers, the estimated pressure head can be quite
				//inaccurate, leading to a completely wrong estimation of these fluxes!
				if(Se[i]>convergencecriterionthreshold || i==nsoillayers_richardssolver-1 || i==nsoillayers_richardssolver || i==lowernode || i==lowernode-1) {
					if ((i!=lowernode || aBottomBC==NEUMANN) && (i!=uppernode || aTopBC==NEUMANN)) {
						//First check general accuarcy:
						if(fabs(delta_h[memstate%nmemstates][i])>track_accuracy_h) {
							track_trigger_layer_accuracy=i;
							track_accuracy_h=fabs(delta_h[memstate%nmemstates][i]);
						}
						//Now check against convergence criterion:
						if(fabs(delta_h[memstate%nmemstates][i])>REQUIRED_ACCURACY_H) {
							trigger_layer_accuracy=i;
							accuracy=fabs(delta_h[memstate%nmemstates][i]);
						}
					}
				}

				//Absolute accuracy in theta. This doesn't behave stable, especially during saturated soil conditions, when the accuracy is set too low.
				//See Huang (1996), which proposes this, and also discusses the need for a higher accuracy:
				//if(not(Se[i]>convergencecriterionthreshold || i==uppernode || i==lowernode || i==nsoillayers_richardssolver-1 || i==nsoillayers_richardssolver)) {
				if(not(Se[i]>convergencecriterionthreshold || i==nsoillayers_richardssolver-1 || i==nsoillayers_richardssolver || i==lowernode || i==lowernode-1)) {
					if ((i!=lowernode || aBottomBC==NEUMANN) && (i!=uppernode || aTopBC==NEUMANN)) {
						//First check general accuarcy:
						if(fabs(delta_theta[i]+delta_theta_i[i]*(Constants::density_ice/Constants::density_water))>track_accuracy_theta) {
							track_trigger_layer_accuracy=i;
							track_accuracy_theta=fabs(delta_theta[i]+delta_theta_i[i]*(Constants::density_ice/Constants::density_water));
						}
						//Now check against convergence criterion:
						if( fabs(delta_theta[i]+delta_theta_i[i]*(Constants::density_ice/Constants::density_water)) > REQUIRED_ACCURACY_THETA ) {
							trigger_layer_accuracy=i;
							accuracy=fabs(delta_theta[i]+delta_theta_i[i]*(Constants::density_ice/Constants::density_water));
						}
					}
				}
			}

			//Apply boundary conditions
			if(aTopBC==DIRICHLET) {
				h_np1_mp1[uppernode]=htop;		//Dirichlet
			}
			if(aBottomBC==DIRICHLET) {
				h_np1_mp1[lowernode]=hbottom;		//Dirichlet
			}

			//Check mass balance:
			//-- Calculate mass change:
			massbalanceerror=mass1-mass2;
			//-- Determine top and bottom flux:
			double tmp_mb_topflux=0.;
			double tmp_mb_bottomflux=0.;
			if(aTopBC==NEUMANN) {		//If we use Neumann, the massbalance should incorporate the applied TopFluxRate:
				tmp_mb_topflux=TopFluxRate*dt;
			} else {			//Else when using Dirichlet, we should estimate the influx: (Note that basically with Dirichlet, the change of theta in the element is 0., so the influx in the model domain is equal to the flux from the upper element to the one below.)
				tmp_mb_topflux=((theta_np1_mp1[uppernode]+theta_i_np1_mp1[uppernode]*(Constants::density_ice/Constants::density_water))-(theta_n[uppernode] + theta_i_n[uppernode]*(Constants::density_ice/Constants::density_water)))*dz[uppernode] + ((((h_np1_mp1[uppernode]-h_np1_mp1[uppernode-1])/dz_down[uppernode])+cos_sl)*k_np1_m_im12[uppernode]*dt);
			}
			if(aBottomBC==NEUMANN) {	//If we use Neumann, the massbalance should incorporate the applied BottomFluxRate:
				tmp_mb_bottomflux=BottomFluxRate*dt;
			} else {			//Else when using Dirichlet, we should estimate the outflux: (Note that basically with Dirichlet, the change of theta in the element is 0., so the outflux in the model domain is equal to the flux from the element above the lowest one to the lowest one.)
				tmp_mb_bottomflux=-1.*(((theta_np1_mp1[lowernode]+theta_i_np1_mp1[lowernode]*(Constants::density_ice/Constants::density_water))-(theta_n[lowernode] + theta_i_n[lowernode]*(Constants::density_ice/Constants::density_water)))*dz[lowernode]-((((h_np1_mp1[lowernode+1]-h_np1_mp1[lowernode])/dz_up[lowernode])+cos_sl)*k_np1_m_ip12[lowernode]*dt));
			}
			massbalanceerror+=tmp_mb_topflux;		//Add topflux (note: topflux>0. means influx)
			massbalanceerror-=tmp_mb_bottomflux;		//Substract bottomflufx (note: bottomflux>0. means outflux)
			massbalanceerror+=totalsourcetermflux*dt;	//Add the sink/source term flux.
			if(WriteOutNumerics_Level2==true) printf("MASSBALANCETEST: mass1 %.8E    mass2 %.8E    topflux %.8E (%.8E)  bottomflux %.8E (%.8E) sourceflux %.8E    delta %.8E\n", mass1, mass2, tmp_mb_topflux, ((theta_np1_mp1[uppernode]+theta_i_np1_mp1[uppernode]*(Constants::density_ice/Constants::density_water))-(theta_n[uppernode] + theta_i_n[uppernode]*(Constants::density_ice/Constants::density_water)))*dz[uppernode] + ((((h_np1_m[uppernode]-h_np1_m[uppernode-1])/dz_down[uppernode])+1.)*k_np1_m_im12[uppernode]*dt), tmp_mb_bottomflux, -1.*(((theta_np1_mp1[lowernode]+theta_i_np1_mp1[lowernode]*(Constants::density_ice/Constants::density_water))-(theta_n[lowernode] + theta_i_n[lowernode]*(Constants::density_ice/Constants::density_water)))*dz[lowernode]-((((h_np1_m[lowernode+1]-h_np1_m[lowernode])/dz_up[lowernode])+cos_sl)*k_np1_m_ip12[lowernode]*dt)), totalsourcetermflux*dt, massbalanceerror);

			//Make sure to trigger a rewind by making max_delta_h very large in case the mass balance is violated or change in head are too large.
			if(fabs(massbalanceerror)>1E-1 || max_delta_h>MAX_ALLOWED_DELTA_H) {
				max_delta_h=2.*MAX_ALLOWED_DELTA_H;
			}

			if (accuracy > 1E-20 || fabs(massbalanceerror)>maxallowedmassbalanceerror) {		//Check whether we converged. Note that accuracy is only assigned when the layer exceeds the prescribed required accuracy. This is because we want to have more than one convergence criterion (both h and theta based), we say accuracy=0 is a sign of convergence in the whole domain.
				boolConvergence=false;
			}


			if(WriteOutNumerics_Level1==true) printf("CONVERGENCE:  layer: %d/%d --- acc_h: %.10f acc_theta: %.10f acc: %.10f def_norm: %f converged? %s\n", track_trigger_layer_accuracy, trigger_layer_accuracy, track_accuracy_h, track_accuracy_theta, accuracy, deficit_vector_norm, (boolConvergence)?"yes":"no");

			//Copy solution, to prepare for next iteration
			for (i = uppernode; i >= lowernode; i--) {
				h_np1_m[i]=h_np1_mp1[i];
				theta_np1_m[i]=theta_np1_mp1[i];
				theta_i_np1_m[i]=theta_i_np1_mp1[i];
			}

			//Rewind time step control: retry finding solution with smaller time step when MAX_ITER is exceeded. Also when max_delta_h is too large we do a rewind, else the model is starting to blow up.
			if((niter>MAX_ITER || max_delta_h>MAX_ALLOWED_DELTA_H) && dt > MIN_VAL_TIMESTEP && boolConvergence==false) {

				if(WriteOutNumerics_Level1==true) std::cout << "REWIND: timestep " << std::setprecision(20) << dt << std::setprecision(6) << " ---> ";

				niter_seqrewinds++;				//We increase the sequential rewinds counter. For this first rewind, this will give a power of 1 for the new time step, etc.
										//When we find a good solution again, we reset this counter to 0.

				TimeAdvance-=dt;				//We do a rewind of the time step, so adjust TimeAdvance with the time step.
				dt*=pow(0.333, double(niter_seqrewinds));	//Now make the time step smaller, we use niter_seqrewinds to speed up when we encounter multiple consecutive rewinds.
										//The value of 0.333 is taken from the HYDRUS-manual, where they do this in case a rewind is necessary.

				if(WriteOutNumerics_Level1==true) std::cout << std::setprecision(20) << dt << " (trigger layer: " << trigger_layer_blowup << "  accuracy: " << std::setprecision(10) << accuracy << " max_delta_h: " << max_delta_h << ")\n" << std::setprecision(6);

				niter=0;					//Because of the rewind, we start again with the iterations.
				niter_nrewinds++;				//Increase rewind counter
				stats_nrewinds++;				//Increase the statistics rewind counter
				boolConvergence=false;				//Of course, when we need to rewind, we have had no convergence.
				DoRewindFlag=true;				//Set DoRewindFlag to true, to quit the iteration loop.
				StopLoop=false;					//In case StopLoop was set true (last time step), we set it back to false. It might be that the smaller time step won't match the SNOWPACK time step any longer.
				for (i = uppernode; i >= lowernode; i--) {	//We have to reset the whole domain, because we do the time step for the whole domain.
					h_np1_m[i]=h_n[i];
					theta_np1_m[i]=theta_n[i];
					theta_i_np1_m[i]=theta_i_n[i];		//Set back ice content due to soil freezing/thawing
					delta_Te_i[i]=0.;			//Reset temperature change due to soil freezing/thawing
					delta_Te_adv_i[i]=0.;			//Reset temperature change due to heat advection by water flowing
				}
			}

			if(niter>500) {
				//Print latest state for debugging:
				if(WriteOutNumerics_Level0==true) {
					for (i = uppernode; i >= lowernode; i--) {
						printf("ITER: %d i: %d active? %s; h_n: %.15f (h_np1: %.15f) theta: %.15f (%f-%f) ice: %f/%f (vg_params: %f %f %f)\n", niter, i, (activelayer[i])?"true":"false", h_n[i], h_np1_m[i], theta_n[i], theta_r[i], theta_s[i], EMS[SnowpackElement[i]].theta[ICE], theta_i_n[i], alpha[i], m[i], n[i]);
					}
					printf("BOUNDARYTOPFLUX: [ BC: %d ] %.15f %.15f\n", TopBC, TopFluxRate, surfacefluxrate);
				}
				if(SafeMode==false) {
					if(niter>500) {
						prn_msg(__FILE__, __LINE__, "err", Date(), "Richards-Equation solver did not converge: reached maximum number of iterations (500), with a time step: %G\n", dt);
					}
					std::cout << "  POSSIBLE SOLUTIONS:\n  =============================================================================\n";
					if(snowpack_dt>900) std::cout << "    - SNOWPACK time step is larger than 15 minutes. This numerical problem\n      may be resolved by using a time step of 15 minutes.\n";
#ifndef CLAPACK
					std::cout << "    - SNOWPACK was not compiled with BLAS and CLAPACK libraries.\n      Try installing libraries BLAS and CLAPACK and use solver TGSV (default).\n";
#endif
					std::cout << "    - Verify that the soil is not initialized in a very dry or a very\n      wet state.\n";
					if(BottomBC!=FREEDRAINAGE) std::cout << "    - If the soil is saturated, try setting LB_COND_WATERFLUX = FREEDRAINAGE\n      in the [SnowpackAdvanced] section of the ini-file.\n";
					if(BottomBC!=WATERTABLE) std::cout << "    - If the soil is dry, try setting LB_COND_WATERFLUX = WATERTABLE in the\n      [SnowpackAdvanced] section of the ini-file.\n";
					std::cout << "    - Try bucket scheme, by setting WATERTRANSPORTMODEL_SNOW = BUCKET and\n      WATERTRANSPORTMODEL_SOIL = BUCKET in the [SnowpackAdvanced] section\n      of the ini-file.\n";
					std::cout << "    - When using Canopy module, there is a known issue with transpiration.\n      Please see issue 471 (http://models.slf.ch/p/snowpack/issues/471/).\n";
					throw;		//We are lost. We cannot do another rewind and decrease time step (if we could, niter is reset).
				} else {
					if(seq_safemode>3) {
						std::cout << "ERROR in Richards-Equation solver: SafeMode was not able to rescue simulation!\n";
						throw;		//We are lost. We cannot do another rewind and decrease time step (if we could, niter is reset).
					}
					std::cout << "WARNING in Richards-Equation solver: SafeMode was needed to rescue simulation! ";
					//We rescue the simulation ...
					double SafeMode_MBE=0.;		// Mass balance error (kg/m^2) due to SafeMode

					//Do a rewind
					TimeAdvance-=dt;
					dt=1E-3;

					niter=0;
					niter_nrewinds++;				//Increase rewind counter
					stats_nrewinds++;				//Increase the statistics rewind counter
					seq_safemode++;					//Increase counter for sequential safemode
					boolConvergence=false;				//Of course, when we need to rewind, we have had no convergence.
					DoRewindFlag=true;				//Set DoRewindFlag to true, to quit the iteration loop.
					StopLoop=false;					//In case StopLoop was set true (last time step), we set it back to false. It might be that the smaller time step won't match the SNOWPACK time step any longer.
					mass1=0.;					//Because we fiddle around with theta, we should update mass1 (mass at beginning of time step)
					for (i = uppernode; i >= lowernode; i--) {	//We have to reset the whole domain, because we do the time step for the whole domain.
						// Update the SafeMode mass balance error tracking variable by "removing" all water
						SafeMode_MBE-=(theta_n[i]+theta_i_n[i])*dz[i]*Constants::density_water;
						// Make sure pressure head is in secure limits:
						h_n[i]=MAX(h_d, MIN(h_e[i], h_n[i]));
						theta_n[i]=fromHtoTHETAforICE(h_n[i], theta_r[i], theta_s[i], alpha[i], m[i], n[i], Sc[i], h_e[i], theta_i_n[i]);
						//Deal with dry layers
						if(theta_n[i]+theta_i_n[i] < theta_r[i]+(REQUIRED_ACCURACY_THETA/1000.)) {
							theta_n[i]=theta_r[i]+(REQUIRED_ACCURACY_THETA/1000.);
							h_n[i]=fromTHETAtoHforICE(theta_n[i], theta_r[i], theta_s[i], alpha[i], m[i], n[i], Sc[i], h_e[i], h_d, theta_i_n[i]);
						}
						//Deal with wet layers
						if(theta_n[i]+theta_i_n[i] > theta_s[i]-(REQUIRED_ACCURACY_THETA/1000.)) {
							theta_i_n[i]*=0.90;
							theta_n[i]=((theta_n[i]-theta_r[i])*0.9)+theta_r[i];
							h_n[i]=fromTHETAtoHforICE(theta_n[i], theta_r[i], theta_s[i], alpha[i], m[i], n[i], Sc[i], h_e[i], h_d, theta_i_n[i]);
						}
						// Update the SafeMode mass balance error tracking variable by "adding" the water again
						SafeMode_MBE+=(theta_n[i]+theta_i_n[i])*dz[i]*Constants::density_water;

						h_np1_m[i]=h_n[i];			//Reset initial guess for next iteration
						theta_np1_m[i]=theta_n[i];		//Reset initial guess for next iteration
						theta_i_np1_m[i]=theta_i_n[i];		//Set back ice content due to soil freezing/thawing

						delta_Te_i[i]=0.;			//Reset temperature change due to soil freezing/thawing
						delta_Te_adv_i[i]=0.;			//Reset temperature change due to heat advection by water flowing

						//The real rescue is to throw away the sink/source terms:
						SafeMode_MBE+=s[i]*(sn_dt-TimeAdvance)*Constants::density_water*dz[i];
						s[i]=0.;
						//And to halve the TopFluxRate:
						SafeMode_MBE+=(surfacefluxrate/2.)*(sn_dt-TimeAdvance)*Constants::density_water;
						surfacefluxrate/=2.;
						//Now update mass1, which may have changed due to SafeMode throwing some water away, or introducing some water:
						mass1+=(theta_n[i]+(theta_i_n[i]*(Constants::density_ice/Constants::density_water)))*dz[i];
					}
					std::cout << "Estimated mass balance error due to SafeMode: " << std::scientific << SafeMode_MBE << std::fixed << " kg/m^2\n";
				}
			}
		}	//End of iteration loop


		//If the solver wants to do a rewind, it exits the previous loop. However, in this case the time step is already adjusted and the matrices are put back to initial state. So we should skip the preparation for the next time step.
		if(DoRewindFlag==false) {
			niter_seqrewinds=0;		//We found a good solution (again), so we reset this niter_seqrewind counter.
			seq_safemode=0;
			stats_niters+=niter;		//Update the statistics with the number of iterations necessary to find the solution, ignoring possible iterations done before a rewind.
			stats_nsteps++;			//Update the statistics of the number of steps.

			//Prepare for next time step:
			for (i = uppernode; i >= lowernode; i--) {				//Cycle through all Richards solver domain layers.
				//Apply change in temperature due to soil freezing or thawing and heat advection by flowing water:
				if(SnowpackElement[i]<int(Xdata.SoilNode)) {			//HACK, TODO: remove type inconstency in comparison
					//Freezing and thawing
					if(fabs(delta_Te_i[i]) > 0.) {				//Check if phase change did occur in soil
						delta_Te[i]+=delta_Te_i[i];
						EMS[SnowpackElement[i]].QIntmf+=(Constants::density_ice*(theta_i_np1_mp1[i]-theta_i_n[i])*(Constants::specific_heat_water-Constants::specific_heat_ice)*(T_melt[i]-Constants::melting_tk))/dt;
						EMS[SnowpackElement[i]].melting_tk=EMS[SnowpackElement[i]].freezing_tk=T_melt[i];
						// Now that we have performed a phase change, we should correct the nodal temperatures too. This will be done later in PhaseChange,
						// by using Qmf to determine amount of phase change that occurred.
					}

					delta_Te_adv[i]+=delta_Te_adv_i[i];
				}

				//We adapted the elements and nodes to the temperature change, so set it to 0.
				delta_Te_i[i]=0.;
				delta_Te_adv_i[i]=0.;

				//Set initial solution for next iteration
				if(activelayer[i]==true) {
					theta_np1_mp1[i]=fromHtoTHETAforICE(h_np1_m[i], theta_r[i], theta_s[i], alpha[i], m[i], n[i], Sc[i], h_e[i], theta_i_np1_mp1[i]);
				} else {
					theta_np1_mp1[i]=theta_n[i];
				}

				delta_h_dt[i]=(h_np1_mp1[i]-h_n[i])/dt;				//We make delta_h_dt relative to time step. If time step is allowed to change, we can use this delta_h_dt (actually a derivative) to better estimate solution in next time step.
				delta_theta_dt[i]=(theta_np1_mp1[i]-theta_n[i])/dt;		//We make delta_theta_dt relative to time step. If time step is allowed to change, we can use this delta_h_dt (actually a derivative) to better estimate solution in next time step.
				delta_theta_i_dt[i]=(theta_i_np1_mp1[i]-theta_i_n[i])/dt;	//We make delta_theta_i_dt relative to time step. If time step is allowed to change, we can use this delta_h_dt (actually a derivative) to better estimate solution in next time step.

				//Copy current state:
				h_n[i]=h_np1_mp1[i];
				theta_n[i]=theta_np1_mp1[i];
				theta_i_n[i]=theta_i_np1_mp1[i];
			}



			//Determine (estimate) flux across boundaries (downward ==> positive flux):
			//This is an additional check for the boundaries.
			actualtopfluxcheck+=((delta_theta_dt[uppernode]*dt)*dz[uppernode])+(((h_n[uppernode]-h_n[uppernode-1])/dz_down[uppernode])+cos_sl)*k_np1_m_im12[uppernode]*dt;
			actualtopflux+=TopFluxRate*dt;
			refusedtopflux+=(surfacefluxrate-TopFluxRate)*dt;
			if(aBottomBC==DIRICHLET) {
				actualbottomflux+=-1.*((delta_theta_dt[lowernode]+(delta_theta_i_dt[lowernode]*(Constants::density_ice/Constants::density_water))*dt)-((((h_np1_mp1[lowernode+1]-h_np1_mp1[lowernode])/dz_up[lowernode])+cos_sl)*k_np1_m_ip12[lowernode]*dt));
			} else {
				actualbottomflux+=BottomFluxRate*dt;
			}

			massbalanceerror_sum+=massbalanceerror;
			if(WriteOutNumerics_Level1==true) printf("MASSBALANCE: mass1 %.8f    mass2 %.8f    delta %.8f\n", mass1, mass2, massbalanceerror);


			//Determine flux at soil snow interface (note: postive=flux upward, negative=flux downward):
			if (int(nsoillayers_snowpack)<int(nE)) {	//We have snow layers
				if(toplayer==nsoillayers_snowpack) {	//We run RE-solver only for soil layers AND have snow layers in the model (meaning TopFluxRate is coming from snow)
					//These lines are not active, as this particular case (RE for soil, other for snow), is dealt with in compTransportMass.
					//snowsoilinterfaceflux1+=surfacefluxrate*dt;
					//snowsoilinterfaceflux2+=surfacefluxrate*dt;
				} else {
					// See McCord (1996). Note: I think there is a minus sign missing there.
					// In any case: snowsoilinterfaceflux > 0 means influx!
					snowsoilinterfaceflux_after=((((h_n[nsoillayers_richardssolver]-h_n[nsoillayers_richardssolver-1])/dz_up[nsoillayers_richardssolver-1])+cos_sl)*k_np1_m_ip12[nsoillayers_richardssolver-1]*dt);
					snowsoilinterfaceflux1+=snowsoilinterfaceflux_after;
					// Other method to estimate soil snow interface flux (based on average before and end of time step).
					snowsoilinterfaceflux2+=0.5*(snowsoilinterfaceflux_before+snowsoilinterfaceflux_after);
				}
			} else {
				//Make the commented lines active if you whish to add the TopFluxRate to the snowsoilinterfaceflux even when no snow is present.
				//snowsoilinterfaceflux1+=TopFluxRate*dt;
				//snowsoilinterfaceflux2+=TopFluxRate*dt;
			}

			if(WriteOutNumerics_Level2==true) printf("CONTROL: %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %f\n", surfacefluxrate, TopFluxRate, actualtopflux, actualtopfluxcheck, BottomFluxRate, actualbottomflux, snowsoilinterfaceflux1, snowsoilinterfaceflux2, dt);

			//Time step control
			//This time step control increases the time step when niter is below a certain value. When rewinds occurred in the time step, no change is done (dt already adapted by the rewind-mechanim), if too many iterations, time step is decreased.
			if((niter+(MAX_ITER*niter_nrewinds))<INCR_ITER) {
				dt*=1.25;
			} else {
				if(niter_nrewinds==0 && niter>DECR_ITER) {
					dt*=0.5;
				}
			}

			//Limit time steps:
			if ( dt < MIN_VAL_TIMESTEP) dt=MIN_VAL_TIMESTEP;
			if ( dt > MAX_VAL_TIMESTEP) dt=MAX_VAL_TIMESTEP;

			//Special limit in case of snow:
			if(int(nsoillayers_snowpack)<int(nE) && dt>MAX_VAL_TIMESTEP_FOR_SNOW) {
				dt=MAX_VAL_TIMESTEP_FOR_SNOW;
			}

			//Time step statistics
			if(stats_min_dt>dt) stats_min_dt=dt;
			if(stats_max_dt<dt) stats_max_dt=dt;

			//Update mass balance status variable (mass1 becomes mass2 to serve as reference for the next iteration)
			mass1=mass2;
			//And snowsoilinterfaceflux_before becomes snowsoilinterfaceflux_after for the next time step.
			snowsoilinterfaceflux_before=snowsoilinterfaceflux_after;

			if(WriteOutNumerics_Level1==true) printf("NSTEPS: %d, TIME ADVANCE: %f, ITERS NEEDED: %d [%d], ACTUALTOPFLUX: %.10f     ---> new dt: %.15f (Problematic layer: %d of %d)\n", nsteps, TimeAdvance, niter, niter_nrewinds, actualtopflux, dt, track_trigger_layer_accuracy, uppernode);
		}	//END DoRewindFlag==false
	}
	while(StopLoop==false);							//This is the main loop to perform 1 SNOWPACK time step

	//Because the Richards solver domain and the snowpack domain does not necessarily match (Richards solver domain can have more layers), we do the following trick:
	//1) We first empty all the layers in the snowpack domain
	for (i = toplayer-1; i >= 0; i--) {							//We loop over all SNOWPACK layers ...
		EMS[i].theta[WATER]=0.;								//... and set water content to 0
		EMS[i].theta_r=0.;								// and set residual water content to 0
		if(EMS[i].theta[SOIL]>Constants::eps2) {					//We are in soil
			EMS[i].theta[ICE]=0.;							//... set ice content to 0.
		}
	}
	//2) Now we fill them, scaling with the Richards solver domain layer heights. The "dictionary" makes that all the Richards solver domain layers belonging to a single SNOWPACK layer are summed together.
	for (i = uppernode; i >= lowernode; i--) {						//We loop over all Richards solver domain layers
		if(EMS[SnowpackElement[i]].theta[SOIL]>Constants::eps2) {			//We are in soil
			if(activelayer[i]==true) {						//If we have more water than theta_r, we copy the end state from the solver
				//EMS[SnowpackElement[i]].theta[WATER]+=dz[i]*fromHtoTHETA(h_n[i], theta_r[i], theta_s[i], alpha[i], m[i], n[i], Sc[i], h_e[i], theta_d[i]);
				EMS[SnowpackElement[i]].theta[WATER]+=dz[i]*fromHtoTHETAforICE(h_n[i], theta_r[i], theta_s[i], alpha[i], m[i], n[i], Sc[i], h_e[i], theta_i_n[i]);
				EMS[i].theta[ICE]+=dz[i]*theta_i_n[i];
				/////EMS[i].theta[ICE]+=dz[i]*theta_i_n[i]*(Constants::density_water/Constants::density_ice);
			} else {								//We are in "dry" conditions, and we just copy the initial state.
				EMS[SnowpackElement[i]].theta[WATER]+=dz[i]*theta_n[i];
			}
		} else {									//We are in snow
			if(activelayer[i]==true) {
				//EMS[SnowpackElement[i]].theta[WATER]+=dz[i]*fromHtoTHETA(h_n[i], theta_r[i], theta_s[i], alpha[i], m[i], n[i], Sc[i], h_e[i], theta_d[i]);
				EMS[SnowpackElement[i]].theta[WATER]+=dz[i]*fromHtoTHETAforICE(h_n[i], theta_r[i], theta_s[i], alpha[i], m[i], n[i], Sc[i], h_e[i], theta_i_n[i]);
			} else {
				EMS[SnowpackElement[i]].theta[WATER]+=dz[i]*theta_n[i];
			}
		}
		EMS[SnowpackElement[i]].theta_r+=dz[i]*theta_r[i];
		//EMS[SnowpackElement[i]].head+=dz[i]*h_n[i];
	}
	//3) Now scale them with the SNOWPACK domain layer heights and adjust properties accordingly
	const double MIN_VAL_THETA_SNOWPACK=0.0;						//Minimum water content that is allowed to be passed on to the rest of SNOWPACK. When theta is below this value, it is truncated to 0, before leaving the Richards solver.
	for (i = toplayer-1; i >= 0; i--) {							//We loop over all SNOWPACK layers
		//Do the actual scaling
		if( (i>nsoillayers_snowpack-1 && EMS[i].theta[WATER]/EMS[i].L>MIN_VAL_THETA_SNOWPACK) || (i<=nsoillayers_snowpack-1)) {	//If we are in snow and have enough water to pass the water on to the rest of SNOWPACK, OR if we are in soil:
			EMS[i].theta[WATER]/=EMS[i].L;							//Scale each layer with SNOWPACK layer height
			if(EMS[i].theta[SOIL]>Constants::eps2) {					//We are in soil ...
				EMS[i].theta[ICE]/=EMS[i].L;						//... scale ice content.
			}
			EMS[i].theta_r/=EMS[i].L;							//Scale each layer with SNOWPACK layer height
			//EMS[i].head/=EMS[i].L;								//Scale each layer with SNOWPACK layer height

			//In case we had to melt ice to get theta_r, we have to adjust the temperature:
			EMS[i].Te -= dT[i];
			if(i==int(nE)-1 && i>=0) {							//HACK, TODO: remove type inconstency in comparison
				NDS[i+1].T-=dT[i];
				NDS[i].T-=dT[i];
			}

			//And adjust all the properties accordingly
			EMS[i].theta[AIR]=1.-EMS[i].theta[WATER]-EMS[i].theta[ICE]-EMS[i].theta[SOIL];
			//Now we have checked everything, we make it fit between [0, 1]: to get rid off all round-off errors
			EMS[i].theta[AIR]=MAX(0, MIN(1., EMS[i].theta[AIR]));
			EMS[i].theta[WATER]=MAX(0, MIN(1., EMS[i].theta[WATER]));
			EMS[i].theta[ICE]=MAX(0, MIN(1., EMS[i].theta[ICE]));
			EMS[i].Rho = (EMS[i].theta[ICE] * Constants::density_ice) + (EMS[i].theta[WATER] * Constants::density_water) + (EMS[i].theta[SOIL] * EMS[i].soil[SOIL_RHO]);
			EMS[i].M=EMS[i].L*EMS[i].Rho;
			EMS[i].heatCapacity();

			//Every change in ice content in a specific layer must be associated with phase changes. Store the associated energy accordingly.
			EMS[i].Qmf += ((EMS[i].theta[ICE]-snowpackBACKUPTHETAICE[i]) * Constants::density_ice * Constants::lh_fusion) / snowpack_dt;	// Units: [W m-3]
			//We transferred the temperature change of the element due to soil freezing/thawing in Qmf, so reset delta_Te:
			delta_Te[i]=0.;
		} else {										//We are in snow and don't have enough water, snow should be dry, so set back to initial values.
			//NOTE: there is an issue to be solved here when Richard domain does not match snowpack domain (use of sublayers)!!
			wateroverflow[i]+=(EMS[i].theta[WATER]-theta_d[i]);				//This is water which stays or is taken out from the domain by this layer.

			EMS[i].theta[ICE]=snowpackBACKUPTHETAICE[i];
			EMS[i].theta[WATER]=snowpackBACKUPTHETAWATER[i];
			//Now we have checked everything, we make it fit between [0, 1]: to get rid off all round-off errors
			EMS[i].theta[AIR]=MAX(0, MIN(1., EMS[i].theta[AIR]));
			EMS[i].theta[WATER]=MAX(0, MIN(1., EMS[i].theta[WATER]));
			EMS[i].theta[ICE]=MAX(0, MIN(1., EMS[i].theta[ICE]));
		}

		//Then check the volumetric contents. This we do, to make a crash at this place, and we have information about the Richards solver available in the core file.
		//Do some checks on volumetric contents:
		const double sum=EMS[i].theta[AIR] + EMS[i].theta[WATER] + EMS[i].theta[ICE] + EMS[i].theta[SOIL];
		if(EMS[i].theta[WATER]<0.-Constants::eps2 || EMS[i].theta[AIR]<0.-Constants::eps2 || EMS[i].theta[AIR] > 1.+Constants::eps2 || EMS[i].theta[ICE]<0.-Constants::eps2 || EMS[i].theta[ICE] > 1.+Constants::eps2) {
			printf("ERROR at layer %d: sum=%f air=%f ice=%f soil=%f water=%f\n", i, sum, EMS[i].theta[AIR], EMS[i].theta[ICE], EMS[i].theta[SOIL], EMS[i].theta[WATER]);
			printf("   -- if this happens and ice<0, check theta_d. Maybe there was so much water created, that it was more than there was ice. This is not accounted for.\n");
			throw;
		}
		if(sum > 1.+Constants::eps2) {
			printf("ERROR at layer %d: sum=%f air=%f ice=%f soil=%f water=%f\n", i, sum, EMS[i].theta[AIR], EMS[i].theta[ICE], EMS[i].theta[SOIL], EMS[i].theta[WATER]);
			throw;
		}
		if(sum < 1.-Constants::eps2) {
			printf("ERROR at layer %d: sum=%f air=%f ice=%f soil=%f water=%f\n", i, sum, EMS[i].theta[AIR], EMS[i].theta[ICE], EMS[i].theta[SOIL], EMS[i].theta[WATER]);
			throw;
		}
	}

	for (i = toplayer-1; i >= 0; i--) {							//We loop over all SNOWPACK layers ...
		//Heat advection by water flow
		double deltaN=0.;
		if(i == int(nE)-1) {								//HACK, TODO: remove type inconstency in comparison
			deltaN=(delta_Te_adv[i] * (EMS[i].c[TEMPERATURE]*EMS[i].Rho*EMS[i].L)) / (EMS[i].c[TEMPERATURE]*EMS[i].Rho*EMS[i].L + 0.5*EMS[i-1].c[TEMPERATURE]*EMS[i-1].Rho*EMS[i-1].L);
		} else {
			if(i==0) {
				deltaN=(delta_Te_adv[i] * (EMS[i].c[TEMPERATURE]*EMS[i].Rho*EMS[i].L)) / (0.5*EMS[i+1].c[TEMPERATURE]*EMS[i+1].Rho*EMS[i+1].L + EMS[i].c[TEMPERATURE]*EMS[i].Rho*EMS[i].L);
			} else {
				deltaN=(delta_Te_adv[i] * (EMS[i].c[TEMPERATURE]*EMS[i].Rho*EMS[i].L)) / (0.5*EMS[i+1].c[TEMPERATURE]*EMS[i+1].Rho*EMS[i+1].L + EMS[i].c[TEMPERATURE]*EMS[i].Rho*EMS[i].L + 0.5*EMS[i-1].c[TEMPERATURE]*EMS[i-1].Rho*EMS[i-1].L);
			}
		}
		NDS[i+1].T+=deltaN;
		NDS[i].T+=deltaN;
		if(fabs(deltaN)>0.) {
			if(i < int(nE)-1) EMS[i+1].Te=0.5*(NDS[i+2].T+NDS[i+1].T);		//HACK, TODO: remove type inconstency in comparison
			EMS[i].Te=0.5*(NDS[i+1].T+NDS[i].T);
			if(i > 0) EMS[i-1].Te=0.5*(NDS[i].T+NDS[i-1].T);
		}
		if (WriteOutNumerics_Level2==true) printf("SENDING at layer %d: sum=%f air=%.15f ice=%.15f soil=%.15f water=%.15f Te=%.15f\n", i, EMS[i].theta[AIR]+EMS[i].theta[ICE]+EMS[i].theta[SOIL]+EMS[i].theta[WATER], EMS[i].theta[AIR], EMS[i].theta[ICE], EMS[i].theta[SOIL], EMS[i].theta[WATER], EMS[i].Te);
	}

	double totalwateroverflow=0.;					//Total water outflow due to numerical issues (requiring minimum theta_r, maximum theta_s, etc), in m^3/m^2
	for (i = uppernode; i>=lowernode; i--) {
		totalwateroverflow+=wateroverflow[i]*dz[i];
		if(i==nsoillayers_richardssolver) {
			// I decided to put all wateroverflow from snow directly in the snowsoilinterfaceflux, although the wateroverflow may occur somewhere in the snowpack.
			snowsoilinterfaceflux1+=totalwateroverflow;
		}
	}


	if(WriteOutNumerics_Level1==true) {
		printf("ACTUALTOPFLUX: [ BC: %d ] %.15f %.15f %.15f CHK: %.15f %f\n", TopBC, actualtopflux/snowpack_dt, refusedtopflux/snowpack_dt, surfacefluxrate, actualtopfluxcheck/snowpack_dt, (surfacefluxrate!=0.)?(actualtopflux/snowpack_dt)/surfacefluxrate:0.);
		printf("ACTUALBOTTOMFLUX: [ BC: %d ] %.15f %.15f %.15f %f    K_ip1=%.15f\n", BottomBC, actualbottomflux, actualbottomflux/snowpack_dt, BottomFluxRate, (BottomFluxRate!=0.)?(actualbottomflux/snowpack_dt)/BottomFluxRate:0., k_np1_m_ip12[lowernode]);
		// This is more or less for testing only. This snowsoilinterfaceflux should anyway be stored in MS_SNOWPACK_RUNOFF and found in the met file
		printf("SNOWSOILINTERFACEFLUX: %.15f %.15f\n", snowsoilinterfaceflux1/snowpack_dt, snowsoilinterfaceflux2/snowpack_dt);
	}


	if(WriteOutNumerics_Level0==true) printf("WATERBALANCE: %.15f %.15f %.15f CHK1: %.15f CHK2: %.15f  WATEROVERFLOW: %.15f MB_ERROR: %.15f\n", actualtopflux/snowpack_dt, refusedtopflux/snowpack_dt, surfacefluxrate, (surfacefluxrate!=0.)?(actualtopflux/snowpack_dt)/surfacefluxrate:0., actualtopfluxcheck/snowpack_dt, totalwateroverflow, massbalanceerror_sum);

	//Update soil runoff (mass[MS_SOIL_RUNOFF] = kg/m^2). Note: it does not matter whether SNOWPACK is run with soil or not. MS_SOIL_RUNOFF is always runoff from lower boundary.
	Sdata.mass[SurfaceFluxes::MS_SOIL_RUNOFF] += actualbottomflux*Constants::density_water;

	// Update snow pack runoff (mass[MS_SNOWPACK_RUNOFF] = kg/m^2 (almost equal to mm/m^2), surfacefluxrate=m^3/m^2/s and snowsoilinterfaceflux = m^3/m^2):
	// NOTE: snowsoilinterfaceflux will only be non-zero IF there is a snowpack AND we solve the richards equation also for snow! Else, snowpack runoff is calculated in the original WaterTransport functions.
	Sdata.mass[SurfaceFluxes::MS_SNOWPACK_RUNOFF] += snowsoilinterfaceflux1*Constants::density_water;

	//Deal with the situation that evaporation flux was limited in case of snow. Then, sublimate ice matrix.
	if (refusedtopflux<0. && toplayer>nsoillayers_snowpack) {
		//Be careful: refusedtopflux = m^3/m^2 and not m^3/m^2/s!!!
		//Now invert the calculation of ql, using refusedtopflux. This amount of ql should be used for sublimation.
		double ql=(refusedtopflux/sn_dt)*Constants::density_water*Constants::lh_vaporization;

		double dL=0.;
		std::vector<double> M_Solutes(Xdata.number_of_solutes, 0.); // Mass of solutes from disappearing phases
		size_t e = nE-1;
		while ((e >= Xdata.SoilNode) && (ql < -Constants::eps2)) {  // While energy is available and we are in snow
			const double L0 = EMS[e].L;
			// If there is no water or if there was not enough water ...
			// Note: as we do not pass through mergeElements anymore, we must assure that elements do not disappear here.
			// By specifying a minimum value just below the Snowpack::min_ice_content, we make sure the element gets removed the next time it passes mergeElements.
			const double theta_i0 = MAX(0., EMS[e].theta[ICE] - (0.99*Snowpack::min_ice_content));
			double M = theta_i0*Constants::density_ice*L0;
			double dM = ql*sn_dt/Constants::lh_sublimation;
			if (-dM > M) dM = -M;

			dL = dM/(EMS[e].Rho);
			if (e < Xdata.SoilNode) {
				dL = 0.;
			}
			NDS[e+1].z += dL; EMS[e].L0 = EMS[e].L = L0 + dL;
			NDS[e+1].z += NDS[e+1].u; NDS[e+1].u = 0.0;

			EMS[e].E = EMS[e].dE = EMS[e].Ee = EMS[e].Ev = EMS[e].S = 0.0;
			EMS[e].theta[ICE] *= L0/EMS[e].L;
			EMS[e].theta[ICE] += dM/(Constants::density_ice*EMS[e].L);
			EMS[e].theta[WATER] *= L0/EMS[e].L;
			for (size_t ii = 0; ii < Xdata.number_of_solutes; ii++) {
				EMS[e].conc[ICE][ii] *= L0*theta_i0/(EMS[e].theta[ICE]*EMS[e].L);
			}

			EMS[e].M += dM;
			// Instead of evaporating, we sublimate the ice matrix:
			Sdata.mass[SurfaceFluxes::MS_EVAPORATION] -= dM*(Constants::lh_sublimation/Constants::lh_vaporization);	//Correct evaporation for sublimated mass
			Sdata.mass[SurfaceFluxes::MS_SUBLIMATION] += dM;							//Add mass to sublimation
			ql -= dM*Constants::lh_sublimation/sn_dt;     // Update the energy used

			//Update volumetric contents
			EMS[e].theta[AIR]=1.-EMS[e].theta[ICE]-EMS[e].theta[WATER]-EMS[e].theta[SOIL];
			EMS[e].Rho = (EMS[e].theta[ICE] * Constants::density_ice) + (EMS[e].theta[WATER] * Constants::density_water) + (EMS[e].theta[SOIL] * EMS[e].soil[SOIL_RHO]);
			EMS[e].heatCapacity();

			e--;
		}
		//Remaining energy should go back again into refusedtopflux and also should not be counted as evaporation
		Sdata.mass[SurfaceFluxes::MS_EVAPORATION]-=ql*sn_dt/Constants::lh_vaporization;
		refusedtopflux=MIN(0., (ql*sn_dt)/(Constants::density_water*Constants::lh_vaporization));
	}
	if(refusedtopflux<0. && toplayer==nsoillayers_snowpack) {
		//Be careful: refusedtopflux = m^3/m^2 and not m^3/m^2/s!!!
		//Now invert the calculation of ql, using refusedtopflux. This amount of ql should be used for sublimation.
		double ql=(refusedtopflux/sn_dt)*Constants::density_water*Constants::lh_vaporization;
		refusedtopflux=0.;
		//Remaining energy should not be counted as evaporation
		Sdata.mass[SurfaceFluxes::MS_EVAPORATION]-=ql*sn_dt/Constants::lh_vaporization;
		//The energy is substracted from the top element
		const double tmp_delta_Te = ql / (EMS[nsoillayers_snowpack-1].c[TEMPERATURE] * EMS[nsoillayers_snowpack-1].Rho);
		NDS[nsoillayers_snowpack].T += 2.*tmp_delta_Te;
		EMS[nsoillayers_snowpack-1].Te += tmp_delta_Te;
	}

	//If we could not handle all incoming water at top boundary AND we have snow AND we solve RE for snow:
	if(refusedtopflux>0. && toplayer>int(Xdata.SoilNode)) {
		Sdata.mass[SurfaceFluxes::MS_SNOWPACK_RUNOFF] += refusedtopflux*Constants::density_water;
	}
	//If we could not handle all snowpack runoff when not modelling snow with RE:
	if(refusedtopflux>0. && toplayer==int(Xdata.SoilNode) && Xdata.getNumberOfElements()>Xdata.SoilNode ){
		Sdata.mass[SurfaceFluxes::MS_SNOWPACK_RUNOFF] += refusedtopflux*Constants::density_water;
	}
	//We want wateroverflow in the snow to be a source/sink term. Therefore, these lines are inactive.
	//if(totalwateroverflow>0. && toplayer>Xdata.SoilNode) {
	//	Sdata.mass[SurfaceFluxes::MS_SNOWPACK_RUNOFF] += totalwateroverflow*Constants::density_water;
	//}

	surfacefluxrate=0.;			//As we now have used the rate for the current time step, reset the value.


	//Now calculate freezing point depression:
	for (i = 0; i<=toplayer-1; i++) {
		if(EMS[i].theta[SOIL]<Constants::eps2) {
			EMS[i].melting_tk=T_0;
		} else {
			if(AllowSoilFreezing == true) {
				//For soil layers solved with Richards Equation, everything (water transport and phase change) is done in this routine, except calculating the heat equation.
				//To suppress phase changes in PhaseChange.cc, set the melting and freezing temperature equal to the element temperature:
				EMS[i].freezing_tk=EMS[i].melting_tk=EMS[i].Te;
			} else {
				EMS[i].freezing_tk=EMS[i].melting_tk=T_melt[i];
				//This is a trick. Now that we deal with phase change in soil right here, we set the melting and freezing temperatures equal to the current Element temperature, so that
				//in CompTemperatures, the element temperature will not be adjusted to freezing temperature just because there is water in it!
				EMS[i].Te=MAX(EMS[i].Te, T_0);	//Because we don't allow soil freezing, soil remains 0 degC.
				NDS[i].T=MAX(NDS[i].T, T_0);
				NDS[i+1].T=MAX(NDS[i+1].T, T_0);
				EMS[i].melting_tk=EMS[i].Te;
				EMS[i].freezing_tk=EMS[i].Te;
			}
		}
		if(WriteOutNumerics_Level2==true) 
			std::cout << "EMS[" << i << "].melting_tk = " << EMS[i].melting_tk << ", EMS[" << i << "].freezing_tk = " << EMS[i].freezing_tk << " (ice: " << EMS[i].theta[ICE] << ")\n";
	}

	//print solver statistics
	if(WriteOutNumerics_Level0==true) {
		std::cout << "SOLVERSTATISTICS: max_dt= " << std::setprecision(5) << stats_max_dt << "   min_dt= " << std::setprecision(20) << stats_min_dt << std::setprecision(6) << "   nsteps_total= " << stats_nsteps << "   niter_total= " << stats_niters << "   nrewinds_total= " << stats_nrewinds << "  Last active solver: ";
		switch (ActiveSolver) {
			case DGESVD:
				std::cout << "DGESVD/DGESDD.";
				break;
			case DGTSV:
				std::cout << "DGTSV.";
				break;
			case TDMA:
				std::cout << "TDMA.";
				break;
		}
		std::cout << " BS_avg_iter: " << double(double(bs_stats_totiter)/double(stats_niters*nsoillayers_richardssolver)) << " BS_max_iter: " << bs_stats_maxiter << "\n";
	}
}

#ifdef __clang__
#pragma clang diagnostic pop
#endif
