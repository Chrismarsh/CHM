//---------------------------------------------------------------------------

#ifndef g_soilH
#define g_soilH
//---------------------------------------------------------------------------



using namespace std;

double g_soil(
	double	rho,	/* snow layer's density (kg/m^3)	     */
	double	tsno,	/* snow layer's temperature (K)		     */
	double	tg,	/* soil temperature (K)			     */
	double	ds,	/* snow layer's thickness (m)		     */
	double	dg,	/* dpeth of soil temperature measurement (m) */
	double	pa);	/* air pressure (Pa)			     */

        #endif
 