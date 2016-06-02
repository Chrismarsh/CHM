#ifndef _PHYSUNITS_H
#define _PHYSUNITS_H

/*
 *  Units of physical measurement
 */

/* ------------------------------------------------------------------------ */

/*
 *  convert Celsius to Kelvin
 */
#define C_TO_K(c)		((c) + FREEZE)

/*
 *  convert Kelvin to Celsius
 */
#define K_TO_C(k)		((k) - FREEZE)

/*
 *  Convert kilograms to grams.
 */
#define KG_TO_G(kg)		((kg) * 1000.0)

/*
 *  Convert grams to kilograms.
 */
#define G_TO_KG(g)		((g) * 0.001)

/*
 *  Convert Joules to calories.
 */
#define J_TO_CAL(j)		((j) * 0.238846)

/*
 *  Convert calories to Joules
 */
#define CAL_TO_J(c)		((c) * 4.186798188)

/*
 *  convert wavelength (um) to wave number (1/cm)
 */
#define WAVENO(x)                       ( (int)(10000. / (x) + .5) )

/*
 *  convert wave number (1/cm) to wavelength (um)
 */
#define WAVELEN(nu)                     ( 10000. / (nu) )


/* ------------------------------------------------------------------------ */

#endif  /* _PHYSUNITS_H */
