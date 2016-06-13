/*
* Copyright 2016 naehrwert
* Licensed under the terms of the GNU GPL, version 2
* http://www.gnu.org/licenses/old-licenses/gpl-2.0.txt
*/

#ifndef _PAIRING_H_
#define _PAIRING_H_

#include "ec_pqr.h"


/*! Comparison results. */
/*! Less. */
#define POINT_CMP_L (-1)
/*! Greater. */
#define POINT_CMP_G (1)
/*! Equal. */
#define POINT_CMP_E (0)

/*!
* \brief Compute line.
*/
poly_t *_pairing_line(poly_t *d, ec_pqr_point_t *P, ec_pqr_point_t *R, ec_pqr_point_t *Q, ec_pqr_group_t *ecg);

/*!
* \brief Miller's algorithm.
*/
poly_t *_pairing_miller(poly_t *d, ec_pqr_point_t *P, ec_pqr_point_t *Q, bn_t *n, ec_pqr_group_t *ecg);

/*!
* \brief Weil pairing.
*/
poly_t *pairing_weil(poly_t *d, ec_pqr_point_t *P, ec_pqr_point_t *Q, bn_t *n, ec_pqr_group_t *ecg);

#endif // _PAIRING_H_
