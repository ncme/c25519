/* Conversion functions using curve isomorphisms
 *
 * These functions rely on the birational equivalence of the curves
 *    Wei25519 <-> Curve25519 <-> Ed25519
 *
 * Daniel Beer <dlbeer@gmail.com>, 18 Jan 2014
 * Nikolas Rösener <nroesener@uni-bremen.de> 1 Aug 2018
 *
 */

#ifndef MORPH25519_H_
#define MORPH25519_H_

#include <stdint.h>

/* Convert an Edwards Y to a Montgomery X (Edwards X is not used).
 * Resulting coordinate is normalized.
 */
void morph25519_e2m(uint8_t *mx, const uint8_t *ey);

/* Return a parity bit for the Edwards X coordinate */
static inline int morph25519_eparity(const uint8_t *edwards_x)
{
	return edwards_x[0] & 1;
}

void morph25519_montgomery_recovery(
			uint8_t *xQ, uint8_t *yQ, uint8_t *zQ,
			const uint8_t *xP, const uint8_t *yP,
			const uint8_t *XQ, const uint8_t *ZQ,
			const uint8_t *xD, const uint8_t *zD);

/* Convert a Montgomery X and a parity bit to an Edwards X/Y. Returns
 * non-zero if successful.
 */
uint8_t morph25519_m2e(uint8_t *ex, uint8_t *ey,
			   const uint8_t *mx, int parity);

/*
 * Transforms an affine point on the Montgomery curve Curve25519
 * to an affine point on the short Weierstrass curve Wei25519.
 */
void morph25519_m2w(uint8_t* wx, uint8_t* wy, const uint8_t* mx, const uint8_t* my);

/*
 * Transforms an affine point on the short Weierstrass curve Wei25519
 * to an affine point on the Montgomery curve Curve25519.
 */
void morph25519_w2m(uint8_t* mx, uint8_t* my, const uint8_t* wx, const uint8_t* wy);

/*
 * Transforms an affine point on the Edwards curve Ed25519
 * to an affine point on the short Weierstrass curve Wei25519.
 */
void morph25519_e2w(uint8_t* wx, uint8_t* wy, const uint8_t* ex, const uint8_t* ey);

/*
 * Transforms an affine point on the short Weierstrass curve Wei25519
 * to an affine point on the Edwards curve Ed25519.
 */
void morph25519_w2e(uint8_t* ex, uint8_t* ey, const uint8_t* wx, const uint8_t* wy);

#endif
