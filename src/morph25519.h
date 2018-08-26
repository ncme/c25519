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

/*
 * Transforms the y-coordinate of a point on the Edwards curve Ed25519
 * to the x-coordinate of a point on the Montgomery curve Curve25519.
 */
void morph25519_ey2mx(uint8_t *mx, const uint8_t *ey);

/*
 * Transforms the x-coordinate of a point on the Montgomery curve Curve25519
 * to the y-coordinate of a point on the Edwards curve Ed25519.
 *
 * Return non-zero if successful.
 */
void morph25519_mx2ey(uint8_t *ey, const uint8_t *mx);

/*
 * Recovers the x-coordinate of the Edwards curve Ed25519
 * from the y-coordinate and a parity bit.
 *
 * Return non-zero if successful.
 */
uint8_t morph25519_ey2ex(uint8_t *x, const uint8_t *y, int parity);

/*
 * Recovers the x-coordinate of the short Weierstrass curve Wei25519
 * from the y-coordinate and a sign bit.
 *
 * Return non-zero if successful.
 */
uint8_t morph25519_wx2wy(uint8_t *wy, const uint8_t *wx, int sign);

/* Return a parity bit for the Edwards X coordinate */
static inline int morph25519_eparity(const uint8_t *edwards_x)
{
	return edwards_x[0] & 1;
}

/*
 * Okeya–Sakurai y-coordinate recovery
 * Input:
 * 	(xP : yP : 1) = P,
 * 	(XQ : ZQ) = x(Q),
 * 	(X⊕ : Z⊕) = x(P ⊕ Q) for P and Q in E(A,B)(Fq)
 * 	 with P /∈ E(A,B)[2] and Q /∈ {P, ⊖P, O}.
 * Output:
 * 	(X′ : Y′ : Z′) = Q
 * Cost:
 * 	10M + 1S + 2c + 3a + 3s
 */
void morph25519_montgomery_recovery(
			uint8_t *xQ, uint8_t *yQ, uint8_t *zQ,
			const uint8_t *xP, const uint8_t *yP,
			const uint8_t *XQ, const uint8_t *ZQ,
			const uint8_t *xD, const uint8_t *zD);

/*
 * Transforms the x-coordinate of a point on the Montgomery curve Curve25519
 * to the x-coordinate of a point on the short Weierstrass curve Wei25519.
 */
void morph25519_mx2wx(uint8_t* wx, const uint8_t* mx);

/*
 * Transforms the x-coordinate of a point on the short Weierstrass curve Wei25519
 * to the x-coordinate of a point on the Montgomery curve Curve25519.
 */
void morph25519_wx2mx(uint8_t* mx, const uint8_t* wx);

/*
 * Transforms the x-coordinate of a point on the Montgomery curve Curve25519
 * to x- and y-coordinate of the Edwards curve Ed25519.
 *
 * Return non-zero if successful.
 */
uint8_t morph25519_mx2e(uint8_t *ex, uint8_t *ey,
			   const uint8_t *mx, int parity);

/*
 * Transforms an affine point on the Montgomery curve Curve25519
 * to an affine point on the short Weierstrass curve Wei25519.
 *
 * Input:
 * 	(MX, MY) on Curve25519 - not (0,0) of order two or the point at infinity!
 * Output:
 *  (WX, WY) on Wei25519
 */
void morph25519_m2w(uint8_t* wx, uint8_t* wy, const uint8_t* mx, const uint8_t* my);

/*
 * Transforms an affine point on the short Weierstrass curve Wei25519
 * to an affine point on the Montgomery curve Curve25519.
 *
 * Input:
 *  (WX, WY) on Wei25519 - not (A/3,0) or the point at infinity!
 * Output:
 *  (MX, MY) on Curve25519
 */
void morph25519_w2m(uint8_t* mx, uint8_t* my, const uint8_t* wx, const uint8_t* wy);

/*
 * Transforms an affine point on the Edwards curve Ed25519
 * to an affine point on the short Weierstrass curve Wei25519.
 *
 * Input:
 * 	(EX, EY) on Ed25519 - not (0,1) or (0,-1)
 * Output:
 *  (WX, WY) on Wei25519
 */
void morph25519_e2w(uint8_t* wx, uint8_t* wy, const uint8_t* ex, const uint8_t* ey);

/*
 * Transforms an affine point on the short Weierstrass curve Wei25519
 * to an affine point on the Edwards curve Ed25519.
 *
 * Input:
 * 	(WX, WY) on Wei25519 - not (A/3,0) or the point at infinity!
 * Output:
 *  (EX, EY) on Ed25519
 */
void morph25519_w2e(uint8_t* ex, uint8_t* ey, const uint8_t* mx, const uint8_t* my);

/*
 * Transforms an affine point on the Edwards curve Ed25519
 * to an affine point on the short Weierstrass curve Wei25519.
 *
 * Input:
 * 	(EX, EY) on Ed25519 - not the neutral point (0,1) or (0, -1)
 * Output:
 *  (MX, MY) on Curve25519
 */
void morph25519_e2m(uint8_t* mx, uint8_t* my, const uint8_t* ex, const uint8_t* ey);

/*
 * Transforms an affine point on the short Weierstrass curve Wei25519
 * to an affine point on the Edwards curve Ed25519.
 *
 * Input:
 * 	(MX, MY) on Curve25519 - not (0,0) of order two or the point at infinity!
 * Output:
 *  (WX, WY) on Wei25519
 */
void morph25519_m2e(uint8_t* ex, uint8_t* ey, const uint8_t* mx, const uint8_t* my);

#endif
