/* Elliptic Curve Digital Signature Algorithm
 * Nikolas Rösener <nroesener@uni-bremen.de>, 01.08.2018
 *
 */

#ifndef ECDSA_H_
#define ECDSA_H_

#include <stdint.h>
#include <stddef.h>
#include "fprime.h"

/*
 * Generate an ecdsa public key (x,y of an affine point on Wei25519)
 *
 * input:
 *	secret: a secret in the interval [1, n-1]
 *
 * output:
 * 	wx: x-coordinate of the point
 *  wy: y-coordinate of the point
 */
void ecdsa_pubkey(uint8_t *wx, uint8_t *wy, const uint8_t *secret);

/**
 * Calculate the ecdsa signature.
 *
 * For a description of this algorithm see
 * https://en.wikipedia.org/wiki/Elliptic_Curve_DSA#Signature_generation_algorithm
 *
 * input:
 *  d: private key on the curve wei25519 (32 bytes)
 *  e: hash to sign (32 bytes)
 *  k: random data, this must be changed for every signature (32 bytes)
 *
 * output:
 *  r: r value of the signature (36 bytes)
 *  s: s value of the signature (36 bytes)
 *
 * return:
 *   1: everything is ok
 *   0: can not create signature, try again with different k.
 */
uint8_t ecdsa_sign(uint8_t *r, uint8_t *s, const uint8_t *d,
		 const uint8_t *e, const uint8_t *k);

/**
 * Verifies a ecdsa signature.
 *
 * For a description of this algorithm see
 * https://en.wikipedia.org/wiki/Elliptic_Curve_DSA#Signature_verification_algorithm
 *
 * input:
 *  x: x coordinate of the public key (32 bytes)
 *  y: y coordinate of the public key (32 bytes)
 *  e: hash to verify the signature of (32 bytes)
 *  r: r value of the signature (32 bytes)
 *  s: s value of the signature (32 bytes)
 *
 * return:
 *  1: signature is ok
 *  0: signature check failed the signature is invalid
 */
uint8_t ecdsa_verify(const uint8_t *x, const uint8_t *y,
		      const uint8_t *e, const uint8_t *r, const uint8_t *s);

#endif
