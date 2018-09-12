/* Conversion functions using curve isomorphisms
 *
 * These functions rely on the birational equivalence of the curves
 *    Wei25519 <-> Curve25519 <-> Ed25519
 *
 * Daniel Beer <dlbeer@gmail.com>, 18 Jan 2014
 * Nikolas Rösener <nroesener@uni-bremen.de> 1 Aug 2018
 *
 */

#include "morph25519.h"
#include "f25519.h"

const uint8_t f25519_A[F25519_SIZE] = {             // = 486662 mod 255^19
	0x06, 0x6d, 0x07, 0x00, 0x00, 0x00, 0x00, 0x00,
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
};
const uint8_t f25519_three[F25519_SIZE] = {3};      // =  3 mod 255^19
const uint8_t f25519_delta[F25519_SIZE] = {         // = (255^19 + A) / 3 mod 255^19
	0x51, 0x24, 0xad, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa,
	0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa,
	0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa,
	0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0x2a
};
const uint8_t f25519_c[F25519_SIZE] = {             // = sqrt(-(A + 2)) mod 255^19
	0xe7, 0x81, 0xba, 0x00, 0x55, 0xfb, 0x91, 0x33,
	0x7d, 0xe5, 0x82, 0xb4, 0x2e, 0x2c, 0x5e, 0x3a,
	0x81, 0xb0, 0x03, 0xfc, 0x23, 0xf7, 0x84, 0x2d,
	0x44, 0xf9, 0x5f, 0x9f, 0x0b, 0x12, 0xd9, 0x70
};

void morph25519_ey2mx(uint8_t *mx, const uint8_t *ey)
{
	uint8_t yplus[F25519_SIZE];
	uint8_t yminus[F25519_SIZE];

	f25519_sub(yplus, f25519_one, ey);      // y+ =           1 - ey
	f25519_inv__distinct(yminus, yplus);    // y- =          (1 - ey)^-1
	f25519_add(yplus, f25519_one, ey);      // y+ =  1 + ey
	f25519_mul__distinct(mx, yplus, yminus);// mx = (1 + ey) * (1 - ey)^-1
	f25519_normalize(mx);                   // mx = (1 + ey) * (1 - ey)^-1 (mod p)
}

void morph25519_mx2ey(uint8_t *ey, const uint8_t *mx)
{
	uint8_t n[F25519_SIZE];
	uint8_t d[F25519_SIZE];

	f25519_add(n, mx, f25519_one);  //  n =             mx + 1
	f25519_inv__distinct(d, n);     //  d =            (mx + 1)^-1
	f25519_sub(n, mx, f25519_one);  //  n =  mx - 1
	f25519_mul__distinct(ey, n, d); // ey = (mx - 1) * (mx + 1)^-1
}

uint8_t morph25519_ey2ex(uint8_t *x, const uint8_t *y, int parity)
{
	static const uint8_t d[F25519_SIZE] = {
		0xa3, 0x78, 0x59, 0x13, 0xca, 0x4d, 0xeb, 0x75,
		0xab, 0xd8, 0x41, 0x41, 0x4d, 0x0a, 0x70, 0x00,
		0x98, 0xe8, 0x79, 0x77, 0x79, 0x40, 0xc7, 0x8c,
		0x73, 0xfe, 0x6f, 0x2b, 0xee, 0x6c, 0x03, 0x52
	};

	uint8_t a[F25519_SIZE];
	uint8_t b[F25519_SIZE];
	uint8_t c[F25519_SIZE];

	/* Compute c = y^2 */
	f25519_mul__distinct(c, y, y);

	/* Compute b = (1+dy^2)^-1 */
	f25519_mul__distinct(b, c, d);
	f25519_add(a, b, f25519_one);
	f25519_inv__distinct(b, a);

	/* Compute a = y^2-1 */
	f25519_sub(a, c, f25519_one);

	/* Compute c = a*b = (y^2-1)/(1+dy^2) */
	f25519_mul__distinct(c, a, b);

	/* Compute a, b = +/-sqrt(c), if c is square */
	f25519_sqrt(a, c);
	f25519_neg(b, a);

	/* Select one of them, based on the parity bit */
	f25519_select(x, a, b, (a[0] ^ parity) & 1);

	/* Verify that x^2 = c */
	f25519_mul__distinct(a, x, x);
	f25519_normalize(a);
	f25519_normalize(c);

	return f25519_eq(a, c);
}

uint8_t morph25519_wx2wy(uint8_t *wy, const uint8_t *wx, int sign)
{
	static const uint8_t a[F25519_SIZE] = {
		0x44, 0xa1, 0x14, 0x49, 0x98, 0xaa, 0xaa, 0xaa,
		0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa,
		0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa,
		0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0x2a
	};  // curve parameter a_4 = a
	static const uint8_t b[F25519_SIZE] = {
		0x64, 0xc8, 0x10, 0x77, 0x9c, 0x5e, 0x0b, 0x26,
		0xb4, 0x97, 0xd0, 0x5e, 0x42, 0x7b, 0x09, 0xed,
		0x25, 0xb4, 0x97, 0xd0, 0x5e, 0x42, 0x7b, 0x09,
		0xed, 0x25, 0xb4, 0x97, 0xd0, 0x5e, 0x42, 0x7b
	};  // curve parameter a_6 = b

	// y = +/- SQRT(x³ + ax + b)
	uint8_t T1[F25519_SIZE];
	uint8_t T2[F25519_SIZE];
	uint8_t T3[F25519_SIZE];

	/* Compute T2 = x^3 */
	f25519_mul__distinct(T1, wx, wx);
	f25519_mul__distinct(T2, T1, wx);

	/* Compute T1 = ax */
	f25519_mul__distinct(T1, a, wx);

	/* Compute T3 = x^3 + ax */
	f25519_add(T3, T2, T1);

	/* Compute T1 = x^3 + ax + b */
	f25519_add(T1, T3, b);

	/* Compute T2, T3 = +/-sqrt(T1), if T1 is square */
	f25519_sqrt(T2, T1);
	f25519_neg(T3, T2);

	/* Select one of them, based on the sign bit */
	f25519_copy(wy, f25519_zero);
	f25519_select(wy, T2, T3, sign);

	/* Verify that T2 = wy^2 == T1 */
	f25519_mul__distinct(T2, wy, wy);
	f25519_normalize(T1);
	f25519_normalize(T2);

	return f25519_eq(T1, T2);
}

void morph25519_montgomery_recovery(
			uint8_t *xQ, uint8_t *yQ, uint8_t *zQ,
			const uint8_t *xP, const uint8_t *yP,
			const uint8_t *XQ, const uint8_t *ZQ,
			const uint8_t *xD, const uint8_t *zD)
{
	static const uint32_t A2 = 973324; 	// 2 * A
	static const uint32_t B2 = 2;		// 2 * B

	uint8_t v1[F25519_SIZE], v2[F25519_SIZE];
	uint8_t	v3[F25519_SIZE], v4[F25519_SIZE];

	f25519_mul(v1, xP, ZQ);   // 1 v1 ← xP · ZQ 	1M
	f25519_add(v2, XQ, v1);	  // 2 v2 ← XQ + v1 	1a
	f25519_sub(v3, XQ, v1);	  // 3 v3 ← XQ − v1 	1s
	f25519_mul(v3, v3, v3);	  // 4 v3 ← v3^2 		1S
	f25519_mul(v3, v3, xD);	  // 5 v3 ← v3 · X⊕ 	1M
	f25519_mul_c(v1, ZQ, A2); // 6 v1 ← 2A · ZQ 	1c
	f25519_add(v2, v2, v1);	  // 7 v2 ← v2 + v1 	1a
	f25519_mul(v4, xP, XQ);	  // 8 v4 ← xP · XQ 	1M
	f25519_add(v4, v4, ZQ);	  // 9 v4 ← v4 + ZQ 	1a
	f25519_mul(v2, v2, v4);	  // 10 v2 ← v2 · v4 	1M
	f25519_mul(v1, v1, ZQ);	  // 11 v1 ← v1 · ZQ 	1M
	f25519_sub(v2, v2, v1);	  // 12 v2 ← v2 − v1 	1s
	f25519_mul(v2, v2, zD);	  // 13 v2 ← v2 · Z⊕ 	1M
	f25519_sub(yQ, v2, v3);	  // 14 Y′ ← v2 − v3 	1s
	f25519_mul_c(v1, yP, B2); // 15 v1 ← 2B · yP 	1c
	f25519_mul(v1, v1, ZQ);	  // 16 v1 ← v1 · ZQ 	1M
	f25519_mul(v1, v1, zD);	  // 17 v1 ← v1 · Z⊕ 	1M
	f25519_mul(xQ, v1, XQ);	  // 18 X′ ← v1 · XQ 	1M
	f25519_mul(zQ, v1, ZQ);	  // 19 Z′ ← v1 · ZQ 	1M
	// 20 return (X′ : Y′ : Z′)
}

uint8_t morph25519_mx2e(uint8_t *ex, uint8_t *ey,
			   const uint8_t *mx, int parity)
{
	uint8_t ok;

	morph25519_mx2ey(ey, mx);
	ok = morph25519_ey2ex(ex, ey, parity);

	f25519_normalize(ex);
	f25519_normalize(ey);

	return ok;
}

void morph25519_wx2mx(uint8_t* mx, const uint8_t* wx)
{
	/*
		The following code calculates:
		wx == 0 ? 0 : (wx - delta)
	*/
	uint8_t  tmp[F25519_SIZE];
	f25519_sub(tmp, wx, f25519_delta);
	f25519_normalize(tmp);
	f25519_select(mx, tmp, f25519_zero, f25519_eq(wx, f25519_zero));
}

void morph25519_mx2wx(uint8_t* wx, const uint8_t* mx)
{
	/*
		The following code calculates:
		mx == 0 ? 0 : (mx + delta)
	*/
	uint8_t  tmp[F25519_SIZE];
	f25519_add(tmp, mx, f25519_delta);
	f25519_normalize(tmp);
	f25519_select(wx, tmp, f25519_zero, f25519_eq(mx, f25519_zero));
}

void morph25519_m2w(uint8_t* wx, uint8_t* wy, const uint8_t* mx, const uint8_t* my)
{
	morph25519_mx2wx(wx, mx);
	f25519_copy(wy, my);
}

void morph25519_w2m(uint8_t* mx, uint8_t* my, const uint8_t* wx, const uint8_t* wy)
{
	morph25519_wx2mx(mx, wx);
	f25519_copy(my, wy);
}

void morph25519_e2w(uint8_t* wx, uint8_t* wy, const uint8_t* ex, const uint8_t* ey)
{
	/*
		The following code calculates:
		wx = (1 + ey) / ((1 - ey) + delta)   (mod p)
		wy = (c * (1 + ey)) / (1 - ey) * ex  (mod p)
	*/
	uint8_t nom[F25519_SIZE];   // nominator
	uint8_t den[F25519_SIZE];   // denominator
	uint8_t inv[F25519_SIZE];   // inversion result
	uint8_t mul[F25519_SIZE];   // multiplication result

	f25519_add(nom, f25519_one, ey);    // nom =   1 + ey
	f25519_sub(den, f25519_one, ey);    // den =              1 - ey
	f25519_inv__distinct(inv, den);     // inv =             (1 - ey)^-1
	f25519_mul__distinct(mul, nom, inv);// mul =  (1 + ey) * (1 - ey)^-1
	f25519_add(wx, mul, f25519_delta);  //  wx = ((1 + ey) * (1 - ey)^-1) + delta
	f25519_normalize(wx);        		//  wx = ((1 + ey) * (1 - ey)^-1) + delta  (mod p)

	f25519_mul__distinct(mul, f25519_c, nom);	// mul =  c * (1 + ey)
	f25519_mul__distinct(inv, den, ex);			// inv =                   (1 - ey) * ex
	f25519_inv__distinct(den, inv);          	// den =                  ((1 - ey) * ex)^-1
	f25519_mul__distinct(wy, mul, den);      	//  wy = (c * (1 + ey)) * ((1 - ey) * ex)^-1
	f25519_normalize(wy);        				//  wy = (c * (1 + ey)) * ((1 - ey) * ex)^-1  (mod p)
}

void morph25519_w2e(uint8_t* ex, uint8_t* ey, const uint8_t* mx, const uint8_t* my)
{
	/*
		The following code calculates:
		pa = 3 * p.x - A
		ex = (c * pa) / (3 * my)
		ey = (pa - 3) / (pa + 3)
	*/
	uint8_t  pa[F25519_SIZE]; // intermediate result
	uint8_t nom[F25519_SIZE]; // nominator
	uint8_t den[F25519_SIZE]; // denominator
	uint8_t inv[F25519_SIZE]; // inverted denominator

	f25519_mul_c(inv, mx, 3);					// inv = 3 * mx
	f25519_sub(pa, inv, f25519_A);              // pa  = 3 * mx - A

	f25519_mul__distinct(nom, f25519_c, pa);    // nom =  c * pa
	f25519_mul_c(den, my, 3);					// den =             3 * my
	f25519_inv__distinct(inv, den);             // inv =            (3 * my)^-1
	f25519_mul__distinct(ex, nom, inv);         // ex  = (c * pa) * (3 * my)^-1
	f25519_normalize(ex);                       // ex  = (c * pa) * (3 * my)^-1 (mod p)

	f25519_sub(nom, pa, f25519_three);          // nom =  pa - 3
	f25519_add(den, pa, f25519_three);          // den =             pa + 3
	f25519_inv__distinct(inv, den);             // inv =            (pa + 3)^-1
	f25519_mul__distinct(ey, nom, inv);         //  ey = (pa - 3) * (pa + 3)^-1
	f25519_normalize(ey);                       //  ey = (pa - 3) * (pa + 3)^-1 (mod p)
}

void morph25519_e2m(uint8_t* mx, uint8_t* my, const uint8_t* ex, const uint8_t* ey)
{
	/*
		The following code calculates:
		mx = (1 + ey) / (1 - ey) (mod p)
		my = c * (1 + ey) / (1 - ey) * ex  (mod p)
	*/
	uint8_t nom[F25519_SIZE];   // nominator
	uint8_t den[F25519_SIZE];   // denominator
	uint8_t inv[F25519_SIZE];   // inversion result

	f25519_add(nom, f25519_one, ey);    // nom =   1 + ey
	f25519_sub(den, f25519_one, ey);    // den =              1 - ey
	f25519_inv__distinct(inv, den);     // inv =             (1 - ey)^-1
	f25519_mul__distinct(mx, nom, inv); // mul =  (1 + ey) * (1 - ey)^-1
	f25519_normalize(mx);        		//  mx = ((1 + ey) * (1 - ey)^-1) (mod p)

	f25519_mul__distinct(inv, den, ex);			// inv =                 (1 - ey) * ex
	f25519_inv__distinct(den, inv);     		// den =                ((1 - ey) * ex)^-1
	f25519_mul__distinct(inv, f25519_c, nom); 	// inv = c * (1 + ey)
	f25519_mul__distinct(my, inv, den); 		//  my = c * (1 + ey) * ((1 - ey) * ex)^-1
	f25519_normalize(my);        				//  my = c * (1 + ey) * ((1 - ey) * ex)^-1  (mod p)
}

void morph25519_m2e(uint8_t* ex, uint8_t* ey, const uint8_t* mx, const uint8_t* my)
{
	/*
		The following code calculates:
		ex = (c * mx) * my^-1
		ey = mx-1 * mx+1
	*/
	uint8_t nom[F25519_SIZE]; // nominator
	uint8_t den[F25519_SIZE]; // denominator
	uint8_t inv[F25519_SIZE]; // inverted denominator

	f25519_mul__distinct(nom, f25519_c, mx);    // nom =  c * mx
	f25519_inv__distinct(inv, my);              // inv =      my^-1
	f25519_mul__distinct(ex, nom, inv);         // ex  = mx * my^-1
	f25519_normalize(ex);                       // ex  = mx * my^-1 (mod p)

	f25519_sub(nom, mx, f25519_one);            // nom =  mx - 1
	f25519_add(den, mx, f25519_one);            // den =             mx + 1
	f25519_inv__distinct(inv, den);             // inv =            (mx + 1)^-1
	f25519_mul__distinct(ey, nom, inv);         //  ey = (mx - 1) * (mx + 1)^-1
	f25519_normalize(ey);                       //  ey = (mx - 1) * (mx + 1)^-1 (mod p)
}