/* Elliptic Curve Digital Signature Algorithm
 * Nikolas Rösener <nroesener@uni-bremen.de>, 01.08.2018
 *
 */
#include "ed25519.h"
#include "fprime.h"
#include "ecdsa.h"
#include "morph25519.h"

static const uint8_t n[FPRIME_SIZE] = {
	0xed, 0xd3, 0xf5, 0x5c, 0x1a, 0x63, 0x12, 0x58,
	0xd6, 0x9c, 0xf7, 0xa2, 0xde, 0xf9, 0xde, 0x14,
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x10
};

static void rshift(uint8_t* A, int t){
	for (int j = 0; j < t; j++) {
		int c = 0;
		uint8_t cOld = 0;
		for (int i = FPRIME_SIZE; i--;)
		{
			c = A[i]&0x1;
			A[i] = A[i]>>1 | cOld<<7;
			cOld = c;
		}
	}
}

void ecdsa_pubkey(uint8_t *wx, uint8_t *wy, const uint8_t *secret)
{
	struct ed25519_pt p1;
	uint8_t ex[F25519_SIZE], ey[F25519_SIZE];
	ed25519_smult(&p1, &ed25519_base, secret);
	ed25519_unproject(ex, ey, &p1);
	morph25519_e2w(wx, wy, ex, ey);
}

uint8_t ecdsa_sign(uint8_t *r, uint8_t *s, const uint8_t *d,
		 const uint8_t *e, const uint8_t *k)
{
	uint8_t t[FPRIME_SIZE], z[FPRIME_SIZE];

	if (fprime_eq(k, fprime_zero))
		return 0;

	struct ed25519_pt p1;
	uint8_t ex[F25519_SIZE], ey[F25519_SIZE];
	uint8_t wx[F25519_SIZE], wy[F25519_SIZE];
	// 4. Calculate the curve point (x_1, y_1) = k * G.
	ed25519_smult(&p1, &ed25519_base, k);
	ed25519_unproject(ex, ey, &p1);
	morph25519_e2w(wx, wy, ex, ey);

	// 5. Calculate r = x_1 \pmod{n}.
	fprime_from_bytes(r, wx, F25519_SIZE, n);

	// 5. If r = 0, go back to step 3.
	if (fprime_eq(r, fprime_zero))
		return 0;

	// 6. Calculate s = k^{-1}(z + r d) \pmod{n}:
	// 6. r * d
	fprime_mul(t, r, d, n);
	//fprime_normalize(t,n);

	// 6. z + (r d)
	fprime_copy(z, e);
	rshift(z,3);
	fprime_add(z, t, n);

	// 6. k^{-1}
	fprime_inv(t, k, n);

	// 6. s = (k^{-1}) (z + (r d))
	fprime_mul(s, t, z, n);
	fprime_normalize(s, n);

	// 6. If s = 0, go back to step 3.
	if (fprime_eq(s, fprime_zero))
		return 0;

	return 1;
}

uint8_t ecdsa_verify(const uint8_t *x, const uint8_t *y,
		      const uint8_t *e, const uint8_t *r, const uint8_t *s)
{
	struct ed25519_pt p1;
	struct ed25519_pt p2;
	struct ed25519_pt Q;
	uint8_t w[FPRIME_SIZE], z[FPRIME_SIZE];
	uint8_t u1[FPRIME_SIZE], u2[FPRIME_SIZE];
	uint8_t ex[F25519_SIZE], ey[F25519_SIZE];
	uint8_t wx[F25519_SIZE], wy[F25519_SIZE];

	// 3. Let z be the L_n leftmost bits of e
	fprime_copy(z, e);
	rshift(z,3);

	// 4. Calculate w = s^-1 mod n
	fprime_inv(w, s, n);

	// 5. Calculate u_1 = zw mod n
	fprime_mul(u1, z, w, n);

	// and  u_2 = r*w mod n
	fprime_mul(u2, r, w, n);

	// 5. Calculate the curve point (x_1, y_1) = u_1 * G + u_2 * Q_A.
	// tmp1 = u_1 * G
	morph25519_w2e(ex, ey, x, y);
	ed25519_project(&Q, ex, ey);
	ed25519_smult(&p1, &ed25519_base, u1);
	ed25519_smult(&p2, &Q, u2);
	ed25519_add(&Q, &p1, &p2);
	ed25519_unproject(ex, ey, &Q);
	morph25519_e2w(wx, wy, ex, ey);

	// 7. The signature is valid of r == x1 mod n
	fprime_normalize(wx, n);
	return f25519_eq(wx, r);
}
