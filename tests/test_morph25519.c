/* Montgomery <-> Edwards isomorphism
 * Daniel Beer <dlbeer@gmail.com>, 18 Jan 2014
 *
 * This file is in the public domain.
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "morph25519.h"
#include "c25519.h"
#include "ed25519.h"

static void print_elem(const uint8_t *e)
{
	int i;

	for (i = 0; i < F25519_SIZE; i++)
		printf("%02x", e[i]);
}

static void test_morph_wx2wy()
{
	static const uint8_t Gx[F25519_SIZE] = {
		0x5a, 0x24, 0xad, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa,
		0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa,
		0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa,
		0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0x2a
	};  // the x coordinate of the base point
	static const uint8_t Gy[F25519_SIZE] = {
		0xd9, 0xd3, 0xce, 0x7e, 0xa2, 0xc5, 0xe9, 0x29,
		0xb2, 0x61, 0x7c, 0x6d, 0x7e, 0x4d, 0x3d, 0x92,
		0x4c, 0xd1, 0x48, 0x77, 0x2c, 0xdd, 0x1e, 0xe0,
		0xb4, 0x86, 0xa0, 0xb8, 0xa1, 0x19, 0xae, 0x20
	};  // the y coordinate of the base point

	uint8_t Gy2[F25519_SIZE];
	assert(morph25519_wx2wy(Gy2, Gx, morph25519_eparity(Gy)));
	printf("  ");
	print_elem(Gy);
	printf(" [%d] ~ \n  ", morph25519_eparity(Gy));
	print_elem(Gy2);
	printf("\n");
	assert(f25519_eq(Gy, Gy2));
}

static void test_morph_e2w(const uint8_t *ex, const uint8_t *ey)
{
	uint8_t w1x[F25519_SIZE], w2x[F25519_SIZE];
	uint8_t w1y[F25519_SIZE], w2y[F25519_SIZE], w3y[F25519_SIZE];

	uint8_t e1x[F25519_SIZE], e1y[F25519_SIZE];
	morph25519_e2w(w1x, w1y, ex, ey);
	morph25519_w2e(e1x, e1y, w1x, w1y);
	assert(f25519_eq(e1x, ex));
	assert(f25519_eq(e1y, ey));

	uint8_t mx[F25519_SIZE], my[F25519_SIZE];
	morph25519_ey2mx(mx, ey);
	morph25519_m2w(w2x, w2y, mx, my);

	assert(morph25519_wx2wy(w3y, w1x, !morph25519_eparity(w1y)));

	assert(f25519_eq(w1x, w2x));
	//assert(f25519_eq(w1y, w2y));
	printf("  ");
	print_elem(w1y);
	printf(" [%d] ~ \n  ", !morph25519_eparity(w1y));
	print_elem(w3y);
	printf("\n");
	assert(f25519_eq(w1y, w3y));
}

static void test_morph(const uint8_t *mx,
		       const uint8_t *ex, const uint8_t *ey)
{
	const int parity = morph25519_eparity(ex);
	uint8_t mx_test[F25519_SIZE];
	uint8_t ex_test[F25519_SIZE];
	uint8_t ey_test[F25519_SIZE];

	printf("  ");
	print_elem(mx);
	printf(" [%d]\n    ~ (", parity);
	print_elem(ex);
	printf(", ");
	print_elem(ey);
	printf(")\n");

	morph25519_ey2mx(mx_test, ey);
	morph25519_mx2e(ex_test, ey_test, mx, parity);

	assert(f25519_eq(mx_test, mx));
	assert(f25519_eq(ex_test, ex));
	assert(f25519_eq(ey_test, ey));
}

static void test_sm(void)
{
	uint8_t e[C25519_EXPONENT_SIZE];
	uint8_t mx[F25519_SIZE];
	uint8_t ex[F25519_SIZE];
	uint8_t ey[F25519_SIZE];
	struct ed25519_pt p;
	unsigned int i;

	for (i = 0; i < sizeof(e); i++)
		e[i] = random();

	c25519_prepare(e);
	c25519_smult(mx, c25519_base_x, e);

	ed25519_smult(&p, &ed25519_base, e);
	ed25519_unproject(ex, ey, &p);

	test_morph(mx, ex, ey);
}

int main(void)
{
	int i;

	srandom(0);

	printf("test_base\n");
	test_morph(c25519_base_x, ed25519_base.x, ed25519_base.y);

	printf("test_sm\n");
	for (i = 0; i < 32; i++)
		test_sm();

	printf("test_morph_wx2wy\n");
	test_morph_wx2wy();

	printf("test_morph_e2w\n");
	test_morph_e2w(ed25519_base.x, ed25519_base.y);

	return 0;
}
