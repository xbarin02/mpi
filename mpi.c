#include "mpi.h"
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>

void mpi_init(mpi_t rop)
{
	rop->nmemb = 0;
	rop->data = NULL;
}

void mpi_clear(mpi_t rop)
{
	free(rop->data);
}

void mpi_enlarge(mpi_t rop, size_t nmemb)
{
	if (nmemb > rop->nmemb) {
		size_t min = rop->nmemb;

		rop->nmemb = nmemb;

		rop->data = realloc(rop->data, nmemb * sizeof(uint32_t));

		if (rop->data == NULL && nmemb != 0) {
			fprintf(stderr, "Out of memory (%zu words requested)\n", nmemb);
			abort();
		}

		for (size_t n = min; n < nmemb; ++n) {
			rop->data[n] = 0;
		}
	}
}

void mpi_compact(mpi_t rop)
{
	size_t nmemb = rop->nmemb;

	if (rop->nmemb == 0) {
		return;
	}

	for (size_t i = rop->nmemb - 1; i != (size_t)-1; --i) {
		if (rop->data[i] == 0) {
			nmemb--;
		} else {
			break;
		}
	}

	assert(nmemb != (size_t)-1);

	rop->data = realloc(rop->data, nmemb * sizeof(uint32_t));
	rop->nmemb = nmemb;

	if (rop->data == NULL && nmemb != 0) {
		fprintf(stderr, "Out of memory (%zu words requested)\n", nmemb);
		abort();
	}
}

static size_t ceil_div(size_t n, size_t d)
{
	return (n + d) / d;
}

void mpi_set_u64(mpi_t rop, uint64_t op)
{
	size_t nmemb = ceil_div(64, 31);

	mpi_enlarge(rop, nmemb);

	for (size_t n = 0; n < nmemb; ++n) {
		rop->data[n] = op & 0x7fffffff;
		op >>= 31;
	}

	for (size_t n = nmemb; n < rop->nmemb; ++n) {
		rop->data[n] = 0;
	}
}

void mpi_set_u32(mpi_t rop, uint32_t op)
{
	size_t nmemb = ceil_div(32, 31);

	mpi_enlarge(rop, nmemb);

	for (size_t n = 0; n < nmemb; ++n) {
		rop->data[n] = op & 0x7fffffff;
		op >>= 31;
	}

	for (size_t n = nmemb; n < rop->nmemb; ++n) {
		rop->data[n] = 0;
	}
}

uint64_t mpi_get_u64(const mpi_t op)
{
	size_t nmemb = op->nmemb;

	if (nmemb > ceil_div(64, 31)) {
		nmemb = ceil_div(64, 31);
	}

	uint64_t r = 0;

	assert(nmemb > 0);

	for (size_t n = nmemb - 1; n != (size_t)-1; --n) {
		r <<= 31;
		r |= op->data[n];
	}

	return r;
}

uint32_t mpi_get_u32(const mpi_t op)
{
	size_t nmemb = op->nmemb;

	if (nmemb > ceil_div(32, 31)) {
		nmemb = ceil_div(32, 31);
	}

	uint32_t r = 0;

	assert(nmemb > 0);

	for (size_t n = nmemb - 1; n != (size_t)-1; --n) {
		r <<= 31;
		r |= op->data[n];
	}

	return r;
}

void mpi_add(mpi_t rop, const mpi_t op1, const mpi_t op2)
{
	size_t nmemb = op1->nmemb > op2->nmemb ? op1->nmemb : op2->nmemb;

	mpi_enlarge(rop, nmemb);

	uint32_t c = 0;

	/* op1 + op2 */
	for (size_t n = 0; n < rop->nmemb; ++n) {
		uint32_t r1 = (n < op1->nmemb) ? op1->data[n] : 0;
		uint32_t r2 = (n < op2->nmemb) ? op2->data[n] : 0;
		rop->data[n] = r1 + r2 + c;
		c = rop->data[n] >> 31;
		rop->data[n] &= 0x7fffffff;
	}

	if (c != 0) {
		mpi_enlarge(rop, nmemb + 1);
		rop->data[nmemb] = 0 + c;
	}
}

void mpi_sub(mpi_t rop, const mpi_t op1, const mpi_t op2)
{
	size_t nmemb = op1->nmemb > op2->nmemb ? op1->nmemb : op2->nmemb;

	mpi_enlarge(rop, nmemb);

	uint32_t c = 0;

	/* op1 - op2 */
	for (size_t n = 0; n < rop->nmemb; ++n) {
		uint32_t r1 = (n < op1->nmemb) ? op1->data[n] : 0;
		uint32_t r2 = (n < op2->nmemb) ? op2->data[n] : 0;
		rop->data[n] = r1 - r2 - c;
		c = rop->data[n] >> 31;
		rop->data[n] &= 0x7fffffff;
	}

	if (c != 0) {
		abort();
	}
}

void mpi_add_u64(mpi_t rop, const mpi_t op1, uint64_t op2)
{
	size_t nmemb = op1->nmemb > ceil_div(64, 31) ? op1->nmemb : ceil_div(64, 31);

	mpi_enlarge(rop, nmemb);

	uint32_t c = 0;

	/* op1 + op2 */
	for (size_t n = 0; n < rop->nmemb; ++n) {
		uint32_t r1 = (n < op1->nmemb) ? op1->data[n] : 0;
		uint32_t r2 = op2 & 0x7fffffff;
		op2 >>= 31;
		rop->data[n] = r1 + r2 + c;
		c = rop->data[n] >> 31;
		rop->data[n] &= 0x7fffffff;
	}

	if (c != 0) {
		mpi_enlarge(rop, nmemb + 1);
		rop->data[nmemb] = 0 + c;
	}
}

void mpi_sub_u64(mpi_t rop, const mpi_t op1, uint64_t op2)
{
	size_t nmemb = op1->nmemb > ceil_div(64, 31) ? op1->nmemb : ceil_div(64, 31);

	mpi_enlarge(rop, nmemb);

	uint32_t c = 0;

	/* op1 + op2 */
	for (size_t n = 0; n < rop->nmemb; ++n) {
		uint32_t r1 = (n < op1->nmemb) ? op1->data[n] : 0;
		uint32_t r2 = op2 & 0x7fffffff;
		op2 >>= 31;
		rop->data[n] = r1 - r2 - c;
		c = rop->data[n] >> 31;
		rop->data[n] &= 0x7fffffff;
	}

	if (c != 0) {
		abort();
	}
}

void mpi_add_u32(mpi_t rop, const mpi_t op1, uint32_t op2)
{
	size_t nmemb = op1->nmemb > ceil_div(32, 31) ? op1->nmemb : ceil_div(32, 31);

	mpi_enlarge(rop, nmemb);

	uint32_t c = 0;

	/* op1 + op2 */
	for (size_t n = 0; n < rop->nmemb; ++n) {
		uint32_t r1 = (n < op1->nmemb) ? op1->data[n] : 0;
		uint32_t r2 = op2 & 0x7fffffff;
		op2 >>= 31;
		rop->data[n] = r1 + r2 + c;
		c = rop->data[n] >> 31;
		rop->data[n] &= 0x7fffffff;
	}

	if (c != 0) {
		mpi_enlarge(rop, nmemb + 1);
		rop->data[nmemb] = 0 + c;
	}
}

void mpi_sub_u32(mpi_t rop, const mpi_t op1, uint32_t op2)
{
	size_t nmemb = op1->nmemb > ceil_div(32, 31) ? op1->nmemb : ceil_div(32, 31);

	mpi_enlarge(rop, nmemb);

	uint32_t c = 0;

	/* op1 + op2 */
	for (size_t n = 0; n < rop->nmemb; ++n) {
		uint32_t r1 = (n < op1->nmemb) ? op1->data[n] : 0;
		uint32_t r2 = op2 & 0x7fffffff;
		op2 >>= 31;
		rop->data[n] = r1 - r2 - c;
		c = rop->data[n] >> 31;
		rop->data[n] &= 0x7fffffff;
	}

	if (c != 0) {
		abort();
	}
}

void mpi_mul_u32(mpi_t rop, const mpi_t op1, uint32_t op2)
{
	size_t nmemb = op1->nmemb + 1;

	mpi_enlarge(rop, nmemb);

	uint32_t c = 0;

	/* op1 * op2 */
	for (size_t n = 0; n < op1->nmemb; ++n) {
		assert(op1->data[n] <= (UINT64_MAX - c) / op2);
		uint64_t r = (uint64_t)op1->data[n] * op2 + c;
		rop->data[n] = r & 0x7fffffff;
		c = r >> 31;
	}

	while (c != 0) {
		mpi_enlarge(rop, nmemb + 1);
		rop->data[nmemb] = c & 0x7fffffff;
		c >>= 31;
	}
}

int mpi_set_str(mpi_t rop, const char *str, int base)
{
	assert(base == 10);

	size_t len = strlen(str);

	mpi_set_u32(rop, (uint32_t)0);

	for (size_t i = 0; i < len; ++i) {
		mpi_mul_u32(rop, rop, (uint32_t)10);
		assert(str[i] >= '0' && str[i] <= '9');
		mpi_add_u32(rop, rop, (uint32_t)(str[i] - '0'));
	}

	return 0;
}

void mpi_swap(mpi_t rop1, mpi_t rop2)
{
	mpi_t t;

	*t = *rop1;
	*rop1 = *rop2;
	*rop2 = *t;
}

void mpi_set(mpi_t rop, const mpi_t op)
{
	mpi_enlarge(rop, op->nmemb);

	for (size_t n = 0; n < op->nmemb; ++n) {
		rop->data[n] = op->data[n];
	}

	for (size_t n = op->nmemb; n < rop->nmemb; ++n) {
		rop->data[n] = 0;
	}
}

void mpi_mul_karatsuba(mpi_t rop, const mpi_t op1, const mpi_t op2)
{
	/* end recursion */
	if (op1->nmemb < 4 || op2->nmemb < 4) {
		mpi_mul(rop, op1, op2);
		return;
	}

	size_t nmemb = op1->nmemb > op2->nmemb ? op1->nmemb : op2->nmemb;

	size_t m = nmemb / 2;

	mpi_t x0, x1, y0, y1;

	mpi_init(x0);
	mpi_init(x1);
	mpi_init(y0);
	mpi_init(y1);

	/* x = op1 */
	/* y = op2 */
	mpi_fdiv_r_2exp(x0, op1, 31 * m);
	mpi_fdiv_q_2exp(x1, op1, 31 * m);
	mpi_fdiv_r_2exp(y0, op2, 31 * m);
	mpi_fdiv_q_2exp(y1, op2, 31 * m);

	mpi_t z0, z1, z2;

	mpi_init(z0);
	mpi_init(z1);
	mpi_init(z2);

	mpi_mul(z2, x1, y1);
	mpi_mul(z0, x0, y0);

	mpi_t w0, w1;

	mpi_init(w0);
	mpi_init(w1);

	mpi_add(w0, x0, x1);
	mpi_add(w1, y0, y1);

	mpi_mul(z1, w0, w1);
	mpi_sub(z1, z1, z2);
	mpi_sub(z1, z1, z0);

	mpi_mul_2exp(z2, z2, 31 * 2 * m);
	mpi_mul_2exp(z1, z1, 31 * m);

	mpi_add(rop, z0, z1);
	mpi_add(rop, rop, z2);

	mpi_clear(w0);
	mpi_clear(w1);

	mpi_clear(z0);
	mpi_clear(z1);
	mpi_clear(z2);

	mpi_clear(x0);
	mpi_clear(x1);
	mpi_clear(y0);
	mpi_clear(y1);

	mpi_compact(rop);
}

void mpi_mul(mpi_t rop, const mpi_t op1, const mpi_t op2)
{
	size_t nmemb = op1->nmemb + op2->nmemb;

	mpi_t tmp;

	mpi_init(tmp);

	mpi_enlarge(tmp, nmemb);

	for (size_t n = 0; n < tmp->nmemb; ++n) {
		tmp->data[n] = 0;
	}

	for (size_t n = 0; n < op1->nmemb; ++n) {
		for (size_t m = 0; m < op2->nmemb; ++m) {
			uint64_t r = (uint64_t)op1->data[n] * op2->data[m];
			uint64_t c = 0;
			for (size_t k = m + n; c != 0 || r != 0; ++k) {
				if (k >= tmp->nmemb) {
					mpi_enlarge(tmp, tmp->nmemb + 1);
				}
				tmp->data[k] += (r & 0x7fffffff) + c;
				r >>= 31;
				c = tmp->data[k] >> 31;
				tmp->data[k] &= 0x7fffffff;
			}
		}
	}

	mpi_set(rop, tmp);

	mpi_compact(rop);

	mpi_clear(tmp);
}

int mpi_cmp(const mpi_t op1, const mpi_t op2)
{
	size_t nmemb = op1->nmemb > op2->nmemb ? op1->nmemb : op2->nmemb;

	if (nmemb == 0) {
		return 0;
	}

	for (size_t n = nmemb - 1; n != (size_t)-1; --n) {
		uint32_t r1 = (n < op1->nmemb) ? op1->data[n] : 0;
		uint32_t r2 = (n < op2->nmemb) ? op2->data[n] : 0;

		if (r1 < r2) {
			return -1;
		}

		if (r1 > r2) {
			return +1;
		}
	}

	return 0;
}

int mpi_cmp_u32(const mpi_t op1, uint32_t op2)
{
	size_t nmemb = op1->nmemb;

	if (nmemb == 0) {
		return 0;
	}

	for (size_t n = nmemb - 1; n != (size_t)-1; --n) {
		uint32_t r1 = (n < op1->nmemb) ? op1->data[n] : 0;
		uint32_t r2 = 0;

		switch (n) {
			case 0: r2 = (op2      ) & 0x7fffffff; break;
			case 1: r2 = (op2 >> 31) & 0x7fffffff; break;
		}

		if (r1 < r2) {
			return -1;
		}

		if (r1 > r2) {
			return +1;
		}
	}

	return 0;
}

int mpi_odd_p(const mpi_t op)
{
	if (op->nmemb == 0) {
		return 0;
	}

	return op->data[0] & 1;
}

int mpi_even_p(const mpi_t op)
{
	if (op->nmemb == 0) {
		return 1;
	}

	return !(op->data[0] & 1);
}

uint64_t mpi_get_word_u64(const mpi_t op, size_t n)
{
	uint64_t r = 0;

	if (n + 0 < op->nmemb) {
		r |= (uint64_t)op->data[n + 0];
	}

	if (n + 1 < op->nmemb) {
		r |= (uint64_t)op->data[n + 1] << 31;
	}

	if (n + 2 < op->nmemb) {
		r |= (uint64_t)op->data[n + 2] << 62;
	}

	return r;
}

void mpi_fdiv_q_2exp(mpi_t q, const mpi_t n, mp_bitcnt_t b)
{
	size_t words = b / 31; /* shift by whole words/limbs */
	size_t bits = b % 31; /* and shift by bits */

	size_t nmemb = n->nmemb >= words ? n->nmemb - words : 0;

	mpi_t tmp;

	mpi_init(tmp);

	mpi_enlarge(tmp, nmemb);

	for (size_t i = 0; i < tmp->nmemb; ++i) {
		uint32_t r = (uint32_t)((mpi_get_word_u64(n, i + words) >> bits) & 0x7fffffff);

		tmp->data[i] = r;
	}

	mpi_set(q, tmp);

	mpi_clear(tmp);
}

void mpi_fdiv_r_2exp(mpi_t r, const mpi_t n, mp_bitcnt_t b)
{
	size_t words = b / 31; /* shift by whole words/limbs */
	size_t bits = b % 31; /* and shift by bits */

	size_t nmemb = words + 1;

	mpi_t tmp;

	mpi_init(tmp);

	mpi_enlarge(tmp, nmemb);

	for (size_t i = 0; i < words; ++i) {
		tmp->data[i] = n->data[i];
	}

	tmp->data[words] = n->data[words] & ((UINT32_C(1) << bits) - 1);

	mpi_set(r, tmp);

	mpi_clear(tmp);
}

void mpi_mul_2exp(mpi_t rop, const mpi_t op1, mp_bitcnt_t op2)
{
	size_t words = ceil_div(op2, 31); /* shift left by whole words/limbs */
	size_t bits = (31 - op2 % 31); /* and shift right by bits */

	size_t nmemb = op1->nmemb + words;

	mpi_t tmp;

	mpi_init(tmp);

	mpi_enlarge(tmp, nmemb);

	for (size_t i = 0; i < words; ++i) {
		tmp->data[i] = 0;
	}

	for (size_t i = words; i < tmp->nmemb; ++i) {
		tmp->data[i] = op1->data[i - words];
	}

	mpi_fdiv_q_2exp(tmp, tmp, bits);

	mpi_set(rop, tmp);

	mpi_clear(tmp);
}

int mpi_get_bit(const mpi_t op, mp_bitcnt_t b)
{
	size_t word = b / 31;
	size_t bit = b % 31;

	int r = 0;

	if (word < op->nmemb) {
		r = (op->data[word] >> bit) & 1;
	}

	return r;
}

mp_bitcnt_t mpi_scan1(const mpi_t op, mp_bitcnt_t starting_bit)
{
	size_t bits = 31 * op->nmemb;

	for (size_t i = starting_bit; i < bits; ++i) {
		if (mpi_get_bit(op, i) == 1) {
			return i;
		}
	}

	return (mp_bitcnt_t)-1;
}

void mpi_ui_pow_ui(mpi_t rop, uint32_t base, uint32_t exp)
{
	mpi_t b;

	mpi_init(b);

	mpi_set_u32(b, base);

	mpi_set_u32(rop, 1);

	while (exp != 0) {
		if (exp & 1) {
			mpi_mul(rop, rop, b);
		}
		mpi_mul(b, b, b);
		exp >>= 1;
	}

	mpi_clear(b);
}
