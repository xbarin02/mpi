#include "mpx.h"
#include <stdlib.h>
#include <assert.h>

void mpx_init(mpx_t rop)
{
	rop->nmemb = 0;
	rop->data = NULL;
}

void mpx_clear(mpx_t rop)
{
	free(rop->data);
}

void mpx_enlarge(mpx_t rop, size_t nmemb)
{
	if (nmemb > rop->nmemb) {
		size_t min = rop->nmemb;

		rop->nmemb = nmemb;

		rop->data = realloc(rop->data, nmemb * sizeof(uint32_t));

		if (rop->data == NULL) {
			abort();
		}

		for (size_t n = min; n < nmemb; ++n) {
			rop->data[n] = 0;
		}
	}
}

static size_t ceil_div(size_t n, size_t d)
{
	return (n + d) / d;
}

void mpx_set_u64(mpx_t rop, uint64_t op)
{
	size_t nmemb = ceil_div(64, 31);

	mpx_enlarge(rop, nmemb);

	for (size_t n = 0; n < nmemb; ++n) {
		rop->data[n] = op & 0x7fffffff;
		op >>= 31;
	}

	for (size_t n = nmemb; n < rop->nmemb; ++n) {
		rop->data[n] = 0;
	}
}

void mpx_set_u32(mpx_t rop, uint32_t op)
{
	size_t nmemb = ceil_div(32, 31);

	mpx_enlarge(rop, nmemb);

	for (size_t n = 0; n < nmemb; ++n) {
		rop->data[n] = op & 0x7fffffff;
		op >>= 31;
	}

	for (size_t n = nmemb; n < rop->nmemb; ++n) {
		rop->data[n] = 0;
	}
}

uint64_t mpx_get_u64(const mpx_t op)
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

uint32_t mpx_get_u32(const mpx_t op)
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

void mpx_add(mpx_t rop, const mpx_t op1, const mpx_t op2)
{
	mpx_t l, r;

	*l = (op1->nmemb < op2->nmemb) ? *op1 : *op2;
	*r = (op1->nmemb < op2->nmemb) ? *op2 : *op1;

	assert(l->nmemb <= r->nmemb);

	mpx_enlarge(rop, r->nmemb);

	uint64_t c = 0;

	/* l + r */
	for (size_t n = 0; n < l->nmemb; ++n) {
		rop->data[n] = l->data[n] + r->data[n] + c;
		c = rop->data[n] >> 31;
		rop->data[n] &= 0x7fffffff;
	}

	/* r */
	for (size_t n = l->nmemb; n < r->nmemb; ++n) {
		rop->data[n] = r->data[n] + c;
		c = rop->data[n] >> 31;
		rop->data[n] &= 0x7fffffff;
	}

	/* carry */
	if (c != 0) {
		mpx_enlarge(rop, r->nmemb + 1);
		rop->data[r->nmemb] = c;
	}
}
