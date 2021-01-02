#include "mpi.h"
#include <stdlib.h>
#include <assert.h>

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
	size_t nmemb = op1->nmemb;

	if (nmemb < ceil_div(64, 31)) {
		nmemb = ceil_div(64, 31);
	}

	mpi_enlarge(rop, nmemb);

	uint32_t c = 0;

	/* l + r */
	for (size_t n = 0; n < nmemb; ++n) {
		uint32_t acc = c;
		acc = c;
		acc += op2 & 0x7fffffff;
		op2 >>= 31;
		if (n < op1->nmemb) {
			acc += op1->data[n];
		}
		c = acc >> 31;
		rop->data[n] = acc & 0x7fffffff;
	}

	for (size_t n = nmemb; n < rop->nmemb; ++n) {
		rop->data[n] = 0;
	}

	/* carry */
	if (c != 0) {
		mpi_enlarge(rop, nmemb + 1);
		rop->data[nmemb] = c;
	}
}

void mpi_add_u32(mpi_t rop, const mpi_t op1, uint32_t op2)
{
	size_t nmemb = op1->nmemb;

	if (nmemb < ceil_div(32, 31)) {
		nmemb = ceil_div(32, 31);
	}

	mpi_enlarge(rop, nmemb);

	uint32_t c = 0;

	/* l + r */
	for (size_t n = 0; n < nmemb; ++n) {
		uint32_t acc = c;
		acc = c;
		acc += op2 & 0x7fffffff;
		op2 >>= 31;
		if (n < op1->nmemb) {
			acc += op1->data[n];
		}
		c = acc >> 31;
		rop->data[n] = acc & 0x7fffffff;
	}

	for (size_t n = nmemb; n < rop->nmemb; ++n) {
		rop->data[n] = 0;
	}

	/* carry */
	if (c != 0) {
		mpi_enlarge(rop, nmemb + 1);
		rop->data[nmemb] = c;
	}
}
