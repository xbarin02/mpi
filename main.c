#include "mpi.h"
#include <stdlib.h>
#include <assert.h>
#include <inttypes.h>

uint32_t rand_u32()
{
	 return (uint64_t)rand() << 24 ^ (uint64_t)rand();
}

uint64_t rand_u64()
{
	 return (uint64_t)rand() << 48 ^ (uint64_t)rand() << 24 ^ (uint64_t)rand();
}

void get_max(mpi_t max, mpi_t n)
{
	mpi_set(max, n);

	while (0 != mpi_cmp_u32(n, 1)) {
		if (mpi_odd_p(n)) {
			mpi_mul_u32(n, n, 3);
			mpi_add_u32(n, n, 1);
		}
		mpi_fdiv_q_2exp(n, n, 1);

		if (mpi_cmp(n, max) > 0) {
			mpi_set(max, n);
		}
	}
}

int main()
{
	srand(42);

	{
		mpi_t s1, s2;

		mpi_init(s1);
		mpi_init(s2);

		uint32_t residue = rand_u32();

		mpi_set_u32(s1, UINT32_C(0));
		mpi_set_u32(s2, residue);

		for (int j = 0; j < 10000; ++j) {
			uint32_t rnd = rand_u32();

			mpi_t t;
			mpi_init(t);
			mpi_set_u32(t, rnd);
			mpi_add(s1, s1, t);
			mpi_add_u32(s2, s2, rnd);
			mpi_clear(t);
		}

		mpi_sub(s2, s2, s1);
		assert(residue == mpi_get_u64(s2));

		mpi_clear(s1);
		mpi_clear(s2);
	}

	{
		mpi_t s;

		mpi_init(s);

		mpi_set_str(s, "1234567890", 10);
		assert(UINT64_C(1234567890) == mpi_get_u64(s));

		mpi_set_str(s, "18446744073709551615", 10);
		assert(UINT64_C(18446744073709551615) == mpi_get_u64(s));

		mpi_set_str(s, "0", 10);
		assert(UINT64_C(0) == mpi_get_u64(s));

		mpi_clear(s);
	}

	{
		mpi_t r, s, t;

		mpi_init(s);
		mpi_init(r);
		mpi_init(t);

		mpi_set_str(s, "1853020188851841", 10);
		mpi_set_str(r, "22876792454961", 10);

		mpi_set_str(t, "42391158275216203514294433201", 10);

		mpi_mul(r, r, s);

		assert(0 == mpi_cmp(r, t));

		mpi_clear(s);
		mpi_clear(r);
		mpi_clear(t);
	}

	{
		mpi_t r, s, t;

		mpi_init(s);
		mpi_init(r);
		mpi_init(t);

		mpi_set_str(s, "1853020188851841", 10);
		mpi_set_str(r, "22876792454961", 10);
		mpi_set_str(t, "22876792454961", 10);

		mpi_swap(s, r);

		assert(0 == mpi_cmp(s, t));

		mpi_clear(s);
		mpi_clear(r);
		mpi_clear(t);
	}

	{
		mpi_t r, s;

		mpi_init(s);
		mpi_init(r);

		mpi_set_str(s, "42391158275216203514294433201", 10);
		mpi_fdiv_q_2exp(s, s, 23);
		mpi_set_str(r, "5053419861223245085989", 10);
		assert(0 == mpi_cmp(s, r));

		mpi_clear(s);
		mpi_clear(r);
	}

	{
		mpi_t r;
		mpi_init(r);
		mpi_set_str(r, "505341612", 10);
		assert(0 == mpi_cmp_u32(r, UINT32_C(505341612)));
		mpi_clear(r);
	}

	{
		mpi_t n, max, r;

		mpi_init(n);
		mpi_init(max);
		mpi_init(r);

		mpi_set_str(n, "274133054632352106267", 10);

		get_max(max, n);

		mpi_set_str(r, "56649062372194325899121269007146717645316", 10);

		assert(0 == mpi_cmp(max, r));

		mpi_clear(n);
		mpi_clear(max);
		mpi_clear(r);
	}

	return 0;
}
