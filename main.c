#include "mpi.h"
#include <stdlib.h>
#include <assert.h>

uint32_t rand_u32()
{
	 return (uint64_t)rand() << 24 ^ (uint64_t)rand();
}

uint64_t rand_u64()
{
	 return (uint64_t)rand() << 48 ^ (uint64_t)rand() << 24 ^ (uint64_t)rand();
}

int main()
{
	srand(42);

/******************************************************************************/

	for (int i = 0; i < 10; ++i) {
		mpi_t r;
		mpi_init(r);

		mpi_t s;
		mpi_init(s);
		mpi_set_u64(s, (uint64_t)0);

		mpi_set_u64(r, rand_u64());
		mpi_set_u64(r, (uint64_t)0);

		uint64_t acc = 0;

		for (int j = 0; j < 10000; ++j) {
			uint32_t b = rand_u32();

			acc += b;
			mpi_add_u64(r, r, (uint64_t)b);

			acc += b;
			mpi_add_u32(r, r, b);

			acc += b;
			mpi_t t;
			mpi_init(t);
			mpi_set_u32(t, b);
			mpi_add(r, r, t);
			mpi_sub(s, s, t);
			mpi_sub(s, s, t);
			mpi_clear(t);
		}

		assert(acc == mpi_get_u64(r));

		mpi_clear(r);
		mpi_clear(s);
	}

/******************************************************************************/

	mpi_t r;

	mpi_init(r);

	for (int i = 0; i < 10000; ++i) {
		uint64_t a = rand_u64();

		mpi_set_u64(r, a);

		uint64_t b = mpi_get_u64(r);

		assert(a == b);
	}

	mpi_clear(r);

/******************************************************************************/

	mpi_t Rnd, Acc;

	mpi_init(Rnd);
	mpi_init(Acc);

	for (int i = 0; i < 10000; ++i) {
		uint64_t acc = 0;
		mpi_set_u32(Acc, (uint32_t)0);

		for (int j = 0; j < 10000; ++j) {
			uint32_t rnd = rand_u32();
			acc += rnd;

			mpi_set_u32(Rnd, rnd);
			mpi_add(Acc, Acc, Rnd);

			uint64_t r = mpi_get_u64(Acc);

			assert(acc == r);
		}
	}

	mpi_clear(Rnd);
	mpi_clear(Acc);

	return 0;
}
