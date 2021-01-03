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
		mpi_t r, s;

		/* ... */
	}

	return 0;
}
