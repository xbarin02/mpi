#include "mpx.h"
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

	mpx_t r;

	mpx_init(r);

	for (int i = 0; i < 10000; ++i) {
		uint64_t a = rand_u64();

		mpx_set_u64(r, a);

		uint64_t b = mpx_get_u64(r);

		assert(a == b);
	}

	mpx_clear(r);

/******************************************************************************/

	mpx_t Rnd, Acc;

	mpx_init(Rnd);
	mpx_init(Acc);

	for (int i = 0; i < 10000; ++i) {
		uint64_t acc = 0;
		mpx_set_u32(Acc, (uint32_t)0);

		for (int j = 0; j < 1000; ++j) {
			uint32_t rnd = rand_u32();
			acc += rnd;

			mpx_set_u32(Rnd, rnd);
			mpx_add(Acc, Acc, Rnd);

			uint64_t r = mpx_get_u64(Acc);

			assert(acc == r);
		}
	}

	mpx_clear(Rnd);
	mpx_clear(Acc);

	return 0;
}
