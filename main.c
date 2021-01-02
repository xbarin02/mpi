#include "mpx.h"
#include <stdlib.h>
#include <assert.h>

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

		mpz_set_u64(r, a);

		uint64_t b = mpz_get_u64(r);

		assert(a == b);
	}

	mpx_clear(r);

	return 0;
}
