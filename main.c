#include "mpi.h"
#include <stdlib.h>
#include <assert.h>
#include <inttypes.h>
#include <stdio.h>

uint32_t rand_u32()
{
	 return (uint64_t)rand() << 24 ^ (uint64_t)rand();
}

uint64_t rand_u64()
{
	 return (uint64_t)rand() << 48 ^ (uint64_t)rand() << 24 ^ (uint64_t)rand();
}

mp_bitcnt_t mpi_ctz(const mpi_t n)
{
	return mpi_scan1(n, 0);
}

void mpi_pow3(mpi_t r, uint32_t n)
{
	mpi_ui_pow_ui(r, 3, n);
}

int collatz_max(const char *n_str, const char *max_str)
{
	mpi_t n, max, r;

	mpi_init(n);
	mpi_init(max);
	mpi_init(r);

	mpi_set_str(n, n_str, 10);

	mpi_set(max, n);

#if 0
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
#else
	while (0 != mpi_cmp_u32(n, 1)) {
		mp_bitcnt_t alpha, beta;

		mpi_add_u32(n, n, 1); /* n++ */

		alpha = mpi_ctz(n);

		mpi_fdiv_q_2exp(n, n, alpha); /* n >>= alpha */

		if (alpha > UINT32_MAX) {
			alpha = UINT32_MAX;
		}

		mpi_t a;
		mpi_init(a);
		mpi_pow3(a, (uint32_t)alpha); /* n *= lut[alpha] */
		mpi_mul(n, n, a);
		mpi_clear(a);

		mpi_sub_u32(n, n, 1);

		if (mpi_cmp(n, max) > 0) {
			mpi_set(max, n);
		}

		beta = mpi_ctz(n);

		mpi_fdiv_q_2exp(n, n, beta); /* n >>= ctz(n) */
	}
#endif

	mpi_set_str(r, max_str, 10);

	int ret = mpi_cmp(max, r) == 0;

	mpi_clear(n);
	mpi_clear(max);
	mpi_clear(r);

	return ret;
}

int llt(mp_bitcnt_t p)
{
	mpi_t s;
	mpi_t m;

	mpi_init(s);
	mpi_init(m);

	mpi_set_u32(m, 1);
	mpi_mul_2exp(m, m, p);
	mpi_sub_u32(m, m, 1);

	mpi_set_u32(s, 4);

	for (size_t i = 0; i < p - 2; ++i) {
		mpi_mul(s, s, s); /* s = s^2 */
		mpi_add(s, s, m); /* s = s + m */
		mpi_sub_u32(s, s, 2); /* s = s - 2 */

		mpi_t q;

		mpi_init(q);

		mpi_fdiv_q_2exp(q, s, p);
		mpi_fdiv_r_2exp(s, s, p);
		assert(mpi_cmp(s, m) <= 0 );
		mpi_add(s, s, q);

		mpi_clear(q);

		while (mpi_cmp(s, m) >= 0) {
			mpi_sub(s, s, m);
		}

		assert(mpi_cmp(s, m) <= 0 );
	}

	int ret = mpi_cmp_u32(s, 0) == 0;

	mpi_clear(s);
	mpi_clear(m);

	return ret;
}

int main()
{
	srand(42);

	printf("mpi_init, mpi_clear\n");
	{
		mpi_t r;
		mpi_init(r);
		mpi_clear(r);
	}

	printf("mpi_set_u32, mpi_get_u32\n");
	{
		mpi_t r;
		mpi_init(r);
		mpi_set_u32(r, UINT32_C(4294967295));
		assert(mpi_get_u32(r) == UINT32_C(4294967295));
		mpi_clear(r);
	}

	printf("mpi_set_u64, mpi_get_u64\n");
	{
		mpi_t r;
		mpi_init(r);
		mpi_set_u64(r, UINT64_C(0xFFFFFFFFFFFFFFFF));
		assert(mpi_get_u64(r) == UINT64_C(0xFFFFFFFFFFFFFFFF));
		mpi_clear(r);
	}

	printf("mpi_set_str\n");
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

	printf("mpi_swap\n");
	{
		mpi_t r, s, t;

		mpi_init(s);
		mpi_init(r);
		mpi_init(t);

		mpi_set_str(s, "1853020188851841", 10);
		mpi_set_str(r, "22876792454961", 10);
		mpi_set_str(t, "22876792454961", 10);
		mpi_swap(s, r);
		assert(mpi_cmp(s, t) == 0);

		mpi_clear(s);
		mpi_clear(r);
		mpi_clear(t);
	}

	printf("mpi_cmp\n");
	{
		mpi_t r, s;
		mpi_init(r);
		mpi_init(s);

		mpi_set_str(r, "3433683820292512484657849089280", 10);
		mpi_set_str(s, "3433683820292512484657849089279", 10);

		assert(mpi_cmp(r, r) == 0);
		assert(mpi_cmp(r, s) > 0);
		assert(mpi_cmp(s, r) < 0);

		mpi_clear(r);
		mpi_clear(s);
	}

	printf("mpi_cmp_u32\n");
	{
		mpi_t r;
		mpi_init(r);

		mpi_set_str(r, "123456", 10);
		assert(mpi_cmp_u32(r, UINT32_C(123456)) == 0);
		assert(mpi_cmp_u32(r, UINT32_C(123455)) > 0);
		assert(mpi_cmp_u32(r, UINT32_C(123457)) < 0);

		mpi_clear(r);
	}

	printf("mpi_add_u32\n");
	{
		mpi_t r, s;
		mpi_init(r);
		mpi_init(s);

		mpi_set_str(r, "3433683820292512484657849089280", 10);
		mpi_add_u32(r, r, UINT32_C(2172748161));
		mpi_set_str(s, "3433683820292512484660021837441", 10);
		assert(mpi_cmp(r, s) == 0);

		mpi_clear(r);
		mpi_clear(s);
	}

	printf("mpi_add_u64\n");
	{
		mpi_t r, s;
		mpi_init(r);
		mpi_init(s);

		mpi_set_str(r, "3433683820292512484657849089280", 10);
		mpi_add_u64(r, r, UINT64_C(142393223512449));
		mpi_set_str(s, "3433683820292512627051072601729", 10);
		assert(mpi_cmp(r, s) == 0);

		mpi_clear(r);
		mpi_clear(s);
	}

	printf("mpi_add\n");
	{
		mpi_t r, s, t;
		mpi_init(r);
		mpi_init(s);
		mpi_init(t);

		mpi_set_str(r, "3433683820292512484657849089280", 10);
		mpi_set_str(t, "1144561273430837494885949696424", 10);
		mpi_add(r, r, t);
		mpi_set_str(s, "4578245093723349979543798785704", 10);
		assert(mpi_cmp(r, s) == 0);

		mpi_set_str(r, "3433683820292512484657849089280", 10);
		mpi_set_str(t, "42391158275216203514294433201", 10);
		mpi_add(r, r, t);
		mpi_set_str(s, "3476074978567728688172143522481", 10);
		assert(mpi_cmp(r, s) == 0);

		mpi_clear(r);
		mpi_clear(s);
		mpi_clear(t);
	}

	printf("mpi_sub_u32\n");
	{
		mpi_t r, s;
		mpi_init(r);
		mpi_init(s);

		mpi_set_str(r, "3433683820292512484657849089280", 10);
		mpi_sub_u32(r, r, 2);
		mpi_set_str(s, "3433683820292512484657849089278", 10);
		assert(mpi_cmp(r, s) == 0);

		mpi_set_str(r, "18446744073709551616", 10);
		mpi_sub_u32(r, r, 2);
		mpi_set_str(s, "18446744073709551614", 10);
		assert(mpi_cmp(r, s) == 0);

		mpi_clear(r);
		mpi_clear(s);
	}

	printf("mpi_sub\n");
	{
		mpi_t r, s, t;
		mpi_init(r);
		mpi_init(s);
		mpi_init(t);

		mpi_set_str(r, "3433683820292512484657849089280", 10);
		mpi_set_str(t, "1144561273430837494885949696424", 10);
		mpi_sub(r, r, t);
		mpi_set_str(s, "2289122546861674989771899392856", 10);
		assert(mpi_cmp(r, s) == 0);

		mpi_set_str(r, "423911582752162035142944332014", 10);
		mpi_set_str(t, "11445612734308374948859496924", 10);
		mpi_sub(r, r, t);
		mpi_set_str(s, "412465970017853660194084835090", 10);
		assert(mpi_cmp(r, s) == 0);

		mpi_clear(r);
		mpi_clear(s);
		mpi_clear(t);
	}

	printf("mpi_mul\n");
	{
		mpi_t r, s, t;

		mpi_init(s);
		mpi_init(r);
		mpi_init(t);

		mpi_set_str(s, "1853020188851841", 10);
		mpi_set_str(r, "22876792454961", 10);
		mpi_set_str(t, "42391158275216203514294433201", 10);
		mpi_mul(r, r, s);
		assert(mpi_cmp(r, t) == 0);

		mpi_set_str(s, "1797010299914431210413179829509605039731475627537851106400", 10);
		mpi_set_str(r, "42391158275216203514294433201", 10);
		mpi_set_str(t, "76177348045866392339289727720615561750424801402395196723959174586681921139518743586400", 10);
		mpi_mul(r, r, s);
		assert(mpi_cmp(r, t) == 0);

		mpi_set_str(s, "2147483648", 10);
		mpi_mul(s, s, s);
		mpi_set_str(t, "4611686018427387904", 10);
		assert(mpi_cmp(s, t) == 0);

		mpi_clear(s);
		mpi_clear(r);
		mpi_clear(t);
	}

	printf("mpi_fdiv_q_2exp\n");
	{
		mpi_t r, s;
		mpi_init(s);
		mpi_init(r);

		mpi_set_str(s, "42391158275216203514294433201", 10);
		mpi_fdiv_q_2exp(s, s, 23);
		mpi_set_str(r, "5053419861223245085989", 10);
		assert(mpi_cmp(s, r) == 0);

		mpi_set_str(s, "42391158275216203514294433201", 10);
		mpi_fdiv_q_2exp(s, s, 31);
		mpi_set_str(r, "19739921332903301117", 10);
		assert(mpi_cmp(s, r) == 0);

		mpi_set_str(s, "42391158275216203514294433201", 10);
		mpi_fdiv_q_2exp(s, s, 35);
		mpi_set_str(r, "1233745083306456319", 10);
		assert(mpi_cmp(s, r) == 0);

		mpi_set_str(s, "1797010299914431210413179829509605039731475627537851106400", 10);
		mpi_fdiv_q_2exp(s, s, 31);
		mpi_set_str(r, "836798129563420643291054214122521243864426215895", 10);
		assert(mpi_cmp(s, r) == 0);

		mpi_set_str(s, "4611686018427387903", 10);
		mpi_fdiv_q_2exp(s, s, 31);
		mpi_set_str(r, "2147483647", 10);
		assert(mpi_cmp(s, r) == 0);

		mpi_set_str(s, "9223372036854775807", 10);
		mpi_fdiv_q_2exp(s, s, 31);
		mpi_set_str(r, "4294967295", 10);
		assert(mpi_cmp(s, r) == 0);

		mpi_set_str(s, "1144561273430837494885949696425", 10);
		mpi_fdiv_q_2exp(s, s, 31);
		mpi_set_str(r, "532977875988389130162", 10);
		assert(mpi_cmp(s, r) == 0);

		mpi_clear(s);
		mpi_clear(r);
	}

	printf("mpi_fdiv_r_2exp\n");
	{
		mpi_t r, s;
		mpi_init(s);
		mpi_init(r);

		mpi_set_str(s, "42391158275216203514294433201", 10);
		mpi_fdiv_r_2exp(s, s, 23);
		mpi_set_str(r, "6419889", 10);
		assert(mpi_cmp(s, r) == 0);

		mpi_set_str(s, "42391158275216203514294433201", 10);
		mpi_fdiv_r_2exp(s, s, 31);
		mpi_set_str(r, "316798385", 10);
		assert(mpi_cmp(s, r) == 0);

		mpi_set_str(s, "42391158275216203514294433201", 10);
		mpi_fdiv_r_2exp(s, s, 35);
		mpi_set_str(r, "28234085809", 10);
		assert(mpi_cmp(s, r) == 0);

		mpi_set_str(s, "1797010299914431210413179829509605039731475627537851106400", 10);
		mpi_fdiv_r_2exp(s, s, 31);
		mpi_set_str(r, "820921440", 10);
		assert(mpi_cmp(s, r) == 0);

		mpi_set_str(s, "1144561273430837494885949696425", 10);
		mpi_fdiv_r_2exp(s, s, 31);
		mpi_set_str(r, "2111105449", 10);
		assert(mpi_cmp(s, r) == 0);

		mpi_clear(s);
		mpi_clear(r);
	}

	printf("mpi_mul_2exp\n");
	{
		mpi_t r, s;
		mpi_init(r);
		mpi_init(s);

		mpi_set_u32(r, 123456);
		mpi_mul_2exp(r, r, 89);
		mpi_set_str(s, "76415562745007953608973140099072", 10);
		assert(mpi_cmp(r, s) == 0);

		mpi_set_str(r, "532977875988389130162", 10);
		mpi_mul_2exp(r, r, 31);
		mpi_set_str(s, "1144561273430837494883838590976", 10);
		assert(mpi_cmp(r, s) == 0);

		mpi_clear(r);
		mpi_clear(s);
	}

	printf("mpi_fdiv_q_2exp, mpi_fdiv_r_2exp, mpi_mul_2exp\n");
	{
		mpi_t s, q, r;
		mpi_init(s);
		mpi_init(q);
		mpi_init(r);

		mpi_set_str(s, "1144561273430837494885949696425", 10);
		mpi_fdiv_q_2exp(q, s, 31);
		mpi_fdiv_r_2exp(r, s, 31);
		mpi_mul_2exp(q, q, 31);
		mpi_add(q, q, r);
		assert(mpi_cmp(q, s) == 0);

		mpi_clear(s);
		mpi_clear(q);
		mpi_clear(r);
	}

	printf("mpi_scan1\n");
	{
		mpi_t s;
		mpi_init(s);

		mpi_set_str(s, "173056", 10);
		assert(mpi_scan1(s, 0) == 10);
		assert(mpi_scan1(s, 11) == 13);

		mpi_set_str(s, "18446744073709551616", 10);
		assert(mpi_scan1(s, 0) == 64);

		mpi_clear(s);
	}

	printf("mpi_ui_pow_ui\n");
	{
		mpi_t s, r;
		mpi_init(s);
		mpi_init(r);

		mpi_ui_pow_ui(s, 3, 63);
		mpi_set_str(r, "1144561273430837494885949696427", 10);
		assert(mpi_cmp(s, r) == 0);

		mpi_clear(s);
		mpi_clear(r);
	}

	printf("Collatz problem\n");
	{
		assert(collatz_max("1980976057694848447", "32012333661096566765082938647132369010"));
		assert(collatz_max("35136221158664800255", "92102545196486820634779928066830214436"));
		assert(collatz_max("48503373501652785087", "296696710908147364747230298439288489642"));
		assert(collatz_max("55247846101001863167", "482192631091346876742345874501396316228"));
		assert(collatz_max("71149323674102624415", "4527691962113372170289733115168874698466"));
		assert(collatz_max("274133054632352106267", "56649062372194325899121269007146717645316"));
	}

	printf("Lucas-Lehmer test\n");
	{
		assert(llt(3) == 1);
		assert(llt(5) == 1);
		assert(llt(7) == 1);
		assert(llt(9) == 0);
		assert(llt(11) == 0);
		assert(llt(13) == 1);
		assert(llt(15) == 0);
		assert(llt(17) == 1);
		assert(llt(19) == 1);
		assert(llt(31) == 1);
		assert(llt(61) == 1);
		assert(llt(89) == 1);
		assert(llt(107) == 1);
		assert(llt(127) == 1);
	}

	return 0;
}
