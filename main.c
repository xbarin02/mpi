#include "mpi.h"
#include <stdlib.h>
#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
#include <string.h>

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
	mpi_ui_pow_u32(r, 3, n);
}

int collatz_max(const char *n_str, const char *max_str)
{
	mpi_t n, max, r;

	mpi_init(n);
	mpi_init(max);
	mpi_init(r);

	mpi_set_str(n, n_str, 10);

	mpi_set(max, n);

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

		mpi_set_str(s, "1144561273430837494885949696425", 10);
		mpi_fdiv_q_2exp(s, s, 100);
		mpi_set_str(r, "0", 10);
		assert(mpi_cmp(s, r) == 0);

		mpi_set_str(s, "1144561273430837494885949696425", 10);
		mpi_fdiv_q_2exp(s, s, 200);
		mpi_set_str(r, "0", 10);
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

		mpi_set_str(r, "3", 10);
		mpi_mul_2exp(r, r, 1);
		mpi_set_str(s, "6", 10);
		assert(mpi_cmp(r, s) == 0);

		mpi_set_str(r, "3", 10);
		mpi_mul_2exp(r, r, 32);
		mpi_set_str(s, "12884901888", 10);
		assert(mpi_cmp(r, s) == 0);

		mpi_set_str(r, "2147483647", 10);
		mpi_mul_2exp(r, r, 1);
		mpi_set_str(s, "4294967294", 10);
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

		mpi_set_str(s, "75212820489", 10);
		assert(mpi_scan1(s, 0) == 0);
		assert(mpi_scan1(s, 1) == 3);
		assert(mpi_scan1(s, 4) == 12);
		assert(mpi_scan1(s, 33) == 36);

		mpi_set_str(s, "72057594037927936", 10);
		assert(mpi_scan1(s, 0) == 56);
		assert(mpi_scan1(s, 10) == 56);
		assert(mpi_scan1(s, 56) == 56);

		mpi_clear(s);
	}

	printf("mpi_ui_pow_u32\n");
	{
		mpi_t s, r;
		mpi_init(s);
		mpi_init(r);

		mpi_ui_pow_u32(s, 3, 63);
		mpi_set_str(r, "1144561273430837494885949696427", 10);
		assert(mpi_cmp(s, r) == 0);

		mpi_ui_pow_u32(s, 5, 22);
		mpi_set_str(r, "2384185791015625", 10);
		assert(mpi_cmp(s, r) == 0);

		mpi_ui_pow_u32(s, 7, 31);
		mpi_set_str(r, "157775382034845806615042743", 10);
		assert(mpi_cmp(s, r) == 0);

		mpi_ui_pow_u32(s, 10, 10);
		mpi_set_str(r, "10000000000", 10);
		assert(mpi_cmp(s, r) == 0);

		mpi_ui_pow_u32(s, 51, 51);
		mpi_set_str(r, "1219211305094648479473193481872927834667576992593770717189298225284399541977208231315051", 10);
		assert(mpi_cmp(s, r) == 0);

		mpi_clear(s);
		mpi_clear(r);
	}

	printf("mpi_divisible_u32_p\n");
	{
		mpi_t s;
		mpi_init(s);

		mpi_set_str(s, "2432902008176640000", 10);
		assert(mpi_divisible_u32_p(s, 19) == 1);
		assert(mpi_divisible_u32_p(s, 20) == 1);
		assert(mpi_divisible_u32_p(s, 23) == 0);
		assert(mpi_divisible_u32_p(s, 31) == 0);

		mpi_clear(s);
	}

	printf("mpz_fdiv_u32\n");
	{
		mpi_t s;
		mpi_init(s);

		mpi_set_str(s, "123456789", 10);
		assert(mpz_fdiv_u32(s, 97) == 39);
		assert(mpz_fdiv_u32(s, 23) == 11);
		assert(mpz_fdiv_u32(s, 1000) == 789);

		mpi_clear(s);
	}

	printf("mpi_tstbit\n");
	{
		mpi_t s;
		mpi_init(s);

		mpi_set_str(s, "4886718345", 10);
		assert(mpi_tstbit(s, 0) == 1);
		assert(mpi_tstbit(s, 10) == 1);
		assert(mpi_tstbit(s, 31) == 0);
		assert(mpi_tstbit(s, 32) == 1);
		assert(mpi_tstbit(s, 33) == 0);
		assert(mpi_tstbit(s, 100) == 0);

		mpi_clear(s);
	}

	printf("mpi_setbit\n");
	{
		mpi_t s;
		mpi_init(s);

		mpi_set_str(s, "0", 10);

		mpi_setbit(s, 1);
		assert(mpi_cmp_u32(s, 2) == 0);
		mpi_setbit(s, 0);
		assert(mpi_cmp_u32(s, 3) == 0);
		mpi_setbit(s, 31);
		assert(mpi_cmp_u32(s, UINT32_C(2147483651)) == 0);

		mpi_clear(s);
	}

	printf("mpi_sizeinbase\n");
	{
		mpi_t s;
		mpi_init(s);

		mpi_set_str(s, "49152", 10);
		assert(mpi_sizeinbase(s, 2) == 16);
		mpi_set_str(s, "4295016448", 10);
		assert(mpi_sizeinbase(s, 2) == 33);

		mpi_clear(s);
	}

	printf("mpi_fdiv_qr\n");
	{
		mpi_t n, d, q, r;
		mpi_init(n);
		mpi_init(d);
		mpi_init(q);
		mpi_init(r);

		mpi_set_str(n, "549755813889", 10);
		mpi_set_str(d, "1234", 10);
		mpi_fdiv_qr(q, r, n, d);
		assert(mpi_cmp_u32(q, 445507142) == 0);
		assert(mpi_cmp_u32(r, 661) == 0);

		mpi_clear(n);
		mpi_clear(d);
		mpi_clear(q);
		mpi_clear(r);
	}

	printf("mpi_fdiv_qr_u32\n");
	{
		mpi_t n, q, r;
		mpi_init(n);
		mpi_init(q);
		mpi_init(r);

		mpi_set_str(n, "549755813889", 10);
		assert(mpi_fdiv_qr_u32(q, r, n, 1234) == 661);
		assert(mpi_cmp_u32(q, 445507142) == 0);
		assert(mpi_cmp_u32(r, 661) == 0);

		mpi_clear(n);
		mpi_clear(q);
		mpi_clear(r);
	}

	printf("gmp_sprintf\n");
	{
		char buffer[4096];

		int i = 42;

		gmp_sprintf(buffer, "i = %i\n", i);

		assert(strcmp(buffer, "i = 42\n") == 0);

		mpi_t n;
		mpi_init(n);

		mpi_set_str(n, "1234567890", 10);

		gmp_sprintf(buffer, "i = %i, n = %Zi\n", i, n);
		assert(strcmp(buffer, "i = 42, n = 1234567890\n") == 0);

		mpi_clear(n);

		gmp_sprintf(buffer, "l = %li\n", (long int)1234);
		assert(strcmp(buffer, "l = 1234\n") == 0);
	}

	printf("Collatz problem\n");
	{
		assert(collatz_max("212581558780141311", "2176718166004315761101410771585688"));
		assert(collatz_max("255875336134000063", "2415428612584587115646993931986234"));
		assert(collatz_max("484549993128097215", "4332751846533208436890106993276834"));
		assert(collatz_max("562380758422254271", "6718947974962862349115040884231904"));
		assert(collatz_max("628226286374752923", "31268160888027375005205169043314754"));
		assert(collatz_max("891563131061253151", "140246903347442029303138585287425762"));
		assert(collatz_max("1038743969413717663", "159695671984678120932209599662553676"));
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
		assert(llt(23) == 0);
		assert(llt(25) == 0);
		assert(llt(31) == 1);
		assert(llt(61) == 1);
		assert(llt(89) == 1);
		assert(llt(107) == 1);
		assert(llt(127) == 1);
		assert(llt(521) == 1);
		assert(llt(607) == 1);
		assert(llt(1279) == 1);
		assert(llt(2203) == 1);
		assert(llt(2281) == 1);
		assert(llt(3217) == 1);
		assert(llt(4253) == 1);
		assert(llt(4423) == 1);
		assert(llt(9689) == 1);
	}

	return 0;
}
