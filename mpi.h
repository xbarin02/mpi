#ifndef MPI_H
#define MPI_H

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

struct mpi {
	uint32_t *data;
	size_t nmemb;
};

typedef struct mpi mpi_t[1];

typedef size_t mp_bitcnt_t;

/* Initialization Functions */

void mpi_init(mpi_t rop);
void mpi_clear(mpi_t rop);

/* Assignment Functions */

void mpi_set(mpi_t rop, const mpi_t op);
void mpi_set_u32(mpi_t rop, uint32_t op);
void mpi_set_u64(mpi_t rop, uint64_t op);

int mpi_set_str(mpi_t rop, const char *str, int base);

void mpi_swap(mpi_t rop1, mpi_t rop2);

/* Conversion Functions */

uint32_t mpi_get_u32(const mpi_t op);
uint64_t mpi_get_u64(const mpi_t op);

/* Arithmetic Functions */

void mpi_add(mpi_t rop, const mpi_t op1, const mpi_t op2);
void mpi_add_u32(mpi_t rop, const mpi_t op1, uint32_t op2);
void mpi_add_u64(mpi_t rop, const mpi_t op1, uint64_t op2);

void mpi_sub(mpi_t rop, const mpi_t op1, const mpi_t op2);
void mpi_sub_u32(mpi_t rop, const mpi_t op1, uint32_t op2);
void mpi_sub_u64(mpi_t rop, const mpi_t op1, uint64_t op2);

void mpi_mul(mpi_t rop, const mpi_t op1, const mpi_t op2);
void mpi_mul_u32(mpi_t rop, const mpi_t op1, uint32_t op2);
void mpi_mul_2exp(mpi_t rop, const mpi_t op1, mp_bitcnt_t op2);

/* Division Functions */

void mpi_fdiv_qr(mpi_t q, mpi_t r, const mpi_t n, const mpi_t d);
uint32_t mpi_fdiv_qr_u32(mpi_t q, mpi_t r, const mpi_t n, uint32_t d);

void mpi_fdiv_q_2exp(mpi_t q, const mpi_t n, mp_bitcnt_t b);
void mpi_fdiv_r_2exp(mpi_t r, const mpi_t n, mp_bitcnt_t b);

uint32_t mpz_fdiv_u32(const mpi_t n, uint32_t d);

int mpi_divisible_u32_p(const mpi_t n, unsigned long int d);

/* Integer Exponentiation */

void mpi_ui_pow_u32(mpi_t rop, uint32_t base, uint32_t exp);

/* Number Theoretic Functions */

void mpi_gcd(mpi_t rop, const mpi_t op1, const mpi_t op2);

/* Comparison Functions */

int mpi_cmp(const mpi_t op1, const mpi_t op2);
int mpi_cmp_u32(const mpi_t op1, uint32_t op2);

/* Integer Logic and Bit Fiddling */

mp_bitcnt_t mpi_scan1(const mpi_t op, mp_bitcnt_t starting_bit);

int mpi_tstbit(const mpi_t op, mp_bitcnt_t bit_index);
void mpi_setbit(mpi_t rop, mp_bitcnt_t bit_index);

/* I/O of Integers */

size_t mpi_out_str(FILE *stream, int base, const mpi_t op);

int gmp_vfprintf(FILE *fp, const char *fmt, va_list ap);
int gmp_fprintf(FILE *fp, const char *fmt, ...);
int gmp_vsprintf(char *buf, const char *fmt, va_list ap);
int gmp_sprintf(char *buf, const char *fmt, ...);

/* Miscellaneous Functions */

int mpi_odd_p(const mpi_t op);
int mpi_even_p(const mpi_t op);

size_t mpi_sizeinbase(const mpi_t op, int base);

#endif
