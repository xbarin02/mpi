#ifndef MPX_H
#define MPX_H

#include <stddef.h>
#include <stdint.h>

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
void mpi_set_u64(mpi_t rop, uint64_t op);
void mpi_set_u32(mpi_t rop, uint32_t op);

int mpi_set_str(mpi_t rop, const char *str, int base);

void mpi_swap(mpi_t rop1, mpi_t rop2);

/* Conversion Functions */

uint64_t mpi_get_u64(const mpi_t op);
uint32_t mpi_get_u32(const mpi_t op);

/* Arithmetic Functions */

void mpi_add(mpi_t rop, const mpi_t op1, const mpi_t op2);
void mpi_add_u64(mpi_t rop, const mpi_t op1, uint64_t op2);
void mpi_add_u32(mpi_t rop, const mpi_t op1, uint32_t op2);

void mpi_sub(mpi_t rop, const mpi_t op1, const mpi_t op2);
void mpi_sub_u64(mpi_t rop, const mpi_t op1, uint64_t op2);
void mpi_sub_u32(mpi_t rop, const mpi_t op1, uint32_t op2);

void mpi_mul(mpi_t rop, const mpi_t op1, const mpi_t op2);
void mpi_mul_u32(mpi_t rop, const mpi_t op1, uint32_t op2);

/* Division Functions */

void mpi_fdiv_q_2exp(mpi_t q, const mpi_t n, mp_bitcnt_t b);
void mpi_fdiv_r_2exp(mpi_t r, const mpi_t n, mp_bitcnt_t b);

/* Comparison Functions */

int mpi_cmp(const mpi_t op1, const mpi_t op2);
int mpi_cmp_u32(const mpi_t op1, uint32_t op2);

/* Miscellaneous Functions */

int mpi_odd_p(const mpi_t op);
int mpi_even_p(const mpi_t op);

#endif
