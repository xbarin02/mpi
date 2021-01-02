#ifndef MPX_H
#define MPX_H

#include <stddef.h>
#include <stdint.h>

struct mpi {
	uint32_t *data;
	size_t nmemb;
};

typedef struct mpi mpi_t[1];

void mpi_init(mpi_t rop);

void mpi_clear(mpi_t rop);

void mpi_set_u64(mpi_t rop, uint64_t op);
void mpi_set_u32(mpi_t rop, uint32_t op);

uint64_t mpi_get_u64(const mpi_t op);
uint32_t mpi_get_u32(const mpi_t op);

void mpi_add(mpi_t rop, const mpi_t op1, const mpi_t op2);
void mpi_add_u64(mpi_t rop, const mpi_t op1, uint64_t op2);
void mpi_add_u32(mpi_t rop, const mpi_t op1, uint32_t op2);

void mpi_sub(mpi_t rop, const mpi_t op1, const mpi_t op2);

#endif
