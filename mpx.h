#ifndef MPX_H
#define MPX_H

#include <stddef.h>
#include <stdint.h>

struct mpx {
	uint32_t *data;
	size_t nmemb;
};

typedef struct mpx mpx_t[1];

void mpx_init(mpx_t rop);

void mpx_clear(mpx_t rop);

void mpx_set_u64(mpx_t rop, uint64_t op);
void mpx_set_u32(mpx_t rop, uint32_t op);

uint64_t mpx_get_u64(const mpx_t op);

void mpx_add(mpx_t rop, const mpx_t op1, const mpx_t op2);

#endif
