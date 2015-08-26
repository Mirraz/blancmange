#include <gmp.h>

// op++
static inline void mpz_inc(mpz_t op) {
	mpz_t res; mpz_init(res);
	mpz_add_ui(res, op, 1);
	mpz_set(op, res);
	mpz_clear(res);
}

// op--
static inline void mpz_dec(mpz_t op) {
	mpz_t res; mpz_init(res);
	mpz_sub_ui(res, op, 1);
	mpz_set(op, res);
	mpz_clear(res);
}

// rop += op_ui
static inline void mpz_add_set_ui(mpz_t rop, unsigned long int op) {
	mpz_t res; mpz_init(res);
	mpz_add_ui(res, rop, op);
	mpz_set(rop, res);
	mpz_clear(res);
}

// rop += op
static inline void mpz_add_set(mpz_t rop, mpz_t op) {
	mpz_t res; mpz_init(res);
	mpz_add(res, rop, op);
	mpz_set(rop, res);
	mpz_clear(res);
}

// rop -= op_ui
static inline void mpz_sub_set_ui(mpz_t rop, unsigned long int op) {
	mpz_t res; mpz_init(res);
	mpz_sub_ui(res, rop, op);
	mpz_set(rop, res);
	mpz_clear(res);
}

// rop -= op
static inline void mpz_sub_set(mpz_t rop, mpz_t op) {
	mpz_t res; mpz_init(res);
	mpz_sub(res, rop, op);
	mpz_set(rop, res);
	mpz_clear(res);
}

// rop *= op_ui
static inline void mpz_mul_set_ui(mpz_t rop, unsigned long int op) {
	mpz_t res; mpz_init(res);
	mpz_mul_ui(res, rop, op);
	mpz_set(rop, res);
	mpz_clear(res);
}

// rop *= op
static inline void mpz_mul_set(mpz_t rop, mpz_t op) {
	mpz_t res; mpz_init(res);
	mpz_mul(res, rop, op);
	mpz_set(rop, res);
	mpz_clear(res);
}

