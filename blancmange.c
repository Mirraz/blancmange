#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <gmp.h>
#include <limits.h>

#include "mul_group_mod_element_order.c"

static inline mp_bitcnt_t mpz_highest_bit(const mpz_t op) {
	assert(sizeof(mp_bitcnt_t) == sizeof(unsigned long int));
	unsigned long int bit_num_old, bit_num_new = 0;
	do {
		bit_num_old = bit_num_new;
		bit_num_new = mpz_scan1(op, bit_num_old+1);
	} while (bit_num_new < ULONG_MAX);
	return bit_num_old;
}

/*
// mod >= 2, elem >= 1, gcd(elem, mod) == 1
void mul_group_mod_element_order(mpz_t ord, const mpz_t mod, const mpz_t elem) {
	mpz_t pow; mpz_init(pow);
	for (mpz_set_ui(ord, 1); mpz_cmp(ord, mod) < 0; mpz_inc(ord)) {
		mpz_powm(pow, elem, ord, mod);
		if (mpz_cmp_ui(pow, 1) == 0) break;
	}
	mpz_clear(pow);
}

void mul_group_mod_element_order_test() {
	mpz_t elem; mpz_init(elem);
	mpz_set_ui(elem, 2);

	mpz_t mod; mpz_init(mod);
	for (mpz_set_ui(mod, 3); mpz_cmp_ui(mod, 4096) < 0; mpz_add_set_ui(mod, 2)) {
		mpz_out_str(stdout, 10, mod); printf("\t");
		mpz_t ord; mpz_init(ord);
		mul_group_mod_element_order(ord, mod, elem);
		mpz_out_str(stdout, 10, ord); printf("\n");
		mpz_clear(ord);
	}
	mpz_clear(mod);
	
	mpz_clear(elem);
}
*/

long int get_bits(const mpz_t val, unsigned long int len) {
	assert(sizeof(mp_bitcnt_t) <= sizeof(unsigned long int));
	assert(len <= LONG_MAX);
	
	unsigned long int ones = mpz_popcount(val);
	assert(ones <= len);
	return (long int)(len - ones) - (long int)ones;
}

void get_bits_test() {
	unsigned long len = 10;
	mpz_t val; mpz_init(val);
	for (mpz_set_ui(val, 1); mpz_cmp_ui(val, 1024) < 0; mpz_inc(val)) {
		mpz_out_str(stdout, 10, val); printf("\t");
		long int bits = get_bits(val, len);
		printf("%li\n", bits);
	}
	mpz_clear(val);
}

void get_coef(mpz_t coef, const mpz_t val, unsigned long int len) {
	assert(sizeof(mp_bitcnt_t) <= sizeof(unsigned long int));
	
	mpz_set(coef, val);
	
	if (len < 2) return;
	
	int prev_bit = mpz_tstbit(val, 0);
	
	mpz_t add; mpz_init(add);
	mpz_tdiv_q_2exp(add, val, 1);

	mp_bitcnt_t bit_num_last = len - 2;		// mpz_highest_bit(add)
	mp_bitcnt_t bit_num;
	
	mpz_t tmp; mpz_init(tmp);
	for (bit_num = 0; bit_num <= bit_num_last; ++bit_num) {
		if (prev_bit) mpz_add(tmp, coef, add);
		else          mpz_sub(tmp, coef, add);
		mpz_set(coef, tmp);
		prev_bit = mpz_tstbit(add, bit_num);
		mpz_clrbit(add, bit_num);
	}
	mpz_clear(tmp);
	
	mpz_clear(add);
}

void get_coef_test() {
	mpz_t val; mpz_init(val);
	for (mpz_set_ui(val, 1); mpz_cmp_ui(val, 4096) < 0; mpz_inc(val)) {
		mpz_out_str(stdout, 10, val); printf("\t");
		mpz_t coef; mpz_init(coef);
		get_coef(coef, val, 12);
		mpz_out_str(stdout, 10, coef); printf("\n");
		mpz_clear(coef);
	}
	mpz_clear(val);
}

static mpz_t cached_n;
static unsigned long int cached_period_len;

unsigned long int get_period_len(const mpz_t n) {
	if (mpz_cmp(n, cached_n) != 0) {
		mpz_t period_len_z; mpz_init(period_len_z);
		mpz_t base; mpz_init(base); mpz_set_ui(base, 2);
		mul_group_mod_element_order(period_len_z, n, base);
		assert(mpz_fits_ulong_p(period_len_z));
		cached_period_len = mpz_get_ui(period_len_z);
		mpz_clear(period_len_z);
		mpz_clear(base);
		mpz_set(cached_n, n);
	}
	return cached_period_len;
}

#define MAX_PERIOD_LEN 1048576ul

static unsigned long int dbg_period_len;
static long int dbg_bits;
static mpz_t dbg_coef;

// 0 < m < n, n - odd, gcd(m, n) == 1
void blanc_period(mpq_t res, const mpz_t m, const mpz_t n) {
	unsigned long int period_len = get_period_len(n);
dbg_period_len = period_len;
/*if ((period_len & 1) == 1) {
	mpq_set_ui(res, 0, 1);
	return;
}*/
	
	if (period_len > MAX_PERIOD_LEN) {fprintf(stderr, "period_len = %lu (too big)\n", period_len); mpq_set_ui(res, 0, 1); return;}

	mpz_t period_pow; mpz_init(period_pow);
	{
		mpz_t const1; mpz_init(const1); mpz_set_ui(const1, 1);
		mpz_mul_2exp(period_pow, const1, period_len);
		mpz_clear(const1);
	}

	
	mpz_t period_pow_1; mpz_init(period_pow_1);
	mpz_sub_ui(period_pow_1, period_pow, 1);
	
	mpz_t period_val; mpz_init(period_val);
	{
		mpz_t tmp; mpz_init(tmp);
		mpz_divexact(tmp, period_pow_1, n);
		mpz_mul(period_val, tmp, m);
		mpz_clear(tmp);
	}
	
	// ( bits(pv,pl) * 2^pl * m / n + 2 * coef(pv) ) / (2^pl - 1)
	
	// res
	{
		mpq_t div1_q; mpq_init(div1_q);
		{
			mpq_t add1_q; mpq_init(add1_q);
			{
				mpq_t dv1_q; mpq_init(dv1_q);
				{
					mpz_t dv1_z; mpz_init(dv1_z);
					{
						mpz_t mul1_z; mpz_init(mul1_z);
						{
							long int bits = get_bits(period_val, period_len);
dbg_bits = bits;
/*if (bits != 0) {
	mpz_clear(mul1_z);
	mpz_clear(dv1_z);
	mpq_clear(dv1_q);
	mpq_clear(add1_q);
	mpq_clear(div1_q);
	mpz_clear(period_pow);
	mpz_clear(period_pow_1);
	mpz_clear(period_val);
	mpq_set_ui(res, 0, 1);
	return;
}*/
							mpz_mul_si(mul1_z, period_pow, bits);
						}
						mpz_mul(dv1_z, mul1_z, m);
						mpz_clear(mul1_z);
					}
					mpq_set_z(dv1_q, dv1_z);
					mpz_clear(dv1_z);
				}
				mpq_t dv2_q; mpq_init(dv2_q);
				mpq_set_z(dv2_q, n);
				mpq_div(add1_q, dv1_q, dv2_q);
				mpq_clear(dv1_q);
				mpq_clear(dv2_q);
			}
			mpq_t add2_q; mpq_init(add2_q);
			{
				mpz_t add2_z; mpz_init(add2_z);
				{
					mpz_t coef; mpz_init(coef);
					get_coef(coef, period_val, period_len);
mpz_set(dbg_coef, coef);
					mpz_mul_2exp(add2_z, coef, 1);
					mpz_clear(coef);
				}
				mpq_set_z(add2_q, add2_z);
				mpz_clear(add2_z);
			}
			mpq_add(div1_q, add1_q, add2_q);
			mpq_clear(add1_q);
			mpq_clear(add2_q);
		}
		mpq_t div2_q; mpq_init(div2_q);
		mpq_set_z(div2_q, period_pow_1);
		mpq_div(res, div1_q, div2_q);
		mpq_clear(div1_q);
		mpq_clear(div2_q);
	}
	
	mpz_clear(period_pow);
	mpz_clear(period_pow_1);
	mpz_clear(period_val);
}

void blanc_period_test() {
	mpz_t m; mpz_init(m);
	mpz_t n; mpz_init(n);
	for (mpz_set_ui(n, 3); mpz_cmp_ui(n, 64) < 0; mpz_add_set_ui(n, 2)) {
		for (mpz_set_ui(m, 1); mpz_cmp(m, n) < 0; mpz_inc(m)) {
			mpz_t gcd; mpz_init(gcd);
			mpz_gcd(gcd, m, n);
			unsigned int flag = (mpz_cmp_ui(gcd, 1) > 0);
			mpz_clear(gcd);
			if (flag) continue;
			
			mpq_t blanc; mpq_init(blanc);
			blanc_period(blanc, m, n);
			mpz_out_str(stdout, 10, m); printf(" / ");
			mpz_out_str(stdout, 10, n); printf("\t");
			mpq_out_str(stdout, 10, blanc); printf("\n");
			mpq_clear(blanc);
		}
	}
	mpz_clear(m);
	mpz_clear(n);
}

// 0 < m < n, gcd(m, n) == 1
void blanc_prop_irred(mpq_t res, const mpz_t m, const mpz_t n) {
	assert(sizeof(mp_bitcnt_t) == sizeof(unsigned long int));
	unsigned long int shift_len = mpz_scan1(n, 0);
	
	mpz_t n_sh; mpz_init(n_sh);
	mpz_tdiv_q_2exp(n_sh, n, shift_len);
	
	mpz_t shift_val; mpz_init(shift_val);
	mpz_t m_sh; mpz_init(m_sh);
	mpz_fdiv_qr(shift_val, m_sh, m, n_sh);
	
	// bits(shv,shl) * m / n + ( 2 * coef(shv) + blanc_period(m_sh,n_sh) ) / 2^shl
	// res
	{
		mpq_t add1_q; mpq_init(add1_q);
		{
			mpq_t div1_q; mpq_init(div1_q);
			{
				long int bits = get_bits(shift_val, shift_len);
				mpz_t mul; mpz_init(mul);
				mpz_mul_si(mul, m, bits);
				mpq_set_z(div1_q, mul);
				mpz_clear(mul);
			}
			mpq_t div2_q; mpq_init(div2_q);
			mpq_set_z(div2_q, n);
			mpq_div(add1_q, div1_q, div2_q);
			mpq_clear(div1_q);
			mpq_clear(div2_q);
		}
		mpq_t add2_q; mpq_init(add2_q);
		{
			mpq_t div1_q; mpq_init(div1_q);
			{
				mpq_t ad1_q; mpq_init(ad1_q);
				{
					mpz_t ad1_z; mpz_init(ad1_z);
					{
						mpz_t coef; mpz_init(coef);
						get_coef(coef, shift_val, shift_len);
						mpz_mul_2exp(ad1_z, coef, 1);
						mpz_clear(coef);
					}
					mpq_set_z(ad1_q, ad1_z);
					mpz_clear(ad1_z);
				}
				mpq_t ad2_q; mpq_init(ad2_q);
				if (mpz_cmp_ui(n_sh, 1) > 0) {
					blanc_period(ad2_q, m_sh, n_sh);
				} else {
					assert(mpz_sgn (m_sh) == 0);
					mpq_set_ui(ad2_q, 0, 1);
				}
				mpq_add(div1_q, ad1_q, ad2_q);
				mpq_clear(ad1_q);
				mpq_clear(ad2_q);
			}
			mpq_t div2_q; mpq_init(div2_q);
			{
				mpz_t shift_pow; mpz_init(shift_pow);
				{
					mpz_t const1; mpz_init(const1); mpz_set_ui(const1, 1);
					mpz_mul_2exp(shift_pow, const1, shift_len);
					mpz_clear(const1);
				}
				mpq_set_z(div2_q, shift_pow);
				mpz_clear(shift_pow);
			}
			mpq_div(add2_q, div1_q, div2_q);
			mpq_clear(div1_q);
			mpq_clear(div2_q);
		}
		mpq_add(res, add1_q, add2_q);
		mpq_clear(add1_q);
		mpq_clear(add2_q);
	}
	
	mpz_clear(n_sh);
	mpz_clear(shift_val);
	mpz_clear(m_sh);
	
	assert(mpq_sgn(res) >= 0);
	
#ifndef NDEBUG
	mpq_t const23; mpq_init(const23);
	mpq_set_ui(const23, 2, 3);
	assert(mpq_cmp(res, const23) <= 0);
	mpq_clear(const23);
#endif
}

void blanc_prop_irred_test() {
	mpz_t m; mpz_init(m);
	mpz_t n; mpz_init(n);
	for (mpz_set_ui(n, 2); mpz_cmp_ui(n, 1024) < 0; mpz_inc(n)) {
		for (mpz_set_ui(m, 1); mpz_cmp(m, n) < 0; mpz_inc(m)) {
			mpz_t gcd; mpz_init(gcd);
			mpz_gcd(gcd, m, n);
			unsigned int flag = (mpz_cmp_ui(gcd, 1) > 0);
			mpz_clear(gcd);
			if (flag) continue;
			mpz_out_str(stdout, 10, m); printf("/");
			mpz_out_str(stdout, 10, n); printf("\t");
			
			mpq_t blanc; mpq_init(blanc);
			blanc_prop_irred(blanc, m, n);
			mpq_out_str(stdout, 10, blanc); printf("\n");
			mpq_clear(blanc);
		}
	}
	mpz_clear(m);
	mpz_clear(n);
}

void blanc_prop_irred_test_run(unsigned long int m_ui, unsigned long int n_ui) {
	assert(m_ui > 0 && n_ui >= 2 && m_ui < n_ui);
	mpz_t m; mpz_init(m); mpz_set_ui(m, m_ui);
	mpz_t n; mpz_init(n); mpz_set_ui(n, n_ui);
	{
		mpz_t gcd; mpz_init(gcd);
		mpz_gcd(gcd, m, n);
		assert(mpz_cmp_ui(gcd, 1) == 0);
		mpz_clear(gcd);
	}
	mpz_out_str(stdout, 10, m); printf("/");
	mpz_out_str(stdout, 10, n); printf("\t");
/*int n_probab_prime = mpz_probab_prime_p (n, 25);
assert(n_probab_prime == 2 || n_probab_prime == 0);
if (n_probab_prime == 2) {printf("n is prime\n"); return;}*/

	mpq_t blanc; mpq_init(blanc);
	blanc_prop_irred(blanc, m, n);
	mpq_out_str(stdout, 10, blanc); printf("\n");
	mpq_clear(blanc);

	mpz_clear(m);
	mpz_clear(n);
}

void blanc_prop_irred_dbg_test() {
	mpz_t m; mpz_init(m);
	mpz_t n; mpz_init(n);
	for (mpz_set_ui(n, 45); mpz_cmp_ui(n, 45) <= 0; mpz_add_set_ui(n, 2)) {
		for (mpz_set_ui(m, 1); mpz_cmp(m, n) < 0; mpz_inc(m)) {
			mpz_t gcd; mpz_init(gcd);
			mpz_gcd(gcd, m, n);
			unsigned int flag_gcd = (mpz_cmp_ui(gcd, 1) > 0);
			mpz_clear(gcd);
			if (flag_gcd) continue;
			
			mpq_t blanc; mpq_init(blanc);
			blanc_prop_irred(blanc, m, n);
			printf("%lu", dbg_period_len); printf("\t");
			printf("%li", dbg_bits); printf("\t");
			mpz_out_str(stdout, 10, m); printf(" / ");
			mpz_out_str(stdout, 10, n); printf("\t");
			mpz_out_str(stdout, 10, dbg_coef); printf("\t");
			mpq_out_str(stdout, 10, blanc); printf("\n");
			mpq_clear(blanc);
		}
	}
	mpz_clear(m);
	mpz_clear(n);
}

void blanc_23() {
	mpq_t const23; mpq_init(const23);
	mpq_set_ui(const23, 2, 3);

	mpz_t m; mpz_init(m);
	mpz_t n; mpz_init(n);
	for (mpz_set_ui(n, 3); mpz_cmp_ui(n, 4096+1) <= 0; mpz_add_set_ui(n, 2)) {
		for (mpz_set_ui(m, 1); mpz_cmp(m, n) < 0; mpz_inc(m)) {
			mpz_t gcd; mpz_init(gcd);
			mpz_gcd(gcd, m, n);
			unsigned int flag = (mpz_cmp_ui(gcd, 1) > 0);
			mpz_clear(gcd);
			if (flag) continue;
			
			mpq_t blanc; mpq_init(blanc);
			blanc_prop_irred(blanc, m, n);
			unsigned int flag_23 = mpq_equal(blanc, const23);
			mpq_clear(blanc);
			if (flag_23) {
				printf("%lu", dbg_period_len); printf("\t");
				printf("%li", dbg_bits); printf("\t");
				mpz_out_str(stdout, 10, dbg_coef); printf("\t");
				mpz_out_str(stdout, 10, m); printf(" / ");
				mpz_out_str(stdout, 10, n); printf("\n");
				break;
			}
		}
		//fprintf(stderr, "\t\t\t"); mpz_out_str(stderr, 10, n); fprintf(stderr, "\n");
	}
	mpz_clear(m);
	mpz_clear(n);
	
	mpq_clear(const23);
}

bool blanc_23_check(mpz_t m, const mpz_t n) {
	mpq_t const23; mpq_init(const23); mpq_set_ui(const23, 2, 3);
	mpz_t gcd; mpz_init(gcd);
	mpq_t blanc; mpq_init(blanc);
	bool m_found = false;
	mpz_t max; mpz_init(max); mpz_mul_ui(m, n, 2); mpz_cdiv_q_ui(max, m, 5);
	for (mpz_fdiv_q_ui(m, n, 3); mpz_cmp(m, max) <= 0; mpz_inc(m)) {
		mpz_gcd(gcd, m, n);
		if (mpz_cmp_ui(gcd, 1) > 0) continue;
		blanc_prop_irred(blanc, m, n);
		if (mpq_equal(blanc, const23)) {
			m_found = true;
			break;
		}
	}
	/*{
		mpz_t n3; mpz_init(n3);
		mpz_fdiv_q_ui(n3, n, 3);
		if (mpz_cmp(m, n3) <= 0) {
			printf("\n!!! m <= n/3\n");
		}
		mpz_clear(n3);
	}*/
	mpq_clear(const23);
	mpz_clear(gcd);
	mpq_clear(blanc);
	mpz_clear(max);
	return m_found;
}

void blanc_23_check_test() {
	mpz_t m; mpz_init(m);
	mpz_t n; mpz_init(n); mpz_set_ui(n, 1025);
	if (blanc_23_check(m, n)) {
		mpz_out_str(stdout, 10, m); printf(" / ");
		mpz_out_str(stdout, 10, n); printf("\n");
	}
	mpz_clear(m);
	mpz_clear(n);
}

static unsigned long int opt_p_l;

bool blanc_23_divisors_cb(const mpz_t n, void *data) {
	(void)data;
	if (mpz_cmp_ui(n, 1) == 0) return true;
	unsigned long int period_len = get_period_len(n);
	if (period_len != opt_p_l) return true;
	
	mpz_t m; mpz_init(m);
	if (blanc_23_check(m, n)) {
		printf("%lu", dbg_period_len); printf("\t");
		printf("%li", dbg_bits); printf("\t");
		mpz_out_str(stdout, 10, dbg_coef); printf("\t");
		mpz_out_str(stdout, 10, m); printf(" / ");
		mpz_out_str(stdout, 10, n); printf("\n");
	}
	mpz_clear(m);
	return true;
}

void blanc_23_divisors() {
	mpz_t const1; mpz_init(const1); mpz_set_ui(const1, 1);
	mpz_t n; mpz_init(n);
	unsigned int i;
	for (i=1; i<=17; ++i) {
		printf("4^%u\n", i);
		opt_p_l = i*2;		// [+1] opt_p_l = i*4;
		mpz_mul_2exp(n, const1, i*2);
		mpz_dec(n);			// [+1] mpz_inc(n);
		iterate_divisors(n, blanc_23_divisors_cb, NULL);
		printf("\n");
		fflush(stdout);
	}
	mpz_clear(const1);
	mpz_clear(n);
}

static mpz_t opt_period_pow_1;
static mpz_t opt_period_pow_1_3;

bool blanc_23_divs_coef_cb(const mpz_t n, void *data) {
	(void)data;
	if (mpz_cmp_ui(n, 1) == 0) return true;
	unsigned long int period_len = get_period_len(n);
	if (period_len != opt_p_l) return true;
	
	mpz_t con; mpz_init(con);
	mpz_divexact(con, opt_period_pow_1, n);
	
	mpz_t gcd; mpz_init(gcd);
	mpz_t val; mpz_init(val);
	mpz_t coef; mpz_init(coef);
	mpz_t m; mpz_init(m);
	mpz_t max; mpz_init(max);
	//if (mpz_cmp_ui(n, 17) < 0) {
		mpz_mul_ui(m, n, 2); mpz_cdiv_q_ui(max, m, 5);
	//} else {
	//	mpz_mul_ui(m, n, 6); mpz_cdiv_q_ui(max, m, 17);
	//}
	for (mpz_fdiv_q_ui(m, n, 3); mpz_cmp(m, max) <= 0; mpz_inc(m)) {
	//for (mpz_set_ui(m, 1); mpz_cmp(m, opt_period_pow_1) < 0; mpz_inc(m)) {
		mpz_gcd(gcd, m, n);
		if (mpz_cmp_ui(gcd, 1) > 0) continue;
		mpz_mul(val, con, m);
		if (get_bits(val, period_len) != 0) continue;
		get_coef(coef, val, period_len);
		if (mpz_cmp(coef, opt_period_pow_1_3) == 0) {
			long double d = (long double)mpz_get_ui(m)/(long double)mpz_get_ui(n) - 1.0L/3.0L;
			printf("%.40Lf", -log2l(d)); printf("\t");
			printf("%lu", period_len); printf("\t");
			mpz_out_str(stdout, 10, n); printf("\t");
			mpz_out_str(stdout, 10, con); printf("\t");
			mpz_out_str(stdout, 10, m); printf("\t");
			mpz_out_str(stdout, 10, val); printf("\n");
			break;
		}
	}
	
	mpz_clear(con);
	mpz_clear(gcd);
	mpz_clear(val);
	mpz_clear(coef);
	mpz_clear(m);
	mpz_clear(max);

	return true;
}

bool blanc_23_divs_coef_desc_cb(const mpz_t n, void *data) {
	(void)data;
	if (mpz_cmp_ui(n, 1) == 0) return true;
	unsigned long int period_len = get_period_len(n);
	if (period_len != opt_p_l) return true;
	
	mpz_t con; mpz_init(con);
	mpz_divexact(con, opt_period_pow_1, n);
	
	mpz_t gcd; mpz_init(gcd);
	mpz_t val; mpz_init(val);
	mpz_t coef; mpz_init(coef);
	mpz_t m; mpz_init(m);
	mpz_t min; mpz_init(min);
	mpz_fdiv_q_ui(min, n, 3);
	for (mpz_fdiv_q_ui(m, n, 2); mpz_cmp(m, min) >= 0; mpz_dec(m)) {
		mpz_gcd(gcd, m, n);
		if (mpz_cmp_ui(gcd, 1) > 0) continue;
		mpz_mul(val, con, m);
		if (get_bits(val, period_len) != 0) continue;
		get_coef(coef, val, period_len);
		if (mpz_cmp(coef, opt_period_pow_1_3) == 0) {
			long double d = 5.0L/12.0L - (long double)mpz_get_ui(m)/(long double)mpz_get_ui(n);
			printf("%.40Lf", -log2l(d)); printf("\t");
			printf("%lu", period_len); printf("\t");
			mpz_out_str(stdout, 10, n); printf("\t");
			mpz_out_str(stdout, 10, con); printf("\t");
			mpz_out_str(stdout, 10, m); printf("\t");
			mpz_out_str(stdout, 10, val); printf("\n");
			break;
		}
	}
	
	mpz_clear(con);
	mpz_clear(gcd);
	mpz_clear(val);
	mpz_clear(coef);
	mpz_clear(m);
	mpz_clear(min);

	return true;
}

//define POW4INC
void blanc_23_divs_coef() {
	mpz_t const1; mpz_init(const1); mpz_set_ui(const1, 1);
	mpz_init(opt_period_pow_1);
	mpz_init(opt_period_pow_1_3);
	mpz_t base; mpz_init(base);
	
	unsigned int i;
	for (i=1; i<=13/*17*/; ++i) {
#ifndef POW4INC
		opt_p_l = i*2;
#else
		opt_p_l = i*4;
#endif
		mpz_mul_2exp(opt_period_pow_1, const1, opt_p_l);
		mpz_dec(opt_period_pow_1);
		mpz_divexact_ui(opt_period_pow_1_3, opt_period_pow_1, 3);
#ifndef POW4INC
		mpz_set(base, opt_period_pow_1);
#else
		mpz_mul_2exp(base, const1, i*2);
		mpz_inc(base);
#endif
		iterate_divisors(base, blanc_23_divs_coef_cb, NULL);
		fflush(stdout);
	}
	
	mpz_clear(const1);
	mpz_clear(opt_period_pow_1);
	mpz_clear(opt_period_pow_1_3);
	mpz_clear(base);
}

void period_value_mix() {
	mpz_t const1; mpz_init(const1); mpz_set_ui(const1, 1);
	mpz_t pow_1; mpz_init(pow_1);
	mpz_t pow_1_3; mpz_init(pow_1_3);
	mpz_t coef; mpz_init(coef);
	mpz_t val; mpz_init(val);
	mpz_t coef_max; mpz_init(coef_max);
	
	unsigned long int len;
	for (len=2; len<=24; len+=2) {
		mpz_mul_2exp(pow_1, const1, len);
		mpz_dec(pow_1);
		// [max] mpz_set_ui(coef_max, 0);
		mpz_divexact_ui(pow_1_3, pow_1, 3);
		
		for (mpz_set_ui(val, 1); mpz_cmp(val, pow_1) < 0; mpz_inc(val)) {
			if (get_bits(val, len) != 0) continue;
			get_coef(coef, val, len);
			// [max] if (mpz_cmp(coef_max, coef) < 0) mpz_set(coef_max, coef);
			if (mpz_cmp(coef, pow_1_3) >= 0) {
				printf("%lu", len); printf("\t");
				mpz_out_str(stdout, 10, val); printf("\t");
				long double part = (long double)mpz_get_ui(pow_1)/(long double)mpz_get_ui(val);
				printf("%.40Lf", part); printf("\n");
				break;
			}
		}
		/* [max]
		printf("%lu", len); printf("\t");
		mpz_out_str(stdout, 10, coef_max); printf("\t");
		long double part = (long double)mpz_get_ui(pow_1)/(long double)mpz_get_ui(coef_max);
		printf("%.40Lf", part); printf("\n");
		*/
	}
	
	mpz_clear(val);
	mpz_clear(coef);
	mpz_clear(pow_1_3);
	mpz_clear(pow_1);
	mpz_clear(const1);
}

void period_value_max() {
	mpz_t const1; mpz_init(const1); mpz_set_ui(const1, 1);
	mpz_t pow_1; mpz_init(pow_1);
	mpz_t pow_1_3; mpz_init(pow_1_3);
	mpz_t coef; mpz_init(coef);
	mpz_t val; mpz_init(val);
	mpz_t coef_max; mpz_init(coef_max);
	
	unsigned long int len;
	for (len=2; len<=30; len+=2) {
		mpz_mul_2exp(pow_1, const1, len);
		mpz_dec(pow_1);
		// [max] mpz_set_ui(coef_max, 0);
		mpz_divexact_ui(pow_1_3, pow_1, 3);
		
		for (mpz_fdiv_q_ui(val, pow_1, 2); mpz_cmp_ui(val, 1) >= 0; mpz_dec(val)) {
			if (get_bits(val, len) != 0) continue;
			get_coef(coef, val, len);
			// [max] if (mpz_cmp(coef_max, coef) < 0) mpz_set(coef_max, coef);
			if (mpz_cmp(coef, pow_1_3) >= 0) {
				printf("%lu", len); printf("\t");
				mpz_out_str(stdout, 10, val); printf("\t");
				long double part = (long double)mpz_get_ui(val)/(long double)mpz_get_ui(pow_1);
				printf("%.40Lf", 5.0L/12.0L - part); printf("\n");
				break;
			}
		}
		/* [max]
		printf("%lu", len); printf("\t");
		mpz_out_str(stdout, 10, coef_max); printf("\t");
		long double part = (long double)mpz_get_ui(pow_1)/(long double)mpz_get_ui(coef_max);
		printf("%.40Lf", part); printf("\n");
		*/
	}
	
	mpz_clear(val);
	mpz_clear(coef);
	mpz_clear(pow_1_3);
	mpz_clear(pow_1);
	mpz_clear(const1);
}

bool blanc_23_coef_cb(const mpz_t n, void *data) {
	(void)data;
	
	if (mpz_cmp_ui(n, 1) == 0) return true;
	unsigned long int len = get_period_len(n);
	if (len != opt_p_l) return true;
	printf("# "); mpz_out_str(stdout, 10, n); printf("\n");
	
	mpz_t const1; mpz_init(const1); mpz_set_ui(const1, 1);
	mpz_t pow_1; mpz_init(pow_1);
	mpz_mul_2exp(pow_1, const1, len);
	mpz_dec(pow_1);
	mpz_t pow_1_3; mpz_init(pow_1_3);
	mpz_divexact_ui(pow_1_3, pow_1, 3);
	
	mpz_t con; mpz_init(con);
	mpz_divexact(con, pow_1, n);
	
	mpz_t gcd; mpz_init(gcd);
	mpz_t val; mpz_init(val);
	mpz_t coef; mpz_init(coef);
	mpz_t d; mpz_init(d);
	mpz_t m; mpz_init(m);
	mpz_t max; mpz_init(max);
	mpz_cdiv_q_ui(max, n, 2);
	
	for (mpz_fdiv_q_ui(m, n, 3); mpz_cmp(m, max) <= 0; mpz_inc(m)) {
	//for (mpz_set_ui(m, 1); mpz_cmp(m, opt_period_pow_1) < 0; mpz_inc(m)) {
		mpz_gcd(gcd, m, n);
		if (mpz_cmp_ui(gcd, 1) > 0) continue;
		mpz_mul(val, con, m);
		if (get_bits(val, len) != 0) continue;
		get_coef(coef, val, len);
		
		mpz_sub(d, pow_1_3, coef);
		mpz_out_str(stdout, 10, d); printf("\t");
		
		long double part = (long double)mpz_get_ui(m)/(long double)mpz_get_ui(n);
		printf("%.40Lf", part); printf("\n");
	}
	
	mpz_clear(con);
	mpz_clear(gcd);
	mpz_clear(val);
	mpz_clear(coef);
	mpz_clear(d);
	mpz_clear(m);
	mpz_clear(max);
	mpz_clear(pow_1);
	mpz_clear(pow_1_3);
	mpz_clear(const1);
	
	printf("\n");
	return true;
}

void blanc_23_coef(unsigned long int len) {
	mpz_t const1; mpz_init(const1); mpz_set_ui(const1, 1);
	
	opt_p_l = len;
	mpz_t pow_1; mpz_init(pow_1);
	mpz_mul_2exp(pow_1, const1, len);
	mpz_dec(pow_1);
	iterate_divisors(pow_1, blanc_23_coef_cb, NULL);
	
	mpz_clear(pow_1);
	mpz_clear(const1);
}

void blanc_23_all_coefs(unsigned long int len) {
	mpz_t const1; mpz_init(const1); mpz_set_ui(const1, 1);
	mpz_t pow_1; mpz_init(pow_1);
	mpz_mul_2exp(pow_1, const1, len);
	mpz_dec(pow_1);
	mpz_t pow_1_3; mpz_init(pow_1_3);
	mpz_divexact_ui(pow_1_3, pow_1, 3);
	
	mpz_t coef; mpz_init(coef);
	mpz_t d; mpz_init(d);
	mpz_t gcd; mpz_init(gcd);
	mpz_t n; mpz_init(n);
	mpz_t max; mpz_init(max); mpz_set(max, pow_1); //mpz_fdiv_q_ui(max, pow_1, 2);
	mpz_t val; mpz_init(val);
	for (mpz_set_ui(val, 1); mpz_cmp(val, max) < 0; mpz_inc(val)) {
		long int bits = get_bits(val, len);
		if (bits != 0) continue;
		get_coef(coef, val, len);
		//if (mpz_cmp(coef, pow_1_3) != 0) continue;
		/*mpz_sub(d, pow_1_3, coef);
		mpz_gcd(gcd, pow_1, val);
		mpz_divexact(n, pow_1, gcd);*/
		long double part = (long double)mpz_get_ui(val)/(long double)mpz_get_ui(pow_1);
		long double c_pt = (long double)mpz_get_ui(coef)/(long double)mpz_get_ui(pow_1);
		
		//printf("%li", bits); printf("\t");
		/*mpz_out_str(stdout, 10, val); printf("\t");
		mpz_out_str(stdout, 10, coef); printf("\t");
		mpz_out_str(stdout, 10, d); printf("\t");
		mpz_out_str(stdout, 10, n); printf("\t");
		printf("%.40Lf", part); printf("\n");*/
		
		printf("%.40Lf", part); printf("\t");
		printf("%.40Lf", c_pt); printf("\n");
	}
	
	mpz_clear(val);
	mpz_clear(max);
	mpz_clear(n);
	mpz_clear(gcd);
	mpz_clear(d);
	mpz_clear(coef);
	
	mpz_clear(pow_1_3);
	mpz_clear(pow_1);
	mpz_clear(const1);
}

void run_blanc_23_all_coefs() {
	//mpz_init(dbg_coef);
	mpz_init(cached_n); mpz_set_ui(cached_n, 0);
	
	blanc_23_all_coefs(16);
	
	//mpz_clear(dbg_coef);
	mpz_clear(cached_n);
}

void print_err_msg_and_exit(const char *msg) {
	fprintf(stderr, "%s\n", msg);
	exit(EXIT_FAILURE);
}

int main() {
	mpq_t input; mpq_init(input);
	size_t read_bytes = mpq_inp_str(input, stdin, 10);
	if (read_bytes == 0) print_err_msg_and_exit("input error");
	if (!(mpq_cmp_si(input, 0, 1) > 0 && mpq_cmp_si(input, 1, 1) < 0)) print_err_msg_and_exit("incorrect input");
	mpz_t numerator;   mpz_init(numerator);
	mpz_t denominator; mpz_init(denominator);
	mpq_get_num(numerator, input);
	mpq_get_den(denominator, input);
	assert(mpz_cmp_ui(numerator, 0) > 0 && mpz_cmp(numerator, denominator) < 0);
	mpq_t result; mpq_init(result);
	blanc_prop_irred(result, numerator, denominator);
	size_t write_bytes = mpq_out_str(stdout, 10, result);
	if (write_bytes == 0) print_err_msg_and_exit("output error");
	printf("\n");
	return 0;
}
