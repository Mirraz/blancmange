#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <memory.h>
#include <assert.h>
#include <limits.h>
#include <gmp.h>
#include <math.h>		// sqrt for factorise_simple

#include "gmp_ext.c"

typedef unsigned int factors_list_size_type;

typedef struct {
	mpz_t prime;
	unsigned long int pow;
} prime_pow_struct;

typedef struct {
	prime_pow_struct *list;
	factors_list_size_type size;
} factors_list_struct;

void factors_list_clear(factors_list_struct factors) {
	factors_list_size_type i;
	for (i=0; i<factors.size; ++i) {
		mpz_clear(factors.list[i].prime);
	}
	free(factors.list);
}

void factors_list_print(factors_list_struct factors) {
	factors_list_size_type i;
	for (i=0; i<factors.size; ++i) {
		mpz_out_str(stdout, 10, factors.list[i].prime);
		printf("^");
		printf("%lu", factors.list[i].pow);
		if (i+1<factors.size) printf(" ");
	}
}

void factorise_simple_print(unsigned long int n) {
	printf("%lu", n); printf(":");
	unsigned long int n_sqrt = sqrt(n);
	unsigned long int p;
	for (p=2; p<=n_sqrt+1; ++p) {
		if ((n % p) == 0) {
			do {
				n /= p;
				printf(" %lu", p);
			} while ((n % p) == 0);
			n_sqrt = sqrt(n);
		}
	}
	if (n > 1) printf(" %lu", n);
	printf("\n");
}

#define primorials_size 16
static unsigned long int primorials[primorials_size] = {
	1ul, 2ul, 6ul, 30ul,
	210ul, 2310ul, 30030ul, 510510ul,
	9699690ul, 223092870ul, 6469693230ul, 200560490130ul,
	7420738134810ul, 304250263527210ul, 13082761331670030ul, 614889782588491410ul
};

factors_list_struct factorise_simple(unsigned long int n) {
	assert(n >= 2);
	
	assert(sizeof(unsigned long int) <= 8);			// => primorials_size == 16
	factors_list_size_type i;
	for (i=2; i<primorials_size && n>=primorials[i]; ++i);
	factors_list_struct res;
	res.size = i-1;
	res.list = malloc(sizeof(prime_pow_struct)*res.size);
	assert(res.size == 0 || res.list != NULL);

	i = 0;
	unsigned long int n_sqrt = (unsigned long int)sqrt(n) + 1;
	unsigned long int p;
	for (p=2; p<=n_sqrt; ++p) {
		if ((n % p) == 0) {
			assert(i<res.size);
			mpz_init(res.list[i].prime);
			mpz_set_ui(res.list[i].prime, p);
			unsigned long int pow = 0;
			do {
				n /= p;
				++pow;
			} while ((n % p) == 0);
			res.list[i].pow = pow;
			++i;
			n_sqrt = (unsigned long int)sqrt(n) + 1;
		}
	}
	if (n > 1) {
		assert(i<res.size);
		mpz_init(res.list[i].prime);
		mpz_set_ui(res.list[i].prime, n);
		res.list[i].pow = 1;
		++i;
	}
	res.size = i;
	return res;
}

void factorise_simple_test() {
	unsigned long int n;
	for(n = 2; n<=1024; ++n) {
		printf("%lu", n); printf("\t");
		factors_list_struct factors = factorise_simple(n);
		factors_list_print(factors);
		factors_list_clear(factors);
		printf("\n");
	}
}

factors_list_struct factorise(const mpz_t n) {
	assert(mpz_fits_ulong_p(n));
	return factorise_simple(mpz_get_ui(n));
}

factors_list_struct factors_lists_merge(const factors_list_struct factors1, const factors_list_struct factors2) {
	factors_list_size_type i = 0, j = 0, uniq_count = 0;
	while (i < factors1.size || j < factors2.size) {
		if (i == factors1.size) {
			++j;
		} else if (j == factors2.size) {
			++i;
		} else {
			int cmp = mpz_cmp(factors1.list[i].prime, factors2.list[j].prime);
			if (cmp < 0) {
				++i;
			} else if (cmp > 0) {
				++j;
			} else {
				++i;
				++j;
			}
		}
		++uniq_count;
	}
	
	factors_list_struct res;
	res.size = uniq_count;
	res.list = malloc(sizeof(prime_pow_struct)*res.size);
	assert(res.size == 0 || res.list != NULL);
	
	i = 0; j = 0;
	factors_list_size_type k = 0;
	while (i < factors1.size || j < factors2.size) {
		assert(k < res.size);
		mpz_init(res.list[k].prime);
		if (i == factors1.size) {
			mpz_set(res.list[k].prime, factors2.list[j].prime);
			res.list[k].pow = factors2.list[j].pow;
			++j;
		} else if (j == factors2.size) {
			mpz_set(res.list[k].prime, factors1.list[i].prime);
			res.list[k].pow = factors1.list[i].pow;
			++i;
		} else {
			int cmp = mpz_cmp(factors1.list[i].prime, factors2.list[j].prime);
			if (cmp < 0) {
				mpz_set(res.list[k].prime, factors1.list[i].prime);
				res.list[k].pow = factors1.list[i].pow;
				++i;
			} else if (cmp > 0) {
				mpz_set(res.list[k].prime, factors2.list[j].prime);
				res.list[k].pow = factors2.list[j].pow;
				++j;
			} else {
				mpz_set(res.list[k].prime, factors1.list[i].prime);
				res.list[k].pow = factors1.list[i].pow + factors2.list[j].pow;
				++i;
				++j;
			}
		}
		++k;
	}
	
	return res;
}

factors_list_struct phi_factors(const mpz_t n) {
	factors_list_struct n_factors = factorise(n);
//printf("n_factors = "); factors_list_print(n_factors); printf("\n");
	
	factors_list_size_type not_single_count = 0;
	factors_list_size_type i;
	for (i=0; i<n_factors.size; ++i) {
		assert(n_factors.list[i].pow > 0);
		if (n_factors.list[i].pow > 1) ++not_single_count;
	}
//printf("single_count = %u\n", n_factors.size - not_single_count);
//printf("not_single_count = %u\n", not_single_count);
	
	factors_list_struct res;
	res.size = not_single_count;
	res.list = malloc(sizeof(prime_pow_struct)*res.size);
	assert(res.size == 0 || res.list != NULL);
	
	if (not_single_count > 0) {
		factors_list_size_type j = 0;
		for (i=0; i<n_factors.size; ++i) {
			if (n_factors.list[i].pow == 1) continue;
			mpz_init(res.list[j].prime);
			mpz_set(res.list[j].prime, n_factors.list[i].prime);
			res.list[j].pow = n_factors.list[i].pow - 1;
			++j;
		}
	}
	
	if (mpz_cmp_ui(n_factors.list[0].prime, 2) == 0) i = 1;
	else i = 0;
	for (; i<n_factors.size; ++i) {
		mpz_t p_1; mpz_init(p_1);
		mpz_sub_ui(p_1, n_factors.list[i].prime, 1);
		factors_list_struct p_1_factors = factorise(p_1);
		mpz_clear(p_1);
//printf("res = "); factors_list_print(res); printf("\n");
//printf("p_1_factors = "); factors_list_print(p_1_factors); printf("\n");
		factors_list_struct res_mul = factors_lists_merge(res, p_1_factors);
//printf("res_mul = "); factors_list_print(res_mul); printf("\n");
		factors_list_clear(res);
		res.size = res_mul.size; res.list = res_mul.list;
		factors_list_clear(p_1_factors);
	}
	
	factors_list_clear(n_factors);
	return res;
}



bool check_order(const mpz_t mod, const mpz_t elem, const mpz_t ord) {
	mpz_t pow; mpz_init(pow);
	mpz_powm(pow, elem, ord, mod);
	bool res = (mpz_cmp_ui(pow, 1) == 0);
	mpz_clear(pow);
	return res;
}

void convolute_factors_list(mpz_t res, const factors_list_struct factors) {
	mpz_set_ui(res, 1);
	factors_list_size_type i;
	for (i=0; i<factors.size; ++i) {
		mpz_t pow; mpz_init(pow);
		mpz_pow_ui(pow, factors.list[i].prime, factors.list[i].pow);
		mpz_mul_set(res, pow);
		mpz_clear(pow);
	}
}

// mod >= 2, elem >= 1, gcd(elem, mod) == 1
void mul_group_mod_element_order(mpz_t ord, const mpz_t mod, const mpz_t elem) {
	factors_list_struct factors = phi_factors(mod);
//printf("phi_factors = "); factors_list_print(factors); printf("\n");
	//convolute_factors_list(ord, factors);
	//assert(check_order(mod, elem, ord));
	
	if (factors.size > 0) {
		bool check_flag;
		factors_list_size_type i;
		do {
			check_flag = false;
			i = factors.size - 1;
			while (true) {
				if (factors.list[i].pow != 0) {
					--factors.list[i].pow;
					convolute_factors_list(ord, factors);
	//printf("convoluted ord = "); mpz_out_str(stdout, 10, ord); printf("\n");
					if (check_order(mod, elem, ord)) {check_flag = true; break;}
					++factors.list[i].pow;
				}
				if (i == 0) break;
				--i;
			}
		} while (check_flag);
	}
	
	convolute_factors_list(ord, factors);
	assert(check_order(mod, elem, ord));
	factors_list_clear(factors);
}

// mod >= 2, elem >= 1, gcd(elem, mod) == 1
void mul_group_mod_element_order__naive(mpz_t ord, const mpz_t mod, const mpz_t elem) {
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
	//for (mpz_set_ui(mod, 3); mpz_cmp_ui(mod, 8192) <= 0; mpz_add_set_ui(mod, 2)) {
	for (mpz_set_ui(mod, ULONG_MAX-64); mpz_cmp_ui(mod, ULONG_MAX) <= 0; mpz_add_set_ui(mod, 2)) {
		mpz_out_str(stdout, 10, mod); printf("\t");
		mpz_t ord; mpz_init(ord);
		mul_group_mod_element_order(ord, mod, elem);
		mpz_out_str(stdout, 10, ord); printf("\n");
		/*{
			mpz_t ord_; mpz_init(ord_);
			mul_group_mod_element_order__naive(ord_, mod, elem);
			if (mpz_cmp(ord, ord_) != 0) {
				mpz_out_str(stdout, 10, ord_); printf("\n");
				assert(0);
			}
			mpz_clear(ord_);
		}*/
		mpz_clear(ord);
	}
	mpz_clear(mod);
	
	mpz_clear(elem);
}

void mul_group_mod_element_order_test2() {
	mpz_t elem; mpz_init(elem);
	mpz_t mod; mpz_init(mod);
	
	for (mpz_set_ui(elem, 2); mpz_cmp_ui(elem, 1024) <= 0; mpz_inc(elem)) {
		for (mpz_set_ui(mod, 2); mpz_cmp_ui(mod, 1024) <= 0; mpz_inc(mod)) {
			mpz_t gcd; mpz_init(gcd);
			mpz_gcd(gcd, elem, mod);
			unsigned int flag = (mpz_cmp_ui(gcd, 1) > 0);
			mpz_clear(gcd);
			if (flag) continue;
	
			mpz_out_str(stdout, 10, elem); printf("\t");
			mpz_out_str(stdout, 10, mod); printf("\t");
			mpz_t ord; mpz_init(ord);
			mul_group_mod_element_order(ord, mod, elem);
			mpz_out_str(stdout, 10, ord); printf("\n");
			/*{
				mpz_t ord_; mpz_init(ord_);
				mul_group_mod_element_order__naive(ord_, mod, elem);
				if (mpz_cmp(ord, ord_) != 0) {
					mpz_out_str(stdout, 10, ord_); printf("\n");
					assert(0);
				}
				mpz_clear(ord_);
			}*/
			mpz_clear(ord);
		}
	}
	mpz_clear(mod);
	mpz_clear(elem);
}

// -------------------

typedef bool (*iterate_divisors_cb)(const mpz_t divisor, void *data);

// ~descending
void iterate_divisors(const mpz_t n, iterate_divisors_cb cb, void *data) {
	factors_list_struct factors = factorise(n);
	unsigned long int pows[factors.size];
	factors_list_size_type i;
	for (i=0; i<factors.size; ++i) pows[i] = factors.list[i].pow;
	mpz_t divisor; mpz_init(divisor);
	
	while (true) {
		convolute_factors_list(divisor, factors);
		if (!cb(divisor, data)) break;
		
		for (i=0; i<factors.size && factors.list[i].pow == 0; ++i) factors.list[i].pow = pows[i];
		if (i == factors.size) break;
		--factors.list[i].pow;
	}
	
	mpz_clear(divisor);
	factors_list_clear(factors);
}

// ~ascending
void iterate_divisors_asc(const mpz_t n, iterate_divisors_cb cb, void *data) {
	factors_list_struct factors = factorise(n);
	unsigned long int pows[factors.size];
	factors_list_size_type i;
	for (i=0; i<factors.size; ++i) {
		pows[i] = factors.list[i].pow;
		factors.list[i].pow = 0;
	}
	mpz_t divisor; mpz_init(divisor);
	
	while (true) {
		convolute_factors_list(divisor, factors);
		if (!cb(divisor, data)) break;
		
		for (i=0; i<factors.size && factors.list[i].pow == pows[i]; ++i) factors.list[i].pow = 0;
		if (i == factors.size) break;
		++factors.list[i].pow;
	}
	
	mpz_clear(divisor);
	factors_list_clear(factors);
}

bool iterate_divisors_test_cb(const mpz_t divisor, void *data) {
	(void)data;
	mpz_out_str(stdout, 10, divisor); printf("\n");
	return true;
}

void iterate_divisors_test() {
	mpz_t n; mpz_init(n); mpz_set_ui(n, 4095);
	iterate_divisors(n, iterate_divisors_test_cb, NULL);
	mpz_clear(n);
}

