#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
//#include <stdbool.h>

/* *************************** */

static unsigned int phi[64] = {
	 0,  1,  1,  2,  2,  4,  2,  6,		//  0 -  7
	 4,  6,  4, 10,  4, 12,  6,  8,		//  8 - 15
	 8, 16,  6, 18,  8, 12, 10, 22,		// 16 - 23
	 8, 20, 12, 18, 12, 28,  8, 30,		// 24 - 31
	16, 20, 16, 24, 12, 36, 18, 24,		// 32 - 39
	16, 40, 12, 42, 20, 24, 22, 46,		// 40 - 47
	16, 42, 20, 32, 24, 52, 18, 40,		// 48 - 55
	24, 36, 28, 58, 16, 60, 30, 36		// 56 - 63
};

unsigned int period_len_f(unsigned int n) {
	assert(n > 0 && n < 64);
	while ((n & 1) == 0) n = n >> 1;
	unsigned long long int n_l = n;
	unsigned int p;
	for (p=1; p<=phi[n]/2; ++p) {
		if ((phi[n] % p) == 0) {
			if ( ( ( (1llu << p) - 1 ) % n_l ) == 0) return p;
		}
	}
	return phi[n];
}

void period_len_f_test() {
	unsigned int n;
	for (n=1; n<64; ++n) printf("%u\t%u\n", n, period_len_f(n));
}

static unsigned int period_len[64] = {
	 0,  1,  1,  2,  1,  4,  2,  3,		//  0 -  7
	 1,  6,  4, 10,  2, 12,  3,  4,		//  8 - 15
	 1,  8,  6, 18,  4,  6, 10, 11,		// 16 - 23
	 2, 20, 12, 18,  3, 28,  4,  5,		// 24 - 31
	 1, 10,  8, 12,  6, 36, 18, 12,		// 32 - 39
	 4, 20,  6, 14, 10, 12, 11, 23,		// 40 - 47
	 2, 21, 20,  8, 12, 52, 18, 20,		// 48 - 55
	 3, 18, 28, 58,  4, 60,  5,  6		// 56 - 63
};

typedef struct {
	unsigned long long int int_val;
	unsigned long long int shift_val;
	unsigned long long int period_val;
	unsigned int shift_len;
	unsigned int period_len;
} repeating_binary_struct;

repeating_binary_struct get_repeating_binary(unsigned int m, unsigned int n) {
	assert(n != 0);
	repeating_binary_struct res;
	res.shift_len = 0;
	res.shift_val = 0;
	res.period_len = 1;
	res.period_val = 0;
	
	res.int_val = m / n;
	m %= n;
	if (m == 0) return res;
	
	unsigned int k, nod = 1;
	for (k=2; k<=n; ++k) if ((m % k) == 0 && (n % k) == 0) nod = k;
	m /= nod; n /= nod;
	
	while ((n & 1) == 0) {
		n = n >> 1;
		++res.shift_len;
	}
	res.shift_val = m / n;
	m %= n;
	if (m == 0) return res;
	
	assert(n < 64);
	res.period_len = period_len[n];		// phi[n]
	res.period_val = ((1llu << res.period_len) - 1) / n * m;
	return res;
}

void get_repeating_binary_test() {
	unsigned int m, n;
	for (n=1; n<64; ++n) {
		for (m=0; m<128; ++m) {
			repeating_binary_struct frac = get_repeating_binary(m, n);
			printf("%u / %u = %llu + (%llu + %llu / (2^%u - 1)) / 2^%u\n",
					m, n, frac.int_val, frac.shift_val, frac.period_val, frac.period_len, frac.shift_len);
		}
	}
}

unsigned int bits_count(unsigned long long int x) {
	unsigned int count = 0;
	while (x > 0) {
		if (x & 1) ++count;
		x >>= 1;
	}
	return count;
}

int get_bits(unsigned long long int val, unsigned int len) {
	if (len == 0) return 0;
	int bits = 0;
	unsigned long long int mask;
	for (mask = (1llu << (len-1)); mask > 0ll; mask >>= 1) {
		bits += ((val & mask) ? -1 : 1);
	}
	return bits;
}

long long int get_coef(unsigned long long int val) {
	long long int coef = 0;
	unsigned long long int p, pow;
	for (p = val >> 1, pow = 1llu; p > 0llu; p >>= 1, pow <<= 1) {
		coef += ((val & pow) ? 1ll : -1ll) * (long long int)(p * pow);
	}
	coef += (long long int)val;
	return coef;
}

/* *************************** */

#include <math.h>

typedef long double DOUBLE_TYPE;
#define FORMAT "%.40Lf"

DOUBLE_TYPE sum(DOUBLE_TYPE a) {
	DOUBLE_TYPE s = 0.0L;
	DOUBLE_TYPE pow = 1.0L;
	while (1) {
		DOUBLE_TYPE ap = a * pow;
		if (!isnormal(ap)) break;
		DOUBLE_TYPE l = floorl(ap);
		DOUBLE_TYPE dl = ap-l, dh = l+1-ap;
		if (dl < 0.0L || dh < 0.0L) break;
		DOUBLE_TYPE d = (dl < dh ? dl : dh)/pow;
		if (!isnormal(d)) break;
		s += d;
		//printf("%Lf/%Lf" "\t" FORMAT "\t" FORMAT "\n", (dl < dh ? l : (l+1)), pow, d, s);
		pow *= 2.0L;
		if (!isnormal(pow)) break;
	}
	return s;
}

void sum_test() {
	unsigned int m, n;
	for (n=1; n<64; ++n) {
		for (m=0; m<64; ++m) {
			DOUBLE_TYPE s = sum((DOUBLE_TYPE)m/(DOUBLE_TYPE)n);
			printf("%u / %u" "\t" FORMAT "\n", m, n, s);
		}
	}
}

void run01() {
	unsigned int m, n;
	for (n=1; n<64; ++n) {
		if (n == 53 || n == 59 || n == 61) continue;
		DOUBLE_TYPE f = n * (powl(2.0, period_len[n]) - 1);
		for (m=0; m<64; ++m) {
			DOUBLE_TYPE a = (DOUBLE_TYPE)m/(DOUBLE_TYPE)n;
			DOUBLE_TYPE s = sum(a);
			DOUBLE_TYPE mul = s * f;
			DOUBLE_TYPE r = roundl(mul);
			if (fabsl(mul - r) > 0.0000006L) {
				printf("%u / %u" "\t" FORMAT "\t" FORMAT "\n", m, n, s, mul);
			}
		}
	}
}

void run02() {
	unsigned int m, n;
	for (n=1; n<64; ++n) {
		if (n == 53 || n == 59 || n == 61) continue;
		for (m=0; m<64; ++m) {
			DOUBLE_TYPE a = (DOUBLE_TYPE)m/(DOUBLE_TYPE)n;
			DOUBLE_TYPE s = sum(a);
			
			unsigned long long int f = 1;
			repeating_binary_struct frac = get_repeating_binary(m, n);
			if (frac.period_len > 1 && (frac.period_len & 1 || bits_count(frac.period_val)*2 != frac.period_len)) {
					// odd period: 7, 23, 31, 47
					// even period, count(1) != priod/2: 15=3*5, 39=3*13, 51=3*17, 55=5*11
					f *= n;
			} else {
				if (frac.shift_len > 0) f *= (1llu << frac.shift_len);
			}
			f *= ((1llu << frac.period_len) - 1);
			//printf("%u / %u" "\t" FORMAT "\t" "%llu" "\n", m, n, s, f);
			
			DOUBLE_TYPE mul = s * (DOUBLE_TYPE)f;
			DOUBLE_TYPE r = roundl(mul);
			if (fabsl(mul - r) > 0.0000001L) {
				printf("%u / %u" "\t" FORMAT "\t" FORMAT "\n", m, n, s, mul);
			}
		}
	}
}

DOUBLE_TYPE exact_algo(unsigned int m, unsigned int n) {
	repeating_binary_struct frac = get_repeating_binary(m, n);
	
	unsigned long long int p_pow = 1llu << frac.period_len;
	unsigned long long int p_pow_1 = p_pow - 1llu;
	unsigned long long int h_pow = 1llu << frac.shift_len;
	
	DOUBLE_TYPE s =
		(DOUBLE_TYPE)get_bits(frac.period_val, frac.period_len) *
			(DOUBLE_TYPE)frac.period_val *
			(DOUBLE_TYPE)p_pow / ((DOUBLE_TYPE)p_pow_1 * (DOUBLE_TYPE)p_pow_1)
		+
		(DOUBLE_TYPE)(get_coef(frac.period_val) * 2ll) / (DOUBLE_TYPE)p_pow_1;

	DOUBLE_TYPE p = s / (DOUBLE_TYPE)h_pow;
	
	DOUBLE_TYPE fr = (
			(DOUBLE_TYPE)frac.shift_val +
			(DOUBLE_TYPE)frac.period_val / (DOUBLE_TYPE)p_pow_1
		) / (DOUBLE_TYPE)h_pow;
	DOUBLE_TYPE h = (DOUBLE_TYPE)get_bits(frac.shift_val, frac.shift_len) * fr +
		(DOUBLE_TYPE)(get_coef(frac.shift_val) * 2ll) / (DOUBLE_TYPE)h_pow;
	
	return h + p;
}

void exact_algo_test() {
	unsigned int m, n;
	for (n=1; n<64; ++n) {
		for (m=0; m<64; ++m) {
			DOUBLE_TYPE s = exact_algo(m, n);
			DOUBLE_TYPE s_old = sum((DOUBLE_TYPE)m/(DOUBLE_TYPE)n);
			DOUBLE_TYPE diff = fabsl(s - s_old);
			if (diff > 0.000000000000000007L) {
				printf("%u / %u" "\t" FORMAT "\t" FORMAT "\t" FORMAT "\n", m, n, s, s_old, diff);
			}
		}
	}
}

/* *************************** */

int main() {
	assert(sizeof(unsigned int) >= 1);
	assert(sizeof(unsigned long long int) >= 8);
	
	/*unsigned int m, n;
	for (n=1; n<64; ++n) {
		for (m=0; m<64; ++m) {
			if (m == 0 || m >= n) continue;
			unsigned int k, nod = 1;
			for (k=2; k<=n; ++k) if ((m % k) == 0 && (n % k) == 0) nod = k;
			if (nod != 1) continue;
			DOUBLE_TYPE s = exact_algo(m, n);
			printf("%u / %u" "\t" FORMAT "\n", m, n, s);
		}
	}*/
	
	DOUBLE_TYPE s = sum(33.0L/67.0L);
	printf(FORMAT "\n", s);
	
	return 0;
}














/*

1/19 = 0.(000011010111100101)
	               0 	* 2^1
	           -1    	* 2^4
	           0     	* 2^5
	         +1      	* 2^6
	        +2       	* 2^7
	      +2         	* 2^9
	    +2           	* 2^11
	   +3            	* 2^12

0-	    5544333332110	(righter "0"s)
1+	    8776654322211	(righter "1"s)

+	00001101011110010
-	0000110101111001
+	000011010111100
-	00001101011110
-	0000110101111
+	000011010111
+	00001101011
+	0000110101
+	000011010
-	00001101
+	0000110
-	000011
+	00001
+	0000

   000011010111100101
+	                1
+	              1
+	           1
+	          1
+	         1
+	        1
+	      1
+	    1
+	   1

*/
