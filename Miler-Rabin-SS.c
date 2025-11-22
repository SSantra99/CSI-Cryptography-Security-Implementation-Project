#include <stdio.h>
#include <gmp.h>
#include <time.h>
#include <stdlib.h>
#include <x86intrin.h>
#include <stdint.h>

int miller_rabin_round(const mpz_t n, const mpz_t a) {
    mpz_t nm1, d, x;
    mpz_init(nm1); mpz_sub_ui(nm1, n, 1);
    mpz_init_set(d, nm1);
    unsigned int s = 0;
    while (mpz_even_p(d)) {
        mpz_divexact_ui(d, d, 2);
        s++;
    }
    mpz_init(x);
    mpz_powm(x, a, d, n);
    if (mpz_cmp_ui(x, 1) == 0 || mpz_cmp(x, nm1) == 0) {
        mpz_clear(x); mpz_clear(d); mpz_clear(nm1);
        return 1;  // Pass (probably prime)
    }
    for (unsigned int r = 1; r < s; r++) {
        mpz_powm_ui(x, x, 2, n);
        if (mpz_cmp(x, nm1) == 0) {
            mpz_clear(x); mpz_clear(d); mpz_clear(nm1);
            return 1;
        }
    }
    mpz_clear(x); mpz_clear(d); mpz_clear(nm1);
    return 0;  // Fail (composite)
}

int is_probably_prime(const mpz_t n, unsigned int k, gmp_randstate_t state) {
    if (mpz_cmp_ui(n, 2) < 0) return 0;
    if (mpz_cmp_ui(n, 2) == 0 || mpz_cmp_ui(n, 3) == 0) return 1;
    if (mpz_even_p(n)) return 0;
    mpz_t range;
    mpz_init(range); mpz_sub_ui(range, n, 3);
    for (unsigned int i = 0; i < k; i++) {
        mpz_t a;
        mpz_init(a);
        mpz_urandomm(a, state, range);
        mpz_add_ui(a, a, 2);
        if (miller_rabin_round(n, a) == 0) {
            mpz_clear(a); mpz_clear(range);
            return 0;
        }
        mpz_clear(a);
    }
    mpz_clear(range);
    return 1;
}

void generate_random_prime(mpz_t p, unsigned int bits, gmp_randstate_t state, unsigned int k) {
    do {
        mpz_urandomb(p, state, bits);
        mpz_setbit(p, bits - 1);  // Set MSB for full bit length
        mpz_setbit(p, 0);  // Make odd
    } while (!is_probably_prime(p, k, state));
}

int main() {
    srand(time(NULL));
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, time(NULL) * rand());
    mpz_t p, q, n;
    mpz_init(p); mpz_init(q); mpz_init(n);
    unsigned int bits = 256;
    unsigned int k_prime = 40;
    generate_random_prime(p, bits, state, k_prime);
    generate_random_prime(q, bits, state, k_prime);
    mpz_mul(n, p, q);
    gmp_printf("p = %Zd\nq = %Zd\nn = %Zd\n", p, q, n);
    unsigned int num_tests = 100000;
    unsigned int false_positives = 0;
    mpz_t range;
    mpz_init(range); mpz_sub_ui(range, n, 3);
    for (unsigned int i = 0; i < num_tests; i++) {
        mpz_t a;
        mpz_init(a);
        mpz_urandomm(a, state, range);
        mpz_add_ui(a, a, 2);
        if (miller_rabin_round(n, a)) {
            false_positives++;
        }
        mpz_clear(a);
    }
    printf("False positives: %u out of %u\n", false_positives, num_tests);
    double observed_prob = (double)false_positives / num_tests;
    printf("Observed probability: %.4f\n", observed_prob);
    printf("Theoretical upper bound: 0.2500\n");
    mpz_clear(p); mpz_clear(q); mpz_clear(n); mpz_clear(range);
    gmp_randclear(state);
    return 0;
}
