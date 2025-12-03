#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

#define PRIME_BITS 256

int inverse_demo() {
    // Big integers
    mpz_t p, q, n, e, d, g, tmp;
    // RNG state
    gmp_randstate_t state;

    mpz_inits(p, q, n, e, d, g, tmp, NULL);

    // Initialize RNG (default or MT). Then seed it.
    gmp_randinit_default(state);
    gmp_randseed_ui(state, (unsigned long) time(NULL));

    // --- Generate random 256-bit prime p ---
    mpz_urandomb(p, state, PRIME_BITS);
    mpz_setbit(p, PRIME_BITS - 1);  // force size
    mpz_setbit(p, 0);               // make odd
    mpz_nextprime(p, p);            // advance to next prime

    // --- Generate random 256-bit prime q (≠ p) ---
    do {
        mpz_urandomb(q, state, PRIME_BITS);
        mpz_setbit(q, PRIME_BITS - 1);
        mpz_setbit(q, 0);
        mpz_nextprime(q, q);
    } while (mpz_cmp(p, q) == 0);

    // n = p * q
    mpz_mul(n, p, q);

    // If you want true RSA-style inverse, use phi = (p-1)(q-1):
    // mpz_t phi;
    // mpz_init(phi);
    // mpz_sub_ui(tmp, p, 1);
    // mpz_sub_ui(d, q, 1);      // reusing d as a temp here would be confusing; keep tmp/d separate if you enable this
    // mpz_mul(phi, tmp, d);
    // ... then use 'phi' instead of 'n' in the gcd/invert below.

    // --- Choose random e with gcd(e, n) = 1 (your request) ---
    // e in [2, n-1], gcd(e, n) == 1
    do {
        mpz_urandomm(e, state, n);       // 0 <= e < n
        if (mpz_cmp_ui(e, 2) < 0)        // ensure e >= 2
            continue;
        mpz_gcd(g, e, n);
    } while (mpz_cmp_ui(g, 1) != 0);

    gmp_printf("p: %Zd\nq: %Zd\nn = p*q: %Zd\n", p, q, n);
    gmp_printf("Chosen e (gcd(e, n)=1): %Zd\n", e);

    // --- Compute d = e^{-1} mod n (as you asked) ---
    if (mpz_invert(d, e, n) == 0) {
        fprintf(stderr, "Error: inverse of e modulo n does not exist (shouldn't happen since gcd=1).\n");
        // Cleanup
        mpz_clears(p, q, n, e, d, g, tmp, NULL);
        gmp_randclear(state);
        return 1;
    }

    gmp_printf("d = e^{-1} mod n: %Zd\n", d);

    // Cleanup
    mpz_clears(p, q, n, e, d, g, tmp, NULL);
    gmp_randclear(state);
    return 0;
}

int main(void) {
    return inverse_demo();
}