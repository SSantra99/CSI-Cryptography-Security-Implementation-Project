/*
 * Solovay–Strassen 512-bit prime generator using GMP
 * with CPU-cycle benchmarking (min / max / avg) over RUNS runs.
 *
 * Build:
 *   gcc -O2 -Wall -Wextra -std=c11 ss_512prime_bench.c -lgmp -o ss_512prime_bench
 *
 * Note: this uses x86 __rdtsc / __rdtscp and thus is for x86/x86_64 platforms.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <unistd.h>
#include <gmp.h>
#include <x86intrin.h>   // for __rdtsc and __rdtscp

/* ----------------------------- Tunable params ----------------------------- */
/* Number of Solovay–Strassen rounds (higher => smaller error prob). */
#define SS_ROUNDS 64

/* Bit size of the prime to generate */
#define PRIME_BITS 512

/* How many iterations to benchmark */
#define RUNS 10000

/* ----------------------------- Small primes ------------------------------- */
/* Quick trial division by a handful of small primes makes testing faster. */
static const unsigned small_primes[] = {
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31,
    37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
    79, 83, 89, 97, 101, 103, 107, 109, 113
};
static const size_t small_primes_count = sizeof(small_primes)/sizeof(small_primes[0]);

/* ------------------------------ rdtsc helpers ------------------------------ */
/* Use __rdtsc / __rdtscp and cpuid for serialization.
   This pattern is simpler and less error-prone than writing raw asm outputs. */
static inline uint64_t rdtsc_start(void) {
    /* serialize */
    asm volatile("cpuid" ::: "rax", "rbx", "rcx", "rdx");
    return __rdtsc();
}

static inline uint64_t rdtsc_end(void) {
    unsigned aux;
    uint64_t t = __rdtscp(&aux);       /* rdtscp is ordered wrt later instructions */
    asm volatile("cpuid" ::: "rax", "rbx", "rcx", "rdx"); /* serialize again */
    return t;
}

/* ----------------------------- RNG seeding -------------------------------- */
/* Initialize a GMP MT RNG. We mix time, pid, and clock ticks into a 64-bit seed.
   For stronger seeding, you could use /dev/urandom. */
static void init_rng(gmp_randstate_t st) {
    gmp_randinit_mt(st);

    uint64_t seed = 0;
    seed ^= (uint64_t)time(NULL);
    seed ^= (uint64_t)clock() << 32;
    seed ^= (uint64_t)getpid() * 0x9e3779b97f4a7c15ULL;

    mpz_t s;
    mpz_init(s);
    mpz_import(s, 1, 1, sizeof(seed), 0, 0, &seed);
    gmp_randseed(st, s);
    mpz_clear(s);
}

/* --------------------- 512-bit odd candidate generation ------------------- */
/* Generate a random 512-bit integer with MSB=1 (exact size) and LSB=1 (odd). */
static void random_odd_candidate_512(mpz_t n, gmp_randstate_t st) {
    mpz_urandomb(n, st, PRIME_BITS);       /* n in [0, 2^512 - 1] */
    mpz_setbit(n, PRIME_BITS - 1);         /* Ensure MSB=1 => exactly 512 bits */
    mpz_setbit(n, 0);                      /* Ensure odd */
}

/* -------------------------- Small-prime screening ------------------------- */
/* Quickly reject numbers divisible by any small prime. */
static int divisible_by_small_prime(const mpz_t n) {
    if (mpz_cmp_ui(n, 2) < 0) return 1;        /* n < 2: treat as composite for our purpose */
    if (mpz_cmp_ui(n, 2) == 0) return 0;       /* 2 is prime (won't occur because we force odd) */
    if (mpz_even_p(n)) return 1;               /* even => composite */

    for (size_t i = 0; i < small_primes_count; ++i) {
        unsigned p = small_primes[i];
        unsigned long rem = mpz_fdiv_ui(n, p);
        if (rem == 0) return 1;                /* divisible => composite */
    }
    return 0;                                   /* passed quick screen */
}

/* ---------------------- Solovay–Strassen primality ------------------------ */
/*
   Return 1 if n is a probable prime by k rounds of Solovay–Strassen, else 0.
*/
static int is_probable_prime_ss(const mpz_t n, int k, gmp_randstate_t st) {
    if (mpz_cmp_ui(n, 2) < 0) return 0;
    if (mpz_cmp_ui(n, 2) == 0) return 1;
    if (mpz_even_p(n)) return 0;

    mpz_t a, g, exp, p, jac_mod, n_minus_1, n_minus_3;
    mpz_inits(a, g, exp, p, jac_mod, n_minus_1, n_minus_3, NULL);

    mpz_sub_ui(n_minus_1, n, 1);   /* n - 1 */
    mpz_sub_ui(n_minus_3, n, 3);   /* n - 3 (upper bound for a) */

    for (int i = 0; i < k; ++i) {
        /* a ∈ [2, n-2]  -> create uniform a in [0, n-4], then add 2 */
        mpz_urandomm(a, st, n_minus_3);   /* [0, n-4] */
        mpz_add_ui(a, a, 2);              /* [2, n-2] */

        /* g = gcd(a, n) > 1 => composite */
        mpz_gcd(g, a, n);
        if (mpz_cmp_ui(g, 1) != 0) {
            mpz_clears(a,g,exp,p,jac_mod,n_minus_1,n_minus_3,NULL);
            return 0;
        }

        /* jac = Jacobi(a, n)  (GMP provides this as an int) */
        int jac = mpz_jacobi(a, n);
        if (jac == 0) {
            mpz_clears(a,g,exp,p,jac_mod,n_minus_1,n_minus_3,NULL);
            return 0;
        }

        /* p = a^((n-1)/2) mod n */
        mpz_set(exp, n_minus_1);
        mpz_fdiv_q_2exp(exp, exp, 1);   /* exp = (n-1)/2 */
        mpz_powm(p, a, exp, n);

        /* Convert jac to residue mod n: -1 -> n-1, +1 -> 1 */
        if (jac == -1) {
            mpz_sub_ui(jac_mod, n, 1);
        } else {
            mpz_set_ui(jac_mod, 1);
        }

        /* If p != jac_mod (mod n), n is composite */
        if (mpz_cmp(p, jac_mod) != 0) {
            mpz_clears(a,g,exp,p,jac_mod,n_minus_1,n_minus_3,NULL);
            return 0;
        }
    }

    mpz_clears(a,g,exp,p,jac_mod,n_minus_1,n_minus_3,NULL);
    return 1;  /* Passed all rounds => probable prime */
}

/* ------------------------- 512-bit prime generator ------------------------ */
/* Keep drawing random 512-bit odd candidates until one passes:
   1) small-prime screen
   2) SS_ROUNDS rounds of Solovay–Strassen
*/
static void generate_prime_512(mpz_t prime, gmp_randstate_t st) {
    for (;;) {
        random_odd_candidate_512(prime, st);
        if (divisible_by_small_prime(prime)) continue;
        if (is_probable_prime_ss(prime, SS_ROUNDS, st)) return;  /* found probable prime */
    }
}

/* ---------------------------------- main ---------------------------------- */
int main(void) {
    /* Initialize RNG (Mersenne Twister in GMP) */
    gmp_randstate_t st;
    init_rng(st);

    mpz_t prime;
    mpz_init(prime);

    uint64_t total = 0;
    uint64_t min_cycles = (uint64_t)-1;  /* initialize to max */
    uint64_t max_cycles = 0;

    /* Run the benchmark RUNS times */
    for (int i = 0; i < RUNS; ++i) {
        uint64_t start = rdtsc_start();
        generate_prime_512(prime, st);
        uint64_t end = rdtsc_end();

        uint64_t cycles = end - start;
        total += cycles;
        if (cycles < min_cycles) min_cycles = cycles;
        if (cycles > max_cycles) max_cycles = cycles;

        /* Optional: print progress every 1000 runs (comment out to avoid clutter) */
        if ((i + 1) % 1000 == 0) {
            fprintf(stderr, "Completed %d/%d runs\n", i + 1, RUNS);
        }
    }

    double avg = (double)total / (double)RUNS;

    printf("Ran %d prime generations (512-bit, Solovay-Strassen).\n", RUNS);
    printf("Min cycles : %llu\n", (unsigned long long)min_cycles);
    printf("Max cycles : %llu\n", (unsigned long long)max_cycles);
    printf("Avg cycles : %.2f\n", avg);

    /* Optionally show the last generated prime (hex) */
    gmp_printf("Last generated prime (hex):\n%Zx\n", prime);

    mpz_clear(prime);
    gmp_randclear(st);
    return 0;
}