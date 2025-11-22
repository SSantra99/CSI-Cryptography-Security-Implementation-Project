#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <unistd.h>
#include <x86intrin.h> // for __rdtsc and __rdtscp
#include <math.h>

#define PRIME_BITS 512
#define MSG_BITS 1023
#define TRIALS 50000

// rdtsc serialization: begin uses cpuid + rdtsc, end uses rdtscp + cpuid
static inline uint64_t rdtsc_serialized_begin(void) {
    unsigned int eax, ebx, ecx, edx;
    asm volatile("cpuid"
                 : "=a"(eax), "=b"(ebx), "=c"(ecx), "=d"(edx)
                 : "a"(0)
                 : "memory");
    return __rdtsc();
}

static inline uint64_t rdtsc_serialized_end(void) {
    unsigned int aux;
    uint64_t t = __rdtscp(&aux);
    asm volatile("cpuid" : : : "rax", "rbx", "rcx", "rdx", "memory");
    return t;
}

#define INIT_STATS(minv, maxv, totalv) \
    minv = UINT64_MAX; maxv = 0; totalv = 0

#define UPDATE_STATS(val, minv, maxv, totalv) \
    if (val < minv) minv = val; \
    if (val > maxv) maxv = val; \
    totalv += (uint64_t)(val)

// convert 128-bit unsigned to long double safely
static long double u128_to_ld(__uint128_t v) {
    unsigned long long low = (unsigned long long)v;
    unsigned long long high = (unsigned long long)(v >> 64);
    return (long double)high * powl(2.0L, 64) + (long double)low;
}

int main(void) {
    gmp_randstate_t state;

    // Better entropy seed: try /dev/urandom first
    unsigned long seed = 0;
    FILE *ur = fopen("/dev/urandom", "rb");
    if (ur) {
        if (fread(&seed, sizeof(seed), 1, ur) != 1) {
            seed = (unsigned long)time(NULL) ^ (unsigned long)getpid();
        }
        fclose(ur);
    } else {
        seed = (unsigned long)time(NULL) ^ (unsigned long)getpid();
    }

    gmp_randinit_mt(state);
    gmp_randseed_ui(state, seed);

    // stats variables
    uint64_t p_min, p_max, q_min, q_max, n_min, n_max, phi_min, phi_max;
    __uint128_t p_total, q_total, n_total, phi_total;

    INIT_STATS(p_min, p_max, p_total);
    INIT_STATS(q_min, q_max, q_total);
    INIT_STATS(n_min, n_max, n_total);
    INIT_STATS(phi_min, phi_max, phi_total);

    for (int t = 0; t < TRIALS; t++) {
        mpz_t p, q, N, phi, e, d, tmp, p1, q1;
        mpz_inits(p, q, N, phi, e, d, tmp, p1, q1, NULL);

        //STEP-1
        // --- generate the first prime p ---
        uint64_t start = 0, end = 0;
        do {
            start = rdtsc_serialized_begin();
            mpz_urandomb(tmp, state, PRIME_BITS);
            mpz_setbit(tmp, PRIME_BITS - 1); // guarantee bit-length
            mpz_setbit(tmp, 0); // Setting the LSB to 1
            mpz_nextprime(p, tmp);
            end = rdtsc_serialized_end();

            // check if (p-1) divisible by e (65537); if so, reject and loop
            mpz_sub_ui(tmp, p, 1);
        } while (mpz_divisible_ui_p(tmp, 65537));
        UPDATE_STATS(end - start, p_min, p_max, p_total);

        //STEP-2
        // --- generate the second prime q ---
        uint64_t q_cycles = 0;
        do {
            start = rdtsc_serialized_begin();
            mpz_urandomb(tmp, state, PRIME_BITS);
            mpz_setbit(tmp, PRIME_BITS - 1); // guarantee bit-length
            mpz_setbit(tmp, 0); // Setting the LSB to 1
            mpz_nextprime(q, tmp);
            end = rdtsc_serialized_end();

            mpz_sub_ui(tmp, q, 1);
            q_cycles = end - start;
        } while (mpz_cmp(p, q) == 0 || mpz_divisible_ui_p(tmp, 65537));
        UPDATE_STATS(q_cycles, q_min, q_max, q_total);

        //STEP- 3
        // --- compute N ---
        start = rdtsc_serialized_begin();
        mpz_mul(N, p, q);
        end = rdtsc_serialized_end();
        UPDATE_STATS(end - start, n_min, n_max, n_total);

        //STEP-4(Computing the Euler totient function)
        // --- compute phi ---
        start = rdtsc_serialized_begin();
        mpz_sub_ui(p1, p, 1);              //p1=(p-1)
        mpz_sub_ui(q1, q, 1);              //q1=(q-1)    
        mpz_mul(phi, p1, q1);              //phi(N)=(p-1)(q-1)
        end = rdtsc_serialized_end();
        UPDATE_STATS(end - start, phi_min, phi_max, phi_total);

        mpz_set_ui(e, 65537);          //One of the commonly used key e=2^16+1
        if (mpz_invert(d, e, phi) == 0) {
            // This should no longer happen because we reject primes with (p-1) or (q-1)
            // divisible by 65537. But keep a diagnostic in case of unexpected issues.
            fprintf(stderr, "mpz_invert failed on trial %d -- skipping\n", t);
            mpz_clears(p, q, N, phi, e, d, tmp, p1, q1, NULL);
            continue;
        }

        mpz_clears(p, q, N, phi, e, d, tmp, p1, q1, NULL);
    }

    // print averages (convert 128-bit totals to long double first)
    long double p_avg_ld   = u128_to_ld(p_total) / (long double)TRIALS;
    long double q_avg_ld   = u128_to_ld(q_total) / (long double)TRIALS;
    long double n_avg_ld   = u128_to_ld(n_total) / (long double)TRIALS;
    long double phi_avg_ld = u128_to_ld(phi_total) / (long double)TRIALS;

    printf("Over %d many trials (PRIME_BITS=%d):\n\n", TRIALS, PRIME_BITS);
    printf("STEP-1 for generating the two prime numbers:\n");
    printf("First Prime p generation: min=%llu, max=%llu, avg=%.2Lf\n",
           (unsigned long long)p_min, (unsigned long long)p_max, p_avg_ld);
    printf("Second Prime q generation: min=%llu, max=%llu, avg=%.2Lf\n\n",
           (unsigned long long)q_min, (unsigned long long)q_max, q_avg_ld);
    printf("STEP-2 for computing the big number 'N' :\n");
    printf("N computation:      min=%llu, max=%llu, avg=%.2Lf\n",
           (unsigned long long)n_min, (unsigned long long)n_max, n_avg_ld);
    printf("STEP-3 for computing phi(N)");
    printf("phi(N) computation:    min=%llu, max=%llu, avg=%.2Lf\n\n",
           (unsigned long long)phi_min, (unsigned long long)phi_max, phi_avg_ld);

    // Let's observe the whole thing with a proper example
    mpz_t p, q, N, phi, e, d, tmp, p1, q1, msg, encrypted, rec;
    mpz_inits(p, q, N, phi, e, d, tmp, p1, q1, msg, encrypted, rec, NULL);

    // Generate two distinct primes p and q with guaranteed bit-length and rejecting e-factor
    do {
        mpz_urandomb(tmp, state, PRIME_BITS);
        mpz_setbit(tmp, PRIME_BITS - 1);
        mpz_setbit(tmp, 0);
        mpz_nextprime(p, tmp);
        mpz_sub_ui(tmp, p, 1);
    } while (mpz_divisible_ui_p(tmp, 65537));

    do {
        mpz_urandomb(tmp, state, PRIME_BITS);
        mpz_setbit(tmp, PRIME_BITS - 1);
        mpz_setbit(tmp, 0);
        mpz_nextprime(q, tmp);
        mpz_sub_ui(tmp, q, 1);
    } while (mpz_cmp(p, q) == 0 || mpz_divisible_ui_p(tmp, 65537));

    mpz_mul(N, p, q);
    mpz_sub_ui(p1, p, 1);
    mpz_sub_ui(q1, q, 1);
    mpz_mul(phi, p1, q1);
    mpz_set_ui(e, 65537);

    // time mpz_invert (computing d)
    uint64_t start = rdtsc_serialized_begin();
    if (mpz_invert(d, e, phi) == 0) {
        fprintf(stderr, "mpz_invert failed: e has no inverse mod phi\n");
        return 1;
    }
    uint64_t end   = rdtsc_serialized_end();
    printf("STEP-4 for finding clock cycle when computing 'd':\n");
    printf("Number of clock cycles for computing d: %llu\n\n",
           (unsigned long long)(end - start));

    // generate message (ensure MSG_BITS sized and < n)
    do {
        mpz_urandomb(msg, state, MSG_BITS);
        mpz_setbit(msg, MSG_BITS - 1);
    } while (mpz_cmp(msg, N) >= 0);

    // encryption timing (small exponent)
    printf("STEP-5 for findinf clock cycle during ENCRYPTION:\n");
    printf("Starting encryption...\n"); fflush(stdout);
    start = rdtsc_serialized_begin();
    mpz_powm(encrypted, msg, e, N);
    end   = rdtsc_serialized_end();
    printf("Here is the EnNCRYPTION . Number of clock cycles for encryption:  %llu\n",
           (unsigned long long)(end - start));

    // CRT-based decryption and timing it
    mpz_t dP, dQ, qInv, m1, m2, h;
    mpz_inits(dP, dQ, qInv, m1, m2, h, NULL);

    // dP = d mod (p-1)
    mpz_sub_ui(tmp, p, 1);
    mpz_mod(dP, d, tmp);
    // dQ = d mod (q-1)
    mpz_sub_ui(tmp, q, 1);
    mpz_mod(dQ, d, tmp);
    // qInv = q^{-1} mod p
    if (mpz_invert(qInv, q, p) == 0) {
        fprintf(stderr, "mpz_invert failed computing qInv\n");
        return 1;
    }

    printf("Starting decryption...\n"); fflush(stdout);
    start = rdtsc_serialized_begin();
    mpz_powm(m1, encrypted, dP, p);
    mpz_powm(m2, encrypted, dQ, q);

    // h = (m1 - m2) mod p
    mpz_sub(h, m1, m2);
    mpz_mod(h, h, p);

    // h = h * qInv mod p
    mpz_mul(h, h, qInv);
    mpz_mod(h, h, p);

    // rec = m2 + h * q
    mpz_mul(tmp, h, q);
    mpz_add(rec, tmp, m2);
    end = rdtsc_serialized_end();
    printf("Here is the DECRYPTION. Number of clock clock cycles for CRT decryption: %llu\n",
           (unsigned long long)(end - start));

    // Verify CRT result
    if (mpz_cmp(msg, rec) != 0) {
        fprintf(stderr, " Decryption did NOT recover the original message!\n");
    } else {
        printf("DECRYPTION has been done succesfully, matching with the original plain text message.\n");
    }

    // print example values (hex)
    gmp_printf("\nLet's see an example\n EXAMPLE:\n");
    gmp_printf("PLAIN TEXT is given by (hex):     %Zx\n\n", msg);
    gmp_printf("The encrypted message is  (hex):   %Zx\n\n", encrypted);
    gmp_printf("The decrypted message is (hex):   %Zx\n", rec);
    gmp_printf(" Observe that the decrypted message matches with the original plain text.");

    mpz_clears(p, q, N, phi, e, d, tmp, p1, q1, msg, encrypted, rec, NULL);
    mpz_clears(dP, dQ, qInv, m1, m2, h, NULL);
    gmp_randclear(state);
    return 0;
}