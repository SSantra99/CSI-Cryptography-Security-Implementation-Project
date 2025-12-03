#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <x86intrin.h>  // for __rdtsc()

#define ROUNDS 20  // ChaCha20

#define ROTL(a,b) (((a) << (b)) | ((a) >> (32 - (b))))

#define QR(a, b, c, d)        \
    a += b; d ^= a; d = ROTL(d,16); \
    c += d; b ^= c; b = ROTL(b,12); \
    a += b; d ^= a; d = ROTL(d, 8); \
    c += d; b ^= c; b = ROTL(b, 7);

void chacha20_block(uint32_t out[16], const uint32_t in[16]) {
    int i;
    uint32_t x[16];
    memcpy(x, in, sizeof(x));

    for (i = 0; i < ROUNDS; i += 2) {
        // Odd round
        QR(x[0], x[4], x[8], x[12])
        QR(x[1], x[5], x[9], x[13])
        QR(x[2], x[6], x[10], x[14])
        QR(x[3], x[7], x[11], x[15])
        // Even round
        QR(x[0], x[5], x[10], x[15])
        QR(x[1], x[6], x[11], x[12])
        QR(x[2], x[7], x[8], x[13])
        QR(x[3], x[4], x[9], x[14])
    }

    for (i = 0; i < 16; ++i) {
        out[i] = x[i] + in[i];
    }
}

void print_block(uint32_t block[16]) {
    for (int i = 0; i < 16; i++) {
        printf("%08x ", block[i]);
        if ((i + 1) % 4 == 0) printf("\n");
    }
}

int main() {
    uint32_t input[16] = {
        0x61707865, 0x3320646e, 0x79622d32, 0x6b206574, // constant
        0x03020100, 0x07060504, 0x0b0a0908, 0x0f0e0d0c, // key part 1
        0x13121110, 0x17161514, 0x1b1a1918, 0x1f1e1d1c, // key part 2
        0x00000001, 0x00000000,                         // counter, nonce
        0x4a000000, 0x00000000
    };
    uint32_t output[16];

    unsigned long long start = __rdtsc();  // Start clock count
    chacha20_block(output, input);
    unsigned long long end = __rdtsc();    // End clock count

    printf("Output block:\n");
    print_block(output);

    printf("\nCPU clock cycles: %llu\n", end - start);
    return 0;
}