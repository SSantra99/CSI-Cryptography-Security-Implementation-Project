#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>

#define TRIALS 1000000
#define KEY_LENGTH 16
#define SBOX_LENGTH 256

// Function to get time in microseconds
double get_time_microseconds(clock_t start, clock_t end) {
    return ((double)(end - start)) / CLOCKS_PER_SEC * 1e6;
}

// RC4 Key Scheduling Algorithm (KSA)
void rc4_ksa(unsigned char *key, int key_len, double *key_time, double *swap_time) {
    unsigned char S[SBOX_LENGTH];
    int i, j = 0;

    clock_t key_start = clock();

    // Initialize S-box
    for (i = 0; i < SBOX_LENGTH; i++) {
        S[i] = i;
    }

    // Fill temporary key array
    unsigned char K[SBOX_LENGTH];
    for (i = 0; i < SBOX_LENGTH; i++) {
        K[i] = key[i % key_len];
    }

    clock_t key_end = clock();
    *key_time = get_time_microseconds(key_start, key_end);

    // Perform the actual key-scheduling (swapping)
    clock_t swap_start = clock();

    for (i = 0; i < SBOX_LENGTH; i++) {
        j = (j + S[i] + K[i]) % SBOX_LENGTH;

        // Swap S[i] and S[j]
        unsigned char temp = S[i];
        S[i] = S[j];
        S[j] = temp;
    }

    clock_t swap_end = clock();
    *swap_time = get_time_microseconds(swap_start, swap_end);
}

int main() {
    srand(time(NULL));

    // Statistics
    double min_key_time = 0, max_key_time = 0, total_key_time = 0;
    double min_swap_time = 0, max_swap_time = 0, total_swap_time = 0;

    for (int t = 0; t < TRIALS; t++) {
        unsigned char key[KEY_LENGTH];

        // Generate random key
        for (int i = 0; i < KEY_LENGTH; i++) {
            key[i] = rand() % 256;
        }

        double key_time = 0, swap_time = 0;
        rc4_ksa(key, KEY_LENGTH, &key_time, &swap_time);

        // Update stats for key time
        if (key_time < min_key_time) min_key_time = key_time;
        if (key_time > max_key_time) max_key_time = key_time;
        total_key_time += key_time;

        // Update stats for swap time
        if (swap_time < min_swap_time) min_swap_time = swap_time;
        if (swap_time > max_swap_time) max_swap_time = swap_time;
        total_swap_time += swap_time;
    }

    // Display results
    printf("RC4 KSA Performance over %d trials:\n\n", TRIALS);
    printf("Key Generation Time (microseconds):\n");
    printf("  Min:    %.3f µs\n", min_key_time);
    printf("  Max:    %.3f µs\n", max_key_time);
    printf("  Avg:    %.3f µs\n", total_key_time / TRIALS);

    printf("\nSwapping Time (microseconds):\n");
    printf("  Min:    %.3f µs\n", min_swap_time);
    printf("  Max:    %.3f µs\n", max_swap_time);
    printf("  Avg:    %.3f µs\n", total_swap_time / TRIALS);

    return 0;
}
