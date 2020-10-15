#include <assert.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include "montgomery.h"

void montgomery_dump(montgomery_state_t *state) {
    printf("montgomery {\n");
    printf("  modulo      = %zu\n", state->modulus);
    printf("  reducerbits = %zu\n", state->reducerbits);
    printf("  reducer     = %zu\n", state->reducer);
    printf("  mask        = %zu\n", state->mask);
    printf("  reciprocal  = %zu\n", state->reciprocal);
    printf("  factor      = %zu\n", state->factor);
    printf("}\n");
}

uint64_t gcd(uint64_t a, uint64_t b) {
    uint64_t c;
    while ( a != 0 ) {
        c = a;
        a = b%a;
        b = c;
    }
    return b;
}

/// Fails for all moduli that are not odd.
void test_even() {
    montgomery_state_t state;

    for(size_t i = 0; i < 10000; i++) {
        assert(0 != montgomery_init(&state, 2 * i));
    }
}

/// Succeeds for all odd moduli that are at least 3.
void test_odd() {
    montgomery_state_t state;

    for(size_t i = 1; i < 10000; i++) {
        assert(0 == montgomery_init(&state, 2 * i + 1));
    }
}

/// Check invariants of the state after init.
void test_state(montgomery_state_t *state) {
    assert(state->reducer > state->modulus);
    assert(gcd(state->reducer, state->modulus) == 1);
    assert(((state->reciprocal * state->reducer) % state->modulus) == 1);
}

void test_init() {
    montgomery_state_t state;

    for(size_t i = 1; i < 10000; i++) {
        assert(0 == montgomery_init(&state, 2 * i + 1));
        test_state(&state);
    }
}

/// Make sure converting into montgomery form and out works.
void test_convert() {
    montgomery_state_t state;

    for(size_t i = 1; i < 10000; i++) {
        size_t modulus = 2 * i + 1;
        assert(0 == montgomery_init(&state, modulus));

        for(size_t i = 1; i < modulus; i++) {
            uint64_t num = montgomery_in(&state, i);
            uint64_t orig = montgomery_out(&state, num);
            assert(i == orig);
        }
    }
}

/// Test multiplications for correctness.
void test_multiply() {
    montgomery_state_t state;

    for(size_t i = 1; i < 10000; i *= 3) {
        size_t modulus = 2 * i + 1;
        assert(0 == montgomery_init(&state, modulus));

        for(size_t a = 0; a < modulus; a++) {
            for(size_t b = 0; b <= a; b++) {
                uint64_t a_conv = montgomery_in(&state, a);
                uint64_t b_conv = montgomery_in(&state, b);
                uint64_t product_conv = montgomery_multiply(&state, a_conv, b_conv);
                uint64_t product = montgomery_out(&state, product_conv);
                assert(product == ((a * b) % modulus));
            }
        }
    }
}

int main(int argc, char *argv[]) {
    test_even();
    test_odd();
    test_init();
    test_convert();
    test_multiply();

    return 0;
}
