#include "montgomery.h"
#include <limits.h>
#include <assert.h>

#define IS_ODD(num) (num & 1)
#define IS_EVEN(num) (!(num & 1))
#define BIT_LENGTH(num) ((CHAR_BIT * sizeof(unsigned long long)) - __builtin_clzll(num))

/// Initialize montgomery state with a given modulo.
int montgomery_init(montgomery_state_t *state, uint64_t modulus) {
    if(IS_EVEN(modulus)) {
        return 1;
    }

    state->modulus = modulus;
    state->reducerbits = (BIT_LENGTH(modulus) / CHAR_BIT + 1) * CHAR_BIT;
    state->reducer = 1 << state->reducerbits;
    state->mask = state->reducer - 1;

    state->reciprocal = montgomery_multinvmod(state->reducer % modulus, modulus);
    state->factor = (state->reducer * state->reciprocal - 1) / modulus;
    state->convertedone = state->reducer % modulus;

    return 0;
}

uint64_t montgomery_in(montgomery_state_t *state, uint64_t number) {
    return (number << state->reducerbits) % state->modulus;
}

uint64_t montgomery_out(montgomery_state_t *state, uint64_t number) {
    return (number * state->reciprocal) % state->modulus;
}

uint64_t montgomery_multiply(montgomery_state_t *state, uint64_t a, uint64_t b) {
    assert(a < state->modulus);
    assert(b < state->modulus);

    uint64_t product = a * b;
    uint64_t temp = ((product & state->mask) * state->factor) & state->mask;
    uint64_t reduced = (product + temp * state->modulus) >> state->reducerbits;
    uint64_t result = (reduced < state->modulus) ? reduced : (reduced - state->modulus);
    assert(result < state->modulus);

    return result;
}

/// Calculate the multiplicative inverse of a number with a given modulo using
/// the extended euclidian algorithm.
uint64_t montgomery_multinvmod(uint64_t num, uint64_t mod) {
    assert(num < mod);
    assert(mod > 0);

    int64_t y = num;
    int64_t x = mod;
    int64_t a = 0;
    int64_t b = 1;
    int64_t tmp;

    while(y != 0) {
        tmp = b;
        b = a - x / y * b;
        a = tmp;

        tmp = y;
        y = x % tmp;
        x = tmp;
    }

    if(x == 1) {
        // FIXME inefficient?
        while(a < 0) a += mod;
        return a % mod;
    } else {
        assert(0);
    }
}
