#pragma once
#include <stdint.h>

typedef struct montgomery_state {
    uint64_t modulus;
    uint64_t reducerbits;
    uint64_t reducer;
    uint64_t mask;
    uint64_t reciprocal;
    uint64_t factor;
    uint64_t convertedone;
} montgomery_state_t;

int montgomery_init(montgomery_state_t *state, uint64_t modulus);
uint64_t montgomery_in(montgomery_state_t *state, uint64_t number);
uint64_t montgomery_out(montgomery_state_t *state, uint64_t number);
uint64_t montgomery_multiply(montgomery_state_t *state, uint64_t a, uint64_t b);
uint64_t montgomery_multinvmod(uint64_t num, uint64_t mod);
