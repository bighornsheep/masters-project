#pragma once
#ifndef MPZ_QFORM__INCLUDED
#define MPZ_QFORM__INCLUDED

#include <gmp.h>
#include <stdint.h>

#include "math_mpz.h"
#include "mpz_xgcd.h"
#include "qform_group.h"

extern const group_cost_t mpz_qform_costs;

typedef struct {
  mpz_t a;
  mpz_t b;
  mpz_t c;
} mpz_qform_t;

typedef struct {
  mpz_t x, q, r;
} mpz_qform_reduce_t;

typedef struct {
  mpz_t a1, b1, c1;
  mpz_t a2, b2, c2;
  mpz_t ss, m;
  mpz_t SP, S;
  mpz_t v1, u2, v2;
  mpz_t K, T;
  mpz_t temp, temp2;
  mpz_t R1, R2, C1, C2, M1, M2;
  mpz_t N, L, B;
  mpz_xgcd_t gcd;
} mpz_qform_compose_t;
  
typedef struct {
  // must be the first member of struct
  qform_group_t desc;
  
  mpz_t D;    // the discriminant of the form
  mpz_t S;    // S = floor(|D|^{1/2})
  mpz_t L;    // L = floor(|D|^{1/4})
  
  // temporaries (be careful not to use these temporaries across functions)
  mpz_t t;
  mpz_t m;
  mpz_t n;
  
  mpz_qform_reduce_t reduce;
  mpz_qform_compose_t compose; // nucomp, nudupl, nucube
} mpz_qform_group_t;

void mpz_qform_group_init(mpz_qform_group_t* group);
void mpz_qform_group_clear(mpz_qform_group_t* group);

static inline void mpz_qform_init(mpz_qform_group_t* group, mpz_qform_t* form);
static inline void mpz_qform_clear(mpz_qform_group_t* group, mpz_qform_t* form);
static inline      uint32_t mpz_qform_hash32(mpz_qform_group_t* group, const mpz_qform_t* form);
void               mpz_qform_set_id(mpz_qform_group_t* group, mpz_qform_t* form);
static inline int  mpz_qform_is_id(mpz_qform_group_t* group, const mpz_qform_t* form);
static inline int  mpz_qform_is_ambiguous(mpz_qform_group_t* group, const mpz_qform_t* form);
static inline int  mpz_qform_equal(mpz_qform_group_t* group, const mpz_qform_t* A, const mpz_qform_t* B);
static inline void mpz_qform_set(mpz_qform_group_t* group, mpz_qform_t* R, const mpz_qform_t* A);
static inline void mpz_qform_inverse(mpz_qform_group_t* group, mpz_qform_t* form);
static inline void mpz_qform_print(mpz_qform_group_t* group, const mpz_qform_t* form);


void mpz_qform_compose(mpz_qform_group_t* group, mpz_qform_t* R, const mpz_qform_t* A, const mpz_qform_t* B);
void mpz_qform_square(mpz_qform_group_t* group, mpz_qform_t* R, const mpz_qform_t* A);
void mpz_qform_cube(mpz_qform_group_t* group, mpz_qform_t* R, const mpz_qform_t* A);
void mpz_qform_group_set_discriminant(mpz_qform_group_t* group, const mpz_t D);
void mpz_qform_c(mpz_qform_group_t* group, mpz_t c, mpz_t a, mpz_t b);
void mpz_qform_c1(mpz_qform_group_t* group, mpz_t c, mpz_t a, mpz_t b);
void mpz_qform_reduce(mpz_qform_group_t* group, mpz_qform_t* form);

int mpz_qform_split_ambiguous(mpz_qform_group_t* group, mpz_t d, const mpz_t N, const mpz_qform_t* form);
int mpz_qform_is_primeform(mpz_qform_group_t* group, mpz_qform_t* form, const int p);


/**
 * Inline methods
 */

static inline void mpz_qform_init(mpz_qform_group_t* group, mpz_qform_t* form) {
  mpz_init(form->a);
  mpz_init(form->b);
  mpz_init(form->c);
}

static inline void mpz_qform_clear(mpz_qform_group_t* group, mpz_qform_t* form) {
  mpz_clear(form->a);
  mpz_clear(form->b);
  mpz_clear(form->c);
}

static inline uint32_t mpz_qform_hash32(mpz_qform_group_t* group, const mpz_qform_t* form) {
  // Magic number is smallest prime larger than (1+sqrt(5))/2 * 2^31,
  // which is the golden ratio.
  // NOTE: It's important that we only use the a and c components, since
  //       b is redundant, but more importantly, this means that
  //       the form and its inverse hash to the same value
  //       and will get picked up by collision detection routines.
  const uint32_t magic = 3474701543UL;
  return (((mpz_get_u32(form->a) * magic) + mpz_get_u32(form->c)) * magic);
}


static inline int mpz_qform_is_id(mpz_qform_group_t* group, const mpz_qform_t* form) {
  return (mpz_cmp_ui(form->a, 1) == 0);
}

static inline void mpz_qform_set(mpz_qform_group_t* group, mpz_qform_t* R, const mpz_qform_t* A) {
  mpz_set(R->a, A->a);
  mpz_set(R->b, A->b);
  mpz_set(R->c, A->c);
}

static inline int mpz_qform_equal(mpz_qform_group_t* group, const mpz_qform_t* A, const mpz_qform_t* B) {
  return ((mpz_cmp(A->a, B->a) == 0) && (mpz_cmp(A->b, B->b) == 0) && (mpz_cmp(A->c, B->c) == 0)); // c is computed from a and b
}

static inline void mpz_qform_inverse(mpz_qform_group_t* group, mpz_qform_t* form) {
  if (mpz_cmp(form->a, form->b) != 0 && mpz_cmp(form->a, form->c) != 0) {
    mpz_neg(form->b, form->b);
  }
}

static inline void mpz_qform_print(mpz_qform_group_t* group, const mpz_qform_t* form) {
  gmp_printf("gp.Qfb(%Zd, %Zd, %Zd)\n%Zd\n", form->a, form->b, form->c, group->D);
}

static inline int mpz_qform_is_ambiguous(mpz_qform_group_t* group, const mpz_qform_t* form) {
  return (mpz_cmp_ui(form->a,1) > 0 && (mpz_sgn(form->b) == 0 || mpz_cmp(form->a, form->b) == 0 || mpz_cmp(form->a, form->c) == 0));
}

#endif // MPZ_QFORM__INCLUDED


