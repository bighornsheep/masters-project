/**
 * 128bit negative discriminant qforms
 */

#include <assert.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "gcd_binary_l2r.h"
#include "gcd_brent.h"
#include "gcd_lehmer.h"
#include "gcd_mpz128.h"

#include "math64.h"
#include "math_mpz.h"
#include "primes.h"
#include "s128_t.h"
#include "sqrtmodp_list.h"
#include "s128_qform.h"
#include "u128_t.h"
#include "s128_pow_reps.h"
#include "mpz_qform.h"

// 0 - NUCUBE only
// 1 - Multiply with Square only
// 2 - NUCUBE <= 69 bits.  Multiply with Square otherwise.
#define S128_QFORM_CUBING_STYLE 2

// The number of bits for which Brent's XGCD partial is faster.
#define  XGCD_CROSSOVER 96
#define PXGCD_CROSSOVER 88

// The number of bits in the discriminant for the cubing cross over.
// <= CUBING_CROSSOVER uses genuine cubing.
// > CUBING_CROSSOVER uses multiplication with squaring.
#define CUBING_CROSSOVER 69

/// Average cost to compose, square, and cube in nanoseconds
/// a form with a 118-bit discriminant.
const group_cost_t s128_qform_costs = {
  736.07672,
  623.03309,
#if (S128_QFORM_CUBING_STYLE != 0)
  1359.10981  // This is the cost to multiply with its square.
#else
  1589.80961  // This is the cost to cube.
#endif
};

// Actual GCD methods to use.
#define xgcd_s64(s, t, u, v) xgcd_binary_l2r_s64(s, t, u, v)
#define xgcd_left_s64(s, u, v) xgcd_left_binary_l2r_s64(s, u, v)
#define xgcd_partial_s64(R1, R0, C1, C0, bound) xgcd_partial_brent_s64(R1, R0, C1, C0, bound)

typedef
void xgcd_s128_f(s128_t* d,
		 s128_t* s, s128_t* t,
		 const s128_t* a, const s128_t* b);

typedef
void xgcd_shortpartial_s128_f(s128_t* R1, s128_t* R0,
			      int64_t* C1, int64_t* C0,
			      const int64_t bound);

static xgcd_s128_f* xgcd_s128;
static xgcd_shortpartial_s128_f* xgcd_shortpartial_s128;

//#define xgcd_s128(g, s, t, u, v) xgcd_lehmer_s128_s64l2r(g, s, t, u, v)
//#define xgcd_s128(g, s, t, u, v) xgcd_binary_l2r_s128(g, s, t, u, v)
//#define xgcd_shortpartial_s128(R1, R0, C1, C0, bound) xgcd_shortpartial_brent_s128(R1, R0, C1, C0, bound)
//#define xgcd_shortpartial_s128(R1, R0, C1, C0, bound) xgcd_shortpartial_lehmer_s128_s64l2r(R1, R0, C1, C0, bound)

/*
// Dynamically choose the XGCD.
static void xgcd_s128(s128_t* g,
		      s128_t* s, s128_t* t,
		      const s128_t* u, const s128_t* v) {
  int k1 = numbits_s128(s);
  int k0 = numbits_s128(t);
  int k = k1 > k0 ? k1 : k0;
  if (k < 63) {
    int64_t gg, ss, tt, uu, vv;
    uu = get_s64_from_s128(u);
    vv = get_s64_from_s128(v);
    gg = xgcd_binary_l2r_s64(&ss, &tt, uu, vv);
    set_s128_s64(g, gg);
    set_s128_s64(s, ss);
    set_s128_s64(t, tt);
  } else {
    xgcd_lehmer_s128_s64l2r(g, s, t, u, v);
  }
}

/// Dynamically chose the XGCD.
static void xgcd_shortpartial_s128(s128_t* R1, s128_t* R0,
				   int64_t* C1, int64_t* C0,
				   int64_t bound) {
  int k1 = numbits_s128(R1);
  int k0 = numbits_s128(R0);
  int k = k1 > k0 ? k1 : k0;
  if (k < 63) {
    int64_t r1 = get_s64_from_s128(R1);
    int64_t r0 = get_s64_from_s128(R0);
    xgcd_partial_brent_s64(&r1, &r0, C1, C0, bound);
    set_s128_s64(R1, r1);
    set_s128_s64(R0, r0);
  } else {
    xgcd_shortpartial_lehmer_s128_s64l2r(R1, R0, C1, C0, bound);
  }
}
*/

#ifdef _DEBUG
#define assert64(t, s) \
  if (!s128_is_s64(t)) {						\
    gmp_printf("%s overflows at %d in %s\n", s, __LINE__, __FILE__);	\
    gmp_printf("%"PRId64"<<64 + %"PRIu64"\n", (t)->v1, (t)->v0);	\
    exit(-1);								\
  }
#else
#define assert64(t, s)
#endif

#ifdef _DEBUG
#define assert64_mpz(t, s) \
  if (mpz_sizeinbase(t, 2) > 63) {					\
    gmp_printf("%s overflows at %d in %s\n", s, __LINE__, __FILE__);	\
    gmp_printf("%Zd\n", t);						\
    exit(-1);								\
  }
#else
#define assert64_mpz(t, s)
#endif


static inline int64_t avg_s64(const int64_t a, const int64_t b) {
  // This can't be optimized into an addq/rcrq.
  // Suppose a = 10000000b and b = 01111111b.
  // Then the result is 01111111b when it should be 11111111b.
  s128_t t;
  set_s128_s64(&t, a);
  add_s128_s64(&t, b);
  shr_s128(&t);
  return get_s64_from_s128(&t);
}

// c=(b^2-D)/(4a)
static inline void s128_qform_c(const s128_qform_group_t* group, s128_t* c, int64_t a, int64_t b) {
  s128_t t;
  mul_s128_s64_s64(&t, b, b);
  sub_s128_s128(&t, &group->D);
  shr_s128(&t);
  shr_s128(&t);
  div_s128_s128_s64(c, &t, a);
}

/**
 * r = a >> (size(a)/2) and a > 0
 * r ~ SquareRoot(a)
 */
static inline uint64_t half_rshift_u64(uint64_t a) {
  uint64_t r = a;
#if defined(__x86_64)
  asm("xorq %%rcx, %%rcx\n\t"
      "bsrq %0, %%rcx\n\t"
      "addq $1, %%rcx\n\t"
      "shrq $1, %%rcx\n\t"
      "shrq %%cl, %0"
      : "=rm"(r)
      : "0"(r)
      : "cc", "rcx");
#else
  r >>= ((msb_u64(r)+1)>>1);
#endif
  // should not return zero since this causes weirdness with partial gcd
  if (r == 0) return 1;
  return r;
}

#ifdef _DEBUG
/**
 * Partial Euclidean algorithm.
 * Short partial is preferred since the output arguments tend to fit into a machine word.
 * (for Book's version of NUCOMP, NUDUPL, and NUCUBE algorithm).
 *
 * Input:  R2 = R_{-1} , R1 = R_{0}, bound
 *  - R_i is the R - sequence from "Solving the Pell Equation"
 *    ( R_i = R_{i-2}-q_i R_{i-1} )
 * Output: R2 = R_{i-1}, R1 = R_i, C2 = C_{i-1}, C1 = C_i,
 *  - R_i = 0 or R_i <= bound < R_{i-1}
 *  - C_i sequence from "Solving the Pell Equation" defined as
 *     C_{-1}=0, C_{1}=-1  C_i=C_{i-2}-q_i C_{i-1}
 */
static void xgcd_partial_s128(s128_t* R1, s128_t* R0, s128_t* C1, s128_t* C0, const int64_t bound) {
  s128_t q;
  s128_t t128;
  int64_t t;
  int64_t r0, r1;
  
  set_s128_s64(C0, -1);
  set_s128_s64(C1, 0);

  while (cmp_s128_s64(R0, bound) > 0 && (!s128_is_s64(R1) || !s128_is_s64(R0))) {
    divrem_s128_s128_s128_s128(&q, R1, R1, R0);
    swap_s128_s128(R1, R0);
    
    mul_s128_s128_s128(&t128, &q, C0);
    sub_s128_s128(C1, &t128);
    t128 = (*C0);
    (*C0) = (*C1);
    (*C1) = t128;
  }
  
  if (s128_is_s64(R1) && s128_is_s64(R0)) {
    r0 = R0->v0;
    r1 = R1->v0;
    while (r0 > bound) {
      divrem_s64(&q.v1, &r1, r1, r0);
      t = r1;
      r1 = r0;
      r0 = t;
      
      mul_s128_s128_s64(&t128, C0, q.v1);
      sub_s128_s128(C1, &t128);
      t128 = (*C0);
      (*C0) = (*C1);
      (*C1) = t128;
    }
    
    set_s128_s64(R1, r1);
    set_s128_s64(R0, r0);
  }
}
#endif

/**
 * Partial Euclidean algorithm.
 * (for Book's version of NUCOMP, NUDUPL, and NUCUBE algorithm).
 *
 * Input:  R2 = R_{-1} , R1 = R_{0}, bound
 *  - R_i is the R - sequence from "Solving the Pell Equation"
 *    ( R_i = R_{i-2}-q_i R_{i-1} )
 * Output: R2 = R_{i-1}, R1 = R_i, C2 = C_{i-1}, C1 = C_i,
 *  - R_i = 0 or R_i <= bound < R_{i-1}
 *  - C_i sequence from "Solving the Pell Equation" defined as
 *     C_{-1}=0, C_{1}=-1  C_i=C_{i-2}-q_i C_{i-1}
 */
/*
static void xgcd_shortpartial_divrem_s128(s128_t* R1, s128_t* R0, int64_t* C1, int64_t* C0, const int64_t bound) {
    s128_t q;
    int64_t t;
    int64_t r0, r1;

    (*C0) = -1;
    (*C1) = 0;

    // 128bit div/rem
    while (cmp_s128_s64(R0, bound) > 0 && (!s128_is_s64(R1) || !s128_is_s64(R0))) {
        divrem_s128_s128_s128_s128(&q, R1, R1, R0);
        swap_s128_s128(R1, R0);

        (*C1) -= ((int64_t)q.v0) * (*C0);
        t = (*C0);
        (*C0) = (*C1);
        (*C1) = t;
    }

    if (s128_is_s64(R0) && s128_is_s64(R1)) {
        r0 = R0->v0;
        r1 = R1->v0;
        // 64bit div/rem
        while (r0 > bound) {
            divrem_s64(&q.v1, &r1, r1, r0);
            t = r1;
            r1 = r0;
            r0 = t;

            (*C1) -= q.v1 * (*C0);
            t = (*C0);
            (*C0) = (*C1);
            (*C1) = t;
        }
        set_s128_s64(R1, r1);
        set_s128_s64(R0, r0);
    }
}
*/

void s128_qform_set_id(s128_qform_group_t* group, s128_qform_t* form) {
  form->a = 1;
  form->b = (group->D.v0 & 3) != 0;
  s128_qform_c(group, &form->c, form->a, form->b);
}

/**
 * Check for a primeform where a=p.
 * Tests b = +/- sqrt(D) mod p.
 * Note: -0 mod p is p.  This way we try ambiguous forms.
 */
int s128_qform_is_primeform(s128_qform_group_t* group,
			    s128_qform_t* form,
			    const int p) {
  int Dmodp;
  const short* sqrtp;
  
  if (p > sqrtmodp_maxp)
    return 0; // p is too large for the table
  
  sqrtp = sqrtmodp[p];
  if (sqrtp == 0)
    return 0; // p was not prime
  
  // a = p
  form->a = p;
  
  // b = sqrt(D) (mod p)
  Dmodp = mod_s64_s128_s64(&group->D, p);
  if (Dmodp < 0) {
    Dmodp += p;
  }
  if (Dmodp == 0) {
    // prime divides discriminant
    return 0;
  }
  form->b = sqrtp[Dmodp];
  if (form->b == -1)
    return 0; // D does not have a sqrt (mod p)
  
  // We know that p | b^2-D for +/- b.
  // Compute c = (b^2 - d) / (4a) if possible.
  if (p == 2) {
    // special case if p == 2
    mul_s128_s64_s64(&form->c, form->b, form->b);
    sub_s128_s128(&form->c, &group->D);
    if ((form->c.v0 & 7) == 0) {
      // 4a | b^2-D, a=2, so we divide by 8
      shr_s128(&form->c);
      shr_s128(&form->c);
      shr_s128(&form->c);
      return 1;
    }
    if (form->b == 1)
      return 0;
    
    // try -b = 2
    form->b = 2;
    set_s128_s64(&form->c, 4);
    sub_s128_s128(&form->c, &group->D);
    if ((form->c.v0 & 7) != 0)
      return 0;
    
    // 4a | b^2-D, a=2, so we divide by 8
    shr_s128(&form->c);
    shr_s128(&form->c);
    shr_s128(&form->c);
    return 1;
  }
  
  // p != 2
  mul_s128_s64_s64(&form->c, form->b, form->b);
  sub_s128_s128(&form->c, &group->D);
  
  if ((form->c.v0 & 3) == 0) {
    // 4a | b^2-D
    // divide by 4
    shr_s128(&form->c);
    shr_s128(&form->c);
    // divide by a
    div_s128_s128_s64(&form->c, &form->c, form->a);
    return 1;
  }
  
  // try -b
  form->b = p-form->b;
  
  // make sure that 4 | b^2-D
  mul_s128_s64_s64(&form->c, form->b, form->b);
  sub_s128_s128(&form->c, &group->D);
  if ((form->c.v0 & 3) != 0)
    return 0;
  
  // divide by 4
  shr_s128(&form->c);
  shr_s128(&form->c);
  
  // divide by a
  div_s128_s128_s64(&form->c, &form->c, form->a);
  return 1;
}

/**
 * If N is split by the form (a,b,C), then d is a factor and 1 is returned
 * otherwise 0 is returned.
 * Assumes that form is an ambiguous form.
 */
int s128_qform_split_ambiguous(s128_qform_group_t* group, mpz_t out_d, const mpz_t in_N, const s128_qform_t* form) {
  s128_t m;
  s128_t d;
  s128_t N;
  mpz_get_s128(&N, in_N);
  
  // ambiguous forms are of three types
  // (a, 0, c), (a, a, c), and (a, b, a)
  if (form->b == 0 || form->a == form->b) {
    // first two forms (a, 0, c) and (a, a, c)
    // (a, 0, c) => -4ac = D => a | D
    // (a, a, c) => a^2-4ac = a(a-4c) = D => a | d
    set_s128_s64(&m, form->a);
  } else {
    // final form (a, b, a)
    // (a, b, a) => b^2 - 4a^2 = (b+2a)(b-2a) = D
    // pick m = min |b +/- 2a|
    set_s128_s64(&m, form->a);
    shl_s128(&m);
    if (form->b < 0) {
      add_s128_s64(&m, form->b);
    } else {
      sub_s128_s64(&m, form->b);
    }
  }
  abs_s128_s128(&m, &m);
  gcd_binary_l2r_u128((u128_t*)&d, (const u128_t*)&m, (const u128_t*)&N);
  if (cmp_s128_s64(&d, 1) > 0 && cmp_s128_s128(&d, &N) < 0) {
    mpz_set_s128(out_d, &d);
    return 1;
  }
  return 0;
}

/**
 * initialize the group
 */
void s128_qform_group_init(s128_qform_group_t* group) {
  group->desc.group.elem_init = (group_elem_init_f*)&s128_qform_init;
  group->desc.group.elem_clear = (group_elem_clear_f*)&s128_qform_clear;
  group->desc.group.elem_size = sizeof(s128_qform_t);
  
  group->desc.group.hash32 = (group_hash32_f*)&s128_qform_hash32;
  group->desc.group.set_id = (group_set_id_f*)&s128_qform_set_id;
  group->desc.group.is_id = (group_is_id_f*)&s128_qform_is_id;
  group->desc.group.set = (group_set_f*)&s128_qform_set;
  group->desc.group.equal = (group_equal_f*)&s128_qform_equal;
  group->desc.group.inverse = (group_inverse_f*)&s128_qform_inverse;
  group->desc.group.compose = (group_compose_f*)&s128_qform_compose;
  group->desc.group.square = (group_square_f*)&s128_qform_square;
  group->desc.group.cube = (group_cube_f*)&s128_qform_cube;
  group->desc.group.print = (group_print_f*)&s128_qform_print;
  
  group->desc.discriminant_max_bits = s128_qform_group_max_bits;
  group->desc.clear = (qform_group_clear_f*)&s128_qform_group_clear;
  group->desc.set_discriminant = (qform_group_set_discriminant_f*)&s128_qform_group_set_discriminant;
  group->desc.reduce = (qform_reduce_f*)&s128_qform_reduce;
  group->desc.is_primeform = (qform_is_primeform_f*)&s128_qform_is_primeform;
  group->desc.is_ambiguous = (qform_is_ambiguous_f*)&s128_qform_is_ambiguous;
  group->desc.split_ambiguous = (qform_split_ambiguous_f*)&s128_qform_split_ambiguous;
  
  group->desc.pow_rep_sizes = s128_pow_rep_sizes;
  group->desc.pow_reps = s128_pow_reps;

  mpz_init(group->tmp);
  mpz_init(group->tmp2);
  mpz_init(group->tmp3);
}

/**
 * release the group data
 */
void s128_qform_group_clear(s128_qform_group_t* group) {
  mpz_clear(group->tmp);
  mpz_clear(group->tmp2);
  mpz_clear(group->tmp3);
}

/**
 * Saves the discriminant and computes the fourth root.
 * Assumes that D is negative.
 */
void s128_qform_group_set_discriminant(s128_qform_group_t* group,
				       const mpz_t D) {
  s128_t root;
  mpz_get_s128(&group->D, D);
  // compute roots
  abs_s128_s128(&root, &group->D);
  sqrt_s128_s128(&root, &root);
  group->S = get_u64_from_u128((u128_t*)&root);
  sqrt_s128_s128(&root, &root);
  group->L = get_u64_from_u128((u128_t*)&root);

  // Setup 128-bit XGCD based on size of discriminant.
  // 96-bits is the changing point.  This came from timing both
  // XGCD methods in the context of ideal arithmetic.
  long n = mpz_sizeinbase(D, 2);
  if (n <= XGCD_CROSSOVER) {
    xgcd_s128 = &xgcd_binary_l2r_s128;
  } else {
    xgcd_s128 = &xgcd_lehmer_s128_s64l2r;
  }
  if (n <= PXGCD_CROSSOVER) {
    xgcd_shortpartial_s128 = &xgcd_shortpartial_brent_s128;
  } else {
    xgcd_shortpartial_s128 = &xgcd_shortpartial_lehmer_s128_brent64;
  }
}

/**
 * Saves the discriminant and computes the fourth root.
 * Assumes that D is negative.
 */
void s128_qform_group_set_discriminant_s128(s128_qform_group_t* group, const s128_t* D) {
  s128_t root;
  set_s128_s128(&group->D, D);

  // compute roots
  abs_s128_s128(&root, &group->D);
  sqrt_s128_s128(&root, &root);
  group->S = get_u64_from_u128((u128_t*)&root);
  sqrt_s128_s128(&root, &root);
  group->L = get_u64_from_u128((u128_t*)&root);
}

void s128_qform_reduce(s128_qform_group_t* group, s128_qform_t* form) {
#if !defined(__x86_64)
  s128_t t;
#endif
  int64_t q;
  int64_t r;
  
  // First make sure form is positive definite
  while (1) {
    if (cmp_s64_s128(form->a, &form->c) > 0) {
      // Swap a and c and invert b.
      r = form->a;
      form->a = (int64_t)form->c.v0;
      set_s128_s64(&form->c, r);
      form->b = -form->b;
    }
    if ((form->b > form->a) || (form->b <= -form->a)) {
      // Find r such that -a < r <= a
      // and r = b (mod 2a).
      // q = b/2a = (b/a)/2
      // r = b%2a = (b%a) +- a*((b/a)%2)
      // where '/' is divide with floor.
      // NOTE: $a$ is always positive. $b$ may be negative, however,
      //       so $q$ and $r$ may be negative.
      divrem_s64(&q, &r, form->b, form->a);
      // Add/sub 1 to q while bringing r closer to 0.
      int64_t qm = -(q & 1);            // qm = q & 1 ? -1 : 0;
      int64_t rm = (r <= 0) - 1;        // rm = r <= 0 ? 0 : -1;
      r += ((form->a ^ rm) - rm) & qm;  // move r towards 0
      q -= (rm | 1) & qm;               // make q even
      q >>= 1;
      
      // c -= (q * (b+r)) / 2
#if defined(__x86_64)
      asm("movq %5, %%rax\n"
	  "addq %6, %%rax\n"
	  "imulq %4\n"
	  "sarq $1, %%rdx\n"
	  "rcrq $1, %%rax\n"
	  "subq %%rax, %1\n"
	  "sbbq %%rdx, %0\n"
	  : "=r"(form->c.v1), "=r"(form->c.v0)
	  : "0"(form->c.v1), "1"(form->c.v0), "r"(q), "r"(form->b), "r"(r)
	  : "cc", "rax", "rdx");
#else
      mul_s128_s64_s64(&t, q, form->b+r);
      shr_s128(&t);
      sub_s128_s128(&form->c, &t);
#endif
            
      form->b = r;
    } else {
      break;
    }
  }
  
  // account for special case
  if (is_equal_s128_s64(&form->c, form->a) && form->b < 0) {
    form->b = -form->b;
  }
}

/**
 * NUCOMP algorithm. Adapted from "Solving the Pell Equation"
 * by Michael J. Jacobson, Jr. and Hugh C. Williams.
 * http://www.springer.com/mathematics/numbers/book/978-0-387-84922-5
 */
void s128_qform_compose(s128_qform_group_t* group,
			s128_qform_t* C,
			const s128_qform_t* A,
			const s128_qform_t* B) {
  int64_t a1, a2, b1, b2;
  s128_t c2;
  int64_t g, s, x, y, z;
  int64_t m12, p12, u;
  int64_t C1, C0, r1, r0, m1, m2;
  int64_t bound;
  s128_t tmp;
  s128_t tmp2;
  
#ifdef _DEBUG
  mul_s128_s64_s64(&tmp, A->b, A->b);
  mul_s128_s128_s64(&tmp2, &A->c, A->a);
  shl_s128_int(&tmp2, 2);
  sub_s128_s128(&tmp, &tmp2);
  assert(cmp_s128_s128(&group->D, &tmp) == 0);
  
  mul_s128_s64_s64(&tmp, B->b, B->b);
  mul_s128_s128_s64(&tmp2, &B->c, B->a);
  shl_s128_int(&tmp2, 2);
  sub_s128_s128(&tmp, &tmp2);
  assert(cmp_s128_s128(&group->D, &tmp) == 0);
#endif  // _DEBUG
    
  // if A is the identity, return B
  if (s128_qform_is_id(group, A)) {
    s128_qform_set(group, C, B);
    s128_qform_reduce(group, C);
    return;
  }
  
  // if B is the identity, return A
  if (s128_qform_is_id(group, B)) {
    s128_qform_set(group, C, A);
    s128_qform_reduce(group, C);
    return;
  }
  
  // Make sure Norm(A) > Norm(B).
  if (A->a > B->a) {
    a1 = A->a;
    b1 = A->b;
    a2 = B->a;
    b2 = B->b;
    c2 = B->c;
  } else {
    a1 = B->a;
    b1 = B->b;
    a2 = A->a;
    b2 = A->b;
    c2 = A->c;
  }
  assert(a2 != 0);
  
  // Compute d1=gcd(a1, a2) and u1 such that  $a2 | (d1 - a2 * u1)$
  g = xgcd_left_s64(&x, a2, a1);
  
  // Compute gcd((b1+b2)/2, g) = s = y * (b1+b2)/2 + z * g
  p12 = avg_s64(b1, b2);
  m12 = p12 - b2;
  
  // Compute u = x*z*(b1-b2)/2 - y*c mod a1
  u = mulmod_s64(x, m12, a1);
  s = 1;
  if (g != 1) {
    s = xgcd_s64(&y, &z, p12, g);
    if (s != 1) {
      a1 /= s;
      a2 /= s;
    }
    
    u = mulmod_s64(u, z, a1);
    int64_t t = mulmod_s64(y, mod_s64_s128_s64(&c2, a1), a1);
    u = submod_s64(u, t, a1);
  }
  u += a1 & (u >> 63);  // if (u < 0) u += a1;
  
  // compute bounds
  bound = half_rshift_u64(a1 / a2) * group->L;
  //  bound = sqrt((double)group->S * (double)a1/(double)a2);
  //  bound = sqrt_u64(a1 / a2) * group->L;
  if (a1 <= bound) {
    // normal composition
    C->a = a1*a2;
    
    // C->b = 2*a2*u + b2 (mod 2C->a)
    mul_s128_s64_s64(&tmp, a2, u);
    shl_s128(&tmp);
    add_s128_s64(&tmp, b2);
    C->b = mod_s64_s128_u64(&tmp, C->a << 1);
  } else {
    // nucomp steps
    r1 = a1;
    r0 = u;
    
    // partial xgcd
    xgcd_partial_s64(&r1, &r0, &C1, &C0, bound);
    
    // m1 = (a2 * r0 + m12 * C0)/a1
    m1 = muladdmuldiv_s64(a2, r0, m12, C0, a1);
    
    // m2 = (p12*r0 - s*C0*c2)/a1
    mul_s128_s64_s64(&tmp, p12, r0);
    mul_s128_s64_s64(&tmp2, s, C0);
    mul_s128_s128_s128(&tmp2, &tmp2, &c2);
    sub_s128_s128(&tmp, &tmp2);
    div_s128_s128_s64(&tmp, &tmp, a1);
    assert64(&tmp, "m2");
    m2 = get_s64_from_s128(&tmp);
    
    // a_{i+1} = r0*m1 - C0*m2
    muladdmul_s128_4s64(&tmp, r0, m1, -C0, m2);
    assert64(&tmp, "C->a");
    C->a = get_s64_from_s128(&tmp);
    
    // b_{i+1}= 2(a2*r0 - a*|C1|) / C0 - b2  (mod 2a)
    muladdmul_s128_4s64(&tmp, a2, r0, -C->a, abs_s64(C1));
    div_s128_s128_s64(&tmp, &tmp, C0);
    u = mod_s64_s128_u64(&tmp, C->a << 1);
    C->b = (u - b2) + u; // = 2u -b2; possibly avoids overflow
  }
  s128_qform_c(group, &C->c, C->a, C->b);
  s128_qform_reduce(group, C);
}

/**
 * NUDUPL. Simplified from compose above.
 */
void s128_qform_square(s128_qform_group_t* group, s128_qform_t* C, const s128_qform_t* A) {
  int64_t a1, b1;
  s128_t c1;
  int64_t s, y;
  int64_t u;
  int64_t C1, C0, r1, r0, m2;
  s128_t tmp;
  s128_t tmp2;
  
  // if A is the identity, return the identity
  if (s128_qform_is_id(group, A)) {
    s128_qform_set_id(group, C);
    return;
  }
  
  // Local copies
  a1 = A->a;
  b1 = A->b;
  c1 = A->c;
  
  // Compute d1=gcd(a1, a2) and u1 such that  $a2 | (d1 - a2 * u1)$
  // Compute gcd((b1+b2)/2, g) = s = y * (b1+b2)/2 + z * g
  // Compute u = x*z*(b1-b2)/2 - y*c mod a1
  s = xgcd_left_s64(&y, b1, a1);
  a1 /= s;
  
  r0 = mod_s64_s128_s64(&c1, a1);
  r1 = mulmod_s64(y, r0, a1);
  u = -r1;
  u += a1 & (u >> 63);  // if (u < 0) u += a1;
  
  if (a1 <= group->L) {
    // normal composition
    C->a = a1*a1;
    
    // C->b = 2*a1*u + b1 (mod 2C->a)
    mul_s128_s64_s64(&tmp, a1, u);
    shl_s128(&tmp);
    add_s128_s64(&tmp, b1);
    C->b = mod_s64_s128_u64(&tmp, C->a << 1);
  } else {
    // nucomp steps
    r1 = a1;
    r0 = u;
    
    // partial xgcd
    xgcd_partial_s64(&r1, &r0, &C1, &C0, group->L);
    
    // m2 = (b1 * r0 - s*C0*c1) / a1
    mul_s128_s64_s64(&tmp, b1, r0);
    mul_s128_s64_s64(&tmp2, s, C0);
    mul_s128_s128_s128(&tmp2, &tmp2, &c1);
    sub_s128_s128(&tmp, &tmp2);
    div_s128_s128_s64(&tmp, &tmp, a1);
    assert64(&tmp, "m2");
    m2 = get_s64_from_s128(&tmp);
    
    // a_{i+1} = r0^2 - C0*m2
    muladdmul_s128_4s64(&tmp, r0, r0, -C0, m2);
    assert64(&tmp, "C->a");
    C->a = get_s64_from_s128(&tmp);
    
    // b_{i+1} = (a1 * r0 - a * |C1|)/C0  (mod 2a)
    muladdmul_s128_4s64(&tmp, a1, r0, -C->a, abs_s64(C1));
    div_s128_s128_s64(&tmp, &tmp, C0);
    u = mod_s64_s128_u64(&tmp, C->a << 1);
    C->b = (u - b1) + u; // = 2u -b2; possibly avoids overflow
  }
  
  s128_qform_c(group, &C->c, C->a, C->b);
  s128_qform_reduce(group, C);
}

/**
 * Computes a reduced ideal equivalent to the cube of an ideal.
 * Adapted from "Fast Ideal Cubing in Imaginary Quadratic Number
 * and Function Fields" by Laurent Imbert, Michael J. Jacobson, Jr. and
 * Arthur Schmidt.
 * www.lirmm.fr/~imbert/pdfs/cubing_amc_2010.pdf
 */
#if (S128_QFORM_CUBING_STYLE != 1)

// r = (a * b) % m
// mz should be filled out already.  az and bz will be filled out only if needed.
static inline void mulmod_mixed(s128_t* r, s128_t* a, s128_t* b, s128_t* m,
				mpz_t az, mpz_t bz, mpz_t mz, int mbits) {
  int n = numbits_s128(a) + numbits_s128(b);
  if (n > 127) {
    s128_to_mpz(a, az);
    s128_to_mpz(b, bz);
    mpz_mul(az, az, bz);
    if (n >= mbits) {
      mpz_fdiv_r(az, az, mz);
    }
    s128_from_mpz(r, az);
  } else {
    mul_s128_s128_s128(r, a, b);
    if (n >= mbits) {
      mod_s128_s128_s128(r, r, m);
    }
  }
}

// r = (a * b) % m
// mz should be filled out already.  az and bz will be filled out only if needed.
static inline void mulmod_mixed64(s128_t* r, s128_t* a, int64_t b, s128_t* m,
				  mpz_t az, mpz_t bz, mpz_t mz, int mbits) {
  int n = numbits_s128(a) + numbits_s64(b);
  if (n > 127) {
    mpz_set_s128(az, a);
    mpz_set_s64(bz, b);
    mpz_mul(az, az, bz);
    if (n >= mbits) {
      mpz_fdiv_r(az, az, mz);
    }
    s128_from_mpz(r, az);
  } else {
    mul_s128_s128_s64(r, a, b);
    if (n >= mbits) {
      mod_s128_s128_s128(r, r, m);
    }
  }
}

static void s128_qform_genuine_cube(s128_qform_group_t* group,
				    s128_qform_t* R,
				    const s128_qform_t* A) {
  int64_t a1;
  int64_t b1;
  s128_t c1;
  int64_t SP;
  int64_t v1;
  int64_t N;
  int64_t B;
  int64_t T;
  int64_t S, u2, v2;
  int64_t R1, C1, C2;
  int64_t M1, M2;
  int64_t L_64;
  int64_t t1_64, t2_64;
  s128_t L;
  s128_t K;
  s128_t temp;
  s128_t temp2;
  s128_t S_128, u2_128, v2_128;
  s128_t R1_128, R2_128;
  s128_t T_128;
  
  // if A is the identity, return the identity
  if (s128_qform_is_id(group, A)) {
    s128_qform_set_id(group, R);
    return;
  }
  
  // cubing an ambiguous form results in the original form
  if (s128_qform_is_ambiguous(group, A)) {
    s128_qform_set(group, R, A);
    s128_qform_reduce(group, R);
    return;
  }
  
  a1 = A->a;
  b1 = A->b;
  c1 = A->c;
  
  // solve SP = v1 b + u1 a (only need v1)
  SP = xgcd_left_s64(&v1, b1, a1);
  
  if (SP == 1) {
    // N = a
    N = a1;
    
    // L = a^2 (this is always positive)
    mul_s128_s64_s64(&L, a1, a1);
    
    // K = c v1 (v1(b - a c v1) - 2) (mod L)
    if (s128_is_s64(&L)) {
      // 64bit modulus
      L_64 = get_s64_from_s128(&L);
      
      // use the remainder of v1 nearest to zero (L >= 0)
      v1 %= L_64;
      v1 += cond_negate_s64((L_64>>1)-v1, L_64);
      t2_64 = mod_s64_s128_s64(&c1, L_64);
      t1_64 = mulmod_s64(t2_64, a1, L_64);
      t1_64 = mulmod_s64(t1_64, v1, L_64);
      t1_64 = mulmod_s64(b1 - t1_64, v1, L_64) - 2;
      t1_64 = mulmod_s64(t1_64, v1, L_64);
      t1_64 = mulmod_s64(t1_64, t2_64, L_64);
      // use a positive remainder
      t1_64 += L_64 & (t1_64 >> 63);  // Add L if t1 < 0.
      set_s128_s64(&K, t1_64);
    } else {
      // 128bit modulus (use GMP)
      mpz_set_s128(group->tmp3, &c1);
      mpz_mul_s64(group->tmp, group->tmp3, a1);
      mpz_mul_s64(group->tmp, group->tmp, v1);
      mpz_set_s64(group->tmp2, b1);
      mpz_sub(group->tmp, group->tmp2, group->tmp);
      mpz_mul_s64(group->tmp, group->tmp, v1);
      mpz_sub_ui(group->tmp, group->tmp, 2);
      mpz_mul_s64(group->tmp, group->tmp, v1);
      mpz_mul(group->tmp, group->tmp, group->tmp3);
      mpz_set_s128(group->tmp2, &L);
      mpz_fdiv_r(group->tmp, group->tmp, group->tmp2);
      mpz_get_s128(&K, group->tmp);
      if (cmpzero_s128(&K) < 0) {
	add_s128_s128(&K, &L);
      }
      /* // Use mixed mode arithmetic. */
      /* int Lbits = numbits_s128(&L); */
      /* s128_to_mpz(&L, group->tmp3); */
      /* s128_t c_mod_L; */
      /* mod_s128_s128_s128(&c_mod_L, &c1, &L); */
      /* mulmod_mixed64(&K, &c_mod_L, a1, &L, group->tmp, group->tmp2, group->tmp3, Lbits); */
      /* mulmod_mixed64(&K, &K, v1, &L, group->tmp, group->tmp2, group->tmp3, Lbits); */
      /* sub_s128_s64(&K, b1); */
      /* neg_s128_s128(&K, &K); */
      /* mulmod_mixed64(&K, &K, v1, &L, group->tmp, group->tmp2, group->tmp3, Lbits); */
      /* sub_s128_s64(&K, 2); */
      /* mulmod_mixed64(&K, &K, v1, &L, group->tmp, group->tmp2, group->tmp3, Lbits); */
      /* mulmod_mixed(&K, &K, &c_mod_L, &L, group->tmp, group->tmp2, group->tmp3, Lbits); */
      /* if (cmpzero_s128(&K) < 0) { */
      /* 	add_s128_s128(&K, &L); */
      /* } */
    }
  } else {
    // S = u2 (a SP) + v2 (b^2 - ac)
    // (this has been verified up to 121bit discriminants
    //  to produce the correct results)
    mul_s128_s64_s64(&temp2, b1, b1);
    mul_s128_s128_s64(&temp, &c1, a1);
    sub_s128_s128(&temp2, &temp);
    mul_s128_s64_s64(&temp, SP, a1);
    if (s128_is_s64(&temp) && s128_is_s64(&temp2)) {
      // (a*SP) and (b^2-ac) fit into 64bit each
      S = xgcd_s64(&u2, &v2, get_s64_from_s128(&temp), get_s64_from_s128(&temp2));
      set_s128_s64(&u2_128, u2);
    } else {
      // 128bit xgcd required
      xgcd_s128(&S_128, &u2_128, &v2_128, &temp, &temp2);
      assert64(&S_128, "S");
      assert64(&v2_128, "v2");
      S = get_s64_from_s128(&S_128);
      v2 = get_s64_from_s128(&v2_128);
    }
    
    // N = a/S
    N = a1 / S;
    
    // L = N a
    mul_s128_s64_s64(&L, N, a1);
    
    // K = -c(v1 u2 a + v2 b) (mod L)
    if (s128_is_s64(&L)) {
      // 64bit modulus
      L_64 = get_s64_from_s128(&L);
      t1_64 = mod_s64_s128_s64(&u2_128, L_64);
      t1_64 = mulmod_s64(t1_64, v1, L_64);
      t1_64 = mulmod_s64(t1_64, a1, L_64);
      t2_64 = mulmod_s64(v2, b1, L_64);
      t1_64 = addmod_s64(t1_64, t2_64, L_64);
      t2_64 = mod_s64_s128_s64(&c1, L_64);
      t1_64 = mulmod_s64(t1_64, -t2_64, L_64);
      t1_64 += L_64 & (t1_64 >> 63);  // Add L if t1 < 0.
      set_s128_s64(&K, t1_64);
    } else {
      // 128bit modulus (use GMP)
      mpz_set_s128(group->tmp, &u2_128);
      mpz_mul_s64(group->tmp, group->tmp, v1);
      mpz_mul_s64(group->tmp, group->tmp, a1);
      mul_s128_s64_s64(&temp, v2, b1);
      mpz_set_s128(group->tmp2, &temp);
      mpz_add(group->tmp, group->tmp, group->tmp2);
      mpz_set_s128(group->tmp2, &c1);
      mpz_mul(group->tmp, group->tmp, group->tmp2);
      mpz_neg(group->tmp, group->tmp);
      mpz_set_s128(group->tmp2, &L);
      mpz_fdiv_r(group->tmp, group->tmp, group->tmp2);
      mpz_get_s128(&K, group->tmp);
      /* // Use mixed mode arithmetic. */
      /* s128_t c_mod_L; */
      /* mod_s128_s128_s128(&c_mod_L, &c1, &L); */
      /* mpz_set_s128(group->tmp3, &L); */
      /* int Lbits = numbits_s128(&L); */
      /* mulmod_mixed64(&K, &u2_128, v1, &L, group->tmp, group->tmp2, group->tmp3, Lbits); */
      /* mulmod_mixed64(&K, &K, a1, &L, group->tmp, group->tmp2, group->tmp3, Lbits); */
      /* mul_s128_s64_s64(&temp, v2, b1); */
      /* mod_s128_s128_s128(&temp, &temp, &L); */
      /* add_s128_s128(&K, &temp); */
      /* mulmod_mixed(&K, &K, &c_mod_L, &L, group->tmp, group->tmp2, group->tmp3, Lbits); */
      /* neg_s128_s128(&K, &K); */
      if (cmpzero_s128(&K) < 0) {
	add_s128_s128(&K, &L);
      }
    }
    // C = Sc
    mul_s128_s128_s64(&c1, &c1, S);
  }
  
  // Compute NUCOMP termination bound,
  // sqrt(a1)*root(4, group->D/4)
  // This is roughly (a1*group->S/2) shifted right by half the number
  // of bits.
  mul_s128_s64_s64(&temp, group->S, a1);
  shr_s128_int(&temp, (msb_s128(&temp) + 1) >> 1); // approximate sqrt
  //  mul_s128_s64_s64(&temp, group->S, a1);
  //  shr_s128(&temp);
  //  sqrt_s128_s128(&temp, &temp);
  assert64(&temp, "B");
  B = get_s64_from_s128(&temp);
  B |= !B;
  if (cmp_s128_s64(&L, B) < 0) {
    // compute with regular cubing formula (result will be reduced)
    // T = NK
    assert64(&K, "K");
    T = N * get_s64_from_s128(&K);
    // R.a = NL
    assert64(&L, "L");
    R->a = N * get_s64_from_s128(&L);
    // C.b = b + 2 T
    R->b = b1 + (T << 1);
    // C.c = (S c + K (T + b)) / L
    // K is apparently 64bit
    mul_s128_s64_s64(&R->c, get_s64_from_s128(&K), T + b1);
    add_s128_s128(&R->c, &c1);
    div_s128_s128_s128(&R->c, &R->c, &L);
  } else {
    // use NUCOMP formulas
    
    // Execute partial reduction
    // (This has been verified upto 121bit discriminants
    //  to produce the correct answer while only needing
    //  64bit C1,C2).
    set_s128_s128(&R2_128, &L);
    set_s128_s128(&R1_128, &K);
    if (s128_is_s64(&R2_128) && s128_is_s64(&R1_128)) {
      int64_t R2 = get_s64_from_s128(&R2_128);
      R1 = get_s64_from_s128(&R1_128);
      xgcd_partial_s64(&R2, &R1, &C2, &C1, B);
      set_s128_s64(&R2_128, R2);
      set_s128_s64(&R1_128, R1);
    } else {
      xgcd_shortpartial_s128(&R2_128, &R1_128, &C2, &C1, B);
      assert64(&R1_128, "R1");
      R1 = get_s64_from_s128(&R1_128);
    }
    
#ifdef _DEBUG
    // verify that shortpartial gives the right answer
    set_s128_s128(&R2_128, &L);
    set_s128_s128(&R1_128, &K);
    xgcd_partial_s128(&R2_128, &R1_128, &temp2, &temp, B);
    assert(cmp_s128_s64(&temp2, C2) == 0);
    assert(cmp_s128_s64(&temp,  C1) == 0);
    assert(cmp_s128_s64(&R1_128, R1) == 0);
#endif
    
    // T = N K (mod L)
    if (numbits_s128(&K) + numbits_s64(N) > 126) {
      mpz_set_s128(group->tmp, &K);
      mpz_mul_s64(group->tmp, group->tmp, N);
      mpz_set_s128(group->tmp2, &L);
      mpz_fdiv_r(group->tmp, group->tmp, group->tmp2);
      mpz_get_s128(&T_128, group->tmp);
    } else {
      mul_s128_s128_s64(&T_128, &K, N);
      mod_s128_s128_s128(&T_128, &T_128, &L);
    }
    if (cmpzero_s128(&T_128) < 0) {
      add_s128_s128(&T_128, &L);
    }
    
    // M1 = (N R1 + T C1) / L
    if (numbits_s128(&T_128) + numbits_s64(C1) <= 127) {
      mul_s128_s64_s64(&temp, R1, N);
      mul_s128_s128_s64(&temp2, &T_128, C1);
      add_s128_s128(&temp, &temp2);
      div_s128_s128_s128(&temp, &temp, &L);
      assert64(&temp, "M1");
      M1 = get_s64_from_s128(&temp);
    } else {
      // use GMP
      mul_s128_s64_s64(&temp, R1, N);
      mpz_set_s128(group->tmp, &T_128);
      mpz_mul_s64(group->tmp, group->tmp, C1);
      mpz_set_s128(group->tmp2, &temp);
      mpz_add(group->tmp, group->tmp, group->tmp2);
      mpz_set_s128(group->tmp2, &L);
      mpz_divexact(group->tmp, group->tmp, group->tmp2);
      assert64_mpz(group->tmp, "M1");
      M1 = mpz_get_s64(group->tmp);
    }
    
    // M2 = (R1(b + T) - c S C1) / L
    if (s128_is_s64(&T_128) && s128_is_s64(&c1)) {
      mul_s128_s64_s64(&temp, b1 + get_s64_from_s128(&T_128), R1);
      mul_s128_s64_s64(&temp2, get_s64_from_s128(&c1), C1);
      sub_s128_s128(&temp, &temp2);
      div_s128_s128_s128(&temp, &temp, &L);
      assert64(&temp, "M2");
      M2 = get_s64_from_s128(&temp);
    } else {
      // use GMP
      mpz_set_s128(group->tmp, &T_128);
      mpz_set_s64(group->tmp2, b1);
      mpz_add(group->tmp, group->tmp, group->tmp2);
      mpz_mul_s64(group->tmp, group->tmp, R1);
      mpz_set_s128(group->tmp2, &c1);
      mpz_mul_s64(group->tmp2, group->tmp2, C1);
      mpz_sub(group->tmp, group->tmp, group->tmp2);
      mpz_set_s128(group->tmp2, &L);
      mpz_divexact(group->tmp, group->tmp, group->tmp2);
      assert64_mpz(group->tmp, "M2");
      M2 = mpz_get_s64(group->tmp);
    }
    
    // C.a = (-1)^(i-1) (R1 M1 - C1 M2)
    R->a = R1*M1 - C1*M2;
    R->a = cond_negate_s64(-C1, R->a);
    
    // C.b = 2 (N R1 - C.a C2) / C1 - b (mod 2C.a)
    muladdmul_s128_4s64(&temp2, N, R1, -R->a, C2);
    shl_s128(&temp2);
    div_s128_s128_s64(&temp2, &temp2, C1);
    sub_s128_s64(&temp2, b1);
    set_s128_s64(&temp, R->a);
    shl_s128(&temp);
    mod_s128_s128_s128(&temp2, &temp2, &temp);
    if (cmpzero_s128(&temp2) < 0) {
      add_s128_s128(&temp2, &temp);
    }
    assert64(&temp2, "R->b");
    R->b = get_s64_from_s128(&temp2);
    
    R->a = abs_s64(R->a);
    
    // C.c = (C.b^2 - D) / 4 C.a
    s128_qform_c(group, &R->c, R->a, R->b);
  }
  
  // normalize and reduce
  s128_qform_reduce(group, R);
}
#endif

#if (S128_QFORM_CUBING_STYLE != 0)
static void s128_qform_multiply_and_square(s128_qform_group_t* group,
					   s128_qform_t* R,
					   const s128_qform_t* A) {
  s128_qform_t tmp;
  s128_qform_square(group, &tmp, A);
  s128_qform_compose(group, R, &tmp, A);
}
#endif

#if (S128_QFORM_CUBING_STYLE == 2)
static void s128_qform_dynamic_cube(s128_qform_group_t* group,
				    s128_qform_t* R,
				    const s128_qform_t* A) {
  int k = numbits_s128(&group->D);
  if (k <= CUBING_CROSSOVER) s128_qform_genuine_cube(group, R, A);
  else s128_qform_multiply_and_square(group, R, A);
}
#endif

void s128_qform_cube(s128_qform_group_t* group,
		     s128_qform_t* R,
		     const s128_qform_t* A) {
#if (S128_QFORM_CUBING_STYLE == 0)
  s128_qform_genuine_cube(group, R, A);
#elif (S128_QFORM_CUBING_STYLE == 1)
  s128_qform_multiply_and_square(group, R, A);
#elif (S128_QFORM_CUBING_STYLE == 2)
  s128_qform_dynamic_cube(group, R, A);
#else
#error "Unrecognized cubing style."
#endif
}


void s128_qform_print(s128_qform_group_t* group, const s128_qform_t* form) {
  char cbuffer[41];
  to_decstr_s128(cbuffer, 40, &form->c);
  printf("Qfb(%"PRId64", %"PRId64", %s)", form->a, form->b, cbuffer);
}
