/**
 * @file math64.h
 * Fast 64 bit arithmetic functions.
 */
#pragma once
#ifndef MATH64__INCLUDED
#define MATH64__INCLUDED

#include <assert.h>
#include <stdint.h>

#include "math32.h"
#include "s128_t.h"

#ifndef INT32_MIN
#define INT32_MIN              (-2147483647-1)
#endif

#ifndef INT32_MAX
#define INT32_MAX              (2147483647)
#endif

/// Random value [0,2^64-1]
static inline uint64_t rand_u64(void) {
  uint64_t res;
  res = rand_u32();
  res <<= 32;
  res |= rand_u32();
  return res;
}

/// Compute a - b and let m = -1 if a < b and 0 otherwise.
static inline int64_t sub_with_mask_u64(uint64_t* m,
					const uint64_t a,
					const uint64_t b) {
  //#if defined(__x86_64)
  //  int64_t r;
  //  asm("subq %3, %0\n\t"
  //      "sbbq %1, %1\n\t"  // %1 is either 0 or -1
  //      : "=r"(r), "=&r"(*m)
  //      : "0"(a), "r"(b)
  //      : "cc");
  //  return r;
  //#else
  //*m = a < b ? -1 : 0;
  *m = (a >= b) - 1;  // a < b ? -1 : 0;
  return a - b;
  //#endif
}

/// Compute a - b and let m = -1 if a < b and 0 otherwise.
static inline int64_t sub_with_mask_s64(uint64_t* m,
					const int64_t a,
					const int64_t b) {
#if defined(__x86_64)
  int64_t r;
  asm("subq %3, %0\n\t"
      "sbbq %1, %1\n\t"  // %1 is either 0 or -1
      : "=r"(r), "=&r"(*m)
      : "0"(a), "r"(b)
      : "cc");
  return r;
#else
  *m = a < b ? -1 : 0;
  return a - b;
#endif
}

static inline void swap_s64(int64_t* u, int64_t* v) {
  *u ^= *v;
  *v ^= *u;
  *u ^= *v;
}

/// Conditionally swap u with v if u < v.
static inline void cond_swap_s64(int64_t* u, int64_t* v) {
  uint64_t m;
  int64_t d = sub_with_mask_s64(&m, *u, *v);
  d &= m;
  *u -= d;
  *v += d;
}

/// Conditionally swap u with v if u2 < v2.
static inline void cond_swap2_s64(int64_t* u1, int64_t* u2,
				  int64_t* v1, int64_t* v2) {
  uint64_t m;
  int64_t d2 = sub_with_mask_s64(&m, *u2, *v2);
  int64_t d1 = (*u1 - *v1) & m;
  d2 &= m;
  *u1 -= d1;
  *u2 -= d2;
  *v1 += d1;
  *v2 += d2;
}

/// Conditionally swap u with v if u3 < v3.
static inline uint64_t cond_swap3_s64(int64_t* u1,
				      int64_t* u2,
				      int64_t* u3,
				      int64_t* v1,
				      int64_t* v2,
				      int64_t* v3) {
  uint64_t m;
  int64_t d3 = sub_with_mask_s64(&m, *u3, *v3);
  int64_t d1 = (*u1 - *v1) & m;
  int64_t d2 = (*u2 - *v2) & m;
  d3 &= m;
  *u1 -= d1;
  *u2 -= d2;
  *u3 -= d3;
  *v1 += d1;
  *v2 += d2;
  *v3 += d3;
  return m;
}

/// Conditionally swap u with v if u3 < v3.
static inline void cond_swap3_s64_mixed(int64_t*  u1,
					int64_t*  u2,
					uint64_t* u3,
					int64_t*  v1,
					int64_t*  v2,
					uint64_t* v3) {
  uint64_t m;
  int64_t d3 = sub_with_mask_u64(&m, *u3, *v3);
  int64_t d1 = (*u1 - *v1) & m;
  int64_t d2 = (*u2 - *v2) & m;
  d3 &= m;
  *u1 -= d1;
  *u2 -= d2;
  *u3 -= d3;
  *v1 += d1;
  *v2 += d2;
  *v3 += d3;
}

/// Negate using a mask. m must be either -1 or 0.
static inline int64_t negate_using_mask_s64(const uint64_t m,
					    const int64_t x) {
  assert(m == 0 || m == (uint64_t)(-1));
  return (x ^ m) - m;
}

/// Negate x when c < 0.
static inline int64_t cond_negate_s64(const int64_t c,
				      const int64_t x) {
  return negate_using_mask_s64(c >> 63, x);
}

/// Absolute.
static inline uint64_t abs_s64(const int64_t x) {
  return cond_negate_s64(x, x);
}

/// true if x fits in a signed 32-bit integer
static inline int s64_is_s32(int64_t x) {
  return (x >= INT32_MIN && x <= INT32_MAX);
}

/// most significant bit
static inline int msb_u64(uint64_t x) {
#if defined(__x86_64)
  int64_t k = -1;
  asm("bsrq %1, %0\n\t"
      : "=r"(k)
      : "r"(x), "0"(k)
      : "cc");
  return k;
#else
  // a binary search approach to finding the most significant set bit
  int n = 0;
  if (x == 0) return -1;
  if (x > 0xFFFFFFFFULL) { n += 32; x >>= 32; }
  if (x > 0xFFFF) { n += 16; x >>= 16; }
  if (x > 0xFF) { n += 8;  x >>= 8; }
  if (x > 0xF) { n += 4;  x >>= 4; }
  if (x > 0x7) { n += 2;  x >>= 2; }
  if (x > 0x3) { n += 1;  x >>= 1; }
  if (x > 0x1) { n ++; }
  return n;
#endif
}

/// least significant set bit
static inline int lsb_u64(uint64_t x) {
#if defined(__x86_64)
  int64_t k = -1;
  asm("bsfq %1, %0\n\t"
      : "=r"(k)
      : "r"(x), "0"(k)
      : "cc");
  return k;
#else
  // a binary search approach to finding the least significant set bit
  int k = 1;
  if (x == 0) return -1;
  if ((x & 0xFFFFFFFF) == 0) { k += 32;  x >>= 32; }
  if ((x & 0xFFFF) == 0) { k += 16;  x >>= 16; }
  if ((x & 0xFF) == 0) { k += 8;  x >>= 8; }
  if ((x & 0xF) == 0) { k += 4;  x >>= 4; }
  if ((x & 0x7) == 0) { k += 2;  x >>= 2; }
  if ((x & 0x3) == 0) { k += 1;  x >>= 1; }
  k -= x&1;
  return k;
#endif
}

/// least significant set bit
static inline int lsb_s64(int64_t x) {
#if defined(__x86_64)
  int64_t k = -1;
  asm("bsfq %1, %0\n\t"
      : "=r"(k)
      : "r"(x), "0"(k)
      : "cc");
  return k;
#else
  // a binary search approach to finding the least significant set bit
  int k = 1;
  if (x == 0) return -1;
  if ((x & 0xFFFFFFFF) == 0) { k += 32;  x >>= 32; }
  if ((x & 0xFFFF) == 0) { k += 16;  x >>= 16; }
  if ((x & 0xFF) == 0) { k += 8;  x >>= 8; }
  if ((x & 0xF) == 0) { k += 4;  x >>= 4; }
  if ((x & 0x7) == 0) { k += 2;  x >>= 2; }
  if ((x & 0x3) == 0) { k += 1;  x >>= 1; }
  k -= x&1;
  return k;
#endif
}

/// Set the i^th bit.
static inline uint64_t setbit_u64(const uint64_t x, const int i) {
  return x | (1 << i);
}

/// Clear the i^th bit.
static inline uint64_t clrbit_u64(const uint64_t x, const int i) {
  return x & ~(1 << i);
}


/// the number of bits in x, i.e. the smallest k such that 2^k > x
static inline int numbits_s64(int64_t x) {
  return msb_u64((uint64_t)abs_s64(x)) + 1;
}

/// res = n % m
static inline uint32_t mod_u32_u64_u32(const uint64_t in_n,
				       const uint32_t in_m) {
  return in_n % in_m;
  // We used to have this fancy branching technique, but in
  // practice, it's faster letting the compiler do it.
  /*
#if defined(__x86_64)
  if (in_n < in_m) {
    // zero divides
    return in_n;
  }

  uint32_t r;
  const uint32_t nhi = in_n >> 32;
  const uint32_t nlo = in_n & 0xFFFFFFFF;
  if (in_m > nhi) {
    // 32-bit divide, since result will fit.
    asm("movl %1, %%edx\n\t"
	"movl %2, %%eax\n\t"
	"divl %3\n\t"
	: "=&d"(r)
	: "rm"(nhi), "rm"(nlo), "rm"(in_m)
	: "cc", "eax");
  } else {
    // 64-bit divide.
    asm("xorq %%rdx, %%rdx\n\t"
	"divq %2\n\t"
	: "=&d"(r)
	: "a"(in_n), "r"((uint64_t)in_m)
	: "cc");
  }
  return r;
#else
  return in_n % in_m;
#endif
  */
}

/// Computes the signed remainder (n % m) that is closest to 0.
static inline int32_t mod_s32_s64_u32(const int64_t in_n,
				      const uint32_t in_m) {
  uint64_t m1 = in_n >> 63;
  uint32_t a1 = mod_u32_u64_u32(negate_using_mask_s64(m1, in_n), in_m);
  int32_t r1  = negate_using_mask_s32(m1, a1);
  int32_t r2  = r1 - negate_using_mask_s32(m1, in_m);
  uint32_t a2 = negate_using_mask_s32(~m1, r2);
  uint32_t m  = (a1 < a2) - 1;
  int32_t r   = (r1 & ~m) | (r2 & m);
  return r;
}

/// r = s1+s2 (mod m)
static inline int64_t addmod_s64(const int64_t s1,
				 const int64_t s2,
				 const int64_t m) {
  assert(m > 0);
  int64_t r;
#if defined(__x86_64)
  asm("movq %2, %%r11\n\t"
      "movq %1, %%rdx\n\t"
      "movq %2, %%r10\n\t"
      "movq %1, %%rax\n\t"
      "sarq $63, %%r11\n\t"
      "sarq $63, %%rdx\n\t"
      "addq %%r10, %%rax\n\t"
      "adcq %%r11, %%rdx\n\t"
      "idivq %3\n\t"
      : "=&d"(r)
      : "rm"(s1), "rm"(s2), "rm"(m)
      : "cc", "rax", "r10", "r11");
#else
  s128_t a;
  s128_t b;
  set_s128_s64(&a, s1);
  set_s128_s64(&b, s2);
  add_s128_s128(&a, &b);
  r = mod_s64_s128_s64(&a, m);
#endif
  return r;
}

/// r = s1-s2 (mod m)
static inline int64_t submod_s64(const int64_t s1,
				 const int64_t s2,
				 const int64_t m) {
  int64_t r;
#if defined(__x86_64)
  asm("movq %2, %%r11\n\t"
      "movq %1, %%rdx\n\t"
      "movq %2, %%r10\n\t"
      "movq %1, %%rax\n\t"
      "sarq $63, %%r11\n\t"
      "sarq $63, %%rdx\n\t"
      "subq %%r10, %%rax\n\t"
      "sbbq %%r11, %%rdx\n\t"
      "idivq %3\n\t"
      : "=&d"(r)
      : "rm"(s1), "rm"(s2), "rm"(m)
      : "cc", "rax", "r10", "r11");
#else
  s128_t a;
  s128_t b;
  set_s128_s64(&a, s1);
  set_s128_s64(&b, s2);
  sub_s128_s128(&a, &b);
  r = mod_s64_s128_s64(&a, m);
#endif
  return r;
}

/// n = qd+r
static inline void divrem_u64(uint64_t* out_q,
			      uint64_t* out_r,
			      const uint64_t in_n,
			      const uint64_t in_d) {
/*
#if defined(__x86_64)
  asm("movq %2, %%rax\n\t"
      "xorq %%rdx, %%rdx\n\t"
      "divq %3\n\t"
      : "=&a"(*out_q), "=&d"(*out_r)
      : "rm"(in_n), "rm"(in_d)
      : "cc");
#else
*/
  // the compiler appears to be smart enough
  // to optimize this into a single divide
  *out_q = in_n / in_d;
  *out_r = in_n % in_d;
//#endif
}

/// n = qd+r
static inline void divrem_s64(int64_t* out_q,
			      int64_t* out_r,
			      const int64_t in_n,
			      const int64_t in_d) {
/*
#if defined(__x86_64)
  asm("movq %2, %%rdx\n\t"
      "movq %2, %%rax\n\t"
      "sarq $63, %%rdx\n\t"
      "idivq %3\n\t"
      : "=&a"(*out_q), "=&d"(*out_r)
      : "rm"(in_n), "rm"(in_d)
      : "cc");
#else
*/
  // the compiler appears to be smart enough
  // to optimize this into a single divide
  *out_q = in_n / in_d;
  *out_r = in_n % in_d;
//#endif
}


/// res = (x*y) % m
static inline uint64_t mulmod_u64(const uint64_t x,
				  const uint64_t y,
				  const uint64_t m) {
#if defined(__x86_64)
  uint64_t res;
  uint64_t tmp;
  asm("movq %2, %%rax\n\t"
      "mulq %3\n\t"
      "cmpq %4, %%rdx\n\t"
      "jb 1f\n\t"
      "movq %%rax, %1\n\t"
      "movq %%rdx, %%rax\n\t"
      "xorq %%rdx, %%rdx\n\t"
      "divq %4\n\t"
      "movq %1, %%rax\n\t"
      "1:\n\t"
      "divq %4\n\t"
      : "=&d"(res), "=&r"(tmp)
      : "rm"(x), "rm"(y), "r"(m)
      : "cc", "rax");
  return res;
#else
  u128_t t;
  uint64_t r;
  mul_u128_u64_u64(&t, x, y);
  mod_u64_u128_u64(&r, &t, m);
  return r;
#endif
}

/// res = x*y (mod m)
static inline int64_t mulmod_s64(const int64_t x,
				 const int64_t y,
				 const int64_t m) {
  int64_t m2 = (int64_t)abs_s64(m);
  
  // Make sure x and y are positive and compute the result sign.
  // This is a non-branching trick for the following:
  //  if (x < 0) {
  //    s = 1;
  //    x2 = -x;
  //  }
  //  if (y < 0) {
  //    s = 1-s;
  //    y2 = -y;
  //  }
  int64_t xt = x >> 63;  // xt is either 0 or -1
  int64_t yt = y >> 63;
  int64_t x2 = (x^xt) - xt;  // negate x if xt==-1
  int64_t y2 = (y^yt) - yt;
  int64_t s = (xt^yt);  // s is either all 0s or all 1s
  
  // perform multiply with remainder
  int64_t r = (int64_t)mulmod_u64(x2, y2, m2);
  
  // use the remainder that is closest to 0
  uint64_t mask;
  sub_with_mask_s64(&mask, m2 >> 1, r);
  r -= m2 & mask;
  
  // Correct the sign of the remainder
  return (r^s) - s;  // negates r is s is -1.
}

/// res = f1*f2+f3*f4
static inline int64_t muladdmul_s64_4s32(const int32_t f1,
					 const int32_t f2,
					 const int32_t f3,
					 const int32_t f4) {
#if defined(__x86_64)
  int64_t res;
  asm("movl %1, %%eax\n\t"
      "imull %2\n\t"
      "movl %%eax, %%ebx\n\t"
      "movl %3, %%eax\n\t"
      "movl %%edx, %%ecx\n\t"
      "imull %4\n\t"
      "addl %%ebx, %%eax\n\t"
      "adcl %%ecx, %%edx\n\t"
      "shlq $32, %%rdx\n\t"
      "orq %%rdx, %%rax\n\t"
      : "=&a"(res)
      : "rm"(f1), "rm"(f2), "rm"(f3), "rm"(f4)
      : "rbx", "rcx", "rdx", "cc");
  return res;
#else
  return ((int64_t)f1 * (int64_t)f2) + ((int64_t)f3 * (int64_t)f4);
#endif 
}

// res = (f1*f2+f3*f4)/d
// NOTE: possible overflow if the result does not fit within 64bits signed
static inline int64_t muladdmuldiv_s64(const int64_t f1,
				       const int64_t f2,
				       const int64_t f3,
				       const int64_t f4,
				       const int64_t d) {
#if defined(__x86_64)
  int64_t res, t1, t2;
  asm("movq %3, %%rax\n\t"
      "imulq %4\n\t"
      "movq %%rax, %1\n\t"
      "movq %5, %%rax\n\t"
      "movq %%rdx, %2\n\t"
      "imulq %6\n\t"
      "addq %1, %%rax\n\t"
      "adcq %2, %%rdx\n\t"
      "idivq %7\n\t"
      : "=&a"(res), "=&r"(t1), "=&r"(t2)
      : "r"(f1), "r"(f2), "r"(f3), "r"(f4), "r"(d)
      : "rdx", "cc");
  return res;
#else
  s128_t t1;
  s128_t t2;
  s128_t q;
  mul_s128_s64_s64(&t1, f1, f2);
  mul_s128_s64_s64(&t2, f3, f4);
  add_s128_s128(&t1, &t2);
  div_s128_s128_s64(&q, &t1, d);
  return get_s64_from_s128(&q);
#endif
}

/// the largest s such that s^2 <= x
uint64_t sqrt_u64(const uint64_t x);

/// true if there exists an integer s such that s^2 == x
static inline int is_square_u64(const uint64_t x) {
  uint64_t s = sqrt_u64(x);
  return s * s == x;
}

/// true if there exists an integer s such that s^2 == x
static inline int is_square_s64(const int64_t x) {
  return is_square_u64(abs_s64(x));
}

/// compute a^e mod m using binary exponentiation
uint64_t expmod_u64(uint64_t a, uint64_t e, uint64_t m);

#endif  // MATH64__INCLUDED


