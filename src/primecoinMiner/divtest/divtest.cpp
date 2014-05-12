#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <gmpxx.h>

#include "cldefs.h"
void printbn(uint64_t *d, uint32_t sz);
#include "bn_div.cl"

void copy_mpz(mpz_t n,mpzcl_t *ret){
	int i;
	int size = n->_mp_size;
	printf("Size: %d\n", size);
	for(i=0; i < size; i++){
		ret->d[i] = n->_mp_d[i];
	}
	ret->size[0] = size;
}

void printbn(uint64_t *d, uint32_t sz){
	int i;
	for(i=0; i < sz; i++){
		printf("%16.16lX", d[sz-1-i]);	
	}
	printf("\n");
}

uint32_t one_count(uint64_t v){
	uint32_t ret=0,i;
	for(i=0; i < 64; i++){
		if(v & (1ULL<<i))
			ret++;
	}
	return ret;
}

#define sub_ddmmss(sh, sl, ah, al, bh, bl) \
  do {									\
    uint64_t __x;								\
    __x = (al) - (bl);							\
    (sh) = (ah) - (bh) - ((al) < (bl));                                 \
    (sl) = __x;								\
  } while (0)


#define add_ssaaaa(sh, sl, ah, al, bh, bl) \
  do {									\
    uint64_t __x;								\
    __x = (al) + (bl);							\
    (sh) = (ah) + (bh) + (__x < (al));					\
    (sl) = __x;								\
  } while (0)

#define W_TYPE_SIZE 64

#define __ll_B ((uint64_t) 1 << (W_TYPE_SIZE / 2))
#define __ll_lowpart(t) ((uint64_t) (t) & (__ll_B - 1))
#define __ll_highpart(t) ((uint64_t) (t) >> (W_TYPE_SIZE / 2))

/* Define this unconditionally, so it can be used for debugging.  */
#define __udiv_qrnnd_c(q, r, n1, n0, d) \
  do {									\
    uint64_t __d1, __d0, __q1, __q0, __r1, __r0, __m;			\
									\									
    __d1 = __ll_highpart (d);						\
    __d0 = __ll_lowpart (d);						\
									\
    __q1 = (n1) / __d1;							\
    __r1 = (n1) - __q1 * __d1;						\
    __m = __q1 * __d0;							\
    __r1 = __r1 * __ll_B | __ll_highpart (n0);				\
    if (__r1 < __m)							\
      {									\
	__q1--, __r1 += (d);						\
	if (__r1 >= (d)) 						\
	  if (__r1 < __m)						\
	    __q1--, __r1 += (d);					\
      }									\
    __r1 -= __m;							\
									\
    __q0 = __r1 / __d1;							\
    __r0 = __r1  - __q0 * __d1;						\
    __m = __q0 * __d0;							\
    __r0 = __r0 * __ll_B | __ll_lowpart (n0);				\
    if (__r0 < __m)							\
      {									\
	__q0--, __r0 += (d);						\
	if (__r0 >= (d))						\
	  if (__r0 < __m)						\
	    __q0--, __r0 += (d);					\
      }									\
    __r0 -= __m;							\
									\
    (q) = __q1 * __ll_B | __q0;						\
    (r) = __r0;								\
  } while (0)

#define WORD_BITS 64
typedef unsigned long word_t;
typedef unsigned int dword_t __attribute__((mode(TI)));

typedef unsigned int qword_t __attribute__((mode(QI)));

static inline
void precompute_inverse1(word_t *dinv, word_t *norm, word_t d)
{
   dword_t t;
   word_t normi = clz(d)+32;

   d <<= normi;
   t = (~(dword_t) 0) - (((dword_t) d) << WORD_BITS);
   
   *dinv = t / d;
   *norm = normi;
}

uint32_t mclz(uint64_t x){
	int32_t i;
	uint32_t ret=0;
	for(i=63; i>=0; i--){
		if(! (x & (1ULL<<i)))
			ret++;
		else
			break;
	}
	return ret;
}

uint64_t precompute_inverse2(uint64_t d){
	uint64_t q,r;
   	uint32_t norm = mclz(d);

//	printf("D: %X, Norm: %d\n",norm);
   	d <<= norm;

	__udiv_qrnnd_c(q,r,~d,(uint64_t)-1LL,d);
	return q;
#if 0
	uint32_t norm = clz(d)+32;//gcc clz is fubar
	d <<= norm;
	
	uint32_t d9 = d >> (64-9);
	uint32_t d32 = d >> (64-32);

	printf("Norm: %d, d9: %X, d32: %X\n", norm, d9, d32); 

	uint32_t v0 = (((1<<18)-(1<<8))/d9)<<6;

	uint32_t v1 = (v0 << 17) - (( (uint64_t)(v0*v0)*d32) << 1 );

	printf("v0: %X, v1: %X\n", v0, v1);

	uint64_t v2 = ((uint64_t)v1<<33) - (mul_hi((uint64_t)v1*v1,d)<<1);

	uint64_t v2m_lo = v2<<2;
	uint64_t v2m_hi = v2>>62;

	__uint128_t foo = (__uint128_t)v2*v2;
	foo >>= 64;
	foo *= d;
	foo >>= 62;

//	uint64_t v22_hi = mul_hi(v2,v2) >> 62;
//	uint64_t v22_lo = mul_lo(v2,d);

	uint64_t v3_lo = v2m_lo - ((uint64_t)(foo)) - 1;

	uint64_t v3 = v3_lo;

	printf("v2: %lx, v3: %lx\n", v2, v3);

	uint64_t v4 = v3 - mul_hi(v3 + 1, d);

	return v4;

#endif
}

/*
Given a double word u, a normalised divisor d and a precomputed
inverse dinv of d, computes the quotient and remainder of u by d.
*/
#define divrem21_preinv1(q, r, u, d, dinv) \
do { \
dword_t __q = ((u)>>WORD_BITS) * (dword_t) (dinv); \
printf("NT__q: %16.16lx%16.16lx\n", (uint64_t)(__q>>64), (uint64_t)__q); \
__q += u; \
printf("NT__q: %16.16lx%16.16lx\n", (uint64_t)(__q>>64), (uint64_t)__q); \
word_t __q1 = (word_t)(__q >> WORD_BITS) + 1; \
word_t __q0 = (word_t) __q; \
word_t __r1 = (word_t)(u) - __q1*(d); \
printf("NT__q1r1: %16.16lx%16.16lx\n", __q1, __r1);\
if (__r1 >= __q0) \
{ \
__q1--; \
__r1 += (d); \
} \
if (__r1 >= (d)) \
{ \
(q) = __q1 + 1; \
(r) = __r1 - (d); \
} else \
{ \
(q) = __q1; \
(r) = __r1; \
} \
} while (0)

word_t nn_divrem1_preinv_c(word_t *q, word_t *a, uint32_t m,
                            word_t d, uint64_t dinv, uint32_t norm)
{	
   word_t ci=0;
   dword_t t;
   long i;
   
   d <<= norm;

   for (i = m - 1; i >= 0; i--)
   {
      t = (((dword_t) ci) << WORD_BITS) + (((dword_t) a[i]) << norm);
	printf("NTt: %16.16lx%16.16lx\n", (uint64_t)(t>>64), (uint64_t)t);
       divrem21_preinv1(q[i], ci, t, d, dinv);
   }

   return (ci >> norm);
}

uint64_t divrem1_c(	uint32_t tid,
			__global uint64_t *restrict q, 
			__global uint64_t *restrict a, 
			__global uint32_t *restrict m,
                        uint64_t d, uint64_t dinv,
			uint32_t norm)
{
	uint64_t ci=0;
   	int32_t i;
	uint64_t t_lo, t_hi;

   	for (i = m[tid] - 1; i >= 0; i--)
   	{
		t_lo = D_REF(a,i,tid) << norm;
		t_hi = ci + (D_REF(a,i,tid) >> (64 - norm));


		printf("MEt: %16.16lx%16.16lx\n", t_hi, t_lo);
		{ 
			uint64_t __q_lo = t_hi * dinv;
			uint64_t __q_hi = mul_hi(t_hi,dinv);
		printf("ME__q: %16.16lx%16.16lx\n", __q_hi, __q_lo);
			if(__q_lo + t_lo < __q_lo)
				__q_hi++;
			__q_lo+=t_lo;
			__q_hi+=t_hi;
		printf("ME__q: %16.16lx%16.16lx\n", __q_hi, __q_lo);
			uint64_t __q1 = __q_hi+1; 
			uint64_t __r1 = t_lo - __q1*(d); 

		printf("ME__q1r1: %16.16lx%16.16lx\n", __q1, __r1);

			if (__r1 >= __q_lo) 
			{ 
				__q1--; 
				__r1 += d; 
			} 
			if (__r1 >= (d)) 
			{ 
				D_REF(q,i,tid) = __q1 + 1; 
				ci = __r1 - d; 
			} else { 
				D_REF(q,i,tid) = __q1; 
				ci = __r1; 
			} 
		}
   	} 
 
   	return (ci >> norm);
}


int main(){
  	mpz_t n,d; 
  	mpz_init(n);
  	mpz_init(d);
	
	mpzcl_t ncl, dcl, qcl;

	srand(time(0));
  	mpz_set_ui(n,1);
  	mpz_set_ui(d,1);

	int i;
	for(i=0; i < 16; i++){ 
		mpz_mul_ui(n,n,rand());
		mpz_mul_ui(d,d,rand());
	}

	copy_mpz(n,&ncl);
	copy_mpz(d,&dcl);

	printf("Starting n: %X\n");
  	mpz_out_str(stdout,16,n);

	printf("\nStarting d: %X\n");
  	mpz_out_str(stdout,16,d);

	mpz_divexact(n,n,d);

	printf("\nResult: \n");
	mpz_out_str(stdout,16,n);

	printf("\nExtracted n:");
	printbn(ncl.d,ncl.size[0]);

	printf("\nExtracted d:");
	printbn(dcl.d,dcl.size[0]);

/////////
#if 0
	for(i=0; i < 8; i++){
		ncl.d[i] = 0;
		dcl.d[i] = 0;	
	}

	ncl.size[0] = 2;
	dcl.size[0] = 1;
	ncl.d[1] = 0x333333330;
	ncl.d[0] = 0x0;
	dcl.d[0] = 3;
	dcl.d[1] = 0;
#endif
/////////

	divBN(0,&qcl,&ncl,&dcl);

	printf("CLQuot: ");
	printbn(qcl.d,qcl.size[0]);

	word_t inv;
	word_t norm;

	uint64_t thenum[2], q[2];
	thenum[1] = 3;
	thenum[0] = 0xfed0;

	for(i=0; i < 1; i++){
		uint64_t den = i;

		uint64_t inv2 = precompute_inverse2(den);
		uint32_t norm = mclz(den);
		uint64_t normed = den << norm;
		uint32_t size = 2;

		uint64_t r1 = nn_divrem1_preinv_c(q,thenum,2,den,inv2,norm);
		uint64_t r2 = divrem1_c(0,q,thenum,&size,normed,inv2,norm);

		printf("R1: %lx, R2: %lx\n", r1, r2);
	}	

	return 0;
}
