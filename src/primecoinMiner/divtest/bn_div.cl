
#define __private
#define __global
#define __constant
#define restrict 

#include "assert.h"

uint64_t mul_hi(uint64_t u, uint64_t v)
{
 	return  ((__int128_t)u * v) >> 64;
}

uint64_t mul_hi32(uint32_t u, uint32_t v)
{
 	return  ((uint64_t)u * v) >> 32;
}

uint64_t mul_lo(uint64_t u, uint64_t v){
	return u * v;
}

#define umul_ppmm(h,l,m,n) \
	h = mul_hi(m,n); \
	l = m * n;

__constant uint32_t  invert_table[128] = {
  0x01, 0xAB, 0xCD, 0xB7, 0x39, 0xA3, 0xC5, 0xEF,
  0xF1, 0x1B, 0x3D, 0xA7, 0x29, 0x13, 0x35, 0xDF,
  0xE1, 0x8B, 0xAD, 0x97, 0x19, 0x83, 0xA5, 0xCF,
  0xD1, 0xFB, 0x1D, 0x87, 0x09, 0xF3, 0x15, 0xBF,
  0xC1, 0x6B, 0x8D, 0x77, 0xF9, 0x63, 0x85, 0xAF,
  0xB1, 0xDB, 0xFD, 0x67, 0xE9, 0xD3, 0xF5, 0x9F,
  0xA1, 0x4B, 0x6D, 0x57, 0xD9, 0x43, 0x65, 0x8F,
  0x91, 0xBB, 0xDD, 0x47, 0xC9, 0xB3, 0xD5, 0x7F,
  0x81, 0x2B, 0x4D, 0x37, 0xB9, 0x23, 0x45, 0x6F,
  0x71, 0x9B, 0xBD, 0x27, 0xA9, 0x93, 0xB5, 0x5F,
  0x61, 0x0B, 0x2D, 0x17, 0x99, 0x03, 0x25, 0x4F,
  0x51, 0x7B, 0x9D, 0x07, 0x89, 0x73, 0x95, 0x3F,
  0x41, 0xEB, 0x0D, 0xF7, 0x79, 0xE3, 0x05, 0x2F,
  0x31, 0x5B, 0x7D, 0xE7, 0x69, 0x53, 0x75, 0x1F,
  0x21, 0xCB, 0xED, 0xD7, 0x59, 0xC3, 0xE5, 0x0F,
  0x11, 0x3B, 0x5D, 0xC7, 0x49, 0x33, 0x55, 0xFF
};

//Requires n be odd
uint64_t invert64(uint64_t n){
    	uint64_t inv = invert_table[(n/2) & 0x7F]; /*  8 */		
	inv = 2 * inv - inv * inv * n;
	inv = 2 * inv - inv * inv * n;
	inv = 2 * inv - inv * inv * n;
//	inv = 2 * inv - inv * inv * n;

	return inv;
}

#define clz(x) __builtin_clz(x)

uint32_t count_trailing_zeros(uint32_t d){
	uint32_t ret = 31-clz(d&-d);
	printf("clz of %x is %d\n", d, ret);
	return ret;
}



void mpn_rshift(uint32_t tid, __private uint64_t *restrict d, __global uint64_t *restrict src, uint32_t ofst, uint32_t length, uint32_t shift){
	uint32_t i;
	printf("Shift %d length: %d, ofst %d\n", shift, length, ofst);
	for(i=0; i < length-1; i++){
//		printf("Tick\n");
		d[i] = (D_REF(src,i+ofst,tid) >> shift) | (D_REF(src,i+ofst+1,tid) << (64 - shift)); 
	}
	d[i] = (D_REF(src, i+ofst,tid) >> shift);
} 

void mpn_copy(uint32_t tid, __private uint64_t *restrict d, __global uint64_t *restrict src, uint32_t ofst, uint32_t length){
	uint32_t i;
	for(i=0; i < length; i++){
		d[i] = D_REF(src,i+ofst,tid); 
	}
}

uint64_t mpn_addmul_1 (__private uint64_t *restrict sum, __private uint64_t *restrict x, uint64_t xsz, uint64_t a){
	uint64_t carry=0;
	uint32_t i;
	uint64_t ul,lpl,hpl,rl;

	for(i=0; i < xsz; i++){
		
      		ul = x[i];
      		umul_ppmm (hpl, lpl, ul, a);

      		lpl += carry;
      		carry = (lpl < carry) + hpl;

      		rl = sum[i];
      		lpl = rl + lpl;
      		carry += lpl < rl;
      		sum[i] = lpl;
    	}

  	return carry;
}


uint64_t mpn_add (__private uint64_t *restrict sum, __private uint64_t *restrict a, uint64_t asz,  __private uint64_t *restrict b, uint64_t bsz){
	uint64_t carry = 0;
	uint32_t i = 0;
	uint64_t len = asz < bsz ? bsz : asz;
	for(; i < len; ++i) {
		uint64_t a_i = i < asz ? a[i] : 0;
		uint64_t b_i = i < bsz ? b[i] : 0;
		sum[i] = a_i + b_i + carry;
		carry = sum[i] < carry ? 1 : 0;
	}
	return carry;
}

uint64_t mpn_add_1 (__private uint64_t *restrict sum, __private uint64_t *restrict a, uint64_t sz, uint64_t b){
	uint64_t carry = b;
	uint32_t i = 0;
	for (; i < sz; i++) {
		sum[i] = a[i] + carry;
		carry = sum[i] < carry;
	}
	return carry;
}

uint64_t mpn_add_n (__private uint64_t *restrict sum, __private uint64_t *restrict a, __private uint64_t *restrict b, uint64_t sz){
	uint64_t carry = 0;
	uint32_t i = 0;
	for (; i < sz; i++) {
		sum[i] = a[i] + b[i] + carry;
		carry = sum[i] < carry ? 1 : 0;
	}
	return carry;
}

uint64_t mpn_sub_n (__private uint64_t *restrict diff, __private uint64_t *restrict a, __private uint64_t *restrict b, uint64_t sz){
	uint64_t carry = 0;
	uint32_t i = 0;
	for (; i < sz; ++i) {
		diff[i] = a[i] - b[i] - carry;
		carry = (diff[i] < a[i]) | (diff[i] - carry > diff[i])? 1 : 0;
	}
	return carry;
}

void mpn_bdiv (__private uint64_t *restrict qp, 
		   __private uint64_t *restrict np, uint32_t nsize,
		   __private uint64_t *restrict dp, uint32_t dsize, uint64_t dinv)
{
//  	uint32_t qsize;
  	uint32_t i;
	uint32_t qofst=0,nofst=0;
	uint64_t q,cy;
//  	uint64_t rh;
//  	uint64_t ql;

	assert(dp[0] & 1 != 0);

//printf("nsize: %d, dsize: %d, dinv %lX\n", nsize, dsize, dinv);

  	for (i = nsize - dsize; i > 0; i--)
    	{
		printf("Tick foo1111111\n");
      		q = dinv * np[nofst];
      		qp[qofst++] = ~q;
      		cy = mpn_addmul_1 (&np[nofst], dp, dsize, q);
      		mpn_add_1 (&np[dsize+nofst], &np[dsize+nofst], i, cy);
		assert (np[nofst] == 0);
      		nofst++;
    	}


//printbn(np,8);
  	for (i = dsize; i > 1; i--)
    	{
//s		printf("Tick2\n");
      		q = dinv * np[nofst];
      		qp[qofst++] = ~q;
      		mpn_addmul_1 (&np[nofst], dp, i, q);
		assert (np[nofst] == 0);
      		nofst++;
//printbn(np,8);
    	}

  	/* Final limb */
  	q = dinv * np[nofst];
  	qp[qofst] = ~q;

//	printbn(qp,8);
  	mpn_add_1 (&qp[qofst - nsize + 1], &qp[qofst - nsize + 1], nsize, 1);
//	printbn(qp,8);
//	printbn(np,8);
}

//TODO: Not debugged
void divBN(uint32_t tid, __global mpzcl_t *restrict q, __global mpzcl_t *restrict n, __global mpzcl_t *restrict d){

	uint32_t nsize = n->size[tid];
	uint32_t dsize = d->size[tid];

	uint64_t qp[8],np[8],dp[8];
	uint64_t di;

	//Catch this early to avoid crash
	if(nsize < dsize){
		q->size[tid] = 0;
		D_REF(q->d,0,tid) = 0;
		printf("nsize < dsize\n");
		return;
	}	


//  	if (quot == num || quot == den)
//    		qp = TMP_ALLOC_LIMBS (qn);

	uint32_t ofst=0;
  	while (D_REF(d->d,ofst,tid) == 0 && ofst < d->size[tid])
    	{
		printf("trimming zero\n");
      		ofst++;
      		nsize--;
      		dsize--;
    	}

	uint32_t qsize = nsize - dsize + 1;

	uint32_t shift = count_trailing_zeros (D_REF(d->d,ofst,tid));
	printf("shift: %d\n", shift);
    	{
#if 1
      		uint32_t ss = (dsize > qsize) ? qsize + 1 : dsize;

		printf("qsize: %d nsize: %d, dsize %d\n", qsize, nsize, dsize);
      		mpn_rshift (tid,np, n->d, ofst, qsize + 1, shift);
      		mpn_rshift (tid,dp, d->d, ofst, ss, shift);

#else
      		mpn_rshift (tid,np, n->d, ofst, nsize, shift);
      		mpn_rshift (tid,dp, d->d, ofst, dsize, shift);
#endif
    	}

	printf("After rshift d: ");
	printbn(dp,dsize);
	printf("\n");

  	di = invert64(dp[0]);  
	di = -di;

  	mpn_bdiv (qp, np,nsize, dp, dsize, di);
	uint32_t foo = qsize;
	printf("q size: %d\n", foo);
	int i;
	int sz=-1;
	for(i=0; i <qsize; i++){
		D_REF(q->d,i,tid) = qp[i];
		if(qp[i])
			sz=i;
	}
	q->size[tid] = sz+1;
}

