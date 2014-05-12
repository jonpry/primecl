/*OpenCL primecoin miner
Copyright (C) 2014 Jon Pry <jonpry@gmail.com>

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.*/

void mpn_rshift(uint32_t tid, __private uint64_t *restrict d, __global uint64_t *restrict src, uint32_t ofst, uint32_t length, uint32_t shift){
	uint32_t i;
//	printf("Shift %d length: %d, ofst %d\n", shift, length, ofst);
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

//	assert(dp[0] & 1 != 0);

//printf("nsize: %d, dsize: %d, dinv %lX\n", nsize, dsize, dinv);

  	for (i = nsize - dsize; i > 0; i--)
    	{
//		printf("Tick foo1111111\n");
      		q = dinv * np[nofst];
      		qp[qofst++] = ~q;
      		cy = mpn_addmul_1 (&np[nofst], dp, dsize, q);
      		mpn_add_1 (&np[dsize+nofst], &np[dsize+nofst], i, cy);
//		assert (np[nofst] == 0);
      		nofst++;
    	}


//printbn(np,8);
  	for (i = dsize; i > 1; i--)
    	{
//s		printf("Tick2\n");
      		q = dinv * np[nofst];
      		qp[qofst++] = ~q;
      		mpn_addmul_1 (&np[nofst], dp, i, q);
//		assert (np[nofst] == 0);
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
//		printf("nsize < dsize\n");
		return;
	}	


//  	if (quot == num || quot == den)
//    		qp = TMP_ALLOC_LIMBS (qn);

	uint32_t ofst=0;
  	while (D_REF(d->d,ofst,tid) == 0 && ofst < d->size[tid])
    	{
//		printf("trimming zero\n");
      		ofst++;
      		nsize--;
      		dsize--;
    	}

	uint32_t qsize = nsize - dsize + 1;

	uint32_t shift = count_trailing_zeros (D_REF(d->d,ofst,tid));
//	printf("shift: %d\n", shift);
    	{
#if 1
      		uint32_t ss = (dsize > qsize) ? qsize + 1 : dsize;

//		printf("qsize: %d nsize: %d, dsize %d\n", qsize, nsize, dsize);
      		mpn_rshift (tid,np, n->d, ofst, qsize + 1, shift);
      		mpn_rshift (tid,dp, d->d, ofst, ss, shift);

#else
      		mpn_rshift (tid,np, n->d, ofst, nsize, shift);
      		mpn_rshift (tid,dp, d->d, ofst, dsize, shift);
#endif
    	}

//	printf("After rshift d: ");
//	printbn(dp,dsize);
//	printf("\n");

  	di = invert64(dp[0]);  
	di = -di;

  	mpn_bdiv (qp, np,nsize, dp, dsize, di);
	uint32_t foo = qsize;
//	printf("q size: %d\n", foo);
	int i;
	int sz=-1;
	for(i=0; i <qsize; i++){
		D_REF(q->d,i,tid) = qp[i];
		if(qp[i])
			sz=i;
	}
	q->size[tid] = sz+1;
}

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

uint64_t invert(uint64_t d){
	uint64_t dummy;
   	uint64_t dinv;//TODO: invert d;
   	uint64_t negone = -1;
	__udiv_qrnnd_c(dinv,dummy,~d,negone,d);
	return dinv;
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

		{ 
			uint64_t __q_lo = t_hi * dinv;
			uint64_t __q_hi = mul_hi(t_hi,dinv);
			if(__q_lo + t_lo < __q_lo)
				__q_hi++;
			__q_lo+=t_lo;
			__q_hi+=t_hi;

			uint64_t __q1 = __q_hi+1; 
			uint64_t __r1 = t_lo - __q1*(d); 
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

//TODO: does not set size on quotient
uint64_t modBN(uint32_t tid, __global mpzcl_t *restrict q, __global mpzcl_t *restrict a, uint64_t d){
   	uint32_t norm = clz(d);
	d <<= norm;

	uint64_t dinv = invert(d);
 
	return divrem1_c(tid,q->d,a->d,a->size,d,dinv,norm);
} 

