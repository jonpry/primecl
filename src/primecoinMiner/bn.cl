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

void setBN(uint32_t tid, __global mpzcl_t *restrict mpzPrimorial, uint64_t v){
	D_REF(mpzPrimorial->d,0,tid) = v;
	mpzPrimorial->size[tid] = 1;
}

uint32_t count_trailing_zeros(uint32_t d){
	return 31-clz(d&-d);
}


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

#define SUBC(cout, w, x, y)        \
  do {                                  \
    uint64_t  __x = (x);               \
    uint64_t  __y = (y);               \
    uint64_t  __w = __x - __y;         \
    (w) = __w;                          \
    (cout) = __w > __x;                 \
  } while (0)

#define umul_ppmm(h,l,m,n) \
	h = mul_hi(m,n); \
	l = m * n;

uint32_t BNModExactOdd(uint32_t tid, __global mpzcl_t *restrict mpz, uint32_t d){
  	uint64_t inverse = invert64(d);
	uint32_t i;
	uint64_t s,c=0, h=0, l, dummy;

//	bn[0][tid] = 105*41+1;
  	for (i = 0; i < mpz->size[tid]; i++)
    	{
      		s = D_REF(mpz->d,i,tid);
      		SUBC (c, l, s, c);

		//c = l;
      		l = l * inverse;
		//c=l;

      		umul_ppmm (h, dummy, l, (uint64_t)d);
      		c += h;
    	}
 
	return c;
}

#define LOW_ZEROS_MASK(n)  (((n) & -(n)) - 1)

//debugged
uint32_t BN256DivisibleBy(uint32_t tid, __global mpzcl_t *restrict mpz, uint32_t d){
	uint32_t twos; 
	/*  asize = SIZ(a);
  	if (UNLIKELY (d == 0))
    	return (asize == 0);

  	if (asize == 0)  // 0 divisible by any d 
    		return 1; */

  	if (! (d & 1))
    	{
     	 	/* Strip low zero bits to get odd d required by modexact.  If d==e*2^n
         		and a is divisible by 2^n and by e, then it's divisible by d. */

      		if ((D_REF(mpz->d,0,tid) & LOW_ZEROS_MASK (d)) != 0)
        		return 0;

      		twos = count_trailing_zeros (d);
      		d >>= twos;
    	}

	return BNModExactOdd(tid,mpz, d) == 0;
//	return BNModExactOdd(tid,bn, d);

}

uint64_t mpn_addmul_1g (uint32_t tid, __global uint64_t *restrict sum, uint32_t sofst, __global uint64_t *restrict x, uint64_t xsz, uint64_t a){
	uint64_t carry=0;
	uint32_t i;
	uint64_t ul,lpl,hpl,rl;

	for(i=0; i < xsz; i++){
		
      		ul = D_REF(x,i,tid);
      		umul_ppmm (hpl, lpl, ul, a);

      		lpl += carry;
      		carry = (lpl < carry) + hpl;

      		rl = D_REF(sum,i+sofst,tid);
      		lpl = rl + lpl;
      		carry += lpl < rl;
      		D_REF(sum,i+sofst,tid) = lpl;
    	}

  	return carry;
}

//Debugged
void mulBN(uint32_t tid, __global mpzcl_t *restrict prod, __global mpzcl_t *restrict m1, uint64_t sml){
  uint64_t ul, cl, hpl, lpl;
  uint32_t i;
  cl = 0;
  for(i=0; i < m1->size[tid]; i++) {
      ul = D_REF(m1->d,i,tid);
      umul_ppmm (hpl, lpl, ul, sml);

      lpl += cl;
      cl = (lpl < cl) + hpl;

      D_REF(prod->d,i,tid) = lpl;
    }

    D_REF(prod->d,m1->size[tid],tid) = cl;
    prod->size[tid] = m1->size[tid] + (cl != (uint64_t)0);
}

void mulBN2(uint32_t dtid, uint32_t stid, __global mpzcl_t *restrict prod, __global mpzcl_t *restrict m1, uint64_t sml){
  uint64_t ul, cl, hpl, lpl;
  uint32_t i;
  cl = 0;
  for(i=0; i < m1->size[stid]; i++) {
      ul = D_REF(m1->d,i,stid);
      umul_ppmm (hpl, lpl, ul, sml);

      lpl += cl;
      cl = (lpl < cl) + hpl;

      D_REF(prod->d,i,dtid) = lpl;
    }

    D_REF(prod->d,m1->size[stid],dtid) = cl;
    prod->size[dtid] = m1->size[stid] + (cl != 0);
}

void addBN (uint64_t tid, __global mpzcl_t *restrict sum, __global mpzcl_t *restrict a, uint64_t b){
	uint64_t carry = b;
	uint32_t i = 0;
	for (; i < a->size[tid]; i++) {
		uint64_t psum = D_REF(sum->d,i,tid);
		psum = D_REF(a->d,i,tid) + carry;
		D_REF(sum->d,i,tid) = psum;
		carry = psum - carry > psum;
	}
	D_REF(sum->d,a->size[tid],tid) = carry;
	sum->size[tid] = a->size[tid] + (carry?1:0);
}

void subBN (uint64_t tid, __global mpzcl_t *restrict diff, __global mpzcl_t *restrict a, uint64_t b){
	uint64_t borrow = b;
	uint32_t i = 0;
	for (; i < a->size[tid]; i++) {
		uint64_t pdiff = D_REF(a->d,i,tid) - borrow;
		D_REF(diff->d,i,tid) = pdiff;
		borrow = pdiff + borrow < pdiff;
	}
	uint32_t sz = a->size[tid];
	//if(D_REF(diff->d,a->size[tid]-1,tid) == (uint64_t)0)
	//	sz--;
	diff->size[tid] = sz; //No support negative
}

void subBNBN(uint32_t tid, __global mpzcl_t *restrict diff, __global mpzcl_t *restrict a, __global mpzcl_t *restrict b){
  	uint64_t cy = 0, al, bl, sl, rl;
	uint32_t cy1, cy2;
  	uint32_t i=0;
	uint32_t asz = a->size[tid];
	uint32_t bsz = b->size[tid];
  	for(i=0; i < asz; i++)
    	{
      		al = D_REF(a->d,i,tid);
      		bl = i < bsz ? D_REF(b->d,i,tid) : 0;
      		sl = al - bl;
      		cy1 = sl > al;
      		rl = sl - cy;
      		cy2 = rl > sl;
      		cy = cy1 | cy2;
      		D_REF(diff->d,i,tid) = rl;
    	}
	diff->size[tid] = asz;
}
	 
//Not debugged 
void mulBNBNS(uint32_t tid, __global mpzcl_t *restrict rp, __global mpzcl_t *restrict up, __global mpzcls_t *restrict vp)
{
	uint32_t vofst=1,rofst=1;
	mulBN(tid,rp,up,vp->d[0]);   
 	//rp += 1, vp += 1, vn -= 1;
  	while (vofst < vp->size) {
		//clear high word //TODO: right 
	    	D_REF(rp->d,rp->size[tid]+1,tid) = 0;

            	D_REF(rp->d,up->size[tid]+rofst,tid) = mpn_addmul_1g (tid, rp->d,rofst , up->d, up->size[tid], vp->d[vofst]);
	    	vofst++; rofst++;
	    	rp->size[tid]++;
        }
}

void mulBNBN(uint32_t tid, __global mpzcl_t *restrict rp, __global mpzcl_t *restrict up, __global mpzcl_t *restrict vp)
{
	if(up->size[tid] < vp->size[tid]){
		__global mpzcl_t* temp = up;
		up = vp;
		vp = temp;
	}

	uint32_t vofst=1,rofst=1;
	mulBN(tid,rp,up,D_REF(vp->d,0,tid));   
 	//rp += 1, vp += 1, vn -= 1;
	   // 	rp->size[tid]++;

  	while (vofst < vp->size[tid]) {
		//clear high word //TODO: right 
	//	printf("Size: %d\n", rp->size[tid]);
	    	D_REF(rp->d,rp->size[tid]+1,tid) = 0;

            	D_REF(rp->d,up->size[tid]+rofst,tid) = mpn_addmul_1g (tid, rp->d,rofst , up->d, up->size[tid], D_REF(vp->d,vofst,tid));
	    	vofst++; rofst++;
	    	rp->size[tid]++;
        }

	if(D_REF(rp->d,up->size[tid] + vp->size[tid] - 1,tid) != (uint64_t)0)
		rp->size[tid]++;
}

void rshiftBN(uint32_t tid, __global mpzcl_t *restrict d, __global mpzcl_t *restrict src, uint32_t shift){
	uint32_t i;
	uint32_t ofst = shift/64;
	//printf("rshift: %d\n" ,shift);
	shift = shift % 64;
	//printf("rshift: %d\n" ,shift);
	int32_t sz = src->size[tid] - ofst;

	//If there is nothing left, can't run next code
	if(sz<=0){
		d->size[tid]=0;
		D_REF(d->d,0,tid)=0;
		return;
	}
	for(i=0; i < sz-1; i++){
		D_REF(d->d,i,tid) = (D_REF(src->d,i+ofst,tid) >> shift) | ((shift==0)?0:(D_REF(src->d,i+ofst+1,tid) << (64 - shift))); 
	}
	D_REF(d->d,i,tid) = (D_REF(src->d, i+ofst,tid) >> shift);
	//printf("%16.16LX\n", D_REF(d->d,i,tid));
	
	if(D_REF(d->d,i,tid)==(uint64_t)0)
		sz--;

	d->size[tid] = sz;
} 

void lshiftBN(uint32_t tid, __global mpzcl_t *restrict d, __global mpzcl_t *restrict src, uint32_t shift){
	uint32_t i;

	//printf("shift %d\n", shift);
	uint32_t ofst = shift/64;
	shift = shift % 64;

	//printf("shift %d\n", shift);

	for(i=0; i < ofst; i++)
		D_REF(d->d,i,tid) = 0;

	D_REF(d->d,i,tid) = D_REF(src->d,0,tid) << shift; 

	for(i=1; i < src->size[tid]; i++){
		D_REF(d->d,i+ofst,tid) = (D_REF(src->d,i-1,tid) >> (64-shift)) | (D_REF(src->d,i,tid) << shift); 
	}
	uint64_t v = D_REF(src->d, i-1,tid) >> (64-shift);
	D_REF(d->d,i+ofst,tid) = shift!=0?v:0;
	//printf("src size: %d, ofst: %d\n", src->size[tid], ofst);
	d->size[tid] = src->size[tid] + ofst+ 1;//(shift != 0 && v != 0)?1:0;
} 

void copyBN(uint32_t tid, __global mpzcl_t *restrict d, __global mpzcl_t *restrict s){
	uint32_t sz = s->size[tid];
	d->size[tid] = sz;
	uint32_t i;
	for(i=0; i < sz; i++){
		D_REF(d->d,i,tid) = D_REF(s->d,i,tid);
	}
}

void normalizeBN(uint32_t tid, __global mpzcl_t *d){
	uint32_t sz = d->size[tid];
	while(sz && D_REF(d->d,sz-1,tid)==0){
		sz--;
	}
	d->size[tid] = sz;
}

void xorBNBN(uint32_t tid, __global mpzcl_t *restrict mpzD, __global mpzcl_t *restrict mpzA, __global mpzcl_t *restrict mpzB){
	uint32_t i;
	uint32_t asz = mpzA->size[tid];
	uint32_t bsz = mpzB->size[tid];

	for(i=0; i < max(asz,bsz); i++){
		uint64_t a = (i < asz) ? D_REF(mpzA->d,i,tid) : 0;
		uint64_t b = (i < bsz) ? D_REF(mpzB->d,i,tid) : 0;
		D_REF(mpzD->d,i,tid) =  a ^ b;
	}
	mpzD->size[tid] = max(asz,bsz);
}

void andBNBN(uint32_t tid, __global mpzcl_t *restrict mpzD, __global mpzcl_t *restrict mpzA, __global mpzcl_t *restrict mpzB){
	uint32_t i;
	uint32_t asz = mpzA->size[tid];
	uint32_t bsz = mpzB->size[tid];

	for(i=0; i < min(asz,bsz); i++){
		D_REF(mpzD->d,i,tid) = D_REF(mpzA->d,i,tid) & D_REF(mpzB->d,i,tid);
	}
	mpzD->size[tid] = min(asz,bsz);
}

void addBNBN (uint32_t tid, __global mpzcl_t *restrict sum, __global mpzcl_t *restrict a, __global mpzcl_t *restrict b){
	uint64_t carry = 0;
	uint32_t i = 0;
	uint32_t asz = a->size[tid];
	uint32_t bsz = b->size[tid];
	for (; i < max(asz,bsz); i++) {
		uint64_t asum = (i < asz)?D_REF(a->d,i,tid):0;
		uint64_t bsum = (i < bsz)?D_REF(b->d,i,tid):0;
		uint64_t accum = asum + bsum;
		uint32_t cy1 = accum < asum;
		uint64_t accum2 = accum + carry;
		carry = cy1 || (accum2 < accum);	

		D_REF(sum->d,i,tid) = accum2;;
	}
	D_REF(sum->d,max(asz,bsz),tid) = carry;
	sum->size[tid] = max(asz,bsz) + (carry?1:0);
}

void toBN32(uint32_t tid, __global mpzcl32_t* mpzD, __global mpzcl_t* mpzS){
	uint32_t i;
	uint32_t len=0;
	for(i=0; i < mpzS->size[tid]; i++){
		uint64_t v = D_REF(mpzS->d,i,tid);
		D_REF(mpzD->d,i*2,tid) = v;
		D_REF(mpzD->d,i*2+1,tid) = v >> 32;
		if(v>>32 != 0)
			len = 2*(i+1);
		else if(v!=0)
			len = 2*i + 1;
	}
	mpzD->size[tid] = len;
}

void fromBN32(uint32_t tid, __global mpzcl_t* mpzD, __global mpzcl32_t* mpzS){
	uint32_t i;
	uint32_t vlast;
	for(i=0; i < mpzS->size[tid]; i++){
		uint32_t v = D_REF(mpzS->d,i,tid);
		if(i%2==0)
			vlast = v;
		else
			D_REF(mpzD->d,i/2,tid) = (uint64_t)v << 32 | vlast;
	}
	if(i%2==1)
		D_REF(mpzD->d,i/2,tid) = vlast;
	mpzD->size[tid] = i/2 + ((i%2)?1:0);
}

void setBN32(uint32_t tid, __global mpzcl32_t* mpzD, uint32_t v){
	mpzD->size[tid] = 1;
	D_REF(mpzD->d,0,tid) = v;
}

void copyBN32(uint32_t tid, __global mpzcl32_t *mpzD, __global mpzcl32_t *mpzS){
	uint32_t i;
	for(i=0; i < mpzS->size[tid]; i++){
		D_REF(mpzD->d,i,tid) = D_REF(mpzS->d,i,tid);
	}
	mpzD->size[tid] = mpzS->size[tid];
}

void copyBN2(uint32_t dtid, __global mpzcl_t *mpzD, uint32_t stid, __global mpzcl32_t *mpzS){
	uint32_t i;
	for(i=0; i < mpzS->size[stid]; i++){
		D_REF(mpzD->d,i,dtid) = D_REF(mpzS->d,i,stid);
	}
	mpzD->size[dtid] = mpzS->size[stid];
}

void clearBN32(uint32_t tid, __global mpzcl32_t *mpzD, uint32_t s){
	uint32_t i;
	for(i=0; i < s; i++){
		D_REF(mpzD->d,i,tid) = 0;
	}
}

void clearBN(uint32_t tid, __global mpzcl_t *mpzD, uint32_t s){
	uint32_t i;
	for(i=0; i < s; i++){
		D_REF(mpzD->d,i,tid) = 0;
	}
}

uint32_t bitIsSet(uint32_t tid, __global mpzcl_t *mpz, uint32_t bit){
	uint32_t word = bit / 64;
	bit = bit % 64;
	return (D_REF(mpz->d,word,tid)>>bit)&1;
}

uint32_t cmpBNEQ(uint32_t tid, __global mpzcl_t *restrict mpz, uint64_t v){
	uint32_t i;
	for(i=mpz->size[tid]-1; i>0; i--){
		if(D_REF(mpz->d,i,tid))
			return 0;
	}
	return (D_REF(mpz->d,0,tid)==v)?1:0;
}
