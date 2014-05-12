
uint32_t norm(uint32_t tid, __global mpzcl_t *restrict mpz){
	int32_t i,j;
	for(i=mpz->size[tid]-1; i>=0; i--){
		uint64_t v = D_REF(mpz->d,i,tid);
		if(v!=0){
//			printf("%16.16lx\n",v);
			for(j=63; j>= 0; j--)
				if(v&(1ULL<<j))
					return (j + i * 64);
		}
	}
	return 0;
}

uint32_t norm32(uint32_t v){
	int32_t i;
	uint32_t ret=0;
	for(i=31; i>=0; i--){
		if(v & (1<<i))
			return i;
	}
	return 0;
}

uint64_t mul_hi(uint64_t a, uint64_t b){
	return ((__uint128_t)a * b) >> 64;
}

#define umul_ppmm(h,l,m,n) \
	h = mul_hi(m,n); \
	l = m * n;

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
    prod->size[tid] = m1->size[tid] + (cl != 0);
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

void setBN(uint32_t tid, __global mpzcl_t *restrict mpzPrimorial, uint64_t v){
	D_REF(mpzPrimorial->d,0,tid) = v;
	mpzPrimorial->size[tid] = 1;
}

void copyBN(uint32_t tid, __global mpzcl_t *restrict d, __global mpzcl_t *restrict s){
	uint32_t sz = s->size[tid];
	d->size[tid] = sz;
	uint32_t i;
	for(i=0; i < sz; i++){
		D_REF(d->d,i,tid) = D_REF(s->d,i,tid);
	}
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

#define MAX(x,y) ((x>y)?x:y)
#define MIN(x,y) ((x>y)?y:x)

void xorBNBN(uint32_t tid, __global mpzcl_t *restrict mpzD, __global mpzcl_t *restrict mpzA, __global mpzcl_t *restrict mpzB){
	uint32_t i;
	uint32_t asz = mpzA->size[tid];
	uint32_t bsz = mpzB->size[tid];

	for(i=0; i < MAX(asz,bsz); i++){
		uint64_t a = (i < asz) ? D_REF(mpzA->d,i,tid) : 0;
		uint64_t b = (i < bsz) ? D_REF(mpzB->d,i,tid) : 0;
		D_REF(mpzD->d,i,tid) =  a ^ b;
	}
	mpzD->size[tid] = MAX(asz,bsz);
}

void andBNBN(uint32_t tid, __global mpzcl_t *restrict mpzD, __global mpzcl_t *restrict mpzA, __global mpzcl_t *restrict mpzB){
	uint32_t i;
	uint32_t asz = mpzA->size[tid];
	uint32_t bsz = mpzB->size[tid];

	for(i=0; i < MIN(asz,bsz); i++){
		D_REF(mpzD->d,i,tid) = D_REF(mpzA->d,i,tid) & D_REF(mpzB->d,i,tid);
	}
	mpzD->size[tid] = MIN(asz,bsz);
}

void addBNBN (uint32_t tid, __global mpzcl_t *restrict sum, __global mpzcl_t *restrict a, __global mpzcl_t *restrict b){
	uint64_t carry = 0;
	uint32_t i = 0;
	uint32_t asz = a->size[tid];
	uint32_t bsz = b->size[tid];
	for (; i < MAX(asz,bsz); i++) {
		uint64_t asum = (i < asz)?D_REF(a->d,i,tid):0;
		uint64_t bsum = (i < bsz)?D_REF(b->d,i,tid):0;
		uint64_t accum = asum + bsum;
		uint32_t cy1 = accum < asum;
		uint64_t accum2 = accum + carry;
		carry = cy1 || (accum2 < accum);	

		D_REF(sum->d,i,tid) = accum2;;
	}
	D_REF(sum->d,MAX(asz,bsz),tid) = carry;
	sum->size[tid] = MAX(asz,bsz) + (carry?1:0);
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
	
	if(D_REF(d->d,i,tid)==0ULL)
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

void newton_invert(uint32_t tid, __global mpzcl_t *restrict mpzInv, __global mpzcl_t *restrict mpzDen, __global fermatTemp_t *temp){
    	lshiftBN(tid,&temp->mpzNewtonDen,mpzDen,7);
	//printbn(tid,&temp->mpzNewtonDen);
    	uint32_t e = norm(tid,&temp->mpzNewtonDen);
	//printf("e: %d\n", e);
	setBN(tid,&temp->mpzNewtonX1,0xB4B4B4B4);
	//printf("X1:\n");
	//printbn(tid,&temp->mpzNewtonX1);
	lshiftBN(tid,&temp->mpzNewtonX1s,&temp->mpzNewtonX1,e-(29));
	//printf("X1s:\n");
	//printbn(tid,&temp->mpzNewtonX1s);
	mulBN(tid,&temp->mpzNewtonX2,&temp->mpzNewtonDen,0xF0F0F0F0);
	//printf("X2:\n");
	//printbn(tid,&temp->mpzNewtonX2);
	rshiftBN(tid,&temp->mpzNewtonX2s,&temp->mpzNewtonX2,31);

	subBNBN(tid,&temp->mpzNewtonX,&temp->mpzNewtonX1s,&temp->mpzNewtonX2s);
	//printf("X:: \n");
	//printbn(tid,&temp->mpzNewtonX);
//    x = (0xB4B4B4B4 << (e-29)) - ((0xF0F0F0F0 * den_) >> 31)
    	uint32_t p = 2*e+2;
    	p = 1+norm32((p+1)/norm32(17));
	//p=10;
	uint32_t i;

	//printf("p: %d\n", p);
	for(i=0; i < p; i++){
		mulBNBN(tid,(mpzcl_t*)&temp->mpzNewtonDenP,&temp->mpzNewtonDen,&temp->mpzNewtonX);
		//printf("DenP:: \n");
		//printbn(tid,(mpzcl_t*)&temp->mpzNewtonDenP);
		rshiftBN(tid,&temp->mpzNewtonDenPS,(mpzcl_t*)&temp->mpzNewtonDenP,e+1);
		//printf("DenPS:: \n");
		//printbn(tid,(mpzcl_t*)&temp->mpzNewtonDenPS);
		setBN(tid,&temp->mpzNewton2,2);
		lshiftBN(tid,&temp->mpzNewton2s,&temp->mpzNewton2,e+1);
		//printf("2S:: \n");
		//printbn(tid,(mpzcl_t*)&temp->mpzNewton2s);
		subBNBN(tid,&temp->mpzNewtonDiff,&temp->mpzNewton2s,&temp->mpzNewtonDenPS);
		mulBNBN(tid,(mpzcl_t*)&temp->mpzNewtonProd,&temp->mpzNewtonDiff,&temp->mpzNewtonX);
		rshiftBN(tid,&temp->mpzNewtonX,(mpzcl_t*)&temp->mpzNewtonProd,e+1);	
		//printf("X:: \n");
		//printbn(tid,&temp->mpzNewtonX);	
//        x = x*((2<<(e+1))-((den_*x)>>(e+1)))>>(e+1)
	}
	rshiftBN(tid,mpzInv,&temp->mpzNewtonX,7);
//    return x>>2 */

}

void xbinGCD(uint32_t tid, __global mpzcl_t* mpzA, __global mpzcl_t* mpzB, __global fermatTemp_t *temp){
    	setBN(tid,&temp->mpzXbinU,1);
    	setBN(tid,&temp->mpzXbinV,0);
	uint32_t k = norm(tid,mpzA);
	printf("k: %d\n", k);
	uint32_t i;
	for(i=0; i <= k; i++){
		if((D_REF(temp->mpzXbinU.d,0,tid) & 1) == 0){
			rshiftBN(tid,&temp->mpzXbinTemp,&temp->mpzXbinU,1);
			//printf("zero:\n");
			//printbn(tid,&temp->mpzXbinTemp);

			copyBN(tid,&temp->mpzXbinU,&temp->mpzXbinTemp);
			//printbn(tid,&temp->mpzXbinU);

			rshiftBN(tid,&temp->mpzXbinTemp,&temp->mpzXbinV,1);
			copyBN(tid,&temp->mpzXbinV,&temp->mpzXbinTemp);


			//printbn(tid,&temp->mpzXbinV);
		}else{
            		//u = ((u ^ beta) >> 1) + (u & beta)
			//printf("Not zero:%d\n",i);
			//printbn(tid,&temp->mpzXbinU);
			//printbn(tid,mpzB);
			//printbn(tid,mpzA);

			xorBNBN(tid,&temp->mpzXbinTemp,&temp->mpzXbinU,mpzB);
			rshiftBN(tid,&temp->mpzXbinTemp2,&temp->mpzXbinTemp,1); 
			//printbn(tid,&temp->mpzXbinTemp);


			andBNBN(tid,&temp->mpzXbinTemp3,&temp->mpzXbinU,mpzB);
			addBNBN(tid,&temp->mpzXbinU,&temp->mpzXbinTemp2,&temp->mpzXbinTemp3);
			//printbn(tid,&temp->mpzXbinU);

            		//v = (v>>1) + alpha
			rshiftBN(tid,&temp->mpzXbinTemp,&temp->mpzXbinV,1);
			addBNBN(tid,&temp->mpzXbinV,&temp->mpzXbinTemp,mpzA);

			//printf("Not zero:%d\n",i);
			//printbn(tid,mpzA);
			//printbn(tid,&temp->mpzXbinV);

		}
	}
}

void normalizeBN(uint32_t tid, __global mpzcl_t *d){
	uint32_t sz = d->size[tid];
	while(sz && D_REF(d->d,sz-1,tid)==0){
		sz--;
	}
	d->size[tid] = sz;
}

bool cmpBNBNGTE(uint32_t tid, __global mpzcl_t *restrict a, __global mpzcl_t *restrict b, __global fermatTemp_t *temp){
	int32_t i;
	uint32_t asz = a->size[tid];
	uint32_t bsz = b->size[tid];
	for(i = MAX(asz,bsz) - 1; i>=0; i--){

		uint64_t al = D_REF(a->d,i,tid);
		if(i >= asz)
			al = 0;
		uint64_t bl = D_REF(b->d,i,tid);
		if(i >= bsz)
			bl = 0;

//		temp->al[tid] = al; temp->bl[tid] = bl;
//		temp->il[tid] = i;

		if(al>bl){
			return true;
		}
		if(bl>al)
			return false;
	}
	return true;
}

void wiki_barrett(uint32_t tid,
		__global mpzcl_t *restrict r, 
		__global mpzcl_t *restrict a, 
		__global mpzcl_t *restrict n, 
		__global mpzcl_t *restrict m, 
		__global fermatTemp_t *restrict temp) {
	__global mpzcl_t *restrict t1 = (mpzcl_t*)&temp->mpzBarrettT1;
	__global mpzcl_t *restrict t2 = (mpzcl_t*)&temp->mpzBarrettT2;
	__global mpzcl_t *restrict t3 = (mpzcl_t*)&temp->mpzBarrettT3;
	__global mpzcl_t *restrict t4 = (mpzcl_t*)&temp->mpzBarrettT4;

    	uint32_t k = norm(tid,n)+1;
	//copyBN(tid,&temp->mpzBarrettN,n);
	//copyBN(tid,&temp->mpzBarrettM,m);
	//copyBN(tid,&temp->mpzBarrettA,a);

    	mulBNBN(tid,t1,m,a);
     	rshiftBN(tid,t2,t1,2*k);

 	mulBNBN(tid,t3,t2,n);
	
	//t3->size[tid]++;
	printf("T3:\n");
	printbn(tid,t3);

 	subBNBN(tid,t4,a,t3);
	normalizeBN(tid,t4);	

	if(cmpBNBNGTE(tid,t4,n,temp)){
		subBNBN(tid,r,t4,n);
	}else
		copyBN(tid,r,t4);  
}

void barrett_square(uint32_t tid, __global mpzcl_t* mpzR, __global mpzcl_t* mpzH, __global mpzcl_t* mpzM, __global mpzcl_t* mpzInv, __global fermatTemp_t *restrict temp){
	//printf("H:\n");	
	//printbn(tid,mpzH);
	wiki_barrett(tid,&temp->mpzSqMod,mpzH,mpzM,mpzInv,temp);
	printf("Sqmod:\n");
	printbn(tid,&temp->mpzSqMod);
	mulBNBN(tid,(mpzcl_t*)&temp->mpzSq,&temp->mpzSqMod,&temp->mpzSqMod);
	//printbn(tid,(mpzcl_t*)&temp->mpzSq);
	wiki_barrett(tid,mpzR,(mpzcl_t*)&temp->mpzSq,mpzM,mpzInv,temp);
	printf("Barrett:\n");
	printbn(tid,mpzR);
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

void monPro(uint32_t tid, __global mpzcl32_t* mpzR, __global mpzcl32_t* mpzA, 
		__global mpzcl32_t* mpzB, uint32_t k, __global mpzcl32_t* mpzM, 
		uint32_t m0, __global fermatTemp_t *temp){
    	uint32_t s = k/32;
	uint32_t i,j;
	mpzcl32_t *mpzT = &temp->mpzMonProT;
	for(i=0; i < s+2; i++){
		D_REF(mpzT->d,i,tid) = 0;
	}
    
	for(i=0; i < s; i++){
		uint32_t c = 0;
		uint64_t cs;
		for(j=0; j < s; j++){
	    		cs = D_REF(mpzT->d,j,tid) + ((uint64_t)D_REF(mpzA->d,j,tid) * D_REF(mpzB->d,i,tid)) + c;
	    		c = cs >> 32;
	    		D_REF(mpzT->d,j,tid) = cs & 0xFFFFFFFF;
		}
		cs = D_REF(mpzT->d,s,tid) + (uint64_t)c;
		D_REF(mpzT->d,s,tid) =  cs & 0xFFFFFFFF;
		D_REF(mpzT->d,s+1,tid) = cs >> 32;

		uint32_t mf = (D_REF(mpzT->d,0,tid) * m0) & 0xFFFFFFFF;
		cs = D_REF(mpzT->d,0,tid) + (uint64_t)mf * D_REF(mpzM->d,0,tid);
		c = cs >> 32;

		for(j=1; j < s; j++) {
			cs = D_REF(mpzT->d,j,tid) + (uint64_t)mf*D_REF(mpzM->d,j,tid) + c;
			c = cs >> 32;
			D_REF(mpzT->d,j-1,tid) = cs & 0xFFFFFFFF;
		}
		cs = D_REF(mpzT->d,s,tid) + (uint64_t)c;
		D_REF(mpzT->d,s-1,tid) =  cs & 0xFFFFFFFF;
		D_REF(mpzT->d,s,tid) = D_REF(mpzT->d,s+1,tid) + (cs >> 32);
	}
	mpzT->size[tid] = s;
	copyBN32(tid,mpzR,mpzT);
}

void powmBN(uint32_t tid, __global fermatTemp_t *restrict temp){
	newton_invert(tid,&temp->mpzInv,&temp->mpzM,temp);

   	uint32_t maxbit = norm(tid,&temp->mpzM);
    uint32_t k = maxbit + 32 - (maxbit%32);

	k=352;
	clearBN(tid,&temp->mpzOne,(k+63)/64);
	clearBN(tid,&temp->mpzHalfR,(k+63)/64);
	clearBN(tid,&temp->mpzR,(k+63)/64);

//////////////
	D_REF(temp->mpzM.d,temp->mpzM.size[tid],tid)=0;
//////////////


	setBN(tid,&temp->mpzOne,1);
    	//r = 1 << k
	lshiftBN(tid,&temp->mpzR,&temp->mpzOne,k);
	rshiftBN(tid,&temp->mpzHalfR,&temp->mpzR,1);
	//printbn(tid,&temp->mpzHalfR);

   	//(r_, m_) = xbinGCD(r>>1, m)
	xbinGCD(tid,&temp->mpzHalfR,&temp->mpzM,temp);
	
	clearBN32(tid,&temp->mpzResult,k/32);
	clearBN32(tid,&temp->mpzBase,k/32);
	clearBN32(tid,&temp->mpzM32,k/32);
	clearBN32(tid,&temp->mpzH32,k/32);
	clearBN32(tid,&temp->mpzV32,k/32);

/////////////
	printbn(tid,&temp->mpzXbinU);
	printbn(tid,&temp->mpzXbinV);

	printbn(tid,&temp->mpzInv);

	//printbn(tid,&temp->mpzR);

/////////////

	//h = barrett_square(r,m,minv)
	//printbn(tid,&temp->mpzR);
	barrett_square(tid,&temp->mpzH,&temp->mpzR,&temp->mpzM,&temp->mpzInv,temp);

	printbn(tid,&temp->mpzR);
	printbn(tid,&temp->mpzM);
	printbn(tid,&temp->mpzInv);
	//printbn(tid,&temp->mpzH);

	toBN32(tid,&temp->mpzH32,&temp->mpzH);
	toBN32(tid,&temp->mpzM32,&temp->mpzM);
	toBN32(tid,&temp->mpzV32,&temp->mpzXbinV);

	setBN32(tid,&temp->mpzResult,1);
	setBN32(tid,&temp->mpzBase,2);

	uint32_t mi0 = D_REF(temp->mpzV32.d,0,tid);

	//To montgomery space
	monPro(tid,&temp->mpzResult,&temp->mpzResult,&temp->mpzH32,k,&temp->mpzM32,mi0,temp);
	monPro(tid,&temp->mpzBase,&temp->mpzBase,&temp->mpzH32,k,&temp->mpzM32,mi0,temp);

	printf("Montgomery space:\n");

    printbn32(tid,&temp->mpzH32);
    printbn32(tid,&temp->mpzResult);
	printf("Base:\n");
    printbn32(tid,&temp->mpzBase);

	//return;

	//printbn(tid,&temp->mpzE);

	uint32_t i;
	for(i=0; i <= maxbit; i++){
		if(bitIsSet(tid,&temp->mpzE,i)){
			//printf("b: %d\n", i);
			//printbn32(tid,&temp->mpzResult);
			//printbn32(tid,&temp->mpzBase);
			monPro(tid,&temp->mpzResult,&temp->mpzResult,&temp->mpzBase,k,&temp->mpzM32,mi0,temp);
			//printbn32(tid,&temp->mpzResult);
		}
		monPro(tid,&temp->mpzBase,&temp->mpzBase,&temp->mpzBase,k,&temp->mpzM32,mi0,temp);
	printf("Base:\n");
    printbn32(tid,&temp->mpzBase);
	}

	//From montgomery space
	clearBN32(tid,&temp->mpzOne32,k/32);
	setBN32(tid,&temp->mpzOne32,1);
	monPro(tid,&temp->mpzResult,&temp->mpzResult,&temp->mpzOne32,k,&temp->mpzM32,mi0,temp);
	fromBN32(tid,&temp->mpzR,&temp->mpzResult);
}
