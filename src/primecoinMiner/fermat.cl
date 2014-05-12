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

bool cmpBNGT(uint32_t tid, __global mpzcl_t *restrict mpz, uint64_t l){
	int32_t i;
	for(i = mpz->size[tid]-1; i>0; i--){
		if(D_REF(mpz->d,i,tid))
			return true;
	}
	if(D_REF(mpz->d,0,tid) > l)
		return true;
	return false;
}

bool cmpBNBNGTE(uint32_t tid, __global mpzcl_t *restrict a, __global mpzcl_t *restrict b, __global fermatTemp_t *temp){
	int32_t i;
	uint32_t asz = a->size[tid];
	uint32_t bsz = b->size[tid];
	for(i = max(asz,bsz) - 1; i>=0; i--){

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

uint32_t norm(uint32_t tid, __global mpzcl_t *restrict mpz){
	int32_t i;
	for(i=mpz->size[tid]-1; i>=0; i--){
		uint32_t v = clz(D_REF(mpz->d,i,tid));
		if(v!=64)
			return (64 - v + i * 64)-1;
	}
	return 0;
}

uint32_t norm32(uint32_t v){
	return 32-clz(v)-1;
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
		mulBNBN(tid,(__global mpzcl_t*)&temp->mpzNewtonDenP,&temp->mpzNewtonDen,&temp->mpzNewtonX);
		//printf("DenP:: \n");
		//printbn(tid,(mpzcl_t*)&temp->mpzNewtonDenP);
		rshiftBN(tid,&temp->mpzNewtonDenPS,(__global mpzcl_t*)&temp->mpzNewtonDenP,e+1);
		//printf("DenPS:: \n");
		//printbn(tid,(mpzcl_t*)&temp->mpzNewtonDenPS);
		setBN(tid,&temp->mpzNewton2,2);
		lshiftBN(tid,&temp->mpzNewton2s,&temp->mpzNewton2,e+1);
		//printf("2S:: \n");
		//printbn(tid,(mpzcl_t*)&temp->mpzNewton2s);
		subBNBN(tid,&temp->mpzNewtonDiff,&temp->mpzNewton2s,&temp->mpzNewtonDenPS);
		mulBNBN(tid,(__global mpzcl_t*)&temp->mpzNewtonProd,&temp->mpzNewtonDiff,&temp->mpzNewtonX);
		rshiftBN(tid,&temp->mpzNewtonX,(__global mpzcl_t*)&temp->mpzNewtonProd,e+1);	
		//printf("X:: \n");
		//printbn(tid,&temp->mpzNewtonX);	
//        x = x*((2<<(e+1))-((den_*x)>>(e+1)))>>(e+1)
	}
	rshiftBN(tid,mpzInv,&temp->mpzNewtonX,7);
//    return x>>2 */

}

void xbinGCD(uint32_t tid, __global mpzcl_t* mpzA, __global mpzcl_t* mpzB, __global fermatTemp_t *temp){
	uint32_t k = 351;
	uint32_t i;
	uint64_t a_0 = D_REF(mpzA->d,0,tid);
	uint64_t a_1 = D_REF(mpzA->d,1,tid);
	uint64_t a_2 = D_REF(mpzA->d,2,tid);
	uint64_t a_3 = D_REF(mpzA->d,3,tid);
	uint64_t a_4 = D_REF(mpzA->d,4,tid);
	uint64_t a_5 = D_REF(mpzA->d,5,tid);
	uint64_t b_0 = D_REF(mpzB->d,0,tid);
	uint64_t b_1 = D_REF(mpzB->d,1,tid);
	uint64_t b_2 = D_REF(mpzB->d,2,tid);
	uint64_t b_3 = D_REF(mpzB->d,3,tid);
	uint64_t b_4 = D_REF(mpzB->d,4,tid);
	uint64_t b_5 = D_REF(mpzB->d,5,tid);
	uint64_t v_0 = 0;
	uint64_t v_1 = 0;
	uint64_t v_2 = 0;
	uint64_t v_3 = 0;
	uint64_t v_4 = 0;
	uint64_t v_5 = 0;
	uint64_t u_0 = 1;
	uint64_t u_1 = 0;
	uint64_t u_2 = 0;
	uint64_t u_3 = 0;
	uint64_t u_4 = 0;
	uint64_t u_5 = 0;
	for(i=0; i <= k; i++) {
		uint32_t odd=(u_0 & 1)?0:1;
		if(odd){
			{ //Rshift BEGIN
				u_0 = (u_0 >> 1) | (u_1 << 63);
				u_1 = (u_1 >> 1) | (u_2 << 63);
				u_2 = (u_2 >> 1) | (u_3 << 63);
				u_3 = (u_3 >> 1) | (u_4 << 63);
				u_4 = (u_4 >> 1) | (u_5 << 63);
				u_5 = (u_5 >> 1);
			} //Rshift END
		}else{
			{ //And Xor Add BEGIN
				 uint32_t cy=0;
				{
					uint64_t a1 = u_0 ^ b_0;
					uint64_t a2 = u_1 ^ b_1;
					uint64_t a = (a1 >> 1) | (a2 << 63);
					uint64_t b = u_0 & b_0;
					uint64_t sum = a + b;
					uint32_t c0 = (sum < a)?1:0;
					u_0 = sum + cy;
					if(u_0 < sum)
						c0=1;
					cy = c0;
				}
				{
					uint64_t a1 = u_1 ^ b_1;
					uint64_t a2 = u_2 ^ b_2;
					uint64_t a = (a1 >> 1) | (a2 << 63);
					uint64_t b = u_1 & b_1;
					uint64_t sum = a + b;
					uint32_t c0 = (sum < a)?1:0;
					u_1 = sum + cy;
					if(u_1 < sum)
						c0=1;
					cy = c0;
				}
				{
					uint64_t a1 = u_2 ^ b_2;
					uint64_t a2 = u_3 ^ b_3;
					uint64_t a = (a1 >> 1) | (a2 << 63);
					uint64_t b = u_2 & b_2;
					uint64_t sum = a + b;
					uint32_t c0 = (sum < a)?1:0;
					u_2 = sum + cy;
					if(u_2 < sum)
						c0=1;
					cy = c0;
				}
				{
					uint64_t a1 = u_3 ^ b_3;
					uint64_t a2 = u_4 ^ b_4;
					uint64_t a = (a1 >> 1) | (a2 << 63);
					uint64_t b = u_3 & b_3;
					uint64_t sum = a + b;
					uint32_t c0 = (sum < a)?1:0;
					u_3 = sum + cy;
					if(u_3 < sum)
						c0=1;
					cy = c0;
				}
				{
					uint64_t a1 = u_4 ^ b_4;
					uint64_t a2 = u_5 ^ b_5;
					uint64_t a = (a1 >> 1) | (a2 << 63);
					uint64_t b = u_4 & b_4;
					uint64_t sum = a + b;
					uint32_t c0 = (sum < a)?1:0;
					u_4 = sum + cy;
					if(u_4 < sum)
						c0=1;
					cy = c0;
				}
				{
					uint64_t a1 = u_5 ^ b_5;
					uint64_t a = (a1 >> 1);
					uint64_t b = u_5 & b_5;
					uint64_t sum = a + b;
					uint32_t c0 = (sum < a)?1:0;
					u_5 = sum + cy;
					if(u_5 < sum)
						c0=1;
					cy = c0;
				}
			} //And Xor Add END
		}
			{ //Rshift BEGIN
				v_0 = (v_0 >> 1) | (v_1 << 63);
				v_1 = (v_1 >> 1) | (v_2 << 63);
				v_2 = (v_2 >> 1) | (v_3 << 63);
				v_3 = (v_3 >> 1) | (v_4 << 63);
				v_4 = (v_4 >> 1) | (v_5 << 63);
				v_5 = (v_5 >> 1);
			} //Rshift END
			if(!odd)
			{ //Add BEGIN
				 uint32_t cy=0;
				{
					uint64_t sum = v_0 + a_0;
					uint32_t c0 = (sum < v_0)?1:0;
					v_0 = sum + cy;
					if(v_0 < sum)
						c0=1;
					cy = c0;
				}
				{
					uint64_t sum = v_1 + a_1;
					uint32_t c0 = (sum < v_1)?1:0;
					v_1 = sum + cy;
					if(v_1 < sum)
						c0=1;
					cy = c0;
				}
				{
					uint64_t sum = v_2 + a_2;
					uint32_t c0 = (sum < v_2)?1:0;
					v_2 = sum + cy;
					if(v_2 < sum)
						c0=1;
					cy = c0;
				}
				{
					uint64_t sum = v_3 + a_3;
					uint32_t c0 = (sum < v_3)?1:0;
					v_3 = sum + cy;
					if(v_3 < sum)
						c0=1;
					cy = c0;
				}
				{
					uint64_t sum = v_4 + a_4;
					uint32_t c0 = (sum < v_4)?1:0;
					v_4 = sum + cy;
					if(v_4 < sum)
						c0=1;
					cy = c0;
				}
				{
					uint64_t sum = v_5 + a_5;
					uint32_t c0 = (sum < v_5)?1:0;
					v_5 = sum + cy;
					if(v_5 < sum)
						c0=1;
					cy = c0;
				}
			} //Add END
	}
	D_REF(temp->mpzXbinV.d,0,tid) = v_0;
	D_REF(temp->mpzXbinV.d,1,tid) = v_1;
	D_REF(temp->mpzXbinV.d,2,tid) = v_2;
	D_REF(temp->mpzXbinV.d,3,tid) = v_3;
	D_REF(temp->mpzXbinV.d,4,tid) = v_4;
	D_REF(temp->mpzXbinV.d,5,tid) = v_5;
	D_REF(temp->mpzXbinU.d,0,tid) = u_0;
	D_REF(temp->mpzXbinU.d,1,tid) = u_1;
	D_REF(temp->mpzXbinU.d,2,tid) = u_2;
	D_REF(temp->mpzXbinU.d,3,tid) = u_3;
	D_REF(temp->mpzXbinU.d,4,tid) = u_4;
	D_REF(temp->mpzXbinU.d,5,tid) = u_5;
	temp->mpzXbinV.size[tid] = 6;
	temp->mpzXbinU.size[tid] = 6;
}






void wiki_barrett(uint32_t tid,
		__global mpzcl_t *restrict r, 
		__global mpzcl_t *restrict a, 
		__global mpzcl_t *restrict n, 
		__global mpzcl_t *restrict m, 
		__global fermatTemp_t *restrict temp) {
	__global mpzcl_t *restrict t1 = (__global mpzcl_t*)&temp->mpzBarrettT1;
	__global mpzcl_t *restrict t2 = (__global mpzcl_t*)&temp->mpzBarrettT2;
	__global mpzcl_t *restrict t3 = (__global mpzcl_t*)&temp->mpzBarrettT3;
	__global mpzcl_t *restrict t4 = (__global mpzcl_t*)&temp->mpzBarrettT4;

    	uint32_t k = norm(tid,n)+1;
	//copyBN(tid,&temp->mpzBarrettN,n);
	//copyBN(tid,&temp->mpzBarrettM,m);
	//copyBN(tid,&temp->mpzBarrettA,a);

    	mulBNBN(tid,t1,m,a);
     	rshiftBN(tid,t2,t1,2*k);

 	mulBNBN(tid,t3,t2,n);

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
	//printf("Sqmod:\n");
	//printbn(tid,&temp->mpzSqMod);
	mulBNBN(tid,(__global mpzcl_t*)&temp->mpzSq,&temp->mpzSqMod,&temp->mpzSqMod);
	//printbn(tid,(mpzcl_t*)&temp->mpzSq);
	wiki_barrett(tid,mpzR,(__global mpzcl_t*)&temp->mpzSq,mpzM,mpzInv,temp);
	//printbn(tid,mpzR);
}

void monPro(uint32_t tid, __global mpzcl32_t* mpzR, __global mpzcl32_t* mpzA, 
		__global mpzcl32_t* mpzB, uint32_t k, __global mpzcl32_t* mpzM, 
		uint32_t m0, __global fermatTemp_t *temp){
    	uint32_t size = k/32;
	uint32_t i,j;
	__global mpzcl32_t *mpzT = &temp->mpzMonProT;
	for(i=0; i < size+2; i++){
		D_REF(mpzT->d,i,tid) = 0;
	}

	uint32_t c, c2,mf;
	uint32_t s,x,y,z,s2;
    
	for(i=0; i < size; i++){
		c = 0;
		for(j=0; j < size; j++){
			x = D_REF(mpzA->d,j,tid);
			y = D_REF(mpzB->d,i,tid);
			z = D_REF(mpzT->d,j,tid);
			s = x * y;
			c2 = mul_hi(x,y);
			s2 = s + z;
			if(s2 < s)
				c2++;
			s = s2 + c;
			if(s < s2)
				c2++;
			c = c2;
	    		D_REF(mpzT->d,j,tid) = s;
		}
		x = D_REF(mpzT->d,size,tid);
		s = x + c;
		c = (s < x)?1:0;
		D_REF(mpzT->d,size,tid) =  s;
		D_REF(mpzT->d,size+1,tid) = c;

		mf = D_REF(mpzT->d,0,tid) * m0;
		x = D_REF(mpzM->d,0,tid);
		y = D_REF(mpzT->d,0,tid);
		s2 = x * mf;
		c = mul_hi(x,mf);
		s = s2 + y;
		if(s < s2)
			c++;

		for(j=1; j < size; j++) {
			x = D_REF(mpzM->d,j,tid);
			y =  D_REF(mpzT->d,j,tid);
			s = x * mf;
			c2 = mul_hi(x,mf);
			s2 = s + y;
			if(s2 < s)
				c2++;
			s = s2 + c;
			if(s < s2)
				c2++;
			c = c2;
			D_REF(mpzT->d,j-1,tid) = s;
		}
		x = D_REF(mpzT->d,size,tid);
		s = x + c;
		c = (s < x)?1:0;
		D_REF(mpzT->d,size-1,tid) =  s;
		D_REF(mpzT->d,size,tid) = D_REF(mpzT->d,size+1,tid) + c;
	}
	mpzT->size[tid] = size;
	copyBN32(tid,mpzR,mpzT);
}

void powmBNu(uint32_t tid, __global fermatTemp_t *restrict temp);
void powmBN(uint32_t tid, __global fermatTemp_t *restrict temp){
	newton_invert(tid,&temp->mpzInv,&temp->mpzM,temp);

    	uint32_t maxbit = norm(tid,&temp->mpzM);
    	uint32_t k = maxbit + 32 - (maxbit%32);
	setBN(tid,&temp->mpzOne,1);
    	//r = 1 << k
	lshiftBN(tid,&temp->mpzR,&temp->mpzOne,k);
	rshiftBN(tid,&temp->mpzHalfR,&temp->mpzR,1);
	//printbn(tid,&temp->mpzHalfR);

   	//(r_, m_) = xbinGCD(r>>1, m)
	xbinGCD(tid,&temp->mpzHalfR,&temp->mpzM,temp);
	
	clearBN32(tid,&temp->mpzResult,k/32);
	clearBN32(tid,&temp->mpzBase,k/32);

	//h = barrett_square(r,m,minv)
	//printbn(tid,&temp->mpzR);
	barrett_square(tid,&temp->mpzH,&temp->mpzR,&temp->mpzM,&temp->mpzInv,temp);
	toBN32(tid,&temp->mpzH32,&temp->mpzH);
	toBN32(tid,&temp->mpzM32,&temp->mpzM);
	toBN32(tid,&temp->mpzV32,&temp->mpzXbinV);

	setBN32(tid,&temp->mpzResult,1);
	setBN32(tid,&temp->mpzBase,2);


	uint32_t mi0 = D_REF(temp->mpzV32.d,0,tid);

	//To montgomery space
	monPro(tid,&temp->mpzResult,&temp->mpzResult,&temp->mpzH32,k,&temp->mpzM32,mi0,temp);
	monPro(tid,&temp->mpzBase,&temp->mpzBase,&temp->mpzH32,k,&temp->mpzM32,mi0,temp);

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
	}

	//From montgomery space
	clearBN32(tid,&temp->mpzOne32,k/32);
	setBN32(tid,&temp->mpzOne32,1);
	monPro(tid,&temp->mpzResult,&temp->mpzResult,&temp->mpzOne32,k,&temp->mpzM32,mi0,temp);
	fromBN32(tid,&temp->mpzR,&temp->mpzResult);
}

uint32_t FermatProbablePrimalityTestFast(uint32_t tid, __global mpzcl_t *restrict origin,  __global fermatTemp_t *restrict temp){
	subBN(tid,&temp->mpzE,origin,1);
	copyBN(tid,&temp->mpzM,origin);
	powmBNu(tid,temp);
	return cmpBNEQ(tid,&temp->mpzR,1);
}

uint32_t ProbablePrimeChainTestFast(uint32_t tid, __global fermatTemp_t *restrict temp, uint32_t type, uint32_t length){
	//Type is going to be the same per warp, so divergence is not a problem. Some blocks will have padded data though
	if(type==-1)
		return 0;

	__global mpzcl_t *origin;
	if(type==BITWIN_TYPE){ //Bitwin
		if(length/2 != 0)
			lshiftBN(tid,&temp->mpzOriginShift,&temp->mpzChainOrigin,length/2);
		else
			copyBN(tid,&temp->mpzOriginShift,&temp->mpzChainOrigin);
		if(length%2){
			addBN(tid,&temp->mpzOriginPlusOne,&temp->mpzOriginShift,1);
			origin = &temp->mpzOriginPlusOne;
		}else{
			subBN(tid,&temp->mpzOriginMinusOne,&temp->mpzOriginShift,1);
			origin = &temp->mpzOriginMinusOne;
		}
	}else{
		if(length != 0)
			lshiftBN(tid,&temp->mpzOriginShift,&temp->mpzChainOrigin,length);
		else
			copyBN(tid,&temp->mpzOriginShift,&temp->mpzChainOrigin);

		if(type==CC1_TYPE){
			subBN(tid,&temp->mpzOriginMinusOne,&temp->mpzOriginShift,1);
			origin = &temp->mpzOriginMinusOne;
		}else if(type==CC2_TYPE){
			addBN(tid,&temp->mpzOriginPlusOne,&temp->mpzOriginShift,1);
			origin = &temp->mpzOriginPlusOne;
		}
	}
	normalizeBN(tid,origin);
	return FermatProbablePrimalityTestFast(tid,origin,temp);
}
