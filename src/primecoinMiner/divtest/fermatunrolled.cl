#define mul_hi(a,b) (uint32_t)(((uint64_t)a*b)>>32)

void xbinGCDu(uint32_t tid, __global mpzcl_t* mpzA, __global mpzcl_t* mpzB, __global fermatTemp_t *temp){
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


void powmBNu(uint32_t tid, __global fermatTemp_t *restrict temp) {
  newton_invert(tid, &temp->mpzInv, &temp->mpzM, temp);

  uint32_t maxbit = norm(tid, &temp->mpzM);
  uint32_t k = maxbit + 32 - (maxbit % 32);
  k=352;

  clearBN(tid,&temp->mpzOne,(k+63)/64);
  clearBN(tid,&temp->mpzHalfR,(k+63)/64);
  clearBN(tid,&temp->mpzR,(k+63)/64);	

//////////////
	D_REF(temp->mpzM.d,temp->mpzM.size[tid],tid)=0;
//////////////

  setBN(tid, &temp->mpzOne, 1);
  // r = 1 << k
  lshiftBN(tid, &temp->mpzR, &temp->mpzOne, k);
  rshiftBN(tid, &temp->mpzHalfR, &temp->mpzR, 1);
  // printbn(tid,&temp->mpzHalfR);

  //(r_, m_) = xbinGCD(r>>1, m)
  xbinGCDu(tid, &temp->mpzHalfR, &temp->mpzM, temp);

  clearBN32(tid,&temp->mpzResult,k/32);
  clearBN32(tid,&temp->mpzBase,k/32);
  clearBN32(tid,&temp->mpzM32,k/32);
  clearBN32(tid,&temp->mpzH32,k/32);
  clearBN32(tid,&temp->mpzV32,k/32);

  // h = barrett_square(r,m,minv)
  // printbn(tid,&temp->mpzR);
  barrett_square(tid, &temp->mpzH, &temp->mpzR, &temp->mpzM, &temp->mpzInv,
                 temp);
  toBN32(tid, &temp->mpzH32, &temp->mpzH);
  toBN32(tid, &temp->mpzM32, &temp->mpzM);
  toBN32(tid, &temp->mpzV32, &temp->mpzXbinV);
  setBN32(tid, &temp->mpzResult, 1);
  setBN32(tid, &temp->mpzBase, 2);

  uint32_t mi0 = D_REF(temp->mpzV32.d, 0, tid);

  // To montgomery space
  monPro(tid, &temp->mpzResult, &temp->mpzResult, &temp->mpzH32, k,
         &temp->mpzM32, mi0, temp);
  monPro(tid, &temp->mpzBase, &temp->mpzBase, &temp->mpzH32, k, &temp->mpzM32,
         mi0, temp);

	printf("Montgomery space:\n");

    printbn32(tid,&temp->mpzH32);
    printbn32(tid,&temp->mpzResult);
    printbn32(tid,&temp->mpzBase);

  // printbn(tid,&temp->mpzE);

 	uint32_t base0;
	uint32_t base1;
	uint32_t base2;
	uint32_t base3;
	uint32_t base4;
	uint32_t base5;
	uint32_t base6;
	uint32_t base7;
	uint32_t base8;
	uint32_t base9;
	uint32_t base10;
	uint32_t result0;
	uint32_t result1;
	uint32_t result2;
	uint32_t result3;
	uint32_t result4;
	uint32_t result5;
	uint32_t result6;
	uint32_t result7;
	uint32_t result8;
	uint32_t result9;
	uint32_t result10;
	uint32_t m0;
	uint32_t m1;
	uint32_t m2;
	uint32_t m3;
	uint32_t m4;
	uint32_t m5;
	uint32_t m6;
	uint32_t m7;
	uint32_t m8;
	uint32_t m9;
	uint32_t m10;
	uint32_t t0;
	uint32_t t1;
	uint32_t t2;
	uint32_t t3;
	uint32_t t4;
	uint32_t t5;
	uint32_t t6;
	uint32_t t7;
	uint32_t t8;
	uint32_t t9;
	uint32_t t10;
	uint32_t t11;
	uint32_t t12;
	uint32_t c,c2,s,s2,mf,x,y,z;
	uint32_t size = k/32;
	{ //Load mpzcl32 to regs BEGIN
		m0 = D_REF(temp->mpzM32.d,0,tid);
		m1 = D_REF(temp->mpzM32.d,1,tid);
		m2 = D_REF(temp->mpzM32.d,2,tid);
		m3 = D_REF(temp->mpzM32.d,3,tid);
		m4 = D_REF(temp->mpzM32.d,4,tid);
		m5 = D_REF(temp->mpzM32.d,5,tid);
		m6 = D_REF(temp->mpzM32.d,6,tid);
		m7 = D_REF(temp->mpzM32.d,7,tid);
		m8 = D_REF(temp->mpzM32.d,8,tid);
		m9 = D_REF(temp->mpzM32.d,9,tid);
		m10 = D_REF(temp->mpzM32.d,10,tid);
	} //Load END
	{ //Load mpzcl32 to regs BEGIN
		base0 = D_REF(temp->mpzBase.d,0,tid);
		base1 = D_REF(temp->mpzBase.d,1,tid);
		base2 = D_REF(temp->mpzBase.d,2,tid);
		base3 = D_REF(temp->mpzBase.d,3,tid);
		base4 = D_REF(temp->mpzBase.d,4,tid);
		base5 = D_REF(temp->mpzBase.d,5,tid);
		base6 = D_REF(temp->mpzBase.d,6,tid);
		base7 = D_REF(temp->mpzBase.d,7,tid);
		base8 = D_REF(temp->mpzBase.d,8,tid);
		base9 = D_REF(temp->mpzBase.d,9,tid);
		base10 = D_REF(temp->mpzBase.d,10,tid);
	} //Load END
	{ //Load mpzcl32 to regs BEGIN
		result0 = D_REF(temp->mpzResult.d,0,tid);
		result1 = D_REF(temp->mpzResult.d,1,tid);
		result2 = D_REF(temp->mpzResult.d,2,tid);
		result3 = D_REF(temp->mpzResult.d,3,tid);
		result4 = D_REF(temp->mpzResult.d,4,tid);
		result5 = D_REF(temp->mpzResult.d,5,tid);
		result6 = D_REF(temp->mpzResult.d,6,tid);
		result7 = D_REF(temp->mpzResult.d,7,tid);
		result8 = D_REF(temp->mpzResult.d,8,tid);
		result9 = D_REF(temp->mpzResult.d,9,tid);
		result10 = D_REF(temp->mpzResult.d,10,tid);
	} //Load END
	uint32_t l;
	for(l=0; l <= maxbit; l++){
		if(bitIsSet(tid,&temp->mpzE,l)){
			{ //monpro BEGIN
				t0 = 0;
				t1 = 0;
				t2 = 0;
				t3 = 0;
				t4 = 0;
				t5 = 0;
				t6 = 0;
				t7 = 0;
				t8 = 0;
				t9 = 0;
				t10 = 0;
				t11 = 0;
				t12 = 0;

				{
					c = 0;
					{
						x = result0;
						y = base0;
						z = t0;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t0 = s;
					}
					{
						x = result1;
						y = base0;
						z = t1;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t1 = s;
					}
					{
						x = result2;
						y = base0;
						z = t2;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t2 = s;
					}
					{
						x = result3;
						y = base0;
						z = t3;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t3 = s;
					}
					{
						x = result4;
						y = base0;
						z = t4;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t4 = s;
					}
					{
						x = result5;
						y = base0;
						z = t5;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t5 = s;
					}
					{
						x = result6;
						y = base0;
						z = t6;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t6 = s;
					}
					{
						x = result7;
						y = base0;
						z = t7;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t7 = s;
					}
					{
						x = result8;
						y = base0;
						z = t8;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t8 = s;
					}
					{
						x = result9;
						y = base0;
						z = t9;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t9 = s;
					}
					{
						x = result10;
						y = base0;
						z = t10;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t10 = s;
					}
					x = t11;
					s = x + c;
					c = (s < x)?1:0;
					t11 =  s;
					t12 = c;
					mf = (t0 * mi0);
					x = m0;
					y = t0;
					s2 = x * mf;
					c = mul_hi(x,mf);
					s = s2 + y;
					if(s < s2)
						c++;
					{
						x = m1;
						y = t1;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t0 = s;
					}
					{
						x = m2;
						y = t2;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t1 = s;
					}
					{
						x = m3;
						y = t3;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t2 = s;
					}
					{
						x = m4;
						y = t4;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t3 = s;
					}
					{
						x = m5;
						y = t5;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t4 = s;
					}
					{
						x = m6;
						y = t6;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t5 = s;
					}
					{
						x = m7;
						y = t7;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t6 = s;
					}
					{
						x = m8;
						y = t8;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t7 = s;
					}
					{
						x = m9;
						y = t9;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t8 = s;
					}
					{
						x = m10;
						y = t10;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t9 = s;
					}
					x = t11;
					s = x + c;
					c = (s < x)?1:0;
					t10 =  s;
					t11 = t12 + c;
				}
				{
					c = 0;
					{
						x = result0;
						y = base1;
						z = t0;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t0 = s;
					}
					{
						x = result1;
						y = base1;
						z = t1;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t1 = s;
					}
					{
						x = result2;
						y = base1;
						z = t2;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t2 = s;
					}
					{
						x = result3;
						y = base1;
						z = t3;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t3 = s;
					}
					{
						x = result4;
						y = base1;
						z = t4;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t4 = s;
					}
					{
						x = result5;
						y = base1;
						z = t5;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t5 = s;
					}
					{
						x = result6;
						y = base1;
						z = t6;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t6 = s;
					}
					{
						x = result7;
						y = base1;
						z = t7;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t7 = s;
					}
					{
						x = result8;
						y = base1;
						z = t8;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t8 = s;
					}
					{
						x = result9;
						y = base1;
						z = t9;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t9 = s;
					}
					{
						x = result10;
						y = base1;
						z = t10;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t10 = s;
					}
					x = t11;
					s = x + c;
					c = (s < x)?1:0;
					t11 =  s;
					t12 = c;
					mf = (t0 * mi0);
					x = m0;
					y = t0;
					s2 = x * mf;
					c = mul_hi(x,mf);
					s = s2 + y;
					if(s < s2)
						c++;
					{
						x = m1;
						y = t1;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t0 = s;
					}
					{
						x = m2;
						y = t2;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t1 = s;
					}
					{
						x = m3;
						y = t3;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t2 = s;
					}
					{
						x = m4;
						y = t4;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t3 = s;
					}
					{
						x = m5;
						y = t5;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t4 = s;
					}
					{
						x = m6;
						y = t6;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t5 = s;
					}
					{
						x = m7;
						y = t7;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t6 = s;
					}
					{
						x = m8;
						y = t8;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t7 = s;
					}
					{
						x = m9;
						y = t9;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t8 = s;
					}
					{
						x = m10;
						y = t10;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t9 = s;
					}
					x = t11;
					s = x + c;
					c = (s < x)?1:0;
					t10 =  s;
					t11 = t12 + c;
				}
				{
					c = 0;
					{
						x = result0;
						y = base2;
						z = t0;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t0 = s;
					}
					{
						x = result1;
						y = base2;
						z = t1;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t1 = s;
					}
					{
						x = result2;
						y = base2;
						z = t2;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t2 = s;
					}
					{
						x = result3;
						y = base2;
						z = t3;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t3 = s;
					}
					{
						x = result4;
						y = base2;
						z = t4;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t4 = s;
					}
					{
						x = result5;
						y = base2;
						z = t5;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t5 = s;
					}
					{
						x = result6;
						y = base2;
						z = t6;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t6 = s;
					}
					{
						x = result7;
						y = base2;
						z = t7;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t7 = s;
					}
					{
						x = result8;
						y = base2;
						z = t8;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t8 = s;
					}
					{
						x = result9;
						y = base2;
						z = t9;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t9 = s;
					}
					{
						x = result10;
						y = base2;
						z = t10;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t10 = s;
					}
					x = t11;
					s = x + c;
					c = (s < x)?1:0;
					t11 =  s;
					t12 = c;
					mf = (t0 * mi0);
					x = m0;
					y = t0;
					s2 = x * mf;
					c = mul_hi(x,mf);
					s = s2 + y;
					if(s < s2)
						c++;
					{
						x = m1;
						y = t1;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t0 = s;
					}
					{
						x = m2;
						y = t2;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t1 = s;
					}
					{
						x = m3;
						y = t3;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t2 = s;
					}
					{
						x = m4;
						y = t4;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t3 = s;
					}
					{
						x = m5;
						y = t5;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t4 = s;
					}
					{
						x = m6;
						y = t6;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t5 = s;
					}
					{
						x = m7;
						y = t7;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t6 = s;
					}
					{
						x = m8;
						y = t8;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t7 = s;
					}
					{
						x = m9;
						y = t9;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t8 = s;
					}
					{
						x = m10;
						y = t10;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t9 = s;
					}
					x = t11;
					s = x + c;
					c = (s < x)?1:0;
					t10 =  s;
					t11 = t12 + c;
				}
				{
					c = 0;
					{
						x = result0;
						y = base3;
						z = t0;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t0 = s;
					}
					{
						x = result1;
						y = base3;
						z = t1;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t1 = s;
					}
					{
						x = result2;
						y = base3;
						z = t2;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t2 = s;
					}
					{
						x = result3;
						y = base3;
						z = t3;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t3 = s;
					}
					{
						x = result4;
						y = base3;
						z = t4;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t4 = s;
					}
					{
						x = result5;
						y = base3;
						z = t5;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t5 = s;
					}
					{
						x = result6;
						y = base3;
						z = t6;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t6 = s;
					}
					{
						x = result7;
						y = base3;
						z = t7;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t7 = s;
					}
					{
						x = result8;
						y = base3;
						z = t8;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t8 = s;
					}
					{
						x = result9;
						y = base3;
						z = t9;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t9 = s;
					}
					{
						x = result10;
						y = base3;
						z = t10;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t10 = s;
					}
					x = t11;
					s = x + c;
					c = (s < x)?1:0;
					t11 =  s;
					t12 = c;
					mf = (t0 * mi0);
					x = m0;
					y = t0;
					s2 = x * mf;
					c = mul_hi(x,mf);
					s = s2 + y;
					if(s < s2)
						c++;
					{
						x = m1;
						y = t1;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t0 = s;
					}
					{
						x = m2;
						y = t2;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t1 = s;
					}
					{
						x = m3;
						y = t3;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t2 = s;
					}
					{
						x = m4;
						y = t4;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t3 = s;
					}
					{
						x = m5;
						y = t5;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t4 = s;
					}
					{
						x = m6;
						y = t6;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t5 = s;
					}
					{
						x = m7;
						y = t7;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t6 = s;
					}
					{
						x = m8;
						y = t8;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t7 = s;
					}
					{
						x = m9;
						y = t9;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t8 = s;
					}
					{
						x = m10;
						y = t10;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t9 = s;
					}
					x = t11;
					s = x + c;
					c = (s < x)?1:0;
					t10 =  s;
					t11 = t12 + c;
				}
				{
					c = 0;
					{
						x = result0;
						y = base4;
						z = t0;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t0 = s;
					}
					{
						x = result1;
						y = base4;
						z = t1;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t1 = s;
					}
					{
						x = result2;
						y = base4;
						z = t2;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t2 = s;
					}
					{
						x = result3;
						y = base4;
						z = t3;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t3 = s;
					}
					{
						x = result4;
						y = base4;
						z = t4;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t4 = s;
					}
					{
						x = result5;
						y = base4;
						z = t5;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t5 = s;
					}
					{
						x = result6;
						y = base4;
						z = t6;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t6 = s;
					}
					{
						x = result7;
						y = base4;
						z = t7;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t7 = s;
					}
					{
						x = result8;
						y = base4;
						z = t8;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t8 = s;
					}
					{
						x = result9;
						y = base4;
						z = t9;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t9 = s;
					}
					{
						x = result10;
						y = base4;
						z = t10;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t10 = s;
					}
					x = t11;
					s = x + c;
					c = (s < x)?1:0;
					t11 =  s;
					t12 = c;
					mf = (t0 * mi0);
					x = m0;
					y = t0;
					s2 = x * mf;
					c = mul_hi(x,mf);
					s = s2 + y;
					if(s < s2)
						c++;
					{
						x = m1;
						y = t1;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t0 = s;
					}
					{
						x = m2;
						y = t2;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t1 = s;
					}
					{
						x = m3;
						y = t3;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t2 = s;
					}
					{
						x = m4;
						y = t4;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t3 = s;
					}
					{
						x = m5;
						y = t5;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t4 = s;
					}
					{
						x = m6;
						y = t6;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t5 = s;
					}
					{
						x = m7;
						y = t7;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t6 = s;
					}
					{
						x = m8;
						y = t8;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t7 = s;
					}
					{
						x = m9;
						y = t9;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t8 = s;
					}
					{
						x = m10;
						y = t10;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t9 = s;
					}
					x = t11;
					s = x + c;
					c = (s < x)?1:0;
					t10 =  s;
					t11 = t12 + c;
				}
				{
					c = 0;
					{
						x = result0;
						y = base5;
						z = t0;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t0 = s;
					}
					{
						x = result1;
						y = base5;
						z = t1;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t1 = s;
					}
					{
						x = result2;
						y = base5;
						z = t2;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t2 = s;
					}
					{
						x = result3;
						y = base5;
						z = t3;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t3 = s;
					}
					{
						x = result4;
						y = base5;
						z = t4;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t4 = s;
					}
					{
						x = result5;
						y = base5;
						z = t5;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t5 = s;
					}
					{
						x = result6;
						y = base5;
						z = t6;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t6 = s;
					}
					{
						x = result7;
						y = base5;
						z = t7;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t7 = s;
					}
					{
						x = result8;
						y = base5;
						z = t8;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t8 = s;
					}
					{
						x = result9;
						y = base5;
						z = t9;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t9 = s;
					}
					{
						x = result10;
						y = base5;
						z = t10;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t10 = s;
					}
					x = t11;
					s = x + c;
					c = (s < x)?1:0;
					t11 =  s;
					t12 = c;
					mf = (t0 * mi0);
					x = m0;
					y = t0;
					s2 = x * mf;
					c = mul_hi(x,mf);
					s = s2 + y;
					if(s < s2)
						c++;
					{
						x = m1;
						y = t1;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t0 = s;
					}
					{
						x = m2;
						y = t2;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t1 = s;
					}
					{
						x = m3;
						y = t3;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t2 = s;
					}
					{
						x = m4;
						y = t4;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t3 = s;
					}
					{
						x = m5;
						y = t5;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t4 = s;
					}
					{
						x = m6;
						y = t6;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t5 = s;
					}
					{
						x = m7;
						y = t7;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t6 = s;
					}
					{
						x = m8;
						y = t8;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t7 = s;
					}
					{
						x = m9;
						y = t9;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t8 = s;
					}
					{
						x = m10;
						y = t10;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t9 = s;
					}
					x = t11;
					s = x + c;
					c = (s < x)?1:0;
					t10 =  s;
					t11 = t12 + c;
				}
				{
					c = 0;
					{
						x = result0;
						y = base6;
						z = t0;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t0 = s;
					}
					{
						x = result1;
						y = base6;
						z = t1;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t1 = s;
					}
					{
						x = result2;
						y = base6;
						z = t2;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t2 = s;
					}
					{
						x = result3;
						y = base6;
						z = t3;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t3 = s;
					}
					{
						x = result4;
						y = base6;
						z = t4;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t4 = s;
					}
					{
						x = result5;
						y = base6;
						z = t5;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t5 = s;
					}
					{
						x = result6;
						y = base6;
						z = t6;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t6 = s;
					}
					{
						x = result7;
						y = base6;
						z = t7;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t7 = s;
					}
					{
						x = result8;
						y = base6;
						z = t8;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t8 = s;
					}
					{
						x = result9;
						y = base6;
						z = t9;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t9 = s;
					}
					{
						x = result10;
						y = base6;
						z = t10;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t10 = s;
					}
					x = t11;
					s = x + c;
					c = (s < x)?1:0;
					t11 =  s;
					t12 = c;
					mf = (t0 * mi0);
					x = m0;
					y = t0;
					s2 = x * mf;
					c = mul_hi(x,mf);
					s = s2 + y;
					if(s < s2)
						c++;
					{
						x = m1;
						y = t1;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t0 = s;
					}
					{
						x = m2;
						y = t2;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t1 = s;
					}
					{
						x = m3;
						y = t3;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t2 = s;
					}
					{
						x = m4;
						y = t4;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t3 = s;
					}
					{
						x = m5;
						y = t5;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t4 = s;
					}
					{
						x = m6;
						y = t6;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t5 = s;
					}
					{
						x = m7;
						y = t7;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t6 = s;
					}
					{
						x = m8;
						y = t8;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t7 = s;
					}
					{
						x = m9;
						y = t9;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t8 = s;
					}
					{
						x = m10;
						y = t10;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t9 = s;
					}
					x = t11;
					s = x + c;
					c = (s < x)?1:0;
					t10 =  s;
					t11 = t12 + c;
				}
				{
					c = 0;
					{
						x = result0;
						y = base7;
						z = t0;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t0 = s;
					}
					{
						x = result1;
						y = base7;
						z = t1;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t1 = s;
					}
					{
						x = result2;
						y = base7;
						z = t2;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t2 = s;
					}
					{
						x = result3;
						y = base7;
						z = t3;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t3 = s;
					}
					{
						x = result4;
						y = base7;
						z = t4;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t4 = s;
					}
					{
						x = result5;
						y = base7;
						z = t5;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t5 = s;
					}
					{
						x = result6;
						y = base7;
						z = t6;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t6 = s;
					}
					{
						x = result7;
						y = base7;
						z = t7;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t7 = s;
					}
					{
						x = result8;
						y = base7;
						z = t8;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t8 = s;
					}
					{
						x = result9;
						y = base7;
						z = t9;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t9 = s;
					}
					{
						x = result10;
						y = base7;
						z = t10;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t10 = s;
					}
					x = t11;
					s = x + c;
					c = (s < x)?1:0;
					t11 =  s;
					t12 = c;
					mf = (t0 * mi0);
					x = m0;
					y = t0;
					s2 = x * mf;
					c = mul_hi(x,mf);
					s = s2 + y;
					if(s < s2)
						c++;
					{
						x = m1;
						y = t1;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t0 = s;
					}
					{
						x = m2;
						y = t2;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t1 = s;
					}
					{
						x = m3;
						y = t3;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t2 = s;
					}
					{
						x = m4;
						y = t4;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t3 = s;
					}
					{
						x = m5;
						y = t5;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t4 = s;
					}
					{
						x = m6;
						y = t6;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t5 = s;
					}
					{
						x = m7;
						y = t7;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t6 = s;
					}
					{
						x = m8;
						y = t8;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t7 = s;
					}
					{
						x = m9;
						y = t9;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t8 = s;
					}
					{
						x = m10;
						y = t10;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t9 = s;
					}
					x = t11;
					s = x + c;
					c = (s < x)?1:0;
					t10 =  s;
					t11 = t12 + c;
				}
				{
					c = 0;
					{
						x = result0;
						y = base8;
						z = t0;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t0 = s;
					}
					{
						x = result1;
						y = base8;
						z = t1;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t1 = s;
					}
					{
						x = result2;
						y = base8;
						z = t2;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t2 = s;
					}
					{
						x = result3;
						y = base8;
						z = t3;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t3 = s;
					}
					{
						x = result4;
						y = base8;
						z = t4;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t4 = s;
					}
					{
						x = result5;
						y = base8;
						z = t5;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t5 = s;
					}
					{
						x = result6;
						y = base8;
						z = t6;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t6 = s;
					}
					{
						x = result7;
						y = base8;
						z = t7;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t7 = s;
					}
					{
						x = result8;
						y = base8;
						z = t8;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t8 = s;
					}
					{
						x = result9;
						y = base8;
						z = t9;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t9 = s;
					}
					{
						x = result10;
						y = base8;
						z = t10;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t10 = s;
					}
					x = t11;
					s = x + c;
					c = (s < x)?1:0;
					t11 =  s;
					t12 = c;
					mf = (t0 * mi0);
					x = m0;
					y = t0;
					s2 = x * mf;
					c = mul_hi(x,mf);
					s = s2 + y;
					if(s < s2)
						c++;
					{
						x = m1;
						y = t1;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t0 = s;
					}
					{
						x = m2;
						y = t2;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t1 = s;
					}
					{
						x = m3;
						y = t3;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t2 = s;
					}
					{
						x = m4;
						y = t4;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t3 = s;
					}
					{
						x = m5;
						y = t5;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t4 = s;
					}
					{
						x = m6;
						y = t6;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t5 = s;
					}
					{
						x = m7;
						y = t7;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t6 = s;
					}
					{
						x = m8;
						y = t8;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t7 = s;
					}
					{
						x = m9;
						y = t9;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t8 = s;
					}
					{
						x = m10;
						y = t10;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t9 = s;
					}
					x = t11;
					s = x + c;
					c = (s < x)?1:0;
					t10 =  s;
					t11 = t12 + c;
				}
				{
					c = 0;
					{
						x = result0;
						y = base9;
						z = t0;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t0 = s;
					}
					{
						x = result1;
						y = base9;
						z = t1;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t1 = s;
					}
					{
						x = result2;
						y = base9;
						z = t2;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t2 = s;
					}
					{
						x = result3;
						y = base9;
						z = t3;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t3 = s;
					}
					{
						x = result4;
						y = base9;
						z = t4;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t4 = s;
					}
					{
						x = result5;
						y = base9;
						z = t5;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t5 = s;
					}
					{
						x = result6;
						y = base9;
						z = t6;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t6 = s;
					}
					{
						x = result7;
						y = base9;
						z = t7;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t7 = s;
					}
					{
						x = result8;
						y = base9;
						z = t8;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t8 = s;
					}
					{
						x = result9;
						y = base9;
						z = t9;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t9 = s;
					}
					{
						x = result10;
						y = base9;
						z = t10;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t10 = s;
					}
					x = t11;
					s = x + c;
					c = (s < x)?1:0;
					t11 =  s;
					t12 = c;
					mf = (t0 * mi0);
					x = m0;
					y = t0;
					s2 = x * mf;
					c = mul_hi(x,mf);
					s = s2 + y;
					if(s < s2)
						c++;
					{
						x = m1;
						y = t1;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t0 = s;
					}
					{
						x = m2;
						y = t2;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t1 = s;
					}
					{
						x = m3;
						y = t3;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t2 = s;
					}
					{
						x = m4;
						y = t4;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t3 = s;
					}
					{
						x = m5;
						y = t5;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t4 = s;
					}
					{
						x = m6;
						y = t6;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t5 = s;
					}
					{
						x = m7;
						y = t7;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t6 = s;
					}
					{
						x = m8;
						y = t8;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t7 = s;
					}
					{
						x = m9;
						y = t9;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t8 = s;
					}
					{
						x = m10;
						y = t10;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t9 = s;
					}
					x = t11;
					s = x + c;
					c = (s < x)?1:0;
					t10 =  s;
					t11 = t12 + c;
				}
				{
					c = 0;
					{
						x = result0;
						y = base10;
						z = t0;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t0 = s;
					}
					{
						x = result1;
						y = base10;
						z = t1;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t1 = s;
					}
					{
						x = result2;
						y = base10;
						z = t2;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t2 = s;
					}
					{
						x = result3;
						y = base10;
						z = t3;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t3 = s;
					}
					{
						x = result4;
						y = base10;
						z = t4;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t4 = s;
					}
					{
						x = result5;
						y = base10;
						z = t5;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t5 = s;
					}
					{
						x = result6;
						y = base10;
						z = t6;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t6 = s;
					}
					{
						x = result7;
						y = base10;
						z = t7;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t7 = s;
					}
					{
						x = result8;
						y = base10;
						z = t8;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t8 = s;
					}
					{
						x = result9;
						y = base10;
						z = t9;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t9 = s;
					}
					{
						x = result10;
						y = base10;
						z = t10;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t10 = s;
					}
					x = t11;
					s = x + c;
					c = (s < x)?1:0;
					t11 =  s;
					t12 = c;
					mf = (t0 * mi0);
					x = m0;
					y = t0;
					s2 = x * mf;
					c = mul_hi(x,mf);
					s = s2 + y;
					if(s < s2)
						c++;
					{
						x = m1;
						y = t1;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t0 = s;
					}
					{
						x = m2;
						y = t2;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t1 = s;
					}
					{
						x = m3;
						y = t3;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t2 = s;
					}
					{
						x = m4;
						y = t4;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t3 = s;
					}
					{
						x = m5;
						y = t5;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t4 = s;
					}
					{
						x = m6;
						y = t6;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t5 = s;
					}
					{
						x = m7;
						y = t7;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t6 = s;
					}
					{
						x = m8;
						y = t8;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t7 = s;
					}
					{
						x = m9;
						y = t9;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t8 = s;
					}
					{
						x = m10;
						y = t10;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t9 = s;
					}
					x = t11;
					s = x + c;
					c = (s < x)?1:0;
					t10 =  s;
					t11 = t12 + c;
				}
				{
					result0 = t0;
					result1 = t1;
					result2 = t2;
					result3 = t3;
					result4 = t4;
					result5 = t5;
					result6 = t6;
					result7 = t7;
					result8 = t8;
					result9 = t9;
					result10 = t10;
				}
			} //END
		}
			{ //monpro BEGIN
				t0 = 0;
				t1 = 0;
				t2 = 0;
				t3 = 0;
				t4 = 0;
				t5 = 0;
				t6 = 0;
				t7 = 0;
				t8 = 0;
				t9 = 0;
				t10 = 0;
				t11 = 0;
				t12 = 0;

				{
					c = 0;
					{
						x = base0;
						y = base0;
						z = t0;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t0 = s;
					}
					{
						x = base1;
						y = base0;
						z = t1;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t1 = s;
					}
					{
						x = base2;
						y = base0;
						z = t2;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t2 = s;
					}
					{
						x = base3;
						y = base0;
						z = t3;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t3 = s;
					}
					{
						x = base4;
						y = base0;
						z = t4;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t4 = s;
					}
					{
						x = base5;
						y = base0;
						z = t5;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t5 = s;
					}
					{
						x = base6;
						y = base0;
						z = t6;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t6 = s;
					}
					{
						x = base7;
						y = base0;
						z = t7;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t7 = s;
					}
					{
						x = base8;
						y = base0;
						z = t8;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t8 = s;
					}
					{
						x = base9;
						y = base0;
						z = t9;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t9 = s;
					}
					{
						x = base10;
						y = base0;
						z = t10;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t10 = s;
					}
					x = t11;
					s = x + c;
					c = (s < x)?1:0;
					t11 =  s;
					t12 = c;
					mf = (t0 * mi0);
					x = m0;
					y = t0;
					s2 = x * mf;
					c = mul_hi(x,mf);
					s = s2 + y;
					if(s < s2)
						c++;
					{
						x = m1;
						y = t1;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t0 = s;
					}
					{
						x = m2;
						y = t2;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t1 = s;
					}
					{
						x = m3;
						y = t3;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t2 = s;
					}
					{
						x = m4;
						y = t4;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t3 = s;
					}
					{
						x = m5;
						y = t5;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t4 = s;
					}
					{
						x = m6;
						y = t6;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t5 = s;
					}
					{
						x = m7;
						y = t7;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t6 = s;
					}
					{
						x = m8;
						y = t8;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t7 = s;
					}
					{
						x = m9;
						y = t9;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t8 = s;
					}
					{
						x = m10;
						y = t10;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t9 = s;
					}
					x = t11;
					s = x + c;
					c = (s < x)?1:0;
					t10 =  s;
					t11 = t12 + c;
				}
				{
					c = 0;
					{
						x = base0;
						y = base1;
						z = t0;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t0 = s;
					}
					{
						x = base1;
						y = base1;
						z = t1;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t1 = s;
					}
					{
						x = base2;
						y = base1;
						z = t2;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t2 = s;
					}
					{
						x = base3;
						y = base1;
						z = t3;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t3 = s;
					}
					{
						x = base4;
						y = base1;
						z = t4;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t4 = s;
					}
					{
						x = base5;
						y = base1;
						z = t5;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t5 = s;
					}
					{
						x = base6;
						y = base1;
						z = t6;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t6 = s;
					}
					{
						x = base7;
						y = base1;
						z = t7;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t7 = s;
					}
					{
						x = base8;
						y = base1;
						z = t8;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t8 = s;
					}
					{
						x = base9;
						y = base1;
						z = t9;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t9 = s;
					}
					{
						x = base10;
						y = base1;
						z = t10;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t10 = s;
					}
					x = t11;
					s = x + c;
					c = (s < x)?1:0;
					t11 =  s;
					t12 = c;
					mf = (t0 * mi0);
					x = m0;
					y = t0;
					s2 = x * mf;
					c = mul_hi(x,mf);
					s = s2 + y;
					if(s < s2)
						c++;
					{
						x = m1;
						y = t1;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t0 = s;
					}
					{
						x = m2;
						y = t2;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t1 = s;
					}
					{
						x = m3;
						y = t3;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t2 = s;
					}
					{
						x = m4;
						y = t4;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t3 = s;
					}
					{
						x = m5;
						y = t5;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t4 = s;
					}
					{
						x = m6;
						y = t6;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t5 = s;
					}
					{
						x = m7;
						y = t7;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t6 = s;
					}
					{
						x = m8;
						y = t8;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t7 = s;
					}
					{
						x = m9;
						y = t9;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t8 = s;
					}
					{
						x = m10;
						y = t10;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t9 = s;
					}
					x = t11;
					s = x + c;
					c = (s < x)?1:0;
					t10 =  s;
					t11 = t12 + c;
				}
				{
					c = 0;
					{
						x = base0;
						y = base2;
						z = t0;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t0 = s;
					}
					{
						x = base1;
						y = base2;
						z = t1;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t1 = s;
					}
					{
						x = base2;
						y = base2;
						z = t2;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t2 = s;
					}
					{
						x = base3;
						y = base2;
						z = t3;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t3 = s;
					}
					{
						x = base4;
						y = base2;
						z = t4;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t4 = s;
					}
					{
						x = base5;
						y = base2;
						z = t5;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t5 = s;
					}
					{
						x = base6;
						y = base2;
						z = t6;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t6 = s;
					}
					{
						x = base7;
						y = base2;
						z = t7;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t7 = s;
					}
					{
						x = base8;
						y = base2;
						z = t8;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t8 = s;
					}
					{
						x = base9;
						y = base2;
						z = t9;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t9 = s;
					}
					{
						x = base10;
						y = base2;
						z = t10;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t10 = s;
					}
					x = t11;
					s = x + c;
					c = (s < x)?1:0;
					t11 =  s;
					t12 = c;
					mf = (t0 * mi0);
					x = m0;
					y = t0;
					s2 = x * mf;
					c = mul_hi(x,mf);
					s = s2 + y;
					if(s < s2)
						c++;
					{
						x = m1;
						y = t1;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t0 = s;
					}
					{
						x = m2;
						y = t2;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t1 = s;
					}
					{
						x = m3;
						y = t3;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t2 = s;
					}
					{
						x = m4;
						y = t4;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t3 = s;
					}
					{
						x = m5;
						y = t5;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t4 = s;
					}
					{
						x = m6;
						y = t6;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t5 = s;
					}
					{
						x = m7;
						y = t7;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t6 = s;
					}
					{
						x = m8;
						y = t8;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t7 = s;
					}
					{
						x = m9;
						y = t9;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t8 = s;
					}
					{
						x = m10;
						y = t10;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t9 = s;
					}
					x = t11;
					s = x + c;
					c = (s < x)?1:0;
					t10 =  s;
					t11 = t12 + c;
				}
				{
					c = 0;
					{
						x = base0;
						y = base3;
						z = t0;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t0 = s;
					}
					{
						x = base1;
						y = base3;
						z = t1;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t1 = s;
					}
					{
						x = base2;
						y = base3;
						z = t2;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t2 = s;
					}
					{
						x = base3;
						y = base3;
						z = t3;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t3 = s;
					}
					{
						x = base4;
						y = base3;
						z = t4;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t4 = s;
					}
					{
						x = base5;
						y = base3;
						z = t5;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t5 = s;
					}
					{
						x = base6;
						y = base3;
						z = t6;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t6 = s;
					}
					{
						x = base7;
						y = base3;
						z = t7;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t7 = s;
					}
					{
						x = base8;
						y = base3;
						z = t8;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t8 = s;
					}
					{
						x = base9;
						y = base3;
						z = t9;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t9 = s;
					}
					{
						x = base10;
						y = base3;
						z = t10;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t10 = s;
					}
					x = t11;
					s = x + c;
					c = (s < x)?1:0;
					t11 =  s;
					t12 = c;
					mf = (t0 * mi0);
					x = m0;
					y = t0;
					s2 = x * mf;
					c = mul_hi(x,mf);
					s = s2 + y;
					if(s < s2)
						c++;
					{
						x = m1;
						y = t1;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t0 = s;
					}
					{
						x = m2;
						y = t2;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t1 = s;
					}
					{
						x = m3;
						y = t3;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t2 = s;
					}
					{
						x = m4;
						y = t4;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t3 = s;
					}
					{
						x = m5;
						y = t5;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t4 = s;
					}
					{
						x = m6;
						y = t6;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t5 = s;
					}
					{
						x = m7;
						y = t7;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t6 = s;
					}
					{
						x = m8;
						y = t8;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t7 = s;
					}
					{
						x = m9;
						y = t9;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t8 = s;
					}
					{
						x = m10;
						y = t10;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t9 = s;
					}
					x = t11;
					s = x + c;
					c = (s < x)?1:0;
					t10 =  s;
					t11 = t12 + c;
				}
				{
					c = 0;
					{
						x = base0;
						y = base4;
						z = t0;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t0 = s;
					}
					{
						x = base1;
						y = base4;
						z = t1;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t1 = s;
					}
					{
						x = base2;
						y = base4;
						z = t2;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t2 = s;
					}
					{
						x = base3;
						y = base4;
						z = t3;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t3 = s;
					}
					{
						x = base4;
						y = base4;
						z = t4;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t4 = s;
					}
					{
						x = base5;
						y = base4;
						z = t5;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t5 = s;
					}
					{
						x = base6;
						y = base4;
						z = t6;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t6 = s;
					}
					{
						x = base7;
						y = base4;
						z = t7;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t7 = s;
					}
					{
						x = base8;
						y = base4;
						z = t8;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t8 = s;
					}
					{
						x = base9;
						y = base4;
						z = t9;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t9 = s;
					}
					{
						x = base10;
						y = base4;
						z = t10;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t10 = s;
					}
					x = t11;
					s = x + c;
					c = (s < x)?1:0;
					t11 =  s;
					t12 = c;
					mf = (t0 * mi0);
					x = m0;
					y = t0;
					s2 = x * mf;
					c = mul_hi(x,mf);
					s = s2 + y;
					if(s < s2)
						c++;
					{
						x = m1;
						y = t1;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t0 = s;
					}
					{
						x = m2;
						y = t2;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t1 = s;
					}
					{
						x = m3;
						y = t3;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t2 = s;
					}
					{
						x = m4;
						y = t4;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t3 = s;
					}
					{
						x = m5;
						y = t5;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t4 = s;
					}
					{
						x = m6;
						y = t6;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t5 = s;
					}
					{
						x = m7;
						y = t7;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t6 = s;
					}
					{
						x = m8;
						y = t8;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t7 = s;
					}
					{
						x = m9;
						y = t9;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t8 = s;
					}
					{
						x = m10;
						y = t10;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t9 = s;
					}
					x = t11;
					s = x + c;
					c = (s < x)?1:0;
					t10 =  s;
					t11 = t12 + c;
				}
				{
					c = 0;
					{
						x = base0;
						y = base5;
						z = t0;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t0 = s;
					}
					{
						x = base1;
						y = base5;
						z = t1;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t1 = s;
					}
					{
						x = base2;
						y = base5;
						z = t2;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t2 = s;
					}
					{
						x = base3;
						y = base5;
						z = t3;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t3 = s;
					}
					{
						x = base4;
						y = base5;
						z = t4;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t4 = s;
					}
					{
						x = base5;
						y = base5;
						z = t5;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t5 = s;
					}
					{
						x = base6;
						y = base5;
						z = t6;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t6 = s;
					}
					{
						x = base7;
						y = base5;
						z = t7;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t7 = s;
					}
					{
						x = base8;
						y = base5;
						z = t8;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t8 = s;
					}
					{
						x = base9;
						y = base5;
						z = t9;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t9 = s;
					}
					{
						x = base10;
						y = base5;
						z = t10;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t10 = s;
					}
					x = t11;
					s = x + c;
					c = (s < x)?1:0;
					t11 =  s;
					t12 = c;
					mf = (t0 * mi0);
					x = m0;
					y = t0;
					s2 = x * mf;
					c = mul_hi(x,mf);
					s = s2 + y;
					if(s < s2)
						c++;
					{
						x = m1;
						y = t1;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t0 = s;
					}
					{
						x = m2;
						y = t2;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t1 = s;
					}
					{
						x = m3;
						y = t3;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t2 = s;
					}
					{
						x = m4;
						y = t4;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t3 = s;
					}
					{
						x = m5;
						y = t5;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t4 = s;
					}
					{
						x = m6;
						y = t6;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t5 = s;
					}
					{
						x = m7;
						y = t7;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t6 = s;
					}
					{
						x = m8;
						y = t8;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t7 = s;
					}
					{
						x = m9;
						y = t9;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t8 = s;
					}
					{
						x = m10;
						y = t10;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t9 = s;
					}
					x = t11;
					s = x + c;
					c = (s < x)?1:0;
					t10 =  s;
					t11 = t12 + c;
				}
				{
					c = 0;
					{
						x = base0;
						y = base6;
						z = t0;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t0 = s;
					}
					{
						x = base1;
						y = base6;
						z = t1;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t1 = s;
					}
					{
						x = base2;
						y = base6;
						z = t2;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t2 = s;
					}
					{
						x = base3;
						y = base6;
						z = t3;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t3 = s;
					}
					{
						x = base4;
						y = base6;
						z = t4;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t4 = s;
					}
					{
						x = base5;
						y = base6;
						z = t5;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t5 = s;
					}
					{
						x = base6;
						y = base6;
						z = t6;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t6 = s;
					}
					{
						x = base7;
						y = base6;
						z = t7;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t7 = s;
					}
					{
						x = base8;
						y = base6;
						z = t8;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t8 = s;
					}
					{
						x = base9;
						y = base6;
						z = t9;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t9 = s;
					}
					{
						x = base10;
						y = base6;
						z = t10;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t10 = s;
					}
					x = t11;
					s = x + c;
					c = (s < x)?1:0;
					t11 =  s;
					t12 = c;
					mf = (t0 * mi0);
					x = m0;
					y = t0;
					s2 = x * mf;
					c = mul_hi(x,mf);
					s = s2 + y;
					if(s < s2)
						c++;
					{
						x = m1;
						y = t1;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t0 = s;
					}
					{
						x = m2;
						y = t2;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t1 = s;
					}
					{
						x = m3;
						y = t3;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t2 = s;
					}
					{
						x = m4;
						y = t4;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t3 = s;
					}
					{
						x = m5;
						y = t5;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t4 = s;
					}
					{
						x = m6;
						y = t6;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t5 = s;
					}
					{
						x = m7;
						y = t7;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t6 = s;
					}
					{
						x = m8;
						y = t8;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t7 = s;
					}
					{
						x = m9;
						y = t9;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t8 = s;
					}
					{
						x = m10;
						y = t10;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t9 = s;
					}
					x = t11;
					s = x + c;
					c = (s < x)?1:0;
					t10 =  s;
					t11 = t12 + c;
				}
				{
					c = 0;
					{
						x = base0;
						y = base7;
						z = t0;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t0 = s;
					}
					{
						x = base1;
						y = base7;
						z = t1;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t1 = s;
					}
					{
						x = base2;
						y = base7;
						z = t2;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t2 = s;
					}
					{
						x = base3;
						y = base7;
						z = t3;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t3 = s;
					}
					{
						x = base4;
						y = base7;
						z = t4;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t4 = s;
					}
					{
						x = base5;
						y = base7;
						z = t5;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t5 = s;
					}
					{
						x = base6;
						y = base7;
						z = t6;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t6 = s;
					}
					{
						x = base7;
						y = base7;
						z = t7;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t7 = s;
					}
					{
						x = base8;
						y = base7;
						z = t8;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t8 = s;
					}
					{
						x = base9;
						y = base7;
						z = t9;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t9 = s;
					}
					{
						x = base10;
						y = base7;
						z = t10;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t10 = s;
					}
					x = t11;
					s = x + c;
					c = (s < x)?1:0;
					t11 =  s;
					t12 = c;
					mf = (t0 * mi0);
					x = m0;
					y = t0;
					s2 = x * mf;
					c = mul_hi(x,mf);
					s = s2 + y;
					if(s < s2)
						c++;
					{
						x = m1;
						y = t1;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t0 = s;
					}
					{
						x = m2;
						y = t2;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t1 = s;
					}
					{
						x = m3;
						y = t3;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t2 = s;
					}
					{
						x = m4;
						y = t4;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t3 = s;
					}
					{
						x = m5;
						y = t5;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t4 = s;
					}
					{
						x = m6;
						y = t6;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t5 = s;
					}
					{
						x = m7;
						y = t7;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t6 = s;
					}
					{
						x = m8;
						y = t8;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t7 = s;
					}
					{
						x = m9;
						y = t9;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t8 = s;
					}
					{
						x = m10;
						y = t10;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t9 = s;
					}
					x = t11;
					s = x + c;
					c = (s < x)?1:0;
					t10 =  s;
					t11 = t12 + c;
				}
				{
					c = 0;
					{
						x = base0;
						y = base8;
						z = t0;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t0 = s;
					}
					{
						x = base1;
						y = base8;
						z = t1;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t1 = s;
					}
					{
						x = base2;
						y = base8;
						z = t2;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t2 = s;
					}
					{
						x = base3;
						y = base8;
						z = t3;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t3 = s;
					}
					{
						x = base4;
						y = base8;
						z = t4;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t4 = s;
					}
					{
						x = base5;
						y = base8;
						z = t5;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t5 = s;
					}
					{
						x = base6;
						y = base8;
						z = t6;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t6 = s;
					}
					{
						x = base7;
						y = base8;
						z = t7;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t7 = s;
					}
					{
						x = base8;
						y = base8;
						z = t8;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t8 = s;
					}
					{
						x = base9;
						y = base8;
						z = t9;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t9 = s;
					}
					{
						x = base10;
						y = base8;
						z = t10;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t10 = s;
					}
					x = t11;
					s = x + c;
					c = (s < x)?1:0;
					t11 =  s;
					t12 = c;
					mf = (t0 * mi0);
					x = m0;
					y = t0;
					s2 = x * mf;
					c = mul_hi(x,mf);
					s = s2 + y;
					if(s < s2)
						c++;
					{
						x = m1;
						y = t1;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t0 = s;
					}
					{
						x = m2;
						y = t2;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t1 = s;
					}
					{
						x = m3;
						y = t3;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t2 = s;
					}
					{
						x = m4;
						y = t4;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t3 = s;
					}
					{
						x = m5;
						y = t5;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t4 = s;
					}
					{
						x = m6;
						y = t6;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t5 = s;
					}
					{
						x = m7;
						y = t7;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t6 = s;
					}
					{
						x = m8;
						y = t8;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t7 = s;
					}
					{
						x = m9;
						y = t9;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t8 = s;
					}
					{
						x = m10;
						y = t10;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t9 = s;
					}
					x = t11;
					s = x + c;
					c = (s < x)?1:0;
					t10 =  s;
					t11 = t12 + c;
				}
				{
					c = 0;
					{
						x = base0;
						y = base9;
						z = t0;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t0 = s;
					}
					{
						x = base1;
						y = base9;
						z = t1;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t1 = s;
					}
					{
						x = base2;
						y = base9;
						z = t2;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t2 = s;
					}
					{
						x = base3;
						y = base9;
						z = t3;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t3 = s;
					}
					{
						x = base4;
						y = base9;
						z = t4;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t4 = s;
					}
					{
						x = base5;
						y = base9;
						z = t5;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t5 = s;
					}
					{
						x = base6;
						y = base9;
						z = t6;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t6 = s;
					}
					{
						x = base7;
						y = base9;
						z = t7;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t7 = s;
					}
					{
						x = base8;
						y = base9;
						z = t8;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t8 = s;
					}
					{
						x = base9;
						y = base9;
						z = t9;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t9 = s;
					}
					{
						x = base10;
						y = base9;
						z = t10;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t10 = s;
					}
					x = t11;
					s = x + c;
					c = (s < x)?1:0;
					t11 =  s;
					t12 = c;
					mf = (t0 * mi0);
					x = m0;
					y = t0;
					s2 = x * mf;
					c = mul_hi(x,mf);
					s = s2 + y;
					if(s < s2)
						c++;
					{
						x = m1;
						y = t1;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t0 = s;
					}
					{
						x = m2;
						y = t2;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t1 = s;
					}
					{
						x = m3;
						y = t3;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t2 = s;
					}
					{
						x = m4;
						y = t4;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t3 = s;
					}
					{
						x = m5;
						y = t5;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t4 = s;
					}
					{
						x = m6;
						y = t6;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t5 = s;
					}
					{
						x = m7;
						y = t7;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t6 = s;
					}
					{
						x = m8;
						y = t8;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t7 = s;
					}
					{
						x = m9;
						y = t9;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t8 = s;
					}
					{
						x = m10;
						y = t10;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t9 = s;
					}
					x = t11;
					s = x + c;
					c = (s < x)?1:0;
					t10 =  s;
					t11 = t12 + c;
				}
				{
					c = 0;
					{
						x = base0;
						y = base10;
						z = t0;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t0 = s;
					}
					{
						x = base1;
						y = base10;
						z = t1;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t1 = s;
					}
					{
						x = base2;
						y = base10;
						z = t2;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t2 = s;
					}
					{
						x = base3;
						y = base10;
						z = t3;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t3 = s;
					}
					{
						x = base4;
						y = base10;
						z = t4;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t4 = s;
					}
					{
						x = base5;
						y = base10;
						z = t5;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t5 = s;
					}
					{
						x = base6;
						y = base10;
						z = t6;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t6 = s;
					}
					{
						x = base7;
						y = base10;
						z = t7;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t7 = s;
					}
					{
						x = base8;
						y = base10;
						z = t8;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t8 = s;
					}
					{
						x = base9;
						y = base10;
						z = t9;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t9 = s;
					}
					{
						x = base10;
						y = base10;
						z = t10;
						s = x * y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t10 = s;
					}
					x = t11;
					s = x + c;
					c = (s < x)?1:0;
					t11 =  s;
					t12 = c;
					mf = (t0 * mi0);
					x = m0;
					y = t0;
					s2 = x * mf;
					c = mul_hi(x,mf);
					s = s2 + y;
					if(s < s2)
						c++;
					{
						x = m1;
						y = t1;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t0 = s;
					}
					{
						x = m2;
						y = t2;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t1 = s;
					}
					{
						x = m3;
						y = t3;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t2 = s;
					}
					{
						x = m4;
						y = t4;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t3 = s;
					}
					{
						x = m5;
						y = t5;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t4 = s;
					}
					{
						x = m6;
						y = t6;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t5 = s;
					}
					{
						x = m7;
						y = t7;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t6 = s;
					}
					{
						x = m8;
						y = t8;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t7 = s;
					}
					{
						x = m9;
						y = t9;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t8 = s;
					}
					{
						x = m10;
						y = t10;
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
						t9 = s;
					}
					x = t11;
					s = x + c;
					c = (s < x)?1:0;
					t10 =  s;
					t11 = t12 + c;
				}
				{
					base0 = t0;
					base1 = t1;
					base2 = t2;
					base3 = t3;
					base4 = t4;
					base5 = t5;
					base6 = t6;
					base7 = t7;
					base8 = t8;
					base9 = t9;
					base10 = t10;
				}
			} //END
	}
	{ //Save mpzcl32 from regs BEGIN
		D_REF(temp->mpzResult.d,0,tid) = result0;
		D_REF(temp->mpzResult.d,1,tid) = result1;
		D_REF(temp->mpzResult.d,2,tid) = result2;
		D_REF(temp->mpzResult.d,3,tid) = result3;
		D_REF(temp->mpzResult.d,4,tid) = result4;
		D_REF(temp->mpzResult.d,5,tid) = result5;
		D_REF(temp->mpzResult.d,6,tid) = result6;
		D_REF(temp->mpzResult.d,7,tid) = result7;
		D_REF(temp->mpzResult.d,8,tid) = result8;
		D_REF(temp->mpzResult.d,9,tid) = result9;
		D_REF(temp->mpzResult.d,10,tid) = result10;
	} //Save END

  // From montgomery space
  clearBN32(tid,&temp->mpzOne32,k/32);

  setBN32(tid, &temp->mpzOne32, 1);
  monPro(tid, &temp->mpzResult, &temp->mpzResult, &temp->mpzOne32, k,
         &temp->mpzM32, mi0, temp);
  fromBN32(tid, &temp->mpzR, &temp->mpzResult);
}

