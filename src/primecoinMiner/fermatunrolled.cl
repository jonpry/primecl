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

//#define mul_hi(a,b) (uint32_t)(((uint64_t)a*b)>>32)

void powmBNu(uint32_t tid, __global fermatTemp_t *restrict temp) {
  newton_invert(tid, &temp->mpzInv, &temp->mpzM, temp); //9ms

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
  normalizeBN(tid,&temp->mpzR);

  rshiftBN(tid, &temp->mpzHalfR, &temp->mpzR, 1);
  // printbn(tid,&temp->mpzHalfR);

  //(r_, m_) = xbinGCD(r>>1, m)
  xbinGCD(tid, &temp->mpzHalfR, &temp->mpzM, temp); //110ms

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
  uint32_t j=0,i=0;	

  // To montgomery space
  monPro(tid, &temp->mpzResult, &temp->mpzResult, &temp->mpzH32, k,
         &temp->mpzM32, mi0, temp);
  monPro(tid, &temp->mpzBase, &temp->mpzBase, &temp->mpzH32, k, &temp->mpzM32,
         mi0, temp);

	//Totoal = 26ms with below disabled
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
#pragma unroll 0
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

					for(i=0; i < 11; i++)
				{
					c = 0;
						switch(i){
						case 0 :
							y = base0; break;
						case 1 :
							y = base1; break;
						case 2 :
							y = base2; break;
						case 3 :
							y = base3; break;
						case 4 :
							y = base4; break;
						case 5 :
							y = base5; break;
						case 6 :
							y = base6; break;
						case 7 :
							y = base7; break;
						case 8 :
							y = base8; break;
						case 9 :
							y = base9; break;
						case 10 :
							y = base10; break;
						}
					for(j=0; j < 11; j++){
                                        if(j<4){
                                                if(j<2){
							x = result1;
							z = t1;
                                                        if(j<1){
								x = result0;
								z = t0;
							}
                                                }else{
							x = result3;
							z = t3;
                                                        if(j<3){
								x = result2;
								z = t2;
							}
                                                }
                                        }else{
                                                if(j<8){
                                                        if(j<6){
								x = result5;
								z = t5;
                                                                if(j<5){
									x = result4;
									z = t4;
								}
                                                        }else{
								x = result7;
								z = t7;
                                                                if(j<7){
									x = result6;
									z = t6;
								}
                                                        }
                                                }else{
                                                        if(j<10){
								x = result9;
								z = t9;
								if(j<9){
									x = result8;
									z = t8;
								}
                                                        }else{
								x = result10;
								z = t10;
            						}
                                                }
                                        }

					{
						s = x*y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
					}
                                        if(j<4){
                                                if(j<2){
                                                        if(j<1)
                                                                t0 = s;
                                                        else
                                                                t1=s;
                                                }else{
                                                        if(j<3)
                                                                t2=s;
                                                        else
                                                                t3=s;
                                                }
                                        }else{
                                                if(j<8){
                                                        if(j<6){
                                                                if(j<5)
                                                                        t4=s;
                                                                else
                                                                        t5=s;
                                                        }else{
                                                                if(j<7)
                                                                        t6=s;
                                                                else
                                                                        t7=s;
                                                        }
                                                }else{
                                                        if(j<10)
								if(j<9)
									t8=s;
								else
	                                                                t9=s;
                                                        else
                                                                t10=s;
                                                }
                                        }
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
					for(j=1; j < 11; j++){
						switch(j){
					case 1:{
						x = m1;
						y = t1;
					}break;
					case 2:{
						x = m2;
						y = t2;
					}break;
					case 3:{
						x = m3;
						y = t3;
					}break;
					case 4:{
						x = m4;
						y = t4;
					}break;
					case 5:{
						x = m5;
						y = t5;
					}break;
					case 6:{
						x = m6;
						y = t6;
					}break;
					case 7:{
						x = m7;
						y = t7;
					}break;
					case 8:{
						x = m8;
						y = t8;
					}break;
					case 9:{
						x = m9;
						y = t9;
					}break;
					case 10:{
						x = m10;
						y = t10;
					}break;
						}
					{
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
					}
						switch(j){
					case 1:{
						t0 = s;
					}break;
					case 2:{
						t1 = s;
					}break;
					case 3:{
						t2 = s;
					}break;
					case 4:{
						t3 = s;
					}break;
					case 5:{
						t4 = s;
					}break;
					case 6:{
						t5 = s;
					}break;
					case 7:{
						t6 = s;
					}break;
					case 8:{
						t7 = s;
					}break;
					case 9:{
						t8 = s;
					}break;
					case 10:{
						t9 = s;
					}break;
						}
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

					for(i=0; i < 11; i++)
				{
					c = 0;
						switch(i){
						case 0 :
							y = base0; break;
						case 1 :
							y = base1; break;
						case 2 :
							y = base2; break;
						case 3 :
							y = base3; break;
						case 4 :
							y = base4; break;
						case 5 :
							y = base5; break;
						case 6 :
							y = base6; break;
						case 7 :
							y = base7; break;
						case 8 :
							y = base8; break;
						case 9 :
							y = base9; break;
						case 10 :
							y = base10; break;
						}
					for(j=0; j < 11; j++){
						switch(j){
					case 0:{
						x = base0;
						z = t0;
					}break;
					case 1:{
						x = base1;
						z = t1;
					}break;
					case 2:{
						x = base2;
						z = t2;
					}break;
					case 3:{
						x = base3;
						z = t3;
					}break;
					case 4:{
						x = base4;
						z = t4;
					}break;
					case 5:{
						x = base5;
						z = t5;
					}break;
					case 6:{
						x = base6;
						z = t6;
					}break;
					case 7:{
						x = base7;
						z = t7;
					}break;
					case 8:{
						x = base8;
						z = t8;
					}break;
					case 9:{
						x = base9;
						z = t9;
					}break;
					case 10:{
						x = base10;
						z = t10;
					}break;
						}
					{
						s = x*y;
						c2 = mul_hi(x,y);
						s2 = s + z;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
					}
                                        if(j<4){
                                                if(j<2){
                                                        if(j<1)
                                                                t0 = s;
                                                        else
                                                                t1=s;
                                                }else{
                                                        if(j<3)
                                                                t2=s;
                                                        else
                                                                t3=s;
                                                }
                                        }else{
                                                if(j<8){
                                                        if(j<6){
                                                                if(j<5)
                                                                        t4=s;
                                                                else
                                                                        t5=s;
                                                        }else{
                                                                if(j<7)
                                                                        t6=s;
                                                                else
                                                                        t7=s;
                                                        }
                                                }else{
                                                        if(j<10)
								if(j<9)
									t8=s;
								else
	                                                                t9=s;
                                                        else
                                                                t10=s;
                                                }
                                        }
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
					for(j=1; j < 11; j++){
						switch(j){
					case 1:{
						x = m1;
						y = t1;
					}break;
					case 2:{
						x = m2;
						y = t2;
					}break;
					case 3:{
						x = m3;
						y = t3;
					}break;
					case 4:{
						x = m4;
						y = t4;
					}break;
					case 5:{
						x = m5;
						y = t5;
					}break;
					case 6:{
						x = m6;
						y = t6;
					}break;
					case 7:{
						x = m7;
						y = t7;
					}break;
					case 8:{
						x = m8;
						y = t8;
					}break;
					case 9:{
						x = m9;
						y = t9;
					}break;
					case 10:{
						x = m10;
						y = t10;
					}break;
						}
					{
						s = x * mf;
						c2 = mul_hi(x,mf);
						s2 = s + y;
						if(s2 < s)
							c2++;
						s = s2 + c;
						if(s < s2)
							c2++;
						c = c2;
					}
						switch(j){
					case 1:{
						t0 = s;
					}break;
					case 2:{
						t1 = s;
					}break;
					case 3:{
						t2 = s;
					}break;
					case 4:{
						t3 = s;
					}break;
					case 5:{
						t4 = s;
					}break;
					case 6:{
						t5 = s;
					}break;
					case 7:{
						t6 = s;
					}break;
					case 8:{
						t7 = s;
					}break;
					case 9:{
						t8 = s;
					}break;
					case 10:{
						t9 = s;
					}break;
						}
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

