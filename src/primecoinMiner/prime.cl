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


//TODO: this is awful on gpu. but maybe it is not a hot path
static unsigned int int_invert(unsigned int a, unsigned int nPrime)
{
   // Extended Euclidean algorithm to calculate the inverse of a in finite field defined by nPrime
   int rem0 = nPrime, rem1 = a % nPrime, rem2;
   int aux0 = 0, aux1 = 1, aux2;
   int quotient, inverse;

   while (1)
   {
      if (rem1 <= 1)
      {
         inverse = aux1;
         break;
      }

      rem2 = rem0 % rem1;
      quotient = rem0 / rem1;
      aux2 = -quotient * aux1 + aux0;

      if (rem2 <= 1)
      {
         inverse = aux2;
         break;
      }

      rem0 = rem1 % rem2;
      quotient = rem1 / rem2;
      aux0 = -quotient * aux2 + aux1;

      if (rem0 <= 1)
      {
         inverse = aux0;
         break;
      }

      rem1 = rem2 % rem0;
      quotient = rem2 / rem0;
      aux1 = -quotient * aux0 + aux2;
   }

   return (inverse + nPrime) % nPrime;
}


void SieveInit(uint32_t tid, uint32_t lSieveTarget, uint32_t lSieveBTTarget, __global sieveOutput_t *restrict output, __global sieveTemp_t *restrict temp){
	uint32_t i;

	for(i=0; i < MULT_WORDS; i++){
		D_REF(temp->nFixedInverse,i,tid) = 0;
		D_REF(temp->flatFixedInverse,i,tid) = 0;
	}

	for(i=tid; i < MULT_WORDS; i+=STRIDE){
		temp->flatCC1Multipliers[i] = 0xFFFFFFFFUL;
		temp->flatCC2Multipliers[i] = 0xFFFFFFFFUL;
	}

	if(tid==0){
		output->nBitwins = 0;
		output->nCC1s = 0;
		output->nCC2s = 0;
		temp->nSievesComplete=0;
	}

	barrier(CLK_GLOBAL_MEM_FENCE);

/*		this->nAllocatedSieveSize = nSieveSize;
      this->nSievePercentage = nSievePercentage;
      this->nSieveExtensions = nSieveExtensions;
        this->mpzHash = mpzHash;
        this->mpzFixedMultiplier = mpzFixedMultiplier;
      this->nChainLength = nTargetChainLength;
      this->nBTChainLength = nTargetBTLength;
      this->nSieveLayers = nChainLength + nSieveExtensions;
      this->nTotalPrimes = vPrimes.size();
      this->nPrimes = (uint64)nTotalPrimes * nSievePercentage / 100;
      this->nPrimeSeq = nMinPrimeSeq;
      this->nCandidateCount = 0;
      this->nCandidateMultiplier = 0;
      this->nCandidateIndex = 0;
      this->fCandidateIsExtended = false;
      this->nCandidateActiveExtension = 0;
      this->nCandidatesWords = (nSieveSize + nWordBits - 1) / nWordBits;
      this->nCandidatesBytes = nCandidatesWords * sizeof(sieve_word_t); */
}


void atomicOr(volatile __local uint32_t *mem, uint32_t bits){
	//TODO: must resolve the conflicts amongs the warp first?
#if 0	
	while((*mem & bits) != bits){
		*mem |= bits;
		while((*mem & bits) != bits){
			*mem |= bits;
		}
	}
#else
	atomic_or(mem,bits);
#endif
}

void ProcessMultiplier2(	uint32_t tid, __local uint32_t * restrict vfOut, 
				__global uint32_t *restrict vPrimes, uint32_t nPrimes,
				__global uint32_t *restrict multipliers, 
				uint32_t minPrime, uint32_t sieveBegin,
				uint32_t sieveEnd, uint32_t start){
	//Main sieve loop;
	uint32_t i,j;
	local uint32_t sem;
	local uint32_t idx[STRIDE/32];
	if(tid==0)
		sem=start;
	barrier(CLK_LOCAL_MEM_FENCE);

	//TODO: not so good here because we want warps to get new primes as they become available
	while(1){
		if(tid%32==0){
		   	i = atomic_add(&sem,32);
			idx[tid/32] = i;
		}
		mem_fence(CLK_LOCAL_MEM_FENCE);

		//barrier?
		i = idx[tid/32] + tid%32;
 		if(i >= nPrimes-1)
			break;

		if(i < minPrime) 
			continue;

	//	if(i >= 32)
	//		break;

		uint32_t prime = vPrimes[i];
		j=multipliers[i];

		//Happens during switching layers because we only do half of the sieve
	//	if (j < sieveBegin)
	//		j += (sieveBegin - j + prime - 1) / prime * prime;

		for(; j < sieveEnd; j+=prime){
		//	if(j > sieveBegin)
			atomic_or(&vfOut[(j-sieveBegin)/32],1<<(j%32));
		}
		multipliers[i] = j; //Save this for next round of sieve
	}
}

void ProcessMultiplierFast(	uint32_t tid, __local uint32_t * restrict vfOut, 
				uint32_t prime, __global uint32_t *restrict multiplier, 
				uint32_t sieveBegin, uint32_t sieveEnd){
	//Main sieve loop;
	uint32_t i,j;
	int32_t oprime = prime;

	//while(1)
	{
		j=*multiplier;
		if(j==0xFFFFFFFF)
			return;

		//TODO: When we switch to new layer, this multiplier will be wrong. It should be possible to rebase the multipliers
		//if(j<sieveBegin)
		//	j += (sieveBegin - j + prime - 1) / prime * prime;

		j+=tid*prime;
		prime*=STRIDE;

		for(; j < sieveEnd; j+=prime){
			//if(j >= sieveBegin)
//			vfOut[(j-sieveBegin)/32] |= 1<<(j%32);
			atomic_or(&vfOut[(j-sieveBegin)/32], 1<<(j%32));
		}

		//*multiplier needs to be set to the first index of prime in next sieve. not index of prime*STRIDE. 
		//some thread will have the right one
		int32_t diff = (int32_t)j - (int32_t)sieveEnd;
		if(diff < oprime){
			*multiplier = j;
		}
	}
}

void ProcessMultiplier2Fast(	uint32_t tid, __local uint32_t * restrict vfOut, 
				__global uint32_t *vPrime, uint32_t index, __global uint32_t *restrict vMultiplier, 
				uint32_t sieveBegin, uint32_t sieveEnd){
	//Main sieve loop;
	uint32_t i,j;
	index += tid/(STRIDE/2);
	uint32_t prime = vPrime[index];
	int32_t oprime = prime;

	//while(1)
	{
		j=vMultiplier[index];
		if(j==0xFFFFFFFF)
			return;

		j+=(tid%(STRIDE/2))*prime;
		prime*=STRIDE/2;

		for(; j < sieveEnd; j+=prime){
			//if(j >= sieveBegin)
//			vfOut[(j-sieveBegin)/32] |= 1<<(j%32);
			atomic_or(&vfOut[(j-sieveBegin)/32], 1<<(j%32));
		}

		//*multiplier needs to be set to the first index of prime in next sieve. not index of prime*STRIDE. 
		//some thread will have the right one
		int32_t diff = (int32_t)j - (int32_t)sieveEnd;
		if(diff < oprime){
			vMultiplier[index] = j;
		}
	}
}

void ProcessMultiplier4Fast(	uint32_t tid, __local uint32_t * restrict vfOut, 
				__global uint32_t *vPrime, uint32_t index, __global uint32_t *restrict vMultiplier, 
				uint32_t sieveBegin, uint32_t sieveEnd){
	//Main sieve loop;
	uint32_t i,j;
	index += tid/(STRIDE/4);
	uint32_t prime = vPrime[index];
	int32_t oprime = prime;

	//while(1)
	{
		j=vMultiplier[index];
		if(j==0xFFFFFFFF)
			return;

		j+=(tid%(STRIDE/4))*prime;
		prime*=STRIDE/4;

		for(; j < sieveEnd; j+=prime){
			//if(j >= sieveBegin)
//			vfOut[(j-sieveBegin)/32] |= 1<<(j%32);
			atomic_or(&vfOut[(j-sieveBegin)/32], 1<<(j%32));
		}

		//*multiplier needs to be set to the first index of prime in next sieve. not index of prime*STRIDE. 
		//some thread will have the right one
		int32_t diff = (int32_t)j - (int32_t)sieveEnd;
		if(diff < oprime){
			vMultiplier[index] = j;
		}

	}
}

void ProcessMultiplier8Fast(	uint32_t tid, __local uint32_t * restrict vfOut, 
				__global uint32_t *vPrime, uint32_t index, __global uint32_t *restrict vMultiplier, 
				uint32_t sieveBegin, uint32_t sieveEnd){
	//Main sieve loop;
	uint32_t i,j;
	index += tid/(STRIDE/8);
	uint32_t prime = vPrime[index];
	int32_t oprime = prime;

	//while(1)
	{
		j=vMultiplier[index];
		if(j==0xFFFFFFFF)
			return;

		j+=(tid%(STRIDE/8))*prime;
		prime*=STRIDE/8;

		for(; j < sieveEnd; j+=prime){
			//if(j >= sieveBegin)
//			vfOut[(j-sieveBegin)/32] |= 1<<(j%32);
			atomic_or(&vfOut[(j-sieveBegin)/32], 1<<(j%32));
		}

		//*multiplier needs to be set to the first index of prime in next sieve. not index of prime*STRIDE. 
		//some thread will have the right one
		int32_t diff = (int32_t)j - (int32_t)sieveEnd;
		if(diff < oprime){
			vMultiplier[index] = j;
		}

	}
}

void ProcessMultiplier16Fast(	uint32_t tid, __local uint32_t * restrict vfOut, 
				__global uint32_t *vPrime, uint32_t index, __global uint32_t *restrict vMultiplier, 
				uint32_t sieveBegin, uint32_t sieveEnd){
	//Main sieve loop;
	uint32_t i,j;
	index += tid/(STRIDE/16);
	uint32_t prime = vPrime[index];
	int32_t oprime = prime;

	//while(1)
	{
		j=vMultiplier[index];
		if(j==0xFFFFFFFF)
			return;

		j+=(tid%(STRIDE/16))*prime;
		prime*=STRIDE/16;

		for(; j < sieveEnd; j+=prime){
			//if(j >= sieveBegin)
//			vfOut[(j-sieveBegin)/32] |= 1<<(j%32);
			atomic_or(&vfOut[(j-sieveBegin)/32], 1<<(j%32));
		}

		//*multiplier needs to be set to the first index of prime in next sieve. not index of prime*STRIDE. 
		//some thread will have the right one
		int32_t diff = (int32_t)j - (int32_t)sieveEnd;
		if(diff < oprime){
			vMultiplier[index] = j;
		}
	}
}

void ProcessMultiplier32Fast(	uint32_t tid, __local uint32_t * restrict vfOut, 
				__global uint32_t *vPrime, uint32_t index, __global uint32_t *restrict vMultiplier, 
				uint32_t sieveBegin, uint32_t sieveEnd){
	//Main sieve loop;
	uint32_t i,j;
	index += tid/(STRIDE/32);
	uint32_t prime = vPrime[index];
	int32_t oprime = prime;

	//while(1)
	{
		j=vMultiplier[index];
		if(j==0xFFFFFFFF)
			return;

		j+=(tid%(STRIDE/32))*prime;
		prime*=STRIDE/32;

		for(; j < sieveEnd; j+=prime){
			//if(j >= sieveBegin)
//			vfOut[(j-sieveBegin)/32] |= 1<<(j%32);
			atomic_or(&vfOut[(j-sieveBegin)/32], 1<<(j%32));
		}

		//*multiplier needs to be set to the first index of prime in next sieve. not index of prime*STRIDE. 
		//some thread will have the right one
		int32_t diff = (int32_t)j - (int32_t)sieveEnd;
		if(diff < oprime){
			vMultiplier[index] = j;
		}
	}
}

void ProcessMultiplier64Fast(	uint32_t tid, __local uint32_t * restrict vfOut, 
				__global uint32_t *vPrime, uint32_t index, __global uint32_t *restrict vMultiplier, 
				uint32_t sieveBegin, uint32_t sieveEnd){
	//Main sieve loop;
	uint32_t i,j;
	index += tid/(STRIDE/64);
	uint32_t prime = vPrime[index];
	int32_t oprime = prime;

	//while(1)
	{
		j=vMultiplier[index];
		if(j==0xFFFFFFFF)
			return;

		j+=(tid%(STRIDE/64))*prime;
		prime*=STRIDE/64;

		for(; j < sieveEnd; j+=prime){
			//if(j >= sieveBegin)
//			vfOut[(j-sieveBegin)/32] |= 1<<(j%32);
			atomic_or(&vfOut[(j-sieveBegin)/32], 1<<(j%32));
		}

		//*multiplier needs to be set to the first index of prime in next sieve. not index of prime*STRIDE. 
		//some thread will have the right one
		int32_t diff = (int32_t)j - (int32_t)sieveEnd;
		if(diff < oprime){
			vMultiplier[index] = j;
		}
	}
}

void ProcessMultiplier128Fast(	uint32_t tid, __local uint32_t * restrict vfOut, 
				__global uint32_t *vPrime, uint32_t index, __global uint32_t *restrict vMultiplier, 
				uint32_t sieveBegin, uint32_t sieveEnd){
	//Main sieve loop;
	uint32_t i,j;
	index += tid/(STRIDE/128);
	uint32_t prime = vPrime[index];
	int32_t oprime = prime;

	//while(1)
	{
		j=vMultiplier[index];
		if(j==0xFFFFFFFF)
			return;

		j+=(tid%(STRIDE/128))*prime;
		prime*=STRIDE/128;

		for(; j < sieveEnd; j+=prime){
			//if(j >= sieveBegin)
//			vfOut[(j-sieveBegin)/32] |= 1<<(j%32);
			atomic_or(&vfOut[(j-sieveBegin)/32], 1<<(j%32));
		}

		//*multiplier needs to be set to the first index of prime in next sieve. not index of prime*STRIDE. 
		//some thread will have the right one
		int32_t diff = (int32_t)j - (int32_t)sieveEnd;
		if(diff < oprime){
			vMultiplier[index] = j;
		}
	}
}


uint32_t morton(uint32_t x)
{
    x = x & 0x55555555;
    x = (x | (x >> 1)) & 0x33333333;
    x = (x | (x >> 2)) & 0x0F0F0F0F;
    x = (x | (x >> 4)) & 0x00FF00FF;
    x = (x | (x >> 8)) & 0x0000FFFF;
    return x;
}

__constant uint32_t fPrimes[] = {1,2,3,5,7,11,13,17,19,23,29,31,37};
 
void ProcessMultiplier3(uint32_t tid, __local uint32_t* temp, __global uint32_t * restrict vfBase, __global uint32_t *restrict vPrimes, uint32_t nPrimes, __global uint32_t *restrict multipliers, uint32_t minPrime, uint32_t layerSeq){
	uint32_t j,k,l,i=0;
#if 1	
	barrier( CLK_GLOBAL_MEM_FENCE );

	__global uint32_t *vfOut = vfBase + layerSeq*SIEVE_WORDS;

#if 0
	if(layerSeq!=0 && layerSeq < NCHAIN_LENGTH){
		__global uint32_t *vfIn = vfBase + (layerSeq-1)*SIEVE_WORDS;
		for(i=0; i < SIEVE_WORDS/2; i+=32*1024/4){
			for(j=0; j < 32*1024/4; j+=STRIDE){
				uint32_t v1 = vfIn[(i+j+tid) * 2];
				uint32_t v2 = vfIn[(i+j+tid) * 2 + 1];
				v1 = morton(v1);
				v1 |= morton(v2) << 16;
				vfOut[i+j+tid] = v1;
			}
		}
	}
	if(layerSeq>=NCHAIN_LENGTH){
		i=SIEVE_WORDS/2;
	}

#endif
	for(; i < SIEVE_WORDS; i+=32*1024/4){
		//Clear the memory
		for(j=0; j < 32*1024/4; j+=STRIDE){
			temp[j+tid] = 0;
		}
		barrier( CLK_LOCAL_MEM_FENCE ); //Make sure the data exists before we start touching OPM. 

		//TODO: why is it that multiplier table starts at 13? - Aha, sieve cannot work if origin mod p == 0. 
		//so at least all primes in the primorial factor are out. 

		//Primes that need for speed - Split up the sieve a little. Different primes work
		//best with different number of threads dedicated to each. Eventually the density is
		//so low we move to 1 thread per prime
#if 1
		k=12; //***************************************Hard coded for PRIMORIAL 41!!!!!! **************************************
		for(l=0; l < FAST_PRIMES; l++){
			ProcessMultiplierFast(tid, temp, fPrimes[k], multipliers+k,i*32,i*32+32*1024*8);
			k++;
		//	barrier(CLK_LOCAL_MEM_FENCE);
		}
		for(l=0; l < FAST_PRIMES2; l++){
			ProcessMultiplier2Fast(tid, temp, vPrimes, k, multipliers,i*32,i*32+32*1024*8);
			k+=2;
			//barrier(CLK_LOCAL_MEM_FENCE);
		}
#endif
#if 1
		for(l=0; l < FAST_PRIMES4; l++){
			ProcessMultiplier4Fast(tid, temp, vPrimes, k, multipliers,i*32,i*32+32*1024*8);
			k+=4;
			//barrier(CLK_LOCAL_MEM_FENCE);
		}

		for(l=0; l < FAST_PRIMES8; l++){
			ProcessMultiplier8Fast(tid, temp, vPrimes, k, multipliers,i*32,i*32+32*1024*8);
			k+=8;
			//barrier(CLK_LOCAL_MEM_FENCE);
		}

		for(l=0; l < FAST_PRIMES16; l++){
			ProcessMultiplier16Fast(tid, temp, vPrimes, k, multipliers,i*32,i*32+32*1024*8);
			k+=16;
			//barrier(CLK_LOCAL_MEM_FENCE);
		}

		for(l=0; l < FAST_PRIMES32; l++){
			ProcessMultiplier32Fast(tid, temp, vPrimes, k, multipliers,i*32,i*32+32*1024*8);
			k+=32;
			//barrier(CLK_LOCAL_MEM_FENCE);
		}

		for(l=0; l < FAST_PRIMES64; l++){
			ProcessMultiplier64Fast(tid, temp, vPrimes, k, multipliers,i*32,i*32+32*1024*8);
			k+=64;
			//barrier(CLK_LOCAL_MEM_FENCE);
		}

		for(l=0; l < FAST_PRIMES128; l++){
			ProcessMultiplier128Fast(tid, temp, vPrimes, k, multipliers,i*32,i*32+32*1024*8);
			k+=128;
			//barrier(CLK_LOCAL_MEM_FENCE);
		}


		ProcessMultiplier2(tid, temp, vPrimes, nPrimes, multipliers,minPrime,i*32,i*32+32*1024*8,k);
#endif
		barrier(CLK_LOCAL_MEM_FENCE);
		for(j=tid; j < 32*1024/4; j+=STRIDE){
			vfOut[i+j] = temp[j];
		}
	}
#endif
}

void ProcessMultiplier4(uint32_t tid, __local uint32_t* temp, __global uint32_t * restrict vfOut, __global uint32_t *restrict vPrimes, uint32_t nPrimes, __global uint32_t *restrict multipliers, uint32_t minPrime, uint32_t layerSeq, uint32_t nSieve){
	ProcessMultiplier3(tid,temp,vfOut,vPrimes,nPrimes,multipliers,minPrime, layerSeq);
}
 
void transpose(uint32_t tid, __local uint32_t *temp, __global uint32_t *in, __global uint32_t *out, uint32_t nSize){
	//Use a stride of 33 in local memory because it stops bank conflicts. This causes all kinds of nasty
	//math for the read index, but oh well. better than 32 cycles for local memory read

	uint32_t i,j;
	

	//Disable half the threads. Sucks, but not enough local memory

	for(i=0; i < nSize; i+=32){ 
		for(j=0; j < 32 && tid < 256; j++){
			temp[tid*33 + j] = D_REF(in,i+j,tid); //32 words for each of 256 threads
		}

		barrier( CLK_LOCAL_MEM_FENCE ); //Make sure the data exists before we start touching OPM. 

		for(j=0; j < 32 && tid < 256; j++){
			uint32_t rdidx = j*8+(tid/32);
			out[rdidx*nSize+(tid%32)+i] = temp[rdidx*33+(tid%32)];
		}
	}
}

void flattenMultipliers(uint32_t tid, __local uint32_t *temp, __global uint32_t *in, __global uint32_t *out, int nPrimes){
	//flattenMultipliers can only handle half the array, so we pump it up twice
	uint32_t i;
	uint32_t rounds = STRIDE/256;
	for(i=0; i < rounds; i++){
		transpose(tid,temp,&in[i*STRIDE/rounds],&out[i*STRIDE*MULT_SIZE/rounds],MULT_SIZE);
	}
}

void SpinSieve(uint32_t tid, __global sieveTemp_t *restrict temp, __global uint32_t *restrict vPrimes, uint32_t nPrimes, uint32_t minPrime, uint32_t nLayer, uint32_t nSieve){
	uint32_t i;
	for(i=minPrime+tid; i < nPrimes; i+=STRIDE){
		uint32_t nPrime = vPrimes[i];
		uint32_t findex = i+nSieve*MULT_WORDS;
		uint32_t nFixedInverse = temp->flatFixedInverse[findex];
		if(nFixedInverse==0){
			temp->flatCC1Multipliers[i] = 0xFFFFFFFF;
			temp->flatCC2Multipliers[i] = 0xFFFFFFFF;
			continue;
		}
		uint32_t nTwoInverse = (nPrime + 1) / 2;

		// Weave the sieve for the prime
		{
			// Find the first number that's divisible by this prime
			uint32_t cc1 = nFixedInverse;
		 	uint32_t cc2 = nPrime - nFixedInverse;

#if 0 //doing all whole layers now
			if(nLayer){
				if(cc1 < SIEVE_SIZE/2)
					cc1 += ((SIEVE_SIZE/2) - cc1 + nPrime - 1) / nPrime * nPrime;
				if(cc2 < SIEVE_SIZE/2)
					cc2 += ((SIEVE_SIZE/2) - cc2 + nPrime - 1) / nPrime * nPrime;
			}
#endif

			// Make sure they are even
			if(cc1%2==1) cc1 += nPrime;
			if(cc2%2==1) cc2 += nPrime;

			temp->flatCC1Multipliers[i] = cc1;
		 	temp->flatCC2Multipliers[i] = cc2;

			// For next number in chain
			nFixedInverse = (uint64_t) cc1 *nTwoInverse % nPrime;
			temp->flatFixedInverse[findex] = nFixedInverse;
		}
	}
}

//This code takes 40ms
void createList(uint32_t tid, __global uint32_t *restrict input, __local volatile uint32_t *atom, 
				__global result_t *restrict output, uint32_t type, uint32_t hash){
	uint32_t i,j;
	for(i=1; i < EXTENSION_LAYERS+1; i++){
		j=tid;
		for(j = SIEVE_WORDS/2; j < SIEVE_WORDS; j+=STRIDE){
			uint32_t v = ~input[j+i*SIEVE_WORDS];
			uint32_t idx;		
			while((idx = clz(v)) != 32){
				idx = 31-idx;
				v &= ~(1<<(idx));
				if(idx==0 && j==0)
					continue;
				uint32_t loc = atomic_add(atom,1);
				//write the prime somewhere	
				output->mults[loc] = ( (2*(j*32+idx))+1)<<i;
				output->hash[loc] = hash;
				output->type[loc] = type;			
			}
		}
	}

	barrier( CLK_LOCAL_MEM_FENCE );
}

void combineResults(uint32_t tid, __local volatile uint32_t *atom, 
				__local volatile uint32_t *ninputsptr, 
				__global bresult_t *restrict output,
				__global result_t *restrict input){
	barrier( CLK_LOCAL_MEM_FENCE );
	uint32_t ninput = *ninputsptr;
	uint32_t offset = *atom;
	uint32_t i;
	for(i=tid; i < ninput; i+=STRIDE){
			output->mults[i+offset] = input->mults[i];
			output->hash[i+offset] = input->hash[i];
			output->type[i+offset] = input->type[i];
	}

	if(ninput%32==0){
		*atom = ninput+offset;
		barrier( CLK_LOCAL_MEM_FENCE );
		return;
	}

	for(i=tid+ninput; i < ninput+32-(ninput%32); i+=STRIDE){
		output->type[i+offset] = -1;
	}
	if(tid==0)
		*atom = ninput+offset+32 - (ninput%32);
	barrier( CLK_LOCAL_MEM_FENCE );
}

//This code takes 65ms to execute. 
void combineSieves(uint32_t tid, uint32_t layer, __global sieveTemp_t *temp, 
							uint32_t nBiTwinCC1Layers, uint32_t nBiTwinCC2Layers){
		uint32_t j,k;

		//TODO: this would have to be regenerated for different numbers of layers and extensions
		for(k=tid; k < SIEVE_WORDS; k+=STRIDE){
			uint32_t cc1_0 = temp->vfCC1[0*SIEVE_WORDS+k];
			uint32_t cc1_1 = temp->vfCC1[1*SIEVE_WORDS+k];
			uint32_t cc1_2 = temp->vfCC1[2*SIEVE_WORDS+k];
			uint32_t cc1_3 = temp->vfCC1[3*SIEVE_WORDS+k];
			uint32_t cc1_4 = temp->vfCC1[4*SIEVE_WORDS+k];
			uint32_t cc1_5 = temp->vfCC1[5*SIEVE_WORDS+k];
			uint32_t cc1_6 = temp->vfCC1[6*SIEVE_WORDS+k];
			uint32_t cc1_7 = temp->vfCC1[7*SIEVE_WORDS+k];
			uint32_t cc1_8 = temp->vfCC1[8*SIEVE_WORDS+k];
			uint32_t cc1_9 = temp->vfCC1[9*SIEVE_WORDS+k];
			uint32_t cc1_10 = temp->vfCC1[10*SIEVE_WORDS+k];
			uint32_t cc1_11 = temp->vfCC1[11*SIEVE_WORDS+k];
			uint32_t cc1_12 = temp->vfCC1[12*SIEVE_WORDS+k];
			uint32_t cc1_13 = temp->vfCC1[13*SIEVE_WORDS+k];
			uint32_t cc1_14 = temp->vfCC1[14*SIEVE_WORDS+k];
			uint32_t cc1_15 = temp->vfCC1[15*SIEVE_WORDS+k];
			uint32_t cc1_16 = temp->vfCC1[16*SIEVE_WORDS+k];
			uint32_t cc1_17 = temp->vfCC1[17*SIEVE_WORDS+k];
			uint32_t cc1_18 = temp->vfCC1[18*SIEVE_WORDS+k];
			uint32_t cc1_19 = temp->vfCC1[19*SIEVE_WORDS+k];

			uint32_t cc2_0 = temp->vfCC2[0*SIEVE_WORDS+k];
			uint32_t cc2_1 = temp->vfCC2[1*SIEVE_WORDS+k];
			uint32_t cc2_2 = temp->vfCC2[2*SIEVE_WORDS+k];
			uint32_t cc2_3 = temp->vfCC2[3*SIEVE_WORDS+k];
			uint32_t cc2_4 = temp->vfCC2[4*SIEVE_WORDS+k];
			uint32_t cc2_5 = temp->vfCC2[5*SIEVE_WORDS+k];
			uint32_t cc2_6 = temp->vfCC2[6*SIEVE_WORDS+k];
			uint32_t cc2_7 = temp->vfCC2[7*SIEVE_WORDS+k];
			uint32_t cc2_8 = temp->vfCC2[8*SIEVE_WORDS+k];
			uint32_t cc2_9 = temp->vfCC2[9*SIEVE_WORDS+k];
			uint32_t cc2_10 = temp->vfCC2[10*SIEVE_WORDS+k];
			uint32_t cc2_11 = temp->vfCC2[11*SIEVE_WORDS+k];
			uint32_t cc2_12 = temp->vfCC2[12*SIEVE_WORDS+k];
			uint32_t cc2_13 = temp->vfCC2[13*SIEVE_WORDS+k];
			uint32_t cc2_14 = temp->vfCC2[14*SIEVE_WORDS+k];
			uint32_t cc2_15 = temp->vfCC2[15*SIEVE_WORDS+k];
			uint32_t cc2_16 = temp->vfCC2[16*SIEVE_WORDS+k];
			uint32_t cc2_17 = temp->vfCC2[17*SIEVE_WORDS+k];
			uint32_t cc2_18 = temp->vfCC2[18*SIEVE_WORDS+k];
			uint32_t cc2_19 = temp->vfCC2[19*SIEVE_WORDS+k];

			{
				uint32_t v=0;
				v |= cc1_0;
				v |= cc1_1;
				v |= cc1_2;
				v |= cc1_3;
				v |= cc1_4;
				v |= cc1_5;
				v |= cc1_6;
				v |= cc1_7;
				v |= cc1_8;
				v |= cc1_9;
				temp->vfExtCC1[k+0*SIEVE_WORDS] = v;
			}

			{
				uint32_t v=0;
				v |= cc2_0;
				v |= cc2_1;
				v |= cc2_2;
				v |= cc2_3;
				v |= cc2_4;
				v |= cc2_5;
				v |= cc2_6;
				v |= cc2_7;
				v |= cc2_8;
				v |= cc2_9;
				temp->vfExtCC2[k+0*SIEVE_WORDS] = v;
			}

			{
				uint32_t v=0;
				v |= cc1_0;
				v |= cc1_1;
				v |= cc1_2;
				v |= cc1_3;
				v |= cc1_4;
				v |= cc2_0;
				v |= cc2_1;
				v |= cc2_2;
				v |= cc2_3;
				v |= cc2_4;
				temp->vfExtBitwin[k+0*SIEVE_WORDS] = v;
			}

			{
				uint32_t v=0;
				v |= cc1_1;
				v |= cc1_2;
				v |= cc1_3;
				v |= cc1_4;
				v |= cc1_5;
				v |= cc1_6;
				v |= cc1_7;
				v |= cc1_8;
				v |= cc1_9;
				v |= cc1_10;
				temp->vfExtCC1[k+1*SIEVE_WORDS] = v;
			}

			{
				uint32_t v=0;
				v |= cc2_1;
				v |= cc2_2;
				v |= cc2_3;
				v |= cc2_4;
				v |= cc2_5;
				v |= cc2_6;
				v |= cc2_7;
				v |= cc2_8;
				v |= cc2_9;
				v |= cc2_10;
				temp->vfExtCC2[k+1*SIEVE_WORDS] = v;
			}

			{
				uint32_t v=0;
				v |= cc1_1;
				v |= cc1_2;
				v |= cc1_3;
				v |= cc1_4;
				v |= cc1_5;
				v |= cc2_1;
				v |= cc2_2;
				v |= cc2_3;
				v |= cc2_4;
				v |= cc2_5;
				temp->vfExtBitwin[k+1*SIEVE_WORDS] = v;
			}

			{
				uint32_t v=0;
				v |= cc1_2;
				v |= cc1_3;
				v |= cc1_4;
				v |= cc1_5;
				v |= cc1_6;
				v |= cc1_7;
				v |= cc1_8;
				v |= cc1_9;
				v |= cc1_10;
				v |= cc1_11;
				temp->vfExtCC1[k+2*SIEVE_WORDS] = v;
			}

			{
				uint32_t v=0;
				v |= cc2_2;
				v |= cc2_3;
				v |= cc2_4;
				v |= cc2_5;
				v |= cc2_6;
				v |= cc2_7;
				v |= cc2_8;
				v |= cc2_9;
				v |= cc2_10;
				v |= cc2_11;
				temp->vfExtCC2[k+2*SIEVE_WORDS] = v;
			}

			{
				uint32_t v=0;
				v |= cc1_2;
				v |= cc1_3;
				v |= cc1_4;
				v |= cc1_5;
				v |= cc1_6;
				v |= cc2_2;
				v |= cc2_3;
				v |= cc2_4;
				v |= cc2_5;
				v |= cc2_6;
				temp->vfExtBitwin[k+2*SIEVE_WORDS] = v;
			}

			{
				uint32_t v=0;
				v |= cc1_3;
				v |= cc1_4;
				v |= cc1_5;
				v |= cc1_6;
				v |= cc1_7;
				v |= cc1_8;
				v |= cc1_9;
				v |= cc1_10;
				v |= cc1_11;
				v |= cc1_12;
				temp->vfExtCC1[k+3*SIEVE_WORDS] = v;
			}

			{
				uint32_t v=0;
				v |= cc2_3;
				v |= cc2_4;
				v |= cc2_5;
				v |= cc2_6;
				v |= cc2_7;
				v |= cc2_8;
				v |= cc2_9;
				v |= cc2_10;
				v |= cc2_11;
				v |= cc2_12;
				temp->vfExtCC2[k+3*SIEVE_WORDS] = v;
			}

			{
				uint32_t v=0;
				v |= cc1_3;
				v |= cc1_4;
				v |= cc1_5;
				v |= cc1_6;
				v |= cc1_7;
				v |= cc2_3;
				v |= cc2_4;
				v |= cc2_5;
				v |= cc2_6;
				v |= cc2_7;
				temp->vfExtBitwin[k+3*SIEVE_WORDS] = v;
			}

			{
				uint32_t v=0;
				v |= cc1_4;
				v |= cc1_5;
				v |= cc1_6;
				v |= cc1_7;
				v |= cc1_8;
				v |= cc1_9;
				v |= cc1_10;
				v |= cc1_11;
				v |= cc1_12;
				v |= cc1_13;
				temp->vfExtCC1[k+4*SIEVE_WORDS] = v;
			}

			{
				uint32_t v=0;
				v |= cc2_4;
				v |= cc2_5;
				v |= cc2_6;
				v |= cc2_7;
				v |= cc2_8;
				v |= cc2_9;
				v |= cc2_10;
				v |= cc2_11;
				v |= cc2_12;
				v |= cc2_13;
				temp->vfExtCC2[k+4*SIEVE_WORDS] = v;
			}

			{
				uint32_t v=0;
				v |= cc1_4;
				v |= cc1_5;
				v |= cc1_6;
				v |= cc1_7;
				v |= cc1_8;
				v |= cc2_4;
				v |= cc2_5;
				v |= cc2_6;
				v |= cc2_7;
				v |= cc2_8;
				temp->vfExtBitwin[k+4*SIEVE_WORDS] = v;
			}

			{
				uint32_t v=0;
				v |= cc1_5;
				v |= cc1_6;
				v |= cc1_7;
				v |= cc1_8;
				v |= cc1_9;
				v |= cc1_10;
				v |= cc1_11;
				v |= cc1_12;
				v |= cc1_13;
				v |= cc1_14;
				temp->vfExtCC1[k+5*SIEVE_WORDS] = v;
			}

			{
				uint32_t v=0;
				v |= cc2_5;
				v |= cc2_6;
				v |= cc2_7;
				v |= cc2_8;
				v |= cc2_9;
				v |= cc2_10;
				v |= cc2_11;
				v |= cc2_12;
				v |= cc2_13;
				v |= cc2_14;
				temp->vfExtCC2[k+5*SIEVE_WORDS] = v;
			}

			{
				uint32_t v=0;
				v |= cc1_5;
				v |= cc1_6;
				v |= cc1_7;
				v |= cc1_8;
				v |= cc1_9;
				v |= cc2_5;
				v |= cc2_6;
				v |= cc2_7;
				v |= cc2_8;
				v |= cc2_9;
				temp->vfExtBitwin[k+5*SIEVE_WORDS] = v;
			}

			{
				uint32_t v=0;
				v |= cc1_6;
				v |= cc1_7;
				v |= cc1_8;
				v |= cc1_9;
				v |= cc1_10;
				v |= cc1_11;
				v |= cc1_12;
				v |= cc1_13;
				v |= cc1_14;
				v |= cc1_15;
				temp->vfExtCC1[k+6*SIEVE_WORDS] = v;
			}

			{
				uint32_t v=0;
				v |= cc2_6;
				v |= cc2_7;
				v |= cc2_8;
				v |= cc2_9;
				v |= cc2_10;
				v |= cc2_11;
				v |= cc2_12;
				v |= cc2_13;
				v |= cc2_14;
				v |= cc2_15;
				temp->vfExtCC2[k+6*SIEVE_WORDS] = v;
			}

			{
				uint32_t v=0;
				v |= cc1_6;
				v |= cc1_7;
				v |= cc1_8;
				v |= cc1_9;
				v |= cc1_10;
				v |= cc2_6;
				v |= cc2_7;
				v |= cc2_8;
				v |= cc2_9;
				v |= cc2_10;
				temp->vfExtBitwin[k+6*SIEVE_WORDS] = v;
			}

			{
				uint32_t v=0;
				v |= cc1_7;
				v |= cc1_8;
				v |= cc1_9;
				v |= cc1_10;
				v |= cc1_11;
				v |= cc1_12;
				v |= cc1_13;
				v |= cc1_14;
				v |= cc1_15;
				v |= cc1_16;
				temp->vfExtCC1[k+7*SIEVE_WORDS] = v;
			}

			{
				uint32_t v=0;
				v |= cc2_7;
				v |= cc2_8;
				v |= cc2_9;
				v |= cc2_10;
				v |= cc2_11;
				v |= cc2_12;
				v |= cc2_13;
				v |= cc2_14;
				v |= cc2_15;
				v |= cc2_16;
				temp->vfExtCC2[k+7*SIEVE_WORDS] = v;
			}

			{
				uint32_t v=0;
				v |= cc1_7;
				v |= cc1_8;
				v |= cc1_9;
				v |= cc1_10;
				v |= cc1_11;
				v |= cc2_7;
				v |= cc2_8;
				v |= cc2_9;
				v |= cc2_10;
				v |= cc2_11;
				temp->vfExtBitwin[k+7*SIEVE_WORDS] = v;
			}

			{
				uint32_t v=0;
				v |= cc1_8;
				v |= cc1_9;
				v |= cc1_10;
				v |= cc1_11;
				v |= cc1_12;
				v |= cc1_13;
				v |= cc1_14;
				v |= cc1_15;
				v |= cc1_16;
				v |= cc1_17;
				temp->vfExtCC1[k+8*SIEVE_WORDS] = v;
			}

			{
				uint32_t v=0;
				v |= cc2_8;
				v |= cc2_9;
				v |= cc2_10;
				v |= cc2_11;
				v |= cc2_12;
				v |= cc2_13;
				v |= cc2_14;
				v |= cc2_15;
				v |= cc2_16;
				v |= cc2_17;
				temp->vfExtCC2[k+8*SIEVE_WORDS] = v;
			}

			{
				uint32_t v=0;
				v |= cc1_8;
				v |= cc1_9;
				v |= cc1_10;
				v |= cc1_11;
				v |= cc1_12;
				v |= cc2_8;
				v |= cc2_9;
				v |= cc2_10;
				v |= cc2_11;
				v |= cc2_12;
				temp->vfExtBitwin[k+8*SIEVE_WORDS] = v;
			}

			{
				uint32_t v=0;
				v |= cc1_9;
				v |= cc1_10;
				v |= cc1_11;
				v |= cc1_12;
				v |= cc1_13;
				v |= cc1_14;
				v |= cc1_15;
				v |= cc1_16;
				v |= cc1_17;
				v |= cc1_18;
				temp->vfExtCC1[k+9*SIEVE_WORDS] = v;
			}

			{
				uint32_t v=0;
				v |= cc2_9;
				v |= cc2_10;
				v |= cc2_11;
				v |= cc2_12;
				v |= cc2_13;
				v |= cc2_14;
				v |= cc2_15;
				v |= cc2_16;
				v |= cc2_17;
				v |= cc2_18;
				temp->vfExtCC2[k+9*SIEVE_WORDS] = v;
			}

			{
				uint32_t v=0;
				v |= cc1_9;
				v |= cc1_10;
				v |= cc1_11;
				v |= cc1_12;
				v |= cc1_13;
				v |= cc2_9;
				v |= cc2_10;
				v |= cc2_11;
				v |= cc2_12;
				v |= cc2_13;
				temp->vfExtBitwin[k+9*SIEVE_WORDS] = v;
			}

			{
				uint32_t v=0;
				v |= cc1_10;
				v |= cc1_11;
				v |= cc1_12;
				v |= cc1_13;
				v |= cc1_14;
				v |= cc1_15;
				v |= cc1_16;
				v |= cc1_17;
				v |= cc1_18;
				v |= cc1_19;
				temp->vfExtCC1[k+10*SIEVE_WORDS] = v;
			}

			{
				uint32_t v=0;
				v |= cc2_10;
				v |= cc2_11;
				v |= cc2_12;
				v |= cc2_13;
				v |= cc2_14;
				v |= cc2_15;
				v |= cc2_16;
				v |= cc2_17;
				v |= cc2_18;
				v |= cc2_19;
				temp->vfExtCC2[k+10*SIEVE_WORDS] = v;
			}

			{
				uint32_t v=0;
				v |= cc1_10;
				v |= cc1_11;
				v |= cc1_12;
				v |= cc1_13;
				v |= cc1_14;
				v |= cc2_10;
				v |= cc2_11;
				v |= cc2_12;
				v |= cc2_13;
				v |= cc2_14;
				temp->vfExtBitwin[k+10*SIEVE_WORDS] = v;
			}

		}
}
 
void WeaveInit(uint32_t tid, __global primecoinInput_t *restrict input, __global uint32_t *restrict vPrimes, __global sieveOutput_t *restrict output,
			 __global sieveTemp_t *restrict temp, __global shaOutput_t *restrict sha){
	mulBNBNS(tid,&temp->mpzFixedFactor,&sha->mpzHash,&input->mpzFixedMultiplier);


	local uint32_t ltemp[256*33];

	uint32_t nCombinedEndSeq = input->primeSeq;
   	uint32_t nFixedFactorCombinedMod = 0;

	//Generate the modulo multipliers table
	uint32_t i;
#if 1
	for(i=input->primeSeq; i < input->nPrimes-1; i++){

      		uint32_t nPrime = vPrimes[i];
      		if (i >= nCombinedEndSeq)
      		{
         		// Combine multiple primes to produce a big divisor
         		uint32_t nPrimeCombined = 1;
         		while (nPrimeCombined < UINT_MAX / vPrimes[nCombinedEndSeq]  )
         		{
            			nPrimeCombined *= vPrimes[nCombinedEndSeq];
            			nCombinedEndSeq++;
         		}
 
            		nFixedFactorCombinedMod = modBN(tid,&temp->mpzQuot,&temp->mpzFixedFactor, nPrimeCombined);
			temp->nPrimeCombined[tid] = nPrimeCombined;
		}
		temp->nFixedFactorCombinedMod[tid] = nFixedFactorCombinedMod;
	
#if 1 
		uint32_t nFixedFactorMod = nFixedFactorCombinedMod % nPrime;

		if (nFixedFactorMod == 0) {
			// Nothing in the sieve is divisible by this prime
			continue;
		}

		// Find the modulo inverse of fixed factor
		uint32_t nFixedInverse = int_invert (nFixedFactorMod, nPrime);
		if (!nFixedInverse){
			//TODO: error
			//("CSieveOfEratosthenes::Weave(): int_invert of fixed factor failed for prime #%u=%u",
		}

		D_REF(temp->nFixedInverse,i,tid) = nFixedInverse;
#endif
  	}
#endif  
	flattenMultipliers(tid, ltemp, temp->nFixedInverse, temp->flatFixedInverse, input->nPrimes);

	barrier( CLK_GLOBAL_MEM_FENCE );
}

void WeavePart(uint32_t tid, __global primecoinInput_t *restrict input, __global uint32_t *restrict vPrimes, __global sieveOutput_t *restrict output,
			 __global sieveTemp_t *restrict temp, __global shaOutput_t *restrict sha){

	local uint32_t ltemp[256*33];
	uint32_t i;

	// Calculate the number of CC1 and CC2 layers needed for BiTwin candidates
	const unsigned int nBiTwinCC1Layers = (input->lSieveBTTarget + 1) / 2;
	const unsigned int nBiTwinCC2Layers = input->lSieveBTTarget / 2;

	volatile local uint32_t* ncc1s = &ltemp[256*33-1];
	volatile local uint32_t* ncc2s = &ltemp[256*33-2];
	volatile local uint32_t* nbitwins = &ltemp[256*33-3];

	//Clear the accumulators
	output->nCC1s = 0;
	output->nCC2s = 0;
	output->nBitwins = 0;

	if(tid==0){
		*ncc1s = 0;//output->nCC1s;
		*ncc2s = 0;//output->nCC2s;
		*nbitwins = 0;//output->nBitwins;
	}
	barrier( CLK_LOCAL_MEM_FENCE );

	for(i=temp->nSievesComplete; i < (STRIDE/ROUNDS) + temp->nSievesComplete; i++){//STRIDE; i++){ 
		uint32_t j,k;

		//Business time -- Calculate all layers of this particular sieve at once
		for(j=0; j < NCHAIN_LENGTH+EXTENSION_LAYERS; j++){
			SpinSieve(tid,temp,vPrimes,input->nPrimes,input->primeSeq,j,i);	

			ProcessMultiplier4(tid, ltemp, temp->vfCC1, vPrimes, input->nPrimes, temp->flatCC1Multipliers, input->primeSeq, j, i); 
			ProcessMultiplier4(tid, ltemp, temp->vfCC2, vPrimes, input->nPrimes, temp->flatCC2Multipliers, input->primeSeq, j, i); 
		}
		barrier( CLK_GLOBAL_MEM_FENCE );
#if 1	

#if 1
		combineSieves(tid,0,temp,nBiTwinCC1Layers,nBiTwinCC2Layers);

#endif
#if 1
		//Compute the results
		barrier( CLK_GLOBAL_MEM_FENCE );
		createList(tid,temp->vfExtCC1,ncc1s,&temp->cc1s,CC1_TYPE, i);
		createList(tid,temp->vfExtCC2,ncc2s,&temp->cc2s,CC2_TYPE, i);
		createList(tid,temp->vfExtBitwin,nbitwins,&temp->bitwins,BITWIN_TYPE, i);
#endif
#endif
	}

	if(tid==0){
		output->nCC1s = *ncc1s;
		output->nCC2s = *ncc2s;
		output->nBitwins = *nbitwins;
		output->nLayer = 0;
		output->nTested = 0;
		temp->nSievesComplete += (STRIDE/ROUNDS);
	}
}

void WeaveComplete(uint32_t tid, __global primecoinInput_t *restrict input, __global uint32_t *restrict vPrimes, __global sieveOutput_t *restrict output,
			 __global sieveTemp_t *restrict temp, __global shaOutput_t *restrict sha){

	local uint32_t ltemp[4];
	volatile local uint32_t* ncc1s = &ltemp[0];
	volatile local uint32_t* ncc2s = &ltemp[1];
	volatile local uint32_t* nbitwins = &ltemp[2];
	volatile local uint32_t* ncombined = &ltemp[3];

	if(tid==0){
		*ncc1s = output->nCC1s;
		*ncc2s = output->nCC2s;
		*nbitwins = output->nBitwins;
		*ncombined = 0;
	}
	barrier( CLK_LOCAL_MEM_FENCE );

	//Write the counts into results
	if(tid==0){
		if(*ncc1s%32==0)
			output->nCC1s = *ncc1s;
		else
			output->nCC1s = *ncc1s + 32 - *ncc1s%32;
		if(*ncc2s%32==0)
			output->nCC2s = *ncc2s + output->nCC1s;
		else
			output->nCC2s = *ncc2s + 32 - *ncc2s%32 + output->nCC1s;
		if(*nbitwins%32==0)
			output->nBitwins = *nbitwins + output->nCC2s;
		else
			output->nBitwins = *nbitwins + 32 - *nbitwins%32 + output->nCC2s;
	}

	barrier( CLK_GLOBAL_MEM_FENCE );
	combineResults(tid,ncombined,ncc1s,&output->results,&temp->cc1s);
	combineResults(tid,ncombined,ncc2s,&output->results,&temp->cc2s);
	combineResults(tid,ncombined,nbitwins,&output->results,&temp->bitwins);
}     

void MineProbablePrimeChain(uint32_t tid, __global primecoinInput_t *restrict input, __global sieveOutput_t *restrict output, 
		__global sieveTemp_t *restrict temp, __global uint32_t *restrict vPrimes, __global shaOutput_t *restrict sha ){
	SieveInit(tid,input->lSieveTarget,input->lSieveBTTarget, output, temp);

	WeaveInit(tid,input,vPrimes,output,temp,sha);
}           

void MineProbablePrimeChainPart(uint32_t tid, __global primecoinInput_t *restrict input, __global sieveOutput_t *restrict output, 
		__global sieveTemp_t *restrict temp, __global uint32_t *restrict vPrimes, __global shaOutput_t *restrict sha ){
	WeavePart(tid,input,vPrimes,output,temp,sha);
}           

void MineProbablePrimeChainComplete(uint32_t tid, __global primecoinInput_t *restrict input, __global sieveOutput_t *restrict output, 
		__global sieveTemp_t *restrict temp, __global uint32_t *restrict vPrimes, __global shaOutput_t *restrict sha ){
	WeaveComplete(tid,input,vPrimes,output,temp,sha);
}           
