#pragma OPENCL EXTENSION cl_khr_byte_addressable_store : enable


typedef unsigned long uint64_t;
typedef unsigned uint32_t;
typedef unsigned char uint8_t;

#include "cldefs.h"
 
#include "sha256.cl" 
#include "bn.cl"
//#include "bn_div.cl" 
//#include "prime.cl"    
    

void memcpy_progblock(uint32_t tid, __global uint32_t dst[], __global primecoinBlock_t *src){
	int i; 

	D_REF(dst,0,tid) = D_REF(src->version,i,tid); 
	for(i=0; i < 8; i++){
		D_REF(dst,i+1,tid) = D_REF(src->prevBlockHash,i,tid); 
		D_REF(dst,i+9,tid) = D_REF(src->merkleRoot,i,tid);
	}

	D_REF(dst,17,tid) = D_REF(src->timestamp,i,tid); 
	D_REF(dst,18,tid) = D_REF(src->nBits,i,tid); 
	D_REF(dst,19,tid) = D_REF(src->nonce,i,tid); 
}

// Compute Primorial number p#
// Fast 32-bit version assuming that p <= 23

uint32_t PrimorialFast(uint32_t p, __global uint32_t *vPrimes, uint32_t vPrimesSize){
   unsigned int nPrimorial = 1;
   for(uint32_t i=0; i<vPrimesSize; i++) 
   {
      unsigned int nPrime = vPrimes[i];
      if (nPrime > p) break;
      nPrimorial *= nPrime;
   }
   return nPrimorial;
}

void Primorial(uint32_t tid, uint32_t multiplier, __global uint32_t *vPrimes, uint32_t vPrimesSize, __global mpz_t *mpzPrimorial){
	setBN(tid,mpzPrimorial,1);
#if 1
	for(uint32_t i=0; i < vPrimesSize; i++){
		uint32_t nPrime = vPrimes[i];
		if(nPrime > multiplier)
			break;
		mulBN(tid,mpzPrimorial,mpzPrimorial,nPrime);		
	} 
#endif
}


void generateHeaderHash(uint32_t tid, __global    sha256_context*    ctx,
			__global primecoinBlock_t* blocks, __global uint32_t temp1[20*STRIDE],
			 __global uint32_t hashOutput[8*STRIDE])
{
    	sha256_starts(tid,ctx);  

//    	memcpy_progblock(tid,temp1,blocks);
//    	sha256_update(tid,ctx, temp1, 80);

    	sha256_finish(tid,ctx, hashOutput);	
//    	sha256_starts(tid,ctx); // is this line needed?

//    	sha256_update(tid,ctx, hashOutput, 32);
//    	sha256_finish(tid,ctx, hashOutput);		 
} 

void hashToBN(uint32_t tid, __global mpz_t *mpz, __global uint32_t hash[8*STRIDE]){
	int i,j;
	for(i=0; i < 4; i++){
		uint64_t v = 0;
		for(j=0; j < 2; j++){
			v |= (uint64_t)D_REF(hash,i*2 + j,tid) << (j*32);
		}
		D_REF(mpz->d,i,tid) = v;
	}
	mpz->size[tid] = 4;
} 

void initPrimeMin(uint32_t tid, __global mpz_t *mpz){
	mpz->size[tid] = 4;
	uint32_t i;
	for(i=0; i < 3; i++){
		D_REF(mpz->d,i,tid) = 0;
	}
	D_REF(mpz->d,3,tid) = 1UL << 63;
}

__kernel void rot13				 	
    (   __global    sha256_context*    ctx
    ,   __global primecoinInput_t *input  
    ,   __global primecoinTemp_t *temp
    ,   __global uint32_t *vPrimes
    ,   uint32_t vPrimesSize	 	
    ,   unsigned length		
    ,   __global primecoinOutput_t *output
    )			 	 			
{							
    	const uint tid = get_global_id(0);	 

#if 0	
    	if(tid > length) 
     	 	return;	  
#endif
   	sha256_init(tid,ctx);  
	barrier(0);
	initPrimeMin(tid,&temp->mpzPrimeMin);	

	//TODO: do this on cpu 
	uint32_t nHashFactor = PrimorialFast(input->nPrimorialHashFactor, vPrimes, vPrimesSize);	
	//just for fun  
	output->factor = nHashFactor;    
	
	input->blocks.nonce[tid] = 0x00400000 * tid;
  	uint32_t maxNonce = (0x00400000 * tid) | (0x00FFFFFF >> 2); 
	  
 
    	//TODO: temporary fudge the block hash. in reality whoever starts up gpu miner will already have block setup
    	input->blocks.timestamp[tid] = tid;

	// Primecoin: try to find hash divisible by primorial
	generateHeaderHash(tid,ctx,&input->blocks,temp->hashTemp,output->hash);
	hashToBN(tid,&temp->mpzHash,output->hash);

	uint32_t isdiv = BN256DivisibleBy(tid,&temp->mpzHash,nHashFactor);
 #if 0
	while ((D_REF(output->hash,7,tid) < 0x80000000U || !isdiv) && input->blocks.nonce[tid] < maxNonce){ 
		input->blocks.nonce[tid]++;
		generateHeaderHash(tid,ctx,&input->blocks,temp->hashTemp,output->hash);
		hashToBN(tid,&temp->mpzHash,output->hash);
		isdiv = BN256DivisibleBy(tid,&temp->mpzHash,nHashFactor);
		output->mod[tid] = isdiv;  
	}       
 #endif     
	//Some hashes are probably not good. But in statistical sense it is fine to just continue. Network will just reject the bad ones. 
	Primorial(tid,input->nPrimorialMultiplier, vPrimes, vPrimesSize, &output->mpzPrimorial);
	  
	//mpz_class mpzMultiplierMin = mpzPrimeMin * nHashFactor / mpzHash + 1;   
	mulBN(tid,&output->mpzMultiplierMin,&temp->mpzPrimeMin,nHashFactor);       
//	divBN(tid,&output->mpzMultiplierMin,&output->mpzMultiplierMin,&temp->mpzHash);  
 	//TODO: fix the division so this can work  
#if 0 
	while (mpzPrimorial < mpzMultiplierMin )
	{      
		if (!PrimeTableGetNextPrime(nPrimorialMultiplier))
			error("PrimecoinMiner() : primorial minimum overflow");
		Primorial(nPrimorialMultiplier, mpzPrimorial); 
	}
#endif	    
        
#if 0
      	mpz_class mpzFixedMultiplier;
      	if (mpzPrimorial > nHashFactor) {
            	mpzFixedMultiplier = mpzPrimorial / nHashFactor;
        } else {
            	mpzFixedMultiplier = 1;
        }
#endif
 
	// Primecoin: mine for prime chain
	//pretty much go time
	//MineProbablePrimeChain(psieve, primecoinBlock, mpzFixedMultiplier, fNewBlock, nTriedMultiplier, nProbableChainLength, nTests, nPrimesHit, threadIndex, mpzHash, nPrimorialMultiplier);
//	MineProbablePrimeChain(tid,input,output,temp);			
}		         					


