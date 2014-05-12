/********************************************************/
/* NV OpenCL XPM Miner */
/* Copyright 2014 Jon Pry */
/********************************************************/
 

#pragma OPENCL EXTENSION cl_khr_byte_addressable_store : enable
                                                                                            
               
typedef unsigned long uint64_t;
typedef unsigned uint32_t;
typedef int int32_t;
typedef unsigned char uint8_t;
   
#include "cldefs.h"
 
#include "sha256.cl" 
#include "bn.cl"
#include "bn_div.cl" 
#include "prime.cl"   
#include "fermat.cl"
#include "fermatunrolled.cl"  
                                  
                                             
void memcpy_progblock(uint32_t tid, __global uint32_t *restrict dst, __global primecoinBlockcl_t *restrict src){
	int i;  
#if 1
	D_REF(dst,0,tid) = src->version[tid];  
	for(i=0; i < 8; i++){ 
		D_REF(dst,i+1,tid) = D_REF(src->prevBlockHash,i,tid); 
		D_REF(dst,i+9,tid) = D_REF(src->merkleRoot,i,tid);
	}  
   
	D_REF(dst,17,tid) = src->timestamp[tid]; 
	D_REF(dst,18,tid) = src->nBits[tid]; 
	D_REF(dst,19,tid) = src->nonce[tid];         
#else 
	for(i=0; i < 20; i++){
		D_REF(dst,i,tid) = 0;   
	}
//	D_REF(dst,17,tid) = src->timestamp[tid];  
//	D_REF(dst,18,tid) = D_REF(src->nBits,i,tid); 
//	D_REF(dst,19,tid) = D_REF(src->nonce,i,tid);      

#endif            
}                                         		   
                                                          
void generateHeaderHash(uint32_t tid, __global    sha256cl_context *restrict ctx,
			__global primecoinBlockcl_t *restrict blocks, __global uint32_t *restrict temp1,
			 __global uint32_t *restrict hashOutput)
{
    	sha256_starts(tid,ctx);   

    	memcpy_progblock(tid,temp1,blocks);

    	sha256_update(tid,ctx, temp1, 80);   

    	sha256_finish(tid,ctx, hashOutput); 	
    	sha256_starts(tid,ctx); // is this line needed?

    	sha256_update(tid,ctx, hashOutput, 32); 
    	sha256_finish(tid,ctx, hashOutput);	  	 
}               
               
void hashToBN(uint32_t tid, __global mpzcl_t *restrict mpz, __global uint32_t *restrict hash){
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
                  
void initPrimeMin(uint32_t tid, __global mpzcl_t *restrict mpz){
	mpz->size[tid] = 4;
	uint32_t i;
	for(i=0; i < 3; i++){
		D_REF(mpz->d,i,tid) = 0;
	}
	D_REF(mpz->d,3,tid) = 1UL << 63; 
}                           
                             
__kernel void sha				 	
    (   __global    sha256cl_context *restrict ctx
    ,   __global primecoinInput_t *restrict input  
    ,   __global shaTemp_t *restrict temp
    ,   unsigned length		
    ,   __global shaOutput_t *restrict output
    )	 		 	 			
{	 	  	       				 
    	const uint32_t tid = get_local_id(0);	 
	const uint32_t gid = get_group_id(0);
	ctx += gid;
	input += gid;
	temp += gid;
	output += gid;
 
#if 0	  
    	if(tid > length)   
     	 	return;	  
#endif 
//	return;
   	sha256_init(tid,ctx);      
//	return;	  
	initPrimeMin(tid,&temp->mpzPrimeMin);	
//	return;
	input->blocks.nonce[tid] = 0x00400000 * tid;
  	uint32_t maxNonce = (0x00400000 * tid) | (0x00FFFFFF >> 2);  
  
	// Primecoin: try to find hash divisible by primorial
 #if 1
	uint32_t isdiv=0;
	while ((D_REF(temp->hash,7,tid) < 0x80000000U || !isdiv) && input->blocks.nonce[tid] < maxNonce){
		input->blocks.nonce[tid]++;
		generateHeaderHash(tid,ctx,&input->blocks,temp->hashTemp,temp->hash);
		hashToBN(tid,&output->mpzHash,temp->hash);
		isdiv = BN256DivisibleBy(tid,&output->mpzHash,input->nHashFactor);
		output->mod[tid] = isdiv;  
	}          
 #endif  
	output->nonce[tid] = input->blocks.nonce[tid];     
        
	//Some hashes are probably not good. But in statistical sense it is fine to just continue. Network will just reject the bad ones. 
}                 
                                                                   
                                      
__kernel void sieve				 	
    (   __global primecoinInput_t *restrict input  
    ,   __global sieveTemp_t *restrict temp
    ,   unsigned length		
    ,   __global sieveOutput_t *restrict output
    ,   __global uint32_t *restrict vPrimes
    ,	__global shaOutput_t *restrict sha
    )			 	 			
{	    	   		 			  
    	const uint32_t tid = get_local_id(0);	  
	const uint32_t gid = get_group_id(0);
	input += gid;
	temp += gid;
	output += gid;                 
	sha += gid; 
                                 
	MineProbablePrimeChain(tid,input,output,temp,vPrimes,sha);	 	 	
}    
     
__kernel void sieve_part				 	  
    (   __global primecoinInput_t *restrict input  
    ,   __global sieveTemp_t *restrict temp
    ,   unsigned length		
    ,   __global sieveOutput_t *restrict output
    ,   __global uint32_t *restrict vPrimes
    ,	__global shaOutput_t *restrict sha
    )			 	 			
{	    	   		 			  
    	const uint32_t tid = get_local_id(0);	                   
	const uint32_t gid = get_group_id(0);
	input += gid;
	temp += gid;
	output += gid;                 
	sha += gid;
                                                          
	MineProbablePrimeChainPart(tid,input,output,temp,vPrimes,sha);	 	 	
}       
      
__kernel void sieve_complete		 	
    (   __global primecoinInput_t *restrict input  
    ,   __global sieveTemp_t *restrict temp
    ,   unsigned length		
    ,   __global sieveOutput_t *restrict output
    ,   __global uint32_t *restrict vPrimes
    ,	__global shaOutput_t *restrict sha
    )			 	 			
{	    	   		 			  
    	const uint32_t tid = get_local_id(0);	                   
	const uint32_t gid = get_group_id(0);
	input += gid;
	temp += gid;
	output += gid;                 
	sha += gid;
                                                          
	MineProbablePrimeChainComplete(tid,input,output,temp,vPrimes,sha);	 	 	
}    
	                 
uint32_t compact(uint32_t tid, __global bresult_t *restrict output, __global bresult_t *restrict input, 
				__local volatile uint32_t *sem, uint32_t begin, uint32_t end){
#if 1
	uint32_t i;
	for(i=begin+tid; i < end; i+=STRIDE){
		if(!input->prime[i])
			continue;
		uint32_t loc = atomic_add(sem,1); 
		output->prime[loc] = input->prime[i];		
		output->mults[loc] = input->mults[i];		
		output->hash[loc] = input->hash[i];
		output->type[loc] = input->type[i];				
	}

	barrier( CLK_LOCAL_MEM_FENCE );
	uint32_t val = *sem;
	//Pad the array with some -1 types to warp alignment
	if((val % 32) && (tid < 32 - (val % 32)))
		output->type[val + tid] = -1;
	if(tid==0){
		val = (val % 32)?val + 32 - (val % 32) : val; 
		*sem = val;
	}
	barrier( CLK_LOCAL_MEM_FENCE );
	return val;
#else
	uint32_t i;
	for(i=begin+tid; i < end; i+= STRIDE){
		output->prime[i] = input->prime[i];		
		output->mults[i] = input->mults[i];		
		output->hash[i] = input->hash[i];
		output->type[i] = input->type[i];		
	}
	*sem = end;
	return end;
#endif
} 	            					
                    
__kernel void fermat				 	
    (   __global primecoinInput_t *restrict input  	
    ,   __global sieveOutput_t *restrict sieve
    ,   __global sieveOutput_t *restrict output
    ,   __global uint32_t *restrict vPrimes
    ,	__global shaOutput_t *restrict sha
    ,   __global fermatTemp_t *restrict temp
    ,   uint32_t length
    )	{
	const uint32_t tid = get_local_id(0);	
	const uint32_t gid = get_group_id(0);
	input += gid;
	sieve += gid;
	output += gid;                 
	sha += gid;
	temp += gid;

	if((length>>16) != 0){
	//	output += (length>>16)*SM_COUNT;
		length&=0xFF;
	}
                         
#if 1   
	//Find chain origin
	mulBNBNS(tid,&temp->mpzFixedFactor,&sha->mpzHash,&input->mpzFixedMultiplier);
	barrier( CLK_GLOBAL_MEM_FENCE );

	uint32_t i,nCandidates = min(sieve->nBitwins,sieve->nTested+1024*32); 
	volatile local uint32_t sem;
	volatile local uint32_t idx[STRIDE/32]; 
	if(tid==0)
		sem=sieve->nTested;
	barrier(CLK_LOCAL_MEM_FENCE);

	//Mine the primes  
	while(1){
		if(tid%32==0){
		   	i = atomic_add(&sem,32);
			idx[tid/32] = i;
		}
		mem_fence(CLK_LOCAL_MEM_FENCE);

		i = idx[tid/32] + tid%32;
 		if(i >= nCandidates)
			break;
		int32_t type = sieve->results.type[i];
		if(type >= 0)
		{
			//TODO: make sure there is a way to correlate results with specific hashes and multipliers
			temp->i[tid] = i;
		
			mulBN2(tid,sieve->results.hash[i],&temp->mpzChainOrigin,&temp->mpzFixedFactor,sieve->results.mults[i]);
			
			sieve->results.prime[i] = ProbablePrimeChainTestFast(tid,temp,sieve->results.type[i],length);
		}else{
			sieve->results.prime[i] = 0;
		}		
	}

	sieve->nTested = nCandidates; 

	barrier( CLK_GLOBAL_MEM_FENCE );

	volatile local uint32_t sem2;
	if(tid==0)
		sem2=0;
	barrier(CLK_GLOBAL_MEM_FENCE);

	uint32_t cc1s,cc2s,bitwins;
	//Compact the list
	cc1s = compact(tid,&output->results,&sieve->results,&sem2,0,sieve->nCC1s);

	cc2s = compact(tid,&output->results,&sieve->results,&sem2,sieve->nCC1s,sieve->nCC2s);
	
	bitwins = compact(tid,&output->results,&sieve->results,&sem2,sieve->nCC2s,sieve->nBitwins);
	
	if(tid==0){
		output->nCC1s = cc1s;
		output->nCC2s = cc2s;
		output->nBitwins = bitwins;
		output->nLayer = sieve->nLayer+1;
		output->nTested=0;
	}	
}
     
          
__kernel void fermat_finish				 	
    (   __global primecoinInput_t *restrict input  	
    ,   __global sieveOutput_t *restrict sieve
    ,   __global sieveOutput_t *restrict output
    ,   __global uint32_t *restrict vPrimes
    ,	__global shaOutput_t *restrict sha
    ,   __global fermatTemp_t *restrict temp
    ,   uint32_t length
    )	{
	const uint32_t tid = get_local_id(0);
	const uint32_t gid = get_group_id(0);
	input += gid;
	sieve += gid;
	output += gid;                 
	sha += gid;
	temp += gid;

	if((length>>16) != 0){
	//	output += (length>>16)*SM_COUNT;
		length&=0xFF;
	}
 
	volatile local uint32_t sem2;
	if(tid==0)
		sem2=0;
	barrier(CLK_GLOBAL_MEM_FENCE);

	uint32_t cc1s,cc2s,bitwins;
	//Compact the list
	cc1s = compact(tid,&output->results,&sieve->results,&sem2,0,sieve->nCC1s);

	cc2s = compact(tid,&output->results,&sieve->results,&sem2,sieve->nCC1s,sieve->nCC2s);
	
	bitwins = compact(tid,&output->results,&sieve->results,&sem2,sieve->nCC2s,sieve->nBitwins);
	
	if(tid==0){
		output->nCC1s = cc1s;
		output->nCC2s = cc2s;
		output->nBitwins = bitwins;
		output->nLayer = sieve->nLayer+1;
		output->nTested=0;
	}	
#endif   
}                           
                                                                         
                             
