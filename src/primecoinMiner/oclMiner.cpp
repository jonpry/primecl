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
#include <stdio.h>
#include <string.h>

#include <CL/cl.h>
#include <vector>
#include <iostream>
#include <gmp.h>

#include "global.h"

#include "cldefs.h"

#include "clprime.h"
#include "ptx.h"

using namespace std;

static std::vector<unsigned int> vPrimes2;
static uint32_t vPrimesSize;
static struct primecoinBlock *tblock;

extern commandlineInput_t commandlineInput;

static void GeneratePrimeTable2(unsigned int nSieveSize)
{
   const unsigned int nPrimeTableLimit = nSieveSize ;
   vPrimes2.clear();
   // Generate prime table using sieve of Eratosthenes
   std::vector<bool> vfComposite (nPrimeTableLimit, false);
   for (unsigned int nFactor = 2; nFactor * nFactor < nPrimeTableLimit; nFactor++)
   {
      if (vfComposite[nFactor])
         continue;
      for (unsigned int nComposite = nFactor * nFactor; nComposite < nPrimeTableLimit; nComposite += nFactor)
         vfComposite[nComposite] = true;
   }
   for (unsigned int n = 2; n < nPrimeTableLimit; n++)
      if (!vfComposite[n])
         vPrimes2.push_back(n);
	   std::cout << "GeneratePrimeTable2() : prime table [1, " << nPrimeTableLimit << "] generated with " << vPrimes2.size() << " primes" << std::endl;
   vPrimesSize = vPrimes2.size(); 
	std::cout << "Last Prime: " << vPrimes2[vPrimesSize-1] << "\n"; 
}

static uint32_t PrimorialFast(uint32_t p, uint32_t *vPrimes, uint32_t vPrimesSize){
   unsigned int nPrimorial = 1;
   for(uint32_t i=0; i<vPrimesSize; i++) 
   {
      unsigned int nPrime = vPrimes[i];
      if (nPrime > p) break;
      nPrimorial *= nPrime;
   }
   return nPrimorial;
}

static void Primorial(uint32_t multiplier, uint32_t * vPrimes, uint32_t vPrimesSize, mpz_t primorial){
	mpz_set_ui(primorial,1);
	for(uint32_t i=0; i < vPrimesSize; i++){
		uint32_t nPrime = vPrimes[i];
		if(nPrime > multiplier)
			break;
		mpz_mul_ui(primorial,primorial,nPrime);
	} 
}

void copy_mpz(mpz_t n,mpzcls_t *ret){
	int i;
	int size = n->_mp_size;
	for(i=0; i < size; i++){
		ret->d[i] = n->_mp_d[i];
	}
	ret->size = size;
}

void copyBN(uint32_t idx, mpz_t n,mpzcl_t *src){
	uint64_t t[8];
	int i;
	int size = src->size[idx];
	if(size < 0 || size > 100){
		printf("BN too large! %s\n", __func__);
		exit(0);
	}


	for(i=0; i < size; i++){
		t[i] = D_REF(src->d,i,idx);
	}

	mpz_import(n, size, -1, 8, -1, 0, t);
}

double get_exec_time(cl_event* event){
	cl_ulong time_start, time_end;
	double total_time;

	clGetEventProfilingInfo(*event, CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL);
	clGetEventProfilingInfo(*event, CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL);
	total_time = time_end - time_start;

	return total_time / 1000000.0;
}

void print_exec_time(const char* name, cl_event *event){

	printf("\n%s: Execution time in milliseconds = %0.3f ms\n", name,  get_exec_time(event));
}

void print_bn(mpzcl_t* mpz, const char *name, uint32_t tid){
	uint32_t sz = mpz->size[tid];
	printf("%s: %d ", name, sz);
	int32_t i;
	for(i=sz-1; i>=0 && i < 50; i--){
		printf("%16.16lX", D_REF(mpz->d,i,tid));
	}
	printf("\n");
}

void initsieve(primecoinInput_t* input){
   	if (input->nOverrideTargetValue > 0)
      		input->lSieveTarget = input->nOverrideTargetValue;
   	else
      		input->lSieveTarget = TargetGetLength(input->blocks.nBits[0]);

  	// If Difficulty gets within 1/36th of next target length, its actually more efficent to
   	// increase the target length.. While technically worse for share val/hr - this should
   	// help block rate.
   	// Discussions with jh00 revealed this is non-linear, and graphs show that 0.1 diff is enough
   	// to warrant a switch
   
	if (GetChainDifficulty(input->blocks.nBits[0]) >= ((input->lSieveTarget + 1) - 0.1f))
      		input->lSieveTarget++;

  	if (input->nOverrideBTTargetValue > 0)
      		input->lSieveBTTarget = input->nOverrideBTTargetValue;
   	else
      		input->lSieveBTTarget = input->lSieveTarget; // Set to same as target
}


char *loadfile(const char* fname, size_t *sz_ptr){
	char *ret;
	FILE *fil=fopen(fname,"r");
	if(!fil)
		return 0;
	fseek(fil, 0L, SEEK_END);
	size_t sz = ftell(fil);

	fseek(fil, 0L, SEEK_SET);
	ret = new char[sz];
	fread(ret, sz, 1, fil);
	fclose(fil);
	
	*sz_ptr = sz;
	
	return ret;
}

void createFermat(cl_program prog, cl_kernel *k_fermat, cl_kernel *k_fermat_finish,
		cl_mem *input_mem, cl_mem *sieveoutput_mem, cl_mem *fermatoutput_mem,
		cl_mem *prime_mem, cl_mem *shaoutput_mem, cl_mem *fermattemp_mem,
		uint32_t *fermat_depth){
	int error=0;
	// get a handle and map parameters for the fermat kernel
	*k_fermat=clCreateKernel(prog, "fermat", &error);
	CheckErr("k_fermat");
	clSetKernelArg(*k_fermat, 0, sizeof(*input_mem), input_mem);
	clSetKernelArg(*k_fermat, 1, sizeof(*sieveoutput_mem),sieveoutput_mem);
	clSetKernelArg(*k_fermat, 2, sizeof(*fermatoutput_mem),fermatoutput_mem);
	clSetKernelArg(*k_fermat, 3, sizeof(*prime_mem),prime_mem);
	clSetKernelArg(*k_fermat, 4, sizeof(*shaoutput_mem), shaoutput_mem);
	clSetKernelArg(*k_fermat, 5, sizeof(*fermattemp_mem), fermattemp_mem);
	clSetKernelArg(*k_fermat, 6, sizeof(*fermat_depth), fermat_depth);

	// get a handle and map parameters for the fermat_finish kernel

	*k_fermat_finish=clCreateKernel(prog, "fermat_finish", &error);
	CheckErr("fermat_finish");
	clSetKernelArg(*k_fermat_finish, 0, sizeof(*input_mem), input_mem);
	clSetKernelArg(*k_fermat_finish, 1, sizeof(*sieveoutput_mem),sieveoutput_mem);
	clSetKernelArg(*k_fermat_finish, 2, sizeof(*fermatoutput_mem),fermatoutput_mem);
	clSetKernelArg(*k_fermat_finish, 3, sizeof(*prime_mem),prime_mem);
	clSetKernelArg(*k_fermat_finish, 4, sizeof(*shaoutput_mem), shaoutput_mem);
	clSetKernelArg(*k_fermat_finish, 5, sizeof(*fermattemp_mem), fermattemp_mem);
	clSetKernelArg(*k_fermat_finish, 6, sizeof(*fermat_depth), fermat_depth); 
}

void createSieve(cl_program prog, cl_kernel *k_sieve, const char *name, 
		cl_mem *input_mem, cl_mem *sievetemp_mem, uint32_t* length,
		cl_mem *sieveoutput_mem, cl_mem *prime_mem, cl_mem *shaoutput_mem){
	// get a handle and map parameters for the sieve kernel
	int error=0;
	*k_sieve=clCreateKernel(prog, name, &error);
	CheckErr("k_sieve");
	clSetKernelArg(*k_sieve, 0, sizeof(*input_mem), input_mem);
	clSetKernelArg(*k_sieve, 1, sizeof(*sievetemp_mem), sievetemp_mem);
	clSetKernelArg(*k_sieve, 2, sizeof(*length), length);
	clSetKernelArg(*k_sieve, 3, sizeof(*sieveoutput_mem),sieveoutput_mem);
	clSetKernelArg(*k_sieve, 4, sizeof(*prime_mem),prime_mem);
	clSetKernelArg(*k_sieve, 5, sizeof(*shaoutput_mem), shaoutput_mem);
}


cl_program load_program(const char *fname, cl_context context, cl_device_id device){
	size_t srcsize=0;
	cl_int error=0;
	const char *src=loadfile(fname,&srcsize);

	printf("%p\n",src);

	const char *srcptr[]={src};
	// Submit the source code of the rot13 kernel to OpenCL
	cl_program prog;
	if(src){
		prog = clCreateProgramWithSource(context,
		1, srcptr, &srcsize, &error);
	}else{
		unsigned char *pptx_code = ptx_code;
		size_t binsize = sizeof(ptx_code);
		int binstat=0;
		prog = clCreateProgramWithBinary(context,
		1, &device, &binsize, (const unsigned char**)&pptx_code, &binstat, &error);
	}
	CheckErr("create_program");

	// and compile it (after this we could extract the compiled version)
	error=clBuildProgram(prog, 0, NULL, "-I ./ -I ./src/primecoinMiner/ -cl-nv-verbose -cl-nv-maxrregcount=64", NULL, NULL);
	delete src;

	if (1){ //error < 0) {
    		// Determine the size of the log
    		size_t log_size;
    		clGetProgramBuildInfo(prog, device, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);

    		// Allocate memory for the log
    		char *log = (char *) malloc(log_size);

    		// Get the log
    		clGetProgramBuildInfo(prog, device, CL_PROGRAM_BUILD_LOG, log_size, log, NULL);

    		// Print the log
    		printf("%s\n", log);

			if(error<0)
				exit(0);
	}{
		size_t sizes[16];
		clGetProgramInfo(prog, CL_PROGRAM_BINARY_SIZES, 8, &sizes, 0);
		printf("%ld\n", sizes[0]);
		char **progs = new char *;
		progs[0] = new char[sizes[0]];
		clGetProgramInfo(prog, CL_PROGRAM_BINARIES, 8, progs, 0);
		
		if(commandlineInput.ptx){
			FILE *out = fopen("out.ptx","w");
			fwrite(progs[0],1,sizes[0],out);
			fclose(out);
		}
	//	printf("%s\n", progs[0]);
		//clLogPtx(prog, device, "oclConvolution.ptx");
	}

	CheckErr("build_program");

	return prog;
}

#define mpz_shift mpz_mul_2exp	

int isPrime(mpz_t m, int depth, int type, bool fraction){	
	mpz_t mshift,p,e,two,r;
	mpz_init(mshift);
	mpz_init(p);
	mpz_init(e);
	mpz_init(two);
	mpz_init(r);
	if(type == BITWIN_TYPE){
		mpz_shift(mshift,m,depth/2);
		if(depth % 2)
	    		mpz_add_ui(p,mshift,1);
		else
	    		mpz_sub_ui(p,mshift,1);
    	}else{
		mpz_shift(mshift,m,depth);
    		if(type == CC1_TYPE)
	    		mpz_sub_ui(p,mshift,1);
        	if(type == CC2_TYPE)
       	    		mpz_add_ui(p,mshift,1);
	}
	mpz_set_ui(two,2);
	mpz_sub_ui(e,p,1);
	mpz_powm(r,two,e,p);
	if(!fraction){
	  	int ret = mpz_cmp_ui(r,1) == 0;
		mpz_clear(mshift);
		mpz_clear(p);
		mpz_clear(e);
		mpz_clear(two);
		mpz_clear(r);
		return ret;
	}

   	// Failed Fermat test, calculate fractional length
   	mpz_sub(e, p, r);
   	mpz_mul_2exp(r, e, nFractionalBits);
   	mpz_tdiv_q(e, r, p);
   	unsigned int nFractionalLength = mpz_get_ui(e);
	if (nFractionalLength >= (1 << nFractionalBits))
      		printf("FermatProbablePrimalityTest() : fractional assert");


	mpz_clear(mshift);
	mpz_clear(p);
	mpz_clear(e);
	mpz_clear(two);
	mpz_clear(r);

   	return nFractionalLength & ((1<<24) - 1);
}

int getChainLength(mpz_t m, int knownLength, int type){
	while(isPrime(m,knownLength,type,0)){
		knownLength++;
	}
	return (knownLength << 24) | isPrime(m,knownLength,type,1);
}

bool checkBlock(void){
	return tblock->serverData.blockHeight != jhMiner_getCurrentWorkBlockHeight(0);
}

bool checkBlock(primecoinBlock_t *block){
	return block->serverData.blockHeight != jhMiner_getCurrentWorkBlockHeight(0);
}

void _fermatDone(void* d, void* target){
	//TODO: d can be used to hit clMiner class

	clMiner* miner = (clMiner*)d;
	miner->fermatDone((sieveOutput_t*)target);

}

void clMiner::fermatDone(sieveOutput_t* sout){
	uint32_t i,j;
	//printf("finish\n");

	// read some data back from fermat

	for(j=0; j < SM_COUNT; j++){
		for(i=0; i < sout[j].nBitwins && !checkBlock() && (sout[j].nBitwins < 1000U); i++){ //TODO: sometimes sieve crashes and spews huge number of chains
			if(sout[j].results.type[i]<0){
				//printf("Nothing at: %d Type: %d Hash: %d\n", i, sout->results.type[i], sout->results.hash[i]);
				continue;
			}
			//printf("Prime\n");
			primeStats.primeChainsFound++;

			clPrime *prime = new clPrime();
			prime->block = new primecoinBlock_t;
			*prime->block = *tblock;
			prime->sm = j;
			prime->tid = sout[j].results.hash[i];
			mpz_init(prime->mpzOrigin);

			uint32_t mult = sout[j].results.mults[i];
			prime->mult = mult;
			mpz_t m,hash;
			mpz_init(m);
			mpz_init(hash);
			mpz_mul_ui(m,fixedMultiplier,mult);

			if(sout[j].results.hash[i] >= STRIDE){
				printf("Opps: %d, %d\n", i, sout[j].results.hash[i]);
				exit(0);
			}

			copyBN(sout[j].results.hash[i],hash,&output[j].mpzHash);
			mpz_mul(prime->mpzOrigin,m,hash);
			prime->block->nonce = output[j].nonce[sout[j].results.hash[i]];
			prime->block->timestamp = linput[j].blocks.timestamp[0]; //nonce[sout[j].results.hash[i]];
			prime->type = sout[j].results.type[i];
			con_queue.push(prime);

			mpz_clear(m);
			mpz_clear(hash);


		}
	}

	primeStats.nSieveRounds+=STRIDE*(SIEVE_SIZE/(1024*1024))/ROUNDS*SM_COUNT;
}

void * chain_thread_func(void* d){
	clMiner* miner = (clMiner*)d;
	miner->chainFunc();
	return 0;
}

void clMiner::chainFunc(){
	while(1){
		clPrime *prime = (clPrime*)con_queue.pop();
		//printf("Work\n");

		if(!checkBlock(prime->block)){
			
	
			//Three is already known length
			int length = getChainLength(prime->mpzOrigin,3,prime->type);
			if(length < 0x05000000)
				continue;

			//printf("Prime: %d, len %X\n", i, length);
	
		      	if( length >= 0x06000000 )
				printf("Found 6-chain\n");

			int shareDifficultyMajor = 0;
		      	if( length >= 0x06000000 ){
			         shareDifficultyMajor = (int)(length>>24);
	      		}else{
	         		continue;
	    		}

			cout << "SM#:" << prime->sm << " Thread: " << prime->tid << " Origin: " << hex << prime->mpzOrigin << dec << "\n";

			if( length > primeStats.bestPrimeChainDifficulty )
	         		primeStats.bestPrimeChainDifficulty = length;

	         	// Update Stats
	         	primeStats.chainCounter[0][std::min(shareDifficultyMajor,12)]++;
	         	primeStats.chainCounter[prime->type+1][std::min(shareDifficultyMajor,12)]++;
	         	primeStats.nChainHit++;
	
	    		if(length < prime->block->serverData.nBitsForShare)
				continue;

			if(mpz_cmp(lastBlockOrigin,prime->mpzOrigin)==0){
				printf("Trying to submit duplicate block!");
				exit(0);
			}
			mpz_set(lastBlockOrigin,prime->mpzOrigin);

			mpz_mul_ui(prime->block->mpzPrimeChainMultiplier.get_mpz_t(),fixedMultiplier,prime->mult);
      
     			// generate block raw data;
			printf("Mult: %d, nonce: %x\n", prime->mult, prime->block->nonce);

	    		float shareDiff = GetChainDifficulty(length);
		        std::cout << "SHARE FOUND! - (SM#:" << prime->sm << ", Thread: " << prime->tid << ") - DIFF:" << shareDiff << " - TYPE:" << prime->type << std::endl;
		        jhMiner_pushShare_primecoin(0, prime->block);
	        	primeStats.foundShareCount ++;
		}
		delete prime->block;
		mpz_clear(prime->mpzOrigin);
		delete prime;
	}
}

void clMiner::oclInit(){
	uint32_t i;
	cl_int error;
	cl_platform_id platform;
	cl_device_id devicea[10];
	cl_uint platforms, devices;

  	pthread_t chain_thread;  
  	pthread_create(&chain_thread, NULL, chain_thread_func, this);


	mpz_init(lastBlockOrigin);
	mpz_init(lastOrigin);

	// Fetch the Platform and Device IDs; we only want one.
	error=clGetPlatformIDs(1, &platform, &platforms);
	CheckErr("platform_ids");
	error=clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 10, devicea, &devices);
	CheckErr("device_ids");
	printf("Found %d GPUs, selecting %d\n", devices, commandlineInput.gpu);
	device = devicea[commandlineInput.gpu];
	cl_context_properties properties[]={
		CL_CONTEXT_PLATFORM, (cl_context_properties)platform,
		0};
	// Note that nVidia's OpenCL requires the platform property
	context=clCreateContext(properties, 1, &device, NULL, NULL, &error);
	CheckErr("create_context");

	cq = clCreateCommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &error);
	CheckErr("create_cq");	

	prog = load_program("src/primecoinMiner/rot13.cl", context, device);
	GeneratePrimeTable2(160000);
	uint32_t prime_size = vPrimesSize;

	// Allocate memory for the kernel to work with
	sha_mem=clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(sha256cl_context)*SM_COUNT, NULL, &error);
	CheckErr("sha_mem");
	input_mem=clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(primecoinInput_t)*SM_COUNT, NULL, &error);
	CheckErr("input_mem");
	input_mem_host = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR, sizeof(primecoinInput_t)*SM_COUNT, NULL, &error);
	CheckErr("input_mem_host");

	shatemp_mem=clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(shaTemp_t)*SM_COUNT, NULL, &error);
	CheckErr("shatemp_mem");
	shaoutput_mem=clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(shaOutput_t)*SM_COUNT, NULL, &error);
	CheckErr("shaoutput_mem");

	//Sieve stuff
	sievetemp_mem=clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(sieveTemp_t)*SM_COUNT, NULL, &error);
	CheckErr("sievetemp_mem");
	sieveoutput_mem=clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(sieveOutput_t)*SM_COUNT, NULL, &error);
	CheckErr("sieveoutput_mem");
	prime_mem=clCreateBuffer(context,CL_MEM_READ_WRITE,prime_size*sizeof(uint32_t),NULL,&error);
	CheckErr("prime_mem");

	//Fermat stuff
	fermattemp_mem=clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(fermatTemp_t)*SM_COUNT, NULL, &error);
	CheckErr("fermattemp_mem");
	fermatoutput_mem=clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(sieveOutput_t)*SM_COUNT, NULL, &error);
	CheckErr("fermatoutput_mem");
	fermatoutputr_mem=clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(sieveOutput_t)*SM_COUNT, NULL, &error);
	CheckErr("fermatoutputr_mem");

	fermatoutputr_mem_host=clCreateBuffer(context, CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR, sizeof(sieveOutput_t)*SM_COUNT, NULL, &error);
	CheckErr("fermatoutputr_mem_host");

	printf("Data bytes : %d\n", sizeof(sha256cl_context) + sizeof(primecoinInput_t) + 
					sizeof(shaTemp_t) + sizeof(shaOutput_t) + 
					sizeof(sieveTemp_t) + sizeof(sieveOutput_t) + 
					sizeof(fermatTemp_t) + sizeof(sieveOutput_t));

	length = STRIDE;
	fermat_depth0 = 0;
	fermat_depth1 = 1;
	fermat_depth2 = 2;

	// get a handle and map parameters for the sha kernel
	k_sha=clCreateKernel(prog, "sha", &error);
	CheckErr("k_sha");
	clSetKernelArg(k_sha, 0, sizeof(sha_mem), &sha_mem);
	clSetKernelArg(k_sha, 1, sizeof(input_mem), &input_mem);
	clSetKernelArg(k_sha, 2, sizeof(shatemp_mem), &shatemp_mem);
	clSetKernelArg(k_sha, 3, sizeof(length), &length);
	clSetKernelArg(k_sha, 4, sizeof(shaoutput_mem), &shaoutput_mem);

	createSieve(prog,&k_sieve,"sieve",&input_mem,&sievetemp_mem,&length,
		&sieveoutput_mem,&prime_mem,&shaoutput_mem);
	createSieve(prog,&k_sieve_part,"sieve_part",&input_mem,&sievetemp_mem,&length,
		&sieveoutput_mem,&prime_mem,&shaoutput_mem);
	createSieve(prog,&k_sieve_complete,"sieve_complete",&input_mem,&sievetemp_mem,&length,
		&sieveoutput_mem,&prime_mem,&shaoutput_mem);

	createFermat(prog,&k_fermat,&k_fermat_finish,&input_mem,&sieveoutput_mem,&fermatoutput_mem,
			&prime_mem, &shaoutput_mem, &fermattemp_mem, &fermat_depth0);
	createFermat(prog,&k_fermat2,&k_fermat_finish2,&input_mem,&fermatoutput_mem,&sieveoutput_mem,
			&prime_mem, &shaoutput_mem, &fermattemp_mem, &fermat_depth1);
	createFermat(prog,&k_fermato,&k_fermat_finisho,&input_mem,&sieveoutput_mem,&fermatoutputr_mem,
			&prime_mem, &shaoutput_mem, &fermattemp_mem, &fermat_depth2);


	//Allocate host memories
	ftemp = (fermatTemp_t*)malloc(sizeof(fermatTemp_t)*SM_COUNT);
	stemp = (sieveTemp_t*)malloc(sizeof(sieveTemp_t)*SM_COUNT);
	sout = (sieveOutput_t*)malloc(sizeof(sieveOutput_t)*SM_COUNT);
	output = (shaOutput_t*)malloc(sizeof(shaOutput_t)*SM_COUNT);
	fermatresults = (sieveOutput_t*)clEnqueueMapBuffer(cq,fermatoutputr_mem_host, CL_TRUE, CL_MAP_WRITE, 0, sizeof(sieveOutput_t)*SM_COUNT, 0, NULL, NULL, NULL);

	//Assemble the stoppable opencl queue
	for(i=0; i < ROUNDS; i++){
		fermat_deptho[i] = (i << 16) + 2;
	}

	//Assemble the stoppable opencl queue
	mq.add(k_sha,0,new clTime("sha")); 
	mq.add(shaoutput_mem,0,sizeof(shaOutput_t)*SM_COUNT,output); 

	mq.add(k_sieve,1,new clTime("sieve init"));
	
	clTime* sieve_part_time = new clTime("sieve_part");
	clTime* sieve_complete_time = new clTime("sieve_complete");
	clTime* fermat_time = new clTime("fermat");
	clTime* fermat2_time = new clTime("fermat2");
	clTime* fermato_time = new clTime("fermato");

	for(i=0; i < ROUNDS; i++){
		mq.add(k_sieve_part,1,sieve_part_time);
		mq.add(k_sieve_complete,1,sieve_complete_time);
		mq.add(k_fermat,2,fermat_time);
		//mq.add(k_fermat_finish,2);
		mq.add(k_fermat2,2,fermat2_time);
		//mq.add(k_fermat_finish2,2);

		//mq.add(k_fermato,6,&fermat_deptho[i]);
		//mq.add(k_fermat_finisho,6,&fermat_deptho[i]);
		mq.add(k_fermato,2,fermato_time);
		//mq.add(k_fermat_finisho,2);
		mq.add(fermatoutputr_mem,0,sizeof(sieveOutput_t)*SM_COUNT, fermatresults, _fermatDone, (void*)this);
	}

	primes = new uint32_t[prime_size];
	for(i=0; i < prime_size; i++)
		primes[i] = vPrimes2[i];

	//Input data common to all kernels
	linput = (primecoinInput_t*)clEnqueueMapBuffer(cq,input_mem_host, CL_TRUE, CL_MAP_READ, 0, sizeof(primecoinInput_t)*SM_COUNT, 0, NULL, NULL, NULL);
	memset(linput,0,sizeof(primecoinInput_t)*SM_COUNT);

	for(i=0; i < SM_COUNT; i++){
		linput[i].nPrimorialHashFactor = 7;
		linput[i].nPrimorialMultiplier = 41;
		linput[i].nHashFactor = PrimorialFast(linput[i].nPrimorialHashFactor, primes, vPrimesSize);
		linput[i].nPrimes = vPrimesSize;
		linput[i].primeSeq = 4;
	}
	mpz_init(primorial);

	Primorial(linput->nPrimorialMultiplier,primes,vPrimesSize,primorial);

	for(i=0; i < SM_COUNT; i++)
		copy_mpz(primorial,&linput[i].mpzPrimorial);
	//TODO: It appears that mpzPrimorial < mpzMultiplierMin is never true 
	mpz_init(fixedMultiplier);
	mpz_div_ui(fixedMultiplier,primorial,linput->nHashFactor);
//	mpz_out_str(stdout,16,fixedMultiplier);

	for(i=0; i < SM_COUNT; i++)
		copy_mpz(fixedMultiplier,&linput[i].mpzFixedMultiplier);
}

void clMiner::mine(struct primecoinBlock *block) {
	tblock = block;
	int i,j,k;
	
	printf("Mine\n");

	size_t worksize = STRIDE;

	cl_int error;

	uint32_t prime_size = vPrimesSize;

	//Input data common to all kernels

#if 0
	//Hash testing codes
	for(i=0; i < 32; i++){
	//	block->prevBlockHash[i] = 0;
	//	block->merkleRoot[i] = 0;
	}
	//block->timestamp = 0;
	//block->version = 0;
	//block->nBits = 0x9FEF014;
#endif
	//this is where we pipein data from primecoin
	for(i=0; i < STRIDE; i++){
		for(j=0; j < SM_COUNT; j++){
			linput[j].blocks.nBits[i] = block->nBits;
			linput[j].blocks.timestamp[i] = block->timestamp+j; //Fudge timestamp to get different work per SM
			linput[j].blocks.version[i] = block->version;

			for(k=0; k < 8; k++){
				D_REF(linput[j].blocks.prevBlockHash,k,i) = *(uint32_t*)&block->prevBlockHash[k*4];
				D_REF(linput[j].blocks.merkleRoot,k,i) = *(uint32_t*)&block->merkleRoot[k*4];
			}
		}
	}

//	printf("Mine2\n");


	//Init sieve uses nBits	
	for(i=0; i < SM_COUNT; i++)
		initsieve(&linput[i]);

 
//	printf("Mine3\n");

	//Launch sha kernel
	// Send input data to OpenCL (async, don't alter the buffer!)
	error=clEnqueueWriteBuffer(cq, input_mem, CL_FALSE, 0, sizeof(primecoinInput_t)*SM_COUNT, linput, 0, NULL, NULL);
	CheckErr("wb input_mem");

	// Send input data to OpenCL (async, don't alter the buffer!)
	error=clEnqueueWriteBuffer(cq, prime_mem, CL_FALSE, 0, prime_size * sizeof(uint32_t), primes, 0, NULL, NULL);
	CheckErr("wb prime_mem");

//	printf("Mine4\n");

	// Perform the operation
	if(mq.run(cq,checkBlock))
		return;

	//Pretty much can just return
	return;

	error=clEnqueueReadBuffer(cq, shaoutput_mem, CL_FALSE, 0, sizeof(shaOutput_t)*SM_COUNT, output, 0, NULL, NULL);
	//printf("Err: %d\n", error);

	// read some data back from sieve
	error=clEnqueueReadBuffer(cq, sievetemp_mem, CL_FALSE, 0, sizeof(sieveTemp_t), stemp, 0, NULL, NULL);
	error=clEnqueueReadBuffer(cq, sieveoutput_mem, CL_FALSE, 0, sizeof(sieveOutput_t)*SM_COUNT, sout, 0, NULL, NULL);

	// read some data back from fermat
	error=clEnqueueReadBuffer(cq, fermattemp_mem, CL_FALSE, 0, sizeof(fermatTemp_t)*SM_COUNT, ftemp, 0, NULL, NULL);

	//printf("Err: %d\n", error);

	// Await completion of all the above
	error=clFinish(cq);

	//printf("Err: %d\n", error);
	if(checkBlock())
		return;
#if 0
	uint64_t v;
	for(i=0; i < 32; i++){
		printf("%d: ", i);
		for(j=0; j < 4; j++){
			v = D_REF(output.mpzHash.d,j,i);
			for(k=0; k < 8; k++){	
				printf("%2.2X", ((uint8_t*)&v)[k]);
			}
		}
		printf("\n");
	}

	//printf("Factor %u\n", output.factor);
	int hashes=0;
	for(i=0; i < worksize; i++){
		if(output.mod[i])
			hashes++;
	}

	printf("Generated: %d hashes\n", hashes);
	printf("Fixed factor: %d ", stemp->mpzFixedFactor.size[0]);
	for(i=0; i < stemp->mpzFixedFactor.size[0]; i++){
		printf("%16.16lx",stemp->mpzFixedFactor.d[(stemp->mpzFixedFactor.size[0]-i-1)*STRIDE]);
	}
	printf("\n");
#endif

#if 0
	for(i=0; i < STRIDE; i++){
		printf("%d\n", stemp->mpzFixedFactor.size[i]);
	}
#endif

#if 0
	printf("nFixedFactorCombinedMod: %d\n", stemp->nFixedFactorCombinedMod[0]);
	printf("nPrimeCombined: %d\n", stemp->nPrimeCombined[0]);
	printf("nFixedInverse: %X\n", stemp->nFixedInverse[0]);

	uint32_t accum1=0,accum2=0,accum3=0;
	for(i=linput->primeSeq; i < prime_size-1; i++){
		uint32_t multiplierPos = i;
		uint32_t j;
		for(j=0; j < 1; j++){
			accum1 += D_REF(stemp->vCunningham1Multipliers,multiplierPos,0);
			accum2 += D_REF(stemp->vCunningham2Multipliers,multiplierPos,0);
			accum3 += stemp->flatCC1Multipliers[multiplierPos];
		}
	}

	printf("Accum1: %X Accum2: %X Accum3: %X\n", accum1, accum2, accum3);

	printf("Elem10: %X Elem10: %X\n", stemp->flatCC1Multipliers[256], D_REF(stemp->vCunningham1Multipliers,256,0));


	accum1=0;
	for(i=0; i < SIEVE_WORDS; i++){
		for(j=0; j < 1; j++){
			accum1 += stemp->vfExtBitwin[i];
		}
	}

	printf("Accum1: %X\n", accum1);

	for(i=0; i < 16; i++){
		uint32_t j;
		printf("VCC1:%4.4X ",i*16);
		for(j=0; j < 16; j++)
			printf("%8.8X ", stemp->vfExtCC1[i*16+j]);
		printf("\n");
	}

	bool pass=1;
	for(i=0; i < 512; i++){
		uint32_t j;
		uint32_t err_cnt=0;
		for(j=0; j < MULT_SIZE; j++){
			if(D_REF(stemp->vCunningham1Multipliers,j,i) != stemp->flatCC1Multipliers[i*MULT_WORDS+j]){
				if(err_cnt++==0)
					pass=0;
			}
		}
		if(err_cnt){
	//		printf("Mismatch at %d %d\n", i, err_cnt);
		}
	}

	printf("Flatten test: %d\n", pass);

	printf("Bitwins: %d, CC1: %d, CC2: %d\n", sout->nBitwins, sout->nCC1s, sout->nCC2s);

	printf("Bitwins: %d, CC1: %d, CC2: %d\n", sout->nBitwins-sout->nCC2s, sout->nCC1s, sout->nCC2s-sout->nCC1s);

	int32_t hashidx[512];
	for(i=0; i < STRIDE; i++){
//  		printf("Prime: %d\n", sout->prime[i]);
//		hashidx[ftemp->i[i]] = i;
	}

	for(i=0; i < STRIDE; i++){
	//	print_bn((mpzcl_t*)&ftemp->mpzOriginMinusOne,"O", hashidx[i]);
	//	print_bn((mpzcl_t*)&ftemp->mpzR,"R", hashidx[i]);
	}
#endif
	printf("WTF!!!!!!!!!!!!!!!!!\n");
	uint32_t fprimes=0;
	uint32_t candidates=0;
	for(j=0; j < SM_COUNT; j++){
		for(i=0; i < sout[j].nBitwins && !checkBlock() && (sout[j].nBitwins < 1000U); i++){ //TODO: sometimes sieve crashes and spews huge number of chains
			if(sout[j].results.type[i]<0){
				//printf("Nothing at: %d Type: %d Hash: %d\n", i, sout->results.type[i], sout->results.hash[i]);
				continue;
			}
			primeStats.primeChainsFound++;

			uint32_t mult = sout[j].results.mults[i];
			mpz_t m,hash;
			mpz_init(m);
			mpz_init(hash);
			mpz_mul_ui(m,fixedMultiplier,mult);

			if(sout[j].results.hash[i] >= STRIDE){
				printf("Opps: %d, %d\n", i, sout[j].results.hash[i]);
				exit(0);
			}

			copyBN(sout[j].results.hash[i],hash,&output[j].mpzHash);
			mpz_mul(m,m,hash);

			//Four is already known length
			int length = getChainLength(m,4,sout[j].results.type[i]);
			if(length < 0x05000000)
				continue;

			//printf("Prime: %d, len %X\n", i, length);
	
	      		if( length >= 0x06000000 )
				printf("Found 6-chain\n");

			int shareDifficultyMajor = 0;
	      		if( length >= 0x06000000 )
      			{
			         shareDifficultyMajor = (int)(length>>24);
      			}else{
        	 		continue;
	      		}

			cout << "Origin: " << hex << m << dec << "\n";

			if( length > primeStats.bestPrimeChainDifficulty )
	         		primeStats.bestPrimeChainDifficulty = length;
	
	    		if(length < block->serverData.nBitsForShare)
				continue;

			if(mpz_cmp(lastBlockOrigin,m)==0){
				printf("Trying to submit duplicate block!");
				exit(0);
			}
			mpz_set(lastBlockOrigin,m);

	         	// Update Stats
	         	primeStats.chainCounter[0][std::min(shareDifficultyMajor,12)]++;
	         	primeStats.chainCounter[sout[j].results.type[i]+1][std::min(shareDifficultyMajor,12)]++;
	         	primeStats.nChainHit++;

			mpz_mul_ui(block->mpzPrimeChainMultiplier.get_mpz_t(),fixedMultiplier,mult);
      
     			// generate block raw data
			block->nonce = output[j].nonce[sout[j].results.hash[i]];
			block->timestamp = linput[j].blocks.timestamp[0]; //nonce[sout[j].results.hash[i]];
			printf("Hash: %d, nonce: %x\n", sout[j].results.hash[i], block->nonce);

///////////////

	     		if(0){
				primecoinBlock_generateHeaderHash(block, block->blockHeaderHash.begin());
				uint256 phash = block->blockHeaderHash;
		        	mpz_class mpzHash;
		        	mpz_set_uint256(mpzHash.get_mpz_t(), phash);
	
				cout << hex << mpzHash << "\n";
				cout << hex << hash << "\n";
//				exit(0);
			}
///////////////

	    		float shareDiff = GetChainDifficulty(length);
		        std::cout << "SHARE FOUND! - (SM#:" << j << ") - DIFF:" << shareDiff << " - TYPE:" << sout[j].results.type[i] << std::endl;
		        jhMiner_pushShare_primecoin(0, block);
	        	primeStats.foundShareCount ++;
		}
	}
//   	primeStats.nWaveTime += 14000;
//   	primeStats.nWaveRound +=512;
   	primeStats.nSieveRounds+=STRIDE*SM_COUNT;
   	primeStats.nCandidateCount += 150000*SM_COUNT;

//	printf("Found %d primes out of %d\n", fprimes, candidates);

#if 0
	hashidx+=11; //The crashing hash

//	exit(0);


	print_bn((mpzcl_t*)&ftemp->mpzE,"E", hashidx);
	print_bn((mpzcl_t*)&ftemp->mpzB,"B", hashidx);
	print_bn(&ftemp->mpzChainOrigin,"Origin", hashidx);
	print_bn(&ftemp->mpzOriginPlusOne,"OriginP1", hashidx);
	print_bn(&ftemp->mpzOriginMinusOne,"OriginM1", hashidx);
	print_bn(&ftemp->mpzNewtonDen,"NewtonDen", hashidx);
	print_bn(&ftemp->mpzNewtonX1s,"NewtonX1s", hashidx);
	print_bn(&ftemp->mpzNewtonX2,"NewtonX2", hashidx);
	print_bn(&ftemp->mpzNewtonX2s,"NewtonX2s", hashidx);
	print_bn(&ftemp->mpzNewtonX,"NewtonX", hashidx);
	print_bn((mpzcl_t*)&ftemp->mpzNewtonDenP,"NewtonDenP", hashidx);
	print_bn(&ftemp->mpzNewtonDenPS,"NewtonDenPS", hashidx);
	print_bn(&ftemp->mpzNewton2s,"Newton2s", hashidx);
	print_bn(&ftemp->mpzNewtonDiff,"NewtonDiff", hashidx);
	print_bn((mpzcl_t*)&ftemp->mpzNewtonProd,"NewtonProd", hashidx);
	print_bn((mpzcl_t*)&ftemp->mpzInv,"Inv", hashidx);
	print_bn((mpzcl_t*)&ftemp->mpzBarrettT1,"BarrettT1", hashidx);
	print_bn((mpzcl_t*)&ftemp->mpzBarrettT2,"BarrettT2", hashidx);
	print_bn((mpzcl_t*)&ftemp->mpzBarrettT3,"BarrettT3", hashidx);
	print_bn((mpzcl_t*)&ftemp->mpzBarrettT4,"BarrettT4", hashidx);
	print_bn((mpzcl_t*)&ftemp->mpzBarrettN,"BarrettN", hashidx);
	print_bn((mpzcl_t*)&ftemp->mpzBarrettM,"BarrettM", hashidx);
	print_bn((mpzcl_t*)&ftemp->mpzBarrettA,"BarrettA", hashidx);
	print_bn((mpzcl_t*)&ftemp->mpzPowTemp,"PowTemp", hashidx);
	printf("E: %d\n", ftemp->e[hashidx]);
	printf("P: %d\n", ftemp->p[hashidx]);
	printf("al: %8.8X\n", ftemp->al[hashidx]);
	printf("bl: %8.8X\n", ftemp->bl[hashidx]);
	printf("il: %8.8X\n", ftemp->il[hashidx]);

	for(i=0; i < 50; i++){
		if((i+4)%5==0)
			printf("\n");
		print_bn((mpzcl_t*)&ftemp->mpzR[i],"R", hashidx);
	}
#endif
	//printf("I: %d, Hash: %d\n", ftemp->i, ftemp->hash); 

	//printf("primorial: %lx %d\n", D_REF(output.mpzPrimorial.d,0,0), output.mpzPrimorial.size[0]); 
	//printf("multiplier min: %x %d\n", D_REF(output.mpzMultiplierMin.d,0,0), output.mpzMultiplierMin.size[0]); 

	//printf("difficulty: %d\n", output.difficulty[0]); 

	//print_exec_time("SHA",&event_sha);
	//print_exec_time("Sieve",&event_sieve);
	//print_exec_time("Fermat",&event_fermat);
}
