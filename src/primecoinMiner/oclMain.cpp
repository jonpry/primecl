/* Simple Hello World for OpenCL, written in C.
 * For real code, check for errors. The error code is stored in all calls here,
 * but no checking is done, which is certainly bad. It'll likely simply crash
 * right after a failing call.
 *
 * On GNU/Linux with nVidia OpenCL, program builds with -lOpenCL.
 * Not sure about other platforms.
 */

#include <stdio.h>
#include <string.h>

#include <CL/cl.h>
#include <vector>
#include <iostream>
#include <gmp.h>
#include "cldefs.h"

#include "clprime.h"

using namespace std;

std::vector<unsigned int> vPrimes;
uint32_t vPrimesSize;

void GeneratePrimeTable(unsigned int nSieveSize)
{
   const unsigned int nPrimeTableLimit = nSieveSize ;
   vPrimes.clear();
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
         vPrimes.push_back(n);
	   std::cout << "GeneratePrimeTable() : prime table [1, " << nPrimeTableLimit << "] generated with " << vPrimes.size() << " primes" << std::endl;
   vPrimesSize = vPrimes.size(); 
	std::cout << "Last Prime: " << vPrimes[vPrimesSize-1] << "\n"; 
}

uint32_t PrimorialFast(uint32_t p, uint32_t *vPrimes, uint32_t vPrimesSize){
   unsigned int nPrimorial = 1;
   for(uint32_t i=0; i<vPrimesSize; i++) 
   {
      unsigned int nPrime = vPrimes[i];
      if (nPrime > p) break;
      nPrimorial *= nPrime;
   }
   return nPrimorial;
}

void Primorial(uint32_t multiplier, uint32_t * vPrimes, uint32_t vPrimesSize, mpz_t primorial){
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
	for(i=0; i < size; i++){
		t[i] = D_REF(src->d,i,idx);
	}

	mpz_import(n, size, -1, 8, -1, 0, t);
}

uint32_t TargetGetLength(unsigned int nBits)
{
   	return ((nBits & TARGET_LENGTH_MASK) >> nFractionalBits);
}

float GetChainDifficulty(unsigned int nChainLength)
{
	return (float)nChainLength / 16777216.0f;
}

void print_exec_time(const char* name, cl_event *event){
	cl_ulong time_start, time_end;
	double total_time;

	clGetEventProfilingInfo(*event, CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL);
	clGetEventProfilingInfo(*event, CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL);
	total_time = time_end - time_start;
	printf("\n%s: Execution time in milliseconds = %0.3f ms\n", name, (total_time / 1000000.0) );
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
	// get a handle and map parameters for the fermat kernel
	*k_fermat=clCreateKernel(prog, "fermat", 0);
	clSetKernelArg(*k_fermat, 0, sizeof(*input_mem), input_mem);
	clSetKernelArg(*k_fermat, 1, sizeof(*sieveoutput_mem),sieveoutput_mem);
	clSetKernelArg(*k_fermat, 2, sizeof(*fermatoutput_mem),fermatoutput_mem);
	clSetKernelArg(*k_fermat, 3, sizeof(*prime_mem),prime_mem);
	clSetKernelArg(*k_fermat, 4, sizeof(*shaoutput_mem), shaoutput_mem);
	clSetKernelArg(*k_fermat, 5, sizeof(*fermattemp_mem), fermattemp_mem);
	clSetKernelArg(*k_fermat, 6, sizeof(*fermat_depth), fermat_depth);

	// get a handle and map parameters for the fermat_finish kernel
	*k_fermat_finish=clCreateKernel(prog, "fermat_finish", 0);
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
	*k_sieve=clCreateKernel(prog, name, 0);
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
	cl_program prog=clCreateProgramWithSource(context,
		1, srcptr, &srcsize, &error);
	printf("Err: %d\n", error);

	// and compile it (after this we could extract the compiled version)
	error=clBuildProgram(prog, 0, NULL, "-I ./ -cl-nv-verbose", NULL, NULL);
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
		
		FILE *out = fopen("out.ptx","w");
		fwrite(progs[0],1,sizes[0],out);
		fclose(out);
	//	printf("%s\n", progs[0]);
		//clLogPtx(prog, device, "oclConvolution.ptx");
	}

	printf("Err: %d\n", error);

	return prog;
}

bool checkBlock(){
	return false;
}

concurrent_queue con_queue;

void fermatDone(void* d, void* target){
	//TODO: d can be used to hit clMiner class
	uint32_t i,j;

	printf("finish\n");

	// read some data back from fermat
	sieveOutput_t *sout = (sieveOutput_t*)target;

	for(j=0; j < SM_COUNT; j++){
		printf("nprimes: %d\n", sout[j].nBitwins);
		for(i=0; i < sout[j].nBitwins; i++){
			if(sout[j].results.type[i]<0){
				//printf("Nothing at: %d Type: %d Hash: %d\n", i, sout->results.type[i], sout->results.hash[i]);
				continue;
			}
			clPrime *prime = new clPrime();
			mpz_init(prime->mpzOrigin);
			con_queue.push(prime);
		}
	}
}

int main() {
	int i,j,k;
	char buf[]="Hello, World!";
	size_t srcsize, worksize=strlen(buf);
	
	cl_int error;
	cl_platform_id platform;
	cl_device_id device, devicea[4];
	cl_uint platforms, devices;

	// Fetch the Platform and Device IDs; we only want one.
	error=clGetPlatformIDs(1, &platform, &platforms);
printf("Err1: %d\n", error);
	error=clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 4, devicea, &devices);
printf("Err2: %d %d\n", error, devices);
	device = devicea[devices-1];
	cl_context_properties properties[]={
		CL_CONTEXT_PLATFORM, (cl_context_properties)platform,
		0};
	// Note that nVidia's OpenCL requires the platform property
	cl_context context=clCreateContext(properties, 1, &device, NULL, NULL, &error);
printf("Err3: %d\n", error);
	cl_command_queue cq = clCreateCommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &error);
	
	cl_program prog = load_program("rot13.cl", context, device);

/*
    (   __global    sha256_context*    ctx
    ,   __global uint8_t hashOutput[32][STRIDE]	
    ,   __global primecoinBlock_t *blocks
    ,   __global uint8_t temp1[][STRIDE]		
    ,   unsigned length	 */

	GeneratePrimeTable(100000);
	uint32_t prime_size = vPrimesSize;

	// Allocate memory for the kernel to work with
	cl_mem sha_mem, input_mem, shatemp_mem, shaoutput_mem, sievetemp_mem, sieveoutput_mem, prime_mem, fermattemp_mem, fermatoutput_mem;
	cl_mem fermatoutputr_mem;
	sha_mem=clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(sha256cl_context)*SM_COUNT, NULL, &error);
	input_mem=clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(primecoinInput_t)*SM_COUNT, NULL, &error);
	shatemp_mem=clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(shaTemp_t)*SM_COUNT, NULL, &error);
	shaoutput_mem=clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(shaOutput_t)*SM_COUNT, NULL, &error);

	//Sieve stuff
	sievetemp_mem=clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(sieveTemp_t)*SM_COUNT, NULL, &error);
	printf("Err: %d\n", error);
	sieveoutput_mem=clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(sieveOutput_t)*SM_COUNT, NULL, &error);
	prime_mem=clCreateBuffer(context,CL_MEM_READ_WRITE,prime_size*sizeof(uint32_t),NULL,&error);

	//Fermat stuff
	fermattemp_mem=clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(fermatTemp_t)*SM_COUNT, NULL, &error);
	fermatoutput_mem=clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(sieveOutput_t)*SM_COUNT, NULL, &error);
	
	fermatoutputr_mem=clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(sieveOutput_t)*SM_COUNT*ROUNDS, NULL, &error);


	printf("Data bytes sha_ctx: %d\n", sizeof(sha256cl_context));
	printf("Data bytes primecoinInput: %d\n",sizeof(primecoinInput_t));
	printf("Data bytes sha_temp: %d\n", sizeof(shaTemp_t));
	printf("Data bytes sha_output: %d\n", sizeof(shaOutput_t));
	printf("Data bytes sieve_temp: %d\n", sizeof(sieveTemp_t));
	printf("Data bytes sieve_output x2: %d\n", sizeof(sieveOutput_t)*18);
	printf("Data bytes fermat_temp: %d\n", sizeof(fermatTemp_t));

	
	printf("Data bytes total: %d\n", sizeof(sha256cl_context) + sizeof(primecoinInput_t) + 
					sizeof(shaTemp_t) + sizeof(shaOutput_t) + 
					sizeof(sieveTemp_t) + 2*sizeof(sieveOutput_t) + 
					sizeof(fermatTemp_t) + ROUNDS*sizeof(sieveOutput_t));


	unsigned length = STRIDE;
	unsigned fermat_depth0 = 0;
	unsigned fermat_depth1 = 1;
	unsigned fermat_depth2 = 2;

	// get a handle and map parameters for the sha kernel	unsigned fermat_depth2 = 2;

	cl_kernel k_sha=clCreateKernel(prog, "sha", &error);
	clSetKernelArg(k_sha, 0, sizeof(sha_mem), &sha_mem);
	clSetKernelArg(k_sha, 1, sizeof(input_mem), &input_mem);
	clSetKernelArg(k_sha, 2, sizeof(shatemp_mem), &shatemp_mem);
	clSetKernelArg(k_sha, 3, sizeof(length), &length);
	clSetKernelArg(k_sha, 4, sizeof(shaoutput_mem), &shaoutput_mem);

	cl_kernel k_sieve,k_sieve_part,k_sieve_complete;
	createSieve(prog,&k_sieve,"sieve",&input_mem,&sievetemp_mem,&length,
		&sieveoutput_mem,&prime_mem,&shaoutput_mem);
	createSieve(prog,&k_sieve_part,"sieve_part",&input_mem,&sievetemp_mem,&length,
		&sieveoutput_mem,&prime_mem,&shaoutput_mem);
	createSieve(prog,&k_sieve_complete,"sieve_complete",&input_mem,&sievetemp_mem,&length,
		&sieveoutput_mem,&prime_mem,&shaoutput_mem);

	cl_kernel k_fermat, k_fermat_finish, k_fermat2, k_fermat_finish2;
 	cl_kernel k_fermato, k_fermat_finisho;

	createFermat(prog,&k_fermat,&k_fermat_finish,&input_mem,&sieveoutput_mem,&fermatoutput_mem,
			&prime_mem, &shaoutput_mem, &fermattemp_mem, &fermat_depth0);
	createFermat(prog,&k_fermat2,&k_fermat_finish2,&input_mem,&fermatoutput_mem,&sieveoutput_mem,
			&prime_mem, &shaoutput_mem, &fermattemp_mem, &fermat_depth1);

	createFermat(prog,&k_fermato,&k_fermat_finisho,&input_mem,&sieveoutput_mem,&fermatoutputr_mem,
			&prime_mem, &shaoutput_mem, &fermattemp_mem, &fermat_depth2);

	int fermat_deptho[ROUNDS];
	for(i=0; i < ROUNDS; i++){
		fermat_deptho[i] = (i << 16) + 2;
	}

	void *fermatresults = malloc(sizeof(sieveOutput_t)*SM_COUNT);
	shaOutput_t *output = new shaOutput_t[SM_COUNT];

	//Assemble the stoppable opencl queue
	clQueue mq;
	mq.add(k_sha,0); 
	mq.add(shaoutput_mem,0,sizeof(shaOutput_t)*SM_COUNT,output); 
#if 1
	mq.add(k_sieve,1);
	for(i=0; i < 1; i++){
		mq.add(k_sieve_part,1);
		mq.add(k_sieve_complete,1);
		mq.add(k_fermat,2);
		mq.add(k_fermat_finish,2);
		mq.add(k_fermat2,2);
		mq.add(k_fermat_finish2,2);

		mq.add(k_fermato,6,&fermat_deptho[i]);
		mq.add(k_fermat_finisho,6,&fermat_deptho[i]);
		mq.add(k_fermato,2);
		mq.add(k_fermat_finisho,2);
		mq.add(fermatoutputr_mem,sizeof(sieveOutput_t)*SM_COUNT*i,sizeof(sieveOutput_t)*SM_COUNT, fermatresults, fermatDone, fermatresults);
	}
#endif
	//prime_size = 100;

	uint32_t *primes = new uint32_t[prime_size];
	for(i=0; i < prime_size; i++)
		primes[i] = vPrimes[i];
	
	//Input data common to all kernels
	primecoinInput_t* linput = new primecoinInput_t[SM_COUNT];
	memset(linput,0,sizeof(primecoinInput_t));

	//TODO: this is where we pipein data from primecoin
	for(i=0; i < STRIDE; i++){
		for(j=0; j < SM_COUNT; j++){
			linput[j].blocks.nBits[i] = 0x9FEF014;
			//Fudge the stamp for moar entropy
			linput[j].blocks.timestamp[i] = j;
		}
	}

	for(i=0; i < SM_COUNT; i++){
		linput[i].nPrimorialHashFactor = 7;
		linput[i].nPrimorialMultiplier = 41;
		linput[i].nHashFactor = PrimorialFast(linput->nPrimorialHashFactor, primes, vPrimesSize);
		linput[i].nPrimes = vPrimesSize;
		linput[i].primeSeq = 4;
	}

	for(i=0; i < SM_COUNT; i++)
		initsieve(&linput[i]);

	mpz_t primorial;
	mpz_init(primorial);

	Primorial(linput->nPrimorialMultiplier,primes,vPrimesSize,primorial);
//	mpz_out_str(stdout,16,primorial);

	for(i=0; i < SM_COUNT; i++)
		copy_mpz(primorial,&linput[i].mpzPrimorial);
	//TODO: It appears that mpzPrimorial < mpzMultiplierMin is never true 

	mpz_t fixedMultiplier;
	mpz_init(fixedMultiplier);
	mpz_div_ui(fixedMultiplier,primorial,linput->nHashFactor);
//	mpz_out_str(stdout,16,fixedMultiplier);

	for(i=0; i < SM_COUNT; i++)
		copy_mpz(fixedMultiplier,&linput[i].mpzFixedMultiplier);

	worksize = STRIDE;
	cl_event event_sha, event_sieve, event_fermat;
	 
	//Launch sha kernel
	// Send input data to OpenCL (async, don't alter the buffer!)
	error=clEnqueueWriteBuffer(cq, input_mem, CL_FALSE, 0, sizeof(primecoinInput_t)*SM_COUNT, linput, 0, NULL, NULL);

	// Send input data to OpenCL (async, don't alter the buffer!)
	error=clEnqueueWriteBuffer(cq, prime_mem, CL_FALSE, 0, prime_size * sizeof(uint32_t), primes, 0, NULL, NULL);

	// Perform the operation
	mq.run(cq,checkBlock);


	// read some data back from sieve
	sieveTemp_t *stemp = (sieveTemp_t*)malloc(sizeof(sieveTemp_t));
	error=clEnqueueReadBuffer(cq, sievetemp_mem, CL_FALSE, 0, sizeof(sieveTemp_t), stemp, 0, NULL, NULL);
	sieveOutput_t *sout = (sieveOutput_t*)malloc(sizeof(sieveOutput_t)*SM_COUNT);
	error=clEnqueueReadBuffer(cq, fermatoutputr_mem, CL_FALSE, 0, sizeof(sieveOutput_t)*SM_COUNT, sout, 0, NULL, NULL);

	// read some data back from fermat
	fermatTemp_t *ftemp = (fermatTemp_t*)malloc(sizeof(fermatTemp_t)*SM_COUNT);
	error=clEnqueueReadBuffer(cq, fermattemp_mem, CL_FALSE, 0, sizeof(fermatTemp_t)*SM_COUNT, ftemp, 0, NULL, NULL);

	printf("Err: %d\n", error);

	// Await completion of all the above
	error=clFinish(cq);

	printf("Err: %d\n", error);

	// Finally, output out happy message.
//	for(i=0; i < 32; i++){
//		printf("%2.2X", buf3[i*STRIDE]);
//	}
//	printf("\n");

	for(i=0; i < 32; i++){
		printf("%d: ", i);
		for(j=0; j < 4; j++){
			uint64_t v = D_REF(output[0].mpzHash.d,j,i);
			printf("%16.16lX",v);
		}
		printf("\n");
	}

	//printf("Factor %u\n", output.factor);
	int hashes=0;
	for(i=0; i < worksize; i++){
		if(output[0].mod[i])
			hashes++;
	}

	printf("Generated: %d hashes\n", hashes);

	printf("Fixed factor: %d ", stemp->mpzFixedFactor.size[0]);
	for(i=0; i < stemp->mpzFixedFactor.size[0]; i++){
		printf("%16.16lx",stemp->mpzFixedFactor.d[(stemp->mpzFixedFactor.size[0]-i-1)*STRIDE]);
	}
	printf("\n");

#if 0
	for(i=0; i < STRIDE; i++){
		printf("%d\n", stemp->mpzFixedFactor.size[i]);
	}
#endif
	printf("nFixedFactorCombinedMod: %d\n", stemp->nFixedFactorCombinedMod[0]);
	printf("nPrimeCombined: %d\n", stemp->nPrimeCombined[0]);
	printf("nFixedInverse: %X\n", stemp->nFixedInverse[0]);

	uint32_t accum1=0,accum2=0,accum3=0;
	for(i=linput->primeSeq; i < prime_size-1; i++){
		uint32_t multiplierPos = i;
		uint32_t j;
		for(j=0; j < 1; j++){
		//	accum1 += D_REF(stemp->vCunningham1Multipliers,multiplierPos,0);
		//	accum2 += D_REF(stemp->vCunningham2Multipliers,multiplierPos,0);
			accum3 += stemp->flatCC1Multipliers[multiplierPos];
		}
	}

	printf("Accum1: %X Accum2: %X Accum3: %X\n", accum1, accum2, accum3);

	//printf("Elem10: %X Elem10: %X\n", stemp->flatCC1Multipliers[256], D_REF(stemp->vCunningham1Multipliers,256,0));


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
#if 0
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
#endif
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

	uint32_t fprimes=0;
	uint32_t candidates=0;
	for(j=0; j < SM_COUNT; j++){
		for(i=0; i < sout[j].nBitwins; i++){
			if(sout[j].results.type[i]<0){
				//printf("Nothing at: %d Type: %d Hash: %d\n", i, sout->results.type[i], sout->results.hash[i]);
				continue;
			}

			uint32_t mult = sout[j].results.mults[i];
			mpz_t m,hash;
			mpz_init(m);
			mpz_init(hash);
			mpz_mul_ui(m,fixedMultiplier,mult);

			if(sout[j].results.hash[i] >= 512){
				printf("Opps: %d, %d\n", i, sout[j].results.hash[i]);
			}

			copyBN(sout[j].results.hash[i],hash,&output[j].mpzHash);
			mpz_mul(m,m,hash);

			mpz_t p,e,two,r;
			mpz_init(p);
			mpz_init(e);
			mpz_init(two);
			mpz_init(r);
			if(sout[j].results.type[i] == CC2_TYPE)
				mpz_add_ui(p,m,1);
			else
				mpz_sub_ui(p,m,1);

			mpz_sub_ui(e,p,1);

			mpz_set_ui(two,2);

			mpz_powm(r,two,e,p);
		
			if(0){//(mpz_cmp_ui(r,1)==0) != sout->results.prime[i]){
				printf("System failure at %d Type: %d Mult: %d Hash: %d Reported as: %d\n", i, 
							sout[j].results.type[i], 
							sout[j].results.mults[i], sout[j].results.hash[i],
							sout[j].results.prime[i]);

				cout << sout[j].results.type[i] << "," << sout[j].results.mults[i] << "," << hex << p << dec << "\n";
#if 0
				print_bn((mpzcl_t*)&ftemp->mpzM,"M", hashidx[i]);
				print_bn((mpzcl_t*)&ftemp->mpzR,"R", hashidx[i]);
				print_bn((mpzcl_t*)&ftemp->mpzInv,"I", hashidx[i]);
				print_bn((mpzcl_t*)&ftemp->mpzXbinU,"U", hashidx[i]);
				print_bn((mpzcl_t*)&ftemp->mpzXbinV,"V", hashidx[i]);
				print_bn((mpzcl_t*)&ftemp->mpzH,"H", hashidx[i]);
				print_bn((mpzcl_t*)&ftemp->mpzSqMod,"SqMod", hashidx[i]);
				print_bn((mpzcl_t*)&ftemp->mpzSq,"Sq", hashidx[i]);
				print_bn((mpzcl_t*)&ftemp->mpzBarrettT2,"T2", hashidx[i]);
				print_bn((mpzcl_t*)&ftemp->mpzBarrettT3,"T3", hashidx[i]);
				print_bn((mpzcl_t*)&ftemp->mpzBarrettT4,"T4", hashidx[i]);
#endif
			}

			if(mpz_cmp_ui(r,1)==0){
				//	printf("%d,",sout->results.mults[i]);
				//cout << sout->results.type[i] << "," << sout->results.mults[i] << "," << hex << m << dec << "\n";
				fprimes++;
			}
			candidates++;

		 	//printf("%d,",mpz_sizeinbase (p, 2));
		
			//	fprimes++;

			//cout << sout->results.type[i] << "," << hex << m << dec << "\n";

			//TODO: put this shit in a list

			//mulBN2(tid,sieve->hash[i],&temp->mpzChainOrigin,&temp->mpzFixedFactor,sieve->mults[i]);
//			if(sout->hash[i]!=0 && sout->type[i] != -1)
//				continue;
//			printf("Multiplier: %d, Hash: %d, Type: %d, Idx: %d\n", sout->mults[i], sout->hash[i], sout->type[i], i);
		}
	}

	printf("Found %d primes out of %d\n", fprimes, candidates);

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
