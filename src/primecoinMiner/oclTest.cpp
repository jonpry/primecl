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

int main() {
	char buf[]="Hello, World!";
	size_t srcsize, worksize=strlen(buf);
	
	cl_int error;
	cl_platform_id platform;
	cl_device_id device, devicea[4];
	cl_uint platforms, devices;

	// Fetch the Platform and Device IDs; we only want one.
	error=clGetPlatformIDs(1, &platform, &platforms);
printf("Err1: %d\n", error);
	error=clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, 4, devicea, &devices);
printf("Err2: %d %d\n", error, devices);
	device = devicea[0];
	cl_context_properties properties[]={
		CL_CONTEXT_PLATFORM, (cl_context_properties)platform,
		0};
	// Note that nVidia's OpenCL requires the platform property
	cl_context context=clCreateContext(properties, 1, &device, NULL, NULL, &error);
printf("Err3: %d\n", error);
	cl_command_queue cq = clCreateCommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &error);
	
	const char *src=loadfile("foo.amdil",&srcsize);

	printf("%p\n",src);

	const char *srcptr[]={src};
	// Submit the source code of the rot13 kernel to OpenCL
	cl_program prog=clCreateProgramWithBinary(context,
		 1, devicea, &srcsize, (const unsigned char**)srcptr, 0, &error);
	printf("Err: %d\n", error);

	// and compile it (after this we could extract the compiled version)
	error=clBuildProgram(prog, 1, devicea, "-I ./ ", NULL, NULL);
	delete src;

	if (error < 0) {
    		// Determine the size of the log
    		size_t log_size;
    		clGetProgramBuildInfo(prog, device, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);

    		// Allocate memory for the log
    		char *log = (char *) malloc(log_size);

    		// Get the log
    		clGetProgramBuildInfo(prog, device, CL_PROGRAM_BUILD_LOG, log_size, log, NULL);

    		// Print the log
    		printf("%s\n", log);
	}else{
		size_t sizes[16];
		clGetProgramInfo(prog, CL_PROGRAM_BINARY_SIZES, 8, &sizes, 0);
		printf("%ld\n", sizes[0]);
		char **progs = new char *;
		progs[0] = new char[sizes[0]];
		clGetProgramInfo(prog, CL_PROGRAM_BINARIES, 8, progs, 0);
		
	//	printf("%s\n", progs[0]);
		//clLogPtx(prog, device, "oclConvolution.ptx");
	}

	printf("Err: %d\n", error);

/*
    (   __global    sha256_context*    ctx
    ,   __global uint8_t hashOutput[32][STRIDE]	
    ,   __global primecoinBlock_t *blocks
    ,   __global uint8_t temp1[][STRIDE]		
    ,   unsigned length	 */

	// Allocate memory for the kernel to work with
	cl_mem temp_mem;
	temp_mem=clCreateBuffer(context, CL_MEM_READ_WRITE, 64, NULL, &error);


	unsigned length = 512;

	int a = 5;
	int b = 6;

	// get a handle and map parameters for the kernel
	cl_kernel k_rot13=clCreateKernel(prog, "if_eq", &error);
	clSetKernelArg(k_rot13, 0, sizeof(temp_mem), &temp_mem);
	clSetKernelArg(k_rot13, 1, sizeof(a), &a);
	clSetKernelArg(k_rot13, 2, sizeof(b), &b);	



	// Target buffer just so we show we got the data from OpenCL
	unsigned int buf2[64] = {5};
	
	// Perform the operation
	worksize = 1;
	cl_event event;
	error=clEnqueueNDRangeKernel(cq, k_rot13, 1, NULL, &worksize, &worksize, 0, NULL,  &event);
	// Read the result back into buf2
	printf("Err: %d\n", error);
	
	error=clEnqueueReadBuffer(cq, temp_mem, CL_FALSE, 0, 64, &buf2, 0, NULL, NULL);
	printf("Err: %d\n", error);

	// Await completion of all the above
	error=clFinish(cq);

	printf("Err: %d\n", error);
	
	printf("Generated: %d\n", buf2[0]);
}
