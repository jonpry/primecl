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

#ifndef CLPRIME_H
#define CLPRIME_H

#endif //CLPRIME_H

#include <CL/cl.h>
#include <gmp.h>
#include <iostream>
#include <queue>
#include <pthread.h>
#include <time.h>
#include <cstdlib>

using namespace std;

void oclInit();

#define TYPE_ARG 1
#define TYPE_KERN 2
#define TYPE_MEM 3

#define CheckErr(x) do { \
		if(error){ \
			printf("OpenCL error at %s: %d\n", x, error); \
			exit(-1); \
		} \
	} while(0)

class clTime{
public:
	clTime(const char *n) {name = n;time=0;count=0;}
	const char* name;
	double time;
	uint64_t count;
};

class clWork{
public:
	int type;
	cl_kernel kernel;
	int evt;
	void(*complete)(void*,void*);
	void *data;
	int arg;
	int* argval;
	cl_mem mem;
	uint32_t ofst, size;
	void *target;
	clTime* time;
};

class clQueue{
public:
	void add(cl_kernel k, int evt, clTime*,void(*complete)(void*,void*)=0,void *data=0);
	void add(cl_kernel k, int arg, int *argval, clTime* t);
	void add(cl_mem mem, uint32_t ofst, uint32_t size, void* target, void(*complete)(void*,void*)=0,void *data=0);
	bool run(cl_command_queue, bool(*pred)(void));

	vector<clWork*> works;
};

class concurrent_queue
{
    private:
        queue<void*> _queue_;
        pthread_mutex_t push_mutex;
        pthread_mutex_t pop_mutex;
        pthread_cond_t cond;
    public:
        concurrent_queue()
        {
            pthread_mutex_init(&push_mutex, NULL);
            pthread_mutex_init(&pop_mutex, NULL);
            pthread_cond_init(&cond, NULL);
        }
 
        void push(void* data)
        {
            pthread_mutex_lock(&push_mutex);
            _queue_.push(data);
            pthread_cond_signal(&cond);
            pthread_mutex_unlock(&push_mutex);
        }
 
        void* pop()
        {
            pthread_mutex_lock(&pop_mutex);
            while(_queue_.empty())
                pthread_cond_wait(&cond, &pop_mutex);
            void *popped_value = _queue_.front();
            _queue_.pop();
            pthread_mutex_unlock(&pop_mutex);
	    return popped_value;
        }

	bool empty(){
            pthread_mutex_lock(&pop_mutex);
	    bool ret = _queue_.empty();
            pthread_mutex_unlock(&pop_mutex);
	    return ret;
	}
};

void print_exec_time(const char* name, cl_event *event);
double get_exec_time(cl_event* event);

struct primecoinInput;
struct primecoinBlock;
struct fermatTemp;
struct sieveTemp;
struct sieveOutput;
struct shaOutput;

class clMiner{
public:
	void oclInit();	
	void mine(struct primecoinBlock *);
	void fermatDone(struct sieveOutput*);
	void chainFunc();
	cl_context context;
	cl_device_id device;


	clQueue mq;
	cl_program prog;
	cl_mem sha_mem, input_mem, shatemp_mem, shaoutput_mem;
	cl_mem sievetemp_mem, sieveoutput_mem, prime_mem;
	cl_mem fermattemp_mem, fermatoutput_mem, fermatoutputr_mem;
	cl_mem input_mem_host, fermatoutputr_mem_host;
	cl_kernel k_sha;
	cl_kernel k_sieve,k_sieve_part,k_sieve_complete;
	cl_kernel k_fermat, k_fermat_finish, k_fermat2;
	cl_kernel k_fermat_finish2;
	cl_kernel k_fermato, k_fermat_finisho;
	cl_command_queue cq;
	uint32_t *primes;
	primecoinInput* linput;
	uint32_t prime_size;
	mpz_t fixedMultiplier,primorial;
	mpz_t lastBlockOrigin;
	mpz_t lastOrigin;

	concurrent_queue con_queue;

	struct fermatTemp *ftemp;
	struct sieveTemp *stemp;
	struct sieveOutput *sout;
	struct shaOutput *output;
	void * fermatresults;

	int fermat_deptho[16];

	unsigned length;
	unsigned fermat_depth0;
	unsigned fermat_depth1;
	unsigned fermat_depth2;
	unsigned fermat_depth3;
};

class clPrime {
public:
	struct primecoinBlock *block;
	mpz_t mpzOrigin;	
	int type;
	uint32_t mult;
	uint32_t sm;
	uint32_t tid;
};
 
