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

#include "global.h"
#include "clprime.h"

#include "cldefs.h"

extern commandlineInput_t commandlineInput;	

void clQueue::add(cl_kernel k, int evt, clTime* t, void (*complete)(void*,void*), void *data){
	clWork *work = new clWork();
	work->type = TYPE_KERN;
	work->kernel = k;
	work->evt = evt;
	work->complete = complete;
	work->data = data;
	work->target = 0;
	work->time = t;

	works.push_back(work);
}

void clQueue::add(cl_kernel k, int arg, int *argval, clTime* t){
	clWork *work = new clWork;
	work->type = TYPE_ARG;
	work->kernel = k;
	work->arg = arg;
	work->argval = argval;
	work->complete = 0;
	work->time = t;

	works.push_back(work);
}

void clQueue::add(cl_mem mem, uint32_t ofst, uint32_t size, void* target, void (*complete)(void*,void*), void *data){
	clWork *work = new clWork;
	work->type = TYPE_MEM;
	work->mem = mem;
	work->ofst = ofst;
	work->size = size;
	work->complete = complete;
	work->data = data;
	work->target = target;
	work->time = 0;

	works.push_back(work);
}

//Instead of putting everything in queue and running, we try to keep one item in async queue to give good performance
//but still allow the queue to be cleaned and done away with in short order. 
bool clQueue::run(cl_command_queue cq, bool (*pred)(void)){
	size_t worksize = STRIDE*SM_COUNT;
	size_t localsize = STRIDE;
	cl_event event_odd, event_even;
	int idx=0;
	int i=0;
	int last_evt=0;
	int error=0;
	//TODO: check predicate at each iteration
	for(vector<clWork*>::iterator it = works.begin(); it!= works.end() && !pred(); it++){
		clWork* work = *it;
		if(work->type==TYPE_KERN||work->type==TYPE_MEM){
			if(work->type==TYPE_KERN){
				error = clEnqueueNDRangeKernel(cq, work->kernel, 1, NULL, &worksize, &localsize, 0, NULL,  idx%2?&event_odd:&event_even);
				CheckErr("enqueue_kernel");
			}else{
				error = clEnqueueReadBuffer(cq, work->mem, CL_FALSE, work->ofst, work->size, work->target, 0, NULL, idx%2?&event_odd:&event_even);
				CheckErr("enqueue_read");
			}
			error = clFlush(cq);
			CheckErr("clQueue_flush");
			if(idx>0){
				error = clWaitForEvents(1,idx%2?&event_even:&event_odd);
				if(error){
					printf("Error: %d\n", error);	
					printf("Type: %d\n", works[last_evt]->type);
					if(works[last_evt]->time)
						printf("Name: %s\n", works[last_evt]->time->name);
				}else if(works[last_evt]->time && commandlineInput.profile) {
					works[last_evt]->time->count++;
					works[last_evt]->time->time+=get_exec_time(idx%2?&event_even:&event_odd);
					printf("%s: Exec time %fms, %f, %d\n", works[last_evt]->time->name, works[last_evt]->time->time/works[last_evt]->time->count, works[last_evt]->time->time,works[last_evt]->time->count);
				}
				CheckErr("clQueue_wait");

				clReleaseEvent(idx%2?event_even:event_odd);
				if(works[last_evt]->complete)
					works[last_evt]->complete(works[last_evt]->data,works[last_evt]->target);
			}
			last_evt = i;
			idx++;
		}
		if(work->type==TYPE_ARG){
			clSetKernelArg(work->kernel,work->arg,sizeof(*work->argval),work->argval);
		}
		i++;
	}
	//printf("Cleanup\n");
		//printf("Type:%d\n", works[last_evt]->type);
	if(idx>0){		
		clReleaseEvent(idx%2?event_even:event_odd);
		clFinish(cq);
		if(works[last_evt]->complete)
			works[last_evt]->complete(works[last_evt]->data,works[last_evt]->target);
	}
	return pred();
}
