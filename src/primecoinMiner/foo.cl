
typedef unsigned int uint32_t;
typedef unsigned long uint64_t;

#define STRIDE 512

#define D_REF(a,b,t) a[(b)*STRIDE+t]

#if 0
__kernel void memcpy_32(uint32_t tid, __global uint32_t *restrict dst, uint32_t dofst, __global uint32_t *restrict src, uint32_t sofst, uint32_t length){
	uint32_t i;
	uint32_t save  = D_REF(dst,dofst/4,tid) & (uint32_t)(0xFFFFFFFFUL >> ((4-(dofst%4)) * 8) );
	uint32_t ssave = (uint64_t)D_REF(src,sofst/4,tid) >> ((sofst%4) * 8);

//	printf("save: %8.8X %8.8X\n", save, ssave);

	for(i=0; i < length/4 + (length%4?1:0); i++){
		uint32_t sword = D_REF(src,1+i+sofst/4,tid);
		uint32_t tword = (sword << ((4 - (sofst%4))*8)) | ssave;
		D_REF(dst,i+dofst/4,tid) = save | (tword << ((dofst % 4) * 8));
		save = ((uint64_t)tword >> ((4 - (dofst % 4)) * 8)) ;
		ssave = (uint64_t)sword >> ((sofst%4)*8);
//		printf("tick: %8.8X %8.8X %8.8X %8.8X\n",sword, tword, save,ssave);
	}
}
#endif

#if 0
__kernel void if_eq(__global unsigned long* out, unsigned  arg0, unsigned long arg1)
{
#if 1
	unsigned i=0;
	for(i = 0; i < len; i++){
   		 out[i] =  (arg0 + arg1) >> args;
   	}
#else

	unsigned id = get_global_id(0);
  	out[id] = id*2;
#endif
}
#endif

#if 1
__kernel void if_eq(__global unsigned * out, unsigned  arg0, unsigned arg1)
{
	*out = arg0 == arg1 ? 1 : 0;
}
#endif
