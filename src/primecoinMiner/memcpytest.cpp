#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define STRIDE 1
#define D_REF(a,b,t) a[(b)*STRIDE+t]


void memcpy_foo(uint32_t tid, uint32_t dst[], uint32_t dofst, uint32_t src[], uint32_t sofst, uint32_t length){
	uint32_t i;
	uint32_t save  = D_REF(dst,dofst/4,tid) & (uint32_t)(0xFFFFFFFFUL >> ((4-(dofst%4)) * 8) );
	uint32_t ssave = (uint64_t)D_REF(src,sofst/4,tid) >> ((sofst%4) * 8);

	printf("save: %8.8X %8.8X\n", save, ssave);

	for(i=0; i < length/4 + (length%4?1:0); i++){
		uint32_t sword = D_REF(src,1+i+sofst/4,tid);
		uint32_t tword = ((uint64_t)sword << ((4 - (sofst%4))*8)) | ssave;
		D_REF(dst,i+dofst/4,tid) = save | (tword << ((dofst % 4) * 8));
		save = ((uint64_t)tword >> ((4 - (dofst % 4)) * 8)) ;
		ssave = (uint64_t)sword >> ((sofst%4)*8);
		printf("tick: %8.8X %8.8X %8.8X %8.8X\n",sword, tword, save,ssave);
	}
}

int main(){
	char *input = new char[256];
	int i;
	for(i = 0; i < 256; i++){
		input[i] = i%52 > 25 ? 'A' + i % 26 : 'a' + i%26;
	}

	uint32_t *inputint = (uint32_t*)input;

	uint32_t *outputint = new uint32_t[256/4];

	for(i=0; i < 5; i++){
		memcpy_foo(0,outputint,i*32,inputint,1,8);
	}

	char *output = (char*)outputint;	
	
	for(i=0; i < 256; i++){
		if(output[i] == 0)
			printf("0");
		else			
			printf("%c",output[i]);
	}

	printf("\n");

	return 0;
}
