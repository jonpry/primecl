#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>


#include "sha256.h"

#include "cldefs.h"

#define __global	
#define __constant
#define restrict
#include "sha256.cl"

//d8835751e651bba83e5a93e5c002ff1b7ff946a4687d76395947fa618eb56f5b - sha of 80 zeros
//14508459B221041EAB257D2BAAA7459775BA748246C8403609EB708F0E57E74B - primecoin hash of all zeros

uint8_t do_sha(){
	uint8_t block[80] = {0};
	uint8_t digest[32];
	uint8_t digest2[32];

	block[72] = 1;

	//*(uint32_t*)&block[72] = 0x9FEF014;
	//*(uint32_t*)&block[76] = 0x2a400100;

	//2a400100

	sha256_context ctx;

	sha256_starts(&ctx);
	sha256_update(&ctx,block,80);
	sha256_finish(&ctx,digest);

	sha256_starts(&ctx);
	sha256_update(&ctx,digest,32);
	sha256_finish(&ctx,digest2);

#if 1
	int i;
	for(i=0; i < 32; i++){
		printf("%2.2X",digest[31-i]);
	}	
	printf("\n");

	for(i=0; i < 32; i++){
		printf("%2.2X",digest2[31-i]);
	}	
	printf("\n");
#endif
	return digest[0];

}

void do_gpu_sha(){
	int i;
	uint32_t digest[8*STRIDE];
	uint32_t block[20*STRIDE] = {0};

	block[(72/4)*STRIDE] = 0x9FEF014;
	block[(76/4)*STRIDE] = 0x2a400100;

	sha256cl_context ctx;
	sha256_init(0,&ctx);

	sha256_starts(0,&ctx);
	sha256_update(0,&ctx,block,80);
	sha256_finish(0,&ctx,digest);

	for(i=0; i < 8; i++)
		printf("%8.8X",D_REF(digest,7-i,0));
	printf("\n");

	sha256_starts(0,&ctx);
	sha256_update(0,&ctx,digest,32);
	sha256_finish(0,&ctx,digest);

	for(i=0; i < 8; i++)
		printf("%8.8X",D_REF(digest,7-i,0));
	printf("\n");
}

int main(){
	int i;
	for(i=0; i < 1; i++){
		do_sha();
		do_gpu_sha();
	}

	return 0;
}
