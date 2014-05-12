
#define STRIDE 512

#define SIEVE_SIZE (1024*1024)
#define SIEVE_WORDS (SIEVE_SIZE / 32)

#define PRIME_PCT 100
#define NTOTAL_PRIMES 9592
#define NPRIMES ((NTOTAL_PRIMES * PRIME_PCT) / 100)
#define NCHAIN_LENGTH 10

#define LAYERS NCHAIN_LENGTH
#define MULT_SIZE (10240)
#define MULT_WORDS (MULT_SIZE)

#define D_REF(a,b,t) a[(b)*STRIDE+t]

#define nFractionalBits 24
#define TARGET_FRACTIONAL_MASK  ((1u<<nFractionalBits) - 1)
#define TARGET_LENGTH_MASK  (~TARGET_FRACTIONAL_MASK)

#define RESULT_SIZE (64*1024)

#define BITWIN_TYPE 0
#define CC1_TYPE 1
#define CC2_TYPE 2

typedef struct {
 	uint64_t v[4]; 
} uint256_t;

typedef struct  
{
	uint32_t total[2*STRIDE];
	uint32_t state[8*STRIDE];
	uint32_t buffer[16*STRIDE];
	uint32_t W[64*STRIDE];
	uint32_t msglen[2*STRIDE];
	uint32_t padding[16*STRIDE];  
} sha256_context;
 

typedef struct  
{
	uint32_t	version[STRIDE];
	uint32_t	prevBlockHash[8*STRIDE];
	uint32_t	merkleRoot[8*STRIDE];
	uint32_t	timestamp[STRIDE];
	uint32_t	nBits[STRIDE];
	uint32_t	nonce[STRIDE];
	// GetHeaderHash() goes up to this offset (4+32+32+4+4+4=80 bytes)
	uint256_t blockHeaderHash;
	//CBigNum bnPrimeChainMultiplierBN; unused
/*	mpz_class mpzPrimeChainMultiplier;
	// other
	serverData_t serverData;
	uint32 threadIndex; // the index of the miner thread
	bool xptMode;*/
}primecoinBlock_t;

typedef struct {
	uint32_t size[STRIDE];
	uint64_t d[8*STRIDE];
} mpzcl_t;

typedef struct {
	uint32_t size[STRIDE];
	uint64_t d[16*STRIDE];
} mpzcll_t;

typedef struct {
	uint32_t size[STRIDE];
	uint64_t d[24*STRIDE];
} mpzclel_t;

typedef struct {
	uint64_t d[8];
	uint32_t size;
} mpzcls_t;

typedef struct {
	uint32_t size[STRIDE];
	uint32_t d[16*STRIDE];
} mpzcl32_t;

typedef struct {
	uint32_t vfBitwin[SIEVE_WORDS*STRIDE];
	uint32_t vfCC1[SIEVE_WORDS*STRIDE];
	uint32_t vfCC2[SIEVE_WORDS*STRIDE];

	uint32_t vCunningham1Multipliers[MULT_WORDS*STRIDE];
	uint32_t vCunningham2Multipliers[MULT_WORDS*STRIDE];

	mpzcl_t mpzFixedFactor;
	mpzcl_t mpzQuot;

	uint32_t nFixedFactorCombinedMod[STRIDE];
	uint32_t nPrimeCombined[STRIDE];
	uint32_t nFixedInverse[MULT_WORDS*STRIDE];

	uint32_t flatMultipliers[MULT_WORDS*STRIDE];
} sieveTemp_t;

typedef struct {
	primecoinBlock_t blocks;

	uint32_t nPrimorialHashFactor;
	uint32_t nHashFactor;
	uint32_t nPrimorialMultiplier;
	uint32_t nOverrideTargetValue;
	uint32_t nOverrideBTTargetValue;
	uint32_t lSieveTarget;
	uint32_t lSieveBTTarget;

	uint32_t nPrimes;
	uint32_t primeSeq; 
	
	mpzcls_t mpzPrimorial;
	mpzcls_t mpzFixedMultiplier;
} primecoinInput_t;

typedef struct { 
	uint32_t mults[RESULT_SIZE];
	uint32_t hash[RESULT_SIZE];
	uint32_t type[RESULT_SIZE];

	uint32_t nBitwins;
	uint32_t nCC1s;
	uint32_t nCC2s;
} sieveOutput_t;

typedef struct {
	uint32_t mod[STRIDE];
	mpzcl_t mpzHash;   
} shaOutput_t;

typedef struct {
	uint32_t hash[8*STRIDE];
	uint32_t hashTemp[20*STRIDE];
	mpzcl_t mpzPrimeMin;
} shaTemp_t;

typedef struct {
	mpzcl_t mpzM;
	mpzcl_t mpzE;
	mpzcl_t mpzR;
	mpzcl_t mpzHalfR;
	mpzcl_t mpzInv;
	mpzcl_t mpzOne;
	mpzcl_t mpzH;

	mpzcl32_t mpzResult;
	mpzcl32_t mpzBase;
	mpzcl32_t mpzH32;
	mpzcl32_t mpzM32;
	mpzcl32_t mpzV32;
	mpzcl32_t mpzMonProT;
	mpzcl32_t mpzOne32;

	mpzcl_t mpzXbinU;
	mpzcl_t mpzXbinV;
	mpzcl_t mpzXbinTemp;
	mpzcl_t mpzXbinTemp2;
	mpzcl_t mpzXbinTemp3;

	mpzcl_t mpzSqMod;
	mpzcll_t mpzSq;

	mpzclel_t mpzBarrettT1;
	mpzclel_t mpzBarrettT2;
	mpzclel_t mpzBarrettT3;
	mpzclel_t mpzBarrettT4;
	mpzcl_t mpzBarrettN;
	mpzcl_t mpzBarrettA;
	mpzcl_t mpzBarrettM;

	mpzcl_t mpzNewtonDen;
	mpzcl_t mpzNewtonX1;
	mpzcl_t mpzNewtonX2;
	mpzcl_t mpzNewtonX1s;
	mpzcl_t mpzNewtonX2s;
	mpzcl_t mpzNewtonX;
	mpzclel_t mpzNewtonDenP;
	mpzcl_t mpzNewtonDenPS;
	mpzcl_t mpzNewton2;
	mpzcl_t mpzNewton2s;
	mpzcl_t mpzNewtonDiff;
	mpzclel_t mpzNewtonProd;
} fermatTemp_t;

