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

#define GET_UINT8(a,b,t) ((D_REF(a,b>>2,t) >> ((b%4)*8)) & 0xFF)

#define GET_UINT32(n,w,b,i,t,o)                       \
{                                               \
    D_REF(n,w,t) = (GET_UINT8(b,i +   o,t) << 24 )       \
                 | (GET_UINT8(b,i + 1+o,t) << 16 )       \
                 | (GET_UINT8(b,i + 2+o,t) <<  8 )       \
                 | (GET_UINT8(b,i + 3+o,t) );      \
}

#define PUT_UINT32(n,b,i,t)                       \
{                                               \
    D_REF(b,i/4,t) = ((((n) >> 24 ) & 0xFF)      )  \
    		   | ((((n) >> 16 ) & 0xFF) << 8 )    \
    		   | ((((n) >>  8 ) & 0xFF) << 16)     \
    		   | ((((n)       ) & 0xFF) << 24);      \
}


__constant static uint32_t sha256_padding[16] =
{
 0x80UL, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

void sha256_init( uint32_t tid, __global sha256cl_context *restrict ctx){
	int i;
	for(i=0; i < 16; i++)
		D_REF(ctx->padding,i,tid) = sha256_padding[i];
}

void sha256_starts( unsigned tid, __global sha256cl_context *restrict volatile ctx )
{
    D_REF(ctx->total,0,tid) = 0;
    D_REF(ctx->total,1,tid) = 0;

    D_REF(ctx->msglen,0,tid) = 0;
    D_REF(ctx->msglen,1,tid) = 0;

    D_REF(ctx->state,0,tid) = 0x6A09E667;
    D_REF(ctx->state,1,tid) = 0xBB67AE85;
    D_REF(ctx->state,2,tid) = 0x3C6EF372;
    D_REF(ctx->state,3,tid) = 0xA54FF53A;
    D_REF(ctx->state,4,tid) = 0x510E527F;
    D_REF(ctx->state,5,tid) = 0x9B05688C;
    D_REF(ctx->state,6,tid) = 0x1F83D9AB;
    D_REF(ctx->state,7,tid) = 0x5BE0CD19;
}


void sha256_process( uint32_t tid, __global sha256cl_context *restrict ctx, __global uint32_t *restrict data, uint32_t ofst )
{
    uint32_t temp1, temp2;
    uint32_t A, B, C, D, E, F, G, H;

//TODO: this is a really shitty way to do this. should do something like memcpy that does auto be32 conversion
    GET_UINT32( ctx->W, 0,  data,  0 , tid, ofst);
    GET_UINT32( ctx->W, 1,  data,  4 , tid, ofst);
    GET_UINT32( ctx->W, 2,  data,  8 , tid, ofst);
    GET_UINT32( ctx->W, 3,  data, 12 , tid, ofst);
    GET_UINT32( ctx->W, 4,  data, 16 , tid, ofst);
    GET_UINT32( ctx->W, 5,  data, 20 , tid, ofst);
    GET_UINT32( ctx->W, 6,  data, 24 , tid, ofst);
    GET_UINT32( ctx->W, 7,  data, 28 , tid, ofst);
    GET_UINT32( ctx->W, 8,  data, 32 , tid, ofst);
    GET_UINT32( ctx->W, 9,  data, 36 , tid, ofst);
    GET_UINT32( ctx->W, 10, data, 40 , tid, ofst);
    GET_UINT32( ctx->W, 11, data, 44 , tid, ofst);
    GET_UINT32( ctx->W, 12, data, 48 , tid, ofst);
    GET_UINT32( ctx->W, 13, data, 52 , tid, ofst);
    GET_UINT32( ctx->W, 14, data, 56 , tid, ofst);
    GET_UINT32( ctx->W, 15, data, 60 , tid, ofst);

#define  SHR(x,n) ((x & 0xFFFFFFFF) >> n)
#define ROTR(x,n) (SHR(x,n) | (x << (32 - n)))

#define S0(x) (ROTR(x, 7) ^ ROTR(x,18) ^  SHR(x, 3))
#define S1(x) (ROTR(x,17) ^ ROTR(x,19) ^  SHR(x,10))

#define S2(x) (ROTR(x, 2) ^ ROTR(x,13) ^ ROTR(x,22))
#define S3(x) (ROTR(x, 6) ^ ROTR(x,11) ^ ROTR(x,25))

#define F0(x,y,z) ((x & y) | (z & (x | y)))
#define F1(x,y,z) (z ^ (x & (y ^ z)))

#define R(t)                                    \
(                                               \
    D_REF(ctx->W,t,tid) = S1(D_REF(ctx->W,t-2,tid)) + D_REF(ctx->W,t-7,tid) +          \
           S0(D_REF(ctx->W,t-15,tid)) + D_REF(ctx->W,t-16,tid)            \
)

#define P(a,b,c,d,e,f,g,h,x,K)                  \
{                                               \
    temp1 = h + S3(e) + F1(e,f,g) + K + x;      \
    temp2 = S2(a) + F0(a,b,c);                  \
    d += temp1; h = temp1 + temp2;              \
}

    A = D_REF(ctx->state,0,tid);
    B = D_REF(ctx->state,1,tid);
    C = D_REF(ctx->state,2,tid);
    D = D_REF(ctx->state,3,tid);
    E = D_REF(ctx->state,4,tid);
    F = D_REF(ctx->state,5,tid);
    G = D_REF(ctx->state,6,tid);
    H = D_REF(ctx->state,7,tid);

    P( A, B, C, D, E, F, G, H, D_REF(ctx->W, 0,tid), 0x428A2F98 );
    P( H, A, B, C, D, E, F, G, D_REF(ctx->W, 1,tid), 0x71374491 );
    P( G, H, A, B, C, D, E, F, D_REF(ctx->W, 2,tid), 0xB5C0FBCF );
    P( F, G, H, A, B, C, D, E, D_REF(ctx->W, 3,tid), 0xE9B5DBA5 );
    P( E, F, G, H, A, B, C, D, D_REF(ctx->W, 4,tid), 0x3956C25B );
    P( D, E, F, G, H, A, B, C, D_REF(ctx->W, 5,tid), 0x59F111F1 );
    P( C, D, E, F, G, H, A, B, D_REF(ctx->W, 6,tid), 0x923F82A4 );
    P( B, C, D, E, F, G, H, A, D_REF(ctx->W, 7,tid), 0xAB1C5ED5 );
    P( A, B, C, D, E, F, G, H, D_REF(ctx->W, 8,tid), 0xD807AA98 );
    P( H, A, B, C, D, E, F, G, D_REF(ctx->W, 9,tid), 0x12835B01 );
    P( G, H, A, B, C, D, E, F, D_REF(ctx->W,10,tid), 0x243185BE );
    P( F, G, H, A, B, C, D, E, D_REF(ctx->W,11,tid), 0x550C7DC3 );
    P( E, F, G, H, A, B, C, D, D_REF(ctx->W,12,tid), 0x72BE5D74 );
    P( D, E, F, G, H, A, B, C, D_REF(ctx->W,13,tid), 0x80DEB1FE );
    P( C, D, E, F, G, H, A, B, D_REF(ctx->W,14,tid), 0x9BDC06A7 );
    P( B, C, D, E, F, G, H, A, D_REF(ctx->W,15,tid), 0xC19BF174 );
    P( A, B, C, D, E, F, G, H, R(16), 0xE49B69C1 );
    P( H, A, B, C, D, E, F, G, R(17), 0xEFBE4786 );
    P( G, H, A, B, C, D, E, F, R(18), 0x0FC19DC6 );
    P( F, G, H, A, B, C, D, E, R(19), 0x240CA1CC );
    P( E, F, G, H, A, B, C, D, R(20), 0x2DE92C6F );
    P( D, E, F, G, H, A, B, C, R(21), 0x4A7484AA );
    P( C, D, E, F, G, H, A, B, R(22), 0x5CB0A9DC );
    P( B, C, D, E, F, G, H, A, R(23), 0x76F988DA );
    P( A, B, C, D, E, F, G, H, R(24), 0x983E5152 );
    P( H, A, B, C, D, E, F, G, R(25), 0xA831C66D );
    P( G, H, A, B, C, D, E, F, R(26), 0xB00327C8 );
    P( F, G, H, A, B, C, D, E, R(27), 0xBF597FC7 );
    P( E, F, G, H, A, B, C, D, R(28), 0xC6E00BF3 );
    P( D, E, F, G, H, A, B, C, R(29), 0xD5A79147 );
    P( C, D, E, F, G, H, A, B, R(30), 0x06CA6351 );
    P( B, C, D, E, F, G, H, A, R(31), 0x14292967 );
    P( A, B, C, D, E, F, G, H, R(32), 0x27B70A85 );
    P( H, A, B, C, D, E, F, G, R(33), 0x2E1B2138 );
    P( G, H, A, B, C, D, E, F, R(34), 0x4D2C6DFC );
    P( F, G, H, A, B, C, D, E, R(35), 0x53380D13 );
    P( E, F, G, H, A, B, C, D, R(36), 0x650A7354 );
    P( D, E, F, G, H, A, B, C, R(37), 0x766A0ABB );
    P( C, D, E, F, G, H, A, B, R(38), 0x81C2C92E );
    P( B, C, D, E, F, G, H, A, R(39), 0x92722C85 );
    P( A, B, C, D, E, F, G, H, R(40), 0xA2BFE8A1 );
    P( H, A, B, C, D, E, F, G, R(41), 0xA81A664B );
    P( G, H, A, B, C, D, E, F, R(42), 0xC24B8B70 );
    P( F, G, H, A, B, C, D, E, R(43), 0xC76C51A3 );
    P( E, F, G, H, A, B, C, D, R(44), 0xD192E819 );
    P( D, E, F, G, H, A, B, C, R(45), 0xD6990624 );
    P( C, D, E, F, G, H, A, B, R(46), 0xF40E3585 );
    P( B, C, D, E, F, G, H, A, R(47), 0x106AA070 );
    P( A, B, C, D, E, F, G, H, R(48), 0x19A4C116 );
    P( H, A, B, C, D, E, F, G, R(49), 0x1E376C08 );
    P( G, H, A, B, C, D, E, F, R(50), 0x2748774C );
    P( F, G, H, A, B, C, D, E, R(51), 0x34B0BCB5 );
    P( E, F, G, H, A, B, C, D, R(52), 0x391C0CB3 );
    P( D, E, F, G, H, A, B, C, R(53), 0x4ED8AA4A );
    P( C, D, E, F, G, H, A, B, R(54), 0x5B9CCA4F );
    P( B, C, D, E, F, G, H, A, R(55), 0x682E6FF3 );
    P( A, B, C, D, E, F, G, H, R(56), 0x748F82EE );
    P( H, A, B, C, D, E, F, G, R(57), 0x78A5636F );
    P( G, H, A, B, C, D, E, F, R(58), 0x84C87814 );
    P( F, G, H, A, B, C, D, E, R(59), 0x8CC70208 );
    P( E, F, G, H, A, B, C, D, R(60), 0x90BEFFFA );
    P( D, E, F, G, H, A, B, C, R(61), 0xA4506CEB );
    P( C, D, E, F, G, H, A, B, R(62), 0xBEF9A3F7 );
    P( B, C, D, E, F, G, H, A, R(63), 0xC67178F2 );

    D_REF(ctx->state,0,tid) += A;
    D_REF(ctx->state,1,tid) += B;
    D_REF(ctx->state,2,tid) += C;
    D_REF(ctx->state,3,tid) += D;
    D_REF(ctx->state,4,tid) += E;
    D_REF(ctx->state,5,tid) += F;
    D_REF(ctx->state,6,tid) += G;
    D_REF(ctx->state,7,tid) += H; 
}

void memcpy(uint32_t tid, __global uint32_t *restrict dst, uint32_t dofst, __global uint32_t *restrict src, uint32_t sofst, uint32_t length){
	uint32_t i;
	uint32_t save  = D_REF(dst,dofst/4,tid) & (uint32_t)(0xFFFFFFFFUL >> ((4-(dofst%4)) * 8) );
	uint32_t ssave = (uint64_t)D_REF(src,sofst/4,tid) >> ((sofst%4) * 8);

//	printf("save: %8.8X %8.8X\n", save, ssave);

	for(i=0; i < length/4 + (length%4?1:0); i++){
		uint32_t sword = D_REF(src,1+i+sofst/4,tid);
		uint32_t tword = ((uint64_t)sword << ((4 - (sofst%4))*8)) | ssave;
		D_REF(dst,i+dofst/4,tid) = save | ((uint64_t)tword << ((dofst % 4) * 8));
		save = ((uint64_t)tword >> ((4 - (dofst % 4)) * 8)) ;
		ssave = (uint64_t)sword >> ((sofst%4)*8);
//		printf("tick: %8.8X %8.8X %8.8X %8.8X\n",sword, tword, save,ssave);
	}
}

void sha256_update( uint32_t tid, __global sha256cl_context *restrict ctx, __global uint32_t *restrict input, uint32_t length )
{
    uint32_t left, fill, ofst;

    if( !length ) return;

    ofst=0;

    left = D_REF(ctx->total,0,tid) & 0x3F;
    fill = 64 - left;

    D_REF(ctx->total,0,tid) += length;
    D_REF(ctx->total,0,tid) &= 0xFFFFFFFF;

    if( D_REF(ctx->total,0,tid) < length )
        D_REF(ctx->total,1,tid)++;

    if( left && length >= fill )
    {
        memcpy(tid, ctx->buffer, left, input, ofst, fill );
        sha256_process( tid, ctx, ctx->buffer, 0 );
        length -= fill;
        ofst  += fill;
        left = 0;
    }

    while( length >= 64 )
    {
        sha256_process( tid, ctx, input, ofst );
        length -= 64;
        ofst  += 64;
    }

    if( length )
    {
        memcpy( tid, ctx->buffer,left,input, ofst, length );
    }
}


void sha256_finish(uint32_t tid, __global sha256cl_context *restrict ctx, __global uint32_t *restrict digest )
{
    uint32_t last, padn;
    uint32_t high, low;

    high = ( D_REF(ctx->total,0,tid) >> 29 )
         | ( D_REF(ctx->total,1,tid) <<  3 );
    low  = ( D_REF(ctx->total,0,tid) <<  3 );

    PUT_UINT32( high, ctx->msglen, 0 ,tid);
    PUT_UINT32( low,  ctx->msglen, 4 ,tid);

    last = D_REF(ctx->total,0,tid) & 0x3F;
    padn = ( last < 56 ) ? ( 56 - last ) : ( 120 - last );

    sha256_update( tid, ctx, ctx->padding, padn );
    sha256_update( tid, ctx, ctx->msglen, 8 );

    PUT_UINT32( D_REF(ctx->state,0,tid), digest,  0 , tid);
    PUT_UINT32( D_REF(ctx->state,1,tid), digest,  4 , tid);
    PUT_UINT32( D_REF(ctx->state,2,tid), digest,  8 , tid);
    PUT_UINT32( D_REF(ctx->state,3,tid), digest, 12 , tid);
    PUT_UINT32( D_REF(ctx->state,4,tid), digest, 16 , tid);
    PUT_UINT32( D_REF(ctx->state,5,tid), digest, 20 , tid);
    PUT_UINT32( D_REF(ctx->state,6,tid), digest, 24 , tid);
    PUT_UINT32( D_REF(ctx->state,7,tid), digest, 28 , tid); 
}
