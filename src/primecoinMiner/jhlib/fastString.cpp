#include"./JHLib.h"
#include<stdio.h>
#include <cstdarg>
#include <iostream>

#ifndef _WIN32
#include <cstdlib>
#endif

// local variables
extern fStr_format_t fStr_formatInfo_ASCII;
extern fStr_format_t fStr_formatInfo_UTF8;

fStr_t* fStr_alloc(uint32 bufferSize, uint32 format)
{
	fStr_t *fStr = (fStr_t*)malloc(sizeof(fStr_t)+bufferSize);
	fStr->str = (uint8*)(fStr+1);
	fStr->allocated = true;
	fStr->length = 0;
	fStr->limit = bufferSize;
	if( format == FSTR_FORMAT_ASCII )
	{
		fStr->format = &fStr_formatInfo_ASCII;
	}
	else
	{
#ifdef _WIN32
		__debugbreak(); // unknown format
#else
	    raise(SIGTRAP);
#endif
	}
	return fStr;
}

fStr_t* fStr_alloc(uint32 bufferSize)
{
	return fStr_alloc(bufferSize, FSTR_FORMAT_ASCII);
}

void _fStr_allocForBuffer(fStr_t* fStr, uint8* buffer, uint32 bufferSize, uint32 format)
{
	fStr->str = buffer;
	fStr->allocated = false;
	fStr->length = 0;
	fStr->limit = bufferSize;
	if( format == FSTR_FORMAT_ASCII )
	{
		fStr->format = &fStr_formatInfo_ASCII;
	}
	else if( format == FSTR_FORMAT_UTF8 )
	{
		fStr->format = &fStr_formatInfo_UTF8;
	}
	else
#ifdef _WIN32
		__debugbreak(); // unknown format
#else
	    raise(SIGTRAP);
#endif
}

/*
 * Fast allocate a stack fStr object + an assigned buffer of 128 bytes
 */
fStr_t* fStr_alloc(fStr_buffer128b_t* fStrBuffer, uint32 format)
{
	_fStr_allocForBuffer(&fStrBuffer->fStrObject, fStrBuffer->bufferData, sizeof(fStrBuffer->bufferData), format);
	return &fStrBuffer->fStrObject;
}

/*
 * Fast allocate a stack fStr object + an assigned buffer of 256 bytes
 */
fStr_t* fStr_alloc(fStr_buffer256b_t* fStrBuffer, uint32 format)
{
	_fStr_allocForBuffer(&fStrBuffer->fStrObject, fStrBuffer->bufferData, sizeof(fStrBuffer->bufferData), format);
	return &fStrBuffer->fStrObject;
}

/*
 * Fast allocate a stack fStr object + an assigned buffer of 1024 bytes
 */
fStr_t* fStr_alloc(fStr_buffer1kb_t* fStrBuffer, uint32 format)
{
	_fStr_allocForBuffer(&fStrBuffer->fStrObject, fStrBuffer->bufferData, sizeof(fStrBuffer->bufferData), format);
	return &fStrBuffer->fStrObject;
}

/*
 * Fast allocate a stack fStr object + an assigned buffer of 4096 bytes
 */
fStr_t* fStr_alloc(fStr_buffer4kb_t* fStrBuffer, uint32 format)
{
	_fStr_allocForBuffer(&fStrBuffer->fStrObject, fStrBuffer->bufferData, sizeof(fStrBuffer->bufferData), format);
	return &fStrBuffer->fStrObject;
}


void fStr_free(fStr_t* fStr)
{
	free(fStr);
}

void fStr_reset(fStr_t* fStr)
{
	fStr->length = 0;
}

char* fStr_get(fStr_t* fStr)
{
	fStr->str[fStr->length] = '\0';
	return (char*)fStr->str;
}

uint32 fStr_getLimit(fStr_t* fStr)
{
	return fStr->limit;
}

void fStr_setLength(fStr_t* fStr, uint32 length)
{
	fStr->length = length;
}

int fStr_len(fStr_t* fStr)
{
	return fStr->length;
}

/* default string operations */

void fStr_copy(fStr_t* fStr, char *sourceASCII)
{
	fStr->format->fstr_copyASCII(fStr, sourceASCII);
	if( fStr->length >= fStr->limit )
	{
#ifdef _WIN32
		OutputDebugString("fStr: Bufferoverflow detected");
#endif
		exit(-32002);
	}
}

void fStr_append(fStr_t* fStr, char *sourceASCII)
{
	fStr->format->fStr_appendASCII(fStr, sourceASCII);
	if( fStr->length >= fStr->limit )
	{
#ifdef _WIN32
		OutputDebugString("fStr: Bufferoverflow detected");
#endif
		exit(-32002);
	}
}

void fStr_copy(fStr_t* fStr, fStr_t* source)
{
	fStr->format->fstr_copy(fStr, source);
	if( fStr->length >= fStr->limit )
	{
#ifdef _WIN32
		OutputDebugString("fStr: Bufferoverflow detected");
#endif
		exit(-32001);
	}
}

void fStr_append(fStr_t* fStr, fStr_t* source)
{
	fStr->format->fStr_append(fStr, source);
	if( fStr->length >= fStr->limit )
	{
#ifdef _WIN32
		OutputDebugString("fStr: Bufferoverflow detected");
#endif
		exit(-32001);
	}
}

/* extended string operations */

/*
 * Formats and appends a string
 */
int fStr_appendFormatted(fStr_t* fStr, char *format, ...)
{
	// use some dirty trick to access varying arguments
	unsigned int formattedLength = 0;

	va_list arguments;
	va_start ( arguments, format );           // Initializing arguments to store all values after format

	//Do parsing
	char *p = format;
	char *out = (char*)(fStr->str + fStr->length);
	char *o = out;
	while(*p)
	{
		if( *p == '%' )
		{
			p++;
			if( *p == '%' )
			{
				*o = *p;
				p++;
				o++;
			}
			else
			{
				//Parse format
				//%[-][#][0][width][.precision]type
				int PadRight = 0;
				int AutoPrefix = 0;
				int PadZero = 0;
				int Width = -1;
				int Precision = 1337;
				//PadRight
				if( *p == '-' )
				{
					PadRight = 1;
					p++;
				}
				//AutoPrefix
				if( *p == '#' )
				{
					AutoPrefix = 1;
					p++;
				}
				//ZeroPad
				if( *p == '0' )
				{
					PadZero = 1;
					p++;
				}
				//Width (2 digits at max)
				if( *p >= '0' && *p <= '9' )
				{
					Width = *p-'0';
					p++;
					if( *p >= '0' && *p <= '9' )
					{
						Width = Width*10 + (*p-'0');
						p++;
					}
				}
				//Precision
				if( *p == '.' )
				{
					p++;
					Precision = *p-'0';
					p++;
					if( *p >= '0' && *p <= '9' )
					{
						Precision = Precision*10 + (*p-'0');
						p++;
					}
				}
				//Parse type
				int LongMode = 0;
				if( *p == 'l' )
				{
					LongMode = 1;
					p++;
				}
				//Now check case

				if( *p == 'd' ) //signed integer
				{
					o += esprintf_d(o, va_arg(arguments,sint64), PadRight, PadZero, Width);
				}
				else if( p[0] == 'u' && p[1] == 't' && p[2] == 'f' && p[3] == '8' ) //utf8 string
				{
					
					o += esprintf_utf8(o, va_arg(arguments,char*), PadRight, PadZero, Width);
					p += 3;
				}
				else if( p[0] == 'x' && p[1] == 'u' && p[2] == 't' && p[3] == 'f' && p[4] == '8' ) //hex-encoded utf8 string
				{
					o += esprintf_xutf8(o, va_arg(arguments,char*), PadRight, PadZero, Width);
					p += 4;
				}
				else if( *p == 'u' ) //signed integer
				{
					o += esprintf_u(o, va_arg(arguments,unsigned int), PadRight, PadZero, Width);
				}
				else if( *p == 'c' ) //signed ascii char
				{
					o += esprintf_c(o, va_arg(arguments,int), PadRight, PadZero, Width);
				}
				else if( *p == 'b' ) //signed long integer
				{
					o += esprintf_b(o, va_arg(arguments,signed long long), PadRight, PadZero, Width);
				}
				else if( *p == 'B' ) //boolean
				{
					o += esprintf_B(o, va_arg(arguments,int), PadRight, PadZero, Width);
				}
				else if( *p == 's' ) //string
				{
					o += esprintf_s(o, va_arg(arguments,char*), PadRight, PadZero, Width);
				}
				else if( *p == 'X' ) //unsigned integer as hex
				{
					o += esprintf_X(o, va_arg(arguments,unsigned int), PadRight, PadZero, Width, 1);
				}
				else if( *p == 'x' ) //unsigned integer as hex
				{
					o += esprintf_X(o, va_arg(arguments,unsigned int), PadRight, PadZero, Width, 0);
				}
				else if( p[0] == 'h' && p[1] == 'f' ) // helper float
				{
					o += esprintf_hf(o, va_arg(arguments,double), PadRight, PadZero, Width);
					p += 1;
				}
				p++;

			}

		}
		else
		{
			*o = *p;
			p++;
			o++;
		}
	}
	*o = '\0';
	va_end(arguments);
	formattedLength = (unsigned int)(o-out);
	fStr->length += formattedLength;
	return formattedLength;

}

/* format specific - ASCII */

void fstr_ASCII_copyASCII(fStr_t* fStr, char *str)
{
	uint32 c = 0;
	while( *str )
	{
		fStr->str[c] = *str;
		str++;
		c++;
	}
	fStr->length = c;
}

void fstr_ASCII_appendASCII(fStr_t* fStr, char *source)
{
	uint32 c = fStr->length;
	while( *source )
	{
		fStr->str[c] = *source;
		source++;
		c++;
	}
	fStr->length = c;
}

void fstr_ASCII_copy(fStr_t* dest, fStr_t* src)
{
	for(uint32 c=0; c<src->length; c++)
	{
		dest->str[c] = src->str[c];
	}
	dest->length = src->length;
}

void fstr_ASCII_append(fStr_t* dest, fStr_t* src)
{
	for(uint32 c=0; c<src->length; c++)
	{
		dest->str[c+dest->length] = src->str[c];
	}
	dest->length += src->length;
}


fStr_format_t fStr_formatInfo_ASCII =
{
	fstr_ASCII_copyASCII,
	fstr_ASCII_appendASCII,
	fstr_ASCII_copy,
	fstr_ASCII_append,
};

/* format specific - UTF8 */

void fstr_UTF8_copyASCII(fStr_t* fStr, char *str)
{
	uint32 c = 0;
	while( *str )
	{
		fStr->str[c] = *str;
		str++;
		c++;
	}
	fStr->length = c;
}

void fstr_UTF8_appendASCII(fStr_t* fStr, char *source)
{
	uint32 c = fStr->length;
	while( *source )
	{
		fStr->str[c] = *source;
		source++;
		c++;
	}
	fStr->length = c;
}

void fstr_UTF8_copy(fStr_t* dest, fStr_t* src)
{
	for(uint32 c=0; c<src->length; c++)
	{
		dest->str[c] = src->str[c];
	}
	dest->length = src->length;
}

void fstr_UTF8_append(fStr_t* dest, fStr_t* src)
{
	for(uint32 c=0; c<src->length; c++)
	{
		dest->str[c+dest->length] = src->str[c];
	}
	dest->length += src->length;
}


fStr_format_t fStr_formatInfo_UTF8 =
{
	fstr_UTF8_copyASCII,
	fstr_UTF8_appendASCII,
	fstr_UTF8_copy,
	fstr_UTF8_append,
};

/* general (non-safe) string operations */

char* fStrDup(char *src)
{
#ifdef _WIN32
	return _strdup(src);
#else
	return strdup(src);
#endif
}

char* fStrDup(char *src, sint32 length)
{
	char* ns = (char*)malloc(length+1);
	memcpy(ns, src, length);
	ns[length] = '\0';
	return ns;
}


void fStrCpy(char *dst, char *src, unsigned int limit)
{
	limit--;
	while( *src && limit )
	{
		*dst = *src;
		src++; dst++;
		limit--;
	}
	*dst = '\0';
}

int fStrLen(char *src)
{
	int c = 0;
	while( *src )
	{
		src++;
		c++;
	}
	return c;
}

char** fStrTokenize(char* src, char* tokens)
{
	// create token lookup table
	bool lookup[256];
	memset(lookup, 0, sizeof(lookup));
	while( *tokens )
	{
		lookup[(uint8)(*tokens)] = true;
		tokens++;
	}
	// tokenize string
	char** params = (char**)malloc(sizeof(char*) * 64);
	memset(params, 0, sizeof(char*) * 64);
	uint32 indexStart = 0;
	uint32 paramCounter = 0;
	for(uint32 i=0; i<0x7FFFFFFF; i++)
	{
		if( lookup[src[i]] || src[i] == '\0' )
		{
			// split param
			uint32 pLen = i-indexStart;
			if( pLen > 0 )
			{
				params[paramCounter] = fStrDup(src+indexStart, pLen);
				paramCounter++;
			}
			indexStart = i+1;
			if( src[i] == '\0' || paramCounter >= 64 )
				break;
		}
	}
	if( paramCounter == 0 )
	{
		free(params);
		return NULL;
	}
	return params;
}

void fStrTokenizeClean(char** values)
{
	if( values == NULL )
		return;
	for(uint32 i=0; i<64; i++)
	{
		if( values[i] )
			free(values[i]);
	}
	free(values);
}

/*
 * Calculates the fast hash of an ASCII string (Method 1)
 */
uint32 fStrGenHashA(char* str)
{
	uint32 h = 0xef4cb8f0;
	while( *str )
	{
		uint32 c = (uint32)*str;
		h ^= c;
		h = (h<<3) | (h>>29);
		str++;
	}
	return h;
}

/*
 * Calculates the fast hash of an ASCII string (Method 2)
 */
uint32 fStrGenHashB(char* str)
{
	uint32 h1 = 0x77b192d3;
	uint32 h2 = 0xab9941bb;
	while( *str )
	{
		uint32 c = (uint32)*str;
		h1 ^= c;
		h2 += c;
		h1 = (h1>>3) | (h1<<29);
		h2 = (h2<<5) | (h2>>27);
		h1 += h2;
		str++;
	}
	return h1;
}

/*
 * Converts the string to lower case
 */
void fStrConvertToLowercase(char* str)
{
	while( *str )
	{
		uint8 c = *str;
		if( c >= 'A' && c <= 'Z' )
		{
			c -= ('A'-'a');
			*str = c;
		}
		str++;
	}

}

/*
 * Adds an array a lower case hex string (without a preleading '0x')
 */
void fStr_addHexString(fStr_t* fStr, uint8* data, uint32 dataLength)
{
	uint8* ptr = (fStr->str+fStr->length);
	for(uint32 i=0; i<dataLength; i++)
	{
		sprintf((char*)ptr, "%02x", data[i]);
		ptr += 2;
	}
	fStr->length += dataLength*2;
}


/*
 * Returns the difference of the un-equal character
 * The string size is expected to be at least length
 * The NT character is ignored
 */
sint32 fStrCmpCaseInsensitive(uint8* str1, uint8* str2, uint32 length)
{
	for(uint32 i=0; i<length; i++)
	{
		uint8 c1 = str1[i];
		uint8 c2 = str2[i];
		// convert both chars to lower case
		if( c1 >= 'A' && c1 <= 'Z' )
			c1 -= ('A'-'a');
		if( c2 >= 'A' && c2 <= 'Z' )
			c2 -= ('A'-'a');
		if( c1 != c2 )
			return c2-c1;
	}
	return 0;
}

