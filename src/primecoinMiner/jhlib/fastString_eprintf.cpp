#include"./JHLib.h"
#include <cstdarg>
#include <iostream>

// #include"fastString.h"

static char esprintf_HEX_LOWER[16] = {'0','1','2','3','4','5','6','7','8','9','a','b','c','d','e','f'};
static char esprintf_HEX_UPPER[16] = {'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F'};


int esprintf_s(char *out, char *s, int padRight, int padZero, int width)
{
	int c = 0;
	while(*s)
	{
		*out = *s;
		c++;
		out++;
		s++;
	}
	return c;
}

int esprintf_utf8(char *out, char *s, int padRight, int padZero, int width)
{
	int c = 0;
	int len = fStrLen(s);
	while(len)
	{
		len--;
		*out = *s;
		c++;
		out++;
		s++;
	}
	return c;
}

int esprintf_xutf8(char *out, char *s, int padRight, int padZero, int width)
{
	int c = 0;
	int len = fStrLen(s);
	while(len)
	{
		len--;
		unsigned int v = ((*s)>>4)&0xF;
		*out = esprintf_HEX_UPPER[v]; out++;
		v = ((*s))&0xF;
		*out = esprintf_HEX_UPPER[v]; out++;
		c += 2;
		s++;
	}
	return c;
}

int esprintf_d(char *out, int value, int padRight, int padZero, int width)
{
	//int c = 0; unused
	char DezimalStr[32];
	int dl = 0; //DezimalLength
	int negative = 0;
	if( value < 0 )
	{
		negative = 1;
		value *= -1;
	}
	if( value == 0 )
	{
		DezimalStr[dl] = '0';
		dl++;		
	}
	else
	{
		while(value)
		{
			DezimalStr[dl] = '0'+(value%10);
			dl++;			
			value /= 10;
		}
	}

	if( width >= 0 )
	{
		int r = 0;
		int totalDigits = negative + dl;
		if( width < totalDigits )
			totalDigits = width;
		if( dl >= totalDigits )
			negative = 0;
		int NumberDigits = totalDigits - negative;
		int PadDigits = width - totalDigits;
		char PadDigit = ' ';
		if( padZero )
			PadDigit = '0';
		//Negative sign first
		if( negative )
		{
			out[r] = '-'; r++;
		}
		//Pad chars
		for(int i=0; i<PadDigits; i++)
		{
			out[r] = PadDigit; r++;
		}
		//The number
		for(int i=0; i<NumberDigits; i++)
		{
			out[r] = DezimalStr[NumberDigits-1-i]; r++;
		}
		return r;
	}
	else
	{
		//Copy dl to output
		if( negative )
		{
			*out = '-';
			out++;		
		}
		for(int i=0; i<dl; i++)
		{
			out[i] = DezimalStr[dl-1-i];
		}
	}
	return dl+negative;
}

int esprintf_hf(char *out, float valueF, int padRight, int padZero, int width)
{
	char* outInitial = out;
	sint32 value = (sint32)floor(valueF);
//	int c = 0;   unused?
	char DezimalStr[32];
	int dl = 0; //DezimalLength
	int negative = 0;
	if( value < 0 )
	{
		negative = 1;
		value *= -1;
	}
	if( value == 0 )
	{
		DezimalStr[dl] = '0';
		dl++;		
	}
	else
	{
		while(value)
		{
			DezimalStr[dl] = '0'+(value%10);
			dl++;			
			value /= 10;
		}
	}

	if( width >= 0 )
	{
		int r = 0;
		int totalDigits = negative + dl;
		if( width < totalDigits )
			totalDigits = width;
		if( dl >= totalDigits )
			negative = 0;
		int NumberDigits = totalDigits - negative;
		int PadDigits = width - totalDigits;
		char PadDigit = ' ';
		if( padZero )
			PadDigit = '0';
		//Negative sign first
		if( negative )
		{
			out[r] = '-'; r++;
		}
		//Pad chars
		for(int i=0; i<PadDigits; i++)
		{
			out[r] = PadDigit; r++;
		}
		//The number
		for(int i=0; i<NumberDigits; i++)
		{
			out[r] = DezimalStr[NumberDigits-1-i]; r++;
		}
		return r;
	}
	else
	{
		//Copy dl to output
		if( negative )
		{
			*out = '-';
			out++;		
		}
		for(int i=0; i<dl; i++)
		{
			*out = DezimalStr[dl-1-i];
			out++;
		}
	}
	// coma
	*out = '.';
	out++;		

	//int dl2 = dl; unused

	value = (sint32)((valueF - floor(valueF))*100000.0f);
	//c = 0; unused
	dl = 0; //DezimalLength
	padZero = 1;
	width = 6;
	//padRight = 0; unused
	negative = 0;
	if( value < 0 )
	{
		negative = 1;
		value *= -1;
	}
	if( value == 0 )
	{
		DezimalStr[dl] = '0';
		dl++;		
	}
	else
	{
		while(value)
		{
			DezimalStr[dl] = '0'+(value%10);
			dl++;			
			value /= 10;
		}
	}

	if( width >= 0 )
	{
		int r = 0;
		int totalDigits = negative + dl;
		if( width < totalDigits )
			totalDigits = width;
		if( dl >= totalDigits )
			negative = 0;
		int NumberDigits = totalDigits - negative;
		int PadDigits = width - totalDigits;
		char PadDigit = ' ';
		if( padZero )
			PadDigit = '0';
		//Negative sign first
		if( negative )
		{
			out[r] = '-'; r++;
		}
		//Pad chars
		for(int i=0; i<PadDigits; i++)
		{
			out[r] = PadDigit; r++;
		}
		//The number
		for(int i=0; i<NumberDigits; i++)
		{
			out[r] = DezimalStr[NumberDigits-1-i]; r++;
		}
		out += r;
	}
	else
	{
		//Copy dl to output
		if( negative )
		{
			*out = '-';
			out++;		
		}
		for(int i=0; i<dl; i++)
		{
			*out = DezimalStr[dl-1-i];
			out++;
		}
	}



	return out-outInitial;//dl2+1+dl+((valueF<0.0f)?-1:0);
}


int esprintf_c(char *out, char value, int padRight, int padZero, int width)
{
	// todo: support for padRight etc.
	//int c = 0; unused
	out[0] = value;
	return 1;
}


int esprintf_B(char *out, bool value, int padRight, int padZero, int width) //Boolean
{
	if(value){
		*out = '1';
	}else{
		*out = '0';
	}
	return 1;
}

int esprintf_b(char *out, signed long long value, int padRight, int padZero, int width) //"Big" - signed long long
{
	//int c = 0; unused
	char DezimalStr[32];
	int dl = 0; //DezimalLength
	int negative = 0;
	if( value < 0 )
	{
		negative = 1;
		value *= -1;
	}
	if( value == 0 )
	{
		DezimalStr[dl] = '0';
		dl++;		
	}
	else
	{
		while(value)
		{
			DezimalStr[dl] = '0'+(value%10);
			dl++;			
			value /= 10;
		}
	}

	if( width >= 0 )
	{
		int r = 0;
		int totalDigits = negative + dl;
		totalDigits = std::min(totalDigits, width);
		if( dl >= totalDigits )
			negative = 0;
		int NumberDigits = totalDigits - negative;
		int PadDigits = width - totalDigits;
		char PadDigit = ' ';
		if( padZero )
			PadDigit = '0';
		//Negative sign first
		if( negative )
		{
			out[r] = '-'; r++;
		}
		//Pad chars
		for(int i=0; i<PadDigits; i++)
		{
			out[r] = PadDigit; r++;
		}
		//The number
		for(int i=0; i<NumberDigits; i++)
		{
			out[r] = DezimalStr[NumberDigits-1-i]; r++;
		}
		return r;
	}
	else
	{
		//Copy dl to output
		if( negative )
		{
			*out = '-';
			out++;		
		}
		for(int i=0; i<dl; i++)
		{
			out[i] = DezimalStr[dl-1-i];
		}
	}
	return dl+negative;
}


int esprintf_u(char *out, unsigned int value, int padRight, int padZero, int width)
{
	//int c = 0; unused
	char DezimalStr[32];
	int dl = 0; //DezimalLength
	if( value == 0 )
	{
		DezimalStr[dl] = '0';
		dl++;		
	}
	else
	{
		while(value)
		{
			DezimalStr[dl] = '0'+(value%10);
			dl++;			
			value /= 10;
		}
	}

	if( width >= 0 )
	{
		int r = 0;
		int totalDigits = dl;
		totalDigits = std::min(totalDigits, width);
		int NumberDigits = totalDigits;
		int PadDigits = width - totalDigits;
		char PadDigit = ' ';
		if( padZero )
			PadDigit = '0';
		//Pad chars
		for(int i=0; i<PadDigits; i++)
		{
			out[r] = PadDigit; r++;
		}
		//The number
		for(int i=0; i<NumberDigits; i++)
		{
			out[r] = DezimalStr[NumberDigits-1-i]; r++;
		}
		return r;
	}
	else
	{
		//Copy dl to output
		for(int i=0; i<dl; i++)
		{
			out[i] = DezimalStr[dl-1-i];
		}
	}
	return dl;
}

int esprintf_X(char *out, unsigned int value, int padRight, int padZero, int width, int UpperCase)
{
	//int c = 0; unused
	char DezimalStr[32];
	int dl = 0; //HexLength
	while(value)
	{
		if( UpperCase )
			DezimalStr[dl] = esprintf_HEX_UPPER[value&0xF];
		else
			DezimalStr[dl] = esprintf_HEX_LOWER[value&0xF];
		dl++;			
		value >>= 4;
	}

	if( width >= 0 )
	{
		int r = 0;
		int totalDigits = dl;
		totalDigits = std::min(totalDigits, width);
		int NumberDigits = totalDigits;
		int PadDigits = width - totalDigits;
		char PadDigit = ' ';
		if( padZero )
			PadDigit = '0';
		//Pad chars
		for(int i=0; i<PadDigits; i++)
		{
			out[r] = PadDigit; r++;
		}
		//The number
		for(int i=0; i<NumberDigits; i++)
		{
			out[r] = DezimalStr[NumberDigits-1-i]; r++;
		}
		return r;
	}
	else
	{
		//Copy dl to output
		for(int i=0; i<dl; i++)
		{
			out[i] = DezimalStr[dl-1-i];
		}
	}
	return dl;
}

/*
 * Param = pointer to parameters used for format insertion
 */
#ifdef _WIN64
void __cdecl _esprintf(char *out, char *format, unsigned int *lengthOut)
#elif defined (_WIN32)
void __cdecl _esprintf(char *out, char *format, unsigned int *lengthOut)
#else // gcc
void __attribute__((cdecl)) _esprintf(char *out, char *format, unsigned int *lengthOut) 
#endif
{

}

void esprintf(char *out, char *format, ...)
{
	// use some dirty trick to access varying arguments
	//unsigned int *param = (unsigned int*)_ADDRESSOF(format);
	//param++; // skip first parameter
	//_esprintf(out, format, param, NULL);
	#ifdef _WIN64
		uint64 *param = (uint64*)_ADDRESSOF(format);
		param++; // skip first parameter
		unsigned int formattedLength = 0;
		_esprintf(out, format, param, &formattedLength);
	#else
	va_list arguments;
	va_start ( arguments, format );           // Initializing arguments to store all values after format

	unsigned int *lengthOut = 0;
	//Do parsing
	char *p = format;
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
					o += esprintf_c(o, va_arg(arguments,char), PadRight, PadZero, Width);
				}
				else if( *p == 'b' ) //signed long integer
				{
					o += esprintf_b(o, va_arg(arguments,signed long long), PadRight, PadZero, Width);
				}
				else if( *p == 'B' ) //boolean
				{
					o += esprintf_B(o, va_arg(arguments,bool), PadRight, PadZero, Width);
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
					o += esprintf_hf(o, va_arg(arguments,float), PadRight, PadZero, Width);
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
	if( lengthOut )
		*lengthOut = (unsigned int)(o-out);
	va_end(arguments);

#endif

}

