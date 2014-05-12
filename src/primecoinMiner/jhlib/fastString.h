

// fast string general

#define FSTR_FORMAT_ASCII	1
#define FSTR_FORMAT_UTF8	2

typedef struct _fStr_t fStr_t;

typedef struct _fStr_format_t
{
	void (*fstr_copyASCII)(fStr_t* dest, char *str);
	void (*fStr_appendASCII)(fStr_t* dest, char *str);
	void (*fstr_copy)(fStr_t* dest, fStr_t* src);
	void (*fStr_append)(fStr_t* dest, fStr_t* src);
}fStr_format_t;

typedef struct _fStr_t 
{
	uint8*				str;
	fStr_format_t*		format; /* pointer to table of functions */
	uint32				length; /* in bytes */
	uint32				limit; /* in bytes */
	bool				allocated; /* if true, the string was allocated by the fastString lib */
}fStr_t;

fStr_t* fStr_alloc(uint32 bufferSize, uint32 format);
fStr_t* fStr_alloc(uint32 bufferSize);

void fStr_free(fStr_t* fStr);
char* fStr_get(fStr_t* fStr);
int fStr_len(fStr_t* fStr);

void fStr_reset(fStr_t* fStr);

// direct access
uint32 fStr_getLimit(fStr_t* fStr);
void fStr_setLength(fStr_t* fStr, uint32 length);

// simple
void fStr_copy(fStr_t* fStr, char *sourceASCII);
void fStr_append(fStr_t* fStr, char *sourceASCII);

void fStr_copy(fStr_t* fStr, fStr_t* source);
void fStr_append(fStr_t* fStr, fStr_t* source);

// extended
int fStr_appendFormatted(fStr_t* fStr, char *format, ...);

// non-instanced
char* fStrDup(char *src);
char* fStrDup(char *src, sint32 length);

void fStrCpy(char *dst, char *src, unsigned int limit);
int fStrLen(char *src);

// other
#ifdef _WIN64
void __cdecl _esprintf(char *out, char *format, uint64 *param, unsigned int *lengthOut);
#elif _WIN32
void __cdecl _esprintf(char *out, char *format, unsigned int *lengthOut);
#else
void  _esprintf(char *out, char *format, unsigned int *lengthOut);
#endif
//void _esprintf(char *out, char *format, unsigned int *param, unsigned int *lengthOut); // used only internally
void esprintf(char *out, char *format, ...);


int esprintf_X(char *out, unsigned int value, int padRight, int padZero, int width, int UpperCase);
int esprintf_u(char *out, unsigned int value, int padRight, int padZero, int width);
int esprintf_b(char *out, signed long long value, int padRight, int padZero, int width); //"Big" - signed long long
int esprintf_B(char *out, bool value, int padRight, int padZero, int width); //Boolean
int esprintf_c(char *out, char value, int padRight, int padZero, int width);
int esprintf_hf(char *out, float valueF, int padRight, int padZero, int width);
int esprintf_d(char *out, int value, int padRight, int padZero, int width);
int esprintf_xutf8(char *out, char *s, int padRight, int padZero, int width);
int esprintf_utf8(char *out, char *s, int padRight, int padZero, int width);
int esprintf_s(char *out, char *s, int padRight, int padZero, int width);











// tokenizer
char** fStrTokenize(char* src, char* tokens);
void fStrTokenizeClean(char** values);

// hashing
uint32 fStrGenHashA(char* str);
uint32 fStrGenHashB(char* str);

// case
void fStrConvertToLowercase(char* str);
sint32 fStrCmpCaseInsensitive(uint8* str1, uint8* str2, uint32 length);

// misc
void fStr_addHexString(fStr_t* fStr, uint8* data, uint32 dataLength);

// fast string (stack) buffers
typedef struct  
{
	fStr_t fStrObject;
	uint8 bufferData[128];
}fStr_buffer128b_t;

typedef struct  
{
	fStr_t fStrObject;
	uint8 bufferData[256];
}fStr_buffer256b_t;

typedef struct  
{
	fStr_t fStrObject;
	uint8 bufferData[512];
}fStr_buffer512b_t;

typedef struct  
{
	fStr_t fStrObject;
	uint8 bufferData[1024];
}fStr_buffer1kb_t;

typedef struct  
{
	fStr_t fStrObject;
	uint8 bufferData[4096];
}fStr_buffer4kb_t;

fStr_t* fStr_alloc(fStr_buffer128b_t* fStrBuffer, uint32 format=FSTR_FORMAT_UTF8);
fStr_t* fStr_alloc(fStr_buffer256b_t* fStrBuffer, uint32 format=FSTR_FORMAT_UTF8);
fStr_t* fStr_alloc(fStr_buffer1kb_t* fStrBuffer, uint32 format=FSTR_FORMAT_UTF8);
fStr_t* fStr_alloc(fStr_buffer4kb_t* fStrBuffer, uint32 format=FSTR_FORMAT_UTF8);



