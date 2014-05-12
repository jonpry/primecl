

#ifndef __JHSYSTEMLIB
#define __JHSYSTEMLIB

#ifdef _WIN32
#define NOMINMAX
#include<Windows.h>
#else
#include <signal.h>
#include <stdint.h>
#endif
#include <cstring> // for memcpy/memset
#include<math.h>
#include <algorithm>


#ifndef _WIN32 // temporary replacement for _ADDRESSOF, replace with boost
template< class T >
T* tmp_addressof(T& arg) {
    return reinterpret_cast<T*>(
               &const_cast<char&>(
                  reinterpret_cast<const volatile char&>(arg)));
}
#endif
/*
typedef unsigned long long 	uint64;
typedef signed long long	sint64;

typedef unsigned int 	uint32;
typedef signed int 		sint32;

typedef unsigned short 	uint16;
typedef signed short 	sint16;

typedef unsigned char 	uint8;
typedef signed char 	sint8;*/

#ifdef _WIN32
typedef __int64           sint64;
typedef unsigned __int64  uint64;
typedef __int32           sint32;
typedef unsigned __int32  uint32;
typedef __int16           sint16;
typedef unsigned __int16  uint16;
typedef __int8            sint8;
typedef unsigned __int8   uint8;
#else
typedef int64_t sint64;
typedef uint64_t uint64;
typedef int32_t sint32;
typedef uint32_t uint32;
typedef int16_t sint16;
typedef uint16_t uint16;
typedef int8_t sint8t;
typedef uint8_t uint8;
#endif


#define JHCALLBACK	__fastcall


void* _ex1_malloc(int size);
void* _ex1_realloc(void* old, int size);
void _ex1_free(void* p);

void _ex2_initialize();

void* _ex2_malloc(int size, char* file, sint32 line);
void* _ex2_realloc(void* old, int size, char* file, sint32 line);
void _ex2_free(void* p, char* file, sint32 line);
void _ex2_analyzeMemoryLog();

// memory validator
//#define malloc(x) _ex1_malloc(x)
//#define realloc(x,y) _ex1_realloc(x,y)
//#define free(x) _ex1_free(x)

// memory logger
//#define MEMORY_LOGGER_ACTIVE			1
//#define MEMORY_LOGGER_ANALYZE_ACTIVE	1

#ifdef MEMORY_LOGGER_ACTIVE
#define malloc(x) _ex2_malloc(x,__FILE__,__LINE__)
#define realloc(x,y) _ex2_realloc(x,y,__FILE__,__LINE__)
#define free(x) _ex2_free(x,__FILE__,__LINE__)
#endif

/*#include".\streamWrapper.h"
#include".\fastString.h"
#include".\hashTable.h"
#include".\fastSorter.h"
#include".\fileMgr.h"
#include".\sData.h"
#include".\bmp.h"
#include".\tgaLib.h"
#include".\fMath.h"
#include".\packetBuffer.h"
#include".\msgQueue.h"
#include".\simpleList.h"
#include".\customBuffer.h"*/

#include"streamWrapper.h"
#include"fastString.h"
#include"hashTable.h"
#ifdef _WIN32
#include"fMath.h"
#include"sData.h"
#endif
#include"fastSorter.h"
#include"fileMgr.h"
#include"bmp.h"
#include"tgaLib.h"
#include"msgQueue.h"
#include"packetBuffer.h"
#include"simpleList.h"
#include"customBuffer.h"


/* error */
#define assertFatal(condition, message, errorCode) if( condition ) { OutputDebugString(message);  ExitProcess(errorCode); } 


#endif

