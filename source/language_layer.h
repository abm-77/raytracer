#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <memory.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>

/// NOTE(bryson): Context Cracking
#pragma region CONTEXT
#if defined (__clang__)
    #define COMPILER_CLANG 1
    #if defined (_WIN32)
        #define OS_WINDOWS 1
    #elif defined (__gnu_linux__)
        #define OS_LINUX 1
    #elif defined (__APPLE__) && defined (__MACH__)
        #define OS_MAC 1
    #else
        #error Missing OS detection
    #endif

    #if defined(__amd64__)
        #define ARCH_X64 1
    #elif defined(__i386__)
        #define ARCH_X86 1
    #elif defined(__arm__)
        #define ARCH_ARM 1
    #elif defined(__arch64__)
        #define ARCH_ARM64 1
    #else
        #error Missing ARCH detection
    #endif

#elif defined (_MSC_VER)
    #define COMPILER_CL 1
    #if defined (_WIN32)
        #define OS_WINDOWS 1
    #else
        #error Missing OS detection
    #endif

    #if defined(_M_AMD64)
        #define ARCH_X64 1
    #elif defined (_M_IX86)
        #define ARCH_X86 1
    #elif defined (_M_ARM)
        #define ARCH_ARM 1
    #else
        #error Missing ARCH detection
    #endif

#elif defined (__GNUC__)
    #define COMPILER_GCC 1
    #if defined (_WIN32)
        #define OS_WINDOWS 1
    #elif defined (__gnu_linux__)
        #define OS_LINUX 1
    #elif defined (__APPLE__) && defined (__MACH__)
        #define OS_MAC 1
    #else
        #error Missing OS detection
    #endif

    #if defined(__amd64__)
        #define ARCH_X64 1
    #elif defined(__i386__)
        #define ARCH_X86 1
    #elif defined(__arm__)
        #define ARCH_ARM 1
    #elif defined(__arch64__)
        #define ARCH_ARM64 1
    #else
        #error Missing ARCH detection
    #endif
#else
    #error Could not determine compiler!
#endif

#pragma endregion

// NOTE(bryson): Helper Macros
#pragma region HELPERS

#define _Stringify(S) #S
#define Stringify(S) _Stringify(S)
#define _Glue(A,B) A##B
#define Glue(A,B) _Glue(A,B)

#define ArrayLengthStatic(a) (sizeof((a)) / sizeof((a)[0]))

#define IntFromPtr(p) (unsigned long long) ((char*)p - (char*)0)
#define PtrFromInt(i) (void*)((char*)0 + i)

#define Member(T, m) (((T*)0)->m)
#define OffsetOfMember(T,m) IntFromPtr(&Member(T,m))

#define Min(a,b)(((a)<(b))?(a):(b))
#define Max(a,b)(((a)>(b))?(a):(b))
#define Clamp(x, l, u) (((x)<(l))?(l):((u)<(x))?(u):(x))
#define ClampTop(a,b) Min(a,b)
#define ClampBot(a,b) Max(a,b)

#define AlignUpPow2(x,p) ((x) + (p) - 1) & ~((p) - 1)
#define AlignDown(x,p) ((x) & ~((p) - 1)

#define global static
#define local_persist static
#define internal static  

#define CLinkageBegin extern "C" {
#define CLinkageEnd }
#define CLinkage extern "C"

#define MemoryMove memmove
#define MemorySet memset
#define MemoryCompare memcmp
#define MemoryReallocate realloc

#define MemoryZero(p,z) MemorySet((p), 0, (z))
#define MemoryZeroStruct(p) MemoryZero(p, sizeof(*(p)))
#define MemoryZeroArray(p) MemoryZero((p), sizeof(p))
#define MemoryZeroTyped(p,c) MemoryZero((p), sizeof(*(p)) * (c))

#define MemoryMatch(a,b,z) (MemoryCompare((a), (b), (z)) == 0)

#define MemoryCopy(d, s, z) MemoryMove((d), (s), (z))
#define MemoryCopyStruct(d, s) MemoryCopy((d), (s),\
                                                Min(sizeof(*(d)), sizeof(*(s))))
#define MemoryCopyArray(d, s) MemoryCopy((d), (s), Min(sizeof(s), sizeof(d))) 
#define MemoryCopyTyped(d, s, c) MemoryCopy((d), (s),\
                                                Min(sizeof(*(d)), sizeof(*(s)) * (c)))
#define ApproxLess(x, y) ((x) < (y + MACHINE_EPSILON_F32))
#define ApproxGreater(x, y) ((x) > (y - MACHINE_EPSILON_F32))
#define ApproxEqual(x, y) ApproxLess(x,y) && ApproxGreater(x,y)

#define SquareRoot sqrt
#define AbsoluteValue fabsf
#pragma endregion


// NOTE(bryson): Symbolic Constants
#pragma region SYMBOLIC_CONSTANTS
#define Bytes(n) ((u64)(n))
#define Kilobytes(n) (Bytes(n) * 1024)
#define Megabytes(n) (Kilobytes(n) * 1024)
#define Gigabytes(n) (Megabytes(n) * 1024)
#define Terabytes(n) (Gigabytes(n) * 1024) 

enum Side {
    SIDE_MIN,
    SIDE_MAX
};


enum OperatingSystem {
    OPERATINGSYSTEM_NULL,
    OPERATINGSYSTEM_WINDOWS,
    OPERATINGSYSTEM_MAC,
    OPERATINGSYSTEM_LINUX,
    OPERATINGSYSTEM_COUNT,
};

enum Architecture {
    ARCHITECTURE_NULL,
    ARCHITECTURE_X86,
    ARCHITECTURE_X64,
    ARCHITECTURE_ARM,
    ARCHITECTURE_ARM64,
    ARCHITECTURE_COUNT,
};

enum Month {
    MONTH_JAN,
    MONTH_FEB,
    MONTH_MAR,
    MONTH_APR,
    MONTH_MAY,
    MONTH_JUN,
    MONTH_JUL,
    MONTH_AUG,
    MONTH_SEP,
    MONTH_OCT,
    MONTH_NOV,
    MONTH_DEC,
};

enum Weekday {
    WEEKDAY_MON,
    WEEKDAY_TUE,
    WEEKDAY_WED,
    WEEKDAY_THU,
    WEEKDAY_FRI,
    WEEKDAY_SAT,
    WEEKDAY_SUN,
};
/*
internal OperatingSystem OperatingSystemFromContext(void);
internal Architecture ArchitectureFromContext(void);
internal const char* StringFromOperatingSystem(OperatingSystem os);
internal const char* StringFromArchitecture(Architecture arch);
internal const char* StringFromMonth(Month month);
internal const char* StringFromWeekday(Weekday day);
*/
#pragma endregion


// NOTE(bryson): Types
#pragma region TYPES
typedef int8_t   i8;
typedef int16_t  i16;
typedef int32_t  i32;
typedef int64_t  i64;
typedef uint8_t  u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;
typedef i8       s8;
typedef i16      s16;
typedef i32      s32;
typedef i64      s64;
typedef i8       b8;
typedef i16      b16;
typedef i32      b32;
typedef i64      b64;
typedef float    f32;
typedef double   f64;
typedef unsigned char byte;

typedef void VoidFunc(void);

// Basic Constants
global const i8 MIN_I8 = (i8) 0x80;
global const i16 MIN_I16 = (i16) 0x8000;
global const i32 MIN_I32 = (i32) 0x80000000;
global const i64 MIN_I64 = (i64) 0x8000000000000000llu;

global const i8 MAX_I8 = (i8) 0x7f;
global const i16 MAX_I16 = (i16) 0x7fff;
global const i32 MAX_I32 = (i32) 0x7fffffff;
global const i64 MAX_I64 = (i64) 0x7fffffffffffffffllu;

global const i8 MAX_U8 = (i8) 0xff;
global const i16 MAX_U16 = (i16) 0xffff;
global const i32 MAX_U32 = (i32) 0xffffffff;
global const i64 MAX_U64 = (i64) 0xffffffffffffffffllu;

global const f32 MACHINE_EPSILON_F32 = 1.1920929e-07f;
global const f32 PI_F32 = 3.14159265359f;
global const f32 TAU_F32 = 6.28318530718f;
global const f32 e_F32 = 2.71828182846f;
global const f32 GOLD_BIG_F32 = 1.61803398875f;
global const f32 GOLD_SMALL_F32 = 0.61803398875f;

global const f64 MACHINE_EPSILON_F64 = 1.1920929e-07;
global const f64 PI_F64 = 3.14159265359;
global const f64 TAU_F64 = 6.28318530718;
global const f32 e_F64 = 2.71828182846;
global const f64 GOLD_BIG_F64 = 1.61803398875;
global const f64 GOLD_SMALL_F64 = 0.61803398875;

global const f64 INF_F64 = std::numeric_limits<f64>::infinity();
#pragma endregion


// Debug
#define Statement(S) do{S}while(0)

#define BreakDebugger() _DebugBreakInternal()

#if ENABLE_ASSERT 
    #define Assert HardAssert
    #define HardAssert(c) Statement( if (!(c)) { _AssertFailure(#c, __FILE__, __LINE__, true); })
    #define SoftAssert(c) Statement( if (!(c)) { _AssertFailure(#c, __FILE__, __LINE__, false); })
#else
    #define Assert(c)
#endif

enum LogType {
    LOG_INFO,
    LOG_WARNING,
    LOG_ERROR
};

#define Log(...)            _DebugLog(LOG_INFO, __FILE__, __LINE__, __VA_ARGS__)
#define LogWarning(...)     _DebugLog(LOG_WARNING, __FILE__, __LINE__, __VA_ARGS__)
#define LogError(...)       _DebugLog(LOG_ERROR, __FILE__, __LINE__, __VA_ARGS__)

void _DebugLog(LogType log_type, const char* file, u32 line, const char* format, ...);
void _AssertFailure(const char* expression, const char* file, u32 line, b32 crash);
void _DebugBreakInternal(void);

// Note(bryson): Zero fill missing context macros
#pragma region DEFAULTS
#if !defined (COMPILER_CL) 
    #define COMPILER_CL 0
#endif
#if !defined (COMPILER_CLANG)
    #define COMPILER_CLANG 0 
#endif
#if !defined (COMPILER_GCC)
    #define COMPILER_GCC 0 
#endif
#if !defined (OS_WINDOWS)
    #define OS_WINDOWS 0
#endif
#if !defined (OS_LINUX)
    #define OS_LINUX 0 
#endif
#if !defined (OS_MAC)
    #define OS_MAC 0 
#endif
#if !defined (ARCH_X64)
    #define ARCH_X64 0
#endif
#if !defined (ARCH_X86)
    #define ARCH_X86 0
#endif
#if !defined (ARCH_ARM)
    #define ARCH_ARM 0
#endif
#if !defined (ARCH_ARM64)
    #define ARCH_ARM64 0
#endif
#if !defined (ARCH_ARM64)
    #define ARCH_ARM64 0
#endif
#if !defined (ENABLE_ASSERT)
    #define ENABLE_ASSERT 0
#endif

#pragma endregion

#define EvalPrint(x) printf("%s = %d\n", #x, (i32) (x))
#define EvalPrintLL(x) printf("%s = %lld\n", #x, (i64) (x))
#define EvalPrintU(x) printf("%s = %u\n", #x, (u32) (x))
#define EvalPrintULL(x) printf("%s = %llu\n", #x, (u64) (x))
#define EvalPrintF(x) printf("%s = %e\n", #x, (f64) (x))
#define EvalPrintV2F(v) printf("%s = \tx:%f, \ty:%f, \tz:%f\n", #v, v.x, v.y)
#define EvalPrintV3F(v) printf("%s = \tx:%f, \ty:%f, \tz:%f\n", #v, v.x, v.y, v.z)
#define EvalPrintV4F(v) printf("%s = \tx:%f, \ty:%f, \tz:%f\n", #v, v.x, v.y, v.z, v.w)
#define EvalPrintMat4(m) printf("%s = \t%f \t%f \t%f \t%f \n \t%f \t%f \t%f \t%f \n \t%f \t%f \t%f \t%f \n \t%f \t%f \t%f \t%f", #m, m.v[0], m.v[1], m.v[2], m.v[3], m.v[4], m.v[5], m.v[6], m.v[7], m.v[8], m.v[9], m.v[10], m.v[11], m.v[12], m.v[13], m.v[14], m.v[15])
#define EvalPrintQuaternion(q) printf("%s = \ts: %f,\tx: %f,\ty: %f,\tz: %f\n", #q, q.s, q.x, q.y, q.z)