#pragma once

#include <new>

#define M_ARENA_DEFAULT_COMMIT_SIZE Kilobytes(4)
#define M_ARENA_DEFAULT_RESERVE_SIZE Megabytes(64)

typedef struct M_Arena {
	byte* base;
	u64 capacity;
	u64 previous_alloc_pos;
	u64 alloc_pos;
	u64 commit_pos;
} Arena;

typedef struct M_TempArenaSavePoint {
	M_Arena* arena;
	u64 position;
} M_TempArenaSavepoint;
