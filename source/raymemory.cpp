internal void* Reserve (u64 size) {
	void* memory = malloc(size);
	return memory;
}

internal void Commit (void* memory, u64 size) {

}

internal void Decommit (void* memory, u64 size) {

}

internal void Release (void* memory) {
	free (memory);
}

internal M_Arena M_ArenaMakeReserve (u64 reserve_size) {
	M_Arena arena = {0};
	arena.base = (byte*) Reserve(reserve_size);
	arena.capacity = reserve_size;
	return arena;
}

internal M_Arena M_ArenaMake (void) {
	return M_ArenaMakeReserve(M_ARENA_DEFAULT_RESERVE_SIZE);
}

internal void M_ArenaRelease(Arena* arena) {
    Release(arena->base);
}

internal void* M_ArenaPush(M_Arena* arena, u64 size) {
    void* memory = 0;
    if (arena->alloc_pos + size > arena->commit_pos) {
		// round up commit size to a power of two on the granularity of the default commit size
        u64 commit_size = size;
        commit_size += M_ARENA_DEFAULT_COMMIT_SIZE - 1;
        commit_size -= commit_size % M_ARENA_DEFAULT_COMMIT_SIZE;

		// commit the memory
        Commit((u8*) arena->base + arena->commit_pos, commit_size);
        arena->commit_pos += commit_size;
    }

	// allocate memory by moving pointer
    memory = (u8*) arena->base + arena->alloc_pos;
	arena->previous_alloc_pos = arena->alloc_pos;
    arena->alloc_pos += size;

    return memory;
}

template<typename T, typename... Args>
internal T* M_ArenaPush(M_Arena* arena, Args&&... args) {
	void* memory = M_ArenaPush(arena, sizeof(T));
	return new (memory) T(args...);
}

internal void* M_ArenaPushZero(M_Arena* arena, u64 size) {
    void* result = M_ArenaPush(arena, size);
    MemoryZero(result, size);
    return result;
}

internal void M_ArenaPop(M_Arena* arena, u64 size) {
    if (size > arena->alloc_pos) {
        size = arena->alloc_pos;
    }

    arena->alloc_pos -= size;
}

internal void M_ArenaPopTo(M_Arena* arena, u64 alloc_position) {
    M_ArenaPop(arena, arena->alloc_pos - alloc_position);
}

void* M_ArenaRealloc (M_Arena* arena, void* old_memory, u64 old_size, u64 new_size) {
	void* memory = NULL;

	if (old_memory == NULL || old_size == 0) {
		// this is essentially an alloc call
		return M_ArenaPush(arena, new_size);
	}
	// if the old memory is within the arena's allocated bounds
	else if (arena->base <= old_memory && old_memory < arena->base + arena->alloc_pos) {
		// if the old memory was from the last allocation
		if (arena->base + arena->previous_alloc_pos == old_memory) {
			arena->alloc_pos = arena->previous_alloc_pos + new_size;
			if (new_size > old_size) {
				MemoryZero(arena->base + arena->alloc_pos, new_size - old_size);
			}
			memory = old_memory;
		}
		else {
			void* new_memory = M_ArenaPush(arena, new_size);
			u64 copy_size = Min(old_size, new_size);
			// copy old memory to new memory
			MemoryMove(new_memory, old_memory, copy_size);
			memory = new_memory;
		}
	}
	else {
		//LogError("Memory is out of bounds of the buffer in this arena.");
		Assert(0);
	}

	return memory;
}


internal M_TempArenaSavePoint M_TempArenaBegin (M_Arena* arena) {
	M_TempArenaSavePoint temp;
	temp.arena = arena;
	temp.position = arena->alloc_pos;
	return temp;
}

internal void M_TempArenaEnd (M_TempArenaSavePoint temp) {
	M_ArenaPopTo(temp.arena, temp.position);
}