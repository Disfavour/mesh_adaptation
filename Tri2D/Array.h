// Array.h : BlockList container
//

#pragma once


template<typename T, size_t N> class BlockList
{
	struct Record
	{
		T arr[N];
		Record *next;
	};

public:
	class iterator
	{
		T *ptr, *last;

		iterator(T *arr, size_t i) : ptr(arr + i), last(arr + N)
		{
		}

		Record *start()
		{
			return reinterpret_cast<Record *>(last - N);
		}

	public:
		iterator(iterator &it) : ptr(it.ptr), last(it.last)
		{
		}

		iterator &operator = (iterator &it)
		{
			ptr = it.ptr;  last = it.last;  return *this;
		}

		operator T *()
		{
			return ptr;
		}

		T *operator -> ()
		{
			return ptr;
		}

		iterator operator ++ ()
		{
			if(++ptr < last)return *this;  ptr = start()->next->arr;
			last = ptr + N;  return *this;
		}

		iterator operator ++ (int)
		{
			iterator old = *this;  ++(*this);  return old;
		}

		bool operator == (iterator &it) const
		{
			return ptr == it.ptr;
		}

		bool operator != (iterator &it) const
		{
			return ptr != it.ptr;
		}

		friend class BlockList;
	};


private:
	Record *first, **next;
	iterator last;

	BlockList(BlockList &);  // TODO: NoCopy
	BlockList &operator = (BlockList &);

public:
	BlockList() : first(0), next(&first), last(0, 0)
	{
	}

	~BlockList()
	{
		while(first)
		{
			Record *del = first;  first = first->next;  delete del;
		}
	}


	iterator begin()
	{
		return iterator(first->arr, 0);
	}

	iterator end()
	{
		return last;
	}

	iterator append()
	{
		if(last)return last++;  Record *ptr = new Record;
		ptr->next = 0;  *next = ptr;  next = &ptr->next;
		last = iterator(ptr->arr, 1);  return iterator(ptr->arr, 0);
	}

	void resize(iterator end)
	{
		last = end;
	}

	void free()
	{
		if(!last)return;  next = &last.start()->next;
		for(Record *rec = *next; rec;)
		{
			Record *del = rec;  rec = rec->next;  delete del;
		}
		*next = 0;
	}
};


template<typename T, size_t N> class BlockListAllocator
{
	union MemoryBlock
	{
		char mem[sizeof(T)];
		MemoryBlock *next;
	};

	struct Record
	{
		MemoryBlock arr[N];
		Record *next;

		Record(MemoryBlock **prev) : next(0)
		{
			for(size_t i = 0; i < N; i++)
			{
				*prev = arr + i;  prev = &arr[i].next;
			}
			*prev = 0;
		}
	};


	Record *first, **next;
	MemoryBlock *free;

	BlockListAllocator(BlockListAllocator &);  // TODO: NoCopy
	BlockListAllocator &operator = (BlockListAllocator &);

	T *AllocMemory()
	{
		if(!free)
		{
			Record *ptr = new Record(&free);  *next = ptr;  next = &ptr->next;
		}
		T *ptr = reinterpret_cast<T *>(free);  free = free->next;  return ptr;
	}


public:
	struct Freer
	{
		BlockListAllocator &alloc;

		Freer(BlockListAllocator &a) : alloc(a)
		{
		}

		void operator () (T *ptr)
		{
			alloc.Free(ptr);
		}
	};


	BlockListAllocator() : first(0), next(&first), free(0)
	{
	}

	~BlockListAllocator()  // nested objects mush be already freed
	{
		while(first)
		{
			Record *del = first;  first = first->next;  delete del;
		}
	}

	T *Alloc()
	{
		T *ptr = AllocMemory();  new(ptr) T();  return ptr;
	}

	template<typename T1> T *Alloc(T1 p1)
	{
		T *ptr = AllocMemory();  new(ptr) T(p1);  return ptr;
	}

	template<typename T1, typename T2> T *Alloc(T1 p1, T2 p2)
	{
		T *ptr = AllocMemory();  new(ptr) T(p1, p2);  return ptr;
	}

	template<typename T1, typename T2, typename T3> T *Alloc(T1 p1, T2 p2, T3 p3)
	{
		T *ptr = AllocMemory();  new(ptr) T(p1, p2, p3);  return ptr;
	}

	template<typename T1, typename T2, typename T3, typename T4>
		T *Alloc(T1 p1, T2 p2, T3 p3, T4 p4)
	{
		T *ptr = AllocMemory();  new(ptr) T(p1, p2, p3, p4);  return ptr;
	}

	template<typename T1, typename T2, typename T3, typename T4, typename T5>
		T *Alloc(T1 p1, T2 p2, T3 p3, T4 p4, T5 p5)
	{
		T *ptr = AllocMemory();  new(ptr) T(p1, p2, p3, p5);  return ptr;
	}

	void Free(T *ptr)
	{
		if(!ptr)return;  ptr->~T();
		MemoryBlock *mem = reinterpret_cast<MemoryBlock *>(ptr);
		mem->next = free;  free = mem;
	}

	Freer GetFreer()
	{
		return Freer(*this);
	}
};
