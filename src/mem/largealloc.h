#pragma once

#include "../ds/flaglock.h"
#include "../ds/helpers.h"
#include "../ds/mpmcstack.h"
#include "../pal/pal.h"
#include "allocstats.h"
#include "baseslab.h"
#include "sizeclass.h"

#include <new>

namespace snmalloc
{
  template<class PAL>
  class MemoryProviderStateMixin;

  class Largeslab : public Baseslab
  {
    // This is the view of a contiguous memory area when it is being kept
    // in the global size-classed caches of available contiguous memory areas.
  private:
    template<class a, Construction c>
    friend class MPMCStack;
    template<class PAL>
    friend class MemoryProviderStateMixin;
    std::atomic<Largeslab*> next;

  public:
    void init()
    {
      kind = Large;
    }
  };

  /**
   * A slab that has been decommitted.  The first page remains committed and
   * the only fields that are guaranteed to exist are the kind and next
   * pointer from the superclass.
   */
  struct Decommittedslab : public Largeslab
  {
    /**
     * Constructor.  Expected to be called via placement new into some memory
     * that was formerly a superslab or large allocation and is now just some
     * spare address space.
     */
    Decommittedslab()
    {
      kind = Decommitted;
    }
  };

  // This represents the state that the large allcoator needs to add to the
  // global state of the allocator.  This is currently stored in the memory
  // provider, so we add this in.
  template<class PAL>
  class MemoryProviderStateMixin : public PAL
  {
    std::atomic_flag lock = ATOMIC_FLAG_INIT;
    void* bump;
    size_t remaining;

    void new_block()
    {
      size_t size = SUPERSLAB_SIZE;
      void* r = reserve<false>(&size, SUPERSLAB_SIZE);

      if (size < SUPERSLAB_SIZE)
        error("out of memory");

      PAL::template notify_using<NoZero>(r, OS_PAGE_SIZE);

      bump = r;
      remaining = size;
    }

    /**
     * The last time we saw a low memory notification.
     */
    std::atomic<uint64_t> last_low_memory_epoch = 0;
    std::atomic_flag lazy_decommit_guard;
    void lazy_decommit()
    {
      // If another thread is try to do lazy decommit, let it continue.  If
      // we try to parallelise this, we'll most likely end up waiting on the
      // same page table locks.
      if (!lazy_decommit_guard.test_and_set())
      {
        return;
      }
      // When we hit low memory, iterate over size classes and decommit all of
      // the memory that we can.  Start with the small size classes so that we
      // hit cached superslabs first.
      // FIXME: We probably shouldn't do this all at once.
      for (size_t large_class = 0; large_class < NUM_LARGE_CLASSES;
           large_class++)
      {
        if (!PAL::expensive_low_memory_check())
        {
          break;
        }
        size_t rsize = bits::one_at_bit(SUPERSLAB_BITS) << large_class;
        size_t decommit_size = rsize - OS_PAGE_SIZE;
        // Grab all of the chunks of this size class.
        auto* slab = large_stack[large_class].pop_all();
        while (slab)
        {
          // Decommit all except for the first page and then put it back on
          // the stack.
          if (slab->get_kind() != Decommitted)
          {
            PAL::notify_not_using(
              pointer_offset(slab, OS_PAGE_SIZE), decommit_size);
          }
          // Once we've removed these from the stack, there will be no
          // concurrent accesses and removal should have established a
          // happens-before relationship, so it's safe to use relaxed loads
          // here.
          auto next = slab->next.load(std::memory_order_relaxed);
          large_stack[large_class].push(new (slab) Decommittedslab());
          slab = next;
        }
      }
      lazy_decommit_guard.clear();
    }

  public:
    /**
     * Stack of large allocations that have been returned for reuse.
     */
    ModArray<NUM_LARGE_CLASSES, MPMCStack<Largeslab, PreZeroed>> large_stack;

    /**
     * Primitive allocator for structure that are required before
     * the allocator can be running.
     */
    template<typename T, size_t alignment, typename... Args>
    T* alloc_chunk(Args&&... args)
    {
      // Cache line align
      size_t size = bits::align_up(sizeof(T), 64);

      void* p;
      {
        FlagLock f(lock);

        if constexpr (alignment != 0)
        {
          uint8_t* aligned_bump = pointer_align_up<alignment, uint8_t>(bump);

          size_t bump_delta =
            static_cast<size_t>(aligned_bump - static_cast<uint8_t*>(bump));

          if (bump_delta > remaining)
          {
            new_block();
          }
          else
          {
            remaining -= bump_delta;
            bump = aligned_bump;
          }
        }

        if (remaining < size)
        {
          new_block();
        }

        p = bump;
        bump = pointer_offset(bump, size);
        remaining -= size;
      }

      auto page_start = pointer_align_down<OS_PAGE_SIZE, uint8_t>(p);
      auto page_end =
        pointer_align_up<OS_PAGE_SIZE, uint8_t>(pointer_offset(p, size));

      PAL::template notify_using<NoZero>(
        page_start, static_cast<size_t>(page_end - page_start));

      return new (p) T(std::forward<Args...>(args)...);
    }

    /**
     * Returns the number of low memory notifications that have been received
     * (over the lifetime of this process).  If the underlying system does not
     * support low memory notifications, this will return 0.
     */
    SNMALLOC_FAST_PATH
    uint64_t low_memory_epoch()
    {
      if constexpr (pal_supports<LowMemoryNotification, PAL>())
      {
        return PAL::low_memory_epoch();
      }
      else
      {
        return 0;
      }
    }

    template<bool committed>
    void* reserve(size_t* size, size_t align) noexcept
    {
      if constexpr (pal_supports<AlignedAllocation, PAL>())
      {
        return PAL::template reserve<committed>(size, align);
      }
      else
      {
        size_t request = *size;
        // Add align, so we can guarantee to provide at least size.
        request += align;
        // Alignment must be a power of 2.
        assert(align == bits::next_pow2(align));

        void* p = PAL::template reserve<committed>(&request);

        *size = request;
        auto p0 = address_cast(p);
        auto start = bits::align_up(p0, align);

        if (start > p0)
        {
          uintptr_t end = bits::align_down(p0 + request, align);
          *size = end - start;
          PAL::notify_not_using(p, start - p0);
          PAL::notify_not_using(pointer_cast<void>(end), (p0 + request) - end);
          p = pointer_cast<void>(start);
        }
        return p;
      }
    }

    SNMALLOC_FAST_PATH void lazy_decommit_if_needed()
    {
#ifdef TEST_LAZY_DECOMMIT
      static_assert(
        TEST_LAZY_DECOMMIT > 0,
        "TEST_LAZY_DECOMMIT must be a positive integer value.");
      static std::atomic<uint64_t> counter;
      auto c = counter++;
      if (c % TEST_LAZY_DECOMMIT == 0)
      {
        lazy_decommit();
      }
#else
      if constexpr (decommit_strategy == DecommitSuperLazy)
      {
        auto new_epoch = low_memory_epoch();
        auto old_epoch = last_low_memory_epoch.load(std::memory_order_acquire);
        if (new_epoch > old_epoch)
        {
          // Try to update the epoch to the value that we've seen.  If
          // another thread has seen a newer epoch than us (or done the same
          // update) let them win.
          do
          {
            if (last_low_memory_epoch.compare_exchange_strong(
                  old_epoch, new_epoch))
            {
              lazy_decommit();
            }
          } while (old_epoch <= new_epoch);
        }
      }
#endif
    }
  };

  using Stats = AllocStats<NUM_SIZECLASSES, NUM_LARGE_CLASSES>;

  enum AllowReserve
  {
    NoReserve,
    YesReserve
  };

  template<class MemoryProvider>
  class LargeAlloc
  {
    void* reserved_start = nullptr;
    void* reserved_end = nullptr;

  public:
    // This will be a zero-size structure if stats are not enabled.
    Stats stats;

    MemoryProvider& memory_provider;

    LargeAlloc(MemoryProvider& mp) : memory_provider(mp) {}

    template<AllowReserve allow_reserve>
    bool reserve_memory(size_t need, size_t add)
    {
      assert(reserved_start <= reserved_end);

      /*
       * Spell this comparison in terms of pointer subtraction like this,
       * rather than "reserved_start + need < reserved_end" becuase the
       * sum might not be representable on CHERI.
       */
      if (
        static_cast<size_t>(
          static_cast<char*>(reserved_end) -
          static_cast<char*>(reserved_start)) < need)
      {
        if constexpr (allow_reserve == YesReserve)
        {
          stats.segment_create();
          reserved_start =
            memory_provider.template reserve<false>(&add, SUPERSLAB_SIZE);
          reserved_end = pointer_offset(reserved_start, add);
          reserved_start = pointer_align_up<SUPERSLAB_SIZE>(reserved_start);

          if (add < need)
            return false;
        }
        else
        {
          return false;
        }
      }

      return true;
    }

    template<ZeroMem zero_mem = NoZero, AllowReserve allow_reserve = YesReserve>
    void* alloc(size_t large_class, size_t size)
    {
      size_t rsize = bits::one_at_bit(SUPERSLAB_BITS) << large_class;
      if (size == 0)
        size = rsize;

      void* p = memory_provider.large_stack[large_class].pop();
      memory_provider.lazy_decommit_if_needed();

      if (p == nullptr)
      {
        assert(reserved_start <= reserved_end);
        size_t add;

        if ((rsize + SUPERSLAB_SIZE) < RESERVE_SIZE)
          add = RESERVE_SIZE;
        else
          add = rsize + SUPERSLAB_SIZE;

        if (!reserve_memory<allow_reserve>(rsize, add))
          return nullptr;

        p = reserved_start;
        reserved_start = pointer_offset(p, rsize);

        stats.superslab_fresh();
        // All memory is zeroed since it comes from reserved space.
        memory_provider.template notify_using<NoZero>(p, size);
      }
      else
      {
        stats.superslab_pop();

        if constexpr (decommit_strategy == DecommitSuperLazy)
        {
          if (static_cast<Baseslab*>(p)->get_kind() == Decommitted)
          {
            // The first page is already in "use" for the stack element,
            // this will need zeroing for a YesZero call.
            if constexpr (zero_mem == YesZero)
              memory_provider.template zero<true>(p, OS_PAGE_SIZE);

            // Notify we are using the rest of the allocation.
            // Passing zero_mem ensures the PAL provides zeroed pages if
            // required.
            memory_provider.template notify_using<zero_mem>(
              pointer_offset(p, OS_PAGE_SIZE),
              bits::align_up(size, OS_PAGE_SIZE) - OS_PAGE_SIZE);
          }
          else
          {
            if constexpr (zero_mem == YesZero)
              memory_provider.template zero<true>(
                p, bits::align_up(size, OS_PAGE_SIZE));
          }
        }
        else if ((decommit_strategy != DecommitNone) || (large_class > 0))
        {
          // The first page is already in "use" for the stack element,
          // this will need zeroing for a YesZero call.
          if constexpr (zero_mem == YesZero)
            memory_provider.template zero<true>(p, OS_PAGE_SIZE);

          // Notify we are using the rest of the allocation.
          // Passing zero_mem ensures the PAL provides zeroed pages if required.
          memory_provider.template notify_using<zero_mem>(
            pointer_offset(p, OS_PAGE_SIZE),
            bits::align_up(size, OS_PAGE_SIZE) - OS_PAGE_SIZE);
        }
        else
        {
          // This is a superslab that has not been decommitted.
          if constexpr (zero_mem == YesZero)
            memory_provider.template zero<true>(
              p, bits::align_up(size, OS_PAGE_SIZE));
        }
      }

      return p;
    }

    void
    decommit_most(void* p, size_t large_class, MemoryProvider& memory_provider)
    {
      /* Decommit all but the first page of this object */
      memory_provider.notify_not_using(
        pointer_offset(p, OS_PAGE_SIZE),
        large_sizeclass_to_size(large_class) - OS_PAGE_SIZE);
    }

    void dealloc(void* p, size_t large_class)
    {
      stats.superslab_push();
      memory_provider.large_stack[large_class].push(static_cast<Largeslab*>(p));
      memory_provider.lazy_decommit_if_needed();
    }
  };

  using GlobalVirtual = MemoryProviderStateMixin<Pal>;
  /**
   * The memory provider that will be used if no other provider is explicitly
   * passed as an argument.
   */
  HEADER_GLOBAL GlobalVirtual default_memory_provider;
} // namespace snmalloc
