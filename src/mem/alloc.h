#pragma once

#ifdef _MSC_VER
#  define ALLOCATOR __declspec(allocator)
#else
#  define ALLOCATOR
#endif

#include "../test/histogram.h"
#include "allocstats.h"
#include "largealloc.h"
#include "mediumslab.h"
#include "pagemap.h"
#include "pooled.h"
#include "remoteallocator.h"
#include "sizeclasstable.h"
#include "slab.h"

#if SNMALLOC_REVOKE_QUARANTINE == 1
extern "C"
{
#  if defined(__CHERI__)
#    include <cheri/cheric.h>
#  endif
#  include <cheri/caprevoke.h>
#  include <inttypes.h>
#  include <sys/caprevoke.h>
}
#endif

namespace snmalloc
{
  enum Boundary
  {
    /**
     * The location of the first byte of this allocation.
     */
    Start,
    /**
     * The location of the last byte of the allocation.
     */
    End,
    /**
     * The location one past the end of the allocation.  This is mostly useful
     * for bounds checking, where anything less than this value is safe.
     */
    OnePastEnd
  };

  enum PageMapSuperslabKind
  {
    PMNotOurs = 0,
    PMSuperslab = 1,
    PMMediumslab = 2

    /*
     * Values 3 (inclusive) through SUPERSLAB_BITS (exclusive) are as yet
     * unused.
     *
     * Values SUPERSLAB_BITS (inclusive) through 64 (exclusive, as it would
     * represent the entire address space) are used for log2(size) at the
     * heads of large allocations.  See SuperslabMap::set_large_size.
     *
     * Values 64 (inclusive) through 128 (exclusive) are used for entries
     * within a large allocation.  A value of x at pagemap entry p indicates
     * that there are at least 2^(x-64) (inclusive) and at most 2^(x+1-64)
     * (exclusive) page map entries between p and the start of the
     * allocation.  See SuperslabMap::set_large_size and external_address's
     * handling of large reallocation redirections.
     *
     * Values 128 (inclusive) through 255 (inclusive) are as yet unused.
     */

  };

  /* Ensure that PageMapSuperslabKind values are actually disjoint */
  static_assert(SUPERSLAB_BITS > 2, "Large allocations possibly too small");

  /*
   * In CHERI, we have to be able to rederive pointers to headers and
   * metadata given the address of the allocation, since the capabilities we
   * give out have bounds narrowed to the allocation itself.  Since snmalloc
   * already holds a map of the address space, here's a great place to do
   * that.  Rather than store sizes per each SUPERSLAB_SIZE sized piece of
   * memory, we store a capability.
   *
   * We have a lot of "metadata" bits at the least-significant end of the
   * address of a capability in this map, since its bounds must, at least,
   * cover a SUPERSLAB_SIZE sized object (or a large allocation).  To
   * minimize churn, we stash the existing enum PageMapSuperslabKind values
   * in the bottom 8 bits of the address.
   *
   * We could cut that down to 6 bits by reclaiming all values above 64; we
   * can test that the capability given to us to free is has address equal
   * to the base of the capability stored here in the page map.
   */
  static constexpr void* PMNotOursPtr = NULL;
  static constexpr int PAGEMAP_PTR_ALIGN = 0x100;

  using PAGEMAP_VALUE_T =
    std::conditional_t<SNMALLOC_PAGEMAP_POINTERS, void*, uint8_t>;
  static constexpr PAGEMAP_VALUE_T PAGEMAP_VALUE_Z = 0;

#ifndef SNMALLOC_MAX_FLATPAGEMAP_SIZE
// Use flat map is under a single node.
#  define SNMALLOC_MAX_FLATPAGEMAP_SIZE PAGEMAP_NODE_SIZE
#endif
  static constexpr bool USE_FLATPAGEMAP =
    (pal_supports<LazyCommit>() &&
     ((1 << 24) >= sizeof(FlatPagemap<SUPERSLAB_BITS, PAGEMAP_VALUE_T>))) ||
    (SNMALLOC_MAX_FLATPAGEMAP_SIZE >=
     sizeof(FlatPagemap<SUPERSLAB_BITS, PAGEMAP_VALUE_T>));

  using SuperslabPagemap = std::conditional_t<
    USE_FLATPAGEMAP,
    FlatPagemap<SUPERSLAB_BITS, PAGEMAP_VALUE_T>,
    Pagemap<SUPERSLAB_BITS, PAGEMAP_VALUE_T, PAGEMAP_VALUE_Z>>;

  HEADER_GLOBAL SuperslabPagemap global_pagemap;

  /**
   * Mixin used by `SuperslabMap` to directly access the pagemap via a global
   * variable.  This should be used from within the library or program that
   * owns the pagemap.
   */
  struct GlobalPagemap
  {
    /**
     * Returns the pagemap.
     */
    SuperslabPagemap& pagemap()
    {
      return global_pagemap;
    }
  };

  /**
   * Optionally exported function that accesses the global pagemap provided by
   * a shared library.
   */
  extern "C" void* snmalloc_pagemap_global_get(snmalloc::PagemapConfig const**);

  /**
   * Mixin used by `SuperslabMap` to access the global pagemap via a
   * type-checked C interface.  This should be used when another library (e.g.
   * your C standard library) uses snmalloc and you wish to use a different
   * configuration in your program or library, but wish to share a pagemap so
   * that either version can deallocate memory.
   */
  class ExternalGlobalPagemap
  {
    /**
     * A pointer to the pagemap.
     */
    SuperslabPagemap* external_pagemap;

  public:
    /**
     * Constructor.  Accesses the pagemap via the C ABI accessor and casts it to
     * the expected type, failing in cases of ABI mismatch.
     */
    ExternalGlobalPagemap()
    {
      const snmalloc::PagemapConfig* c;
      external_pagemap =
        SuperslabPagemap::cast_to_pagemap(snmalloc_pagemap_global_get(&c), c);
      // FIXME: Report an error somehow in non-debug builds.
      assert(external_pagemap);
    }

    /**
     * Returns the exported pagemap.
     */
    SuperslabPagemap& pagemap()
    {
      return *external_pagemap;
    }
  };

  /**
   * Class that defines an interface to the pagemap.  This is provided to
   * `Allocator` as a template argument and so can be replaced by a compatible
   * implementation (for example, to move pagemap updates to a different
   * protection domain).
   */
  template<typename PagemapProvider = GlobalPagemap>
  struct SuperslabMap : public PagemapProvider
  {
    using PagemapProvider::PagemapProvider;
    /**
     * Get the pagemap entry corresponding to a specific address.
     *
     * Despite the type, the return value is an enum PageMapSuperslabKind
     * or one of the reserved values described therewith.
     */
    uint8_t get(address_t p)
    {
      if constexpr (SNMALLOC_PAGEMAP_POINTERS)
      {
        return (uintptr_t)(PagemapProvider::pagemap().get(p));
      }
      else
      {
        return PagemapProvider::pagemap().get(p);
      }
    }

    /**
     * Get the pagemap entry corresponding to a specific address.
     */
    uint8_t get(void* p)
    {
      return get(address_cast(p));
    }

    /**
     * Fetch a pointer through the pagemap.  The return value and the
     * argument point at the same location, if the template parameter offset
     * is true, but may have different metadata associated with them, such
     * as CHERI bounds.  If the pagemap is storing pointers, then the
     * template parameter offset may be false, in which case, the return
     * value is exactly the pagemap entry, which will point at the bottom of
     * the pagemap granule and include the argument in its bounds.
     */
    template<bool offset = true>
    void* getp(void* p)
    {
      if constexpr (SNMALLOC_PAGEMAP_POINTERS)
      {
        void* pmp = pointer_align_down<PAGEMAP_PTR_ALIGN, void>(
          PagemapProvider::pagemap().get(address_cast(p)));
        if constexpr (offset)
        {
#if defined(__CHERI__)
          size_t delta = pointer_diff_without_provenance(pmp, p);
#else
          size_t delta = pointer_diff(pmp, p);
#endif
          return (void*)(pointer_offset(pmp, delta));
        }
        else
        {
          return pmp;
        }
      }
      else
      {
        static_assert(offset);
        return p;
      }
    }

    /**
     * Set a pagemap entry indicating that there is a superslab at the
     * specified index.
     */
    void set_slab(Superslab* slab)
    {
      if constexpr (SNMALLOC_PAGEMAP_POINTERS)
      {
        set(slab, pointer_offset(slab, PMSuperslab));
      }
      else
      {
        set(slab, static_cast<size_t>(PMSuperslab));
      }
    }
    /**
     * Add a pagemap entry indicating that a medium slab has been allocated.
     */
    void set_slab(Mediumslab* slab)
    {
      if constexpr (SNMALLOC_PAGEMAP_POINTERS)
      {
        set(slab, pointer_offset(slab, PMMediumslab));
      }
      else
      {
        set(slab, static_cast<size_t>(PMMediumslab));
      }
    }
    /**
     * Remove an entry from the pagemap corresponding to a superslab.
     */
    void clear_slab(Superslab* slab)
    {
      assert(get(slab) == PMSuperslab);
      if constexpr (SNMALLOC_PAGEMAP_POINTERS)
      {
        set(slab, PMNotOursPtr);
      }
      else
      {
        set(slab, static_cast<size_t>(PMNotOurs));
      }
    }
    /**
     * Remove an entry corresponding to a medium slab.
     */
    void clear_slab(Mediumslab* slab)
    {
      assert(get(slab) == PMMediumslab);
      if constexpr (SNMALLOC_PAGEMAP_POINTERS)
      {
        set(slab, PMNotOursPtr);
      }
      else
      {
        set(slab, static_cast<size_t>(PMNotOurs));
      }
    }
    /**
     * Update the pagemap to reflect a large allocation, of `size` bytes from
     * address `p`.
     */
    void set_large_size(void* p, size_t size)
    {
      size_t size_bits = bits::next_pow2_bits(size);
      if constexpr (SNMALLOC_PAGEMAP_POINTERS)
      {
        set(p, pointer_offset(p, size_bits));
        // Set redirect slide
        char* ss = reinterpret_cast<char*>(pointer_offset(p, SUPERSLAB_SIZE));
        for (size_t i = 0; i < size_bits - SUPERSLAB_BITS; i++)
        {
          size_t run = 1ULL << i;
          PagemapProvider::pagemap().set_range(
            address_cast(ss), pointer_offset(p, 64 + i + SUPERSLAB_BITS), run);
          ss += SUPERSLAB_SIZE * run;
        }
      }
      else
      {
        set(p, static_cast<uint8_t>(size_bits));
        // Set redirect slide
        auto ss = address_cast(p) + SUPERSLAB_SIZE;
        for (size_t i = 0; i < size_bits - SUPERSLAB_BITS; i++)
        {
          size_t run = 1ULL << i;
          PagemapProvider::pagemap().set_range(
            ss, static_cast<uint8_t>(64 + i + SUPERSLAB_BITS), run);
          ss = ss + SUPERSLAB_SIZE * run;
        }
      }
    }
    /**
     * Update the pagemap to remove a large allocation, of `size` bytes from
     * address `p`.
     */
    void clear_large_size(void* vp, size_t size)
    {
      auto p = address_cast(vp);
      size_t rounded_size = bits::next_pow2(size);
      assert(get(p) == bits::next_pow2_bits(size));
      auto count = rounded_size >> SUPERSLAB_BITS;
      if constexpr (SNMALLOC_PAGEMAP_POINTERS)
      {
        PagemapProvider::pagemap().set_range(p, PMNotOursPtr, count);
      }
      else
      {
        PagemapProvider::pagemap().set_range(p, PMNotOurs, count);
      }
    }

  private:
    /**
     * Helper function to set a pagemap entry.  This is not part of the public
     * interface and exists to make it easy to reuse the code in the public
     * methods in other pagemap adaptors.
     */
    void set(void* p, PAGEMAP_VALUE_T x)
    {
      PagemapProvider::pagemap().set(address_cast(p), x);
    }
  };

#ifndef SNMALLOC_DEFAULT_PAGEMAP
#  define SNMALLOC_DEFAULT_PAGEMAP snmalloc::SuperslabMap<>
#endif

  // This class is just used so that the free lists are the first entry
  // in the allocator and hence has better code gen.
  // It contains a free list per small size class.  These are used for
  // allocation on the fast path. This part of the code is inspired by mimalloc.
  class FastFreeLists
  {
  protected:
    FreeListHead small_fast_free_lists[NUM_SMALL_CLASSES];

  public:
    FastFreeLists() : small_fast_free_lists() {}
  };

  SNMALLOC_FAST_PATH void* no_replacement(void*)
  {
    return nullptr;
  }

  /**
   * Allocator.  This class is parameterised on three template parameters.  The
   * `MemoryProvider` defines the source of memory for this allocator.
   * Allocators try to reuse address space by allocating from existing slabs or
   * reusing freed large allocations.  When they need to allocate a new chunk
   * of memory they request space from the `MemoryProvider`.
   *
   * The `PageMap` parameter provides the adaptor to the pagemap.  This is used
   * to associate metadata with large (16MiB, by default) regions, allowing an
   * allocator to find the allocator responsible for that region.
   *
   * The next template parameter, `IsQueueInline`, defines whether the
   * message queue for this allocator should be stored as a field of the
   * allocator (`true`) or provided externally, allowing it to be anywhere else
   * in the address space (`false`).
   *
   * The final template parameter provides a hook to allow the allocator in use
   * to be dynamically modified.  This is used to implement a trick from
   * mimalloc that avoids a conditional branch on the fast path.  We initialise
   * the thread-local allocator pointer with the address of a global allocator,
   * which never owns any memory.  When we try to allocate memory, we call the
   * replacement function.
   */
  template<
    class MemoryProvider = GlobalVirtual,
    class PageMap = SNMALLOC_DEFAULT_PAGEMAP,
    bool IsQueueInline = true,
    void* (*Replacement)(void*) = no_replacement>
  class Allocator
  : public FastFreeLists,
    public Pooled<
      Allocator<MemoryProvider, PageMap, IsQueueInline, Replacement>>
  {
    LargeAlloc<MemoryProvider> large_allocator;
    PageMap page_map;

#if SNMALLOC_QUARANTINE_DEALLOC == 1

    /*
     * In our design for a temporally-safe allocator, objects (more
     * properly, virtual addresses) transition through a QUARANTINED state
     * between ALLOCATED and FREE.  They linger here until we are sure that
     * there are no pointers *to* them surviving outside the TCB (i.e., held
     * by clients of the allocator; naturally, perhaps, the kernel and we
     * both retain access to these addresses, but we trust both to not
     * mishandle their authority).  In particular, we ensure this by
     * forcibly removing such pointers from the system, a process we call
     * "revocation".
     *
     * In CHERI, we have the ability to generate pointers which cannot
     * reference outside an allocation (see SNMALLOC_CHERI_SETBOUNDS); we
     * refer to these as "external" pointers and will refer to the
     * whole-slab pointers we use internally as "high" pointers (we avoid
     * the use of "internal", as that might mean "internal to an
     * allocation", and of "privileged" as all pointers, in some sense, are
     * privilege-bearing; naming is hard).  Occasionally, we have
     * "high-external" pointers which are intermediate forms: they have the
     * privileges of high pointers but bounds of an allocation.
     *
     * These two aspects interact in quarantining memory:
     *
     *   dealloc() now transitions memory from ALLOCATED to QUARANTINED.  It
     *   takes an "external" pointer from the client and rederives a "high"
     *   pointer to the slab/large allocation in use.  It stores this "high"
     *   pointer to quarantine, as that will survive revocation.  There is
     *   no need to store the object size, as it can be easily recovered
     *   from actions already necessary after revocation.
     *
     *   After revocation, the "high" pointers of quarantine are fed to the
     *   dealloc_real() function.
     *
     *   This separation unfortunately means that static size information is
     *   of little use to dealloc().  We preserve the increasingly-static
     *   dealloc_real() functions for comparison, but they will not be
     *   called if a quarantine is in use.
     */

    /*
     * A piece of the quarantine region for this allocator.  Objects in
     * quarantine have had their revocation requested but are waiting
     * for notification that this has happened.
     *
     * A quarantine node undergoes the following lifecycle:
     *
     *   FILLING: free()'d objects are written here while there is room.
     *
     *   FULL: The revocation epoch has been measured and this has been
     *         threaded onto the linked list of waiting nodes.  This list
     *         is used as a queue, with older entries are at the front, so
     *         that the allocator can quickly find the collection of nodes
     *         that can advance to FREE.
     *
     *   FREE: Sufficient revocation epochs have passed that this node
     *         can be freed back to general circulation.
     *
     * XXX What do we do about epoch rollover?  It'd be really nice if we
     * could force other threads to advance, somehow... can we take a
     * lightweight mutex on entry, maybe?  (We might need such a thing to
     * properly defend against untrusted code anyway??)  A better idea would
     * be to use a capability itself to mark that an epoch had definitely
     * passed.  That seems like a post-ASPLOS idea.
     *
     * XXX Right now, we're using the quarantine as double-free protection,
     * because that was the least invasive change.  There are other schemes
     * with deeper quarantine that could be used instead, but they have
     * additional subtlety and/or complexity.  (If one is naive about
     * quarantine, double-free is especially insidious: it's possible to put
     * something into quarantine in two epochs, such that an earlier epoch
     * causes the shoot-down of the re-allocated object in a later epoch!
     * That seems messy and worth commenting on in the paper!)
     *
     */
    struct QuarantineNode
    {
      /* DLL of this Allocator's FULL quarantine nodes */
      QuarantineNode* prev;
      QuarantineNode* next;

      uint64_t full_epoch; /* When did this fill? */
      size_t footprint; /* How many bytes of the heap do we represent? */
      uint32_t self_foot; /* What's the size of this object? */
      uint16_t n_ents; /* How many pointers are there in this node? */
      uint16_t first_ent; /* How full did we get?  0: completely */

      /*
       * A QuarantineNode is followed by n_ents struct QuarantineEntry-s.
       * We omit the field here and just access past the end of this
       * structure.
       */
    };

    struct QuarantineEntry
    {
      void* privp; /* Internal pointer to object, with slab permissions */
      union
      {
#  if (SNMALLOC_REVOKE_QUARANTINE == 1)
        uint8_t* revbitmap; /* Direct access to large object's shadow */
#  endif
        uint8_t pmsc; /* Page map size class */
      } addl;
#  if (SNMALLOC_REVOKE_PARANOIA == 1)
      void* origp;
#  endif
    };

    class QuarantineState
    {
      DLList<struct QuarantineNode> waiting;
      uint8_t n_waiting;
      uint64_t waiting_footprint;

      struct QuarantineNode* filling;
      uint16_t filling_left;

      static struct QuarantineNode* newqn(Allocator* a, uint8_t sc)
      {
        struct QuarantineNode* qn;
        size_t size = sizeclass_to_size(sc);

        qn = static_cast<struct QuarantineNode*>(
          a->medium_alloc<YesZero, YesReserve>(sc, size, size));

#  if defined(__CHERI_PURE_CAPABILITY__) && (SNMALLOC_PAGEMAP_REDERIVE == 1)
        qn = static_cast<struct QuarantineNode*>(cheri_csetbounds(qn, size));
#  endif

        if (qn != nullptr)
        {
          qn->n_ents = (sizeclass_to_size(sc) - sizeof(struct QuarantineNode)) /
            sizeof(struct QuarantineEntry);
          qn->self_foot = sizeclass_to_size(sc);
          qn->footprint = qn->self_foot;
          qn->full_epoch = UINT64_MAX;
        }

        return qn;
      }

      static void
      deqqn(Allocator* a, struct QuarantineNode* qn, uint16_t initpos)
      {
        uint16_t qix = initpos;
        struct QuarantineEntry* q = reinterpret_cast<struct QuarantineEntry*>(
          pointer_offset(qn, sizeof(struct QuarantineNode)));
#  if SNMALLOC_REVOKE_CHATTY == 1
        uint64_t bitmap_cycles = 0;
        uint64_t cyc_start;
#  endif

        q = pointer_offset(q, initpos * sizeof(struct QuarantineEntry));

        for (; qix < qn->n_ents; qix++, q++)
        {
          void* privp = q->privp;
          q->privp = nullptr;

#  if defined(__CHERI_PURE_CAPABILITY__)
          assert(cheri_gettag(privp));
#  endif

#  if SNMALLOC_REVOKE_QUARANTINE == 1
          uint8_t* revbitmap;
          void* p;
#  endif

#  if SNMALLOC_REVOKE_QUARANTINE == 1
          if (cheri_gettag(q->addl.revbitmap))
#  else
          if (q->addl.pmsc >= SUPERSLAB_BITS)
#  endif
          {
            /*
             * This is a large object and we have been handed its shadow
             * bitmap access.  privp can be used directly, as its bounds are
             * precisely those of the large object.
             */
#  if SNMALLOC_REVOKE_QUARANTINE == 1
            revbitmap = q->addl.revbitmap;
            p = privp;
            a->large_dealloc(privp, cheri_getlen(privp));
#  else
            a->large_dealloc(
              privp, large_sizeclass_to_size(q->addl.pmsc - SUPERSLAB_BITS));
#  endif
          }
          else if (q->addl.pmsc == PMSuperslab)
          {
#  if SNMALLOC_REVOKE_QUARANTINE == 1
            Superslab* super = Superslab::get(privp);
            Metaslab& meta = super->get_meta(Slab::get(privp));
            sizeclass_t sizeclass = meta.sizeclass;

            revbitmap = super->get_revbitmap();
            p = cheri_csetbounds(privp, sizeclass_to_size(sizeclass));
#  endif
            a->dealloc_real_small(privp);
          }
          else
          {
            assert(q->addl.pmsc == PMMediumslab);
#  if SNMALLOC_REVOKE_QUARANTINE == 1
            Mediumslab* slab = Mediumslab::get(privp);
            sizeclass_t sizeclass = slab->get_sizeclass();

            revbitmap = slab->get_revbitmap();
            p = cheri_csetbounds(privp, sizeclass_to_size(sizeclass));
#  endif
            a->dealloc_real_medium(privp);
          }

#  if SNMALLOC_REVOKE_QUARANTINE == 1
#    if SNMALLOC_REVOKE_CHATTY == 1
          cyc_start = AAL::tick();
#    endif
          caprev_shadow_nomap_clear(reinterpret_cast<uint64_t*>(revbitmap), p);
#    if SNMALLOC_REVOKE_CHATTY == 1
          bitmap_cycles += AAL::tick() - cyc_start;
#    endif

#    if SNMALLOC_REVOKE_PARANOIA == 1
          /* Verify that the original pointer has had its tag cleared */
          assert(!cheri_gettag(q->origp));
#    endif
#  endif
        }

#  if SNMALLOC_REVOKE_CHATTY == 1
        fprintf(
          stderr,
          "dequar: a=%p foot=0x%zx bmcyc=0x%" PRIx64 "\n",
          a,
          qn->footprint,
          bitmap_cycles);
#  endif
      }

#  if SNMALLOC_REVOKE_CHATTY == 1
      void print_revoke_stats(
        FILE* f,
        const char* what,
        Allocator* a,
        struct caprevoke_stats* crst,
        uint64_t ccount)
      {
#    ifdef __mips__
        //__asm__ __volatile__("li $0, 0xdead"); // trace off
#    endif
        fprintf(
          f,
          "revoke: %s"
          " a=%p"
          " einit=0x%03" PRIx64 " efini=0x%03" PRIx64 " scand=0x%" PRIx64
          " pfro=0x%" PRIx64 " pfrw=0x%" PRIx64 " pfsk=0x%" PRIx64
          " clook=0x%" PRIx64 " cnuke=0x%" PRIx64 " tpage=0x%" PRIx64
          " tcrev=0x%" PRIx64 "\n",
          what,
          a,
          crst->epoch_init,
          crst->epoch_fini,
          crst->pages_scanned,
          crst->pages_faulted_ro,
          crst->pages_faulted_rw,
          crst->pages_fault_skip,
          crst->caps_found,
          crst->caps_cleared,
          crst->page_scan_cycles,
          ccount);
#    ifdef __mips__
        //__asm__ __volatile__("li $0, 0xbeef"); // trace on
#    endif
      }
#  endif

      static bool epoch_clears(uint64_t now, uint64_t test)
      {
#  if SNMALLOC_REVOKE_QUARANTINE == 1
        return caprevoke_epoch_clears(now, test);
#  else
        (void)now;
        (void)test;
        return 1;
#  endif
      }

      /*
       * Given the current epoch clock, pull off sufficiently old quarantine
       * nodes and push all the pointers therein back to free.
       */
      void drain(Allocator* a, uint64_t epoch)
      {
        struct QuarantineNode* qn = waiting.get_head();
        while ((qn != nullptr) && epoch_clears(epoch, qn->full_epoch))
        {
          assert(waiting_footprint >= qn->footprint);
          waiting_footprint -= qn->footprint;

          deqqn(a, qn, qn->first_ent);

          waiting.pop();
          n_waiting--;
          struct QuarantineNode* qnext = waiting.get_head();

          if (filling == nullptr)
          {
            qn->prev = nullptr;
            qn->next = nullptr;
            qn->full_epoch = UINT64_MAX;
            qn->footprint = qn->self_foot;
            filling = qn;
            filling_left = filling->n_ents;
          }
          else
          {
#  if defined(__CHERI_PURE_CAPABILITY__) && (SNMALLOC_PAGEMAP_REDERIVE == 1)
            qn = static_cast<struct QuarantineNode*>(a->pagemap().getp(qn));
#  endif
            a->dealloc_real(qn);
          }

          qn = qnext;
        }
      }

      /*
       * Alright, we tried, but we are facing pressure; finish the
       * epoch of our oldest waiting node.
       */
      void quarantine_step_drain(Allocator* a, const char* cause)
      {
        uint64_t epoch;

        (void)cause;
#  if SNMALLOC_REVOKE_QUARANTINE == 1
        struct caprevoke_stats crst;
        QuarantineNode* qn = waiting.get_head();
        do
        {
#    if SNMALLOC_REVOKE_CHATTY == 1
          uint64_t cyc_init = AAL::tick();
#    endif
#    if SNMALLOC_REVOKE_DRY_RUN == 0
          int res = caprevoke(CAPREVOKE_LAST_PASS, qn->full_epoch, &crst);
          UNUSED(res);
          assert(res == 0);
#    endif
#    if SNMALLOC_REVOKE_CHATTY == 1
          uint64_t cyc_fini = AAL::tick();
          print_revoke_stats(stderr, "stdr", a, &crst, cyc_fini - cyc_init);
#    endif
        } while (!epoch_clears(crst.epoch_fini, qn->full_epoch));
        epoch = crst.epoch_fini;
#  else
        epoch = 4;
#  endif
        assert(n_waiting > 0);
        drain(a, epoch);

        a->handle_message_queue();
      }

      void quarantine_step(Allocator* a)
      {
        uint64_t epoch;

#  if SNMALLOC_REVOKE_QUARANTINE == 1
        {
          struct caprevoke_stats crst;
#    if SNMALLOC_REVOKE_CHATTY == 1
          uint64_t cyc_init = AAL::tick();
#    endif
#    if SNMALLOC_REVOKE_DRY_RUN == 0
          int res = caprevoke(
            CAPREVOKE_IGNORE_START | CAPREVOKE_ONLY_IF_OPEN, 0, &crst);
          UNUSED(res);
          assert(res == 0);
#    endif
#    if SNMALLOC_REVOKE_CHATTY == 1
          uint64_t cyc_fini = AAL::tick();
          print_revoke_stats(stderr, "step", a, &crst, cyc_fini - cyc_init);
#    endif
          filling->full_epoch = crst.epoch_init;
          epoch = crst.epoch_fini;
        }
#  else
        epoch = 4;
#  endif

#  if SNMALLOC_REVOKE_CHATTY == 1
        fprintf(stderr, "enquar: foot=0x%zx\n", filling->footprint);
#  endif

        filling->first_ent = filling_left;
        waiting.insert_back(filling);
        waiting_footprint += filling->footprint;
        filling = nullptr;
        n_waiting++;

#  if SNMALLOC_REVOKE_QUARANTINE == 1
        /*
         * This one might be a no-op, if no nodes are sufficiently old.
         * Do this only if we are revoking, so that we exercise the
         * multiple node paths on non-revoking platforms.
         */

        drain(a, epoch);
#  endif

        if ((filling == nullptr) && (n_waiting > 1)) // XXX n_waiting vs. ?
        {
          quarantine_step_drain(a, "full");
        }
        else
        {
          QuarantineNode* qn = newqn(a, NUM_SMALL_CLASSES);
          if (qn == nullptr)
          {
            quarantine_step_drain(a, "noqn");
          }
          else
          {
            filling = qn;
            filling_left = qn->n_ents;
          }
        }

        assert(filling != nullptr);
      }

    public:
      void init(Allocator* a)
      {
        n_waiting = 0;
        waiting_footprint = 0;

        filling = newqn(a, NUM_SMALL_CLASSES);
        assert(filling != nullptr);

        filling_left = filling->n_ents;
      }

      void quarantine(
        Allocator* a,
        uint8_t* revbitmap,
        void* privp,
        void* privpred,
        void* p,
        size_t psize,
        uint8_t pmsc)
      {
#  if SNMALLOC_REVOKE_QUARANTINE == 1
        assert(revbitmap != nullptr);
        /*
         * Paint the revocation bitmap using the external pointer so that we
         * correctly fence against concurrent double-free and frees from
         * earlier revocation epochs.
         */
        if (
          caprev_shadow_nomap_set(
            reinterpret_cast<uint64_t*>(revbitmap), privpred, p) != 0)
        {
          return;
        }
#  else
        UNUSED(revbitmap);
        UNUSED(privpred);
        UNUSED(p);
        UNUSED(privp);
        UNUSED(pmsc);
#  endif

        struct QuarantineEntry* q = reinterpret_cast<struct QuarantineEntry*>(
          pointer_offset(filling, sizeof(struct QuarantineNode)));

        filling_left--;
        if ((pmsc == PMSuperslab) || (pmsc == PMMediumslab))
        {
          /*
           * Small objects store their pmsc to avoid another lookup and
           * a sufficiently privileged pointer that we can get back to
           * the containing slab's headers.
           */
          q[filling_left].addl.pmsc = pmsc;
          q[filling_left].privp = privp;
        }
        else
        {
          /*
           * Large objects store their revocation shadow capability directly
           * and the reduced-bounds privileged pointer which encapsulates
           * the entire object.  There is no need to reach outside this
           * region going forward: snmalloc does not coalesce large objects
           * nor manage them in aggregation.
           */
#  if SNMALLOC_REVOKE_QUARANTINE == 1
          q[filling_left].addl.revbitmap = revbitmap;
#  else
          q[filling_left].addl.pmsc = pmsc;
#  endif
          q[filling_left].privp = privpred;
        }
#  if (SNMALLOC_REVOKE_PARANOIA == 1)
        q[filling_left].origp = p;
#  endif

        filling->footprint += psize;

        if (filling_left == 0)
        {
          quarantine_step(a);
        }
      }

      /*
       * Like drain(), but unconditional on epoch and also all the pointers
       * from the filling qn.
       *
       * Like it says on the tin, only for debugging use.
       */
      void debug_drain_all(Allocator* a)
      {
        struct QuarantineNode* qn = waiting.get_head();

#  if SNMALLOC_REVOKE_QUARANTINE == 1
        struct caprevoke_stats crst;
        uint64_t start_epoch;
        bool did_revoke = false;

#    if SNMALLOC_REVOKE_DRY_RUN == 0
        int res = caprevoke(
          CAPREVOKE_NO_WAIT_OK | CAPREVOKE_IGNORE_START |
            CAPREVOKE_LAST_NO_EARLY,
          0,
          &crst);
        UNUSED(res);
        assert(res == 0);
        start_epoch = crst.epoch_init;
#    else
        start_epoch = 0;
        crst.epoch_fini = 4;
#    endif
#  endif

        /* Drain queue, forcing revocations as we go */
        while (qn != nullptr)
        {
          waiting.pop();
          n_waiting--;
          struct QuarantineNode* qnext = waiting.get_head();

#  if SNMALLOC_REVOKE_QUARANTINE == 1
          while (!epoch_clears(crst.epoch_fini, qn->full_epoch))
          {
#    if SNMALLOC_REVOKE_CHATTY == 1
            uint64_t cyc_init = AAL::tick();
#    endif
#    if SNMALLOC_REVOKE_DRY_RUN == 0
            int res = caprevoke(CAPREVOKE_LAST_PASS, qn->full_epoch, &crst);
            UNUSED(res);
            assert(res == 0);
#    endif
#    if SNMALLOC_REVOKE_CHATTY == 1
            uint64_t cyc_fini = AAL::tick();
            print_revoke_stats(stderr, "dbgw", a, &crst, cyc_fini - cyc_init);
#    endif
            did_revoke = true;
          }
#  endif

          assert(waiting_footprint >= qn->footprint);
          waiting_footprint -= qn->footprint;
          deqqn(a, qn, qn->first_ent);

#  if defined(__CHERI_PURE_CAPABILITY__) && (SNMALLOC_PAGEMAP_REDERIVE == 1)
          qn = static_cast<struct QuarantineNode*>(a->pagemap().getp(qn));
#  endif
          a->dealloc_real(qn);

          qn = qnext;
        }
        assert(n_waiting == 0);
        assert(waiting_footprint == 0);

        if (filling_left != filling->n_ents)
        {
#  if SNMALLOC_REVOKE_QUARANTINE == 1
          while (!did_revoke && !epoch_clears(crst.epoch_fini, start_epoch))
          {
#    if SNMALLOC_REVOKE_CHATTY == 1
            uint64_t cyc_init = AAL::tick();
#    endif
#    if SNMALLOC_REVOKE_DRY_RUN == 0
            int res = caprevoke(
              CAPREVOKE_LAST_PASS | CAPREVOKE_IGNORE_START, start_epoch, &crst);
            UNUSED(res);
            assert(res == 0);
#    endif
#    if SNMALLOC_REVOKE_CHATTY == 1
            uint64_t cyc_fini = AAL::tick();
            print_revoke_stats(stderr, "dbgf", a, &crst, cyc_fini - cyc_init);
#    endif
          }
#  endif
          deqqn(a, filling, filling_left);
          filling_left = filling->n_ents;
        }
      }
    };

    QuarantineState quarantine;
#endif

  public:
    Stats& stats()
    {
      return large_allocator.stats;
    }

    template<class MP>
    friend class AllocPool;

    // Allocate memory of a statically known size.
    template<
      size_t ssize,
      ZeroMem zero_mem = NoZero,
      AllowReserve allow_reserve = YesReserve>
    SNMALLOC_FAST_PATH ALLOCATOR void* alloc()
    {
      static_assert(ssize != 0, "Size must not be zero.");
#ifdef USE_MALLOC
      static_assert(
        allow_reserve == YesReserve,
        "When passing to malloc, cannot require NoResereve");
      if constexpr (zero_mem == NoZero)
        return malloc(ssize);
      else
        return calloc(1, ssize);
#else

#  if SNMALLOC_CHERI_ALIGN == 1
      constexpr size_t size =
        bits::align_up_const(ssize, 1 << CHERI_ALIGN_SHIFT(ssize));
#  else
      constexpr size_t size = ssize;
#  endif

      constexpr sizeclass_t sizeclass = size_to_sizeclass_const(size);
      void* ret;

      stats().alloc_request(size);

      if constexpr (sizeclass < NUM_SMALL_CLASSES)
      {
        ret = small_alloc<zero_mem, allow_reserve>(size);
      }
      else if constexpr (sizeclass < NUM_SIZECLASSES)
      {
        handle_message_queue();
        constexpr size_t rsize = sizeclass_to_size(sizeclass);
        ret = medium_alloc<zero_mem, allow_reserve>(sizeclass, rsize, size);
      }
      else
      {
        handle_message_queue();
        ret = large_alloc<zero_mem, allow_reserve>(size);
      }

#  if SNMALLOC_CHERI_SETBOUNDS == 1
      ret = (ret == NULL) ?
        ret :
        cheri_andperm(
          cheri_csetbounds(ret, size),
          CHERI_PERMS_USERSPACE_DATA & ~CHERI_PERM_CHERIABI_VMMAP);
#  endif

      return ret;

#endif
    }

    // Allocate memory of a dynamically known size.
    template<ZeroMem zero_mem = NoZero, AllowReserve allow_reserve = YesReserve>
    SNMALLOC_FAST_PATH ALLOCATOR void* alloc(size_t dsize)
    {
#ifdef USE_MALLOC
      static_assert(
        allow_reserve == YesReserve,
        "When passing to malloc, cannot require NoResereve");
      if constexpr (zero_mem == NoZero)
        return malloc(dsize);
      else
        return calloc(1, dsize);
#else

#  if SNMALLOC_CHERI_ALIGN == 1
      int shift = CHERI_ALIGN_SHIFT(dsize);
      size_t size = bits::align_up(dsize, 1 << shift);
#  else
      size_t size = dsize;
#  endif

      void* ret;

      stats().alloc_request(size);

      // Perform the - 1 on size, so that zero wraps around and ends up on
      // slow path.
      if (likely((size - 1) <= (sizeclass_to_size(NUM_SMALL_CLASSES - 1) - 1)))
      {
        // Allocations smaller than the slab size are more likely. Improve
        // branch prediction by placing this case first.
        ret = small_alloc<zero_mem, allow_reserve>(size);
      }
      else
      {
        ret = alloc_not_small<zero_mem, allow_reserve>(size);
      }

#  if SNMALLOC_CHERI_SETBOUNDS == 1
      ret = (ret == NULL) ?
        ret :
        cheri_andperm(
          cheri_csetbounds(ret, size),
          CHERI_PERMS_USERSPACE_DATA & ~CHERI_PERM_CHERIABI_VMMAP);
#  endif

      return ret;
    }

    template<ZeroMem zero_mem = NoZero, AllowReserve allow_reserve = YesReserve>
    SNMALLOC_SLOW_PATH ALLOCATOR void* alloc_not_small(size_t size)
    {
      handle_message_queue();

      if (size == 0)
      {
        return small_alloc<zero_mem, allow_reserve>(1);
      }

      sizeclass_t sizeclass = size_to_sizeclass(size);
      if (sizeclass < NUM_SIZECLASSES)
      {
        size_t rsize = sizeclass_to_size(sizeclass);
        return medium_alloc<zero_mem, allow_reserve>(sizeclass, rsize, size);
      }
      else
      {
        return large_alloc<zero_mem, allow_reserve>(size);
      }

#endif
    }

#if defined(CHECK_CLIENT) && (SNMALLOC_UNSAFE_FREES_CHECK == 1)
  private:
    void verify_size(void* p, size_t claimed)
    {
      void* privp = p;
      if constexpr (SNMALLOC_PAGEMAP_REDERIVE)
      {
        privp = page_map.getp(p);
      }

      uint8_t pmsc = pagemap().get(address_cast(privp));
      if (pmsc == PMSuperslab)
      {
        Superslab* super = Superslab::get(privp);
        Slab* slab = Slab::get(privp);
        Metaslab& meta = super->get_meta(slab);
        uint8_t sc = meta.sizeclass;
        assert(sc == size_to_sizeclass(claimed));
      }
      else if (pmsc == PMMediumslab)
      {
        Mediumslab* slab = Mediumslab::get(privp);
        size_t sc = slab->get_sizeclass();
        assert(sc == size_to_sizeclass(claimed));
      }
      else
      {
        assert(pmsc < 64);
        assert(pmsc == bits::next_pow2_bits(claimed));
      }
    }
#endif

#if SNMALLOC_UNSAFE_FREES == 1
  public:
    /*
     * Free memory of a statically known size. Must be called with an
     * external pointer.
     */
    template<size_t size>
    void dealloc(void* p)
    {
#  if defined(CHECK_CLIENT) && (SNMALLOC_UNSAFE_FREES_CHECK == 1)
      verify_size(p, size);
#  endif

#  ifdef USE_MALLOC
      UNUSED(size);
      return free(p);
#  else
#    if SNMALLOC_QUARANTINE_DEALLOC == 0
      handle_message_queue();
#    endif

      void* privp = p;
      if constexpr (SNMALLOC_PAGEMAP_REDERIVE)
      {
        privp = page_map.getp(p);
      }

#    if SNMALLOC_QUARANTINE_DEALLOC == 1

      uint8_t pmsc;
      constexpr uint8_t sizeclass = size_to_sizeclass_const(size);
      uint8_t* revbitmap = nullptr;
      void* privpred = privp;

      if constexpr (sizeclass < NUM_SMALL_CLASSES)
      {
#      if SNMALLOC_REVOKE_QUARANTINE == 1
        Superslab* super = Superslab::get(privp);
        revbitmap = super->get_revbitmap();
#      endif
        pmsc = PMSuperslab;
      }
      else if constexpr (sizeclass < NUM_SIZECLASSES)
      {
#      if SNMALLOC_REVOKE_QUARANTINE == 1
        Mediumslab* slab = Mediumslab::get(privp);
        revbitmap = slab->get_revbitmap();
#      endif
        pmsc = PMMediumslab;
      }
      else
      {
#      if SNMALLOC_REVOKE_QUARANTINE == 1
        // XXX System call!  Would rather fuse into madvise(), if possible.
        int res = caprevoke_shadow(
          CAPREVOKE_SHADOW_NOVMMAP,
          privp,
          reinterpret_cast<void**>(&revbitmap));
        (void)res; /* quiet NDEBUG builds */
        assert(res == 0);
        privpred = cheri_csetbounds(privp, size);
#      endif
        pmsc = page_map.get(p);
      }

      quarantine.quarantine(this, revbitmap, privp, privpred, p, size, pmsc);

#    else
      dealloc_real<size>(privp);
#    endif
#  endif
    }

  private:
    template<size_t size>
    void dealloc_real(void* p)
    {
      // p is a "high" pointer, not "external"; xref predealloc

      constexpr sizeclass_t sizeclass = size_to_sizeclass_const(size);

      if constexpr (sizeclass < NUM_SMALL_CLASSES)
      {
        Superslab* super = Superslab::get(p);
        RemoteAllocator* target = super->get_allocator();

        if (target == public_state())
          small_dealloc(super, p, sizeclass);
        else
          remote_dealloc(target, p, sizeclass);
      }
      else if constexpr (sizeclass < NUM_SIZECLASSES)
      {
        Mediumslab* slab = Mediumslab::get(p);
        RemoteAllocator* target = slab->get_allocator();

        if (target == public_state())
          medium_dealloc(slab, p, sizeclass);
        else
          remote_dealloc(target, p, sizeclass);
      }
      else
      {
        large_dealloc(p, size);
      }
    }

  public:
    /*
     * Free memory of a dynamically known size. Must be called with an
     * external pointer.
     */
    void dealloc(void* p, size_t size)
    {
#  if defined(CHECK_CLIENT) && (SNMALLOC_UNSAFE_FREES_CHECK == 1)
      verify_size(p, size);
#  endif

#  ifdef USE_MALLOC
      UNUSED(size);
      return free(p);
#  else
#    if SNMALLOC_QUARANTINE_DEALLOC == 0
      handle_message_queue();
#    endif

      void* privp = p;
      if constexpr (SNMALLOC_PAGEMAP_REDERIVE)
      {
        privp = page_map.getp(p);
      }

#    if SNMALLOC_QUARANTINE_DEALLOC == 1

      uint8_t pmsc;
      uint8_t sizeclass = size_to_sizeclass(size);
      uint8_t* revbitmap = nullptr;
      void* privpred = privp;

      if (sizeclass < NUM_SMALL_CLASSES)
      {
#      if SNMALLOC_REVOKE_QUARANTINE == 1
        Superslab* super = Superslab::get(privp);
        revbitmap = super->get_revbitmap();
#      endif
        pmsc = PMSuperslab;
      }
      else if (sizeclass < NUM_SIZECLASSES)
      {
#      if SNMALLOC_REVOKE_QUARANTINE == 1
        Mediumslab* slab = Mediumslab::get(privp);
        revbitmap = slab->get_revbitmap();
#      endif
        pmsc = PMMediumslab;
      }
      else
      {
#      if SNMALLOC_REVOKE_QUARANTINE == 1
        // XXX System call!  Would rather fuse into madvise(), if possible.
        int res = caprevoke_shadow(
          CAPREVOKE_SHADOW_NOVMMAP,
          privp,
          reinterpret_cast<void**>(&revbitmap));
        (void)res; /* quiet NDEBUG builds */
        assert(res == 0);
        privpred = cheri_csetbounds(privp, size);
#      endif
        pmsc = page_map.get(p);
      }

      quarantine.quarantine(this, revbitmap, privp, privpred, p, size, pmsc);
#    else
      dealloc_real(privp, size);
#    endif
#  endif
    }

  private:
    void dealloc_real(void* p, size_t size)
    {
      sizeclass_t sizeclass = size_to_sizeclass(size);

      if (sizeclass < NUM_SMALL_CLASSES)
      {
        Superslab* super = Superslab::get(p);
        RemoteAllocator* target = super->get_allocator();

        if (target == public_state())
          small_dealloc(super, p, sizeclass);
        else
          remote_dealloc(target, p, sizeclass);
      }
      else if (sizeclass < NUM_SIZECLASSES)
      {
        Mediumslab* slab = Mediumslab::get(p);
        RemoteAllocator* target = slab->get_allocator();

        if (target == public_state())
          medium_dealloc(slab, p, sizeclass);
        else
          remote_dealloc(target, p, sizeclass);
      }
      else
      {
        large_dealloc(p, size);
      }
    }

#else /* SNMALLOC_UNSAFE_FREES */

  public:
    template<size_t size>
    void dealloc(void* p)
    {

#  if defined(CHECK_CLIENT) && (SNMALLOC_UNSAFE_FREES_CHECK == 1)
      verify_sizeclass(p, size_to_sizeclass_const(size));
#  endif

      dealloc(p);
    }

    void dealloc(void* p, size_t size)
    {

#  if defined(CHECK_CLIENT) && (SNMALLOC_UNSAFE_FREES_CHECK == 1)
      verify_sizeclass(p, size_to_sizeclass_const(size));
#  endif

      (void)size;
      dealloc(p);
    }

#endif /* SNMALLOC_UNSAFE_FREES */

  public:
    /*
     * Free memory of an unknown size. Must be called with an external
     * pointer.
     */
    SNMALLOC_FAST_PATH void dealloc(void* p)
    {
#ifdef USE_MALLOC
      return free(p);
#else

      void* privp = p;
      if constexpr (SNMALLOC_PAGEMAP_REDERIVE)
      {
        privp = page_map.getp(p);
      }

#  if SNMALLOC_QUARANTINE_DEALLOC == 1
      size_t size;
      uint8_t pmsc = pagemap().get(address_cast(privp));
      if (pmsc == PMNotOurs)
      {
        error("Not allocated by this allocator");
      }

      uint8_t* revbitmap = nullptr;
      void* privpred = privp;

      if (pmsc == PMSuperslab)
      {
        Superslab* super = Superslab::get(privp);
        size = sizeclass_to_size(super->get_meta(Slab::get(privp)).sizeclass);
#    if SNMALLOC_REVOKE_QUARANTINE == 1
        revbitmap = super->get_revbitmap();
#    endif
      }
      else if (pmsc == PMMediumslab)
      {
        Mediumslab* slab = Mediumslab::get(privp);
        size = sizeclass_to_size(slab->get_sizeclass());
#    if SNMALLOC_REVOKE_QUARANTINE == 1
        revbitmap = slab->get_revbitmap();
#    endif
      }
      else
      {
        size = large_sizeclass_to_size(pmsc - SUPERSLAB_BITS);
#    if SNMALLOC_REVOKE_QUARANTINE == 1
        // XXX System call!  Would rather fuse into madvise(), if possible.
        int res = caprevoke_shadow(
          CAPREVOKE_SHADOW_NOVMMAP,
          privp,
          reinterpret_cast<void**>(&revbitmap));
        (void)res; /* quiet NDEBUG builds */
        assert(res == 0);
        privpred = cheri_csetbounds(privp, size);
#    endif
      }

      quarantine.quarantine(this, revbitmap, privp, privpred, p, size, pmsc);
#  else
      dealloc_real(privp);
#  endif
#endif
    }

  private:
    SNMALLOC_FAST_PATH void dealloc_real(void* p)
    {
      uint8_t size = pagemap().get(address_cast(p));

      if (likely(size == PMSuperslab))
      {
        dealloc_real_small(p);
        return;
      }
      dealloc_not_small(p, size);
    }

    SNMALLOC_FAST_PATH void dealloc_real_small(void* p)
    {
      Superslab* super = Superslab::get(p);
      RemoteAllocator* target = super->get_allocator();
      Slab* slab = Slab::get(p);
      Metaslab& meta = super->get_meta(slab);

      // Reading a remote sizeclass won't fail, since the other allocator
      // can't reuse the slab, as we have not yet deallocated this
      // pointer.
      sizeclass_t sizeclass = meta.sizeclass;

      if (likely(target == public_state()))
        small_dealloc(super, p, sizeclass);
      else
        remote_dealloc(target, p, sizeclass);
    }

    SNMALLOC_SLOW_PATH void dealloc_not_small(void* p, uint8_t size)
    {
#if SNMALLOC_QUARANTINE_DEALLOC == 0
      /*
       * If we're quarantining, we'll process the message queue more
       * sporadically, per quarantine drain, rather than per call to
       * dealloc_real.  See quarantine_step_drain, above.
       */
      handle_message_queue();
#endif

      if (p == nullptr)
        return;

      if (size == PMMediumslab)
      {
        dealloc_real_medium(p);
        return;
      }

      if (size == 0)
      {
        error("Not allocated by this allocator");
      }

#ifdef CHECK_CLIENT
      Superslab* super = Superslab::get(p);
      if (size > 64 || address_cast(super) != address_cast(p))
      {
        error("Not deallocating start of an object");
      }
#endif

      large_dealloc(p, 1ULL << size);
    }

    inline void dealloc_real_medium(void* p)
    {
      Mediumslab* slab = Mediumslab::get(p);
      RemoteAllocator* target = slab->get_allocator();

      // Reading a remote sizeclass won't fail, since the other allocator
      // can't reuse the slab, as we have not yet deallocated this pointer.
      sizeclass_t sizeclass = slab->get_sizeclass();

      if (target == public_state())
        medium_dealloc(slab, p, sizeclass);
      else
        remote_dealloc(target, p, sizeclass);
    }

  public:
    template<Boundary location = Start>
    address_t external_address(void* p)
    {
#ifdef USE_MALLOC
      error("Unsupported");
      UNUSED(p);
#else
      uint8_t size = pagemap().get(address_cast(p));

      void* privp = p;
      if constexpr (SNMALLOC_PAGEMAP_REDERIVE)
      {
        privp = page_map.getp(p);
      }

#  if defined(__CHERI__)
      /*
       * It is possible that the user-provided capability p is untagged; it
       * might, for example, have been revoked.  If that's true now, bail
       * out.  Otherwise, it's possible that p might become revoked while
       * we're working, but we work with privp when accessing internal
       * structures and tread carefully, ensuring that our answer is in
       * terms of the original p, thereby preventing aliasing due to races
       * with this function.  It's possible that we will reveal something
       * about our internal state, but no user data and no capabilities
       * will flow to the caller.
       *
       * We might data-race with another thread here if p is freed during
       * operation, but nothing that thread does should cause us to fault;
       * the headers we're after are in the first page of a superslab,
       * which are never decommitted.
       */

      if (!cheri_gettag(p))
      {
        return address_cast(static_cast<void*>(nullptr));
      }
#  endif

      Superslab* super = Superslab::get(privp);
      if (size == PMSuperslab)
      {
        Slab* slab = Slab::get(privp);
        Metaslab& meta = super->get_meta(slab);

        sizeclass_t sc = meta.sizeclass;
        void* slab_end = pointer_offset(slab, SLAB_SIZE);

        return external_pointer<location>(p, sc, slab_end);
      }
      if (size == PMMediumslab)
      {
        Mediumslab* slab = Mediumslab::get(privp);

        sizeclass_t sc = slab->get_sizeclass();
        void* slab_end = pointer_offset(slab, SUPERSLAB_SIZE);

        return external_pointer<location>(p, sc, slab_end);
      }

      address_t ss;

      if constexpr (1 && SNMALLOC_PAGEMAP_POINTERS)
      {
        /*
         * There's no reason to do anything logarithmic; we're directly
         * storing the pointer to the start here.  Just use that.
         */

        ss = address_cast(pagemap().template getp<false>(super));
        size = pagemap().get(ss);
      }
      else
      {
        ss = address_cast(super);
        while (size > 64)
        {
          // This is a large alloc redirect.
          ss = ss - (1ULL << (size - 64));
          size = pagemap().get(ss);
        }
      }

      if constexpr (0 && SNMALLOC_PAGEMAP_POINTERS)
      {
        assert(
          (size == 0) ||
          (ss == address_cast(pagemap().template getp<false>(super))));
      }

      if (size == 0)
      {
        if constexpr ((location == End) || (location == OnePastEnd))
          // We don't know the End, so return MAX_PTR
          return UINTPTR_MAX;
        else
          // We don't know the Start, so return MIN_PTR
          return 0;
      }

#  if defined(__CHERI__)
      ss = address_cast(cheri_setaddress(p, ss));
#  endif

      // This is a large alloc, mask off to the slab size.
      if constexpr (location == Start)
        return ss;
      else if constexpr (location == End)
        return (ss + (1ULL << size) - 1ULL);
      else
        return (ss + (1ULL << size));
#endif
    }

    template<Boundary location = Start>
    void* external_pointer(void* p)
    {
      return pointer_cast<void>(external_address<location>(p));
    }

    size_t alloc_size(void* p)
    {
      // This must be called on an external pointer.
      size_t size = pagemap().get(address_cast(p));

      if constexpr (SNMALLOC_PAGEMAP_REDERIVE)
      {
        p = page_map.getp(p);
      }

      if (size == 0)
      {
        error("Not allocated by this allocator");
      }
      else if (size == PMSuperslab)
      {
        Superslab* super = Superslab::get(p);

        // Reading a remote sizeclass won't fail, since the other allocator
        // can't reuse the slab, as we have no yet deallocated this pointer.
        Slab* slab = Slab::get(p);
        Metaslab& meta = super->get_meta(slab);

        return sizeclass_to_size(meta.sizeclass);
      }
      else if (size == PMMediumslab)
      {
        Mediumslab* slab = Mediumslab::get(p);
        // Reading a remote sizeclass won't fail, since the other allocator
        // can't reuse the slab, as we have no yet deallocated this pointer.
        return sizeclass_to_size(slab->get_sizeclass());
      }

      return 1ULL << size;
    }

    size_t get_id()
    {
      return id();
    }

  private:
    using alloc_id_t = typename Remote::alloc_id_t;

    /*
     * A singly-linked list of Remote objects, supporting append and
     * take-all operations.  Intended only for the private use of this
     * allocator; the Remote objects here will later be taken and pushed
     * to the inter-thread message queues.
     */
    struct RemoteList
    {
      /*
       * A stub Remote object that will always be the head of this list;
       * never taken for further processing.
       */
      Remote head;

      Remote* last;

      RemoteList()
      {
        clear();
      }

      void clear()
      {
        last = &head;
      }

      bool empty()
      {
        return last == &head;
      }
    };

    struct RemoteCache
    {
      /**
       * The total amount of memory stored awaiting dispatch to other
       * allocators.  This is initialised to the maximum size that we use
       * before caching so that, when we hit the slow path and need to dispatch
       * everything, we can check if we are a real allocator and lazily provide
       * a real allocator.
       */
      size_t size = REMOTE_CACHE;
      RemoteList list[REMOTE_SLOTS];

      /// Used to find the index into the array of queues for remote
      /// deallocation
      /// r is used for which round of sending this is.
      inline size_t get_slot(size_t id, size_t r)
      {
        constexpr size_t allocator_size = sizeof(
          Allocator<MemoryProvider, PageMap, IsQueueInline, Replacement>);
        constexpr size_t initial_shift =
          bits::next_pow2_bits_const(allocator_size);
        assert((initial_shift - (r * REMOTE_SLOT_BITS)) < 64);
        return (id >> (initial_shift + (r * REMOTE_SLOT_BITS))) & REMOTE_MASK;
      }

      SNMALLOC_FAST_PATH void
      dealloc_sized(alloc_id_t target_id, void* p, size_t objectsize)
      {
        this->size += objectsize;

        Remote* r = static_cast<Remote*>(p);
        r->set_target_id(target_id);
        assert(r->target_id() == target_id);

        RemoteList* l = &list[get_slot(target_id, 0)];
        l->last->non_atomic_next = r;
        l->last = r;
      }

      SNMALLOC_FAST_PATH void
      dealloc(alloc_id_t target_id, void* p, sizeclass_t sizeclass)
      {
        dealloc_sized(target_id, p, sizeclass_to_size(sizeclass));
      }

      void post(alloc_id_t id)
      {
        // When the cache gets big, post lists to their target allocators.
        size = 0;

        size_t post_round = 0;

        while (true)
        {
          auto my_slot = get_slot(id, post_round);

          for (size_t i = 0; i < REMOTE_SLOTS; i++)
          {
            if (i == my_slot)
              continue;

            RemoteList* l = &list[i];
            Remote* first = l->head.non_atomic_next;

            if (!l->empty())
            {
              // Send all slots to the target at the head of the list.
              Superslab* super = Superslab::get(first);
              super->get_allocator()->message_queue.enqueue(first, l->last);
              l->clear();
            }
          }

          RemoteList* resend = &list[my_slot];
          if (resend->empty())
            break;

          // Entries could map back onto the "resend" list,
          // so take copy of the head, mark the last element,
          // and clear the original list.
          Remote* r = resend->head.non_atomic_next;
          resend->last->non_atomic_next = nullptr;
          resend->clear();

          post_round++;

          while (r != nullptr)
          {
            // Use the next N bits to spread out remote deallocs in our own
            // slot.
            size_t slot = get_slot(r->target_id(), post_round);
            RemoteList* l = &list[slot];
            l->last->non_atomic_next = r;
            l->last = r;

            r = r->non_atomic_next;
          }
        }
      }
    };

    SlabList small_classes[NUM_SMALL_CLASSES];
    DLList<Mediumslab> medium_classes[NUM_MEDIUM_CLASSES];

    DLList<Superslab> super_available;
    DLList<Superslab> super_only_short_available;

    RemoteCache remote;

    std::conditional_t<IsQueueInline, RemoteAllocator, RemoteAllocator*>
      remote_alloc;

#ifdef CACHE_FRIENDLY_OFFSET
    size_t remote_offset = 0;

    void* apply_cache_friendly_offset(void* p, sizeclass_t sizeclass)
    {
      size_t mask = sizeclass_to_cache_friendly_mask(sizeclass);

      size_t offset = remote_offset & mask;
      remote_offset += CACHE_FRIENDLY_OFFSET;

      return (void*)((uintptr_t)p + offset);
    }
#else
    void* apply_cache_friendly_offset(void* p, sizeclass_t sizeclass)
    {
      UNUSED(sizeclass);
      return p;
    }
#endif

    auto* public_state()
    {
      if constexpr (IsQueueInline)
      {
        return &remote_alloc;
      }
      else
      {
        return remote_alloc;
      }
    }

    alloc_id_t id()
    {
      return public_state()->id();
    }

    auto& message_queue()
    {
      return public_state()->message_queue;
    }

    template<class A, class MemProvider>
    friend class Pool;

  public:
    Allocator(
      MemoryProvider& m,
      PageMap&& p = PageMap(),
      RemoteAllocator* r = nullptr,
      bool isFake = false)
    : large_allocator(m), page_map(p)
    {
      if constexpr (IsQueueInline)
      {
        assert(r == nullptr);
        (void)r;
      }
      else
      {
        remote_alloc = r;
      }

      if (id() >= static_cast<alloc_id_t>(-1))
        error("Id should not be -1");

      // If this is fake, don't do any of the bits of initialisation that may
      // allocate memory.
      if (isFake)
        return;

      init_message_queue();
      message_queue().invariant();

#if SNMALLOC_QUARANTINE_DEALLOC == 1
      quarantine.init(this);
#endif

#ifndef NDEBUG
      for (sizeclass_t i = 0; i < NUM_SIZECLASSES; i++)
      {
        size_t size = sizeclass_to_size(i);
        sizeclass_t sc1 = size_to_sizeclass(size);
        sizeclass_t sc2 = size_to_sizeclass_const(size);
        size_t size1 = sizeclass_to_size(sc1);
        size_t size2 = sizeclass_to_size(sc2);

        // All medium size classes are page aligned.
        if (i > NUM_SMALL_CLASSES)
        {
          assert(bits::is_aligned_block<OS_PAGE_SIZE>(nullptr, size1));
        }

        assert(sc1 == i);
        assert(sc1 == sc2);
        assert(size1 == size);
        assert(size1 == size2);
      }
#endif
    }

    template<Boundary location>
    static uintptr_t
    external_pointer(void* p, sizeclass_t sizeclass, void* end_point)
    {
      size_t rsize = sizeclass_to_size(sizeclass);

      void* end_point_correction = location == End ?
        (static_cast<uint8_t*>(end_point) - 1) :
        (location == OnePastEnd ? end_point :
                                  (static_cast<uint8_t*>(end_point) - rsize));

      ptrdiff_t offset_from_end =
        (static_cast<uint8_t*>(end_point) - 1) - static_cast<uint8_t*>(p);

      size_t end_to_end =
        round_by_sizeclass(rsize, static_cast<size_t>(offset_from_end));

      uint8_t* priv_res =
        static_cast<uint8_t*>(end_point_correction) - end_to_end;

#if defined(__CHERI__)
      uint8_t* res =
        static_cast<uint8_t*>(cheri_setaddress(p, cheri_getaddress(priv_res)));
#else
      uint8_t* res = priv_res;
#endif

      return address_cast<uint8_t>(res);
    }

    void init_message_queue()
    {
      // Manufacture an allocation to prime the queue
      // Using an actual allocation removes a conditional of a critical path.
      //
      // Bypass the alloc() wrapper and go straight for the internal
      // small_alloc() so that CHERI bounds are not set on the resulting
      // internal object, so that it is suitable to be passed directly to
      // the internal deallocation routines which expect to be able to align
      // to find Superslab headers.
      Remote* dummy = reinterpret_cast<Remote*>(
        small_alloc<YesZero, YesReserve>(MIN_ALLOC_SIZE));
      dummy->set_target_id(id());
      message_queue().init(dummy);
    }

    SNMALLOC_FAST_PATH void handle_dealloc_remote(Remote* p)
    {
      Superslab* super = Superslab::get(p);

#ifdef CHECK_CLIENT
      if (p->target_id() != super->get_allocator()->id())
        error("Detected memory corruption.  Potential use-after-free");
#endif
      if (likely(super->get_kind() == Super))
      {
        Slab* slab = Slab::get(p);
        Metaslab& meta = super->get_meta(slab);
        if (likely(p->target_id() == id()))
        {
          small_dealloc_offseted(super, p, meta.sizeclass);
          return;
        }
      }
      handle_dealloc_remote_slow(p);
    }

    SNMALLOC_SLOW_PATH void handle_dealloc_remote_slow(Remote* p)
    {
      Superslab* super = Superslab::get(p);
      if (likely(super->get_kind() == Medium))
      {
        Mediumslab* slab = Mediumslab::get(p);
        if (p->target_id() == id())
        {
          sizeclass_t sizeclass = slab->get_sizeclass();
          void* start = remove_cache_friendly_offset(p, sizeclass);
          medium_dealloc(slab, start, sizeclass);
        }
        else
        {
          // Queue for remote dealloc elsewhere.
          remote.dealloc(p->target_id(), p, slab->get_sizeclass());
        }
      }
      else
      {
        assert(likely(p->target_id() != id()));
        Slab* slab = Slab::get(p);
        Metaslab& meta = super->get_meta(slab);
        // Queue for remote dealloc elsewhere.
        remote.dealloc(p->target_id(), p, meta.sizeclass);
      }
    }

    SNMALLOC_SLOW_PATH void handle_message_queue_inner()
    {
      for (size_t i = 0; i < REMOTE_BATCH; i++)
      {
        auto r = message_queue().dequeue();

        if (unlikely(!r.second))
          break;

        handle_dealloc_remote(r.first);
      }

      // Our remote queues may be larger due to forwarding remote frees.
      if (likely(remote.size < REMOTE_CACHE))
        return;

      stats().remote_post();
      remote.post(id());
    }

    SNMALLOC_FAST_PATH void handle_message_queue()
    {
      // Inline the empty check, but not necessarily the full queue handling.
      if (likely(message_queue().is_empty()))
        return;

      handle_message_queue_inner();
    }

    template<AllowReserve allow_reserve>
    Superslab* get_superslab()
    {
      Superslab* super = super_available.get_head();

      if (super != nullptr)
        return super;

      super = reinterpret_cast<Superslab*>(
        large_allocator.template alloc<NoZero, allow_reserve>(
          0, SUPERSLAB_SIZE));

      if ((allow_reserve == NoReserve) && (super == nullptr))
        return super;

      uint8_t* revbitmap = nullptr;
#if SNMALLOC_REVOKE_QUARANTINE == 1
      int res = caprevoke_shadow(
        CAPREVOKE_SHADOW_NOVMMAP, super, reinterpret_cast<void**>(&revbitmap));
      (void)res; /* quiet NDEBUG builds */
      assert(res == 0);
#endif

      super->init(public_state(), revbitmap);
      pagemap().set_slab(super);
      super_available.insert(super);
      return super;
    }

    void reposition_superslab(Superslab* super)
    {
      switch (super->get_status())
      {
        case Superslab::Full:
        {
          // Remove from the list of superslabs that have available slabs.
          super_available.remove(super);
          break;
        }

        case Superslab::Available:
        {
          // Do nothing.
          break;
        }

        case Superslab::OnlyShortSlabAvailable:
        {
          // Move from the general list to the short slab only list.
          super_available.remove(super);
          super_only_short_available.insert(super);
          break;
        }

        case Superslab::Empty:
        {
          // Can't be empty since we just allocated.
          error("Unreachable");
          break;
        }
      }
    }

    template<AllowReserve allow_reserve>
    Slab* alloc_slab(sizeclass_t sizeclass)
    {
      stats().sizeclass_alloc_slab(sizeclass);
      if (Superslab::is_short_sizeclass(sizeclass))
      {
        // Pull a short slab from the list of superslabs that have only the
        // short slab available.
        Superslab* super = super_only_short_available.pop();

        if (super != nullptr)
        {
          Slab* slab =
            super->alloc_short_slab(sizeclass, large_allocator.memory_provider);
          assert(super->is_full());
          return slab;
        }

        super = get_superslab<allow_reserve>();

        if ((allow_reserve == NoReserve) && (super == nullptr))
          return nullptr;

        Slab* slab =
          super->alloc_short_slab(sizeclass, large_allocator.memory_provider);
        reposition_superslab(super);
        return slab;
      }

      Superslab* super = get_superslab<allow_reserve>();

      if ((allow_reserve == NoReserve) && (super == nullptr))
        return nullptr;

      Slab* slab =
        super->alloc_slab(sizeclass, large_allocator.memory_provider);
      reposition_superslab(super);
      return slab;
    }

    template<ZeroMem zero_mem, AllowReserve allow_reserve>
    SNMALLOC_FAST_PATH void* small_alloc(size_t size)
    {
      MEASURE_TIME_MARKERS(
        small_alloc,
        4,
        16,
        MARKERS(
          zero_mem == YesZero ? "zeromem" : "nozeromem",
          allow_reserve == NoReserve ? "noreserve" : "reserve"));

      SNMALLOC_ASSUME(size <= SLAB_SIZE);
      sizeclass_t sizeclass = size_to_sizeclass(size);

      assert(sizeclass < NUM_SMALL_CLASSES);
      auto& fl = small_fast_free_lists[sizeclass];
      void* head = fl.value;
      if (likely(head != nullptr))
      {
        // Read the next slot from the memory that's about to be allocated.
        fl.value = Metaslab::follow_next(head);

        void* p = remove_cache_friendly_offset(head, sizeclass);
        if constexpr (zero_mem == YesZero)
        {
          large_allocator.memory_provider.zero(p, size);
        }
        stats().sizeclass_alloc(sizeclass);
        return p;
      }

      return small_alloc_slow<zero_mem, allow_reserve>(sizeclass);
    }

    template<ZeroMem zero_mem, AllowReserve allow_reserve>
    SNMALLOC_SLOW_PATH void* small_alloc_slow(sizeclass_t sizeclass)
    {
      if (void* replacement = Replacement(this))
      {
        return reinterpret_cast<Allocator*>(replacement)
          ->template small_alloc_slow<zero_mem, allow_reserve>(sizeclass);
      }
      handle_message_queue();
      size_t rsize = sizeclass_to_size(sizeclass);
      auto& sl = small_classes[sizeclass];
      stats().sizeclass_alloc(sizeclass);

      Slab* slab;

      if (!sl.is_empty())
      {
        SlabLink* link = sl.get_head();
        slab = link->get_slab();
      }
      else
      {
        slab = alloc_slab<allow_reserve>(sizeclass);

        if ((allow_reserve == NoReserve) && (slab == nullptr))
          return nullptr;

        sl.insert_back(slab->get_link());
      }
      auto& ffl = small_fast_free_lists[sizeclass];
      return slab->alloc<zero_mem>(
        sl, ffl, rsize, large_allocator.memory_provider);
    }

    SNMALLOC_FAST_PATH void
    small_dealloc(Superslab* super, void* p, sizeclass_t sizeclass)
    {
#ifdef CHECK_CLIENT
      Slab* slab = Slab::get(p);
      if (!slab->is_start_of_object(super, p))
      {
        error("Not deallocating start of an object");
      }
#endif

      void* offseted = apply_cache_friendly_offset(p, sizeclass);
      small_dealloc_offseted(super, offseted, sizeclass);
    }

    SNMALLOC_FAST_PATH void
    small_dealloc_offseted(Superslab* super, void* p, sizeclass_t sizeclass)
    {
      MEASURE_TIME(small_dealloc, 4, 16);
      stats().sizeclass_dealloc(sizeclass);

      Slab* slab = Slab::get(p);
      if (likely(slab->dealloc_fast(super, p)))
        return;

      small_dealloc_offseted_slow(super, p, sizeclass);
    }

    SNMALLOC_SLOW_PATH void small_dealloc_offseted_slow(
      Superslab* super, void* p, sizeclass_t sizeclass)
    {
      bool was_full = super->is_full();
      SlabList* sl = &small_classes[sizeclass];
      Slab* slab = Slab::get(p);
      Superslab::Action a =
        slab->dealloc_slow(sl, super, p, large_allocator.memory_provider);
      if (likely(a == Superslab::NoSlabReturn))
        return;
      stats().sizeclass_dealloc_slab(sizeclass);

      if (a == Superslab::NoStatusChange)
        return;

      switch (super->get_status())
      {
        case Superslab::Full:
        {
          error("Unreachable");
          break;
        }

        case Superslab::Available:
        {
          if (was_full)
          {
            super_available.insert(super);
          }
          else
          {
            super_only_short_available.remove(super);
            super_available.insert(super);
          }
          break;
        }

        case Superslab::OnlyShortSlabAvailable:
        {
          super_only_short_available.insert(super);
          break;
        }

        case Superslab::Empty:
        {
          super_available.remove(super);

          if constexpr (decommit_strategy == DecommitSuper)
          {
            large_allocator.memory_provider.notify_not_using(
              pointer_offset(super, OS_PAGE_SIZE),
              SUPERSLAB_SIZE - OS_PAGE_SIZE);
          }
          else if constexpr (decommit_strategy == DecommitSuperLazy)
          {
            static_assert(
              pal_supports<LowMemoryNotification, MemoryProvider>(),
              "A lazy decommit strategy cannot be implemented on platforms "
              "without low memory notifications");
          }

          pagemap().clear_slab(super);
          large_allocator.dealloc(super, 0);
          stats().superslab_push();
          break;
        }
      }
    }

    template<ZeroMem zero_mem, AllowReserve allow_reserve>
    void* medium_alloc(sizeclass_t sizeclass, size_t rsize, size_t size)
    {
      MEASURE_TIME_MARKERS(
        medium_alloc,
        4,
        16,
        MARKERS(
          zero_mem == YesZero ? "zeromem" : "nozeromem",
          allow_reserve == NoReserve ? "noreserve" : "reserve"));

      sizeclass_t medium_class = sizeclass - NUM_SMALL_CLASSES;

      DLList<Mediumslab>* sc = &medium_classes[medium_class];
      Mediumslab* slab = sc->get_head();
      void* p;

      if (slab != nullptr)
      {
        p = slab->alloc<zero_mem>(size, large_allocator.memory_provider);

        if (slab->full())
          sc->pop();
      }
      else
      {
        if (void* replacement = Replacement(this))
        {
          return reinterpret_cast<Allocator*>(replacement)
            ->template medium_alloc<zero_mem, allow_reserve>(
              sizeclass, rsize, size);
        }
        slab = reinterpret_cast<Mediumslab*>(
          large_allocator.template alloc<NoZero, allow_reserve>(
            0, SUPERSLAB_SIZE));

        if ((allow_reserve == NoReserve) && (slab == nullptr))
          return nullptr;

        uint8_t* revbitmap = nullptr;
#if SNMALLOC_REVOKE_QUARANTINE == 1
        int res = caprevoke_shadow(
          CAPREVOKE_SHADOW_NOVMMAP, slab, reinterpret_cast<void**>(&revbitmap));
        (void)res; /* quiet NDEBUG builds */
        assert(res == 0);
#endif

        slab->init(public_state(), revbitmap, sizeclass, rsize);
        pagemap().set_slab(slab);
        p = slab->alloc<zero_mem>(size, large_allocator.memory_provider);

        if (!slab->full())
          sc->insert(slab);
      }

      stats().sizeclass_alloc(sizeclass);
      return p;
    }

    void medium_dealloc(Mediumslab* slab, void* p, sizeclass_t sizeclass)
    {
      MEASURE_TIME(medium_dealloc, 4, 16);
      stats().sizeclass_dealloc(sizeclass);
      slab->decommit(p, large_allocator.memory_provider);
      bool was_full = slab->dealloc(p);

#ifdef CHECK_CLIENT
      if (!is_multiple_of_sizeclass(
            sizeclass_to_size(sizeclass),
            pointer_diff(p, pointer_offset(slab, SUPERSLAB_SIZE))))
      {
        error("Not deallocating start of an object");
      }
#endif

      if (slab->empty())
      {
        if (!was_full)
        {
          sizeclass_t medium_class = sizeclass - NUM_SMALL_CLASSES;
          DLList<Mediumslab>* sc = &medium_classes[medium_class];
          sc->remove(slab);
        }

        if constexpr (decommit_strategy == DecommitSuper)
        {
          large_allocator.memory_provider.notify_not_using(
            pointer_offset(slab, OS_PAGE_SIZE), SUPERSLAB_SIZE - OS_PAGE_SIZE);
        }

        pagemap().clear_slab(slab);
        large_allocator.dealloc(slab, 0);
        stats().superslab_push();
      }
      else if (was_full)
      {
        sizeclass_t medium_class = sizeclass - NUM_SMALL_CLASSES;
        DLList<Mediumslab>* sc = &medium_classes[medium_class];
        sc->insert(slab);
      }
    }

    template<ZeroMem zero_mem, AllowReserve allow_reserve>
    void* large_alloc(size_t size)
    {
      MEASURE_TIME_MARKERS(
        large_alloc,
        4,
        16,
        MARKERS(
          zero_mem == YesZero ? "zeromem" : "nozeromem",
          allow_reserve == NoReserve ? "noreserve" : "reserve"));

      if (void* replacement = Replacement(this))
      {
        return reinterpret_cast<Allocator*>(replacement)
          ->template large_alloc<zero_mem, allow_reserve>(size);
      }

      size_t size_bits = bits::next_pow2_bits(size);
      size_t large_class = size_bits - SUPERSLAB_BITS;
      assert(large_class < NUM_LARGE_CLASSES);

      void* p = large_allocator.template alloc<zero_mem, allow_reserve>(
        large_class, size);

      pagemap().set_large_size(p, size);

      stats().large_alloc(large_class);
      return p;
    }

    void large_dealloc(void* p, size_t size)
    {
      MEASURE_TIME(large_dealloc, 4, 16);

      size_t size_bits = bits::next_pow2_bits(size);
      size_t rsize = bits::one_at_bit(size_bits);
      assert(rsize >= SUPERSLAB_SIZE);
      size_t large_class = size_bits - SUPERSLAB_BITS;

      pagemap().clear_large_size(p, size);

      stats().large_dealloc(large_class);

      if ((decommit_strategy != DecommitNone) || (large_class > 0))
        large_allocator.memory_provider.notify_not_using(
          pointer_offset(p, OS_PAGE_SIZE), rsize - OS_PAGE_SIZE);

      // Initialise in order to set the correct SlabKind.
      Largeslab* slab = static_cast<Largeslab*>(p);
      slab->init();
      large_allocator.dealloc(slab, large_class);
    }

    // Note that this is on the slow path as it lead to better code.
    // As it is tail, not inlining means that it is jumped to, so has no perf
    // impact on the producer consumer scenarios, and doesn't require register
    // spills in the fast path for local deallocation.
    SNMALLOC_SLOW_PATH
    void remote_dealloc(RemoteAllocator* target, void* p, sizeclass_t sizeclass)
    {
      MEASURE_TIME(remote_dealloc, 4, 16);
      assert(target->id() != id());

#if SNMALLOC_QUARANTINE_DEALLOC == 0
      handle_message_queue();
#endif

      void* offseted = apply_cache_friendly_offset(p, sizeclass);

      // Check whether this will overflow the cache first.  If we are a fake
      // allocator, then our cache will always be full and so we will never hit
      // this path.
      size_t sz = sizeclass_to_size(sizeclass);
      if ((remote.size + sz) < REMOTE_CACHE)
      {
        stats().remote_free(sizeclass);
        remote.dealloc_sized(target->id(), offseted, sz);
        return;
      }
      // Now that we've established that we're in the slow path (if we're a
      // real allocator, we will have to empty our cache now), check if we are
      // a real allocator and construct one if we aren't.
      if (void* replacement = Replacement(this))
      {
        // We have to do a dealloc, not a remote_dealloc here because this may
        // have been allocated with the allocator that we've just had returned.
        reinterpret_cast<Allocator*>(replacement)->dealloc(p);
        return;
      }

      stats().remote_free(sizeclass);
      remote.dealloc(target->id(), offseted, sizeclass);

      stats().remote_post();
      remote.post(id());
    }

    PageMap& pagemap()
    {
      return page_map;
    }
  };
} // namespace snmalloc
