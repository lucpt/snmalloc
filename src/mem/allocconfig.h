#pragma once

#include "../ds/bits.h"

namespace snmalloc
{
#ifndef SNMALLOC_CHERI_ALIGN
#  define SNMALLOC_CHERI_ALIGN 0
#endif

#ifndef SNMALLOC_PAGEMAP_POINTERS
#  define SNMALLOC_PAGEMAP_POINTERS 0
#endif

#ifndef SNMALLOC_PAGEMAP_REDERIVE
#  define SNMALLOC_PAGEMAP_REDERIVE 0
#endif

#if (SNMALLOC_PAGEMAP_REDERIVE == 1) && (SNMALLOC_PAGEMAP_POINTERS == 0)
#  error Need pointers in pagemap for rederivation
#endif

#ifndef SNMALLOC_CHERI_SETBOUNDS
#  define SNMALLOC_CHERI_SETBOUNDS 0
#endif

#if SNMALLOC_CHERI_SETBOUNDS == 1
#  if SNMALLOC_PAGEMAP_REDERIVE == 0
#    error Must rederive pointers if setting bounds
#  endif
#  ifndef __CHERI__
#    error Unable to set bounds without CHERI
#  endif
#endif

#if (SNMALLOC_CHERI_SETBOUNDS == 1) && (SNMALLOC_CHERI_ALIGN == 0)
#  error CHERI cannot safely bound objects with misaligned sizes
#endif

#ifndef SNMALLOC_QUARANTINE_DEALLOC
#  define SNMALLOC_QUARANTINE_DEALLOC 0
#endif

#if SNMALLOC_QUARANTINE_DEALLOC == 1
/* Quarantine policy knob(s). */

/*
 * Each allocator is only willing to pin some footprint in quarantine, but
 * note that allocators are per-thread
 */
#  ifndef SNMALLOC_QUARANTINE_PER_ALLOC_THRESHOLD
#    define SNMALLOC_QUARANTINE_PER_ALLOC_THRESHOLD (32 * 1024 * 1024)
#  endif

/*
 * Optionally, prevent the quarantine queue from growing too big, even
 * if not pinning too much physmem.
 */
#  ifdef SNMALLOC_QUARANTINE_PER_ALLOC_CHUNK_THRESHOLD
#    if SNMALLOC_QUARANTINE_PER_ALLOC_CHUNK_THRESHOLD < 0
#      error Quarantine must be allowed at least one chunk
#    endif
#  endif

/*
 * The size of the "chunk"s holding quarantine might matter, so make it a
 * knob here.  This is a sizeclass parameter and has only been tested with
 * medium classes, * somewhere between NUM_SMALL_CLASSES and
 * NUM_SIZECLASSES-1.
 *
 * For ease of refering to the endpoints of the spectrum from outside the
 * allocator's C source we also map some implausible numbers to those, and
 * reserve the right to extend this table:
 *
 *   -1024 : NUM_SMALL_CLASSES
 *   -2048 : NUM_SIZECLASSES-1
 */
#  if !defined(SNMALLOC_QUARANTINE_CHUNK_SIZECLASS)
#    define SNMALLOC_QUARANTINE_CHUNK_SIZECLASS (NUM_SMALL_CLASSES)
#  elif SNMALLOC_QUARANTINE_CHUNK_SIZECLASS == -1024
#    undef SNMALLOC_QUARANTINE_CHUNK_SIZECLASS
#    define SNMALLOC_QUARANTINE_CHUNK_SIZECLASS (NUM_SMALL_CLASSES)
#  elif SNMALLOC_QUARANTINE_CHUNK_SIZECLASS == -2048
#    undef SNMALLOC_QUARANTINE_CHUNK_SIZECLASS
#    define SNMALLOC_QUARANTINE_CHUNK_SIZECLASS (NUM_SIZECLASSES-1)
#  endif

#  if !defined(SNMALLOC_QUARANTINE_POLICY_STRICT)
#    define SNMALLOC_QUARANTINE_POLICY_STRICT 1
#  endif

#endif

#ifndef SNMALLOC_QUARANTINE_CHATTY
#  define SNMALLOC_QUARANTINE_CHATTY 0
#endif

#ifndef SNMALLOC_REVOKE_QUARANTINE
#  define SNMALLOC_REVOKE_QUARANTINE 0
#endif

#if SNMALLOC_REVOKE_QUARANTINE == 1

#  ifndef SNMALLOC_REVOKE_PARANOIA
#    define SNMALLOC_REVOKE_PARANOIA 0
#  endif

#  ifndef SNMALLOC_REVOKE_DRY_RUN
#    define SNMALLOC_REVOKE_DRY_RUN 0
#  elif (SNMALLOC_REVOKE_PARANOIA == 1) && (SNMALLOC_REVOKE_DRY_RUN == 1)
#    error Doing nothing while being paranoid will not work out well.
#  endif

#  if (SNMALLOC_REVOKE_DRY_RUN == 1) && (SNMALLOC_QUARANTINE_CHATTY == 1)
#    error Refusing to raise much ado about nothing
#  endif

#  if SNMALLOC_QUARANTINE_DEALLOC == 0
#    error Revocation depends upon quarantine.
#  endif
#  if SNMALLOC_CHERI_SETBOUNDS == 0
#    error Bounds information is used in revocation; this will not work.
#  endif
#endif

/*
 * Deallocation routines which take a size from the caller are unsafe in
 * that they can introduce type confusion within the allocator: we could be
 * convinced to treat user memory as a superslab, for example.  Set
 * SNMALLOC_UNSAFE_FREES to 0 to stub them out and redirect them to the
 * generic deallocation routine which just takes a pointer.  Instead, one
 * can set SNMALLOC_UNSAFE_FREES_CHECK (with CHECK_CLIENT) to verify that
 * the size indicated is compatible with the actual pointer size obtained by
 * metadata lookup.
 */
#ifndef SNMALLOC_UNSAFE_FREES
#  define SNMALLOC_UNSAFE_FREES 1
#endif

// The CHECK_CLIENT macro is used to turn on minimal checking of the client
// calling the API correctly.
#if !defined(NDEBUG) && !defined(CHECK_CLIENT)
#  define CHECK_CLIENT
#endif

#if (SNMALLOC_UNSAFE_FREES == 1) && !defined(SNMALLOC_UNSAFE_FREES_CHECK)
#  ifdef CHECK_CLIENT
#    define SNMALLOC_UNSAFE_FREES_CHECK 1
#  else
#    define SNMALLOC_UNSAFE_FREES_CHECK 0
#  endif
#endif

  // 0 intermediate bits results in power of 2 small allocs. 1 intermediate
  // bit gives additional sizeclasses at the midpoint between each power of 2.
  // 2 intermediate bits gives 3 intermediate sizeclasses, etc.
  static constexpr size_t INTERMEDIATE_BITS =
#ifdef USE_INTERMEDIATE_BITS
    USE_INTERMEDIATE_BITS
#else
    2
#endif
    ;

  // Return remote small allocs when the local cache reaches this size.
  static constexpr size_t REMOTE_CACHE =
#ifdef USE_REMOTE_CACHE
    USE_REMOTE_CACHE
#else
    1 << 20
#endif
    ;

  // Handle at most this many object from the remote dealloc queue at a time.
  static constexpr size_t REMOTE_BATCH =
#ifdef USE_REMOTE_BATCH
    REMOTE_BATCH
#else
    64
#endif
    ;

  // Specifies smaller slab and super slab sizes for address space
  // constrained scenarios.
  static constexpr size_t ADDRESS_SPACE_CONSTRAINED =
#ifdef IS_ADDRESS_SPACE_CONSTRAINED
    true
#else
    // In 32 bit uses smaller superslab.
    (!bits::is64())
#endif
    ;

  static constexpr size_t RESERVE_MULTIPLE =
#ifdef USE_RESERVE_MULTIPLE
    USE_RESERVE_MULTIPLE
#else
    bits::is64() ? 16 : 2
#endif
    ;

  enum DecommitStrategy
  {
    /**
     * Never decommit memory.
     */
    DecommitNone,
    /**
     * Decommit superslabs when they are entirely empty.
     */
    DecommitSuper,
    /**
     * Decommit all slabs once they are empty.
     */
    DecommitAll,
    /**
     * Decommit superslabs only when we are informed of memory pressure by the
     * OS, do not decommit anything in normal operation.
     */
    DecommitSuperLazy,
    /**
     * Decommit medium and large objects on entry to quarantine
     */
    DecommitQuarantine,
    /**
     * Decommit medium and large objects on entry to quarantine
     * and superslabs when they become empty.
     */
    DecommitQuarantineSuper,
    /* XXX DecommitQuarantineSuperLazy, if ever quarantine on Windows? */

    /* XXX DecommitQuarantineMP, DecommitQuarantineSuperMP: if we get
     * MPROT_QUARANTINE landed, we need to know not do to anything when
     * reissuing.
     */
  };

  static constexpr DecommitStrategy decommit_strategy =
#ifdef USE_DECOMMIT_STRATEGY
    USE_DECOMMIT_STRATEGY
#elif defined(_WIN32) && !defined(OPEN_ENCLAVE)
    DecommitSuperLazy
#elif SNMALLOC_QUARANTINE_DEALLOC == 1
    DecommitQuarantineSuper
#else
    DecommitSuper
#endif
    ;

  // The DecommitQuarantine* options only make sense if quarantining.
  static_assert(
    ((decommit_strategy != DecommitQuarantine) &&
     (decommit_strategy != DecommitQuarantineSuper)) ||
    (SNMALLOC_QUARANTINE_DEALLOC == 1));

  // The remaining values are derived, not configurable.
  static constexpr size_t POINTER_BITS =
    bits::next_pow2_bits_const(sizeof(uintptr_t));

  // Used to isolate values on cache lines to prevent false sharing.
  static constexpr size_t CACHELINE_SIZE = 64;

  // Used to keep Superslab metadata committed.
  static constexpr size_t OS_PAGE_SIZE = 0x1000;
  static constexpr size_t PAGE_ALIGNED_SIZE = OS_PAGE_SIZE << INTERMEDIATE_BITS;
  // Some system headers (e.g. Linux' sys/user.h, FreeBSD's machine/param.h)
  // define `PAGE_SIZE` as a macro.  We don't use `PAGE_SIZE` as our variable
  // name, to avoid conflicts, but if we do see a macro definition then check
  // that our value matches the platform's expected value.
#ifdef PAGE_SIZE
  static_assert(
    PAGE_SIZE == OS_PAGE_SIZE,
    "Page size from system header does not match snmalloc config page size.");
#endif

  // Minimum allocation size is space for two pointers.
#ifdef __CHERI_PURE_CAPABILITY__
  static constexpr size_t MIN_ALLOC_BITS = 5;
#else
  static constexpr size_t MIN_ALLOC_BITS = bits::is64() ? 4 : 3;
#endif
  static constexpr size_t MIN_ALLOC_SIZE = 1 << MIN_ALLOC_BITS;

  // Slabs are 64 KiB unless constrained to 16 KiB.
  static constexpr size_t SLAB_BITS = ADDRESS_SPACE_CONSTRAINED ? 14 : 16;
  static constexpr size_t SLAB_SIZE = 1 << SLAB_BITS;
  static constexpr size_t SLAB_MASK = ~(SLAB_SIZE - 1);

  // Superslabs are composed of this many slabs. Slab offsets are encoded as
  // a byte, so the maximum count is 256. This must be a power of two to
  // allow fast masking to find a superslab start address.
  static constexpr size_t SLAB_COUNT_BITS = ADDRESS_SPACE_CONSTRAINED ? 6 : 8;
  static constexpr size_t SLAB_COUNT = 1 << SLAB_COUNT_BITS;
  static constexpr size_t SUPERSLAB_SIZE = SLAB_SIZE * SLAB_COUNT;
  static constexpr size_t SUPERSLAB_MASK = ~(SUPERSLAB_SIZE - 1);
  static constexpr size_t SUPERSLAB_BITS = SLAB_BITS + SLAB_COUNT_BITS;
  static constexpr size_t RESERVE_SIZE = SUPERSLAB_SIZE * RESERVE_MULTIPLE;

  static_assert((1ULL << SUPERSLAB_BITS) == SUPERSLAB_SIZE, "Sanity check");

  // Number of slots for remote deallocation.
  static constexpr size_t REMOTE_SLOT_BITS = 6;
  static constexpr size_t REMOTE_SLOTS = 1 << REMOTE_SLOT_BITS;
  static constexpr size_t REMOTE_MASK = REMOTE_SLOTS - 1;

  static_assert(
    INTERMEDIATE_BITS < MIN_ALLOC_BITS,
    "INTERMEDIATE_BITS must be less than MIN_ALLOC_BITS");
  static_assert(
    MIN_ALLOC_SIZE >= (sizeof(void*) * 2),
    "MIN_ALLOC_SIZE must be sufficient for two pointers");
  static_assert(
    SLAB_BITS <= (sizeof(uint16_t) * 8),
    "SLAB_BITS must not be more than the bits in a uint16_t");
  static_assert(
    SLAB_COUNT == bits::next_pow2_const(SLAB_COUNT),
    "SLAB_COUNT must be a power of 2");
  static_assert(
    SLAB_COUNT <= (UINT8_MAX + 1), "SLAB_COUNT must fit in a uint8_t");
} // namespace snmalloc
