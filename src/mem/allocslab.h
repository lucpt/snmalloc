#pragma once

#include "../mem/baseslab.h"
#include "remoteallocator.h"

namespace snmalloc
{
  class Allocslab : public Baseslab
  {
  protected:
    RemoteAllocator* allocator;

#if SNMALLOC_REVOKE_QUARANTINE == 1
    uint8_t* revbitmap;
#endif

  public:
    RemoteAllocator* get_allocator()
    {
      return allocator;
    }

#if SNMALLOC_REVOKE_QUARANTINE == 1
    uint8_t* get_revbitmap()
    {
      return revbitmap;
    }
#endif
  };
} // namespace snmalloc
