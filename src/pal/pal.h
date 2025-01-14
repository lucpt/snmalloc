#pragma once

#include "pal_consts.h"

namespace snmalloc
{
  void error(const char* const str);
} // namespace snmalloc

// If simultating OE, then we need the underlying platform
#if !defined(OPEN_ENCLAVE) || defined(OPEN_ENCLAVE_SIMULATION)
#  include "pal_apple.h"
#  include "pal_cheribsd.h"
#  include "pal_free_bsd_kernel.h"
#  include "pal_freebsd.h"
#  include "pal_linux.h"
#  include "pal_windows.h"
#endif
#include "pal_open_enclave.h"
#include "pal_plain.h"

namespace snmalloc
{
#if !defined(OPEN_ENCLAVE) || defined(OPEN_ENCLAVE_SIMULATION)
  using DefaultPal =
#  if defined(_WIN32)
    PALWindows;
#  elif defined(__APPLE__)
    PALApple;
#  elif defined(__linux__)
    PALLinux;
#  elif defined(FreeBSD_KERNEL)
    PALFreeBSDKernel;
#  elif defined(__FreeBSD__)
#    if defined(__CHERI_PURE_CAPABILITY__)
    PALCHERIBSD;
#    else
    PALFBSD;
#    endif
#  else
#    error Unsupported platform
#  endif
#endif

  using Pal =
#if defined(SNMALLOC_MEMORY_PROVIDER)
    PALPlainMixin<SNMALLOC_MEMORY_PROVIDER>;
#elif defined(OPEN_ENCLAVE)
    PALPlainMixin<PALOpenEnclave>;
#else
    DefaultPal;
#endif

  inline void error(const char* const str)
  {
    Pal::error(str);
  }

  /**
   * Query whether the PAL supports a specific feature.
   */
  template<PalFeatures F, typename PAL = Pal>
  constexpr static bool pal_supports()
  {
    return (PAL::pal_features & F) == F;
  }
} // namespace snmalloc
