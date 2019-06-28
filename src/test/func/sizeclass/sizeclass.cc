#include <iostream>
#include <snmalloc.h>

NOINLINE
snmalloc::sizeclass_t size_to_sizeclass(size_t size)
{
  return snmalloc::size_to_sizeclass(size);
}

int main(int, char**)
{
  bool failed = false;
  size_t size_low = 0;
  size_t zsc = (size_t)snmalloc::size_to_sizeclass(0);

  std::cout << "INTERMEDIATE: " << snmalloc::INTERMEDIATE_BITS
            << " MIN_ALLOC: " << snmalloc::MIN_ALLOC_BITS << std::endl;

  std::cout << "0 has sizeclass: " << zsc << std::endl;

  std::cout << "size(sizeclass(0)) = "
            << (size_t)snmalloc::sizeclass_to_size(zsc) << std::endl;

  std::cout << "sizeclass |-> [size_low, size_high] " << std::endl;

  for (snmalloc::sizeclass_t sz = 0; sz < snmalloc::NUM_SIZECLASSES; sz++)
  {
    // Separate printing for small and medium sizeclasses
    if (sz == snmalloc::NUM_SMALL_CLASSES)
      std::cout << std::endl;

    size_t size = snmalloc::sizeclass_to_size(sz);
    std::cout << (size_t)sz << " |-> "
              << "[" << size_low + 1 << ", " << size << "]";

#ifdef __CHERI_PURE_CAPABILITY__
    size_t cheri_low = snmalloc::bits::align_up(
      size_low + 1, 1 << CHERI_ALIGN_SHIFT(size_low + 1));
    size_t cheri_size =
      snmalloc::bits::align_up(size, 1 << CHERI_ALIGN_SHIFT(size));
    std::cout << " CHERI [" << cheri_low << ", " << cheri_size << "] "
              << (cheri_size == size ? "ok" : "bad");
#endif

    std::cout << std::endl;

    if (size < size_low)
    {
      std::cout << "Sizeclass " << (size_t)sz << " is " << size
                << " which is less than " << size_low << std::endl;
      failed = true;
    }

    for (size_t i = size_low + 1; i <= size; i++)
    {
      if (size_to_sizeclass(i) != sz)
      {
        std::cout << "Size " << i << " has sizeclass "
                  << (size_t)size_to_sizeclass(i) << " but expected sizeclass "
                  << (size_t)sz << std::endl;
        failed = true;
      }
    }

    size_low = size;
  }

  if (failed)
    abort();
}