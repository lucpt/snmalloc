/*
 * This is a somewhat, OK, very, artificial benchmark, intended to highlight
 * the potential gains of segmented quarantines in snmalloc+revocation.
 *
 * The basic design is as follows.  We create a series of threads which hand
 * malloc()'d messages around and free() them; the goal is that we should be
 * able to control the relative rates at which threads free(), so that we
 * can, indeed, see that the fastest thread is usually the one doing
 * revocation sweeps and that the rest of the system continues to function.
 *
 */

#include <mutex>
#include <thread>
#include <snmalloc.h>
#include <test/measuretime.h>

using namespace snmalloc;
using namespace std::chrono_literals;

static const int SCALE_FACTOR = 6; // Pass more messages?
static const int NTHR = 4;

class Msg
{
public:
  std::atomic<Msg*> next;
  int destbits;
};

std::mutex diag_mtx; // serialize printouts
std::atomic<bool> shutdown;
std::atomic<size_t> global_msgs;
MPSCQ<Msg> queues[NTHR];

int distribute(Msg* m)
{
  int ix, ixb;
  for (ix = 0, ixb = 1; ix < NTHR; ix++, ixb <<= 1)
  {
    if (m->destbits & ixb) {
      m->destbits &= ~ixb;
      queues[ix].enqueue(m, m);
      return 1;
    }
  }
  return 0;
}

int make(int dest)
{
  Msg* m = new Msg();
  m->destbits = (dest & ((1 << NTHR) - 1));
  if (m->destbits == 0)
    m->destbits = ((1 << NTHR) - 1);
  if (!distribute(m))
  {
    delete m;
    return 0;
  }
  else
  {
    return 1;
  }
}

void worker(int ix)
{
  size_t c_yield = 0;
  size_t c_msg = 0;
  size_t c_new = 0;
  size_t c_del = 0;

  while(true)
  {
    auto mp = queues[ix].dequeue();
    if (!mp.second)
    {
      c_yield++;
      std::this_thread::sleep_for(1ms);
      if (shutdown)
        break;
      else
        continue;
    }
    Msg* m = mp.first;

    /* With some probability, distribute a new message */
    if ((c_msg % (1 << (ix+2))) == 0)
    {
      global_msgs++;
      make(~m->destbits);
      c_new++;
    }

    /* Distribute the message dequeued */
    if (!distribute(m))
    {
      size_t om = global_msgs--;
      if (om == 0)
        abort();
      c_del++;
      delete m;
    }

    c_msg++;
  }

  {
    std::lock_guard<std::mutex> guard(diag_mtx);
    std::cout << "Worker ix=" << ix
              << " a="        << ThreadAlloc::get()
              << " c_yield="  << c_yield
              << " c_msg="    << c_msg
              << " c_new="    << c_new
              << " c_del="    << c_del
              << std::endl;
  }
}

int main(int, char **)
{
  int ix;
  std::thread* thrs[NTHR];

  shutdown = false;
  global_msgs = 0;
  for (ix = 0; ix < NTHR; ix++)
  {
    global_msgs++;
    queues[ix].init(new Msg());
  }

  const int initial_batch = 1 << (10 + SCALE_FACTOR);
  global_msgs += initial_batch;
  for (ix = 0; ix < initial_batch; ix++)
  {
    make(ix);
  }
  
  for (ix = 0; ix < NTHR; ix++)
  {
    thrs[ix] = new std::thread(worker, ix);
  }

  for (ix = 0; ix < (1<<8); ix++)
  {
    const int wait_thresh = 1 << (10 + SCALE_FACTOR);
    while(global_msgs > wait_thresh)
    {
      std::this_thread::sleep_for(1ms);
    }

    const int batch = 1 << (8 + SCALE_FACTOR);
    global_msgs += batch;
    for (size_t jx = 0; jx < batch; jx++)
    {
      make(jx);
    }
  }

  while(global_msgs > (1<<10))
  {
    ThreadAlloc::get()->handle_message_queue();
    std::this_thread::sleep_for(1ms);
  }

  shutdown = true;
  for (ix = 1; ix < NTHR; ix++)
  {
    thrs[ix]->join();
    delete thrs[ix];
  }

  return 0;
}
