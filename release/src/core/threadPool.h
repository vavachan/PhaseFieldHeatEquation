/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Adrian Kummerlaender
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
 */

#ifndef THREAD_POOL_H
#define THREAD_POOL_H

#include <thread>
#include <mutex>
#include <atomic>
#include <condition_variable>
#include <queue>
#include <functional>
#include <future>
#include <vector>

#include "core/olbDebug.h"
#include "io/ostreamManager.h"

namespace olb {

/// Pool of threads for CPU-based background processing
/**
 * Not intended for parallelization of performance critical paths but
 * for background processing of e.g. VTK output to reduce blocking of
 * simulation process.
 **/
#ifdef FEATURE_EMSCRIPTEN
class ThreadPool { // Dummy implementation for Emscripten
public:
  void init(int nThreads, bool verbose) {}

  unsigned size() const { return 0; }

  template <typename F>
  void scheduleAndForget(F&& f) {
    throw std::runtime_error("ThreadPool::scheduleAndForget not implemented for Emscripten");
    f();
  }

  template <typename F, typename R = std::invoke_result_t<std::decay_t<F>>>
  std::future<R> schedule(F&& f) {
    throw std::runtime_error("ThreadPool::schedule not implemented for Emscripten");
    return std::async(std::launch::async, f);
  }

  void wait() {}

  template <typename T>
  void waitFor(std::vector<std::future<T>>& futures) {
    for (auto& future : futures) {
      future.wait();
    }
  }
};

#else

class ThreadPool {
private:
  void work(unsigned iThread)
  {
    while (_active) {
      std::unique_lock lock(_mutex);
      _available.wait(lock, [&]() {
        return !_queue.empty()
            || !_active;
      });
      if (_active) {
        std::function<void()> task = std::move(_queue.front());
        _queue.pop();
        lock.unlock();
        task();
        --_taskCount;
        if (_waiting) {
          _done.notify_one();
        }
      }
    }
  }

  mutable std::mutex _mutex;

  std::atomic<bool> _waiting;
  std::atomic<bool> _active;

  std::atomic<std::size_t> _taskCount;

  std::condition_variable _available;
  std::condition_variable _done;

  std::vector<std::thread> _threads;
  std::queue<std::function<void()>> _queue;

  bool _initialized;

public:
  ThreadPool():
    _mutex{},
    _waiting{false},
    _active{true},
    _taskCount{0},
    _available{},
    _done{},
    _threads{1},
    _queue{},
    _initialized{false}
  { }

  /// Initialization to be called by olbInit
  void init(int nThreads, bool verbose)
  {
    OLB_PRECONDITION(!_initialized);
    OstreamManager clout(std::cout, "ThreadPool");
    if (nThreads > 1) {
      _threads.resize(nThreads);
    }
    for (unsigned iThread=0; iThread < _threads.size(); ++iThread) {
      _threads[iThread] = std::thread(&ThreadPool::work, this, iThread);
    }
    if (verbose) {
      clout << "Sucessfully initialized, numThreads=" << std::to_string(_threads.size()) << std::endl;
    }
    _initialized = true;
  }

  ~ThreadPool()
  {
    _active = false;
    _available.notify_all();
    for (std::thread& thread : _threads) {
      thread.join();
    }
  }

  /// Returns number of threads
  unsigned size() const
  {
    return _threads.size();
  }

  /// Schedule F, tracking neither its return value nor completion
  template <typename F>
  void scheduleAndForget(F&& f)
  {
    OLB_PRECONDITION(_initialized);
    {
      const std::scoped_lock lock(_mutex);
      _queue.push(f);
    }
    ++_taskCount;
    _available.notify_one();
  }

  /// Schedule F and return future of its return value
  template <typename F, typename R = std::invoke_result_t<std::decay_t<F>>>
  std::future<R> schedule(F&& f)
  {
    OLB_PRECONDITION(_initialized);
    auto packagedF = std::make_shared<std::packaged_task<R()>>(f);
    {
      const std::scoped_lock lock(_mutex);
      _queue.push([packagedF]() {
        (*packagedF)();
      });
    }
    ++_taskCount;
    _available.notify_one();
    return packagedF->get_future();
  }

  /// Blocks until all tasks are completed
  void wait()
  {
    OLB_PRECONDITION(_initialized);
    _waiting = true;
    std::unique_lock lock(_mutex);
    _done.wait(lock, [&]() {
      return _taskCount == 0;
    });
    _waiting = false;
  }

  /// Blocks until all tasks producing the given futures are completed
  template <typename T>
  void waitFor(std::vector<std::future<T>>& futures)
  {
    for (auto& future : futures) {
      future.wait();
    }
  }

};
#endif

}

#endif
