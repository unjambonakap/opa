/* Copyright 2015 The TensorFlow Authors. All Rights Reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
==============================================================================*/

#include "glib/core/blocking_counter.h"
#include "glib/core/threadpool.h"
#include "glib/platform/test.h"
#include "glib/platform/test_benchmark.h"

namespace glib {
namespace {

TEST(BlockingCounterTest, TestZero) {
  BlockingCounter bc(0);
  bc.Wait();
}

TEST(BlockingCounterTest, TestSingleThread) {
  BlockingCounter bc(2);
  bc.DecrementCount();
  bc.DecrementCount();
  bc.Wait();
}

TEST(BlockingCounterTest, TestMultipleThread) {
  int N = 3;
  thread::ThreadPool* thread_pool =
      new thread::ThreadPool(Env::Default(), "test", N);

  BlockingCounter bc(N);
  for (int i = 0; i < N; ++i) {
    thread_pool->Schedule([&bc] { bc.DecrementCount(); });
  }

  bc.Wait();
  delete thread_pool;
}

}  // namespace

static void BM_BlockingCounter(int iters, int num_threads,
                               int shards_per_thread) {
  testing::StopTiming();
  std::unique_ptr<thread::ThreadPool> thread_pool(
      new thread::ThreadPool(Env::Default(), "test", num_threads));
  const int num_shards = num_threads * shards_per_thread;
  testing::StartTiming();
  for (int i = 0; i < iters; ++i) {
    BlockingCounter bc(num_shards);
    for (int j = 0; j < num_threads; ++j) {
      thread_pool->Schedule([&bc, shards_per_thread] {
        for (int k = 0; k < shards_per_thread; ++k) {
          bc.DecrementCount();
        }
      });
    }
    bc.Wait();
  }
  testing::StopTiming();
}

BENCHMARK(BM_BlockingCounter)->RangePair(1, 12, 1, 1000);

}  // namespace glib
