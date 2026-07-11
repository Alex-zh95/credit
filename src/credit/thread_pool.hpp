#ifndef THREAD_POOL_HPP
#define THREAD_POOL_HPP

#include <condition_variable>
#include <cstddef>
#include <functional>
#include <future>
#include <memory>
#include <mutex>
#include <queue>
#include <thread>
#include <utility>
#include <vector>

/* Description
 * -----------
 * Minimal fixed-size thread pool.
 *
 * Why: NLopt's derivative-free optimizers evaluate the objective function
 * thousands of times per fit. Spawning hardware_concurrency() threads inside
 * every evaluation (the previous design) pays thread creation/teardown cost on
 * each of those calls. Here the workers are created once, live for the duration
 * of the fit, and repeatedly pick jobs off a queue.
 *
 * submit() returns a std::future<void> so callers can (a) wait for a batch of
 * jobs to complete and (b) observe any exception a job threw - with raw
 * std::thread an escaping exception would call std::terminate instead.
 */
class ThreadPool {
  public:
    explicit ThreadPool(std::size_t n_threads) {
        workers.reserve(n_threads);
        for (std::size_t i = 0; i < n_threads; ++i)
            workers.emplace_back([this] { worker_loop(); });
    }

    ~ThreadPool() {
        {
            std::lock_guard lock(mtx);
            stopping = true;
        }
        cv.notify_all();
        for (auto& worker : workers) worker.join();
    }

    ThreadPool(const ThreadPool&) = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;

    template <typename F>
    std::future<void> submit(F&& job) {
        // packaged_task is move-only but std::function requires copyable
        // callables, hence the shared_ptr wrapper
        auto task = std::make_shared<std::packaged_task<void()>>(std::forward<F>(job));
        auto result = task->get_future();
        {
            std::lock_guard lock(mtx);
            jobs.emplace([task] { (*task)(); });
        }
        cv.notify_one();
        return result;
    }

  private:
    void worker_loop() {
        for (;;) {
            std::function<void()> job;
            {
                std::unique_lock lock(mtx);
                cv.wait(lock, [this] { return stopping || !jobs.empty(); });
                if (stopping && jobs.empty()) return;
                job = std::move(jobs.front());
                jobs.pop();
            }
            job();
        }
    }

    std::vector<std::thread> workers;
    std::queue<std::function<void()>> jobs;
    std::mutex mtx;
    std::condition_variable cv;
    bool stopping = false;
};

#endif // THREAD_POOL_HPP
