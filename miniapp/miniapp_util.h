#pragma once

#include "SNLS_config.h"
#if defined(SNLS_RAJA_PORT_SUITE)
#include "SNLS_memory_manager.h"
#endif

#include "SNLS_unused.h"

// We're going to use this to determine what RAJA code to run for our
// kernels.
// The HIP backend won't be able to run on AMD GPGPUs
// until device function pointers are supported.
enum class ExecutionStrategy { CPU, GPU, OPENMP };

template<class T>
class memoryManager {
public:
    memoryManager() = delete;

    memoryManager(const size_t num_items) : total_items(num_items) {
        assert(num_items > 0 && "num_items must be greater than 0...");
#if defined(SNLS_RAJA_PORT_SUITE)
        auto mm = snls::memoryManager::getInstance();
        buffer = mm.allocManagedArray<T>(num_items);
#else
    buffer = new T[num_items];
#endif
    }

    ~memoryManager() {
#if defined(SNLS_RAJA_PORT_SUITE)
        buffer.free();
#else
        if (buffer) {
            delete buffer;
        }
#endif
    }

    T* getNew(const size_t num_items, const ecmech::ExecutionStrategy UNUSED_GPU(strat)) {
        assert((num_items + offset) <= total_items && "Requested too large of an allocation");
        const size_t old_offset = offset;
        offset += num_items;
#if defined(SNLS_RAJA_PORT_SUITE)
        chai::ExecutionSpace ses;
        switch (strat) {
            case ecmech::ExecutionStrategy::GPU: {
                ses = chai::ExecutionSpace::GPU;
                break;
            }
            case ecmech::ExecutionStrategy::OPENMP:
            case ecmech::ExecutionStrategy::CPU:
            default: {
                ses = chai::ExecutionSpace::CPU;
                break;
            }
        }
        return &(buffer.data(ses)[old_offset]);
#else
        return &(buffer[old_offset]);
#endif
    }

private:
#if defined(SNLS_RAJA_PORT_SUITE)
    chai::ManagedArray<T> buffer;
#else
    T* buffer = nullptr;
#endif
    size_t offset = 0;
    const size_t total_items;
};

