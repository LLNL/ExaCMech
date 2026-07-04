/**
 * @file miniapp_util.h
 *
 * @brief Declares `memoryManager<T>`, the miniapp's single-allocation sub-buffer
 * carver: one big backing array is allocated up front and every per-quantity array used
 * in `orientation_evolution.cxx`'s main loop (`state_vars`, `velocity_grad`,
 * `cauchy_stress_array`, ...) is handed out as a non-overlapping slice of it via
 * `getNew`, so the whole simulation's working set lives in one CHAI-managed (or plain
 * host) allocation that is freed automatically when the `memoryManager` goes out of
 * scope.
 */

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
//
// @note This local `ExecutionStrategy` is currently unused/dead: every call site in the
// miniapp (orientation_evolution.cxx, and memoryManager::getNew below) uses
// `ecmech::ExecutionStrategy` from the core library instead. Left in place rather than
// removed as part of this documentation-only pass -- see DOCUMENTATION_TODO.md.
enum class ExecutionStrategy { CPU, GPU, OPENMP };

/**
 * @brief Bump-allocator view over one big backing array, handing out non-overlapping
 * `T*` sub-buffers by size.
 *
 * Construct once with the total number of `T` elements every sub-buffer will need
 * combined (`orientation_evolution.cxx` computes this as `nqpts * (num_state_vars +
 * num_var_variables)`), then call `getNew()` repeatedly to carve off each named array in
 * turn. There's no way to free an individual sub-buffer -- the whole backing allocation
 * is freed together when the `memoryManager` itself is destroyed.
 *
 * @tparam T Element type of the backing array (this miniapp always uses `double`).
 */
template<class T>
class memoryManager {
public:
    memoryManager() = delete;

    /**
     * @brief Allocate the single backing array.
     * @param num_items Total number of `T` elements needed across every sub-buffer that
     * will be requested via `getNew()`.
     */
    memoryManager(const size_t num_items) : total_items(num_items) {
        assert(num_items > 0 && "num_items must be greater than 0...");
#if defined(SNLS_RAJA_PORT_SUITE)
        auto mm = snls::memoryManager::getInstance();
        buffer = mm.allocManagedArray<T>(num_items);
#else
    buffer = new T[num_items];
#endif
    }

    /** @brief Frees the backing array (all sub-buffers handed out by `getNew` become invalid). */
    ~memoryManager() {
#if defined(SNLS_RAJA_PORT_SUITE)
        buffer.free();
#else
        if (buffer) {
            delete buffer;
        }
#endif
    }

    /**
     * @brief Carve off the next `num_items`-sized slice of the backing array.
     *
     * Slices are handed out contiguously in call order (a simple bump allocator, no
     * reuse/freeing of individual slices) -- so callers must request them in the same
     * order every time and must not request more total elements than were reserved in
     * the constructor (checked by the `assert` below).
     *
     * @param num_items Number of `T` elements this slice should contain.
     * @param strat Execution strategy the returned pointer will be used under; only
     * meaningful in the CHAI-backed (`SNLS_RAJA_PORT_SUITE`) build, where it selects
     * which execution space's copy of the data is returned (host build ignores it, hence
     * `UNUSED_GPU`).
     * @return Pointer to the start of the newly carved slice.
     */
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
    /** @brief Backing storage (CHAI build): one managed array shared by every sub-buffer. */
    chai::ManagedArray<T> buffer;
#else
    /** @brief Backing storage (plain host build): one heap array shared by every sub-buffer. */
    T* buffer = nullptr;
#endif
    /** @brief Number of elements already handed out by `getNew`; the start offset of the next slice. */
    size_t offset = 0;
    /** @brief Total capacity reserved in the constructor. */
    const size_t total_items;
};

