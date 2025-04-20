
/*
    Copyright (c) 2005-2022 Intel Corporation

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

        http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.
*/

// Adapted and modified by Bernard Gingold on 20/04/2025

/*MIT License
Copyright (c) 2020 Bernard Gingold
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/


#ifndef __GMS_CACHE_ALIGNED_ALLOCATOR_H__
#define __GMS_CACHE_ALIGNED_ALLOCATOR_H__


#if __cplusplus >= 201703L
#include <memory_resource>
#else
#error "Required a C++ version >= 201703l"
#endif
#include <cstdlib>
#include <cassert>
#include "GMS_common.h"


namespace gms {

          class cache_aligned_allocator : public std::pmr::memory_resource {
                
                std::pmr::memory_resource * __restrict m_mem_resource;

                public:

                cache_aligned_allocator()
                : cache_aligned_allocator(std::pmr::get_default_resource()) {}

                explicit cache_aligned_allocator(std::pmr::memory_resource * __restrict mem_resource)
                :
                m_mem_resource(mem_resource) {}

                std::pmr::memory_resource * __restrict get_mem_resource() const noexcept(true) 
                {
                    return this->m_mem_resource;
                }

                private:

                void * do_allocate(std::size_t nbytes, std::size_t alignment) override 
                {
                    const std::size_t cache_line_alignment = alignment_correction(alignment);
                    const std::size_t mem_space  = alignment_correction(nbytes)+cache_line_alignment;
                    std::uintptr_t    mem_base   = reinterpret_cast<std::uintptr_t>(this->m_mem_resource->allocate(mem_space));
                    assert(mem_base != 0, "NULL-pointer returned by the allocate!!");
                    const std::uintptr_t mem_result = (mem_base+cache_line_alignment) & ~(cache_line_alignment-1ULL);
                    assert((mem_result-mem_base) >= sizeof(std::uintptr_t), "Size <= to minimal allowable!!");
                    assert(mem_space-(mem_result-mem_base) >= bytes, "Running out of available memory for this request!!");

                    (reinterpret_cast<std::uintptr_t*>(mem_result))[-1ULL] = mem_base;
                    return reinterpret_cast<void*>(mem_result);
                }

                void do_deallocate(void * ptr, std::size_t nbytes, std::size_t alignment) override 
                {
                     if(ptr)
                     {
                        std::uintptr_t mem_base = (reinterpret_cast<std::uintptr_t*>(ptr))[-1];
                        this->m_mem_resource->deallocate(reinterpret_cast<void*>(mem_base),
                                                         size_correction(nbytes)+alignment_correction(alignment));
                     }
                }

                std::size_t alignment_correction(std::size_t alignment) noexcept(true)
                {
                    const bool is_pow2 = alignment && !(alignment & (alignment-1ULL));
                    assert(is_pow2, "ALignment must be a power of 2!!");
#ifdef __cpp_lib_hardware_interference_size
                    std::size_t cache_line_size = std::hardware_destructive_interference_size;
#else 
                    std::size_t cache_line_size = 64ULL;
#endif
                    return alignment < cache_line_size ? cache_line_size : alignment;
                }

                std::size_t size_correction(std::size_t nbytes) noexcept(true)
                {
                    return nbytes < sizeof(std::uintptr_t) ? sizeof(std::uintptr_t) : nbytes;
                }
          };  
}













#endif /*__GMS_CACHE_ALIGNED_ALLOCATOR__*/