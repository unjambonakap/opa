/*
    Copyright (c) 2011-2014 University of Zurich
    
    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:
    
    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.
    
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
    THE SOFTWARE.
*/

#ifndef PLLL_INCLUDE_GUARD__KECCAK_HPP
#define PLLL_INCLUDE_GUARD__KECCAK_HPP

#include <limits>
#include <stdint.h>
#include <deque>
#include <plll/helper.hpp>

namespace plll
{
    namespace helper
    {
        namespace keccak
        {
            namespace implementation
            {
                template<unsigned ell>
                struct Keccak_helper
                {
                    enum { bitmask = Keccak_helper<ell - 1>::bitmask | (static_cast<uint64_t>(1) << (ell - 1)) };
                };
                
                template<>
                struct Keccak_helper<0>
                {
                    enum { bitmask = 0 };
                };
                
                template<>
                struct Keccak_helper<1>
                {
                    enum { bitmask = 1 };
                };
                
                template<unsigned ell>
                struct Keccak_register
                // An `ell`-bit register.
                {
                    typedef typename plll::helper::SelectFirstType<static_cast<unsigned>(std::numeric_limits<unsigned>::digits) >= ell,
                                                                   unsigned,
                                                                   typename plll::helper::SelectFirstType<static_cast<unsigned>(std::numeric_limits<unsigned long>::digits) >= ell,
                                                                                                          unsigned long,
                                                                                                          uint64_t>::result >::result register_type;
                    
                    static register_type rotate(register_type input, unsigned offset)
                    // Assumes `0 <= offset < ell`.
                    {
                        return (offset == 0)
                            ? input
                            : ( ((input << offset) & Keccak_helper<ell>::bitmask) | (input >> (ell - offset)) );
                    }
                };
                
                template<unsigned n>
                struct Keccak_power_of_two;
                
                template<>
                struct Keccak_power_of_two<25 * 1>
                {
                    enum { ell = 0 };
                    enum { ttell = 1, b = 25 * 1 };
                };
                
                template<>
                struct Keccak_power_of_two<25 * 2>
                {
                    enum { ell = 1 };
                    enum { ttell = 2, b = 25 * 2 };
                };
                
                template<>
                struct Keccak_power_of_two<25 * 4>
                {
                    enum { ell = 2 };
                    enum { ttell = 4, b = 25 * 4 };
                };
                
                template<>
                struct Keccak_power_of_two<25 * 8>
                {
                    enum { ell = 3 };
                    enum { ttell = 8, b = 25 * 8 };
                };
                
                template<>
                struct Keccak_power_of_two<25 * 16>
                {
                    enum { ell = 4 };
                    enum { ttell = 16, b = 25 * 16 };
                };
                
                template<>
                struct Keccak_power_of_two<25 * 32>
                {
                    enum { ell = 5 };
                    enum { ttell = 32, b = 25 * 32 };
                };
                
                template<>
                struct Keccak_power_of_two<25 * 64>
                {
                    enum { ell = 6 };
                    enum { ttell = 64, b = 25 * 64 };
                };
                
                template<unsigned b>
                struct Keccak_f
                {
                    typedef Keccak_power_of_two<b> Power_of_two;
                    typedef Keccak_helper<Power_of_two::ttell> Helper;
                    typedef Keccak_register<Power_of_two::ttell> Register;
                    typedef typename Register::register_type register_type;
                    
                    enum { packed_size = (b + std::numeric_limits<unsigned>::digits - 1)
                                       / std::numeric_limits<unsigned>::digits,
                           last_packed_leftover = packed_size * std::numeric_limits<unsigned>::digits - b,
                           number_of_rounds = 12 + 2 * Power_of_two::ell };
                    
                    typedef register_type state_type[5][5];
                    typedef unsigned packed_state_type[packed_size];
                    
                    static void unpack_state(state_type &, const packed_state_type &);
                    static void unpack_state_xor(state_type &, const packed_state_type &);
                    static void pack_state(packed_state_type &, const state_type &);
                    
                    static void f(state_type &);
                };
            }
            
            template<unsigned b, unsigned r>
            class Sponge
            {
            private:
                typedef implementation::Keccak_f<b> f;
                enum { buffer_fill_size = r / f::Power_of_two::ttell };
                
                typename f::state_type d_state;
                
                typename f::packed_state_type d_buffer;
                unsigned d_buffer_index;
                
                inline void flush_buffer_to_state()
                {
                    f::unpack_state_xor(d_state, d_buffer);
                    f::f(d_state);
                    d_buffer_index = 0;
                }
                
                inline void flush_state_to_buffer()
                {
                    f::f(d_state);
                    f::pack_state(d_buffer, d_state);
                    d_buffer_index = 0;
                }
                
            public:
                inline Sponge()
                    : d_buffer_index(0)
                {
                    PLLL_INTERNAL_STATIC_CHECK((r % f::Power_of_two::ttell) == 0, RMustBeMultipleOfLaneLength);
                    PLLL_INTERNAL_STATIC_CHECK(r > 0, RCannotBeZero);
                    PLLL_INTERNAL_STATIC_CHECK(r <= b, RCannotNotBeLargerThanB);
                    for (unsigned x = 0; x < 5; ++x)
                        for (unsigned y = 0; y < 5; ++y)
                            d_state[x][y] = 0;
                    for (unsigned i = buffer_fill_size; i < f::packed_size; ++i)
                        d_buffer[i] = 0;
                }
                
                inline void feed(unsigned value)
                // Put data into sponge. Must only be called before calling `squeeze()`.
                {
                    d_buffer[d_buffer_index] = value;
                    if (++d_buffer_index == buffer_fill_size)
                        flush_buffer_to_state();
                }
                
                template<typename InputIt>
                inline void feed(InputIt begin, InputIt end)
                // Put data into sponge. Iterators should point to a range of `unsigned`s. Must only be called
                // before calling `squeeze()`.
                {
                    while (begin != end)
                    {
                        feed(*begin);
                        ++begin;
                    }
                }
                
                void squeeze()
                // Squeeze sponge. Must be called precisely once.
                {
                    // Still something in the buffer?
                    if (d_buffer_index > 0)
                    {
                        // Add padding
                        d_buffer[d_buffer_index] |= 1;
                        d_buffer[buffer_fill_size - 1] |=
                            static_cast<typename f::register_type>(1)
                            << (std::numeric_limits<unsigned>::digits - f::last_packed_leftover - 1);
                        // Flush buffer
                        flush_buffer_to_state();
                    }
                    // Now unpack the state
                    f::pack_state(d_buffer, d_state);
                    d_buffer_index = 0;
                }
                
                inline unsigned read()
                // Reads one `unsigned` from the squeezed sponge. Must only be called after `squeeze()`.
                {
                    unsigned value = d_buffer[d_buffer_index];
                    if (++d_buffer_index == buffer_fill_size)
                        flush_state_to_buffer();
                    return value;
                }
                
                template<typename OutputIt>
                inline void read(OutputIt begin, OutputIt end)
                // Read a range of `unsigned`s from the squeezed sponge. Must only be called after `squeeze()`.
                {
                    while (begin != end)
                    {
                        *begin = read();
                        ++begin;
                    }
                }
            };
            
            template<unsigned b, unsigned r>
            class DuplexSponge
            {
            private:
                typedef implementation::Keccak_f<b> f;
                enum { buffer_fill_size = r / f::Power_of_two::ttell };
                
                typename f::state_type d_state;
                
                typename f::packed_state_type d_buffer;
                unsigned d_buffer_index;
                
                std::deque<unsigned> d_input_buffer;
                
                inline void populate_buffer()
                {
                    // Fill d_buffer with input
                    if (!d_input_buffer.empty())
                    {
                        // Fill buffer with input data (as much as available)
                        if (d_input_buffer.size() >= buffer_fill_size)
                        {
                            // We can fill the buffer completely
                            for (unsigned i = 0; i < buffer_fill_size; ++i)
                                d_buffer[i] = d_input_buffer[i];
                            d_input_buffer.erase(d_input_buffer.begin(), d_input_buffer.begin() + buffer_fill_size);
                        }
                        else
                        {
                            // We can fill the buffer partially
                            unsigned i = 0;
                            for (unsigned count = d_input_buffer.size(); count > 0; --count, ++i)
                                d_buffer[i] = d_input_buffer[i];
                            d_input_buffer.clear();
                            // Add padding
                            d_buffer[i] = 1;
                            for (++i; i < buffer_fill_size; ++i)
                                d_buffer[i] = 0;
                            d_buffer[buffer_fill_size - 1] |=
                                static_cast<typename f::register_type>(1)
                                << (std::numeric_limits<unsigned>::digits - f::last_packed_leftover - 1);
                        }
                        // Set rest to zero
                        for (unsigned i = buffer_fill_size; i < f::packed_size; ++i)
                            d_buffer[i] = 0;
                        // XOR to state
                        f::unpack_state_xor(d_state, d_buffer);
                    }
                    // Apply f function
                    f::f(d_state);
                    // Now unpack the state
                    f::pack_state(d_buffer, d_state);
                    d_buffer_index = 0;
                }
                
            public:
                inline DuplexSponge()
                    : d_buffer_index(0)
                {
                    PLLL_INTERNAL_STATIC_CHECK((r % f::Power_of_two::ttell) == 0, RMustBeMultipleOfLaneLength);
                    PLLL_INTERNAL_STATIC_CHECK(r > 0, RCannotBeZero);
                    PLLL_INTERNAL_STATIC_CHECK(r <= b, RCannotNotBeLargerThanB);
                    for (unsigned x = 0; x < 5; ++x)
                        for (unsigned y = 0; y < 5; ++y)
                            d_state[x][y] = 0;
                    for (unsigned i = buffer_fill_size; i < f::packed_size; ++i)
                        d_buffer[i] = 0;
                }
                
                inline void feed(unsigned value)
                // Put data into sponge. Must only be called before calling `squeeze()`.
                {
                    d_input_buffer.push_back(value);
                }
                
                template<typename InputIt>
                inline void feed(InputIt begin, InputIt end)
                // Put data into sponge. Iterators should point to a range of `unsigned`s. Must only be called
                // before calling `squeeze()`.
                {
                    d_input_buffer.insert(d_input_buffer.end(), begin, end);
                }
                
                void squeeze()
                // Squeeze sponge. Must be called precisely once.
                {
                    populate_buffer();
                }
                
                inline unsigned read()
                // Reads one `unsigned` from the squeezed sponge. Must only be called after `squeeze()`.
                {
                    unsigned value = d_buffer[d_buffer_index];
                    if (++d_buffer_index == buffer_fill_size)
                        populate_buffer();
                    return value;
                }
                
                template<typename OutputIt>
                inline void read(OutputIt begin, OutputIt end)
                // Read a range of `unsigned`s from the squeezed sponge. Must only be called after `squeeze()`.
                {
                    while (begin != end)
                    {
                        *begin = read();
                        ++begin;
                    }
                }
            };
        }
    }
}

#endif
