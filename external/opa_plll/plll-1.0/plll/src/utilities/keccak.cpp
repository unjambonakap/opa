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

#include "keccak.hpp"
#include <iostream>
#include <cassert>

namespace plll
{
    namespace helper
    {
        namespace keccak
        {
            namespace implementation
            {
                template<unsigned b>
                void Keccak_f<b>::unpack_state(typename Keccak_f<b>::state_type & state, const typename Keccak_f<b>::packed_state_type & packed_state)
                {
                    unsigned index = 0, bitindex = 0;
                    for (unsigned y = 0; y < 5; ++y)
                        for (unsigned x = 0; x < 5; ++x)
                        {
                            register_type v = 0;
                            unsigned bits_left = Power_of_two::ttell;
                            while (bits_left)
                            {
                                unsigned bits_find = std::numeric_limits<unsigned>::digits - bitindex;
                                if (bits_find > bits_left)
                                    bits_find = bits_left;
                                // Take the correct number of bits from packed_state[index] and put them into v
                                register_type ww = packed_state[index] >> bitindex;
                                v |= ww << (Power_of_two::ttell - bits_left);
                                // Increase counters
                                bitindex += bits_find;
                                bits_left -= bits_find;
                                if (bitindex == static_cast<unsigned>(std::numeric_limits<unsigned>::digits))
                                {
                                    ++index;
                                    bitindex = 0;
                                }
                            }
                            state[x][y] = v & Helper::bitmask;
                        }
                }
                
                template<unsigned b>
                void Keccak_f<b>::unpack_state_xor(typename Keccak_f<b>::state_type & state, const typename Keccak_f<b>::packed_state_type & packed_state)
                {
                    unsigned index = 0, bitindex = 0;
                    for (unsigned y = 0; y < 5; ++y)
                        for (unsigned x = 0; x < 5; ++x)
                        {
                            register_type v = 0;
                            unsigned bits_left = Power_of_two::ttell;
                            while (bits_left)
                            {
                                unsigned bits_find = std::numeric_limits<unsigned>::digits - bitindex;
                                if (bits_find > bits_left)
                                    bits_find = bits_left;
                                // Take the correct number of bits from packed_state[index] and put them into v
                                register_type ww = packed_state[index] >> bitindex;
                                v |= ww << (Power_of_two::ttell - bits_left);
                                // Increase counters
                                bitindex += bits_find;
                                bits_left -= bits_find;
                                if (bitindex == static_cast<unsigned>(std::numeric_limits<unsigned>::digits))
                                {
                                    ++index;
                                    bitindex = 0;
                                }
                            }
                            state[x][y] ^= v & Helper::bitmask;
                        }
                }
                
                template<unsigned b>
                void Keccak_f<b>::pack_state(typename Keccak_f<b>::packed_state_type & packed_state, const typename Keccak_f<b>::state_type & state)
                {
                    unsigned index = 0, bitindex = 0;
                    unsigned current_ps = 0;
                    for (unsigned y = 0; y < 5; ++y)
                        for (unsigned x = 0; x < 5; ++x)
                        {
                            register_type v = state[x][y];
                            unsigned bits_left = Power_of_two::ttell;
                            while (bits_left)
                            {
                                unsigned bits_find = std::numeric_limits<unsigned>::digits - bitindex;
                                if (bits_find > bits_left)
                                    bits_find = bits_left;
                                // Take the correct number of bits from v and put them into current_ps
                                current_ps |= v << bitindex;
                                // Increase counters
                                bitindex += bits_find;
                                bits_left -= bits_find;
                                v >>= bits_find;
                                if (bitindex == static_cast<unsigned>(std::numeric_limits<unsigned>::digits))
                                {
                                    packed_state[index] = current_ps;
                                    ++index;
                                    bitindex = 0;
                                    current_ps = 0;
                                }
                            }
                        }
                }
                
                static uint64_t Keccak_round_constants[24] =
                {
                    0x0000000000000001ULL,
                    0x0000000000008082ULL,
                    0x800000000000808AULL,
                    0x8000000080008000ULL,
                    0x000000000000808BULL,
                    0x0000000080000001ULL,
                    0x8000000080008081ULL,
                    0x8000000000008009ULL,
                    0x000000000000008AULL,
                    0x0000000000000088ULL,
                    0x0000000080008009ULL,
                    0x000000008000000AULL,
                    0x000000008000808BULL,
                    0x800000000000008BULL,
                    0x8000000000008089ULL,
                    0x8000000000008003ULL,
                    0x8000000000008002ULL,
                    0x8000000000000080ULL,
                    0x000000000000800AULL,
                    0x800000008000000AULL,
                    0x8000000080008081ULL,
                    0x8000000000008080ULL,
                    0x0000000080000001ULL,
                    0x8000000080008008ULL,
                };
                
                static unsigned rotation_offsets[5][5] =
                { { 3, 36, 0, 210, 105 }, { 10, 300, 1, 66, 45 }, { 171, 6, 190, 253, 15 }, { 153, 55, 28, 120, 21 }, { 231, 276, 91, 78, 136 } };
                
                template<unsigned b>
                void Keccak_f<b>::f(typename Keccak_f<b>::state_type & state)
                {
                    for (unsigned round = 0; round < number_of_rounds; ++round)
                    {
                        // Theta step
                        register_type C[5], D[5];
                        for (unsigned x = 0; x < 5; ++x)
                            C[x] = state[x][0] ^ state[x][1] ^ state[x][2] ^ state[x][3] ^ state[x][4];
                        for (unsigned x = 0; x < 5; ++x)
                        {
                            D[x] = C[(x - 1) % 5] ^ Register::rotate(C[(x + 1) % 5], 1);
                            for (unsigned y = 0; y < 5; ++y)
                                state[x][y] ^= D[x];
                        }
                        
                        // Rho and pi steps
                        register_type B[5][5];
                        for (unsigned x = 0; x < 5; ++x)
                            for (unsigned y = 0; y < 5; ++y)
                                B[y][(2 * x + 3 * y) % 5] = Register::rotate(state[x][y], rotation_offsets[x][y] & (Power_of_two::ttell - 1));
                        
                        // Chi step
                        for (unsigned x = 0; x < 5; ++x)
                            for (unsigned y = 0; y < 5; ++y)
                                state[x][y] = B[x][y] ^ ((~B[(x + 5) % 5][y]) & B[(x + 2) % 5][y]);
                        
                        // Iota step
                        state[0][0] ^= Keccak_round_constants[round] & Helper::bitmask;
                    }
                }
                
                // Explicitly instantiate Keccak f functions
                template struct Keccak_f<25>;
                template struct Keccak_f<50>;
                template struct Keccak_f<100>;
                template struct Keccak_f<200>;
                template struct Keccak_f<400>;
                template struct Keccak_f<800>;
                template struct Keccak_f<1600>;
            }
        }
    }
}
