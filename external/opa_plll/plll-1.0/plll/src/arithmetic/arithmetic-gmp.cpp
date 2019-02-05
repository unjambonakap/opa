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

#include <plll/config.hpp>
#include <plll/arithmetic.hpp>
#include <plll/helper.hpp>
#include "myalloc.hpp"
#include "buffer.hpp"
#include "keccak.hpp"
#include <cstring>
#include <iostream>
#include <fstream>
#include <boost/thread/tss.hpp>
#include <boost/thread.hpp>
#include <boost/random/mersenne_twister.hpp>
#if (BOOST_VERSION >= 104700)
  #include <boost/random/seed_seq.hpp>
#endif

#if defined(PLLL_INTERNAL_HR_POSIX)
  #include <time.h>
#endif

namespace plll
{
    namespace arithmetic
    {
#ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
        namespace internal
        {
            bool Real_precision_check_enabled = true;
        };
#endif
        
        void initArithmeticThreadAllocators()
        {
            // Sets the Thread-Local Secure Allocator as the standard allocator for GMP and MPFR.
            
            initTLAlloc(); // Make sure that the allocator is initialized for this thread
            
            mp_set_memory_functions(TLSalloc, TLSrealloc, TLSfree); // Set the GMP/MPFR allocation
                                                                    // functions
        }
        
        namespace
        {
            // This is not visible to other modules
            boost::thread_specific_ptr<RealContext> g_realcontexts;
            
            template<typename T>
            class ArrayStore
            {
            private:
                T * d_storage;
                
            public:
                ArrayStore(std::size_t n = 0)
                    : d_storage(n ? new T[n] : NULL)
                {
                }
                
                ~ArrayStore()
                {
                    delete [] d_storage;
                }
                
                void reset(std::size_t n)
                {
                    delete [] d_storage;
                    d_storage = new T[n];
                }
                
                void release()
                {
                    delete [] d_storage;
                    d_storage = NULL;
                }
                
                T * get()
                {
                    return d_storage;
                }
                
                const T * get() const
                {
                    return d_storage;
                }
                
                T & operator [] (std::size_t i)
                {
                    return d_storage[i];
                }
                
                T & operator [] (std::size_t i) const
                {
                    return d_storage[i];
                }
            };
        }
        
        RealContext & getThreadRealContext()
        // returns context for current thread
        {
            if (g_realcontexts.get() == NULL)
                g_realcontexts.reset(new RealContext(true));
            return *g_realcontexts;
        }
        
        template<typename Sponge>
        class SpongeSeeder
        {
        private:
            Sponge & d_sponge;
            
        public:
            SpongeSeeder(Sponge & sponge)
                : d_sponge(sponge)
            {
            }
            
            template<typename It>
            void generate(It begin, It end) const
            // For boost 1.47.0 or newer
            {
                while (begin != end)
                {
                    *begin = d_sponge.read();
                    ++begin;
                }
            }
            
            unsigned operator () () const
            // For boost 1.46.1 or older
            {
                return d_sponge.read();
            }
        };
        
        class RNGImpl
        {
#if (BOOST_VERSION >= 104700)
        public:
            enum { statelen = 312 * 8 }; // number of bytes per state
            
        private:
            boost::random::mt19937_64 d_rng;
            typedef uint64_t OutputT; // this is what d_rng() returns
            enum { outchunks = sizeof(OutputT) }; // number of chars in OutputT
#else
        public:
            enum { statelen = 624 * 4 }; // number of bytes per state
            
        private:
            boost::mt19937 d_rng; // note the missing namespace random::, which is apparently not
                                  // used here in older version of boost...
            typedef uint32_t OutputT; // this is what d_rng() returns
            enum { outchunks = sizeof(OutputT) }; // number of chars in OutputT
#endif
            
        public:
            inline void seed(const mpz_t & x)
            {
#if (BOOST_VERSION >= 104700)
                // First, get uint32_t array
                size_t size = mpz_size(x) * mp_bits_per_limb / 8; // first, find out size in bytes
                size = (size + sizeof(uint32_t) - 1) / sizeof(uint32_t); // convert to number of uint32_t's
                ArrayStore<uint32_t> buffer(size);
                size_t s;
                mpz_export(buffer.get(), &s, 1, sizeof(uint32_t), 0, 0, x);
                assert(s == size);
                // Next, generate sequence seeder
                boost::random::seed_seq seedseq(buffer.get(), buffer.get() + size);
                buffer.release();
                // Seed RNG with sequence seeder
                d_rng.seed(seedseq);
#else
                // First, get unsigned array
                size_t size = mpz_size(x) * mp_bits_per_limb / 8; // first, find out size in bytes
                size = (size + sizeof(unsigned) - 1) / sizeof(unsigned); // convert to number of unsigned's
                ArrayStore<unsigned> buffer(size);
                size_t s;
                mpz_export(buffer.get(), &s, 1, sizeof(unsigned), 0, 0, x);
                assert(s == size);
                // Next, generate and feed sponge
                typedef helper::keccak::DuplexSponge<1600, 1024> Sponge;
                Sponge sponge;
                sponge.feed(buffer.get(), buffer.get() + size);
                buffer.release();
                sponge.squeeze();
                // Seed RNG with sponge seeder
                SpongeSeeder<Sponge> seeder(sponge);
                d_rng.seed(seeder);
#endif
            }
            
            inline void genLimbs(mpz_t & x, size_t limbs)
            {
                unsigned count = (limbs * mp_bits_per_limb / 8 + sizeof(OutputT) - 1) / sizeof(OutputT);
                ArrayStore<OutputT> data(count);
                for (unsigned i = 0; i < count; ++i)
                    data[i] = d_rng();
                mpz_import(x, limbs, 1, sizeof(mp_limb_t), 0, 0, data.get());
            }
            
            inline void genBytes(char * dest, size_t count)
            {
                // "Giant steps"
                while (count >= sizeof(OutputT))
                {
                    *((OutputT*)dest) = d_rng();
                    dest += sizeof(OutputT);
                    count -= sizeof(OutputT);
                }
                // "Baby steps"
                if (count)
                {
                    OutputT r = d_rng();
                    while (count)
                    {
                        *dest = r;
                        ++dest;
                        r >>= 64 / sizeof(OutputT);
                        --count;
                    }
                }
            }
            
            inline void genBits(mpz_t & x, size_t bits)
            {
                size_t limbs = (bits + mp_bits_per_limb - 1) / mp_bits_per_limb;
                genLimbs(x, limbs);
                // Make sure we have the right size
                if (limbs * mp_bits_per_limb > bits)
                    mpz_tdiv_q_2exp(x, x, limbs * mp_bits_per_limb - bits);
            }
            
            inline void genRand(mpz_t & x, const mpz_t & bound)
            {
                assert(mpz_sgn(bound) > 0);
                
                size_t limbs = mpz_size(x) + 1;
                mpz_t xx;
                mpz_init(xx);
                do
                {
                    // Generate random number (using full length of bound)
                    genLimbs(xx, limbs);
                    mpz_mod(x, xx, bound);
                    if (!mpz_tstbit(xx, limbs * mp_bits_per_limb - 1))
                        break;
                    mpz_sub(xx, xx, x);
                    mpz_add(xx, xx, bound);
                }
                while (mpz_size(xx) > limbs);
                mpz_clear(xx);
            }
            
        private:
            static inline bool isPOT(unsigned long ul)
            // C++ standard enforces two's complement for unsigned integer types; therefore, this
            // works on any platform.
            {
                return (ul & (-ul)) == ul;
            }
            
            inline unsigned long genUL()
            {
                if (sizeof(unsigned long) > sizeof(OutputT))
                {
                    unsigned long output = d_rng();
                    for (unsigned i = 1; i < sizeof(unsigned long) / sizeof(OutputT); ++i)
                    {
                        output <<= (std::numeric_limits<OutputT>::digits % std::numeric_limits<unsigned long>::digits);
                        // The only reason of the modulo is to stop GCC from producing a warning
                        // ("left shift count >= width of type"), which is absurd, since in case the
                        // warning is true, this part of the if clause will never be executed...
                        output |= d_rng();
                    }
                    return output;
                }
                else
                    return d_rng();
            }
            
        public:
            inline unsigned long genRand(unsigned long bound)
            {
                assert(bound > 0);
                // Special case: power of 2
                if (isPOT(bound))
                    return genUL() & (bound - 1);
                unsigned long v; 
                // Large?
                if (bound > std::numeric_limits<unsigned long>::max() / 2)
                {
                    // General case: try until done
                    // Expected: at most two iterations in average
                    do
                    {
                        v = genUL();
                    }
                    while (v >= bound);
                }
                else
                {
                    // General case: try until done
                    // Know: 2 * bound < maximal value of unsigned long
                    unsigned long r, end = std::numeric_limits<unsigned long>::max() - bound;
                    do
                    {
                        // Generate random number (using full length of bound)
                        r = genUL();
                        v = r % bound;
                    }
                    while (r - v > end);
                }
                return v;
            }
        };
        
        RandomNumberGenerator::RandomNumberGenerator()
        {
            d_state = new RNGImpl();
        }
        
        RandomNumberGenerator::RandomNumberGenerator(const RandomNumberGenerator & rg)
            : d_seed(rg.d_seed)
        {
            RNGImpl * other_state = reinterpret_cast<RNGImpl*>(rg.d_state);
            d_state = new RNGImpl(*other_state);
        }
        
        RandomNumberGenerator::~RandomNumberGenerator()
        {
            RNGImpl * state = reinterpret_cast<RNGImpl*>(d_state);
            delete state;
        }
        
        RandomNumberGenerator & RandomNumberGenerator::operator = (const RandomNumberGenerator & rg)
        {
            RNGImpl * state = reinterpret_cast<RNGImpl*>(d_state);
            if (this != &rg)
            {
                RNGImpl * other_state = reinterpret_cast<RNGImpl*>(rg.d_state);
                *state = *other_state;
                d_seed = rg.d_seed;
            }
            return *this;
        }
        
        Integer RandomNumberGenerator::getSeed() const
        {
            return d_seed;
        }
        
        void RandomNumberGenerator::setSeed(const Integer & seed)
        {
            RNGImpl * state = reinterpret_cast<RNGImpl*>(d_state);
            d_seed = seed;
            state->seed(seed.d_value);
        }
        
        namespace
        {
            boost::mutex s_rng_mutex;
            RandomNumberGenerator * s_rng = NULL;
        }
        
        static void initSeeder(double goodness)
        {
            if (!s_rng)
            {
                // Warning! This could very well be blocking, depending on the size of the system's entropy cache!
                s_rng = new RandomNumberGenerator();
                std::cout << "Initializing seeder using /dev/(u)random (might be blocking)..." << std::flush;
                s_rng->randomizeDevRandom(false/*true*/);
                std::cout << " done.\n";
            }
        }
        
        Integer RandomNumberGenerator::createSeed()
        // create seed from internal RNG (initialized using /dev/random)
        {
            Integer res;
            s_rng_mutex.lock();
            initSeeder(0.1);
            res = s_rng->randomBits(RNGImpl::statelen * 8);
            s_rng_mutex.unlock();
            return res;
        }
        
        void RandomNumberGenerator::randomizeTime()
        {
#if defined(PLLL_INTERNAL_HR_POSIX)
            timespec t;
            clock_gettime(CLOCK_REALTIME, &t);
#else
            time_t t = time(0);
#endif
            Integer seed;
            mpz_import(seed.d_value, sizeof(t), 1, 1, 0, 0, &t);
            setSeed(seed);
            seed = randomBits(RNGImpl::statelen * 8);
            setSeed(seed);
        }
        
        void RandomNumberGenerator::randomizeDevRandom(double goodness)
        // Parameter is in range [0, 1]. Decides how much is used from /dev/random and how much from /dev/urandom
        {
            if (goodness < 0)
                goodness = 0;
            else if (goodness > 1)
                goodness = 1;
            // Use /dev/random and /dev/urandom to initialize RNG
            const int size = RNGImpl::statelen;
            // Map [0,1] monotonically increasing to [0, size] such that 0.01 is mapped to 2
            int size1 = 0.5 + ((goodness <= 0.01) ? (2 * (goodness / 0.01)) : (2 + (goodness - 0.01) / (1 - 0.01) * (size - 2)));
            if (size1 < 0)
                size1 = 0;
            else if (size1 > size)
                size1 = size;
            // Fill the bytes
            char bytes[size];
            std::ifstream f;
            if (size1 > 0)
            {
                f.open("/dev/random");
                f.read(bytes, size1);
                f.close();
            }
            if (size1 < size)
            {
                f.open("/dev/urandom");
                f.read(bytes + size1, size - size1);
                f.close();
            }
            // Create seed
            Integer seed;
            mpz_import(seed.d_value, size, 1, 1, 0, 0, bytes);
            setSeed(seed);
        }
        
        void RandomNumberGenerator::initializeSeeder(double goodness)
        // Can (optionally) be called so that the seeder is ready when randomizeSeed() is called
        {
            if (!s_rng)
            {
                s_rng_mutex.lock();
                initSeeder(goodness);
                s_rng_mutex.unlock();
            }
        }
        
        void RandomNumberGenerator::randomizeSeed()
        {
            setSeed(createSeed());
        }
        
        void RandomNumberGenerator::random(Integer & r, const Integer & bound)
        // returns a random integer x with 0 <= x < bound
        {
            RNGImpl * state = reinterpret_cast<RNGImpl*>(d_state);
            state->genRand(r.d_value, bound.d_value);
        }
        
        Integer RandomNumberGenerator::random(const Integer & bound)
        // returns a random integer x with 0 <= x < bound
        {
            RNGImpl * state = reinterpret_cast<RNGImpl*>(d_state);
            Integer r;
            state->genRand(r.d_value, bound.d_value);
            return r;
        }
        
        unsigned long RandomNumberGenerator::random(unsigned long bound)
        // returns a random integer x with 0 <= x < bound
        {
            RNGImpl * state = reinterpret_cast<RNGImpl*>(d_state);
            return state->genRand(bound);
        }
        
        void RandomNumberGenerator::randomBits(void * dest, unsigned long size)
        {
            RNGImpl * state = reinterpret_cast<RNGImpl*>(d_state);
            state->genBytes(reinterpret_cast<char*>(dest), size);
        }
        
        void RandomNumberGenerator::randomBits(Integer & r, unsigned long bits)
        // returns a random integer x with 0 <= x < 2^bits
        {
            RNGImpl * state = reinterpret_cast<RNGImpl*>(d_state);
            state->genBits(r.d_value, bits);
        }
        
        Integer RandomNumberGenerator::randomBits(unsigned long bits)
        // returns a random integer x with 0 <= x < 2^bits
        {
            RNGImpl * state = reinterpret_cast<RNGImpl*>(d_state);
            Integer r;
            state->genBits(r.d_value, bits);
            return r;
        }
        
        void RandomNumberGenerator::randomLen(Integer & r, unsigned long bits)
        // returns a random integer x with 2^(bits-1) <= x < 2^bits
        {
            RNGImpl * state = reinterpret_cast<RNGImpl*>(d_state);
            if (bits)
            {
                if (bits > 1)
                    state->genBits(r.d_value, bits - 1);
                setbit(r, bits - 1, true);
            }
            else
                setZero(r);
        }
        
        Integer RandomNumberGenerator::randomLen(unsigned long bits)
        // returns a random integer x with 2^(bits-1) <= x < 2^bits
        {
            RNGImpl * state = reinterpret_cast<RNGImpl*>(d_state);
            Integer r;
            if (bits)
            {
                if (bits > 1)
                    state->genBits(r.d_value, bits - 1);
                setbit(r, bits - 1, true);
            }
            return r;
        }
        
        void RandomNumberGenerator::randomUniform(Real & r)
        // generates a uniform Real number in [0, 1)
        {
            unsigned long b = mpfr_get_default_prec();
            mpfr_set_z(r.d_value, randomBits(b).d_value, MPFR_RNDN);
            shr(r, r, b);
        }
        
        void RandomNumberGenerator::randomUniform(Real & r, const RealContext & rc)
        // generates a uniform Real number in [0, 1)
        {
            unsigned long b = rc.getRealPrecision();
            mpfr_set_z(r.d_value, randomBits(b).d_value, MPFR_RNDN);
            shr(r, r, b);
        }
        
        void randomUniform(float & r, RandomNumberGenerator & rng)
        {
            // Assume that mantissa is at most 64 bits
            uint64_t x;
            rng.randomBits(&x, sizeof(x));
            float rr = x;
            r = ldexp(rr, -64);
        }
        
        void randomUniform(double & r, RandomNumberGenerator & rng)
        {
            // Assume that mantissa is at most 64 bits
            uint64_t x;
            rng.randomBits(&x, sizeof(x));
            double rr = x;
            r = ldexp(rr, -64);
        }
        
        void randomUniform(long double & r, RandomNumberGenerator & rng)
        {
    #if PLLL_INTERNAL_LONGDOUBLE == PLLL_INTERNAL_LONGDOUBLE_EP
            uint64_t x;
            rng.randomBits(&x, sizeof(x));
            long double rr = x;
            r = ldexp(rr, -64);
    #elif PLLL_INTERNAL_LONGDOUBLE == PLLL_INTERNAL_LONGDOUBLE_QP
            uint64_t x;
            rng.randomBits(&x, sizeof(x));
            long double rr = x;
            r = ldexp(rr, -128);
            rng.randomBits(&x, sizeof(x));
            rr = x;
            r += ldexp(rr, -64);
    #endif
        }
        
        static inline bool isPOT(unsigned long long ul)
        // C++ standard enforces two's complement for unsigned integer types; therefore, this works
        // on any platform.
        {
            return (ul & (-ul)) == ul;
        }
        
        namespace implementation
        {
            void randomInt(int & l, long b, RandomNumberGenerator & rng)
            {
                l = rng.random(b);
            }
            
            void randomInt(long & l, long b, RandomNumberGenerator & rng)
            {
                l = rng.random(b);
            }
            
            void randomInt(long long & v, long long bound, RandomNumberGenerator & rng)
            {
                if (sizeof(long long) <= sizeof(long))
                {
                    v = rng.random((long)bound);
                }
                else
                {
                    assert(bound > 0);
                    // Special case: power of 2
                    if (isPOT(bound))
                    {
                        rng.randomBits(&v, sizeof(v));
                        v &= bound - 1;
                        return;
                    }
                    // Large?
                    if (bound * 2 < bound)
                    {
                        // General case: try until done
                        // Expected: at most two iterations in average
                        do
                        {
                            rng.randomBits(&v, sizeof(v));
                        }
                        while (v >= bound);
                    }
                    else
                    {
                        // General case: try until done
                        // Know: 2 * bound < maximal value of unsigned long
                        unsigned long long r;
                        do
                        {
                            // Generate random number (using full length of bound)
                            rng.randomBits(&r, sizeof(r));
                            v = r % bound;
                        }
                        while (r - v + (unsigned long long)bound < (unsigned long long)bound);
                    }
                }
            }
            
            void randomBits(int & l, unsigned long b, RandomNumberGenerator & rng)
            {
                PLLL_INTERNAL_STATIC_CHECK(std::numeric_limits<unsigned int>::radix == 2, Implementation_requires_radix_to_be_2);
                unsigned int ll;
                rng.randomBits(&ll, sizeof(l));
                l = (int)(ll >> (std::numeric_limits<unsigned int>::digits - b));
            }
            
            void randomBits(long & l, unsigned long b, RandomNumberGenerator & rng)
            {
                PLLL_INTERNAL_STATIC_CHECK(std::numeric_limits<unsigned long>::radix == 2, Implementation_requires_radix_to_be_2);
                unsigned long ll;
                rng.randomBits(&ll, sizeof(l));
                l = (long)(ll >> (std::numeric_limits<unsigned long>::digits - b));
            }
            
            void randomBits(long long & l, unsigned long b, RandomNumberGenerator & rng)
            {
                PLLL_INTERNAL_STATIC_CHECK(std::numeric_limits<unsigned long long>::radix == 2, Implementation_requires_radix_to_be_2);
                unsigned long long ll;
                rng.randomBits(&ll, sizeof(l));
                l = (long long)(ll >> (std::numeric_limits<unsigned long long>::digits - b));
            }
        }
        
        void Integer::mpz_init_set_ld(mpz_t & v, long double d)
        {
            mpz_init(v);
            Integer::mpz_set_ld(v, d);
        }
        
        struct LDExtractEP
        // This only works if long double is in Extended Precision (80 bit) format!
        {
            unsigned long mantissa : 64;
            unsigned long exp : 15;
            unsigned long sign : 1;
            unsigned long rest : 48;
        };
        
        unsigned find_msb_64(unsigned long x)
        // Assume that unsigned long has precisely 64 bits. x is garuanteed to be non-zero.
        {
            // Adapted from http://www.hackersdelight.org/HDcode/nlz.c.txt
            unsigned n = 0;
            if (x <= 0x00000000FFFFFFFFul) { n = n + 32; x = x << 32; }
            if (x <= 0x0000FFFFFFFFFFFFul) { n = n + 16; x = x << 16; }
            if (x <= 0x00FFFFFFFFFFFFFFul) { n = n +  8; x = x <<  8; }
            if (x <= 0x0FFFFFFFFFFFFFFFul) { n = n +  4; x = x <<  4; }
            if (x <= 0x3FFFFFFFFFFFFFFFul) { n = n +  2; x = x <<  2; }
            if (x <= 0x7FFFFFFFFFFFFFFFul) { n = n +  1; }
            return 63 - n;
        }
        
        namespace
        {
            long double anon_ld_infty = std::numeric_limits<long double>::infinity();
            long double anon_ld_minfty = -std::numeric_limits<long double>::infinity();
        }
        
        void Integer::mpz_set_ld(mpz_t & v, long double d)
        // As mpz_set_d, we truncate (and do not round or something like that).
        {
#if GMP_LIMB_BITS == 64
            PLLL_INTERNAL_STATIC_CHECK(sizeof(long double) == 16, LongDoubleShouldHave16Bytes);
            PLLL_INTERNAL_STATIC_CHECK(std::numeric_limits<long double>::has_infinity, LongDoubleShouldHaveInfinity);
            // Check for NaN, +-inf, +-0. This can easily be done platform independent.
            if ((d != d) || // NaN
                (d == anon_ld_infty) || (d == anon_ld_minfty) || // +-inf
                (d == 0.0l)) // +-0
            {
                mpz_set_ui(v, 0);
                return;
            }
  #define PLLL_INTERNAL_LONGDOUBLE_EP 1
  #define PLLL_INTERNAL_LONGDOUBLE_QP 2
  #if PLLL_INTERNAL_LONGDOUBLE == PLLL_INTERNAL_LONGDOUBLE_EP
            union {
                long double ld;
                LDExtractEP lde;
            } ld;
            ld.ld = d;
            
            // Get mantissa
            mpz_set_ui(v, ld.lde.mantissa);
            
            // Shift correctly
            signed exp = (signed)ld.lde.exp - (signed)((1ul << 14) - 1);
            exp -= 63;
            if (exp < 0)
                mpz_tdiv_q_2exp(v, v, -exp);
            else
                mpz_mul_2exp(v, v, exp);
            
            // Consider sign
            if (ld.lde.sign)
                mpz_neg(v, v);
  #elif PLLL_INTERNAL_LONGDOUBLE == PLLL_INTERNAL_LONGDOUBLE_QP
            union
            {
                long double ld;
                unsigned long l[2];
            } ld;
            ld.ld = d;
            
            // Get mantissa
            mpz_set_ui(v, (ld.l[0] & 0x0000FFFFFFFFFFFF) | 0x0001000000000000);
            mpz_mul_2exp(v, v, 64);
            mpz_add_ui(v, v, ld.l[1]);
            
            // Get exponent
            signed exp = (signed)((ld.l[0] & 0x7FFF000000000000) >> 48) - (signed)((1ul << 14) - 1);
            exp -= 112;
            if (exp < 0)
                mpz_tdiv_q_2exp(v, v, -exp);
            else
                mpz_mul_2exp(v, v, exp);
            
            // Consider sign
            if (ld.l[0] & 0x8000000000000000)
                mpz_neg(v, v);
  #else // #if PLLL_INTERNAL_LONGDOUBLE
    #error Long double type not supported yet!
  #endif
#else
            // Fallback: use MPFR. This is *SLOW*
            mpfr_t tmp;
            mpfr_init2(tmp, std::numeric_limits<long double>::digits);
            mpfr_set_ld(tmp, d, MPFR_RNDN);
            mpfr_get_z(v, tmp, MPFR_RNDN);
            mpfr_clear(tmp);
#endif
        }
        
        long double Integer::mpz_get_ld(const mpz_t & v)
        {
#if GMP_LIMB_BITS == 64
            PLLL_INTERNAL_STATIC_CHECK(sizeof(long double) == 16, LongDoubleShouldHave16Bytes);
            PLLL_INTERNAL_STATIC_CHECK(std::numeric_limits<long double>::has_infinity, LongDoubleShouldHaveInfinity);
  #define PLLL_INTERNAL_LONGDOUBLE_EP 1
  #define PLLL_INTERNAL_LONGDOUBLE_QP 2
  #if PLLL_INTERNAL_LONGDOUBLE == PLLL_INTERNAL_LONGDOUBLE_EP
            union
            {
                long double ld;
                LDExtractEP lde;
            } ld;
            ld.ld = 0.0l;
            if (mpz_sgn(v) == 0)
                return ld.ld;
            if (mpz_sgn(v) < 0)
                ld.lde.sign = 1;
            size_t size = mpz_size(v);
            mp_limb_t msl = mpz_getlimbn(v, size - 1);
            mp_limb_t smsl = size > 1 ? mpz_getlimbn(v, size - 2) : 0;
            size_t msb = find_msb_64(msl);
            if (msb < 64-1)
            {
                msl <<= (64-1) - msb;
                msl |= (mp_limb_t)smsl >> (msb + 1);
            }
            signed long exp = ((1ul << 14) - 1) + msb + (size - 1) * 64;
            // round?
            if ((smsl >> msb) & 1)
            {
                // overflow?
                if (msl == 0xFFFFFFFFFFFFFFFF)
                {
                    ++exp;
                    msl = 0x8000000000000000;
                }
                else
                    ++msl;
            }
            if (exp >= 0x7FFF)
            {
                // Too large: this is infinity!
                return mpz_sgn(v) < 0 ? -std::numeric_limits<long double>::infinity() : std::numeric_limits<long double>::infinity();
            }
            if (exp <= 0)
            {
                // Too small: ...
                return mpz_sgn(v) < 0 ? -0.0l : 0.0l; // one should use denormalized values instead!!! ??? ...
            }
            // now the MSB of msl is 1
            ld.lde.mantissa = msl;
            ld.lde.exp = (unsigned long)exp;
            return ld.ld;
  #elif PLLL_INTERNAL_LONGDOUBLE == PLLL_INTERNAL_LONGDOUBLE_QP
            union
            {
                long double ld;
                unsigned long l[2];
            } ld;
            ld.ld = 0.0l;
            if (mpz_sgn(v) == 0)
                return ld.ld;
            if (mpz_sgn(v) < 0)
                ld.l[0] |= (1ul << (64 - 1));
            size_t size = mpz_size(v);
            mp_limb_t msl = mpz_getlimbn(v, size - 1);
            mp_limb_t smsl = size > 1 ? mpz_getlimbn(v, size - 2) : 0;
            mp_limb_t tmsl = size > 2 ? mpz_getlimbn(v, size - 3) : 0;
            size_t msb = find_msb_64(msl);
            if (msb < 64-1)
            {
                msl <<= (64-1) - msb;
                msl |= (unsigned long)smsl >> (msb + 1);
                smsl <<= (64-1) - msb;
                smsl |= (unsigned long)tmsl >> (msb + 1);
            }
            bool round = (smsl >> 14) & 1;
            smsl >>= 15;
            smsl |= msl << (64 - 15);
            msl >>= 15;
            signed long exp = ((1ul << 14) - 1) + msb + (size - 1) * 64;
            // round?
            if (round)
            {
                // overflow?
                if (smsl == 0xFFFFFFFFFFFFFFFF)
                {
                    smsl = 0;
                    ++msl;
                    if (msl == 0x0001FFFFFFFFFFFF)
                    {
                        ++exp;
                        msl = 0x0001000000000000;
                    }
                }
                else
                    ++msl;
            }
            if (exp >= 0x7FFF)
            {
                // Too large: this is infinity!
                return mpz_sgn(v) < 0 ? -std::numeric_limits<long double>::infinity() : std::numeric_limits<long double>::infinity();
            }
            if (exp <= 0)
            {
                // Too small: ...
                return mpz_sgn(v) < 0 ? -0.0l : 0.0l; // one should use denormalized values instead!!! ??? ...
            }
            // The leading 1 of the mantissa is thrown away...
            ld.l[0] |= ((unsigned long)exp << 48) | (msl & 0x0000FFFFFFFFFFFF);
            ld.l[1] = smsl;
            return ld.ld;
  #else // #if PLLL_INTERNAL_LONGDOUBLE
    #error Long double type not supported yet!
  #endif
#else // #if GMP_LIMB_BITS
            // Fallback: use MPFR. This is *SLOW*
            mpfr_t tmp;
            mpfr_init2(tmp, std::numeric_limits<long double>::digits);
            mpfr_set_z(tmp, v, MPFR_RNDN);
            long double r = mpfr_get_ld(tmp, MPFR_RNDN);
            mpfr_clear(tmp);
            return r;
#endif
        }
        
        std::ostream & operator << (std::ostream & s, const Integer & i)
        // Output to stream
        {
            char * str = mpz_get_str(NULL, 10, i.d_value);
            s << str;
            mpfr_free_str(str); // There is no nice way to directly call the free() method of
            // GMP. But since MPFR is using the same allocators as GMP, we can
            // just use mpfr_free_str() instead...
            // std::free(str);
            return s;
        }
        
        std::istream & operator >> (std::istream & s, Integer & i)
        // Input from stream
        {
            Buffer<char> buf;
            char c;
            c = s.get();
            if (!s)
            {
                mpz_set_ui(i.d_value, 0);
                s.setstate(std::ios_base::failbit);
                return s;
            }
            if ((c == '-') || (c == '+'))
            {
                buf.add(c);
                c = s.get();
                if (!s)
                {
                    mpz_set_ui(i.d_value, 0);
                    s.setstate(std::ios_base::failbit);
                    return s;
                }
            }
            if (!isdigit(c))
            {
                mpz_set_ui(i.d_value, 0);
                s.setstate(std::ios_base::failbit);
                return s;
            }
            while (isdigit(c))
            {
                buf.add(c);
                c = s.get();
                if (!s)
                {
                    if (s.eof())
                        break;
                    mpz_set_ui(i.d_value, 0);
                    s.setstate(std::ios_base::failbit);
                    return s;
                }
            }
            if (!s.eof())
                s.unget();
            else
            {
                if (!s.bad())
                    s.clear(std::ios_base::eofbit);
            }
            buf.add(0);
            if (mpz_set_str(i.d_value, buf.data(), 10))
            {
                // Not a valid integer, set to zero
                mpz_set_ui(i.d_value, 0);
                // Mark stream as bad
                s.setstate(std::ios_base::failbit);
            }
            return s;
        }
        
        namespace implementation
        {
            bool from_string_conversion<IntegerContext>::convert(Integer & res, const std::string & s, const IntegerContext & c)
            {
                return mpz_set_str(res.d_value, s.c_str(), 10) == 0;
            }
            
            bool from_string_conversion<IntegerContext>::convert(Integer & res, const char * s, const IntegerContext & c)
            {
                return mpz_set_str(res.d_value, s, 10) == 0;
            }
            
            std::string to_string_conversion<Integer>::convert(const Integer & val)
            {
                char * str = mpz_get_str(NULL, 10, val.d_value);
                std::string res(str);
                mpfr_free_str(str); // There is no nice way to directly call the free() method of
                                    // GMP. But since MPFR is using the same allocators as GMP, we can
                                    // just use mpfr_free_str() instead...
                return res;
            }
            
            std::string to_string_conversion<Integer>::convert(const Integer & val, unsigned base)
            {
                char * str = mpz_get_str(NULL, base, val.d_value);
                std::string res(str);
                mpfr_free_str(str); // There is no nice way to directly call the free() method of
                                    // GMP. But since MPFR is using the same allocators as GMP, we can
                                    // just use mpfr_free_str() instead...
                return res;
            }
            
            void to_string_conversion<Integer>::convert(std::string & res, const Integer & val)
            {
                char * str = mpz_get_str(NULL, 10, val.d_value);
                res = str;
                mpfr_free_str(str); // There is no nice way to directly call the free() method of
                                    // GMP. But since MPFR is using the same allocators as GMP, we can
                                    // just use mpfr_free_str() instead...
            }
            
            void to_string_conversion<Integer>::convert(std::string & res, const Integer & val, unsigned base)
            {
                char * str = mpz_get_str(NULL, base, val.d_value);
                res = str;
                mpfr_free_str(str); // There is no nice way to directly call the free() method of
                                    // GMP. But since MPFR is using the same allocators as GMP, we can
                                    // just use mpfr_free_str() instead...
            }
        }
        
        enum OutputType { OT_Default, OT_Fixed, OT_Scientific };
        enum BaseType { BT_Dec, BT_Hex, BT_Oct };
        enum AdjustType { AT_Internal, AT_Left, AT_Right };
        
        namespace
        {
            class CallMPFRFreeStr
            {
            private:
                char * d_str;
                
            public:
                CallMPFRFreeStr(char * str)
                    : d_str(str)
                {
                }
                
                ~CallMPFRFreeStr()
                {
                    mpfr_free_str(d_str);
                }
            };
        }
        
        void printToStream(std::ostream & s, const mpfr_t & r,
                           OutputType ot, BaseType bt, AdjustType at, unsigned width, unsigned precision,
                           char fill, bool showbase, bool showpoint, bool showpos, bool uppercase)
        // Output to stream
        {
            // Retrieve string and exponent
            mpfr_exp_t exp;
            int b = 10;
            switch (bt)
            {
            case BT_Dec: b = 10; break;
            case BT_Hex: b = 16; break;
            case BT_Oct: b =  8; break;
            }
            int digits = precision;
            if (ot == OT_Scientific)
                ++digits;
            if ((ot == OT_Fixed) && !mpfr_zero_p(r))
            {
                char * str = mpfr_get_str(NULL, &exp, b, 2, r, MPFR_RNDN);
                CallMPFRFreeStr strrel(str);
                digits = precision + exp;
                if (digits < 1)
                    digits = 1;
            }
            digits += 3; // just in case
            char * str = mpfr_get_str(NULL, &exp, b, digits, r, MPFR_RNDN);
            CallMPFRFreeStr strrel(str);
            char * rstr = str;
            ArrayStore<char> rstr_storage;
            if (str == NULL)
                // Error!
                return;
            // Handle sign and base
            char sign = 0;
            const char * base = NULL;
            if ((mpfr_sgn(r) < 0) || (*rstr == '-'))
            {
                if (*rstr == '-') // safe guard
                    ++rstr;
                sign = '-';
            }
            else
                if (showpos && (mpfr_sgn(r) > 0))
                    sign = '+';
            if (showbase)
                switch (bt)
                {
                case BT_Dec: break;
                case BT_Hex: base = uppercase ? "0X" : "0x"; break;
                case BT_Oct: base = "0"; break;
                }
            // Compute number of leading digits and length
            int leading = 1, len = std::strlen(rstr);
            if (mpfr_zero_p(r))
                exp = leading;
            // Shorten if necessary
            if (len > digits - 3)
            {
                len = digits - 3;
                rstr[len] = 0;
            }
            // Convert digits to uppercase if necessary
            if (uppercase && (bt == BT_Hex))
                for (int i = 0; i < len; ++i)
                    if (rstr[i] >= 'a')
                        rstr[i] -= 'a' - 'A';
            // Decide how many leading digits to display
            switch (ot)
            {
            case OT_Default:
                if ((exp > 0) && (exp <= len))
                    leading = exp;
                // Remove trailing zeros from the fractional part
                while (len > leading)
                    if (rstr[len - 1] == '0')
                        --len;
                    else
                        break;
                rstr[len] = 0;
                // Prepend zeros?
                if ((leading == 1) && (exp <= 0) && (len - exp - 3 < precision))
                {
                    // Insert -exp+1 zeros
                    for (int i = len; i >= 0; --i)
                        rstr[i + 1 - exp] = rstr[i];
                    for (unsigned i = 0; i < 1 - exp; ++i)
                        rstr[i] = '0';
                
                    // Increase len
                    len += 1-exp;
                
                    // Set exp to leading == 1
                    exp = 1;
                }
                break;
            case OT_Scientific:
                // Nothing to do.
                break;
            case OT_Fixed:
                // Shorten if necessary
                if ((exp > 0) && (exp + precision < len))
                    len = exp + precision;
            
                // Insert zeros if necessary
                if (exp <= 0)
                {
                    rstr_storage.reset(precision + 2);
                    for (unsigned i = 0; i <= precision; ++i)
                        rstr_storage[i] = '0';
                    rstr_storage[precision + 1] = 0;
                    for (unsigned i = -exp + 1; i <= precision; ++i)
                        if (i + exp <= len)
                            rstr_storage[i] = rstr[i + exp - 1];
                    rstr = rstr_storage.get();
                    len = precision + 1;
                    leading = 1;
                    exp = 1;
                }
            
                if ((exp > 0) && (exp < len))
                    leading = exp;
                break;
            }
            // Determine output length
            unsigned outlen = 0;
            if (sign) ++outlen;
            if (base) outlen += strlen(base);
            outlen += len;
            if ((leading < len) || (showpoint)) ++outlen;
            if ((ot == OT_Fixed) && ((unsigned)len < precision + (unsigned)leading))
                outlen += precision + leading - len;
            unsigned pot = 1;
            if ((exp - leading) || (ot == OT_Scientific))
            {
                outlen += 2;
                unsigned t = std::abs(exp - leading);
                do
                {
                    ++outlen;
                    t /= 10;
                    pot *= 10;
                } while (t);
                pot /= 10;
            }
            // Pad in front?
            if ((outlen < width) && (at == AT_Right))
            {
                for (unsigned i = outlen; i < width; ++i)
                    s.put(fill);
                outlen = width;
            }
            // Output sign/base
            if (sign)
            {
                s << sign;
            }
            if (base)
                s << base;
            // Pad in middle?
            if ((outlen < width) && (at == AT_Internal))
            {
                for (unsigned i = outlen; i < width; ++i)
                    s.put(fill);
                outlen = width;
            }
            // Output leading digits
            for (int i = 0; i < leading; ++i)
                s.put(rstr[i]);
            // Print decimal dot
            if ((leading < len) || (showpoint))
                s.put('.');
            // Is something left to show?
            for (int i = leading; i < len; ++i)
                s.put(rstr[i]);
            if (ot == OT_Fixed)
                for (unsigned i = len - leading; i < precision; ++i)
                    s.put('0');
            // Output exponent if needed
            if ((exp - leading) || (ot == OT_Scientific))
            {
                s.put(uppercase ? 'E' : 'e');
                s.put((exp >= leading) ? '+' : '-');
                unsigned e = std::abs(exp - leading);
                do
                {
                    s.put('0' + (e / pot));
                    e %= pot;
                    pot /= 10;
                } while (pot);
            }
            // Pad in end?
            if ((outlen < width) && (at == AT_Left))
            {
                for (unsigned i = outlen; i < width; ++i)
                    s.put(fill);
                outlen = width;
            }
        }
    
        std::ostream & operator << (std::ostream & s, const Real & r)
        // Output to stream
        {
            char fill = s.fill(); // fill character
            unsigned width = s.width(); // field width
            unsigned prec = s.precision(); // output precision (maximum number of digits)
            OutputType ot = OT_Default;
            switch (s.flags() & std::ios::floatfield)
            {
            case std::ios::scientific: ot = OT_Scientific; break;
            case std::ios::fixed: ot = OT_Fixed; break;
            default: ;
            }
            bool showbase = s.flags() & std::ios::showbase;
            bool showpoint = s.flags() & std::ios::showpoint;
            bool showpos = s.flags() & std::ios::showpos;
            bool uppercase = s.flags() & std::ios::uppercase;
            BaseType bt = BT_Dec;
            switch (s.flags() & std::ios::basefield)
            {
            case std::ios::dec: bt = BT_Dec; break;
            case std::ios::hex: bt = BT_Hex; break;
            case std::ios::oct: bt = BT_Oct; break;
            default: ;
            }
            AdjustType at = AT_Right;
            switch (s.flags() & std::ios::adjustfield)
            {
            case std::ios::internal: at = AT_Internal; break;
            case std::ios::left: at = AT_Left; break;
            case std::ios::right: at = AT_Right; break;
            default: ;
            }
            s.width(0); // reset width
            printToStream(s, r.d_value, ot, bt, at, width, prec, fill, showbase, showpoint, showpos, uppercase);
            return s;
        }
    
        std::istream & operator >> (std::istream & s, Real & r)
        // Input from stream
        {
            Buffer<char> buf;
            char c;
            c = s.get();
            if (!s)
            {
                s.setstate(std::ios_base::failbit);
                return s;
            }
            if ((c == '-') || (c == '+'))
            {
                buf.add(c);
                c = s.get();
                if (!s)
                {
                    s.setstate(std::ios_base::failbit);
                    return s;
                }
            }
            if (!isdigit(c) && (c != '.'))
            {
                s.setstate(std::ios_base::failbit);
                return s;
            }
            while (isdigit(c))
            {
                buf.add(c);
                c = s.get();
                if (!s)
                {
                    if (s.eof())
                    {
                        c = ' ';
                        break;
                    }
                    s.setstate(std::ios_base::failbit);
                    return s;
                }
            }
            if (c == '.')
            {
                buf.add(c);
                c = s.get();
                if (!s)
                {
                    s.setstate(std::ios_base::failbit);
                    return s;
                }
                while (isdigit(c))
                {
                    buf.add(c);
                    c = s.get();
                    if (!s)
                    {
                        if (s.eof())
                        {
                            c = ' ';
                            break;
                        }
                        s.setstate(std::ios_base::failbit);
                        return s;
                    }
                }
            }
            if ((c == 'e') || (c == 'E'))
            {
                buf.add(c);
                c = s.get();
                if (!s)
                {
                    s.setstate(std::ios_base::failbit);
                    return s;
                }
                if ((c == '-') || (c == '+'))
                {
                    buf.add(c);
                    c = s.get();
                    if (!s)
                    {
                        s.setstate(std::ios_base::failbit);
                        return s;
                    }
                }
                if (!isdigit(c))
                {
                    s.setstate(std::ios_base::failbit);
                    return s;
                }
                while (isdigit(c))
                {
                    buf.add(c);
                    c = s.get();
                    if (!s)
                    {
                        if (s.eof())
                        {
                            c = ' ';
                            break;
                        }
                        s.setstate(std::ios_base::failbit);
                        return s;
                    }
                }
            }
            if (!s.eof())
                s.unget();
            else
            {
                if (!s.bad())
                    s.clear(std::ios_base::eofbit);
            }
            buf.add(0);
            char * s_end;
            mpfr_strtofr(r.d_value, buf.data(), &s_end, 10, MPFR_RNDN);
            if (*s_end != 0)
                s.setstate(std::ios_base::failbit);
            return s;
        }
        
        namespace implementation
        {
            bool from_string_conversion<RealContext>::convert(Real & res, const std::string & s, const RealContext & c)
            {
                return convert(res, s.c_str(), c);
            }
            
            bool from_string_conversion<RealContext>::convert(Real & res, const char * s, const RealContext & c)
            {
                if (*s == 0)
                    return false;
                char * s_end;
                res.setContext(c);
                mpfr_strtofr(res.d_value, s, &s_end, 10, MPFR_RNDN);
                return s + strlen(s) == s_end;
            }
            
            std::string to_string_conversion<Real>::convert(const Real & val)
            {
                mpfr_exp_t exp;
                char * str = mpfr_get_str(NULL, &exp, 10, 0, val.d_value, MPFR_RNDN);
                std::string s;
                if (sign(val) < 0)
                    s += '-';
                s += "0.";
                s += (sign(val) < 0 ? str + 1 : str);
                mpfr_free_str(str);
                const int expstrsize = 20;
                char expstr[expstrsize];
                snprintf(expstr, expstrsize, "e%li", exp);
                s += expstr;
                return s;
            }
            
            std::string to_string_conversion<Real>::convert(const Real & val, unsigned base)
            {
                mpfr_exp_t exp;
                char * str = mpfr_get_str(NULL, &exp, base, 0, val.d_value, MPFR_RNDN);
                std::string s;
                if (sign(val) < 0)
                    s += '-';
                s += "0.";
                s += (sign(val) < 0 ? str + 1 : str);
                mpfr_free_str(str);
                const int expstrsize = 20;
                char expstr[expstrsize];
                snprintf(expstr, expstrsize, "e%li", exp);
                s += expstr;
                return s;
            }
            
            void to_string_conversion<Real>::convert(std::string & res, const Real & val)
            {
                mpfr_exp_t exp;
                char * str = mpfr_get_str(NULL, &exp, 10, 0, val.d_value, MPFR_RNDN);
                res.clear();
                if (sign(val) < 0)
                    res += '-';
                res += "0.";
                res += (sign(val) < 0 ? str + 1 : str);
                mpfr_free_str(str);
                const int expstrsize = 20;
                char expstr[expstrsize];
                snprintf(expstr, expstrsize, "e%li", exp);
                res += expstr;
            }
            
            void to_string_conversion<Real>::convert(std::string & res, const Real & val, unsigned base)
            {
                mpfr_exp_t exp;
                char * str = mpfr_get_str(NULL, &exp, base, 0, val.d_value, MPFR_RNDN);
                res.clear();
                if (sign(val) < 0)
                    res += '-';
                res += "0.";
                res += (sign(val) < 0 ? str + 1 : str);
                mpfr_free_str(str);
                const int expstrsize = 20;
                char expstr[expstrsize];
                snprintf(expstr, expstrsize, "e%li", exp);
                res += expstr;
            }
        }
        
        void Integer::mpz_set_ll(mpz_t & value, long long i)
        {
            PLLL_INTERNAL_STATIC_CHECK(sizeof(long long) == 8, Implementation_requires_longlong_to_be_64bits);
        
            bool neg = i < 0;
            unsigned long long u = neg ? -i : i;
            mpz_set_ui(value, (unsigned long)(u >> 32));
            mpz_mul_2exp(value, value, 32);
            mpz_add_ui(value, value, (unsigned long)(u & 0xFFFFFFFFull));
        }
        
        long long Integer::mpz_get_ll(const mpz_t & value)
        {
            unsigned long long u;
#if GMP_LIMB_BITS == 64
            u = mpz_getlimbn(value, 0);
#else
  #error Implementation currently only supports limb size of 64!
#endif
            return (mpz_sgn(value) < 0) ? -u : u;
        }
        
        void Real::mpfr_set_ll(mpfr_t & v, long long i, mpfr_rnd_t rnd)
        {
#if _MPFR_H_HAVE_INTMAX_T
            PLLL_INTERNAL_STATIC_CHECK(sizeof(long long) <= sizeof(intmax_t), Implementation_requires_longlong_to_be_64bits);
            mpfr_set_sj(v, i, rnd);
#else
            PLLL_INTERNAL_STATIC_CHECK(sizeof(long long) == 8, Implementation_requires_longlong_to_be_64bits);
            bool neg = i < 0;
            unsigned long long u = neg ? -i : i;
            mpfr_set_ui(v, (unsigned long)(u >> 32), rnd);
            mpfr_mul_2ui(v, v, 32, rnd);
            mpfr_add_ui(v, v, (unsigned long)(u & 0xFFFFFFFFull), rnd);
#endif
        }
        
        long long Real::mpfr_get_ll(const mpfr_t & v, mpfr_rnd_t rnd)
        {
#if _MPFR_H_HAVE_INTMAX_T
            PLLL_INTERNAL_STATIC_CHECK(sizeof(long long) <= sizeof(intmax_t), Implementation_requires_longlong_to_be_64bits);
            return mpfr_get_sj(v, rnd);
#else
            bool neg = mpfr_signbit(v);
            PLLL_INTERNAL_STATIC_CHECK(sizeof(long long) == 8, Implementation_requires_longlong_to_be_64bits);
            mpfr_t t;
            mpfr_prec_t p = mpfr_get_prec(v);
            mpfr_init2(t, p);
            mpfr_rint(t, v, rnd);
            mpfr_abs(t, t, MPFR_RNDN);
            unsigned long long res;
            unsigned long tt;
            mpfr_div_2ui(t, t, 32, MPFR_RNDN);
            tt = mpfr_get_ui(t, MPFR_RNDZ);
            mpfr_sub_ui(t, t, tt, MPFR_RNDN);
            mpfr_mul_2ui(t, t, 32, MPFR_RNDN);
            res = (((unsigned long long)tt) << 32) + (unsigned long long)mpfr_get_ui(t, MPFR_RNDZ);
            mpfr_clear(t);
            return neg ? -res : res;
#endif
        }
        
        long long Real::mpfr_get_ll(const mpfr_t & v, mpfr_rnd_t rnd, bool & roundUp)
        {
#if _MPFR_H_HAVE_INTMAX_T
            PLLL_INTERNAL_STATIC_CHECK(sizeof(long long) <= sizeof(intmax_t), Implementation_requires_longlong_to_be_64bits);
            long long w = mpfr_get_sj(v, rnd);
            roundUp = (rnd == MPFR_RNDD) ? false : (w > mpfr_get_sj(v, MPFR_RNDD));
            return w;
#else
            bool neg = mpfr_signbit(v);
            PLLL_INTERNAL_STATIC_CHECK(sizeof(long long) == 8, Implementation_requires_longlong_to_be_64bits);
            mpfr_t t;
            mpfr_prec_t p = mpfr_get_prec(v);
            mpfr_init2(t, p);
            mpfr_rint(t, v, rnd);
            roundUp = mpfr_cmp(v, t) < 0;
            mpfr_abs(t, t, MPFR_RNDN);
            unsigned long long res;
            unsigned long tt;
            mpfr_div_2ui(t, t, 32, MPFR_RNDN);
            tt = mpfr_get_ui(t, MPFR_RNDZ);
            mpfr_sub_ui(t, t, tt, MPFR_RNDN);
            mpfr_mul_2ui(t, t, 32, MPFR_RNDN);
            res = (((unsigned long long)tt) << 32) + (unsigned long long)mpfr_get_ui(t, MPFR_RNDZ);
            mpfr_clear(t);
            return neg ? -res : res;
#endif
        }
    }
}
