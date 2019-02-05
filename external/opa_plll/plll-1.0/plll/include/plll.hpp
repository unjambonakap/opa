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

#ifndef PLLL_INCLUDE_GUARD__PLLL2_LLL_HPP
#define PLLL_INCLUDE_GUARD__PLLL2_LLL_HPP

/**
   \file
   \brief Main header for the `plll` library.
   
   This is the main header for the `plll` lattice reduction library. It is usually the only header
   you have to include to use the library. It provides the main namespace, [plll](@ref plll),
   which contains all definitions.
   
   Most prominently, it contains the class [LatticeReduction](@ref plll::LatticeReduction), which
   provides the main interface to all lattice reduction algorithms.
*/

/**
   \dir /plll/include
   \brief Contains the main `plll` header `plll.hpp`.
 */

/**
   \dir /plll/include/plll
   \brief Contains individual public headers for the `plll` library.
 */

#include <plll/config.hpp>
#include <plll/arithmetic.hpp>
#include <plll/arithmetic-nint.hpp>
#include <plll/matrix.hpp>
#include <exception>
#include <iostream>
#include <limits>
#include <boost/function.hpp>

/**
   \brief Contains the `plll` library.
   
   The `plll` namespace contains all definitions of the `plll` library.
*/
namespace plll
{
    class LatticeReductionImpl;
    
    class LatticeReduction
    /**
       \brief Provides an interface to the lattice reduction algorithms of the `plll` library.
       
       This class provides an interface to all lattice reduction algorithms of the `plll`
       library. It is possible to configure the algorithms using the interface and execute them. The
       lattice can be set and retrieved, and also callback functions can be set up.
    */
    {
    private:
        std::auto_ptr<LatticeReductionImpl> d_impl; /// A pointer to the implementation. Used to hide internals
                                                    /// of the implementation from the user.
        
        // Disable copy constructor and operator
        LatticeReduction(const LatticeReduction &);
        LatticeReduction & operator = (const LatticeReduction &);
        
    public:
        /**@{
           \name Exceptions.
         */
        class stop_reduction : public std::exception
        /**
           \brief This exception should be thrown by callback functions to stop reduction and return
                  the current state of reduction to the caller.
           
           \see \ref desc-genconf-callback
         */
        {
        public:
            virtual ~stop_reduction() PLLL_INTERNAL_NOTHROW_POSTFIX_ENFORCE
            {
            }
            
            virtual const char * what() const PLLL_INTERNAL_NOTHROW_POSTFIX_ENFORCE;
        };
        
        class stop_enumeration : public std::exception
        /**
           \brief This exception should be thrown by callback functions to stop enumeration and
                  return the currently found vector to the caller.
           
           \see \ref desc-genconf-callback
         */
        {
        public:
            virtual ~stop_enumeration() PLLL_INTERNAL_NOTHROW_POSTFIX_ENFORCE
            {
            }
            
            virtual const char * what() const PLLL_INTERNAL_NOTHROW_POSTFIX_ENFORCE;
        };
        ///@}
        
        /**@{
           \name Configuration enums.
         */
        enum Arithmetic
        /**
           \brief Specifies which arithmetic to use for Gram-Schmidt orthogonalization.
           
           Specifies which arithmetic should be used by the `plll` library for Gram-Schmidt
           orthogonalizations. The choice of arithmetics might vary depending on the platform, but
           except for "Fast Compile" scenarios, at least rational, real, long double and double
           arithmetic is supported.
           
           \see \ref desc-genconf-arith
        */
        {
            A_Rational,
            /**< Rational arithmetic, i.e. a pair of `arithmetic::Integer` objects is used for numerator
                 and denominator. While this yields infinite precision, it is also quite slow.
            */
            A_Real,
            /**< MPFR-based `arithmetic::Real` arithmetic. This is a floating point type with arbitrary
                 (but fixed) precision. A precision is automatically selected by the library.
            */
            A_LongDouble,
            /**< The system's native `long double` arithmetic. This is a finite-precision floating
                 point type with implementation defined precision. On Intel/AMD/... platforms, it
                 usually provides a mantissa of 80 bits, while on UltraSparc platforms, this is a
                 Float128 type with 112 bit mantissa, whose arithmetic is unfortunately usually
                 implemented in software and not in hardware.
                 
                 This is the default setting for `plll`.
            */
            A_Double,
            /**< The system's native `double` arithmetic. This is a finite-precision floating point
                 type with implementation defined precision. It usually has a mantissa of 53 bits.
            */
            A_DoubleDouble,
            /**< Uses two native `double` floating point numbers to represent floating point numbers
                 with higher precision. While the exponent is still limited to a `double`'s range,
                 the precision is much higher; double double numbers correspond to floating point
                 numbers with a mantissa of 104 bits.
                 
                 Support for this arithmetic depends on the libqd library by Yozo Hida. It
                 apparently does not work on UltraSparc CPUs.
            */
            A_QuadDouble,
            /**< Uses four native `double` floating point numbers to represent floating point
                 numbers with higher precision. While the exponent is still limited to a `double`'s
                 range, the precision is much higher; quad double numbers correspond to floating
                 point numbers with a mantissa of 209 bits.
                 
                 Support for this arithmetic depends on the libqd library by Yozo Hida. It
                 apparently does not work on UltraSparc CPUs.
            */
            A_Default = A_LongDouble
            /**< Alias for the default arithmetic, i.e. for `A_LongDouble`. */
        };
        
        enum Integers
        /**
           \brief Specifies which arithmetic to use for integer operations.
           
           Specifies which arithmetic should be used by the `plll` library for integer
           operations. Besides GMP-based arbitrary precision integers, it also provides support for
           system-dependent `long int` integers and for an automatic decision mode which uses
           arbitrary precision integers if necessary and system integers otherwise.
           
           \see \ref desc-genconf-arith
         */
        {
            I_Auto,
            /**< Decides based on the input whether to use arbitrary precision integer arithmetic or
                 system-dependent integer arithmetic. This decision is analyzed from time to time,
                 and arithmetic is changed when necessary.
                 
                 This is the default setting for `plll`.
            */
            I_ArbitraryPrecision,
            /**< Uses `arithmetic::Integer` arbitrary precision integers. These are provided by GMP and
                 are only limited by the system's available RAM.
            */
            I_LongInt,
            /**< Uses `long int` native CPU arithmetic. On modern CPUs, this yields at least 64 bits
                 of precision.
            */
            I_Default = I_Auto
            /**< Alias for the default integer arithmetic, i.e. for `A_Auto`. */
        };
        
        enum GramSchmidt
        /**
           \brief Specifies which Gram-Schmidt orthogonalization is used.
           
           Specifies which Gram-Schmidt orthogonalization should be used by the `plll`
           library. Always available are classic GS, classic GS with integer arithmetic, and
           numerically stable GS. For floating-point arithmetic, also Givens rotations can be used
           for GS orthogonalization.
           
           \see \ref desc-genconf-gs
         */
        {
            G_Classic,
            /**< "Classical" Gram-Schmidt orthogonalization. Uses formulae to update Gram-Schmidt
                 coefficients on swaps, transformations etc., which can introduce additional error
                 with floating point arithmetic.
            */
            G_ClassicInteger,
            /**< Integer-based Gram-Schmidt orthogonalization. Internally uses arbitrary precision
                 integers respectively rational numbers. Slow, but always accurate.
            */
            G_Givens,
            /**< Uses Givens rotations to compute Gram-Schmidt orthogonalization. Only available for
                 floating-point arithmetic (i.e. everything but `A_Rational`). Uses formulae to
                 update Gram-Schmidt coefficients on swaps, transformations etc., which can
                 introduce additional error with floating point arithmetic.
            */
            G_NumStable,
            /**< Uses numerically stable Gram-Schmidt orthogonalization as described by P. Q. Nguyen
                 and D. Stehle in
                 \cite nguyen-stehle-fplll-revisited. This arithmetic is more resilient against
                 approximation errors than other approaches. Instead of updating GS coefficients for
                 swaps, transformations etc., these are recomputed to ensure maximal precision.
                 
                 This is the default setting for `plll`.
            */
            G_Default = G_NumStable
            /**< Alias for the default Gram-Schmidt orthogonalization, i.e. for `A_NumStable`. */
        };
        
        enum Transform
        /**
           \brief Specifies which transformation matrix is computed.
           
           Specifies which kind of transformation matrix is computed by the `plll` library, if
           computation of transformation matrices is requested.
         */
        {
            T_Normal,
            /**< Computes a transformation matrix \f$T_1\f$ such that if \f$A\f$ is the old and
                 \f$A'\f$ the new matrix representing the lattice before respectively after
                 reduction, then \f$T_1 \cdot A = A'\f$.
            */
            T_Inverse
            /**< Computes a transformation matrix \f$T_2\f$ such that if \f$A\f$ is the old and
                 \f$A'\f$ the new matrix representing the lattice before respectively after
                 reduction, then \f$T_2 \cdot A' = A\f$. Note that if the matrices are square, then
                 both are invertible and \f$T_2 = T_1^{-1}\f$.
            */
        };
        
        enum LLLMode
        /**
           \brief Specifies which LLL condition is used.
           
           Specifies which LLL condition is used by the `plll` library. One can choose between the
           original condition introduced by Lenstra, Lenstra and Lovász, a "unprojected" variant
           comparing the norms in the ambient space, and the Siegel condition.
         */
        {
            LLL_Classic,
            /**< The original condition by Lenstra, Lenstra and Lovász, called the *Lovász
                 condition*, which compares the projected lengths of two adjacent vectors, where the
                 vectors are projected onto the orthogonal complement of all previous vectors.
                 
                 In terms of the Gram-Schmidt orthogonalization \f$\pi_1(b_1), \dots, \pi_k(b_k)\f$
                 of the basis \f$b_1, \dots, b_k\f$, this condition can be expressed as \f$\alpha
                 \cdot \|\pi_k(b_k)\|^2 > \|\pi_k(b_{k+1})\|^2\f$.
                 
                 This is the default setting for `plll`.
                 
                 \see \ref desc-algs-lll-classic
            */
            LLL_Unprojected,
            /**< A simpler condition, which compares the unprojected lengths of two adjacent.

                 In terms of the basis \f$b_1, \dots, b_k\f$, this condition can be expressed as
                 \f$\alpha \cdot \|b_k\|^2 > \|b_{k+1}\|^2\f$.
                 
                 \see \ref desc-algs-lll-unproj
            */
            LLL_Siegel,
            /**< A simpler condition, which allows to prove the same bounds on the output quality,
                 called the *Siegel condition*. Note that the Lovász condition implies the Siegel
                 condition.
                 
                 In terms of the Gram-Schmidt orthogonalization \f$\pi_1(b_1), \dots, \pi_k(b_k)\f$
                 of the basis \f$b_1, \dots, b_k\f$, this condition can be expressed as
                 \f$\frac{3}{4} \cdot \alpha \cdot \|\pi_k(b_k)\|^2 > \|\pi_{k+1}(b_{k+1})\|^2\f$.
                 
                 \see \ref desc-algs-lll-siegel
            */
            LLL_Default = LLL_Classic
            /**< Alias for the default LLL condition, i.e. for `LLL_Classic`. */
        };
        
        enum BKZMode
        /**
           \brief Specifies which BKZ algorithm is used.
           
           Specifies which BKZ algorithm is used by the `plll` library. Note that BKZ is somewhat
           misleading here, since most of these algorithms are known under completely different
           names. What all these algorithms have in common that they compute Korkine-Zolotarev or at
           least SVP bases for certain blocks of (projected) basis vectors, or their dual.
         */
        {
            BKZ_SchnorrEuchner,
            /**< The classical Schnorr-Euchner BKZ algorithm, as described in Section 6 of
                 \cite schnorr-euchner-BKZ by C.-P. Schnorr and M. Euchner.
                 
                 This is the default setting for `plll`.
                 
                 \see \ref dest-algs-bkz-schnorreuchner
            */
            BKZ_Simplified,
            /**< A simplified version of the BKZ algorithm, as described by G. Hanrot, X. Pujol and
                 D. Stehle in
                 \cite hanrot-pujol-stehle-BKZdynamic.
                 
                 \see \ref dest-algs-bkz-simplified
            */
            BKZ_HanrotPujolStehleHKZ,
            /**< A "Terminating BKZ" variant of the simplified version of the BKZ algorithm, as
                 described by G. Hanrot, X. Pujol and D. Stehle in
                 \cite hanrot-pujol-stehle-BKZdynamic. This version computes
                 Hermite-Korkine-Zolotarev (HKZ) bases for every block.
                 
                 \see \ref dest-algs-bkz-terminating
            */
            BKZ_HanrotPujolStehleSVP,
            /**< A "Terminating BKZ" variant of the simplified version of the BKZ algorithm, as
                 described by G. Hanrot, X. Pujol and D. Stehle in
                 \cite hanrot-pujol-stehle-BKZdynamic. This version computes Shortest Vector (SVP)
                 bases for every block.
                 
                 \see \ref dest-algs-bkz-terminating
            */
            BKZ_PrimalDual,
            /**< The primal-dual BKZ algorithm by H. Koy, as described in Section 3 of
                 \cite schnorr-bkzrevisited
                 
                 \see \ref dest-algs-bkz-primaldual
            */
            BKZ_SlideReduction,
            /**< The slide reduction algorithm, as described by N. Gama and P. Q. Nguyen in
                 \cite gama-nguyen-mordellsinequality.
                 
                 \see \ref dest-algs-bkz-slide
            */
            BKZ_ImprovedSlideReduction,
            /**< A improved and accelerated version of Slide Reduction, as described in
                 \cite schnorr-accslide by C.-P. Schnorr (Section 3).
                 
                 \see \ref dest-algs-bkz-improvslide
            */
            BKZ_ImprovedSlideReduction2,
            /**< A improved and accelerated version of Slide Reduction, as described in
                 \cite schnorr-accslide by C.-P. Schnorr (Section 3).
                 
                 This variant uses a larger dual SVP enumeration and a single primal SVP
                 enumeration.
                 
                 \see \ref dest-algs-bkz-improvslide
            */
            BKZ_ImprovedSlideReduction3,
            /**< A improved and accelerated version of Slide Reduction, as described in
                 \cite schnorr-accslide by C.-P. Schnorr (Section 3).
                 
                 This variant uses a single larger primal SVP enumeration and one dual SVP
                 enumeration.
                 
                 \see \ref dest-algs-bkz-improvslide
            */
            BKZ_SemiBlock2k,
            /**< The semi block 2k-reduction, as described by C.-P. Schnorr in Section 2 of
                 \cite schnorr-bkzrevisited.
                 
                 \see \ref dest-algs-bkz-semiblock2k
            */
            BKZ_SamplingReduction,
            /**< The Sampling Reduction, as described by J. A. Buchmann and C. Ludwig in
                 \cite buchmann-ludwig-samplingreduction.
                 
                 \see \ref dest-algs-bkz-sr
            */
            BKZ_Experimental, // An experimental variant which should *not* be used.
            BKZ_Default = BKZ_SchnorrEuchner
            /**< Alias for the default BKZ variant, i.e. for `BKZ_SchnorrEuchner`. */
        };
        
        enum SVPMode
        /**
           \brief Specifies which SVP solver is used.
           
           Specifies which SVP solver is used by the `plll` library. For small enough dimensions,
           always a simple (but efficient) enumeration variant is used, except for the highest level
           reduction, where always the selected solver is used.
         */
        {
            SVP_KannanSchnorrEuchner,
            /**< A variant of the Kannan-Schnorr-Euchner enumeration algorithm
                 \cite schnorr-euchner-BKZ with various improvements; see, for example, Appendix B of
                 \cite gama-nguyen-regev-extremepruning.
                 
                 \see \ref desc-algs-svp-kse
            */
            SVP_ParallelKannanSchnorrEuchner,
            /**< A parallelized variant of the Kannan-Schnorr-Euchner enumeration algorithm
                 \cite schnorr-euchner-BKZ with various improvements. The parallelization uses
                 multiple cores. There is almost no communication between cores, except if a shorter
                 vector is found by one, or if one core runs out of work and asks others to split
                 their workload.
                 
                 This is the default setting for `plll`. Note that in case only one core is
                 available, automatically the solver `SVP_KannanSchnorrEuchner` will be used.
                 
                 \see \ref desc-algs-svp-kse
            */
            SVP_SchnorrFast,
            /**< A (single threaded) version of Schnorr's new SVP solver, as described in
                 \cite schnorr-factoringcvp by C.-P. Schnorr.
                 
                 This implementation is still experimental.
                 
                 \see \ref desc-algs-svp-schnorr
            */
            SVP_VoronoiCellSVP,
            /**< A deterministic SVP solver which first computes the voronoi cell of the
                 lattice. While being deterministic and asymptotically faster than enumeration, in
                 practice this algorithm is only usable for very low dimensions, say at most 10,
                 where enumeration is still much faster. This algorithm was described by
                 D. Micciancio and P. Voulgaris in
                 \cite micciancio-voulgaris-voronoicell.
                 
                 \see \ref desc-algs-svp-voronoi
            */
            SVP_ListSieve,
            /**< The probabilistic list sieve SVP solver. It was described by D. Micciancio and
                 P. Voulgaris in Section 3.1 of
                 \cite micciancio-voulgaris-fastersvp.
                 
                 This implementation is still experimental.
                 
                 \see \ref desc-algs-svp-list
            */
            SVP_ListSieveBirthday,
            /**< The probabilistic list sieve (with birthday paradox exploitation) SVP solver. It
                 was described by X. Pujol and D. Stehle in the paper
                 \cite pujol-stehle-birthdaysieve.
                 
                 This implementation is still experimental.
                 
                 \see \ref desc-algs-svp-lsb
            */
            SVP_GaussSieve,
            /**< The probabilistic Gauss sieve SVP solver. It was described by D. Micciancio and
                 P. Voulgaris in Section 3.2 of
                 \cite micciancio-voulgaris-fastersvp. Our implementation is similar to <a
                 href="http://cseweb.ucsd.edu/~pvoulgar/impl.html">one by P. Voulgaris</a>.
                 
                 \see \ref desc-algs-svp-gauss
            */
            SVP_Default = SVP_ParallelKannanSchnorrEuchner
            /**< Alias for the default SVP solver, i.e. for `SVP_ParallelKannanSchnorrEuchner`. */
        };
        
        enum DIMethod
        /**
           \brief Specifies which Deep Insertion mode is used.
           
           Specifies which Deep Insertion mode is used by the `plll` library. Note that Deep
           Insertions are not supported directly by all LLL and BKZ algorithms, but will be used by
           recursive LLL and BKZ calls.
           
           The idea of Deep Insertions is that while usual LLL only swaps adjacent vectors, Deep
           Insertions also compares the current vectors with more previous vectors to decide where
           it could be inserted.
         */
        {
            DI_None,
            /**< Disable Deep Insertions.
                 
                 This is the default setting for `plll`.
            */
            DI_Classic,
            /**< Uses classic Deep Insertions as described by C.-P. Schnorr and M. Euchner in
                 \cite schnorr-euchner-BKZ.
            */
            DI_MinimizePotential1,
            /**< Uses potential minimizing Deep Insertions (PotLLL) as described by F. Fontein,
                 M. Schneider and U. Wagner in
                 \cite fontein-schneider-wagner-minpot-wcc and
                 \cite fontein-schneider-wagner-potlll.
            */
            DI_MinimizePotential2,
            /**< Uses a variant of the potential minimizing Deep Insertions (PotLLL2) as described
                 by F. Fontein, M. Schneider and U. Wagner in
                 \cite fontein-schneider-wagner-potlll. See "PotLLL2" in Section 3.1.
            */
            DI_Default = DI_None
            /**< Alias for the default Deep Insertions mode, i.e. for `DI_None`. */
        };
        
        enum DIChoice
        /**
           \brief Specifies which Deep Insertion choice is used.
           
           Specifies which vectors are considered for Deep Insertions by the `plll` library.
         */
        {
            DIC_First,
            /**< Deep Insertions considers only the first t vectors of the basis, where t is given.
            */
            DIC_Block,
            /**< Deep Insertions considers only the block of t basis vectors ending at the current
                 vector, where t is given.
            */
            DIC_FirstBlock,
            /**< Deep Insertions considers both the first t basis vectors, and the block of t basis
                 vectors ending at the current vector, where t is given.
                 
                 This is the default setting for the LLL reduction in V. Shoup's Number Theory
                 Library, if Deep Insertions are enabled.
            */
            DIC_All,
            /**< Deep Insertions considers all basis vectors as possible destinations.
                 
                 This is the default setting for `plll`.
            */
            DIC_Default = DIC_All
            /**< Alias for the default Deep Insertions choice, i.e. for `DIC_All`. */
        };
        
        enum DIMode
        /**
           \brief Specifies which Deep Insertions mode is used.
           
           Specifies which kind Deep Insertions mode is used, i.e. whether the insertions are done
           before and/or after size reduction.
         */
        {
            DIM_BeforeSR,
            /**< Do Deep Insertions only before size reduction.
            */
            DIM_AfterSR,
            /**< Do Deep Insertions only after size reduction.
            */
            DIM_Both,
            /**< Do Deep Insertions both before and after size reduction.
                 
                 This is the default setting for `plll`.
            */
            DIM_Default = DIM_Both
            /**< Alias for the default Deep Insertions mode, i.e. for `DIM_Both`. */
        };
        
        enum VerboseOutputLevel
        /**
            \brief Specifies the verbosity output level.
            
            Specifies the verbosity output level of the `plll` library. It can be configured
            between no output, and full output. Note that by outputting, it is meant that the
            specified `VerboseFunction` is called (or its default, which uses `std::cerr`).
            
            \see \ref desc-genconf-verbose
        */
        {
            VOL_None,
            /**< Outputs nothing. */
            VOL_Warnings,
            /**< Outputs only warnings and errors. */
            VOL_Informative,
            /**< Outputs warnings, errors and informative messages. */
            VOL_Full
            /**< Outputs everything. */
        };
        
        enum VerboseLevel
        /**
            \brief Specifies a warning level for a specific output message.
            
            Specifies the level of a specific output message. This is used to inform a
            `VerboseFunction` callback function which level the given message has.
            
            \see \ref desc-genconf-verbose
        */
        {
            VL_Error,
            /**< The given message is an error. */
            VL_Warning,
            /**< The given message is a warning. */
            VL_Information,
            /**< The given message is purely informative. */
            VL_Chatter
            /**< The given message is chatter which can safely be ignored. Can be helpful for
                 debugging purposes or to see more about how the algorithms work. */
        };
        ///@}
        
        struct Statistics
        /**
           \brief Describes lattice basis reduction statistics.
           
           Describes statistics collected while reducting lattice bases. These statistics count the
           numer of operations done.
         */
        {
            unsigned long swaps;
            /**< Counts the number of times two adjacent vectors are swapped.
            */
            unsigned long adds;
            /**< Counts the number of times a non-zero multiple of a basis vector was added to
                 another basis vector.
            */
            unsigned long adds_pm1;
            /**< Counts the number of times a basis vector was added or subtracted from another
                 basis vector. This number is never larger than `adds`.
            */
            unsigned long adds_pm2;
            /**< Counts the number of times a basis vector was added or subtracted twice from
                 another basis vector. This number is never larger than `adds`.
            */
            unsigned long flips;
            /**< Counts the number of times a basis vector was flipped, i.e. all its signs where
                 changed.
            */
            unsigned long trans;
            /**< Counts the number of times a 2 x 2 invertible transformation was applied to two
                 basis vectors.
            */
            unsigned long sizereductions;
            /**< Counts the number of applied size reductions.
            */
            unsigned long deepinsertions;
            /**< Counts the number of Deep Insertions with insertion distance at least 2. (An
                 insertion distance of 1 is counted by `swaps`.)
            */
            unsigned long enumcalls;
            /**< Counts the number of calls to the SVP solver.
            */
            unsigned long enumfails;
            /**< Counts the number of calls to the SVP solver where the solver did not find any
                 vector at all. (This happens if heuristics for the enumeration radius are used
                 which are too small.)
            */
            unsigned long vectorinsertions;
            /**< Counts the number of insertions of (completely) new vectors. This can be a result
                 of enumeration or of Sampling Reduction.
            */
            unsigned long vectorinsertions_rearrange;
            /**< Counts the number of where the basis is rearranged to move one basis vector to an
                 earlier position. This is usually the result of an SVP solver returning a vector
                 which already is part of the basis.
            */
            
            void reset(); ///< \brief Resets the statistics.
            Statistics(); ///< \brief Creates a new statistics object with all values set to zero.
        };
        
        class GramSchmidtInformer
        /**
           \brief Provides information on the Gram-Schmidt coefficients.
           
           An interface for a Gram-Schmidt Informer, which is a little object providing access to
           the Gram-Schmidt coefficients. This is for example provided to annealing callbacks.
        */
        {
        protected:
            GramSchmidtInformer() { }
            virtual ~GramSchmidtInformer() { }
            
        public:
            virtual double getGSCoefficientD(unsigned, unsigned) const = 0;
            /**< \brief Returns the `(i, j)` Gram-Schmidt coefficient as a `double`. */
            virtual long double getGSCoefficientLD(unsigned, unsigned) const = 0;
            /**< \brief Returns the `(i, j)` Gram-Schmidt coefficient as a `long double`. */
            virtual arithmetic::Real getGSCoefficientR(unsigned, unsigned, const arithmetic::RealContext &) const = 0;
            /**< \brief Returns the `(i, j)` Gram-Schmidt coefficient as a `arithmetic::Real`
                        number. */
            
            virtual double getGSSqNormD(unsigned) const = 0;
            /**< \brief Returns the squared norm of the `i`-th Gram-Schmidt orthogonalized basis
                        vector as a `double`.
            */
            virtual long double getGSSqNormLD(unsigned) const = 0;
            /**< \brief Returns the squared norm of the `i`-th Gram-Schmidt orthogonalized basis
                        vector as a `long double`.
            */
            virtual arithmetic::Real getGSSqNormR(unsigned, const arithmetic::RealContext &) const = 0;
            /**< \brief Returns the squared norm of the `i`-th Gram-Schmidt orthogonalized basis
                        vector as a `arithmetic::Real` number.
            */
            
            virtual double computeProjectionLengthD(unsigned k, unsigned b,
                                                    const linalg::math_rowvector<arithmetic::Integer> & vec) const = 0;
            /**< \brief Returns the squared norm of the given vector projected into the orthogonal
                        complement of the first `k`-th basis vectors.
                 
                 The given vector is a linear combination of basis vectors `b`, `b + 1`, ..., where
                 the linear combination is provided in the variable `vec`. The result is returned as
                 a `double`.
            */
            virtual long double computeProjectionLengthLD(unsigned k, unsigned b,
                                                          const linalg::math_rowvector<arithmetic::Integer> & vec) const = 0;
            /**< \brief Returns the squared norm of the given vector projected into the orthogonal
                        complement of the first `k`-th basis vectors.
                 
                 The given vector is a linear combination of basis vectors `b`, `b + 1`, ..., where
                 the linear combination is provided in the variable `vec`. The result is returned as
                 a `long double`.
            */
            virtual arithmetic::Real computeProjectionLengthR(unsigned k, unsigned b,
                                                              const linalg::math_rowvector<arithmetic::Integer> & vec,
                                                              const arithmetic::RealContext &) const = 0;
            /**< \brief Returns the squared norm of the given vector projected into the orthogonal
                        complement of the first `k`-th basis vectors.
                 
                 The given vector is a linear combination of basis vectors `b`, `b + 1`, ..., where
                 the linear combination is provided in the variable `vec`. The result is returned as
                 a `arithmetic::Real` number.
            */
            
            inline arithmetic::Real getGSCoefficientR(unsigned i, unsigned j) const
            /**< \brief Returns the `(i, j)` Gram-Schmidt coefficient as a `arithmetic::Real`
                        number.
                 
                 Uses the default `arithmetic::RealContext`.
            */
            {
                return getGSCoefficientR(i, j, arithmetic::getThreadRealContext());
            }
            
            inline arithmetic::Real getGSSqNormR(unsigned i) const
            /**< \brief Returns the squared norm of the `i`-th Gram-Schmidt orthogonalized basis
                        vector as a `arithmetic::Real` number.
                 
                 Uses the default `arithmetic::RealContext`.
            */
            {
                return getGSSqNormR(i, arithmetic::getThreadRealContext());
            }
            
            inline arithmetic::Real computeProjectionLengthR(unsigned k, unsigned b,
                                                             const linalg::math_rowvector<arithmetic::Integer> & vec) const
            /**< \brief Returns the squared norm of the given vector projected into the orthogonal
                        complement of the first `k`-th basis vectors.
                 
                 The given vector is a linear combination of basis vectors `b`, `b + 1`, ..., where
                 the linear combination is provided in the variable `vec`. The result is returned as
                 a `arithmetic::Real` number.
                 
                 Uses the default `arithmetic::RealContext`.
            */
            {
                return computeProjectionLengthR(k, b, vec, arithmetic::getThreadRealContext());
            }
        };
        
        /**@{
           \name Callback function types.
         */
        typedef boost::function<void(const linalg::math_matrix<arithmetic::Integer> &)> CallbackFunction;
        /**< \brief A generic callback function.
             
             The current lattice will be passed as an argument and shall not be modified. The
             function can throw an exception of type `stop_reduction` to abort reduction.
             
             \see \ref desc-genconf-callback
        */
        typedef boost::function<void(const linalg::math_matrix<arithmetic::NInt<long int> > &)> CallbackFunction_LI;
        /**< \brief A generic callback function for `long int` lattices.
             
             The current lattice will be passed as an argument and shall not be modified. The
             function can throw an exception of type `stop_reduction` to abort reduction.
             
             \see \ref desc-genconf-callback
        */
        
        typedef boost::function<void(const linalg::math_matrix<arithmetic::Integer> &, unsigned,
                                     const arithmetic::Integer &)> MinCallbackFunction;
        /**< \brief A callback function for newly found shortest vectors.
             
             The basis is given as an argument together with an index to the shortest vector in that
             base (encourntered so far). The function can throw an exception of type
             `stop_reduction` to abort reduction.
             
             \see \ref desc-genconf-callback
        */
        typedef boost::function<void(const linalg::math_matrix<arithmetic::NInt<long int> > &, unsigned,
                                     const arithmetic::NInt<long int> &)> MinCallbackFunction_LI;
        /**< \brief A callback function for newly found shortest vectors for `long int` lattices.
             
             The basis is given as an argument together with an index to the shortest vector in that
             base (encourntered so far). The function can throw an exception of type
             `stop_reduction` to abort reduction.
             
             \see \ref desc-genconf-callback
        */
        
        typedef boost::function<void(const linalg::math_matrix<arithmetic::Integer> & basis,
                                     int p,
                                     const linalg::math_rowvector<arithmetic::Integer> & vec)> EnumCallbackFunction;
        /**< \brief A callback function for newly found shortest vectors during enumerations, or
                    more generally SVP solving.
             
             The basis is given as an argument together with an index and a vector with
             coefficients, such that the shortest vector is a linear combination of the basis
             vectors `p`, `p + 1`, ... with the coefficients given in `vec`.
             
             The function can throw an exception of type `stop_enumeration` to stop enumeration (and
             return the currently found vector), or of type `stop_reduction` to abort reduction.
             
             \see \ref desc-genconf-callback
        */
        typedef boost::function<void(const linalg::math_matrix<arithmetic::NInt<long int> > & basis,
                                     int p,
                                     const linalg::math_rowvector<arithmetic::NInt<long int> > & vec)> EnumCallbackFunction_LI;
        /**< \brief A callback function for newly found shortest vectors during enumerations, or
                    more generally SVP solving. This callback is for `long int` lattices.
             
             The basis is given as an argument together with an index and a vector with
             coefficients, such that the shortest vector is a linear combination of the basis
             vectors `p`, `p + 1`, ... with the coefficients given in `vec`.
             
             The function can throw an exception of type `stop_enumeration` to stop enumeration (and
             return the currently found vector), or of type `stop_reduction` to abort reduction.
             
             \see \ref desc-genconf-callback
        */
        
        typedef boost::function<bool(const linalg::math_matrix<arithmetic::Integer> &,
                                     arithmetic::RealContext &,
                                     arithmetic::Real & T)> AnnealCallbackFunction;
        /**< \brief A callback function for Simulated Annealing.
             
             It will be given the current basis and current temperature. It should modify the
             temperature `T`, and return `true` if another annealing round should be started or
             `false` if the algorithm should terminate.
             
             If the given temperature `T` was negative when this function was called, then this call
             was done before the first annealing round and after the first LLL reduction.
        */
        
        typedef boost::function<bool(const linalg::math_matrix<arithmetic::Integer> &,
                                     int k,
                                     arithmetic::RealContext &,
                                     arithmetic::RandomNumberGenerator &,
                                     arithmetic::Real & T,
                                     GramSchmidtInformer *)> LLL_AnnealFunction;
        /**< \brief A annealing function for LLL.
             
             It receives the current basis, the current index, a `arithmetic::Real` context, a
             random number generator, the current temparature `T`, as well as a Gram-Schmidt
             informer object. It can return `true` to indicate that basis vector `k` should be
             swapped with basis vector `k - 1`, or `false` to indicate that nothing should be done.
        */
        
        typedef boost::function<bool (const linalg::math_matrix<arithmetic::Integer> &,
                                      int k, int windowsize,
                                      linalg::math_rowvector<arithmetic::Integer> & lincomb,
                                      arithmetic::RealContext &,
                                      arithmetic::RandomNumberGenerator &,
                                      arithmetic::Real & T,
                                      GramSchmidtInformer *)> BKZ_AnnealFunction;
        /**< \brief A annealing function for BKZ.
             
             It receives the current basis, the current index, the current window size, a linear
             combination `comb` of the basis vectors `k`, `k+1`, ... of the currently shortest
             vector found (by SVP solving), a `arithmetic::Real` context, a random number generator,
             the current temparature `T`, as well as a Gram-Schmidt informer object.
             
             It can return `true` to indicate that the vector given by the linear combination
             `lincomb` should be inserted, or `false` to not insert a vector. Note that the function
             can freely modify the linear combination.
        */
        
        typedef boost::function<void(VerboseLevel, const std::string &)> VerboseFunction;
        /**< \brief A verbose output callback function.
             
             Given a verbose level (of type `VerboseLevel`) and a message, should somehow inform the
             caller about the content of the message.
             
             \see \ref desc-genconf-verbose
        */
        ///@}
        
        LatticeReduction();
        /**< \brief Creates a default `LatticeReduction` object with default settings. */
        LatticeReduction(const linalg::math_matrix<arithmetic::Integer> & lattice);
        /**< \brief Creates a `LatticeReduction` object with default settings, and sets the lattice
                    basis to the given matrix.
             
             \param lattice The matrix whose rows to take as a generating system for the lattice.
             \sa void setLattice(const linalg::math_matrix<arithmetic::Integer> & lattice)
        */
        template<typename IType>
        LatticeReduction(const linalg::math_matrix<arithmetic::NInt<IType> > & lattice);
        /**< \brief Creates a `LatticeReduction` object with default settings, and sets the lattice
                    basis to the given matrix.
             
             \param lattice The matrix whose rows to take as a generating system for the lattice.
             \tparam IType A native CPU integer type. Must be `int`, `long int` or `long long`.
             \sa void setLattice(const linalg::math_matrix<arithmetic::Integer> & lattice)
        */
#if __cplusplus >= 201103L
        LatticeReduction(linalg::math_matrix<arithmetic::Integer> && lattice);
        /**< \brief Creates a `LatticeReduction` object with default settings, and moves the lattice
                    basis from the given matrix.
             
             \param lattice The matrix whose rows to take as a generating system for the lattice.
             \sa void setLattice(linalg::math_matrix<arithmetic::Integer> && lattice)
        */
        LatticeReduction(LatticeReduction &&);
        /**< \brief Moves a given `LatticeReduction` object into a newly constructed one. The old
                    `LatticeReduction` object is left in a state in which only deconstruction and
                    `operator = (LatticeReduction &&)` are valid calls.
        */
        LatticeReduction & operator = (LatticeReduction &&);
        /**< \brief Moves a given `LatticeReduction` object into this one. The old
                    `LatticeReduction` object is left in a state in which only deconstruction and
                    `operator = (LatticeReduction &&)` are valid calls.
        */
#endif
        ~LatticeReduction();
        /**< \brief Destroys the `LatticeReduction` object. */
        
        /**@{
           \name Lattice setting and retrieving functions.
         */
        void setLattice(const linalg::math_matrix<arithmetic::Integer> & lattice);
        /**< \brief Sets the current lattice of the `LatticeReduction` object to the lattice given
                    by the current matrix.
             
             The rows of the matrix `lattice` are taken as the generating set.
        */
        template<typename IType>
        void setLattice(const linalg::math_matrix<arithmetic::NInt<IType> > & lattice);
        /**< \brief Sets the current lattice of the `LatticeReduction` object to the lattice given
                    by the current matrix.
             
             The rows of the matrix `lattice` are taken as the generating set.
             
             \tparam IType A native CPU integer type. Must be `int`, `long int` or `long long`.
        */
#if __cplusplus >= 201103L
        void setLattice(linalg::math_matrix<arithmetic::Integer> && lattice);
        /**< \brief Sets the current lattice of the `LatticeReduction` object to the lattice given
                    by the current matrix. Leaves the matrix in a "moved-away" state.
             
             The rows of the matrix `lattice` are taken as the generating set.
        */
#endif
        const linalg::math_matrix<arithmetic::Integer> & getLattice() const;
        /**< \brief Returns the current lattice generating set as a matrix.
             
             The rows of the matrix are the vectors.
        */
        ///@}
        
        /**@{
           \name Lattice information functions.
         */
        unsigned rank() const;
        /**< \brief Returns the current rank, i.e. the number of generating vectors.
             
             After calling LLL or another reduction algorithm (not `sort` etc.), this equals the
             rank of the lattice since the vectors are in fact a basis.
        */
        unsigned dimension() const;
        /**< \brief Returns the dimension of the ambient space in which the lattice exists.
             
             When there is no linear dependence among the vectors, `dimension()` is never less than
             `rank()`.
        */
        ///@}
        
        /**@{
           \name Arithmetic configuration functions.
         */
        void setArithmetic(Arithmetic);
        /**< \brief Sets which arithmetic will be used for Gram-Schmidt orthogonalizations.
             
             \see \ref desc-genconf-arith
         */
        Arithmetic getArithmetic() const;
        /**< \brief Returns the currently set arithmetic for Gram-Schmidt orthogonalizations.
             
             \see \ref desc-genconf-arith
         */
        bool ensurePrecision(unsigned long);
        /**< \brief Ensures a minimal floating point precision (if possible).
             
             For variable precision floating point arithmetic (i.e. `A_Real`), ensures that the
             given minimum precision is used. Returns `true` if the precision can be set, or `false`
             if it cannot.
             
             \see \ref desc-genconf-arith
        */
        void setIntegers(Integers);
        /**< \brief Sets which integer arithmetic will be used for all lattice operations.
             
             \see \ref desc-genconf-arith
         */
        Integers getIntegers() const;
        /**< \brief Returns the currently set integer arithmetic.
             
             \see \ref desc-genconf-arith
         */
        ///@}
        
        /**@{
           \name Gram-Schmidt orthogonalization configuration functions.
         */
        void setGramSchmidt(GramSchmidt);
        /**< \brief Sets which Gram-Schmidt orthogonalization method will be used.
             
             \see \ref desc-genconf-gs
         */
        void setGramSchmidtRestart(bool);
        /**< \brief Sets that the current Gram-Schmidt method should use restarts.

             This makes only sense for `G_Classic` and `G_Givens`, where `swap()`, `add()`
             etc. modify the Gram-Schmidt coefficients directly and thus might propagate
             approximation errors.
             
             Currently experimental.
             
             \see \ref desc-genconf-gs
        */
        GramSchmidt getGramSchmidt() const;
        /**< \brief Returns the currently used Gram-Schmidt orthogonalization method.
             
             \see \ref desc-genconf-gs
         */
        bool getGramSchmidtRestart() const;
        /**< \brief Returns whether restarts are used for the current Gram-Schmidt orthogonalization
                    method.
             
             \see \ref desc-genconf-gs
        */
        ///@}
        
        /**@{
           \name Gram-Schmidt orthogonalization querying functions.
         */
        void forceGSRebuild(bool makeSureAllComputed = false);
        /**< \brief Forces the Gram-Schmidt coefficients to be rebuild.
             
             If the argument `makeSureAllComputed` is set to `true`, all Gram-Schmidt coefficients
             are build now and not when they are needed.
             
             This can be called when no algorithm is currently running.
        */
        double getGSCoefficientD(unsigned, unsigned) const;
        /**< \brief Returns the (i, j) Gram-Schmidt coefficient as a `double`.
             
             This can be called when no algorithm is currently running.
        */
        long double getGSCoefficientLD(unsigned, unsigned) const;
        /**< \brief Returns the (i, j) Gram-Schmidt coefficient as a `long double`.
             
             This can be called when no algorithm is currently running.
        */
        arithmetic::Real getGSCoefficientR(unsigned, unsigned) const;
        /**< \brief Returns the (i, j) Gram-Schmidt coefficient as a `arithmetic::Real` number.
             
             This can be called when no algorithm is currently running.
        */
        arithmetic::Real getGSCoefficientR(unsigned, unsigned, const arithmetic::RealContext &) const;
        /**< \brief Returns the (i, j) Gram-Schmidt coefficient as a `arithmetic::Real` number.
             
             This can be called when no algorithm is currently running. Uses the default
             `arithmetic::RealContext`.
        */
        double getGSSqNormD(unsigned) const;
        /**< \brief Returns the squared norm of the i-th Gram-Schmidt orthogonalized basis vector as
                    a `double`.
             
             This can be called when no algorithm is currently running.
        */
        long double getGSSqNormLD(unsigned) const;
        /**< \brief Returns the squared norm of the i-th Gram-Schmidt orthogonalized basis vector as
                    a `long double`.
             
             This can be called when no algorithm is currently running.
        */
        arithmetic::Real getGSSqNormR(unsigned) const;
        /**< \brief Returns the squared norm of the i-th Gram-Schmidt orthogonalized basis vector as
                    a `arithmetic::Real` number.
             
             This can be called when no algorithm is currently running.
        */
        arithmetic::Real getGSSqNormR(unsigned, const arithmetic::RealContext &) const;
        /**< \brief Returns the squared norm of the i-th Gram-Schmidt orthogonalized basis vector as
                    a `arithmetic::Real` number.
             
             This can be called when no algorithm is currently running. Uses the default
             `arithmetic::RealContext`.
        */
        ///@}
        
        /**@{
           \name Modifying functions.
         */
        void modFlip(unsigned);
        /**< \brief Flips the signs of the given basis vector.
             
             This can be called when no algorithm is currently running.
        */
        void modSwap(unsigned, unsigned);
        /**< \brief Swaps two given basis vectors.
             
             This can be called when no algorithm is currently running.
        */
        void modAdd(unsigned i, unsigned j, const arithmetic::Integer & m);
        /**< \brief Adds `m` times the `i`-th basis vector to the `j`-th basis vector.
             
             This can be called when no algorithm is currently running.
        */
        ///@}
        
        /**@{
           \name Transformation recording related functions.
         */
        void enableTransform(Transform = T_Normal);
        /**< \brief Enables that from this point on, a transformation matrix is recorded.
             
             The argument allows to chose which transformation matrix is created; by default, a
             normal one is created, such that if \f$A\f$ denotes the basis at this point, \f$A'\f$
             at a later point and \f$T\f$ the transformation matrix at that point, then \f$T \cdot A
             = A'\f$.
             
             \see \ref desc-genconf-trans
        */
        void disableTransform();
        /**< \brief Disables recording of a transformation matrix.
             
             \see \ref desc-genconf-trans
         */
        bool isTransformationRecorded() const;
        /**< \brief Queries whether a transformation matrix is recorded.
             
             \see \ref desc-genconf-trans
         */
        Transform getTransformationMode() const;
        /**< \brief Queries whether transformation matrix (normal or inverse) is recorded.
             
             \see \ref desc-genconf-trans
         */
        const linalg::math_matrix<arithmetic::Integer> * getTransformation() const;
        /**< \brief Queries the current transformation matrix.
             
             If none is recorded, `NULL` is returned.
             
             \see \ref desc-genconf-trans
        */
        ///@}
        
        /**@{
           \name SVP mode related functions.
         */
        void setSVPMode(SVPMode);
        /**< \brief Sets the current SVP solver. */
        SVPMode getSVPMode() const;
        /**< \brief Retrieves the current SVP solver. */
        ///@}
        
        /**@{
           \name Multi-threading related functions.
         */
        void setMaximalCoreUsage(unsigned);
        /**< \brief Setups the maximal number of cores which will be used by the `plll` library.
             
             If 0 is specified, as many as possible (i.e. needed and available) are used. This
             currently only affects parallel enumeration.
             
             The default value is 0.
             
             \see \ref desc-genconf-multithreading
        */
        unsigned getMaximalCoreUsage();
        /**< \brief Returns the current maximal number of cores used.
             
             \see \ref desc-genconf-multithreading
         */
        ///@}
        
        /**@{
           \name Callback related functions.
         */
        void setCallbackFunction(const CallbackFunction &, const CallbackFunction_LI & = CallbackFunction_LI());
        /**< \brief Sets a callback function.
             
             The callback function will be called in regular intervals during calls to reduction
             algorithms. These intervals are garuanteed minimal waiting times between calls; there
             is no upper bound; if for example an enumeration is running, this function will only be
             called when the enumeration eventually finishes.
             
             \warning If only a usual `CallbackFunction` callback function is set and the integer
                      arithmetic is set to `I_LongInt` or `I_Auto` (and `long int` lattices are
                      processed), the lattice must be converted for every call.
             
             \see \ref desc-genconf-callback
        */
        void setCallbackFunction(const CallbackFunction_LI &);
        /**< \brief Sets a callback function.
             
             The callback function will be called in regular intervals during calls to reduction
             algorithms. These intervals are garuanteed minimal waiting times between calls; there
             is no upper bound; if for example an enumeration is running, this function will only be
             called when the enumeration eventually finishes.
             
             \warning If the integer arithmetic is set to `I_ArbitraryPrecision` or `I_Auto` (and
                      `arithmetic::Integer` lattices are processed), the lattice must be converted
                      for every call.
             
             \see \ref desc-genconf-callback
        */
        void setCallbackInterval(double = 60.0*5.0);
        /**< \brief Sets the callback interval for the callback function.
             
             By default, the interval is every 5 minutes. The value can be set in seconds, where
             fractional multiples are allowed.
             
             \see \ref desc-genconf-callback
        */
        std::pair<CallbackFunction, CallbackFunction_LI> getCallbackFunction() const;
        /**< \brief Retrieves the currently set callback function.
             
             \see \ref desc-genconf-callback
         */
        double getCallbackInterval() const;
        /**< \brief Retrieves the current (minimal) callback interval.
             
             \see \ref desc-genconf-callback
         */
        
        void setMinCallbackFunction(const MinCallbackFunction &, const MinCallbackFunction_LI & = MinCallbackFunction_LI());
        /**< \brief Sets a minimum callback function which will be called as soon as a new shortest
                    vector is found during reduction.
             
             It will *not* be called during enumerations (see [setEnumCallbackFunction](@ref
             setEnumCallbackFunction) for this case), but only afterwards if the returned vector is
             shorter than the previous ones.
             
             (Note that the unprojected norms will be used to compare lengths.)
             
             \warning If only a usual `MinCallbackFunction` callback function is set and the integer
                      arithmetic is set to `I_LongInt` or `I_Auto` (and `long int` lattices are
                      processed), the lattice must be converted for every call.
             
             \see \ref desc-genconf-callback
        */
        void setMinCallbackFunction(const MinCallbackFunction_LI &);
        /**< \brief Sets a minimum callback function which will be called as soon as a new shortest
                    vector is found during reduction.
             
             It will *not* be called during enumerations (see [setEnumCallbackFunction](@ref
             setEnumCallbackFunction) for this case), but only afterwards if the returned vector is
             shorter than the previous ones.
             
             (Note that the unprojected norms will be used to compare lengths.)
             
             \warning If the integer arithmetic is set to `I_ArbitraryPrecision` or `I_Auto` (and
                      `arithmetic::Integer` lattices are processed), the lattice must be converted
                      for every call.
             
             \see \ref desc-genconf-callback
        */
        std::pair<MinCallbackFunction, MinCallbackFunction_LI> getMinCallbackFunction() const;
        /**< \brief Retrieves the currently set minimum callback function.
             
             \see \ref desc-genconf-callback
         */
        
        void setEnumCallbackFunction(const EnumCallbackFunction &, const EnumCallbackFunction_LI & = EnumCallbackFunction_LI());
        /**< \brief Sets a enumeration callback function which will be used if a new shortest vector
                    is found during an enumeration.
             
             Note that the found vector might be much longer than vectors found during previous
             enumerations or during previous calls of a `MinCallbackFunction`; it will only be
             shorter than the vectors found during the current enumeration.
             
             \warning If only a usual `EnumCallbackFunction` callback function is set and the
                      integer arithmetic is set to `I_LongInt` or `I_Auto` (and `long int` lattices
                      are processed), the lattice must be converted for every call.
             
             \see \ref desc-genconf-callback
        */
        void setEnumCallbackFunction(const EnumCallbackFunction_LI &);
        /**< \brief Sets a enumeration callback function which will be used if a new shortest vector
                    is found during an enumeration.
             
             Note that the found vector might be much longer than vectors found during previous
             enumerations or during previous calls of a `MinCallbackFunction`; it will only be
             shorter than the vectors found during the current enumeration.
             
             \warning If the integer arithmetic is set to `I_ArbitraryPrecision` or `I_Auto` (and
                      `arithmetic::Integer` lattices are processed), the lattice must be converted
                      for every call.
             
             \see \ref desc-genconf-callback
        */
        std::pair<EnumCallbackFunction, EnumCallbackFunction_LI> getEnumCallbackFunction();
        /**< \brief Retrieves the current enumeration callback function.
             
             \see \ref desc-genconf-callback
         */
        ///@}
        
        /**@{
           \name Simulated Annealing related functions.
         */
        void setDefaultAnnealing();
        /**< \brief Sets the default annealing functions.
             
             This is equivalent to a call both to `setDefaultAnnealingLLL()` and
             `setDefaultAnnealingBKZ()`. */
        void setDefaultAnnealingLLL();
        /**< \brief Sets the defalut LLL annealing function.

             \sa setDefaultAnnealing()
        */
        void setDefaultAnnealingBKZ();
        /**< \brief Sets the defalut BKZ annealing function.

             \sa setDefaultAnnealing()
        */
        void setAnnealing(const AnnealCallbackFunction &, const LLL_AnnealFunction &, const BKZ_AnnealFunction &);
        /**< \brief Sets annealing functions for both LLL and BKZ.
        */
        void setAnnealingLLL(const AnnealCallbackFunction &, const LLL_AnnealFunction &);
        /**< \brief Sets annealing functions for LLL.
             
             \sa setAnnealing()
        */
        void setAnnealingBKZ(const AnnealCallbackFunction &, const BKZ_AnnealFunction &);
        /**< \brief Sets annealing functions for BKZ.
             
             \sa setAnnealing()
        */
        void disableAnnealing();
        /**< \brief Disable annealing for both LLL and BKZ.
             
        */
        void disableAnnealingLLL();
        /**< \brief Disable annealing for LLL.
             
             \sa disableAnnealing()
        */
        void disableAnnealingBKZ();
        /**< \brief Disable annealing for LLL.
             
             \sa disableAnnealing()
        */
        bool isAnnealingLLLEnabled() const;
        /**< \brief Tests whether annealing for LLL is enabled. */
        bool isAnnealingBKZEnabled() const;
        /**< \brief Tests whether annealing for BKZ is enabled. */
        ///@}
        
        /**@{
           \name Deep Insertions related functions.
         */
        void setDeepInsertionMethod(DIMethod, DIMode = DIM_Default);
        /**< \brief Sets the Deep Insertions method and optionally also the mode.
             
             The default for the mode is `DIM_Default`.
             
             \see \ref desc-algconf-di
        */
        void setDeepInsertionMode(DIMode);
        /**< \brief Sets the Deep Insertions mode.
             
             \see \ref desc-algconf-di
         */
        void setDeepInsertionChoice(DIChoice, unsigned = 1);
        /**< \brief Sets the Deep Insertion choice and block size, which by default is 1.
             
             \see \ref desc-algconf-di
        */
        DIMethod getDeepInsertionMethod() const;
        /**< \brief Retrieves the current Deep Insertions method.
             
             \see \ref desc-algconf-di
         */
        DIMode getDeepInsertionMode() const;
        /**< \brief Retrieves the current Deep Insertions mode.
             
             \see \ref desc-algconf-di
         */
        DIChoice getDeepInsertionChoice() const;
        /**< \brief Retrieves the current Deep Insertions choice.
             
             \see \ref desc-algconf-di
         */
        unsigned getDeepInsertionBlocksize() const;
        /**< \brief Retrieves the current Deep Insertions choice block size parameter.
             
             \see \ref desc-algconf-di
         */
        ///@}
        
        /**@{
           \name Range related functions.
         */
        void setRange(unsigned begin, unsigned end = std::numeric_limits<unsigned>::max());
        /**< \brief Retricts all algorithms to the range [begin, end].
             
             If not said otherwise, they will be dealt with after projection into the orthogonal
             complement of the first `begin` basis vectors.
             
             If `end` is larger than the index of the last vector, all basis vectors following index
             `begin` will be used. The default is to handle *all* basis vectors.
             
             \see \ref desc-genconf-range
        */
        std::pair<unsigned, unsigned> getRange() const;
        /**< \brief Retrieves the current range.
             
             \see \ref desc-genconf-range
         */
        ///@}
        
        /**@{
           \name Lattice functions.
         */
        void sort(bool projected = false);
        /**< \brief Sorts the vectors by their norm.
             
             If the parameter `projected` is set to `true`, the Gram-Schmidt projected norms are
             used instead.
        */
        
        void sizereduction();
        /**< \brief Applies only size reduction to the vectors.
             
             \see \ref howitworks
         */
        void lll(double alpha = 0.99, LLLMode mode = LLL_Classic);
        /**< \brief Applies LLL with given reduction parameter `alpha` (\f$\alpha\f$) and mode
                    `mode`.
             
             By default, the reduction parameter is 0.99 and the mode is `LLL_Classic`.
             
             Note that often 0.99 completely suffices as a reduction parameter, and that the default
             `alpha` used by Lenstra, Lenstra and Lovász is 0.75.
             
             \see \ref desc-algs-lll
        */
        void bkz(double alpha = 0.99, unsigned blocksize = 20, BKZMode mode = BKZ_SchnorrEuchner);
        /**< \brief Applies BKZ (or one of its variants) with given reduction parameter `alpha`
                    (\f$\alpha\f$), block size `blocksize` and mode `mode`.
             
             By default, the reduction parameter is 0.99, the block size is 20 and the mode is
             `BKZ_SchnorrEuchner`.
             
             Note that often 0.99 completely suffices as a reduction parameter, and that the default
             `alpha` used by Lenstra, Lenstra and Lovász is 0.75.
             
             \see \ref desc-algs-bkz
        */
        void hkz(bool dual = false);
        /**< \brief Computes a Hermite-Korkine-Zolotarev reduced basis.
             
             This is one of the strongest reduction method known and cannot be computed
             efficiently. It will be computed by \f$k\f$ calls to the SVP solver if a lattice of
             rank \f$k\f$ is given.
             
             In case `dual` is set to `true`, a HKZ basis will be computed of the reversed dual of
             the current lattice.
        */
        void svp(bool make_basis = true, bool extreme = false, bool dual = false);
        /**< \brief Computes a shortest vector of the lattice and inserts it at the beginning.
             
             If `make_basis` is set to `true`, the resulting generating system is reduced to a basis
             using LLL, so that the first vector is the shortest (such bases are also called *SVP
             bases*).
             
             If `extreme` is set to `true`, *Extreme Pruning* as described by N. Gama, P. Q. Nguyen
             and O. Regev in
             \cite gama-nguyen-regev-extremepruning will be used. Note that this is experimental at
             the moment.
             
             In case `dual` is set to `true`, a SVP generating system respectively basis will be
             computed of the reversed dual of the current lattice.
             
             \see \ref desc-algs-svp
        */
        ///@}
        
        /**@{
           \name Lattice reduction testing functions.
         */
        bool isSizeReduced() const;
        /**< \brief Tests whether the given lattice is size reduced. */
        bool isLLLBasis(double alpha = 0.99, LLLMode mode = LLL_Classic) const;
        /**< \brief Tests whether the given lattice basis is LLL reduced with respect to the given
                    reduction parameter `alpha` (\f$\alpha\f$) and mode `mode`.
             
             \sa lll(alpha, mode)
        */
        bool isBKZBasis(double alpha = 0.99, unsigned blocksize = 20, BKZMode mode = BKZ_SchnorrEuchner) const;
        /**< \brief Tests whether the given lattice basis is BKZ reduced with respect to the given
                    reduction parameter `alpha` (\f$\alpha\f$), the given block size `blocksize` and
                    mode `mode`.
             
             \sa bkz(alpha, blocksize, mode)
        */
        bool isHKZBasis(bool dual = false) const;
        /**< \brief Tests whether the given lattice basis (or its reversed dual, in case `dual` is
                    `true`) is HKZ reduced.
             
             \sa hkz(dual)
        */
        bool isSVPBasis(bool dual = false) const;
        /**< \brief Tests whether the given lattice basis (or its reversed dual, in case `dual` is
                    `true`) has a shortest vector at the first index.
             
             \sa svp(*, *, dual)
        */
        ///@}
        
        /**@{
           \name Statistics related functions.
         */
        const Statistics & getStatistics() const;
        /**< \brief Retrieves the current statistics object. */
        void resetStatistics();
        /**< \brief Resets the current statistics object. */
        ///@}
        
        /**@{
           \name Verbosity related functions.
         */
        void setVerbose(VerboseOutputLevel level, const VerboseFunction & = 0);
        /**< \brief Sets the verbose output level to `level`.
             
             If given, the current verbose output function will also be set.
             
             \see \ref desc-genconf-verbose
        */
        VerboseOutputLevel getVerboseOutputLevel();
        /**< \brief Retrieves the current verbose output level.
             
             \see \ref desc-genconf-verbose
         */
        const VerboseFunction & getVerboseFunction();
        /**< \brief Retrieves the current verbose output function.
             
             \see \ref desc-genconf-verbose
         */
        ///@}
    };
}

#endif
