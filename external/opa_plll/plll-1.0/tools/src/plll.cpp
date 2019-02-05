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

#include <plll.hpp>
#include <plll/linalg.hpp>
#include "profiling.hpp"
#include <plll/arguments.hpp>
#include <plll/linalg.hpp>
#include <list>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cctype>
#include <cmath>
#include <math.h> // for C98 functionality; mainly the tgamma() function
#include <unistd.h> // for isatty() (to detect redirection)

using namespace plll;

arithmetic::Integer * determinant = NULL;
bool determinant_is_square;

inline void argsEnum(const std::string & a, const helper::ArgumentParser::Value &)
{
    std::cerr << " " << a;
}

bool loadLattice(LatticeReduction & lr, const std::string & fn)
{
    linalg::math_matrix<arithmetic::Integer> lattice;
    std::ifstream f(fn.c_str());
    f >> lattice;
    if (!f)
        return false;
    while (f)
    {
        // Read additional stuff, like determinant
        std::string line, cmd;
        getline(f, line);
        std::istringstream ss(line);
        ss >> cmd;
        if (ss)
        {
            if (cmd == "DET")
            {
                determinant = new arithmetic::Integer();
                determinant_is_square = false;
                ss >> cmd;
                if (ss && (cmd == "SQRT"))
                {
                    determinant_is_square = true;
                    ss >> cmd;
                }
                if (ss)
                {
                    *determinant = arithmetic::convert<arithmetic::Integer>(cmd.c_str());
                    if (isZero(*determinant))
                        return false;
                }
                else
                    return false;
            }
            else
            {
                std::cerr << "Unknown directive '" << cmd << "'!\n";
                return false;
            }
        }
    }
//    std::cout << lattice << "\n";
    lr.setLattice(lattice);
    return true;
}

template<typename MatrixType>
void writeLatticeNTL(std::ostream & f, const MatrixType & lattice)
{
    f << "[";
    for (unsigned i = 0; i < lattice.rows(); ++i)
    {
        f << "[";
        for (unsigned j = 0; j < lattice.cols(); ++j)
        {
            if (j)
                f << " ";
            f << lattice(i, j);
        }
        f << "]\n";
    }
    f << "]\n";
}

template<typename MatrixType>
bool storeLattice(const MatrixType & lattice, const std::string & fn)
{
    std::ofstream f(fn.c_str());
    writeLatticeNTL(f, lattice);
    if (determinant)
    {
        f << "DET ";
        if (determinant_is_square)
            f << "SQRT ";
        f << *determinant << "\n";
    }
    return true;
}

bool storeLattice(const LatticeReduction & lr, const std::string & fn)
{
    return storeLattice(lr.getLattice(), fn);
}

linalg::math_rowvector<arithmetic::Integer> parseVector(const std::string & str)
{
    linalg::math_rowvector<arithmetic::Integer> output;
    std::istringstream s(str);
    s >> output;
    return output;
}

template<typename IntType>
void outputShortest(const linalg::math_matrix<IntType> & lattice)
{
    int idx = 0, idx2 = 0;
    IntType norm = normSq(lattice.row(0));
    IntType norm2 = normSq(lattice.row(0));
    for (unsigned i = 1; i < lattice.rows(); ++i)
    {
        IntType n = normSq(lattice.row(i));
        if (isZero(n))
            continue;
        if ((n < norm) || isZero(norm))
        {
            idx = i;
            norm = n;
        }
        if (n > norm2)
        {
            idx2 = i;
            norm2 = n;
        }
    }
    std::cout << "=============================================================================\n";
    std::cout << "Shortest vector is #" << idx << " with squared norm " << norm << "\n";
    std::cout << "Longest vector is #" << idx2 << " with squared norm " << norm2 << "\n";
    arithmetic::RealContext rc;
    std::cout << "Norm: " << sqrt(arithmetic::convert(norm, rc)) << ", " << sqrt(arithmetic::convert(norm2, rc)) << "\n";
    std::cout << lattice.row(idx) << "\n";
}

inline void outputShortest(LatticeReduction & lr)
{
    outputShortest<arithmetic::Integer>(lr.getLattice());
}

std::vector<std::string> LLLCFfiles;
unsigned LLLCFindex;

void InitLLLCallbackFunction()
{
    arithmetic::RandomNumberGenerator rng;
    rng.randomizeTime();
    std::string basefn = "svpchallenge-temp-";
    std::string id;
    arithmetic::convert(id, rng.randomBits(64), arithmetic::HexStringContext());
    int add = (64 / 4) - id.size();
    while (add--) basefn += '0';
    basefn += id;
    LLLCFfiles.resize(2);
    LLLCFfiles[0] = basefn + "-1.lattice";
    LLLCFfiles[1] = basefn + "-2.lattice";
    std::cout << "\nWill be using temporary files:";
    for (std::vector<std::string>::iterator i = LLLCFfiles.begin(); i != LLLCFfiles.end(); ++i)
        std::cout << " " << *i;
    std::cout << "\n\n";
    LLLCFindex = 0;
}

std::string timestamp()
{
    time_t ltime = time(NULL);
    struct tm * TM = localtime(&ltime);;
    std::ostringstream ts;
    ts << std::setfill('0');
    ts << std::setw(4) << TM->tm_year + 1900 << "/"
       << std::setw(2) << TM->tm_mon + 1 << "/"
       << std::setw(2) << TM->tm_mday << " "
       << std::setw(2) << TM->tm_hour << ":"
       << std::setw(2) << TM->tm_min << ":"
       << std::setw(2) << TM->tm_sec;
    return ts.str();
}

bool minLogFile = false;
bool dumpLatticeEveryStep = false;
std::fstream minLog;

bool LLLCallbackFunction(const linalg::math_matrix<arithmetic::Integer> & A)
{
    std::cout << timestamp() << " Writing current intermediate result to " << LLLCFfiles[LLLCFindex] << "...\n";
    if (!storeLattice(A, LLLCFfiles[LLLCFindex]))
        std::cerr << "Error while writing " << LLLCFfiles[LLLCFindex] << "!\n";
    outputShortest<arithmetic::Integer>(A);
    if (dumpLatticeEveryStep)
        std::cout << "Basis:\n" << A << "\n";
    LLLCFindex = (LLLCFindex + 1) % LLLCFfiles.size();
    if (minLogFile)
        minLog.flush();
    return true;
}

bool LLLCallbackFunction_LI(const linalg::math_matrix<arithmetic::NInt<long int> > & A)
{
    std::cout << timestamp() << " Writing current intermediate result to " << LLLCFfiles[LLLCFindex] << "...\n";
    if (!storeLattice(A, LLLCFfiles[LLLCFindex]))
        std::cerr << "Error while writing " << LLLCFfiles[LLLCFindex] << "!\n";
    outputShortest<arithmetic::NInt<long int> >(A);
    if (dumpLatticeEveryStep)
        std::cout << "Basis:\n" << A << "\n";
    LLLCFindex = (LLLCFindex + 1) % LLLCFfiles.size();
    if (minLogFile)
        minLog.flush();
    return true;
}

arithmetic::Real * shortestEstimate = NULL;
bool firstmin = true;
arithmetic::Integer * currentmin = NULL;
unsigned latticeRank;
Timer estimateTiming;
CPUTimer estimateTimingCPU;
LatticeReduction * latticered = NULL;
bool mincallback_dumpvectors = false;

void InitMinCallbackFunction(LatticeReduction & lr, const std::string & logfn)
{
    latticered = &lr;
    currentmin = new arithmetic::Integer();
    if (determinant != NULL)
    {
        latticeRank = lr.rank(); // hoping that the vectors generating the lattice are linearly independent
        shortestEstimate = new arithmetic::Real(*determinant);
        if (determinant_is_square)
            *shortestEstimate = sqrt(*shortestEstimate);
        else
            *shortestEstimate = abs(*shortestEstimate);
        arithmetic::RealContext rc;
        *shortestEstimate = power(*shortestEstimate, arithmetic::convert(1.0, rc) / arithmetic::convert(latticeRank, rc));
        if (logfn.size() > 0)
        {
            minLogFile = true;
            minLog.open(logfn.c_str(), std::ios_base::out);
        }
        estimateTimingCPU.start();
        estimateTiming.start();
    }
    else
        shortestEstimate = NULL;
}

template<class Vec, class IntType>
void noteNewShortestVector(const Vec & vec, const IntType & len, bool duringEnum)
{
    *currentmin = arithmetic::convert<arithmetic::Integer>(len);
    std::ostream & s = minLogFile ? minLog : std::cout;
    s << "Current ";
    if (shortestEstimate)
    {
        arithmetic::RealContext rc;
        arithmetic::Real
            norm = sqrt(arithmetic::convert(len, rc)),
            approxfact = norm / *shortestEstimate,
            approxfact1n = power(approxfact, arithmetic::convert(1.0, rc) / arithmetic::convert(latticeRank, rc));
        s << std::setprecision(10)
            // This should be called approximation factor and not Hermite factor!
          << "Hermite factor: " << arithmetic::convert<long double>(approxfact)
          << " (n-th root is " << arithmetic::convert<long double>(approxfact1n) << "); ";
    }
    s << "squared norm: " << *currentmin << "; ";
//    s << "index: " << idx << "; ";
    s << "timestamp: "
      << "wallclock " << (long double)estimateTiming.elapsed() / (long double)Timer::d_one_timeunit << " " << Timer::d_timeunit_name << ", "
      << "CPU " << (long double)estimateTimingCPU.elapsed() / (long double)CPUTimer::d_one_timeunit << " " << CPUTimer::d_timeunit_name << "; ";
    const LatticeReduction::Statistics & stats = latticered->getStatistics();
    s << "[s:" << stats.swaps << " a:" << stats.adds << " a1:" << stats.adds_pm1 << " a2:" << stats.adds_pm2 << " f:" << stats.flips << " t:" << stats.trans
      << " sr:" << stats.sizereductions << " di:" << stats.deepinsertions << " e:" << stats.enumcalls << " ef:" << stats.enumfails << " vi:" << stats.vectorinsertions << " vir:" << stats.vectorinsertions_rearrange << "]";
    if (duringEnum)
        s << " (during enumeration)";
    if (mincallback_dumpvectors)
        s << " shortest vector: " << vec;
    s << "\n";
    firstmin = false;
}

void EnumCallbackFunction(const linalg::math_matrix<arithmetic::Integer> & A, int idx, const linalg::math_rowvector<arithmetic::Integer> & vec)
{
    arithmetic::Integer len;
    if (idx >= 0)
    {
        linalg::math_rowvector<arithmetic::Integer> v(A.cols());
        for (unsigned i = 0; i < vec.size(); ++i)
            if (!isZero(vec[i]))
                v += A.row(idx + i) * vec[i];
        normSq(len, v);
        if (firstmin || (len < *currentmin))
            noteNewShortestVector(v, len, true);
    }
    else
    {
        normSq(len, vec);
        if (firstmin || (len < *currentmin))
            noteNewShortestVector(vec, len, true);
    }
}

void EnumCallbackFunction_LI(const linalg::math_matrix<arithmetic::NInt<long int> > & A, int idx, const linalg::math_rowvector<arithmetic::NInt<long int> > & vec)
{
    arithmetic::NInt<long int> len;
    setZero(len);
    if (idx >= 0)
    {
        linalg::math_rowvector<arithmetic::NInt<long int> > v(A.cols());
        for (unsigned i = 0; i < vec.size(); ++i)
            if (!isZero(vec[i]))
                v += A.row(idx + i) * vec[i];
        normSq(len, v);
        if (firstmin || (arithmetic::convert<arithmetic::Integer>(len) < *currentmin))
            noteNewShortestVector(v, len, true);
    }
    else
    {
        normSq(len, vec);
        if (firstmin || (arithmetic::convert<arithmetic::Integer>(len) < *currentmin))
            noteNewShortestVector(vec, len, true);
    }
}

void MinCallbackFunction(const linalg::math_matrix<arithmetic::Integer> & A, unsigned idx, const arithmetic::Integer & len)
{
    if (firstmin || (len < *currentmin))
        noteNewShortestVector(A.row(idx), len, false);
}

void MinCallbackFunction_LI(const linalg::math_matrix<arithmetic::NInt<long int> > & A, unsigned idx, const arithmetic::NInt<long int> & len)
{
    if (firstmin || (arithmetic::convert<arithmetic::Integer>(len) < *currentmin))
        noteNewShortestVector(A.row(idx), len, false);
}

void MinCallbackLogClose()
{
    if (minLogFile)
    {
        estimateTiming.stop();
        estimateTimingCPU.stop();
        minLog << std::setprecision(10)
               << "Done; timestamp: "
               << "wallclock " << (long double)estimateTiming.elapsed() / (long double)Timer::d_one_timeunit << " " << Timer::d_timeunit_name << ", "
               << "CPU " << (long double)estimateTimingCPU.elapsed() / (long double)CPUTimer::d_one_timeunit << " " << CPUTimer::d_timeunit_name << "; ";
        const LatticeReduction::Statistics & stats = latticered->getStatistics();
        minLog << "[s:" << stats.swaps << " a:" << stats.adds << " a1:" << stats.adds_pm1 << " a2:" << stats.adds_pm2 << " f:" << stats.flips << " t:" << stats.trans
               << " sr:" << stats.sizereductions << " di:" << stats.deepinsertions << " e:" << stats.enumcalls << " ef:" << stats.enumfails << " vi:" << stats.vectorinsertions << " vir:" << stats.vectorinsertions_rearrange << "]";
        minLog << "\n";
        minLog.close();
    }
    delete currentmin;
    currentmin = NULL;
    delete shortestEstimate;
    shortestEstimate = NULL;
}

void EnumCallbackFunctionSilent(const linalg::math_matrix<arithmetic::Integer> & A, int idx, const linalg::math_rowvector<arithmetic::Integer> & vec)
{
}

void EnumCallbackFunctionSilent_LI(const linalg::math_matrix<arithmetic::NInt<long int> > & A, int idx, const linalg::math_rowvector<arithmetic::NInt<long int> > & vec)
{
}

void MinCallbackFunctionSilent(const linalg::math_matrix<arithmetic::Integer> & A, unsigned idx, const arithmetic::Integer & len)
{
}

void MinCallbackFunctionSilent_LI(const linalg::math_matrix<arithmetic::NInt<long int> > & A, unsigned idx, const arithmetic::NInt<long int> & len)
{
}

bool use_colors = true;

const char * RED = "\033[31m";
const char * GREEN = "\033[32m";
const char * YELLOW = "\033[33m";
const char * NORMAL = "\033[0m";

void verify(bool result, const std::string & basis)
{
    std::cout << "VERIFY: ";
    if (result)
    {
        if (use_colors)
            std::cout << GREEN;
        std::cout << "Computed basis is " << basis << ".";
    }
    else
    {
        if (use_colors)
            std::cout << RED;
        std::cout << "Computed basis is not " << basis << "!";
    }
    if (use_colors)
        std::cout << NORMAL;
    std::cout << "\n";
}

void help(const char * cmd = "test-svpchallenge2")
{
    std::cout << cmd << " <lattice name> <options>\n";
    std::cout << "\n";
    std::cout << "Features:\n";
    std::cout << "  -arith=<mode>, -arithmetic (mode = rational/rat, real/r, longdouble/ld [default], double/d, doubledouble/dd, quaddouble/qd)\n";
    std::cout << "  -ints=<mode>, -integers (mode = automatic/auto [default], bigint/big/bi, longint/long/l)\n";
    std::cout << "  -gs=<mode>, -gramschmidt (mode = gramschmidt/gs, integer/int, gsrestart/restart, givens/G, givensrestart/Gr, numstable/numstab/ns [default])\n";
    std::cout << "  -enum=<mode>, -enumeration (mode = kse, parallelkse/pkse [default], schnorrfast, listsieve/ls, listsievebirthday/lsb, gausssieve/gs, voronoicell/vc)\n";
    std::cout << "  -di=<mode>, -deepinsertion (mode = all [default], first, block, both)\n";
    std::cout << "    -dialg=<mode> (mode = classic [default], minimizepotential1/minpot1/minpot, minimizepotential2/minpot2)\n";
    std::cout << "    -dimode=<mode> (mode = before, after [default], both)\n";
    std::cout << "    -dibs=<block size>, -diblocksize\n";
    std::cout << "  -anneal, -annealing\n";
    std::cout << "  -maxthreads=<maximal threads> (0 = as many as CPUs/cores)\n";
    std::cout << "  -begin=<int>, -end=<int>\n";
    std::cout << "  -mintrack, -mintracklog=<filename>, -mintrackincvecs, -mintracksilent\n";
    std::cout << "  -verbose: output more verbose information\n";
    std::cout << "  -verbose=full, -vv: output very detailed verbose information\n";
    std::cout << "Note that deep insertions and annealing only affect LLL and BKZ\n";
    std::cout << "\n";
    std::cout << "Commands w/ options:\n";
    std::cout << "  -info: outputs information on lattice and computes determinant if necessary\n";
    std::cout << "  -gsinfo: outputs Gram-Schmidt coefficients\n";
    std::cout << "  -sort=<mode>: sort vectors by norm (mode = unprojected/unproj [default], projected/proj)\n";
    std::cout << "  -sr: performs size reduction\n";
    std::cout << "  -ssred=<mode>: performs sorting+sizered reduction (mode = unprojected/unproj [default], projected/proj)\n";
    std::cout << "  -lll=<mode>: performs LLL reduction (mode = classic [default], unprojected/unproj/up, siegel/s)\n";
    std::cout << "    -alpha=<alpha> (reduction parameter, default: 0.99)\n";
    std::cout << "  -bkz=<mode>: performs BKZ reduction (mode = schnorreuchner [default], simplified, hanrotpujolstehle-hkz/hps-hkz, hanrotpujolstehle-svp/hps-svp, primaldual/pd, slide/s, improvedslide/is, improvedslide[23]/is[23], semiblock2k/sb, samplingreduction/sr, experimental/ex)\n";
    std::cout << "    -alpha=<alpha> (reduction parameter, default: 0.99)\n";
    std::cout << "    -window=<window size>, -windowsize, -bs, -blocksize\n";
    std::cout << "  -hkz: performs HKZ reduction\n";
    std::cout << "  -dhkz: performs HKZ reduction on the dual basis\n";
    std::cout << "  -svp=<mode>: computes a shortest vector (mode = normal [default], extreme/ex)\n";
    std::cout << "    -nosvpbasis: inserts shortest vector, but does not turn resulting generating system into basis\n";
    std::cout << "  -dsvp=<mode>: computes a shortest vector in the dual basis (mode = normal [default], extreme/ex)\n";
    std::cout << "    -nosvpbasis: inserts shortest vector, but does not turn resulting generating system into basis\n";
    std::cout << "  -test=<vector>: tests whether the vector is in the lattice; if yes, writes it in terms of the basis\n";
    std::cout << "  -hnf, -hnf=right: compute basis in Hermite Normal Form (diagonal on the left [default] or right)\n";
    std::cout << "  -pp, -preprocessing: applies a round of LLL without any Deep Insertion / Annealing settings before applying any other reduction; reduction factor can be given as value, default is 0.75\n";
    std::cout << "\n";
    std::cout << "Output:\n";
    std::cout << "  -write=<filename>\n";
    std::cout << "  -dump, -dump=ntl: dump lattice to standard output (in linalg::math_matrix or NTL format)\n";
    std::cout << "\n";
}

void computeDeterminant(const LatticeReduction & lr)
{
    if (determinant == NULL)
        determinant = new arithmetic::Integer();
    const linalg::math_matrix<arithmetic::Integer> & lattice = lr.getLattice();
    std::cout << "Computing determinant..." << std::flush;
    if (lattice.rows() == lattice.cols())
    {
        // Check for special form: lower/upper triangular matrix
        bool is_upper = true, is_lower = true;
        for (unsigned i = 0; i < lattice.rows(); ++i)
            for (unsigned j = 0; j < lattice.cols(); ++j)
                if (!isZero(lattice(i, j)))
                {
                    if (j < i)
                        is_lower = false;
                    if (j > i)
                        is_upper = false;
                }
        // Is it a triangular matrix?
        if (is_upper || is_lower)
        {
            // Yes: determinant is product of diagonal elements
            setOne(*determinant);
            for (unsigned i = 0; i < lattice.rows(); ++i)
                *determinant *= lattice(i, i);
        }
        else
            // No: use default approach to compute determinant
            *determinant = det(lattice);
        // This is the non-squared determinant
        determinant_is_square = false;
    }
    else
    {
        // Compute determinant of Gram matrix and take square root
        *determinant = det(lattice * lattice.transpose());
        determinant_is_square = true;
    }
    std::cout << "\n";
}

void MyVerboseFunction(LatticeReduction::VerboseLevel level, const std::string & s)
{
    if (use_colors)
    {
        switch (level)
        {
        case LatticeReduction::VL_Chatter: std::cerr << NORMAL; break;
        case LatticeReduction::VL_Information: std::cerr << GREEN; break;
        case LatticeReduction::VL_Warning: std::cerr << YELLOW; break;
        case LatticeReduction::VL_Error: std::cerr << RED; break;
        }
    }
    std::cerr << s;
    if (use_colors)
        std::cerr << NORMAL;
    std::cerr << "\n";
}

int main(int argc, char **argv)
{
    arithmetic::initArithmeticThreadAllocators();
    
    helper::ArgumentParser args(argc, argv);
    
    if ((args.getValue("help") != NULL) || (args.getValue("h") != NULL))
    {
        help(argv[0]);
        return -1;
    }
    
    if (!isatty(fileno(stdout)) || !isatty(fileno(stderr)))
        use_colors = false;
    
    std::string fn = "";
    
    // Go through names (i.e. anything not beginning with - or --)
    const std::list<std::string> & names = args.getNames();
    std::list<std::string>::const_iterator namesit = names.begin();
    if (namesit != names.end())
    {
        fn = *namesit;
        ++namesit;
    }
    if (namesit != names.end())
    {
        std::cerr << "Too many arguments!\n";
        return -1;
    }
    if (fn == "")
    {
        std::cerr << "No lattice filename specified!\n";
        return -1;
    }
    
    // Options
    const helper::ArgumentParser::Value * v, * vv;
    // Verbose level
    LatticeReduction::VerboseOutputLevel verb = LatticeReduction::VOL_Warnings; // VOL_Informative, VOL_Full
    if ((v = args.getValue("v")) == NULL)
        v = args.getValue("verbose");
    if (v != NULL)
    {
        verb = LatticeReduction::VOL_Informative;
        if (v->getText() == "full")
            verb = LatticeReduction::VOL_Full;
    }
    if ((v = args.getValue("vv")) != NULL)
        verb = LatticeReduction::VOL_Full;
    LatticeReduction lr;
    lr.setVerbose(verb, MyVerboseFunction);
    
    // Load lattice
    std::cout << "Loading '" << fn << "'...\n";
    if (!loadLattice(lr, fn))
    {
        std::cerr << "Could not load '" << fn << "'!\n";
        return -1;
    }
    outputShortest(lr);
    
    // Options
    std::ostringstream config;
    // RNG seed
    if ((v = args.getValue("seed")) != NULL)
    {
        double goodness = 0.01;
        if (v->getText() == "seed")
            ; // default value
        else if (v->isInteger() || (v->getInteger() == 1) || (v->getInteger() == 0))
            goodness = 1.0;
        else if (v->isFloat() && (v->getFloat() > 0) && (v->getFloat() < 1))
            goodness = v->getFloat();
        else
        {
            std::cerr << "seed must be real number in [0, 1]!\n";
            return -1;
        }
        arithmetic::RandomNumberGenerator::initializeSeeder(goodness);
    }
    // Arithmetic
    if ((v = args.getValue("arithmetic")) == NULL)
        v = args.getValue("arith");
    if (v != NULL)
    {
        if ((v->getText() == "rational") || (v->getText() == "rat"))
        {
            lr.setArithmetic(LatticeReduction::A_Rational);
            config << "rational arithmetic";
        }
        else if ((v->getText() == "real") || (v->getText() == "r"))
        {
            lr.setArithmetic(LatticeReduction::A_Real);
            config << "real arithmetic";
        }
        else if ((v->getText() == "longdouble") || (v->getText() == "ld"))
        {
            lr.setArithmetic(LatticeReduction::A_LongDouble);
            config << "long double arithmetic";
        }
        else if ((v->getText() == "double") || (v->getText() == "d"))
        {
            lr.setArithmetic(LatticeReduction::A_Double);
            config << "double arithmetic";
        }
        else if ((v->getText() == "doubledouble") || (v->getText() == "dd"))
        {
            lr.setArithmetic(LatticeReduction::A_DoubleDouble);
            config << "double-double arithmetic";
        }
        else if ((v->getText() == "quaddouble") || (v->getText() == "qd"))
        {
            lr.setArithmetic(LatticeReduction::A_QuadDouble);
            config << "quad-double arithmetic";
        }
        else
        {
            std::cerr << "Unknown arithmetic \"" << v->getText() << "\"\n";
            std::cerr << "Choices: rational/rat, real/r, longdouble/ld [default], double/d, doubledouble/dd, quaddouble/qd\n";
            return -1;
        }
    }
    else
    {
        lr.setArithmetic(LatticeReduction::A_LongDouble);
        config << "long double arithmetic";
    }
    // Integers
    if ((v = args.getValue("integers")) == NULL)
        v = args.getValue("ints");
    if (v != NULL)
    {
        if ((v->getText() == "automatic") || (v->getText() == "auto"))
        {
            lr.setIntegers(LatticeReduction::I_Auto);
            config << ", automatic integer selection";
        }
        else if ((v->getText() == "bigint") || (v->getText() == "big") || (v->getText() == "bi"))
        {
            lr.setIntegers(LatticeReduction::I_ArbitraryPrecision);
            config << ", arbitrary precision integers";
        }
        else if ((v->getText() == "longint") || (v->getText() == "long") || (v->getText() == "l"))
        {
            lr.setIntegers(LatticeReduction::I_LongInt);
            config << ", long integers";
        }
        else
        {
            std::cerr << "Unknown integers selection\"" << v->getText() << "\"\n";
            std::cerr << "Choices: automatic/auto [default], bigint/big/bi, longint/long/l\n";
            return -1;
        }
    }
    else
    {
        lr.setIntegers(LatticeReduction::I_Auto);
        config << ", automatic integer selection";
    }
    // Gram-Schmidt
    if ((v = args.getValue("gramschmidt")) == NULL)
        v = args.getValue("gs");
    if (v != NULL)
    {
        if ((v->getText() == "gramschmidt") || (v->getText() == "gs"))
        {
            lr.setGramSchmidt(LatticeReduction::G_Classic);
            config << ", classical Gram-Schmidt";
        }
        else if ((v->getText() == "integer") || (v->getText() == "int"))
        {
            lr.setGramSchmidt(LatticeReduction::G_ClassicInteger);
            config << ", integer Gram-Schmidt";
        }
        else if ((v->getText() == "gsrestart") || (v->getText() == "restart"))
        {
            lr.setGramSchmidt(LatticeReduction::G_Classic);
            lr.setGramSchmidtRestart(true);
            config << ", classical Gram-Schmidt w/ restarts";
        }
        else if ((v->getText() == "givens") || (v->getText() == "G"))
        {
            lr.setGramSchmidt(LatticeReduction::G_Givens);
            config << ", Givens rotations";
        }
        else if ((v->getText() == "givensrestart") || (v->getText() == "Gr"))
        {
            lr.setGramSchmidt(LatticeReduction::G_Givens);
            lr.setGramSchmidtRestart(true);
            config << ", Givens rotations w/ restarts";
        }
        else if ((v->getText() == "numstable") || (v->getText() == "numstab") || (v->getText() == "ns"))
        {
            lr.setGramSchmidt(LatticeReduction::G_NumStable);
            config << ", numerically stable GS";
        }
        else
        {
            std::cerr << "Unknown Gram-Schmidt mode \"" << v->getText() << "\"\n";
            std::cerr << "Choices: gramschmidt/gs, integer/int, gsrestart/restart, givens/G, givensrestart/Gr, numstable/numstab/ns [default]\n";
            return -1;
        }
    }
    else
    {
        lr.setGramSchmidt(LatticeReduction::G_NumStable);
        config << ", numerically stable GS";
    }
    // Enumeration
    if ((v = args.getValue("enumeration")) == NULL)
        v = args.getValue("enum");
    if (v != NULL)
    {
        if (v->getText() == "kse")
        {
            lr.setSVPMode(LatticeReduction::SVP_KannanSchnorrEuchner);
            config << ", Kannan-Schnorr-Euchner enumeration";
        }
        else if ((v->getText() == "parallelkse") || (v->getText() == "pkse"))
        {
            lr.setSVPMode(LatticeReduction::SVP_ParallelKannanSchnorrEuchner);
            config << ", parallel Kannan-Schnorr-Euchner enumeration";
        }
        else if (v->getText() == "schnorrfast")
        {
            lr.setSVPMode(LatticeReduction::SVP_SchnorrFast);
            config << ", fast Schnorr enumeration";
        }
        else if ((v->getText() == "listsieve") || (v->getText() == "ls"))
        {
            lr.setSVPMode(LatticeReduction::SVP_ListSieve);
            config << ", List-Sieve";
        }
        else if ((v->getText() == "listsievebirthday") || (v->getText() == "lsb"))
        {
            lr.setSVPMode(LatticeReduction::SVP_ListSieveBirthday);
            config << ", List-Sieve with Birthday";
        }
        else if ((v->getText() == "gausssieve") || (v->getText() == "gs"))
        {
            lr.setSVPMode(LatticeReduction::SVP_GaussSieve);
            config << ", Gauss-Sieve";
        }
        else if ((v->getText() == "voronoicell") || (v->getText() == "vc"))
        {
            lr.setSVPMode(LatticeReduction::SVP_VoronoiCellSVP);
            config << ", Voronoi Cell computation";
        }
        else
        {
            std::cerr << "Unknown enumeration mode \"" << v->getText() << "\"\n";
            std::cerr << "Choices: kse, parallelkse/pkse [default], schnorrfast, listsieve/ls, listsievebirthday/lsb, gausssieve/gs, voronoicell/vc\n";
            return -1;
        }
    }
    else
    {
        lr.setSVPMode(LatticeReduction::SVP_ParallelKannanSchnorrEuchner);
        config << ", parallel Kannan-Schnorr-Euchner enumeration";
    }
    // Threads
    unsigned maxthreads = 0;
    if ((v = args.getValue("maxthreads")) != NULL)
    {
        if (!v->isInteger() || (v->getInteger() < 0))
        {
            std::cerr << "Maximal thread number must be non-negative integer!\n";
            return -1;
        }
        maxthreads = v->getInteger();
        lr.setMaximalCoreUsage(maxthreads);
    }
    // Range
    signed begin = 0, end = -1;
    if ((v = args.getValue("begin")) != NULL)
    {
        if (!v->isInteger() || (v->getInteger() < 0))
        {
            std::cerr << "Begin should be non-negative!\n";
            return -1;
        }
        else
        {
            begin = v->getInteger();
            if (begin < 0)
                begin = 0;
            if ((unsigned)begin >= lr.rank())
                begin = lr.rank() - 1;
        }
    }
    if ((v = args.getValue("end")) != NULL)
    {
        if (!v->isInteger() || (v->getInteger() < begin))
        {
            std::cerr << "End should be at least " << begin << "!\n";
            return -1;
        }
        else
            end = v->getInteger();
    }
    if ((begin > 0) || (end >= 0))
    {
        if ((end < 0) || ((unsigned)end >= lr.rank()))
            end = lr.rank() - 1;
        if (end < begin)
            end = begin;
        lr.setRange(begin, end);
        config << ", range [" << begin << ", " << end << "]";
    }
    // MinTrack
    if ((v = args.getValue("mintracksilent")) != NULL)
    {
        lr.setMinCallbackFunction(MinCallbackFunctionSilent, MinCallbackFunctionSilent_LI);
        lr.setEnumCallbackFunction(EnumCallbackFunctionSilent, EnumCallbackFunctionSilent_LI);
    }
    if ((v = args.getValue("mintrack")) != NULL)
    {
        vv = args.getValue("mintracklog");
        InitMinCallbackFunction(lr, vv ? vv->getText() : "");
        lr.setMinCallbackFunction(MinCallbackFunction, MinCallbackFunction_LI);
        lr.setEnumCallbackFunction(EnumCallbackFunction, EnumCallbackFunction_LI);
        vv = args.getValue("mintrackincvecs");
        if (vv)
            mincallback_dumpvectors = true;
    }
    // Undocumented debug mode
    if ((v = args.getValue("fulldebug")) != NULL) // undocumented: does _A LOT_ of output
    {
        lr.setCallbackInterval(0);
        dumpLatticeEveryStep = true;
    }
    // Init timeres
    CPUTimer cputime;
    Timer wctime;
    InitLLLCallbackFunction();
    lr.setCallbackFunction(LLLCallbackFunction, LLLCallbackFunction_LI);
    // Preprocessing
    if ((v = args.getValue("preprocessing")) == NULL)
        v = args.getValue("pp");
    if (v != NULL)
    {
        std::ostringstream localmode;
        double alpha = 0.75;
        if (v->isInteger() || (v->getInteger() == 1))
            alpha = 1.0;
        else
            if (v->isFloat() && (v->getFloat() >= 0.5) && (v->getFloat() <= 1))
                alpha = v->getFloat();
        localmode << ", alpha=" << alpha;
        
        std::cout << "Running preprocessing LLL [" << config.str() << localmode.str() << "]...\n";
        lr.lll(alpha, LatticeReduction::LLL_Classic);
        outputShortest(lr); 
    }
    // Deep Insertion
    if ((v = args.getValue("deepinsertion")) == NULL)
        v = args.getValue("di");
    if (v != NULL)
    {
        lr.setDeepInsertionMethod(LatticeReduction::DI_Classic);
        std::string modestring = "Schnorr-Euchner";
        if ((vv = args.getValue("dialg")) != NULL)
        {
            if (vv->getText() == "classic")
            {
                lr.setDeepInsertionMethod(LatticeReduction::DI_Classic);
                modestring = "Schnorr-Euchner";
            }
            else if ((vv->getText() == "minimizepotential1") || (vv->getText() == "minpot1") || (vv->getText() == "minpot"))
            {
                lr.setDeepInsertionMethod(LatticeReduction::DI_MinimizePotential1);
                modestring = "potential minizming I";
            }
            else if ((vv->getText() == "minimizepotential2") || (vv->getText() == "minpot2"))
            {
                lr.setDeepInsertionMethod(LatticeReduction::DI_MinimizePotential2);
                modestring = "potential minizming II";
            }
            else
            {
                std::cerr << "Unknown deep insertion method \"" << vv->getText() << "\"\n";
                std::cerr << "Choices: classic [default], minimizepotential1/minpot1/minpot, minimizepotential2/minpot2\n";
                return -1;
            }
        }
        if ((v->getText() == "all") || (v->getText() == "di"))
        {
            lr.setDeepInsertionChoice(LatticeReduction::DIC_All);
            config << ", full " << modestring << " deep insertions";
        }
        else
        {
            unsigned bs = 1;
            if ((vv = args.getValue("diblocksize")) == NULL)
                vv = args.getValue("dibs");
            if (vv != NULL)
            {
                if (!vv->isInteger() || (vv->getInteger() <= 0))
                {
                    std::cerr << "Deep insertion blocksize be positive!\n";
                    return -1;
                }
                else
                    bs = vv->getInteger();
            }
            if (v->getText() == "first")
            {
                lr.setDeepInsertionChoice(LatticeReduction::DIC_First, bs);
                config << ", " << modestring << " deep insertions for first " << bs;
            }
            else if (v->getText() == "block")
            {
                lr.setDeepInsertionChoice(LatticeReduction::DIC_Block, bs);
                config << ", " << modestring << " deep insertions w/ blocksize " << bs;
            }
            else if (v->getText() == "both")
            {
                lr.setDeepInsertionChoice(LatticeReduction::DIC_FirstBlock, bs);
                config << ", " << modestring << " deep insertions for first and w/ blocksize " << bs;
            }
            else
            {
                std::cerr << "Unknown deep insertion choice \"" << v->getText() << "\"\n";
                std::cerr << "Choices: all [default], first, block, both\n";
                return -1;
            }
        }
        if ((vv = args.getValue("dimode")) != NULL)
        {
            if (vv->getText() == "before")
            {
                lr.setDeepInsertionMode(LatticeReduction::DIM_BeforeSR);
                config << ", do DI before SR";
            }
            else if (vv->getText() == "after")
            {
                lr.setDeepInsertionMode(LatticeReduction::DIM_AfterSR);
                config << ", do DI after SR";
            }
            else if (vv->getText() == "both")
            {
                lr.setDeepInsertionMode(LatticeReduction::DIM_Both);
                config << ", do DI before+after SR";
            }
            else
            {
                std::cerr << "Unknown deep insertion mode \"" << vv->getText() << "\"\n";
                std::cerr << "Choices: before, after [default], both\n";
                return -1;
            }
        }
        else
        {
            lr.setDeepInsertionMode(LatticeReduction::DIM_AfterSR);
            config << ", do DI after SR";
        }
    }
    // Annealing
    if ((v = args.getValue("annealing")) == NULL)
        v = args.getValue("anneal");
    if (v != NULL)
    {
        lr.setDefaultAnnealing();
        config << ", simulated annealing";
    }
    cputime.start();
    wctime.start();
    
    // Commands
    if ((v = args.getValue("sr")) != NULL)
    {
        std::cout << "Running size reduction [" << config.str() << "]...\n";
        {
            lr.sizereduction();
            if ((v = args.getValue("verify")) != NULL)
            {
              #ifndef PLLL_CONFIG_NO_GS_CLASSICINT
                lr.setGramSchmidt(LatticeReduction::G_ClassicInteger);
              #endif
              #ifndef PLLL_CONFIG_NO_ARITHMETIC_RATIONAL
                lr.setArithmetic(LatticeReduction::A_Rational);
              #endif
                verify(lr.isSizeReduced(), "size reduced basis");
            }
        }
        outputShortest(lr); 
    }
    if ((v = args.getValue("sort")) != NULL)
    {
        std::ostringstream localmode;
        bool proj = false;
        if ((v->getText() == "unprojected") || (v->getText() == "unproj") || (v->getText() == "sort"))
        {
            proj = false;
            localmode << ", unprojected";
        }
        else if ((v->getText() == "projected") || (v->getText() == "proj"))
        {
            proj = true;
            localmode << ", projected";
        }
        else
        {
            std::cerr << "Unknown sorting mode \"" << v->getText() << "\"\n";
            std::cerr << "Choices: unprojected/unproj [default], projected/proj\n";
            return -1;
        }
        
        std::cout << "Running sorting [" << config.str() << localmode.str() << "]...\n";
        {
            lr.sort(proj);
        }
        outputShortest(lr); 
    }
    if ((v = args.getValue("ssred")) != NULL)
    {
        std::ostringstream localmode;
        bool proj = false;
        if ((v->getText() == "unprojected") || (v->getText() == "unproj") || (v->getText() == "ssred"))
        {
            proj = false;
            localmode << ", unprojected";
        }
        else if ((v->getText() == "projected") || (v->getText() == "proj"))
        {
            proj = true;
            localmode << ", projected";
        }
        else
        {
            std::cerr << "Unknown sorting mode \"" << v->getText() << "\"\n";
            std::cerr << "Choices: unprojected/unproj [default], projected/proj\n";
            return -1;
        }
        
        std::cout << "Running sorting+sizered reduction [" << config.str() << localmode.str() << "]...\n";
        {
            while (true)
            {
                bool modified = false;
                linalg::math_matrix<arithmetic::Integer> old = lr.getLattice();
                std::cout << "Size reduction...\n";
                lr.sizereduction();
                if (old != lr.getLattice())
                    modified = true;
                std::cout << "Sorting...\n";
                lr.sort(proj);
                if (!modified)
                    if (old == lr.getLattice())
                        break;
                LLLCallbackFunction(lr.getLattice());
            }
        }
        outputShortest(lr); 
    }
    if ((v = args.getValue("lll")) != NULL)
    {
        std::ostringstream localmode;
        LatticeReduction::LLLMode lllmode; 
        if ((v->getText() == "classic") || (v->getText() == "lll"))
        {
            lllmode = LatticeReduction::LLL_Classic;
            localmode << ", classic LLL";
        }
        else if ((v->getText() == "unprojected") || (v->getText() == "unproj") || (v->getText() == "up"))
        {
            lllmode = LatticeReduction::LLL_Unprojected;
            localmode << ", unprojected LLL";
        }
        else if ((v->getText() == "siegel") || (v->getText() == "s"))
        {
            lllmode = LatticeReduction::LLL_Siegel;
            localmode << ", LLL with Siegel condition";
        }
        else
        {
            std::cerr << "Unknown LLL mode \"" << v->getText() << "\"\n";
            std::cerr << "Choices: classic [default], unprojected/unproj/up, siegel/s\n";
            return -1;
        }
        
        double alpha = 0.99;
        if ((v = args.getValue("alpha")) != NULL)
        {
            if (v->isInteger() || (v->getInteger() == 1))
                alpha = 1.0;
            else
            {
                if (!v->isFloat() || (v->getFloat() < 0.25) || (v->getFloat() > 1))
                {
                    std::cerr << "alpha must be real number in [0.25, 1]!\n";
                    return -1;
                }
                alpha = v->getFloat();
            }
        }
        localmode << ", alpha=" << alpha;
        
        std::cout << "Running LLL [" << config.str() << localmode.str() << "]...\n";
        {
            lr.lll(alpha, lllmode);
            if ((v = args.getValue("verify")) != NULL)
            {
              #ifndef PLLL_CONFIG_NO_GS_CLASSICINT
                lr.setGramSchmidt(LatticeReduction::G_ClassicInteger);
              #endif
              #ifndef PLLL_CONFIG_NO_ARITHMETIC_RATIONAL
                lr.setArithmetic(LatticeReduction::A_Rational);
              #endif
                verify(lr.isLLLBasis(alpha, lllmode), "LLL basis");
            }
        }
        outputShortest(lr); 
    }
    if ((v = args.getValue("bkz")) != NULL)
    {
        std::ostringstream localmode;
        std::string name;
        LatticeReduction::BKZMode bkzmode; 
        if ((v->getText() == "schnorreuchner") || (v->getText() == "bkz"))
        {
            bkzmode = LatticeReduction::BKZ_SchnorrEuchner;
            name = "Schnorr-Euchner BKZ";
        }
        else if (v->getText() == "simplified")
        {
            bkzmode = LatticeReduction::BKZ_Simplified;
            name = "Simplified BKZ";
        }
        else if ((v->getText() == "hps-hkz") || (v->getText() == "hanrotpujolstehle-hkz"))
        {
            bkzmode = LatticeReduction::BKZ_HanrotPujolStehleHKZ;
            name = "Hanrot-Pujol-Stehle BKZ (HKZ)";
        }
        else if ((v->getText() == "hps-svp") || (v->getText() == "hanrotpujolstehle-svp"))
        {
            bkzmode = LatticeReduction::BKZ_HanrotPujolStehleSVP;
            name = "Hanrot-Pujol-Stehle BKZ (SVP)";
        }
        else if ((v->getText() == "pd") || (v->getText() == "primaldual"))
        {
            bkzmode = LatticeReduction::BKZ_PrimalDual;
            name = "Koy's Primal-Dual BKZ";
        }
        else if ((v->getText() == "s") || (v->getText() == "slide"))
        {
            bkzmode = LatticeReduction::BKZ_SlideReduction;
            name = "Slide Reduction";
        }
        else if ((v->getText() == "is") || (v->getText() == "improvedslide"))
        {
            bkzmode = LatticeReduction::BKZ_ImprovedSlideReduction;
            name = "Improved Slide Reduction";
        }
        else if ((v->getText() == "is2") || (v->getText() == "improvedslide2"))
        {
            bkzmode = LatticeReduction::BKZ_ImprovedSlideReduction2;
            name = "Improved Slide Reduction (larger DSVP)";
        }
        else if ((v->getText() == "is3") || (v->getText() == "improvedslide3"))
        {
            bkzmode = LatticeReduction::BKZ_ImprovedSlideReduction3;
            name = "Improved Slide Reduction (larger SVP)";
        }
        else if ((v->getText() == "sb") || (v->getText() == "semiblock2k"))
        {
            bkzmode = LatticeReduction::BKZ_SemiBlock2k;
            name = "Semi Block 2k Reduction";
        }
        else if ((v->getText() == "sr") || (v->getText() == "samplingreduction"))
        {
            bkzmode = LatticeReduction::BKZ_SamplingReduction;
            name = "Sampling Reduction";
        }
        else if ((v->getText() == "ex") || (v->getText() == "experimental"))
        {
            bkzmode = LatticeReduction::BKZ_Experimental;
            name = "Experimental Block Reduction";
        }
        else
        {
            std::cerr << "Unknown BKZ mode \"" << v->getText() << "\"\n";
            std::cerr << "Choices: schnorreuchner [default], simplified, hanrotpujolstehle-hkz/hps-hkz, hanrotpujolstehle-svp/hps-svp, primaldual/pd, slide/s, improvedslide/is, improvedslide[23]/is[23], semiblock2k/sb, samplingreduction/sr, experimental/ex\n";
            return -1;
        }
        localmode << ", " << name;
        
        double alpha = 0.99;
        if ((v = args.getValue("alpha")) != NULL)
        {
            if (v->isInteger() || (v->getInteger() == 1))
                alpha = 1.0;
            else
            {
                if (!v->isFloat() || (v->getFloat() < 0.25) || (v->getFloat() > 1))
                {
                    std::cerr << "alpha must be real number in [0.25, 1]!\n";
                    return -1;
                }
                alpha = v->getFloat();
            }
        }
        localmode << ", alpha=" << alpha;
        
        int bkz_window_size = 10;
        if ((v = args.getValue("window")) == NULL)
            if ((v = args.getValue("windowsize")) == NULL)
                if ((v = args.getValue("bs")) == NULL)
                    v = args.getValue("blocksize");
        if (v != NULL)
        {
            if (!v->isInteger() || (v->getInteger() < 1))
            {
                std::cerr << "BKZ window size must be positive integer!\n";
                return -1;
            }
            bkz_window_size = v->getInteger();
            if (bkz_window_size < 2)
                bkz_window_size = 2;
        }
        localmode << ", window size " << bkz_window_size;
        
        std::cout << "Running BKZ [" << config.str() << localmode.str() << "]...\n";
        {
            lr.bkz(alpha, bkz_window_size, bkzmode);
            if ((v = args.getValue("verify")) != NULL)
            {
              #ifndef PLLL_CONFIG_NO_GS_CLASSICINT
                lr.setGramSchmidt(LatticeReduction::G_ClassicInteger);
              #endif
              #ifndef PLLL_CONFIG_NO_ARITHMETIC_RATIONAL
                lr.setArithmetic(LatticeReduction::A_Rational);
              #endif
                verify(lr.isBKZBasis(alpha, bkz_window_size, bkzmode), name + " basis");
            }
        }
        outputShortest(lr); 
    }
    if ((v = args.getValue("hkz")) != NULL)
    {
        std::cout << "Running HKZ [" << config.str() << "]...\n";
        {
            lr.hkz();
            if ((v = args.getValue("verify")) != NULL)
            {
              #ifndef PLLL_CONFIG_NO_GS_CLASSICINT
                lr.setGramSchmidt(LatticeReduction::G_ClassicInteger);
              #endif
              #ifndef PLLL_CONFIG_NO_ARITHMETIC_RATIONAL
                lr.setArithmetic(LatticeReduction::A_Rational);
              #endif
                verify(lr.isHKZBasis(), "HKZ basis");
            }
        }
        outputShortest(lr); 
    }
    if ((v = args.getValue("dhkz")) != NULL)
    {
        std::cout << "Running Dual HKZ [" << config.str() << "]...\n";
        {
            lr.hkz(true);
            if ((v = args.getValue("verify")) != NULL)
            {
              #ifndef PLLL_CONFIG_NO_GS_CLASSICINT
                lr.setGramSchmidt(LatticeReduction::G_ClassicInteger);
              #endif
              #ifndef PLLL_CONFIG_NO_ARITHMETIC_RATIONAL
                lr.setArithmetic(LatticeReduction::A_Rational);
              #endif
                verify(lr.isHKZBasis(true), "HKZ basis");
            }
        }
        outputShortest(lr); 
    }
    if ((v = args.getValue("svp")) != NULL)
    {
        std::cout << "Running SVP [" << config.str() << "]...\n";
        {
            bool basis = true, extreme = false;
            if ((v->getText() == "svp") || (v->getText() == "normal"))
                extreme = false;
            else if ((v->getText() == "ex") || (v->getText() == "extreme"))
                extreme = true;
            else
            {
                std::cerr << "Unknown SVP command \"" << v->getText() << "\"\n";
                std::cerr << "Choices: normal [default], extreme/ex\n";
                return -1;
            }
            if ((v = args.getValue("nosvpbasis")) != NULL)
                basis = false;
            lr.svp(basis, extreme);
            if ((v = args.getValue("verify")) != NULL)
            {
              #ifndef PLLL_CONFIG_NO_GS_CLASSICINT
                lr.setGramSchmidt(LatticeReduction::G_ClassicInteger);
              #endif
              #ifndef PLLL_CONFIG_NO_ARITHMETIC_RATIONAL
                lr.setArithmetic(LatticeReduction::A_Rational);
              #endif
                verify(lr.isSVPBasis(), "SVP basis");
            }
        }
        outputShortest(lr); 
    }
    if ((v = args.getValue("dsvp")) != NULL)
    {
        std::cout << "Running Dual-SVP [" << config.str() << "]...\n";
        {
            bool basis = true, extreme = false;
            if ((v->getText() == "svp") || (v->getText() == "normal"))
                extreme = false;
            else if ((v->getText() == "ex") || (v->getText() == "extreme"))
                extreme = true;
            else
            {
                std::cerr << "Unknown SVP command \"" << v->getText() << "\"\n";
                std::cerr << "Choices: normal [default], extreme/ex\n";
                return -1;
            }
            if ((v = args.getValue("nosvpbasis")) != NULL)
                basis = false;
            lr.svp(basis, extreme, true);
            if ((v = args.getValue("verify")) != NULL)
            {
              #ifndef PLLL_CONFIG_NO_GS_CLASSICINT
                lr.setGramSchmidt(LatticeReduction::G_ClassicInteger);
              #endif
              #ifndef PLLL_CONFIG_NO_ARITHMETIC_RATIONAL
                lr.setArithmetic(LatticeReduction::A_Rational);
              #endif
                verify(lr.isSVPBasis(true), "SVP basis");
            }
        }
        outputShortest(lr); 
    }
    if ((v = args.getValue("computedeterminant")) != NULL)
    {
        computeDeterminant(lr);
        std::cout << "Lattice has " << (determinant_is_square ? "square " : "") << "determinant " << *determinant << "\n";
    }
    if ((v = args.getValue("info")) != NULL)
    {
        lr.forceGSRebuild(true);
        linalg::math_rowvector<arithmetic::Real> q(lr.rank() ? lr.rank() - 1 : 0);
        arithmetic::Real min, max;
        for (unsigned i = 0; i < q.size(); ++i)
        {
            q[i] = lr.getGSSqNormR(i + 1) / lr.getGSSqNormR(i);
            if (i == 0)
                min = max = q[i];
            else
            {
                if (min > q[i])
                    min = q[i];
                if (max < q[i])
                    max = q[i];
            }
        }
        std::cout << "q = " << q << "\n";
        std::cout << "minimal q: " << arithmetic::convert<double>(min)
                  << ", maximal q: " << arithmetic::convert<double>(max) << "\n";
        
        // Compute the "estimate" for the length of the shortest vector
        // 1.05 * Gamma(n/2 + 1)^(1/n) / sqrt(pi) * (det(lattice))^(1/n)
        if ((determinant == NULL) || ((v = args.getValue("force")) != NULL))
            computeDeterminant(lr);
        std::cout << "+===========================================================================+\n";
        std::cout << "| Lattice has " << (determinant_is_square ? "square " : "") << "determinant " << *determinant << "\n";
        arithmetic::RealContext rc;
        arithmetic::Real est = arithmetic::convert(1.05l / std::sqrt(std::atan(1.0l) * 4.0l), rc);
        arithmetic::Real r = arithmetic::convert(*determinant, rc);
        if (determinant_is_square)
            r = sqrt(r);
        else
            r = abs(r);
        r *= arithmetic::convert(tgamma(((long double)lr.rank()) / 2.0l + 1.0l), rc);
        est *= power(r, arithmetic::convert(1, rc) / arithmetic::convert(lr.rank(), rc));
        std::cout << "| Shortest vector has estimated length <= " << est << "\n";
        std::cout << "| Shortest vector has estimated square length <= " << est * est << "\n";
        std::cout << "+===========================================================================+\n";
    }
    if ((v = args.getValue("gsinfo")) != NULL)
    {
        lr.forceGSRebuild(true);
        std::cout << "Gram-Schmidt squared norms:\n";
        for (unsigned i = 0; i < lr.rank(); ++i)
            std::cout << "  " << std::setw(3) << i << ": " << lr.getGSSqNormR(i) << "\n";
        std::cout << "Gram-Schmidt coefficients:\n";
        for (unsigned i = 0; i < lr.rank(); ++i)
        {
            std::cout << "  " << std::setw(3) << i << ":";
            for (unsigned j = 0; j < i; ++j)
                std::cout << " " << lr.getGSCoefficientR(i, j);
            std::cout << "\n";
        }
    }
    if ((v = args.getValue("test")) != NULL)
    {
        linalg::math_rowvector<arithmetic::Integer> vec = parseVector(v->getText());
        if (vec.size() != lr.dimension())
        {
            std::cerr << "Vector has invalid dimension (" << vec.size() << ")!\n";
            return -1;
        }
        else
        {
            arithmetic::RealContext rc;
            std::cout << vec << " (squared norm " << normSq(vec) << "; norm " << sqrt(arithmetic::convert(normSq(vec), rc)) << ")\n";
            linalg::math_colvector<arithmetic::Integer> res = solveUniqInt(lr.getLattice().transpose(), vec.transpose());
            if (res.size() == 0)
                std::cout << "No solution found!\n";
            else
            {            
                std::cout << "Coefficient vector: " << res << "\n";
                assert(res * lr.getLattice() == vec);
            }
        }
    }
    if ((v = args.getValue("hnf")) != NULL)
    {
        std::cout << "Computing Hermite Normal Form of basis...\n";
        linalg::math_matrix<arithmetic::Integer> A = lr.getLattice();
        if (v->getText() == "right")
            for (unsigned i = 0; i < A.cols()/2; ++i)
                linalg::swap(A.col(i), A.col(A.cols() - 1 - i));
        hnf(A);
        if (v->getText() == "right")
        {
            for (unsigned i = 0; i < A.cols()/2; ++i)
                linalg::swap(A.col(i), A.col(A.cols() - 1 - i));
            for (unsigned i = 0; i < A.rows()/2; ++i)
                linalg::swap(A.row(i), A.row(A.rows() - 1 - i));
        }
        lr.setLattice(A);
        outputShortest(lr); 
    }
    
    // Close minimal vector length log (if open) and stop timers
    cputime.stop();
    wctime.stop();
    MinCallbackLogClose();
    
    // Write output
    if ((v = args.getValue("write")) != NULL)
    {
        if (v->isEmpty())
        {
            std::cerr << "Error: argument -write needs output file name (-write=blah)\n";
            return -1;
        }
        std::cout << "Storing as '" << v->getText() << "'...\n";
        if (!storeLattice(lr, v->getText()))
            std::cerr << "Error while storing lattice!\n";
        else
            std::cout << "  done.\n";
    }
    
    std::cout << "\n==LATTICE====================================================\n";
    std::cout << "Dimension: " << lr.dimension() << "\n";
    std::cout << "Rank: " << lr.rank() << "\n";
    if ((v = args.getValue("dump")) != NULL)
    {
        if (v->getText() == "ntl")
            writeLatticeNTL(std::cout, lr.getLattice());
        else
            std::cout << lr.getLattice() << "\n";
    }
    
    // Anything left?
    if (args.hasUnprocessedArguments())
    {
        std::cerr << "\nWARNING: Ignoring the following arguments:";
        args.enumerateUnprocessedArguments(argsEnum);
        std::cerr << "\n";
    }
    
    // Output statistics
    const LatticeReduction::Statistics & stats = lr.getStatistics();
    std::cout << "\n==STATISTICS=================================================\n";
    std::cout << "Swaps:                         " << stats.swaps << "\n";
    std::cout << "Adds:                          " << stats.adds << " (+-1: " << stats.adds_pm1;
    if (stats.adds)
        std::cout << " [" << std::setprecision(4) << 100.0 * stats.adds_pm1 / stats.adds << "%]";
    std::cout << ", +-2: " << stats.adds_pm2;
    if (stats.adds)
        std::cout << " [" << std::setprecision(4) << 100.0 * stats.adds_pm2 / stats.adds << "%]";
    std::cout << ")\n";
    std::cout << "Flips:                         " << stats.flips << "\n";
    std::cout << "Generic 2x2 transformations:   " << stats.trans << "\n";
    std::cout << "Size reductions:               " << stats.sizereductions << "\n";
    std::cout << "Deep insertions:               " << stats.deepinsertions << "\n";
    std::cout << "Calls to enumeration:          " << stats.enumcalls << "\n";
    std::cout << "Failed calls to enumeration:   " << stats.enumfails;
    if (stats.enumcalls)
        std::cout << " [" << std::setprecision(4) << 100.0 * stats.enumfails / stats.enumcalls << "%]";
    std::cout << "\n";
    std::cout << "Vector insertions:             " << stats.vectorinsertions;
    if (stats.enumcalls)
        std::cout << " [" << std::setprecision(4) << 100.0 * stats.vectorinsertions / stats.enumcalls << "%]";
    std::cout << "\n";
    std::cout << "Vector insertions (rearrange): " << stats.vectorinsertions_rearrange;
    if (stats.enumcalls)
        std::cout << " [" << std::setprecision(4) << 100.0 * stats.vectorinsertions_rearrange / stats.enumcalls << "%]";
    std::cout << "\n";
    std::cout << "Total wallclock time:          " << std::setprecision(6) << (long double)wctime.elapsed() / (long double)Timer::d_one_timeunit << " " << Timer::d_timeunit_name << "\n";
    std::cout << "Total CPU time:                " << std::setprecision(6) << (long double)cputime.elapsed() / (long double)CPUTimer::d_one_timeunit << " " << CPUTimer::d_timeunit_name << "\n";
    
    if (determinant)
    {
        delete determinant;
        determinant = NULL;
    }
}
