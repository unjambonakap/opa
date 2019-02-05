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
#include <plll/arguments.hpp>
#include <cassert>
#include <cstring>
#include <cstdlib>
#include <cctype>
#include <cmath>
#include <limits>
#include <iostream>

namespace plll
{
    namespace helper
    {
        ArgumentParser::Value::Value(const std::string & value, std::string & id)
            : d_value(value), d_type(VT_Empty), d_int(0), d_float(0.0), d_id_overflow(false)
              // value is of format "x=y" (VT_Text, VT_Integer, VT_Float) or just "x" (VT_Empty)
        {
            size_t i = value.find('=');
            if (i == std::string::npos)
            {
                // No value
                id = value;
            }
            else
            {
                // First: treat as text
                id = value.substr(0, i);
                d_value = value.substr(i + 1);
                d_type = VT_Text;
                
                // Next: identify numerical data
                bool couldBeInteger = d_value.size() > 0, couldBeFloat = d_value.size() > 0;
                unsigned i = 0;
                if ((d_value[i] == '-') || (d_value[i] == '+'))
                {
                    if (d_value.size() == 1)
                        couldBeInteger = couldBeFloat = false;
                    else
                        ++i;
                }
                int mantissadigits = 0; // number of (decimal) digits in mantissa
                int mantissaexp = 0; // absolute value m of mantissa satisfies 10^{mantissaexp-1} <= m < 10^{mantissaexp}
                int exponent = 0; // in base 10
                bool decdot = false;
                while ((i < d_value.size()) && (couldBeInteger || couldBeFloat))
                    // Test for mantissa
                {
                    if (std::isdigit(d_value[i]))
                    {
                        if (mantissaexp || (d_value[i] != '0'))
                            ++mantissaexp;
                        ++mantissadigits;
                        ++i;
                    }
                    else
                    {
                        couldBeInteger = false;
                        if ((d_value[i] == '.') && !decdot)
                        {
                            ++i;
                            decdot = true;
                        }
                        else if ((d_value[i] == 'e') || (d_value[i] == 'E'))
                        {
                            ++i;
                            break;
                        }
                        else
                            couldBeFloat = false;
                    }
                }
                if (couldBeFloat && !couldBeInteger)
                {
                    // Test for exponent
                    bool exppos = true;
                    if ((d_value[i] == '-') || (d_value[i] == '+'))
                    {
                        if (d_value.size() == i + 1)
                            couldBeFloat = false;
                        if (d_value[i] == '-')
                            exppos = false;
                    }
                    while ((i < d_value.size()) && couldBeFloat)
                    {
                        if (std::isdigit(d_value[i]))
                        {
                            exponent = 10 * exponent + (int)(d_value[i] - 10);
                            ++i;
                        }
                        else
                            couldBeFloat = false;
                    }
                    if (!exppos)
                        exponent = -exponent;
                }
                
                // Mark & convert
                if (couldBeInteger)
                {
                    d_type = VT_Integer;
                    d_int = std::atoi(d_value.c_str());
                    // Test for overflow; very basic test, not completely trustable
                    d_id_overflow = (mantissaexp > std::numeric_limits<int>::digits10) ||
                        (d_int == std::numeric_limits<int>::max()) || (d_int == std::numeric_limits<int>::min()); // correct way according to documentation on std::atoi
                }
                else if (couldBeFloat)
                {
                    d_type = VT_Float;
                    d_float = std::atof(d_value.c_str());
                    // Test for overflow; very basic test, not completely trustable
                    d_id_overflow = (exponent + mantissaexp > std::numeric_limits<double>::max_exponent10) || (exponent + mantissaexp - 1 < std::numeric_limits<double>::min_exponent10) ||
                        (d_float == HUGE_VAL); // correct way according to documentation of std::atof
                }
            }
        }

        ArgumentParser::ArgumentParser(int argc, char **argv)
        {
            for (int i = 1; i < argc; ++i)
            {
                assert(std::strlen(argv[i]) > 0);
        
                if ((std::strcmp(argv[i], "--") == 0) && (i < argc - 1))
                    d_names.push_back(argv[++i]);
                else if (argv[i][0] == '-')
                {
                    std::string id;
                    Value val(argv[i] + ((argv[i][1] != '-') ? 1 : 2), id);
                    std::pair<std::map<std::string, Value>::iterator, bool> i = d_args.insert(std::make_pair(id, val));
                    if (i.second)
                        d_unprocessed_args.insert(i.first);
                    else
                        std::cerr << "WARNING: Argument \"" << id << "\" specified twice. Ignoring second value.\n";
                }
                else
                    d_names.push_back(argv[i]);
            }
        }
    }
}
