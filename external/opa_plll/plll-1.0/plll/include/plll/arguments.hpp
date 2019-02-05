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

#ifndef PLLL_INCLUDE_GUARD__ARGUMENTS_HPP
#define PLLL_INCLUDE_GUARD__ARGUMENTS_HPP

#include <map>
#include <set>
#include <list>
#include <string>

/**
   \file
   \brief Command line argument parsing.
   
   This header provides a simple command line argument parsing class used for the `plll` main
   program.
*/
namespace plll
{
    namespace helper
    {
        class ArgumentParser
        /**
           \brief Parses command line arguments and makes them accessible.
           
           Typically, the command line of a command line program is an array of strings:
           
               "<executable name>" "-value=1" "-R" "name1" "-h=" "name2" "--" "-name3", "--debug=all"
               
           The `ArgumentParser` class' constructor accepts such a list of C strings as provided to
           the `main()` function, and splits it up in a list of names, in this case "name1", "name2"
           and "-name3" (because "-name3" is preceeded with "--", it is not treated as a value), and
           a list of values with lables. Here, the values with their labels are:
           
           - `"value"` -> `"1"`,
           - `"R"` -> (empty),
           - `"h"` -> `""`,
           - `"-debug"` -> `"all"`.
           
           Integers and floating point numbers are detected, and it is tried to parse them (as
           `long` or `double`), with a warning system (`Value::hasOverflow()` returns `true` if the
           value apparently does not fit in these types).
        */
        {
        public:
            class Value
            /**
               \brief Stores the value of one command line argument.
            */
            {
                friend class ArgumentParser;
                
            private:
                enum Type { VT_Empty, VT_Text, VT_Integer, VT_Float };
                
                // VT_Empty is used for arguments of type "--blah". The argument "--blah=" is
                // treated as VT_Text with getText().size() == 0.
                
                std::string d_value;
                Type d_type;
                long d_int;
                double d_float;
                bool d_id_overflow; // d_type is VT_Integer or VT_Float, but value is too large,
                                    // whence d_int/d_float is not meaningful and d_value should be
                                    // used
                
                Value(const std::string & value, std::string & id);
                
            public:
                /**
                   \brief Tests whether value is empty.
                   
                   Allows to test whether the command line argument had a value or not. An argument
                   of the type `--blah` will have `isEmpty() == true`, while an argument of the type
                   `--blah=` will have `isEmpty() == false`.
                */
                bool isEmpty() const
                {
                    return d_type == VT_Empty;
                }
                
                /**
                   \brief Returns whether the value is text.
                   
                   Informs whether the command line argument's value is not an integer, floating
                   point number or anything else which could be parsed.
                */
                bool isText() const
                {
                    return d_type == VT_Text;
                }
                
                /**
                   \brief Returns the raw text value.
                   
                   Returns the raw text value of the command line argument, i.e. the part after `=`
                   in `--blah=value`. In case no `=` was found, `-blah` (for `--blah`) or `blah`
                   (for `-blah`) itself is returned.
                */
                const std::string & getText() const
                {
                    return d_value;
                }
                
                /**
                   \brief Checks for overflows.
                   
                   Informs in case of an integer or floating point number whether an overflow
                   occured while casting to a native type. In case this returns `true`, one should
                   use `getText()` and parse the result using `plll::arithmetic::Integer`
                   respectively `plll::arithmetic::Real`.
                */
                bool hasOverflow() const
                {
                    return d_id_overflow;
                }
                
                /**
                   \brief Informs whether the result is an integer.
                */
                bool isInteger() const
                {
                    return d_type == VT_Integer;
                }
                
                /**
                   \brief Return integer value.
                   
                   Returns the (native) integer value, or 0 in case `isInteger()` returns
                   `false`. Undefined in case of `hasOverflow() == true`.
                */
                long getInteger() const
                {
                    return d_int;
                }
                
                /**
                   \brief Informs whether the result is a floating point number.
                   
                   Informs whether the result is a floating point number.
                */
                bool isFloat() const
                {
                    return d_type == VT_Float;
                }
                
                /**
                   \brief Return floating point value.
                   
                   Returns the (native) floating point value, or 0 in case `isFloat()` returns
                   `false`. Undefined in case of `hasOverflow() == true`.
                */
                double getFloat() const
                {
                    return d_float;
                }
            };
            
        private:
            class MapIteratorComparator
            {
            public:
                bool operator() (const std::map<std::string, Value>::const_iterator & A,
                                 const std::map<std::string, Value>::const_iterator & B) const
                {
                    std::less<const Value*> less;
                    return less(&(A->second), &(B->second));
                }
            };
            
            std::map<std::string, Value> d_args;
            std::list<std::string> d_names;
            mutable std::set<std::map<std::string, Value>::const_iterator, MapIteratorComparator> d_unprocessed_args;
            
        public:
            /**
               \brief Parses the command line arguments. Ignores the name of the current program.
            */
            ArgumentParser(int argc, char **argv);
            
            /**
               \brief Quick tests whether arguments were given.
            */
            bool hasNoArguments() const
            {
                return d_args.empty() && d_names.empty();
            }
            
            /**
               \brief Retrieves a list of names, i.e. arguments not starting with `-` or preceeded with `--`.
            */
            const std::list<std::string> & getNames() const
            {
                return d_names;
            }
            
            /**
               \brief Looks up the value of a label.
               
               Retrieves a `plll::ArgumentParser::Value` object for the given argument label `arg`
               if it was found, or `NULL` in case the argument was not found.
            */
            const Value * getValue(const std::string & arg) const
            {
                std::map<std::string, Value>::const_iterator i = d_args.find(arg);
                if (i == d_args.end())
                    return NULL;
                else
                {
                    std::set<std::map<std::string, Value>::const_iterator, MapIteratorComparator>::iterator it = d_unprocessed_args.find(i);
                    if (it != d_unprocessed_args.end())
                        d_unprocessed_args.erase(it);
                    return &(i->second);
                }
            }
            
            /**
               \brief Tests for unchecked label/values pairs.
               
               Returns `true` if arguments exists whose existence haven't been checked with
               `getValue()`.
            */
            bool hasUnprocessedArguments() const
            {
                return d_unprocessed_args.size();
            }
            
            /**
               \brief Calls a functor for every unchecked label/values pair.
               
               Calls the functor `f` for every argument not checked yet with `getValue()`. The
               functor `f` is called for every such argument with parameters of type `std::string`
               for the argument name and `plll::ArgumentParser::Value *` for its value.
            */
            template<class Fun>
            void enumerateUnprocessedArguments(Fun f) const
            {
                for (std::set<std::map<std::string, Value>::const_iterator, MapIteratorComparator>::iterator
                         i = d_unprocessed_args.begin(); i != d_unprocessed_args.end(); ++i)
                    f((*i)->first, (*i)->second);
            }
        };
    }
}

#endif
