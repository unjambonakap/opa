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

#ifndef PLLL_INCLUDE_GUARD__BUFFER_HPP
#define PLLL_INCLUDE_GUARD__BUFFER_HPP

namespace plll
{
    template<class T>
    class Buffer
    {
    private:
        T * d_data;
        unsigned d_size, d_used;
        
        void reserve(unsigned minadd)
        {
            unsigned minsize = d_used + minadd;
            // Make multiple of 256 = 2^8
            minsize += 0xFF;
            minsize &= ~0xFF;
            // Reserve enough memory
            T * d = new T[minsize];
            for (unsigned i = 0; i < d_used; ++i)
                d[i] = d_data[i];
            d_size = minsize;
            delete[] d_data;
            d_data = d;
        }
        
    public:
        Buffer()
            : d_data(NULL), d_size(0), d_used(0)
        {
        }
        
        ~Buffer()
        {
            delete[] d_data;
        }
        
        const T * data() const
        {
            return d_data;
        }
        
        unsigned size() const
        {
            return d_used;
        }
        
        void clear()
        {
            d_used = 0;
        }
        
        void add(const T & c)
        {
            if (d_used + 1 >= d_size)
                reserve(1);
            d_data[d_used++] = c;
        }
        
        void add(const T * c, unsigned n)
        {
            if (d_used + n >= d_size)
                reserve(n);
            while (n)
            {
                d_data[d_used++] = *c++;
                --n;
            }
        }
    };
}

#endif
