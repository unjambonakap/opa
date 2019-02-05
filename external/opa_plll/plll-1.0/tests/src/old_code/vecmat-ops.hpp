#ifndef PLLL_INCLUDE_GUARD__VECMAT_OPS_HPP
#define PLLL_INCLUDE_GUARD__VECMAT_OPS_HPP

namespace plll
{
    namespace linalgold
    {
        /*
          Operations between vectors (and vectors) yielding vectors:
              add, sub, neg
          Operations between vectors and scalars yielding vectors:
              multiply, divide, modulo
        */
        
        template<class Container, class ScT>
        class VecWrapper
        {
        private:
            const Container & d_vec;
    
        public:
            enum { isVector = true };
            typedef ScT ScalarType;
    
            VecWrapper(const Container & vec)
                : d_vec(vec)
            {
            }
    
            inline operator math_vector<ScalarType> () const
            {
                return d_vec;
            }

            template<class V>
            inline void assignTo(V & x) const
            {
                x = d_vec;
            }
    
            inline math_vector<ScalarType> evaluate() const
            {
                return d_vec;
            }
    
            inline unsigned size() const
            {
                return d_vec.size();
            }
    
            inline const ScalarType & evalComponent(unsigned i) const
            {
                return d_vec[i];
            }
    
            inline const void evalComponentTo(ScalarType & dest, unsigned i) const
            {
                dest = d_vec[i];
            }
    
            inline const void evalComponentAddTo(ScalarType & dest, unsigned i) const
            {
                dest += d_vec[i];
            }
    
            inline const void evalComponentSubFrom(ScalarType & dest, unsigned i) const
            {
                dest -= d_vec[i];
            }

            template<class V>
            inline bool usesVec(const V & v) const
            {
                return (void*)&v == (void*)&d_vec;
            }
    
            template<class V>
            inline bool usesVecNF(const V & v) const
            {
                return false;
            }
        };

        template<class Op, class Data, class ScT>
        class VecExpression
        {
        private:
            Op d_op;
            Data d_data;
    
        public:
            enum { isVector = true };
            typedef ScT ScalarType;
    
            VecExpression(const Op & op, const Data & data)
                : d_op(op), d_data(data)
            {
            }
    
            inline operator math_vector<ScalarType> () const
            {
                math_vector<ScalarType> x(size());
                d_op.assignTo(x, d_data);
                return x;
            }
    
            template<class V>
            inline void assignTo(V & x) const
            {
                d_op.assignTo(x, d_data);
            }
    
            inline math_vector<ScalarType> evaluate() const
            {
                return d_op.evaluate(d_data);
            }
    
            unsigned size() const
            {
                return d_op.size(d_data);
            }
    
            ScalarType evalComponent(unsigned i) const
            {
                return d_op.evalComponent(d_data, i);
            }
    
            inline const void evalComponentTo(ScalarType & dest, unsigned i) const
            {
                d_op.evalComponentTo(dest, d_data, i);
            }
    
            inline const void evalComponentAddTo(ScalarType & dest, unsigned i) const
            {
                d_op.evalComponentAddTo(dest, d_data, i);
            }
    
            inline const void evalComponentSubFrom(ScalarType & dest, unsigned i) const
            {
                d_op.evalComponentSubFrom(dest, d_data, i);
            }
    
            inline ScT operator[] (unsigned i) const
            {
                return d_op.evalComponent(d_data, i);
            }
    
            template<class V>
            inline bool usesVec(const V & v) const
            {
                return d_op.usesVec(v, d_data);
            }
    
            template<class V>
            inline bool usesVecNF(const V & v) const
            {
                return d_op.usesVecNF(v, d_data);
            }

            inline ScT normSq() const
            {
                ScT r;
                setZero(r);
                for (unsigned i = 0; i < size(); ++i)
                    r += square(d_op.evalComponent(d_data, i));
                return r;
            }
        };

        template<class ScT>
        class VecAddOp
        {
        public:
            template<class A, class B, class C>
            inline void assignTo(A & d, const std::pair<B, C> & data) const
            {
                unsigned s = size(data);
                for (unsigned i = 0; i < s; ++i)
                    d[i] = data.first.evalComponent(i) + data.second.evalComponent(i);
            }
    
            template<class B, class C>
            inline math_vector<ScT> evaluate(const std::pair<B, C> & data) const
            {
                math_vector<ScT> res(size());
                assignTo(res, data);
            }
    
            template<class B, class C>
            inline unsigned size(const std::pair<B, C> & data) const
            {
                return data.first.size();
            }
    
            template<class B, class C>
            inline ScT evalComponent(const std::pair<B, C> & data, unsigned i) const
            {
                return data.first.evalComponent(i) + data.second.evalComponent(i);
            }
    
            template<class B, class C>
            inline void evalComponentTo(ScT & dest, const std::pair<B, C> & data, unsigned i) const
            {
                data.first.evalComponentTo(dest, i); 
                data.second.evalComponentAddTo(dest, i);
            }
    
            template<class B, class C>
            inline void evalComponentAddTo(ScT & dest, const std::pair<B, C> & data, unsigned i) const
            {
                data.first.evalComponentAddTo(dest, i); 
                data.second.evalComponentAddTo(dest, i);
            }
    
            template<class B, class C>
            inline void evalComponentSubFrom(ScT & dest, const std::pair<B, C> & data, unsigned i) const
            {
                data.first.evalComponentSubFrom(dest, i);
                data.second.evalComponentSubFrom(dest, i);
            }
    
            template<class V, class B, class C>
            inline bool usesVec(const V & v, const std::pair<B, C> & data) const
            {
                return data.first.usesVec(v) || data.second.usesVec(v);
            }
    
            template<class V, class B, class C>
            inline bool usesVecNF(const V & v, const std::pair<B, C> & data) const
            {
                return data.first.usesVecNF(v) || data.second.usesVec(v);
            }
        };

        template<class ScT>
        class VecSubOp
        {
        public:
            template<class A, class B, class C>
            inline void assignTo(A & d, const std::pair<B, C> & data) const
            {
                unsigned s = size(data);
                for (unsigned i = 0; i < s; ++i)
                    d[i] = data.first.evalComponent(i) - data.second.evalComponent(i);
            }
    
            template<class B, class C>
            inline math_vector<ScT> evaluate(const std::pair<B, C> & data) const
            {
                math_vector<ScT> res(size());
                assignTo(res, data);
            }
    
            template<class B, class C>
            inline unsigned size(const std::pair<B, C> & data) const
            {
                return data.first.size();
            }
    
            template<class B, class C>
            inline ScT evalComponent(const std::pair<B, C> & data, unsigned i) const
            {
                return data.first.evalComponent(i) - data.second.evalComponent(i);
            }
    
            template<class B, class C>
            inline void evalComponentTo(ScT & dest, const std::pair<B, C> & data, unsigned i) const
            {
                data.first.evalComponentTo(dest, i); 
                data.second.evalComponentSubFrom(dest, i);
            }
    
            template<class B, class C>
            inline void evalComponentAddTo(ScT & dest, const std::pair<B, C> & data, unsigned i) const
            {
                data.first.evalComponentAddTo(dest, i); 
                data.second.evalComponentSubFrom(dest, i);
            }
    
            template<class B, class C>
            inline void evalComponentSubFrom(ScT & dest, const std::pair<B, C> & data, unsigned i) const
            {
                data.first.evalComponentSubFrom(dest, i);
                data.second.evalComponentAddTo(dest, i);
            }
    
            template<class V, class B, class C>
            inline bool usesVec(const V & v, const std::pair<B, C> & data) const
            {
                return data.first.usesVec(v) || data.second.usesVec(v);
            }
    
            template<class V, class B, class C>
            inline bool usesVecNF(const V & v, const std::pair<B, C> & data) const
            {
                return data.first.usesVecNF(v) || data.second.usesVec(v);
            }
        };

        template<class ScT>
        class VecNegOp
        {
        public:
            template<class A, class B>
            inline void assignTo(A & d, const B & data) const
            {
                unsigned s = size(data);
                for (unsigned i = 0; i < s; ++i)
                    d[i] = -data.evalComponent(i);
            }
    
            template<class B, class C>
            inline math_vector<ScT> evaluate(const std::pair<B, C> & data) const
            {
                math_vector<ScT> res(size());
                assignTo(res, data);
            }
    
            template<class B>
            inline unsigned size(const B & data) const
            {
                return data.size();
            }
    
            template<class B>
            inline ScT evalComponent(const B & data, unsigned i) const
            {
                return -data.evalComponent(i);
            }
    
            template<class B>
            inline void evalComponentTo(ScT & dest, const B & data, unsigned i) const
            {
                dest = -data.evalComponent(i); 
            }
    
            template<class B>
            inline void evalComponentAddTo(ScT & dest, const B & data, unsigned i) const
            {
                data.evalComponentSubFrom(dest, i); 
            }
    
            template<class B>
            inline void evalComponentSubFrom(ScT & dest, const B & data, unsigned i) const
            {
                data.evalComponentAddTo(dest, i);
            }
    
            template<class V, class B>
            inline bool usesVec(const V & v, const B & data) const
            {
                return data.usesVec(v);
            }
    
            template<class V, class B>
            inline bool usesVecNF(const V & v, const B & data) const
            {
                return false;
            }
        };

        template<class ScT>
        class VecSMulOp
        {
        private:
            const ScT & d_val;
    
        public:
            VecSMulOp(const ScT & val)
                : d_val(val)
            {
            }
    
            template<class A, class B>
            inline void assignTo(A & d, const B & data) const
            {
                unsigned s = size(data);
                for (unsigned i = 0; i < s; ++i)
                    d[i] = data.evalComponent(i) * d_val;
            }
    
            template<class B>
            inline math_vector<ScT> evaluate(const B & data) const
            {
                math_vector<ScT> res(size());
                assignTo(res, data);
            }
    
            template<class B>
            inline unsigned size(const B & data) const
            {
                return data.size();
            }
    
            template<class B>
            inline ScT evalComponent(const B & data, unsigned i) const
            {
                return data.evalComponent(i) * d_val;
            }
    
            template<class B>
            inline void evalComponentTo(ScT & dest, const B & data, unsigned i) const
            {
                dest = data.evalComponent(i) * d_val;
            }
    
            template<class B>
            inline void evalComponentAddTo(ScT & dest, const B & data, unsigned i) const
            {
                dest += data.evalComponent(i) * d_val;
            }
    
            template<class B>
            inline void evalComponentSubFrom(ScT & dest, const B & data, unsigned i) const
            {
                dest -= data.evalComponent(i) * d_val;
            }
    
            template<class V, class B>
            inline bool usesVec(const V & v, const B & data) const
            {
                return data.usesVec(v);
            }
    
            template<class V, class B>
            inline bool usesVecNF(const V & v, const B & data) const
            {
                return false;
            }
        };

        template<class ScT>
        class VecSDivOp
        {
        private:
            const ScT & d_val;
    
        public:
            VecSDivOp(const ScT & val)
                : d_val(val)
            {
            }
    
            template<class A, class B>
            inline void assignTo(A & d, const B & data) const
            {
                unsigned s = size(data);
                for (unsigned i = 0; i < s; ++i)
                    d[i] = data.evalComponent(i) / d_val;
            }
    
            template<class B>
            inline math_vector<ScT> evaluate(const B & data) const
            {
                math_vector<ScT> res(size());
                assignTo(res, data);
            }
    
            template<class B>
            inline unsigned size(const B & data) const
            {
                return data.size();
            }
    
            template<class B>
            inline ScT evalComponent(const B & data, unsigned i) const
            {
                return data.evalComponent(i) / d_val;
            }
    
            template<class B>
            inline void evalComponentTo(ScT & dest, const B & data, unsigned i) const
            {
                dest = data.evalComponent(i) / d_val;
            }
    
            template<class B>
            inline void evalComponentAddTo(ScT & dest, const B & data, unsigned i) const
            {
                dest += data.evalComponent(i) / d_val;
            }
    
            template<class B>
            inline void evalComponentSubFrom(ScT & dest, const B & data, unsigned i) const
            {
                dest -= data.evalComponent(i) / d_val;
            }
    
            template<class V, class B>
            inline bool usesVec(const V & v, const B & data) const
            {
                return data.usesVec(v);
            }
    
            template<class V, class B>
            inline bool usesVecNF(const V & v, const B & data) const
            {
                return false;
            }
        };

        template<class ScT>
        class VecSModOp
        {
        private:
            const ScT & d_val;
    
        public:
            VecSModOp(const ScT & val)
                : d_val(val)
            {
            }
    
            template<class A, class B>
            inline void assignTo(A & d, const B & data) const
            {
                unsigned s = size(data);
                for (unsigned i = 0; i < s; ++i)
                    d[i] = data.evalComponent(i) % d_val;
            }
    
            template<class B>
            inline math_vector<ScT> evaluate(const B & data) const
            {
                math_vector<ScT> res(size());
                assignTo(res, data);
            }
    
            template<class B>
            inline unsigned size(const B & data) const
            {
                return data.size();
            }
    
            template<class B>
            inline ScT evalComponent(const B & data, unsigned i) const
            {
                return data.evalComponent(i) % d_val;
            }
    
            template<class B>
            inline void evalComponentTo(ScT & dest, const B & data, unsigned i) const
            {
                dest = data.evalComponent(i) % d_val;
            }
    
            template<class B>
            inline void evalComponentAddTo(ScT & dest, const B & data, unsigned i) const
            {
                dest += data.evalComponent(i) % d_val;
            }
    
            template<class B>
            inline void evalComponentSubFrom(ScT & dest, const B & data, unsigned i) const
            {
                dest -= data.evalComponent(i) % d_val;
            }
    
            template<class V, class B>
            inline bool usesVec(const V & v, const B & data) const
            {
                return data.usesVec(v);
            }
    
            template<class V, class B>
            inline bool usesVecNF(const V & v, const B & data) const
            {
                return false;
            }
        };

        template<class ScT, class ST, class V>
        inline VecExpression<VecAddOp<ScT>, std::pair<VecWrapper<math_vector<ScT, ST>, ScT>, VecWrapper<V, ScT> >, ScT> operator + (const math_vector<ScT, ST> & v, const V & w)
        { return VecExpression<VecAddOp<ScT>, std::pair<VecWrapper<math_vector<ScT, ST>, ScT>, VecWrapper<V, ScT> >, ScT>(VecAddOp<ScT>(), std::make_pair(VecWrapper<math_vector<ScT, ST>, ScT>(v), VecWrapper<V, ScT>(w))); }
        template<class ScT, class ST, class V>
        inline VecExpression<VecSubOp<ScT>, std::pair<VecWrapper<math_vector<ScT, ST>, ScT>, VecWrapper<V, ScT> >, ScT> operator - (const math_vector<ScT, ST> & v, const V & w)
        { return VecExpression<VecSubOp<ScT>, std::pair<VecWrapper<math_vector<ScT, ST>, ScT>, VecWrapper<V, ScT> >, ScT>(VecSubOp<ScT>(), std::make_pair(VecWrapper<math_vector<ScT, ST>, ScT>(v), VecWrapper<V, ScT>(w))); }
        template<class ScT, class ST, class V>
        inline VecExpression<VecAddOp<ScT>, std::pair<VecWrapper<MatrixRow<ScT, ST>, ScT>, VecWrapper<V, ScT> >, ScT> operator + (const MatrixRow<ScT, ST> & v, const V & w)
        { return VecExpression<VecAddOp<ScT>, std::pair<VecWrapper<MatrixRow<ScT, ST>, ScT>, VecWrapper<V, ScT> >, ScT>(VecAddOp<ScT>(), std::make_pair(VecWrapper<MatrixRow<ScT, ST>, ScT>(v), VecWrapper<V, ScT>(w))); }
        template<class ScT, class ST, class V>
        inline VecExpression<VecSubOp<ScT>, std::pair<VecWrapper<MatrixRow<ScT, ST>, ScT>, VecWrapper<V, ScT> >, ScT> operator - (const MatrixRow<ScT, ST> & v, const V & w)
        { return VecExpression<VecSubOp<ScT>, std::pair<VecWrapper<MatrixRow<ScT, ST>, ScT>, VecWrapper<V, ScT> >, ScT>(VecSubOp<ScT>(), std::make_pair(VecWrapper<MatrixRow<ScT, ST>, ScT>(v), VecWrapper<V, ScT>(w))); }
        template<class ScT, class ST, class V>
        inline VecExpression<VecAddOp<ScT>, std::pair<VecWrapper<MatrixCRow<ScT, ST>, ScT>, VecWrapper<V, ScT> >, ScT> operator + (const MatrixCRow<ScT, ST> & v, const V & w)
        { return VecExpression<VecAddOp<ScT>, std::pair<VecWrapper<MatrixCRow<ScT, ST>, ScT>, VecWrapper<V, ScT> >, ScT>(VecAddOp<ScT>(), std::make_pair(VecWrapper<MatrixCRow<ScT, ST>, ScT>(v), VecWrapper<V, ScT>(w))); }
        template<class ScT, class ST, class V>
        inline VecExpression<VecSubOp<ScT>, std::pair<VecWrapper<MatrixCRow<ScT, ST>, ScT>, VecWrapper<V, ScT> >, ScT> operator - (const MatrixCRow<ScT, ST> & v, const V & w)
        { return VecExpression<VecSubOp<ScT>, std::pair<VecWrapper<MatrixCRow<ScT, ST>, ScT>, VecWrapper<V, ScT> >, ScT>(VecSubOp<ScT>(), std::make_pair(VecWrapper<MatrixCRow<ScT, ST>, ScT>(v), VecWrapper<V, ScT>(w))); }
        template<class ScT, class ST, class V>
        inline VecExpression<VecAddOp<ScT>, std::pair<VecWrapper<MatrixColumn<ScT, ST>, ScT>, VecWrapper<V, ScT> >, ScT> operator + (const MatrixColumn<ScT, ST> & v, const V & w)
        { return VecExpression<VecAddOp<ScT>, std::pair<VecWrapper<MatrixColumn<ScT, ST>, ScT>, VecWrapper<V, ScT> >, ScT>(VecAddOp<ScT>(), std::make_pair(VecWrapper<MatrixColumn<ScT, ST>, ScT>(v), VecWrapper<V, ScT>(w))); }
        template<class ScT, class ST, class V>
        inline VecExpression<VecSubOp<ScT>, std::pair<VecWrapper<MatrixColumn<ScT, ST>, ScT>, VecWrapper<V, ScT> >, ScT> operator - (const MatrixColumn<ScT, ST> & v, const V & w)
        { return VecExpression<VecSubOp<ScT>, std::pair<VecWrapper<MatrixColumn<ScT, ST>, ScT>, VecWrapper<V, ScT> >, ScT>(VecSubOp<ScT>(), std::make_pair(VecWrapper<MatrixColumn<ScT, ST>, ScT>(v), VecWrapper<V, ScT>(w))); }
        template<class ScT, class ST, class V>
        inline VecExpression<VecAddOp<ScT>, std::pair<VecWrapper<MatrixCColumn<ScT, ST>, ScT>, VecWrapper<V, ScT> >, ScT> operator + (const MatrixCColumn<ScT, ST> & v, const V & w)
        { return VecExpression<VecAddOp<ScT>, std::pair<VecWrapper<MatrixCColumn<ScT, ST>, ScT>, VecWrapper<V, ScT> >, ScT>(VecAddOp<ScT>(), std::make_pair(VecWrapper<MatrixCColumn<ScT, ST>, ScT>(v), VecWrapper<V, ScT>(w))); }
        template<class ScT, class ST, class V>
        inline VecExpression<VecSubOp<ScT>, std::pair<VecWrapper<MatrixCColumn<ScT, ST>, ScT>, VecWrapper<V, ScT> >, ScT> operator - (const MatrixCColumn<ScT, ST> & v, const V & w)
        { return VecExpression<VecSubOp<ScT>, std::pair<VecWrapper<MatrixCColumn<ScT, ST>, ScT>, VecWrapper<V, ScT> >, ScT>(VecSubOp<ScT>(), std::make_pair(VecWrapper<MatrixCColumn<ScT, ST>, ScT>(v), VecWrapper<V, ScT>(w))); }

        template<class ScT, class ST, class O2, class D2, class T2>
        inline VecExpression<VecAddOp<ScT>, std::pair<VecWrapper<math_vector<ScT, ST>, ScT>, VecExpression<O2, D2, T2> >, ScT> operator + (const math_vector<ScT, ST> & v, const VecExpression<O2, D2, T2> & w)
        { return VecExpression<VecAddOp<ScT>, std::pair<VecWrapper<math_vector<ScT, ST>, ScT>, VecExpression<O2, D2, T2> >, ScT>(VecAddOp<ScT>(), std::make_pair(VecWrapper<math_vector<ScT, ST>, ScT>(v), w)); }
        template<class ScT, class ST, class O2, class D2, class T2>
        inline VecExpression<VecSubOp<ScT>, std::pair<VecWrapper<math_vector<ScT, ST>, ScT>, VecExpression<O2, D2, T2> >, ScT> operator - (const math_vector<ScT, ST> & v, const VecExpression<O2, D2, T2> & w)
        { return VecExpression<VecSubOp<ScT>, std::pair<VecWrapper<math_vector<ScT, ST>, ScT>, VecExpression<O2, D2, T2> >, ScT>(VecSubOp<ScT>(), std::make_pair(VecWrapper<math_vector<ScT, ST>, ScT>(v), w)); }
        template<class ScT, class ST, class O2, class D2, class T2>
        inline VecExpression<VecAddOp<ScT>, std::pair<VecWrapper<MatrixRow<ScT, ST>, ScT>, VecExpression<O2, D2, T2> >, ScT> operator + (const MatrixRow<ScT, ST> & v, const VecExpression<O2, D2, T2> & w)
        { return VecExpression<VecAddOp<ScT>, std::pair<VecWrapper<MatrixRow<ScT, ST>, ScT>, VecExpression<O2, D2, T2> >, ScT>(VecAddOp<ScT>(), std::make_pair(VecWrapper<MatrixRow<ScT, ST>, ScT>(v), w)); }
        template<class ScT, class ST, class O2, class D2, class T2>
        inline VecExpression<VecSubOp<ScT>, std::pair<VecWrapper<MatrixRow<ScT, ST>, ScT>, VecExpression<O2, D2, T2> >, ScT> operator - (const MatrixRow<ScT, ST> & v, const VecExpression<O2, D2, T2> & w)
        { return VecExpression<VecSubOp<ScT>, std::pair<VecWrapper<MatrixRow<ScT, ST>, ScT>, VecExpression<O2, D2, T2> >, ScT>(VecSubOp<ScT>(), std::make_pair(VecWrapper<MatrixRow<ScT, ST>, ScT>(v), w)); }
        template<class ScT, class ST, class O2, class D2, class T2>
        inline VecExpression<VecAddOp<ScT>, std::pair<VecWrapper<MatrixCRow<ScT, ST>, ScT>, VecExpression<O2, D2, T2> >, ScT> operator + (const MatrixCRow<ScT, ST> & v, const VecExpression<O2, D2, T2> & w)
        { return VecExpression<VecAddOp<ScT>, std::pair<VecWrapper<MatrixCRow<ScT, ST>, ScT>, VecExpression<O2, D2, T2> >, ScT>(VecAddOp<ScT>(), std::make_pair(VecWrapper<MatrixCRow<ScT, ST>, ScT>(v), w)); }
        template<class ScT, class ST, class O2, class D2, class T2>
        inline VecExpression<VecSubOp<ScT>, std::pair<VecWrapper<MatrixCRow<ScT, ST>, ScT>, VecExpression<O2, D2, T2> >, ScT> operator - (const MatrixCRow<ScT, ST> & v, const VecExpression<O2, D2, T2> & w)
        { return VecExpression<VecSubOp<ScT>, std::pair<VecWrapper<MatrixCRow<ScT, ST>, ScT>, VecExpression<O2, D2, T2> >, ScT>(VecSubOp<ScT>(), std::make_pair(VecWrapper<MatrixCRow<ScT, ST>, ScT>(v), w)); }
        template<class ScT, class ST, class O2, class D2, class T2>
        inline VecExpression<VecAddOp<ScT>, std::pair<VecWrapper<MatrixColumn<ScT, ST>, ScT>, VecExpression<O2, D2, T2> >, ScT> operator + (const MatrixColumn<ScT, ST> & v, const VecExpression<O2, D2, T2> & w)
        { return VecExpression<VecAddOp<ScT>, std::pair<VecWrapper<MatrixColumn<ScT, ST>, ScT>, VecExpression<O2, D2, T2> >, ScT>(VecAddOp<ScT>(), std::make_pair(VecWrapper<MatrixColumn<ScT, ST>, ScT>(v), w)); }
        template<class ScT, class ST, class O2, class D2, class T2>
        inline VecExpression<VecSubOp<ScT>, std::pair<VecWrapper<MatrixColumn<ScT, ST>, ScT>, VecExpression<O2, D2, T2> >, ScT> operator - (const MatrixColumn<ScT, ST> & v, const VecExpression<O2, D2, T2> & w)
        { return VecExpression<VecSubOp<ScT>, std::pair<VecWrapper<MatrixColumn<ScT, ST>, ScT>, VecExpression<O2, D2, T2> >, ScT>(VecSubOp<ScT>(), std::make_pair(VecWrapper<MatrixColumn<ScT, ST>, ScT>(v), w)); }
        template<class ScT, class ST, class O2, class D2, class T2>
        inline VecExpression<VecAddOp<ScT>, std::pair<VecWrapper<MatrixCColumn<ScT, ST>, ScT>, VecExpression<O2, D2, T2> >, ScT> operator + (const MatrixCColumn<ScT, ST> & v, const VecExpression<O2, D2, T2> & w)
        { return VecExpression<VecAddOp<ScT>, std::pair<VecWrapper<MatrixCColumn<ScT, ST>, ScT>, VecExpression<O2, D2, T2> >, ScT>(VecAddOp<ScT>(), std::make_pair(VecWrapper<MatrixCColumn<ScT, ST>, ScT>(v), w)); }
        template<class ScT, class ST, class O2, class D2, class T2>
        inline VecExpression<VecSubOp<ScT>, std::pair<VecWrapper<MatrixCColumn<ScT, ST>, ScT>, VecExpression<O2, D2, T2> >, ScT> operator - (const MatrixCColumn<ScT, ST> & v, const VecExpression<O2, D2, T2> & w)
        { return VecExpression<VecSubOp<ScT>, std::pair<VecWrapper<MatrixCColumn<ScT, ST>, ScT>, VecExpression<O2, D2, T2> >, ScT>(VecSubOp<ScT>(), std::make_pair(VecWrapper<MatrixCColumn<ScT, ST>, ScT>(v), w)); }

        template<class T, class ST>
        inline VecExpression<VecNegOp<T>, VecWrapper<math_vector<T, ST>, T>, T> operator - (const math_vector<T, ST> & v)
        { return VecExpression<VecNegOp<T>, VecWrapper<math_vector<T, ST>, T>, T>(VecNegOp<T>(), VecWrapper<math_vector<T, ST>, T>(v)); }
        template<class ScT, class T, class ST>
        inline VecExpression<VecSMulOp<ScT>, VecWrapper<math_vector<T, ST>, T>, ScT> operator * (const ScT & s, const math_vector<T, ST> & v)
        { return VecExpression<VecSMulOp<ScT>, VecWrapper<math_vector<T, ST>, T>, ScT>(VecSMulOp<ScT>(s), VecWrapper<math_vector<T, ST>, T>(v)); }
        template<class ScT, class T, class ST>
        inline VecExpression<VecSMulOp<T>, VecWrapper<math_vector<T, ST>, T>, T> operator * (const math_vector<T, ST> & v, const ScT & s)
        { return VecExpression<VecSMulOp<T>, VecWrapper<math_vector<T, ST>, T>, T>(VecSMulOp<T>(s), VecWrapper<math_vector<T, ST>, T>(v)); }
        template<class ScT, class T, class ST>
        inline VecExpression<VecSDivOp<T>, VecWrapper<math_vector<T, ST>, T>, ScT> operator / (const math_vector<T, ST> & v, const ScT & s)
        { return VecExpression<VecSDivOp<T>, VecWrapper<math_vector<T, ST>, T>, T>(VecSDivOp<T>(s), VecWrapper<math_vector<T, ST>, T>(v)); }
        template<class ScT, class T, class ST>
        inline VecExpression<VecSModOp<T>, VecWrapper<math_vector<T, ST>, T>, ScT> operator % (const math_vector<T, ST> & v, const ScT & s)
        { return VecExpression<VecSModOp<T>, VecWrapper<math_vector<T, ST>, T>, T>(VecSModOp<T>(s), VecWrapper<math_vector<T, ST>, T>(v)); }
        template<class T, class ST>
        inline VecExpression<VecNegOp<T>, VecWrapper<MatrixRow<T, ST>, T>, T> operator - (const MatrixRow<T, ST> & v)
        { return VecExpression<VecNegOp<T>, VecWrapper<MatrixRow<T, ST>, T>, T>(VecNegOp<T>(), VecWrapper<MatrixRow<T, ST>, T>(v)); }
        template<class ScT, class T, class ST>
        inline VecExpression<VecSMulOp<ScT>, VecWrapper<MatrixRow<T, ST>, T>, ScT> operator * (const ScT & s, const MatrixRow<T, ST> & v)
        { return VecExpression<VecSMulOp<ScT>, VecWrapper<MatrixRow<T, ST>, T>, ScT>(VecSMulOp<ScT>(s), VecWrapper<MatrixRow<T, ST>, T>(v)); }
        template<class ScT, class T, class ST>
        inline VecExpression<VecSMulOp<T>, VecWrapper<MatrixRow<T, ST>, T>, T> operator * (const MatrixRow<T, ST> & v, const ScT & s)
        { return VecExpression<VecSMulOp<T>, VecWrapper<MatrixRow<T, ST>, T>, T>(VecSMulOp<T>(s), VecWrapper<MatrixRow<T, ST>, T>(v)); }
        template<class ScT, class T, class ST>
        inline VecExpression<VecSDivOp<T>, VecWrapper<MatrixRow<T, ST>, T>, ScT> operator / (const MatrixRow<T, ST> & v, const ScT & s)
        { return VecExpression<VecSDivOp<T>, VecWrapper<MatrixRow<T, ST>, T>, T>(VecSDivOp<T>(s), VecWrapper<MatrixRow<T, ST>, T>(v)); }
        template<class ScT, class T, class ST>
        inline VecExpression<VecSModOp<T>, VecWrapper<MatrixRow<T, ST>, T>, ScT> operator % (const MatrixRow<T, ST> & v, const ScT & s)
        { return VecExpression<VecSModOp<T>, VecWrapper<MatrixRow<T, ST>, T>, T>(VecSModOp<T>(s), VecWrapper<MatrixRow<T, ST>, T>(v)); }
        template<class T, class ST>
        inline VecExpression<VecNegOp<T>, VecWrapper<MatrixCRow<T, ST>, T>, T> operator - (const MatrixCRow<T, ST> & v)
        { return VecExpression<VecNegOp<T>, VecWrapper<MatrixCRow<T, ST>, T>, T>(VecNegOp<T>(), VecWrapper<MatrixCRow<T, ST>, T>(v)); }
        template<class ScT, class T, class ST>
        inline VecExpression<VecSMulOp<ScT>, VecWrapper<MatrixCRow<T, ST>, T>, ScT> operator * (const ScT & s, const MatrixCRow<T, ST> & v)
        { return VecExpression<VecSMulOp<ScT>, VecWrapper<MatrixCRow<T, ST>, T>, ScT>(VecSMulOp<ScT>(s), VecWrapper<MatrixCRow<T, ST>, T>(v)); }
        template<class ScT, class T, class ST>
        inline VecExpression<VecSMulOp<T>, VecWrapper<MatrixCRow<T, ST>, T>, T> operator * (const MatrixCRow<T, ST> & v, const ScT & s)
        { return VecExpression<VecSMulOp<T>, VecWrapper<MatrixCRow<T, ST>, T>, T>(VecSMulOp<T>(s), VecWrapper<MatrixCRow<T, ST>, T>(v)); }
        template<class ScT, class T, class ST>
        inline VecExpression<VecSDivOp<T>, VecWrapper<MatrixCRow<T, ST>, T>, ScT> operator / (const MatrixCRow<T, ST> & v, const ScT & s)
        { return VecExpression<VecSDivOp<T>, VecWrapper<MatrixCRow<T, ST>, T>, T>(VecSDivOp<T>(s), VecWrapper<MatrixCRow<T, ST>, T>(v)); }
        template<class ScT, class T, class ST>
        inline VecExpression<VecSModOp<T>, VecWrapper<MatrixCRow<T, ST>, T>, ScT> operator % (const MatrixCRow<T, ST> & v, const ScT & s)
        { return VecExpression<VecSModOp<T>, VecWrapper<MatrixCRow<T, ST>, T>, T>(VecSModOp<T>(s), VecWrapper<MatrixCRow<T, ST>, T>(v)); }
        template<class T, class ST>
        inline VecExpression<VecNegOp<T>, VecWrapper<MatrixColumn<T, ST>, T>, T> operator - (const MatrixColumn<T, ST> & v)
        { return VecExpression<VecNegOp<T>, VecWrapper<MatrixColumn<T, ST>, T>, T>(VecNegOp<T>(), VecWrapper<MatrixColumn<T, ST>, T>(v)); }
        template<class ScT, class T, class ST>
        inline VecExpression<VecSMulOp<ScT>, VecWrapper<MatrixColumn<T, ST>, T>, ScT> operator * (const ScT & s, const MatrixColumn<T, ST> & v)
        { return VecExpression<VecSMulOp<ScT>, VecWrapper<MatrixColumn<T, ST>, T>, ScT>(VecSMulOp<ScT>(s), VecWrapper<MatrixColumn<T, ST>, T>(v)); }
        template<class ScT, class T, class ST>
        inline VecExpression<VecSMulOp<T>, VecWrapper<MatrixColumn<T, ST>, T>, T> operator * (const MatrixColumn<T, ST> & v, const ScT & s)
        { return VecExpression<VecSMulOp<T>, VecWrapper<MatrixColumn<T, ST>, T>, T>(VecSMulOp<T>(s), VecWrapper<MatrixColumn<T, ST>, T>(v)); }
        template<class ScT, class T, class ST>
        inline VecExpression<VecSDivOp<T>, VecWrapper<MatrixColumn<T, ST>, T>, ScT> operator / (const MatrixColumn<T, ST> & v, const ScT & s)
        { return VecExpression<VecSDivOp<T>, VecWrapper<MatrixColumn<T, ST>, T>, T>(VecSDivOp<T>(s), VecWrapper<MatrixColumn<T, ST>, T>(v)); }
        template<class ScT, class T, class ST>
        inline VecExpression<VecSModOp<T>, VecWrapper<MatrixColumn<T, ST>, T>, ScT> operator % (const MatrixColumn<T, ST> & v, const ScT & s)
        { return VecExpression<VecSModOp<T>, VecWrapper<MatrixColumn<T, ST>, T>, T>(VecSModOp<T>(s), VecWrapper<MatrixColumn<T, ST>, T>(v)); }
        template<class T, class ST>
        inline VecExpression<VecNegOp<T>, VecWrapper<MatrixCColumn<T, ST>, T>, T> operator - (const MatrixCColumn<T, ST> & v)
        { return VecExpression<VecNegOp<T>, VecWrapper<MatrixCColumn<T, ST>, T>, T>(VecNegOp<T>(), VecWrapper<MatrixCColumn<T, ST>, T>(v)); }
        template<class ScT, class T, class ST>
        inline VecExpression<VecSMulOp<ScT>, VecWrapper<MatrixCColumn<T, ST>, T>, ScT> operator * (const ScT & s, const MatrixCColumn<T, ST> & v)
        { return VecExpression<VecSMulOp<ScT>, VecWrapper<MatrixCColumn<T, ST>, T>, ScT>(VecSMulOp<ScT>(s), VecWrapper<MatrixCColumn<T, ST>, T>(v)); }
        template<class ScT, class T, class ST>
        inline VecExpression<VecSMulOp<T>, VecWrapper<MatrixCColumn<T, ST>, T>, T> operator * (const MatrixCColumn<T, ST> & v, const ScT & s)
        { return VecExpression<VecSMulOp<T>, VecWrapper<MatrixCColumn<T, ST>, T>, T>(VecSMulOp<T>(s), VecWrapper<MatrixCColumn<T, ST>, T>(v)); }
        template<class ScT, class T, class ST>
        inline VecExpression<VecSDivOp<T>, VecWrapper<MatrixCColumn<T, ST>, T>, ScT> operator / (const MatrixCColumn<T, ST> & v, const ScT & s)
        { return VecExpression<VecSDivOp<T>, VecWrapper<MatrixCColumn<T, ST>, T>, T>(VecSDivOp<T>(s), VecWrapper<MatrixCColumn<T, ST>, T>(v)); }
        template<class ScT, class T, class ST>
        inline VecExpression<VecSModOp<T>, VecWrapper<MatrixCColumn<T, ST>, T>, ScT> operator % (const MatrixCColumn<T, ST> & v, const ScT & s)
        { return VecExpression<VecSModOp<T>, VecWrapper<MatrixCColumn<T, ST>, T>, T>(VecSModOp<T>(s), VecWrapper<MatrixCColumn<T, ST>, T>(v)); }

        template<class O1, class D1, class T1, class V>
        inline VecExpression<VecAddOp<T1>, std::pair<VecExpression<O1, D1, T1>, VecWrapper<V, T1> >, T1> operator + (const VecExpression<O1, D1, T1> & v, const V & w)
        { return VecExpression<VecAddOp<T1>, std::pair<VecExpression<O1, D1, T1>, VecWrapper<V, T1> >, T1>(VecAddOp<T1>(), std::make_pair(v, VecWrapper<V, T1>(w))); }
        template<class O1, class D1, class T1, class V>
        inline VecExpression<VecSubOp<T1>, std::pair<VecExpression<O1, D1, T1>, VecWrapper<V, T1> >, T1> operator - (const VecExpression<O1, D1, T1> & v, const V & w)
        { return VecExpression<VecSubOp<T1>, std::pair<VecExpression<O1, D1, T1>, VecWrapper<V, T1> >, T1>(VecSubOp<T1>(), std::make_pair(v, VecWrapper<V, T1>(w))); }

        template<class O1, class D1, class T1, class O2, class D2, class T2>
        inline VecExpression<VecAddOp<T1>, std::pair<VecExpression<O1, D1, T1>, VecExpression<O2, D2, T2> >, T1> operator + (const VecExpression<O1, D1, T1> & v, const VecExpression<O2, D2, T2> & w)
        { return VecExpression<VecAddOp<T1>, std::pair<VecExpression<O1, D1, T1>, VecExpression<O2, D2, T2> >, T1>(VecAddOp<T1>(), std::make_pair(v, w)); }
        template<class O1, class D1, class T1, class O2, class D2, class T2>
        inline VecExpression<VecSubOp<T1>, std::pair<VecExpression<O1, D1, T1>, VecExpression<O2, D2, T2> >, T1> operator - (const VecExpression<O1, D1, T1> & v, const VecExpression<O2, D2, T2> & w)
        { return VecExpression<VecSubOp<T1>, std::pair<VecExpression<O1, D1, T1>, VecExpression<O2, D2, T2> >, T1>(VecSubOp<T1>(), std::make_pair(v, w)); }

        template<class O1, class D1, class T1>
        inline VecExpression<VecNegOp<T1>, VecExpression<O1, D1, T1>, T1> operator - (const VecExpression<O1, D1, T1> & v)
        { return VecExpression<VecNegOp<T1>, VecExpression<O1, D1, T1>, T1>(VecNegOp<T1>(), v); }
        template<class O1, class D1, class T1, class ScT>
        inline VecExpression<VecSMulOp<ScT>, VecExpression<O1, D1, T1>, ScT> operator * (const ScT & s, const VecExpression<O1, D1, T1> & v)
        { return VecExpression<VecSMulOp<ScT>, VecExpression<O1, D1, T1>, ScT>(VecSMulOp<ScT>(s), v); }
        template<class O1, class D1, class T1, class ScT>
        inline VecExpression<VecSMulOp<T1>, VecExpression<O1, D1, T1>, T1> operator * (const VecExpression<O1, D1, T1> & v, const ScT & s)
        { return VecExpression<VecSMulOp<T1>, VecExpression<O1, D1, T1>, T1>(VecSMulOp<T1>(s), v); }
        template<class O1, class D1, class T1, class ScT>
        inline VecExpression<VecSDivOp<T1>, VecExpression<O1, D1, T1>, T1> operator / (const VecExpression<O1, D1, T1> & v, const ScT & s)
        { return VecExpression<VecSDivOp<T1>, VecExpression<O1, D1, T1>, T1>(VecSDivOp<T1>(s), v); }
        template<class O1, class D1, class T1, class ScT>
        inline VecExpression<VecSModOp<T1>, VecExpression<O1, D1, T1>, T1> operator % (const VecExpression<O1, D1, T1> & v, const ScT & s)
        { return VecExpression<VecSModOp<T1>, VecExpression<O1, D1, T1>, T1>(VecSModOp<T1>(s), v); }
    }
}

#endif
