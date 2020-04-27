
#ifndef _SYMBA_FIELDS_VECTOR_HPP_
#define _SYMBA_FIELDS_VECTOR_HPP_

#include<cmath>
#include<algorithm>
#include<array>
#include<vector>
#include<numeric>
#include<sstream>
#include<iterator>
#include <initializer_list>

#include "util.h"

namespace symba {
namespace fields {

template<class Field, std::size_t dimension1, std::size_t dimension2> class matrix;

template<class Field, std::size_t dimension> class vector {
    private:
        typedef vector<Field, dimension> vector_type;
        typedef std::array<Field, dimension> storage_type;
        typedef typename std::array<Field, dimension>::iterator storage_type_iterator;
        storage_type values;
    public:
        vector(const Field &_value) : values() { values.fill(_value); }
        vector() : vector(0) { }
        vector(const storage_type &_values) : values(_values) { }
        template<class Iterator>vector(Iterator _begin, Iterator _end) : values() { storage_type_iterator _cur=values.begin(); for(Iterator _v=_begin; (_v!=_end) && (_cur!=values.end()); ++_v, ++_cur) *_cur = *_v; }
        vector(std::initializer_list<Field> const &_values) : vector(_values.begin(), _values.end()) { }
        const storage_type &get_values() const { return values; }
        matrix<Field, 1, dimension> get_row_matrix() const { return matrix<Field, 1, dimension>(values); }
        matrix<Field, dimension, 1> get_column_matrix() const { return matrix<Field, dimension, 1>(values); }
        const Field& operator[](std::size_t i) const { return values[i]; }
        Field& operator[](std::size_t i) { return values[i]; }
        vector_type operator*=(const Field &c) {
            if (c == 0) {
                std::fill(values.begin(), values.end(), 0);
            } else {
                for (auto value=values.begin(); value != values.end(); ++value) (*value) *= c;
            }
            return *this;
        }
        vector_type operator/=(const Field &c) {
            for (auto value=values.begin(); value != values.end(); ++value) (*value) /= c;
            return *this;
        }
        vector_type operator+=(const vector_type &b) {
            for (std::size_t i=0; i<values.size(); ++i) values[i] += b[i];
            return *this;
        }
        vector_type operator-=(const vector_type &b) {
            for (std::size_t i=0; i<values.size(); ++i) values[i] -= b[i];
            return *this;
        }
        Field dot(const vector_type& b) const {
            Field res(0);
            for (std::size_t i=0; i<values.size(); ++i) res += values[i] * b[i];
            return res;
        }
        Field inner_product(const vector_type& b) const {
            return dot(b);
        }
        Field l2_norm_sq() const {
            return dot(*this);
        }
        friend vector_type operator+(vector_type a, const vector_type &b) {
            a += b;
            return a;
        }
        friend vector_type operator-(vector_type a, const vector_type &b) {
            a -= b;
            return a;
        }
        friend vector_type operator+(vector_type a, const Field &c) {
            a += c;
            return a;
        }
        friend vector_type operator-(vector_type a, const Field &b) {
            a -= b;
            return a;
        }
        friend vector_type operator*(const Field &c, vector_type a) {
            a *= c;
            return a;
        }
        friend vector_type operator*(vector_type a, const Field &c) {
            a *= c;
            return a;
        }
        friend vector_type operator/(const Field &b, vector_type a) {
            a /= b;
            return a;
        }
        friend vector_type operator/(vector_type a, const Field &b) {
            a /= b;
            return a;
        }
        friend Field dot(vector_type a, const vector_type &b) {
            a.dot(b);
            return a;
        }
        string to_string() const {
            ostringstream stream;
            stream << "(";
            for (auto value=values.begin(); value != values.end(); ++value) {
                if (value != values.begin()) stream << ", ";
                stream << *value;
            }
            stream << ")";
            return stream.str();
        }
        friend string to_string(vector_type a) {
            return a.to_string();
        }
        friend ostream& operator<<(ostream& os, const vector_type& a) {
            os << a.to_string();
            return os;
        }
        friend inline int cmp(const vector_type& a, const vector_type& b) {
            return cmp(a.l2_norm()-b.l2_norm(),0);
        }
        friend inline bool operator==(const vector_type& a, const vector_type& b){ return cmp(a,b) == 0; }
        friend inline bool operator!=(const vector_type& a, const vector_type& b){ return cmp(a,b) != 0; }
        friend inline bool operator< (const vector_type& a, const vector_type& b){ return cmp(a,b) <  0; }
        friend inline bool operator> (const vector_type& a, const vector_type& b){ return cmp(a,b) >  0; }
        friend inline bool operator<=(const vector_type& a, const vector_type& b){ return cmp(a,b) <= 0; }
        friend inline bool operator>=(const vector_type& a, const vector_type& b){ return cmp(a,b) >= 0; }        
        friend inline bool operator==(const vector_type& a, const Field& b){ return cmp(a,b) == 0; }
        friend inline bool operator!=(const vector_type& a, const Field& b){ return cmp(a,b) != 0; }
        friend inline bool operator< (const vector_type& a, const Field& b){ return cmp(a,b) <  0; }
        friend inline bool operator> (const vector_type& a, const Field& b){ return cmp(a,b) >  0; }
        friend inline bool operator<=(const vector_type& a, const Field& b){ return cmp(a,b) <= 0; }
        friend inline bool operator>=(const vector_type& a, const Field& b){ return cmp(a,b) >= 0; }
        static inline vector_type zeros() {
            return vector_type();
        }
        static inline vector_type ones() {
            return vector_type(Field(1));
        }
        static inline vector_type e(std::size_t i) {
            vector_type v(0);
            v[i-1] = Field(1);
            return v;
        }
        static inline vector_type e1() {
            return e(1);
        }
        static inline vector_type e2() {
            return e(2);
        }
        static inline vector_type e3() {
            return e(3);
        }
        inline pair<vector_type, std::vector<vector_type> > factor(const std::vector<vector_type> &values, util::side_type side) {
            std::vector<vector_type> ret(values.size());
            vector_type common_factor = values[0];
            for(std::size_t i=0; i < dimension; i++) {
                std::vector<Field> tmp_input(values.size());
                for(size_t j=0; j < values.size(); j++) tmp_input[j] = values[j][i];
                pair<Field, std::vector<Field> > tmp = Field::factor(tmp_input);
                common_factor[i] = get<Field>(tmp);
                for(size_t j=0; j < ret.size(); j++) ret[j][i] = get<std::vector<Field> >(tmp)[j];
            }
            return pair<vector_type, std::vector<vector_type> >(common_factor, ret);
        }
        friend inline vector_type abs(const vector_type &_v) {
            return _v.l2_norm();
        }
        friend inline vector_type exp(const vector_type &_v) {
            vector_type v = _v;
            for(auto value=v.values.begin(); value != v.values.end(); ++value) *value = Field::exp(*value);
            return v;
        }
        friend inline vector_type log(const vector_type &_v) {
            vector_type v = _v;
            for(auto value=v.values.begin(); value != v.values.end(); ++value) *value = Field::log(*value);
            return v;
        }
}; // class vector

template<std::size_t dimension> class VectorIntField {
    public:
        using value_type=vector<int, dimension>;
        using power_type=unsigned int;
        using coefficient_type=vector<int, dimension>;
        struct power_type_max {
            enum { value = 10000 };
        };
        struct power_type_fill {
            enum { value = 8 };
        };
        struct is_value_addition_commutative {
            enum { value = 1 };
        };
        struct is_coefficient_multiplication_commutative {
            enum { value = 0 };
        };
        struct is_value_multiplication_commutative {
            enum { value = 0 };
        };
        inline static pair<value_type, std::vector<value_type> > factor(const std::vector<value_type> &values, const util::side_type &side=util::side_type::both) {
            return factor(values, side);
        }
        inline static value_type gcd(const value_type& val1, const value_type& val2) {
            return value_type(1);
        }
        inline static value_type abs(const value_type& val) {
            return abs(val);
        }
        inline static value_type exp(const value_type& val) {
            return exp(val);
        }
        inline static value_type log(const value_type& val) {
            return log(val);
        }
};

template<class Field, class ForwardIterator, std::size_t dimension> class vector_view {
    private:
        typedef fields::vector<Field, dimension> vector_type;
        typedef vector_view<Field, ForwardIterator, dimension> vector_view_type;
        typedef ForwardIterator iterator_type;
        typedef util::strided_iterator<ForwardIterator> vector_view_iterator;
        iterator_type _begin;
        std::size_t stride;
    public:
        vector_view(const iterator_type& __begin, std::size_t _stride=1) : _begin(__begin), stride(_stride) { }
        Field& operator[](std::size_t i) { return *(_begin+i*stride); }
        const Field& operator[](std::size_t i) const { return *(_begin+i*stride); }
        vector_view_iterator begin() const { return vector_view_iterator(_begin, stride); }
        vector_view_iterator end() const { return vector_view_iterator(_begin+dimension*stride, stride); }
        vector_type get_vector() { return vector_type(begin(), end()); }
        std::size_t get_stride() const { return stride; }
        void assign(const vector_type &other) {
            for(std::size_t i=0; i < dimension; i++) {
                *(_begin + i*stride) = other[i];
            }
        }
        template<class InputIterator>void assign(InputIterator __begin) {
            for(std::size_t i=0; i < dimension; i++) {
                *(_begin + i*stride) = *(__begin+i);
            }
        }
        Field dot(const vector_type &other) const {
            Field ret(0);
            for(std::size_t i=0; i < dimension; i++) {
                ret += *(_begin + i*stride) * other[i];
            }
            return ret;
        }
        Field l2_norm_sq() const {
            Field ret(0);
            for(std::size_t i=0; i < dimension; i++) {
                ret += (*(_begin + i*stride)) * (*(_begin + i*stride));
            }
            return ret;
        }
        string to_string() const {
            ostringstream stream;
            stream << "(";
            for (std::size_t i=0; i < dimension; i++) {
                if (i > 0) stream << ", ";
                stream << *(_begin+i*stride);
            }
            stream << ")";
            return stream.str();
        }
        friend string to_string(vector_view_type a) {
            return a.to_string();
        }
        friend ostream& operator<<(ostream& os, const vector_view_type& a) {
            os << a.to_string();
            return os;
        }
        inline bool operator==(const vector_view_type& rhs) const {return (_begin==rhs._begin) && (stride==rhs.get_stride()); }
        inline bool operator!=(const vector_view_type& rhs) const {return !(*this==rhs); }
};

}; // namespace fields
}; // namespace symba
#endif // _SYMBA_FIELDS_VECTOR_HPP_
