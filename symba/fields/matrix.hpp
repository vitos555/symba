
#ifndef _SYMBA_FIELDS_MATRIX_HPP_
#define _SYMBA_FIELDS_MATRIX_HPP_

#include<cmath>
#include<algorithm>
#include<array>
#include<vector>
#include<numeric>
#include<sstream>
#include<iterator>
#include <initializer_list>

#include "util.hpp"

namespace symba {
namespace fields {

enum orientation {
    type_defined=0,
    rows=1,
    columns=2
};

template<class Field, std::size_t dimension> class vector;
template<class Field, class ForwardIterator, std::size_t dimension> class vector_view; 

template<class Field, std::size_t dimension1, std::size_t dimension2> class _matrix {
    protected:
        typedef _matrix<Field, dimension2, dimension1> transposed_matrix_type;
        typedef _matrix<Field, dimension1, dimension1> dimension1_square_matrix_type;
        typedef _matrix<Field, dimension2, dimension2> dimension2_square_matrix_type;
        typedef _matrix<Field, dimension1, dimension2> matrix_type;
        typedef array<Field, dimension1*dimension2> storage_type;
        typedef typename array<Field, dimension1*dimension2>::iterator storage_type_iterator;
        typedef typename array<Field, dimension1*dimension2>::const_iterator storage_type_const_iterator;
        typedef vector<Field, dimension1> column_vector;
        typedef vector_view<Field, storage_type_iterator, dimension1> column_vector_view;
        typedef vector_view<Field, storage_type_const_iterator, dimension1> column_const_vector_view;
        typedef vector<Field, dimension2> row_vector;
        typedef vector_view<Field, storage_type_iterator, dimension2> row_vector_view;
        typedef vector_view<Field, storage_type_const_iterator, dimension2> row_const_vector_view;
        storage_type values;
    public:
        _matrix(const Field &_value) : values() { values.fill(_value); }
        _matrix() : _matrix(0) { }
        _matrix(const storage_type &_values) : values(_values) { }
        template<class InputIterator>_matrix(InputIterator _begin, InputIterator _end) : values() { storage_type_iterator _cur=values.begin(); for(InputIterator _v=_begin; (_v!=_end) && (_cur!=values.end()); ++_v, ++_cur ) *_cur = *_v; }
        _matrix(std::initializer_list<Field> const &_values) : _matrix(_values.begin(), _values.end()) { }
        template<class _vector>_matrix(const std::vector<_vector> &_vectors, orientation _orientation=orientation::type_defined) : values(0) {
            if constexpr(is_same_v<row_vector, column_vector> && is_same_v<_vector, row_vector>) {
                if (_orientation==orientation::columns)
                    for(std::size_t i=0; i < dimension2; i++) get_column(i).assign(_vectors[i]);
                else
                    for(std::size_t i=0; i < dimension1; i++) get_row(i).assign(_vectors[i]);
            } else if constexpr(is_same_v<_vector, row_vector>) {
                if ((_orientation==orientation::type_defined)||(_orientation==orientation::rows))
                    for(std::size_t i=0; i < dimension1; i++) get_row(i).assign(_vectors[i]);
            } else if constexpr(is_same_v<_vector, column_vector>) {
                if ((_orientation==orientation::type_defined)||(_orientation==orientation::columns))
                    for(std::size_t i=0; i < dimension2; i++) get_column(i).assign(_vectors[i]);
            }
        }
        row_vector_view operator[](std::size_t i) { return row_vector_view(values.begin() + i*dimension2, 1); }
        row_vector_view get_row(std::size_t i) { return row_vector_view(values.begin() + i*dimension2, 1); }
        row_const_vector_view const_row(std::size_t i) const { return row_const_vector_view(values.begin() + i*dimension2, 1); }
        column_vector_view get_column(std::size_t i) { return column_vector_view(values.begin()+i, dimension2); }
        column_const_vector_view const_column(std::size_t i) const { return column_const_vector_view(values.begin()+i, dimension2); }
        const Field& at(std::size_t i, std::size_t j) const { return values[i*dimension2 + j]; }
        const Field& at(std::size_t i) const { return values[i]; }
        Field& at(std::size_t i, std::size_t j) { return values[i*dimension2 + j]; }
        Field& at(std::size_t i) { return values[i]; }
        matrix_type operator*=(const Field &c) {
            for(auto _v=values.begin(); _v!=values.end(); ++_v) *_v *= c;
            return *this;
        }
        matrix_type operator/=(const Field &c) {
            for(auto _v=values.begin(); _v!=values.end(); ++_v) *_v /= c;
            return *this;
        }
        matrix_type operator+=(const matrix_type &b) {
            for(std::size_t i=0; i<dimension1+dimension2; i++) values[i] += b.at(i);
            return *this;
        }
        matrix_type operator-=(const matrix_type &b) {
            for(std::size_t i=0; i<dimension1+dimension2; i++) values[i] -= b.at(i);
            return *this;
        }
        friend column_vector operator*(const matrix_type &a, const row_vector &_row) {
            column_vector ret = row_vector(0);
            for(std::size_t i=0; i<dimension1; i++) ret[i] = a.const_row(i).dot(_row);
            return ret;
        }
        friend row_vector operator*(const column_vector &_column, const matrix_type &a) {
            column_vector ret = column_vector(0);
            for(std::size_t i=0; i<dimension2; i++) ret[i] = a.const_column(i).dot(_column);
            return ret;
        }
        friend matrix_type operator+(matrix_type a, const matrix_type &b) {
            a += b;
            return a;
        }
        friend matrix_type operator-(matrix_type a, const matrix_type &b) {
            a -= b;
            return a;
        }
        friend matrix_type operator*(matrix_type a, const Field &c) {
            a *= c;
            return a;
        }
        friend matrix_type operator/(matrix_type a, const Field &b) {
            a /= b;
            return a;
        }
        friend matrix_type operator*(const Field &c, matrix_type a) {
            a *= c;
            return a;
        }
        friend matrix_type operator/(const Field &b, matrix_type a) {
            a /= b;
            return a;
        }
        template<std::size_t dimension3>friend _matrix<Field, dimension1, dimension3> operator*(matrix_type a, _matrix<Field, dimension2, dimension3> b) {
            _matrix<Field, dimension1, dimension3> ret = _matrix<Field, dimension1, dimension3>(0);
            for(std::size_t i=0; i < dimension1; i++)
                for(std::size_t j=0; j < dimension3; j++)
                    ret.at(i,j) = a.const_row(i).get_vector().dot(b.const_column(j).get_vector());
            return ret;
        }
        transposed_matrix_type transpose() const {
            transposed_matrix_type ret=transposed_matrix_type(0);
            for(size_t i=0; i<dimension1; i++)
                ret.get_column(i).assign(const_row(i).begin());
            return ret;
        }
        transposed_matrix_type T() const {
            return transpose();
        }
        Field frobenius_norm() const {
            Field ret(0);
            for(std::size_t i=0; i < dimension1*dimension2; i++) ret += values[i]*values[i];
            return ret;
        }
        string to_string() const {
            ostringstream stream;
            for(std::size_t i=0; i < dimension1*dimension2; i++) {
                if (i % dimension2 == 0) stream << endl;
                else stream << " ";
                stream << values[i];
            }
            return stream.str();
        }
        friend string to_string(matrix_type a) {
            return a.to_string();
        }
        friend ostream& operator<<(ostream& os, const matrix_type& a) {
            os << a.to_string();
            return os;
        }
        friend inline int cmp(const matrix_type& a, const matrix_type& b) {
            return cmp(a.frobenius_norm()-b.frobenius_norm(), 0);
        }
        friend inline bool operator==(const matrix_type& a, const matrix_type& b){
            bool ret = true;
            for(size_t i=0; i < dimension1*dimension2; i++) {
                if(a.at(i)!=b.at(i)) {
                    ret = false;
                    break;
                }
            }
            return ret; 
        }
        friend inline bool operator!=(const matrix_type& a, const matrix_type& b){ return !(a==b); }
        friend inline bool operator< (const matrix_type& a, const matrix_type& b){ return cmp(a,b) <  0; }
        friend inline bool operator> (const matrix_type& a, const matrix_type& b){ return cmp(a,b) >  0; }
        friend inline bool operator<=(const matrix_type& a, const matrix_type& b){ return cmp(a,b) <= 0; }
        friend inline bool operator>=(const matrix_type& a, const matrix_type& b){ return cmp(a,b) >= 0; }        
        friend inline bool operator==(const matrix_type& a, const Field& b){ return cmp(a,b) == 0; }
        friend inline bool operator!=(const matrix_type& a, const Field& b){ return cmp(a,b) != 0; }
        friend inline bool operator< (const matrix_type& a, const Field& b){ return cmp(a,b) <  0; }
        friend inline bool operator> (const matrix_type& a, const Field& b){ return cmp(a,b) >  0; }
        friend inline bool operator<=(const matrix_type& a, const Field& b){ return cmp(a,b) <= 0; }
        friend inline bool operator>=(const matrix_type& a, const Field& b){ return cmp(a,b) >= 0; }
        friend inline matrix_type abs(const matrix_type& val) {
            return val.frobenius_norm();
        }
        friend inline pair<matrix_type, std::vector<matrix_type> > factor(const std::vector<matrix_type> &values, util::side_type side) {
            std::vector<matrix_type> ret(values.size());
            matrix_type common_factor = values[0];
            for(std::size_t i=0; i < dimension1*dimension2; i++) {
                std::vector<Field> tmp_input(values.size());
                for(size_t j=0; j < values.size(); j++) tmp_input[j] = values[j][i];
                pair<Field, std::vector<Field> > tmp = Field::factor(tmp_input);
                common_factor[i] = get<Field>(tmp);
                for(size_t j=0; j < ret.size(); j++) ret[j][i] = get<std::vector<Field> >(tmp)[j];
            }
            return pair<matrix_type, std::vector<matrix_type> >(common_factor, ret);
        }
        friend inline matrix_type exp(const matrix_type &_v) {
            matrix_type v = _v;
            for(auto value=v.values.begin(); value != v.values.end(); ++value) *value = Field::exp(*value);
            return v;
        }
        friend inline matrix_type log(const matrix_type &_v) {
            matrix_type v = _v;
            for(auto value=v.values.begin(); value != v.values.end(); ++value) *value = Field::log(*value);
            return v;
        }
}; // class _matrix

template<class Field, std::size_t dimension1, std::size_t dimension2=dimension1> class matrix : public _matrix<Field, dimension1, dimension2> {
    public:
        using _matrix<Field, dimension1, dimension2>::_matrix;
}; // class matrix - general case

template<class Field, std::size_t dimension> class matrix<Field, dimension, dimension> : public _matrix<Field, dimension, dimension> {
    protected:
        typedef vector<Field, dimension> vector_type;
        typedef _matrix<Field, dimension, dimension> parent_type;
        typedef matrix<Field, dimension, dimension> matrix_type;
    public:
        using _matrix<Field, dimension, dimension>::_matrix;
        matrix(const vector_type &_diagonal) : parent_type(0) { for(size_t i=0; i<dimension;i++) parent_type::at(i,i) = _diagonal[i]; }
        friend inline Field det(const matrix_type &a) {
            Field ret(0);
            int sign = 1;
            std::vector<std::size_t> permutation(dimension);
            for(std::size_t i=0; i<dimension; i++) permutation[i] = i;
            do {
                Field prod(sign);
                for(std::size_t i=0; i<dimension; i++) prod *= a.at(i, permutation[i]);
                ret += prod;
                sign *= -1;
            } while(std::next_permutation(permutation.begin(), permutation.end()));
            return ret;
        }
        static inline matrix_type I() {
            return matrix(vector_type::ones());
        }
        static inline matrix_type zeros() {
            return matrix_type();
        }
        static inline matrix_type ones() {
            return matrix_type(Field(1));
        }
}; // class matrix - square specialization

class Matrix33IntField {
    public:
        using value_type=matrix<int, 3, 3>;
        using power_type=unsigned int;
        using coefficient_type=int;
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
            enum { value = 1 };
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

}; // namespace fields
}; // namespace symba
#endif // _SYMBA_FIELDS_MATRIX_HPP_
