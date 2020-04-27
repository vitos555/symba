
#ifndef _SYMBA_MATRIX_HPP_
#define _SYMBA_MATRIX_HPP_

#include "core.hpp"
#include "linalg.hpp"
#include "fields.hpp"
#include "vector.hpp"

namespace symba {
namespace linalg {
using namespace std;
using namespace util;
using namespace polynomial;
using namespace core;

template<
    class FieldClass,
    std::size_t dimension1,
    std::size_t dimension2,
    class AllocatorsClass=Allocators<FieldClass>,
    class MatrixAllocator=typename MatrixAllocators<FieldClass, dimension1, dimension2, AllocatorsClass>::matrix_allocator,
    class TransposedMatrixAllocator=typename MatrixAllocators<FieldClass, dimension2, dimension1, AllocatorsClass>::matrix_allocator, 
    class ColumnVectorAllocator=typename VectorAllocators<FieldClass, dimension1, AllocatorsClass>::vector_allocator,
    class RowVectorAllocator=typename VectorAllocators<FieldClass, dimension2, AllocatorsClass>::vector_allocator> 
class Matrix {
    private:
        typedef typename FieldClass::value_type ValueType;
        typedef typename FieldClass::coefficient_type CoefficientType;
        typedef typename FieldClass::power_type PowerType;
        typedef typename AllocatorsClass::coefficient_allocator CoefficientAllocator;
        typedef typename AllocatorsClass::constant_allocator ConstantAllocator;
        typedef typename AllocatorsClass::variable_allocator VariableAllocator;
        typedef typename AllocatorsClass::term_allocator TermAllocator;
        typedef typename AllocatorsClass::monomial_allocator MonomialAllocator;
        typedef typename AllocatorsClass::polynomial_allocator PolynomialAllocator;
        typedef typename AllocatorsClass::rational_allocator RationalAllocator;
        typedef typename AllocatorsClass::entity_allocator EntityAllocator;
        typedef Coefficient<FieldClass> CoefficientClass;
        typedef Constant<FieldClass> ConstantClass;
        typedef Variable<FieldClass, AllocatorsClass> VariableClass;
        typedef Term<FieldClass, AllocatorsClass> TermClass;
        typedef Entity<FieldClass, AllocatorsClass> EntityClass;
        typedef Monomial<FieldClass, AllocatorsClass> MonomialClass;
        typedef Polynomial<FieldClass, AllocatorsClass> PolynomialClass;
        typedef Rational<FieldClass, AllocatorsClass> RationalClass;
        typedef Vector<FieldClass, dimension1, AllocatorsClass, ColumnVectorAllocator> ColumnVectorClass;
        typedef Vector<FieldClass, dimension2, AllocatorsClass, RowVectorAllocator> RowVectorClass;
        typedef Matrix<FieldClass, dimension1, dimension2, AllocatorsClass, MatrixAllocator, TransposedMatrixAllocator, ColumnVectorAllocator, RowVectorAllocator> MatrixClass;
        typedef Matrix<FieldClass, dimension2, dimension1, AllocatorsClass, TransposedMatrixAllocator, MatrixAllocator, RowVectorAllocator, ColumnVectorAllocator> TransposedMatrixClass;
        typedef array<EntityClass, dimension1*dimension2> StorageClass;
        typedef typename StorageClass::iterator StorageIterator;
        typedef typename StorageClass::const_iterator StorageConstIterator;
        typedef VectorView<FieldClass, StorageIterator, dimension1, AllocatorsClass, ColumnVectorAllocator> ColumnVectorViewClass;
        typedef VectorView<FieldClass, StorageIterator, dimension2, AllocatorsClass, RowVectorAllocator> RowVectorViewClass;
        typedef VectorView<FieldClass, StorageConstIterator, dimension1, AllocatorsClass, ColumnVectorAllocator> ColumnConstVectorViewClass;
        typedef VectorView<FieldClass, StorageConstIterator, dimension2, AllocatorsClass, RowVectorAllocator> RowConstVectorViewClass;
        CoefficientAllocator oa;
        TermAllocator ta;
        MonomialAllocator ma;
        PolynomialAllocator pa;
        EntityAllocator ea;
        MatrixAllocator lma;
        ColumnVectorAllocator lcva;
        RowVectorAllocator lrva;
        StorageClass values;
    public:
        Matrix() : values(), ta(), ma(), pa(), ea(), lma(), lcva(), lrva() { values.fill(EntityClass(ConstantClass(0))); }
        Matrix(const EntityClass &_value) : values(), ta(), ma(), pa(), ea(), lma(), lcva(), lrva() { values.fill(_value); }
        Matrix(const StorageClass &_values) : values(_values), ta(), ma(), pa(), ea(), lma(), lcva(), lrva() { }
        template<class InputIterator>Matrix(InputIterator _begin, InputIterator _end) : values(), ta(), ma(), pa(), ea(), lma(), lcva(), lrva() { StorageIterator _cur=values.begin(); for(InputIterator _value=_begin; (_value!=_end()) && (_cur!=values.end()); ++_value, ++_cur) *_cur = *_value; }

        const StorageClass& get_values() const { return values; }

        EntityClass& at(std::size_t i) { return values[i]; }
        EntityClass& at(std::size_t i, std::size_t j) { return values[i*dimension1 + j]; }
        EntityClass& operator[](std::size_t i) { return values[i]; }
        const EntityClass& at(std::size_t i) const { return values[i]; }
        const EntityClass& at(std::size_t i, std::size_t j) const { return values[i*dimension1 + j]; }
        const EntityClass& operator[](std::size_t i) const { return values[i]; }

        bool is_zero() const {
            bool ret = true;
            for(auto value=values.begin(); value != values.end(); ++value) {
                if (!value->is_zero()) {
                    ret = false;
                    break;
                }
            }
            return ret;
        }

        string signature() const {
            ostringstream stream;
            stream << "99LM|";
            for(auto value=values.begin(); value != values.end(); ++value) stream << value->signature() << "|";
            return stream.str();
        }

        ValueType evaluate(EvaluationMap<FieldClass> &_values) const {
            return ValueType(0);
        }

        fields::matrix<ValueType, dimension1, dimension2> evaluate_matrix(EvaluationMap<FieldClass> &_values) const {
            fields::matrix<ValueType, dimension1, dimension2> ret = fields::matrix<ValueType, dimension1, dimension2>(0);
            for(size_t i=0; i < dimension1*dimension2; i++) ret.at(i) = at(i).evaluate(_values);
            return ret;
        }

        void simplify() {
            for(auto value=values.begin(); value != values.end(); ++value) value->simplify();
        }

        void expand() {
            for(auto value=values.begin(); value != values.end(); ++value) value->expand();
        }

        void factor(side_type side=side_type::both) {
            for(auto value=values.begin(); value != values.end(); ++value) value->factor();
        }

        void flatten() {
            for(auto value=values.begin(); value != values.end(); ++value) value->flatten();
        }

        void substitute(SubstitutionMap<FieldClass> &_values) {
            for(auto value=values.begin(); value != values.end(); ++value) value->substitute(_values);
        }

        RowVectorViewClass get_row(std::size_t i) {
            return RowVectorViewClass(values.begin() + i*dimension2, 1);
        }
        RowConstVectorViewClass const_row(std::size_t i) const {
            return RowConstVectorViewClass(values.begin() + i*dimension2, 1);
        }
        ColumnVectorViewClass get_column(std::size_t i) {
            return ColumnVectorViewClass(values.begin() + i, dimension2);
        }
        ColumnConstVectorViewClass const_column(std::size_t i) const {
            return ColumnConstVectorViewClass(values.begin() + i, dimension2);
        }

        TransposedMatrixClass transpose() const {
            TransposedMatrixClass ret = TransposedMatrixClass();
            for(size_t i=0; i<dimension1; i++) {
                ret.get_column(i).assign(const_row(i).begin());
            }
            return ret;            
        }
        TransposedMatrixClass T() const {
            return transpose();
        }

        Matrix<FieldClass, dimension1, dimension1, AllocatorsClass, MatrixAllocator, TransposedMatrixAllocator, ColumnVectorAllocator, RowVectorAllocator> multiply(const TransposedMatrixClass &other) {
            Matrix<FieldClass, dimension1, dimension1, AllocatorsClass, MatrixAllocator, TransposedMatrixAllocator, ColumnVectorAllocator, RowVectorAllocator> ret = Matrix<FieldClass, dimension1, dimension1, AllocatorsClass, MatrixAllocator, TransposedMatrixAllocator, ColumnVectorAllocator, RowVectorAllocator>();
            for(std::size_t i=0; i < dimension1; i++)
                for(std::size_t j=0; j < dimension1; j++)
                    ret.at(i,j) = const_row(i).get_vector().dot(other.const_column(j).get_vector());
            return ret;            
        }

        void add(const MatrixClass& other) {
            for(std::size_t i=0; i<dimension1*dimension2; i++) values[i].add(other[i]);
        }

        void multiply(const EntityClass &other) {
            for(std::size_t i=0; i<dimension1*dimension2; i++) values[i].multiply(other[i]);
        }

        string to_string() const {
            ostringstream stream;
            for (std::size_t i=0; i < dimension1*dimension2; i++) {
                if (i % dimension2 == 0) stream << endl;
                else stream << " | ";
                stream << at(i).to_string();
            }
            stream << endl;
            return  stream.str();
        }        
};

// Square matrix
template<
    class FieldClass,
    std::size_t dimension,
    class AllocatorsClass,
    class MatrixAllocator, 
    class VectorAllocator
    > 
class Matrix<FieldClass, dimension, dimension, AllocatorsClass, MatrixAllocator, MatrixAllocator, VectorAllocator, VectorAllocator>  {
    private:
        typedef typename FieldClass::value_type ValueType;
        typedef typename FieldClass::coefficient_type CoefficientType;
        typedef typename FieldClass::power_type PowerType;
        typedef typename AllocatorsClass::coefficient_allocator CoefficientAllocator;
        typedef typename AllocatorsClass::constant_allocator ConstantAllocator;
        typedef typename AllocatorsClass::variable_allocator VariableAllocator;
        typedef typename AllocatorsClass::term_allocator TermAllocator;
        typedef typename AllocatorsClass::monomial_allocator MonomialAllocator;
        typedef typename AllocatorsClass::polynomial_allocator PolynomialAllocator;
        typedef typename AllocatorsClass::rational_allocator RationalAllocator;
        typedef typename AllocatorsClass::entity_allocator EntityAllocator;
        typedef Coefficient<FieldClass> CoefficientClass;
        typedef Constant<FieldClass> ConstantClass;
        typedef Variable<FieldClass, AllocatorsClass> VariableClass;
        typedef Term<FieldClass, AllocatorsClass> TermClass;
        typedef Entity<FieldClass, AllocatorsClass> EntityClass;
        typedef Monomial<FieldClass, AllocatorsClass> MonomialClass;
        typedef Polynomial<FieldClass, AllocatorsClass> PolynomialClass;
        typedef Rational<FieldClass, AllocatorsClass> RationalClass;
        typedef Vector<FieldClass, dimension, AllocatorsClass, VectorAllocator> VectorClass;
        typedef Matrix<FieldClass, dimension, dimension, AllocatorsClass, MatrixAllocator, MatrixAllocator, VectorAllocator, VectorAllocator> MatrixClass;
        typedef array<EntityClass, dimension*dimension> StorageClass;
        typedef typename StorageClass::iterator StorageIterator;
        typedef typename StorageClass::const_iterator StorageConstIterator;
        typedef VectorView<FieldClass, StorageIterator, dimension, AllocatorsClass, VectorAllocator> VectorViewClass;
        typedef VectorView<FieldClass, StorageConstIterator, dimension, AllocatorsClass, VectorAllocator> ConstVectorViewClass;
        CoefficientAllocator oa;
        TermAllocator ta;
        MonomialAllocator ma;
        PolynomialAllocator pa;
        EntityAllocator ea;
        MatrixAllocator lma;
        VectorAllocator lva;
        StorageClass values;
    public:
        Matrix() : values(), ta(), ma(), pa(), ea(), lma(), lva() { values.fill(EntityClass(ConstantClass(0))); }
        Matrix(const EntityClass &_value) : values(), ta(), ma(), pa(), ea(), lma(), lva() { values.fill(_value); }
        Matrix(const StorageClass &_values) : values(_values), ta(), ma(), pa(), ea(), lma(), lva() { }
        Matrix(const VectorClass &_values) : Matrix() { for (std::size_t i=0; i < dimension; i++) at(i, i) = _values[i]; }
        template<class InputIterator>Matrix(InputIterator _begin, InputIterator _end) : values(), ta(), ma(), pa(), ea(), lma(), lva() { StorageIterator _cur=values.begin(); for(InputIterator _value=_begin; (_value!=_end()) && (_cur!=values.end()); ++_value, ++_cur) *_cur = *_value; }

        const StorageClass& get_values() const { return values; }

        EntityClass& at(std::size_t i) { return values[i]; }
        EntityClass& at(std::size_t i, std::size_t j) { return values[i*dimension + j]; }
        EntityClass& operator[](std::size_t i) { return values[i]; }
        const EntityClass& at(std::size_t i) const { return values[i]; }
        const EntityClass& at(std::size_t i, std::size_t j) const { return values[i*dimension + j]; }
        const EntityClass& operator[](std::size_t i) const { return values[i]; }

        bool is_zero() const {
            bool ret = true;
            for(auto value=values.begin(); value != values.end(); ++value) {
                if (!value->is_zero()) {
                    ret = false;
                    break;
                }
            }
            return ret;
        }

        string signature() const {
            ostringstream stream;
            stream << "99LM2|";
            for(auto value=values.begin(); value != values.end(); ++value) stream << value->signature() << "|";
            return stream.str();
        }

        ValueType evaluate(EvaluationMap<FieldClass> &_values) const {
            return ValueType(0);
        }

        fields::matrix<ValueType, dimension, dimension> evaluate_matrix(EvaluationMap<FieldClass> &_values) const {
            fields::matrix<ValueType, dimension, dimension> ret = fields::matrix<ValueType, dimension, dimension>(0);
            for(size_t i=0; i < dimension*dimension; i++) ret.at(i) = at(i).evaluate(_values);
            return ret;
        }

        void simplify() {
            for(auto value=values.begin(); value != values.end(); ++value) value->simplify();
        }

        void expand() {
            for(auto value=values.begin(); value != values.end(); ++value) value->expand();
        }

        void factor(side_type side=side_type::both) {
            for(auto value=values.begin(); value != values.end(); ++value) value->factor();
        }

        void flatten() {
            for(auto value=values.begin(); value != values.end(); ++value) value->flatten();
        }

        void substitute(SubstitutionMap<FieldClass> &_values) {
            for(auto value=values.begin(); value != values.end(); ++value) value->substitute(_values);
        }

        VectorViewClass get_row(std::size_t i) {
            return VectorViewClass(values.begin() + i*dimension, 1);
        }
        ConstVectorViewClass const_row(std::size_t i) const {
            return ConstVectorViewClass(values.begin() + i*dimension, 1);
        }
        VectorViewClass get_column(std::size_t i) {
            return VectorViewClass(values.begin() + i, dimension);
        }
        ConstVectorViewClass const_column(std::size_t i) const {
            return ConstVectorViewClass(values.begin() + i, dimension);
        }

        MatrixClass transpose() const {
            MatrixClass ret = MatrixClass();
            for(size_t i=0; i<dimension; i++) {
                ret.get_column(i).assign(const_row(i).begin());
            }
            return ret;            
        }
        MatrixClass T() const {
            return transpose();
        }

        MatrixClass multiply(const MatrixClass &other) {
            MatrixClass ret = MatrixClass();
            for(std::size_t i=0; i < dimension; i++)
                for(std::size_t j=0; j < dimension; j++)
                    ret.at(i,j) = const_row(i).get_vector().dot(other.const_column(j).get_vector());
            return ret;            
        }

        void add(const MatrixClass& other) {
            for(std::size_t i=0; i<dimension*dimension; i++) values[i].add(other[i]);
        }

        void multiply(const EntityClass &other) {
            for(std::size_t i=0; i<dimension*dimension; i++) values[i].multiply(other[i]);
        }

        EntityClass det() {
            EntityClass ret;
            int sign = 1;
            std::size_t k=0;
            std::vector<std::size_t> permutation(dimension);
            std::size_t N = util::factorial<std::size_t>(dimension);
            MonomialClass *mobj = ma.allocate(N);
            for(std::size_t i=0; i<dimension; i++) permutation[i] = i;
            do {
                TermClass *tobj = ta.allocate(dimension+1);
                allocator_traits<TermAllocator>::construct(ta, tobj, ConstantClass(sign));
                for(std::size_t i=0; i<dimension; i++) allocator_traits<TermAllocator>::construct(ta, tobj+i+1, at(i, permutation[i]));
                allocator_traits<MonomialAllocator>::construct(ma, mobj+k, std::vector<TermClass>(tobj, tobj+dimension+1));
                for(std::size_t i=0; i < dimension + 1; i++) allocator_traits<TermAllocator>::destroy(ta, tobj+i);
                k++;
                sign *= -1;
                ta.deallocate(tobj, dimension+1);
            } while(std::next_permutation(permutation.begin(), permutation.end()));
            PolynomialClass *pobj = pa.allocate(1);
            allocator_traits<PolynomialAllocator>::construct(pa, pobj, std::vector<MonomialClass>(mobj, mobj+N));
            ret = EntityClass(pobj[0]);
            allocator_traits<PolynomialAllocator>::destroy(pa, pobj);
            for(std::size_t i=0; i<N; i++) allocator_traits<MonomialAllocator>::destroy(ma, mobj+i);
            pa.deallocate(pobj, 1);
            ma.deallocate(mobj, N);
            return ret;
        }

        string to_string() const {
            ostringstream stream;
            for (std::size_t i=0; i < dimension*dimension; i++) {
                if (i % dimension == 0) stream << endl;
                else stream << " | ";
                stream << at(i).to_string();
            }
            stream << endl;
            return  stream.str();
        }        
};





}; // namespace linalg
}; // namespace symba
#endif // _SYMBA_MATRIX_HPP_