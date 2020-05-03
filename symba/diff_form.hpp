
#ifndef _SYMBA_DIFF_FORM_HPP_
#define _SYMBA_DIFF_FORM_HPP_

#include "core.hpp"
#include "linalg.hpp"
#include "fields.hpp"
#include "util.hpp"
#include "vector.hpp"

namespace symba {
namespace linalg {
using namespace std;
using namespace util;
using namespace polynomial;
using namespace core;

template<class FieldClass,
		 std::size_t dimension,
		 std::size_t order,
		 class AllocatorsClass=Allocators<FieldClass>,
		 class VectorAllocator=class VectorAllocators<FieldClass, dimension, AllocatorsClass>::vector_allocator,
		 class _DiffFormAllocators=class DiffFormAllocators<FieldClass, dimension, order, AllocatorsClass> > 
class DiffForm {
    private:
    	enum terms { cnt = cexpr_limited_factorial(dimension, dimension-order) / cexpr_factorial(order) };
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
        typedef DiffForm<FieldClass, dimension, order, AllocatorsClass, VectorAllocator, _DiffFormAllocators> DiffFormClass;
        typedef DiffForm<FieldClass, dimension+1, order, AllocatorsClass, VectorAllocator, _DiffFormAllocators> DiffFormIncClass;
        typedef DiffForm<FieldClass, dimension, order, AllocatorsClass, VectorAllocator, _DiffFormAllocators> DiffFormDecClass;
        typedef array<EntityClass, terms::cnt> StorageClass;
        typedef typename StorageClass::iterator StorageIterator;
        typedef typename StorageClass::const_iterator StorageConstIterator;
        CoefficientAllocator oa;
        TermAllocator ta;
        MonomialAllocator ma;
        PolynomialAllocator pa;
        EntityAllocator ea;
        VectorAllocator lva;
        _DiffFormAllocators ldfa;
        StorageClass values;
    public:
        DiffForm() : values(), ta(), ma(), pa(), ea(), lva() { values.fill(EntityClass(ConstantClass(0))); }
        DiffForm(const EntityClass &_value) : values(), ta(), ma(), pa(), ea(), lva() { values.fill(_value); }
        DiffForm(const StorageClass &_values) : values(_values), ta(), ma(), pa(), ea(), lva() { }
        template<class InputIterator>DiffForm(InputIterator _begin, InputIterator _end) : values(), ta(), ma(), pa(), ea(), lva() { StorageIterator _cur=values.begin(); for(InputIterator _value=_begin; (_value!=_end) && (_cur!=values.end()); ++_value, ++_cur) *_cur = *_value; }

        const StorageClass& get_values() const { return values; }

        EntityClass& at(std::size_t i) { return values[i]; }
        EntityClass& operator[](std::size_t i) { return values[i]; }
        const EntityClass& at(std::size_t i) const { return values[i]; }
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
            stream << "99LDF|";
            for(auto value=values.begin(); value != values.end(); ++value) stream << value->signature() << "|";
            return stream.str();
        }

        ValueType evaluate(EvaluationMap<FieldClass> &_values) const {
            return ValueType(0);
        }

        // fields::vector<ValueType, dimension> evaluate_diff_form(EvaluationMap<FieldClass> &_values) const {
        //     fields::vector<ValueType, dimension> ret = fields::vector<ValueType, dimension>(0);
        //     for(size_t i=0; i < dimension; i++) ret[i] = at(i).evaluate(_values);
        //     return ret;
        // }

         d() {

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

        void add(const VectorClass& other) {
            for(std::size_t i=0; i<terms::cnt; i++) {
                values[i].add(other[i]);
            }
        }

        void multiply(const EntityClass &other) {
            for(std::size_t i=0; i<terms::cnt; i++) {
                values[i].multiply(other[i]);
            }
        }

        string to_string() const {
            ostringstream stream;
            stream << "(";
            for (auto value=values.begin(); value!=values.end(); ++value) {
                if (value != values.begin()) stream << ", ";
                stream << value->to_string();
            }
            stream << ")";
            return  stream.str();
        }        
};

template<class FieldClass,
		 std::size_t dimension,
		 class AllocatorsClass,
		 class VectorAllocator,
		 class DiffFormAllocator> 
class DiffForm<FieldClass, dimension, 1, AllocatorsClass, VectorAllocator, DiffFormAllocator> {
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
        typedef DiffForm<FieldClass, dimension, 1, AllocatorsClass, VectorAllocator, DiffFormAllocator> DiffFormClass;
        typedef array<EntityClass, dimension> StorageClass;
        typedef typename StorageClass::iterator StorageIterator;
        typedef typename StorageClass::const_iterator StorageConstIterator;
        CoefficientAllocator oa;
        TermAllocator ta;
        MonomialAllocator ma;
        PolynomialAllocator pa;
        EntityAllocator ea;
        VectorAllocator lva;
        DiffFormAllocator ldfa;
        StorageClass values;
    public:
        DiffForm() : values(), ta(), ma(), pa(), ea(), lva() { values.fill(EntityClass(ConstantClass(0))); }
        DiffForm(const EntityClass &_value) : values(), ta(), ma(), pa(), ea(), lva() { values.fill(_value); }
        DiffForm(const StorageClass &_values) : values(_values), ta(), ma(), pa(), ea(), lva() { }
        template<class InputIterator>DiffForm(InputIterator _begin, InputIterator _end) : values(), ta(), ma(), pa(), ea(), lva() { StorageIterator _cur=values.begin(); for(InputIterator _value=_begin; (_value!=_end) && (_cur!=values.end()); ++_value, ++_cur) *_cur = *_value; }

        const StorageClass& get_values() const { return values; }

        EntityClass& at(std::size_t i) { return values[i]; }
        EntityClass& operator[](std::size_t i) { return values[i]; }
        const EntityClass& at(std::size_t i) const { return values[i]; }
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
            stream << "99LDF|";
            for(auto value=values.begin(); value != values.end(); ++value) stream << value->signature() << "|";
            return stream.str();
        }

        ValueType evaluate(EvaluationMap<FieldClass> &_values) const {
            return ValueType(0);
        }

        // fields::vector<ValueType, dimension> evaluate_diff_form(EvaluationMap<FieldClass> &_values) const {
        //     fields::vector<ValueType, dimension> ret = fields::vector<ValueType, dimension>(0);
        //     for(size_t i=0; i < dimension; i++) ret[i] = at(i).evaluate(_values);
        //     return ret;
        // }

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

        void add(const VectorClass& other) {
            for(std::size_t i=0; i<dimension; i++) {
                values[i].add(other[i]);
            }
        }

        void multiply(const EntityClass &other) {
            for(std::size_t i=0; i<dimension; i++) {
                values[i].multiply(other[i]);
            }
        }

        EntityClass dot(const VectorClass& other) {
            EntityClass ret;
            auto _cur = values.begin();
            auto _other = other.values.begin();
            MonomialClass *mobj = ma.allocate(dimension);
            for(size_t i=0;(_cur != values.end())&&(_other != other.values.end()); ++_cur, ++_other, ++i) {
                TermClass *tobj = ta.allocate(2);
                allocator_traits<TermAllocator>::construct(ta, tobj, *_cur);
                allocator_traits<TermAllocator>::construct(ta, tobj+1, *_other);
                allocator_traits<MonomialAllocator>::construct(ma, mobj+i, vector<TermClass>({tobj[0], tobj[1]}));
                allocator_traits<TermAllocator>::destroy(ta, tobj);
                allocator_traits<TermAllocator>::destroy(ta, tobj+1);
                ta.deallocate(tobj, 2);
            }
            PolynomialClass *pobj = pa.allocate(1);
            allocator_traits<PolynomialAllocator>::construct(pa, pobj, vector<MonomialClass>(mobj, mobj+dimension));
            ret = EntityClass(pobj[0]);
            for(size_t i=0; i<dimension; i++) allocator_traits<MonomialAllocator>::destroy(ma, mobj + i);
            allocator_traits<PolynomialAllocator>::destroy(pa, pobj);
            ma.deallocate(mobj, dimension);
            pa.deallocate(pobj, 1);
            return ret;
        }

        string to_string() const {
            ostringstream stream;
            stream << "(";
            for (auto value=values.begin(); value!=values.end(); ++value) {
                if (value != values.begin()) stream << ", ";
                stream << value->to_string();
            }
            stream << ")";
            return  stream.str();
        }        
};



}; // namespace linalg
}; // namespace symba
#endif // _SYMBA_DIFF_FORM_HPP_