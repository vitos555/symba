
#ifndef _SYMBA_CORE_HPP_
#define _SYMBA_CORE_HPP_

#include<cmath>
#include<algorithm>
#include<variant>
#include<string>
#include<vector>
#include<iostream>
#include<map>
#include<numeric>
#include<sstream>
#include<type_traits>
#include<iomanip>

#include "util.hpp"

namespace symba {

namespace polynomial {

template<class FieldClass> class Coefficient;
template<class FieldClass> class Constant;
template<class FieldClass, class Allocators> class Variable;
template<class FieldClass, class Allocators> class Term;
template<class FieldClass, class Allocators> class Monomial;
template<class FieldClass, class Allocators> class Polynomial;

}; // namespace polynomial

namespace rational {

template<class FieldClass, class Allocators> class Rational;

}; // namespace rational

namespace function {

template<class FieldClass, class Allocators> class FunctionAllocators;
template<class FieldClass, class Allocators, class FunctionAllocatorsClass> class Function;

}; // namespace function

namespace core {
using namespace std;
using namespace util;
using namespace polynomial;
using namespace rational;
using namespace function;

template<class FieldClass> class Allocators;
template<class FieldClass, class Allocators> class Entity;

template<class FieldClass> class Allocators {
    public:
        typedef allocator<Coefficient<FieldClass> > coefficient_allocator;
        typedef allocator<Constant<FieldClass> > constant_allocator;
        typedef allocator<Variable<FieldClass, Allocators<FieldClass> > > variable_allocator;
        typedef allocator<Term<FieldClass, Allocators<FieldClass> > > term_allocator;
        typedef allocator<Monomial<FieldClass, Allocators<FieldClass> > > monomial_allocator;
        typedef allocator<Polynomial<FieldClass, Allocators<FieldClass> > > polynomial_allocator;
        typedef allocator<Rational<FieldClass, Allocators<FieldClass> > > rational_allocator;
        typedef allocator<Function<FieldClass, Allocators<FieldClass>, FunctionAllocators<FieldClass, Allocators<FieldClass> > > > function_allocator;
        typedef allocator<Entity<FieldClass, Allocators<FieldClass> > > entity_allocator;
};

template<class FieldClass, class AllocatorsClass=Allocators<FieldClass> > class SubstitutionMap {
    private:
        typedef typename FieldClass::value_type ValueType;
        typedef Entity<FieldClass, AllocatorsClass> EntityClass;
        map<string, EntityClass> values;
    public:
        explicit SubstitutionMap(const map<string, EntityClass> &_values) : values(_values) { }
        size_t count(string name) {
            return values.count(name);
        }
        const EntityClass &get(string name) {
            return values[name];
        }
};

template<class FieldClass, class AllocatorsClass=Allocators<FieldClass> > class Entity {
    private:
        typedef typename FieldClass::value_type ValueType;
        typedef Coefficient<FieldClass> CoefficientClass;
        typedef Constant<FieldClass> ConstantClass;
        typedef Variable<FieldClass, AllocatorsClass> VariableClass;
        typedef Term<FieldClass, AllocatorsClass> TermClass;
        typedef Monomial<FieldClass, AllocatorsClass> MonomialClass;
        typedef Polynomial<FieldClass, AllocatorsClass> PolynomialClass;
        typedef Rational<FieldClass, AllocatorsClass> RationalClass;
        typedef Function<FieldClass, AllocatorsClass, FunctionAllocators<FieldClass, AllocatorsClass> > FunctionClass;
        typedef Entity<FieldClass, AllocatorsClass> EntityClass;
        typedef typename AllocatorsClass::term_allocator TermAllocator;
        typedef typename AllocatorsClass::monomial_allocator MonomialAllocator;
        typedef typename AllocatorsClass::polynomial_allocator PolynomialAllocator;
        typedef variant<CoefficientClass, ConstantClass, VariableClass, MonomialClass, PolynomialClass, RationalClass, FunctionClass> variant_class;
        variant_class obj;
        TermAllocator ta;
        MonomialAllocator ma;
        PolynomialAllocator pa;
    public:
        Entity() : obj(ConstantClass(0)), ta(), ma(), pa() {}
        explicit Entity(const variant_class &_obj) : obj(_obj), ta(), ma(), pa() { }
        Entity(const CoefficientClass &_o) : obj(_o), ta(), ma(), pa() { }
        Entity(const ConstantClass &_c) : obj(_c), ta(), ma(), pa() { }
        Entity(const VariableClass &_v) : obj(_v), ta(), ma(), pa() { }
        Entity(const MonomialClass &_p) : obj(_p), ta(), ma(), pa() { }
        Entity(const PolynomialClass &_p) : obj(_p), ta(), ma(), pa() { }
        Entity(const RationalClass &_r) : obj(_r), ta(), ma(), pa() { }
        Entity(const FunctionClass &_f) : obj(_f), ta(), ma(), pa() { }

        bool is_coefficient() const {
            return holds_alternative<ConstantClass>(obj);
        }
        bool is_constant() const {
            return holds_alternative<ConstantClass>(obj);
        }
        bool is_variable() const {
            return holds_alternative<VariableClass>(obj);
        }
        bool is_monomial() const {
            return holds_alternative<MonomialClass>(obj);
        }
        bool is_polynomial() const {
            return holds_alternative<PolynomialClass>(obj);
        }
        bool is_rational() const {
            return holds_alternative<RationalClass>(obj);
        }
        bool is_function() const {
            return holds_alternative<FunctionClass>(obj);
        }

        string signature() const {
            return visit([](auto &&arg) -> string { return arg.signature(); }, obj);
        }
        ValueType evaluate(EvaluationMap<FieldClass> &values) const {
            return visit([&values](auto &&arg) -> ValueType { return arg.evaluate(values); }, obj);
        }
        void substitute(SubstitutionMap<FieldClass> &values) {
            if (holds_alternative<ConstantClass>(obj) || holds_alternative<CoefficientClass>(obj)) {
                // do nothing
            } else if (holds_alternative<VariableClass>(obj)) {
                // do substitution
                const VariableClass &var = get<VariableClass>(obj);
                if (values.count(var.get_name()) > 0) {
                    obj = values.get(var.get_name()).get_obj();
                }
            } else {
                // do recursive substitution
                visit([&values](auto &&arg) { arg.substitute(values); }, obj);
            }
        }
        void simplify() {
            visit([](auto &&arg) { arg.simplify(); }, obj);
        }
        void expand() {
            visit([](auto &&arg) { arg.expand(); }, obj);
        }
        void flatten() {
            visit([](auto &&arg) { arg.flatten(); }, obj);
        }
        void factor(side_type side=side_type::both) {
            visit([side](auto &&arg) { arg.factor(side); }, obj);
        }
        string to_string() const {
            return visit([](auto &&arg) -> string { return arg.to_string(); }, obj);
        }
        void replace(const variant_class &_obj) {
            obj = _obj;
        }
        const variant_class &get_obj() const {
            return obj;
        }
        void add(const EntityClass& other, const side_type &side=side_type::right) {
            if (holds_alternative<PolynomialClass>(obj)) {
                PolynomialClass &pobj = get<PolynomialClass>(obj);
                pobj.add(other, side);
            } else {
                TermClass *tobj = ta.allocate(1);
                MonomialClass *mobj = ma.allocate(1);
                PolynomialClass *pobj = pa.allocate(1);
                allocator_traits<TermAllocator>::construct(ta, tobj, *this);
                allocator_traits<MonomialAllocator>::construct(ma, mobj, vector<TermClass>({tobj[0]}));
                allocator_traits<PolynomialAllocator>::construct(pa, pobj, vector<MonomialClass>({mobj[0]}));
                pobj[0].add(other, side);
                obj = pobj[0];
                allocator_traits<PolynomialAllocator>::destroy(pa, pobj);
                allocator_traits<MonomialAllocator>::destroy(ma, mobj);
                allocator_traits<TermAllocator>::destroy(ta, tobj);
                pa.deallocate(pobj, 1);
                ma.deallocate(mobj, 1);
                ta.deallocate(tobj, 1);
            }
        }
        void multiply(const EntityClass &other, const side_type &side=side_type::right) {
            if (holds_alternative<MonomialClass>(obj)) {
                const MonomialClass &mobj = get<MonomialClass>(obj);
                mobj.multiply(other, side);
            } else {
                TermClass *tobj = ta.allocate(1);
                MonomialClass *mobj = ma.allocate(1);
                allocator_traits<TermAllocator>::construct(ta, tobj, *this);
                allocator_traits<MonomialAllocator>::construct(ma, mobj, vector<TermClass>({tobj[0]}));
                mobj[0].multiply(other, side);
                obj = mobj[0];
                allocator_traits<MonomialAllocator>::destroy(ma, mobj);
                allocator_traits<TermAllocator>::destroy(ta, tobj);
                ma.deallocate(mobj, 1);
                ta.deallocate(tobj, 1);
            }
        }
};

}; // namespace core
}; // namespace symba

#endif // _SYMBA_CORE_HPP_