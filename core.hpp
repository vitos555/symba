
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
        using coefficient_allocator=allocator<Coefficient<FieldClass> >;
        using constant_allocator=allocator<Constant<FieldClass> >;
        using variable_allocator=allocator<Variable<FieldClass, Allocators<FieldClass> > >;
        using term_allocator=allocator<Term<FieldClass, Allocators<FieldClass> > >;
        using monomial_allocator=allocator<Monomial<FieldClass, Allocators<FieldClass> > >;
        using polynomial_allocator=allocator<Polynomial<FieldClass, Allocators<FieldClass> > >;
        using rational_allocator=allocator<Rational<FieldClass, Allocators<FieldClass> > >;
        using function_allocator=allocator<Function<FieldClass, Allocators<FieldClass>, FunctionAllocators<FieldClass, Allocators<FieldClass> > > >;
        using entity_allocator=allocator<Entity<FieldClass, Allocators<FieldClass> > >;
};

template<class FieldClass, class AllocatorsClass=Allocators<FieldClass> > class SubstitutionMap {
    private:
        using ValueType = typename FieldClass::value_type;
        using EntityClass = Entity<FieldClass, AllocatorsClass>;
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
        using ValueType = typename FieldClass::value_type;
        using CoefficientClass = Coefficient<FieldClass>;
        using ConstantClass = Constant<FieldClass>;
        using VariableClass = Variable<FieldClass, AllocatorsClass>;
        using MonomialClass = Monomial<FieldClass, AllocatorsClass>;
        using PolynomialClass = Polynomial<FieldClass, AllocatorsClass>;
        using RationalClass = Rational<FieldClass, AllocatorsClass>;
        using FunctionClass = Function<FieldClass, AllocatorsClass, FunctionAllocators<FieldClass, AllocatorsClass> >;
        using EntityClass = Entity<FieldClass, AllocatorsClass>;
        using variant_class = variant<CoefficientClass, ConstantClass, VariableClass, MonomialClass, PolynomialClass, RationalClass, FunctionClass>;
        variant_class obj;
    public:
        Entity() : obj(ConstantClass(0)) {}
        explicit Entity(const variant_class &_obj) : obj(_obj) { }
        Entity(const CoefficientClass &_o) : obj(_o) { }
        Entity(const ConstantClass &_c) : obj(_c) { }
        Entity(const VariableClass &_v) : obj(_v) { }
        Entity(const MonomialClass &_p) : obj(_p) { }
        Entity(const PolynomialClass &_p) : obj(_p) { }
        Entity(const RationalClass &_r) : obj(_r) { }
        Entity(const FunctionClass &_f) : obj(_f) { }

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
};

}; // namespace core
}; // namespace symba

#endif // _SYMBA_CORE_HPP_