
#ifndef _SYMBA_FUNCTION_HPP_
#define _SYMBA_FUNCTION_HPP_

#include "core.hpp"

namespace symba {
namespace function {
using namespace std;
using namespace util;
using namespace core;

template<class FieldClass, class AllocatorsClass, class FunctionAllocatorsClass> class FunctionImplementation;
template<class FieldClass, class AllocatorsClass, class FunctionAllocatorsClass> class IdentityImplementation;
template<class FieldClass, class AllocatorsClass, class FunctionAllocatorsClass> class ExpImplementation;
template<class FieldClass, class AllocatorsClass, class FunctionAllocatorsClass> class LogImplementation;

template<class FieldClass, class AllocatorsClass> class FunctionAllocators {
    public:
        using function_implementation_allocator=allocator<FunctionImplementation<FieldClass, Allocators<FieldClass>, FunctionAllocators<FieldClass, AllocatorsClass> > >;
        using identity_implementation_allocator=allocator<IdentityImplementation<FieldClass, AllocatorsClass, FunctionAllocators<FieldClass, AllocatorsClass> > >;    
        using exp_implementation_allocator=allocator<ExpImplementation<FieldClass, AllocatorsClass, FunctionAllocators<FieldClass, AllocatorsClass> > >;
        using log_implementation_allocator=allocator<LogImplementation<FieldClass, AllocatorsClass, FunctionAllocators<FieldClass, AllocatorsClass> > >;    
};

template<class FieldClass, class AllocatorsClass=Allocators<FieldClass>, class FunctionAllocatorsClass=FunctionAllocators<FieldClass, Allocators<FieldClass> > > class IdentityImplementation {
    protected:
        using ValueType = typename FieldClass::value_type;
        using EntityClass = Entity<FieldClass, AllocatorsClass>;
        using VariableClass = Variable<FieldClass, AllocatorsClass>;
        using IdentityImpClass = IdentityImplementation<FieldClass, AllocatorsClass>;
        EntityClass argument;
    public:
        explicit IdentityImplementation(const EntityClass &_argument) : argument(_argument) { }
        const EntityClass &get_argument() const {
            return argument;
        }
        string signature() const {
            return "7F0001";
        }
        string get_name() const {
         return "identity";
        }
        string get_inverse_name() const {
         return "";
        }
        void simplify() {
            argument.simplify();
        }
        ValueType evaluate(EvaluationMap<FieldClass> &values) const { 
         return argument.evaluate(values);
        }
        void substitute(SubstitutionMap<FieldClass> &values) {
            argument.substitute(values);
        }
        void expand() {
            argument.expand();
        }
        void flatten() {
            argument.flatten();
        }
        void factor() {
            argument.factor();
        }
        void derivative(const VariableClass &by_variable) {
            argument.derivative(by_variable);
        }
        string to_string() const {
            return "I(" + argument.to_string() + ")";
        }
        friend ostream& operator<<(ostream& os, const IdentityImpClass& i) {
            os << i.to_string();
            return os;
        }
};

template<class FieldClass, class AllocatorsClass=Allocators<FieldClass>, class FunctionAllocatorsClass=FunctionAllocators<FieldClass, Allocators<FieldClass> > > class ExpImplementation {
    protected:
        using ValueType = typename FieldClass::value_type;
        using EntityClass = Entity<FieldClass, AllocatorsClass>;
        using VariableClass = Variable<FieldClass, AllocatorsClass>;
        using ExpImpClass = ExpImplementation<FieldClass, AllocatorsClass>;
        EntityClass argument;
    public:
        explicit ExpImplementation(const EntityClass &_argument) : argument(_argument) { }
        const EntityClass &get_argument() const {
            return argument;
        }
        string signature() const {
            return "7F0002";
        }
        string get_name() const {
         return "exp";
        }
        string get_inverse_name() const {
            return "log";
        }
        void simplify() {
            argument.simplify();
        }
        ValueType evaluate(EvaluationMap<FieldClass> &values) const { 
         return FieldClass::exp(argument.evaluate(values));
        }
        void substitute(SubstitutionMap<FieldClass> &values) {
            argument.substitute(values);
        }
        void expand() {
            argument.expand();
        }
        void flatten() {
            argument.flatten();
        }
        void factor() {
            argument.factor();
        }
        void derivative(const VariableClass &by_variable) { }
        string to_string() const {
         return "Exp(" + argument.to_string() + ")";
        }
        friend ostream& operator<<(ostream& os, const ExpImpClass& e) {
            os << e.to_string();
            return os;
        }
};

template<class FieldClass, class AllocatorsClass=Allocators<FieldClass>, class FunctionAllocatorsClass=FunctionAllocators<FieldClass, AllocatorsClass> > class LogImplementation {
    protected:
        using ValueType = typename FieldClass::value_type;
        using EntityClass = Entity<FieldClass, AllocatorsClass>;
        using VariableClass = Variable<FieldClass, AllocatorsClass>;
        using LogImpClass = LogImplementation<FieldClass, AllocatorsClass>;
        EntityClass argument;
    public:
        explicit LogImplementation(const EntityClass &_argument) : argument(_argument) { }
        const EntityClass &get_argument() const {
            return argument;
        }
        string signature() const {
            return "7F0003";
        }
        string get_name() const {
            return "log";
        }
        string get_inverse_name() const {
            return "exp";
        }
        void simplify() {
            argument.simplify();
        }
        ValueType evaluate(EvaluationMap<FieldClass> &values) const { 
         return FieldClass::log(argument.evaluate(values));
        }
        void substitute(SubstitutionMap<FieldClass> &values) {
            argument.substitute(values);
        }
        void expand() {
            argument.expand();
        }
        void flatten() {
            argument.flatten();
        }
        void factor() {
            argument.factor();
        }
        void derivative(const VariableClass &by_variable) { }
        string to_string() const {
         return "Log(" + argument.to_string() + ")";
        }
        friend ostream& operator<<(ostream& os, const LogImpClass& l) {
            os << l.to_string();
            return os;
        }
};

template<class FieldClass, class AllocatorsClass=Allocators<FieldClass>, class FunctionAllocatorsClass=FunctionAllocators<FieldClass, AllocatorsClass> > class FunctionImplementation {
    private:
        using ValueType = typename FieldClass::value_type;
        using EntityClass = Entity<FieldClass, AllocatorsClass>;
        using VariableClass = Variable<FieldClass, AllocatorsClass>;
        using FunctionClass = Function<FieldClass, AllocatorsClass, FunctionAllocatorsClass>;
        using FunctionImplementationClass = FunctionImplementation<FieldClass, AllocatorsClass, FunctionAllocatorsClass>;
        using IdImpClass = IdentityImplementation<FieldClass, AllocatorsClass>;
        using ExpImpClass = ExpImplementation<FieldClass, AllocatorsClass>;
        using LogImpClass = LogImplementation<FieldClass, AllocatorsClass>;
        using IdentityAllocator = typename FunctionAllocatorsClass::identity_implementation_allocator;
        using implementation_class = variant< ExpImpClass, LogImpClass, IdImpClass >;
        implementation_class implementation;
        IdentityAllocator ida;
    public:
        explicit FunctionImplementation(implementation_class &_implementation) : implementation(_implementation), ida() { }
        FunctionImplementation(ExpImpClass &_exp) : implementation(_exp), ida() { }
        FunctionImplementation(LogImpClass &_log) : implementation(_log), ida() { }
        FunctionImplementation(IdImpClass &_id) : implementation(_id), ida() { }
        string signature() {
            return visit([](auto &&arg) -> string { return arg.signature(); }, implementation);
        }
        string get_name() const { 
            return visit([](auto &&arg) -> string { return arg.get_name(); }, implementation);
        }
        string get_inverse_name() const { 
            return visit([](auto &&arg) -> string { return arg.get_inverse_name(); }, implementation);
        }
        ValueType evaluate(EvaluationMap<FieldClass> &values) const {
            return visit([&values](auto &&arg) -> ValueType { return arg.evaluate(values); }, implementation);
        }
        void simplify() {
            visit([](auto &&arg) { arg.simplify(); }, implementation);
            if (get_argument().is_function()) {
                const FunctionClass &fn = get<FunctionClass>(get_argument().get_obj());
                if (get_name().compare(fn.get_inverse_name())==0) {
                    IdImpClass *iobj = ida.allocate(1);
                    allocator_traits<IdentityAllocator>::construct(ida, iobj, fn.get_implementation().get_argument());
                    replace(iobj[0]);
                    allocator_traits<IdentityAllocator>::destroy(ida, iobj);
                    ida.deallocate(iobj, 1);
                }
            }
        }
        void factor() {
            visit([](auto &&arg) { arg.factor(); }, implementation);
        }
        void flatten() {
            visit([](auto &&arg) { arg.flatten(); }, implementation);
        }
        void expand() {
            visit([](auto &&arg) { arg.expand(); }, implementation);
        }
        void substitute(SubstitutionMap<FieldClass> &values) {
            visit([&values](auto &&arg) { arg.substitute(values); }, implementation);
        }
        void derivative(const VariableClass& by_variable) {
            visit([by_variable](auto &&arg) { arg.derivative(by_variable); }, implementation);
        }
        string to_string() const {
            return visit([](auto &&arg) -> string { return arg.to_string(); }, implementation);
        }
        const EntityClass &get_argument() const {
            return visit([](auto &&arg) -> const EntityClass& { return arg.get_argument(); }, implementation);
        }
        void replace(const implementation_class &_implementation) {
            implementation = _implementation;
        }
        const implementation_class &get_implementation() const {
            return implementation;
        }
        friend ostream& operator<<(ostream& os, const FunctionImplementationClass& f) {
            os << f.to_string();
            return os;
        }
};

template<class FieldClass, class AllocatorsClass=Allocators<FieldClass>, class FunctionAllocatorsClass=FunctionAllocators<FieldClass, AllocatorsClass> > class Function {
    private:
        using ValueType = typename FieldClass::value_type;
        using EntityClass = Entity<FieldClass, AllocatorsClass>;
        using VariableClass = Variable<FieldClass, AllocatorsClass>;
        using FunctionClass = Function<FieldClass, AllocatorsClass>;
        using ImplementationClass = FunctionImplementation<FieldClass, AllocatorsClass, FunctionAllocatorsClass>;
        using ImplementationAllocator = typename FunctionAllocatorsClass::function_implementation_allocator;
        ImplementationClass *implementation;
        ImplementationAllocator ia;
    public:
        Function(const Function &_function) : ia() {
            implementation = ia.allocate(1);
            allocator_traits<ImplementationAllocator>::construct(ia, implementation, _function.get_implementation());
        }
        ~Function() {
            allocator_traits<ImplementationAllocator>::destroy(ia, implementation);
            ia.deallocate(implementation, 1);
        }
        explicit Function(ImplementationClass &_implementation) : ia() {
            implementation = ia.allocate(1);
            allocator_traits<ImplementationAllocator>::construct(ia, implementation, _implementation);
        }
        string signature() const {
            return "7F_" + implementation->signature();
        }

        string get_name() const {
            return implementation->get_name();
        }

        string get_inverse_name() const {
            return implementation->get_inverse_name();
        }

        ValueType evaluate(EvaluationMap<FieldClass> &values) const {
            return implementation->evaluate(values);
        }

        void simplify() {
            implementation->simplify();
        }

        void expand() {
            implementation->expand();
        }

        void factor(side_type side=side_type::both) { 
            implementation->factor();
        }

        void flatten() {
            implementation->flatten();
        }

        void substitute(SubstitutionMap<FieldClass> &values) {
            implementation->substitute(values);
        }

        void derivative(const VariableClass &by_variable) {
            implementation->derivative(by_variable);
        }

        string to_string() const {
            return implementation->to_string();
        }

        void replace(const ImplementationClass &_implementation) {
            allocator_traits<ImplementationAllocator>::destroy(ia, implementation);
            ia.deallocate(implementation);
            implementation = ia.allocate(1);
            allocator_traits<ImplementationAllocator>::construct(ia, implementation, _implementation);
        }

        const ImplementationClass &get_implementation() const {
            return *implementation;
        }

        friend ostream& operator<<(ostream& os, const FunctionClass& f) {
            os << f.to_string();
            return os;
        }
};

template<class FieldClass, class AllocatorsClass=Allocators<FieldClass>, class FunctionAllocatorsClass=FunctionAllocators<FieldClass, AllocatorsClass> > Function<FieldClass, AllocatorsClass, FunctionAllocatorsClass> exp(const Entity<FieldClass, AllocatorsClass> &arg) {
    ExpImplementation<FieldClass, AllocatorsClass, FunctionAllocatorsClass> _exp(arg);
    FunctionImplementation<FieldClass, AllocatorsClass, FunctionAllocatorsClass> _exp_imp(_exp);
    return Function<FieldClass, AllocatorsClass, FunctionAllocatorsClass>(_exp_imp);
}

template<class FieldClass, class AllocatorsClass=Allocators<FieldClass>, class FunctionAllocatorsClass=FunctionAllocators<FieldClass, AllocatorsClass> > Function<FieldClass, AllocatorsClass, FunctionAllocatorsClass> log(const Entity<FieldClass, AllocatorsClass> &arg) {
    LogImplementation<FieldClass, AllocatorsClass, FunctionAllocatorsClass> _log(arg);
    FunctionImplementation<FieldClass, AllocatorsClass, FunctionAllocatorsClass> _log_imp(_log);
    return Function<FieldClass, AllocatorsClass, FunctionAllocatorsClass>(_log_imp);
}

}; // namespace function
}; // namespace symba
#endif // _SYMBA_FUNCTION_HPP_