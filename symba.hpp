#include<cmath>
#include<algorithm>
#include<variant>
#include<string>
#include<vector>
#include<iostream>
#include<map>
#include<numeric>
#include<sstream>

#include "util.hpp"

#ifndef _SYMBA_HPP_
#define _SYMBA_HPP_

namespace symba {
namespace poly {
using namespace std;
using namespace util;

template<class Field> class Allocators;
template<class Field> class Constant;
template<class Field, class Allocators> class Variable;
template<class Field, class Allocators> class Term;
template<class Field, class Allocators> class Monomial;
template<class Field, class Allocators> class Polynomial;
template<class Field, class Allocators> class Entity;

template<class FieldClass> class Allocators {
    public:
        using constant_allocator=allocator<Constant<FieldClass> >;
        using variable_allocator=allocator<Variable<FieldClass, Allocators<FieldClass> > >;
        using term_allocator=allocator<Term<FieldClass, Allocators<FieldClass> > >;
        using monomial_allocator=allocator<Monomial<FieldClass, Allocators<FieldClass> > >;
        using polynomial_allocator=allocator<Polynomial<FieldClass, Allocators<FieldClass> > >;
        using entity_allocator=allocator<Entity<FieldClass, Allocators<FieldClass> > >;
};

template<class FieldClass, class AllocatorsClass=Allocators<FieldClass> > class SubstitutionMap {
    private:
        using Field = typename FieldClass::value_type;
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

template<class FieldClass> class Constant {
    private:
        using Field = typename FieldClass::value_type;
        using ConstantClass = Constant<FieldClass>;
        Field value;
    public:
        explicit Constant(Field _value) : value(_value) { }
        string signature() const {
            ostringstream stream;
            stream << "0C_"  << value;
            return stream.str();
        }
        Field evaluate(EvaluationMap<FieldClass> &values) const {
            return value;
        }
        Field evaluate() const {
            return value;
        }
        void substitute(SubstitutionMap<FieldClass> &values) { }
        void simplify() { }
        void expand() { }
        string to_string() const {
            ostringstream stream;
            if (value > 0) stream << value;
            else stream << "("  << value << ")";
            return stream.str();
        }
};

template<class FieldClass, class AllocatorsClass=Allocators<FieldClass> > class Variable {
    private:
        using Field = typename FieldClass::value_type;
        using ConstantClass = Constant<FieldClass>;
        using VariableClass = Variable<FieldClass, AllocatorsClass>;
        string name;
    public:
        explicit Variable(const string &_name) : name(_name) { }
        const string& get_name() const {
            return name;
        }
        string signature() const {
            return "1V_" + name;
        }
        Field evaluate(EvaluationMap<FieldClass> &values) const {
            return values.get(name);
        }
        void substitute(SubstitutionMap<FieldClass> &values) { }        
        void simplify() { }
        void expand() { }
        string to_string() const {
            return name;
        }
};

template<class FieldClass, class AllocatorsClass=Allocators<FieldClass> > class Term {
    private:
        using Field = typename FieldClass::value_type;
        using PowerClass = typename FieldClass::power_type;
        using ConstantAllocator = typename AllocatorsClass::constant_allocator;
        using TermAllocator = typename AllocatorsClass::term_allocator;
        using MonomialAllocator = typename AllocatorsClass::monomial_allocator;
        using PolynomialAllocator = typename AllocatorsClass::polynomial_allocator;
        using EntityAllocator = typename AllocatorsClass::entity_allocator;
        using ConstantClass = Constant<FieldClass>;
        using VariableClass = Variable<FieldClass, AllocatorsClass>;
        using TermClass = Term<FieldClass, AllocatorsClass>;
        using EntityClass = Entity<FieldClass, AllocatorsClass>;
        using MonomialClass = Monomial<FieldClass, AllocatorsClass>;
        using PolynomialClass = Polynomial<FieldClass, AllocatorsClass>;
        EntityClass obj;
        PowerClass power;
        ConstantAllocator ca;
        TermAllocator ta;
        MonomialAllocator ma;
        PolynomialAllocator pa;
    public:
        Term(const EntityClass &_entity, const PowerClass &_power) : obj(_entity), power(_power), ca(), ta(), ma(), pa() { }
        Term(const ConstantClass &_entity, const PowerClass &_power=1) : obj(_entity), power(_power), ca(), ta(), ma(), pa() { }
        Term(const VariableClass &_entity, const PowerClass &_power=1) : obj(_entity), power(_power), ca(), ta(), ma(), pa() { }
        Term(const MonomialClass &_entity, const PowerClass &_power=1) : obj(_entity), power(_power), ca(), ta(), ma(), pa() { }
        Term(const PolynomialClass &_entity, const PowerClass &_power=1) : obj(_entity), power(_power), ca() { }

        const PowerClass &get_power() const { return power; }
        PowerClass get_power() { return power; }
        EntityClass &get_obj() { return obj; }
        const EntityClass &get_obj() const { return obj; }

        string signature() const {
            return "2T_" + obj.signature() + "^" + std::to_string(power);
        }

        Field evaluate(EvaluationMap<FieldClass> &values) const {
            return pow<Field, PowerClass>(obj.evaluate(values), power);
        }

        void substitute(SubstitutionMap<FieldClass> &values) {
            if(!is_constant()) {
                obj.substitute(values);
            }
        }

        void simplify() {
            if (is_constant()) {
                const ConstantClass &tmp = get_constant();
                ConstantClass *replacement = ca.allocate(1);
                allocator_traits<ConstantAllocator>::construct(ca, replacement, pow<Field, PowerClass>(tmp.evaluate(), power));
                obj.replace(replacement[0]);
                power = 1;
                allocator_traits<ConstantAllocator>::destroy(ca, replacement);
                ca.deallocate(replacement, 1);
            } else if (is_monomial()) {
                const MonomialClass &tmp = get_monomial();
                MonomialClass *replacement = ma.allocate(1);
                const vector<TermClass> &old_terms = tmp.get_terms();
                TermClass *new_terms = ta.allocate(old_terms.size());
                for(size_t i=0; i<old_terms.size(); ++i) {
                    if (old_terms[i].is_constant()) {
                        const ConstantClass &tmp_c = old_terms[i].get_constant();
                        ConstantClass *cobj = ca.allocate(1);
                        allocator_traits<ConstantAllocator>::construct(ca, cobj, pow<Field, PowerClass>(tmp_c.evaluate(), old_terms[i].get_power() * power));
                        allocator_traits<TermAllocator>::construct(ta, new_terms+i, cobj[0], 1);
                        allocator_traits<ConstantAllocator>::destroy(ca, cobj);
                        ca.deallocate(cobj, 1);
                    } else {
                        allocator_traits<TermAllocator>::construct(ta, new_terms+i, old_terms[i].get_obj(), old_terms[i].get_power() * power);
                    }
                }
                allocator_traits<MonomialAllocator>::construct(ma, replacement, vector<TermClass>(new_terms, new_terms+old_terms.size()));
                obj.replace(replacement[0]);
                power = 1;
                for(size_t i=0; i < old_terms.size(); ++i) {
                    allocator_traits<TermAllocator>::destroy(ta, new_terms+i);
                }
                ta.deallocate(new_terms, old_terms.size());
                allocator_traits<MonomialAllocator>::destroy(ma, replacement);
                ma.deallocate(replacement, 1);
            } else if (is_polynomial()) {
                obj.simplify();
            }
        }

        bool is_zero() const {
            return is_constant() && (get_constant().evaluate() == 0);
        }

        bool is_constant() const {
            return holds_alternative<ConstantClass>(obj.get_obj());
        }

        const ConstantClass &get_constant() const {
            return get<ConstantClass>(obj.get_obj());
        }

        bool is_variable() const {
            return holds_alternative<VariableClass>(obj.get_obj());
        }

        const VariableClass &get_variable() const {
            return get<VariableClass>(obj.get_obj());
        }

        string variable_name() const {
            if (!is_variable()) return "";
            return get_variable().get_name();
        }

        bool same_variable_signature(TermClass &_term) const {
            if (is_variable() && _term.is_variable() && (variable_name().compare(_term.variable_name()) == 0))
                return true;
            if (is_monomial() && _term.is_monomial() && (get_monomial().same_variable_signature(_term.get_monomial())))
                return true;
            if (is_polynomial() && _term.is_polynomial() && (get_polynomial().same_variable_signature(_term.get_polynomial())))
                return true;
            return false;
        }

        bool is_monomial() const {
            return holds_alternative<MonomialClass>(obj.get_obj());
        }

        const MonomialClass &get_monomial() const {
            return get<MonomialClass>(obj.get_obj());
        }

        bool is_polynomial() const {
            return holds_alternative<PolynomialClass>(obj.get_obj());
        }

        const PolynomialClass &get_polynomial() const {
            return get<PolynomialClass>(obj.get_obj());
        }

        string to_string() const {
            string ret = obj.to_string();
            if (power > 1) ret += "^" + std::to_string(power);
            return ret;
        }

        // Generate following combinations
        // of term powers and rotate the terms
        // power==7, polynomial.size()==3:
        // 7 0 0
        // 6 1 0
        // 5 2 0
        // 5 1 1
        // 4 3 0
        // 4 2 1
        // 3 3 1
        // 3 2 2
        void expand() {
            if (is_polynomial() && (power > 1)) {
                const PolynomialClass &polynomial = get_polynomial();
                const vector<MonomialClass> &monomials = polynomial.get_terms();
                vector<MonomialClass> _monomials=vector<MonomialClass>();
                Field stop_power = power / polynomial.size() + ((power % polynomial.size()) > 0);
                vector<PowerClass> powers=vector<PowerClass>(polynomial.size(),0);
                powers[0] = power;
                size_t pointer = 0; // last non-zero power pointer
                while(powers[0] >= stop_power) {
                    Field c = multinomial_coefficient<PowerClass>(power, powers.begin(), powers.begin()+pointer+1);
                    do {
                        size_t idx=0, terms_size = 1;
                        for (size_t i = 0; i < powers.size(); ++i) {
                            if (powers[i] > 0) {
                               terms_size += monomials[i].size();
                            }
                        }
                        MonomialClass *mobj = ma.allocate(1);
                        TermClass *tobj = ta.allocate(terms_size);
                        ConstantClass *cobj = ca.allocate(1);
                        allocator_traits<ConstantAllocator>::construct(ca, cobj, c);
                        allocator_traits<TermAllocator>::construct(ta, tobj, cobj[0]);
                        for (size_t i=0; i < powers.size(); ++i) {
                            if (powers[i] > 0) {
                                const vector<TermClass> &_terms = monomials[i].get_terms();
                                for(size_t j=0; j < _terms.size(); j++) {
                                    allocator_traits<TermAllocator>::construct(ta, tobj+idx+1, _terms[j].get_obj(), _terms[j].get_power()*powers[i]);
                                    idx++;
                                }
                            }
                        }
                        allocator_traits<MonomialAllocator>::construct(ma, mobj, vector<TermClass>(tobj, tobj+terms_size));
                        mobj[0].simplify();
                        _monomials.push_back(mobj[0]);
                        allocator_traits<ConstantAllocator>::destroy(ca, cobj);
                        for (size_t i=0; i < terms_size; ++i) allocator_traits<TermAllocator>::destroy(ta,tobj+i);
                        allocator_traits<MonomialAllocator>::destroy(ma, mobj);
                        ca.deallocate(cobj, 1);
                        ta.deallocate(tobj, terms_size);
                        ma.deallocate(mobj, 1);
                    } while (prev_permutation(powers.begin(), powers.end()));

                    // Generate next combination of powers
                    if ((pointer < powers.size()-1) && ((pointer==0)||((powers[pointer] > powers[pointer+1])&&(powers[pointer]>1)))) {
                        powers[pointer]--;
                        pointer++;
                        powers[pointer]++;
                    } else {
                        Field sum = powers[pointer];
                        Field level = 0;
                        powers[pointer] = 0;
                        pointer--;
                        while((pointer>0) && (sum+1 > (powers[pointer]-1)*(powers.size()-pointer-1))) {
                            sum += powers[pointer];
                            powers[pointer] = 0;
                            pointer--;
                        }
                        if ((pointer==0) && (sum+1 > (powers[pointer]-1)*(powers.size()-pointer-1))) break;
                        level = --powers[pointer];
                        ++sum;
                        while(sum > 0) {
                            if (sum > level) {
                                powers[++pointer] = level;
                                sum -= level;
                            } else {
                                powers[++pointer] = sum;
                                sum = 0;
                            }
                        }
                    }
                }
                PolynomialClass *pobj = pa.allocate(1);
                allocator_traits<PolynomialAllocator>::construct(pa, pobj, _monomials);
                obj.replace(pobj[0]);
                allocator_traits<PolynomialAllocator>::destroy(pa, pobj);
                pa.deallocate(pobj, 1);
            }
        }
};

template<class TermClass> bool compare_terms(TermClass term1, TermClass term2) { return term1.signature() < term2.signature(); }

void increment_multiindex(vector<size_t> &index, const vector<size_t> &sizes) {
    size_t i = 0;
    index[0]++;
    while(i < sizes.size()) {
        if (index[i]==sizes[i]) {
            index[i] = 0;
            index[++i]++;
        } else {
            break;
        }
    }
}

template<class FieldClass, class AllocatorsClass=Allocators<FieldClass> > class Monomial {
    private:
        using Field = typename FieldClass::value_type;
        using PowerClass = typename FieldClass::power_type;
        using ConstantAllocator = typename AllocatorsClass::constant_allocator;
        using VariableAllocator = typename AllocatorsClass::variable_allocator;
        using TermAllocator = typename AllocatorsClass::term_allocator;
        using MonomialAllocator = typename AllocatorsClass::monomial_allocator;
        using PolynomialAllocator = typename AllocatorsClass::polynomial_allocator;
        using ConstantClass = Constant<FieldClass>;
        using VariableClass = Variable<FieldClass, AllocatorsClass>;
        using TermClass = Term<FieldClass, AllocatorsClass>;
        using EntityClass = Entity<FieldClass, AllocatorsClass>;
        using MonomialClass = Monomial<FieldClass, AllocatorsClass>;
        using PolynomialClass = Polynomial<FieldClass, AllocatorsClass>;
        vector<TermClass> terms;
        ConstantAllocator ca;
        VariableAllocator va;
        TermAllocator ta;
        MonomialAllocator ma;
        PolynomialAllocator pa;
    public:
        explicit Monomial(const vector<TermClass> &_terms) : terms(_terms), ca(), va(), ta(), ma(), pa() { sort(terms.begin(), terms.end(), compare_terms<TermClass>); }

        size_t size() const { return terms.size(); }

        const TermClass& first_term() const { return terms[0]; }

        const TermClass& get_term(size_t i) const { return terms[i]; }

        const vector<TermClass> &get_terms() const { return terms; }

        bool is_zero() const {
            return (terms.size() == 0) || ((terms.size() == 1) && (terms[0].is_zero()));
        }

        string variable_signature() const {
            string ret = "3MV_";
            for (auto term=terms.begin(); term!=terms.end(); ++term) {
                if (term->is_variable()) ret += term->signature();
            }
            return ret;
        }

        string signature() const {
            string ret = "3M_";
            for (auto term=terms.begin(); term!=terms.end(); ++term) ret += term->signature();
            return ret;
        }

        Field evaluate(EvaluationMap<FieldClass> &values) const {
            Field ret(1);
            for (auto term=terms.begin(); term!=terms.end(); ++term) ret *= term->evaluate(values);
            return ret;
        }

        void substitute(SubstitutionMap<FieldClass> &values) {
            for(auto term=terms.begin(); term!=terms.end(); ++term) term->substitute(values);
        }

        void simplify() {
            Field c(1);
            PowerClass p = 0;
            size_t constants = 0;
            auto term = terms.begin();
            #pragma omp parallel for
            for (; term!=terms.end(); ++term) term->simplify();
            term = terms.begin();
            while(term!=terms.end()) {
                if (term->is_monomial()) {
                    const MonomialClass &_monomial = term->get_monomial();
                    const vector<TermClass> &_terms = _monomial.get_terms();
                    term = terms.insert(term+1, _terms.begin(), _terms.end());
                    term = terms.erase(term-1);
                } else {
                    ++term;
                }
            }
            sort(terms.begin(), terms.end(), compare_terms<TermClass>);
            term = terms.begin();
            while (term->is_constant() && (term != terms.end())) {
                // Find all constants
                ConstantClass tmp = term->get_constant();
                c *= tmp.evaluate();
                ++term;
                ++constants;
            }
            if (constants > 1) {
                ConstantClass *cobj = ca.allocate(1);
                allocator_traits<ConstantAllocator>::construct(ca, cobj, c);
                TermClass *cterm = ta.allocate(1);
                allocator_traits<TermAllocator>::construct(ta, cterm, cobj[0]);
                terms.erase(terms.begin(), term-1);
                *(terms.begin()) = cterm[0];
                term = terms.begin() + 1;
                allocator_traits<TermAllocator>::destroy(ta, cterm);
                ta.deallocate(cterm, 1);
                allocator_traits<ConstantAllocator>::destroy(ca, cobj);
                ca.deallocate(cobj, 1);
            } else {
                ++term;
            }
            if ((constants > 0) && (terms[0].get_constant().evaluate()==0)) {
                // Remove everything is constant == 0
                terms.erase(terms.begin()+1, terms.end());
            } else {
                // Combine same variables
                size_t processed_variables = 0, same_variables = 0;
                auto first_occurence = term;
                while(term != terms.end()) {
                    while ((term != terms.end()) && ((processed_variables==0)||(term->same_variable_signature(*first_occurence)))) {
                        if (processed_variables==0) {
                            first_occurence = term;
                        }
                        p += term->get_power();
                        ++term;
                        ++same_variables;
                        ++processed_variables;
                    }
                    if (same_variables > 1) {
                        if (first_occurence->is_variable()) {
                            VariableClass *vobj = va.allocate(1);
                            allocator_traits<VariableAllocator>::construct(va, vobj, first_occurence->get_variable().get_name());
                            TermClass *vterm = ta.allocate(1);
                            allocator_traits<TermAllocator>::construct(ta, vterm, vobj[0], p);
                            term = terms.erase(first_occurence, term-1);
                            *(term) = vterm[0];
                            ++term;
                            allocator_traits<TermAllocator>::destroy(ta, vterm);
                            ta.deallocate(vterm, 1);
                            allocator_traits<VariableAllocator>::destroy(va, vobj);
                            va.deallocate(vobj, 1);
                        } else if (first_occurence->is_monomial()) {
                            MonomialClass *mobj = ma.allocate(1);
                            allocator_traits<MonomialAllocator>::construct(ma, mobj, first_occurence->get_monomial().get_terms());
                            TermClass *mterm = ta.allocate(1);
                            allocator_traits<TermAllocator>::construct(ta, mterm, mobj[0], p);
                            term = terms.erase(first_occurence, term-1);
                            *(term) = mterm[0];
                            ++term;
                            allocator_traits<TermAllocator>::destroy(ta, mterm);
                            ta.deallocate(mterm, 1);
                            allocator_traits<MonomialAllocator>::destroy(ma, mobj);
                            ma.deallocate(mobj, 1);
                        } else if (first_occurence->is_polynomial()) {
                            PolynomialClass *pobj = pa.allocate(1);
                            allocator_traits<PolynomialAllocator>::construct(pa, pobj, first_occurence->get_polynomial().get_terms());
                            TermClass *pterm = ta.allocate(1);
                            allocator_traits<TermAllocator>::construct(ta, pterm, pobj[0], p);
                            term = terms.erase(first_occurence, term-1);
                            *(term) = pterm[0];
                            ++term;
                            allocator_traits<TermAllocator>::destroy(ta, pterm);
                            ta.deallocate(pterm, 1);
                            allocator_traits<PolynomialAllocator>::destroy(pa, pobj);
                            pa.deallocate(pobj, 1);
                        }
                    }
                    same_variables = 1;
                    first_occurence = term;
                    if (term != terms.end()) {
                        p = term->get_power();
                        ++term;
                        ++processed_variables;
                    }
                }
            }
        }

        void expand() {
            auto term = terms.begin();
            size_t min_idx=0, total_size = 1, i = 0;
            vector<size_t> polynomial_idxs = {}, polynomial_counters = {},
                polynomial_sizes = {};
            while(term != terms.end()) {
                if (term->is_polynomial()) {
                    if (min_idx==0) min_idx=i;
                    term->expand();
                    polynomial_idxs.push_back(i);
                    polynomial_counters.push_back(0);
                    polynomial_sizes.push_back(term->get_polynomial().size());
                    total_size *= term->get_polynomial().size();
                }
                ++i;
                ++term;
            }
            if (total_size > 1) {
                MonomialClass *mobj = ma.allocate(total_size);
                vector<TermClass> first_terms;
                if (min_idx > 0) {
                    // If we have at least one constant or variable
                    first_terms = vector<TermClass>(terms.begin(),terms.begin()+min_idx);
                } else {
                    // Otherwise, initialize first monomial with constant = 1
                    ConstantClass *cobj = ca.allocate(1);
                    allocator_traits<ConstantAllocator>::construct(ca, cobj, 1);
                    Term<FieldClass, AllocatorsClass> *tobj = ta.allocate(1);
                    allocator_traits<TermAllocator>::construct(ta, tobj, cobj[0]);
                    first_terms = {tobj[0]};
                    allocator_traits<TermAllocator>::destroy(ta, tobj);
                    ta.deallocate(tobj, 1);
                    allocator_traits<ConstantAllocator>::destroy(ca, cobj);
                    ca.deallocate(cobj, 1);
                }
                i = 0;
                while(i < total_size) {
                    vector<TermClass> cur_terms = first_terms;
                    for(size_t j=0; j < polynomial_idxs.size(); j++) {
                        const PolynomialClass &polynomial = terms[polynomial_idxs[j]].get_polynomial();
                        const vector<MonomialClass> &monomials = polynomial.get_terms();
                        const MonomialClass &monomial = monomials[polynomial_counters[j]];
                        const vector<TermClass> &new_terms = monomial.get_terms();
                        cur_terms.insert(cur_terms.end(), new_terms.begin(), new_terms.end());
                    }
                    allocator_traits<MonomialAllocator>::construct(ma, mobj+i, cur_terms);
                    increment_multiindex(polynomial_counters, polynomial_sizes);
                    i++;
                }
                PolynomialClass *pobj = pa.allocate(1);
                vector<MonomialClass> _monomials = vector<MonomialClass>(mobj, mobj+total_size);
                allocator_traits<PolynomialAllocator>::construct(pa, pobj, _monomials);
                terms.clear();
                terms.push_back(pobj[0]);
                for(i=0; i<total_size; i++) allocator_traits<MonomialAllocator>::destroy(ma, mobj+i);
                ma.deallocate(mobj, total_size);
                allocator_traits<PolynomialAllocator>::destroy(pa, pobj);
                pa.deallocate(pobj, 1);
            }
        }

        string to_string() const {
            string ret = "";
            bool first = true;
            for (auto term=terms.begin(); term!=terms.end(); ++term) {
                if (!first) ret += " ";
                else first = false;
                ret += term->to_string();
            }
            return ret;
        }

        Field get_constant_value() {
            if (terms[0].is_constant()) {
                return terms[0].get_constant().evaluate();
            } else {
                return Field(1);
            }
        }

        bool has_constant() {
            return terms[0].is_constant();
        }

        bool same_variable_signature(const MonomialClass &_monomial) const {
            return variable_signature().compare(_monomial.variable_signature()) == 0;
        }

        void replace_constant_value(Field c) {
            if (has_constant()) {
                ConstantClass *cobj = ca.allocate(1);
                allocator_traits<ConstantAllocator>::construct(ca, cobj, c);
                TermClass *cterm = ta.allocate(1);
                allocator_traits<TermAllocator>::construct(ta, cterm, cobj[0]);
                *(terms.begin()) = cterm[0];
                allocator_traits<TermAllocator>::destroy(ta, cterm);
                ta.deallocate(cterm, 1);
                allocator_traits<ConstantAllocator>::destroy(ca, cobj);
                ca.deallocate(cobj, 1);
            } else {
                ConstantClass *cobj = ca.allocate(1);
                allocator_traits<ConstantAllocator>::construct(ca, cobj, c);
                TermClass *cterm = ta.allocate(1);
                allocator_traits<TermAllocator>::construct(ta, cterm, cobj[0]);
                terms.insert(terms.begin(), cterm[0]);
                allocator_traits<TermAllocator>::destroy(ta, cterm);
                ta.deallocate(cterm, 1);
                allocator_traits<ConstantAllocator>::destroy(ca, cobj);
                ca.deallocate(cobj, 1);
            }
        }
};

template<class MonomialClass> bool compare_monomials(MonomialClass term1, MonomialClass term2) { if(term1.variable_signature() == term2.variable_signature()) return term1.signature() < term2.signature(); else return term1.variable_signature() < term2.variable_signature(); }

template<class FieldClass, class AllocatorsClass=Allocators<FieldClass> > class Polynomial {
    private:
        using Field = typename FieldClass::value_type;
        using PowerClass = typename FieldClass::power_type;
        using ConstantAllocator = typename AllocatorsClass::constant_allocator;
        using VariableAllocator = typename AllocatorsClass::variable_allocator;
        using TermAllocator = typename AllocatorsClass::term_allocator;
        using MonomialAllocator = typename AllocatorsClass::monomial_allocator;
        using PolynomialAllocator = typename AllocatorsClass::polynomial_allocator;
        using ConstantClass = Constant<FieldClass>;
        using VariableClass = Variable<FieldClass, AllocatorsClass>;
        using TermClass = Term<FieldClass, AllocatorsClass>;
        using EntityClass = Entity<FieldClass, AllocatorsClass>;
        using MonomialClass = Monomial<FieldClass, AllocatorsClass>;
        using PolynomialClass = Polynomial<FieldClass, AllocatorsClass>;
        vector<MonomialClass> terms;
        ConstantAllocator ca;
        VariableAllocator va;
        TermAllocator ta;
        MonomialAllocator ma;
        PolynomialAllocator pa;
    public:
        explicit Polynomial(const vector<MonomialClass> &_terms) : terms(_terms), ca(), va(), ta(), ma(), pa() { sort(terms.begin(), terms.end(), compare_monomials<MonomialClass>); }

        size_t size() const { return terms.size(); }

        const MonomialClass& first_term() const { return terms[0]; }

        const MonomialClass& get_term(size_t i) const { return terms[i]; }

        const vector<MonomialClass> &get_terms() const { return terms; }

        bool is_zero() const {
            return (terms.size() == 0) || ((terms.size()==1) && (terms[0].is_zero()));
        }

        string signature() const {
            string ret = "4P_";
            for (auto term=terms.begin(); term!=terms.end(); ++term) ret += term->signature();
            return ret;
        }

        string variable_signature() const {
            return signature();
        }

        Field evaluate(EvaluationMap<FieldClass> &values) const {
            Field ret(0);
            for (auto term=terms.begin(); term!=terms.end(); ++term) ret += term->evaluate(values);
            return ret;
        }

        void substitute(SubstitutionMap<FieldClass> &values) {
            for(auto term=terms.begin(); term!=terms.end(); ++term) {
                term->substitute(values);
            }
        }

        void simplify() {
            auto term=terms.begin();
            #pragma omp parallel for
            for(; term!=terms.end(); ++term) term->simplify();
            term=terms.begin();
            while(term!=terms.end()) {
                if (term->is_zero() && (terms.size() > 1)) {
                    term = terms.erase(term);
                } else {
                    ++term;
                }
            }
            term = terms.begin();
            while(term!=terms.end()) {
                // Erase single term monomials with polynomial in them
                if ((term->size()==1) && (term->first_term().is_polynomial())) {
                    const PolynomialClass &_polynomial = term->first_term().get_polynomial();
                    const vector<MonomialClass> &_terms = _polynomial.get_terms();
                    term = terms.insert(term+1, _terms.begin(), _terms.end());
                    term = terms.erase(term-1);
                } else {
                    ++term;
                }
            }
            sort(terms.begin(), terms.end(), compare_monomials<MonomialClass>);
            if (terms.size() > 1) {
                Field c(0);
                size_t processed_monomials=0, same_monomials=0;
                auto first_occurence = term;
                term = terms.begin();
                while(term != terms.end()) {
                    while ((term != terms.end()) && ((processed_monomials==0)||(term->same_variable_signature(*first_occurence)))) {
                        // Capture same variables
                        if (processed_monomials==0) {
                            first_occurence = term;
                        }
                        c += term->get_constant_value();
                        ++term;
                        ++same_monomials;
                        ++processed_monomials;
                    }
                    if (same_monomials > 1) {
                        // Replace multiple same variables with one
                        term = terms.erase(first_occurence, term-1);
                        term->replace_constant_value(c);
                        ++term;
                    }
                    same_monomials = 1;
                    first_occurence = term;
                    if (term != terms.end()) {
                        c = term->get_constant_value();
                        ++term;
                        ++processed_monomials;
                    }
                }
            }
        }

        void expand() {
            #pragma omp parallel for 
            for (auto term=terms.begin(); term!=terms.end(); ++term) term->expand();
            simplify();
        }

        string to_string() const {
            string ret = "(";
            bool first = true;
            for (auto term=terms.begin(); term!=terms.end(); ++term) {
                if (!first) ret += "+"; 
                else first=false;
                ret += term->to_string();
            }
            ret += ")";
            return ret;
        }

        bool same_variable_signature(const PolynomialClass &_polynomial) const {
            return variable_signature().compare(_polynomial.variable_signature()) == 0;
        }        
};


template<class FieldClass, class AllocatorsClass=Allocators<FieldClass> > class Entity {
    private:
        using Field = typename FieldClass::value_type;
        using ConstantClass = Constant<FieldClass>;
        using VariableClass = Variable<FieldClass, AllocatorsClass>;
        using MonomialClass = Monomial<FieldClass, AllocatorsClass>;
        using PolynomialClass = Polynomial<FieldClass, AllocatorsClass>;
        using EntityClass = Entity<FieldClass, AllocatorsClass>;
        using variant_class = variant<ConstantClass, VariableClass, MonomialClass, PolynomialClass>;
        variant_class obj;
    public:
        Entity() : obj(ConstantClass(0)) {}
        explicit Entity(const variant_class &_obj) : obj(_obj) { }
        Entity(const ConstantClass &_c) : obj(_c) { }
        Entity(const VariableClass &_v) : obj(_v) { }
        Entity(const MonomialClass &_p) : obj(_p) { }
        Entity(const PolynomialClass &_p) : obj(_p) { }
        string signature() const {
            return visit([](auto &&arg) -> string { return arg.signature(); }, obj);
        }
        Field evaluate(EvaluationMap<FieldClass> &values) const {
            return visit([&values](auto &&arg) -> Field { return arg.evaluate(values); }, obj);
        }
        void substitute(SubstitutionMap<FieldClass> &values) {
            if (holds_alternative<ConstantClass>(obj)) {
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
        string to_string() const {
            return visit([](auto &&arg) -> string { return arg.to_string(); }, obj);
        }
        void replace(variant_class _obj) {
            obj = _obj;
        }
        const variant_class &get_obj() const {
            return obj;
        }
};

};
};
#endif