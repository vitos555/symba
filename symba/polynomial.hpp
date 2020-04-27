
#ifndef _SYMBA_POLYNOMIAL_HPP_
#define _SYMBA_POLYNOMIAL_HPP_

#include "core.hpp"

namespace symba {
namespace polynomial {
using namespace std;
using namespace util;
using namespace core;

template<class FieldClass> class Coefficient {
    private:
        typedef typename FieldClass::coefficient_type CoefficientType;
        typedef Coefficient<FieldClass> CoefficientClass;
        CoefficientType value;
    public:
        explicit Coefficient(CoefficientType _value) : value(_value) { }
        string signature() const {
            ostringstream stream;
            stream << "0O_"  << value;
            return stream.str();
        }
        CoefficientType evaluate(EvaluationMap<FieldClass> &values) const {
            return value;
        }
        CoefficientType evaluate() const {
            return value;
        }
        void substitute(SubstitutionMap<FieldClass> &values) { }
        void simplify() { }
        void expand() { }
        void flatten() { }
        void factor(side_type side=side_type::both) { }
        string to_string() const {
            ostringstream stream;
            if (value > 0) stream << value;
            else stream << "("  << value << ")";
            return stream.str();
        }
};

template<class FieldClass> class Constant {
    private:
        typedef typename FieldClass::value_type ValueType;
        typedef Constant<FieldClass> ConstantClass;
        ValueType value;
    public:
        explicit Constant(ValueType _value) : value(_value) { }
        string signature() const {
            ostringstream stream;
            stream << "1C_"  << value;
            return stream.str();
        }
        ValueType evaluate(EvaluationMap<FieldClass> &values) const {
            return value;
        }
        ValueType evaluate() const {
            return value;
        }
        void substitute(SubstitutionMap<FieldClass> &values) { }
        void simplify() { }
        void expand() { }
        void flatten() { }
        void factor(side_type side=side_type::both) { }
        string to_string() const {
            ostringstream stream;
            if (value > 0) stream << value;
            else stream << "("  << value << ")";
            return stream.str();
        }
};

template<class FieldClass, class AllocatorsClass=Allocators<FieldClass> > class Variable {
    private:
        typedef typename FieldClass::value_type ValueType;
        typedef Constant<FieldClass> ConstantClass;
        typedef Variable<FieldClass, AllocatorsClass> VariableClass;
        string name;
    public:
        explicit Variable(const string &_name) : name(_name) { }
        const string& get_name() const {
            return name;
        }
        string signature() const {
            return "2V_" + name;
        }
        ValueType evaluate(EvaluationMap<FieldClass> &values) const {
        
            return values.get(name);
        }
        void substitute(SubstitutionMap<FieldClass> &values) { }        
        void simplify() { }
        void expand() { }
        void flatten() { }
        void factor(side_type side=side_type::both) { }
        string to_string() const {
            return name;
        }
};

template<class FieldClass, class AllocatorsClass=Allocators<FieldClass> > class Term {
    private:
        typedef typename FieldClass::value_type ValueType;
        typedef typename FieldClass::coefficient_type CoefficientType;
        typedef typename FieldClass::power_type PowerType;
        typedef typename AllocatorsClass::coefficient_allocator CoefficientAllocator;
        typedef typename AllocatorsClass::constant_allocator ConstantAllocator;
        typedef typename AllocatorsClass::term_allocator TermAllocator;
        typedef typename AllocatorsClass::monomial_allocator MonomialAllocator;
        typedef typename AllocatorsClass::polynomial_allocator PolynomialAllocator;
        typedef typename AllocatorsClass::entity_allocator EntityAllocator;
        typedef Coefficient<FieldClass> CoefficientClass;
        typedef Constant<FieldClass> ConstantClass;
        typedef Variable<FieldClass, AllocatorsClass> VariableClass;
        typedef Term<FieldClass, AllocatorsClass> TermClass;
        typedef Entity<FieldClass, AllocatorsClass> EntityClass;
        typedef Monomial<FieldClass, AllocatorsClass> MonomialClass;
        typedef Polynomial<FieldClass, AllocatorsClass> PolynomialClass;
        EntityClass obj;
        PowerType power;
        CoefficientAllocator oa;
        ConstantAllocator ca;
        TermAllocator ta;
        MonomialAllocator ma;
        PolynomialAllocator pa;
    public:
        explicit Term(const EntityClass &_entity, const PowerType &_power=1) : obj(_entity), power(_power), oa(), ca(), ta(), ma(), pa() { }
        Term(const CoefficientClass &_entity, const PowerType &_power=1) : obj(_entity), power(_power), oa(), ca(), ta(), ma(), pa() { }
        Term(const ConstantClass &_entity, const PowerType &_power=1) : obj(_entity), power(_power), oa(), ca(), ta(), ma(), pa() { }
        Term(const VariableClass &_entity, const PowerType &_power=1) : obj(_entity), power(_power), oa(), ca(), ta(), ma(), pa() { }
        Term(const MonomialClass &_entity, const PowerType &_power=1) : obj(_entity), power(_power), oa(), ca(), ta(), ma(), pa() { }
        Term(const PolynomialClass &_entity, const PowerType &_power=1) : obj(_entity), power(_power), oa(), ca(), ta(), ma(), pa() { }

        const PowerType &get_power() const { return power; }
        PowerType get_power() { return power; }
        const EntityClass &get_obj() const { return obj; }
        EntityClass &get_obj() { return obj; }

        string signature() const {
            ostringstream stream;
            stream << "4T_" << obj.signature() << "^" << setw(FieldClass::power_type_fill::value) << setfill('0') << FieldClass::power_type_max::value - get_power();
            return stream.str();
        }

        string value_signature() const {
            if (is_coefficient()||is_constant()) return "";
            ostringstream stream;            
            stream << "4TV_" << obj.signature();
            return stream.str();
        }

        ValueType evaluate(EvaluationMap<FieldClass> &values = {}) const {
            ValueType ret = 0;
            if (is_coefficient())
                ret = pow<CoefficientType, PowerType>(obj.evaluate(values), power);
            else if (is_constant())
                ret = pow<ValueType, PowerType>(obj.evaluate(values), power);
            else
                ret = pow<ValueType, PowerType>(obj.evaluate(values), power);
            return ret;
        }

        void substitute(SubstitutionMap<FieldClass> &values) {
            if(!is_constant() && !is_coefficient()) {
                obj.substitute(values);
            }
        }

        void simplify() {
            if (is_coefficient()) {
                const CoefficientClass &tmp = get_coefficient();
                CoefficientClass *replacement = oa.allocate(1);
                allocator_traits<CoefficientAllocator>::construct(oa, replacement, pow<CoefficientType, PowerType>(tmp.evaluate(), power));
                obj.replace(replacement[0]);
                power = 1;
                allocator_traits<CoefficientAllocator>::destroy(oa, replacement);
                oa.deallocate(replacement, 1);
            } else if (is_constant()) {
                const ConstantClass &tmp = get_constant();
                ConstantClass *replacement = ca.allocate(1);
                allocator_traits<ConstantAllocator>::construct(ca, replacement, pow<ValueType, PowerType>(tmp.evaluate(), power));
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
                    if (old_terms[i].is_coefficient()) {
                        const CoefficientClass &tmp_o = old_terms[i].get_coefficient();
                        CoefficientClass *oobj = oa.allocate(1);
                        allocator_traits<CoefficientAllocator>::construct(oa, oobj, pow<CoefficientType, PowerType>(tmp_o.evaluate(), old_terms[i].get_power() * power));
                        allocator_traits<TermAllocator>::construct(ta, new_terms+i, oobj[0], 1);
                        allocator_traits<CoefficientAllocator>::destroy(oa, oobj);
                        oa.deallocate(oobj, 1);
                    } else if (old_terms[i].is_constant()) {
                        const ConstantClass &tmp_c = old_terms[i].get_constant();
                        ConstantClass *cobj = ca.allocate(1);
                        allocator_traits<ConstantAllocator>::construct(ca, cobj, pow<ValueType, PowerType>(tmp_c.evaluate(), old_terms[i].get_power() * power));
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
            if (is_coefficient() && (get_coefficient().evaluate() == 0))
                return true;
            if (is_constant() && (get_constant().evaluate() == 0))
                return true;
            return false;
        }

        bool is_one() const {
            if (is_coefficient() && (get_coefficient().evaluate() == 1))
                return true;
            if (is_constant() && (get_constant().evaluate() == 1))
                return true;
            return false;
        }

        bool is_coefficient() const {
            return holds_alternative<CoefficientClass>(obj.get_obj());
        }

        const CoefficientClass &get_coefficient() const {
            return get<CoefficientClass>(obj.get_obj());
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

        bool same_value_signature(TermClass &_term) const {
            bool ret = false;
            if (is_variable() && _term.is_variable()) {
                if (variable_name().compare(_term.variable_name()) == 0)
                    ret = true;
            } else if (is_monomial() && _term.is_monomial()) {
                if (get_monomial().same_value_signature(_term.get_monomial()))
                    ret = true;
            } else if (obj.signature().compare(_term.get_obj().signature()) == 0) {
                ret = true;
            }
            return ret;
        }

        bool same_variable_signature(TermClass &_term) const {
            bool ret = false;
            if (is_variable() && _term.is_variable()) {
                if (variable_name().compare(_term.variable_name()) == 0)
                    ret = true;
            }
            else if (is_monomial() && _term.is_monomial()) {
                if (get_monomial().same_variable_signature(_term.get_monomial()))
                    ret = true;
            } else if (obj.signature().compare(_term.get_obj().signature()) == 0) {
                ret = true;
            }
            return ret;
        }

        bool same_power_signature(TermClass &_term) const {
            bool ret = false;
            if (is_variable() && _term.is_variable()) {
                if ((variable_name().compare(_term.variable_name()) == 0) && (get_power() == _term.get_power()))
                    ret = true;
            } else if (is_monomial() && _term.is_monomial()) {
                if (get_monomial().same_power_signature(_term.get_monomial()))
                    ret = true;
            } else if (obj.power_signature().compare(_term.get_obj().power_signature()) == 0) {
                ret = true;
            }
            return ret;
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
            if (!is_coefficient() && !is_constant() && !is_variable()) {
                ret = "(" + ret + ")";
            }
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
                vector<MonomialClass> _monomials = {};
                PowerType stop_power = power / polynomial.size() + ((power % polynomial.size()) > 0);
                vector<PowerType> powers = vector<PowerType>(polynomial.size(), 0);
                powers[0] = power;
                size_t pointer = 0; // last non-zero power pointer
                while(powers[0] >= stop_power) {
                    CoefficientType c = multinomial_coefficient<PowerType>(power, powers.begin(), powers.begin()+pointer+1);
                    do {
                        size_t idx=0, terms_size = 1;
                        for (size_t i = 0; i < powers.size(); ++i) {
                            if (powers[i] > 0) {
                               terms_size += monomials[i].size();
                            }
                        }
                        MonomialClass *mobj = ma.allocate(1);
                        TermClass *tobj = ta.allocate(terms_size);
                        CoefficientClass *oobj = oa.allocate(1);
                        allocator_traits<CoefficientAllocator>::construct(oa, oobj, c);
                        allocator_traits<TermAllocator>::construct(ta, tobj, oobj[0]);
                        allocator_traits<CoefficientAllocator>::destroy(oa, oobj);
                        oa.deallocate(oobj, 1);
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
                        for (size_t i=0; i < terms_size; ++i) allocator_traits<TermAllocator>::destroy(ta, tobj+i);
                        allocator_traits<MonomialAllocator>::destroy(ma, mobj);
                        ta.deallocate(tobj, terms_size);
                        ma.deallocate(mobj, 1);
                    } while (prev_permutation(powers.begin(), powers.end()));

                    // Generate next combination of powers
                    if ((pointer < powers.size()-1) && ((pointer==0)||((powers[pointer] > powers[pointer+1])&&(powers[pointer]>1)))) {
                        powers[pointer]--;
                        pointer++;
                        powers[pointer]++;
                    } else {
                        PowerType sum = powers[pointer];
                        PowerType level = 0;
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

        void flatten() {
            obj.flatten();
        }

        void factor(side_type side=side_type::both) {
            obj.factor(side);
        }
};

template<class TermClass> bool compare_terms(TermClass term1, TermClass term2) { return term1.signature().compare(term2.signature()) < 0; }

void increment_multi_index(vector<size_t> &index, const vector<size_t> &sizes) {
    if (index.size() > 0) {
        size_t i = 0;
        index[0]++;
        while(i < sizes.size()-1) {
            if (index[i]==sizes[i]) {
                index[i] = 0;
                index[++i]++;
            } else {
                break;
            }
        }
    }
}

template<class FieldClass, class AllocatorsClass=Allocators<FieldClass> > class Monomial {
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
        typedef Coefficient<FieldClass> CoefficientClass;
        typedef Constant<FieldClass> ConstantClass;
        typedef Variable<FieldClass, AllocatorsClass> VariableClass;
        typedef Term<FieldClass, AllocatorsClass> TermClass;
        typedef Entity<FieldClass, AllocatorsClass> EntityClass;
        typedef Monomial<FieldClass, AllocatorsClass> MonomialClass;
        typedef Polynomial<FieldClass, AllocatorsClass> PolynomialClass;
        bool sorted;
        bool combined_similar;
        vector<TermClass> terms;
        CoefficientAllocator oa;
        ConstantAllocator ca;
        VariableAllocator va;
        TermAllocator ta;
        MonomialAllocator ma;
        PolynomialAllocator pa;
    public:
        explicit Monomial(const vector<TermClass> &_terms) : terms(_terms), oa(), ca(), va(), ta(), ma(), pa(), sorted(false), combined_similar(false) {
            sort();
        }

        size_t size() const { return terms.size(); }

        const TermClass& first_term() const { return terms[0]; }

        const TermClass& get_term(size_t i) const { return terms[i]; }

        const vector<TermClass> &get_terms() const { return terms; }

        bool is_zero() const {
            bool ret = false;
            if (terms.size() == 0) {
                ret = true;
            } else {
                if (sorted && terms[0].is_zero()) {
                    ret = true;
                } else {
                    auto term = terms.begin();
                    while(term != terms.end()) {
                        if (term->is_zero()) {
                            ret = true;
                            break;
                        }
                        ++term;
                    }
                }
            }
            return ret;
        }

        bool is_one() const {
            ValueType p(1);
            bool has_non_const = false, ret;
            auto term = terms.begin();
            while(term != terms.end()) {
                if (term->is_coefficient()) {
                    p *= term->get_coefficient().evaluate();
                } else if (term->is_constant()) {
                    p *= term->get_constant().evaluate();
                } else {
                    has_non_const = true;
                    break;
                }
                ++term;
            }
            ret = !has_non_const && (p==1);
            return ret;
        }

        string variable_signature() const {
            ostringstream stream;
            stream << "4MVR_";
            for (auto term=terms.begin(); term!=terms.end(); ++term) 
                if (term->is_variable()) stream << term->signature();
            return stream.str();
        }

        string value_signature() const {
            ostringstream stream;
            stream << "4MVL_";
            for (auto term=terms.begin(); term!=terms.end(); ++term) 
                if (term->is_variable()||term->is_constant()) stream << term->signature();
            return stream.str();
        }

        string power_signature() const {
            ostringstream stream;
            stream << "4MP^"  << setw(FieldClass::power_type_fill::value) << setfill('0') << FieldClass::power_type_max::value - get_power() << "_";
            for (auto term=terms.begin(); term!=terms.end(); ++term) 
                if (term->is_variable()) stream << term->signature();
            return stream.str();
        }

        PowerType get_power() const {
            PowerType ret = 0;
            for (auto term=terms.begin(); term!=terms.end(); ++term) {
                if (!term->is_constant() && !term->is_coefficient()) ret += term->get_power();
            }
            return ret;
        }

        string signature() const {
            ostringstream stream;
            stream << "4M^"  << setw(FieldClass::power_type_fill::value) << setfill('0') << FieldClass::power_type_max::value - get_power() << "_";
            for (auto term=terms.begin(); term!=terms.end(); ++term) stream << term->signature();
            return stream.str();
        }

        bool same_variable_signature(const MonomialClass &_monomial) const {
            return variable_signature().compare(_monomial.variable_signature()) == 0;
        }

        bool same_value_signature(const MonomialClass &_monomial) const {
            return value_signature().compare(_monomial.value_signature()) == 0;
        }

        bool same_power_signature(const MonomialClass &_monomial) const {
            return power_signature().compare(_monomial.power_signature()) == 0;
        }

        ValueType evaluate(EvaluationMap<FieldClass> &values) const {
            ValueType ret(1);
            for (auto term=terms.begin(); term!=terms.end(); ++term) ret *= term->evaluate(values);
            return ret;
        }

        void substitute(SubstitutionMap<FieldClass> &values) {
            for(auto term=terms.begin(); term!=terms.end(); ++term) term->substitute(values);
            sorted = false;
            combined_similar = false;
        }

        inline void sort() {
            if constexpr(FieldClass::is_coefficient_multiplication_commutative::value && FieldClass::is_value_multiplication_commutative::value) {
                if (!sorted) {
                    std::sort(terms.begin(), terms.end(), compare_terms<TermClass>);
                    sorted = true;
                }
            }
        }

        void flatten() {
            auto term = terms.begin();
            while(term!=terms.end()) {
                if (term->is_monomial()) {
                    // Remove monomials in the monomial
                    const MonomialClass &_monomial = term->get_monomial();
                    const vector<TermClass> &_terms = _monomial.get_terms();
                    term = terms.insert(term+1, _terms.begin(), _terms.end());
                    term = terms.erase(term-1);
                } else if ((term->is_polynomial())&&(term->get_polynomial().size()==1)) {
                    // Remove single term polynomials
                    const PolynomialClass &_polynomial = term->get_polynomial();
                    const vector<MonomialClass> &_polynomial_terms = _polynomial.get_terms();
                    const MonomialClass &_monomial = _polynomial_terms[0];
                    const vector<TermClass> &_terms = _monomial.get_terms();                    
                    term = terms.insert(term+1, _terms.begin(), _terms.end());
                    term = terms.erase(term-1);
                } else {
                    term->flatten();
                    ++term;
                }
            }
            sorted = false;
            combined_similar = false;
            sort();
        }

        void combine_similar() {
            if (!combined_similar) {
                CoefficientType o(1);
                ValueType c(1);
                PowerType p = 0;
                size_t coefficients = 0;
                auto term = terms.begin();
                if constexpr(FieldClass::is_coefficient_multiplication_commutative::value && FieldClass::is_value_multiplication_commutative::value) {
                    size_t constants = 0;
                    sort();
                    term = terms.begin();
                    while (term->is_coefficient() && (term != terms.end())) {
                        // Find all coefficients
                        const CoefficientClass &tmp = term->get_coefficient();
                        o *= tmp.evaluate();
                        ++term;
                        ++coefficients;
                    }
                    if constexpr(!is_same_v<CoefficientType, ValueType>) {
                        if (coefficients > 1) {
                            CoefficientClass *oobj = oa.allocate(1);
                            allocator_traits<CoefficientAllocator>::construct(oa, oobj, o);
                            TermClass *oterm = ta.allocate(1);
                            allocator_traits<TermAllocator>::construct(ta, oterm, oobj[0]);
                            terms.erase(terms.begin(), term-1);
                            *(terms.begin()) = oterm[0];
                            term = terms.begin() + 1;
                            allocator_traits<TermAllocator>::destroy(ta, oterm);
                            ta.deallocate(oterm, 1);
                            allocator_traits<CoefficientAllocator>::destroy(oa, oobj);
                            oa.deallocate(oobj, 1);
                        } else {
                            ++term;
                        }
                    } else {
                        if (coefficients > 0) {
                            term = terms.erase(terms.begin(), terms.begin() + coefficients);                        
                        }
                    }
                    while (term->is_constant() && (term != terms.end())) {
                        // Find all constants
                        const ConstantClass &tmp = term->get_constant();
                        c *= tmp.evaluate();
                        ++term;
                        ++constants;
                    }
                    if constexpr(!is_same_v<CoefficientType, ValueType>) {
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
                    } else {
                        if (constants > 0) {
                            term = terms.erase(terms.begin(), terms.begin() + constants);                        
                        }

                    }
                    if constexpr(is_same_v<CoefficientType, ValueType>) {
                        // If CoefficientType is the same as ValueType insert single value
                        o *= c;
                        if ((o != 1) || (size() == 0)) {
                            CoefficientClass *oobj = oa.allocate(1);
                            allocator_traits<CoefficientAllocator>::construct(oa, oobj, o);
                            TermClass *oterm = ta.allocate(1);
                            allocator_traits<TermAllocator>::construct(ta, oterm, oobj[0]);
                            term = terms.insert(terms.begin(), oterm[0]);
                            ++term;
                            allocator_traits<TermAllocator>::destroy(ta, oterm);
                            ta.deallocate(oterm, 1);
                            allocator_traits<CoefficientAllocator>::destroy(oa, oobj);
                            oa.deallocate(oobj, 1);
                        }
                    }
                }
                else if constexpr(FieldClass::is_coefficient_multiplication_commutative::value) {
                    // Remove coefficients first
                    term = terms.begin();
                    while(term!=terms.end()) {
                        if (term->is_coefficient()) {
                            o *= term->get_coefficient().evaluate();
                            term = terms.erase(term);
                            ++coefficients;
                        } else {
                            ++term;
                        }
                    }
                    if ((coefficients > 0) && (o != 1)) {
                        // Insert single coefficient in front
                        CoefficientClass *oobj = oa.allocate(1);
                        allocator_traits<CoefficientAllocator>::construct(oa, oobj, c);
                        TermClass *oterm = ta.allocate(1);
                        allocator_traits<TermAllocator>::construct(ta, oterm, oobj[0]);
                        term = terms.insert(terms.begin(), oterm[0]);
                        ++term;
                        allocator_traits<TermAllocator>::destroy(ta, oterm);
                        ta.deallocate(oterm, 1);
                        allocator_traits<CoefficientAllocator>::destroy(oa, oobj);
                        oa.deallocate(oobj, 1);
                    }
                }

                if ((c==0)||(o==0)) {
                    // Remove everything if there is 0 coefficient or 0 constant
                    terms.clear();
                    ConstantClass *cobj = ca.allocate(1);
                    allocator_traits<ConstantAllocator>::construct(ca, cobj, 0);
                    TermClass *cterm = ta.allocate(1);
                    allocator_traits<TermAllocator>::construct(ta, cterm, cobj[0]);
                    terms.erase(terms.begin(), term-1);
                    terms.push_back(cterm[0]);
                    term = terms.begin() + 1;
                    allocator_traits<TermAllocator>::destroy(ta, cterm);
                    ta.deallocate(cterm, 1);
                    allocator_traits<ConstantAllocator>::destroy(ca, cobj);
                    ca.deallocate(cobj, 1);
                } else {
                    // Combine same objects
                    size_t processed_objects = 0, same_objects = 0;
                    auto first_occurence = term;
                    while(term != terms.end()) {
                        while ((term != terms.end()) && ((processed_objects==0)||(term->same_value_signature(*first_occurence)))) {
                            if (processed_objects==0) {
                                first_occurence = term;
                            }
                            p += term->get_power();
                            ++term;
                            ++same_objects;
                            ++processed_objects;
                        }
                        if (same_objects > 1) {
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
                            } else if (first_occurence->is_constant()) {
                                ConstantClass *cobj = ca.allocate(1);
                                allocator_traits<ConstantAllocator>::construct(ca, cobj, first_occurence->get_constant().evaluate());
                                TermClass *cterm = ta.allocate(1);
                                allocator_traits<TermAllocator>::construct(ta, cterm, cobj[0], p);
                                term = terms.erase(first_occurence, term-1);
                                *(term) = cterm[0];
                                ++term;
                                allocator_traits<TermAllocator>::destroy(ta, cterm);
                                ta.deallocate(cterm, 1);
                                allocator_traits<ConstantAllocator>::destroy(ca, cobj);
                                ca.deallocate(cobj, 1);
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
                        same_objects = 1;
                        first_occurence = term;
                        if (term != terms.end()) {
                            p = term->get_power();
                            ++term;
                            ++processed_objects;
                        }
                    }
                }
                combined_similar = true;
            }
        }        

        void simplify() {
            auto term = terms.begin();
            for (; term!=terms.end(); ++term) term->simplify();
            flatten();
            combine_similar();
        }

        void expand() {
            auto term = terms.begin();
            size_t min_idx=0, total_size = 1, i = 0;
            vector<size_t> polynomial_counters = {},
                polynomial_sizes = {};
            while(term != terms.end()) {
                if (term->is_polynomial()) {
                    if (min_idx==0) min_idx=i;
                    term->expand();
                    polynomial_counters.push_back(0);
                    polynomial_sizes.push_back(term->get_polynomial().size());
                    total_size *= term->get_polynomial().size();
                } else {
                    polynomial_sizes.push_back(1);
                    polynomial_counters.push_back(0);
                }
                ++i;
                ++term;
            }
            if (total_size > 1) {
                MonomialClass *mobj = ma.allocate(total_size);
                vector<TermClass> first_terms = {};
                if (min_idx > 0) {
                    // If we have at least one constant or variable
                    first_terms.insert(first_terms.end(), terms.begin(), terms.begin()+min_idx);
                } else {
                    // Otherwise, initialize first monomial with constant = 1
                    ConstantClass *cobj = ca.allocate(1);
                    allocator_traits<ConstantAllocator>::construct(ca, cobj, 1);
                    Term<FieldClass, AllocatorsClass> *tobj = ta.allocate(1);
                    allocator_traits<TermAllocator>::construct(ta, tobj, cobj[0]);
                    first_terms.push_back(tobj[0]);
                    allocator_traits<TermAllocator>::destroy(ta, tobj);
                    ta.deallocate(tobj, 1);
                    allocator_traits<ConstantAllocator>::destroy(ca, cobj);
                    ca.deallocate(cobj, 1);
                }
                i = 0;
                while(i < total_size) {
                    vector<TermClass> cur_terms = {};
                    for(size_t j=0; j < size(); j++) {
                        if (terms[j].is_polynomial()) {
                            const PolynomialClass &polynomial = terms[j].get_polynomial();
                            const vector<MonomialClass> &monomials = polynomial.get_terms();
                            const MonomialClass &monomial = monomials[polynomial_counters[j]];
                            const vector<TermClass> &new_terms = monomial.get_terms();
                            cur_terms.insert(cur_terms.end(), new_terms.begin(), new_terms.end());
                        } else {
                            cur_terms.push_back(terms[j]);
                        }
                    }
                    allocator_traits<MonomialAllocator>::construct(ma, mobj+i, cur_terms);
                    increment_multi_index(polynomial_counters, polynomial_sizes);
                    i++;
                }
                PolynomialClass *pobj = pa.allocate(1);
                vector<MonomialClass> _monomials = vector<MonomialClass>(mobj, mobj+total_size);
                allocator_traits<PolynomialAllocator>::construct(pa, pobj, _monomials);
                terms.clear();
                terms.push_back(pobj[0]);
                for(i=0; i < total_size; i++) allocator_traits<MonomialAllocator>::destroy(ma, mobj+i);
                ma.deallocate(mobj, total_size);
                allocator_traits<PolynomialAllocator>::destroy(pa, pobj);
                pa.deallocate(pobj, 1);
            }
            sorted = false;
            combined_similar = false;
            sort();
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

        bool has_coefficient() {
            return terms[0].is_coefficient();
        }

        ValueType get_coefficient_value() {
            if (terms[0].is_coefficient()) {
                return terms[0].get_coefficient().evaluate();
            } else {
                return CoefficientType(1);
            }
        }

        void replace_coefficient_value(CoefficientType c) {
            if (has_coefficient()) {
                CoefficientClass *oobj = oa.allocate(1);
                allocator_traits<CoefficientAllocator>::construct(oa, oobj, c);
                TermClass *oterm = ta.allocate(1);
                allocator_traits<TermAllocator>::construct(ta, oterm, oobj[0]);
                *(terms.begin()) = oterm[0];
                allocator_traits<TermAllocator>::destroy(ta, oterm);
                ta.deallocate(oterm, 1);
                allocator_traits<CoefficientAllocator>::destroy(oa, oobj);
                oa.deallocate(oobj, 1);
            } else {
                CoefficientClass *oobj = oa.allocate(1);
                allocator_traits<CoefficientAllocator>::construct(oa, oobj, c);
                TermClass *oterm = ta.allocate(1);
                allocator_traits<TermAllocator>::construct(ta, oterm, oobj[0]);
                terms.insert(terms.begin(), oterm[0]);
                allocator_traits<TermAllocator>::destroy(ta, oterm);
                ta.deallocate(oterm, 1);
                allocator_traits<CoefficientAllocator>::destroy(oa, oobj);
                oa.deallocate(oobj, 1);
            }
        }

        void factor(side_type side=side_type::both) {
            for (auto term=terms.begin(); term!=terms.end(); ++term) term->factor(side);
            sorted = false;
            combined_similar = false;
            flatten();
        }

        // Should be doing following
        // 16x^2z / 6xy = 8 x z, 3 y
        // 16x^2y / 2xy = 8 x, 1
        pair<MonomialClass, MonomialClass> divide(const MonomialClass &divisor, const side_type &side=side_type::both) {
            vector<TermClass> d_terms = vector<TermClass>(), n_terms = vector<TermClass>();
            const vector<TermClass> &_terms = divisor.get_terms();
            auto self_term = terms.begin();
            auto divisor_term = _terms.begin();
            while((self_term != terms.end()) && (divisor_term != _terms.end())){
                if ((self_term->is_coefficient()) || (divisor_term->is_coefficient())) {
                    if (self_term->is_coefficient() && divisor_term->is_coefficient()) {
                        CoefficientType o = FieldClass::gcd(self_term->get_coefficient().evaluate(), divisor_term->get_coefficient().evaluate());
                        if (o > 1) {
                            CoefficientType o1 = self_term->get_coefficient().evaluate() / o;
                            CoefficientType o2 = divisor_term->get_coefficient().evaluate() / o;
                            CoefficientClass *oobj = oa.allocate(2);
                            allocator_traits<CoefficientAllocator>::construct(oa, oobj, o1);
                            allocator_traits<CoefficientAllocator>::construct(oa, oobj+1, o2);
                            n_terms.push_back(oobj[0]);
                            d_terms.push_back(oobj[1]);
                            allocator_traits<CoefficientAllocator>::destroy(oa, oobj);
                            allocator_traits<CoefficientAllocator>::destroy(oa, oobj+1);
                            oa.deallocate(oobj, 2);
                        }
                        ++self_term;
                        ++divisor_term;
                    } else if (self_term->is_coefficient()) {
                        n_terms.push_back(*self_term);
                        ++self_term;
                    } else if (divisor_term->is_coefficient()) {
                        d_terms.push_back(*divisor_term);
                        ++divisor_term;
                    }
                } else {
                    int cmp = self_term->value_signature().compare(divisor_term->value_signature());
                    if (cmp > 0) {
                        n_terms.push_back(*self_term);
                        ++self_term;
                    } else if (cmp < 0) {
                        d_terms.push_back(*divisor_term);
                        ++divisor_term;
                    } else {
                        if (self_term->get_power() > divisor_term->get_power()) {
                            TermClass *tobj = ta.allocate(1);
                            allocator_traits<TermAllocator>::construct(ta, tobj, self_term->get_obj(), self_term->get_power() - divisor_term->get_power());
                            n_terms.push_back(tobj[0]);
                            allocator_traits<TermAllocator>::destroy(ta, tobj);
                            ta.deallocate(tobj, 1);
                        } else if (self_term->get_power() < divisor_term->get_power()) {
                            TermClass *tobj = ta.allocate(1);
                            allocator_traits<TermAllocator>::construct(ta, tobj, self_term->get_obj(), divisor_term->get_power() - self_term->get_power());
                            d_terms.push_back(tobj[0]);
                            allocator_traits<TermAllocator>::destroy(ta, tobj);
                            ta.deallocate(tobj, 1);
                        }
                        ++divisor_term;
                        ++self_term;
                    }
                }
            }
            if (n_terms.size() == 0) {
                CoefficientClass *oobj = oa.allocate(1);
                allocator_traits<CoefficientAllocator>::construct(oa, oobj, 1);
                n_terms.push_back(oobj[0]);
                allocator_traits<CoefficientAllocator>::destroy(oa, oobj);
                oa.deallocate(oobj, 1);
            }
            if (d_terms.size() == 0) {
                CoefficientClass *oobj = oa.allocate(1);
                allocator_traits<CoefficientAllocator>::construct(oa, oobj, 1);
                d_terms.push_back(oobj[0]);
                allocator_traits<CoefficientAllocator>::destroy(oa, oobj);
                oa.deallocate(oobj, 1);
            }
            return pair<MonomialClass, MonomialClass>(MonomialClass(n_terms), MonomialClass(d_terms));
        }

        // Should be doing following
        // 16x^2z & 6xy = 2 x
        void intersect(const MonomialClass &other, const side_type &side=side_type::both) {
            const vector<TermClass> &_terms = other.get_terms();
            auto self_term = terms.begin();
            auto other_term = _terms.begin();
            vector<TermClass> intersection = vector<TermClass>({});
            while((self_term != terms.end()) && (other_term != _terms.end())) {
                if ((self_term->is_coefficient()) || (other_term->is_coefficient())) {
                    if (self_term->is_coefficient() && other_term->is_coefficient()) {
                        CoefficientType o = FieldClass::gcd(self_term->get_coefficient().evaluate(), other_term->get_coefficient().evaluate());
                        if (o > 1) {
                            CoefficientClass *oobj = oa.allocate(1);
                            allocator_traits<CoefficientAllocator>::construct(oa, oobj, o);
                            intersection.push_back(oobj[0]);
                            allocator_traits<CoefficientAllocator>::destroy(oa, oobj);
                            oa.deallocate(oobj, 1);
                        }
                        ++self_term;
                        ++other_term;
                    } else if (self_term->is_coefficient()) {
                        ++self_term;
                    } else if (other_term->is_coefficient()) {
                        ++other_term;
                    }
                } else {
                    int cmp = self_term->value_signature().compare(other_term->value_signature());
                    if (cmp > 0) {
                        ++self_term;
                    } else if (cmp < 0) {
                        ++other_term;
                    } else {
                        if (self_term->get_power() > other_term->get_power()) {
                            intersection.push_back(*other_term);
                        } else {
                            intersection.push_back(*self_term);
                        }
                        ++other_term;
                        ++self_term;
                    }
                }
            }
            terms = intersection;
        }

        void multiply(const EntityClass &other, const side_type &side=side_type::right) {
            TermClass *tobj = ta.allocate(1);
            allocator_traits<TermAllocator>::construct(ta, tobj, other);
            if ((side == side_type::right)||(side == side_type::both)) {
                terms.push_back(tobj[0]);
            }
            if ((side == side_type::left)||(side == side_type::both)) {
                terms.insert(terms.begin(), tobj[0]);
            }
            allocator_traits<TermAllocator>::destroy(ta, tobj);
            ta.deallocate(tobj, 1);
        }
};

template<class MonomialClass> bool compare_monomials(MonomialClass term1, MonomialClass term2) { if(term1.power_signature().compare(term2.power_signature())==0) return term1.signature().compare(term2.signature()) < 0; else return term1.power_signature().compare(term2.power_signature()) < 0; }
template<class MonomialClass> int compare_monomial_variables(MonomialClass term1, MonomialClass term2) { if(term1.variable_signature().compare(term2.variable_signature())==0) return term1.signature().compare(term2.signature()); else return term1.variable_signature.compare(term2.variable_signature()); }

template<class FieldClass, class AllocatorsClass=Allocators<FieldClass> > class Polynomial {
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
        typedef Coefficient<FieldClass> CoefficientClass;
        typedef Constant<FieldClass> ConstantClass;
        typedef Variable<FieldClass, AllocatorsClass> VariableClass;
        typedef Term<FieldClass, AllocatorsClass> TermClass;
        typedef Entity<FieldClass, AllocatorsClass> EntityClass;
        typedef Monomial<FieldClass, AllocatorsClass> MonomialClass;
        typedef Polynomial<FieldClass, AllocatorsClass> PolynomialClass;
        bool combined_similar;
        bool sorted;
        vector<MonomialClass> terms;
        CoefficientAllocator oa;
        ConstantAllocator ca;
        VariableAllocator va;
        TermAllocator ta;
        MonomialAllocator ma;
        PolynomialAllocator pa;
    public:
        explicit Polynomial(const vector<MonomialClass> &_terms) : terms(_terms), oa(), ca(), va(), ta(), ma(), pa(), combined_similar(false), sorted(false) {
            sort();
        }

        size_t size() const { return terms.size(); }

        const MonomialClass& first_term() const { return terms[0]; }

        const MonomialClass& get_term(size_t i) const { return terms[i]; }

        const vector<MonomialClass> &get_terms() const { return terms; }

        bool is_zero() const {
            return (terms.size() == 0) || ((terms.size()==1) && (terms[0].is_zero()));
        }

        PowerType get_power() const {
            PowerType ret = 0;
            for (auto term=terms.begin(); term!=terms.end(); ++term) {
                if (term->get_power() > ret) ret = term->get_power();
            }
            return ret;
        }

        string signature() const {
            ostringstream stream;
            stream << "5P^" << setw(FieldClass::power_type_fill::value) << setfill('0') << FieldClass::power_type_max::value - get_power() << "_";
            for (auto term=terms.begin(); term!=terms.end(); ++term) stream << term->signature();
            return stream.str();
        }

        string power_signature() const {
            ostringstream stream;
            // Factorized comparison
            // if (size() > 1) {
            //     vector<CoefficientType> coefs = vector<CoefficientType>();
            //     for (size_t i=0; i<terms.size(); ++i) coefs[i] = terms[i].get_coefficient_value();
            //     pair<CoefficientType, vector<CoefficientType> > factors = FieldClass::factor(coefs);
            //     if (get<0>(factors) > 1) {
            //         stream << "5PP^" << setw(FieldClass::power_type_fill::value) << setfill('0') << FieldClass::power_type_max::value - get_power() << "_";
            //         for (size_t i=0; i<terms.size(); ++i) {
            //             stream << "O_" << get<1>(factors)[i];
            //             stream << terms[i].value_signature();
            //         }
            //     }
            // } else {
                stream << signature();
            // }
            return stream.str();
        }

        bool same_power_signature(const PolynomialClass &_polynomial) const {
            return power_signature().compare(_polynomial.power_signature()) == 0;
        }

        ValueType evaluate(EvaluationMap<FieldClass> &values) const {
            ValueType ret(0);
            for (auto term=terms.begin(); term!=terms.end(); ++term) ret += term->evaluate(values);
            return ret;
        }

        void substitute(SubstitutionMap<FieldClass> &values) {
            combined_similar = false;
            sorted = false;
            for(auto term=terms.begin(); term!=terms.end(); ++term) {
                term->substitute(values);
            }
        }

        inline void sort() {
            if constexpr(FieldClass::is_value_addition_commutative::value) {
                if (!sorted) {
                    std::sort(terms.begin(), terms.end(), compare_monomials<MonomialClass>);
                    sorted = true;
                }
            }            
        }

        void flatten() {
            auto term = terms.begin();
            while(term!=terms.end()) {
                // Erase single term monomials with polynomial in them
                if ((term->size()==1) && (term->first_term().is_polynomial())) {
                    const PolynomialClass &_polynomial = term->first_term().get_polynomial();
                    const vector<MonomialClass> &_terms = _polynomial.get_terms();
                    term = terms.insert(term+1, _terms.begin(), _terms.end());
                    term = terms.erase(term-1);
                } else {
                    term->flatten();
                    ++term;
                }
            }            
        }

        inline void remove_zeroes() {
            auto term=terms.begin();
            while(term!=terms.end()) {
                if (term->is_zero() && (terms.size() > 1)) {
                    term = terms.erase(term);
                } else {
                    ++term;
                }
            }

        }

        inline void combine_similar() {
            if (!combined_similar) {
                sort();
                if (terms.size() > 1) {
                    CoefficientType c(0);
                    size_t processed_monomials=0, same_monomials=0;
                    auto term = terms.begin();
                    auto first_occurence = term;
                    while(term != terms.end()) {
                        while ((term != terms.end()) && ((processed_monomials==0)||(term->same_power_signature(*first_occurence)))) {
                            // Capture same monomials
                            if (processed_monomials==0) {
                                first_occurence = term;
                            }
                            c += term->get_coefficient_value();
                            ++term;
                            ++same_monomials;
                            ++processed_monomials;
                        }
                        if (same_monomials > 1) {
                            // Replace multiple monimials with one
                            term = terms.erase(first_occurence, term-1);
                            term->replace_coefficient_value(c);
                            ++term;
                        }
                        same_monomials = 1;
                        first_occurence = term;
                        if (term != terms.end()) {
                            c = term->get_coefficient_value();
                            ++term;
                            ++processed_monomials;
                        }
                    }
                }
                combined_similar = true;
            }
        }

        void simplify() {
            auto term=terms.begin();
            for(; term!=terms.end(); ++term) term->simplify();
            remove_zeroes();
            flatten();
            combine_similar();
        }

        void expand() {
            for (auto term=terms.begin(); term!=terms.end(); ++term) term->expand();
            combined_similar = false;
            sorted = false;
            simplify();
        }

        string to_string() const {
            string ret = "";
            bool first = true;
            for (auto term=terms.begin(); term!=terms.end(); ++term) {
                if (!first) ret += "+";
                else first = false;
                ret += term->to_string();
            }
            return ret;
        }

        void factor_coefficients() {
            if (size() > 1) {
                remove_zeroes();
                combine_similar();
                vector<TermClass> _terms = vector<TermClass>();
                vector<CoefficientType> coefs = vector<CoefficientType>(size());
                for (size_t i=0; i < size(); ++i) coefs[i] = terms[i].get_coefficient_value();
                pair<CoefficientType, vector<CoefficientType> > factored_coefs = FieldClass::factor(coefs);
                if (get<0>(factored_coefs) != 1) {
                    PolynomialClass *pobj = pa.allocate(1);
                    CoefficientClass *oobj = oa.allocate(1);
                    vector<MonomialClass> new_terms = terms;
                    for(size_t i=0; i < size(); ++i) {
                        new_terms[i].replace_coefficient_value(get<1>(factored_coefs)[i]);
                    }
                    allocator_traits<PolynomialAllocator>::construct(pa, pobj, new_terms);
                    allocator_traits<CoefficientAllocator>::construct(oa, oobj, get<0>(factored_coefs));
                    _terms.push_back(oobj[0]);
                    _terms.push_back(pobj[0]);
                    allocator_traits<CoefficientAllocator>::destroy(oa, oobj);
                    allocator_traits<PolynomialAllocator>::destroy(pa, pobj);
                    oa.deallocate(oobj, 1);
                    pa.deallocate(pobj, 1);
                }
                if (_terms.size() > 1) {
                    MonomialClass *mobj = ma.allocate(1);
                    allocator_traits<MonomialAllocator>::construct(ma, mobj, _terms);
                    terms.clear();
                    terms.push_back(mobj[0]);
                    allocator_traits<MonomialAllocator>::destroy(ma, mobj);
                    ma.deallocate(mobj, 1);
                }
            }
            combined_similar = false;
            sorted = false;
        }

        void factor(const side_type &side=side_type::both) {
            if (size() > 1) {
                remove_zeroes();
                combine_similar();
                CoefficientType o(1);
                auto cur_term = terms.begin();
                MonomialClass intersection = *cur_term;
                ++cur_term;
                while(cur_term != terms.end()) {
                    intersection.intersect(*cur_term, side);
                    if (intersection.is_one()) {
                        break;
                    }
                    ++cur_term;
                }
                if (!intersection.is_one()) {
                    vector<MonomialClass> new_terms = terms;
                    for(size_t i=0; i < size(); ++i) {
                        pair<MonomialClass, MonomialClass> res = new_terms[i].divide(intersection, side);
                        new_terms[i] = get<0>(res);
                    }
                    MonomialClass *mobj = ma.allocate(1);
                    PolynomialClass *pobj = pa.allocate(2);
                    allocator_traits<PolynomialAllocator>::construct(pa, pobj, vector<MonomialClass>({intersection}));
                    allocator_traits<PolynomialAllocator>::construct(pa, pobj+1, new_terms);
                    allocator_traits<MonomialAllocator>::construct(ma, mobj, vector<TermClass>(pobj, pobj+2));
                    terms.clear();
                    terms.push_back(mobj[0]);
                    allocator_traits<PolynomialAllocator>::destroy(pa, pobj);
                    allocator_traits<PolynomialAllocator>::destroy(pa, pobj + 1);
                    allocator_traits<MonomialAllocator>::destroy(ma, mobj);
                    pa.deallocate(pobj, 2);
                    ma.deallocate(mobj, 1);
                }
                flatten();
            } else {
                terms[0].factor(side);
            }
            combined_similar = false;
            sorted = false;
        }

        // Should be doing following
        // 4 x^2 + 4 x + 1 / (x+1) = 4x, 1
        // 4 x y + 4 y + 4 x / (x+y) = 4x, 4y
        pair<PolynomialClass, PolynomialClass> divide(const PolynomialClass &divisor) {
            sort();

        }
        void add(const EntityClass &other, const side_type &side=side_type::right) {
            TermClass *tobj = ta.allocate(1);
            allocator_traits<TermAllocator>::construct(ta, tobj, other);
            MonomialClass *mobj = ma.allocate(1);
            allocator_traits<MonomialAllocator>::construct(ma, mobj, vector<TermClass>({tobj[0]}));
            if ((side == side_type::right)||(side == side_type::both)) {
                terms.push_back(mobj[0]);
            }
            if ((side == side_type::left)||(side == side_type::both)) {
                terms.insert(terms.begin(), mobj[0]);
            }
            allocator_traits<TermAllocator>::destroy(ta, tobj);
            ta.deallocate(tobj, 1);
            allocator_traits<MonomialAllocator>::destroy(ma, mobj);
            ma.deallocate(mobj, 1);
        }        
};

}; // namespace polynomial
}; // namespace symba
#endif // _SYMBA_POLYNOMIAL_HPP_