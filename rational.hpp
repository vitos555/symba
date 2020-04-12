
#ifndef _SYMBA_RATIONAL_HPP_
#define _SYMBA_RATIONAL_HPP_

#include "core.hpp"

namespace symba {
namespace rational {
using namespace std;
using namespace util;
using namespace polynomial;
using namespace core;

template<class FieldClass, class AllocatorsClass=Allocators<FieldClass> > class Rational {
    private:
        using ValueType = typename FieldClass::value_type;
        using CoefficientType = typename FieldClass::coefficient_type;        
        using PowerType = typename FieldClass::power_type;
        using CoefficientAllocator = typename AllocatorsClass::coefficient_allocator;
        using ConstantAllocator = typename AllocatorsClass::constant_allocator;
        using VariableAllocator = typename AllocatorsClass::variable_allocator;
        using TermAllocator = typename AllocatorsClass::term_allocator;
        using MonomialAllocator = typename AllocatorsClass::monomial_allocator;
        using PolynomialAllocator = typename AllocatorsClass::polynomial_allocator;
        using RationalAllocator = typename AllocatorsClass::rational_allocator;
        using CoefficientClass = Coefficient<FieldClass>;
        using ConstantClass = Constant<FieldClass>;
        using VariableClass = Variable<FieldClass, AllocatorsClass>;
        using TermClass = Term<FieldClass, AllocatorsClass>;
        using EntityClass = Entity<FieldClass, AllocatorsClass>;
        using MonomialClass = Monomial<FieldClass, AllocatorsClass>;
        using PolynomialClass = Polynomial<FieldClass, AllocatorsClass>;
        using RationalClass = Rational<FieldClass, AllocatorsClass>;
        PolynomialClass numerator;
        PolynomialClass denonimator;
        CoefficientAllocator oa;
        ConstantAllocator ca;
        VariableAllocator va;
        TermAllocator ta;
        MonomialAllocator ma;
        PolynomialAllocator pa;
        RationalAllocator ra;
    public:
        Rational(const PolynomialClass &_numerator, const PolynomialClass &_denonimator) : numerator(_numerator), denonimator(_denonimator), oa(), ca(), va(), ta(), ma(), pa(), ra() { }

        const PolynomialClass& get_numerator() const { return numerator; }

        const PolynomialClass& get_denominator() const { return denonimator; }

        bool is_zero() const {
            return numerator.is_zero();
        }

        bool is_infinity() const {
            return denonimator.is_zero();
        }

        string signature() const {
            ostringstream stream;
            stream << "6RD^" << setw(FieldClass::power_type_fill::value) << setfill('0') << FieldClass::power_type_max::value - denonimator.get_power() << denonimator.signature();
            stream << "/N^" << setw(FieldClass::power_type_fill::value) << setfill('0') << FieldClass::power_type_max::value - numerator.get_power() << numerator.signature();
            return stream.str();
        }

        ValueType evaluate(EvaluationMap<FieldClass> &values) const {
            return numerator.evaluate(values)/denonimator.evaluate(values);
        }

        void simplify() {
            numerator.simplify();
            denonimator.simplify();
        }

        void expand() {
            numerator.expand();
            denonimator.expand();
        }

        void factor(side_type side=side_type::both) {
            numerator.factor(side);
            denonimator.factor(side);
        }

        void flatten() {
            numerator.flatten();
            denonimator.flatten();
        }

        void substitute(SubstitutionMap<FieldClass> &values) {
            numerator.substitute(values);
            denonimator.substitute(values);
        }

        void cancel() {
            factor();
            if ((denonimator.get_size()==1) && (numerator.get_size()==1)) {
                if constexpr(FieldClass::is_coefficient_multiplication_commutative::value && FieldClass::is_value_multiplication_commutative::value) {
                    const MonomialClass &d_monomial = denonimator.first_term(),
                        &n_monomial = numerator.first_term();
                    const vector<MonomialClass> &d_terms = d_monomial.get_terms(),
                        &n_terms = n_monomial.get_terms();
                    auto d_term = d_terms.begin();
                    auto n_term = n_terms.begin();
                    vector<MonomialClass> new_d_terms=vector<MonomialClass>(),
                        new_n_terms=vector<MonomialClass>();
                    while((d_term != d_terms.end()) || (n_term != n_terms.end())) {
                        if ((d_term != d_terms.end())&&(n_term != n_terms.end())) {
                            int cmp = compare_monomial_variables<MonomialClass>(*d_term, *n_term);
                            if (cmp > 0) {
                                new_d_terms.push_back(*d_term);
                                ++d_term;
                            } else if (cmp < 0) {
                                new_n_terms.push_back(*n_term);
                                ++n_term;
                            } else {
                                ++d_term;
                                ++n_term;
                            }
                        } else if (d_term != d_terms.end()) {
                            new_d_terms.push_back(*d_term);
                            ++d_term;
                        } else if (n_term != n_terms.end()) {
                            new_n_terms.push_back(*n_term);
                            ++n_term;
                        }
                    }
                    PolynomialClass *pobj = pa.allocate(2);
                    MonomialClass *mobj = ma.allocate(2);
                    allocator_traits<MonomialAllocator>::construct(ma, mobj, new_n_terms);
                    allocator_traits<MonomialAllocator>::construct(ma, mobj+1, new_d_terms);
                    allocator_traits<PolynomialAllocator>::construct(pa, pobj, vector<MonomialClass>(mobj[0]));
                    allocator_traits<PolynomialAllocator>::construct(pa, pobj+1, vector<MonomialClass>(mobj[1]));
                    numerator = pobj[0];
                    denonimator = pobj[1];
                    allocator_traits<PolynomialAllocator>::destroy(pa, pobj);
                    allocator_traits<PolynomialAllocator>::destroy(pa, pobj+1);
                    allocator_traits<MonomialAllocator>::destroy(ma, mobj);
                    allocator_traits<MonomialAllocator>::destroy(ma, mobj+1);
                    ma.deallocate(mobj, 2);
                    pa.deallocate(pobj, 2);
                }
            }
        }

        string to_string() const {
            return "(" + numerator.to_string() + ")/(" + denonimator.to_string() + ")";
        }        
};


}; // namespace rational
}; // namespace symba
#endif // _SYMBA_RATIONAL_HPP_