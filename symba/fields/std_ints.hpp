
#ifndef _SYMBA_FIELDS_STD_INTS_HPP_
#define _SYMBA_FIELDS_STD_INTS_HPP_

#include<cmath>
#include<algorithm>
#include<vector>
#include<numeric>
#include<sstream>

#include "util.hpp"

namespace symba {
namespace fields {
using namespace std;

template<class Field> inline pair<Field, vector<Field> > factor_int(const vector<Field> &values, util::side_type side) {
    vector<Field> ret = values;
    Field common_factor = values[0];
    auto _value = ret.begin();
    for(;_value != ret.end(); ++_value) 
        if (common_factor > *_value) common_factor = *_value;
    _value = ret.begin();
    while((common_factor > 1)&&(_value!=ret.end())) {
        common_factor = gcd<Field, Field>(*_value, common_factor);
        ++_value;
    }
    for(_value=ret.begin(); _value != ret.end(); ++_value) (*_value) /= common_factor;
    return pair<Field, vector<Field> >(common_factor, ret);
}

class IntField {
    public:
        typedef int value_type;
        typedef unsigned int power_type;
        typedef int coefficient_type;
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
            enum { value = 1 };
        };
        inline static pair<value_type, vector<value_type> > factor(const vector<value_type> &values, const util::side_type &side=util::side_type::both) {
            return factor_int<value_type>(values, side);
        }
        inline static value_type gcd(const value_type& val1, const value_type& val2) {
            return std::gcd<value_type, value_type>(val1, val2);
        }
        inline static value_type abs(const value_type& val) {
            return std::abs(val);
        }
        inline static value_type exp(const value_type& val) {
            return std::exp(val);
        }
        inline static value_type log(const value_type& val) {
            return std::log(val);
        }
};

class LongIntField {
    public:
        typedef long int value_type;
        typedef unsigned long int power_type;
        typedef long int coefficient_type;
        struct power_type_max {
            enum { value = 10000 };
        };
        struct power_type_fill {
            enum { value = 16 };
        };
        struct is_value_addition_commutative {
            enum { value = 1 };
        };
        struct is_coefficient_multiplication_commutative {
            enum { value = 1 };
        };
        struct is_value_multiplication_commutative {
            enum { value = 1 };
        };
        inline static pair<value_type, vector<value_type> > factor(const vector<value_type> &values, const util::side_type &side=util::side_type::both) {
            return factor_int<value_type>(values, side);
        }
        inline static value_type gcd(const value_type& val1, const value_type& val2) {
            return std::gcd<value_type, value_type>(val1, val2);
        }
        inline static value_type abs(const value_type& val) {
            return std::abs(val);
        }
        inline static value_type exp(const value_type& val) {
            return std::exp(val);
        }
        inline static value_type log(const value_type& val) {
            return std::log(val);
        }
};

}; // namespace fields
}; // namespace symba
#endif // _SYMBA_FIELDS_STD_INTS_HPP_
