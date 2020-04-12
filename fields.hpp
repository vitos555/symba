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

#ifndef _FIELDS_HPP_
#define _FIELDS_HPP_

namespace symba {
namespace fields {
using namespace std;

template<class Field> class rationals {
    private:
        using rationals_type = rationals<Field>;
        Field numerator;
        Field denominator;
    public:
        rationals() : numerator(1), denominator(0) { }
        rationals(Field _numerator, Field _denominator) : numerator(_numerator), denominator(_denominator) {
            if ((numerator == 0)&&(denominator==0)) { }
            else if ((numerator == 0)&&(denominator != 1)&&(denominator != -1)) {
                if (denominator > 0) denominator = 1; else denominator = -1;
            } else if ((denominator == 0)&&(numerator != 1) && (numerator != -1)) {
                if (numerator > 0) numerator = 1; else numerator = -1;
            }
        }
        rationals(Field _numerator) : numerator(_numerator), denominator(1) { }
        const Field &get_numerator() const { return numerator; }
        const Field &get_denominator() const { return denominator; }
        rationals_type operator+=(const Field &c) {
            numerator += c*denominator;
            return *this;
        }
        rationals_type operator-=(const Field &c) {
            numerator -= c*denominator;
            return *this;
        }
        rationals_type operator*=(const Field &c) {
            if (c == 0) {
                numerator = 0;
                if ((denominator > 0) && (denominator != 1)) denominator = 1;
                if ((denominator < 0) && (denominator != -1)) denominator = -1;
            } else {
                Field divider = gcd<Field, Field>(denominator, c);
                numerator *= c/divider;
                denominator /= divider;
            }
            return *this;
        }
        rationals_type operator/=(const Field &c) {
            if (c == 0) {
                denominator = 0;
                if ((numerator > 0) && (numerator != 1)) numerator = 1;
                if ((numerator < 0) && (numerator != -1)) numerator = -1;
            } else {
                denominator *= c;
                if (denominator < 0) {
                    numerator *= -1;
                    denominator *= -1;
                }
            }
            return *this;
        }
        // rationals_type operator-(const rationals_type &a) {
        //     numerator = (-1)*numerator;
        //     return *this;
        // }
        rationals_type operator+=(const rationals_type &b) {
            if (denominator*numerator*b.denominator*b.numerator==0) {
                if (denominator == 0) {

                } else if (b.denominator == 0) {
                    denominator = 0;
                    if ((numerator > 0) && (numerator != 1)) numerator = 1;
                    if ((numerator < 0) && (numerator != -1)) numerator = -1;                
                } else if (numerator == 0) {
                    numerator = b.numerator;
                    denominator = b.denominator;
                }
            } else {
                Field _denominator = gcd<Field, Field>(denominator, b.denominator);
                if (_denominator>1) {
                    numerator = denominator/_denominator*(b.numerator) + b.denominator/_denominator*(numerator);
                    denominator = denominator / _denominator * b.denominator;
                } else {
                    numerator = b.denominator*numerator + denominator*b.numerator;
                    denominator *= b.denominator;
                }
                Field divider = gcd<Field, Field>(denominator, numerator);
                numerator /= divider;
                denominator /= divider;
            }
            return *this;
        }
        rationals_type operator-=(const rationals_type &b) {
            if (denominator*numerator*b.denominator*b.numerator==0) {
                if (denominator == 0) {

                } else if (b.denominator == 0) {
                    denominator = 0;
                    if ((numerator > 0) && (numerator != 1)) numerator = 1;
                    if ((numerator < 0) && (numerator != -1)) numerator = -1;                
                } else if (numerator == 0) {
                    numerator = b.numerator;
                    denominator = b.denominator;
                }
            } else {
                Field _denominator = gcd<Field, Field>(denominator, b.denominator);
                if (_denominator>1) {
                    numerator = denominator/_denominator*(b.numerator) - b.denominator/_denominator*(numerator);
                    denominator = denominator / _denominator * b.denominator;
                } else {
                    numerator = b.denominator*numerator - denominator*b.numerator;
                    denominator *= b.denominator;
                }
                Field divider = gcd<Field, Field>(denominator, numerator);
                numerator /= divider;
                denominator /= divider;
            }
            return *this;
        }
        rationals_type operator*=(const rationals_type &b) {
            if (denominator*numerator*b.denominator*b.numerator==0) {
                if (denominator == 0) {
                    if (b.numerator == 0) numerator = 0;
                } else if (b.denominator == 0) {
                    denominator = 0;
                    if (b.numerator == 0) numerator = 0;
                    if ((numerator > 0) && (numerator != 1)) numerator = 1;
                    if ((numerator < 0) && (numerator != -1)) numerator = -1;                
                } else if (numerator == 0) {

                } else if (b.numerator == 0) {
                    numerator = 0;
                    if (b.denominator * denominator > 0) denominator = 1; else denominator = -1;
                }
            } else {
                numerator *= b.numerator;
                denominator *= b.denominator;
                Field divider = gcd<Field, Field>(denominator, numerator);
                numerator /= divider;
                denominator /= divider;
                if (denominator < 0) {
                    numerator *= -1;
                    denominator *= -1;
                }
            }
            return *this;
        }
        rationals_type operator/=(const rationals_type &b) {
            if (denominator*numerator*b.denominator*b.numerator==0) {
                if (denominator == 0) {
                    if (b.denominator == 0) numerator = 0;
                } else if (b.numerator == 0) {
                    denominator = 0;
                    if (b.denominator == 0) numerator = 0;
                    if ((numerator > 0) && (numerator != 1)) numerator = 1;
                    if ((numerator < 0) && (numerator != -1)) numerator = -1;                
                } else if (numerator == 0) {

                } else if (b.denominator == 0) {
                    numerator = 0;
                    if (b.numerator * denominator > 0) denominator = 1; else denominator = -1;
                }
            } else {
                numerator *= b.denominator;
                denominator *= b.numerator;
                Field divider = gcd<Field, Field>(denominator, numerator);
                numerator /= divider;
                denominator /= divider;
                if (denominator < 0) {
                    numerator *= -1;
                    denominator *= -1;
                }
                if (numerator == 0) {
                    denominator = 1;
                }
            }
            return *this;
        }
        friend rationals_type operator+(rationals_type a, const rationals_type &b) {
            a += b;
            return a;
        }
        friend rationals_type operator-(rationals_type a, const rationals_type &b) {
            a -= b;
            return a;
        }
        friend rationals_type operator*(rationals_type a, const rationals_type &b) {
            a *= b;
            return a;
        }
        friend rationals_type operator/(rationals_type a, const rationals_type &b) {
            a /= b;
            return a;
        }
        friend rationals_type operator+(rationals_type a, const Field &c) {
            a += c;
            return a;
        }
        friend rationals_type operator-(rationals_type a, const Field &b) {
            a -= b;
            return a;
        }
        friend rationals_type operator*(rationals_type a, const Field &c) {
            a *= c;
            return a;
        }
        friend rationals_type operator/(rationals_type a, const Field &b) {
            a /= b;
            return a;
        }
        string to_string() const {
            ostringstream stream;
            if ((denominator == 1)||((numerator==0)&&(denominator==-1))) {
                stream << numerator;
            } else if (denominator*numerator == 0) {
                if (numerator!=0) {
                    if (numerator==1) stream << "+Infty"; else stream << "-Infty";
                } else {
                    stream << numerator << "/" << denominator;
                }
            } else {
                stream << numerator << "/" << denominator;
            }
            return stream.str();
        }
        friend string to_string(rationals_type a) {
            return a.to_string();
        }
        friend ostream& operator<<(ostream& os, const rationals_type& a) {
            os << a.to_string();
            return os;
        }
        friend inline int cmp(const rationals_type& a, const rationals_type& b) {
            return cmp(a-b,0);
        }
        friend inline int cmp(const rationals_type& a, const Field& b) {
            if ((a.get_numerator() == b) && (a.get_denominator() == 1)) return 0;
            rationals_type c = a-b;
            if (c.get_numerator()*c.get_denominator() > 0) { return 1; } else { return -1; }
        }
        friend inline bool operator==(const rationals_type& a, const rationals_type& b){ return cmp(a,b) == 0; }
        friend inline bool operator!=(const rationals_type& a, const rationals_type& b){ return cmp(a,b) != 0; }
        friend inline bool operator< (const rationals_type& a, const rationals_type& b){ return cmp(a,b) <  0; }
        friend inline bool operator> (const rationals_type& a, const rationals_type& b){ return cmp(a,b) >  0; }
        friend inline bool operator<=(const rationals_type& a, const rationals_type& b){ return cmp(a,b) <= 0; }
        friend inline bool operator>=(const rationals_type& a, const rationals_type& b){ return cmp(a,b) >= 0; }        
        friend inline bool operator==(const rationals_type& a, const Field& b){ return cmp(a,b) == 0; }
        friend inline bool operator!=(const rationals_type& a, const Field& b){ return cmp(a,b) != 0; }
        friend inline bool operator< (const rationals_type& a, const Field& b){ return cmp(a,b) <  0; }
        friend inline bool operator> (const rationals_type& a, const Field& b){ return cmp(a,b) >  0; }
        friend inline bool operator<=(const rationals_type& a, const Field& b){ return cmp(a,b) <= 0; }
        friend inline bool operator>=(const rationals_type& a, const Field& b){ return cmp(a,b) >= 0; }        
};

class RationalsIntField {
    public:
        using value_type=rationals<int>;
        using power_type=unsigned int;
        using coefficient_type=rationals<int>;
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
            return util::factor_rationals<value_type>(values, side);
        }
        inline static value_type gcd(const value_type& val1, const value_type& val2) {
            return value_type(1);
        }
        inline static value_type exp(const value_type& val) {
            return std::exp(val.get_numerator());
        }
        inline static value_type log(const value_type& val) {
            return std::log(val.get_numerator());
        }
};

class IntField {
    public:
        using value_type=int;
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
            enum { value = 1 };
        };
        inline static pair<value_type, vector<value_type> > factor(const vector<value_type> &values, const util::side_type &side=util::side_type::both) {
            return util::factor_int<value_type>(values, side);
        }
        inline static value_type gcd(const value_type& val1, const value_type& val2) {
            return std::gcd<value_type, value_type>(val1, val2);
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
        using value_type=long int;
        using power_type=unsigned long int;
        using coefficient_type=long int;
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
            return util::factor_int<value_type>(values, side);
        }
        inline static value_type gcd(const value_type& val1, const value_type& val2) {
            return std::gcd<value_type, value_type>(val1, val2);
        }
        inline static value_type exp(const value_type& val) {
            return std::exp(val);
        }
        inline static value_type log(const value_type& val) {
            return std::log(val);
        }
};

template<class Field> class Vector {
    public:

};

template<class Field> class LieAlgebra {
    public:

};

template<class Field, class LieAlgebra> class LieGroup {
    private:
        LieAlgebra lie_agebra;
    public:
        LieGroup() : lie_agebra() { }


};

};
};
#endif
