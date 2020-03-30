#include<cmath>
#include<algorithm>
#include<variant>
#include<string>
#include<vector>
#include<iostream>
#include<map>
#include<numeric>
#include<sstream>


#ifndef _FIELDS_HPP_
#define _FIELDS_HPP_

namespace symba {
namespace fields {
using namespace std;

template<class Field> class rational {
    private:
        using rational_type = rational<Field>;
        Field numerator;
        Field denominator;
    public:
        rational() : numerator(1), denominator(0) { }
        rational(Field _numerator, Field _denominator) : numerator(_numerator), denominator(_denominator) {
            if ((numerator == 0)&&(denominator==0)) { }
            else if ((numerator == 0)&&(denominator != 1)&&(denominator != -1)) {
                if (denominator > 0) denominator = 1; else denominator = -1;
            } else if ((denominator == 0)&&(numerator != 1) && (numerator != -1)) {
                if (numerator > 0) numerator = 1; else numerator = -1;
            }
        }
        rational(Field _numerator) : numerator(_numerator), denominator(1) { }
        const Field &get_numerator() const { return numerator; }
        const Field &get_denominator() const { return denominator; }
        rational_type operator+=(const Field &c) {
            numerator += c*denominator;
            return *this;
        }
        rational_type operator-=(const Field &c) {
            numerator -= c*denominator;
            return *this;
        }
        rational_type operator*=(const Field &c) {
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
        rational_type operator/=(const Field &c) {
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
        // rational_type operator-(const rational_type &a) {
        //     numerator = (-1)*numerator;
        //     return *this;
        // }
        rational_type operator+=(const rational_type &b) {
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
                    denominator = _denominator;
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
        rational_type operator-=(const rational_type &b) {
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
                    denominator = _denominator;
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
        rational_type operator*=(const rational_type &b) {
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
        rational_type operator/=(const rational_type &b) {
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
        friend rational_type operator+(rational_type a, const rational_type &b) {
            a += b;
            return a;
        }
        friend rational_type operator-(rational_type a, const rational_type &b) {
            a -= b;
            return a;
        }
        friend rational_type operator*(rational_type a, const rational_type &b) {
            a *= b;
            return a;
        }
        friend rational_type operator/(rational_type a, const rational_type &b) {
            a /= b;
            return a;
        }
        friend rational_type operator+(rational_type a, const Field &c) {
            a += c;
            return a;
        }
        friend rational_type operator-(rational_type a, const Field &b) {
            a -= b;
            return a;
        }
        friend rational_type operator*(rational_type a, const Field &c) {
            a *= c;
            return a;
        }
        friend rational_type operator/(rational_type a, const Field &b) {
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
        friend string to_string(rational_type a) {
            return a.to_string();
        }
        friend ostream& operator<<(ostream& os, const rational_type& a) {
            os << a.to_string();
            return os;
        }
        friend inline int cmp(const rational_type& a, const rational_type& b) {
            return cmp(a-b,0);
        }
        friend inline int cmp(const rational_type& a, const Field& b) {
            if ((a.get_numerator() == b) && (a.get_denominator() == 1)) return 0;
            rational_type c = a-b;
            if (c.get_numerator()*c.get_denominator() > 0) { return 1; } else { return -1; }
        }
        friend inline bool operator==(const rational_type& a, const rational_type& b){ return cmp(a,b) == 0; }
        friend inline bool operator!=(const rational_type& a, const rational_type& b){ return cmp(a,b) != 0; }
        friend inline bool operator< (const rational_type& a, const rational_type& b){ return cmp(a,b) <  0; }
        friend inline bool operator> (const rational_type& a, const rational_type& b){ return cmp(a,b) >  0; }
        friend inline bool operator<=(const rational_type& a, const rational_type& b){ return cmp(a,b) <= 0; }
        friend inline bool operator>=(const rational_type& a, const rational_type& b){ return cmp(a,b) >= 0; }        
        friend inline bool operator==(const rational_type& a, const Field& b){ return cmp(a,b) == 0; }
        friend inline bool operator!=(const rational_type& a, const Field& b){ return cmp(a,b) != 0; }
        friend inline bool operator< (const rational_type& a, const Field& b){ return cmp(a,b) <  0; }
        friend inline bool operator> (const rational_type& a, const Field& b){ return cmp(a,b) >  0; }
        friend inline bool operator<=(const rational_type& a, const Field& b){ return cmp(a,b) <= 0; }
        friend inline bool operator>=(const rational_type& a, const Field& b){ return cmp(a,b) >= 0; }        
};
class RationalIntField {
    public:
        using value_type=rational<int>;
        using power_type=unsigned int;
};

class IntField {
    public:
        using value_type=int;
        using power_type=unsigned int;
};

class LongIntField {
    public:
        using value_type=long int;
        using power_type=unsigned long int;
};
};
};
#endif
