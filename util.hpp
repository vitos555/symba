#include<cmath>
#include<algorithm>
#include<variant>
#include<string>
#include<vector>
#include<iostream>
#include<map>
#include<numeric>
#include<sstream>

#ifndef _UTIL_HPP_
#define _UTIL_HPP_

namespace symba {
namespace util {
using namespace std;

template<class FieldClass> class EvaluationMap {
    private:
        using Field = typename FieldClass::value_type;
        map<string, Field> values;
    public:
        explicit EvaluationMap(const map<string, Field> &_values) : values(_values) { }
        Field get(string name) {
            return values[name];
        }
};

template<class Field> inline Field factorial(Field n) {
    Field factorial_list[] = {Field(1),Field(1),Field(2),Field(6),Field(24),Field(120),Field(720),Field(5040),Field(40320),Field(362880),Field(3628800)};
    if (n <= Field(10)) return factorial_list[n];
    else {
        Field ret = factorial_list[10];
        for(Field i=11; i <= n; i++) {
            ret *= i;
        }
        return ret;
    }
}

template<class Field> Field binomial_coefficient(Field n, Field k) {
    if (n <= Field(10)) return factorial<Field>(n)/factorial<Field>(k)/factorial<Field>(n-k);
    else {
        Field ret=1;
        for(Field i=n; i>=n-k+1; i++) {
            ret *= i;
        }
        return ret/factorial(k);
    }
}

template<class Field> Field multinomial_coefficient(Field n, typename vector<Field>::iterator ks_begin, typename vector<Field>::iterator ks_end) {
    Field ret = factorial<Field>(n);
    for(auto k=ks_begin; k!=ks_end; ++k) {
        ret /= factorial<Field>(*k);
    }
    return ret;
}

template<class Field> Field multinomial_coefficient(Field n, vector<Field> ks) {
    return multinomial_coefficient<Field>(n, ks.begin(), ks.end());
}

template<class Field, class PowerClass> Field pow(const Field &x, const PowerClass &pow) {
    PowerClass tmp_pow=pow, i=PowerClass(0), j;
    Field ret(1), tmp_ret;
    if (pow == PowerClass(0)) return Field(1);
    if (pow == PowerClass(1)) return x;
    while (tmp_pow > 0) {
        if (tmp_pow != ((tmp_pow >> 1) << 1)) {
            tmp_ret = x;
            for(j=0; j<i; j++) tmp_ret *= tmp_ret;
            ret *= tmp_ret;
        }
        tmp_pow = tmp_pow >> 1;
        i++;
    }
    return ret;
}

};
};
#endif
