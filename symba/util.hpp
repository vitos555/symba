
#ifndef _SYMBA_UTIL_HPP_
#define _SYMBA_UTIL_HPP_

#include<cmath>
#include<algorithm>
#include<string>
#include<vector>
#include<iostream>
#include<map>
#include<numeric>
#include<sstream>
#include<utility>
#include<iterator>

namespace symba {
namespace util {
using namespace std;

enum side_type {
    both = 0,
    left = 1,
    right = 2
};

template<class FieldClass> class EvaluationMap {
    private:
        using ValueType = typename FieldClass::value_type;
        map<string, ValueType> values;
    public:
        explicit EvaluationMap(const map<string, ValueType> &_values) : values(_values) { }
        ValueType get(string name) {
            return values[name];
        }
};

constexpr int cexpr_factorial(const int n){
    if (n <= 10) {
        if (n >= 5) {
            if (n==5) return 120;
            else if (n>7) {
                if (n==8) return 40320;
                else if (n==9) return 362880;
                else return 3628800;
            } else {
                if (n==7) return 5040;
                else return 720;
            }
        } else {
            if (n > 2) {
                if (n==3) return 6;
                else return 24;
            } else {
                if (n == 2) return 2;
                else return 1;
            }
        }
    } else {
        return n * cexpr_factorial(n-1);
    }
}

constexpr int cexpr_limited_factorial(const int n, const int limit){
    return (n > limit) ? n*cexpr_limited_factorial(n-1, limit): limit;
}

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

template<typename _Iterator>
    class strided_iterator
    : public iterator<typename iterator_traits<_Iterator>::iterator_category,
              typename iterator_traits<_Iterator>::value_type,
              typename iterator_traits<_Iterator>::difference_type,
              typename iterator_traits<_Iterator>::pointer,
                      typename iterator_traits<_Iterator>::reference>
    {
    protected:
      _Iterator current;

      typedef iterator_traits<_Iterator>        __traits_type;
      std::size_t stride;
    public:
      typedef _Iterator                 iterator_type;
      typedef typename __traits_type::difference_type   difference_type;
      typedef typename __traits_type::pointer       pointer;
      typedef typename __traits_type::reference     reference;

      constexpr strided_iterator(std::size_t _stride=1) : current(), stride(_stride) { }

      explicit constexpr strided_iterator(iterator_type __x, std::size_t _stride=1) : current(__x), stride(_stride)  { }

      constexpr strided_iterator(const strided_iterator& __x)
      : current(__x.current), stride(__x.stride) { }

      strided_iterator& operator=(const strided_iterator&) = default;

      template<typename _Iter>
        constexpr strided_iterator(const strided_iterator<_Iter>& __x)
    : current(__x.base()), stride(__x.get_stride()) { }

      constexpr iterator_type base() const { return current; }
      constexpr std::size_t get_stride() const { return stride; }

      constexpr reference operator*() const
      {
        _Iterator __tmp = current;
        return *__tmp;
      }

      constexpr strided_iterator&
      operator++()
      {
        current += stride;
        return *this;
      }

      constexpr strided_iterator
      operator++(int)
      {
        strided_iterator __tmp = *this;
        current += stride;
        return __tmp;
      }

      constexpr strided_iterator&
      operator--()
      {
        current -= stride;
        return *this;
      }

      constexpr strided_iterator
      operator--(int)
      {
        strided_iterator __tmp = *this;
        current -= stride;
        return __tmp;
      }

      constexpr strided_iterator
      operator+(difference_type _n) const
      { return strided_iterator(current + _n*stride); }

      constexpr strided_iterator&
      operator+=(difference_type _n)
      {
        current += _n*stride;
        return *this;
      }

      constexpr strided_iterator
      operator-(difference_type _n) const
      { return strided_iterator(current - _n*stride); }

      constexpr strided_iterator&
      operator-=(difference_type _n)
      {
        current -= _n*stride;
        return *this;
      }

      constexpr reference
      operator[](difference_type _n) const
      { return *(*this + _n*stride); }      
};

  template<typename _Iterator>
    inline constexpr bool
    operator==(const strided_iterator<_Iterator>& __x,
           const strided_iterator<_Iterator>& __y)
    { return (__x.base() == __y.base()) && (__x.get_stride()==__y.get_stride()); }

  template<typename _Iterator>
    inline constexpr bool
    operator!=(const strided_iterator<_Iterator>& __x,
           const strided_iterator<_Iterator>& __y)
    { return !(__x == __y); }

  template<typename _IteratorL, typename _IteratorR>
    inline constexpr bool
    operator==(const strided_iterator<_IteratorL>& __x,
           const strided_iterator<_IteratorR>& __y)
    { return (__x.base() == __y.base()) && (__x.get_stride() == __y.get_stride()); }

  template<typename _IteratorL, typename _IteratorR>
    inline constexpr bool
    operator!=(const strided_iterator<_IteratorL>& __x,
           const strided_iterator<_IteratorR>& __y)
    { return !(__x == __y); }

  template<typename _Iterator>
    inline constexpr strided_iterator<_Iterator>
    operator+(typename strided_iterator<_Iterator>::difference_type __n,
          const strided_iterator<_Iterator>& __x)
    { return strided_iterator<_Iterator>(__x.base() + __n*__x.get_stride()); }

}; // namespace util
}; // namespace symba
#endif // _SYMBA_UTIL_HPP_
