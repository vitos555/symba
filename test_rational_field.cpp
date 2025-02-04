#include<iostream>
#include<vector>

#include "fields.hpp"
#include "symba.hpp"

using namespace symba;
using namespace symba::util;

int main() {
	fields::rational<int> a(1,2),b(2,3),c(0,5),d(5),e(0,0),f(2,0),g=1,h(1,36),l(4,3);
	cout<<"a: "<<to_string(a)<<endl;
	cout<<"b: "<<to_string(b)<<endl;
	cout<<"c: "<<to_string(c)<<endl;
	cout<<"d: "<<to_string(d)<<endl;
	cout<<"e: "<<to_string(e)<<endl;
	cout<<"f: "<<to_string(f)<<endl;
    cout<<"g: "<<to_string(g)<<endl;
    cout<<"h: "<<to_string(h)<<endl;
    cout<<"l: "<<to_string(l)<<endl;
	cout<<"a+b: "<<to_string(a+b)<<endl;
    cout<<"a: "<<to_string(a)<<endl;
	cout<<"a-b: "<<to_string(a-b)<<endl;
	cout<<"a*b: "<<to_string(a*b)<<endl;
	cout<<"a/b: "<<to_string(a/b)<<endl;
	cout<<"a+c: "<<to_string(a+c)<<endl;
	cout<<"a-d: "<<to_string(a-d)<<endl;
	cout<<"a*f: "<<to_string(a*f)<<endl;
	cout<<"a/f: "<<to_string(a/f)<<endl;
	cout<<"a*e: "<<to_string(a*e)<<endl;
	cout<<"a/e: "<<to_string(a/e)<<endl;
    cout<<"h+l: "<<to_string(h+l)<<endl;
    cout<<"5*a: "<<to_string(5*a)<<endl;
    cout<<"a: "<<to_string(a)<<endl;
    cout<<"5+a: "<<to_string(5+a)<<endl;

    EvaluationMap<fields::RationalIntField> em = EvaluationMap<fields::RationalIntField>({ {"x", a}, {"y", b}, {"z", c} });
    polynomial::Constant<fields::RationalIntField> cc = polynomial::Constant<fields::RationalIntField>(a);
    cout << "polynomial::Constant to_string: " << cc.to_string() << endl;
    cout << "polynomial::Constant evaluation: " << cc.evaluate() << endl;
    polynomial::Variable<fields::RationalIntField> x = polynomial::Variable<fields::RationalIntField>("x");
    cout << "polynomial::Variable to_string: " << x.to_string() << endl;
    cout << "polynomial::Variable name: " << x.get_name() << endl;
    cout << "polynomial::Variable evaluation: " << x.evaluate(em) << endl;
    polynomial::Variable<fields::RationalIntField> y = polynomial::Variable<fields::RationalIntField>("y");
    cout << "polynomial::Variable to_string: " << y.to_string() << endl;
    cout << "polynomial::Variable name: " << y.get_name() << endl;
    cout << "polynomial::Variable evaluation: " << y.evaluate(em) << endl;
    polynomial::Entity<fields::RationalIntField> ec = polynomial::Entity<fields::RationalIntField>(cc);
    cout << "polynomial::Entity to_string: " << ec.to_string() << endl;
    cout << "polynomial::Entity evaluation: " << ec.evaluate(em) << endl;
    polynomial::Entity<fields::RationalIntField> ex = polynomial::Entity<fields::RationalIntField>(x);
    cout << "polynomial::Entity to_string: " << ex.to_string() << endl;
    cout << "polynomial::Entity evaluation: " << ex.evaluate(em) << endl;
    polynomial::Term<fields::RationalIntField> t1 = polynomial::Term<fields::RationalIntField>(x, 2);
    cout << "polynomial::Entity to_string: " << t1.to_string() << endl;
    cout << "polynomial::Entity evaluation: " << t1.evaluate(em) << endl;
    polynomial::Term<fields::RationalIntField> t2 = polynomial::Term<fields::RationalIntField>(cc, 5);
    cout << "polynomial::Entity to_string: " << t2.to_string() << endl;
    cout << "polynomial::Entity evaluation: " << t2.evaluate(em) << endl;
    t2.simplify();
    cout << "polynomial::Entity to_string: " << t2.to_string() << endl;
    cout << "polynomial::Entity evaluation: " << t2.evaluate(em) << endl;
    polynomial::Term<fields::RationalIntField> t3 = polynomial::Term<fields::RationalIntField>(y, 1);
    cout << "polynomial::Entity to_string: " << t3.to_string() << endl;
    cout << "polynomial::Entity evaluation: " << t3.evaluate(em) << endl;
    polynomial::Term<fields::RationalIntField> t4 = polynomial::Term<fields::RationalIntField>(cc, 2);
    cout << "polynomial::Entity to_string: " << t4.to_string() << endl;
    cout << "polynomial::Entity evaluation: " << t4.evaluate(em) << endl;
    polynomial::Monomial<fields::RationalIntField> m1 = polynomial::Monomial<fields::RationalIntField>({t1, t2, t1, t3, t3, t4});
    cout << "polynomial::Monomial to_string: " << m1.to_string() << endl;
    cout << "polynomial::Monomial evaluation: " << m1.evaluate(em) << endl;
    m1.simplify();
    cout << "polynomial::Monomial to_string: " << m1.to_string() << endl;
    cout << "polynomial::Monomial evaluation: " << m1.evaluate(em) << endl;
    polynomial::Monomial<fields::RationalIntField> m2 = polynomial::Monomial<fields::RationalIntField>({t1, t1, t3, t4});
    cout << "polynomial::Monomial to_string: " << m2.to_string() << endl;
    cout << "polynomial::Monomial evaluation: " << m2.evaluate(em) << endl;
    m2.simplify();
    cout << "polynomial::Monomial to_string: " << m2.to_string() << endl;
    cout << "polynomial::Monomial evaluation: " << m2.evaluate(em) << endl;
    polynomial::Polynomial<fields::RationalIntField> p1 = polynomial::Polynomial<fields::RationalIntField>({m1, m1});
    cout << "polynomial::Polynomial to_string: " << p1.to_string() << endl;
    cout << "polynomial::Polynomial evaluation: " << p1.evaluate(em) << endl;
    p1.simplify();
    cout << "polynomial::Polynomial to_string: " << p1.to_string() << endl;
    cout << "polynomial::Polynomial evaluation: " << p1.evaluate(em) << endl;
    polynomial::Polynomial<fields::RationalIntField> p2 = polynomial::Polynomial<fields::RationalIntField>({m1, m2});
    cout << "polynomial::Polynomial to_string: " << p2.to_string() << endl;
    cout << "polynomial::Polynomial evaluation: " << p2.evaluate(em) << endl;
    p2.simplify();
    cout << "polynomial::Polynomial to_string: " << p2.to_string() << endl;
    cout << "polynomial::Polynomial evaluation: " << p2.evaluate(em) << endl;
    p2.factor();
    cout << "polynomial::Polynomial to_string: " << p2.to_string() << endl;
    cout << "polynomial::Polynomial evaluation: " << p2.evaluate(em) << endl;
	return 0;
}