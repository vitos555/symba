#include<iostream>
#include<vector>
#include "symba.hpp"
#include "util.hpp"
#include "fields.hpp"

using namespace symba;
using namespace symba::util;
using namespace symba::fields;

int main() {
	rational<int> a(1,2),b(2,3),c(0,5),d(5),e(0,0),f(2,0),g=1;
	cout<<"a: "<<to_string(a)<<endl;
	cout<<"b: "<<to_string(b)<<endl;
	cout<<"c: "<<to_string(c)<<endl;
	cout<<"d: "<<to_string(d)<<endl;
	cout<<"e: "<<to_string(e)<<endl;
	cout<<"f: "<<to_string(f)<<endl;
    cout<<"g: "<<to_string(g)<<endl;
	cout<<"a+b: "<<to_string(a+b)<<endl;
	cout<<"a-b: "<<to_string(a-b)<<endl;
	cout<<"a*b: "<<to_string(a*b)<<endl;
	cout<<"a/b: "<<to_string(a/b)<<endl;
	cout<<"a+c: "<<to_string(a+c)<<endl;
	cout<<"a-d: "<<to_string(a-d)<<endl;
	cout<<"a*f: "<<to_string(a*f)<<endl;
	cout<<"a/f: "<<to_string(a/f)<<endl;
	cout<<"a*e: "<<to_string(a*e)<<endl;
	cout<<"a/e: "<<to_string(a/e)<<endl;

    EvaluationMap<RationalIntField> em = EvaluationMap<RationalIntField>({ {"x", a}, {"y", b}, {"z", c} });
    poly::Constant<RationalIntField> cc = poly::Constant<RationalIntField>(a);
    cout << "Constant to_string: " << cc.to_string() << endl;
    cout << "Constant evaluation: " << cc.evaluate() << endl;
    poly::Variable<RationalIntField> x = poly::Variable<RationalIntField>("x");
    cout << "Variable to_string: " << x.to_string() << endl;
    cout << "Variable name: " << x.get_name() << endl;
    cout << "Variable evaluation: " << x.evaluate(em) << endl;
    poly::Variable<RationalIntField> y = poly::Variable<RationalIntField>("y");
    cout << "Variable to_string: " << y.to_string() << endl;
    cout << "Variable name: " << y.get_name() << endl;
    cout << "Variable evaluation: " << y.evaluate(em) << endl;
    poly::Entity<RationalIntField> ec = poly::Entity<RationalIntField>(cc);
    cout << "Entity to_string: " << ec.to_string() << endl;
    cout << "Entity evaluation: " << ec.evaluate(em) << endl;
    poly::Entity<RationalIntField> ex = poly::Entity<RationalIntField>(x);
    cout << "Entity to_string: " << ex.to_string() << endl;
    cout << "Entity evaluation: " << ex.evaluate(em) << endl;
    poly::Term<RationalIntField> t1 = poly::Term<RationalIntField>(x, 2);
    cout << "Entity to_string: " << t1.to_string() << endl;
    cout << "Entity evaluation: " << t1.evaluate(em) << endl;
    poly::Term<RationalIntField> t2 = poly::Term<RationalIntField>(cc, 5);
    cout << "Entity to_string: " << t2.to_string() << endl;
    cout << "Entity evaluation: " << t2.evaluate(em) << endl;
    t2.simplify();
    cout << "Entity to_string: " << t2.to_string() << endl;
    cout << "Entity evaluation: " << t2.evaluate(em) << endl;
    poly::Term<RationalIntField> t3 = poly::Term<RationalIntField>(y, 1);
    cout << "Entity to_string: " << t3.to_string() << endl;
    cout << "Entity evaluation: " << t3.evaluate(em) << endl;
    poly::Term<RationalIntField> t4 = poly::Term<RationalIntField>(cc, 2);
    cout << "Entity to_string: " << t4.to_string() << endl;
    cout << "Entity evaluation: " << t4.evaluate(em) << endl;
    poly::Monomial<RationalIntField> m1 = poly::Monomial<RationalIntField>({t1, t2, t1, t3, t3, t4});
    cout << "Monomial to_string: " << m1.to_string() << endl;
    cout << "Monomial evaluation: " << m1.evaluate(em) << endl;
    m1.simplify();
    cout << "Monomial to_string: " << m1.to_string() << endl;
    cout << "Monomial evaluation: " << m1.evaluate(em) << endl;
    poly::Polynomial<RationalIntField> p1 = poly::Polynomial<RationalIntField>({m1, m1});
    cout << "Polynomial to_string: " << p1.to_string() << endl;
    cout << "Polynomial evaluation: " << p1.evaluate(em) << endl;
    p1.simplify();

	return 0;
}