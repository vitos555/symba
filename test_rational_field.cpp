#include<iostream>
#include<vector>

#include "fields.hpp"
#include "symba.hpp"

using namespace symba;
using namespace symba::util;
using namespace symba::fields;
using namespace symba::polynomial;

int main() {
	rationals<int> a(1,2),b(2,3),c(0,5),d(5),e(0,0),f(2,0),g=1,h(1,36),l(4,3);
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

    EvaluationMap<RationalsIntField> em = EvaluationMap<RationalsIntField>({ {"x", a}, {"y", b}, {"z", c} });
    Constant<RationalsIntField> cc = Constant<RationalsIntField>(a);
    cout << "Constant to_string: " << cc.to_string() << endl;
    cout << "Constant evaluation: " << cc.evaluate() << endl;
    Variable<RationalsIntField> x = Variable<RationalsIntField>("x");
    cout << "Variable to_string: " << x.to_string() << endl;
    cout << "Variable name: " << x.get_name() << endl;
    cout << "Variable evaluation: " << x.evaluate(em) << endl;
    Variable<RationalsIntField> y = Variable<RationalsIntField>("y");
    cout << "Variable to_string: " << y.to_string() << endl;
    cout << "Variable name: " << y.get_name() << endl;
    cout << "Variable evaluation: " << y.evaluate(em) << endl;
    Entity<RationalsIntField> ec = Entity<RationalsIntField>(cc);
    cout << "Entity to_string: " << ec.to_string() << endl;
    cout << "Entity evaluation: " << ec.evaluate(em) << endl;
    Entity<RationalsIntField> ex = Entity<RationalsIntField>(x);
    cout << "Entity to_string: " << ex.to_string() << endl;
    cout << "Entity evaluation: " << ex.evaluate(em) << endl;
    Term<RationalsIntField> t1 = Term<RationalsIntField>(x, 2);
    cout << "Entity to_string: " << t1.to_string() << endl;
    cout << "Entity evaluation: " << t1.evaluate(em) << endl;
    Term<RationalsIntField> t2 = Term<RationalsIntField>(cc, 5);
    cout << "Entity to_string: " << t2.to_string() << endl;
    cout << "Entity evaluation: " << t2.evaluate(em) << endl;
    t2.simplify();
    cout << "Entity to_string: " << t2.to_string() << endl;
    cout << "Entity evaluation: " << t2.evaluate(em) << endl;
    Term<RationalsIntField> t3 = Term<RationalsIntField>(y, 1);
    cout << "Entity to_string: " << t3.to_string() << endl;
    cout << "Entity evaluation: " << t3.evaluate(em) << endl;
    Term<RationalsIntField> t4 = Term<RationalsIntField>(cc, 2);
    cout << "Entity to_string: " << t4.to_string() << endl;
    cout << "Entity evaluation: " << t4.evaluate(em) << endl;
    Monomial<RationalsIntField> m1 = Monomial<RationalsIntField>({t1, t2, t1, t3, t3, t4});
    cout << "Monomial to_string: " << m1.to_string() << endl;
    cout << "Monomial evaluation: " << m1.evaluate(em) << endl;
    m1.simplify();
    cout << "Monomial to_string: " << m1.to_string() << endl;
    cout << "Monomial evaluation: " << m1.evaluate(em) << endl;
    Monomial<RationalsIntField> m2 = Monomial<RationalsIntField>({t1, t1, t3, t4});
    cout << "Monomial to_string: " << m2.to_string() << endl;
    cout << "Monomial evaluation: " << m2.evaluate(em) << endl;
    m2.simplify();
    cout << "Monomial to_string: " << m2.to_string() << endl;
    cout << "Monomial evaluation: " << m2.evaluate(em) << endl;
    Polynomial<RationalsIntField> p1 = Polynomial<RationalsIntField>({m1, m1});
    cout << "Polynomial to_string: " << p1.to_string() << endl;
    cout << "Polynomial evaluation: " << p1.evaluate(em) << endl;
    p1.simplify();
    cout << "Polynomial to_string: " << p1.to_string() << endl;
    cout << "Polynomial evaluation: " << p1.evaluate(em) << endl;
    Polynomial<RationalsIntField> p2 = Polynomial<RationalsIntField>({m1, m2});
    cout << "Polynomial to_string: " << p2.to_string() << endl;
    cout << "Polynomial evaluation: " << p2.evaluate(em) << endl;
    p2.simplify();
    cout << "Polynomial to_string: " << p2.to_string() << endl;
    cout << "Polynomial evaluation: " << p2.evaluate(em) << endl;
    p2.factor();
    cout << "Polynomial to_string: " << p2.to_string() << endl;
    cout << "Polynomial evaluation: " << p2.evaluate(em) << endl;
	return 0;
}