#include<iostream>
#include<vector>
#include "symba.hpp"
#include "fields.hpp"

using namespace symba;
using namespace symba::util;
using namespace symba::fields;
using namespace symba::poly;

int main() {
    EvaluationMap<LongIntField> em = EvaluationMap<LongIntField>({ {"x", 1}, {"y", 2} });
    Constant<LongIntField> c = Constant<LongIntField>(2);
    cout << "Constant to_string: " << c.to_string() << endl;
    cout << "Constant evaluation: " << c.evaluate() << endl;
    Variable<LongIntField> x = Variable<LongIntField>("x");
    cout << "Variable to_string: " << x.to_string() << endl;
    cout << "Variable name: " << x.get_name() << endl;
    cout << "Variable evaluation: " << x.evaluate(em) << endl;
    Variable<LongIntField> y = Variable<LongIntField>("y");
    cout << "Variable to_string: " << y.to_string() << endl;
    cout << "Variable name: " << y.get_name() << endl;
    cout << "Variable evaluation: " << y.evaluate(em) << endl;
    Entity<LongIntField> ec = Entity<LongIntField>(c);
    cout << "Entity to_string: " << ec.to_string() << endl;
    cout << "Entity evaluation: " << ec.evaluate(em) << endl;
    Entity<LongIntField> ex = Entity<LongIntField>(x);
    cout << "Entity to_string: " << ex.to_string() << endl;
    cout << "Entity evaluation: " << ex.evaluate(em) << endl;
    Term<LongIntField> t1 = Term<LongIntField>(x, 2);
    cout << "Entity to_string: " << t1.to_string() << endl;
    cout << "Entity evaluation: " << t1.evaluate(em) << endl;
    Term<LongIntField> t2 = Term<LongIntField>(c, 5);
    cout << "Entity to_string: " << t2.to_string() << endl;
    cout << "Entity evaluation: " << t2.evaluate(em) << endl;
    t2.simplify();
    cout << "Entity to_string: " << t2.to_string() << endl;
    cout << "Entity evaluation: " << t2.evaluate(em) << endl;
    Term<LongIntField> t3 = Term<LongIntField>(y, 1);
    cout << "Entity to_string: " << t3.to_string() << endl;
    cout << "Entity evaluation: " << t3.evaluate(em) << endl;
    Term<LongIntField> t4 = Term<LongIntField>(c, 2);
    cout << "Entity to_string: " << t4.to_string() << endl;
    cout << "Entity evaluation: " << t4.evaluate(em) << endl;
    Monomial<LongIntField> m1 = Monomial<LongIntField>({t1, t2, t1, t3, t3, t4});
    cout << "Monomial to_string: " << m1.to_string() << endl;
    cout << "Monomial evaluation: " << m1.evaluate(em) << endl;
    m1.simplify();
    cout << "Monomial to_string: " << m1.to_string() << endl;
    cout << "Monomial evaluation: " << m1.evaluate(em) << endl;
    Polynomial<LongIntField> p1 = Polynomial<LongIntField>({m1, m1});
    cout << "Polynomial to_string: " << p1.to_string() << endl;
    cout << "Polynomial evaluation: " << p1.evaluate(em) << endl;
    p1.simplify();
    cout << "Polynomial to_string: " << p1.to_string() << endl;
    cout << "Polynomial evaluation: " << p1.evaluate(em) << endl;
    Monomial<LongIntField> m2 = Monomial<LongIntField>({t1, t3, t3, t4});
    m2.simplify();
    cout << "Monomial to_string: " << m2.to_string() << endl;
    cout << "Monomial evaluation: " << m2.evaluate(em) << endl;
    Polynomial<LongIntField> p2 = Polynomial<LongIntField>({m1, m1, m2});
    p2.simplify();
    cout << "Polynomial to_string: " << p2.to_string() << endl;
    cout << "Polynomial evaluation: " << p2.evaluate(em) << endl;
    Term<LongIntField> t5 = Term<LongIntField>(p2, 2);
    Monomial<LongIntField> m3 = Monomial<LongIntField>({t1, t3, t5, t5});
    Term<LongIntField> tm = Term<LongIntField>(m1, 2);
    Monomial<LongIntField> m4 = Monomial<LongIntField>({c, x, y, tm, tm, t4});
    m4.simplify();
    cout << "Monomial to_string: " << m4.to_string() << endl;
    cout << "Monomial evaluation: " << m4.evaluate(em) << endl;
    Polynomial<LongIntField> p3 = Polynomial<LongIntField>({m1, m1, m2, m3});
    cout << "Polynomial to_string: " << p3.to_string() << endl;
    cout << "Polynomial evaluation: " << p3.evaluate(em) << endl;
    p3.simplify();
    cout << "Polynomial to_string: " << p3.to_string() << endl;
    cout << "Polynomial evaluation: " << p3.evaluate(em) << endl;
    p3.expand();
    cout << "Polynomial to_string: " << p3.to_string() << endl;
    cout << "Polynomial evaluation: " << p3.evaluate(em) << endl;

    // Substitution tests
    Polynomial<LongIntField> p11 = Polynomial<LongIntField>(p1);
    SubstitutionMap<LongIntField> sem1 = SubstitutionMap<LongIntField>({ {"x", c} });
    cout << "Polynomial11 before substitution1: " << p11.to_string() << endl;
    p11.substitute(sem1);
    cout << "Polynomial11 after substitution1: " << p11.to_string() << endl;
    p11.simplify();
    cout << "Polynomial11 after substitution1(simplified): " << p11.to_string() << endl;
    Polynomial<LongIntField> p12 = Polynomial<LongIntField>(p1);
    SubstitutionMap<LongIntField> sem2 = SubstitutionMap<LongIntField>({ {"y", c} });
    cout << "Polynomial12 before substitution2: " << p12.to_string() << endl;
    p12.substitute(sem2);
    cout << "Polynomial12 after substitution2: " << p12.to_string() << endl;
    p12.simplify();
    cout << "Polynomial12 after substitution2(simplified): " << p12.to_string() << endl;
    Polynomial<LongIntField> p21 = Polynomial<LongIntField>(p2);
    cout << "Polynomial21 before substitution1: " << p21.to_string() << endl;
    p21.substitute(sem1);
    cout << "Polynomial21 after substitution1: " << p21.to_string() << endl;
    p21.simplify();
    cout << "Polynomial21 after substitution1(simplified): " << p21.to_string() << endl;
    Polynomial<LongIntField> p22 = Polynomial<LongIntField>(p2);
    cout << "Polynomial22 before substitution2: " << p22.to_string() << endl;
    p22.substitute(sem2);
    cout << "Polynomial22 after substitution2: " << p22.to_string() << endl;
    p22.simplify();
    cout << "Polynomial22 after substitution2(simplified): " << p22.to_string() << endl;
    Polynomial<LongIntField> p33 = Polynomial<LongIntField>(p2);
    SubstitutionMap<LongIntField> sem3 = SubstitutionMap<LongIntField>({ {"x", p1} });
    cout << "Polynomial33 before substitution2: " << p33.to_string() << endl;
    p33.substitute(sem3);
    cout << "Polynomial33 after substitution2: " << p33.to_string() << endl;
    p33.simplify();
    cout << "Polynomial33 after substitution2(simplified): " << p33.to_string() << endl;
    Polynomial<LongIntField> p44 = Polynomial<LongIntField>(p2);
    SubstitutionMap<LongIntField> sem4 = SubstitutionMap<LongIntField>({ {"x", p2} });
    cout << "Polynomial44 before substitution2: " << p44.to_string() << endl;
    p44.substitute(sem4);
    cout << "Polynomial44 after substitution2: " << p44.to_string() << endl;
    p44.simplify();
    cout << "Polynomial44 after substitution2(simplified): " << p44.to_string() << endl;
    p44.expand();
    cout << "Polynomial44 after substitution2(expanded): " << p44.to_string() << endl;
    return 0;
}