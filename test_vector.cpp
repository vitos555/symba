#include<iostream>
#include<vector>

#include "fields.hpp"
#include "symba.hpp"
#include "vector.hpp"

using namespace symba;
using namespace symba::util;

int main() {
    EvaluationMap<fields::LongIntField> em = EvaluationMap<fields::LongIntField>({ {"x", 1}, {"y", 2}, {"z", 3} });
    polynomial::Constant<fields::LongIntField> cc = polynomial::Constant<fields::LongIntField>(2);
    cout << "polynomial::Constant to_string: " << cc.to_string() << endl;
    cout << "polynomial::Constant evaluation: " << cc.evaluate() << endl;
    polynomial::Variable<fields::LongIntField> x = polynomial::Variable<fields::LongIntField>("x");
    cout << "polynomial::Variable to_string: " << x.to_string() << endl;
    cout << "polynomial::Variable name: " << x.get_name() << endl;
    cout << "polynomial::Variable evaluation: " << x.evaluate(em) << endl;
    polynomial::Variable<fields::LongIntField> y = polynomial::Variable<fields::LongIntField>("y");
    cout << "polynomial::Variable to_string: " << y.to_string() << endl;
    cout << "polynomial::Variable name: " << y.get_name() << endl;
    cout << "polynomial::Variable evaluation: " << y.evaluate(em) << endl;
    polynomial::Entity<fields::LongIntField> ec = polynomial::Entity<fields::LongIntField>(cc);
    cout << "polynomial::Entity to_string: " << ec.to_string() << endl;
    cout << "polynomial::Entity evaluation: " << ec.evaluate(em) << endl;
    polynomial::Entity<fields::LongIntField> ex = polynomial::Entity<fields::LongIntField>(x);
    cout << "polynomial::Entity to_string: " << ex.to_string() << endl;
    cout << "polynomial::Entity evaluation: " << ex.evaluate(em) << endl;
    polynomial::Term<fields::LongIntField> t1 = polynomial::Term<fields::LongIntField>(x, 2);
    cout << "polynomial::Entity to_string: " << t1.to_string() << endl;
    cout << "polynomial::Entity evaluation: " << t1.evaluate(em) << endl;
    polynomial::Term<fields::LongIntField> t2 = polynomial::Term<fields::LongIntField>(cc, 5);
    cout << "polynomial::Entity to_string: " << t2.to_string() << endl;
    cout << "polynomial::Entity evaluation: " << t2.evaluate(em) << endl;
    t2.simplify();
    cout << "polynomial::Entity to_string: " << t2.to_string() << endl;
    cout << "polynomial::Entity evaluation: " << t2.evaluate(em) << endl;
    polynomial::Term<fields::LongIntField> t3 = polynomial::Term<fields::LongIntField>(y, 1);
    cout << "polynomial::Entity to_string: " << t3.to_string() << endl;
    cout << "polynomial::Entity evaluation: " << t3.evaluate(em) << endl;
    polynomial::Term<fields::LongIntField> t4 = polynomial::Term<fields::LongIntField>(cc, 2);
    cout << "polynomial::Entity to_string: " << t4.to_string() << endl;
    cout << "polynomial::Entity evaluation: " << t4.evaluate(em) << endl;
    polynomial::Monomial<fields::LongIntField> m1 = polynomial::Monomial<fields::LongIntField>({t1, t2, t1, t3, t3, t4});
    cout << "polynomial::Monomial to_string: " << m1.to_string() << endl;
    cout << "polynomial::Monomial evaluation: " << m1.evaluate(em) << endl;
    m1.simplify();
    cout << "polynomial::Monomial to_string: " << m1.to_string() << endl;
    cout << "polynomial::Monomial evaluation: " << m1.evaluate(em) << endl;
    polynomial::Monomial<fields::LongIntField> m2 = polynomial::Monomial<fields::LongIntField>({t1, t1, t3, t4});
    cout << "polynomial::Monomial to_string: " << m2.to_string() << endl;
    cout << "polynomial::Monomial evaluation: " << m2.evaluate(em) << endl;
    m2.simplify();
    cout << "polynomial::Monomial to_string: " << m2.to_string() << endl;
    cout << "polynomial::Monomial evaluation: " << m2.evaluate(em) << endl;
    polynomial::Polynomial<fields::LongIntField> p1 = polynomial::Polynomial<fields::LongIntField>({m1, m1});
    cout << "polynomial::Polynomial to_string: " << p1.to_string() << endl;
    cout << "polynomial::Polynomial evaluation: " << p1.evaluate(em) << endl;
    p1.simplify();
    cout << "polynomial::Polynomial to_string: " << p1.to_string() << endl;
    cout << "polynomial::Polynomial evaluation: " << p1.evaluate(em) << endl;
    polynomial::Polynomial<fields::LongIntField> p2 = polynomial::Polynomial<fields::LongIntField>({m1, m2});
    cout << "polynomial::Polynomial to_string: " << p2.to_string() << endl;
    cout << "polynomial::Polynomial evaluation: " << p2.evaluate(em) << endl;
    p2.simplify();
    cout << "polynomial::Polynomial to_string: " << p2.to_string() << endl;
    cout << "polynomial::Polynomial evaluation: " << p2.evaluate(em) << endl;
    p2.factor();
    cout << "polynomial::Polynomial to_string: " << p2.to_string() << endl;
    cout << "polynomial::Polynomial evaluation: " << p2.evaluate(em) << endl;

    // Test vector
    linalg::Vector<fields::LongIntField, 3> v = linalg::Vector<fields::LongIntField, 3>({p1, p2, ec});
    cout << "linalg::Vector to_string: " << v.to_string() << endl;
    cout << "linalg::Vector evaluate: " << v.evaluate(em) << endl;
    cout << "linalg::Vector evaluate_vector: " << v.evaluate_vector(em) << endl;
    v.simplify();
    cout << "linalg::Vector to_string: " << v.to_string() << endl;
    cout << "linalg::Vector evaluate_vector: " << v.evaluate_vector(em) << endl;
    v.expand();
    cout << "linalg::Vector to_string: " << v.to_string() << endl;
    cout << "linalg::Vector evaluate_vector: " << v.evaluate_vector(em) << endl;
    v.factor();
    cout << "linalg::Vector to_string: " << v.to_string() << endl;
    cout << "linalg::Vector evaluate_vector: " << v.evaluate_vector(em) << endl;
    core::Entity<fields::LongIntField> vdotv = v.dot(v);
    cout << "linalg::Vector dot.to_string: " << vdotv.to_string() << endl;
    cout << "linalg::Vector dot.evaluate: " << vdotv.evaluate(em) << endl;
    vdotv.simplify();
    cout << "linalg::Vector dot.to_string: " << vdotv.to_string() << endl;
    cout << "linalg::Vector dot.evaluate: " << vdotv.evaluate(em) << endl;
    return 0;
}