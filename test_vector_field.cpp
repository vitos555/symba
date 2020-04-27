#include<iostream>

#include "symba.hpp"
#include "fields.hpp"

using namespace symba;
using namespace std;
using R2 = fields::vector<int, 2>;
using R3 = fields::vector<int, 3>;

int main() {
	R2 a({1,2}),b({2,3}),c({0,5}),d({5,0}),f({2,0}),g(1);
	R3 aa({1,2,3});
	cout<<"a: "<<to_string(a)<<endl;
	cout<<"aa: "<<to_string(aa)<<endl;
	cout<<"b: "<<to_string(b)<<endl;
	cout<<"c: "<<to_string(c)<<endl;
	cout<<"d: "<<to_string(d)<<endl;
	cout<<"f: "<<to_string(f)<<endl;
	cout<<"g: "<<to_string(g)<<endl;
	cout<<"a+b: "<<to_string(a+b)<<endl;
	cout<<"a-b: "<<to_string(a-b)<<endl;
	cout<<"5*a: " <<to_string(5*a)<<endl;
	cout<<"a.a: " <<to_string(a.dot(a))<<endl;
	cout<<"a.b: " <<to_string(a.dot(b))<<endl;
	cout<<"e1: "<<to_string(R2::e1())<<endl;
	cout<<"e1.e2: "<<to_string(R2::e1().dot(R2::e2()))<<endl;
	cout<<"e1+e2: "<<to_string(R2::e1()+R2::e2())<<endl;
	cout<<"e1[0]: "<<to_string(R2::e1()[0])<<endl;
	cout<<"a innner b: "<<to_string(a.get_row_matrix()*b.get_column_matrix())<<endl;
	cout<<"a outer b: "<<to_string(a.get_column_matrix()*b.get_row_matrix())<<endl;
	return 0;
}