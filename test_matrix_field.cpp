#include<iostream>

#include "fields.hpp"
#include "symba.hpp"

using namespace symba;
using namespace std;
using R33 = fields::matrix<int, 3, 3>;
using R3 = fields::vector<int, 3>;

int main() {
	R33 a({1,2,3,1,2,3,1,2,3}),b({1,1,1,1,1,1,1,1,1});
	R3 c({1,2,3});
	R33 d(c);
	cout<<"a: "<<a<<endl;
	cout<<"b: "<<b<<endl;
	cout<<"c: "<<c<<endl;
	cout<<"d: "<<d<<endl;
	cout<<"a+b: "<<a+b<<endl;
	cout<<"a-b: "<<a-b<<endl;
	cout<<"5*a: "<<5*a<<endl;
	cout<<"a*a: "<<a*a<<endl;
	cout<<"a*b: "<<a*b<<endl;
	cout<<"a.T: "<<a.T()<<endl;
	cout<<"a*c: "<<a*c<<endl;
	cout<<"c*a: "<<c*a<<endl;
	cout<<"ones: "<<R33::ones()<<endl;
	cout<<"I: "<<R33::I()<<endl;
	cout<<"det(I): "<<det(R33::I())<<endl;
	cout<<"det(ones): "<<det(R33::zeros())<<endl;
	return 0;
}