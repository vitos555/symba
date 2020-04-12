#include<algorithm>
#include<iostream>

using namespace std;

template<class Field, class PowerClass> Field pow(const Field &x, const PowerClass &pow) {
    PowerClass tmp_pow=pow, i=0, j;
    Field ret=1, tmp_ret;
    if (pow == 0) return 1;
    if (pow == 1) return x;
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

int main() {
	cout << "pow(0,0): " << pow<int, int>(0, 0) << endl;
	cout << "pow(0,1): " << pow<int, int>(0, 1) << endl;
	cout << "pow(1,0): " << pow<int, int>(1, 0) << endl;
	cout << "pow(1,1): " << pow<int, int>(1, 1) << endl;
	cout << "pow(1,2): " << pow<int, int>(1, 2) << endl;
	cout << "pow(1,3): " << pow<int, int>(1, 3) << endl;
	cout << "pow(1,4): " << pow<int, int>(1, 4) << endl;
	cout << "pow(1,5): " << pow<int, int>(1, 5) << endl;
	cout << "pow(1,6): " << pow<int, int>(1, 6) << endl;
	cout << "pow(2,0): " << pow<int, int>(2, 0) << endl;
	cout << "pow(2,1): " << pow<int, int>(2, 1) << endl;
	cout << "pow(2,2): " << pow<int, int>(2, 2) << endl;
	cout << "pow(2,3): " << pow<int, int>(2, 3) << endl;
	cout << "pow(2,4): " << pow<int, int>(2, 4) << endl;
	cout << "pow(2,5): " << pow<int, int>(2, 5) << endl;
	cout << "pow(2,6): " << pow<int, int>(2, 6) << endl;
	cout << "pow(2,7): " << pow<int, int>(2, 7) << endl;
	cout << "pow(2,8): " << pow<int, int>(2, 8) << endl;
	cout << "pow(2,9): " << pow<int, int>(2, 9) << endl;
	cout << "pow(3,0): " << pow<int, int>(3, 0) << endl;
	cout << "pow(3,1): " << pow<int, int>(3, 1) << endl;
	cout << "pow(3,2): " << pow<int, int>(3, 2) << endl;
	cout << "pow(3,3): " << pow<int, int>(3, 3) << endl;
	cout << "pow(3,4): " << pow<int, int>(3, 4) << endl;
	cout << "pow(3,5): " << pow<int, int>(3, 5) << endl;
	cout << "pow(3,6): " << pow<int, int>(3, 6) << endl;
	return 0;
}