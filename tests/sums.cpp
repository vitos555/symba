#include<iostream>
#include<algorithm>
#include<vector>

using namespace std;

template<class Field> inline Field factorial(Field n) {
    Field factorial_list[] = {1,1,2,6,24,120,720,5040,40320,362880,3628800};
    if (n <= 10) return factorial_list[n];
    else {
        Field ret = factorial_list[10];
        for(Field i=11; i <= n; i++) {
            ret *= i;
        }
        return ret;
    }
}

template<class Field> Field multinomial_coefficient(Field n, typename vector<Field>::iterator ks_begin, typename vector<Field>::iterator ks_end) {
    Field ret = factorial<Field>(n);
    for(auto k=ks_begin; k!=ks_end; ++k) {
        ret /= factorial(*k);
    }
    return ret;
}

template<class Field, class PowerClass> void generate(Field power, const vector<Field> &terms) {
    Field stop_power = power / terms.size() + ((power % terms.size()) > 0);
    vector<PowerClass> powers=vector<PowerClass>(terms.size(),0);
    powers[0] = power;
    size_t pointer = 0; // last non-zero power pointer
    while(powers[0] >= stop_power) {
        Field c = multinomial_coefficient<Field>(power, powers.begin(), powers.begin()+pointer+1);
        do {
            cout << c;
            for (int i = 0; i < powers.size(); ++i) {
                if (powers[i] > 0) {
                    cout << "x" << "_" << i+1 << "^" << powers[i];
                }
            }
            cout << " + ";
        } while (prev_permutation(powers.begin(), powers.end()));

        // Generate next combination of powers
        if ((pointer < powers.size()-1) && ((pointer==0)||((powers[pointer] > powers[pointer+1])&&(powers[pointer]>1)))) {
            powers[pointer]--;
	        pointer++;
            powers[pointer]++;
        } else {
            Field sum = powers[pointer];
            Field level = 0;
            powers[pointer] = 0;
            pointer--;
            while((pointer>0) && (sum+1 > (powers[pointer]-1)*(powers.size()-pointer-1))) {
                sum += powers[pointer];
                powers[pointer] = 0;
                pointer--;
            }
            if ((pointer==0) && (sum+1 > (powers[pointer]-1)*(powers.size()-pointer-1))) break;
            level = --powers[pointer];
            sum++;
            while(sum > 0) {
                if (sum > level) {
                    powers[++pointer] = level;
                    sum -= level;
                } else {
                    powers[++pointer] = sum;
                    sum = 0;
                }
            }
        }
    }
    cout << endl;
}

int main() {
	vector<int> v(3,0), v1(4,0), v2(2,0), v3(5,0);
//	generate<int,int>(7, v);
    generate<int,int>(2, v2);
	generate<int,int>(3, v2);
	generate<int,int>(6, v);
//	generate<int,int>(7, v1);
//	generate<int,int>(6, v1);
//	generate<int,int>(7, v2);
//	generate<int,int>(7, v3);
	return 0;
}