#include<iostream>
#include<algorithm>
#include<vector>

using namespace std;
void generate(int n, const vector<int> &_ks) {
	vector<int> ks=_ks;
	do {
		cout << "Permutation: ";
		for(auto k=ks.begin(); k!=ks.end(); ++k) {
			cout << *k << ":";
		}
		cout << endl;
	} while (prev_permutation(ks.begin(), ks.end()));
}


int main() {
	generate(5, {6,4,4,3,2,1});
	return 0;
}