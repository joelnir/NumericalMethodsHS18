#include<bits/stdc++.h>
#include<algorithm>

using namespace std;

float naive_sum(vector<float> &v) {
	float sum = 0;
	for(float x : v){
		sum += x;
	}
	return sum;
}
/*
 * a)
 *
 * With 1 in the first position the final value of the sum will be 1.
 * As the 1 is added the exponent of the float will be set to 10^0,
 * meaning that all the future additions will have to be done to
 * the mantissa directly. Since the values in [1e-8, 2e-8] are very
 * small they will not change the mantissa at all, basically
 * being rounded to 0.
 *
 * If the 1 is in the end of v, ~10^7 additions have already been done.
 * The sum is then in the interval [1e-1, 2e-1]. Adding 1 to this
 * doesn't change the exponent much, so only a small error is introduced
 * from rounding the last decimals in the mantissa. There is still some error
 * introduced by adding each small number to the sum, since as the sum grows
 * it is much larger than the individual elements of v, resulting in a loss of
 * accuracy in when they are added.
 */

/*
 * b)
 */

// Do some divide and conquer stuff
float sum_help(vector<float> &v){
    // Terminating case
    if(v.size() <= 1){
        return v[0];
    }

    vector<float> v1;
    vector<float> v2;
    for(int i = 0; i < v.size(); ++i){
        if(i < v.size() / 2){
            v1.push_back(v[i]);
        }
        else{
            v2.push_back(v[i]);
        }
    }

    return sum_help(v1) + sum_help(v2);
}

float sum(vector<float> v) {
    float sum = 0;

    // Detect and separate any too big numbers
    for(int i = 0; i < v.size(); ++i){
        if(v[i] > 0.1){
            sum += v[i];
            v[i] = 0;
        }
    }

    sum += sum_help(v);

    return sum;
}

/*
 * c)
 * The naive_ sum would still return 1 if the 1 is in the first place of v,
 * even though the number is now in the range [2, 3]. If the 1 is last in v it still
 * gives a fairly correct value. My sum function performs wel, still giving an
 * accurate value.
 */

/*
 * d)
 */
float acc_sum(vector<float> &v) {
	float sum = 0;
    float acc = 0;

    float oldSum = 0;

    for(float x : v){
        acc += x;

        oldSum = sum;
		sum += acc; // adding x + old error

        // Save cancellation error for next iteration
        acc = acc - (sum - oldSum);
	}

	return sum;
}

int main() {
	srand(time(0));
	cout << setprecision(15);
	int N = 1e7;
	double corr_sum = 0;
	vector<float> v(N);
	for (int i = 0; i < N; i++) {
		double x = 1e-8 + (rand() % 10) * 1e-9;
		corr_sum += x;
		v[i] = (float) x;
	}
	corr_sum += 1;
	v[0]=1;

	cout << "naive_sum(v) = " << naive_sum(v) << endl
		 << "      sum(v) = " << sum(v) << endl
		 << "  acc_sum(v) = " << acc_sum(v) << endl
		 << " correct sum = " << corr_sum << endl;
}
