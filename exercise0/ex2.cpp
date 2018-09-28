#include <iostream>
#include <cmath>
#include <functional>

using namespace std;

/*
 * a)
 */
double power(double a, int b){
    bool negative = false;

    if(b < 0){
        negative = true;
        b *= (-1);
    }

    double res = 1;
    for(int i = 0; i < b; ++i){
        res *= a;
    }

    if(negative){
        return 1.0/res;
    }
    else{
        return res;
    }
}

/*
 * b)
 */
double fast_power(double a, int b){
    if(b < 0){
        return 1.0 / fast_power(a, -b);
    }

    // Stop condition
    if(b == 0){
        return 1;
    }

    if(b % 2 == 0){
        // Even
        double sqrt = fast_power(a, b / 2);
        return sqrt * sqrt;
    }
    else{
        // Odd
        return a * fast_power(a, (b - 1));
    }
}
/*
 * c)
 *
 * a^b = exp(log(a^b)) = exp(b * log(a))
 */
double fast_power2(double a, int b){
    return exp(b * log(a));
}

/*
 * d)
 * No, not for all inputs.
 * It is implemented differently, using a method for minimizing numerical errors.
 */

// Function for manually inputting a and b
void manual_test(){
    while(true){
        double a;
        int b;
        cin >> a;
        cin >> b;

        cout.precision(17);
        cout << "In O(b): " << power(a, b) << endl;
        cout << "in O(log(b)): " << fast_power(a, b) << endl;
        cout << "in O(1): " << fast_power2(a, b) << endl;
        cout << "Using std::pow: " << pow(a, b) << endl;
    }
}

/*
 * e)
 */
double avg_error(function<double(double, int)> pow_func){
    double func_res, true_res, rel_error;
    int n = 0;

    for(double a = 1; a < 100; ++a){
        for(int b = 0; b < 300; ++b){
            func_res = pow_func(a,b);
            true_res = pow(a,b);

            if(!(isinf(func_res) || isinf(true_res))){
                ++n;

                rel_error += abs(func_res - true_res) / true_res;
            }
        }
    }

    return rel_error / n;
}
/*
 * The fast_power implementation minimizes the error.
 * The fast_power2 is completely dependent on the accuracy of the log- and exp-functions.
 * These migh be based on table lookups.
 *
 * Both power and fast_power introduce errors through repeated multiplication, but
 * fast_power perfroms less of those operations, resulting in a somewhat smaller error.
 */

int main(){
    cout << "Error power: " << avg_error(power) << endl;
    cout << "Error fast_power: " << avg_error(fast_power) << endl;
    cout << "Error fast_power2: " << avg_error(fast_power2) << endl;

    /*
     * Results:
     * Error power: 3.00885e-16
     * Error fast_power: 2.57389e-16
     * Error fast_power2: 1.9356e-14
     */
}
