#include <iostream>
#include <cmath>
#include <climits>
#include <limits>
#include <float.h>

using namespace std;

/*
 * a)
 */

template<class T>
T factorial(T x){
    T product = 1;

    while(x >= 1){
        product *= x;
        --x;
    }

    return product;
}

/*
 * b)
 */
template<class T>
double dbl_factorial(T x){
    double x_whole = round(x);

    if(x_whole < 2){
        return 1;
    }

    return x_whole * dbl_factorial(x_whole - 1);
}

/*
 * c)
 *
 * int:
 * Gives correct values up until x = 13.
 *
 * long:
 * Gives correct values up until x = 21.
 *
 * double: (using dbl_factorial)
 * Gives correct values up until x = 171.
 *
 * The output stops being correct wheN the result > INT_MAX ( / LONG_MAX / DBL_MAX ).
 * Greater values can not be represented in the datatype.
 * Double can reach much higher calues because of the way it is structured,
 * with designated bits for representing the exponent.
 * This does however mean that the double will only store the
 * most significant positions of the number, and won't give a complete answer.
 */

int main(){
    while(true){
        double in;
        cin >> in;

        // Print with max precision
        cout.precision(numeric_limits<double>::max_digits10);

        cout << dbl_factorial(in) << endl;
    }
}
