#include <iostream>
#include <stack>

using namespace std;

/*
 * a)
 * Integers are represented in binary two's complement.
 * The first bit is the sign bit, 0 = positive, 1 = negative.
 * For positive numbers this means that a 1 at place n from the right has value 2^n, indexed from 0.
 * Negative numbers are in twos complement. This means they can be calculated as the number that together with it's positive counterpart sums up to 2^N (N = 32 in this case) in binary. The sign bit is included.
 */

/*
 * b)
 */
void int_to_bits(int x){
    int cur_bit;
    stack<int> bits;

    // Push bits to stack
    for(int i = 0; i < 32; ++i){
        cur_bit = x & 1;
        bits.push(cur_bit);

        x = x >> 1;
    };

    // Print bits
    while(!bits.empty()){
        cout << bits.top();
        bits.pop();
    }

    cout << endl;
}

/*
 * c)
 * Floating point numbers are split into mantissa m and exponent e on the form f = m * 2^e.
 * The first bit is also here a sign bit.
 * Next 8 bits give the exponent, but a bias of -127 is added to allow for both positive and negative values.
 * Remaining 23 bits make up the mantissa. An implicit "1." is assumed, and each bit then adds the value 2^-n for position n from the left.
 */

/*
 * d)
 */
void float_to_bits(float x){
    int* xp = (int*)&x;
    int_to_bits(*xp);
}

int main(){
    while(true){
        cout << "Enter a float:" << endl;

        float in;
        cin >> in;

        float_to_bits(in);
    }
}
