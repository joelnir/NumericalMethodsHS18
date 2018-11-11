#include <iostream>
#include <algorithm>
#include <vector>
#include <sstream>
#include <string>
#include <complex>

using namespace std;

const complex<double> I(0,1); // imaginary unit
const double PI = 3.14159265359;

template<class T> struct duplet {
	int ind;
	T val;

	duplet(int p, T v) {
		ind = p;
		val = v;
	}
};

template<class T>
bool dupletSort(duplet<T> d1, duplet<T> d2){
    return d1.ind < d2.ind;
}

template<class T> struct sparse_vec {
	double tol = 1e-6;
	vector<duplet<T> > duplets;
	int len=0;

	sparse_vec(int l) {
		len = l;
	}

	void append(int ind, T val) {
        if(abs(val) >= this->tol){
            duplet<T> newDup(ind, val);
            this->duplets.push_back(newDup);
        }
	}

	void cleanup() {
        // Remove too high indexes and too small values
        vector<duplet<T>> newDuplets;
        for(duplet<T> d : this->duplets){
            if(d.ind < this->len && abs(d.val) >= this->tol){
                // Acceptable
                newDuplets.push_back(d);
            }
        }

        // Sort duplets
        sort(newDuplets.begin(), newDuplets.end(), dupletSort<T>);

        // Add up at same index
        this->duplets.clear();

        if(!newDuplets.empty()){
            duplet<T> curDuplet = duplet<T>(-1, 0);
            for(duplet<T> d : newDuplets){
                if(curDuplet.ind != d.ind){
                    // new index

                    if(curDuplet.ind != -1){
                        // Don't save initial, placeholder value
                        this->duplets.push_back(curDuplet);
                    }

                    curDuplet = d;
                }
                else{
                    // same index, add value
                    curDuplet.val += d.val;
                }
            }

            this->duplets.push_back(curDuplet);
        }
	}

	T get_val(int ind) const {
        auto bound = lower_bound(this->duplets.begin(),
                this-> duplets.end(), ind, dupletSort<T>);

        if(bound == this->duplets.end()){
            // index higher than any duplet
            return 0;
        }
        else if((*bound).ind == ind){
            return *bound;
        }
        else{
            return 0;
        }
	}

	static sparse_vec cwise_mult(const sparse_vec &a, const sparse_vec &b) {
		sparse_vec out(max(a.len,b.len));

        int ia = 0;
        int ib = 0;

        while(ia < a.duplets.size() && ib < b.duplets.size()){
            duplet<T> duplet_a = a.duplets.at(ia);
            duplet<T> duplet_b = b.duplets.at(ib);
            if(duplet_a.ind == duplet_b.ind){
                // Index of non-negligible elemnts coincide
                T res = duplet_a.val * duplet_b.val;

                if(abs(res) >= out.tol){
                    out.append(duplet_a.ind, res);
                }

                ia++;
                ib++;
            }
            else if(duplet_a.ind < duplet_b.ind){
                ia++;
            }
            else{
                // duplet_a.ind > duplet_b.ind
                ib++;
            }
        }

		return out;
	}

	static sparse_vec conv(const sparse_vec &a, const sparse_vec &b) {
		sparse_vec out(a.len + b.len - 1);

        int ia = 0;
        int ib = 0;

        for(int ia = 0; ia < a.duplets.size(); ++ia){
            int a_ind = a.duplets.at(ia).ind;
            T a_val = a.duplets.at(ia).val;

            for(int ib = (b.duplets.size() - 1); ib >= 0; --ib){
                int b_ind = b.duplets.at(ib).ind;
                T b_val = b.duplets.at(ib).val;

                out.append((a_ind + b_ind), (a_val*b_val));
            }
        }

        out.cleanup();
		return out;
	}

	static sparse_vec fft(const sparse_vec &x) {
		int n = x.len;
		sparse_vec tot(n);

        if(n == 1){
            // Base case, single element
            if(x.duplets.size() != 0){
                tot.append(0, x.duplets.at(0).val);
            }
        }
        else{
            // Divide and conquer
            int m = n/2;

            sparse_vec even(m);
            sparse_vec odd(m);

            // Assign to even and odd
            for(duplet<T> d : x.duplets){
                if(d.ind % 2 == 0){
                    even.append(d.ind/2, d.val);
                }
                else{
                    odd.append((d.ind - 1)/2, d.val);
                }
            }

            // Recursive call
            sparse_vec even_f = sparse_vec::fft(even);
            sparse_vec odd_f = sparse_vec::fft(odd);

            // Re-add even components
            for(duplet<T> ec : even_f.duplets){
                tot.append(ec.ind, ec.val);
                tot.append(ec.ind + m, ec.val);
            }

            // Re-add odd components
            for(duplet<T> oc : odd_f.duplets){
                complex<double> coeff = exp((complex<double>(0, -2*PI*oc.ind))/(double)n);

                tot.append(oc.ind, coeff * oc.val);
                tot.append(oc.ind + m, -coeff * oc.val);
            }

            // Perform cleanup to some components added from even and odd
            tot.cleanup();
        }

		return tot;
	}

	static sparse_vec ifft(const sparse_vec &x) {
		double n = x.len;
		sparse_vec out(n);

        // Take complex conjugate of x
        sparse_vec fft_input(n);
        for(duplet<T> d : x.duplets){
            fft_input.append(d.ind, conj(d.val));
        }

        // Forward fft on comp. conj.
        sparse_vec res_fft = sparse_vec::fft(fft_input);

        // Comp. conj. of result as well as (1/n) division
        for(duplet<T> d : res_fft.duplets){
            out.append(d.ind, conj(d.val) * (1/(double)n));
        }

        out.cleanup();

		return out;
	}

	static sparse_vec conv_fft(sparse_vec a, sparse_vec b) {
        // zero pad
        int n = b.len;

        a.len = 2*n;
        b.len = 2*n;

        // Transform
        sparse_vec<T> a_dft = sparse_vec::fft(a);
        sparse_vec<T> b_dft = sparse_vec::fft(b);

        // Mult
        sparse_vec<T> res_dft = sparse_vec::cwise_mult(a_dft, b_dft);

        // Inverse transform
        return sparse_vec::ifft(res_dft);
	}

	std::string to_string() const {
		std::stringstream ss;
		for (auto p : this->duplets) {
			ss << "(" << p.ind << "," << p.val << "),";
		}
		ss << "\n";
		std::string out = ss.str();
		return out;
	}
};




/***** TESTING ******/

int main() {

	sparse_vec<complex<double> > x(5);
	x.append(0,complex<double>(8.2,0));
	x.append(1,complex<double>(1,-2));
	x.append(3,complex<double>(-3,4.66));
	x.append(4,complex<double>(0,4));
	x.cleanup();

    sparse_vec<complex<double> > y(4);
	y.append(1,complex<double>(5,0));
	y.append(2,complex<double>(1.21,-4));
	y.append(3,complex<double>(4,2.4));
	y.cleanup();

	auto m = sparse_vec<complex<double> >::cwise_mult(x,y);
	m.cleanup();
	cout << "TESTS. Correct componentwise multiplication between x and y: (1,(5,-10)),(3,(-23.184,11.44)),\n";
	cout << "cwise_mult(x,y) = " << m.to_string();

	auto c = sparse_vec<complex<double> >::conv(x,y);
	c.cleanup();
	cout << "Correct exact discrete convolution between x and y: (1,(41,0)),(2,(14.922,-42.8)),(3,(26.01,13.26)),(4,(-6.2,17.7)),(5,(15.01,37.6386)),(6,(-7.184,16.28)),(7,(-9.6,16)),\n";
	cout << "conv(x,y) = " << c.to_string();
	auto cf = sparse_vec<complex<double> >::conv_fft(x,y);
	cf.cleanup();
	cout << "conv_fft(x,y) = " << cf.to_string();
}

/*
 * g)
 *
 * Yes, for very sparse x, y the asymptotic complexity of conv is much better.
 * (Consider for example SZa << log(n), SZb << log(n))
 *
 * In particular consider x, y single impulses at some index Xi, Yi.
 * The xomplexity of such signals for conv is O(1),
 * but conv_fft is still O(n*log(n)).
 *
 */

