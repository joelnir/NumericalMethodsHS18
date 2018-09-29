
using namespace std;

struct vec{
private:
    int capacity;
    int size;
    double* data;
public:
    /*
     * a)
     */
    vec(): capacity(10), size(0) {
        data = new double[10];
    }

    /*
     * b)
     */
    void push_back(double x){
        if(size >= capacity){
            // Allocate new, longer, array
            double* new_arr = new double[2*capacity];

            // Copy over existing elements
            for(int i = 0; i < size; ++i){
                new_arr[i] = data[i];
            }

            delete[] data;
            data = new_arr;
            capacity *= 2;
        }

        // Insert value
        data[size] = x;
        ++size;
    }
};

int main(){
    vec my_vec;

    for(int i = 0; i < 100; ++i){
        my_vec.push_back(i);
    }
}
