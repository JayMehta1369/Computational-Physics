#include <iostream>
#include <math.h>
#include <cmath>
#include <random>

using namespace std;

// G_inv function used for imp sampling
double G_inv (double u, int b) {
    // int b stands for the upperbround in our integral
    double kernel = 1 - u*((exp(b) - 1)/exp(b));
    double res = -log(kernel);
    return res;
}

// mc_term returns each terms in imp sampling
double mc_term (double x, int b) {
    double t = G_inv (x, b); // new input
    // function f and g are defined as follows
    double f = exp(-pow(t,2)); 
    // Same as above, int b stands for the upperbround in our integral
    double g = (exp(-t))*(exp(b)/(exp(b) - 1));
    double res = f/g;
    return res;
}

double simpleMC (int N) {
    std::random_device r; // Seed with a real random device
    std::seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
    std::mt19937 engine2(seed);
    std::uniform_real_distribution<> u(0,1);
    double t = u(engine2);
    double f = exp(-pow(t,2));

    double sum = 0;
    for (int n = 0; n < N; ++n) {
        sum += f;
    }
    double res = sum/N;
    return res;
}

int main() { 
    int N;
    cout << "Enter number of terms (N): " << endl;
    cin >> N;

    std::random_device r; // Seed with a real random device
    std::seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
    std::mt19937 engine2(seed);
    std::uniform_real_distribution<> u(0,1);

    // Calculating the MonteCarlo integration and the error terms
    double sum1 = 0;
    double error1 = 0;
    
    double sum2 = 0;
    double error2 = 0;

    for (int n = 0; n < N; ++n) {
        double I1 = mc_term(u(engine2),1);
        double I2 = mc_term(u(engine2),10);
        error1 += pow(I1,2);
        error2 += pow(I2,2);
        sum1 += I1;
        sum2 += I2;
    }
    double res1 = sum1/N; // for b = 1
    double err1 = error1/N - pow(res1,2);
    double res_err1 = sqrt(err1/N);

    double res2 = sum2/N; // for b = 10
    double err2 = error2/N - pow(res2,2);
    double res_err2 = sqrt(err2/N);

    cout << "\n";
    cout << "For N equals to " << N << ":\n";
    cout << "When the upperbound is 1" << endl;
    cout << "The integral using important sampling = " << res1 << endl;
    cout << "The error = " << res_err1 << endl;
    cout << "The simpleMC returns " << simpleMC(N) << endl;
    cout << "\n";
    cout << "When the upperbound is 10" << endl;
    cout << "The integral using important sampling = " << res2 << endl;
    cout << "The error = " << res_err2 << endl;
    return 0;
}