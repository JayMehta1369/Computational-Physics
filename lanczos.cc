#define BIT_SET(a,b) ((a) |= (1U<<(b)))
#define BIT_CLEAR(a,b) ((a) &= ~(1U<<(b)))
#define BIT_FLIP(a,b) ((a) ^= (1U<<(b)))
#define BIT_CHECK(a,b) ((bool)((a) & (1U<<(b))))
// Set on the condition f else clear
//bool f;         // conditional flag
//unsigned int m; // the bit mask
//unsigned int w; // the word to modify:  if (f) w |= m; else w &= ~m; 
#define COND_BIT_SET(a,b,f) ((a) = ((a) & ~(1U<<(b))) | ((-(unsigned int)f) & (1U<<(b))))

#include <random>
#include <chrono>
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>

using namespace std;
using namespace chrono;

typedef vector<double> Row; // One row of the matrix 
typedef vector<Row> Matrix; // Matrix: a vector of rows

inline double dotprod(const vector<double>& y, const vector<double>& x)
{ double sum = 0.0;
  for (int i=0;i<y.size();i++)
    sum += x[i]*y[i];
  return (sum);
}

// Hamiltonian-operated vector
void hv(vector<double>& y, const vector<double>& x,int L)
{ 
  for (vector<double>::iterator it = y.begin() ; it != y.end(); ++it)
                   *it=0.0;
  bool b ;
  unsigned int k;
  for (unsigned int i=0;i<x.size();i++){
    if (abs(x[i])>2.2e-16) {
      int jm = L-1;
      double xov2=x[i]/2.0;
      for (int j=0;j<L;j++){
        k=i;
        COND_BIT_SET(k,jm,BIT_CHECK(i,j));
        COND_BIT_SET(k,j,BIT_CHECK(i,jm));
        y[k] += xov2;
        jm = j;
      }
    }
  }
  for (unsigned int i=0;i<x.size();i++)
    y[i]=y[i]-((double) L)/2.0*x[i]/2.0;
}

// Single Jacobi Rotation of Matrix A, about the element (p,q)
void jac_rotate (Matrix & A, int & p, int & q, double & c, double & s) {
    int n = A.size();

    Row AP(n);
    for (int i = 0; i < n; i++) {
        AP.at(i) = A[p][i];
    }   
    for (int i = 0; i < n; i++) {
        A[p][i] = c*(A[p][i]) - s*(A[q][i]);
        A[q][i] = s*(AP.at(i)) + c*(A[q][i]);
    }  

    Row AP2(n);
    for (int i = 0; i < n; i++) {
        AP2.at(i) = A[i][p];
    }   
    for (int i = 0; i < n; i++) {
        A[i][p] = c*(A[i][p]) - s*(A[i][q]);
        A[i][q] = s*(AP2.at(i)) + c*(A[i][q]);
    }   
    
    A[p][q] = 0.0;
    A[q][p] = 0.0;
}

// Jacobi Diagonalization subroutine
void jacdiag (Matrix & A, vector<double> & d) {
    // Assuming that matrix A is n x n 
    // and hence, d is of size n too. 
    int n = d.size();
    // eps is the limting value of the off-diagonal elements
    double eps = 2.2e-16;

    // Initializing r to the diagonal of A
    Row r(n), v(n);
    // r transforms into d and v is the container of t*a_pq
    for (int p = 0; p < n; p++) {
        r.at(p) = A.at(p).at(p);
        d.at(p) = r.at(p); // updating d here
    }

    for (int sweep = 1; sweep <= 10; sweep++) {
        // Counting the sum of the squares of the remaining elements after each sweep
        double sum = 0.;
        for (int p = 0; p < (n-1); p++) {
            for (int q = (p+1); q < n; q++) {
                sum += (abs(A.at(p).at(q)))*(abs(A.at(p).at(q)));
            }   
        }

        // Applying the Jacobi rotations to all upperdiagonal elements of A       
        for (int p = 0; p < (n-1); p++) {
            for (int q = (p+1); q < n; q++) {
                
                double alpha = (0.5)*(A.at(q).at(q) - A.at(p).at(p))/(A.at(p).at(q));
                double t = 1.0/(abs(alpha) + sqrt(1.0 + (alpha*alpha)));
                if (alpha < 0.0) {
                    t = -t;
                }
                double c = 1.0/sqrt(1 + (t*t));
                double s = t*c;
                
                double diff;
                diff = t*(A.at(p).at(q));
                v[p] -= diff;
                v[q] += diff;
                d[p] -= diff;
                d[q] += diff;

                // A.at(p).at(q) = 0.0;
                // A.at(q).at(p) = 0.0;
                
                // Rotating the matrix A after one Jacobi rotation
                jac_rotate (A, p, q, c, s);   
                
            }   
        }

        // updating the eigenvalues' vector d
        for (int p = 0; p < n; p++) {
            r.at(p) += v.at(p);
            d.at(p) = r.at(p);
            v.at(p) = 0.0; // reinitializing v for each sweep 
        }

        cout << "Sweep = " << sweep << endl;
        cout << "The squared sum of the remaining off-diagonal elements is " << sum << endl;
    }

}

int main() {
    // define the geneator without seeding
    mt19937 generator;
    uniform_real_distribution<double> distribution(0,1.0);
    
    int L,m;
    cout << "Size of system, L: ";
    cin >> L;
    // cout << "Size of Lanczos matrix, m: ";
    // cin >> m;

    int N=pow(2,L);
    vector<double> v1(N),v2(N),f(N),omega(N);

    ofstream file1("E_0_14.txt"), file2("E_1_14.txt"), file3("diff_14.txt");

    vector<double> E_0_14, E_1_14, diff_14;
    double E_const = -4.5154463544;

    for (int lsize = 5; lsize <= 50; lsize += 5) {

        m = lsize;

        Matrix Lan(m,Row(m));

        // Initializing the random vector v1
        for (int i=0;i<N;i++){
            v1[i] = 1.0-2.0*distribution(generator);
            v2[i] = 0.0;
        }

        // Normalizing the vector v1
        double mod_v1 = sqrt(dotprod (v1, v1));
        for (int i = 0; i < N; i++) {
            v1.at(i) = ((v1.at(i))/(mod_v1));    
        }

        // Initializing the vector omega using the function hv
        hv (omega, v1, L);

        // Inserting the alpha_0 value
        double alpha_0 = dotprod (v1, omega);
        Lan.at(0).at(0) = alpha_0;

        // Updating vector f
        for (int i = 0; i < N; i++) {
            f.at(i) = (omega.at(i)) - ((alpha_0)*(v1.at(i)));    
        }

        // Building up the Lanczos matrix
        for (int j = 0; j < (m-1); j++) {

            // Updating beta values (off-diagonal)
            double mod_f = sqrt(dotprod (f, f));
            Lan[j][j+1] = mod_f;
            Lan[j+1][j] = mod_f;

            // Updating v2
            for (int i = 0; i < N; i++) {
                v2.at(i) = ((f.at(i))/(mod_f));    
            }

            // Updating omega
            hv (omega, v2, L);

            for (int i = 0; i < N; i++) {
                omega.at(i) -= (mod_f)*(v1.at(i));
                // Transforming v1 into v2 after each j
                v1.at(i) = v2.at(i);
            }

            // Inserting next alpha value (diagonal)
            Lan[j+1][j+1] = dotprod (v2, omega);

            // Updating vector f
            for (int i = 0; i < N; i++) {
                f.at(i) = (omega.at(i)) - ((Lan[j+1][j+1])*(v2.at(i)));    
            }
        }

        // Printing out the Lanczos matrix
        // cout << "\nLanczos Matrix = " << endl;
        for (Row &r : Lan) {
            for (double &d : r) {
                cout << d << " ";
            }
            cout << "\n";
        }

        Row eigenvalues (m);

        // cout << "\nStarting Diagonalization..." << endl;

        // Diagonalizing here
        jacdiag (Lan, eigenvalues);

        // cout << "\nAnd The diagonalized matrix Lan_D = " << endl;
        // for (Row &j : Lan) {
        //     for (double &k : j) {
        //         cout << k << ' '; 
        //     }
        //     cout << endl;
        // }

        sort (eigenvalues.begin(), eigenvalues.end());

        double lowest = eigenvalues.at(0);
        double excited = eigenvalues.at(1);
        double difference = abs(lowest-E_const);

        E_0_14.push_back(lowest);
        E_1_14.push_back(excited);
        diff_14.push_back(difference);

        file1 << m << "\t" << lowest << "\n";
        file2 << m << "\t" << excited << "\n";
        file3 << m << "\t" << difference << "\n";


        // cout << "\nthe vector containing Eigenvalues is = " << endl;
        // for (double &i : eigenvalues) {
        //     cout << i << ' ';    
        // }

        // cout << "\n\nThe END!" << endl;
    }

    file1.close(), file2.close(), file3.close();

    cout << "\nlowest eigenvalue convergence is as follows: " << endl;

    for (auto &i : E_0_14) {
        cout << i << " ";
    }

    cout << "\n\n";
    return 0;
}
