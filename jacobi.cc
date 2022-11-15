#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace std;

typedef vector<double> Row; // One row of the matrix
typedef vector<Row> Matrix; // Matrix: a vector of rows

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

// Building up the Hamiltonian
void Hamiltonian (Matrix & H, int & M) {
    // Initializing the empty matrix H
    // MxM matrix of zeroes
    for (int r = 0; r < M; r++) {
        Row vec(M);
        H.push_back(vec);
    }

    double epsilon = 2.2e-16;
    for (int i = 0; i < M; i++) {
        int n = i + 1;
        for (int j = 0; j < M; j++) {
            int m = j + 1;
            if (j == i) {
                H.at(i).at(j) = 4*((n*n) + (5.05*n));
            }
            else {
                double term2 = 5.0;
                if ((abs(m-n))%2 != 0) {
                    term2 = -5.0;
                }
                H.at(i).at(j) = 4*(min(n,m))*(0.05 + term2);
            }
            if (abs(H.at(i).at(j)) < epsilon) {
                H.at(i).at(j) = 0;
            }
            // cout << H.at(i).at(j) << " ";
        }
        // cout << "\n";
    }
}

int main () {
    int m;
    cout << "Enter the value of M (size of Hamiltonian): " << endl;
    cin >> m;

    Matrix H;
    Row d(m);
    // Note: we must initialize the size of d (= m) here!

    // Initializing the Hamiltonian here
    Hamiltonian (H, m);
    
    cout << "\nInital Hamiltonian matrix H = " << endl;
    for (Row &j : H) {
        for (double &k : j) {
            cout << k << ' '; 
        }
        cout << endl;
    }

    cout << "\nStarting Diagonalization..." << endl;

    // Diagonalizing here
    jacdiag (H, d);

    cout << "\nAnd The diagonalized matrix H_D = " << endl;
    for (Row &j : H) {
        for (double &k : j) {
            cout << k << ' '; 
        }
        cout << endl;
    }

    sort (d.begin(), d.end());

    cout << "\nthe vector containing Eigenvalues is, d = " << endl;
    for (double &i : d) {
        cout << i << ' ';    
    }

    cout << "\n\nThe END!" << endl;
    return 0;
}
