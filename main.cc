#include <iostream>
#include <math.h>
#include <cmath>
#include <random>
#include <fstream>
#include "MCvar.h" // header file

using namespace std; 

// randomGen() generates a random number b/w 0 and 1
double randomGen () {
    random_device r; // Seed with a real random device
    seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
    mt19937 engine(seed);
    uniform_real_distribution<> u(0,1);

    return u(engine);
}

// Calculating the probability of flipping the spin
double Pflip (double dE, double T) {
    double p = exp(((-1)*dE)/T);
    double res;
    if (p >= 1.0) {
        res = 1.0;
    }
    else {
        res = p;
    }
    return res;
}

// Initializing the Spin Configuration
void initialize (int SpinConf[], int L) {
    int N = L*L;
    // Following routine initializes a random configuration
    for (int n = 1; n <= N; n++) {
        double r = randomGen ();
            if (r < 0.5) {
                SpinConf[n] = 1;
            }
            else {
                SpinConf[n] = -1;
            }   
        // cout << SpinConf[n] << endl;
        }
}

// Runs 'Nmcs' Monte Carlo Sweeps on SpinConf[]
void MCSweeps (int Nmcs, int SpinConf[], int L, double T) {
    int N = L*L;

    // 2D array containing nearest neighbors of each spins
    int neighbor[N+1][4];
    for (int n = 1; n <= N; n++) {
        // up
        if ((n+L) <= N) {
            neighbor[n][0] = n+L;
        }
        else {
            neighbor[n][0] = n-(L*(L-1));
        }
        // right
        if ((n%L) != 0) {
            neighbor[n][1] = n+1;
        }
        else {
            neighbor[n][1] = n-L+1;
        }
        // down
        if ((n-L) >= 1) {
            neighbor[n][2] = n-L;
        }
        else {
            neighbor[n][2] = n+(L*(L-1));
        }
        // left
        if (((n-1)%L) != 0) {
            neighbor[n][3] = n-1;
        }
        else {
            neighbor[n][3] = n+L-1;
        }
    
    }

    // Monte Carlo Algorithm
    int sweep = 0;
    while (sweep < Nmcs) {
        int count = 0;
        while (count < N) {
            // generating a random site
            double r = randomGen ();
            int k = r*N;
            int i = k+1; // i-th site has been chosen

            int nn = 0; 
            // nn stands for the total spin of nearest neighbors
            for (int j = 0; j < 4; j++) {
                int index = neighbor[i][j];
                nn += SpinConf[index];
            }
            // clearly, range of nn = (-4 to 4)

            // dE is the difference in Energy when i-th spin is flipped
            double dE = 2*(SpinConf[i])*(nn);

            // metropolis algorithm
            double r2 = randomGen ();
            if (Pflip(dE, T) > r2) {
                SpinConf[i] = (-1)*(SpinConf[i]); 
            }
            count++;
        }
        sweep++;
    }
}

// Energy() calculates the total energy of the SpinConf[]
double Energy (int SpinConf[], int L) {
    int N = L*L;

    // 2D array containing nearest neighbors of each spins
    int neighbor[N+1][4];
    for (int n = 1; n <= N; n++) {
        // up
        if ((n+L) <= N) {
            neighbor[n][0] = n+L;
        }
        else {
            neighbor[n][0] = n-(L*(L-1));
        }
        // right
        if ((n%L) != 0) {
            neighbor[n][1] = n+1;
        }
        else {
            neighbor[n][1] = n-L+1;
        }
        // down
        if ((n-L) >= 1) {
            neighbor[n][2] = n-L;
        }
        else {
            neighbor[n][2] = n+(L*(L-1));
        }
        // left
        if (((n-1)%L) != 0) {
            neighbor[n][3] = n-1;
        }
        else {
            neighbor[n][3] = n+L-1;
        }
    
    }

    int total_E2 = 0; 
    double E;
    // total_E2 stands for the total energy (counted twice)
    for (int i = 1; i <= N; i++) {
        for (int j = 0; j < 4; j++) {
            int index = neighbor[i][j];
            total_E2 += (SpinConf[i])*(SpinConf[index]);
        }
    }
    E = (-0.5)*(total_E2);
    return E;
}

// Magnetization() calculates the total magnetization of the SpinConf[]
double Magnetization (int SpinConf[], int L) {
    int N = L*L;
    double M = 0; // M is the total magnetization
    for (int i = 1; i <= N; i++) {
        M += SpinConf[i];
    }
    return M;
}

// do_measurement & printout subroutine
// void printout (int Nmeas, double E[], double E2[], double E4[], double M[], double M2[], double M4[], double avrgs[]) {
//     MCvar<double> e, e2, e4, m, m2, m4;
//     for (int meas = 1; meas <= Nmeas; meas++) {
//      // e.push(E[meas]);
        // e2.push(E2[meas]);
        // e4.push(E4[meas]);
        // m.push(M[meas]);
        // m2.push(M2[meas]);
        // m4.push(M4[meas]);
//  }}

// Simulation
void SQIsing (int L, double T, int NWarmup, int NStep, int Nmeas) {

    // Create and open text files
    // ofstream file4("g_8.txt"), file5("g_12.txt"), file6("g_16.txt");

    // L_tilda takes values 8, 12, 16 respectively
    int L_tilda = L;
    for (int l = 0; l < 1; l++) {
        L_tilda += 4*l;

        // Using one-dimensional array, indexing from 1 to L^2
        int N = (L_tilda)*(L_tilda);
        // MySpinConf is the initial configuration here
        int MySpinConf[N+1]; 
        MySpinConf[0] = 1;
        
        ofstream file1("E_avg.txt"), file2("chi.txt"), file3("C_v.txt"); 

        // T_tilda runs from 2.0 to 2.5 
        double T_tilda = T;
        int pts = 50; // loop runs through (pts+1) points
        for (int i = 0; i <= pts; i++) {
            T_tilda += (0.5)*(i/pts);

            // Initializing
            initialize (MySpinConf, L_tilda);

            // performing warmup sweeps
            MCSweeps (NWarmup, MySpinConf, L_tilda, T_tilda);

            // containers for measurement
            double E[Nmeas+1], E2[Nmeas+1], E4[Nmeas+1];
            double M[Nmeas+1], M2[Nmeas+1], M4[Nmeas+1];

            double avrgs[6];

            // performing 'Nmeas' measurements
            for (int meas = 1; meas <= Nmeas; meas++) {
                MCSweeps (NStep, MySpinConf, L_tilda, T_tilda);
                // do_measurement (MySpinConf, L, E[meas], E2[meas], E4[meas], M[meas], M2[meas], M4[meas]);
                E[meas] = Energy (MySpinConf, L_tilda);
                E2[meas] = (E[meas])*(E[meas]);
                E4[meas] = (E2[meas])*(E2[meas]);
                M[meas] = Magnetization (MySpinConf, L_tilda);
                M2[meas] = (M[meas])*(M[meas]);
                M4[meas] = (M2[meas])*(M2[meas]);

                // cout << "The Energy when measurement = " << meas << " is, " << E[meas] << endl;
                // cout << "The magnetization = " << M[meas] << endl;

                // file << M[meas] << "\n";
            }

            // performing printout
            // double E_avrg = printout (Nmeas, E, E2, E4, M, M2, M4);

            // cout << "Printing out results:" << endl;
            double E_avg = 0;
            for (int meas = 1; meas <= Nmeas; meas++) {
                E_avg += E[meas];    
            }
            E_avg = E_avg/Nmeas;

            cout << E_avg << " is the avg energy." << "\n";

            double E2_avg = 0;
            for (int meas = 1; meas <= Nmeas; meas++) {
                E2_avg += E2[meas];    
            }
            E2_avg = E2_avg/Nmeas;

            double M_avg = 0;
            for (int meas = 1; meas <= Nmeas; meas++) {
                M_avg += M[meas];    
            }
            M_avg = M_avg/Nmeas;

            cout << M_avg << " is the avg magnetization." << "\n";

            double M2_avg = 0;
            for (int meas = 1; meas <= Nmeas; meas++) {
                M2_avg += M2[meas];    
            }
            M2_avg = M2_avg/Nmeas;

            double M4_avg = 0;
            for (int meas = 1; meas <= Nmeas; meas++) {
                M4_avg += M4[meas];    
            }
            M4_avg = M4_avg/Nmeas;

            // Defining variable chi
            double chi = (M2_avg)/(T*N); 
            // Defining variable C_v
            double heatCapacity = (E2_avg - (E_avg*E_avg))/((T*T)*N); 
            // Defining variable g
            double g = (0.5)*(3 - (M4_avg/(M2_avg*M2_avg)));

            // Writing to the files 
            file1 << E_avg << "\n";
            file2 << chi << "\n";
            file3 << heatCapacity << "\n";

            // if (L_tilda == 8) {
            //     file4 << g << "\n";
            // }
            // else if (L_tilda == 12) {
            //     file5 << g << "\n";
            // }
            // else {
            //     file6 << g << "\n";
            // }
        }   
        // Closing the files
        file1.close(); file2.close(); file3.close();
    }

    // Closing the files
    // file4.close(), file5.close(), file6.close();

    cout << "The END of SQIsing!" << endl;
} 

// main function 
int main () {
    int L;
    cout << "Enter the value of initial L (size of the lattice): " << endl;
    cin >> L;

    double T;
    cout << "Enter the value of initial T (temperature): " << endl;
    cin >> T;

    SQIsing (L, T, 15, 6, 50);

    cout << "Falling back to main..." << endl;
    return 0;
}