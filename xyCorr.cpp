#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <complex>

using namespace std;

double J = 1.0;

inline double random01() {
    return rand() / (double)RAND_MAX;
}

class SpinSystem {
public:
    vector<double> spins; // Array to store spin angles
    vector<double> horizontal_corr;
    vector<double> vertical_corr;
    vector<double> first_diagonal_corr;
    vector<double> second_diagonal_corr;
    int N;
    int L;
    int connectivity;
    double temp;
    double beta;
    double theta_max;
    double h;
    int tau;
    double phiForH;
    int theta_max_INT;
    double EPR;

    SpinSystem(int arg1, double arg2, double arg3, int arg4, double arg5, int arg6) : L(arg1), temp(arg2), theta_max_INT(arg3), connectivity(arg4), h(arg5), tau(arg6){
        N = L * L;
        beta = 1 / temp; 
        theta_max = 2 * M_PI * theta_max_INT / 360;

        spins.resize(N);
        neighbors.resize(N);
        geometricNeighbors.resize(N);
        angles.resize(connectivity);
        horizontal_pairs.resize(L/2+1);
        vertical_pairs.resize(L/2+1);
        horizontal_corr.resize(L/2);
        vertical_corr.resize(L/2);
        first_diagonal_pairs.resize(L/2+1);
        second_diagonal_pairs.resize(L/2+1);
        first_diagonal_corr.resize(L/2);
        second_diagonal_corr.resize(L/2);
        

        initializeSpins();
        initializeAngles();
        compute_phiForH();
        initialize_pairs();
     
        // Calculate neighbors based on lattice type
        if (connectivity == 4) {
            // Square lattice
            calculateSquareGeometricNeighbors();
        } else if (connectivity == 6) {
            // Hexagonal lattice
            calculateHexagonalGeometricNeighbors();
        }
  
        initializeNeighbors();
       

    }

//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------
//-------------Funzioni per il modello non reciproco-------------------------------------
//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------


  
    void MC_NR_WITH_H(int MCSteps) {
        for (int step = 0; step < MCSteps; step++) {
            for (int i = 0; i < N; i++) {
                // Choose a random spin to update
                int randomSpin = rand() % N;

                // Generate a random angle change between -PI and PI
                double deltaTheta = M_PI * (2.0 * random01() - 1.0);

                // Calculate the change in energy
                double deltaE = energyDiff_NR_WITH_H(randomSpin, deltaTheta, step);
            //  std::cout << "deltaE " << deltaE << endl;
               // spins[randomSpin] += deltaTheta;

                // Metropolis acceptance criterion
                if ( random01() > (1 - tanh((-beta * deltaE) / 2)) / 2  ) {
                    // Reject the change
                //    std::cout << "rejected" << endl;
                    spins[randomSpin] += deltaTheta;
                } else {
                   // std::cout << "Accepted" << endl;
                }
                if (spins[randomSpin] < 0) {
                    spins[randomSpin] += 2 * M_PI;
                }
                if (spins[randomSpin] > 2 * M_PI) {
                    spins[randomSpin] -= 2 * M_PI;
                }
                calculateNeighborForSite(randomSpin, spins[randomSpin], connectivity);
            }
        }
    }  



    void MC_NR(int MCSteps) {
        EPR = 0.0;
        for (int step = 0; step < MCSteps; step++) {
            for (int i = 0; i < N; i++) {
                // Choose a random spin to update
                int randomSpin = rand() % N;

                // Generate a random angle change between -PI and PI
                double deltaTheta = M_PI * (2.0 * random01() - 1.0);

                // Calculate the change in energy
                double deltaE = energyDiff_NR(randomSpin, deltaTheta);
              //  std::cout << "deltaE " << deltaE << endl;
               // spins[randomSpin] += deltaTheta;

                // Metropolis acceptance criterion
                if ( random01() > (1 - tanh((-beta * deltaE) / 2)) / 2  ) {
                    // Reject the change
                //    std::cout << "rejected" << endl;
                    spins[randomSpin] += deltaTheta;
                    EPR += log( (1-tanh(deltaE*beta/2)) / (1-tanh(-deltaE*beta/2)) );
                } else {
                   // std::cout << "Accepted" << endl;
                }
                if (spins[randomSpin] < 0) {
                    spins[randomSpin] += 2 * M_PI;
                }
                if (spins[randomSpin] > 2 * M_PI) {
                    spins[randomSpin] -= 2 * M_PI;
                }
                calculateNeighborForSite(randomSpin, spins[randomSpin], connectivity);
            }
        }
        EPR = EPR / (MCSteps * N);
    }

            double calculateEnergy_NR() {
        double energy = 0.0;
        for (int i = 0; i < N; i++) {
            double spin_i = spins[i];
            for (const int& j : neighbors[i]) {
                double spin_j = spins[j];
                double cos_diff = cos(spin_i - spin_j);
                energy += -J * cos_diff;
            }
        }
        return energy;
    }


    // Function to calculate the energy difference for the XY non reciprocal
    double energyDiff_NR(const int& i, const double& deltaTheta) {
        double energyDiff = 0.0;
        double spin_i = spins[i];
        double spin_iPrime = spin_i + deltaTheta;
        vector<int> neighbors_i = neighbors[i];
        for (const int& j : neighbors_i) {
            double spin_j = spins[j];
            double cos_diff = cos(spin_i - spin_j);
            energyDiff += J * cos_diff; 
        }
        calculateNeighborForSite(i, spin_iPrime, connectivity);
        neighbors_i = neighbors[i];
        for (const int& j : neighbors_i) {
            double spin_j = spins[j];
            double cos_diff = cos(spin_iPrime - spin_j);  
            energyDiff += - J * cos_diff; 
        }
        return energyDiff;
    }

    // Function to calculate the energy difference for the XY non reciprocal
    double energyDiff_NR_WITH_H(const int& i, const double& deltaTheta, int& time) {
        double heff = h * exp(- time / tau);
        //cout << phiForH << endl;
        double energyDiff = 0.0;
        double spin_i = spins[i];
        double spin_iPrime = spin_i + deltaTheta;
        vector<int> neighbors_i = neighbors[i];
        for (const int& j : neighbors_i) {
            double spin_j = spins[j];
            double cos_diff = cos(spin_i - spin_j);
            energyDiff += J * cos_diff; 
            energyDiff += - heff * cos(spin_i - phiForH); 
        }
        calculateNeighborForSite(i, spin_iPrime, connectivity);
        neighbors_i = neighbors[i];
        for (const int& j : neighbors_i) {
            double spin_j = spins[j];
            double cos_diff = cos(spin_iPrime - spin_j);  
            energyDiff += - J * cos_diff; 
            energyDiff += + heff * cos(spin_iPrime - phiForH); 
        }
        return energyDiff;
    }

//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------
//-------------Funzioni per il modello reciproco         ----------------------------
//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------


    // Wolff cluster algorithm for the reciprocal case (what we call M2) with magnetic field
    void wolff_R_WITH_H(int time){
        vector<int> cluster; cluster = {}; //initialize an empty list that has to contain the cluster
        double deltaTheta = M_PI * random01(); // Generate a random angle to rotate the spins from 0 to PI
        int randomSpin = rand() % N; // Choose a random spin from which to start the cluster
        cluster.push_back(randomSpin); // Insert the first spin in the cluster
        int counter = 0; // Dummy counter variable

        // At the first step counter = 0 and cluster.size = 1 so the statement is true. Every time I add a site to the cluster 
        // cluster.size increase by 1. Every time I address a site in order to possibly add its neighbors to the cluster, counter 
        // is increased by one; when I have addressed all the site of the cluster without having added sites to the cluster, the size
        // will be equal to the counter, and the process stops
        while (counter < cluster.size()) {
            int i = cluster[counter];
            double spin_i = spins[i];
            double spin_iPrime = 2 * deltaTheta - spin_i; // it si between 0 and 2PI
            spins[i] = spin_iPrime; // I already update and it actually stays updated forever cause it is in the cluster
            calculateNeighborForSite(i, spin_iPrime, connectivity); // I also change its effective neghbors (the ones it sees)
            // The candidates that can enter the cluster have to be chosen from all the geometric neighbors
            for (int&k : geometricNeighbors[i]) {
                auto it = std::find(cluster.begin(), cluster.end(), k);
                // If the neighbor I'm considering is not in the cluster yet, I consider the possibility to add it
                if (it == cluster.end()){
                    double spin_k = spins[k];
                    double energy = computeEnergyBond_R_WITH_H(spin_iPrime, spin_k, i, k, time); // Compute the first contribution to the energy difference
                    double spin_kPrime = 2 * deltaTheta - spin_k;
                    calculateNeighborForSite(k, spin_kPrime, connectivity); // Upadate of neighbors of k... if the site is rejected I have to undo this operation
                    energy = energy - computeEnergyBond_R_WITH_H(spin_iPrime, spin_kPrime, i, k, time);
                    energy = - energy; // mi torna meglio pensarla così
                    if ( energy < 0 && random01() < (1 - exp((+ beta * energy)))  ) {  //spero qui di non aver sbagliato un segno... al massimo cambio dopo...
                        cluster.push_back(k); // The site is accepted: I add to to the cluster
                    } else {
                        calculateNeighborForSite(k, spin_k, connectivity); // I undo the operation if the site is rejected
                    }
                }           
            }
            counter++; // I add the counter of one unit. Note that cluster.size could have been incresed from 0 to connectivity.
        }

    }

    // Wolff cluster algorithm for the reciprocal case (what we call M2) without magnetic field
    void wolff_R(){
        vector<int> cluster; cluster = {}; //initialize an empty list that has to contain the cluster
        double deltaTheta = M_PI * random01(); // Generate a random angle to rotate the spins from 0 to PI
        int randomSpin = rand() % N; // Choose a random spin from which to start the cluster
        cluster.push_back(randomSpin); // Insert the first spin in the cluster
        int counter = 0; // Dummy counter variable

        // At the first step counter = 0 and cluster.size = 1 so the statement is true. Every time I add a site to the cluster 
        // cluster.size increase by 1. Every time I address a site in order to possibly add its neighbors to the cluster, counter 
        // is increased by one; when I have addressed all the site of the cluster without having added sites to the cluster, the size
        // will be equal to the counter, and the process stops
        while (counter < cluster.size()) {
            int i = cluster[counter];
            double spin_i = spins[i];
            double spin_iPrime = 2 * deltaTheta - spin_i; // it si between 0 and 2PI
            spins[i] = spin_iPrime; // I already update and it actually stays updated forever cause it is in the cluster
            calculateNeighborForSite(i, spin_iPrime, connectivity); // I also change its effective neghbors (the ones it sees)
            // The candidates that can enter the cluster have to be chosen from all the geometrice neighbors
            for (int&k : geometricNeighbors[i]) {
                auto it = std::find(cluster.begin(), cluster.end(), k);
                // If the neighbor I'm considering is not in the cluster yet, I consider the possibility to add it
                if (it == cluster.end()){
                    double spin_k = spins[k];
                    double energy = computeEnergyBond_R(spin_iPrime, spin_k, i, k); // Compute the first contribution to the energy difference
                    double spin_kPrime = 2 * deltaTheta - spin_k;
                    calculateNeighborForSite(k, spin_kPrime, connectivity); // Upadate of neighbors of k... if the site is rejected I have to undo this operation
                    energy = energy - computeEnergyBond_R(spin_iPrime, spin_kPrime, i, k);
                    energy = - energy;
                    if ( energy < 0 && random01() < (1 - exp((+ beta * energy)))  ) {   //spero qui di non aver sbagliato un segno... al massimo cambio dopo...
                        cluster.push_back(k); // The site is accepted: I add to to the cluster
                    } else {
                        calculateNeighborForSite(k, spin_k, connectivity); // I undo the operation if the site is rejected
                    }
                }           
            }
            counter++; // I add the counter of one unit. Note that cluster.size could have been incresed from 0 to connectivity.
        }

    }

    // The energy relatice to a bond in the reciprocal case (M2) with magnetic field
    double computeEnergyBond_R_WITH_H(double spin_i, double spin_k, int i, int k, int time){
        double JR = 0.5;
        vector<int> neighbors_i = neighbors[i];
        vector<int> neighbors_k = neighbors[k];
        double heff = h * exp(- time / tau);
        double energy = 0.0;
        auto it1 = std::find(neighbors_i.begin(), neighbors_i.end(), k);
        if (it1 != neighbors_i.end()){
            double energy =+ -JR * cos(spin_i - spin_k);
        }
        auto it2 = std::find(neighbors_k.begin(), neighbors_k.end(), i);
        if (it2 != neighbors_k.end()){
            double energy =+ -JR * cos(spin_i - spin_k);
        }
        energy =+ + heff * (cos(spin_i - phiForH) + cos(spin_k - phiForH));
        return energy;
    }

     // The energy relatice to a bond in the reciprocal case (M2) without magnetic field
    double computeEnergyBond_R(double spin_i, double spin_k, int i, int k){
        vector<int> neighbors_i = neighbors[i];
        vector<int> neighbors_k = neighbors[k];
        double JR = 0.5;
        // double heff = h * exp(- time / tau);
        double energy = 0.0;
        auto it1 = std::find(neighbors_i.begin(), neighbors_i.end(), k);
        if (it1 != neighbors_i.end()){
            double energy =+ -JR * cos(spin_i - spin_k);
        }
        auto it2 = std::find(neighbors_k.begin(), neighbors_k.end(), i);
        if (it2 != neighbors_k.end()){
            double energy =+ -JR * cos(spin_i - spin_k);
        }
        //energy =+ - heff * (cos(spin_i - phiForH) + cos(spin_k - phiForH));
        return energy;
    }

    void MC_R_WITH_H(int MCSteps) {
        for (int step = 0; step < MCSteps; step++) {
            for (int i = 0; i < N; i++) {
                // Choose a random spin to update
                int randomSpin = rand() % N;

                // Generate a random angle change between -PI and PI
                double deltaTheta = M_PI * (2.0 * random01() - 1.0);

                // Calculate the change in energy
                double deltaE = energyDiff_R_WITH_H(randomSpin, deltaTheta, step);
              //  std::cout << "deltaE " << deltaE << endl;
               // spins[randomSpin] += deltaTheta;

                // Metropolis acceptance criterion
                if ( random01() > (1 - tanh((-beta * deltaE) / 2)) / 2  ) {
                    // Reject the change
                //    std::cout << "rejected" << endl;
                    spins[randomSpin] += deltaTheta;
                } else {
                   // std::cout << "Accepted" << endl;
                }
                if (spins[randomSpin] < 0) {
                    spins[randomSpin] += 2 * M_PI;
                }
                if (spins[randomSpin] > 2 * M_PI) {
                    spins[randomSpin] -= 2 * M_PI;
                }
                calculateNeighborForSite(randomSpin, spins[randomSpin], connectivity);
            }
            wolff_R_WITH_H(step);
        }
    }      

    void MC_R(int MCSteps) {
        EPR = 0.0;
        for (int step = 0; step < MCSteps; step++) {
            for (int i = 0; i < N; i++) {
                // Choose a random spin to update
                int randomSpin = rand() % N;

                // Generate a random angle change between -PI and PI
                double deltaTheta = M_PI * (2.0 * random01() - 1.0);

                // Calculate the change in energy
                double deltaE = energyDiff_R(randomSpin, deltaTheta);
              //  std::cout << "deltaE " << deltaE << endl;
               // spins[randomSpin] += deltaTheta;

                // Metropolis acceptance criterion
                if ( random01() > (1 - tanh((-beta * deltaE) / 2)) / 2  ) {
                    // Reject the change
                //    std::cout << "rejected" << endl;
                    spins[randomSpin] += deltaTheta;
                    EPR += log( (1-tanh(deltaE*beta/2)) / (1-tanh(-deltaE*beta/2)) );
                } else {
                   // std::cout << "Accepted" << endl;
                }
                if (spins[randomSpin] < 0) {
                    spins[randomSpin] += 2 * M_PI;
                }
                if (spins[randomSpin] > 2 * M_PI) {
                    spins[randomSpin] -= 2 * M_PI;
                }
                calculateNeighborForSite(randomSpin, spins[randomSpin], connectivity);
            }
            wolff_R();
        }
        EPR = EPR / (MCSteps * N);
    }  

        double calculateEnergy_R() {
        double JR = 0.5;
        double energy = 0.0;
        for (int i = 0; i < N; i++) {
            double spin_i = spins[i];
            for (const int& j : neighbors[i]) {
                double spin_j = spins[j];
                double cos_diff = cos(spin_i - spin_j);
                energy += -JR * cos_diff;
            }
        }
        return energy;
    }

    // Function to calculate the energy difference for the XY non reciprocal
    double energyDiff_R(const int& i, const double& deltaTheta) {
        double JR = 0.5;
        double energyDiff = 0.0;
        double spin_i = spins[i];
        double spin_iPrime = spin_i + deltaTheta;
        vector<int> neighbors_i = neighbors[i];
        for (const int& j : neighbors_i) {
            double spin_j = spins[j];
            double cos_diff = cos(spin_i - spin_j);
            energyDiff += JR * cos_diff; 
        }
        calculateNeighborForSite(i, spin_iPrime, connectivity);
        neighbors_i = neighbors[i];
        for (const int& j : neighbors_i) {
            double spin_j = spins[j];
            double cos_diff = cos(spin_iPrime - spin_j);  
            energyDiff += - JR * cos_diff; 
        }
        for (const int& k : geometricNeighbors[i]) {
            vector<int> neighbors_k = neighbors[k];
            double spin_k = spins[k];
            auto it = std::find(neighbors_k.begin(), neighbors_k.end(), i);
            if (it != neighbors_k.end()) {
                double cos_diff = cos(spin_k - spin_i) - cos(spin_k - spin_iPrime);
                energyDiff += JR * cos_diff;
            }
        }
        return energyDiff;
    }

    double energyDiff_R_WITH_H(const int& i, const double& deltaTheta, int& time) {
        double JR = 0.5;
        double heff = h * exp(- time / tau);
        double energyDiff = 0.0;
        double spin_i = spins[i];
        double spin_iPrime = spin_i + deltaTheta;
        vector<int> neighbors_i = neighbors[i];
        for (const int& j : neighbors_i) {
            double spin_j = spins[j];
            double cos_diff = cos(spin_i - spin_j);
            energyDiff += JR * cos_diff; 
            energyDiff += - heff * cos(spin_i - phiForH); 
        }
        calculateNeighborForSite(i, spin_iPrime, connectivity);
        neighbors_i = neighbors[i];
        for (const int& j : neighbors_i) {
            double spin_j = spins[j];
            double cos_diff = cos(spin_iPrime - spin_j);  
            energyDiff += - JR * cos_diff; 
            energyDiff += + heff * cos(spin_iPrime - phiForH); 
        }
        for (const int& k : geometricNeighbors[i]) {
            vector<int> neighbors_k = neighbors[k];
            double spin_k = spins[k];
            auto it = std::find(neighbors_k.begin(), neighbors_k.end(), i);
            if (it != neighbors_k.end()) {
                double cos_diff = cos(spin_k - spin_i) - cos(spin_k - spin_iPrime);
                energyDiff += JR * cos_diff;
            }
        }
        return energyDiff;
    }

//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------
//-------------Funzioni per il modello a Lazo         ----------------------------
//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------


    // Wolff cluster algorithm for the reciprocal case (what we call M2) with magnetic field
    void wolff_L_WITH_H(int time){
        vector<int> cluster; cluster = {}; //initialize an empty list that has to contain the cluster
        double deltaTheta = M_PI * random01(); // Generate a random angle to rotate the spins from 0 to PI
        int randomSpin = rand() % N; // Choose a random spin from which to start the cluster
        cluster.push_back(randomSpin); // Insert the first spin in the cluster
        int counter = 0; // Dummy counter variable

        // At the first step counter = 0 and cluster.size = 1 so the statement is true. Every time I add a site to the cluster 
        // cluster.size increase by 1. Every time I address a site in order to possibly add its neighbors to the cluster, counter 
        // is increased by one; when I have addressed all the site of the cluster without having added sites to the cluster, the size
        // will be equal to the counter, and the process stops
        while (counter < cluster.size()) {
            int i = cluster[counter];
            double spin_i = spins[i];
            double spin_iPrime = 2 * deltaTheta - spin_i; // it si between 0 and 2PI
            spins[i] = spin_iPrime; // I already update and it actually stays updated forever cause it is in the cluster
            calculateNeighborForSite(i, spin_iPrime, connectivity); // I also change its effective neghbors (the ones it sees)
            // The candidates that can enter the cluster have to be chosen from all the geometrice neighbors
            for (int&k : geometricNeighbors[i]) {
                auto it = std::find(cluster.begin(), cluster.end(), k);
                // If the neighbor I'm considering is not in the cluster yet, I consider the possibility to add it
                if (it == cluster.end()){
                    double spin_k = spins[k];
                    double energy = computeEnergyBond_L_WITH_H(spin_iPrime, spin_k, i, k, time); // Compute the first contribution to the energy difference
                    double spin_kPrime = 2 * deltaTheta - spin_k;
                    calculateNeighborForSite(k, spin_kPrime, connectivity); // Upadate of neighbors of k... if the site is rejected I have to undo this operation
                    energy = energy - computeEnergyBond_L_WITH_H(spin_iPrime, spin_kPrime, i, k, time);
                    energy = - energy;
                    if ( energy < 0 && random01() < (1 - exp((+ beta * energy)))  ) {   //spero qui di non aver sbagliato un segno... al massimo cambio dopo...
                        cluster.push_back(k); // The site is accepted: I add to to the cluster
                    } else {
                        calculateNeighborForSite(k, spin_k, connectivity); // I undo the operation if the site is rejected
                    }
                }           
            }
            counter++; // I add the counter of one unit. Note that cluster.size could have been incresed from 0 to connectivity.
        }

    }

    // Wolff cluster algorithm for the reciprocal case (what we call M2) without magnetic field
    void wolff_L(){
        vector<int> cluster; cluster = {}; //initialize an empty list that has to contain the cluster
        double deltaTheta = M_PI * random01(); // Generate a random angle to rotate the spins from 0 to PI
        int randomSpin = rand() % N; // Choose a random spin from which to start the cluster
        cluster.push_back(randomSpin); // Insert the first spin in the cluster
        int counter = 0; // Dummy counter variable

        // At the first step counter = 0 and cluster.size = 1 so the statement is true. Every time I add a site to the cluster 
        // cluster.size increase by 1. Every time I address a site in order to possibly add its neighbors to the cluster, counter 
        // is increased by one; when I have addressed all the site of the cluster without having added sites to the cluster, the size
        // will be equal to the counter, and the process stops
        while (counter < cluster.size()) {
            int i = cluster[counter];
            double spin_i = spins[i];
            double spin_iPrime = 2 * deltaTheta - spin_i; // it si between 0 and 2PI
            spins[i] = spin_iPrime; // I already update and it actually stays updated forever cause it is in the cluster
            calculateNeighborForSite(i, spin_iPrime, connectivity); // I also change its effective neghbors (the ones it sees)
            // The candidates that can enter the cluster have to be chosen from all the geometrice neighbors
            for (int&k : geometricNeighbors[i]) {
                auto it = std::find(cluster.begin(), cluster.end(), k);
                // If the neighbor I'm considering is not in the cluster yet, I consider the possibility to add it
                if (it == cluster.end()){
                    double spin_k = spins[k];
                    double energy = computeEnergyBond_L(spin_iPrime, spin_k, i, k); // Compute the first contribution to the energy difference
                    double spin_kPrime = 2 * deltaTheta - spin_k;
                    calculateNeighborForSite(k, spin_kPrime, connectivity); // Upadate of neighbors of k... if the site is rejected I have to undo this operation
                    energy = energy - computeEnergyBond_L(spin_iPrime, spin_kPrime, i, k);
                    energy = -energy;
                    if ( energy < 0 && random01() < (1 - exp((+ beta * energy)))  ) {   //spero qui di non aver sbagliato un segno... al massimo cambio dopo...
                        cluster.push_back(k); // The site is accepted: I add it to to the cluster
                    } else {
                        calculateNeighborForSite(k, spin_k, connectivity); // I undo the operation if the site is rejected
                    }
                }           
            }
            counter++; // I add the counter of one unit. Note that cluster.size could have been incresed from 0 to connectivity.
        }

    }

    // The energy relatice to a bond in the reciprocal case (M2) with magnetic field
    double computeEnergyBond_L_WITH_H(double spin_i, double spin_k, int i, int k, int time){
        //double JR = 0.5;
        vector<int> neighbors_i = neighbors[i];
        vector<int> neighbors_k = neighbors[k];
        double heff = h * exp(- time / tau);
        double energy = 0.0;
        auto it1 = std::find(neighbors_i.begin(), neighbors_i.end(), k);
        auto it2 = std::find(neighbors_k.begin(), neighbors_k.end(), i);
        if (it1 != neighbors_i.end() || it2 != neighbors_k.end()){
            double energy =- J * cos(spin_i - spin_k);
        }
        energy =+ + heff * (cos(spin_i - phiForH) + cos(spin_k - phiForH));
        return energy;
    }

     // The energy relatice to a bond in the reciprocal case (M2) without magnetic field
    double computeEnergyBond_L(double spin_i, double spin_k, int i, int k){
        vector<int> neighbors_i = neighbors[i];
        vector<int> neighbors_k = neighbors[k];
        // double heff = h * exp(- time / tau);
        double energy = 0.0;
        auto it1 = std::find(neighbors_i.begin(), neighbors_i.end(), k);
        auto it2 = std::find(neighbors_k.begin(), neighbors_k.end(), i);
        if (it1 != neighbors_i.end() || it2 != neighbors_k.end()){
            double energy =- J * cos(spin_i - spin_k);
        }
        //energy =+ - heff * (cos(spin_i - phiForH) + cos(spin_k - phiForH));
        return energy;
    }

        void MC_L_WITH_H(int MCSteps) {
        for (int step = 0; step < MCSteps; step++) {
            for (int i = 0; i < N; i++) {
                // Choose a random spin to update
                int randomSpin = rand() % N;

                // Generate a random angle change between -PI and PI
                double deltaTheta = M_PI * (2.0 * random01() - 1.0);

                // Calculate the change in energy
                double deltaE = energyDiff_L_WITH_H(randomSpin, deltaTheta, step);
              //  std::cout << "deltaE " << deltaE << endl;
               // spins[randomSpin] += deltaTheta;

                // Metropolis acceptance criterion
                if ( random01() > (1 - tanh((-beta * deltaE) / 2)) / 2  ) {
                    // Reject the change
                //    std::cout << "rejected" << endl;
                    spins[randomSpin] += deltaTheta;
                } else {
                   // std::cout << "Accepted" << endl;
                }
                if (spins[randomSpin] < 0) {
                    spins[randomSpin] += 2 * M_PI;
                }
                if (spins[randomSpin] > 2 * M_PI) {
                    spins[randomSpin] -= 2 * M_PI;
                }
                calculateNeighborForSite(randomSpin, spins[randomSpin], connectivity);
            }
            wolff_L_WITH_H(step);
        }
    }    

    void MC_L(int MCSteps) {
        EPR = 0.0;
        for (int step = 0; step < MCSteps; step++) {
            for (int i = 0; i < N; i++) {
                // Choose a random spin to update
                int randomSpin = rand() % N;

                // Generate a random angle change between -PI and PI
                double deltaTheta = M_PI * (2.0 * random01() - 1.0);

                // Calculate the change in energy
                double deltaE = energyDiff_L(randomSpin, deltaTheta);
              //  std::cout << "deltaE " << deltaE << endl;
               // spins[randomSpin] += deltaTheta;

                // Metropolis acceptance criterion
                if ( random01() > (1 - tanh((-beta * deltaE) / 2)) / 2  ) {
                    // Reject the change
                //    std::cout << "rejected" << endl;
                    spins[randomSpin] += deltaTheta;
                    EPR += log( (1-tanh(deltaE*beta/2)) / (1-tanh(-deltaE*beta/2)) );
                } else {
                   // std::cout << "Accepted" << endl;
                }
                if (spins[randomSpin] < 0) {
                    spins[randomSpin] += 2 * M_PI;
                }
                if (spins[randomSpin] > 2 * M_PI) {
                    spins[randomSpin] -= 2 * M_PI;
                }
                calculateNeighborForSite(randomSpin, spins[randomSpin], connectivity);
            }
            wolff_L();
        }
        EPR = EPR / (MCSteps * N);
    }   

        double calculateEnergy_L() {
        double energy = 0.0;
        for (int i = 0; i < N; i++) {
            double spin_i = spins[i];
            vector<int> neighbors_i = neighbors[i];
            for (int j = 0; j < connectivity / 2; j++) {
                int k = geometricNeighbors[i][j];
                double spin_k = spins[k];
                vector<int> neighbors_k = neighbors[k];
                auto it1 = std::find(neighbors_k.begin(), neighbors_k.end(), i);
                auto it2 = std::find(neighbors_i.begin(), neighbors_i.end(), k);
                if (it1 != neighbors_k.end() || it2 != neighbors_i.end()) {
                    double cos_diff = cos(spin_i - spin_k);
                    energy += -J * cos_diff; 
                }
            }
        }
        return energy;
    }

    // Function to calculate the energy difference for the XY non reciprocal
    double energyDiff_L(const int& i, const double& deltaTheta) {
        double energyDiff = 0.0;
        double spin_i = spins[i];
        double spin_iPrime = spin_i + deltaTheta;
        vector<int> neighbors_i = neighbors[i];
        for (const int& j : neighbors_i) {
            double spin_j = spins[j];
            double cos_diff = cos(spin_i - spin_j);
            energyDiff += J * cos_diff; 
        }
        calculateNeighborForSite(i, spin_iPrime, connectivity);
        vector<int> neighbors_iPrime = neighbors[i];
        for (const int& j : neighbors_iPrime) {
            double spin_j = spins[j];
            double cos_diff = cos(spin_iPrime - spin_j);  
            energyDiff += - J * cos_diff; 
        }
        for (const int& k : geometricNeighbors[i]) {
            vector<int> neighbors_k = neighbors[k];
            double spin_k = spins[k];
            auto it = std::find(neighbors_k.begin(), neighbors_k.end(), i);
            auto it1 = std::find(neighbors_i.begin(), neighbors_i.end(), k);
            auto it2 = std::find(neighbors_iPrime.begin(), neighbors_iPrime.end(), k);
            if (it != neighbors_k.end() && it1 == neighbors_i.end()) {
                double cos_diff = cos(spin_k - spin_i);
                energyDiff += J * cos_diff;
            }
            if (it != neighbors_k.end() && it2 == neighbors_iPrime.end()) {
                double cos_diff = cos(spin_k - spin_iPrime);
                energyDiff += - J * cos_diff;
            }
        }
        return energyDiff;
    }

    // Function to calculate the energy difference for the XY non reciprocal
    double energyDiff_L_WITH_H(const int& i, const double& deltaTheta, int& time) {
        double heff = h * exp(- time / tau);
        double energyDiff = 0.0;
        double spin_i = spins[i];
        double spin_iPrime = spin_i + deltaTheta;
        vector<int> neighbors_i = neighbors[i];
        for (const int& j : neighbors_i) {
            double spin_j = spins[j];
            double cos_diff = cos(spin_i - spin_j);
            energyDiff += J * cos_diff; 
            energyDiff += - heff * cos(spin_i - phiForH);
        }
        calculateNeighborForSite(i, spin_iPrime, connectivity);
        vector<int> neighbors_iPrime = neighbors[i];
        for (const int& j : neighbors_iPrime) {
            double spin_j = spins[j];
            double cos_diff = cos(spin_iPrime - spin_j);  
            energyDiff += - J * cos_diff; 
            energyDiff += + heff * cos(spin_iPrime - phiForH);
        }
        for (const int& k : geometricNeighbors[i]) {
            vector<int> neighbors_k = neighbors[k];
            double spin_k = spins[k];
            auto it = std::find(neighbors_k.begin(), neighbors_k.end(), i);
            auto it1 = std::find(neighbors_i.begin(), neighbors_i.end(), k);
            auto it2 = std::find(neighbors_iPrime.begin(), neighbors_iPrime.end(), k);
            if (it != neighbors_k.end() && it1 == neighbors_i.end()) {
                double cos_diff = cos(spin_k - spin_i);
                energyDiff += J * cos_diff;
            }
            if (it != neighbors_k.end() && it2 == neighbors_iPrime.end()) {
                double cos_diff = cos(spin_k - spin_iPrime);
                energyDiff += - J * cos_diff;
            }
        }
        return energyDiff;
    }

//-------------------------------------------------------------------------------------------    
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------

    void PrintSpinAngles() const {
        for (int i = 0; i < N; i++) {
            cout << "Spin " << i << ": " << spins[i] << endl;
        }
    }

    // Function to compute the magnetization of the XY model
    double computeMagnetization() {
        double mx = 0.0;
        double my = 0.0;

        for (const double& angle : spins) {
            mx += cos(angle);
            my += sin(angle);
        }

        mx /= N;
        my /= N;

        return sqrt(mx * mx + my * my);
    }

    // Quantity I need to compute the helicity modulus
    double computeE(){
        double eValue = 0.0;
        for (int i = 0; i < N; i++) {
            double spin_i = spins[i];
            double spin_j = spins[geometricNeighbors[i][0]];
            eValue += cos(spin_i - spin_j);
        }
        return eValue / N;
    }
    
    // Quantity I need to compute the helicity modulus
     double computeS(){
        double sValue = 0.0;
        for (int i = 0; i < N; i++) {
            double spin_i = spins[i];
            double spin_j = spins[geometricNeighbors[i][0]];
            sValue += sin(spin_i - spin_j);
        }
        return sValue / N;
    }

    // Quantity I need for the second-moment correlation length
     double compute_M_km_squared(){
        double M_km_squared;
        complex<double> imaginary_unit(0.0, 1.0);
        double km = 2 * M_PI / L; 
        complex<double> mx = 0.0; 
        complex<double> my = 0.0;  
        for (int i = 0; i < N; i++) {
            int xcoordinate = i % L;
            double xc = xcoordinate * 1.0;
            //cout << "questa è la xcoord " << xc << endl;
            double spin_i = spins[i];
            mx += cos(spin_i) * exp(imaginary_unit * km * xc);
            my += sin(spin_i) * exp(imaginary_unit * km * xc);
            //cout << "questo è il mx " << mx << endl; 
        }
        mx /= N; my /= N; 
        M_km_squared = pow(abs(mx),2) + pow(abs(my), 2);
        return M_km_squared;
    }

    // funzione per la funzione di correlazione nel reticolo quadrato
    void correlation_squareLattice() {
        for (int l = 0; l < L/2; l++) {
            horizontal_corr[l] = 0.0;
            vertical_corr[l] = 0.0;
            first_diagonal_corr[l] = 0.0;
            second_diagonal_corr[l] = 0.0;   

            for (int i = 0; i < N; i++) {
                double spin_a = spins[horizontal_pairs[l][i][0]];
                double spin_b = spins[horizontal_pairs[l][i][1]];
                double spin_c = spins[vertical_pairs[l][i][0]];
                double spin_d = spins[vertical_pairs[l][i][1]];
                double spin_e = spins[first_diagonal_pairs[l][i][0]];
                double spin_f = spins[first_diagonal_pairs[l][i][1]];
                double spin_g = spins[second_diagonal_pairs[l][i][0]];
                double spin_h = spins[second_diagonal_pairs[l][i][1]];
                horizontal_corr[l] += cos(spin_a - spin_b);
                vertical_corr[l] += cos(spin_c - spin_d);
                first_diagonal_corr[l] += cos(spin_e - spin_f);
                second_diagonal_corr[l] += cos(spin_g - spin_h);
            }
            horizontal_corr[l] /= N;
            vertical_corr[l] /= N;
            first_diagonal_corr[l] /= N;
            second_diagonal_corr[l] /= N;
        }
    }

//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------        
//------------------------funzioni private           ----------------------------------------
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
   

private:
   // vector<double> spins; // Array to store spin angles
    vector<vector<int>> neighbors; // Array to store neighbors for each spin
    vector<vector<int>> geometricNeighbors; // Array to store geometric neighbors for each spin
    vector<double> angles;
    vector<vector<vector<int>>> vertical_pairs; 
    vector<vector<vector<int>>> horizontal_pairs; 
    vector<vector<vector<int>>> first_diagonal_pairs;
    vector<vector<vector<int>>> second_diagonal_pairs;

    void initialize_pairs() {
        horizontal_pairs = {};
        vertical_pairs = {};
        first_diagonal_pairs = {}; /////  it's this \   //////
        second_diagonal_pairs = {}; /////  it's this /    //////
        for (int l = 0; l <= L/2; l++) {
            cout << l << endl;
            for (int y = 0; y < L; y++) {
                for (int x = 0; x < L; x++) {
                    int siteIndex = y * L + x;
                    int siteIndexHor = y * L + x + l; 
                    if (x + l > L - 1) {
                        siteIndexHor = siteIndexHor - L; 
                    }
                    int siteIndexVer = (y + l) * L + x; 
                    if (y + l > L - 1) {
                        siteIndexVer = siteIndexVer - (L * L); 
                    }

                    int siteIndexFirstDiag = (y + l) * L + x + l;
                    if ( y + l > L - 1 && x + l > L - 1 ) {
                        siteIndexFirstDiag = siteIndexFirstDiag - L - (L * L); 
                    }
                    else if (y + l > L - 1) {
                        siteIndexFirstDiag = siteIndexFirstDiag - (L * L); 
                    }
                    else if (x + l > L - 1) {
                        siteIndexFirstDiag = siteIndexFirstDiag - L; 
                    }

                    int siteIndexSecondDiag = (y + l) * L + x - l;
                    if ( y + l > L - 1 && x - l < 0 ) {
                        siteIndexSecondDiag = siteIndexSecondDiag + L - (L * L); 
                    }
                    else if (y + l > L - 1) {
                        siteIndexSecondDiag = siteIndexSecondDiag - (L * L); 
                    }
                    else if (x + l < 0) {
                        siteIndexSecondDiag = siteIndexSecondDiag + L; 
                    }

                    vertical_pairs[l].push_back({siteIndex, siteIndexVer});
                    horizontal_pairs[l].push_back({siteIndex, siteIndexHor});
                    first_diagonal_pairs[l].push_back({siteIndex, siteIndexFirstDiag});
                    second_diagonal_pairs[l].push_back({siteIndex, siteIndexSecondDiag});
                }
            }
        }
    }

    
    
    

    void compute_phiForH() {
        if (connectivity == 4) {
            phiForH = M_PI / 8;
        }
        if (connectivity == 6) {
            phiForH = M_PI / 12;
        }
        //phiForH = 0.0;
    }

    void initializeAngles(){
        for (int i = 0; i < connectivity; i++){
          //  std::cout << i << std::endl;
            angles[i] = i * 2 * M_PI / connectivity;
            //std::cout << angles[i] << std::endl;
            //std::cout << "h00" << std::endl;

        }
       // std::cout << "finw" << std::endl;
    }

    void initializeSpins() {
        for (int i = 0; i < N; i++) {
            spins[i] = 2 * M_PI * random01();
        }
    }


    // Function to calculate square geometric neighbors with periodic boundary conditions
    void calculateSquareGeometricNeighbors() {
        const int dx[] = {1, 0, -1, 0}; // Displacement in x-direction for neighbors
        const int dy[] = {0, -1, 0, 1}; // Displacement in y-direction for neighbors

        for (int y = 0; y < L; y++) {
            for (int x = 0; x < L; x++) {
                int siteIndex = y * L + x;
                vector<int> siteNeighbors;

                for (int i = 0; i < 4; i++) {
                    int newX = (x + dx[i] + L) % L; // Apply periodic boundary conditions
                    int newY = (y + dy[i] + L) % L; // Apply periodic boundary conditions
                    int neighborIndex = newY * L + newX;
                    siteNeighbors.push_back(neighborIndex);
                }

                geometricNeighbors[siteIndex] = siteNeighbors;
            }
        }
    }


    // Function to calculate hexagonal geometric neighbors with periodic boundary conditions
    void calculateHexagonalGeometricNeighbors() {
        const int dx_even[] = {1, 1, 0, -1, 0, 1}; // Displacement in x-direction for even rows
        const int dx_odd[] = {1, 0, -1, -1, -1, 0}; // Displacement in x-direction for odd rows
        const int dy[] = {0, -1, -1, 0, 1, 1}; // Displacement in y-direction for neighbors

        for (int y = 0; y < L; y++) {
            for (int x = 0; x < L; x++) {
                int siteIndex = y * L + x;
                vector<int> siteNeighbors;

                for (int i = 0; i < 6; i++) {
                    int newX;
                    if (y % 2 == 0){
                        newX = (x + dx_even[i]+L) % L; // Apply periodic boundary conditions
                    }
                    else{
                        newX = (x + dx_odd[i]+L) % L; // Apply periodic boundary conditions
                    }
                    int newY = (y + dy[i]+L) % L; // Apply periodic boundary conditions
                    int neighborIndex = newY * L + newX;
                    siteNeighbors.push_back(neighborIndex);
                }

                geometricNeighbors[siteIndex] = siteNeighbors;
            }
        }
    }

    void calculateNeighborForSite(const int& site, const double& spin, const int& connectivity){
        int siteNumber = site;
        //std::cout << "hui"; 
        double theta = spin;
        neighbors[siteNumber] = {}; 
        for (int i = 0; i < connectivity; i++){
            double angle = angles[i];
        //    std::cout << theta << endl;
          //  std::cout << abs(angle - theta) << endl;
            //std::cout <<  theta_max<< endl;
            if (abs(angle - theta) < theta_max/2 || abs(angle - theta + 2 * M_PI) < theta_max/2 || abs(angle - theta - 2 * M_PI) < theta_max/2){
                neighbors[siteNumber].push_back(geometricNeighbors[siteNumber][i]);
              //  std::cout << "hei" << endl;  
            //    std::cout << neighbors[siteNumber][0] << endl; 
            }   else {
           // std::cout << "h0i" << endl;       
            }
        }
    }
    

    void initializeNeighbors() {
        for (int i = 0; i < N; i++) {
            calculateNeighborForSite(i, spins[i], connectivity);
        }
    }

    

};



   int main(int argc, char *argv[]) {
    srand(42); // Seed the random number generator with current time

    // Check if an input file path is provided as a command-line argument
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <input_file_path>" << endl;
        return 1;
    }

    // Open the input file
    ifstream inputFile(argv[1]);
    if (!inputFile.is_open()) {
        cerr << "Error: Could not open input file." << endl;
        return 1;
    }

    int L;
    double thetaM;
    double temperature;
    int connectivity;
    int MCstep1;
    int MCstep2;
    int modelType;
    string magnetizationName;
    string energyName;
    string latticeName;
    string EPRName;
    string verName;
    string horName;
    string firstDiagName;
    string secondDiagName;
    string eValueName;
    string sValueName;
    string mkValueName;
    
    inputFile >> L;
    inputFile >> thetaM; 
    inputFile >> temperature;   
    inputFile >> modelType;
    inputFile >> connectivity;
    inputFile >> MCstep1;
    inputFile >> MCstep2;
    inputFile >> magnetizationName;
    inputFile >> energyName;
    inputFile >> latticeName;
    inputFile >> EPRName;
    inputFile >> eValueName;
    inputFile >> sValueName;
    inputFile >> mkValueName;

    double theta_max = 1.0 * thetaM / 360 * 2 * M_PI;
    double magnetization;
    double energy;

    ofstream magnFile(magnetizationName);
    ofstream eneFile(energyName);
    ofstream latticeFile(latticeName);
    ofstream EPRFile(EPRName);
    //ofstream horFile(horName);
    //ofstream verFile(verName);
    //ofstream firstDiagFile(firstDiagName);
    //ofstream secondDiagFile(secondDiagName);
    ofstream eValueFile(eValueName);
    ofstream sValueFile(sValueName);
    ofstream mkValueFile(mkValueName);

    int h; int tau;

    SpinSystem system(L, temperature, thetaM, connectivity, h = 1.0, tau=30000);

    if (modelType == 1) {
        if (connectivity == 6) {
            if (thetaM > 61.0){
                system.MC_NR_WITH_H(100000);
            } else {
                system.MC_NR(100000);
            }
        }
        if (connectivity == 4) {
            if (thetaM > 91.0){
                system.MC_NR_WITH_H(100000);
            } else {
                system.MC_NR(100000);
            }
        }

        for (int i=0; i < MCstep1; i++) {
            system.MC_NR(MCstep2);
            energy = system.calculateEnergy_NR() / pow(L, 2);
            magnFile << system.computeMagnetization() << endl;
            eneFile << energy << endl;
            EPRFile << system.EPR << endl;
            eValueFile << system.computeE() << endl;
            sValueFile << system.computeS() << endl;
            mkValueFile << system.compute_M_km_squared() << endl;
         
        }

        
    }

    else if (modelType == 2) {
        if (connectivity == 6) {
            if (thetaM > 61.0){
                system.MC_R_WITH_H(100000);
            } else {
                system.MC_R(100000);
            }
        }
        if (connectivity == 4) {
            if (thetaM > 91.0){
                system.MC_R_WITH_H(100000);
            } else {
                system.MC_R(100000);
            }
        }

        for (int i=0; i < MCstep1; i++) {
            system.MC_R(MCstep2);
            energy = system.calculateEnergy_R() / pow(L, 2);
            magnFile << system.computeMagnetization() << endl;
            eneFile << energy << endl;
            EPRFile << system.EPR << endl;
            eValueFile << system.computeE() << endl;
            sValueFile << system.computeS() << endl;
            mkValueFile << system.compute_M_km_squared() << endl;
            
        }

        
    }

    else if (modelType == 3) {
        if (connectivity == 6) {
            if (thetaM > 61.0){
                system.MC_L_WITH_H(100000);
            } else {
                system.MC_L(100000);
            }
        }
        if (connectivity == 4) {
            if (thetaM > 91.0){
                system.MC_L_WITH_H(100000);
            } else {
                system.MC_L(100000);
            }
        }

        for (int i=0; i < MCstep1; i++) {
            system.MC_L(MCstep2);
            energy = system.calculateEnergy_L() / pow(L, 2);
            magnFile << system.computeMagnetization() << endl;
            eneFile << energy << endl;
            EPRFile << system.EPR << endl;
            system.correlation_squareLattice();
            eValueFile << system.computeE() << endl;
            sValueFile << system.computeS() << endl;
            mkValueFile << system.compute_M_km_squared() << endl;
            //cout << "questa è mkd " << system.compute_M_km_squared() << endl;
            
        }
    }

    return 0;
}
