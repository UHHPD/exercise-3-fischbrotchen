#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>



double poisson(double mu, int k) {
    // computes the Poisson probability for a given mean mu and k
    double poisson_prob = pow(mu, k) * exp(-mu) * 1/(tgamma(k+1));
    return poisson_prob;
}

int main() {
    using namespace std;

    // create vector with 11 entries to count n_i
    vector<int> zaehler(11);

    // open input and output files
    ifstream fin("datensumme.txt");
    ofstream fout("hist.txt");
    ofstream fout_poisson("histpoi.txt");

    // count occurences of n_i and store in vector
    int n_i;
    for(int i = 0 ; i < 234 ; ++i) {
        fin >> n_i;
        zaehler[n_i] += 1;
    }

    double mu = 3.11538;
    int N = 234;

    // print occurences from vector to terminal and write to output file
    for(unsigned int k = 0; k < zaehler.size(); k++) {
        std::cout << k << ":" << zaehler[k] << std::endl;
        fout << k << " " << zaehler[k] << std::endl;
        double poisson_expectation = N * poisson(mu, k);
        fout_poisson << k << " " << zaehler[k] << " " << poisson_expectation << std::endl;
    }

    // close input and output files
    fin.close();
    fout.close();
    fout_poisson.close();
}