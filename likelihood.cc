#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>


double poisson(double mu, int k) {
    // computes the Poisson probability for a given mean mu and k
    double poisson_prob = pow(mu, k) * exp(-mu) * 1/(tgamma(k+1));
    return poisson_prob;
}

double prob(std::vector<int> daten, double mu) {
    double likelihood;
    int k;
    for (int idx = 0; idx < daten.size(); idx++) {
        k = daten[idx];
        if (idx == 0) {
            likelihood = poisson(mu, k);
        }
        else {
            likelihood *= poisson(mu, k);
        }
    }
    return likelihood;
}

double prob_mu_is_k(std::vector<int> daten) {
    double likelihood;
    int k;
    for (int idx = 0; idx < daten.size(); idx++) {
        k = daten[idx];
        if (idx == 0) {
            likelihood = poisson(k, k);
        }
        else {
            likelihood *= poisson(k, k);
        }
    }
    return likelihood;
}

int main() {
    using namespace std;

    // store data from datensumme.txt in vector daten
    vector<int> daten;

    // vectors to store mu uncertainty
    vector<double> estimator_uncertainty;

    ifstream fin("datensumme.txt");

    int n_i;
    for(int i = 0 ; i < 234 ; ++i) {
        fin >> n_i;
        daten.push_back(n_i);
    }
    fin.close();

    // Print out likelihood for mean of sample mu = 3.11538
    double mu = 3.11538;
    double likelihood_mean_of_sample = prob(daten, mu);
    cout<<"Likelihood = "<<likelihood_mean_of_sample<<endl;

    // Scan likelihood for 0 < mu < 6, mu+=0.1
    ofstream fout("likelihood.txt");
    ofstream fout_nll("nll.txt");
    ofstream fout_delta_nll("deltanll.txt");
    double likelihood_mu = 0;
    double min_log_likelihood_diff = 0;
    double min_log_likelihood_mu;


    for (double mu_step = 0; mu_step <= 6; mu_step += 0.1) {
        likelihood_mu = prob(daten, mu_step);
        double nll_likelihood = -2*log(likelihood_mu);
        double delta_likelihood = nll_likelihood - (-2*log(likelihood_mean_of_sample));

        if (delta_likelihood < 1.0) {
            estimator_uncertainty.push_back(mu_step);
        }
        if (mu_step == 0) {
            min_log_likelihood_diff = delta_likelihood;
            min_log_likelihood_mu = mu_step;
        }
        else {
            if (delta_likelihood < min_log_likelihood_diff) {
                min_log_likelihood_diff = delta_likelihood;
                min_log_likelihood_mu = mu_step;
            }
        }

        fout << mu_step << " " << likelihood_mu << endl;
        fout_nll << mu_step << " " <<  likelihood_mu << endl;
        fout_delta_nll << mu_step << " " << delta_likelihood << endl;
    }


    // 4d
    cout << "\n4 d)" << endl;
    cout<< "Minimum of log likelihood diff = " << min_log_likelihood_diff << " for mu = " << min_log_likelihood_mu << endl;
    cout << "This means that " << min_log_likelihood_mu << " is the best estimator" << endl;
    cout << "The uncertainty on that is [" << estimator_uncertainty[0] << ", " << estimator_uncertainty.back() << "]" << endl;
        cout << "The uncertainty on that is " << estimator_uncertainty.back() - estimator_uncertainty[0] << endl;

    // Comparison to uncertainty on sample mean
    double delta_sample_mean = 1.65365/sqrt(234);
    cout << "For comparison: The uncertainty on the sample mean is " << delta_sample_mean << endl;
    fout.close();
    fout_nll.close();
    fout_delta_nll.close();

    // 4e: Compute likelihood ratio
    cout << "\n";
    cout << "4 e) " << endl;
    double likelihood_ratio;
    likelihood_ratio = -2 * log(likelihood_mean_of_sample / prob_mu_is_k(daten));
    cout << "The likelihood ratio -2 ln(Lambda) = " << likelihood_ratio << endl;

    double n_dof = 233;
    double relative_deviation = (likelihood_ratio - n_dof)/(sqrt(2*n_dof));
    cout << "The relative deviation of my likelihood ratio from the mean: z = " << relative_deviation << endl;



}