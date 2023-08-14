//
// Created by andi on 11.08.23.
//

#ifndef LEARNINGPROJECT_SECONDCORRLENGHTHANDLER_H
#define LEARNINGPROJECT_SECONDCORRLENGHTHANDLER_H

#include "../../Header/Helpfunctions and Classes.h"
#include <map>
#include <boost/range/adaptor/reversed.hpp>
#include "calcHandler.h"

class secondCorrLenghtHandler : public calcHandler {
    using calcHandler::calcHandler;     // if i don't need a special constructor i can write this to inherit from base
    vector<double> mu_0_x{};              // vector to hold the m^2 for every thermal realization
    vector<double> mu_2_x{};
    vector<double> mu_0_y{};
    vector<double> mu_2_y{};
    double tau;
public:

    void pre_routine() override {
        calcFile.open(root / "second.corrlength");
        // Generating the vector of Ls
        init_file_header(calcFile);

    }

    void setting_pre_routine(fs::path setting_path) override {
        // we need a map that maps the system size to the vector containing all the m values
        // it has to be reset for every setting / temperature
        fs::path txt_path = findFirstTxtFile(setting_path);
        tau = extractTauValue(txt_path);
    }

    void realization_routine(vector<double> &lat_q, double temp, double t) override {
        // starting with m... we calculate m for every row, that will be mx
        // we iterate over rows(cols)
        print_vector(lat_q);
        exit(0);
        for(int i = 0; i < lat_dim; i++) {
            double mx = 0;
            double my = 0;
            double mu2x = 0;
            double mu2y = 0;
            // we iterate over every element in the row (col) to get mx (my)
            for(int j = 0; j < lat_dim; j++) {
                mx += lat_q[ind_to_1D(i, j, lat_dim)];
                my += lat_q[ind_to_1D(j, i, lat_dim)];      // col fixed
                // for the second moment we need another loop
                for(int k = 0; k < lat_dim; k++) {
                    // we just use every pair?
                    // j and k are the indices that run over the row (col)
                    mu2x += pow((j-k), 2) * lat_q[ind_to_1D(i, j, lat_dim)] * lat_q[ind_to_1D(i, k, lat_dim)];
                    // col stays the same (i)
                    mu2y += pow((j-k), 2) * lat_q[ind_to_1D(j, i, lat_dim)] * lat_q[ind_to_1D(k, i, lat_dim)];

                }
            }
            // normalizing
            mx /= (double)lat_dim;
            my /= (double)lat_dim;

            // we add m^2 to the mu vector
            mu_0_x.push_back(mx * mx);
            mu_0_y.push_back(my * my);
            mu_2_x.push_back(mu2x);
            mu_2_y.push_back(mu2y);

        }
    }

    void setting_post_routine() override {
        // Getting tau into the file, then xi_2nd etc
        calcFile << tau << ",";
        cout << mu_0_x.size() << endl;
        print_vector(mu_0_x);
        print_vector(mu_2_x);
        double chi_x_check = mean(mu_0_x) * (double) lat_dim;
        double chi_x = (double)lat_dim / (double) mu_0_x.size()  *  std::reduce(mu_0_x.begin(), mu_0_x.end(),
                                           0.0, // initial value for the reduction (sum)
                                           std::plus<double>());
        // TODO need to look for error, to tired for that shit
        cout << "chi_x = " << chi_x << " vs chi_x_check =" << chi_x_check << endl;
        double chi_x_err = std::transform_reduce(mu_0_x.begin(), mu_0_x.end(),
                                                 0.0, // initial value for the reduction (sum)
                                                 std::plus<double>(), // transformation (square)
                                                 [&chi_x](double b) { return pow((b * b - chi_x), 2)  ; });
        double chi_y = (double)lat_dim / (double) mu_0_x.size() *  std::reduce(mu_0_y.begin(), mu_0_y.end(),
                                                                0.0, // initial value for the reduction (sum)
                                                                std::plus<double>());
        double mu2_x_check = mean(mu_2_x) / (double) lat_dim;
        double mu2_x = std::reduce(mu_2_x.begin(), mu_2_x.end(),
                                   0.0, // initial value for the reduction (sum)
                                   std::plus<double>()) / (double) lat_dim / (double) mu_0_x.size();
        cout << "mu2_x = " << mu2_x << " vs mu2_x_check =" << mu2_x_check << endl;

        double mu2_y = std::reduce(mu_2_y.begin(), mu_2_y.end(),
                                   0.0, // initial value for the reduction (sum)
                                   std::plus<double>()) / (double) lat_dim / (double) mu_0_x.size();
        cout << endl << "chi_x = " << chi_x << endl;
        cout << endl << "mu_2_x = " << mu2_x << endl;
        double xix = sqrt(mu2_x / (4 * chi_x));
        double xiy = sqrt(mu2_y / (4 * chi_y));
        exit(0);


        calcFile << xix << "," << xiy << "," << endl;
    }

    void post_routine() override {
        cout << "Calling BinderHandler post routine" << endl;
    }

    vector<int> generate_L(int starting_k, int n, int nr_Ls) {
        vector<int> L = {(int)sqrt(n) / (2*starting_k)};
        for (int k = starting_k + 1; k < starting_k + nr_Ls; k++) {
            int L_val = (int)sqrt(n) / (2*k);
            cout << L.back() - L_val << endl;
            if (L.back() - L_val < 2) {
                L_val = L.back() - 2;
            }
            L.push_back(L_val);
        }
        return L;
    }

    void init_file_header(ofstream& file) {
        file << "tau," << "xix," << "xiy," << "xix_err," << "xiy_err," << endl;
    }
    
};


#endif //LEARNINGPROJECT_SECONDCORRLENGHTHANDLER_H
