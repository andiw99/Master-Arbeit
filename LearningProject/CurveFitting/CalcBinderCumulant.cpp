//
// Created by andi on 27.07.23.
//

#include "CalcBinderCumulant.h"

template <class value_type>
value_type mean(vector<value_type> vec) {
    // get lattice dimension for looping and averaging
    auto n = vec.size();
    double m = 0;
    // add everything up
    // for the cumulant it is not important which q value is next to each other
    // so we just leave the vector 1D and only sum over one index
    for(int i=0; i < n; i++) {
        m += vec[i];
    }

    // we have added everything up for one lattice, but we still need to average
    m /= n;
    // this function reduced to returning the mean value of a vector
    return m;
}

template <class value_type>
value_type isotropic_rms(vector<value_type> vec) {
    // get lattice dimension for looping and averaging
    auto n = vec.size();
    double m2 = 0;
    // okay this time we have somehow correlations between the lattice sites,
    // we add up q[i] * q[j]
    for(int i=0; i < n; i++) {
        for (int j=1; j < n; j++) {
            m2 += (vec[i] * vec[j]);
        }
    }

    // we have added everything up for one lattice, but we still need to average
    // this function reduced to returning the mean value of a vector
    return m2;
}


int main(int argc, char* argv[]) {
    // We want to calculate the binder cumulant which is the ratio of <m^4> to <m^2>^2
    // for us <m^k> should just be the sum over all (lattice sites) ^k normalized

    // We just read from the parameter file and compile this when we start the run, but we run when we have finished
    // the run
    fs::path root;
    if(argc >= 2) {
        // if we give some argument, doesnt even matter what argument, we take the parameter file values
        // It is called in the same directory as the run itself, so why would you need the ../ in front
        root = adaptive_tempscan_root;
    } else {
        root = "../../Generated content/Coulomb/Binder2";
    }
    cout << "root in CalcBinderCumulant" << endl;
    cout << root << endl;

    vector<fs::path> temp_directories = list_dir_paths(root);
    print_vector(temp_directories);

    ofstream cumsList;
    cumsList.open(root/"binder.cumulants");
    cumsList << "T,U,m" << endl;
    for(auto path : temp_directories) {
        // We need to read in every csv file and then the last file to calculate Binder and at last we have to average
        // We might should program this with error calculation in the hinterkopf
        vector<fs::path> csv_files = list_csv_files(path);
        // If i want to calculate the error i certainly need the number of realizations i have
        const int nr_csv_files = csv_files.size();
        // we now know the number of csv files which means we know how many spaces our binder
        // comulant vector needs
        // WAIT: We dont need a vector for the cumulant, as the cumulant is the ratio
        // of the expecation/mean values. What we need are vectors for the values of m^2 and m^4.
        // Since we get one value for the magnetization for every realization
        // WAIT: You only need the magnetization m for every lattice you idiot, everything else
        // happens when averaging over the different runs
        vector<double> m = vector<double>(nr_csv_files, 0);
        // checking that binder vec only consits of zeros
        // now we cycle over every csv file, calc <m^4> / <m^2> and fill one space of the binder_vec
        // when cycling we can already add the value of this run to the running value for the average
        double m4 = 0;
        double m2 = 0;
        double m1 = 0;
        double m2_rms = 0;
        double T, t;
        size_t n = 1;
        for(int i = 0; i < nr_csv_files; i++) {
            ifstream file = safe_read(csv_files[i], true);

            auto lat_q = readDoubleValuesAt(file, -1,  T, t);
            if (i == 0) {
                n = lat_q.size();
            }
            // now calcing the magnetization for this run/csv_file
            m[i] = mean(lat_q);
            m4 += pow(m[i], 4);
            m2 += pow(m[i], 2);
            m1 += m[i];
            m2_rms += isotropic_rms(lat_q);
        }
        print_vector(m);
        // okay we have every m, so now we can calc the mean of m^4 and the mean of m^2 squared
        // average over the number of runs
        m4 /= nr_csv_files;
        m2 /= nr_csv_files;
        m1 /= nr_csv_files;
        // ''thermal average''
        m2_rms /= nr_csv_files;
        // root and average over number of lattice sites
        double m_rms = sqrt(m2_rms) / (double)n;

        double cum = m4 / (m2 * m2);

        cout << "Temperature " << T << endl;
        cout << "cumulant value = " << cum << endl;
        cout << "magnetization value = " << m1 << endl;
        cout << endl << endl;

        ofstream ofile;
        ofile.open(path/ "binder.cumulant");
        ofile << cum;
        cumsList << T << "," << cum << "," << m1 << endl;
    }

    return 0;
}