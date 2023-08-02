//
// Created by andi on 28.07.23.
//

#include "OriginalBinderCumulant.h"
//
// Created by andi on 27.07.23.
//

#include "CalcBinderCumulant.h"
#include <boost/range/adaptor/reversed.hpp>
#include "../Header/Helpfunctions and Classes.h"

template <class value_type>
value_type mean(vector<value_type> vec) {
    // get lattice dimension for looping and averaging
    auto n = vec.size();
    double m = 0;
    // add everything up
    // for the cumulant it is not important which q value is next to each other
    // so we just leave the vector 1D and only sum over one index
    for(int i=0; i < n; i++) {
/*        cout << vec[i] <<", ";
        if (i % 20 == 0) {
            cout << endl;
        }*/
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

template <class value_type>
value_type calc_s2(vector<value_type>& cell) {
    int L_sqaured = cell.size();
    double s2 = 0;
    // sum over rows inside the cell
    for(int j = 0; j < L_sqaured; j++) {
        // sum over first row
        for(int i = 0; i < L_sqaured; i++) {
            s2 += cell[i] * cell[j];
        }
    }
    return s2;
}

template <class value_type>
value_type calc_s4(vector<value_type>& cell) {
    // this is not working because of a 4-fach sum
    int L_squared = cell.size();
    double s4 = 0;
    // we have to sum over 4 different indices, we could probably reduce the
    // computational effort by only considering every pair S_i S_j S_k S_j once
    // would this actually be the right way?
    // Lets look what comes out first the double counting should actually just be a factor
    // that would not influence the intersection
    for(int j = 0; j < L_squared; j++) {
        for(int i = j; i < L_squared; i++) {
            for(int k = i; k < L_squared; k++) {
                for(int l = k; l < L_squared; l++) {
                    // if l becomes very large, like 100
                    // this wont be calculateable?
                    // 10000 ** 4
                    s4 += cell[i] * cell[j] * cell[k] * cell[l];
                }
            }
        }
    }
    return s4;
}



void init_s_map(map<int, vector<double>>& map, int n, int nr_Ls, int starting_k) {
    for (int k = starting_k; k < starting_k + nr_Ls; k++) {
        map[sqrt(n) / pow(2, k)] = vector<double>{};
    }
}

void init_file_header(ofstream& file,int n, int nr_Ls, int starting_k) {
    for (int k = starting_k; k < starting_k + nr_Ls; k++) {
        int L = sqrt(n) / pow(2, k);
        file << "," << L;
    }
    file << endl;
}

void init_file_header(ofstream& file, const vector<int>& L_vec) {
    for (int L : L_vec) {
        file << "," << L << "," << L << "_err";
    }
    file << endl;
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
        root = "../../Generated content/New/Coulomb/Critical Exponent";
    }
    cout << "root in CalcBinderCumulant" << endl;
    cout << root << endl;
    int nr_Ls = 11;
    int starting_k = 4;
    vector<fs::path> temp_directories = list_dir_paths(root);
    print_vector(temp_directories);
    cout << endl;

    ofstream cumsList;
    cumsList.open(root/"binder.cumulants");
    cumsList << "T";

    // I think i will use a map that maps the different L's to vectors with the s2, s4, values for
    // the different csv files

    map<int, vector<double>> S2_map = {};
    map<int, vector<double>> S4_map = {};
    int running = 0;
    size_t n = 1;
    vector<int> L_vec;
    for(auto path : temp_directories) {
        map<int, vector<double>> m_map = {};
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

        double T = 0, t;
        for(int i = 0; i < nr_csv_files; i++) {
            ifstream file = safe_read(csv_files[i], false);

            auto lat_q = readDoubleValuesAt(file, -1,  T, t);
            if (running == 0) {
                n = lat_q.size();
                L_vec = generate_L(starting_k, n, nr_Ls);
                running++;
                init_file_header(cumsList, L_vec);
            }
            // now calcing the magnetization for this run/csv_file
            // The original paper to the binder cumulant does something completely different
            // We divide our lattice into blocks of Dimension L and calculate a local magnetization
            // for every block. And the <m^2> is not <(1/ N sum s_i)^2>, but <(1/N² sum s_i s_j)> which
            // is completely different.
            // But I think the important part is that we oenly look at relatively small cells
            // we know n, how do we decide for the different Ls?
            // We want to have L = n/2^k to ensure I think that our blocks don't overlap
            // even though it probably won't make a difference anyway
            // largest L to consider might be n/2 because usually we do not have larger domains
            // At this point it would be nice if our lattice dim was also a power of 2

            // vectors to save the value for this csv file for s2, s4

            for(int L : L_vec) {
                double s2_L = 0;
                double s4_L = 0;
                // okay we now switch to calculating just m for every cell and then summing m_L ** 4
                // calculating the lattice dimension for this
                // the higher k is the more cells we "have to" deal with
                // I think we can calc the local magnetization for every cell and just add
                // them up and at last take a big average over all cells and all realizations
                // since we are just adding up the local magnetization between cells and runs this should
                // be fine
                // the number of cells for the power of k is
                int nr_cells = (int)(n / (L * L));
                // now looping over every cell
                for(int cell_nr = 0; cell_nr < nr_cells; cell_nr++) {
                    double m_L = 0;
                    // now we are in the ith cell and can start calculating the local magnetization powers
                    // it might be easiest and maybe also the most efficient to extract the spins of the cell
                    // out of the large lattice
                    vector<double> q_cell = vector<double>(L * L, 0);
                    extract_cell(lat_q, q_cell, cell_nr, L);
                    m_L = mean(q_cell);
                    m_map[L].push_back(m_L);
                }
                // so what I did now in this function is dividing the lattice and calculating the local
                // magnetization powers for every cell. We could now to the averaging over the lattice site number
/*                s2_L /= pow(L, 2*2);
                s4_L /= pow(L, 4*2);*/
                // we need to save this value for this csv file to do the thermal average later on.
                // this is the s2 value of a single csv file for a single L
/*                S2_map[L].push_back(s2_L);
                S4_map[L].push_back(s4_L);*/
            }
        }
        // okay so i now have added every csv file for one temperature into the m map
        // we now do the thermal average by calcing the powers of m, adding them and averaging
        cumsList << T;

        for(auto pair : boost::adaptors::reverse(m_map)){
            double m_L2 = std::transform_reduce(pair.second.begin(), pair.second.end(),
                                                0.0, // initial value for the reduction (sum)
                                                std::plus<double>(), // transformation (square)
                                                [](double b) { return b * b; });
            m_L2 /= pair.second.size();
            double m_L2_err = std::transform_reduce(pair.second.begin(), pair.second.end(),
                                                    0.0, // initial value for the reduction (sum)
                                                    std::plus<double>(), // transformation (square)
                                                    [&m_L2](double b) { return pow((b * b - m_L2), 2)  ; });
            m_L2_err /= (double)pow(pair.second.size(), 2);
            m_L2_err = sqrt(m_L2_err);
            double m_L4 = std::transform_reduce(pair.second.begin(), pair.second.end(),
                                                0.0, // initial value for the reduction (sum)
                                                std::plus<>(), // transformation (square)
                                                [](double b) { return (pow(b, 4)); });
            m_L4 /= pair.second.size();
            double m_L4_err = std::transform_reduce(pair.second.begin(), pair.second.end(),
                                                    0.0, // initial value for the reduction (sum)
                                                    std::plus<>(), // transformation (square)
                                                    [&m_L4](double b) { return pow(pow(b, 4) - m_L4, 2); });
            m_L4_err /= (double)pow(pair.second.size(), 2);
            m_L4_err = sqrt(m_L4_err);

            // we write it directly to a file?
            // we now have the binder cumulant for one temperature for different sizes.
            // nice would be a large file with the layout
            // T     L1        L2       ...
            // 0.1   1          1       ...
            // ...  ...        ...
            // 1     3          3       ...
            double cum = m_L4 / (m_L2 * m_L2);
            double cum_error = sqrt(pow(1 / m_L2 / m_L2 * m_L4_err, 2) + pow(2 * m_L4 / pow(m_L2, 3) * m_L2_err, 2));
            cumsList << "," << cum << "," << cum_error;
        }

        cumsList << endl;
        cout << "Temperature " << T << endl;
        cout << endl << endl;

        ofstream ofile;
        ofile.open(path/ "binder.cumulant");

        // okay we have every m, so now we can calc the mean of m^4 and the mean of m^2 squared
        // average over the number of runs
        // root and average over number of lattice sites


    }

    return 0;
}