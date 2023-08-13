//
// Created by andi on 06.04.23.
//

#ifndef LEARNINGPROJECT_HELPFUNCTIONS_AND_CLASSES_H
#define LEARNINGPROJECT_HELPFUNCTIONS_AND_CLASSES_H

//
// Created by andi on 28.03.23.
//


#include <iostream>
#include <complex>
#include <random>
#include <fstream>
#include <filesystem>
#include <utility>


// using namespaces!
using namespace std;
namespace fs = std::filesystem;
using namespace fs;
// new state is composed of every
typedef vector<double> entry_type;
typedef  vector<vector<entry_type>> state_type;

template <class map_like>
void write_parameters(ofstream& file, map_like paras) {
    for(const auto& pair : paras) {
        file << pair.first << "," << pair.second << endl;
    }
}

int pymod(int a, int b) {
    return ((b + (a % b)) % b);
}

template <class T>
void print_vector(vector<T> vec) {
    cout << vec.size() << endl;
    for(int i = 0; i < vec.size(); i++) {
        if(i % 20 == 0) {
            cout << "\n";
        }
        cout << vec[i] << ", ";
    }
    cout << endl;
}


void calc_corr(vector<vector<complex<double>>> &f, vector<double> &C_x, vector<double> &C_y) {
    // i guess we could also use just normal doubles, would also give some performance
    // but i actually really don't think that it will be very computationally difficult
    int lat_dim = f.size();
    // I don't think that we have to pay attention to the lattice spacings. If we set them one we are just measuring
    // the distance in multiples of the lattice spacing
    // because of PBC, the maximum distance is half the lat_dim
    int d_max = lat_dim / 2;

    // loop over the rows(cols)
    for(int i = 0; i < lat_dim; i++) {
        // for every row we loop over all possible distances, which are 0, 1, ..., d_max
        for(int d = 0; d <= d_max; d++) {
            // for every distance we loop over every lattice site, so over every col (row)
            for(int j = 0; j < lat_dim; j++) {
            // we only have to add up the +d terms since we have pbc?
            C_x[d] += f[i][j].real() * f[i][(j + d) % lat_dim].real();
            C_y[d] += f[j][i].real() * f[(j + d) % lat_dim][i].real();
            }
        }
    }
    // meaning now for C_x[d] we have had lat_dim(i) * lat_dim(j) additions, so normalization
    for(int d = 0; d <= d_max; d++) {
        C_x[d] /= (lat_dim * lat_dim);
        C_y[d] /= (lat_dim * lat_dim);
    }
}

template <typename T, size_t N>
void print_array(const T (&arr)[N]) {
    for (size_t i = 0; i < N; ++i) {
        std::cout << arr[i] << " ";
    }
    std::cout << std::endl;
}

template <typename T>
void print_array(const T (&arr), size_t L) {
    for (size_t i = 0; i < L; ++i) {
        std::cout << arr[i] << " ";
    }
    std::cout << std::endl;
}

template <class value_type, class result_type>
void calc_corr(vector<vector<value_type>> &f, result_type &C_x, result_type &C_y) {
    // i guess we could also use just normal doubles, would also give some performance
    // but i actually really don't think that it will be very computationally difficult
    int lat_dim = f.size();
    // I don't think that we have to pay attention to the lattice spacings. If we set them one we are just measuring
    // the distance in multiples of the lattice spacing
    // because of PBC, the maximum distance is half the lat_dim
    int d_max = lat_dim / 2;

    // loop over the rows(cols)
    for(int i = 0; i < lat_dim; i++) {
        // for every row we loop over all possible distances, which are 0, 1, ..., d_max
        for(int d = 0; d <= d_max; d++) {
            // for every distance we loop over every lattice site, so over every col (row)
            for(int j = 0; j < lat_dim; j++) {
                // we only have to add up the +d terms since we have pbc?
                C_x[d] += f[i][j] * f[i][(j + d) % lat_dim];
                C_y[d] += f[j][i] * f[(j + d) % lat_dim][i];
            }
        }
    }
    // meaning now for C_x[d] we have had lat_dim(i) * lat_dim(j) additions, so normalization
    for(int d = 0; d <= d_max; d++) {
        C_x[d] /= (lat_dim * lat_dim);
        C_y[d] /= (lat_dim * lat_dim);
    }
}


vector<fs::path> list_dir_paths(const fs::path& root)
{
    vector<fs::path> dir_paths;

    if (fs::is_directory(root))
    {
        for (const auto& entry : fs::directory_iterator(root))
        {
            // if regular file, do nothing?
            // if directory, add
            if (fs::is_directory(entry.path()) && entry.path().filename() != "plots")
            {
                // check for the fucking .ipynb folders
                if(entry.path().filename().string()[0] != '.') {
                    dir_paths.push_back(entry.path());
                    vector<fs::path> sub_dir_paths = list_dir_paths(entry.path());
                    dir_paths.insert(dir_paths.end(), sub_dir_paths.begin(), sub_dir_paths.end());
                }
            }
        }
    }
    return dir_paths;
}


vector<fs::path> list_csv_files(const fs::path& root)
{
    vector<fs::path> csv_files;

    if (fs::is_directory(root))
    {
        for (const auto& entry : fs::directory_iterator(root))
        {
            if (fs::is_regular_file(entry.path()) && entry.path().extension() == ".csv")
            {
                csv_files.push_back(entry.path());
            }
            else if (fs::is_directory(entry.path()))
            {
                vector<fs::path> sub_csv_files = list_csv_files(entry.path());
                csv_files.insert(csv_files.end(), sub_csv_files.begin(), sub_csv_files.end());
            }
        }
    }

    return csv_files;
}

string readLineAt(std::ifstream& file, int i) {
    std::string line;
    std::string running_line;
    int running = 0;
    // Init the line to be first line
    std::getline(file, running_line);
    while (std::getline(file, line) && (running != i))
    {
        running++;
        // overwrite if line number is found
        running_line = line;
    }
    return running_line;
}

int getNrRowsFile(std::ifstream& file) {
    std::string line;
    int running = 0;
    while (std::getline(file, line))
    {
        running++;
    }
    return running;
}

int getNrRows(fs::path csv_path) {
    ifstream file(csv_path);
    // check if file is opened
    if(file.is_open()) {
        // cout << "File successfully opened" << endl;
    } else {
        cout << "Failed to open file" << endl;
        // abort if file is not opened
        exit(0);
    }
    return getNrRowsFile(file);
}

string readLastLine(std::ifstream& file) {
    std::string line;
    std::string lastLine;
    while (std::getline(file, line))
    {
        lastLine = line;
    }
    return lastLine;
}


vector<complex<double>> readValuesAt(ifstream& file, int ind, double& T, double& t) {
    string line = readLineAt(file, ind);
    vector<complex<double>> values;

    stringstream llss(line);
    string token;
    while (std::getline(llss, token, ',')) {
        if(token != "t") {
            complex<double> value = stod(token);
            values.push_back(value);
        }
    }
    // we delete the last value as it is the temperature
    // or do we need the temperature?
    // we might need it but we cannot really give it back
    t = values[0].real();
    values.erase(values.begin());
    T = values[0].real();
    values.erase(values.begin());

    return values;
}

vector<double> readDoubleValuesAt(ifstream& file, int ind, double& T, double& t) {
    string line = readLineAt(file, ind);
    vector<double> values;

    stringstream llss(line);
    string token;
    while (std::getline(llss, token, ',')) {
        if(token != "t") {
            double value = stod(token);

            values.push_back(value);
        }
    }
    // we delete the last value as it is the temperature
    // or do we need the temperature?
    // we might need it but we cannot really give it back
    t = values[0];
    values.erase(values.begin());
    T = values[0];
    values.erase(values.begin());

    return values;
}

vector<complex<double>> readLastValues(ifstream& file, double& T) {
    double t;
    int ind = -1;
    return readValuesAt(file, ind, T, t);
}

template <class value_type, template<class, class> class container>
void chess_trafo(container<value_type, std::allocator<value_type>>& vec) {
    int lat_dim = sqrt(vec.size());
    for (int i = 0; i < lat_dim/2; i++) {
        for (int j = 0; j < lat_dim/2; j++) {
            vec[2*i * lat_dim + 2 * j] *= (-1);
            vec[(2*i+1) * lat_dim + 2 * j + 1] *= (-1);
        }
    }
}

template <class container>
void chess_trafo(container& vec, int lat_dim) {
    for (int i = 0; i < lat_dim/2; i++) {
        for (int j = 0; j < lat_dim/2; j++) {
            vec[2*i * lat_dim + 2 * j] *= (-1);
            vec[(2*i+1) * lat_dim + 2 * j + 1] *= (-1);
        }
    }
}

double stringToDouble(const std::string& str) {
    try {
        return std::stod(str);
    } catch (const std::invalid_argument&) {
        return std::nan(""); // Return NaN if conversion fails
    }
}

template <typename T>
std::string findClosestStem(const std::vector<fs::path>& folders, const T& targetValue) {
    std::string closestFolder = "None";
    double minDifference = std::numeric_limits<double>::max();

    for (const auto& folder : folders) {
        double folderValue = stringToDouble(folder.stem().string());
        double difference = std::abs(folderValue - static_cast<double>(targetValue));
        if (difference < minDifference && difference > 0.0001) {
            minDifference = difference;
            closestFolder = folder.string();
        }
    }

    return closestFolder;
}

template <typename T>
std::string findClosestDir(const std::vector<fs::path>& folders, const T& targetValue, bool lowerTemp = true) {
    std::string closestFolder = "None";
    double minDifference = std::numeric_limits<double>::max();

    for (const auto& folder : folders) {
        double folderValue = stringToDouble(folder.filename().string());
        double difference = std::abs(folderValue - static_cast<double>(targetValue));
        if (difference < minDifference && difference > 0.0001) {
            if(lowerTemp) {
                // if lowerTemp is true, we only accept that this is the closest folder of the folder value is
                // smaller than the target value, this is important for example for binder cum calculations
                if(folderValue < targetValue) {
                    minDifference = difference;
                    closestFolder = folder.string();
                }
            } else {
                // else we set it anyway
                minDifference = difference;
                closestFolder = folder.string();
            }
        }
    }

    return closestFolder;
}


template <class value_type>
vector<vector<value_type>> oneD_to_twoD(vector<value_type> &q) {
    int N = q.size();
    int lat_dim = (int)sqrt(q.size());


    vector<vector<value_type>> q_2D = vector<vector<value_type>>( lat_dim, vector<value_type>(lat_dim, 0));
    // now i fill it with values
    // actually there was probably a function for this in some package...
    int n;
    int m;

    for(int i = 0; i < N; i++) {
        n = (i % lat_dim);
        m = i / lat_dim;
        // Zeile kreuz spalte?
        q_2D[m][n] = q[i];
    }

    return q_2D;
}

template <template<class> class container, class value_type>
void printComplexMatrix(container<value_type>* arr) {
    int size = sizeof(arr) / (sizeof(value_type) * 2);
    int dimension = static_cast<int>(std::sqrt(size));
    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            // Access the real and imaginary parts of the complex number
            double realPart = arr[i * dimension + j][0];
            double imagPart = arr[i * dimension + j][1];

            // Print the complex number as "real + imag i"
            std::cout << realPart << " + " << imagPart << "i\t";
        }
        std::cout << std::endl;
    }
}

template <class value_type>
void printComplexMatrix(value_type* arr, int lat_dim) {
    int dimension = lat_dim;
    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            // Access the real and imaginary parts of the complex number
            double realPart = arr[i * dimension + j][0];
            double imagPart = arr[i * dimension + j][1];

            // Print the complex number as "real + imag i"
            std::cout << realPart << " + " << imagPart << "i\t";
        }
        std::cout << std::endl;
    }
}
template <class value_type, int size>
void printComplexMatrix(value_type(&arr)[size]) {
    cout << size << endl;
    int dimension = static_cast<int>(std::sqrt(size));
    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            // Access the real and imaginary parts of the complex number
            double realPart = arr[i * dimension + j][0];
            double imagPart = arr[i * dimension + j][1];

            // Print the complex number as "real + imag i"
            std::cout << realPart << " + " << imagPart << "i\t";
        }
        std::cout << std::endl;
    }
}

string print_statetype(const state_type x) {
    string str = "";
    int n = x.size();
    str += "x: \n";
    str += "[";
    for (int i = 0; i < n; i++) {
        str += "[";
        for (int j = 0; j < n; j++) {
            if(j == n-1) {
                str += to_string(x[i][j][0]);
            } else {
                str += to_string(x[i][j][0]) + ", ";
            }
        }
        if(i < n-1) {
            str += "]\n ";
        }
    }
    str += "]]\n";
    str += "p: \n";
    str += "[";
    for (int i = 0; i < n; i++) {
        str += "[";
        for (int j = 0; j < n; j++) {
            if(j == n-1) {
                str += to_string(x[i][j][1]);
            } else {
                str += to_string(x[i][j][1]) + ", ";
            }
        }
        if(i < n-1) {
            str += "]\n ";
        }
    }
    str += "]]";
    return str;
}

template <typename T>
void print2DVector(const std::vector<std::vector<T>>& vec) {
    for (const auto& row : vec) {
        for (const auto& element : row) {
            std::cout << element << " ";
        }
        std::cout << std::endl;
    }
}
template <typename Container>
void print2DContainer(const Container& container) {

    for (const auto& row : container) {
        for (const auto& element : row) {
            std::cout << element << " ";
        }
        std::cout << std::endl;
    }
}

template <typename Container>
void print_container(const Container& container) {
    int i = 0;
    for (const auto& element : container) {
        if (i % 30 == 0) {
            cout << endl;
        }
        std::cout << element << ", ";
        i++;

    }
    cout << endl;
}


ifstream safe_read(string readpath, bool verbose=true) {
    if(verbose) {
        cout << "reading: " << readpath << endl;
    }
    ifstream file(readpath);
    // check if file is opened
    if(file.is_open()) {
        if(verbose) {
            cout << "File successfully opened" << endl;
        }
    } else {
        cout << "Failed to open file" << endl;
        // abort if file is not opened
        exit(0);
    }
    return file;
}


string trunc_double(double a, int precision=2) {
    stringstream stream;
    stream << std::fixed << std::setprecision(precision) << a;
    return stream.str();
}


void create_dir(const string dir_name) {
    // check wheter the directory already exists, if not create it
    if(!filesystem::is_directory(dir_name) || !filesystem::exists(dir_name)) {
        filesystem::create_directories(dir_name);
    }
}

string create_tree_name(double eta, double T, double dt, int n, double alpha, double beta, double J, double tau,
                        const string root) {
    string dir_name = root + "eta=" + trunc_double(eta)
                      + "/T=" + trunc_double(T) + "/dt=" +
                      trunc_double(dt, 4) + "/n="+ to_string(n) + "/alpha=" + trunc_double(alpha) + "/beta=" +
                      trunc_double(beta) + "/J=" + trunc_double(J) + "/tau=" + trunc_double(tau);
    create_dir(dir_name);
    return dir_name;
}

void write_parameters(ofstream& file, double eta, double T, double dt, int n, double alpha, double beta, double J,
                      double tau) {
    // insert the parameters
    file << "eta," << eta << ", \n";
    file << "T," << T << ", \n";
    file << "dt," << dt << ", \n";
    file << "n," << n << ", \n";
    file << "alpha," << alpha << ", \n";
    file << "beta," << beta << ", \n";
    file << "J," << J << ", \n";
    file << "tau," << tau << ", \n";
}


int count_files(string dir_name) {
    int i = 0;
    for(const auto & entry : filesystem::directory_iterator(dir_name)) {
        ++i;
    }
    return i;
}

/**
 * lattice_system is the CLASS (not an object) to create instances of
 * Okay i think what i wanted to do is impossible
 * Better Verision of the search grid function
 * @param paras vector of vectors of parameters
 */
template<typename T>
void search_grid_v2(vector<vector<double>> paras) {

    for(vector<double> para_set : paras) {

    }
}


double ind_value(vector<double> paras, int ind) {
    return paras[ind % paras.size()];
}


int ind_value(vector<int> paras, int ind) {
    return paras[ind % paras.size()];
}

template <class T>
vector<T>& operator+ (vector<T> &a, vector<T> &b) {
    if (a.size() != b.size()) {
        throw std::invalid_argument("Vectors must be of equal length.");
    }
    for(int i = 0; i < a.size(); i++) {
        a[i] += b[i];
    }
    return a;
}

template <class value_type>
void extract_cell(vector<value_type>& lattice, vector<value_type>& cell, int cell_nr, int L) {
    // okay so the vector is 1D
    // numbering the cells from left to right, top to bottom, we should be able to extract the
    // indices we need just from the cell nr
    //                   n
    //       ________________________
    //      |            |           |
    //      |     1      |     2     |
    //      |            |           |
    //      |____________|___________|
    //      |            |           |
    //   L  |     3      |      4    |
    //      |            |           |
    //      |____________|___________|
    //             L
    // first cell has the indices 1 : L, n : n + L, 2n : 2n +L, ...., n*(L-1) : n * L
    // second cell has the indices L : 2L, n + L: n + 2L, 2n +L : 2n +2L, ...., n*(L1) : n * L + L
    // third cell has the indices n * (row_nr * L)
    // kth cell has the indices kL : (k +1)L, n + kL
    int n = sqrt(lattice.size());
    int cells_per_row = n / L;
    int cell_in_row = cell_nr % cells_per_row;
    // determine in which rom number we are
    int row = cell_nr / cells_per_row;
    for(int j = 0; j < L; j++) {
        // extract values of j-th row
        for(int i = 0; i < L; i++) {
            int ind = n * (row * L + j) + i + cell_in_row * L;
            cell[j * L + i] = lattice[ind % lattice.size()];
        }
    }
}

fs::path findCsvFile(const fs::path& directory) {
    for (const auto& entry : fs::recursive_directory_iterator(directory)) {
        if (entry.is_regular_file() && entry.path().extension() == ".csv") {
            return entry.path();
        }
    }
    return ""; // Return an empty path if no .csv file is found
}

size_t get_sim_size(const fs::path& root) {
    // TODO its a bit lost that i have to read in a csv file to get the size but it takes like 0.1 seconds
    // so i don't care
    string csv_path = findCsvFile(root);
    ifstream csv_file = safe_read(csv_path);
    double dummy_T, dummy_t;
    vector<double> lattice = readDoubleValuesAt(csv_file, -1,  dummy_T, dummy_t);
    return (size_t) sqrt(lattice.size());
}

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

int ind_to_1D(int row, int col, int lat_dim) {
    return row * lat_dim + col;
}

fs::path findFirstTxtFile(const fs::path& directory) {
    for (const auto& entry : fs::directory_iterator(directory)) {
        if (entry.is_regular_file() && entry.path().extension() == ".txt") {
            return entry.path();
        }
    }

    return {}; // Return an empty path if no .txt file is found
}

double extractTauValue(const fs::path& filePath) {
    ifstream inputFile = safe_read(filePath);

    string line;
    while (std::getline(inputFile, line)) {
        size_t pos = line.find(",");
        if (pos != std::string::npos) {
            std::string parameter = line.substr(0, pos);
            if (parameter == "tau") {
                cout << line.substr(pos + 1);
                double tauValue = std::stod(line.substr(pos + 1));
                return tauValue;
            }
        }
    }

    return 0.0; // Return a default value if 'tau' is not found
}

#endif //LEARNINGPROJECT_HELPFUNCTIONS_AND_CLASSES_H
