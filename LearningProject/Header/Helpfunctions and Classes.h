//
// Created by andi on 06.04.23.
//

#ifndef LEARNINGPROJECT_HELPFUNCTIONS_AND_CLASSES_H
#define LEARNINGPROJECT_HELPFUNCTIONS_AND_CLASSES_H
#define EIGEN_NO_CUDA

//
// Created by andi on 28.03.23.
//


#include <iostream>
#include <complex>
#include <random>
#include <fstream>
#include <filesystem>
#include <utility>
#include <map>
#include <array>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Core>
#include <unsupported/Eigen/NonLinearOptimization>

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

void print_vector_of_vectors(const std::vector<std::vector<double>>& vec) {
    // Loop through each row
    for (const auto& row : vec) {
        // Loop through each element in the row
        for (const auto& element : row) {
            std::cout << element << " ";
        }
        std::cout << std::endl; // Move to the next line after printing each row
    }
}

int pymod(int a, int b) {
    return ((b + (a % b)) % b);
}


template <class T>
void print_vector(vector<T> vec, int breakpoint = 20) {
    cout << vec.size() << endl;
    for(int i = 0; i < vec.size(); i++) {
        if(i % breakpoint == 0) {
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
        std::cout << arr[i] << ", ";
    }
    std::cout << std::endl;
}

template <typename T>
void print_complex_array(const T (&arr), size_t L) {
    for (size_t i = 0; i < L; ++i) {
        std::cout << arr[i][0] << " + " << arr[i][1] << ", ";
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


string get_current_time() {
    auto timepoint = chrono::system_clock::now().time_since_epoch();
    auto hour = chrono::duration_cast<chrono::hours>(timepoint) % 24;
    auto minute = chrono::duration_cast<chrono::minutes>(timepoint) % 60;
    auto second = chrono::duration_cast<chrono::seconds>(timepoint) % 60;
    auto millisecond = chrono::duration_cast<chrono::milliseconds>(timepoint) % 1000;
    auto microsecond = chrono::duration_cast<chrono::microseconds>(timepoint) % 1000;

    string hour_str = to_string(hour.count());
    if (hour.count() < 10) {
        hour_str = "0" + hour_str;
    }
    string minute_str = to_string(minute.count());
    if (minute.count() < 10) {
        minute_str = "0" + minute_str;
    }
    string seconds_str = to_string(second.count());
    if (second.count() < 10) {
        seconds_str = "0" + seconds_str;
    }
    string milliseconds_str = to_string(millisecond.count());
    string microseconds_str = to_string(microsecond.count());
    return hour_str + ":" + minute_str + ":" + seconds_str + ":" + milliseconds_str + ":" + microseconds_str;
}


vector<fs::path> list_dir_paths(const fs::path& root, bool recursive=true)
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
                    if(recursive){
                        vector<fs::path> sub_dir_paths = list_dir_paths(entry.path());
                        dir_paths.insert(dir_paths.end(), sub_dir_paths.begin(), sub_dir_paths.end());
                    }
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
                if(entry.path().filename().string()[0] != '.') {
                    vector<fs::path> sub_csv_files = list_csv_files(entry.path());
                    csv_files.insert(csv_files.end(), sub_csv_files.begin(), sub_csv_files.end());
                }
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
    int nr_errors = 0;
    while (std::getline(llss, token, ',')) {
        if(token != "t") {
            double value;
            try {
                value = stod(token);
            } catch (std::invalid_argument iaex) {
                std::cout << "Caught an error!" << std::endl;
                value = 0;
                nr_errors++;
            }
            values.push_back(value);
        }
    }
    if(nr_errors > 2) {
        cout << "too many stod Errors!" << endl;
        exit(0);
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
/*
template <class value_type, template<class, class> class container>
void chess_trafo_rectangular(container<value_type, std::allocator<value_type>>& vec, size_t dim_size_x) {
    int dim_size_y = vec.size() / dim_size_x;
    cout << "Using this?" << endl;
    int nr_rows = (int)round(dim_size_y / 2.0);
    int nr_cols = (int)round(dim_size_x / 2.0);

    for (int i = 0; i < nr_rows; i++) {            // row
        for (int j = 0; j < nr_cols; j++) {        // col
            vec[2*i * dim_size_x + 2 * j] *= (-1);      // both indices even
            if((nr_cols % 2 == 1) & (j == (nr_cols - 1)) || ((nr_rows % 2 == 1)) & (i == (nr_rows - 1))) {
                continue;
            }
            vec[(2*i+1) * dim_size_x + 2 * j + 1] *= (-1); // both indices uneven
        }
    }
}*/

std::string exec(const char* cmd) {
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    return result;
}
int timeStringToMinutes(const std::string& timeStr) {
    std::vector<int> timeParts;
    std::stringstream ss(timeStr);
    std::string part;
    while (std::getline(ss, part, ':')) {
        timeParts.push_back(std::stoi(part));
    }

    int days = 0, hours = 0, minutes = 0, seconds = 0;
    switch (timeParts.size()) {
        case 4:  // Format is d:hh:mm:ss
            days = timeParts[0];
            hours = timeParts[1];
            minutes = timeParts[2];
            seconds = timeParts[3];
            break;
        case 3:  // Format is hh:mm:ss
            hours = timeParts[0];
            minutes = timeParts[1];
            seconds = timeParts[2];
            break;
        case 2:  // Format is mm:ss
            minutes = timeParts[0];
            seconds = timeParts[1];
            break;
        default:
            std::cout << "Invalid time format!" << std::endl;
            return -1;
    }

    int totalMinutes = days * 24 * 60 + hours * 60 + minutes + seconds / 60;
    return totalMinutes;
}

int get_remaining_minutes() {
    string job_id_str;
    if (const char* job_id = std::getenv("SLURM_JOB_ID")){
        job_id_str = string(job_id);
        string command = "squeue -h -j " + job_id_str  + " -O TimeLeft";
        string timeleft = exec(command.c_str());
        int minutes_left = timeStringToMinutes(timeleft);
        return minutes_left;
    } else {
        return 1000;
    }

}

template <class container>
void chess_trafo_rectangular(container& vec, size_t dim_size_x) {
    int dim_size_y = vec.size() / 2 / dim_size_x;       // TODO important, this only works for the simulation!
    int nr_rows = (int)round(dim_size_y / 2.0);
    int nr_cols = (int)round(dim_size_x / 2.0);

    for (int i = 0; i < nr_rows; i++) {            // row
        for (int j = 0; j < nr_cols; j++) {        // col
            vec[2*i * dim_size_x + 2 * j] *= (-1);      // both indices even
            if((dim_size_x % 2 == 1) & (j == (nr_cols - 1)) || ((dim_size_y % 2 == 1)) & (i == (nr_rows - 1))) {
                continue;
            }
            vec[(2*i+1) * dim_size_x + 2 * j + 1] *= (-1); // both indices uneven
        }
    }
}
template <class container>
void chess_trafo_rectangular(container& vec, size_t dim_size_x, size_t dim_size_y) {
    int nr_rows = (int)round(dim_size_y / 2.0);
    int nr_cols = (int)round(dim_size_x / 2.0);

    for (int i = 0; i < nr_rows; i++) {            // row
        for (int j = 0; j < nr_cols; j++) {        // col
            vec[2*i * dim_size_x + 2 * j] *= (-1);      // both indices even
            if((dim_size_x % 2 == 1) & (j == (nr_cols - 1)) || ((dim_size_y % 2 == 1)) & (i == (nr_rows - 1))) {
                continue;
            }
            vec[(2*i+1) * dim_size_x + 2 * j + 1] *= (-1); // both indices uneven
        }
    }
}

template <class container>
void chess_trafo_rectangular_ptr(container& vec, size_t dim_size_x, size_t dim_size_y) {
    int nr_rows = (int)round(dim_size_y / 2.0);
    int nr_cols = (int)round(dim_size_x / 2.0);

    for (int i = 0; i < nr_rows; i++) {            // row
        for (int j = 0; j < nr_cols; j++) {        // col
            (*(vec[2*i * dim_size_x + 2 * j])) *= (-1);      // both indices even
            if((dim_size_x % 2 == 1) & (j == (nr_cols - 1)) || ((dim_size_y % 2 == 1)) & (i == (nr_rows - 1))) {
                continue;
            }
            (*(vec[(2*i+1) * dim_size_x + 2 * j + 1])) *= (-1); // both indices uneven
        }
    }
}

template <class value_type>
void extract_cell_ptrs(vector<value_type>& lattice, vector<value_type*>& cell, int cell_nr, int Lx, int Ly, int dim_size_x) {
    int cells_per_row = dim_size_x / Lx;
    int cell_in_row = cell_nr % cells_per_row;
    // determine in which rom number we are
    int row = cell_nr / cells_per_row;
    for(int j = 0; j < Ly; j++) {
        // extract values of j-th row
        for(int i = 0; i < Lx; i++) {
            int ind = dim_size_x * (row * Lx + j) + i + cell_in_row * Lx;
            cell[j * Lx + i] = &(*(lattice.begin() + (ind % lattice.size())));
        }
    }
}

template <class container>
void chess_trafo_rectangular_subsystems(container& vec, size_t dim_size_x, size_t dim_size_y, size_t Lx) {
    int nr_subsystems = dim_size_x / Lx;
    for(int i = 0; i < nr_subsystems; i++) {
        vector<double*> cell_ptrs(Lx * dim_size_y);
        extract_cell_ptrs(vec, cell_ptrs, i, Lx, dim_size_y, dim_size_x);
        chess_trafo_rectangular_ptr(cell_ptrs, Lx, dim_size_y);
    }
}


template <class value_type, template<class, class> class container, class Functor>
void trafo_rectangular(container<value_type, std::allocator<value_type>>& vec, Functor functor) {
    // useless, just use transform
    int size = vec.size();
    for (int i = 0; i < size; i++) {            // row
            vec[i] *= functor(vec[i]);
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
        if (difference < minDifference) {
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
template <class value_type1, class value_type2>
void print_pair_vector(vector<pair<value_type1, value_type2>>  &vec){
    cout << "vector size = " << vec.size();
    for(int i = 0; i < vec.size(); i++) {
        if (i % 10 == 0) {
            cout << endl;
        }
        cout << "(" << vec[i].first << ", " << vec[i].second << ")   ";
    }
    cout << endl;
}


vector<double> get_frequencies(int nr_times) {
    vector<double> freqs;
    for(int i = 0; i < nr_times; i++) {
        freqs.push_back(2 * M_PI * (double)(i - nr_times/2) / (double)nr_times);
    }
    return freqs;
}

vector<double> get_frequencies_fftw_order(int nr_times) {
    vector<double> freqs;
    int ind;
    int K = nr_times/ 2;
    for(int i = 0; i < nr_times; i++) {
        //  int ix = (j + Kx < Lx) ? j + Kx : j - Kx;
        ind = (i + K < nr_times) ? i + K : i - K;
        freqs.push_back(2 * M_PI * (double)(ind - nr_times/2) / (double)nr_times);
    }
    return freqs;
}

vector<double> p_to_vec(vector<vector<array<double, 2>>>& p) {
    vector<double> k;
    for(int j = 0; j<p.size(); j++) {
        k.push_back(p[j][j][0]);
    }
    return k;
}

template <class value_type>
void printComplexMatrix(value_type &arr, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            // Access the real and imaginary parts of the complex number
            double realPart = arr[i * cols + j][0];
            double imagPart = arr[i * cols + j][1];

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
void print_container(const Container& container, int total_nr, bool b=false) {
    int i = 0;
    for (const auto& element : container) {
        if(b) {
            if (i % 20 == 0) {
                cout << endl;
            }
        }
        std::cout << element << ", ";
        i++;
        if (i == total_nr) {
            break;
        }
    }
    cout << endl;
}

template <typename Container>
void print_container(const Container& container, bool b=false) {
    print_container(container, 0, b);
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

template <class value_type>
void extract_cell(vector<value_type>& lattice, vector<value_type>& cell, int cell_nr, int Lx, int Ly, int dim_size_x) {
    int cells_per_row = dim_size_x / Lx;
    int cell_in_row = cell_nr % cells_per_row;
    // determine in which rom number we are
    int row = cell_nr / cells_per_row;
    for(int j = 0; j < Ly; j++) {
        // extract values of j-th row
        for(int i = 0; i < Lx; i++) {
            int ind = dim_size_x * (row * Lx + j) + i + cell_in_row * Lx;
            cell[j * Lx + i] = lattice[ind % lattice.size()];
        }
    }
}



template <class container1, class container2>
void extract_cell(container1& lattice, container2& cell, int cell_nr, int Lx, int Ly, int dim_size_x) {
    int cells_per_row = dim_size_x / Lx;
    int cell_in_row = cell_nr % cells_per_row;
    // determine in which rom number we are
    int row = cell_nr / cells_per_row;
    for(int j = 0; j < Ly; j++) {
        // extract values of j-th row
        for(int i = 0; i < Lx; i++) {
            int ind = dim_size_x * (row * Lx + j) + i + cell_in_row * Lx;
            cell[j * Lx + i] = lattice[ind % lattice.size()];
        }
    }
}

bool isCSVFileInCurrentDirectory(path& currentDirectory) {
    // Iterate over the files in the current directory
    for (const auto& entry : std::filesystem::directory_iterator(currentDirectory)) {
        // Check if the file has a .csv extension
        if (entry.is_regular_file() && entry.path().extension() == ".csv") {
            return true;  // Found a .csv file
        }
    }

    return false;  // No .csv file found in the current directory
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
    for (const auto& entry : fs::recursive_directory_iterator(directory)) {
        if (entry.is_regular_file() && entry.path().extension() == ".txt") {
            return entry.path();
        }
    }

    return {}; // Return an empty path if no .txt file is found
}

string extractStringValueFromTxt(const fs::path& filePath, const string& value) {
    ifstream inputFile = safe_read(filePath);

    string line;
    while (std::getline(inputFile, line)) {
        size_t pos = line.find(",");
        if (pos != std::string::npos) {
            std::string parameter = line.substr(0, pos);
            if (parameter == value) {
                string value_str = line.substr(pos + 1);
                return value_str;
            }
        }
    }

    return "0.0"; // Return a default value if 'tau' is not found
}

double extractValueFromTxt(const fs::path& filePath, const string& value) {
    return stod(extractStringValueFromTxt(filePath, value));
}


double extractTauValue(const fs::path& filePath) {
    return (extractValueFromTxt(filePath, "tau"));
}

template <typename T>
void printAsMatrix(const std::vector<T>& vec, int rows, int cols) {
    if (vec.size() != rows * cols) {
        std::cout << "Error: Vector size does not match the specified matrix shape." << std::endl;
        return;
    }

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            std::cout << vec[i * cols + j] << "   "; // Assuming tab-separated values
        }
        std::cout << std::endl;
    }
}



int findHighestCSVNumber(const string& folderPath) {
    if (!std::filesystem::is_directory(folderPath)) {
        return -1;
    }
    int highestNumber = -1;

    for (const auto& entry : fs::directory_iterator(folderPath)) {
        if (entry.is_regular_file()) {
            const string& filename = entry.path().filename().string();

            // Check if the file has the .csv extension
            if (entry.path().extension() == ".csv") {
                // Extract the numeric part of the filename (excluding the extension)
                string numberPart = entry.path().stem().string();
                auto pos = numberPart.find('-');
                if(pos != string::npos) {
                    numberPart = numberPart.substr(0, pos);
                }
                try {
                    // Convert the numeric part to an integer
                    int number = stoi(numberPart);

                    // Update the highest number if the current number is higher
                    if (number > highestNumber) {
                        highestNumber = number;
                    }
                } catch (const invalid_argument&) {
                    // Ignore invalid filenames that cannot be converted to integers
                }
            }
        }
    }

    return highestNumber;
}

template <typename T>
std::vector<T> geomspace(T start, T stop, int num, bool endpoint = true) {
    std::vector<T> result;

    if (num < 1) {
        return result;
    }

    if (endpoint) {
        T factor = std::pow(stop / start, 1.0 / (num - 1));
        for (int i = 0; i < num; ++i) {
            T value = start * std::pow(factor, i);
            result.push_back(value);
        }
    } else {
        T factor = std::pow(stop / start, 1.0 / num);
        for (int i = 0; i < num; ++i) {
            T value = start * std::pow(factor, i);
            result.push_back(value);
        }
    }

    return result;
}


template<typename T>
std::vector<T> linspace(T start_in, T end_in, int num_in)
{

    std::vector<T> linspaced;

    double start = static_cast<double>(start_in);
    double end = static_cast<double>(end_in);
    double num = static_cast<double>(num_in);

    if (num == 0) { return linspaced; }
    if (num == 1)
    {
        linspaced.push_back(start);
        return linspaced;
    }

    double delta = (end - start) / (num - 1);

    for(int i=0; i < num-1; ++i)
    {
        linspaced.push_back((T)(start + delta * i));
    }
    linspaced.push_back(end); // I want to ensure that start and end
    // are exactly the same as the input
    return linspaced;
}

template<typename T, typename S>
std::vector<S> linspace(T start_in, T end_in, int num_in)
{

    std::vector<S> linspaced;

    double start = static_cast<double>(start_in);
    double end = static_cast<double>(end_in);
    double num = static_cast<double>(num_in);

    if (num == 0) { return linspaced; }
    if (num == 1)
    {
        linspaced.push_back((S)start);
        return linspaced;
    }

    double delta = (end - start) / (num - 1);

    for(int i=0; i < num-1; ++i)
    {
        linspaced.push_back((S)(start + delta * i));
    }
    linspaced.push_back((S)end); // I want to ensure that start and end
    // are exactly the same as the input
    return linspaced;
}

template<typename T>
std::vector<T> logspace(T start_in, T end_in, int num_in, T base_in = 2.0)
{
    std::vector<T> logspaced;

    double start = static_cast<double>(start_in);
    double end = static_cast<double>(end_in);
    double base = static_cast<double>(base_in);
    double num = static_cast<double>(num_in);

    if (num == 0) { return logspaced; }
    if (num == 1)
    {
        logspaced.push_back((T)pow(base, start));
        return logspaced;
    }

    double delta = (end - start) / (num - 1);

    for(int i=0; i < num-1; ++i)
    {
        logspaced.push_back((T)pow(base, start + delta * i));
    }
    logspaced.push_back((T)pow(base, end)); // I want to ensure that start and end
    // are exactly the same as the input
    return logspaced;
}

template<typename T, typename S>
std::vector<S> logspace(T start_in, T end_in, int num_in, T base_in = 2.0)
{
    std::vector<S> logspaced;

    double start = static_cast<double>(start_in);
    double end = static_cast<double>(end_in);
    double base = static_cast<double>(base_in);
    double num = static_cast<double>(num_in);

    if (num == 0) { return logspaced; }
    if (num == 1)
    {
        logspaced.push_back((S)pow(base, start));
        return logspaced;
    }

    double delta = (end - start) / (num - 1);

    for(int i=0; i < num-1; ++i)
    {
        logspaced.push_back((S)pow(base, start + delta * i));
    }
    logspaced.push_back((S)pow(base, end)); // I want to ensure that start and end
    // are exactly the same as the input
    return logspaced;
}

template <typename Key, typename Value>
void printMap(const std::map<Key, Value>& myMap) {
    for (const auto& pair : myMap) {
        std::cout << "Key: " << pair.first << ", Value: " << pair.second << std::endl;
    }
}


vector<vector<array<double, 2>>> init_q(size_t lat_dim) {
    // later we can either overload or adjust this function to work with lattice spacings
    vector<vector<array<double, 2>>> q = vector<vector<array<double, 2>>>(
            lat_dim, vector<array<double, 2>>(lat_dim, array<double, 2>()));
    // iterate over the lattice
    for(int m = 0; m < lat_dim; m++) {
        for(int n = 0; n < lat_dim; n++) {
            q[m][n][0]= 1.0 * n;
            q[m][n][1] = 1.0 * m;
        }
    }
    return q;
}
void fill_p(const vector<vector<array<double, 2>>> &q, vector<vector<array<double, 2>>> &p) {
    int lat_dim = q.size();         // f has size of lattice dim
    int N = lat_dim;
    // To center the impulses around the center, we shift by K = N/2
    int K = N / 2;
    // qi corresponds to xn, so i is the col
    // find out the lattice spacings a_x = x_1 - x_0
    double ax = q[0][1][0] - q[0][0][0];
    double ay = q[1][0][1] - q[0][0][1];
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            int i_ft = i < K ? i : i - N;
            int j_ft = j < K ? j : j - N;
            double p_i = 2 * M_PI * (double)i_ft / N / ax;
            double p_j = 2 * M_PI * (double)j_ft / N / ay;
            p[i][j] = array<double, 2>{p_i, p_j};
        }
    }
}

struct sin_functor {
    double p;
    sin_functor(double p): p(p) {}
    sin_functor(): p(1.0) {}
    template <class value_type>
    value_type operator()(value_type x) {
        return sin(p*x);
    }
};

struct cos_functor {
    template <class value_type>
    value_type operator()(value_type x) {
        return cos(x);
    }
};



template<typename scalar, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>
struct Functor
{
    typedef scalar Scalar;
    enum {
        InputsAtCompileTime = NX,
        ValuesAtCompileTime = NY
    };
    typedef Eigen::Matrix<Scalar,InputsAtCompileTime,1> InputType;
    typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,1> ValueType;
    typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;

    int m_inputs, m_values;

    Functor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
    Functor(int inputs, int values) : m_inputs(inputs), m_values(values) {}

    int inputs() const { return m_inputs; }
    int values() const { return m_values; }

};

struct LorentzianPeakFunctor : Functor<double> {
private:
    // Observations for a sample.
    Eigen::MatrixXd X_Y;
public:

    static Eigen::VectorXd get_starting_paras() {
        auto paras = Eigen::VectorXd(2);
        paras << 1.0, 1.0;
        return paras;
    }

    LorentzianPeakFunctor(const Eigen::MatrixXd &X_Y): Functor(X_Y.cols(), X_Y.rows()), X_Y(X_Y) {
    }

    int operator()(const Eigen::VectorXd &paras, Eigen::VectorXd &fvec) const {
        // fvec is the vector with residuals
        for(int i = 0; i < this->values(); i++) {
            fvec(i) = X_Y(i, 1) - paras(0) * 2 * paras(1) /
                                  (1 + X_Y(i, 0) * X_Y(i,0) * paras(1) * paras(1));
        }
        return 0;
    }
};

struct LorentzianPeakOffsetFunctor : Functor<double> {
private:
    // Observations for a sample.
    Eigen::MatrixXd X_Y;
public:

    static Eigen::VectorXd get_starting_paras() {
        auto paras = Eigen::VectorXd(3);
        paras << 1.0, 1.0, 0.0;
        return paras;
    }


    LorentzianPeakOffsetFunctor(const Eigen::MatrixXd &X_Y): Functor(X_Y.cols(), X_Y.rows()), X_Y(X_Y) {
    }

    int operator()(const Eigen::VectorXd &paras, Eigen::VectorXd &fvec) const {
        // fvec is the vector with residuals
        for(int i = 0; i < this->values(); i++) {
            fvec(i) = X_Y(i, 1) - paras(2) + paras(0) * 2 * paras(1) /
                                  (1 + X_Y(i, 0) * X_Y(i,0) * paras(1) * paras(1));
        }
        return 0;
    }
};

template <typename vec>
Eigen::MatrixXd construct_matrix(const vec &a, const vec &b) {
    Eigen::VectorXd a_eigen = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(a.data(), a.size());
    Eigen::VectorXd b_eigen = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(b.data(), b.size());

    Eigen::MatrixXd result_matrix(a.size(), 2);
    result_matrix << a_eigen, b_eigen;

    return result_matrix;
}

template <class value_type>
Eigen::MatrixXd construct_matrix(const value_type* a, const value_type *b, const int size) {
    Eigen::Map<Eigen::VectorXd> a_eigen(a, size);
    Eigen::Map<Eigen::VectorXd> b_eigen(b, size);

    Eigen::MatrixXd result_matrix(size, 2);
    result_matrix << a_eigen, b_eigen;

    return result_matrix;
}

template <class value_type>
Eigen::MatrixXd construct_matrix(vector<double> &a, value_type *b, const int size) {
    Eigen::VectorXd a_eigen = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(a.data(), (long)a.size());
    Eigen::Map<Eigen::VectorXd> b_eigen(b, size);

    Eigen::MatrixXd result_matrix(size, 2);
    result_matrix << a_eigen, b_eigen;

    return result_matrix;
}

Eigen::MatrixXd construct_matrix(vector<double> &a, double *b, const int size) {
    Eigen::VectorXd a_eigen = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(a.data(), (long)a.size());
    Eigen::Map<Eigen::VectorXd> b_eigen(b, size);

    Eigen::MatrixXd result_matrix(size, 2);
    result_matrix << a_eigen, b_eigen;
    return result_matrix;
}


Eigen::MatrixXd construct_matrix(vector<double> &a, vector<double> &b) {
    Eigen::VectorXd a_eigen = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(a.data(), (long)a.size());
    Eigen::VectorXd b_eigen = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(b.data(), (long)b.size());

    Eigen::MatrixXd result_matrix(a.size(), 2);
    result_matrix << a_eigen, b_eigen;

    return result_matrix;
}

void printMatrixXd(const Eigen::MatrixXd& matrix) {
    std::cout << "Matrix values:\n";
    for (int i = 0; i < matrix.rows(); ++i) {
        for (int j = 0; j < matrix.cols(); ++j) {
            std::cout << matrix(i, j) << " ";
        }
        std::cout << std::endl;
    }
}

template <class Functor>
Eigen::VectorXd fit_matrix(Eigen::MatrixXd X_Y_vals) {
    // printMatrixXd(X_Y_vals);

    Eigen::VectorXd params = Functor::get_starting_paras();
    // the params of the fit are the scaling and the correlation length? params(0) = amplitude params(1) = xi
    Functor functor(X_Y_vals);
    Eigen::NumericalDiff<Functor> numericalDiff(functor);
    Eigen::LevenbergMarquardt<Eigen::NumericalDiff<Functor>> lm(numericalDiff);

    Eigen::LevenbergMarquardtSpace::Status status = lm.minimize(params);
    // std::cout << "status: " << status << std::endl;
    return params;
}

Eigen::VectorXd fit_lorentz_peak(vector<double>& k_values, vector<double>& ft_values) {
    // constrtuct the matrix that holds the k and ft values
    Eigen::MatrixXd X_Y_vals(k_values.size(), 2);
    X_Y_vals = construct_matrix(k_values, ft_values);
    return fit_matrix<LorentzianPeakFunctor>(X_Y_vals);
}

Eigen::VectorXd fit_lorentz_peak(vector<double>& k_values, double* ft_values) {
    // constrtuct the matrix that holds the k and ft values
    int L = k_values.size();
    Eigen::MatrixXd X_Y_vals(k_values.size(), 2);
    X_Y_vals = construct_matrix(k_values, ft_values, L);
    return fit_matrix<LorentzianPeakFunctor>(X_Y_vals);
}

Eigen::VectorXd fit_offset_lorentz_peak(vector<double>& k_values, double* ft_values) {
    // constrtuct the matrix that holds the k and ft values
    int L = k_values.size();
    Eigen::MatrixXd X_Y_vals(k_values.size(), 2);
    X_Y_vals = construct_matrix(k_values, ft_values, L);
    return fit_matrix<LorentzianPeakOffsetFunctor>(X_Y_vals);
}

void readXiFromFile(const std::string& filename, std::vector<double>& xix_values, std::vector<double>& xiy_values,
                    std::vector<double>& times) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    std::string line;
    // first line is bullcr..?
    getline(file, line);
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string token;
        std::getline(iss, token, ','); // Read t, but don't store it
        double t = std::stod(token);
        std::getline(iss, token, ','); // Read xix
        double xix = std::stod(token);
        std::getline(iss, token, ','); // Read xiy
        double xiy = std::stod(token);

        times.push_back(t);
        xix_values.push_back(xix);
        xiy_values.push_back(xiy);
    }

    file.close();
}

void readCumFromFile(const std::string& filename, std::vector<double>& U_L_vec, std::vector<double>& times) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    std::string line;
    // the first line is the header
    getline(file, line);
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string token;
        std::getline(iss, token, ','); // Read t, but don't store it
        double t = std::stod(token);
        std::getline(iss, token, ','); // Read U_L
        double U_L = std::stod(token);

        times.push_back(t);
        U_L_vec.push_back(U_L);
    }

    file.close();
}

void readMagFromFile(const std::string& filename, std::vector<double>& m_vec, std::vector<double>& times) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    std::string line;
    // the first line is the header
    getline(file, line);
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string t, mVec;
        if (std::getline(iss, t, ';') && std::getline(iss, mVec, ';')) {
            std::istringstream mIss(mVec);
            std::string mValue;
            while (std::getline(mIss, mValue, ',')) {
                m_vec.push_back(std::stod(mValue));
            }
            times.push_back(stod(t));
        }
    }

    file.close();
}


struct absDiff1DStd
{
    double operator()(double x, double y) const
    {
        return fabs(x - y);
    }
};


double meanAbsDifference(double* f, int f_size) {
    // since what size is the summation on GPU faster?
    std::vector<double> diffs(f_size - 1);
    transform(f, f + f_size - 1, f + 1, diffs.begin(), absDiff1DStd());

/*
    for (int i = 0; i< 20; i++){
        cout << "|" <<f[i] << " - " << f[i+1]<< "|" << " = " <<  diffs[i] << endl;
    }
*/



    double accumulated_abs_diff = reduce(diffs.begin(), diffs.end(), 0.0) / (double)(f_size - 1);

    return accumulated_abs_diff;
}

double getMovingFactor(int nr_values, int min_ind, vector<double>& f, double avg_f) {
    int recent_ind = (int) (nr_values * 0.8);
    // the idea is to check if the differences between the values show a trend to be positive
    // or negative. But I think that comes done to calculating the deviation of the recent
    // average from the total average
    // small differences will
    // Idea would be to compare this deviation to the deviation per step
    // If we are in a high temperature state, the deviation per step will be relatively larger
    // If should not be dependent on the stepsize
    // first of all compute the avergage absolute difference per step
    int recent_size = nr_values - recent_ind;
    double f_start = reduce(f.begin() + min_ind, f.begin() + recent_ind) / (double)(recent_ind - min_ind);
    double f_end = avg_f;

    double* recent_f = &f[recent_ind];      // we need it still for the mean delta?
    double mean_abs_delta = meanAbsDifference(recent_f, recent_size);
    mean_abs_delta = max(1e-7, mean_abs_delta);
/*    cout << "meanAbsDifference: " << mean_abs_delta << endl;
    // and we need the recent mean
    cout << "val_end - val_start = " << fabs(f_end - f_start) << endl;*/
    double moving_factor = fabs(f_end - f_start) / ((double)recent_size * mean_abs_delta);

    return moving_factor;
}

double getMovingFactor(int nr_f_values, int min_ind, vector<double>& f) {
    double avg_f = accumulate(f.begin() + min_ind, f.end(), 0.0) / (double) (f.size() - min_ind);
    return getMovingFactor(nr_f_values, min_ind, f, avg_f);
}

std::vector<int> generateRandomIntegers(int n, int N) {
    std::vector<int> randomIntegers;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, N);

    for (int i = 0; i < n; ++i) {
        randomIntegers.push_back(dis(gen));
    }

    return randomIntegers;
}

double findFWHM(const double* y, const std::vector<double>& x, int size) {
    double halfMax = 0.5 * (*std::max_element(y, y + size)); // Half of the maximum y value
    double x1 = 0.0, x2 = 0.0; // Initialize variables for storing x values

    // Find the x values corresponding to the points where y(x) crosses the half-maximum line
    bool found1 = false, found2 = false;
    for (int i = 0; i < size - 1; ++i) {
        if (!found1 && y[i] >= halfMax) {
            x1 = x[i];
            found1 = true;
        }
        if (!found2 && y[i] < halfMax && y[i + 1] >= halfMax) {
            x2 = x[i];
            found2 = true;
        }
        if (found1 && found2) break; // Break loop if both points are found
    }

    // Calculate FWHM
    double fwhm = fabs(x2 - x1);
    return fwhm;
}

void omitHalfTimes(const string& fileName) {
    ifstream file(fileName);
    ofstream tempFile("temp.txt");

    if (!file.is_open() || !tempFile.is_open()) {
        cerr << "Error opening files." << endl;
        return;
    }

    string line;
    bool firstLine = true;
    int lineCount = 0;

    while (getline(file, line)) {
        if (firstLine) {
            tempFile << line << endl; // Write the first line as it is
            firstLine = false;
        } else {
            if (lineCount % 2 == 0) {
                tempFile << line << endl; // Write lines with even line count
            }
        }
        lineCount++;
    }

    file.close();
    tempFile.close();

    // Rename temp file to original file
    remove(fileName.c_str());
    rename("temp.txt", fileName.c_str());
}

#endif //LEARNINGPROJECT_HELPFUNCTIONS_AND_CLASSES_H
