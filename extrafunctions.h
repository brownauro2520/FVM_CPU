#ifndef EXTRAFUNCTIONS_H_INCLUDED
#define EXTRAFUNCTIONS_H_INCLUDED

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <limits>
#include <fstream>
#include <sstream>
using namespace std;

vector<double> divideVS(const vector<double> &vec, double scalar) {
  vector<double> result(vec.size());
  for (int i = 0; i < vec.size(); i++) {
    result[i] = vec[i] / scalar;
  }
  return result;
}

double sum_V(vector<vector<double>> v) {
  double sum = 0;
  for (int i = 0; i < v.size(); i++) {
    for (int j = 0; j < v[i].size(); j++) {
      sum += abs(v[i][j]);
    }
  }
  return sum;
}

void savefilev(vector<double> vec, string file_name) {
    ofstream outfile(file_name);
    int max = std::numeric_limits<double>::digits10;
    for (int i = 0; i < vec.size(); i++) {
        outfile << setprecision(max) << vec[i] << endl;
    }
    outfile.close();
    cout << "Vector saved to file: " << file_name << endl;
}

void savefilem(vector<vector<double>>& matrix, string fileName) {
    ofstream file;
    file.open(fileName);
    int max = std::numeric_limits<double>::digits10;
    for (int i = 0; i < matrix.size(); i++) {
        for (int j = 0; j < matrix[i].size(); j++) {
            file << setprecision(max) <<  matrix[i][j] << " ";
        }
        file << endl;
    }
    file.close();
    cout << "Matrix saved to file: " << fileName << endl;
}

vector<double> add_v(vector<double> vec1, vector<double> vec2) {
    int size = vec1.size();
    vector<double> result(size);
    for (int i = 0; i < size; i++) {
        result[i] = vec1[i] + vec2[i];
    }
    return result;
}

double dot_product(const vector<double> &v1, const vector<double> &v2) {
    double result = 0;
    if (v1.size() != v2.size()) {
        cout << "Invalid size" << endl;
        return 0;
    }
    for (int i = 0; i < v1.size(); i++) {
        result += v1[i] * v2[i];
    }
    return result;
}

std::vector<double> Scalar(const std::vector<double>& vector, double scalar) {
    std::vector<double> result(vector.size());
    for (int i = 0; i < vector.size(); i++) {
        result[i] = vector[i] * scalar;
    }
    return result;
}

vector<double> subtractVectors(const vector<double>& vector1, const vector<double>& vector2) {
    vector<double> result(vector1.size());
    for (int i = 0; i < vector1.size(); i++) {
        result[i] = vector1[i] - vector2[i];
    }
    return result;
}

vector<double> addVectors(const vector<double>& a, const vector<double>& b) {
    vector<double> result;
    int size = min(a.size(), b.size());
    result.reserve(size);
    for (int i = 0; i < size; i++) {
        result.push_back(a[i] + b[i]);
    }
    return result;
}

vector<vector<double>> readmatrix(string fileName) {
    vector<vector<double>> data;
    ifstream file(fileName);
    if (file.is_open()) {
        string line;
        while (getline(file, line)) {
            vector<double> row;
            stringstream ss(line);
            double value;
            while (ss >> value) {
                row.push_back(value);
            }
            data.push_back(row);
        }
        file.close();
    } else {
        cout << "Unable to open file." << fileName << endl;
    }
    return data;
}

vector<double> readvector(const string& fileName) {
    vector<double> values;
    ifstream file(fileName);
    if (file.is_open()) {
        double value;
        while (file >> value) {
            values.push_back(value);
        }
        file.close();
    } else {
        cerr << "Failed to open file " << fileName << endl;
    }
    return values;
}

#endif // EXTRAFUNCTIONS_H_INCLUDED