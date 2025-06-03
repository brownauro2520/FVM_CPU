#ifndef READFILE_H_INCLUDED
#define READFILE_H_INCLUDED
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

using namespace std;

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

#endif // READFILE_H_INCLUDED