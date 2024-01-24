#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <tuple>
#include <vector>

using namespace std;

int getSuperIndex(int row, int col, const int L) {
  if (row < 0 || row >= L) {
    return -1;
  }
  else if (col < 0 || col == (row*L - 1) || col == (row*L + L)) {
    return -1;
  }
  else {
    return row*L + col;
  }
}

double drawRandomNumber(random_device rd) {
  mt19937 gen(rd());
  uniform_real_distribution<> dis(0.0, 1.0);
  return dis(gen);
}

vector<vector<int>> init(vector<vector<int>> neighbourArray, const int L) {
  for (int row = 0; row < L; row++) {
    for (int col = 0; col < L; col++) {
        int index = getSuperIndex(row, col, L);
        vector<int> neighbours;
        neighbours.push_back(getSuperIndex(row - 1, col, L));
        neighbours.push_back(getSuperIndex(row + 1, col, L));
        neighbours.push_back(getSuperIndex(row, col - 1, L));
        neighbours.push_back(getSuperIndex(row, col + 1, L));
        neighbours.erase(remove(neighbours.begin(), neighbours.end(), -1), neighbours.end());
        neighbourArray.push_back(neighbours);
        cout << "POINT AT (ROW, COL): (" << row << ", " << col << ")" << endl; 
        for (int i: neighbours) {
            cout << i << ", ";
        }
        cout << endl;
        // neighbourArray[index][0] = getSuperIndex(row - 1, col, L);
        // neighbourArray[index][1] = getSuperIndex(row + 1, col, L);
        // neighbourArray[index][2] = getSuperIndex(row, col - 1, L);
        // neighbourArray[index][3] = getSuperIndex(row, col + 1, L);
        // cout << neighbourArray[index][0] << ",";
        // cout << neighbourArray[index][1] << ",";
        // cout << neighbourArray[index][2] << ",";
        // cout << neighbourArray[index][3] << endl;
        // neighbourArray[index].erase(remove(neighbourArray[index].begin(),
        //                                 neighbourArray[index].end(), -1),
        //                                 neighbourArray[index].end());
        // cout << neighbourArray[index][0] << ",";
        // cout << neighbourArray[index][1] << ",";
        // cout << neighbourArray[index][2] << ",";
        // cout << neighbourArray[index][3] << endl;
    }
  }
  return neighbourArray;
}

double delta_H(vector<int> latticeSpin, vector<vector<int>> neighbourArray,
               int latticePoint) {
  double spinProduct = 0;
  for (int neighbour : neighbourArray[latticePoint]) {
    spinProduct += latticeSpin[latticePoint] * latticeSpin[neighbour];
  }
  return -1 * spinProduct;
}

double pSwitch(double E, double beta) {
  return exp(-E * beta);
}

double getMagnetisation(vector<int> latticeSpin) {
  double M = 0.0;
  for (int i : latticeSpin) {
    M += i;
  }
  return M = M / (latticeSpin.size() * latticeSpin.size());
}

vector<int> updateLattice(vector<int> latticeSpin,
                          vector<vector<int>> neighbourArray, double beta) {
  random_device rd;
  mt19937 gen(rd());
  uniform_int_distribution<> dis_int(0, latticeSpin.size() - 1);
  uniform_real_distribution<> dis_real(0, 1);
  // pick position in lattice randomly and flip spin, L**2 times
  for (int i = 0; i < pow(latticeSpin.size(), 2); i++) {
    int latticePoint = dis_int(gen);
    latticeSpin[latticePoint] *= -1;

    // calculate delta H
    int E;
    E = delta_H(latticeSpin, neighbourArray, latticePoint);

    // acceptance step
    if (E >= 0) {
      double acceptance_ratio = pSwitch(E, beta);
      if (dis_real(gen) > acceptance_ratio) {
        latticeSpin[latticePoint] *= -1;
      }
    }
  }
  return latticeSpin;
}

void generate() {
  const int N = 2;
  const int Nskip = 1;
  const int L = 4;
  const double betaLower = 0.1;
  const double betaUpper = 0.9;
  const double betaStep = 0.001;
  // double beta = 0.5;
  const int NThermal = 20;
  const int coldStart = 1;
  // initialize lattice vectors
  vector<int> latticeSpin((L * L), coldStart);

  vector<vector<int>> neighbourArray(L * L, vector<int>(4, 0));

  neighbourArray = init(neighbourArray, L);
  for (double beta = betaLower; beta <= betaUpper; beta += betaStep) {
    cout << "beta: " << beta << endl;
    vector<vector<int>> latticeSpinData;
    vector<double> magnetization;
    for (int i = 0; i <= N; i++) {
      latticeSpin = updateLattice(latticeSpin, neighbourArray, beta);
      if (i % Nskip == 0) {
        double M = getMagnetisation(latticeSpin);
        magnetization.push_back(M);
        latticeSpinData.push_back(latticeSpin);
        cout << "i: " << Nskip << endl;
      }
    }
    // Writing Magnetisation File
    ofstream outFile;
    outFile.open("../data/magnetization" + to_string(beta * 1000) + ".txt");
    if (outFile.is_open()) {
      for (double i : magnetization) {
        outFile << i << ", ";
      }
    } else {
      cout << "Error opening file." << endl;
    }
    outFile.close();
  }
}

int main() {
  generate();
  return 0;
}