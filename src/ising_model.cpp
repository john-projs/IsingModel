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
  } else if (col < 0 || col >= L) {
    return -1;
  } else {
    return row * L + col;
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
      neighbours.erase(remove(neighbours.begin(), neighbours.end(), -1),
                       neighbours.end());
      neighbourArray.push_back(neighbours);
    }
  }
  return neighbourArray;
}

double delta_H(vector<int> latticeSpin, vector<int> neighbours,
               int latticePoint) {
  double spinProduct = 0;
  for (int neighbour : neighbours) {
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
  for (int i = 0; i < latticeSpin.size(); i++) {
    int latticePoint = dis_int(gen);
    latticeSpin[latticePoint] *= -1;
    vector<int> neighbours = neighbourArray[latticePoint];
    // calculate delta H
    int E;
    E = delta_H(latticeSpin, neighbours, latticePoint);

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
  const int N = 10000;
  const int Nskip = 100;
  const int L = 30;
  const double betaLower = 0.1;
  const double betaUpper = 0.9;
  const double betaStep = 0.001;
  const int NThermal = 20;
  const int coldStart = 1;
  // initialize lattice vectors
  vector<int> latticeSpin((L * L), coldStart);

  vector<vector<int>> neighbourArray(L * L, vector<int>(4, 0));

  // Open Magnetisation File
  ofstream outFile;
  outFile.open("../data/magnetization/magnetization.csv");
  outFile << "magnetisation, beta" << endl;
  neighbourArray = init(neighbourArray, L);
  if (outFile.is_open()) {
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
          cout << "i: " << i << endl;
        }
      }

      for (double i : magnetization) {
        outFile << i << ", " << beta << endl;
      }
    }
  }
  else {
    cout << "Error opening file." << endl;
  }
  outFile.close();
}

  int main() {
    generate();
    return 0;
  }