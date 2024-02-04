#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <random>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

using namespace std;

struct data_row {
  float magnetization;
  float beta;
  int i;
};

bool comparebyItr(const data_row &r1, const data_row &r2) {
  if (r1.beta < r2.beta) {
    return true;
  }
  if (r2.beta < r1.beta) {

    return false;
  }
  if (r1.i < r2.i) {

    return true;
  }
  if (r2.i < r1.i) {

    return false;
  }
  return false;
}

int getSuperIndex(int row, int col, const int L) {
  if (row < 0 || row >= L) {
    return -1;
  } else if (col < 0 || col >= L) {
    return -1;
  } else {
    return row * L + col;
  }
}

double drawRandomFloat() {
  static mt19937 gen_real = []() {
    random_device rd;
    return mt19937(rd());
  }();
  uniform_real_distribution<float> dis(0, 1);
  return dis(gen_real);
}

double drawRandomInt() {
  static mt19937 gen_int = []() {
    random_device rd;
    return mt19937(rd());
  }();
  uniform_int_distribution dis(0, 899);
  return dis(gen_int);
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

int delta_H(vector<int> &latticeSpin, int &latticePointSpin,
               vector<int> &neighbours, int &latticePoint) {
  int spinProduct = 0;
  for (int neighbour : neighbours) {
    spinProduct += latticePointSpin * latticeSpin[neighbour];
  }
  return -1 * spinProduct;
}

float pSwitch(int &E, float &beta) { return exp(-E * beta); }

float getMagnetisation(vector<int> &latticeSpin) {
  float M = 0.0;
  for (int i : latticeSpin) {
    M += i;
  }
  return M = M / (latticeSpin.size() * latticeSpin.size());
}

int getLocation(int kbeta, int i, const int N, const int Nskip) {
  return ((kbeta - 100) * N + i) / Nskip;
}

void updateLattice(vector<int> &latticeSpin,
                   vector<vector<int>> &neighbourArray, float &beta,
                   int &latticeSpinLength) {
  // pick position in lattice randomly and flip spin, L**2 times
  for (int i = 0; i < latticeSpinLength; i++) {
      // draw random float and convert to lattice point
    float random_num = drawRandomFloat();
    int latticePoint = random_num * (latticeSpinLength - 1);

    // flip random point in lattice and find its neighbours
    latticeSpin[latticePoint] *= -1;

    // calculate delta H
    int E = delta_H(latticeSpin, latticeSpin[latticePoint],
                    neighbourArray[latticePoint], latticePoint);

    // acceptance step
    if (E >= 0) {
      if (random_num > pSwitch(E, beta)) {
        latticeSpin[latticePoint] *= -1;
      }
    }
  }
  return;
}

void generate() {
  const int N = 10000;
  const int Nskip = 100;
  const int L = 30;
  const int betaLower = 100;
  const int betaUpper = 900;
  const int betaStep = 1;
  const int NThermal = 20;
  const int coldStart = 1;

  // create result vector
  int vector_size = (N / Nskip) * (betaUpper - betaLower) + Nskip + 1;
  cout << vector_size << endl;
  vector<data_row> results(vector_size);

#pragma omp parallel for
  for (int kbeta = betaLower; kbeta <= betaUpper; kbeta += betaStep) {
    // convert kbeta to beta
    float beta = kbeta / 1000.0;
    cout << "beta: " << beta << '\n';

    // initialize lattice vectors
    int latticeSpinLength = L * L;
    vector<int> latticeSpin((latticeSpinLength), coldStart);
    vector<vector<int>> neighbourArray(latticeSpinLength, vector<int>(4, 0));
    neighbourArray = init(neighbourArray, L);

    // main monte carlo simulation
    for (int i = 0; i <= N; i++) {
      updateLattice(latticeSpin, neighbourArray, beta, latticeSpinLength);
      if (i % Nskip == 0) {
        data_row row;
        row.magnetization = getMagnetisation(latticeSpin);
        row.beta = beta;
        row.i = i;
        int vec_loc = getLocation(kbeta, i, N, Nskip);
        results[vec_loc] = row;
      }
    }
  }
  cout << "Ising Model Calculations Complete"
       << "\n";

  // sort beta results vector
  std::sort(results.begin(), results.end(), comparebyItr);

  // Open Magnetisation File
  ofstream outFile;
  outFile.open("data/magnetization/magnetization.csv");

  if (!outFile.is_open()) {
    throw exception("Error opening file");
  } else {
    cout << "Writing file";
    // write contents of results vector to csv
    outFile << "magnetisation, beta, index" << '\n';
    for (data_row row : results) {
      outFile << row.magnetization << ", " << row.beta << ", " << row.i << "\n";
    }
  }
  outFile.close();
}

int main() {
  generate();
  return 0;
}