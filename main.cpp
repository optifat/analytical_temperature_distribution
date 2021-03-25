#include "Temperature_map.h"

#include <iostream>
#include <sstream>
#include <vector>
#include <thread>

int main(int argc, char* argv[]){
    std::stringstream stream;
    stream << argv[1];
    double p;
    stream >> p;

    double dr = 0.5*pow(10, -6);
    double dt = 60*pow(10, -9);
    int Nx = 400;
    int Ny = 100;
    int Nz = 150;
    double absorb = 0.25;
    double D = 5*pow(10, -6);
    double k = 26;
    double cp = 5.2*pow(10, 6);
    int Tm = 1563;
    double a = 25*pow(10, -6);
    double v = 0.1;
    double T0 = 300;

    //double p = h/absorb*(k/D)*Tm*pow(M_PI*M_PI*M_PI*D*v*a*a*a, 0.5);

    Temperature_map A(dr, dt, Nx, Ny, Nz,
                      absorb, D, k, cp,
                      p, a, v, T0, 10000);

    std::vector<std::vector<double>> temp_map;
    for(int i = 0; i < Nx; i++){
        std::vector<double> result;
        for(int j = 0; j < Nz; j++){
            result.push_back(T0);
        }
        temp_map.push_back(result);
    }

    auto calculate = [](std::vector<std::vector<double>> &temp_map, int Nx0, int Nx1, int Nx, int Nz, Temperature_map A) {
        for (int i = Nx0; i < Nx1; i++) {
            for (int j = 0; j < Nz; j++) {
                temp_map[i][j] += A.integrate(i - Nx / 2, 0, j);
            }
        }
    };

    std::thread th1(calculate, std::ref(temp_map), 0, Nx/2, Nx, Nz, A);
    std::thread th2(calculate, std::ref(temp_map), Nx/2, Nx, Nx, Nz, A);

    th1.join();
    th2.join();

    int depth = 0;

    for(int i = 0; i < Nx; i++){
        for(int j = 0; j < Nz; j++){
            if(temp_map[i][j] < Tm){
                if(depth < j){
                    depth = j;
                }
                break;
            }
        }
    }

    std::cout << depth*dr/pow(D*a*pow(2, 0.5)/v, 0.5) << ' ' << depth << std::endl;

    return 0;
}
