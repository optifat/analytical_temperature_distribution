#include "Temperature_map.h"

#include <iostream>
#include <sstream>

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
    double absorb = 0.24;
    double D = 5*pow(10, -6);
    double k = 24;
    int Tm = 1563;
    double a = 50*pow(10, -6);
    double v = 0.8;
    double T0 = 293;

    Temperature_map A(dr, dt, Nx, Ny, Nz,
                      absorb, D, k,
                      p, a, v, T0, 10000);

    auto temp_map = new double*[Nx];
    for(int i = 0; i < Nx; i++){
        temp_map[i] = new double[Nz];
    }

    for(int i = 0; i < Nx; i++){
        for(int j = 0; j < Nz; j++){
                temp_map[i][j] = 293 + A.integrate(i - Nx/2, 0, j);
        }
    }

    int depth = 0;

    for(int i = 0; i < Nx; i++){
        for(int j = 0; j < Nz; j++){
            //std::cout << temp_map[i][j] << " ";
            if(temp_map[i][j] < Tm){
                //std::cout << j << '\n';
                if(depth < j){
                    depth = j;
                }
                break;
            }
        }
        //std::cout << '\n';
    }

    for(int i = 0; i < Nz; i++){
        delete[] temp_map[i];
    }
    delete[] temp_map;

    std::cout << depth*dr/pow(D*a*pow(2, 0.5)/v, 0.5) << ' ' << depth << std::endl;

    return 0;
}
