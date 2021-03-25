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

    //depth
    std::vector<std::vector<double>> temp_map_depth;
    temp_map_depth.reserve(Nx);
    for(int i = 0; i < Nx; i++){
        std::vector<double> result;
        result.reserve(Nx);
        for(int j = 0; j < Nz; j++){
            result.push_back(T0);
        }
        temp_map_depth.push_back(result);
    }

    auto calculate_depth = [](std::vector<std::vector<double>> &temp_map, int Nx0, int Nx1, int Nx, int Nz, Temperature_map A) {
        for (int i = Nx0; i < Nx1; i++) {
            for (int j = 0; j < Nz; j++) {
                temp_map[i][j] += A.integrate(i - Nx / 2, 0, j);
            }
        }
    };

    std::thread th1(calculate_depth, std::ref(temp_map_depth), 0, Nx/2, Nx, Nz, A);
    std::thread th2(calculate_depth, std::ref(temp_map_depth), Nx/2, Nx, Nx, Nz, A);

    th1.join();
    th2.join();

    int depth = 0;
    bool smallNz = false;

    for(int i = 0; i < Nx && !smallNz; i++){
        for(int j = 0; j < Nz; j++){
            if(temp_map_depth[i][j] < Tm){
                if(depth < j){
                    depth = j;
                }
                break;
            }
            else if(temp_map_depth[i][j] >= Tm && j+1 == Nz){
                smallNz = true;
            }
        }
    }

    //width and length
    std::vector<std::vector<double>> temp_map_width_length;
    temp_map_width_length.reserve(Nx);
    for(int i = 0; i < Nx; i++){
        std::vector<double> result;
        result.reserve(Nx);
        for(int j = 0; j < Nz; j++){
            result.push_back(T0);
        }
        temp_map_width_length.push_back(result);
    }

    auto calculate_width_length = [](std::vector<std::vector<double>> &temp_map, int Nx0, int Nx1, int Nx, int Ny, Temperature_map A) {
        for (int i = Nx0; i < Nx1; i++) {
            for (int j = 0; j < Ny; j++) {
                temp_map[i][j] += A.integrate(i - Nx / 2, j, 0);
            }
        }
    };

    std::thread th3(calculate_width_length, std::ref(temp_map_width_length), 0, Nx/2, Nx, Ny, A);
    std::thread th4(calculate_width_length, std::ref(temp_map_width_length), Nx/2, Nx, Nx, Ny, A);

    th3.join();
    th4.join();

    int width = 0;
    int right;
    bool smallNy = false;

    for(int i = 0; i < Nx && !smallNy; i++){
        right = -1;
        for(int j = 0; j < Ny; j++){
            if(temp_map_width_length[i][j] >= Tm && (j == 0 || j+1 == Ny)){
                smallNy = true;
            }
            else if(temp_map_width_length[i][j] >= Tm && right == -1){
                right = j;
            }
            else if(temp_map_width_length[i][j] < Tm && right != -1){
                width = j - right;
                break;
            }
        }
    }

    int length = 0;
    int tail = -1;
    bool smallNx;

    for(int i = 0; i < Nx; i++){
        if(temp_map_width_length[i][0] >= Tm && tail == -1){
            tail = i;
        }
        else if(temp_map_width_length[i][0] < Tm && tail != -1){
            length = i - tail;
            break;
        }
        else if(temp_map_width_length[i][0] >= Tm && (i == 0 || i+1 == Nx)){
            smallNx = true;
            break;
        }
    }


    if(smallNz){
        std::cout << "CAUTION: depth is equal to grid size, restart with bigger Nz\n";
    }
    else{
        std::cout << "Normalized depth: " << depth*dr/pow(D*a*pow(2, 0.5)/v, 0.5)
                  << ", depth in cells: " << depth << std::endl;
    }

    if(smallNy){
        std::cout << "CAUTION: width is equal to grid size, restart with bigger Ny\n";
    }
    else{
        std::cout << "Normalized width: " << width*dr/a
                  << ", width in cells: " << width << std::endl;
    }

    if(smallNx){
        std::cout << "CAUTION: length is equal to grid size, restart with bigger Nx\n";
    }
    else{
        std::cout << "Normalized length: " << length*dr/a
                  << ", length in cells: " << length << std::endl;
    }
    return 0;
}
