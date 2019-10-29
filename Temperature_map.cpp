#include "Temperature_map.h"

#include <iostream>

Temperature_map::Temperature_map(double dr, double dt, int Nx, int Ny, int Nz,
                   double A, double D, double k,
                   double P, double a, double v, double T0, double t_max){

    this->dr = dr;
    this->dt = dt;
    this->Nx = Nx;
    this->Ny = Ny;
    this->Nz = Nz;
    this->current_t = 0;

    this->A = A;
    this->D = D;
    this->k = k;

    this->P = P;
    this->a = a;
    this->v = v;
    this->T0 = T0;
    this->t_max = t_max;

    temp_map = new double**[Nx];
    for(int i = 0; i < Nx; i++){
        temp_map[i] = new double*[Ny];
        for(int j = 0; j < Ny; j++){
            temp_map[i][j] = new double[Nz];
        }
    }

    for(int i = 0; i < Nx; i++){
        for(int j = 0; j < Ny; j++){
            for(int l = 0; l < Nz; l++){
                temp_map[i][j][l] = this->T0;
            }
        }
    }

}

Temperature_map::~Temperature_map(){
    for(int i = 0; i < Nx; i++){
        for(int j = 0; j < Ny; j++){
            delete temp_map[i][j];
        }
        delete temp_map[i];
    }
    delete temp_map;
}

void Temperature_map::update() {
    current_t++;
    for(int i = 0; i < Nx; i++){
        for(int j = 0; j < Ny; j++){
            for(int l = 0; l < Nz; l++){
                temp_map[i][j][l] += calc_dT_2(i, j, l);
            }
        }
    }
}

double Temperature_map::calc_dT(int x, int y, int z) {
    double power = (pow(x*dr+v*current_t*dt, 2)+y*dr*y*dr)/(4*D*current_t*dt+a*a) + z*dr*z*dr/(4*D*current_t*dt);
    double result =  exp(-power);
    result *= dt/(sqrt(current_t*dt)*(current_t*dt+a*a/4/D));
    return A*P/(2*M_PI*k*sqrt(4*M_PI*D))*result;
}

double Temperature_map::calc_dT_2(int x, int y, int z) {
    double power = (pow(x*dr-v*current_t*dt, 2)+y*dr*y*dr)/(4*D*(t_max-current_t*dt)+a*a)
                    + z*dr*z*dr/(4*D*(t_max - current_t*dt));
    double result =  exp(-power);
    result *= dt/(sqrt(t_max - current_t*dt)*(t_max - current_t*dt+a*a/4/D));
    return A*P/(2*M_PI*k*sqrt(4*M_PI*D))*result;
}
