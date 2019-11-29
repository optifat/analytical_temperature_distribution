#include "Temperature_map.h"

#include <iostream>

Temperature_map::Temperature_map(double dr, double dt, int Nx, int Ny, int Nz,
                                 double A, double D, double k,
                                 double P, double a, double v, double T0, double t_max_steps){

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
    this->t_max_steps = t_max_steps;

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
            delete[] temp_map[i][j];
        }
        delete[] temp_map[i];
    }
    delete[] temp_map;
}

void Temperature_map::update() {
    current_t++;
    for(int i = 0; i < Nx; i++){
        for(int j = 0; j < Ny; j++){
            for(int l = 0; l < Nz; l++){
                temp_map[i][j][l] += calc_dT(i, j, l);
            }
        }
    }
}

double Temperature_map::calc_dT(int x, int y, int z) {
    double power = (pow(x*dr-v*current_t*dt, 2)+y*dr*y*dr)/(4*D*(t_max_steps-current_t)*dt+a*a)
                    + z*dr*z*dr/(4*D*(t_max_steps - current_t)*dt);
    double result =  exp(-power);
    result *= dt/(sqrt((t_max_steps - current_t)*dt)*((t_max_steps - current_t)*dt+a*a/4/D));
    return A*P/(2*M_PI*k*sqrt(4*M_PI*D))*result;
}

double Temperature_map::integrate(int x, int y, int z){
    double integral = 0;
    for(int t = 1; t < t_max_steps; t++){
        double power = (pow((x*dr+v*t*dt),2) + y*y*dr*dr)/(4*D*t*dt + 2*a*a) + z*z*dr*dr/(4*D*t*dt);
        double result =  exp(-power);
        result *= dt/sqrt(t*dt)/(2*D*t*dt + a*a) * A*P/(M_PI*k/D*sqrt(4*M_PI*D));
        integral += result;
    }
    return integral;
}

double Temperature_map::integrate_2(int x, int y, int z){
    a = this->a*pow(2, 0.5);
    double p = D/(v*a);

    double integral = 0;
    for(int t = 1; t < t_max_steps; t++){
        double power = (pow((x*dr+t*dt),2) + y*y*dr*dr)/(4*p*t*dt + 1) + z*z*dr*dr/(4*t*dt);
        double result =  exp(-power);
        result *= dt/sqrt(t*dt)/(4*p*t*dt + 1) * A*P/(M_PI*k/D*sqrt(M_PI*D*v*a*a*a));
        integral += result;
    }
    return integral;
}
