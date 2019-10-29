#ifndef TEMPERATURE_MAP_H
#define TEMPERATURE_MAP_H

#include <cmath>

/*
    The solution of heat equation with temperature-independent parameters by the Green function method
    (from Gladush, Smurov "Physics of Laser Material Processing", 2011, par. 2.1)
*/

class Temperature_map{

private:
    double dr;             // space step
    double dt;             // time step
    int Nx;                // number of cell, x dimension
    int Ny;
    int Nz;
    int current_t;         // current time step
    //material properties
    double A;              // absorption coefficient
    double D;              // thermal diffusivity
    double k;              // heat conductivity
    //experiment parameters
    double P;              // laser power
    double a;              // beam size (exp(-r^2/a^2))
    double v;              // velocity
    double T0;             // T(x == y == z == inf)
    double t_max;          // time of experiment

public:
    double*** temp_map;

public:
    Temperature_map(double dr, double dt, int Nx, int Ny, int Nz,
             double A, double D, double k,
             double P, double a, double v, double T0, double t_max);
    ~Temperature_map();


    void update();                         // updates temperature map; I think must be private

private:
    double calc_dT(int x, int y, int z);   // x, y, z -- number in steps (actual value divided by step)
    double calc_dT_2(int x, int y, int z);   // x, y, z -- number in steps (actual value divided by step)
};

#endif //TEMPERATURE_MAP_H