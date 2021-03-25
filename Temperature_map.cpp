#include "Temperature_map.h"


Temperature_map::Temperature_map(double dr, double dt, int Nx, int Ny, int Nz,
                                 double A, double D, double k, double cp,
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
    this->cp = cp;

    this->P = P;
    this->a = a;
    this->v = v;
    this->T0 = T0;
    this->t_max_steps = t_max_steps;

    this->factor = A*P/(M_PI*cp*sqrt(4*M_PI*D));
    this->dr_square = dr*dr;
}

Temperature_map::~Temperature_map() = default;

double Temperature_map::integrate(int x, int y, int z){
    double integral = 0;
    double y_square = y*y*dr_square;
    double z_square = z*z*dr_square;
    for(int t = 1; t < t_max_steps; t++){
        double time = t*dt;
        double divider = 1/(2*D*time + a*a);
        double x_new = x*dr+v*time;
        double power = (x_new*x_new + y_square)*0.5*divider + z_square/(4*D*time);
        double result =  exp(-power);
        result *= dt/sqrt(time) * divider;
        integral += result;
    }
    integral *= factor;
    return integral;
}
