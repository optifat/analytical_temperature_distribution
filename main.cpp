#include "Temperature_map.h"

#include <iostream>

int main(){
    double dt = 1*pow(10, -8);
    Temperature_map A(5*pow(10, -6), dt, 10, 10, 10,
                       0.24, 5*pow(10, -6), 24,
                       100, 25*sqrt(2)*pow(10, -6), 2, 300, 2500*dt);


    for(int i = 0; i<10; i++){
        for(int j = 0; j<10; j++){
            std::cout << A.temp_map[i][j][0] << " ";
        }
        std::cout<<std::endl;
    }
    std::cout << std::endl;
    /*
    for(int i = 0; i<10; i++){
        for(int j = 0; j<10; j++){
            for(int k = 0; k<10; k++){
                A.temp_map[i][j][k] = A.integrate(i, j, k, (1/(1000*dt)));
            }
        }
    }
    */

    for(int i = 0; i<2499; i++){
        A.update();
    }

    for(int i = 0; i<10; i++){
        for(int j = 0; j<10; j++){
            std::cout << A.temp_map[i][j][0] << " ";
        }
        std::cout << std::endl;
    }
}
