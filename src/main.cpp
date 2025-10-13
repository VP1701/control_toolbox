#include <iostream>
#include "pid.h"

int main() {
    Pid Pid(2.0, 1.0, 0.0, 0.01, ControlMethod::Trapezoidal);

    double error = 5.0;

    double control1 = Pid.calculate_control(error);
    double control2 = Pid.calculate_control(error);
    double control3 = Pid.calculate_control(error);

    std::cout << "Calculated control1: " << control1 << std::endl;
    std::cout << "Calculated control2: " << control2 << std::endl;
    std::cout << "Calculated control3: " << control3 << std::endl;
}