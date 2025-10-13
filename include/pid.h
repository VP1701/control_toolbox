// pid.h
#ifndef PID_H
#define PID_H

#include <iostream>
enum class ControlMethod {
    Forward_Euler,
    Backward_Euler,
    Trapezoidal,
    Tustin,
};

class Pid {
    public:
        Pid(double kp, double kd, double ki, double h, ControlMethod method);

        double calculate_control(double e);
        void calculate_pid_coefficients();

    private:
        double K_p;
        double K_d;
        double K_i;
        double h;
        double u;
        double diff_e;
        double summ_e;
        double e_last;
        double e_last_last;
        double u_last;
        double a0;
        double a1;
        double a2;
        ControlMethod method;
        

};


#endif // PID_H