
#include "pid.h"

Pid::Pid(double kp, double kd, double ki, double h, ControlMethod method)
    : K_p(kp), K_d(kd), K_i(ki), h(h), u(0.0), diff_e(0.0), summ_e(0.0),
      e_last(0.0), e_last_last(0.0), u_last(0.0), a0(0.0), a1(0.0), a2(0.0), method(method) {
        calculate_pid_coefficients(); // Initializes coefficients
    }


    double Pid::calculate_control(double e) {
        
                
        u = u_last + a0 * e + a1 * e_last + a2 * e_last_last;
                
        e_last_last = e_last;
        e_last = e;
        u_last = u;

        return u;
    }
    void Pid::calculate_pid_coefficients() {
        switch (method){
            case ControlMethod::Backward_Euler:
                a0 = K_p + K_i * h + K_d / h; 
                a1 = -K_p - 2 * (K_d/h);
                a2 = K_d /h;
                break;

            case ControlMethod::Forward_Euler:
                a0 = K_p+ K_d / h;
                a1 = -K_p + K_i * h - 2 * (K_d/h);
                a2 = K_d /h;
                break;
            
            case ControlMethod::Trapezoidal:
            case ControlMethod::Tustin:
                a0 = K_p + (K_i * h) / 2 + (2 * K_d) / h;
                a1 = -K_p + (K_i * h) / 2 - 2 * (K_d/h);
                a2 = K_d /h;
                break;
        }

    }

    
    