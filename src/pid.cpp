
#include "pid.h"
#include <algorithm>

Pid::Pid(double kp, double kd, double ki, double h, ControlMethod method,
         bool antiwindup, double u_low, double u_high, double Kh)
        : K_p(kp), K_d(kd), K_i(ki), h(h), u(0.0), diff_e(0.0), e_last(0.0),
         e_last_last(0.0), u_last(0.0), a0(0.0), a1(0.0), a2(0.0),
        u_low(u_low), u_high(u_high), aw(0.0), K_h(Kh), method(method), antiwindup(antiwindup) {
            calculate_pid_coefficients(); // Initializes coefficients
        }


    double Pid::calculate_control(double e) {

        if(antiwindup){
            aw = h * K_h * (v_last - u_last_last);
        } else {
            aw = 0.0;
        }
        u = u_last + a0 * e + a1 * e_last + a2 * e_last_last + aw;
     
        e_last_last = e_last;
        e_last = e;
        u_last_last = u_last;
        u_last = u;
        

        if(antiwindup){
            v_last = v;
            
            v = std::clamp(u, u_low, u_high);
        } else {
            v = u;
        }
        
        return v;
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

    
    