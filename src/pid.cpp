
#include "pid.h"

Pid::Pid(double kp, double kd, double ki, double h, ControlMethod method)
    : K_p(kp), K_d(kd), K_i(ki), h(h), u(0.0), diff_e(0.0), summ_e(0.0), e_last(0.0), ulast(0.0), method(method) {

    }


    double Pid::calculate_control(double e) {
        switch (method){
            case ControlMethod::Backward_Euler:
                diff_e = e - e_last;
                
                u = u_last + K_p * diff_e + K_i * h * e + K_d * (e - 2 * e_last + e_last_last) / h
                
                e_last_last = e_last;
                e_last = e;
                u_last = u;

                return u;
                break;

            case ControlMethod::Forward_Euler:
                diff_e = e - e_last;
                
                u = u_last + K_p * diff_e + K_i * h * e_last + K_d * (e - 2 * e_last + e_last_last) / h
                
                e_last_last = e_last;
                e_last = e;
                u_last = u;

                return u;
                break;
            
            case ControlMethod::Trapezoidal:
                diff_e = e - e_last;

                alpha = (2 * K_d)/h;
                beta = (K_i * h)/2;
                a0 = K_p + alpha + beta 
                a1 = -2 * K_p + beta - alpha
                a2 = K_p + alpha - beta 
                u = u_last + K_p * diff_e + K_i * h * e + K_d * (e - 2 * e_last + e_last_last) / h
                
                e_last_last = e_last;
                e_last = e;
                u_last = u;

                return u;
                break;

            default:
                diff_e = (e - e_last) / h;
                summ_e =+ e;
                u = K_p * e + K_d * diff_e + K_i * h * summ_e;

                e_last = e;

                
                break;
            
            }
            return u;
        }

    
    