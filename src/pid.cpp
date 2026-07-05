
#include "pid.h"
#include <algorithm>

PidVelocity::PidVelocity(double K_p, double K_d, double K_i, double h, ControlMethod method,
         bool antiwindup, double u_low, double u_high)
        : K_p(K_p), K_d(K_d), K_i(K_i), h(h),
          u_low(u_low), u_high(u_high), method(method), antiwindup(antiwindup) {
            calculate_coefficients(); // Initializes coefficients
        }


    double PidVelocity::calculate_control(double e) {

        double u = u_last + a0 * e + a1 * e_last + a2 * e_last_last;
        
        // clamp control
        if(antiwindup){
            double v = std::clamp(u, u_low, u_high);
            u_last = v;
        } else {
            double v = u;
            u_last = v;
        }
        
        // store old states
        e_last_last = e_last;
        e_last = e;
        
        
        return u_last;
    }

    void PidVelocity::calculate_coefficients() {
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

    void PidVelocity::reset() {
        e_last_last = 0.0;
        e_last = 0.0;
        u_last = 0.0;
    }

    

PidPositional::PidPositional(double K_p, double K_i, double K_d, double h,
            double u_low, double u_high, double T_f, double T_aw, double beta)
            : K_p(K_p), K_d(K_d), K_i(K_i), h(h),
              u_low(u_low), u_high(u_high),
              T_f(T_f), T_aw(T_aw), beta(beta) {
                calculate_coefficients();
            }

    double PidPositional::calculate_control(double r, double y) {

        if (first_iter) {
            y_last = y;
            first_iter = false;
        }
        
        double P = K_p * (beta * r - y);
        D = ad * D - bd * (y-y_last);

        double u = P + I + D ;
        double v = std::clamp(u, u_low, u_high);

        I = I + bi * (r - y) + br * (v - u);
        y_last = y;

        return v;
    }

    void PidPositional::calculate_coefficients() {
        bi = K_i * h;
        ad = T_f / (T_f + h); 
        bd = K_d / (T_f + h);
        br = h / T_aw;
    }


    void PidPositional::reset() {
        first_iter = true;
        y_last = 0.0;
        D = 0.0;
        I = 0.0;
    }