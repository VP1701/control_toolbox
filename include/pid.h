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
        Pid(double kp, double kd, double ki, double h, ControlMethod method, bool antiwindup = false, double u_low = 0.0, double u_high = 0.0, double Kh = 0.0);

        double calculate_control(double e);     // Method that calculates the next control.
        void calculate_pid_coefficients();      // Method for calculating PID coefficients a0, a1 and a2 based on the discretization method.

    private:
        double K_p; // Proportional gain
        double K_d; // Derivative gain
        double K_i; // integral gain
        double h;   // hold time
        double K_h; // antiwindup gain

        double u;                   // Unclamped control value or just control value if antiwindup = false
        double u_last;              // Unclamped control of time step k-1
        double u_last_last;         // Unclamped control of time step k-2

        double v;                   // Clamped control value 
        double v_last;              // Clamped control value of time step k-1
        double v_last_last;         // Clamped control value of time step k-2

        double diff_e;              // Difference between current and last error
        double e_last;              // Error of time step k-1
        double e_last_last;         // Error of time step k-2
        
        double aw;

        // The values of these gain terms depend on the discretization method.
        double a0;                  // Gain term for e
        double a1;                  // Gain term for e_last
        double a2;                  // Gain term for e_last_last

        bool antiwindup;            // Antiwindup enabled or not
        double u_low;               // Lower bound for antiwindup control clamping
        double u_high;              // Upper bound for antiwindup control clamping

        ControlMethod method;
        

};


#endif // PID_H