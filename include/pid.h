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

class PidVelocity {
    public:
        PidVelocity(double K_p, double K_d, double K_i, double h,
                    ControlMethod method, bool antiwindup = false, double u_low = 0.0, double u_high = 0.0);

        double calculate_control(double e);     // Method that calculates the next control.
        void calculate_coefficients();      // Method for calculating PID coefficients a0, a1 and a2 based on the discretization method.
        void reset(); // Resets the controller state
    private:

        double e; // error



        double K_p; // Proportional gain
        double K_d; // Derivative gain
        double K_i; // integral gain
        double h;   // hold time

        double u_last = 0.0;              // Unclamped control of time step k-1
        double e_last = 0.0;              // Error of time step k-1
        double e_last_last = 0.0;         // Error of time step k-2

        // The values of these gain terms depend on the discretization method.
        double a0;                  // Gain term for e
        double a1;                  // Gain term for e_last
        double a2;                  // Gain term for e_last_last

        bool antiwindup;            // Antiwindup enabled or not
        double u_low;               // Lower bound for antiwindup control clamping
        double u_high;              // Upper bound for antiwindup control clamping

        ControlMethod method;
        

};



class PidPositional {
    // Postional PID algorithm from Feedback systems seond editions. Includes antiwindup and derivative filter
    public:
        PidPositional(double K_p, double K_i, double K_d, double h,
                      double u_low, double u_high, double T_f, double T_aw, double beta = 0.0);

        double calculate_control(double r, double y);     // Method that calculates the next control.
        void calculate_coefficients();      // Method for calculating PID coefficients a0, a1 and a2 based on the discretization method.
        void reset(); // Resets the controller state

    private:
        double y_last;
        bool first_iter = true;
        double beta; // setpoint weight
        double K_p, K_i, K_d; // Pid parematers
        double bi, ad, bd, br; // Controller coefficients
        double h; // Timestep
        double T_f; // 
        double T_aw; // Antiwindup gain  
        double I = 0.0;
        double D = 0.0;
        double u_low;               // Lower bound for antiwindup control clamping
        double u_high;              // Upper bound for antiwindup control clamping
};


#endif // PID_H