//status.hpp

#ifndef STATUS_HPP
#define STATUS_HPP

#include <iostream>

enum class [[nodiscard]] Status {
    Ok,
    DimensionMismatch, 
    NotPosDef, 
    NotInitialized, 
    IterationCountReached, 
    NumericalError,
    Infeasible,
    FileError,    
    ParseError,   
    NotFound      
};


inline void status_check(Status status) {
    if (status != Status::Ok) {
        std::cout << "Status not ok" << "\n";
    }
}
#endif

