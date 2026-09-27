//status.hpp

#ifndef STATUS_HPP
#define STATUS_HPP

enum class [[nodiscard]] Status {
    Ok,
    DimensionMismatch, 
    NotPosDef, 
    NotInitialized, 
    IterationCountReached, 
    NumericalError,
    Infeasible
};

#endif