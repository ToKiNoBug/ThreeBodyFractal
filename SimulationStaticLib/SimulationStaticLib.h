#ifndef SIMULATIONSTATICLIB_H
#define SIMULATIONSTATICLIB_H

#include "defines.h"
class SimulationStaticLib
{
public:
    SimulationStaticLib();

    static void imp_threeBodyFast(
            const Eigen::Map<const BodyVector> & mass,
            const Eigen::TensorMap<const Position> & begPos,
            const Eigen::TensorMap<const Velocity> & begVelocity,
            const TimeSpan & ts,
            const double precison,
            bool * noCollide,
            double * lastTime
            );
};

#endif // SIMULATIONSTATICLIB_H
