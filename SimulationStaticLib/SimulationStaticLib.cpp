#include "SimulationStaticLib.h"
#include "Simulator.h"
SimulationStaticLib::SimulationStaticLib()
{
}

void SimulationStaticLib::imp_threeBodyFast(
            const Eigen::Map<const BodyVector> & mass,
            const Eigen::TensorMap<const Position> & begPos,
            const Eigen::TensorMap<const Velocity> & begVelocity,
            const TimeSpan & ts,
            const double precison,
            bool * noCollide,
            double * lastTime
            )
{
    Simulator s;

    s.setMass(mass*Ms);
    s.simulateRK4Var1<false>(1e-4*year,ts,Statue(begPos*rs,begVelocity*vs),noCollide,precison);
    *lastTime=std::min(s.getResult().back().first,ts.second)/year;
}

