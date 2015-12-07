#include "math/random.h"
#include "potentials/lennardjones.h"
#include "integrators/eulercromer.h"
#include "integrators/velocityverlet.h"
#include "system.h"
#include "statisticssampler.h"
#include "atom.h"
#include "io.h"
#include "unitconverter.h"
#include <iostream>
#include <QElapsedTimer>

using namespace std;

int main(int numberOfArguments, char **argumentList)
{
    int numberOfUnitCells = 10;
    double initialTemperature = UnitConverter::temperatureFromSI(2000); // measured in Kelvin
    double latticeConstant = UnitConverter::lengthFromAngstroms(5.26); // measured in angstroms
    double dt = UnitConverter::timeFromSI(1.0e-15); // measured in seconds

    // If a first argument is provided, it is the number of unit cells
    if(numberOfArguments > 1) numberOfUnitCells = atoi(argumentList[1]);
    // If a second argument is provided, it is the initial temperature (measured in kelvin)
    if(numberOfArguments > 2) initialTemperature = UnitConverter::temperatureFromSI(atof(argumentList[2]));
    // If a third argument is provided, it is the lattice constant determining the density (measured in angstroms)
    if(numberOfArguments > 3) latticeConstant = UnitConverter::lengthFromAngstroms(atof(argumentList[3]));
    // time step in seconds from the command line
    if(numberOfArguments > 5) dt = UnitConverter::timeFromSI(atof(argumentList[5]));

    cout << "One unit of length is " << UnitConverter::lengthToSI(1.0) << " meters" << endl;
    cout << "One unit of velocity is " << UnitConverter::velocityToSI(1.0) << " meters/second" << endl;
    cout << "One unit of time is " << UnitConverter::timeToSI(1.0) << " seconds" << endl;
    cout << "One unit of mass is " << UnitConverter::massToSI(1.0) << " kg" << endl;
    cout << "One unit of temperature is " << UnitConverter::temperatureToSI(1.0) << " K" << endl;

    System system;
    system.createFCCLattice(numberOfUnitCells, latticeConstant, initialTemperature);
    system.setPotential(new LennardJones(3.405, 1.0)); // You must insert correct parameters here
    system.potential()->setRCut(2.5*3.405);
    system.potential()->setRShell(0.3*3.405);

    // can now chose integration method from command line
    int integration_method = 0;
    if(numberOfArguments > 4 ) integration_method = atoi(argumentList[4]);
    switch(integration_method) {
        case 1:
            system.setIntegrator(new EulerCromer);
            break;
        default:
            system.setIntegrator(new VelocityVerlet);
            break;
    }

    system.removeTotalMomentum();

    StatisticsSampler statisticsSampler;
    IO movie; // To write the state to file
    movie.open("movie.xyz");

    cout << "Timestep Time Temperature KineticEnergy PotentialEnergy TotalEnergy" << endl;
    QElapsedTimer elapsedTimer;
    elapsedTimer.start();
    int timesteps = 2000;
    for(int timestep=0; timestep<timesteps; timestep++) {
        system.step(dt);
        statisticsSampler.sample(system);
        if( !(timestep % 100) ) {
            // Print the timestep every 100 timesteps
            cout << system.steps() << "      " << system.time() << "      " << UnitConverter::temperatureToSI(statisticsSampler.temperature()) << "      " << UnitConverter::energyToEv(statisticsSampler.kineticEnergy()) << "      " << UnitConverter::energyToEv(statisticsSampler.potentialEnergy()) << "      " << UnitConverter::energyToEv(statisticsSampler.totalEnergy()) << endl;
        }
        movie.saveState(&system);
    }
    cout << "Simulation finished after " << elapsedTimer.elapsed()/1000.0 << " seconds." << endl;
    cout << 1000.0*timesteps*system.atoms().size() / elapsedTimer.elapsed() << " atom timesteps / sec" << endl;


    movie.close();

    return 0;
}
