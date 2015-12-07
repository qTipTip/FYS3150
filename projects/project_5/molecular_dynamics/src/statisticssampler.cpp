#include "system.h"
#include "statisticssampler.h"
#include "potentials/potential.h"
#include "fstream"
#include "unitconverter.h"
StatisticsSampler::StatisticsSampler()
{

}

void StatisticsSampler::saveToFile(System &system)
{
    // Save the statistical properties for each timestep for plotting etc.

    using namespace std;
    ofstream outfile;
    outfile.open("samples.dat", ios::app);
    outfile << UnitConverter::temperatureToSI(m_temperature) << " "
            << UnitConverter::energyToEv(m_kineticEnergy/system.atoms().size()) << " "
            << UnitConverter::energyToEv(m_potentialEnergy/system.atoms().size()) << " "
            << UnitConverter::energyToEv(totalEnergy()) << " "
            << m_diffusionConstant << " "
            << system.time() << endl;
}

void StatisticsSampler::sample(System &system)
{
    // Here you should measure different kinds of statistical properties and save it to a file.
    sampleKineticEnergy(system);
    samplePotentialEnergy(system);
    sampleTemperature(system);
    sampleDensity(system);
    sampleMeanSquareDisplacement(system);
    sampleDiffusionConstant(system);
    saveToFile(system);
}

void StatisticsSampler::sampleKineticEnergy(System &system)
{
    m_kineticEnergy = 0; // Remember to reset the value from the previous timestep

    // for each atom we compute its mass times velocity squared, half it and add it to the counter
    for(Atom *atom : system.atoms()) {
        m_kineticEnergy += 0.5*atom->mass()*atom->velocity.lengthSquared();
    }
}

void StatisticsSampler::samplePotentialEnergy(System &system)
{
    // the potential is calculated in the Potential-class so
    // we simply extract it from there
    m_potentialEnergy = system.potential()->potentialEnergy();
}

void StatisticsSampler::sampleTemperature(System &system)
{
    //kb = 1
    int n_atoms = system.atoms().size();
    m_temperature = 2.0/3.0 * m_kineticEnergy / (n_atoms);
}



void StatisticsSampler::sampleDensity(System &system)
{

}

void StatisticsSampler::sampleMeanSquareDisplacement(System &system)
{
    // initialize to zero
    m_meanSquareDisplacement = 0.0;

    // for each atom we compute the difference between current position, and initial_position
    // need to make sure that we use true distance, this should be fixed in the applyBoundaryConditions method
    for(Atom *atom : system.atoms()){
        double diff_x = atom->position[0] - atom->initial_position[0];
        double diff_y = atom->position[1] - atom->initial_position[1];
        double diff_z = atom->position[2] - atom->initial_position[2];

        m_meanSquareDisplacement += diff_x * diff_x + diff_y * diff_y + diff_z * diff_z;
    }
    // divide by number of atoms, because we want the mean square displacement for each atom
    m_meanSquareDisplacement /= system.atoms().size();
}

void StatisticsSampler::sampleDiffusionConstant(System &system)
{
    // plotting without divding by time, as this might be more informative
    m_diffusionConstant = m_meanSquareDisplacement / (6*system.time());
}
