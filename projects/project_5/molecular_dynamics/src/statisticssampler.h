#ifndef STATISTICSSAMPLER_H
#define STATISTICSSAMPLER_H

class System;
class StatisticsSampler
{
private:
    double m_kineticEnergy = 0;
    double m_potentialEnergy = 0;
    double m_temperature = 0;
    double m_density = 0;
    double m_meanSquareDisplacement = 0;
    double m_diffusionConstant = 0;
public:
    StatisticsSampler();
    void saveToFile(System &system);
    void sample(System &system);
    void sampleKineticEnergy(System &system);
    void samplePotentialEnergy(System &system);
    void sampleTemperature(System &system);
    void sampleDensity(System &system);
    void sampleMeanSquareDisplacement(System &system);
    void sampleDiffusionConstant(System &system);
    double kineticEnergy() { return m_kineticEnergy; }
    double potentialEnergy() { return m_potentialEnergy; }
    double totalEnergy() { return m_kineticEnergy+m_potentialEnergy; }
    double temperature() { return m_temperature; }
    double density() { return m_density; }
};
#endif
