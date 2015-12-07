#ifndef POTENTIAL_H
#define POTENTIAL_H

#include <string>
#include <vector>
#include "../system.h"

class Potential
{
protected:
    double m_potentialEnergy = 0;
    double m_rCut = -1;
    double m_rShell = -1;
public:
    Potential();
    virtual ~Potential() {}
    virtual void calculateForces(System *system) = 0;
    double potentialEnergy();
    void setPotentialEnergy(double potentialEnergy);
    double rCut() const;
    void setRCut(double rCut);
    double rShell() const;
    void setRShell(double rShell);
    void buildNeighborList(System &system);
};
#endif
