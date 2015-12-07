#include "potential.h"

double Potential::rCut() const
{
    return m_rCut;
}

void Potential::setRCut(double rCut)
{
    m_rCut = rCut;
}

double Potential::rShell() const
{
    return m_rShell;
}

void Potential::setRShell(double rShell)
{
    m_rShell = rShell;
}

Potential::Potential()
{
    
}

double Potential::potentialEnergy()
{
    return m_potentialEnergy;
}

void Potential::setPotentialEnergy(double potentialEnergy)
{
    m_potentialEnergy = potentialEnergy;
}

void Potential::buildNeighborList(System &system) {
    system.cellList().cutoffDistance = m_rCut + m_rShell;
    system.cellList().build(system);
    for(Atom *atom : system.atoms()) {
        atom->neighbors.clear();
        atom->positionLastNeighborlistBuild = atom->position;
    }

    CellList &cellList = system.cellList();
    vec3 systemSize = system.systemSize();
    double radiusSquared = (m_rCut + m_rShell)*(m_rCut + m_rShell);

    for(int cx = 0; cx<cellList.numberOfCellsPerDimension; cx++) {
        for(int cy = 0; cy<cellList.numberOfCellsPerDimension; cy++) {
            for(int cz = 0; cz<cellList.numberOfCellsPerDimension; cz++) {
                vector<Atom*> &cell1 = cellList.cell(cx, cy, cz);

                for(int dx=-1; dx<=1; dx++) {
                    for(int dy=-1; dy<=1; dy++) {
                        for(int dz=-1; dz<=1; dz++) {
                            int i = (cx + dx + cellList.numberOfCellsPerDimension) % cellList.numberOfCellsPerDimension;
                            int j = (cy + dy + cellList.numberOfCellsPerDimension) % cellList.numberOfCellsPerDimension;
                            int k = (cz + dz + cellList.numberOfCellsPerDimension) % cellList.numberOfCellsPerDimension;

                            vector<Atom*> &cell2 = cellList.cell(i, j, k);
                            for(Atom *atom1 : cell1) {
                                for(Atom *atom2 : cell2) {
                                    if(atom1 <= atom2) continue;
                                    double dx = atom1->position[0] - atom2->position[0];
                                    double dy = atom1->position[1] - atom2->position[1];
                                    double dz = atom1->position[2] - atom2->position[2];

                                    if(dx < -0.5*systemSize[0]) dx += systemSize[0];
                                    else if(dx > 0.5*systemSize[0]) dx -= systemSize[0];

                                    if(dy < -0.5*systemSize[1]) dy += systemSize[1];
                                    else if(dy > 0.5*systemSize[1]) dy -= systemSize[1];

                                    if(dz < -0.5*systemSize[2]) dz += systemSize[2];
                                    else if(dz > 0.5*systemSize[2]) dz -= systemSize[2];

                                    double dr2 = dx*dx + dy*dy + dz*dz;
                                    if(dr2 < radiusSquared) {
                                        atom1->neighbors.push_back(atom2);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
