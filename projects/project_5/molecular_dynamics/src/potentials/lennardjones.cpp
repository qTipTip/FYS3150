#include "lennardjones.h"
#include <cmath>
LennardJones::LennardJones(double sigma, double epsilon) :
    m_sigma(sigma),
    m_epsilon(epsilon)
{

}

void LennardJones::calculateForces(System *system) {
    if(system->systemSize().x() < 3*(m_rCut+m_rShell)) {
        calculateForcesN2(system);
    } else {
        calculateForcesNeighborlist(system);
    }
}

void LennardJones::calculateForcesNeighborlist(System *system)
{
    m_potentialEnergy = 0; // Remember to compute this in the loop
    double size_x = system->systemSize().x();
    double size_y = system->systemSize().y();
    double size_z = system->systemSize().z();
    double sigma_six = m_sigma*m_sigma*m_sigma*m_sigma*m_sigma*m_sigma;
    double potentialEnergyAtRcut = 4 * m_epsilon *pow(m_rCut,-6)*sigma_six*(sigma_six*pow(m_rCut,-6)  - 1.0);

    vector<Atom*> &atoms = system->atoms();
    // we sum over all distinct pairs of atoms
    for(Atom *atom_i : atoms) {
        for(Atom *atom_j : atom_i->neighbors) {
            double xij = atom_i->position[0] - atom_j->position[0];
            double yij = atom_i->position[1] - atom_j->position[1];
            double zij = atom_i->position[2] - atom_j->position[2];

            //since our system is periodic, we have to account for the two atoms being closer through the boundaries
            //than through the body of the lattice.

            if(xij > 0.5*size_x) xij -= size_x;
            else if(xij < -0.5*size_x) xij += size_x;

            if(yij > 0.5*size_y) yij -= size_y;
            else if(yij < -0.5*size_y) yij += size_y;

            if(zij > 0.5*size_z) zij -= size_z;
            else if(zij < -0.5*size_z) zij += size_z;

            // to make the program statements less daunting we precompute some values used in the expression for F
            double dr2 = xij*xij + yij*yij + zij*zij;
            if(dr2 >= m_rCut*m_rCut) continue;

            double oneOverDr2 = 1.0/dr2;
            double oneOverDr6 = oneOverDr2*oneOverDr2*oneOverDr2;

            double F =  24*m_epsilon*sigma_six*oneOverDr6*( 2*sigma_six*oneOverDr6 - 1.0)*oneOverDr2;
            double Fx = F*xij;
            double Fy = F*yij;
            double Fz = F*zij;

            atom_i->force[0] += Fx;
            atom_i->force[1] += Fy;
            atom_i->force[2] += Fz;
            atom_j->force[0] -= Fx;
            atom_j->force[1] -= Fy;
            atom_j->force[2] -= Fz;

            m_potentialEnergy += 4 * m_epsilon *oneOverDr6*sigma_six*(sigma_six*oneOverDr6  - 1.0) - potentialEnergyAtRcut;

        }
    }
}

void LennardJones::calculateForcesN2(System *system)
{
    m_potentialEnergy = 0; // Remember to compute this in the loop
    double size_x = system->systemSize().x();
    double size_y = system->systemSize().y();
    double size_z = system->systemSize().z();
    double sigma_six = m_sigma*m_sigma*m_sigma*m_sigma*m_sigma*m_sigma;

    vector<Atom*> &atoms = system->atoms();
    // we sum over all distinct pairs of atoms
    for(int i = 0; i < atoms.size(); i++) {
        Atom *atom_i = atoms[i]; // .at() ) supertregt. Bruk []
        for(int j = i+1; j < atoms.size(); j++) {
            Atom *atom_j = atoms[j];
            // vec3 r_ij = (atom_i->position - atom_j->position); // componentwise mye raskere

            // check distance
            //if(r_ij.length() < distance_correction) { continue; }

            double xij = atom_i->position[0] - atom_j->position[0];
            double yij = atom_i->position[1] - atom_j->position[1];
            double zij = atom_i->position[2] - atom_j->position[2];

            //since our system is periodic, we have to account for the two atoms being closer through the boundaries
            //than through the body of the lattice.

            if(xij > 0.5*size_x) xij -= size_x;
            else if(xij < -0.5*size_x) xij += size_x;

            if(yij > 0.5*size_y) yij -= size_y;
            else if(yij < -0.5*size_y) yij += size_y;

            if(zij > 0.5*size_z) zij -= size_z;
            else if(zij < -0.5*size_z) zij += size_z;

            // to make the program statements less daunting we precompute some values used in the expression for F
            double dr2 = xij*xij + yij*yij + zij*zij;
            double oneOverDr2 = 1.0/dr2;
            double oneOverDr6 = oneOverDr2*oneOverDr2*oneOverDr2;

            double F =  24*m_epsilon*sigma_six*oneOverDr6*( 2*sigma_six*oneOverDr6 - 1.0)*oneOverDr2;
            double Fx = F*xij;
            double Fy = F*yij;
            double Fz = F*zij;

            atom_i->force[0] += Fx;
            atom_i->force[1] += Fy;
            atom_i->force[2] += Fz;
            atom_j->force[0] -= Fx;
            atom_j->force[1] -= Fy;
            atom_j->force[2] -= Fz;

            m_potentialEnergy += 4 * m_epsilon *oneOverDr6*sigma_six*(sigma_six*oneOverDr6  - 1.0);

        }
    }
}
