#include "system.h"
#include "integrators/integrator.h"
#include "potentials/potential.h"
#include "statisticssampler.h"
#include "unitconverter.h"
#include "math/random.h"
#include <cmath>

System::System()
{

}

System::~System()
{
    delete m_potential;
    delete m_integrator;
    for(Atom *atom : m_atoms) {
        delete atom;
    }
    m_atoms.clear();
}

void System::applyPeriodicBoundaryConditions() {
    // Read here: http://en.wikipedia.org/wiki/Periodic_boundary_conditions#Practical_implementation:_continuity_and_the_minimum_image_convention
    // This method also takes into consideration the computation of the diffusion constant, hence
    // we in addition to moving the position, also move the initial position
    // to keep track of the true position of the atom
    for(Atom *atom : m_atoms) {

        if(atom->position[0] > m_systemSize[0]) {
            atom->position[0] -= m_systemSize[0];
            atom->initial_position[0] -= m_systemSize[0];
            atom->positionLastNeighborlistBuild[0] -= m_systemSize[0];
        }
        else if(atom->position[0] < 0) {
            atom->position[0] += m_systemSize[0];
            atom->initial_position[0] += m_systemSize[0];
            atom->positionLastNeighborlistBuild[0] += m_systemSize[0];
        }

        if(atom->position[1] > m_systemSize[1]) {
            atom->position[1] -= m_systemSize[1];
            atom->initial_position[1] -= m_systemSize[1];
            atom->positionLastNeighborlistBuild[1] -= m_systemSize[1];
        }
        else if(atom->position[1] < 0) {
            atom->position[1] += m_systemSize[1];
            atom->initial_position[1] += m_systemSize[1];
            atom->positionLastNeighborlistBuild[1] += m_systemSize[1];
        }
        if(atom->position[2] > m_systemSize[2]) {
            atom->position[2] -= m_systemSize[2];
            atom->initial_position[2] -= m_systemSize[2];
            atom->positionLastNeighborlistBuild[2] -= m_systemSize[2];
        }
        else if(atom->position[2] < 0) {
            atom->position[2] += m_systemSize[2];
            atom->initial_position[2] += m_systemSize[2];
            atom->positionLastNeighborlistBuild[2] += m_systemSize[2];
        }
    }
}

void System::removeTotalMomentum() {
    // Find the total momentum and remove momentum equally on each atom so the total momentum becomes zero.
    vec3 CM_V = vec3(0, 0, 0);
    double total_mass = 0;
    // Compute initial net momentum
    for(Atom *atom : m_atoms) {
        total_mass += atom->mass();
        CM_V += atom->velocity*atom->mass();
    }
    // subtract velocity from each atom
    CM_V /= total_mass; // CM_V now corresponds to the velocity to be subtracted from each atom
    for(Atom *atom : m_atoms) {
        atom->velocity -= CM_V;
    }
}

void System::resetForcesOnAllAtoms() {
    for(Atom *atom : m_atoms) {
        atom->resetForce();
    }
}

void System::createFCCLattice(int numberOfUnitCellsEachDimension, double latticeConstant, double temperature) {
    // You should implement this function properly. Right now, 100 atoms are created uniformly placed in the system of size (10, 10, 10).

    // For each unit cell, we create 4 new atoms, compute the global coordinate R of the unit cell, and then compute the local coordinates
    // r of each atom. The global position of the atom is then the R + r.
    for(int i = 0; i < numberOfUnitCellsEachDimension; i++){
        for(int j = 0; j < numberOfUnitCellsEachDimension; j++){
            for(int k = 0; k < numberOfUnitCellsEachDimension; k++){
                Atom *atom1 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                Atom *atom2 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                Atom *atom3 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                Atom *atom4 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                vec3 R = vec3(i*latticeConstant, j*latticeConstant, k*latticeConstant);
                vec3 r1 = R + vec3(0, 0, 0);
                vec3 r2 = R + vec3(latticeConstant/2, latticeConstant/2, 0);
                vec3 r3 = R + vec3(0, latticeConstant/2 , latticeConstant/2);
                vec3 r4 = R + vec3(latticeConstant/2, 0, latticeConstant/2);
                atom1->position.set(r1.x(), r1.y(), r1.z());
                atom2->position.set(r2.x() ,r2.y() ,r2.z());
                atom3->position.set(r3.x(), r3.y(), r3.z());
                atom4->position.set(r4.x(), r4.y(), r4.z());
                atom1->initial_position.set(r1.x(), r1.y(), r1.z());
                atom2->initial_position.set(r2.x() ,r2.y() ,r2.z());
                atom3->initial_position.set(r3.x(), r3.y(), r3.z());
                atom4->initial_position.set(r4.x(), r4.y(), r4.z());
                atom1->resetVelocityMaxwellian(temperature);
                atom2->resetVelocityMaxwellian(temperature);
                atom3->resetVelocityMaxwellian(temperature);
                atom4->resetVelocityMaxwellian(temperature);
                m_atoms.push_back(atom1);
                m_atoms.push_back(atom2);
                m_atoms.push_back(atom3);
                m_atoms.push_back(atom4);
            }
        }
    }
    double size = numberOfUnitCellsEachDimension*latticeConstant;
    setSystemSize(vec3(size, size, size));
}

void System::calculateForces() {
    resetForcesOnAllAtoms();
    m_potential->calculateForces(this);
}

void System::validateNeighborList() {
    for(Atom *atom : m_atoms) {
        if((atom->position - atom->positionLastNeighborlistBuild).length() > 0.5*m_potential->rShell()) {
            // std::cout << "Building neighbor list..." << std::endl;
            m_potential->buildNeighborList(*this);
            break;
        }
    }
}

void System::step(double dt) {
    if(m_firstStep) {
        m_potential->buildNeighborList(*this);
        m_firstStep = false;
    }

    validateNeighborList();
    m_integrator->integrate(this, dt);
    m_steps++;
    m_time += dt;
}
