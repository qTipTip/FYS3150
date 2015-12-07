#include "velocityverlet.h"
#include "../system.h"
#include "../atom.h"
#define FAST
VelocityVerlet::VelocityVerlet() : m_first_step(true) // sets default value to true at construction
{
}

void VelocityVerlet::integrate(System *system, double dt)
{
    // we first calculate the current forces
    if(m_first_step) {
        system->calculateForces();
        m_first_step = false;
    }
    double dt_half = dt / 2;
    // first we iterate over each atom in the system, and update the positions
    double oneOverMass = 1.0/system->atoms()[0]->mass();
    for(Atom *atom : system->atoms()) {
        atom->velocity[0] += atom->force[0]*dt_half*oneOverMass;
        atom->velocity[1] += atom->force[1]*dt_half*oneOverMass;
        atom->velocity[2] += atom->force[2]*dt_half*oneOverMass;

        atom->position[0] += atom->velocity[0] * dt;
        atom->position[1] += atom->velocity[1] * dt;
        atom->position[2] += atom->velocity[2] * dt;
    }

    // i spent a few good hours wondering why my simulations werent working. This line was missing
    system->applyPeriodicBoundaryConditions();

    // we then recompute the forces on all the atoms in preparation for computing the velocities
    system->calculateForces();

    // then we compute the velocities
    for(Atom *atom : system->atoms()) {
        atom->velocity[0] += atom->force[0]*dt_half*oneOverMass;
        atom->velocity[1] += atom->force[1]*dt_half*oneOverMass;
        atom->velocity[2] += atom->force[2]*dt_half*oneOverMass;
    }


}
