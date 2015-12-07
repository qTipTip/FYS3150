#ifndef VELOCITYVERLET_H
#define VELOCITYVERLET_H
#include "integrators/integrator.h"

class VelocityVerlet : public Integrator
{
public:
    VelocityVerlet();
    virtual void integrate(System *system, double dt) override;
private:
    bool m_first_step;
};
#endif
